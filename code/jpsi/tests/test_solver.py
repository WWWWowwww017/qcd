"""使用 python -B -m unittest discover -s code/jpsi/tests -v 运行。"""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import warnings

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent.parent))
from jpsi import solve
from jpsi import physics as P
from jpsi import raus2024 as R
from jpsi.solver import borel_grid


class SolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        warnings.simplefilter("error", RuntimeWarning)

    def test_historical_nine_windows(self):
        # 只冻结历史网格，用归档的 J/psi 夹具检验当前物理、GSVD 和 alpha 选择。
        with np.load(ROOT / "reference_latest.npz", allow_pickle=False) as ref:
            cases = json.loads(str(ref["cases_json"]))
            for case in cases:
                with self.subTest(**case):
                    # 只冻结离散化，不冻结物理、GSVD 或 alpha 选择代码。
                    legacy_hat = P.sgrid_hats(case["Lambda"], nodes=200)
                    prefix = str(case["id"]) + "_"
                    with patch.object(P, "sgrid_hats", return_value=legacy_hat), patch.object(P, "sgrid_fine", return_value=ref[prefix + "s"]):
                        result = solve(case["M2_min"], case["M2_max"], case["Lambda"], channel=case["channel"])
                    self.assertEqual(result["metadata"]["n_hats"], 200)
                    out = result["channels"][case["channel"]]
                    self.assertAlmostEqual(np.log10(out["alpha_star"]), np.log10(float(ref[prefix + "alpha_star"])), places=10)
                    for key, source in [("residues", "residues"), ("hat", "continuum_hat"), ("rho", "spectra"), ("weighted_residual", "weighted_residual"), ("penalty", "penalty"), ("f_jpsi", "f_jpsi")]:
                        expected = ref[prefix + source]
                        if key == "residues":
                            expected = expected[:len(out["residues"])]
                        np.testing.assert_allclose(out[key], expected, rtol=2e-7, atol=2e-10, err_msg=key)
                    self.assertLess(out["residual_factor_discrepancy"], 2e-12)
                    np.testing.assert_array_equal(out["s"], ref[prefix + "s"])
                    self.assertEqual(out["hat"][0], 0.)
                    self.assertEqual(out["hat"][-1], P.aa_rho_pert(np.array([case["Lambda"]]))[0])

    def test_ope_reference_and_independent_integrals(self):
        # 将向量化 OPE 数值与保存的、独立积分得到的参考分量进行比较。
        records = json.loads((ROOT / "reference_ope.json").read_text())
        for record in records:
            with self.subTest(Lambda=record["Lambda"], M2=record["M2"]):
                values = np.array(P.aa_ope(record["Lambda"], np.array([record["M2"]]), components=True)).ravel()
                np.testing.assert_allclose(values, record["components"], rtol=2e-12, atol=1e-17)
                np.testing.assert_allclose(values, record["independent"], rtol=2e-4, atol=1e-17)

    def test_hat_partition_and_h1_integral(self):
        # 检验帽函数核的闭合恒等式与有限元 H1 二次型。
        for lam in (20., 30., 70.):
            t = np.array([1.5, 15., 80.])
            sums = P.hats(t, P.sgrid_hats(lam)).sum(axis=1)
            np.testing.assert_allclose(sums, np.exp(-P.S_H/t)-np.exp(-lam/t), rtol=2e-9, atol=1e-10)
        s = np.array([0., .3, 1.1, 2.])
        c = np.array([.2, -.4, .7, 1.3])
        self.assertAlmostEqual(float(np.dot(c, np.dot(P.h_mass_stiffness(s), c))), 4.150166666666667, places=13)

    def test_arbitrary_upper_endpoint_and_stationarity(self):
        # 检验精确端点、短谱单元，以及单双极点模型的正规方程驻点条件。
        result = solve(5.003, 20.007, 30.123)
        t = result["M2"]
        self.assertEqual(t[0], 5.003)
        self.assertEqual(t[-1], 20.007)
        self.assertLessEqual(float(np.max(np.diff(t))), .01000000000001)
        np.testing.assert_allclose(result["weights"].sum(), 20.007-5.003, atol=1e-13)
        self.assertEqual(result["metadata"]["n_hats"], 164)
        self.assertEqual(result["metadata"]["ds"], .1)
        self.assertEqual(result["metadata"]["dM2"], .01)
        for out in result["channels"].values():
            np.testing.assert_array_equal(out["s"], out["s_hat"])
            np.testing.assert_array_equal(out["rho"], out["hat"])
            self.assertEqual(out["s_hat"][0], P.S_H)
            self.assertEqual(out["s_hat"][-1], 30.123)
            np.testing.assert_allclose(np.diff(out["s_hat"])[:-1], .1, rtol=0, atol=1e-14)
            self.assertAlmostEqual(np.diff(out["s_hat"])[-1], .0124846976, places=13)
            self.assertEqual(out["factors"]["X"].shape, (162+len(out["masses"]),)*2)
            poles = np.column_stack([np.exp(-m*m/t)/t for m in out["masses"]])
            K = np.column_stack((poles, P.hats(t, out["s_hat"])))
            weighted_res = result["weights"] * out["residual"]
            grad = np.dot(K.T, weighted_res)
            p = len(out["masses"])
            np.testing.assert_allclose(grad[:p], 0., atol=1e-12)
            penalty_grad = np.dot(P.h_mass_stiffness(out["s_hat"]), out["hat"])
            np.testing.assert_allclose(grad[p+1:-1]+out["alpha_star"]*penalty_grad[1:-1], 0., atol=2e-12)
            np.testing.assert_allclose(out["backsub"], out["backsub_poles"]+out["backsub_continuum"], atol=1e-15)
            self.assertLess(out["residual_factor_discrepancy"], 2e-12)

    def test_default_s_grid(self):
        # 验证锚定的 Delta-s 网格，并确认会拒绝过短区间。
        for lam in (20., 30., 30.123, 70., P.S_H+.2, P.S_H+.200001,
                    np.nextafter(P.S_H+.2, -np.inf), P.S_H+.100001):
            with self.subTest(Lambda=lam):
                s = P.sgrid_hats(lam)
                self.assertEqual(s[0], P.S_H)
                self.assertEqual(s[-1], lam)
                self.assertGreaterEqual(len(s), 3)
                self.assertTrue(np.all(np.diff(s) > 0))
                np.testing.assert_allclose(np.diff(s)[:-1], .1, rtol=0, atol=2e-14)
                self.assertLessEqual(s[-1]-s[-2], .1+2e-14)
                np.testing.assert_array_equal(s[:-1], P.S_H+.1*np.arange(len(s)-1))
                np.testing.assert_array_equal(P.sgrid_fine(lam), s)
        self.assertEqual(len(P.sgrid_hats(30.)), 162)
        for lam in (P.S_H, P.S_H+.05, P.S_H+.1, np.nan, np.inf):
            with self.subTest(Lambda=lam), self.assertRaises(ValueError):
                P.sgrid_hats(lam)
        with self.assertRaisesRegex(ValueError, "at least 3 hat nodes"):
            solve(5., 80., P.S_H+.05)
        with patch.object(P, "sgrid_hats", return_value=np.array([P.S_H, 30.])):
            with self.assertRaisesRegex(ValueError, "at least 3 hat nodes"):
                solve(5., 80., 30.)

    def test_actual_free_dimension_validation(self):
        # 只有验证通过后才停止在 OPE 计算之前。
        for lam, nodes in ((30., 162), (30.123, 164), (70., 562)):
            for channel, poles in (("1delta", 1), ("2delta", 2)):
                free = nodes-2+poles
                with self.subTest(Lambda=lam, channel=channel):
                    short_hi = 1.5+.01*(free-2)
                    with self.assertRaisesRegex(ValueError, f"at least {free} Borel points.*got {free-1}"):
                        solve(1.5, short_hi, lam, channel=channel)
                    enough_hi = 1.5+.01*(free-1)
                    with patch("jpsi.solver.lapacke_info", side_effect=RuntimeError("validation passed")):
                        with self.assertRaisesRegex(RuntimeError, "validation passed"):
                            solve(1.5, enough_hi, lam, channel=channel)

    def test_bad_parameters(self):
        # 非法物理窗口和未知极点模型应在数值计算前立即失败。
        for values in [(0,80,30), (-1,80,30), (5,5,30), (80,5,30), (np.nan,80,30), (5,np.inf,30), (5,80,np.nan), (5,80,P.S_H), (True,80,30), ("5",80,30), (5,5.01,30)]:
            with self.subTest(values=values), self.assertRaises(ValueError):
                solve(*values)
        with self.assertRaises(ValueError):
            solve(5,80,30,channel="3delta")

    def test_grid_endpoints_and_old_arithmetic(self):
        # 对齐端点时保留历史的 0.01 Borel 步长算法。
        old = 1.5 + .01*np.arange(7851)
        np.testing.assert_array_equal(borel_grid(5.,80.), old[350:])
        for lo, hi in [(1.49, 4.), (1.501, 5.123), (81., 90.005), (5., 5.0001)]:
            t = borel_grid(lo,hi)
            self.assertEqual(t[0],lo)
            self.assertEqual(t[-1],hi)
            self.assertTrue(np.all(np.diff(t)>0))

    def test_raus_orientation_and_fallback(self):
        # 选出的 alpha 不依赖扫描方向；不存在 Q 局部极小时应暴露 HR 兜底规则。
        alpha = np.logspace(-8,0,801)
        sigma = np.array([1., .03, .00001])
        beta = np.array([1., .1, .001])
        up = R.select(sigma,beta,1e-10,alpha)
        down = R.select(sigma,beta,1e-10,alpha[::-1])
        self.assertEqual(up["alpha"],down["alpha"])
        self.assertEqual(up["index"],len(alpha)-1-down["index"])
        fallback = R.raus2024_rule(np.array([.01,.1,1.]),np.array([1.,2.,3.]),np.array([.01,.2,3.]),.001)
        self.assertEqual(fallback["fallback"],"psiHR_global")


if __name__ == "__main__":
    unittest.main(verbosity=2)
