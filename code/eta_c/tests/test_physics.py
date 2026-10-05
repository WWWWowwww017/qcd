"""ηc 专用的源代码/MATLAB 回归测试与无 pickle 输出测试。"""
from pathlib import Path
import json
import sys
import tempfile
import unittest

import numpy as np
from scipy.io import loadmat

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent.parent))
from eta_c import solve, lapacke_info
from eta_c import physics as P
from eta_c.__main__ import save_result


def replay(out, alpha):
    """在固定 alpha 下回放保存的 GSVD，且不依赖选参器。"""
    f = out["factors"]
    n = len(f["sigma"])
    z = P._dot(f["X"][:, :n], f["sigma"]*f["beta_c"] /
               (f["sigma"]**2+alpha*f["mu"]**2))
    z += P._dot(f["X"][:, n:], f["beta_n"]/f["sigma_unreg"])
    p = len(out["masses"])
    return z[:p], np.r_[out["boundary"][0], z[p:]-f["q"], out["boundary"][1]]


class PhysicsTests(unittest.TestCase):
    def test_eta_constants_and_lo(self):
        # 在测试逆问题求解前，先检查固定常数和微扰谱密度的阈值行为。
        self.assertEqual(P.MC, 1.1334710909)
        self.assertEqual(P.MU, 2.)
        self.assertEqual(P.S_H, 14.99)
        self.assertEqual(P.MASSES, {"1delta": (2.9841,), "2delta": (2.9841, 3.51067)})
        self.assertEqual((P.G2, P.G3), (.038, .013))
        s = np.array([0., 4*P.MC**2, 20., 40.])
        v = np.sqrt(1-4*P.MC**2/s[2:])
        np.testing.assert_array_equal(P.aa_rho_pert(s)[:2], [0., 0.])
        np.testing.assert_allclose(P.aa_rho_pert(s)[2:], 3*v*(1-v*v/3)/(8*np.pi**2), rtol=1e-15)

    def test_source_three_windows(self):
        # 在三组几何参数和多个 alpha 值下复现源夹具，
        # 包括病态的小 alpha 区域。
        ref = loadmat(ROOT.parent / "test_reference.mat", simplify_cells=True)
        for case in ref["cases"]:
            with self.subTest(geometry=case["geometry"]):
                r = solve(*case["geometry"])
                t = case["sample_m2"]
                s = P.sgrid_hats(case["geometry"][2])
                np.testing.assert_allclose(s, case["s"], atol=2e-13, rtol=0)
                np.testing.assert_allclose(P.aa_ope(case["geometry"][2], t), case["g"], rtol=3e-12)
                np.testing.assert_allclose(P.hats(t, s), case["Khat"], rtol=2e-8, atol=1e-13)
                np.testing.assert_allclose(P.h_mass_stiffness(s), case["H"], rtol=3e-12)
                for name, source in (("1delta", "one"), ("2delta", "two")):
                    for j, loga in enumerate(ref["logs"]):
                        residues, rho = replay(r["channels"][name], 10.**loga)
                        expected = np.r_[np.asarray(case[source]["residues"][j]).ravel(), case[source]["rho"][j]]
                        actual = np.r_[residues, rho]
                        error = np.linalg.norm(actual-expected)/np.linalg.norm(expected)
                        self.assertLess(error, 3e-4 if loga < -12 else 5e-8)

    def test_independent_matlab_same_grid(self):
        # 在完全相同的 Borel 与谱网格上，将 Python OPE、GSVD 回放
        # 与独立 MATLAB 参考实现比较。
        r = solve(5., 80., 30.)
        with np.load(ROOT / "reference_matlab.npz", allow_pickle=False) as ref:
            np.testing.assert_allclose(r["M2"], ref["M2"], atol=2e-14, rtol=0)
            np.testing.assert_allclose(r["weights"], ref["weights"], atol=2e-14, rtol=0)
            for key, source in (("pert", "pert"), ("D4", "g2"), ("D6", "g3"), ("total", "g")):
                np.testing.assert_allclose(r["OPE"][key], ref[source], rtol=3e-12, atol=1e-18)
            for name, source in (("1delta", "one"), ("2delta", "two")):
                out = r["channels"][name]
                np.testing.assert_allclose(out["s"], ref["s"], atol=2e-13, rtol=0)
                for j, alpha in enumerate(ref["alpha"]):
                    residues, rho = replay(out, alpha)
                    expected = np.r_[np.asarray(ref[source+"_residues"][j]).ravel(), ref[source+"_rho"][j]]
                    error = np.linalg.norm(np.r_[residues, rho]-expected)/np.linalg.norm(expected)
                    self.assertLess(error, 3e-4 if alpha < 1e-12 else 5e-8)
                index = 6 if name == "1delta" else 7
                self.assertAlmostEqual(np.log10(out["alpha_star"]), np.log10(ref["alpha"][index]), places=10)
                np.testing.assert_allclose(out["residues"], np.atleast_1d(ref[source+"_residues"][index]), rtol=5e-8)
                np.testing.assert_allclose(out["hat"], ref[source+"_rho"][index], rtol=5e-8, atol=1e-10)
                np.testing.assert_allclose(out["penalty"], ref[source+"_penalty"][index], rtol=5e-8)
                np.testing.assert_allclose(out["relative_weighted_residual"], ref[source+"_relative_residual"][index], rtol=5e-8)

    def test_output_and_real_backend(self):
        # 确认使用真实 LAPACKE 后端，并生成不含 pickle 的数值 NPZ 输出。
        self.assertTrue(lapacke_info()["available"])
        r = solve(5., 20., 20., channel="1delta")
        self.assertEqual(set(r["channels"]), {"1delta"})
        self.assertTrue(r["channels"]["1delta"]["factors"]["source"].startswith("LAPACKE_dggsvd3"))
        with tempfile.TemporaryDirectory() as tmp:
            a = save_result(r, tmp)
            b = save_result(r, tmp)
            self.assertNotEqual(a["result_directory"], b["result_directory"])
            dest = Path(a["result_directory"])
            self.assertEqual(json.loads((dest/"summary.json").read_text()), a)
            with np.load(dest/"result.npz", allow_pickle=False) as data:
                for key in data.files:
                    self.assertNotEqual(data[key].dtype, object)
                np.testing.assert_array_equal(data["channels__1delta__hat"], r["channels"]["1delta"]["hat"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
