"""独立的 J/psi n=0 谱反演；不读取缓存物理量或外部工作区。"""
from __future__ import annotations

import math
import numpy as np
from . import physics as P
from .alpha_grid import selection_grid
from .factors import factors_one_window
from .gsvd_backend import _load_dggsvd3
from . import raus2024 as R


def lapacke_info():
    """探测已安装动态库，不执行安装操作。"""
    loaded = _load_dggsvd3()
    return {"available": loaded is not None, "library": loaded[1] if loaded else None}


def _positive(value, name):
    if isinstance(value, (bool, str, bytes)) or np.ndim(value) != 0:
        raise ValueError(f"{name} must be a finite positive scalar")
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite positive scalar") from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive scalar")
    return value


def borel_grid(lo, hi):
    """包含精确输入端点，并保留历史锚定算术。

    当下界位于 1.5 + .01*k 网格上时沿用该表达式，否则从 lo 开始。
    允许最后一个短步；接近端点的样本会被替换，而不会重复添加。
    """
    step = P.M2STEP
    index = round((lo - 1.5) / step)
    anchor = 1.5 + step * index
    tol = 8 * np.finfo(float).eps * max(1., abs(lo), abs(hi))
    n = int(np.floor((hi - lo) / step))
    if abs(anchor - lo) <= tol:
        grid = 1.5 + step * (index + np.arange(n + 1))
    else:
        grid = lo + step * np.arange(n + 1)
    grid[0] = lo
    grid = grid[grid < hi]
    if len(grid) > 1 and hi - grid[-1] <= tol:
        grid[-1] = hi
    else:
        grid = np.r_[grid, hi]
    if not np.all(np.diff(grid) > 0):
        raise ValueError("Borel grid cannot be represented with strictly increasing float64 values")
    return grid


def solve(M2_min, M2_max, Lambda, *, channel=None):
    """返回 {metadata, M2, weights, OPE, channels}；channel=None 表示两个模型。

    所有有量纲输入的单位均为 GeV^2；函数不写入文件。
    channel 可为 '1delta' 或 '2delta'，每个通道独立选择 alpha。
    """
    lo = _positive(M2_min, "M2_min")
    hi = _positive(M2_max, "M2_max")
    lam = _positive(Lambda, "Lambda")
    if hi <= lo:
        raise ValueError("M2_max must exceed M2_min")
    if channel is not None and (not isinstance(channel, str) or channel not in P.MASSES):
        raise ValueError("channel must be None, '1delta', or '2delta'")
    channels = tuple(P.MASSES) if channel is None else (channel,)
    if lam <= P.S_H or any(lam <= m*m for c in channels for m in P.MASSES[c]):
        raise ValueError("Lambda must exceed S_H and all included pole masses squared")
    # 构造物理谱网格，并确认 Borel 窗口包含足够样本以分辨
    # 每个自由极点和连续谱系数。
    s_hat = P.sgrid_hats(lam)
    n_hats = len(s_hat)
    if n_hats < 3:
        raise ValueError("need at least 3 hat nodes")
    t = borel_grid(lo, hi)
    for c in channels:
        free = n_hats - 2 + len(P.MASSES[c])
        if len(t) < free:
            raise ValueError(f"{c} needs at least {free} Borel points (free unknowns); got {len(t)}")
    info = lapacke_info()
    if not info["available"]:
        raise RuntimeError("LAPACKE_dggsvd3 unavailable. Point JPSI_LAPACKE_LIBRARY to an existing LP64 LAPACKE library; no numerical fallback is used.")
    s = P.sgrid_fine(lam)
    # OPE 提供 g(M^2)；K 将极点与连续谱系数映射到同一数据量。
    pert, d4, d6 = P.aa_ope(lam, t, components=True)
    g = pert + d4 + d6
    hat_kernel = P.hats(t, s_hat)
    H = P.h_mass_stiffness(s_hat)
    boundary = np.array([0., float(P.aa_rho_pert(np.array([lam]))[0])])
    weights = P.trapw(t)
    if not all(np.all(np.isfinite(a)) for a in (g, hat_kernel, H)):
        raise RuntimeError("Nonfinite physical operator/OPE for this window")
    results = {}
    for c in channels:
        # 两种极点模型在共享物理网格上独立正则化，
        # 再依据 Raus 诊断量分别选择 alpha。
        masses = np.array(P.MASSES[c]); p = len(masses)
        poles = np.column_stack([np.exp(-m*m/t)/t for m in masses])
        K = np.column_stack((poles, hat_kernel))
        fac = factors_one_window(K, g, t, H, p, boundary, "lapacke")
        alpha, covered = selection_grid(fac["sigma_tilde"], fac["beta_c"], fac["null"])
        sel = R.select(fac["sigma_tilde"], fac["beta_c"], fac["null"], alpha)
        if sel["status"] != "defined":
            raise RuntimeError(f"{c}: alpha selection {sel['status']}")
        a = sel["alpha"]; sig = fac["sigma"]; mu = fac["mu"]; n = len(sig)
        # 对受惩罚模态使用标量 GSVD Tikhonov 滤波；末端未正则化模态
        # 承载极点幅度。
        coef = sig*fac["beta_c"]/(sig*sig+a*mu*mu)
        un = fac["beta_n"]/fac["sigma_unreg"]
        shifted = P._dot(fac["X"][:, :n], coef) + P._dot(fac["X"][:, n:], un)
        # 撤销端点平移；只为展示结果插值连续谱，逆问题求解仍使用帽函数节点向量。
        residues = shifted[:p]
        cont = np.r_[boundary[0], shifted[p:]-fac["q"], boundary[1]]
        rho = np.interp(s, s_hat, cont)
        backsub = P._dot(K, np.r_[residues, cont])
        residual = backsub-g
        wr = float(np.sqrt(np.sum(weights*residual*residual)))
        pen = float(np.sqrt(max(0., float(P._dot(cont, P._dot(H, cont))))))
        f = np.full(p, np.nan)
        positive = residues > 0
        f[positive] = np.sqrt(residues[positive])
        d_md, psi = R.quantities(fac["sigma_tilde"], fac["beta_c"], fac["null"], alpha)
        d_factor = float(np.sqrt(np.sum((a/(fac["sigma_tilde"]**2+a)*fac["beta_c"])**2)+fac["null"]))
        if not all(np.all(np.isfinite(v)) for v in (residues, cont, residual, rho)):
            raise RuntimeError(f"{c}: nonfinite reconstructed solution")
        results[c] = dict(
            masses=masses, s=s, rho=rho, s_hat=s_hat, hat=cont,
            residues=residues, f=f, f_jpsi=float(f[0]), boundary=boundary,
            alpha_star=a, alpha_diagnostics=dict(
                **sel, grid=alpha, d_MD=d_md, psi_Q=psi,
                psi_MQ=R.psi_mq_tikhonov(psi, d_md), psi_HR=d_md/np.sqrt(alpha),
                on_boundary=sel["index"] in (0,len(alpha)-1),
                small_alpha_peak_covered=covered, sigma_min_sq=float(np.min(fac["sigma_tilde"])**2)),
            backsub=backsub, backsub_poles=P._dot(poles, residues),
            backsub_continuum=P._dot(hat_kernel, cont), residual=residual,
            weighted_residual=wr, relative_weighted_residual=wr/max(fac["data_norm"], np.finfo(float).tiny),
            unweighted_residual=float(np.linalg.norm(residual)), penalty=pen,
            objective=wr**2+a*pen**2, residual_factor=d_factor,
            residual_factor_discrepancy=abs(d_factor-wr), min_rho=float(np.min(rho)),
            negative_rho_count=int(np.count_nonzero(rho<0)), factors=fac)
    return dict(metadata=dict(schema="jpsi-standalone-v1", M2_min=lo, M2_max=hi,
                Lambda=lam, n_hats=n_hats, ds=P.DS, dM2=P.M2STEP, mu=P.MU, mc=P.MC,
                alpha_s=P.ALPHA_S, G2=P.G2, G3=P.G3, S_H=P.S_H,
                rule="raus2024_etac_equiv_p0_1", lapacke=info),
                M2=t, weights=weights, OPE=dict(pert=pert, D4=d4, D6=d6, total=g), channels=results)
