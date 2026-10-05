"""Endpoint elimination and boundary cross-term shift, original arithmetic."""
import numpy as np
from .physics import trapw, _dot
from .gsvd_backend import compact_gsvd

def factors_one_window(K, g, m2, H, pab, boundary, method: str):
    """Prepare weighted endpoint-shifted data for compact GSVD replay.

    Pole columns are unpenalized. Only interior continuum hats enter the H1
    penalty, while fixed threshold and Lambda values are moved to the RHS.
    """
    n_s = H.shape[0]
    n_full = K.shape[1]
    # 将 H 划分为内部块和端点块，用于处理边界平移。
    internal = np.arange(1, n_s - 1)
    n_internal = n_s - 2
    b_ii = H[np.ix_(internal, internal)]
    b_if = H[np.ix_(internal, [0, n_s - 1])]
    r_b = np.linalg.cholesky(b_ii).T
    sqrt_w = np.sqrt(trapw(m2))
    # 系数向量先放极点块；最后一个连续谱节点是固定上端点，
    # 因而不属于待求的自由向量。
    fixed = np.array([pab, n_full - 1])
    free = np.r_[np.arange(pab), np.arange(pab + 1, n_full - 1)]
    k0 = sqrt_w[:, None] * K[:, free]
    l0 = np.column_stack((np.zeros((n_internal, pab)), r_b))
    # 经济型 QR 后进行 GSVD，得到 Tikhonov 模态及其滤波因子。
    decomp = compact_gsvd(k0, l0, method=method, qr_reduce=True)
    if decomp.p != n_internal or decomp.k != pab:
        raise RuntimeError(
            f"GSVD ranks p={decomp.p} k={decomp.k}, expected p={n_internal} k={pab}"
        )
    q = np.linalg.solve(b_ii, _dot(b_if, boundary))
    b0 = sqrt_w * (g - _dot(K[:, fixed], boundary)) + _dot(k0[:, pab:], q)
    # 将平移后的数据投影到 GSVD 的左奇异向量上。
    beta_all = _dot(decomp.U.T, b0)
    beta_c = beta_all[:n_internal]
    beta_n = beta_all[n_internal:]
    recon = _dot(decomp.U, beta_all)
    null = float(_dot(b0 - recon, b0 - recon))
    if np.any(decomp.mu <= np.finfo(float).tiny):
        raise RuntimeError("GSVD mu has non-positive compact modes")
    sigma_tilde = decomp.sigma / decomp.mu
    data_norm = float(np.sqrt(np.sum(trapw(m2) * g * g)))
    return {
        "sigma": decomp.sigma,
        "mu": decomp.mu,
        "sigma_unreg": decomp.sigma_unregularized,
        "sigma_tilde": sigma_tilde,
        "beta_c": beta_c,
        "beta_n": beta_n,
        "null": null,
        "X": decomp.X,
        "q": q,
        "data_norm": data_norm,
        "n_data": np.int64(len(g)),
        "source": decomp.source,
        "library": decomp.library,
    }
