"""端点消元与边界交叉项平移，保持原始数值运算。"""
import numpy as np
from .physics import trapw, _dot
from .gsvd_backend import compact_gsvd

def factors_one_window(K, g, m2, H, pab, boundary, method: str):
    """为紧凑 GSVD 回放准备一个加权反问题。

    连续谱端点固定，极点列不受惩罚，只有内部连续谱系数进入 H1 惩罚。
    返回的投影量足以在之后计算任意正 alpha。
    """
    n_s = H.shape[0]
    n_full = K.shape[1]
    # 将 H 分成内部块和端点块，用于处理边界平移。
    internal = np.arange(1, n_s - 1)
    n_internal = n_s - 2
    b_ii = H[np.ix_(internal, internal)]
    b_if = H[np.ix_(internal, [0, n_s - 1])]
    r_b = np.linalg.cholesky(b_ii).T
    sqrt_w = np.sqrt(trapw(m2))
    # 完整系数向量的前 pab 项是极点；最后一个连续谱节点是固定的 Lambda 端点。
    fixed = np.array([pab, n_full - 1])
    free = np.r_[np.arange(pab), np.arange(pab + 1, n_full - 1)]
    k0 = sqrt_w[:, None] * K[:, free]
    l0 = np.column_stack((np.zeros((n_internal, pab)), r_b))
    # QR+GSVD 同时对角化加权数据项和 H1 惩罚项。
    decomp = compact_gsvd(k0, l0, method=method, qr_reduce=True)
    if decomp.p != n_internal or decomp.k != pab:
        raise RuntimeError(
            f"GSVD ranks p={decomp.p} k={decomp.k}, expected p={n_internal} k={pab}"
        )
    q = np.linalg.solve(b_ii, _dot(b_if, boundary))
    b0 = sqrt_w * (g - _dot(K[:, fixed], boundary)) + _dot(k0[:, pab:], q)
    # beta 投影是 Tikhonov 滤波器使用的数据坐标。
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
