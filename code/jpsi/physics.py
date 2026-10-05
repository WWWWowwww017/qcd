"""用于紧凑 GSVD 运行的 J/psi n=0 冻结 OPE 与帽函数核。

Correct n=0 LO*(1+alpha/pi) + D4 + D6. Fixed mu=2, three-loop RG.
The inversion uses 0.1 GeV^2 hat spacing on [s_h, Lambda]; arbitrary Borel windows.
"""
from __future__ import annotations

import numpy as np

MC = 1.0990536546335794
ALPHA_S = 0.3029450988131294
MU = 2.0
G2 = 0.038
G3 = 0.013
S_H = 13.9105153024
DS = 0.1
M2STEP = 0.01
# 极点质量按模型保存；极点留数是不受惩罚的 delta 极点幅度，
# 重构的连续谱位于物理 s 网格上。
MASSES = {
    "1delta": (3.096900,),
    "2delta": (3.096900, 3.686097),
}


def _dot(a, b):
    """1D/2D products without Accelerate matmul's spurious FP warnings.

    dot does not report real overflow either, so reject nonfinite results.
    """
    out = np.dot(a, b)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("matrix product returned NaN or Inf")
    return out


def trapw(x: np.ndarray) -> np.ndarray:
    """返回真实梯形积分权重，包括最后一个短单元。"""
    d = np.diff(x)
    w = np.empty_like(x)
    w[0] = d[0] / 2
    w[-1] = d[-1] / 2
    w[1:-1] = (d[:-1] + d[1:]) / 2
    return w


def aa_rho_pert(s):
    """带固定 NLO 乘法修正的 J/psi 微扰连续谱。"""
    s = np.asarray(s, float)
    out = np.zeros_like(s)
    q = 4 * MC * MC
    m = s > q
    v = np.sqrt(1 - q / s[m])
    out[m] = 3 * v * (1 - v * v / 3) / (8 * np.pi**2) * (1 + ALPHA_S / np.pi)
    return out


def aa_ope(L, m2, nv=16001, nx=4001, components=False):
    """计算 Borel 变换后的冻结 n=0 OPE。

    The perturbative contribution is integrated in the velocity variable;
    D4 and D6 are integrated over the Feynman parameter x. The three terms
    remain separate so the forward model is auditable in the web workbench.
    """
    m2 = np.asarray(m2, float)
    q = 4 * MC * MC
    vel = np.linspace(0, np.sqrt(max(1.0 - q / L, 0.0)), nv)
    denom = np.maximum(1.0 - vel * vel, q / L)
    sg = q / denom
    sg[-1] = L
    ws = trapw(sg)
    rw = aa_rho_pert(sg)
    # 分块处理 Borel 点，控制指数临时数组的大小。
    pert = np.empty_like(m2)
    for a in range(0, len(m2), 256):
        mm = m2[a : a + 256]
        pert[a : a + len(mm)] = _dot(np.exp(-sg[None, :] / mm[:, None]) / mm[:, None] * rw, ws)
    x = np.linspace(1e-8, 1 - 1e-8, nx)
    xb = 1 - x
    wx = trapw(x)
    a3 = -(45 / 8) * x * xb * (x**3 + xb**3) - (69 / 72) * x**2 * xb**2
    b3 = (
        -(3 / 4) * MC**2 * x * xb * (x**4 + xb**4)
        + MC**2 * x**2 * xb**2 * (-(23 / 12) * (x**2 + xb**2) - 11 * x * xb / 3 + 2)
    )
    c3 = (2 / 5) * MC**4 * x * xb * (x**5 + xb**5)
    # 这些系数函数来自固定的源公式，不是逆问题求解器引入的拟合参数。
    g2 = np.empty_like(m2)
    g3 = np.empty_like(m2)
    for a in range(0, len(m2), 128):
        mm = m2[a : a + 128, None]
        E = np.exp(-MC**2 / (x * xb * mm))
        i2 = E * (1 / (2 * mm**2) - MC**2 * (x**3 + xb**3) / (2 * mm**3 * x**2 * xb**2))
        i3 = E * (
            a3 / (2 * mm**3 * x**3 * xb**3)
            + b3 / (6 * mm**4 * x**4 * xb**4)
            + c3 / (24 * mm**5 * x**5 * xb**5)
        )
        g2[a : a + len(mm)] = G2 * _dot(i2, wx) / (6 * np.pi)
        g3[a : a + len(mm)] = G3 * _dot(i3, wx) / (4 * np.pi) ** 2
    if components:
        return pert, g2, g3
    return pert + g2 + g3


def sgrid_fine(L: float) -> np.ndarray:
    """输出网格与默认反演网格相同，不进行加密。"""
    return sgrid_hats(L)


def sgrid_hats(L: float, nodes: int | None = None) -> np.ndarray:
    """从 S_H 锚定并以 DS 步进，保留精确 Lambda 和最后的短区间。

    Explicit nodes selects the historical linspace for regression only.
    Near-coincident endpoints are merged at float64 roundoff, not at DS scale.
    """
    if not np.isfinite(L) or L <= S_H:
        raise ValueError("Lambda must be finite and exceed S_H")
    if nodes is not None:
        if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)) or nodes < 3:
            raise ValueError("need at least 3 hat nodes (integer nodes)")
        return np.linspace(S_H, L, nodes)
    # 从 S_H 锚定网格，并精确保留 Lambda 端点。
    n = int(np.floor((L - S_H) / DS))
    s = S_H + DS * np.arange(n + 1)
    s = s[s < L]
    tol = 8 * np.finfo(float).eps * max(1., abs(S_H), abs(L))
    if len(s) > 1 and L - s[-1] <= tol:
        s[-1] = L
    else:
        s = np.r_[s, L]
    if len(s) < 3:
        raise ValueError("need at least 3 hat nodes on the 0.1 GeV^2 grid; increase Lambda")
    if not np.all(np.diff(s) > 0):
        raise ValueError("hat grid must have strictly increasing float64 values")
    return s


def hats(m2, s):
    """组装分段线性帽函数的解析 Borel 变换。"""
    K = np.zeros((len(m2), len(s)))
    for j, h in enumerate(np.diff(s)):
        q = h / m2
        st = np.exp(-s[j] / m2)
        small = np.abs(q) < 1e-4
        asc = np.empty_like(q)
        des = np.empty_like(q)
        z = q[small]
        asc[small] = 0.5 - z / 3 + z * z / 8 - z**3 / 30 + z**4 / 144
        des[small] = 0.5 - z / 6 + z * z / 24 - z**3 / 120 + z**4 / 720
        z = q[~small]
        e = np.exp(-z)
        asc[~small] = (1 - (1 + z) * e) / z**2
        des[~small] = (z - 1 + e) / z**2
        fac = h / m2
        K[:, j] += st * fac * des
        K[:, j + 1] += st * fac * asc
    return K


def h_mass_stiffness(s: np.ndarray) -> np.ndarray:
    """在物理 s 网格的帽函数基上构造 H = M + M1。"""
    # H 表示有限元质量项与一阶导数惩罚项之和。
    h = np.diff(s)
    n = len(s)
    H = np.zeros((n, n))
    for i, hi in enumerate(h):
        H[i, i] += hi / 3 + 1 / hi
        H[i + 1, i + 1] += hi / 3 + 1 / hi
        H[i, i + 1] += hi / 6 - 1 / hi
        H[i + 1, i] += hi / 6 - 1 / hi
    return H

