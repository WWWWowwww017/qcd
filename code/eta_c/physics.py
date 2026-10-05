"""用于紧凑 GSVD 反演的 eta_c n=0 冻结 OPE 与帽函数核。

A-A n=0 的 LO + D4 + D6，不含 alpha_s 项；mu=2，质量采用领先对数约定。
物理系数来自 2026-09-11 的 eta_aa_physics.py 快照。
反演在 [s_h, Lambda] 上使用 0.1 GeV^2 的帽节点间隔，并允许任意 Borel 窗口。
"""
from __future__ import annotations

import numpy as np

MC = 1.1334710909
MU = 2.0
G2 = 0.038
G3 = 0.013
S_H = 14.99
DS = 0.1
M2STEP = 0.01
# 极点质量按模型保存；极点留数是不受惩罚的 delta 极点幅度，
# 连续谱在 s 网格上重构。
MASSES = {
    "1delta": (2.9841,),
    "2delta": (2.9841, 3.51067),
}


def _dot(a, b):
    """用 dot 完成一维或二维乘法，避免底层矩阵乘法的虚假浮点警告。

    dot 本身也不会自动报告真实溢出，因此这里显式拒绝非有限结果。
    """
    out = np.dot(a, b)
    if not np.all(np.isfinite(out)):
        raise FloatingPointError("matrix product returned NaN or Inf")
    return out


def trapw(x: np.ndarray) -> np.ndarray:
    """返回真实梯形积分权重，并保留最后一个短单元。"""
    d = np.diff(x)
    w = np.empty_like(x)
    w[0] = d[0] / 2
    w[-1] = d[-1] / 2
    w[1:-1] = (d[:-1] + d[1:]) / 2
    return w


def aa_rho_pert(s):
    """计算重夸克阈值以上的领先阶微扰连续谱。"""
    s = np.asarray(s, float)
    out = np.zeros_like(s)
    q = 4 * MC * MC
    m = s > q
    v = np.sqrt(1 - q / s[m])
    out[m] = 3 * v * (1 - v * v / 3) / (8 * np.pi**2)
    return out


def aa_ope(L, m2, nv=16001, nx=4001, components=False):
    """计算 Borel 变换后的冻结 A-A n=0 OPE。

    微扰项使用速度变量积分，D4 和 D6 凝聚项使用 Feynman 参数 x 积分。
    保留三个分项，便于网页绘图和数值检验追溯到 sum rule 输入。
    """
    m2 = np.asarray(m2, float)
    q = 4 * MC * MC
    vel = np.linspace(0, np.sqrt(max(1.0 - q / L, 0.0)), nv)
    denom = np.maximum(1.0 - vel * vel, q / L)
    sg = q / denom
    sg[-1] = L
    ws = trapw(sg)
    rw = aa_rho_pert(sg)
    # 分块处理 Borel 网格，避免指数临时矩阵不必要地增长。
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
        + MC**2 * x**2 * xb**2 * (-(23 / 12) * (x**2 + xb**2) + x * xb / 3 + 2)
    )
    c3 = (2 / 5) * MC**4 * x * xb * (x**5 + xb**5)
    # 以下系数函数是物理快照中的固定 D4/D6 表达式，不引入拟合参数。
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
    """输出网格与默认反演网格相同，不在输出阶段加密。"""
    return sgrid_hats(L)


def sgrid_hats(L: float, nodes: int | None = None) -> np.ndarray:
    """从 S_H 锚定并以 DS 步进，保留精确 Lambda 和最后的短区间。

    显式传入 nodes 只用于历史 linspace 回归。
    只有 float64 舍入造成的近重合端点才会合并，不按 DS 尺度合并。
    """
    if not np.isfinite(L) or L <= S_H:
        raise ValueError("Lambda must be finite and exceed S_H")
    if nodes is not None:
        if isinstance(nodes, bool) or not isinstance(nodes, (int, np.integer)) or nodes < 3:
            raise ValueError("need at least 3 hat nodes (integer nodes)")
        return np.linspace(S_H, L, nodes)
    # 从物理阈值锚定，而不是拉伸到固定节点数，从而精确保留 Lambda。
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
    """组装分片线性帽函数的解析 Borel 变换。"""
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
    # 二次惩罚对应 integral[(rho)^2 + (d rho/ds)^2] ds。
    h = np.diff(s)
    n = len(s)
    H = np.zeros((n, n))
    for i, hi in enumerate(h):
        H[i, i] += hi / 3 + 1 / hi
        H[i + 1, i + 1] += hi / 3 + 1 / hi
        H[i, i + 1] += hi / 6 - 1 / hi
        H[i + 1, i] += hi / 6 - 1 / hi
    return H

