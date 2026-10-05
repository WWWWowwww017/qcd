"""采用最新 ηc 等价的 Raus 2024 规则 r_H、p0=1，而非历史 J/psi 规则。"""
from __future__ import annotations

import numpy as np

TINY = np.finfo(np.float64).tiny
P0_TIKHONOV = 1.0


def quantities(sigma, beta, null, alpha):
    """在单调正则化网格上计算 d_MD 和 psi_Q。"""
    s = np.asarray(sigma, float)
    b = np.asarray(beta, float)
    a = np.asarray(alpha, float)
    if s.shape != b.shape or np.any(s <= 0) or not np.all(np.isfinite(s + b)) or not np.isfinite(null) or null < 0:
        raise ValueError("invalid standard factors")
    if a.ndim != 1 or len(a) < 3 or np.any(a <= 0) or not np.all(np.isfinite(a)):
        raise ValueError("invalid alpha grid")
    den = s[None, :] ** 2 + a[:, None]
    d = np.sqrt(np.sum((a[:, None] ** 1.5 * b[None, :] / den ** 1.5) ** 2, axis=1) + null)
    psi = np.sqrt(np.sum((a[:, None] * s[None, :] * b[None, :] / den ** 2) ** 2, axis=1))
    return d, psi


def local_minima(values):
    """返回严格局部极小值，作为 Q/MQ 候选。"""
    values = np.asarray(values, float)
    return np.array(
        [i for i in range(1, len(values) - 1) if values[i] <= values[i - 1] and values[i] < values[i + 1]],
        dtype=np.int64,
    )


def psi_mq_tikhonov(psi_q, d_md, p0=P0_TIKHONOV):
    """根据差异曲线构造反向 MQ 包络。"""
    psi_q = np.asarray(psi_q, float)
    d_md = np.asarray(d_md, float)
    n = psi_q.size
    mq = np.empty(n)
    mq[-1] = psi_q[-1]
    expn = 2.0 / (2.0 + 1.0 / p0)
    for j in range(n - 2, -1, -1):
        floor = (d_md[j] / max(d_md[j + 1], TINY)) ** expn * mq[j + 1]
        mq[j] = max(psi_q[j], floor)
    return mq


def raus2024_rule(alpha, psi_q, d_md, sigma_min_sq, p0=P0_TIKHONOV):
    """依据归档规则合并 Q、MQ 和 HR 候选。"""
    alpha = np.asarray(alpha, float)
    psi_q = np.asarray(psi_q, float)
    d_md = np.asarray(d_md, float)
    mq = psi_mq_tikhonov(psi_q, d_md, p0)
    psi_hr = d_md / np.sqrt(np.maximum(alpha, TINY))
    in_hq = alpha >= max(float(alpha[0]), float(sigma_min_sq))
    i_hr = int(np.argmin(np.where(in_hq, psi_hr, np.inf)))
    mins = local_minima(psi_q)
    fallback = ""
    if mins.size == 0:
        hr_mins = local_minima(psi_hr)
        if hr_mins.size:
            i_h = int(hr_mins[np.argmin(psi_hr[hr_mins])])
            fallback = "psiHR_local_min"
        else:
            i_h = i_hr
            fallback = "psiHR_global"
        i_q_loc = i_h
        i_mq = i_h
        source = fallback
    else:
        i_q_loc = int(mins[np.argmin(psi_q[mins])])
        i_mq = int(mins[np.argmin(mq[mins])])
        r_h = min(1.0 / alpha[i_q_loc], max(1.0 / alpha[i_mq], 1.0 / alpha[i_hr]))
        i_h = int(np.argmin(np.abs(alpha - 1.0 / r_h)))
        source = "q" if i_h == i_q_loc else ("mq" if i_h == i_mq else ("hr" if i_h == i_hr else "nearest"))
    return {
        "index": i_h,
        "index_q_local": int(i_q_loc),
        "index_mq": int(i_mq),
        "index_hr": int(i_hr),
        "source": source,
        "fallback": fallback,
        "n_psiQ_local_minima": int(mins.size),
    }


def select(sigma, beta, null, alpha):
    """选择 alpha，并保持调用者输入的升序或降序。"""
    a = np.asarray(alpha, float)
    if not (np.all(np.diff(a) > 0) or np.all(np.diff(a) < 0)):
        raise ValueError("alpha must be strictly monotone")
    increasing = bool(len(a) < 2 or np.all(np.diff(a) > 0))
    a_work = a if increasing else a[::-1]
    d, q = quantities(sigma, beta, null, a_work)
    rau = raus2024_rule(a_work, q, d, float(np.min(np.asarray(sigma, float)) ** 2))
    i = int(rau["index"])
    ah = float(a_work[i])
    status = "defined" if np.isfinite(ah) and np.all(np.isfinite(d + q)) else "undefined_nonfinite"
    return dict(
        status=status,
        alpha=ah if status == "defined" else float("nan"),
        alpha_Q=float(a_work[rau["index_q_local"]]),
        alpha_MQ=float(a_work[rau["index_mq"]]),
        alpha_HR=float(a_work[rau["index_hr"]]),
        index=i if increasing else (len(a) - 1 - i),
        n_local_minima=int(rau["n_psiQ_local_minima"]),
        n_alpha=len(a),
        source=str(rau["source"]),
        fallback=str(rau["fallback"]),
    )
