"""实际紧凑 LAPACKE GSVD，不改变分解算术。
这里只保留所需的库探测和分解函数。
"""
from __future__ import annotations
from dataclasses import dataclass
import ctypes
import ctypes.util
from functools import lru_cache
import glob
import os
from pathlib import Path
import sys
import sysconfig
from typing import Any
import numpy as np
from .physics import _dot
_LAPACK_COL_MAJOR = 102
_LAPACK_INT = ctypes.c_int
_DOUBLE_PTR = ctypes.POINTER(ctypes.c_double)
_INT_PTR = ctypes.POINTER(_LAPACK_INT)

@dataclass(frozen=True)
class CompactGSVD:
    """采用 Regularization Tools ``cgsvd`` 排序的紧凑 GSVD。

    前 ``p`` 列是正则化模态，最后 ``k`` 列是未正则化振幅模态。
    因而有 ``A = U C W``、``B = V S W``，并且 ``X = inv(W)``。
    """

    U: np.ndarray
    V: np.ndarray
    X: np.ndarray
    W: np.ndarray
    sigma: np.ndarray
    mu: np.ndarray
    sigma_unregularized: np.ndarray
    source: str
    library: str

    @property
    def p(self) -> int:
        return int(self.sigma.size)

    @property
    def k(self) -> int:
        return int(self.sigma_unregularized.size)


def _candidate_libraries() -> list[str]:
    """返回可能导出 LP64 LAPACKE_dggsvd3 的共享库候选。"""

    candidates: list[str] = []
    explicit = os.environ.get("ETA_C_LAPACKE_LIBRARY")
    if explicit:
        candidates.append(explicit)
    for name in ("lapacke", "lapack", "openblas", "blas", "mkl_rt"):
        found = ctypes.util.find_library(name)
        if found:
            candidates.append(found)
    if os.name == "nt":
        patterns = ("mkl_rt*.dll", "*lapacke*.dll", "*openblas*.dll", "*lapack*.dll", "*blas*.dll")
        directories = [Path(sys.prefix) / "Library" / "bin", Path(sys.prefix) / "DLLs", Path(sys.prefix) / "bin", Path(sys.executable).resolve().parent]
        directories.extend(Path(item) for item in os.environ.get("PATH", "").split(os.pathsep) if item)
    else:
        patterns = ("liblapacke.so*", "liblapack.so*", "libopenblas.so*", "libblas.so*", "liblapacke*.dylib", "liblapack*.dylib", "libopenblas*.dylib")
        directories = [Path(sys.prefix) / "lib", Path(sys.prefix) / "lib64", Path("/usr/lib"), Path("/usr/lib64"), Path("/usr/local/lib"), Path("/usr/local/lib64"), Path("/lib"), Path("/lib64"), Path("/opt/lib"), Path("/opt/homebrew/opt/lapack/lib"), Path("/usr/local/opt/lapack/lib")]
        for root in (Path("/usr/lib"), Path("/lib")):
            directories.extend(root.glob("*-linux-gnu"))
        multiarch = sysconfig.get_config_var("MULTIARCH")
        if multiarch:
            directories.extend((Path("/usr/lib") / multiarch, Path("/lib") / multiarch))
    for directory in directories:
        for pattern in patterns:
            candidates.extend(sorted(glob.glob(str(directory / pattern))))
    numpy_root = Path(np.__file__).resolve().parent
    package_dirs = (numpy_root.parent / "numpy.libs", numpy_root / ".dylibs", numpy_root.parent / ".dylibs")
    for entry in sys.path:
        root = Path(entry)
        package_dirs += (root / "scipy" / ".dylibs", root / "scipy.libs", root / "numpy.libs")
    package_patterns = patterns + ("*lapack*", "*openblas*", "*blas*", "mkl_rt*.dll")
    for directory in package_dirs:
        if directory.is_dir():
            for pattern in package_patterns:
                candidates.extend(sorted(glob.glob(str(directory / pattern))))
    seen: set[str] = set()
    return [item for item in candidates if item and not (item in seen or seen.add(item))]


@lru_cache(maxsize=1)
def _load_dggsvd3() -> tuple[Any, str] | None:
    """不依赖 SciPy 封装，定位 LAPACKE 的 ``dggsvd3``。"""

    # 直接使用 C LAPACKE 接口，使 GSVD 约定匹配 MATLAB 参考实现，
    # 不依赖封装层可能进行的重排序。
    for candidate in _candidate_libraries():
        try:
            library = ctypes.CDLL(candidate)
            routine = library.LAPACKE_dggsvd3
        except (OSError, AttributeError):
            continue
        if "ILP64" in os.environ.get("MKL_INTERFACE_LAYER", "").upper() and "mkl_rt" in Path(candidate).name.lower():
            if candidate == os.environ.get("ETA_C_LAPACKE_LIBRARY"):
                raise RuntimeError("The selected MKL interface is ILP64; this GSVD interface requires LP64.")
            continue
        # The unsuffixed MKL LAPACKE entry is the LP64 C interface; never call *_64.
        config = getattr(library, "openblas_get_config", None)
        if config is not None:
            config.restype = ctypes.c_char_p
            if b"USE64BITINT" in (config() or b""):
                configured = os.environ.get("ETA_C_LAPACKE_LIBRARY")
                same = configured and (candidate == configured or Path(candidate).name == Path(configured).name)
                if same:
                    raise RuntimeError("The specified OpenBLAS uses ILP64; this GSVD interface requires LP64.")
                continue
        routine.restype = _LAPACK_INT
        routine.argtypes = [
            _LAPACK_INT,
            ctypes.c_char,
            ctypes.c_char,
            ctypes.c_char,
            _LAPACK_INT,
            _LAPACK_INT,
            _LAPACK_INT,
            _INT_PTR,
            _INT_PTR,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _DOUBLE_PTR,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _DOUBLE_PTR,
            _LAPACK_INT,
            _INT_PTR,
        ]
        return routine, str(candidate)
    return None


def lapacke_gsvd_available() -> bool:
    """返回当前是否能找到导出 ``dggsvd3`` 的 LAPACKE 实现。"""

    return _load_dggsvd3() is not None


def _compact_from_lapacke_direct(A: np.ndarray, B: np.ndarray) -> CompactGSVD:
    """直接调用 LAPACKE GSVD，保留此路径用于验证 QR 约化。"""
    loaded = _load_dggsvd3()
    if loaded is None:
        raise RuntimeError("LAPACKE_dggsvd3 was not found")
    routine, library_path = loaded
    m, n = A.shape
    p, n_b = B.shape
    if n_b != n or m < n or p >= n:
        raise ValueError("compact M3 GSVD requires m >= n > p")

    a_work = np.array(A, dtype=np.float64, order="F", copy=True)
    b_work = np.array(B, dtype=np.float64, order="F", copy=True)
    alpha = np.empty(n, dtype=np.float64)
    beta = np.empty(n, dtype=np.float64)
    u = np.empty((m, m), dtype=np.float64, order="F")
    v = np.empty((p, p), dtype=np.float64, order="F")
    q = np.empty((n, n), dtype=np.float64, order="F")
    iwork = np.empty(n, dtype=np.intc)
    k = _LAPACK_INT()
    ell = _LAPACK_INT()
    info = routine(
        _LAPACK_COL_MAJOR,
        b"U",
        b"V",
        b"Q",
        m,
        n,
        p,
        ctypes.byref(k),
        ctypes.byref(ell),
        a_work.ctypes.data_as(_DOUBLE_PTR),
        m,
        b_work.ctypes.data_as(_DOUBLE_PTR),
        p,
        alpha.ctypes.data_as(_DOUBLE_PTR),
        beta.ctypes.data_as(_DOUBLE_PTR),
        u.ctypes.data_as(_DOUBLE_PTR),
        m,
        v.ctypes.data_as(_DOUBLE_PTR),
        p,
        q.ctypes.data_as(_DOUBLE_PTR),
        n,
        iwork.ctypes.data_as(_INT_PTR),
    )
    if info != 0:
        raise RuntimeError(f"LAPACKE_dggsvd3 failed with INFO={info}")
    if k.value + ell.value != n or ell.value != p or k.value != n - p:
        raise RuntimeError(
            "unexpected dggsvd3 ranks: "
            f"k={k.value}, l={ell.value}, expected k={n - p}, l={p}"
        )
    if not np.all(np.isfinite(alpha)) or not np.all(np.isfinite(beta)):
        raise RuntimeError("LAPACKE_dggsvd3 returned non-finite generalized singular values")

    # LAPACKE 返回三角因子和广义奇异值；对紧凑正则化模态重排后，
    # 才向上层暴露逆变换。
    raw_r = np.triu(a_work[:n, :n])
    if np.any(np.abs(np.diag(raw_r)) <= np.finfo(float).tiny):
        raise RuntimeError("LAPACKE_dggsvd3 returned a singular triangular factor")
    raw_w = _dot(raw_r, q.T)
    compact_raw = np.arange(k.value, n)
    null_raw = np.arange(k.value)
    ordering = np.r_[compact_raw, null_raw]
    x_raw = _dot(q, np.linalg.solve(raw_r, np.eye(n)))
    return CompactGSVD(
        U=np.asarray(u[:, :n])[:, ordering],
        V=np.asarray(v),
        X=x_raw[:, ordering],
        W=raw_w[ordering, :],
        sigma=np.asarray(alpha[compact_raw]),
        mu=np.asarray(beta[compact_raw]),
        sigma_unregularized=np.asarray(alpha[null_raw]),
        source="LAPACKE_dggsvd3",
        library=library_path,
    )


def _compact_from_lapacke(A: np.ndarray, B: np.ndarray, *, qr_reduce: bool) -> CompactGSVD:
    """执行实际 LAPACKE GSVD，可选精确经济型 QR 约化。

    对物理 M3 系统，``A=K0`` 是高矩阵。若经济型分解给出 ``A=Q_A R_A``，
    最小化 ``||A x-b||`` 与最小化 ``||R_A x-Q_A.T b||`` 只相差与 alpha 无关的
    正交残差。因此 ``(R_A, B)`` 的 GSVD 正好足以求 Tikhonov 解；再用
    ``Q_A @ U_R`` 映回可保留有效投影、重构恒等式和解析导数，同时避免
    LAPACK 分配稠密的 m-by-m U 矩阵。
    """

    # 经济型 QR 只去除高矩阵数据中的 alpha 无关正交残差，
    # 显著减少稠密 U 矩阵的内存分配。
    if not qr_reduce:
        return _compact_from_lapacke_direct(A, B)
    q_active, r_active = np.linalg.qr(A, mode="reduced")
    reduced = _compact_from_lapacke_direct(r_active, B)
    return CompactGSVD(
        U=_dot(q_active, reduced.U),
        V=reduced.V,
        X=reduced.X,
        W=reduced.W,
        sigma=reduced.sigma,
        mu=reduced.mu,
        sigma_unregularized=reduced.sigma_unregularized,
        source="LAPACKE_dggsvd3_economy_qr",
        library=reduced.library,
    )


def compact_gsvd(A, B, *, method="lapacke", qr_reduce=True):
    """暴露唯一支持的物理 GSVD 后端。"""
    if method != "lapacke":
        raise ValueError("Physical solves require LAPACKE; no Gram-matrix fallback")
    return _compact_from_lapacke(np.asarray(A, dtype=np.float64), np.asarray(B, dtype=np.float64), qr_reduce=qr_reduce)
