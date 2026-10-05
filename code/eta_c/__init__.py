"""Independent eta_c n=0 spectral sum-rule solver."""
# 只导出稳定的公开 API；实现模块仍可供复现和单元测试直接导入。
from .solver import solve, lapacke_info

__all__ = ["solve", "lapacke_info"]
