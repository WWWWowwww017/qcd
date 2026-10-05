"""Independent J/psi n=0 spectral sum-rule solver."""
# 保持公开 API 精简；测试仍可直接导入实现模块。
from .solver import solve, lapacke_info

__all__ = ["solve", "lapacke_info"]
