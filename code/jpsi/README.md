# J/ψ 谱函数与衰变常数反演（独立 Python 核心）

本目录独立计算 J/ψ **n=0 谱 sum rule**，给出连续谱、窄极点留数及衰变常数。**不是 LCDA 矩提取，也不是 LCDA 重建。** 运行不读取旧工作区、历史 NPZ、集群文件或 skill；无绘图步骤。

## 最简单的调用

从站点根目录运行：

```python
from jpsi import solve, lapacke_info
print(lapacke_info())
r = solve(5.0, 80.0, 30.0)   # 默认同时求 1delta、2delta
one = r["channels"]["1delta"]
print(one["f"], one["alpha_star"])
s, rho = one["s"], one["rho"]
```

只求一个模型：`solve(5.0, 80.0, 30.0, channel="2delta")`。`solve` 不写文件，单模型也保持相同的嵌套返回结构。

CLI 只需三个位置参数，上下界及 Λ：

```bash
python -B code/jpsi/run.py 5 80 30
```

可选 `--channel 1delta` / `--channel 2delta` 和 `--output ./results/`。从站点根目录进入 `code` 后也可运行 `python -B -m jpsi 5 80 30`。

默认在本目录 `results/` 下独占创建 `jpsi_时间戳_随机后缀/`，输出 `result.npz`、简短 `summary.json`，并向终端打印该 JSON。不会静默覆盖已有结果。JSON 中非有限值写成标准 JSON 的 `null`；NPZ 中保留 NaN。

## 依赖与动态库

Python 3.9+、NumPy，以及提供 **LP64 `LAPACKE_dggsvd3`** 的现有共享库。核心不调用 SciPy 数值接口，但可复用它附带的 OpenBLAS 动态库。程序不安装依赖，不退回 Gram 矩阵特征分解（它不适合本问题的极小广义奇异值）。

`lapacke_info()` 返回 `available` 和探测到的库路径。零环境变量时自动探测 Windows 的 `sys.prefix\Library\bin\mkl_rt*.dll`、NumPy/SciPy 库目录，以及 Ubuntu/Debian 多架构目录（含 `/usr/lib/x86_64-linux-gnu`）；MKL 只调用未加 `_64` 后缀的 LP64 C 入口，OpenBLAS 会拒绝 `USE64BITINT`。标准环境没有可用库时，可先执行 `python -m pip install numpy scipy`，Ubuntu 也可安装 `liblapacke-dev`，或明确指定已有库：

```bash
JPSI_LAPACKE_LIBRARY=<existing-library-path> python -B code/jpsi/run.py 5 80 30
```

兼容原环境变量 `ETA_C_LAPACKE_LIBRARY`，但无需设置它。必须是 LP64 C LAPACKE 接口，不可把只有 Fortran 接口或 ILP64 接口的库当作替代。

## 单位、物理与默认设置

- 输入 `M2_min, M2_max, Lambda` 均为 **GeV²**。输出 `M2,s,s_hat` 同单位；极点质量及 `f` 为 GeV，留数 `residues=f²` 为 GeV²，连续谱 `rho` 无量纲。
- 连续区间 `[S_H, Lambda]`，`S_H=13.9105153024 GeV²`；实际反演 hat 网格与输出 `s` 网格均从 `S_H` 起以 **Δs=0.1 GeV²** 步进，精确保留 Λ 末端，最后不足一格时保留短 interval（浮点舍入级重复端点合并）。不再固定 200 节点，也不再输出 0.01 的细谱；`rho` 与 `hat` 在同一网格上。节点数随 Λ 变化，例如 Λ=30 时 162 个节点、160 个内部连续谱未知数，最后 interval 为 0.0894846976 GeV²。
- `sgrid_hats(L, nodes=200)` 仅保留为显式历史复核接口；`solve` 默认不使用它。
- `1delta` 极点质量为 3.096900 GeV；`2delta` 再加 3.686097 GeV。每个模型独立进行 GSVD 和选参，留数不受正则项惩罚。
- 固定 μ=2 GeV，三圈跑动的冻结数值 `mc=1.0990536546335794 GeV`、`alpha_s=0.3029450988131294`；本程序不根据 Borel 参数重新跑动。
- `G2=0.038 GeV⁴`、`G3=0.013 GeV⁶`。原 `aa_physics.py` 的 n=0 微扰 `LO*(1+alpha_s/pi)`、D4、D6 公式和积分网格原样保留：微扰速度网格 16001 点，Feynman 参数网格 4001 点，x 两端为 `1e-8` 和 `1-1e-8`。
- Borel 步长默认 0.01 GeV²。允许任意有限正上下界，精确包含用户给出的两个 float64 端点，不整窗平移或裁成旧扫描窗。最后不足一步时保留短区间。历史 `1.5+0.01*k` 网格的内部点沿用同一浮点算式；用户端点精确值优先。
- 校验 `max>min`、`Lambda>S_H` 且大于所含极点质量平方；实际网格少于 3 个 hat 节点时清楚拒绝。Borel 点数须至少为实际自由未知数个数 `n_hats-2+n_poles`：Λ=30 时 1delta 要 161 点，2delta 要 162 点。过短窗口报错而不暗中减少 hats。
- 数据范数用 Borel 梯形积分权重；连续谱惩罚为物理 s 网格上的 `H=M+M1`（hat 质量矩阵加刚度矩阵）。两端固定为 `rho(S_H)=0`、`rho(Lambda)=rho_pert(Lambda)`，严格保留端点交叉项平移。
- 算法最小化 `||sqrt(W)(Kz-g)||² + alpha*rho_hat.T@H@rho_hat`，此处 **alpha 不是 alpha²**。GSVD 前用 economy QR 减少高矩阵内存。

输入合法不等于所有极端窗口都能稳定反演；数值非有限或 GSVD 失败会明确报错，不伪造结果。

## 三种历史规则的区别（这里只交付最新一种）

1. **area-rule-2**：旧 Q-curve 的几何面积选参，不是本程序默认，也没有暴露对应开关。
2. **旧 J/ψ `raus2024.py`**：严格局部极小、对平坦/非唯一极小或无内部极小可返回 undefined；Q 候选取全网格最小，HR 未加最新的广义奇异值阈值限制。
3. **本程序：`raus2024_etac_equiv_p0_1`**，原 `raus2024_etac.py`：Q/MQ 均在 ψQ 局部极小候选上选取；MQ 采用 `p0=1` 的 2/3 指数递推包络；HR 全局候选限制 `alpha>=max(alpha_grid[0],sigma_min²)`；合成 `r_H=min(r_Q,max(r_MQ,r_HR))`、`alpha=1/r_H`。无 ψQ 局部极小时沿用原代码的 HR 局部极小/受限全局回退。回退原因和来源均返回。

网格也属于默认规则的一部分：先按旧小 α 峰覆盖逻辑取对数下界 `max(-80,min(-50,2*log10(sigma_min)-8))`，约 0.01 的 log10 步长；未覆盖左侧小 α 峰则每次向左扩 5，止于 -80。随后**严格按最新 `watch_and_select.py`** 取 `10**linspace(log10(first_alpha),0,n_alpha)`，而不是直接换一套网格。覆盖状态作为诊断返回；未覆盖不会暗中改选面积规则。

复制进来的 `raus2024.py` 文件名虽短，其内容是最新 `raus2024_etac.py`，不是上述第 2 种旧实现。

## 返回结构与 NPZ

`solve` 返回字典：

- `metadata`：输入、常数、schema、规则和 LAPACKE 路径；`n_hats` 为实际节点数，`ds=0.1` 为标称 s 步长（末 interval 可短于它），`dM2=0.01` 不变。
- `M2, weights`：本次实际窗口与梯形权重。
- `OPE`：`pert, D4, D6, total`，长度等于 Borel 点数。
- `channels["1delta"]` / `channels["2delta"]`：
  - `masses, residues, f, f_jpsi`（留数/衰变常数只含该模型实际极点，不补假极点）。
  - `s, rho` 是**连续谱**，δ 极点单独保存在 masses/residues，不将 δ 峰画成有限宽度函数；`s_hat, hat, boundary` 为节点解与端点。
  - `alpha_star` 及 `alpha_diagnostics`：选中索引、网格、Q/MQ/HR 候选、`d_MD, psi_Q, psi_MQ, psi_HR`、来源、回退、是否网格边界、小 α 峰覆盖、最小广义奇异值平方等。
  - `backsub` 为总反代，`backsub_poles, backsub_continuum` 分项；`residual=backsub-OPE.total`。
  - `weighted_residual, relative_weighted_residual, unweighted_residual`、`penalty, objective`；`residual_factor` 和 `residual_factor_discrepancy` 检查模态残差与直接反代残差之差。
  - `min_rho, negative_rho_count`；`factors` 包含 sigma/mu、sigma_tilde、未正则化模式、投影 beta、正交残差平方 null、X、边界平移 q 等可复核数据。

不施加正性约束，**负留数和负连续谱原样保留**。正留数给 `f=sqrt(residue)`；非正留数按历史约定给 `f=NaN`（包括恰好 0），不 clip 成 0、不取绝对值。这表示该无约束拟合未给出可报告的正实衰变常数，并非自动物理可接受的结果。

NPZ 用双下划线展开嵌套键，全是普通数值/字符串数组，无 pickle：

```python
import numpy as np
with np.load("结果子目录/result.npz", allow_pickle=False) as z:
    s = z["channels__1delta__s"]
    rho = z["channels__1delta__rho"]
    f = z["channels__1delta__f"]
    ope = z["OPE__total"]
```

## 文件分工

`solver.py` 为统一接口；`physics.py` 保留物理公式；`factors.py` 为端点消元与 GSVD 投影；`gsvd_backend.py` 只保留所需的 LAPACKE 分解/探测代码，不重写分解算法；`alpha_grid.py` 仅复制历史网格所需函数；`raus2024.py` 为最新选参；`__main__.py` 和 `run.py` 提供输出与 CLI。运行时不依赖旧 NPZ、集群调度器、旧目录包装器或绘图工具；`tests/` 仅携带少量历史结果作为独立回归基准。

## 已运行的验证（2026-09-11）

本机 `/usr/bin/python3`（Python 3.9.6、NumPy 2.0.2、SciPy 附带的 LP64 OpenBLAS）实算通过。测试启用 `-W error`，没有屏蔽数值告警。

```bash
python -B -m unittest discover -s code/jpsi/tests -v
```

9 项测试全部通过（`Ran 9 tests in 7.326s`，`OK`），覆盖 15 组 OPE 分项及独立积分基准、帽核与 H1 恒等式、选参升降序及回退、非法参数、默认 0.1 网格端点/短区间/实际自由维数，以及任意窗口 `(5.003,20.007,30.123)` 的端点与最优解一阶条件。该任意 Λ 测试使用新默认 164 节点，最后 interval 为 0.0124846976 GeV²。

历史 9 个单/双 δ 窗口单独通过 mock 显式指定 `sgrid_hats(L,nodes=200)` 及历史输出网格复核，fixture 保持不变；**不要求新离散化硬匹配旧谱数组**。历史结果基准来自 `raus2024_etac_equiv_local`，不是更早的 area-rule-2 或严格 Raus 版本。

### 历史 200 节点实算（不是当前默认）

例 `(M2_min,M2_max,Lambda)=(5,80,30)`，反演为 200 个等距 hat 节点、输出 Δs=0.01：

| 模型 | alpha | f_J/psi (GeV) | f_psi(2S) (GeV) |
|---|---:|---:|---:|
| 1delta | 3.630780547701032e-5 | 0.4793037976586482 | — |
| 2delta | 3.548133892335789e-5 | 0.4549453640611985 | 0.27030329434510975 |

这两项 f_J/psi 与历史结果差分别为 `9.83e-14`、`2.62e-13 GeV`。不同 BLAS/LAPACK 在病态模态上不保证逐位一致；回归检验实际选出的 alpha、留数、谱和回代误差，而不是仅看程序退出码。

### 新默认 Δs=0.1 实算

同一窗口运行 `/usr/bin/python3 -B -W error run.py 5 80 30` 成功（无绘图、自动探测本机 LP64）。实际反演及输出均为 162 节点；末两点为 29.9105153024、30.0 GeV²，最后 interval 为 0.0894846976 GeV²。Borel 仍为 0.01 步长、7501 点。

| 模型 | alpha | f_J/psi (GeV) | f_psi(2S) (GeV) | 加权残差 |
|---|---:|---:|---:|---:|
| 1delta | 3.630780547701032e-5 | 0.4793043028331749 | — | 1.0900467435455031e-4 |
| 2delta | 3.548133892335789e-5 | 0.4549458259873889 | 0.27030113328643024 | 1.7049881179731103e-4 |

两模型均由 Q 候选选中 α，无回退、不在 α 网格边界，小 α 峰覆盖成功。数值变化来自 s 离散化改变，不是物理或 α 算法改动。新输出独立保存在 `results/jpsi_20260911T081159Z_atve8rw8/` 的 `result.npz` 与 `summary.json`，没有覆盖旧实算。

本机 NumPy/Accelerate 的 `@` 路径对有界矩阵会产生虚假的浮点告警。核心统一使用 `np.dot`，并立即检查乘积有限性；真实 Inf/NaN 仍报错，没有忽略 warning。物理算式、积分节点和 GSVD 选参未改变。

**范围限制：** 能计算任意合法窗口，不代表已证明该窗口满足 OPE 收敛或物理稳定性判据；本工具不自动寻找可信 Borel 平台。
