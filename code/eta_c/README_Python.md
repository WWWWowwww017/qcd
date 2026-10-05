# ηc 独立 Python 核心

本目录新增的 Python 包沿用相邻 `jpsi/` 的模块分工、接口、CLI、LAPACKE GSVD 和自动 α 选择。原 MATLAB 文件、`README.md`、`test_reference.mat` 不改动。运行不依赖 jpsi、原工作区、MATLAB、历史结果或网络；无绘图。

## 最简单的调用

从站点根目录运行：

```python
from eta_c import solve, lapacke_info
print(lapacke_info())
r = solve(5.0, 80.0, 30.0)
one = r["channels"]["1delta"]
print(one["f_eta_c"], one["alpha_star"])
s, rho = one["s"], one["rho"]
```

三个必填量为 `M2_min, M2_max, Lambda`，单位均为 GeV²。默认同时计算 `1delta`、`2delta`，两者分别选 α。可用 `solve(5,80,30,channel="2delta")` 只算一个。API 不写文件。

从任何目录运行 CLI：

```bash
python -B code/eta_c/run.py 5 80 30
```

从站点根目录进入 `code` 后也可运行 `python -B -m eta_c 5 80 30`。可选 `--channel 1delta` / `--channel 2delta`、`--output ./results/`。默认在本目录 `results/eta_c_UTC时间戳_随机后缀/` 新建 `result.npz` 和 `summary.json`，并打印 JSON；不覆盖旧结果。JSON 非有限值记为 `null`，NPZ 保留 NaN。

## 步长和物理

“步长还是 0.1”按 Jpsi 当前实际代码解释为 **谱网格 Δs=0.1 GeV²**，反演与输出用同一网格；**Borel ΔM²=0.01 GeV²**，与 Jpsi 和原 ηc MATLAB 一致，不把两个步长混为一谈。谱从 14.99 起固定步进，末端精确保留 Λ，最后一格可短于 0.1；不拉伸成 200 点。Λ=30 时 152 节点、150 个内部连续谱未知数，末格 0.01 GeV²。

物理来自原工作区 `test_experiments/charmonium_4d_20260911/src/vendor/eta_aa_physics.py`，固定 Δs 约定来自 `fixed_ds_global/fixed_physics.py`；这两者只是来源标签，不是运行时依赖。

- 原 **A–A n=0 LO+G2+G3**，不是新推导纯 Πg；不加 Jpsi 的 `1+alpha_s/pi`。
- μ=2 GeV，mc=1.1334710909 GeV（原领先对数质量演化冻结值），G2=0.038 GeV⁴、G3=0.013 GeV⁶，S_H=14.99 GeV²。
- G3 的 b3 多项式保留 ηc 的 `+x*xb/3`，不是 Jpsi 的 `-11*x*xb/3`。
- `1delta` 质量 2.9841 GeV；`2delta` 再加 **χc1(3.51067 GeV)**，不是 ηc(2S)。延续原测试模型，不重新解释谱内容。
- 微扰速度网格 16001 点、凝聚 Feynman 参数网格 4001 点（1e-8 至 1−1e-8）；微扰从 4mc² 积到 Λ，凝聚项不截断于 Λ。
- 解析 hat 核、物理 s 上 H=质量矩阵+刚度矩阵、Borel 真梯形权重；端点 `rho(S_H)=0`、`rho(Lambda)=rho_LO(Lambda)`，保留非零边界交叉项平移。
- 最小化 `||sqrt(W)(Kz-g)||²+alpha*hat.T@H@hat`；alpha **不平方**，极点不惩罚，不约束留数或连续谱正性。非正留数的 f 为 NaN，不取绝对值或 clip。

## 与 Jpsi 一致的 α 支持

只提供 Jpsi 当前的 **`raus2024_etac_equiv_p0_1`**：Q/MQ 在 ψQ 内部局部极小候选中选择；MQ 用 p0=1、2/3 指数递推；HR 全局候选限制 α≥max(网格下界,σmin²)；合成 rH=min(rQ,max(rMQ,rHR))。无 ψQ 内部极小时按原代码回退 HR 局部/受限全局极小，明确返回原因。不是 area-rule-2、旧严格 Raus 或只取 QOC 最小。

α 网格也原样沿用 Jpsi：下界从 max(−80,min(−50,2log10σmin−8)) 起，约 0.01 的 log10 步长；必要时每次左扩 5 至 −80，以覆盖小 α 峰；再按 `10**linspace(log10(first_alpha),0,n_alpha)` 构造最终网格。没有新加手选 α、全 α 扫描或其他选择开关。

**与原 MATLAB 的差别是有意的**：MATLAB 默认输出 −20…0 的整族 α，完全不选参；Python 按 Jpsi 自动选一个 α，保留因子与选择诊断。MATLAB 同网格验证是在明确相同 α 上比较解，不假称 MATLAB 验证了自动规则。

## 输入与返回

与 Jpsi 一致，允许有限正 Borel 上下界，要求 max>min，Λ>S_H 且超过所含极点质量平方；至少 3 个 hat 节点，Borel 点数至少 `n_hats−2+n_poles`。Λ=30 时单/双 δ 分别需 151/152 个 Borel 点。不暗中减少 hats、不采用近似 GSVD 后备。

Borel 两端精确保留；与 `1.5+0.01*k` 对齐时保留原内部算式，否则从用户下界每次加 0.01，最后可有短格。**非格点下界的内部网格沿用 Jpsi，不等同于 MATLAB 永远锚定 1.5 的内点网格**；比较时须先对齐网格。Python 不继承 MATLAB 的 [1.5,100] 硬输入域或欠定窗口支持。域外可算不代表物理可信，GSVD 失败/非有限值明确报错。

返回与 Jpsi 同层级：

- `metadata`：几何、n_hats、ds、dM2、物理常数、`ope_order`、规则与真实 LAPACKE 库路径；无伪造 alpha_s 常数。
- `M2, weights, OPE`，后者包含 `pert,D4,D6,total`。
- `channels[模型]`：`masses,residues,f,f_eta_c,s,rho,s_hat,hat,boundary`；`residues=f²` 单位 GeV²，f/masses 为 GeV，连续谱无量纲。δ 极点单独返回，不伪装成有限宽峰。
- `alpha_star,alpha_diagnostics`：α 网格、Q/MQ/HR 候选、曲线、来源/回退、边界/小峰覆盖检查。
- `backsub,backsub_poles,backsub_continuum,residual`、加权/相对/不加权残差、H1 penalty、objective、模态与直接残差差值、负谱计数，以及可复算的 `factors`。

NPZ 用 `__` 展开层级，例如 `channels__1delta__rho`、`channels__1delta__f_eta_c`、`OPE__total`；可 `np.load(path,allow_pickle=False)`，无对象数组。

## 依赖与验证

Python 3.9+、NumPy、已有 LP64 `LAPACKE_dggsvd3` 库。测试另需 SciPy 读取原 `test_reference.mat`。程序零环境变量自动探测 Windows 的 `sys.prefix\Library\bin\mkl_rt*.dll`、NumPy/SciPy 库目录，以及 Ubuntu/Debian 多架构目录（含 `/usr/lib/x86_64-linux-gnu`）；MKL 只调用未加 `_64` 后缀的 LP64 C 入口，OpenBLAS 会拒绝 `USE64BITINT`。若发行版没有可用库，可先执行 `python -m pip install numpy scipy`，Ubuntu 也可安装 `liblapacke-dev`，或用 `ETA_C_LAPACKE_LIBRARY=<existing-library-path>` 明确指定；只用 eta_c 环境变量，不读取 JPSI 环境变量。QR+实际 GSVD，不用 Gram 特征分解或玩具 NumPy 替代。`np.dot` 后检查有限性，不屏蔽数值告警。

```bash
python -B -m unittest discover -s code/eta_c/tests -v
```

2026-09-11 本机 **11 项测试全部通过（8.198 s）**：ηc 常数与 LO、三组原基准共 36 个固定 α 解、独立 MATLAB 同网格 OPE/谱/留数/残差/H1、hat 闭合与 H1、任意上下界与一阶条件、实际自由维数、非法输入、Raus 升降序/回退、真实 LAPACKE、非覆盖保存及 pickle-free NPZ。

新增 `tests/reference_matlab.npz` 来自 MATLAB R2026a 的 `eta_c_inverse(5,80,30,'Alpha',sort(alpha))`，并用 `eta_c_spectrum` 回放，alpha 为 `[1e-20,1e-16,1e-12,1e-8,1e-4,1,1.0964781961431828e-6,4.073802778041122e-5]`。152 个谱点、7501 个 Borel 点；保留 OPE 分项、算子样本及系数/残差/惩罚作为测试夹具，不用于运行求解。最小 α=1e-20 的系数相对差最大 7.29e-6；α≥1e-12 最大 5.15e-10；各模型选中 α 处分别约 1.14e-13、5.53e-15。极小 α 病态差异如实保留，不宣称逐位相同。

CLI `(5,80,30)` 实测结果：

| 模型 | alpha | f_eta_c (GeV) | 第二极点 f (GeV) |
|---|---:|---:|---:|
| 1delta | 1.0964781961431828e-6 | 0.4009772959677418 | — |
| 2delta | 4.073802778041122e-5 | 0.3219558624907369 | 0.3797736122856722 |

两者均选 Q 候选，无回退，非 α 端点，小峰覆盖成功。只是本窗口无约束反演结果，不代表已完成窗口可信性或正性判定。

文件分工与 Jpsi 相同：`physics.py` 物理，`solver.py` 统一 API，`factors.py` 边界/GSVD 投影，`gsvd_backend.py` 库探测/真实分解，`alpha_grid.py` 网格，`raus2024.py` 选择，`__init__.py` 导出，`__main__.py` 保存/CLI，`run.py` 任意目录入口。
