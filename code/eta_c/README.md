# ηc 反问题：独立 MATLAB 包

只依赖 MATLAB 基础函数（QR、GSVD、Cholesky），不依赖 Python、RegTools、其他工作区或网络。已在 MATLAB R2026a 实测。无绘图、无 UI、无自动 α 选择。

## 使用

把本目录加入 MATLAB 路径，或切换到本目录：

```matlab
out = eta_c_inverse(2,80,30); % M²min, M²max, Λ，均为 GeV²
% 同时返回 out.one_delta 和 out.two_delta；每行对应一个 α。
alpha = out.alpha;
f1 = out.one_delta.f;       % f_eta_c(α)，GeV
f2 = out.two_delta.f;

% 下面的 1e-8 仅演示手动指定，不是推荐值或默认选点。
r = eta_c_spectrum(out.two_delta,1e-8);
% r.s (GeV²)，r.rho (无量纲连续谱)，r.residues (GeV²)
% δ 极点以 r.masses 和 r.residues 单独返回，不把 δ 画成有限宽峰。

% 可任选 α（不必在原扫描网格上），也可查询指定 s：
r = eta_c_spectrum(out.one_delta,[1e-8,1e-6],[0,14.99,20,30,40]);
% r.rho_query：阈值以下连续谱为零；Λ以上接冻结 LO 微扰谱。

% 测试或自定义扫描；Alpha 必须正、有限且严格递增，允许单点。
small = eta_c_inverse(2,80,30,'Alpha',10.^(-12:2:0));
save('my_result.mat','out'); % 用户自行决定是否保存
```

**三个几何输入只确定一族 f(α)，不能唯一确定 f。** 默认完整扫描 `10.^linspace(-20,0,2001)`，来自 2026-09-11 固定 Δs 全局测试；不调用 QOC、首个极小、平台或正性筛选。所有 `selected_alpha` 均为空。没有内部极小、只有一个 α 或留数为负，都不会触发“选择失败”：留数照常返回，仅当第一极点留数严格大于零时令 `f=sqrt(R1)`，否则 `f=NaN`。第二极点留数也不强制为正。

`model` 的输出包括 `alpha`（列）、`residues`（行对应 α）、`f`、`df_dlog10alpha`、`relative_residual`、`penalty`、`minrho`、`negative_count`、`factors`。谱不在完整扫描中全部落盘；独立 helper 从 `factors` 重建任意 α，无需重新装配 OPE。`out.operator` 保留网格、权重、帽子核、H、OPE 及 LO/G2/G3 分项供核验。`penalty` 为完整连续谱的 H1 范数，`relative_residual` 分母为原始 OPE 的 Borel 加权范数。

## 冻结物理与离散化

- 原 **A–A n=0 投影**的 LO+G2+G3，不是新推导的纯 Πg，也不包含额外 αs/NLO 项。
- 1δ：ηc，质量 2.9841 GeV；2δ：ηc + **χc1(3.51067 GeV)**，不是 ηc(2S)。这是原测试模型，不因打包改变谱解释。
- μ=2 GeV；使用原测试冻结值 `mc=1.1334710909 GeV`。原约定为 nf=4 领先对数质量演化 `m(μ)=1.275[αs(μ)/αs(1.275)]^(12/25)`，其中耦合来自项目三圈 RG（αs(mZ)=0.1184、mZ=91.1876、mb=4.18 GeV）。本包不另行重算/混用 J/ψ 的三圈质量，不提供 μ 可选参数。
- G2=0.038 GeV⁴、G3=0.013 GeV⁶，按源代码的凝聚项归一化冻结；s_h=14.99 GeV²。OPE 微扰积分从 4mc² 到 Λ；凝聚积分不作 Λ 截断。
- 模型：`g(M²)=sum_i R_i exp(-mi²/M²)/M² + integral_[s_h,Λ] rho(s) exp(-s/M²)/M² ds`。R1=fηc²；rho 与 g 无量纲。
- 连续谱为物理 s 上的分片线性帽子；Δs=0.1 GeV²，从 s_h 出发，最后不足一步时追加**真实 Λ**，不拉伸成固定 200 节点。端点固定 `rho(s_h)=0`、`rho(Λ)=rho_LO(Λ)`。
- Borel 网格为 1.5 起的 0.01 GeV² 格点，保留任意窗口的精确两端点并直接计算其 OPE；真梯形积分权重，不用单位矩阵，也不是统计协方差。
- 最小化 `||sqrt(W)(Kx-g)||² + α rho' H rho`；H=帽子有限元质量矩阵+刚度矩阵。极点不惩罚、不加正性约束；**α 不平方**。非零边界用 `q=Hii\(Hib*boundary)` 完成平移；使用 MATLAB 上三角 `chol(Hii)`，不能误用未经转置的下三角因子。
- 微扰 OPE 采用速度变量生成 16001 点的 s 梯形积分；凝聚项用 x∈[1e-8,1−1e-8] 的 4001 点梯形积分。帽子积分使用源代码解析式及小参数展开，不用 s 网格 trapz 近似帽子核。

输入域沿用参考 Borel 域 `1.5 ≤ M²min < M²max ≤ 100`；Λ须有至少一个内部帽子节点（Λ>15.09，极近浮点端点可能拒绝）。参考全局计算 Λ=20…70；更大 Λ 虽可输入，未由随包基准覆盖，内存和 GSVD 成本随节点数增长。极窄窗口若无法数值区分不惩罚极点会显式报错。最小 α 处病态放大舍入误差，完整输出不代表每一点具有物理可信度。

## 测试与来源

```matlab
results = runtests('test_eta_c_inverse.m');
assertSuccess(results);
```

`test_reference.mat` 是随包只读数值夹具，运行测试不访问源工作区。其来源根目录为原项目中的 `test_experiments/charmonium_4d_20260911/`，以下均为**来源标签而非运行时路径**：

- `src/physics.py`：A–A 两模型质量、尺度、精确 Borel 端点约定。
- `src/vendor/eta_aa_physics.py`：OPE、帽子核、H；源 SHA256 `f4d702e86133981dfe82afe0d8afda7a96e150b31cd28ffbfacfa8ecfc84c972`。
- `fixed_ds_global/fixed_physics.py`：固定 Δs=0.1、末端短单元。
- `fixed_ds_global/run_fixed_ds.py`：默认 2001 α、六个验证 α。
- `fixed_ds_global/fixed_ds_sweep.m`：QR+GSVD、边界平移、不惩罚极点；`loader.py`：任意 α 留数/谱回放。
- 夹具的三组几何 `(1.5,70,20)`、`(8,85,45)`、`(15,100,70)`，两模型均对照 `chunk_000/full_in_00000.mat`、`chunk_020/full_in_00247.mat`、`chunk_041/full_in_00494.mat`。六个 log10α 为 −20、−16、−12、−8、−4、0；夹具存储源 `loader.from_factors` 回放的留数/连续谱以及 Python 独立算子/OPE。另有 `(2.003,2.027,20.037)` 的非格点端点基准。

验证覆盖两模型共 36 个存储因素回放点，以及 30 个独立增广最小二乘点 `[A;sqrt(α)L] \ [b;0]`、原式残差、H1 范数、短单元、帽子闭合、端点、默认扫描、单 α 无选择、非法输入。R2026a 实测：对源因素回放的系数向量最大相对差约 **7.6e-5**（α=1e-20），α≥1e-12 最大约 **8.0e-9**；增广 LS 最大约 **6.5e-6**（α=1e-20）。小 α 差异如实保留，不宣称机器精度一致。
