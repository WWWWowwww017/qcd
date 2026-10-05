# Inverse Problem Approach for Non-Perturbative QCD

这是一个可直接部署到 GitHub Pages 的静态研究展示页，内容依据：

- `docs/Foundation.pdf`
- `docs/Application.pdf`

## 部署到 GitHub Pages

1. 将本目录内容放入一个 GitHub 仓库的根目录。
2. 推送到 `main` 分支。
3. 在仓库设置中打开 **Settings → Pages**，将构建来源设为 **GitHub Actions**。
4. 工作流完成后，GitHub Pages 会生成公开访问地址。

站点不依赖构建工具；`index.html`、`code-lab.html`、样式、脚本、数据和 `assets/` 均为静态文件，可直接托管。

## 代码与结果工作台

打开 `code-lab.html` 可以：

- 在左侧切换并只读阅读归档中的 MATLAB、Python 和 README 文件。
- 在右侧切换谱函数、OPE、α 诊断和数值摘要，并在 `1δ / 2δ` 模型间切换。
- 在 ηc 项目中比较 MATLAB 参考夹具和 Python 保存结果；J/ψ 项目展示归档中的 Python 结果。

工作台展示的是 `code/results/` 中由归档导出的结果回放，不在浏览器内执行 MATLAB 或重新计算求解器。页面保留原始结果中的负谱节点、δ 极点独立表示以及源文件路径和 SHA-256 信息；浏览器原生的选择、复制、右键和打印行为不再被页面拦截。

完整研究归档通过与发布版本对应的 GitHub Release 附件 `qcd-code-archive.zip` 获取，不再放在 Pages 站点目录。附件包含完整/历史 MATLAB 与 Python 源码、数值夹具和保存结果；本轮只移动原始 ZIP，没有重打包。公开树保留工作台所需的轻量源码、预览 JSON，以及四个测试夹具：`code/eta_c/test_reference.mat`、`code/eta_c/tests/reference_matlab.npz`、`code/jpsi/tests/reference_latest.npz` 和 `code/jpsi/tests/reference_ope.json`。归档中的 `__MACOSX` / `.DS_Store` 条目尚未清理，需作者确认后再生成新的 Release 附件。

工作台运行时只读取公开树中的 `code/lab-data.js`、`code/results/*.json` 和源文件回退路径，不依赖 Release ZIP。页面中的论文图表和四个首页数值是 Application 论文结果；ηc/J/ψ 是同一研究框架下的独立归档回放，不应解读为该论文的重构结果。

本轮性能测量（2026-10-01，localhost，文本类资源按 gzip level 9、二进制资源按原字节计，只统计初始 `script`、`img` 和 stylesheet/icon 资源）显示：首页首屏资源为 238,654 B，工作台首屏资源为 699,791 B，其中 `code/lab-data.js` 为 256,968 B。按项目拆分的静态估算为 ηc 数据 154,223 B、J/ψ 数据 104,186 B；即使按默认 ηc 首屏计算，工作台只减少约 14.7%，低于本轮设定的 30% 门槛，因此保留单一 `lab-data.js`，接受其成本以保持当前回退逻辑和静态结构简单。LCP 未能在本机测得：Playwright/Chromium 启动返回 `browserType.launch: spawn EPERM`，没有用估算值替代实测。

## 内容说明

页面将“理论基础”“数值应用”“可验证结果”和“方法边界”分开呈现。数值结果、作者、日期和限定条件均按两份 PDF 整理；完整公式、图表与参考文献请以 PDF 为准。

## 许可证与引用

本站按学术交流用途公开；论文结果、公式、图表和代码的再使用应遵循相应 PDF、源文件和第三方依赖的许可与引用要求。正式引用请以 `docs/Foundation.pdf` 与 `docs/Application.pdf` 为准，并保留页面中的结果/回放边界说明。
