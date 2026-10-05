const createElement = (iconNode, customAttrs = {}) => {
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  const attributes = {
    xmlns: "http://www.w3.org/2000/svg",
    width: 24,
    height: 24,
    viewBox: "0 0 24 24",
    fill: "none",
    stroke: "currentColor",
    "stroke-width": 2,
    "stroke-linecap": "round",
    "stroke-linejoin": "round",
    ...customAttrs,
  };
  Object.entries(attributes).forEach(([name, value]) => svg.setAttribute(name, String(value)));
  iconNode.forEach(([tag, attrs]) => {
    const child = document.createElementNS("http://www.w3.org/2000/svg", tag);
    Object.entries(attrs).forEach(([name, value]) => child.setAttribute(name, String(value)));
    svg.appendChild(child);
  });
  return svg;
};
const icons = {
  download: [["path", { d: "M12 15V3" }], ["path", { d: "M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4" }], ["path", { d: "m7 10 5 5 5-5" }]],
  copy: [["rect", { width: "14", height: "14", x: "8", y: "8", rx: "2" }], ["path", { d: "M4 16c-1.1 0-2-.9-2-2V4c0-1.1.9-2 2-2h10c1.1 0 2 .9 2 2" }]],
  "file-code": [["path", { d: "M6 22a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h8a2.4 2.4 0 0 1 1.704.706l3.588 3.588A2.4 2.4 0 0 1 20 8v12a2 2 0 0 1-2 2z" }], ["path", { d: "M14 2v5a1 1 0 0 0 1 1h5" }], ["path", { d: "M10 12.5 8 15l2 2.5" }], ["path", { d: "m14 12.5 2 2.5-2 2.5" }]],
  chart: [["path", { d: "M3 3v16a2 2 0 0 0 2 2h16" }], ["path", { d: "M7 16c.5-2 1.5-7 4-7 2 0 2 3 4 3 2.5 0 4.5-5 5-7" }]],
  reset: [["path", { d: "M3 12a9 9 0 1 0 9-9 9.75 9.75 0 0 0-6.74 2.74L3 8" }], ["path", { d: "M3 3v5h5" }]],
  image: [["rect", { width: "18", height: "18", x: "3", y: "3", rx: "2" }], ["circle", { cx: "9", cy: "9", r: "2" }], ["path", { d: "m21 15-3.086-3.086a2 2 0 0 0-2.828 0L6 21" }]],
  maximize: [["path", { d: "M8 3H5a2 2 0 0 0-2 2v3" }], ["path", { d: "M21 8V5a2 2 0 0 0-2-2h-3" }], ["path", { d: "M3 16v3a2 2 0 0 0 2 2h3" }], ["path", { d: "M16 21h3a2 2 0 0 0 2-2v-3" }]],
};

(async () => {
document.querySelectorAll("[data-icon]").forEach((node) => {
  const icon = icons[node.dataset.icon];
  if (icon) node.replaceWith(createElement(icon, { width: 16, height: 16, "aria-hidden": "true", "stroke-width": 1.7 }));
});

const languageButton = document.querySelector(".language-button");
function setLanguage(language) {
  document.body.dataset.language = language;
  document.documentElement.lang = language === "en" ? "en" : "zh-CN";
  document.title = language === "en" ? document.body.dataset.titleEn : document.body.dataset.titleZh;
  document.querySelectorAll("[data-label-zh][data-label-en]").forEach((node) => {
    node.setAttribute("aria-label", language === "en" ? node.dataset.labelEn : node.dataset.labelZh);
  });
  document.querySelectorAll("[data-title-zh][data-title-en]").forEach((node) => {
    node.setAttribute("title", language === "en" ? node.dataset.titleEn : node.dataset.titleZh);
  });
  languageButton?.setAttribute("aria-pressed", language === "en" ? "true" : "false");
  try { window.localStorage.setItem("qcd-inverse-language", language); } catch {}
}
languageButton?.addEventListener("click", () => setLanguage(document.body.dataset.language === "en" ? "zh" : "en"));
try {
  const storedLanguage = window.localStorage.getItem("qcd-inverse-language");
  if (storedLanguage === "en" || storedLanguage === "zh") setLanguage(storedLanguage);
} catch {}

const projects = {
  "eta-c": {
    directory: "eta_c",
    original: "ηc",
    title: "ηc · A–A n=0",
    poles: ["ηc", "χc1"],
    files: [
      "eta_c_inverse.m", "eta_c_operator.m", "eta_c_spectrum.m",
      "test_eta_c_inverse.m", "solver.py", "physics.py", "factors.py",
      "gsvd_backend.py", "alpha_grid.py", "raus2024.py", "run.py",
      "README.md", "README_Python.md",
    ],
  },
  jpsi: {
    directory: "jpsi",
    original: "Jpsi",
    title: "J/ψ · n=0",
    poles: ["J/ψ", "ψ(2S)"],
    files: [
      "solver.py", "physics.py", "factors.py", "gsvd_backend.py",
      "alpha_grid.py", "raus2024.py", "run.py", "README.md",
    ],
  },
};

const $ = (id) => document.getElementById(id);
const state = { project: "eta-c", file: "eta_c_inverse.m", source: "matlab", model: "1delta", view: "spectrum", alphaSlot: 4 };
const sourceCache = new Map();
const dataCache = new Map();
const embeddedLabData = window.QCD_LAB_DATA || null;
let data = null;
let sourceText = "";
let sourceRevision = 0;
let resultRevision = 0;
let plotQueue = Promise.resolve();
let noticeTimer;
let alphaOrder = [];
let chartResizeObserver;

const translate = (zh, en) => document.body.dataset.language === "en" ? en : zh;
const scientific = (value) => Number.isFinite(value) ? value.toExponential(3) : "NaN";
const fixed = (value, digits = 6) => Number.isFinite(value) ? value.toFixed(digits) : "NaN";
const escape = (value) => String(value).replace(/[&<>"']/g, (character) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[character]);

function tokenMarkup(token) {
  if (typeof token === "string") return escape(token);
  if (Array.isArray(token)) return token.map(tokenMarkup).join("");
  const classes = ["token", token.type, ...[token.alias || []].flat()];
  // Preserve whole-file grammar context, but close every token span at a line boundary.
  return tokenMarkup(token.content).split(/\r?\n/).map((part) =>
    `<span class="${escape(classes.join(" "))}">${part}</span>`).join("\n");
}

function notify(message) {
  $("lab-notice").textContent = message;
  clearTimeout(noticeTimer);
  noticeTimer = setTimeout(() => { $("lab-notice").textContent = ""; }, 2400);
}

function fail(error) {
  $("lab-error-text").textContent = translate("资源读取失败：", "Resource loading failed: ") + error.message;
  $("lab-error").hidden = false;
}

function chartStatus(message, kind = "loading") {
  const node = $("lab-chart-status");
  if (!node) return;
  node.textContent = message || "";
  node.dataset.kind = kind;
  node.hidden = !message;
}

function clearChart() {
  const chart = $("lab-chart");
  if (!chart) return;
  if (window.Plotly && (chart.data || chart._fullLayout)) {
    try { window.Plotly.purge(chart); } catch {}
  }
  chart.replaceChildren();
  chart.removeAttribute("data-fallback");
}

function numericSeries(x, y, label) {
  const xValues = Array.from(x ?? [], Number);
  const yValues = Array.from(y ?? [], Number);
  if (xValues.length !== yValues.length || xValues.length < 2) {
    throw new Error(`${label}: invalid series length`);
  }
  const values = { x: xValues, y: yValues };
  if (!values.x.every(Number.isFinite) || !values.y.every(Number.isFinite)) {
    throw new Error(`${label}: non-finite data`);
  }
  return values;
}

function renderSpectrumFallback(traces, error) {
  const chart = $("lab-chart");
  const width = Math.max(chart.clientWidth || 720, 320);
  const height = 326;
  const margin = { left: 58, right: 18, top: 24, bottom: 42 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const points = traces.flatMap((trace) => trace.x.map((x, i) => ({ x, y: trace.y[i] })));
  if (!points.length || points.some((point) => !Number.isFinite(point.x) || !Number.isFinite(point.y))) {
    chartStatus(translate('谱函数数据无效。', 'Spectrum data are invalid.'), 'error');
    return;
  }
  const xMin = Math.min(...points.map((point) => point.x));
  const xMax = Math.max(...points.map((point) => point.x));
  const rawYMin = Math.min(...points.map((point) => point.y));
  const rawYMax = Math.max(...points.map((point) => point.y));
  const yMin = Math.min(0, rawYMin);
  const yMax = Math.max(0, rawYMax);
  const xSpan = Math.max(xMax - xMin, Number.EPSILON);
  const ySpan = Math.max(yMax - yMin, Number.EPSILON);
  const x = (value) => margin.left + (value - xMin) / xSpan * plotWidth;
  const y = (value) => margin.top + (yMax - value) / ySpan * plotHeight;
  const polyline = (trace) => trace.x.map((value, i) => `${x(value).toFixed(2)},${y(trace.y[i]).toFixed(2)}`).join(" ");
  const colors = ["#25506b", "#c79a42"];
  const lines = traces.map((trace, i) =>
    `<polyline points="${polyline(trace)}" fill="none" stroke="${colors[i] || "#55788a"}" stroke-width="${i ? 1.5 : 2.2}"${i ? ' stroke-dasharray="6 4"' : ""}/>`
  ).join("");
  const zero = y(0);
  chart.innerHTML = `<svg class="lab-spectrum-fallback" viewBox="0 0 ${width} ${height}" role="img" aria-label="${escape(translate("谱函数静态回放图", "Static spectrum replay"))}">
    <rect width="${width}" height="${height}" fill="#fff"/>
    <line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${height - margin.bottom}" stroke="#b1c4bd"/>
    <line x1="${margin.left}" y1="${height - margin.bottom}" x2="${width - margin.right}" y2="${height - margin.bottom}" stroke="#b1c4bd"/>
    <line x1="${margin.left}" y1="${zero}" x2="${width - margin.right}" y2="${zero}" stroke="#d7e2de"/>
    <text x="${margin.left}" y="${height - 12}" fill="#61716f" font-size="10">s [GeV²]</text>
    <text x="12" y="${margin.top + 4}" fill="#61716f" font-size="10">ρc(s)</text>
    ${lines}
    <g fill="#61716f" font-size="10"><text x="${margin.left}" y="${height - margin.bottom + 17}">${xMin.toFixed(2)}</text><text x="${width - margin.right - 28}" y="${height - margin.bottom + 17}">${xMax.toFixed(2)}</text><text x="${margin.left - 42}" y="${margin.top + 4}">${yMax.toExponential(2)}</text><text x="${margin.left - 42}" y="${height - margin.bottom}">${yMin.toExponential(2)}</text></g>
  </svg>`;
  chart.dataset.fallback = "true";
  chartStatus(translate(`交互图加载失败，已显示静态谱函数。${error ? `原因：${error.message}` : ""}`, `Interactive plot unavailable; showing a static spectrum. ${error ? `Reason: ${error.message}` : ""}`), "error");
}

async function resource(cache, path, json = false) {
  if (!cache.has(path)) {
    cache.set(path, fetch(path).then((response) => {
      if (!response.ok) throw new Error(`${response.status}: ${path}`);
      return json ? response.json() : response.text();
    }).catch((error) => { cache.delete(path); throw error; }));
  }
  return cache.get(path);
}

function populateFiles() {
  $("lab-file").replaceChildren(...projects[state.project].files.map((filename) => {
    const option = document.createElement("option");
    option.value = filename;
    option.textContent = filename;
    return option;
  }));
  $("lab-file").value = state.file;
  $("lab-file").disabled = false;
}

async function readSource() {
  const revision = ++sourceRevision;
  const project = projects[state.project];
  const file = state.file;
  const path = `./code/${project.directory}/${file}`;
  $("lab-file-name").textContent = file;
  $("lab-language").textContent = file.endsWith(".m") ? "MATLAB" : file.endsWith(".py") ? "PYTHON" : "README";
  $("lab-source-path").textContent = `${project.original}/${file}`;
  $("lab-code").textContent = translate("正在读取归档文件…", "Loading archive source...");
  try {
    const embeddedKey = `${project.directory}/${file}`;
    const text = embeddedLabData?.sources?.[embeddedKey] ?? await resource(sourceCache, path);
    if (typeof text !== "string") throw new Error(`Missing source: ${embeddedKey}`);
    if (revision !== sourceRevision) return;
    sourceText = text;
    const language = file.endsWith(".m") ? "matlab" : file.endsWith(".py") ? "python" : null;
    const highlighted = language && window.Prism
      ? tokenMarkup(window.Prism.tokenize(text, window.Prism.languages[language]))
      : escape(text);
    $("lab-code").innerHTML = highlighted.split(/\r?\n/).map((line) => `<span class="lab-code-line">${line || " "}</span>`).join("");
    $("lab-line-count").textContent = `${text.trimEnd().split(/\r?\n/).length} lines · UTF-8`;
    document.querySelector(".lab-code-viewport").scrollTo(0, 0);
  } catch (error) {
    if (revision === sourceRevision) {
      $("lab-code").textContent = translate("源文件读取失败。", "Source file could not be loaded.");
      fail(error);
    }
  }
}

function setupSource(data) {
  $("lab-source-select").replaceChildren();
  if (data.matlab) {
    const matlab = new Option("MATLAB · reference_matlab.npz", "matlab");
    $("lab-source-select").add(matlab);
    alphaOrder = data.matlab.alpha.map((alpha, index) => ({ alpha, index })).sort((a, b) => a.alpha - b.alpha);
    $("lab-alpha-range").max = String(alphaOrder.length - 1);
    const selected = data.channels[state.model].alpha;
    state.alphaSlot = alphaOrder.reduce((best, item, i) =>
      Math.abs(Math.log(item.alpha / selected)) < Math.abs(Math.log(alphaOrder[best].alpha / selected)) ? i : best, 0);
  } else {
    state.source = "python";
    alphaOrder = [];
  }
  $("lab-source-select").add(new Option("Python · result.npz", "python"));
  $("lab-source-select").value = state.source;
  $("lab-source-select").disabled = !data.matlab;
  $("lab-alpha-range").value = String(state.alphaSlot);
}

function currentResult() {
  const original = data.channels[state.model];
  if (state.source !== "matlab") return { ...original, s: original.s, M2: data.M2, OPE: data.OPE, source: data.source, sha256: data.sha256 };
  const i = alphaOrder[state.alphaSlot].index;
  const ref = data.matlab.channels[state.model];
  return {
    masses: original.masses, alpha: data.matlab.alpha[i],
    residues: ref.residues[i], f: ref.f[i],
    rho: ref.rho[i], s: data.matlab.s,
    relative_residual: ref.relative_residual[i],
    negative_count: ref.negative_count[i], penalty: ref.penalty[i],
    M2: data.matlab.M2, OPE: data.matlab.OPE,
    source: data.matlab.source, sha256: data.matlab.sha256,
  };
}

function layout(xTitle, yTitle) {
  return {
    paper_bgcolor: "#ffffff", plot_bgcolor: "#ffffff",
    font: { family: "Segoe UI, Microsoft YaHei, sans-serif", size: 11, color: "#61716f" },
    margin: { l: 58, r: 20, t: 42, b: 45 },
    xaxis: { title: { text: xTitle, standoff: 12 }, gridcolor: "#e8efec", zerolinecolor: "#b1c4bd", tickfont: { size: 10 }, automargin: false },
    yaxis: { title: { text: yTitle, standoff: 8 }, gridcolor: "#e8efec", zerolinecolor: "#b1c4bd", tickfont: { size: 10 }, automargin: true },
    legend: { orientation: "h", y: 1.15, x: 0, font: { size: 10 } },
    hovermode: "x unified",
    showlegend: true,
  };
}

function line(x, y, name, color, options = {}) {
  return { type: "scatter", mode: "lines", x, y, name, line: { color, width: 2 }, connectgaps: false, ...options };
}

function plotContent(result) {
  let config;
  let traces;
  let caption;
  if (state.view === "spectrum") {
    config = layout("s [GeV²]", "ρc(s) [dimensionless]");
    const spectrum = numericSeries(result.s, result.rho, "spectrum");
    traces = [line(spectrum.x, spectrum.y, "ρc(s)", "#25506b")];
    const m = data.metadata.mc;
    const correction = state.project === "jpsi" ? 1 + data.metadata.alpha_s / Math.PI : 1;
    const perturbative = spectrum.x.map((s) => {
      const v = Math.sqrt(Math.max(0, 1 - 4 * m * m / s));
      return 3 * v * (1 - v * v / 3) / (8 * Math.PI ** 2) * correction;
    });
    traces.push(line(spectrum.x, perturbative, state.project === "jpsi" ? "ρpert · LO+NLO" : "ρLO", "#c79a42", { line: { color: "#c79a42", width: 1.5, dash: "dash" } }));
    caption = translate("连续谱保留全部节点与负值；δ 极点的质量和留数单独列出。", "All continuum nodes and negative values are retained; δ-pole masses and residues are listed separately.");
  } else if (state.view === "ope") {
    config = layout("M² [GeV²]", "g(M²) [dimensionless]");
    traces = [
      line(result.M2, result.OPE.total, "OPE total", "#25506b"),
      line(result.M2, result.OPE.pert, "Perturbative", "#c79a42", { line: { color: "#c79a42", width: 1.5, dash: "dash" } }),
      line(result.M2, result.OPE.D4, "G2 / D4", "#55788a"),
      line(result.M2, result.OPE.D6, "G3 / D6", "#846470"),
    ];
    if (state.source === "python") {
      traces.push(line(result.M2, result.backsub, "Back-substitution", "#173d54", { line: { color: "#173d54", width: 1, dash: "dot" } }));
    }
    caption = translate(`OPE 绘图抽取 ${result.M2.length}/${data.n_borel} 个 Borel 节点；残差使用原始全网格。`, `OPE plot: ${result.M2.length}/${data.n_borel} Borel nodes; residuals use the original full grid.`);
  } else if (state.source === "matlab") {
    config = layout("α", "f [GeV]");
    config.xaxis.type = "log";
    const ref = data.matlab.channels[state.model];
    const alpha = alphaOrder.map((item) => item.alpha);
    traces = result.masses.map((mass, i) => line(alpha, alphaOrder.map((item) => ref.f[item.index][i]), projects[state.project].poles[i], i ? "#c79a42" : "#25506b", { mode: "markers", marker: { size: 6 } }));
    const lower = fixed(result.f[0]);
    caption = translate(`MATLAB 存储的 8 个 α 基准点，未做插值或自动选参；非正留数的 f 为 NaN。当前 f₁ = ${lower}。`, `Eight stored MATLAB α values; no interpolation or automatic selection. Nonpositive residues produce NaN decay constants. Current f₁ = ${lower}.`);
  } else {
    const diagnostic = data.channels[state.model].diagnostics;
    config = layout("α", "ψQ(α)");
    config.xaxis.type = "log";
    config.yaxis.type = "log";
    traces = [
      line(diagnostic.alpha, diagnostic.psi_Q.map((v) => v > 0 ? v : null), "ψQ", "#25506b"),
      { type: "scatter", mode: "markers", x: [result.alpha], y: [diagnostic.selected_psi_Q], name: translate("保存选参", "Saved selection"), marker: { color: "#c79a42", size: 9 } },
    ];
    caption = translate(`Python 保存的选参诊断，原 α 网格 ${diagnostic.n_alpha} 点；绘图保留局部极小值与选中点。`, `Saved Python selection diagnostic; original α grid: ${diagnostic.n_alpha} points. Local minima and the selected point are preserved in the plot.`);
  }
  return { config, traces, caption };
}

function updateFacts(result) {
  $("lab-f").textContent = fixed(result.f[0]);
  $("lab-residual").textContent = scientific(result.relative_residual);
  $("lab-negative").textContent = result.negative_count;
  $("lab-negative").classList.toggle("is-negative", result.negative_count > 0);
  $("lab-model-label").textContent = projects[state.project].poles.slice(0, result.masses.length).join(" + ");
  $("lab-result-path").textContent = result.source;
  $("lab-alpha-control").hidden = state.source !== "matlab";
  $("lab-alpha-output").textContent = scientific(result.alpha);
  $("lab-alpha-range").value = state.alphaSlot;
  $("lab-poles").innerHTML = `<table><thead><tr><th>${translate("δ 极点", "δ pole")}</th><th>m [GeV]</th><th>R [GeV²]</th><th>f [GeV]</th></tr></thead><tbody>${result.masses.map((mass, i) => `<tr><td>${escape(projects[state.project].poles[i])}</td><td>${fixed(mass, 5)}</td><td>${fixed(result.residues[i])}</td><td>${fixed(result.f[i])}</td></tr>`).join("")}</tbody></table>`;
  const facts = [
    [translate("来源", "Source"), state.source === "matlab" ? "MATLAB R2026a reference fixture" : "Python saved solver result"],
    ["M² / Λ [GeV²]", `[${data.metadata.M2_min}, ${data.metadata.M2_max}] / ${data.metadata.Lambda}`],
    ["α", scientific(result.alpha)],
    [translate("α 选择", "α selection"), state.source === "matlab" ? translate("未自动选参", "No automatic selection") : data.metadata.rule],
    ["Δs / ΔM² [GeV²]", `${data.metadata.ds} / ${data.metadata.dM2}`],
    [translate("谱节点 / Borel 节点", "Spectral / Borel nodes"), `${data.metadata.n_hats} / ${data.n_borel}`],
    [translate("相对加权残差", "Relative weighted residual"), scientific(result.relative_residual)],
    ["H¹ penalty", scientific(result.penalty)],
    [translate("最小连续谱", "Minimum continuum"), scientific(Math.min(...result.rho))],
    ["ρ < 0", String(result.negative_count)],
    ["SHA-256", result.sha256],
  ];
  $("lab-summary-panel").innerHTML = `<dl>${facts.map(([key, value]) => `<dt>${escape(key)}</dt><dd>${escape(value)}</dd>`).join("")}</dl>`;
  $("lab-provenance").textContent = state.source === "matlab"
    ? translate("来源：归档中的 MATLAB R2026a 数值夹具，固定窗口 M²∈[5,80]、Λ=30 GeV²。浏览器只回放存储结果，不执行 MATLAB。", "Source: archived MATLAB R2026a fixture at M²∈[5,80], Λ=30 GeV². The browser replays stored results; it does not execute MATLAB.")
    : translate("来源：归档中 Python 求解器的 result.npz（2026-09-11）；只展示保存的窗口和选参结果，未在网页重新计算。", "Source: the archived Python solver result.npz (2026-09-11). Only the stored window and selection are shown; the browser does not recompute them.");
  $("lab-limit-note").textContent = translate(
    `数值示例不等同于论文结论；未施加正性约束，极小 α 可能放大舍入误差。${state.project === "eta-c" ? "双极点模型的第二态为 χc1，不是 ηc(2S)。" : "双极点模型的第二态为 ψ(2S)。"}`,
    `Numerical examples are not paper conclusions. Positivity is not imposed; tiny α may amplify roundoff. ${state.project === "eta-c" ? "The second pole is χc1, not ηc(2S)." : "The second pole is ψ(2S)."}`
  );
}

async function renderResults() {
  if (!data) return;
  const revision = ++resultRevision;
  const result = currentResult();
  updateFacts(result);
  const summary = state.view === "summary";
  $("lab-plot-panel").hidden = summary;
  $("lab-summary-panel").hidden = !summary;
  $("lab-reset").disabled = summary;
  $("lab-plot-panel").setAttribute("aria-labelledby", `tab-${state.view}`);
  document.querySelectorAll("[data-lab-view]").forEach((tab) => {
    const active = tab.dataset.labView === state.view;
    tab.setAttribute("aria-selected", active);
    tab.tabIndex = active ? 0 : -1;
  });
  document.querySelectorAll("[data-lab-model]").forEach((button) => button.setAttribute("aria-pressed", button.dataset.labModel === state.model));
  if (summary) return;
  clearChart();
  chartStatus(translate("正在加载谱函数图…", "Loading spectrum plot..."));
  let plot;
  try {
    plot = plotContent(result);
  } catch (error) {
    if (state.view === "spectrum") chartStatus(translate("谱函数数据无效。", "Spectrum data are invalid."), "error");
    fail(error);
    return;
  }
  $("lab-plot-caption").textContent = plot.caption;
  if (!window.Plotly || typeof window.Plotly.react !== "function") {
    if (state.view === "spectrum") renderSpectrumFallback(plot.traces);
    else chartStatus(translate("交互绘图库未加载。", "Interactive plotting library unavailable."), "error");
    return;
  }
  plotQueue = plotQueue.catch(() => {}).then(async () => {
    if (revision !== resultRevision) return;
    await new Promise((resolve) => requestAnimationFrame(resolve));
    if (!$("lab-chart").clientWidth) throw new Error("chart container has zero width");
    await window.Plotly.react($("lab-chart"), plot.traces, plot.config, {
      responsive: true, displaylogo: false, displayModeBar: false, scrollZoom: false,
    });
    await new Promise((resolve) => requestAnimationFrame(resolve));
    if (state.view === "spectrum" && !$("lab-chart").querySelector(".main-svg")) {
      throw new Error("Plotly returned without a visible SVG");
    }
    window.Plotly.Plots.resize($("lab-chart"));
    chartStatus("");
  });
  try {
    await plotQueue;
  } catch (error) {
    if (state.view === "spectrum") renderSpectrumFallback(plot.traces, error);
    else chartStatus(translate("图表绘制失败。", "The plot could not be rendered."), "error");
  }
}

function observeChartSize() {
  const chart = $("lab-chart");
  if (!chart || !window.ResizeObserver) return;
  chartResizeObserver?.disconnect();
  chartResizeObserver = new ResizeObserver(() => {
    if (chart.data?.length && window.Plotly) window.Plotly.Plots.resize(chart);
  });
  chartResizeObserver.observe(chart);
}

let projectRevision = 0;
async function loadProject() {
  const revision = ++projectRevision;
  const project = state.project;
  document.querySelector(".lab-workspace").setAttribute("aria-busy", "true");
  $("lab-error").hidden = true;
  data = null;
  populateFiles();
  const sourcePromise = readSource();
  try {
    const loaded = embeddedLabData?.results?.[project]
      ?? await resource(dataCache, `./code/results/${project}-preview.json`, true);
    if (revision !== projectRevision) return;
    data = loaded;
    setupSource(data);
    await renderResults();
  } catch (error) { if (revision === projectRevision) fail(error); }
  await sourcePromise;
  if (revision === projectRevision) document.querySelector(".lab-workspace").setAttribute("aria-busy", "false");
}

$("lab-project").addEventListener("change", () => {
  ++resultRevision;
  state.project = $("lab-project").value;
  state.file = projects[state.project].files[0];
  state.model = "1delta";
  state.source = state.project === "eta-c" ? "matlab" : "python";
  loadProject();
});
$("lab-file").addEventListener("change", () => { state.file = $("lab-file").value; readSource(); });
$("lab-source-select").addEventListener("change", () => { state.source = $("lab-source-select").value; renderResults(); });
$("lab-alpha-range").addEventListener("input", () => { state.alphaSlot = Number($("lab-alpha-range").value); renderResults(); });
document.querySelectorAll("[data-lab-model]").forEach((button) => button.addEventListener("click", () => { state.model = button.dataset.labModel; renderResults(); }));
const tabs = [...document.querySelectorAll("[data-lab-view]")];
tabs.forEach((tab, i) => {
  tab.addEventListener("click", () => { state.view = tab.dataset.labView; renderResults(); });
  tab.addEventListener("keydown", (event) => {
    let index;
    if (event.key === "ArrowRight") index = (i + 1) % tabs.length;
    if (event.key === "ArrowLeft") index = (i - 1 + tabs.length) % tabs.length;
    if (event.key === "Home") index = 0;
    if (event.key === "End") index = tabs.length - 1;
    if (index !== undefined) { event.preventDefault(); tabs[index].focus(); tabs[index].click(); }
  });
});
$("lab-wrap").addEventListener("change", () => document.querySelector(".lab-code-viewport").classList.toggle("is-wrapped", $("lab-wrap").checked));
$("lab-reset").addEventListener("click", () => window.Plotly?.relayout($("lab-chart"), { "xaxis.autorange": true, "yaxis.autorange": true }));
$("lab-retry").addEventListener("click", loadProject);
new MutationObserver(() => { if (data) renderResults(); }).observe(document.body, { attributes: true, attributeFilter: ["data-language"] });
window.addEventListener("resize", () => {
  const chart = $("lab-chart");
  if (chart?.data?.length && window.Plotly) window.Plotly.Plots.resize(chart);
});
observeChartSize();
await loadProject();
})();
