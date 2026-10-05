(function () {
  "use strict";

  function renderMath() {
    if (!window.katex) {
      document.querySelectorAll("[data-math]").forEach(function (node) {
        if (node.dataset.mathFallback) node.textContent = node.dataset.mathFallback;
      });
      return;
    }

    document.querySelectorAll("[data-math]").forEach(function (node) {
      const source = node.dataset.math || "";
      const fallback = node.dataset.mathFallback || node.textContent || source;
      try {
        window.katex.render(source, node, {
          displayMode: node.classList.contains("math-display"),
          throwOnError: false,
          output: "htmlAndMathml",
          strict: "warn",
        });
        node.classList.add("math-rendered");
      } catch (error) {
        node.textContent = fallback;
        node.classList.add("math-render-fallback");
      }
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", renderMath, { once: true });
  } else {
    renderMath();
  }
})();
