(function () {
  document.documentElement.classList.add("js");
  const body = document.body;
  const languageButton = document.querySelector(".language-button");
  const languageKey = "qcd-inverse-language";

  function setLanguage(language) {
    body.dataset.language = language;
    document.documentElement.lang = language === "en" ? "en" : "zh-CN";
    document.title = language === "en" ? body.dataset.titleEn : body.dataset.titleZh;
    document.querySelectorAll("[data-label-zh][data-label-en]").forEach(function (node) {
      node.setAttribute("aria-label", language === "en" ? node.dataset.labelEn : node.dataset.labelZh);
    });
    if (languageButton) {
      languageButton.setAttribute("aria-pressed", language === "en" ? "true" : "false");
    }
    try {
      window.localStorage.setItem(languageKey, language);
    } catch (error) {
      // Storage is optional on static pages.
    }
  }

  if (languageButton) {
    languageButton.addEventListener("click", function () {
      setLanguage(body.dataset.language === "en" ? "zh" : "en");
    });
  }

  try {
    const storedLanguage = window.localStorage.getItem(languageKey);
    if (storedLanguage === "en" || storedLanguage === "zh") {
      setLanguage(storedLanguage);
    }
  } catch (error) {
    // Use Chinese as the default when storage is unavailable.
  }

  document.querySelectorAll('a[href="#top"]').forEach(function (link) {
    link.addEventListener("click", function (event) {
      event.preventDefault();
      const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
      window.scrollTo({ top: 0, behavior: reducedMotion ? "auto" : "smooth" });
      window.history.replaceState(null, "", window.location.pathname + window.location.search);
    });
  });

  document.querySelectorAll(".archive-summary[data-open-pdf]").forEach(function (summary) {
    summary.addEventListener("dblclick", function (event) {
      event.preventDefault();
      const pdfUrl = summary.dataset.openPdf;
      if (pdfUrl) {
        window.location.assign(pdfUrl);
      }
    });
  });

  const archiveEntries = Array.from(document.querySelectorAll(".archive-entry"));
  const reducedMotionQuery = window.matchMedia("(prefers-reduced-motion: reduce)");
  const archiveDuration = 380;

  function clearArchiveAnimationStyles(content) {
    content.style.removeProperty("height");
    content.style.removeProperty("opacity");
    content.style.removeProperty("transform");
  }

  function cancelArchiveAnimation(entry, content) {
    if (entry._archiveTimer) {
      window.clearTimeout(entry._archiveTimer);
      entry._archiveTimer = null;
    }
    if (entry._archiveFrame) {
      window.cancelAnimationFrame(entry._archiveFrame);
      entry._archiveFrame = null;
    }
    if (entry._archiveOnTransitionEnd) {
      content.removeEventListener("transitionend", entry._archiveOnTransitionEnd);
      entry._archiveOnTransitionEnd = null;
    }
    if (typeof content.getAnimations === "function") {
      content.getAnimations().forEach(function (animation) { animation.cancel(); });
    }
  }

  function parseTranslateY(transform) {
    if (!transform || transform === "none") return 0;
    const match = transform.match(/^matrix(?:3d)?\((.+)\)$/);
    if (!match) return 0;
    const values = match[1].split(",").map(Number);
    return transform.startsWith("matrix3d") ? (values[13] || 0) : (values[5] || 0);
  }

  function readArchiveVisualState(content) {
    const computed = window.getComputedStyle(content);
    const rect = content.getBoundingClientRect();
    return {
      height: rect.height,
      opacity: Number.parseFloat(computed.opacity) || 0,
      transform: parseTranslateY(computed.transform),
      paddingBottom: Number.parseFloat(computed.paddingBottom) || 0,
    };
  }

  function setArchiveExpanded(entry, expanded) {
    const summary = entry.querySelector(".archive-summary");
    if (summary) {
      summary.setAttribute("aria-expanded", expanded ? "true" : "false");
    }
  }

  function animateArchive(entry, expanded) {
    const content = entry.querySelector(".archive-content");
    if (!content) {
      entry.open = expanded;
      setArchiveExpanded(entry, expanded);
      return;
    }

    const current = readArchiveVisualState(content);
    const wasOpen = entry.open;
    cancelArchiveAnimation(entry, content);
    entry._archiveTarget = expanded;
    setArchiveExpanded(entry, expanded);

    if (reducedMotionQuery.matches) {
      entry.open = expanded;
      entry.removeAttribute("data-animating");
      clearArchiveAnimationStyles(content);
      return;
    }

    entry.setAttribute("data-animating", "");
    entry.open = true;
    const startHeight = wasOpen ? current.height : 0;
    const startOpacity = wasOpen ? current.opacity : 0;
    const startTransform = wasOpen ? current.transform : 8;
    const startPaddingBottom = wasOpen ? current.paddingBottom : 0;
    content.style.transition = "none";
    content.style.height = "auto";
    content.style.removeProperty("padding-bottom");
    const naturalHeight = content.getBoundingClientRect().height;
    const naturalPaddingBottom = window.getComputedStyle(content).paddingBottom;
    content.style.height = startHeight + "px";
    content.style.paddingBottom = startPaddingBottom + "px";
    content.style.opacity = String(startOpacity);
    content.style.transform = "translateY(" + startTransform + "px)";

    // Force the current pixel state to be painted before the next frame targets the new state.
    content.getBoundingClientRect();
    content.style.removeProperty("transition");
    entry._archiveNaturalHeight = naturalHeight;
    entry._archiveNaturalPaddingBottom = naturalPaddingBottom;

    function finish() {
      if (entry._archiveTarget !== expanded) return;
      if (entry._archiveTimer) {
        window.clearTimeout(entry._archiveTimer);
        entry._archiveTimer = null;
      }
      if (entry._archiveFrame) {
        window.cancelAnimationFrame(entry._archiveFrame);
        entry._archiveFrame = null;
      }
      if (entry._archiveOnTransitionEnd) {
        content.removeEventListener("transitionend", entry._archiveOnTransitionEnd);
        entry._archiveOnTransitionEnd = null;
      }

      content.style.height = expanded ? "auto" : "0px";
      content.style.paddingBottom = expanded ? "" : "0px";
      content.style.opacity = expanded ? "1" : "0";
      content.style.transform = expanded ? "translateY(0)" : "translateY(8px)";

      if (!expanded) {
        // The content is already at its closed visual state before details is hidden.
        entry.open = false;
        entry.removeAttribute("data-animating");
        clearArchiveAnimationStyles(content);
        return;
      }

      // Switch the completed open state to natural height without creating a second transition.
      content.style.transition = "none";
      content.getBoundingClientRect();
      content.style.removeProperty("height");
      content.style.removeProperty("opacity");
      content.style.removeProperty("transform");
      content.getBoundingClientRect();
      content.style.removeProperty("transition");
      entry.removeAttribute("data-animating");
    }

    function onTransitionEnd(event) {
      if (event.target === content && event.propertyName === "height") finish();
    }

    entry._archiveOnTransitionEnd = onTransitionEnd;
    content.addEventListener("transitionend", onTransitionEnd);
    entry._archiveFrame = window.requestAnimationFrame(function () {
      entry._archiveFrame = null;
      if (entry._archiveTarget !== expanded) return;
      if (expanded) {
        content.style.removeProperty("padding-bottom");
      }
      const targetHeight = expanded ? entry._archiveNaturalHeight : 0;
      content.style.paddingBottom = expanded ? entry._archiveNaturalPaddingBottom : "0px";
      content.style.height = targetHeight + "px";
      content.style.opacity = expanded ? "1" : "0";
      content.style.transform = expanded ? "translateY(0)" : "translateY(8px)";
    });
    entry._archiveTimer = window.setTimeout(finish, archiveDuration + 160);
  }

  archiveEntries.forEach(function (entry) {
    const summary = entry.querySelector(".archive-summary");
    setArchiveExpanded(entry, entry.open);
    if (!summary) return;
    summary.addEventListener("click", function (event) {
      event.preventDefault();
      const target = entry._archiveTarget === undefined ? entry.open : entry._archiveTarget;
      animateArchive(entry, !target);
    });
  });

  const revealItems = document.querySelectorAll("[data-reveal]");
  if ("IntersectionObserver" in window && !window.matchMedia("(prefers-reduced-motion: reduce)").matches) {
    const observer = new IntersectionObserver(function (entries, currentObserver) {
      entries.forEach(function (entry) {
        if (!entry.isIntersecting) {
          return;
        }
        entry.target.classList.add("is-visible");
        currentObserver.unobserve(entry.target);
      });
    }, { threshold: 0.14 });

    revealItems.forEach(function (item) {
      observer.observe(item);
    });
  } else {
    revealItems.forEach(function (item) {
      item.classList.add("is-visible");
    });
  }
})();
