/* scStudio client-side helpers: language switch + startup animation */
(function () {
  "use strict";

  // ---------------------------------------------------------------------------
  // Language switch — only ONE language visible at a time.
  // Each translatable element has class "i18n" and data-en / data-zh holding the
  // two versions (may be HTML). We swap innerHTML to the active language.
  // ---------------------------------------------------------------------------
  window.scStudioLang = "en";

  function applyLang(root) {
    var lang = window.scStudioLang;
    var scope = root || document;
    var nodes = scope.querySelectorAll(".i18n");
    nodes.forEach(function (el) {
      var val = el.getAttribute("data-" + lang);
      if (val === null) val = el.getAttribute("data-en");
      if (val !== null && el.innerHTML !== val) el.innerHTML = val;
    });
  }

  window.scStudioSetLang = function (lang) {
    window.scStudioLang = (lang === "zh") ? "zh" : "en";
    document.body.setAttribute("data-lang", window.scStudioLang);
    applyLang(document);
    // reflect active state on the segmented control
    document.querySelectorAll(".scstudio-lang-btn").forEach(function (b) {
      b.classList.toggle("active", b.getAttribute("data-lang") === window.scStudioLang);
    });
  };

  // Re-apply language to any freshly rendered (renderUI) content.
  document.addEventListener("shiny:idle", function () { applyLang(document); });
  document.addEventListener("shiny:value", function (e) {
    if (e.target) setTimeout(function () { applyLang(e.target); }, 0);
  });

  // ---------------------------------------------------------------------------
  // Startup animation — scattered "cells" fly in and coalesce into coloured
  // UMAP-like clusters, then the logo fades in; the splash fades out on load.
  // ---------------------------------------------------------------------------
  function runSplash() {
    var splash = document.getElementById("scstudio-splash");
    if (!splash) return;
    var canvas = document.getElementById("scstudio-splash-canvas");
    if (!canvas || !canvas.getContext) { fadeOut(splash); return; }
    var ctx = canvas.getContext("2d");
    var W, H;
    function resize() {
      W = canvas.width = splash.clientWidth;
      H = canvas.height = splash.clientHeight;
    }
    resize();
    window.addEventListener("resize", resize);

    // brand palette (blue-cyan research theme + accents)
    var pal = ["#2f81c7", "#4f9fd0", "#3fb37f", "#e4572e", "#b5179e",
               "#f4a261", "#9c6ade", "#2a9d8f"];
    var K = 6;                     // number of UMAP-like clusters
    var N = 340;                   // number of cells
    var cx = W / 2, cy = H / 2 - 20;
    var R = Math.min(W, H) * 0.26;
    var centers = [];
    for (var k = 0; k < K; k++) {
      var a = (k / K) * Math.PI * 2;
      centers.push({ x: cx + Math.cos(a) * R, y: cy + Math.sin(a) * R,
                     col: pal[k % pal.length] });
    }
    var cells = [];
    for (var i = 0; i < N; i++) {
      var c = centers[i % K];
      cells.push({
        sx: Math.random() * W, sy: Math.random() * H,             // start (scattered)
        tx: c.x + (Math.random() - 0.5) * R * 0.7,                // target (in cluster)
        ty: c.y + (Math.random() - 0.5) * R * 0.7,
        col: c.col, r: 1.5 + Math.random() * 2.2
      });
    }

    var start = null, DUR = 1600;
    function ease(t) { return 1 - Math.pow(1 - t, 3); }          // easeOutCubic
    function frame(ts) {
      if (start === null) start = ts;
      var p = Math.min(1, (ts - start) / DUR);
      var e = ease(p);
      ctx.clearRect(0, 0, W, H);
      for (var i = 0; i < cells.length; i++) {
        var c = cells[i];
        var x = c.sx + (c.tx - c.sx) * e;
        var y = c.sy + (c.ty - c.sy) * e;
        ctx.globalAlpha = 0.35 + 0.55 * e;
        ctx.fillStyle = c.col;
        ctx.beginPath();
        ctx.arc(x, y, c.r, 0, Math.PI * 2);
        ctx.fill();
      }
      ctx.globalAlpha = 1;
      if (p < 1) {
        requestAnimationFrame(frame);
      } else {
        var logo = document.getElementById("scstudio-splash-logo");
        if (logo) logo.classList.add("show");
        setTimeout(function () { fadeOut(splash); }, 750);
      }
    }
    requestAnimationFrame(frame);
    // safety: never let the splash trap the user
    setTimeout(function () { fadeOut(splash); }, 5000);
  }

  function fadeOut(splash) {
    if (!splash || splash.classList.contains("hide")) return;
    splash.classList.add("hide");
    setTimeout(function () { if (splash && splash.parentNode) splash.style.display = "none"; }, 600);
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", runSplash);
  } else {
    runSplash();
  }
})();
