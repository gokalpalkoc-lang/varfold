// Apply the saved theme before first paint to avoid a flash (kept external for a strict CSP).
try {
  var t = localStorage.getItem("varfold-theme");
  if (t === "light" || t === "dark") document.documentElement.dataset.theme = t;
} catch (e) {}
