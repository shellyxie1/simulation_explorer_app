/* ===========================================================================
   Separation conditions: one table per explorer

   Shared by the linear and log-normal explorers, which mount it as their
   "Separation conditions" view. It reads data/meta.json and data/sepcond.json
   (export_sepcond.R) from the explorer's own folder: one value per replicate
   for each functional of the comparison theorem, in the same pct-major slot
   order as coefs.json. The thresholds are formed here from n, m and M'. The
   table has its own replicate slider, so the explorer's global knobs are
   hidden on this view.
   =========================================================================== */

(function () {
const d3 = window.d3;

const state = { mode: "rep", rep: {} };   // rep: replicate index per study (0-based)

/* ---- columns --------------------------------------------------------------
   `get` reads one derived record; `agg` says how 50 of them collapse in the
   all-replicates view: median of the numbers, a count for a yes/no. `kind`
   picks the number format. */
const COLS = [
  { group: "Scenario", head: "$\\epsilon$", csv: "epsilon", get: r => r.eps, kind: "eps", agg: "first" },
  { group: "Scenario", head: "$\\delta$", csv: "delta", get: r => r.shift, kind: "shift", agg: "first" },
  { group: "Scenario", head: "$m$", csv: "m", get: r => r.m, kind: "int", agg: "median" },

  { group: "$\\tilde{\\boldsymbol{\\theta}}$: displaced-start maximizer", head: "$M'(\\tilde{\\boldsymbol{\\theta}})$", csv: "M_prime", get: r => r.Mp, kind: "num", agg: "median" },
  { group: "$\\tilde{\\boldsymbol{\\theta}}$: displaced-start maximizer", head: "$\\dfrac{n\\log n + 2nM'}{m}$", csv: "threshold_m", get: r => r.thr, kind: "num", agg: "median" },

  { group: "$\\boldsymbol{\\theta}$: true-start maximizer", head: "$M^{-}_{\\mathcal{E}}$", csv: "Mneg_truestart", get: r => r.Mneg_t, kind: "num", agg: "median" },
  { group: "$\\boldsymbol{\\theta}$: true-start maximizer", head: "$M^{+}_{\\mathcal{E}}$", csv: "Mpos_truestart", get: r => r.Mpos_t, kind: "num", agg: "median" },
  { group: "$\\boldsymbol{\\theta}$: true-start maximizer", head: "holds", csv: "neg_holds_truestart", get: r => r.hold_t, kind: "bool", agg: "count" },
  { group: "$\\boldsymbol{\\theta}$: true-start maximizer", head: "group $k$", csv: "k_truestart", get: r => r.kneg_t > 0 ? r.kneg_t : null, kind: "int", agg: "kcount" },

  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "$M^{-}_{\\mathcal{E}}$", csv: "Mneg_truth", get: r => r.Mneg_0, kind: "num", agg: "median" },
  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "$M^{+}_{\\mathcal{E}}$", csv: "Mpos_truth", get: r => r.Mpos_0, kind: "num", agg: "median" },
  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "holds", csv: "neg_holds_truth", get: r => r.hold_0, kind: "bool", agg: "count" },
  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "group $k$", csv: "k_truth", get: r => r.kneg_0 > 0 ? r.kneg_0 : null, kind: "int", agg: "kcount" },
  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "$M^{-}_{\\mathcal{E}}$ of $k$", csv: "Mneg_group_truth", get: r => r.Mnegk_0, kind: "num", agg: "median" },
  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "$\\dfrac{n\\log n + 2nM'}{k}$", csv: "threshold_k_truth", get: r => r.thrk_0, kind: "num", agg: "median" },
  { group: "$\\boldsymbol{\\theta} = \\boldsymbol{\\theta}_0$", head: "outliers in $k$", csv: "outliers_in_group_truth", get: r => r.kout_0, kind: "int", agg: "median" },

  { group: "$\\log L_{\\mathrm{ETEL}}$", head: "$\\boldsymbol{\\theta}_0$", csv: "logL_truth", get: r => r.ll_0, kind: "ll", agg: "median" },
  { group: "$\\log L_{\\mathrm{ETEL}}$", head: "true start", csv: "logL_truestart", get: r => r.ll_t, kind: "ll", agg: "median" },
  { group: "$\\log L_{\\mathrm{ETEL}}$", head: "displaced start", csv: "logL_displaced", get: r => r.ll_d, kind: "ll", agg: "median" }
];

/* ---- data ----------------------------------------------------------------- */
const slotOf = st => (ip, ish, isim) => (ip * st.nShift + ish) * st.nSim + isim;

/* One replicate as the theorem sees it. The thresholds use the exported m and
   M'; `thrk_0` is the threshold for the best group of size k at theta_0. */
function record(st, ip, ish, isim) {
  const J = st.sc, s = slotOf(st)(ip, ish, isim), n = J.n;
  const m = J.m[s], Mp = J.Mp[s];
  if (m == null || Mp == null) return null;
  const a = n * Math.log(n) + 2 * n * Mp;
  const ok0 = J.conv_0[s] === 1;
  return {
    eps: +st.meta.pct[ip] / 100, shift: st.meta[st.shiftKey][ish], rep: st.meta.sim[isim],
    m, Mp, thr: a / m,
    Mneg_t: J.Mneg_t[s], Mpos_t: J.Mpos_t[s], hold_t: J.Mneg_t[s] > a / m, kneg_t: J.kneg_t[s], kpos_t: J.kpos_t[s],
    Mneg_0: ok0 ? J.Mneg_0[s] : null, Mpos_0: ok0 ? J.Mpos_0[s] : null,
    hold_0: ok0 ? J.Mneg_0[s] > a / m : null,
    kneg_0: ok0 ? J.kneg_0[s] : null, Mnegk_0: ok0 ? J.Mnegk_0[s] : null,
    thrk_0: ok0 && J.kneg_0[s] > 0 ? a / J.kneg_0[s] : null,
    kout_0: ok0 ? J.kout_0[s] : null, kpos_0: ok0 ? J.kpos_0[s] : null,
    ll_0: ok0 ? J.ll_0[s] : null, ll_t: J.ll_t[s], ll_d: J.ll_d[s]
  };
}

const finite = v => v != null && Number.isFinite(v);

/* ---- formatting ----------------------------------------------------------- */
const fInt = d3.format(",d"), fLL = d3.format(",.0f"), f0 = d3.format(",.0f"), f1 = d3.format(".1f"), f2 = d3.format(".2f");
function fmt(v, kind) {
  if (!finite(v)) return "–";
  switch (kind) {
    case "eps":   return f2(v);
    case "shift": return v > 0 ? `+${d3.format("~g")(v)}` : d3.format("~g")(v);
    case "int":   return fInt(Math.round(v));
    case "ll":    return fLL(v);
    default: {
      const a = Math.abs(v);
      return a >= 1000 ? f0(v) : a >= 100 ? d3.format(".0f")(v) : a >= 10 ? f1(v) : f2(v);
    }
  }
}

/* ---- cells ---------------------------------------------------------------- */
/* Returns {text, cls} for one column of one scenario row, in either view. */
function cell(col, recs) {
  if (state.mode === "rep") {
    const r = recs[0];
    if (!r) return { text: "–", cls: "na" };
    const v = col.get(r);
    if (col.kind === "bool") return v == null ? { text: "–", cls: "na" } : v ? { text: "yes", cls: "yes" } : { text: "no", cls: "no" };
    // A group size is itself a "yes": the condition holds for that group.
    if (col.agg === "kcount") return finite(v) ? { text: fmt(v, col.kind), cls: "yes" } : { text: "none", cls: "no" };
    return { text: fmt(v, col.kind), cls: finite(v) ? "" : "na" };
  }
  const vals = recs.filter(Boolean).map(col.get);
  const N = vals.length;
  if (col.agg === "first") return { text: fmt(vals[0], col.kind), cls: "" };
  if (col.agg === "count") {
    const known = vals.filter(v => v != null), c = known.filter(Boolean).length;
    return { text: known.length ? `${c}/${known.length}` : "–", cls: c > 0 ? "yes" : "no" };
  }
  if (col.agg === "kcount") {
    const c = vals.filter(finite).length;
    return { text: `${c}/${N}`, cls: c > 0 ? "yes" : "no" };
  }
  const med = d3.median(vals.filter(finite));
  return { text: fmt(med, col.kind), cls: finite(med) ? "" : "na" };
}

/* ---- build ---------------------------------------------------------------- */
function build(st) {
  const sec = d3.select(st.host);
  sec.html("");
  sec.append("h2").attr("class", "study-title")
    .html(`${st.title} <span class="dim">· $n = ${st.sc.n}$ · ${st.thetaTex}</span>`);

  const knobs = sec.append("div").attr("class", "knobs sep-knobs");
  const mc = knobs.append("span").attr("class", "ctl");
  mc.append("span").attr("class", "lbl").text("show");
  const modeTabs = mc.append("span").attr("class", "tabs mode-tabs");
  modeTabs.selectAll("button").data([{ key: "rep", label: "one replicate" }, { key: "all", label: `all ${st.nSim} replicates` }])
    .join("button").attr("class", "tab").attr("type", "button").text(d => d.label)
    .on("click", (ev, d) => { state.mode = d.key; renderStudy(st); });

  const rc = knobs.append("label").attr("class", "ctl rep-ctl");
  rc.append("span").attr("class", "lbl").text("replicate");
  rc.append("input").attr("type", "range").attr("min", 0).attr("max", st.nSim - 1).attr("step", 1)
    .property("value", state.rep[st.key])
    .on("input", function () { state.rep[st.key] = +this.value; renderStudy(st); });
  rc.append("span").attr("class", "val");

  const dc = knobs.append("span").attr("class", "ctl");
  dc.append("span").attr("class", "tabs").append("button").attr("class", "tab dl-csv").attr("type", "button").text("download csv")
    .on("click", () => downloadCsv(st));

  const wrap = sec.append("div").attr("class", "tablewrap");
  const table = wrap.append("table").attr("class", "sep");
  const thead = table.append("thead");
  const groups = [];
  COLS.forEach(c => { const g = groups[groups.length - 1]; if (g && g.name === c.group) g.span++; else groups.push({ name: c.group, span: 1 }); });
  thead.append("tr").attr("class", "groups").selectAll("th").data(groups).join("th")
    .attr("colspan", d => d.span).attr("class", (d, i) => i % 2 ? "alt" : "").html(d => d.name);
  thead.append("tr").attr("class", "heads").selectAll("th").data(COLS).join("th")
    .attr("class", d => (d.kind === "bool" ? "c" : "r") + (groupIndex(d) % 2 ? " alt" : "")).html(d => d.head);
  table.append("tbody");

  if (window.renderMathInElement) {
    renderMathInElement(sec.node(), {
      delimiters: [{ left: "$$", right: "$$", display: true }, { left: "$", right: "$", display: false }],
      throwOnError: false
    });
  }
}

function groupIndex(col) {
  let gi = -1, last = null;
  for (const c of COLS) { if (c.group !== last) { gi++; last = c.group; } if (c === col) return gi; }
  return gi;
}

/* ---- rows ----------------------------------------------------------------- */
function scenarioRows(st) {
  const rows = [];
  for (let ip = 0; ip < st.nPct; ip++) for (let ish = 0; ish < st.nShift; ish++) {
    const recs = state.mode === "rep"
      ? [record(st, ip, ish, state.rep[st.key])]
      : d3.range(st.nSim).map(isim => record(st, ip, ish, isim));
    rows.push({ ip, ish, recs });
  }
  return rows;
}

function renderStudy(st) {
  const sec = d3.select(st.host);
  sec.selectAll(".mode-tabs .tab").classed("on", d => d.key === state.mode);
  sec.select(".rep-ctl").style("display", state.mode === "rep" ? "" : "none");
  sec.select(".rep-ctl input").property("value", state.rep[st.key]);
  sec.select(".rep-ctl .val").text(`${st.meta.sim[state.rep[st.key]]} / ${st.nSim}`);

  const rows = scenarioRows(st);
  const tr = sec.select("tbody").selectAll("tr").data(rows).join("tr")
    .attr("class", d => d.ish === 0 ? "first" : "");
  tr.selectAll("td").data(d => COLS.map(c => ({ col: c, ...cell(c, d.recs) })))
    .join("td")
    .attr("class", d => `${d.col.kind === "bool" || d.col.agg === "count" || d.col.agg === "kcount" ? "c" : "r"} ${d.cls}${groupIndex(d.col) % 2 ? " alt" : ""}`)
    .text(d => d.text);
}

/* ---- csv ------------------------------------------------------------------ */
function downloadCsv(st) {
  const rows = scenarioRows(st);
  const head = ["study", "view", ...COLS.map(c => c.csv)];
  const view = state.mode === "rep" ? `replicate ${st.meta.sim[state.rep[st.key]]}` : `median or count over ${st.nSim} replicates`;
  const lines = [head.join(",")];
  for (const row of rows) {
    const vals = COLS.map(c => {
      const { text } = cell(c, row.recs);
      return text === "–" ? "" : text.replace(/,/g, "").replace(/^\+/, "");
    });
    lines.push([st.key, `"${view}"`, ...vals].join(","));
  }
  const blob = new Blob([lines.join("\n") + "\n"], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `separation_${st.key}_${state.mode === "rep" ? "rep" + st.meta.sim[state.rep[st.key]] : "all"}.csv`;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

/* ---- mount ----------------------------------------------------------------
   SepCond.show({ host, key, title, shiftKey, defaultRep, thetaTex }) loads
   the data on the first call and builds the table; later calls do nothing. */
let loading = null;
function show(st) {
  if (!loading) {
    loading = Promise.all([
      fetch("data/meta.json").then(r => r.json()),
      fetch("data/sepcond.json").then(r => { if (!r.ok) throw new Error(r.status); return r.json(); })
    ]).then(([meta, sc]) => {
      st.meta = meta; st.sc = sc;
      st.nPct = meta.pct.length;
      st.nShift = meta[st.shiftKey].length;
      st.nSim = meta.sim.length;
      const def = meta.sim.indexOf(st.defaultRep);
      state.rep[st.key] = def >= 0 ? def : 0;
    }, err => {
      d3.select(st.host).html(`<p class="caption">Could not load data/sepcond.json.</p>`);
      loading = null;   // try again next visit
      throw err;
    }).then(() => { build(st); renderStudy(st); });
  }
  return loading;
}

window.SepCond = { show };
})();
