/* 知识图谱渲染：四个研究方向各一页（页签切换），每页只画该方向的问题、方案与知识闭包。
 * 两种排序：
 *   自上而下（down）：❓控制问题 → 🛠解决方案 → 知识层 ⑦→①（箭头 = 拆解/依赖）
 *   自下而上（up）  ：知识层 ①→⑦ → 🛠解决方案 → ❓控制问题（箭头 = 支撑/通往）
 */
window.Graph = (function () {
  const SVGNS = "http://www.w3.org/2000/svg";
  const NODE_W = 168, NODE_H = 52;
  const COL_GAP = 250, ROW_GAP = 74, TOP = 64;

  const DOMAIN_COLORS = {
    "d-vehicle": "#e67e22",
    "d-ad": "#2f9e44",
    "d-robot": "#0e9db0",
    "d-uav": "#7048e8",
  };
  const SOL_COLOR = "#64748b";
  const TAB_ORDER = ["d-ad", "d-robot", "d-uav", "d-vehicle"];

  let svg, viewport, state, onNodeClickCb;
  let curDomain = "d-ad", curOrder = "up";
  let legendOpen = false;
  const hiddenCats = new Set();

  /* ── 组装单方向子图 ── */
  function buildGraph(domainId, order) {
    const P = window.KG_PROBLEMS, K = window.KG;
    const dom = P.domains.find(d => d.id === domainId);
    const color = DOMAIN_COLORS[domainId] || "#888";
    const problems = P.problems.filter(p => p.domain === domainId);
    const solIds = [...new Set(problems.flatMap(p => p.solutions.map(s => s.ref)))];
    const sols = solIds.map(id => P.solutions.find(s => s.id === id));

    // 知识闭包：方案直接依赖 + 依赖的全部前置
    const preds = {};
    K.nodes.forEach(n => preds[n.id] = []);
    K.edges.forEach(([a, b]) => preds[b].push(a));
    const closure = new Set();
    const stack = sols.flatMap(s => [...s.uses]);
    while (stack.length) {
      const id = stack.pop();
      if (closure.has(id)) continue;
      closure.add(id);
      (preds[id] || []).forEach(p => stack.push(p));
    }
    const kNodes = K.nodes.filter(n => closure.has(n.id));
    const tiersAsc = [...new Set(kNodes.map(n => n.tier))].sort((a, b) => a - b);

    const probCol = {
      label: "❓ 控制问题", prob: true,
      items: problems.map(p => ({ id: p.id, title: p.title, en: dom.icon + " " + dom.name, color, kind: "problem" }))
    };
    const solCol = {
      label: "🛠 解决方案", prob: true,
      items: sols.map(s => ({ id: s.id, title: s.name, en: "方案 · Solution", color: SOL_COLOR, kind: "solution" }))
    };
    const tierCols = (order === "down" ? [...tiersAsc].reverse() : tiersAsc).map(t => ({
      label: K.tiers[t],
      items: kNodes.filter(n => n.tier === t).map(n => ({
        id: n.id, title: n.title, en: n.en,
        color: K.categories[n.cat].color, kind: "knowledge", cat: n.cat
      }))
    }));

    const colDefs = order === "down" ? [probCol, solCol, ...tierCols] : [...tierCols, solCol, probCol];
    const cols = [], nodes = [], edges = [];
    colDefs.forEach((c, ci) => {
      cols.push({ label: c.label, prob: !!c.prob });
      c.items.forEach(it => nodes.push({ ...it, col: ci }));
    });

    // 边统一按绘制方向（左 → 右）给出
    if (order === "down") {
      problems.forEach(p => p.solutions.forEach(s => edges.push([p.id, s.ref])));
      sols.forEach(s => s.uses.forEach(u => edges.push([s.id, u])));
      K.edges.forEach(([a, b]) => { if (closure.has(a) && closure.has(b)) edges.push([b, a]); });
    } else {
      K.edges.forEach(([a, b]) => { if (closure.has(a) && closure.has(b)) edges.push([a, b]); });
      sols.forEach(s => s.uses.forEach(u => edges.push([u, s.id])));
      problems.forEach(p => p.solutions.forEach(s => edges.push([s.ref, p.id])));
    }

    return { cols, nodes, edges, dom };
  }

  function computeLayout(graph) {
    const preds = {}, succs = {};
    graph.nodes.forEach(n => { preds[n.id] = []; succs[n.id] = []; });
    graph.edges.forEach(([a, b]) => {
      if (preds[b] && succs[a]) { preds[b].push(a); succs[a].push(b); }
    });
    let maxRows = 0;
    graph.cols.forEach((_, ci) => {
      maxRows = Math.max(maxRows, graph.nodes.filter(n => n.col === ci).length);
    });
    const H = maxRows * ROW_GAP + TOP + 30;
    const pos = {};
    graph.cols.forEach((_, ci) => {
      const col = graph.nodes.filter(n => n.col === ci);
      const totalH = col.length * ROW_GAP;
      const offset = Math.max(TOP, TOP + (H - TOP - totalH) / 2);
      col.forEach((n, i) => {
        pos[n.id] = { x: 60 + ci * COL_GAP, y: offset + i * ROW_GAP };
      });
    });
    return { pos, preds, succs, H };
  }

  function el(tag, attrs, parent) {
    const e = document.createElementNS(SVGNS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function edgePath(p1, p2) {
    if (Math.abs(p2.x - p1.x) < 10) {
      const x = p1.x, y1 = p1.y + NODE_H / 2, y2 = p2.y + NODE_H / 2;
      return `M ${x} ${y1} C ${x - 48} ${y1}, ${x - 48} ${y2}, ${x} ${y2}`;
    }
    const x1 = p1.x + NODE_W, y1 = p1.y + NODE_H / 2;
    const x2 = p2.x, y2 = p2.y + NODE_H / 2;
    const dx = Math.max(50, (x2 - x1) / 2);
    return `M ${x1} ${y1} C ${x1 + dx} ${y1}, ${x2 - dx} ${y2}, ${x2} ${y2}`;
  }

  /* ── 页签栏（方向切换 + 排序切换） ── */
  function renderTabs() {
    const box = document.getElementById("graph-tabs");
    if (!box) return;
    const P = window.KG_PROBLEMS;
    box.innerHTML = "";

    const row = document.createElement("div");
    row.className = "gt-row";

    const tabsWrap = document.createElement("div");
    tabsWrap.className = "gt-domains";
    TAB_ORDER.forEach(id => {
      const d = P.domains.find(x => x.id === id);
      if (!d) return;
      const btn = document.createElement("button");
      btn.className = "gt-tab" + (id === curDomain ? " active" : "");
      btn.style.setProperty("--c", DOMAIN_COLORS[id]);
      btn.textContent = d.icon + " " + d.name;
      btn.onclick = () => setDomain(id);
      tabsWrap.appendChild(btn);
    });
    row.appendChild(tabsWrap);

    const toggle = document.createElement("div");
    toggle.className = "gt-order";
    [["up", "⬆ 自下而上 · 基础→问题"], ["down", "⬇ 自上而下 · 问题→基础"]].forEach(([o, label]) => {
      const btn = document.createElement("button");
      btn.className = o === curOrder ? "active" : "";
      btn.textContent = label;
      btn.onclick = () => setOrder(o);
      toggle.appendChild(btn);
    });
    row.appendChild(toggle);
    box.appendChild(row);

    const dom = P.domains.find(d => d.id === curDomain);
    const desc = document.createElement("div");
    desc.className = "gt-desc";
    desc.innerHTML = `<span class="gt-desc-link" data-goto="${curDomain}">${dom.icon} ${dom.name} 总览 →</span> ${dom.desc}`;
    box.appendChild(desc);
  }

  function setDomain(id) {
    if (!DOMAIN_COLORS[id] || id === curDomain) return;
    curDomain = id;
    renderTabs();
    renderGraph();
  }
  function setOrder(o) {
    if (o === curOrder) return;
    curOrder = o;
    renderTabs();
    renderGraph();
  }

  /* ── 绘制当前页 ── */
  function renderGraph() {
    svg = document.getElementById("graph-svg");
    svg.innerHTML = "";
    const graph = buildGraph(curDomain, curOrder);
    const layout = computeLayout(graph);
    state = { graph, layout, tx: 20, ty: 0, scale: 1 };

    const defs = el("defs", {}, svg);
    const marker = el("marker", {
      id: "arrow", viewBox: "0 0 10 10", refX: 9, refY: 5,
      markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse"
    }, defs);
    el("path", { d: "M 0 0 L 10 5 L 0 10 z", fill: "#b3c0d3" }, marker);

    viewport = el("g", { id: "viewport" }, svg);
    const bandLayer = el("g", {}, viewport);
    const edgeLayer = el("g", {}, viewport);
    const nodeLayer = el("g", {}, viewport);

    const { pos } = layout;

    // 条带 + 列头
    graph.cols.forEach((c, ci) => {
      const x = 60 + ci * COL_GAP;
      if (ci % 2 === 1) {
        el("rect", {
          x: x - 28, y: 6, width: NODE_W + 56, height: layout.H,
          fill: "rgba(53,101,201,0.045)", rx: 16
        }, bandLayer);
      }
      const cls = c.prob ? "g-tier g-tier-prob" : "g-tier";
      const txt = el("text", { x: x + NODE_W / 2, y: 34, "text-anchor": "middle", class: cls }, bandLayer);
      txt.textContent = c.label;
      if (ci > 0) {
        const ax = x - COL_GAP + NODE_W + (COL_GAP - NODE_W) / 2;
        const t2 = el("text", { x: ax, y: 34, "text-anchor": "middle", class: "g-tier-arrow" }, bandLayer);
        t2.textContent = "➜";
      }
    });

    const edgeEls = [];
    graph.edges.forEach(([a, b]) => {
      if (!pos[a] || !pos[b]) return;
      const p = el("path", {
        class: "g-edge", d: edgePath(pos[a], pos[b]),
        "marker-end": "url(#arrow)", "data-a": a, "data-b": b
      }, edgeLayer);
      edgeEls.push(p);
    });

    graph.nodes.forEach(n => {
      const { x, y } = pos[n.id];
      const g = el("g", { class: "g-node", transform: `translate(${x},${y})`, "data-id": n.id, "data-cat": n.cat || n.kind, "data-kind": n.kind }, nodeLayer);
      el("rect", {
        width: NODE_W, height: NODE_H, rx: 10,
        fill: shade(n.color, 0.9), stroke: n.color
      }, g);
      el("rect", { width: 5, height: NODE_H, rx: 2.5, fill: n.color }, g);
      const t1 = el("text", { x: 14, y: 22 }, g); t1.textContent = n.title;
      fitText(t1, NODE_W - 22);
      const t2 = el("text", { x: 14, y: 40, class: "g-en" }, g); t2.textContent = n.en;
      fitText(t2, NODE_W - 22);

      g.addEventListener("click", () => { if (!state.wasPan) onNodeClickCb(n.id); });
      g.addEventListener("mouseenter", () => highlight(n.id, edgeEls, layout));
      g.addEventListener("mouseleave", () => highlight(null, edgeEls, layout));
    });

    // 重新应用分类过滤
    if (hiddenCats.size) {
      document.querySelectorAll(".g-node").forEach(n =>
        n.classList.toggle("dimmed", hiddenCats.has(n.dataset.cat)));
    }

    setupPanZoom();
    renderLegend();
    renderHint();
    fitToScreen();
  }

  function renderHint() {
    const hint = document.getElementById("graph-hint");
    if (!hint) return;
    hint.textContent = curOrder === "down"
      ? "自上而下：❓问题 → 🛠方案 → 知识层⑦→①（箭头 = 拆解/依赖）· 拖动平移 · 滚轮缩放 · 点击节点查看详情"
      : "自下而上：知识层①→⑦ → 🛠方案 → ❓问题（箭头 = 支撑/通往）· 拖动平移 · 滚轮缩放 · 点击节点查看详情";
  }

  function fitText(textEl, maxW) {
    let size = parseFloat(getComputedStyle(textEl).fontSize) || 14;
    while (textEl.getComputedTextLength && textEl.getComputedTextLength() > maxW && size > 8) {
      size -= 0.5; textEl.style.fontSize = size + "px";
    }
  }

  function shade(hex, amt) {
    const n = parseInt(hex.slice(1), 16);
    let r = (n >> 16), g = (n >> 8) & 255, b = n & 255;
    if (amt < 0) { r *= (1 + amt); g *= (1 + amt); b *= (1 + amt); }
    else { r += (255 - r) * amt; g += (255 - g) * amt; b += (255 - b) * amt; }
    return `rgb(${r | 0},${g | 0},${b | 0})`;
  }

  function highlight(id, edgeEls, layout) {
    const nodes = document.querySelectorAll(".g-node");
    if (!id) {
      nodes.forEach(n => n.classList.remove("dimmed"));
      edgeEls.forEach(e => e.classList.remove("hl"));
      if (hiddenCats.size) {
        nodes.forEach(n => n.classList.toggle("dimmed", hiddenCats.has(n.dataset.cat)));
      }
      return;
    }
    const related = new Set([id, ...(layout.preds[id] || []), ...(layout.succs[id] || [])]);
    nodes.forEach(n => n.classList.toggle("dimmed", !related.has(n.dataset.id)));
    edgeEls.forEach(e => e.classList.toggle("hl", e.dataset.a === id || e.dataset.b === id));
  }

  /* ── 平移缩放 ── */
  function applyTransform() {
    viewport.setAttribute("transform", `translate(${state.tx},${state.ty}) scale(${state.scale})`);
  }
  function setupPanZoom() {
    let dragging = false, lx = 0, ly = 0;
    svg.addEventListener("mousedown", e => {
      e.preventDefault(); // 阻止拖拽时选中节点文字
      dragging = true; state.wasPan = false; lx = e.clientX; ly = e.clientY;
      svg.classList.add("dragging");
    });
    window.addEventListener("mousemove", e => {
      if (!dragging) return;
      const dx = e.clientX - lx, dy = e.clientY - ly;
      if (Math.abs(dx) + Math.abs(dy) > 3) state.wasPan = true;
      state.tx += dx; state.ty += dy; lx = e.clientX; ly = e.clientY;
      applyTransform();
    });
    window.addEventListener("mouseup", () => { dragging = false; svg.classList.remove("dragging"); });
    svg.addEventListener("wheel", e => {
      e.preventDefault();
      const factor = e.deltaY < 0 ? 1.12 : 1 / 1.12;
      const rect = svg.getBoundingClientRect();
      const mx = e.clientX - rect.left, my = e.clientY - rect.top;
      state.tx = mx - (mx - state.tx) * factor;
      state.ty = my - (my - state.ty) * factor;
      state.scale *= factor;
      state.scale = Math.min(3, Math.max(0.15, state.scale));
      applyTransform();
    }, { passive: false });
    let lastTouch = null;
    svg.addEventListener("touchstart", e => { if (e.touches.length === 1) lastTouch = e.touches[0]; });
    svg.addEventListener("touchmove", e => {
      if (e.touches.length === 1 && lastTouch) {
        state.tx += e.touches[0].clientX - lastTouch.clientX;
        state.ty += e.touches[0].clientY - lastTouch.clientY;
        lastTouch = e.touches[0]; applyTransform(); e.preventDefault();
      }
    }, { passive: false });
  }

  function fitToScreen() {
    if (!state) return;
    const rect = svg.getBoundingClientRect();
    if (rect.width < 10 || rect.height < 10) { state.needsFit = true; return; } // 视图隐藏时跳过，显示时再调用
    state.needsFit = false;
    const { pos } = state.layout;
    let maxX = 0, maxY = 0;
    Object.values(pos).forEach(p => { maxX = Math.max(maxX, p.x + NODE_W); maxY = Math.max(maxY, p.y + NODE_H); });
    maxX += 30; maxY += 20;
    const s = Math.min((rect.width - 40) / maxX, (rect.height - 40) / maxY, 1.15);
    state.scale = s;
    state.tx = (rect.width - maxX * s) / 2 + 10;
    state.ty = (rect.height - maxY * s) / 2 + 6;
    applyTransform();
  }

  /* ── 图例（可折叠，默认收起，避免盖住节点） ── */
  function renderLegend() {
    const box = document.getElementById("legend");
    box.classList.toggle("collapsed", !legendOpen);
    box.innerHTML = "";

    const head = document.createElement("div");
    head.className = "legend-head";
    const hiddenNote = hiddenCats.size ? ` · 已隐藏 ${hiddenCats.size} 类` : "";
    head.innerHTML = legendOpen
      ? `🗂 图例与筛选 <span class="legend-toggle">收起 ×</span>`
      : `🗂 图例与筛选${hiddenNote}`;
    head.onclick = () => { legendOpen = !legendOpen; renderLegend(); };
    box.appendChild(head);
    if (!legendOpen) return;

    const h0 = document.createElement("h3");
    h0.textContent = "本页结构";
    box.appendChild(h0);
    const layers = curOrder === "down"
      ? [["❓", "控制问题"], ["🛠", "解决方案"], ["📚", "知识层 ⑦→①（右端最基础）"]]
      : [["📚", "知识层 ①→⑦（左端最基础）"], ["🛠", "解决方案"], ["❓", "控制问题"]];
    layers.forEach(([ic, name]) => {
      const item = document.createElement("div");
      item.className = "legend-item";
      item.style.cursor = "default";
      item.innerHTML = `<span style="width:14px;text-align:center">${ic}</span>${name}`;
      box.appendChild(item);
    });
    const h = document.createElement("h3");
    h.style.marginTop = "10px";
    h.textContent = "知识分类（点击过滤）";
    box.appendChild(h);
    Object.entries(window.KG.categories).forEach(([key, c]) => {
      const item = document.createElement("div");
      item.className = "legend-item";
      if (hiddenCats.has(key)) item.classList.add("off");
      item.innerHTML = `<span class="legend-dot" style="background:${c.color}"></span>${c.name}`;
      item.onclick = () => {
        if (hiddenCats.has(key)) hiddenCats.delete(key);
        else hiddenCats.add(key);
        item.classList.toggle("off", hiddenCats.has(key));
        document.querySelectorAll(".g-node").forEach(n =>
          n.classList.toggle("dimmed", hiddenCats.has(n.dataset.cat)));
      };
      box.appendChild(item);
    });
  }

  function render(onNodeClick) {
    onNodeClickCb = onNodeClick;
    renderTabs();
    renderGraph();
  }

  // 图谱视图重新显示时调用：仅当隐藏期间发生过重绘才重新适配视野
  function shown() {
    if (state && state.needsFit) fitToScreen();
  }

  return { render, fitToScreen, setDomain, setOrder, shown };
})();
