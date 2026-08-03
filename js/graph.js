/* 知识图谱渲染（自上而下）：
 * 列结构：🎯研究方向 → ❓控制问题 → 🛠解决方案 → ⑦鲁棒与先进 → … → ①建模基础
 * 所有箭头统一含义：拆解 / 依赖（左边的大问题被拆成右边更小的知识点）
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

  let svg, viewport, state;

  /* ── 组装完整图：问题层 + 知识层（知识边反向 = 依赖方向） ── */
  function buildGraph() {
    const P = window.KG_PROBLEMS, K = window.KG;
    const cols = [];
    const nodes = [];
    const edges = [];

    // 列 0：研究方向
    cols.push("🎯 研究方向");
    P.domains.forEach(d => nodes.push({
      id: d.id, title: d.icon + " " + d.name, en: "研究方向 · Domain",
      color: DOMAIN_COLORS[d.id] || "#888", col: 0, kind: "domain"
    }));

    // 列 1：控制问题（继承所属方向颜色）
    cols.push("❓ 控制问题");
    P.problems.forEach(p => {
      const dom = P.domains.find(d => d.id === p.domain);
      nodes.push({
        id: p.id, title: p.title, en: dom.icon + " " + dom.name,
        color: DOMAIN_COLORS[p.domain] || "#888", col: 1, kind: "problem"
      });
      edges.push([p.domain, p.id]);
    });

    // 列 2：解决方案
    cols.push("🛠 解决方案");
    P.solutions.forEach(s => nodes.push({
      id: s.id, title: s.name, en: "方案 · Solution",
      color: SOL_COLOR, col: 2, kind: "solution"
    }));
    P.problems.forEach(p => p.solutions.forEach(s => edges.push([p.id, s.ref])));
    P.solutions.forEach(s => s.uses.forEach(u => edges.push([s.id, u])));

    // 列 3+：知识层，⑦ → ①（越右越基础，拆到最小知识点）
    const tierKeys = Object.keys(K.tiers).map(Number).sort((a, b) => b - a);
    tierKeys.forEach((t, i) => {
      cols.push(K.tiers[t]);
      K.nodes.filter(n => n.tier === t).forEach(n => nodes.push({
        id: n.id, title: n.title, en: n.en,
        color: K.categories[n.cat].color, col: 3 + i, kind: "knowledge", cat: n.cat
      }));
    });
    // 知识依赖边反向绘制：B 依赖 A ⇒ B → A（与拆解方向一致）
    K.edges.forEach(([a, b]) => edges.push([b, a]));

    return { cols, nodes, edges };
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

  function render(onNodeClick) {
    svg = document.getElementById("graph-svg");
    svg.innerHTML = "";
    const graph = buildGraph();
    const layout = computeLayout(graph);
    state = { graph, layout, tx: 20, ty: 0, scale: 1, hiddenCats: new Set() };

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
    graph.cols.forEach((label, ci) => {
      const x = 60 + ci * COL_GAP;
      if (ci % 2 === 1) {
        el("rect", {
          x: x - 28, y: 6, width: NODE_W + 56, height: layout.H,
          fill: "rgba(53,101,201,0.045)", rx: 16
        }, bandLayer);
      }
      // 问题层三列用不同色调的列头
      const cls = ci <= 2 ? "g-tier g-tier-prob" : "g-tier";
      const txt = el("text", { x: x + NODE_W / 2, y: 34, "text-anchor": "middle", class: cls }, bandLayer);
      txt.textContent = label;
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

      g.addEventListener("click", () => { if (!state.wasPan) onNodeClick(n.id); });
      g.addEventListener("mouseenter", () => highlight(n.id, edgeEls, layout));
      g.addEventListener("mouseleave", () => highlight(null, edgeEls, layout));
    });

    setupPanZoom();
    renderLegend();
    fitToScreen();
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
    const { pos } = state.layout;
    let maxX = 0, maxY = 0;
    Object.values(pos).forEach(p => { maxX = Math.max(maxX, p.x + NODE_W); maxY = Math.max(maxY, p.y + NODE_H); });
    maxX += 30; maxY += 20;
    const rect = svg.getBoundingClientRect();
    const s = Math.min((rect.width - 40) / maxX, (rect.height - 40) / maxY, 1.15);
    state.scale = s;
    state.tx = (rect.width - maxX * s) / 2 + 10;
    state.ty = (rect.height - maxY * s) / 2 + 6;
    applyTransform();
  }

  /* ── 图例 ── */
  function renderLegend() {
    const box = document.getElementById("legend");
    box.innerHTML = "<h3>图谱结构</h3>";
    const layers = [
      ["🎯", "研究方向"], ["❓", "控制问题"], ["🛠", "解决方案"], ["📚", "知识层 ⑦→①（右端最基础）"],
    ];
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
      item.innerHTML = `<span class="legend-dot" style="background:${c.color}"></span>${c.name}`;
      item.onclick = () => {
        if (state.hiddenCats.has(key)) state.hiddenCats.delete(key);
        else state.hiddenCats.add(key);
        item.classList.toggle("off", state.hiddenCats.has(key));
        document.querySelectorAll(".g-node").forEach(n =>
          n.classList.toggle("dimmed", state.hiddenCats.has(n.dataset.cat)));
      };
      box.appendChild(item);
    });
  }

  return { render, fitToScreen };
})();
