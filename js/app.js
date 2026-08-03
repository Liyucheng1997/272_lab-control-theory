/* 应用主逻辑：视图切换、四类节点详情（知识/问题/方案/方向）、搜索、学习路线 */
(function () {
  const $ = s => document.querySelector(s);
  const views = { graph: $("#graph-view"), detail: $("#detail-view"), path: $("#path-view") };
  const PB = window.KG_PROBLEMS;
  const byId = {}; window.KG.nodes.forEach(n => byId[n.id] = n);
  const probById = {}; PB.problems.forEach(p => probById[p.id] = p);
  const solById = {}; PB.solutions.forEach(s => solById[s.id] = s);
  const domById = {}; PB.domains.forEach(d => domById[d.id] = d);
  const preds = {}, succs = {};
  window.KG.nodes.forEach(n => { preds[n.id] = []; succs[n.id] = []; });
  window.KG.edges.forEach(([a, b]) => { preds[b].push(a); succs[a].push(b); });

  function show(name) {
    Object.entries(views).forEach(([k, v]) => v.hidden = (k !== name));
  }

  /* ── KaTeX 渲染 ── */
  function renderMath(container) {
    if (window.renderMathInElement) {
      renderMathInElement(container, {
        delimiters: [
          { left: "$$", right: "$$", display: true },
          { left: "$", right: "$", display: false },
          { left: "\\[", right: "\\]", display: true },
          { left: "\\(", right: "\\)", display: false },
        ],
        throwOnError: false,
      });
    }
  }

  /* ── Python 代码高亮（轻量） ── */
  function escapeHtml(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  }
  function highlightPython(src) {
    let s = escapeHtml(src);
    const rules = [
      [/(#[^\n]*)/g, '<span class="c">$1</span>'],
      [/(&quot;&quot;&quot;[\s\S]*?&quot;&quot;&quot;|&#39;[^&\n]*?&#39;|&quot;[^&\n]*?&quot;)/g, '<span class="s">$1</span>'],
      [/('[^'\n]*'|"[^"\n]*")/g, '<span class="s">$1</span>'],
      [/\b(def|class|return|if|elif|else|for|while|in|import|from|as|with|lambda|not|and|or|None|True|False|assert|break|continue|pass|raise|try|except|yield)\b/g, '<span class="k">$1</span>'],
      [/\b(\d+\.?\d*(?:e-?\d+)?)\b/g, '<span class="n">$1</span>'],
      [/\b(def\s+)([a-zA-Z_]\w*)/g, '$1<span class="f">$2</span>'],
    ];
    rules.forEach(([re, rep]) => { s = s.replace(re, rep); });
    return s;
  }

  /* ── 公共构件 ── */
  function chip(id) {
    const label = byId[id] ? byId[id].title
      : probById[id] ? probById[id].title
      : solById[id] ? solById[id].name
      : domById[id] ? domById[id].icon + " " + domById[id].name : id;
    return `<span class="d-chip" data-goto="${id}">${label}</span>`;
  }

  function depClosure(ids) {
    const seen = new Set();
    const stack = [...ids];
    while (stack.length) {
      const id = stack.pop();
      if (seen.has(id)) continue;
      seen.add(id);
      (preds[id] || []).forEach(p => stack.push(p));
    }
    return seen;
  }

  // 第一性拆解块：知识点闭包按 ①→⑦ 学习顺序展开
  function closureHtml(uses, open = false) {
    const closure = depClosure([...uses]);
    const byTier = {};
    closure.forEach(id => {
      const n = byId[id];
      if (!n) return;
      (byTier[n.tier] = byTier[n.tier] || []).push(id);
    });
    const tierKeys = Object.keys(byTier).map(Number).sort((a, b) => a - b);
    let html = `<details class="prob-expand"${open ? " open" : ""}><summary>🔬 第一性拆解：共 ${closure.size} 个知识点（按学习顺序 ①→⑦ 展开）</summary><div class="prob-tiers">`;
    tierKeys.forEach(t => {
      html += `<div class="prob-tier-row"><span class="prob-tier-label">${window.KG.tiers[t]}</span><span>${byTier[t].map(chip).join("")}</span></div>`;
    });
    html += `</div></details>`;
    return html;
  }

  /* ── 侧栏：问题层 + 知识层 ── */
  function sidebarHtml(activeId) {
    let html = `<div class="d-side-back" data-goto-graph>← 返回图谱</div>`;
    PB.domains.forEach(d => {
      html += `<div class="d-side-cat"><span data-goto="${d.id}" style="cursor:pointer">${d.icon} ${d.name}</span></div>`;
      d.problems.forEach(pid => {
        html += `<div class="d-side-item${pid === activeId ? " active" : ""}" data-goto="${pid}">
          <span class="legend-dot" style="background:#94a3b8"></span>${probById[pid].title}</div>`;
      });
    });
    html += `<div class="d-side-cat">🛠 解决方案</div>`;
    PB.solutions.forEach(s => {
      html += `<div class="d-side-item${s.id === activeId ? " active" : ""}" data-goto="${s.id}">
        <span class="legend-dot" style="background:#64748b"></span>${s.name}</div>`;
    });
    Object.keys(window.KG.tiers).map(Number).sort((a, b) => a - b).forEach(t => {
      html += `<div class="d-side-cat">${window.KG.tiers[t]}</div>`;
      window.KG.nodes.filter(n => n.tier === t).forEach(n => {
        const c = window.KG.categories[n.cat];
        html += `<div class="d-side-item${n.id === activeId ? " active" : ""}" data-goto="${n.id}">
          <span class="legend-dot" style="background:${c.color}"></span>${n.title}</div>`;
      });
    });
    return `<aside class="d-sidebar">${html}</aside>`;
  }

  function setDetail(activeId, articleHtml) {
    views.detail.innerHTML = sidebarHtml(activeId) +
      `<div class="d-main"><article id="detail-article">${articleHtml}</article></div>`;
    renderMath($("#detail-article"));
    show("detail");
    views.detail.querySelector(".d-main").scrollTop = 0;
  }

  /* ── 知识节点详情（页签式） ── */
  const TAB_GROUPS = [
    { key: "intro", label: "💡 直觉", icons: ["💡"] },
    { key: "def", label: "📐 定义", icons: ["📐"] },
    { key: "proof", label: "🧾 证明", icons: ["🧾", "🔬"] },
  ];
  let currentSim = null;

  function openNode(id) {
    window.Sims && window.Sims.stopAll();
    currentSim = null;
    const n = byId[id];
    const content = (window.KG_CONTENT && window.KG_CONTENT[id]) || null;
    const cat = window.KG.categories[n.cat];

    let head = `
      <div><span class="d-cat-badge" style="background:${cat.color}">${cat.name} · ${window.KG.tiers[n.tier]}</span></div>
      <h1>${n.title}</h1>
      <div class="d-en">${n.en}</div>
      <div class="d-tagline">${n.tagline}</div>`;
    const pre = preds[id], suc = succs[id];
    const usedBy = PB.solutions.filter(s => s.uses.includes(id)).map(s => s.id);
    if (pre.length || suc.length || usedBy.length) {
      head += `<div class="d-links">`;
      if (pre.length) head += `<span class="lk-label">⬅ 前置：</span>` + pre.map(chip).join("");
      if (suc.length) head += `<span class="lk-label" style="margin-left:14px">➡ 通往：</span>` + suc.map(chip).join("");
      head += `</div>`;
      if (usedBy.length) head += `<div class="d-links"><span class="lk-label">🛠 用于方案：</span>` + usedBy.map(chip).join("") + `</div>`;
    }

    const tabs = [];
    if (content) {
      TAB_GROUPS.forEach(g => {
        const secs = content.sections.filter(s => g.icons.includes(s.icon || "🧾"));
        if (secs.length) {
          tabs.push({
            key: g.key, label: g.label,
            html: secs.map(s => `<h2><span class="h-icon">${s.icon}</span>${s.title}</h2>${s.html}`).join("")
          });
        }
      });
      if (content.sim && window.Sims && window.Sims.registry[content.sim]) {
        tabs.push({ key: "sim", label: "🎛️ 实验", html: `<div id="sim-mount"></div>`, sim: content.sim });
      }
      if (content.code) {
        tabs.push({
          key: "code", label: "💻 代码",
          html: `<p class="sim-desc" style="color:var(--text-dim)">${content.code.title || ""}</p>
                 <pre class="code">${highlightPython(content.code.text)}</pre>`
        });
      }
    } else {
      tabs.push({ key: "intro", label: "💡 直觉", html: `<div class="box note"><span class="box-title">内容建设中</span>该节点的详细内容尚未编写。</div>` });
    }

    const tabBar = `<div class="d-tabs">` + tabs.map((t, i) =>
      `<button class="d-tab${i === 0 ? " active" : ""}" data-tab="${t.key}">${t.label}</button>`).join("") + `</div>`;
    const panels = tabs.map((t, i) =>
      `<div class="d-panel${i === 0 ? " active" : ""}" data-panel="${t.key}">${t.html}</div>`).join("");

    const lp = window.KG.learningPath;
    const idx = lp.indexOf(id);
    let nav = "";
    if (idx >= 0) {
      nav = `<div class="d-nav">`;
      nav += idx > 0
        ? `<span class="d-chip" data-goto="${lp[idx - 1]}">← 上一站：${byId[lp[idx - 1]].title}</span>`
        : `<span></span>`;
      nav += idx < lp.length - 1
        ? `<span class="d-chip" data-goto="${lp[idx + 1]}">下一站：${byId[lp[idx + 1]].title} →</span>`
        : `<span></span>`;
      nav += `</div>`;
    }

    setDetail(id, head + tabBar + panels + nav);

    const simTab = tabs.find(t => t.sim);
    views.detail.querySelectorAll(".d-tab").forEach(btn => {
      btn.addEventListener("click", () => {
        views.detail.querySelectorAll(".d-tab").forEach(b => b.classList.toggle("active", b === btn));
        views.detail.querySelectorAll(".d-panel").forEach(p =>
          p.classList.toggle("active", p.dataset.panel === btn.dataset.tab));
        if (simTab && btn.dataset.tab === "sim") {
          const mountEl = $("#sim-mount");
          if (mountEl && (currentSim !== simTab.sim || !mountEl.childElementCount)) {
            mountEl.innerHTML = "";
            window.Sims.mount(simTab.sim, mountEl);
            currentSim = simTab.sim;
          }
        } else {
          window.Sims && window.Sims.stopAll();
          currentSim = null;
        }
        views.detail.querySelector(".d-main").scrollTop = 0;
      });
    });
  }

  /* ── 问题详情 ── */
  function openProblem(id) {
    window.Sims && window.Sims.stopAll();
    const p = probById[id];
    const dom = domById[p.domain];
    const color = { "d-vehicle": "#e67e22", "d-ad": "#2f9e44", "d-robot": "#0e9db0", "d-uav": "#7048e8" }[p.domain] || "#888";
    let html = `
      <div><span class="d-cat-badge" style="background:${color}">${dom.icon} ${dom.name} · 控制问题</span></div>
      <h1>❓ ${p.title}</h1>
      <div class="d-tagline">${p.desc}</div>
      <h2><span class="h-icon">🛠</span>解决方案（点击方案名深入）</h2>`;
    const allUses = new Set();
    p.solutions.forEach(sref => {
      const s = solById[sref.ref];
      s.uses.forEach(u => allUses.add(u));
      html += `<div class="prob-sol">
        <div class="prob-sol-name">${chip(s.id)}</div>
        <div class="prob-sol-note">${sref.note}</div>
        <div class="prob-sol-uses">核心知识：${s.uses.map(chip).join("")}</div>
      </div>`;
    });
    html += `<h2><span class="h-icon">🔬</span>拆到最小知识点</h2>` + closureHtml(allUses, true);
    // 同方向的其它问题
    const siblings = dom.problems.filter(x => x !== id);
    if (siblings.length) {
      html += `<div class="d-links" style="margin-top:26px"><span class="lk-label">同方向问题：</span>${siblings.map(chip).join("")}</div>`;
    }
    setDetail(id, html);
  }

  /* ── 方案详情 ── */
  function openSolution(id) {
    window.Sims && window.Sims.stopAll();
    const s = solById[id];
    const usedIn = PB.problems.filter(p => p.solutions.some(x => x.ref === id));
    let html = `
      <div><span class="d-cat-badge" style="background:#64748b">🛠 解决方案</span></div>
      <h1>${s.name}</h1>
      <div class="d-tagline">${s.desc}</div>
      <h2><span class="h-icon">❓</span>解决哪些问题</h2>
      <div class="d-links">${usedIn.map(p => chip(p.id)).join("")}</div>`;
    usedIn.forEach(p => {
      const note = p.solutions.find(x => x.ref === id).note;
      html += `<div class="prob-sol"><div class="prob-sol-name">${chip(p.id)}</div><div class="prob-sol-note">${note}</div></div>`;
    });
    html += `<h2><span class="h-icon">📚</span>直接依赖的知识</h2>
      <div class="d-links">${s.uses.map(chip).join("")}</div>` +
      closureHtml(s.uses, true);
    setDetail(id, html);
  }

  /* ── 方向详情 ── */
  function openDomain(id) {
    window.Sims && window.Sims.stopAll();
    const d = domById[id];
    const color = { "d-vehicle": "#e67e22", "d-ad": "#2f9e44", "d-robot": "#0e9db0", "d-uav": "#7048e8" }[id] || "#888";
    let html = `
      <div><span class="d-cat-badge" style="background:${color}">🎯 研究方向</span></div>
      <h1>${d.icon} ${d.name}</h1>
      <div class="d-tagline">${d.desc}</div>
      <h2><span class="h-icon">❓</span>核心控制问题</h2>`;
    d.problems.forEach(pid => {
      const p = probById[pid];
      html += `<div class="prob-card" data-goto="${pid}" style="cursor:pointer">
        <div class="prob-title">❓ ${p.title}</div>
        <div class="prob-desc">${p.desc}</div>
        <div class="d-links"><span class="lk-label">方案：</span>${p.solutions.map(x => chip(x.ref)).join("")}</div>
      </div>`;
    });
    setDetail(id, html);
  }

  /* ── 统一分发 ── */
  function open(id) {
    if (byId[id]) return openNode(id);
    if (probById[id]) return openProblem(id);
    if (solById[id]) return openSolution(id);
    if (domById[id]) return openDomain(id);
  }

  /* ── 事件委托 ── */
  document.addEventListener("click", e => {
    const go = e.target.closest("[data-goto]");
    if (go) { open(go.dataset.goto); return; }
    const back = e.target.closest("[data-goto-graph]");
    if (back) { window.Sims && window.Sims.stopAll(); show("graph"); return; }
  });

  $("#brand-home").onclick = () => { window.Sims && window.Sims.stopAll(); show("graph"); };

  /* ── 学习路线 ── */
  $("#btn-path").onclick = () => {
    window.Sims && window.Sims.stopAll();
    const list = $("#path-list");
    let html = `<h1>📜 推荐学习路线</h1>
      <div class="path-sub">从零基础到前沿，共 ${window.KG.learningPath.length} 站，分 ${Object.keys(window.KG.tiers).length} 个层级严格递进——每一站的前置知识都已在之前覆盖。（图谱是自上而下的拆解方向，这里是自下而上的学习方向。）</div>`;
    let lastTier = 0;
    window.KG.learningPath.forEach((id, i) => {
      const n = byId[id];
      const c = window.KG.categories[n.cat];
      if (n.tier !== lastTier) {
        html += `<h2 class="path-tier">${window.KG.tiers[n.tier]}</h2>`;
        lastTier = n.tier;
      }
      html += `<div class="path-item" data-goto="${id}">
        <span class="path-num" style="background:${c.color}">${i + 1}</span>
        <div><div class="pi-title">${n.title} <small style="color:var(--text-dim);font-weight:400">${n.en}</small></div>
        <div class="pi-tag">${n.tagline}</div></div></div>`;
    });
    list.innerHTML = html;
    show("path");
  };

  /* ── 搜索（覆盖四类节点） ── */
  const searchBox = $("#search"), results = $("#search-results");
  let activeIdx = -1;
  const searchIndex = [
    ...window.KG.nodes.map(n => ({ id: n.id, title: n.title, sub: n.en, tag: window.KG.categories[n.cat].name, text: (n.title + n.en + n.id + n.tagline).toLowerCase() })),
    ...PB.problems.map(p => ({ id: p.id, title: "❓ " + p.title, sub: domById[p.domain].name, tag: "控制问题", text: (p.title + p.desc).toLowerCase() })),
    ...PB.solutions.map(s => ({ id: s.id, title: "🛠 " + s.name, sub: "", tag: "解决方案", text: (s.name + s.desc).toLowerCase() })),
    ...PB.domains.map(d => ({ id: d.id, title: d.icon + " " + d.name, sub: "", tag: "研究方向", text: d.name.toLowerCase() })),
  ];
  function doSearch(q) {
    if (!q.trim()) { results.classList.remove("open"); return; }
    const ql = q.toLowerCase();
    const hits = searchIndex.filter(it => it.text.includes(ql)).slice(0, 9);
    if (!hits.length) { results.classList.remove("open"); return; }
    results.innerHTML = hits.map(it =>
      `<div class="sr-item" data-goto="${it.id}"><strong>${it.title}</strong>
       <span>${it.sub}</span><span class="sr-cat">${it.tag}</span></div>`).join("");
    results.classList.add("open");
    activeIdx = -1;
  }
  searchBox.addEventListener("input", e => doSearch(e.target.value));
  searchBox.addEventListener("keydown", e => {
    const items = results.querySelectorAll(".sr-item");
    if (e.key === "ArrowDown") { activeIdx = Math.min(activeIdx + 1, items.length - 1); }
    else if (e.key === "ArrowUp") { activeIdx = Math.max(activeIdx - 1, 0); }
    else if (e.key === "Enter" && items[activeIdx >= 0 ? activeIdx : 0]) {
      items[activeIdx >= 0 ? activeIdx : 0].click();
      results.classList.remove("open"); searchBox.value = ""; return;
    } else if (e.key === "Escape") { results.classList.remove("open"); return; }
    else return;
    items.forEach((it, i) => it.classList.toggle("active", i === activeIdx));
    e.preventDefault();
  });
  document.addEventListener("click", e => {
    if (!e.target.closest("#search") && !e.target.closest("#search-results"))
      results.classList.remove("open");
  });

  /* ── 启动 ── */
  window.Graph.render(open);
  window.addEventListener("resize", () => {
    if (!views.graph.hidden) window.Graph.fitToScreen();
  });
})();
