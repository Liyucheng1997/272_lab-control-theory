/* 控制之美 · 交互课程 —— 路由 / 目录 / 进度 / 内容后处理 */
window.Course = (function () {
  "use strict";
  const lessons = [];
  const byId = {};
  const $ = s => document.querySelector(s);
  const store = {
    get(k, d) { try { const v = localStorage.getItem(k); return v === null ? d : JSON.parse(v); } catch (e) { return d; } },
    set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { } },
  };
  let done = new Set(store.get("cb-done", []));

  const VOLS = {
    1: { name: "卷 1", title: "控制理论从传递函数到状态空间", color: "var(--accent)" },
    2: { name: "卷 2", title: "最优化控制、MPC 与卡尔曼滤波器", color: "var(--c2)" },
  };

  /* 课程注册：{id, vol, ch, chName, no, title, sub, goals[], minutes, html, sims:{name: fn(el)}} */
  function add(L) { lessons.push(L); byId[L.id] = L; }

  function order() { return lessons.slice().sort((a, b) => a.vol - b.vol || a.ord - b.ord); }

  /* ── 侧栏 ── */
  function buildSidebar(cur) {
    const sb = $("#sidebar"); let html = `<div class="sb-item ${!cur ? "active" : ""}" data-go=""><span class="no">⌂</span><span>课程首页 · 学习地图</span></div>`;
    let lastVol = 0, lastCh = "";
    order().forEach(L => {
      if (L.vol !== lastVol) { html += `<div class="sb-vol">${VOLS[L.vol].name}<span>${VOLS[L.vol].title}</span></div>`; lastVol = L.vol; lastCh = ""; }
      if (L.chName !== lastCh) { html += `<div class="sb-ch">${L.chName}</div>`; lastCh = L.chName; }
      html += `<div class="sb-item ${cur === L.id ? "active" : ""}" data-go="${L.id}"><span class="no">${L.no}</span><span>${L.title}</span>${done.has(L.id) ? '<span class="ok">✓</span>' : ""}</div>`;
    });
    sb.innerHTML = html;
    sb.querySelectorAll("[data-go]").forEach(el => el.onclick = () => { location.hash = "#/" + el.dataset.go; document.body.classList.remove("sb-open"); });
    const act = sb.querySelector(".active"); if (act && act.scrollIntoViewIfNeeded) act.scrollIntoViewIfNeeded(); else if (act) act.scrollIntoView({ block: "nearest" });
  }
  function updateProgress() {
    const n = lessons.length, k = lessons.filter(L => done.has(L.id)).length;
    $("#progress-fill").style.width = (n ? k / n * 100 : 0) + "%";
    $("#progress-text").textContent = `已学 ${k} / ${n}`;
  }

  /* ── 内容后处理：公式块 / 测验 / 仿真挂载 ── */
  function enhance(root, L) {
    // 公式块
    root.querySelectorAll(".eq").forEach(eq => {
      if (eq.dataset.done) return; eq.dataset.done = 1;
      const dl = eq.querySelector("dl"); const say = eq.querySelector(".say");
      const main = document.createElement("div"); main.className = "eq-main";
      [...eq.childNodes].forEach(n => { if (n !== dl && n !== say) main.appendChild(n); });
      eq.innerHTML = "";
      if (eq.dataset.title) { const t = document.createElement("div"); t.className = "eq-title"; t.textContent = eq.dataset.title; eq.appendChild(t); }
      eq.appendChild(main);
      if (eq.dataset.no) { const s = document.createElement("span"); s.className = "eq-no"; s.textContent = eq.dataset.no; main.appendChild(s); }
      if (say) eq.appendChild(say);
      if (dl) { dl.className = "legend"; eq.appendChild(dl); }
    });
    // 推导步骤自动编号
    root.querySelectorAll("details.derive .step").forEach((s, i, all) => {
      if (s.querySelector(":scope > i")) return;
      const idx = [...s.parentElement.querySelectorAll(":scope > .step")].indexOf(s) + 1;
      const wrap = document.createElement("div"); while (s.firstChild) wrap.appendChild(s.firstChild);
      const ii = document.createElement("i"); ii.textContent = idx; s.append(ii, wrap);
    });
    CE.renderMath(root);
    // 测验
    root.querySelectorAll(".quiz").forEach(q => {
      q.querySelectorAll(".opt").forEach(b => b.onclick = () => {
        if (q.classList.contains("done")) return;
        q.classList.add("done");
        q.querySelectorAll(".opt").forEach(o => { if (o.hasAttribute("data-ok")) o.classList.add("right"); });
        if (!b.hasAttribute("data-ok")) b.classList.add("wrong");
      });
    });
    // 仿真
    root.querySelectorAll(".sim[data-id]").forEach(el => {
      const fn = L && L.sims && L.sims[el.dataset.id];
      if (!fn) { el.innerHTML = `<div class="warn">仿真 ${el.dataset.id} 未定义</div>`; return; }
      try { fn(el); } catch (e) { console.error(e); el.innerHTML = `<div class="warn">仿真加载失败：${e.message}</div>`; }
    });
    // 图
    root.querySelectorAll("[data-fig]").forEach(el => { const f = window.Figs && Figs[el.dataset.fig]; if (f) { el.innerHTML = f(el.dataset); el.classList.add("diagram"); CE.renderMath(el); } });
    // 仿真卡片里后插入的文字（提示、任务）再渲染一次公式
    root.querySelectorAll(".sim-hint, .tasks").forEach(el => CE.renderMath(el));
  }

  /* ── 首页 ── */
  function renderHome() {
    const page = $("#page");
    const card = L => `<div class="road-card" data-go="${L.id}"><div class="no">${VOLS[L.vol].name} · ${L.no}</div><h4>${L.title}</h4><p>${L.sub || ""}</p>${done.has(L.id) ? '<span class="ok">✓ 已学</span>' : ""}</div>`;
    const vol = v => order().filter(L => L.vol === v).map(card).join("");
    page.innerHTML = `
      <section class="hero">
        <div class="crumb">基于《控制之美》卷 1 · 卷 2（王天威 DR_CAN 等著）的交互式学习页面</div>
        <h1>让一个<em>从没学过控制</em>的人，<br>也能一步步看懂控制理论</h1>
        <p>每一课都按“故事 → 直觉 → 公式 → 亲手调参”的顺序展开。公式下面都有逐个符号的解释，推导过程可以展开，所有图都能拖、能调、能动。你不需要会写代码，只需要会拖滑块。</p>
      </section>
      <div class="tip" data-t="怎么用这套课程">
        <ul>
          <li><b>按顺序学</b>：左侧目录从上到下就是推荐顺序。卷 1 讲“怎样描述一个会动的系统、怎样让它听话”，卷 2 讲“怎样让它以最省力的方式听话，以及看不清时怎么估计”。</li>
          <li><b>先玩再看公式</b>：每个“交互”卡片下方都有“试一试”任务，先动手观察现象，再回头读公式，会容易很多。</li>
          <li><b>推导可以先跳过</b>：带 ▸ 的“推导过程”是给想深究的同学的，第一次读可以只看结论。</li>
          <li>学完一课点底部的 <b>“标记已学完”</b>，顶栏会记录你的进度（只保存在本机浏览器）。</li>
        </ul>
      </div>
      <h2>整条故事线</h2>
      <div class="flow">
        <span class="fstep">① 描述系统<br><small>微分方程·传递函数·状态空间</small></span><span class="farr">→</span>
        <span class="fstep">② 分析系统<br><small>时域响应·稳定性·频率响应</small></span><span class="farr">→</span>
        <span class="fstep">③ 经典设计<br><small>PI/PID·根轨迹</small></span><span class="farr">→</span>
        <span class="fstep">④ 现代设计<br><small>能控能观·极点配置·观测器</small></span><span class="farr">→</span>
        <span class="fstep">⑤ 最优控制<br><small>动态规划·LQR·MPC</small></span><span class="farr">→</span>
        <span class="fstep">⑥ 状态估计<br><small>卡尔曼滤波·EKF</small></span>
      </div>
      <div class="vol-title"><h2>卷 1</h2><span>${VOLS[1].title}</span></div>
      <div class="road">${vol(1)}</div>
      <div class="vol-title"><h2>卷 2</h2><span>${VOLS[2].title}</span></div>
      <div class="road">${vol(2)}</div>
      <p style="margin-top:40px;color:var(--muted);font-size:.85rem">说明：本页面为个人学习用的交互式笔记，文字为对原书内容的重新讲解，例题参数尽量与原书保持一致；部分原书未给出具体数值的仿真参数由本页面补充，均已注明。</p>`;
    page.querySelectorAll("[data-go]").forEach(el => el.onclick = () => location.hash = "#/" + el.dataset.go);
    document.title = "控制之美 · 交互课程";
  }

  /* ── 课程页 ── */
  function renderLesson(L) {
    const page = $("#page"); const list = order(); const i = list.indexOf(L);
    const prev = list[i - 1], next = list[i + 1];
    page.innerHTML = `
      <header class="lesson-head">
        <div class="crumb">${VOLS[L.vol].name} · ${L.chName} · 对应原书 ${L.ref || ""}</div>
        <h1>${L.title}</h1>
        <div class="sub">${L.sub || ""}</div>
        ${L.goals ? `<div class="goals"><b>学完你会</b><ul>${L.goals.map(g => `<li>${g}</li>`).join("")}</ul>${L.pre ? `<b>需要先懂</b><div>${L.pre}</div>` : ""}${L.minutes ? `<b>预计用时</b><div>${L.minutes} 分钟</div>` : ""}</div>` : ""}
      </header>
      <article id="lesson-body">${L.html}</article>
      <div class="lesson-foot">
        ${prev ? `<a class="nav prev" href="#/${prev.id}"><small>← 上一课</small>${prev.title}</a>` : `<a class="nav prev" href="#/"><small>←</small>课程首页</a>`}
        <button class="done-btn ${done.has(L.id) ? "is-done" : ""}" id="done-btn">${done.has(L.id) ? "✓ 已学完" : "标记已学完"}</button>
        ${next ? `<a class="nav next" href="#/${next.id}"><small>下一课 →</small>${next.title}</a>` : `<a class="nav next" href="#/"><small>→</small>回到首页</a>`}
      </div>`;
    enhance(page, L);
    $("#done-btn").onclick = () => {
      if (done.has(L.id)) done.delete(L.id); else done.add(L.id);
      store.set("cb-done", [...done]);
      renderLessonFootOnly(L); buildSidebar(L.id); updateProgress();
    };
    document.title = L.title + " · 控制之美";
  }
  function renderLessonFootOnly(L) { const b = $("#done-btn"); b.classList.toggle("is-done", done.has(L.id)); b.textContent = done.has(L.id) ? "✓ 已学完" : "标记已学完"; }

  function route() {
    CE.stopAll();
    const id = decodeURIComponent((location.hash || "").replace(/^#\/?/, ""));
    const L = byId[id];
    if (L) renderLesson(L); else renderHome();
    buildSidebar(L ? L.id : "");
    updateProgress();
    window.scrollTo(0, 0);
  }

  function start() {
    lessons.forEach((L, i) => { if (L.ord === undefined) L.ord = i; });
    window.addEventListener("hashchange", route);
    $("#menu-btn").onclick = () => document.body.classList.toggle("sb-open");
    $("#theme-btn").onclick = () => {
      // 默认浅色；深色只在用户手动切换后启用
      const r = document.documentElement;
      r.dataset.theme = r.dataset.theme === "dark" ? "light" : "dark";
      try { localStorage.setItem("cb-theme", r.dataset.theme); } catch (e) { }
      CE.refreshColors(); const y = window.scrollY; route(); window.scrollTo(0, y);
    };
    CE.refreshColors();
    route();
  }

  return { add, start, enhance, lessons, byId };
})();
