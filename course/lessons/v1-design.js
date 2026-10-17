/* 卷1 第7–8章：比例积分控制 / 根轨迹 / 补偿器 / PID */
(function () {
  const { M, P, Plot, Scene, card, linspace, logspace, f2, COL, clamp, rk4, cfmt } = CE;

  /* ════════════════ 第7章 比例积分控制 ════════════════ */
  Course.add({
    id: "v1-pi", vol: 1, ord: 120, no: "7", chName: "第7章 比例积分控制", ref: "卷1 第7章",
    title: "比例积分控制：科学地燃烧卡路里",
    sub: "用体重管理当例子，弄清楚比例控制为什么总差一点、积分控制怎么消除误差、PI 为什么是工业界最爱，以及现实中的“饱和”约束。",
    goals: ["理解比例控制 $u=K_Pe$ 及其稳态误差", "会用终值定理 $\\lim_{t\\to\\infty}x=\\lim_{s\\to0}sX(s)$ 求稳态值", "理解积分控制的原理（让 $C(0)=\\infty$）及其“振荡、迟缓”的代价", "会用 PI 控制和饱和函数处理真实约束"],
    minutes: 35,
    html: String.raw`
<div class="story" data-t="引子：燃烧卡路里">回到第 2 章的体重模型。三个 175 cm、20 岁、体重 70 kg 的大学男生（轻体力 α=1.3）：<b>宅男 A</b> 每天吃 2500 kCal 不运动；<b>宅男 B</b> 戒了可乐和烤串，每天 2100 kCal；<b>运动男 C</b> 吃 2500 kCal 但每天慢跑 1 小时（消耗 500 kCal）。几年后他们的体重分别稳定在 <b>92.4 kg、61.7 kg、54 kg</b>。现在宅男 A 已经胖到 90 kg，想减到 65 kg——该怎么安排饮食与运动？这就是一个控制问题。</div>
<p>令 $C=6.25h-5a+S=998.75$（成年后近似不变），定义<b>输入</b> $u=E_i-E_a$（每日净热量），<b>扰动</b> $d=-\alpha C$，体重方程变为：</p>
<div class="eq" data-no="(7.1.3)–(7.1.5)" data-title="体重系统">
$$7000\,\dot x + 10\alpha\,x = u + d,\qquad G(s)=\frac{1}{7000s+10\alpha},\qquad s_p=-\frac{\alpha}{700}$$
<dl><dt>$x$</dt><dd>体重（kg），时间单位：天</dd><dt>$u$</dt><dd>每日净热量输入（kCal）：吃进去的减去额外运动消耗的</dd><dt>$d$</dt><dd>扰动：由基础代谢中与体重无关的部分决定，$d=-\alpha C\approx-1298$</dd></dl>
</div>

<h2><span class="sec">7.2</span>比例控制</h2>
<p>最直觉的策略：“比目标重得越多，就越少吃、越多运动”——控制量与误差成正比：</p>
<div class="eq" data-no="(7.2.1)" data-title="比例控制器">
$$u(t)=K_P\,e(t)=K_P\,(r-x(t)),\qquad C(s)=K_P$$
</div>
<p>闭环解出（$R=r/s$，$D=-\alpha C/s$）：$x(t)\to x(\infty)=\dfrac{K_Pr-\alpha C}{K_P+10\alpha}$，于是</p>
<div class="eq" data-no="(7.2.7)" data-title="稳态误差">
$$e_{ss}=r-x(\infty)=\frac{10\alpha r+\alpha C}{K_P+10\alpha}$$
<div class="say">K_P 越大误差越小，但只有 K_P→∞ 才为 0。现实中每天的热量不可能无限调节——<b>比例控制无法消除稳态误差</b>。</div>
</div>
<p>直观地说：体重要稳定在某处，就需要一个恒定的净热量来“维持”它；而比例控制只有在误差不为零时才会给出非零输入——所以它<b>必须</b>留一点误差来“养活”这个维持量。</p>

<h2><span class="sec">7.3.1</span>终值定理：不解方程就知道终点</h2>
<div class="eq" data-no="(7.3.1)" data-title="终值定理（前提：x(t) 的极限存在，即系统稳定）">
$$\lim_{t\to\infty}x(t)=\lim_{s\to0}s\,X(s)$$
</div>
<p>对一阶对象 $\frac{a}{s+a}$、单位反馈、阶跃参考，用它可以一步算出一般控制器 C(s) 的稳态误差：</p>
<div class="eq" data-no="(7.3.10)">
$$e_{ss}=\lim_{s\to0}s\cdot\frac{(s+a)\frac{r}{s}}{s+a+aC(s)}=\frac{r}{1+\lim_{s\to0}C(s)}$$
<div class="say">想让 $e_{ss}=0$，就需要 $C(0)=\infty$。最简单的选择：$C(s)=\frac{K_I}{s}$——这就是<b>积分控制</b>，是从数学上“倒推”出来的！</div>
</div>

<h2><span class="sec">7.3.2</span>积分控制</h2>
<div class="eq" data-no="(7.3.11)">
$$u(t)=K_I\int_0^t e(\tau)\,\mathrm{d}\tau,\qquad C(s)=\frac{K_I}{s}$$
<div class="say">只要还有误差，积分就会一直累积、控制量一直变化，直到误差彻底为零才停下来。</div>
</div>
<p>代价是什么？对一阶对象 $\frac{a}{s+a}$ 加积分控制后，闭环变成：</p>
<div class="eq" data-no="(7.3.13)–(7.3.14)">
$$\ddot x + a\dot x + aK_I x = aK_I r\quad\Rightarrow\quad \zeta=\frac{a}{2\sqrt{aK_I}},\quad \omega_n=\sqrt{aK_I}$$
<div class="say">一阶系统变成了二阶系统——会振荡、有超调。$K_I$ 越大越快，但 ζ 越小、超调越大：这是积分控制的调参矛盾。</div>
</div>
<div class="sim" data-id="ki"></div>

<h2><span class="sec">7.3.3</span>比例积分（PI）控制</h2>
<div class="eq" data-no="(7.3.16)" data-title="PI 控制器：工业界最常用的控制器">
$$u(t)=K_P\,e(t)+K_I\int_0^t e(\tau)\,\mathrm{d}\tau,\qquad C(s)=K_P+\frac{K_I}{s}$$
<div class="say">比例项负责“快”，积分项负责“准”。</div>
</div>

<h2><span class="sec">7.4</span>含有限制条件的控制器</h2>
<p>PI 控制器（$K_P=200,K_I=1$）一开始给出的 $u=-5000$ kCal——相当于不吃不喝还要跑十个小时，根本做不到。现实中控制量总有上下限（发动机最大扭矩、空调最大风量……）。办法是在控制器后面加一个<span class="kw">饱和函数</span>：</p>
<div class="eq" data-no="(7.4.1)">
$$u_{\text{实际}}=\begin{cases}u_{\max}, & u>u_{\max}\\ u, & u_{\min}\le u\le u_{\max}\\ u_{\min}, & u<u_{\min}\end{cases}\qquad u_{\min}=1000,\;u_{\max}=5000\text{ kCal}$$
<div class="say">$u_{\min}=1000$：每天至少吃 1800 kCal、最多健身 1.5 小时消耗 800 kCal。</div>
</div>
<div class="sim" data-id="weightctl"></div>
<div class="tip">加了饱和后，前约 250 天每日净热量都卡在下限 1000 kCal，体重缓慢地线性下降——<b>减肥是长期工作，急不得</b>。成功之后维持很简单：每天净摄入约 2150 kCal（比如吃 2500、运动消耗 350）。<br>（书中未展开的一点：饱和期间积分项仍在累积，导致体重先“冲过头”再回来，这叫积分饱和/windup，工程上常用“抗饱和”技巧处理。）</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>比例控制 $u=K_Pe$：简单、快，但有稳态误差。</li>
<li>终值定理 $\lim_{t\to\infty}x=\lim_{s\to0}sX(s)$（要求系统稳定）；阶跃误差 $e_{ss}=\frac{r}{1+C(0)}$。</li>
<li>积分控制 $u=K_I\int e$：消除稳态误差，但引入振荡、反应迟缓。</li>
<li>PI 兼顾快与准；有物理约束时用饱和函数限制控制量。</li>
</ol></div>`,
    sims: {
      ki(el) {
        const ui = card(el, { title: "积分增益：快与稳的矛盾", hint: "对象 $G(s)=\\frac{1}{s+1}$，控制器 $C(s)=\\frac{K_I}{s}$，参考 r=1。右图：两个闭环极点随 $K_I$ 移动。" });
        const [a, b] = ui.row("", "");
        const pr = new Plot(a, { height: 260, x: [0, 15], y: [0, 1.6], xlabel: "t", title: "阶跃响应" });
        const ps = new Plot(b, { height: 260, x: [-1.5, 0.5], y: [-2.2, 2.2], equal: true, xlabel: "σ", title: "闭环极点" });
        const sk = ui.slider({ label: "积分增益 $K_I$", min: 0.05, max: 5, step: 0.05, value: 1, onInput: draw });
        const cP = ui.check({ label: "显示书中 K_I = 1, 2, 3 三条曲线", value: true, onChange: draw });
        function draw() {
          const K = sk.get(), ts = linspace(0, 15, 600), a_ = 1;
          const resp = k => CE.step([a_ * k], [1, a_, a_ * k], 15, 0.025);
          const r = resp(K), z = a_ / (2 * Math.sqrt(a_ * K)), wn = Math.sqrt(a_ * K);
          pr.render = () => { pr.clear(); pr.hline(1, COL.c3, [6, 4]); if (cP.get()) [1, 2, 3].forEach((k, i) => { const q = resp(k); pr.line(q.t, q.y, `hsl(${200 + i * 50},45%,65%)`, 1.3); }); pr.line(r.t, r.y, COL.c1, 2.8); };
          const rt = P.roots([1, a_, a_ * K]);
          ps.render = () => { ps.clear(); const ks = logspace(-2, 0.8, 200); ks.forEach(k => P.roots([1, 1, k]).forEach(p => ps.dot(p[0], p[1], COL.axis, 1.4))); rt.forEach(p => ps.cross(p[0], p[1], COL.c4, 8, 3)); };
          pr.render(); ps.render();
          const info = CE.stepInfo(r.t, r.y, 1);
          ui.out(`ζ = <b>${f2(z)}</b>，ωₙ = <b>${f2(wn)}</b><br>超调 Mp ≈ <b>${f2(info.os, 1)}%</b><br>稳态误差 = <b>0</b>（积分的功劳）<br><span class="sim-note">$K_I=0.25$ 时 ζ=1，恰好临界阻尼。</span>`);
        }
        draw();
      },
      weightctl(el) {
        const ui = card(el, { title: "体重控制实验室", hint: "初始 90 kg，目标 65 kg。切换控制器类型，观察体重（上）和每日净热量 u（下，阴影是“人能做到的范围”）。时间单位：天。" });
        const px = new Plot(ui.plotHost(), { height: 230, x: [0, 1000], y: [30, 100], xlabel: "天", ylabel: "体重 kg" });
        const pu = new Plot(ui.plotHost(), { height: 200, x: [0, 1000], y: [-11000, 6000], xlabel: "天", ylabel: "u kCal/天" });
        let mode = "P";
        ui.seg({ label: "控制器", options: [["open", "开环三案例"], ["P", "比例 P"], ["I", "积分 I"], ["PI", "PI"]], value: mode, onChange: v => { mode = v; draw(); } });
        const sP = ui.slider({ label: "$K_P$", min: 0, max: 500, step: 10, value: 200, fmt: v => v.toFixed(0), onInput: draw });
        const sI = ui.slider({ label: "$K_I$", min: 0, max: 5, step: 0.05, value: 1, onInput: draw });
        const sat = ui.check({ label: "启用饱和 [1000, 5000] kCal", value: false, onChange: draw });
        const alpha = 1.3, Cc = 998.75, r = 65;
        function sim(ctrl, x0, T, uSat) {
          const dt = 0.5; let x = x0, I = 0; const ts = [], xs = [], us = [];
          for (let t = 0; t <= T; t += dt) {
            const e = r - x; let u = ctrl(e, I); if (uSat) u = clamp(u, 1000, 5000);
            ts.push(t); xs.push(x); us.push(u);
            I += e * dt; x += dt * (u - alpha * Cc - 10 * alpha * x) / 7000;
          }
          return { ts, xs, us };
        }
        function draw() {
          const Kp = sP.get(), Ki = sI.get(), T = mode === "I" ? 3500 : mode === "open" ? 2500 : mode === "P" ? 250 : 1000;
          px.setRange([0, T], null); pu.setRange([0, T], null);
          let runs = [];
          if (mode === "open") {
            [[2500, 0, "A：2500/0", COL.c4], [2100, 0, "B：2100/0", COL.c2], [2500, 500, "C：2500/500", COL.c3]].forEach(([Ei, Ea, n, c]) => runs.push({ ...sim(() => Ei - Ea, 70, T, false), n, c }));
            px.setRange(null, [45, 100]); pu.setRange(null, [1500, 2700]);
          } else if (mode === "P") {
            [[Kp, COL.c1, "K_P=" + Kp]].concat([200, 300, 400].filter(k => k !== Kp).map((k, i) => [k, `hsl(${30 + i * 40},50%,65%)`, "K_P=" + k])).forEach(([k, c, n]) => runs.push({ ...sim(e => k * e, 90, T, sat.get()), n, c }));
            px.setRange(null, [45, 95]); pu.setRange(null, sat.get() ? [0, 6000] : [-11000, 3000]);
          } else if (mode === "I") {
            runs.push({ ...sim((e, I) => Ki * I, 90, T, sat.get()), n: "K_I=" + Ki, c: COL.c1 });
            px.setRange(null, [20, 110]); pu.setRange(null, [-3000, 6000]);
          } else {
            runs.push({ ...sim((e, I) => Kp * e + Ki * I, 90, T, sat.get()), n: "PI", c: COL.c1 });
            if (!sat.get()) runs.push({ ...sim((e, I) => Kp * e + Ki * I, 90, T, true), n: "PI + 饱和", c: COL.c5 });
            px.setRange(null, [50, 95]); pu.setRange(null, sat.get() ? [0, 6000] : [-6000, 6000]);
          }
          px.render = () => { px.clear(); if (mode !== "open") px.hline(r, COL.c3, [6, 4]); runs.forEach((q, i) => px.line(q.ts, q.xs, q.c, i ? 1.6 : 2.6)); px.legend(runs.map(q => ({ label: q.n, color: q.c })), "tr"); };
          pu.render = () => { pu.clear(); pu.rect(0, 1000, T, 5000, COL.c3 + "18"); runs.forEach((q, i) => pu.line(q.ts, q.us, q.c, i ? 1.6 : 2.4)); };
          px.render(); pu.render();
          const q = runs[0], xf = q.xs[q.xs.length - 1];
          let msg = `终点体重：<b class="big">${f2(xf, 1)} kg</b>`;
          if (mode === "P") msg += `<br>理论 $e_{ss}=\\frac{10\\alpha r+\\alpha C}{K_P+10\\alpha}$ = <b>${f2((13 * r + alpha * Cc) / (Kp + 13), 2)} kg</b><br>开局 u = ${(q.us[0]).toFixed(0)} kCal`;
          if (mode === "I") msg += `<br>ζ = ${f2(13 / (2 * Math.sqrt(7000 * Ki)), 3)}：严重欠阻尼，振荡周期约 ${f2(2 * Math.PI / Math.sqrt(Ki / 7000), 0)} 天！`;
          if (mode === "PI") msg += `<br>最终每日净热量 ≈ <b>${q.us[q.us.length - 1].toFixed(0)}</b> kCal`;
          if (mode === "open") msg = "A → 92.4 kg，B → 61.7 kg，C → 54 kg<br>时间常数约 538 天";
          ui.out(msg);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>P：K_P 从 200 调到 500，体重越来越接近 65 但总差一点；看看开局的 u——负一万千卡！</li><li>I：能到 65，但要振荡好几年。</li><li>PI：几十天就到，且无误差。再勾选“饱和”，对比现实可行的方案。</li></ol></div>`);
      },
    },
  });

  /* ════════════════ 第8章 根轨迹 ════════════════ */
  // 由根列表构造多项式
  const polyOf = rs => { let c = [1]; rs.forEach(r => { c = r[1] !== 0 ? P.mul(c, [1, -2 * r[0], r[0] * r[0] + r[1] * r[1]]) : P.mul(c, [1, -r[0]]); }); return c; };
  // 根轨迹点集：K 从 0 到很大
  function locus(N, D, kmax = 1e3) {
    const pts = []; const ks = [0, ...logspace(-3, Math.log10(kmax), 700)];
    ks.forEach(k => P.roots(P.add(D, P.scale(N, k))).forEach(r => pts.push(r)));
    return pts;
  }
  const RL_PRESETS = {
    e823: { n: "例 8.2.3", z: [[-2, 0], [-4, 0]], p: [[-1, 0], [-3, 0]] },
    e824: { n: "例 8.2.4", z: [[0, 2]], p: [[-1, 1]] },
    e825: { n: "例 8.2.5a", z: [[-2, 0]], p: [[-1, 0], [-3, 0]] },
    e826: { n: "例 8.2.6", z: [], p: [[-1, 0], [-2, 0]] },
    e827: { n: "例 8.2.7", z: [[-1, 0]], p: [[0, 0], [-2, 0], [-3, 0], [-4, 0]] },
    e828: { n: "例 8.2.8 (ωₙ=2)", z: [[0, 0]], p: [[0, 2]], g: 4 },
  };
  // 共轭扩展：存储的复根只存上半，展开时加共轭
  const expand = rs => { const out = []; rs.forEach(r => { out.push([r[0], r[1]]); if (r[1] !== 0) out.push([r[0], -r[1]]); }); return out; };

  Course.add({
    id: "v1-rlocus", vol: 1, ord: 130, no: "8.1–8.3", chName: "第8章 根轨迹法", ref: "卷1 §8.1–8.3",
    title: "根轨迹：增益变化时，闭环极点怎么走",
    sub: "用开环传递函数的零极点，画出闭环极点随增益 K 移动的轨迹。重点不是手绘技巧，而是理解背后的逻辑。",
    goals: ["知道根轨迹研究的对象：$1+KG(s)=0$ 的根随 K 的变化", "理解 6 条基本规则的含义", "掌握幅值条件与相角条件，会判断一个点是否在根轨迹上"],
    minutes: 35,
    html: String.raw`
<p>上一章我们看到：调整控制器参数，闭环极点就会移动，系统表现随之改变。<span class="kw">根轨迹法</span>把这件事画成一张图：“根”是闭环特征方程的根（闭环极点），“轨迹”是它们随增益 $K$ 从 0 到 ∞ 在复平面上走过的路线。作者特别强调：今天软件一两行代码就能画根轨迹，<b>理解规律背后的逻辑比手绘技巧重要得多</b>。</p>
<div data-fig="closedloop"></div>
<div class="eq" data-no="8.1" data-title="研究对象：单位反馈 + 比例增益 K">
$$G_{cl}(s)=\frac{KG(s)}{1+KG(s)},\qquad \text{闭环极点满足}\;1+KG(s)=0\;\Longleftrightarrow\; D(s)+K\,N(s)=0$$
<dl><dt>$G(s)=\frac{N(s)}{D(s)}$</dt><dd>开环传递函数；$N$ 的根是零点（○），$D$ 的根是极点（×）</dd><dt>$n,\ m$</dt><dd>开环极点数、零点数</dd></dl>
<div class="say">研究<b>目标</b>是闭环极点，研究<b>方法</b>是分析开环 G(s)。所以先要把闭环分母整理成 $1+KG(s)$ 的样子。例：$G_{cl}=\frac{Ks}{s^3+3s^2+Ks+1}$ 整理后 $G(s)=\frac{s}{s^3+3s^2+1}$；非单位反馈时用 $G(s)H(s)$。</div>
</div>

<h2><span class="sec">8.2</span>六条基本规则</h2>
<table class="tbl">
<tr><th>规则</th><th>内容</th><th>为什么</th></tr>
<tr><td>1 分支数</td><td>共有 $\max(n,m)$ 条分支</td><td>$D+KN=0$ 的次数</td></tr>
<tr><td>2 起点终点</td><td>K: 0→∞，从开环<b>极点</b>出发，走向开环<b>零点</b></td><td>K=0 时 $D=0$；K→∞ 时 $N\to0$</td></tr>
<tr><td>3 实轴部分</td><td>实轴上，从右往左数第<b>奇数</b>个极点/零点的<b>左边</b>一段有根轨迹</td><td>由相角条件推出</td></tr>
<tr><td>4 对称</td><td>关于实轴对称</td><td>复根总是共轭成对</td></tr>
<tr><td>5 去向无穷</td><td>$n>m$ 时有 $n-m$ 条分支走向无穷</td><td>零点“不够用”</td></tr>
<tr><td>6 渐近线</td><td>交点 $\sigma_a=\dfrac{\sum p-\sum z}{n-m}$，夹角 $\theta=\dfrac{(2q+1)\pi}{n-m}$</td><td>远处看，零极点像集中在一点</td></tr>
</table>
<div class="sim" data-id="rlboard"></div>
<p><b>分离点与汇合点</b>：两条分支在实轴上相遇后同时离开实轴的点叫分离点，从复平面回到实轴的点叫汇合点。书中例 8.2.8 很巧妙：把二阶系统 $\frac{2\zeta\omega_ns}{s^2+2\zeta\omega_ns+\omega_n^2}$ 中的 ζ 当成 K，开环就是 $G=\frac{2\omega_ns}{s^2+\omega_n^2}$。把 $K=-\frac{s^2+\omega_n^2}{2\omega_ns}$ 在实轴上求极值，得汇合点 $\sigma=-\omega_n$、$K=1$——正是临界阻尼 ζ=1！ζ&lt;1 时极点沿半径 ωₙ 的半圆运动（欠阻尼），ζ&gt;1 时在实轴上分开（过阻尼）。<b>分析分离/汇合点时，关键是它对应的 K 值</b>：知道 K 多大时极点进入复平面开始振荡，就能避免或利用振荡。</p>

<h2><span class="sec">8.3</span>几何性质：幅值条件与相角条件</h2>
<p>复数相乘：模相乘、角相加；相除：模相除、角相减。所以 $G(s)$ 在某点的值可以“量”出来：</p>
<div class="eq" data-no="(8.3.2)">
$$|G(s)|=\frac{\prod\text{各零点到 }s\text{ 的距离}}{\prod\text{各极点到 }s\text{ 的距离}},\qquad \angle G(s)=\sum\text{零点到 }s\text{ 的角}-\sum\text{极点到 }s\text{ 的角}$$
</div>
<p>闭环极点满足 $KG(s)=-1$，而 −1 的模是 1、角度是 π，所以：</p>
<div class="eq" data-no="(8.3.6)–(8.3.7)" data-title="根轨迹的两个条件">
$$\angle G(s)=-(2q+1)\pi\quad(\text{相角条件：判断是否在根轨迹上}),\qquad |KG(s)|=1\;\Rightarrow\;K=\frac{1}{|G(s)|}\quad(\text{幅值条件：求对应的 }K)$$
</div>
<div class="sim" data-id="protractor"></div>
<p>例 8.3.2：对 $G=\frac{1}{(s+1)(s+3)}$，点 $-2-j$ 的相角正好是 π，在根轨迹上（K=2）；点 $-3+2\sqrt3j$ 的相角是 $-\frac76\pi$，<b>不在</b>根轨迹上——不管 K 怎么调，闭环极点都到不了那里。想让它到达，就得改造根轨迹本身：加补偿器。下一课继续。</p>

<div class="summary"><h3>本课要点</h3><ol>
<li>根轨迹：K 从 0 到 ∞ 时 $1+KG(s)=0$ 的根的轨迹；从开环极点出发，走向开环零点或无穷远。</li>
<li>六条规则：分支数、起终点、实轴奇数规则、对称、去向无穷、渐近线。</li>
<li>相角条件 $\angle G=-(2q+1)\pi$ 判断点是否在轨迹上；幅值条件 $K=1/|G|$ 求增益。</li>
</ol></div>`,
    sims: {
      rlboard(el) {
        const ui = card(el, { title: "交互式根轨迹画板", hint: "选择书中例题，或<b>拖动</b>极点（×）和零点（○）。拖动 K 滑块看闭环极点（红色大点）的位置。可逐条打开规则提示。" });
        const ps = new Plot(ui.plotHost(), { height: 400, x: [-7, 2], y: [-4.5, 4.5], equal: true, xlabel: "σ", ylabel: "jω" });
        let Z = [], Pp = [], gain = 1;
        const sel = document.createElement("select"); sel.style.cssText = "padding:6px;border-radius:8px;border:1px solid var(--line);background:var(--panel);color:var(--ink);font:inherit;font-size:.85rem";
        sel.innerHTML = Object.entries(RL_PRESETS).map(([k, v]) => `<option value="${k}">${v.n}</option>`).join("");
        sel.onchange = () => load(sel.value); ui.controls.appendChild(sel);
        const sK = ui.slider({ label: "增益 K", min: 0.001, max: 200, value: 1, log: true, fmt: v => v < 1 ? v.toFixed(3) : v.toFixed(2), onInput: draw });
        const rules = { r3: ui.check({ label: "规则3：实轴上的根轨迹", value: true, onChange: draw }), r6: ui.check({ label: "规则6：渐近线", value: true, onChange: draw }) };
        ui.buttons([["+ 极点", () => { Pp.push([-1.5, 0]); draw(); }, "ghost"], ["+ 零点", () => { Z.push([-5, 0]); draw(); }, "ghost"], ["+ 共轭极点", () => { Pp.push([-1, 1.5]); draw(); }, "ghost"]]);
        let drag = null;
        ps.drag({
          down: (x, y) => { drag = null; let bd = 16; [["z", Z], ["p", Pp]].forEach(([t, arr]) => arr.forEach((r, i) => { [r[1], -r[1]].forEach(im => { const d = ps.pxDist(r[0], im, x, y); if (d < bd) { bd = d; drag = [arr, i]; } }); })); return drag ? undefined : false; },
          move: (x, y) => { const [arr, i] = drag; const r = arr[i]; arr[i] = [clamp(x, -6.8, 1.8), r[1] !== 0 ? clamp(Math.abs(y), 0.15, 4.3) : 0]; draw(); },
          hover: (x, y) => [...Z, ...Pp].some(r => ps.pxDist(r[0], r[1], x, y) < 16 || ps.pxDist(r[0], -r[1], x, y) < 16),
        });
        function load(k) { gain = RL_PRESETS[k].g || 1; Z = RL_PRESETS[k].z.map(r => r.slice()); Pp = RL_PRESETS[k].p.map(r => r.slice()); draw(); }
        function draw() {
          const N = P.scale(polyOf(Z), gain), D = polyOf(Pp), K = sK.get(), z = expand(Z), p = expand(Pp), n = p.length, m = z.length;
          const pts = locus(N, D, 1e4), cl = P.roots(P.add(D, P.scale(N, K)));
          ps.render = () => {
            ps.clear();
            if (rules.r3.get()) {
              const reals = [...z, ...p].filter(r => r[1] === 0).map(r => r[0]).sort((a, b) => b - a);
              reals.forEach((v, i) => { if (i % 2 === 0) { const left = i + 1 < reals.length ? reals[i + 1] : -50; ps.line([v, left], [0, 0], COL.c3 + "aa", 7); } });
            }
            if (rules.r6.get() && n !== m) {
              const sa = (p.reduce((s, r) => s + r[0], 0) - z.reduce((s, r) => s + r[0], 0)) / (n - m), k = Math.abs(n - m);
              for (let q = 0; q < k; q++) { const th = (2 * q + 1) * Math.PI / k; ps.line([sa, sa + 20 * Math.cos(th)], [0, 20 * Math.sin(th)], COL.c5, 1.3, [6, 4]); }
              ps.dot(sa, 0, COL.c5, 4); ps.text(sa, -0.45, "σa=" + f2(sa), COL.c5, "center");
            }
            pts.forEach(r => ps.dot(r[0], r[1], COL.c1 + "99", 1.5));
            z.forEach(r => ps.circle(r[0], r[1], COL.c5, 7, 2.6));
            p.forEach(r => ps.cross(r[0], r[1], COL.ink, 7, 2.6));
            cl.forEach(r => ps.dot(r[0], r[1], COL.c4, 6));
          };
          ps.render();
          const mx = Math.max(...cl.map(r => r[0]));
          ui.out(`n = ${n} 个极点，m = ${m} 个零点 → <b>${Math.max(n, m)}</b> 条分支${n > m ? `，${n - m} 条去往无穷` : ""}<br>K = ${f2(K, 3)} 时闭环极点：<br>${cl.map(r => cfmt(r)).join("，")}<br>${mx < 0 ? '<span class="good">稳定</span>' : '<span class="bad">不稳定</span>'}`);
        }
        load("e826");
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>例 8.2.6：两个极点相向而行，在 −1.5 处分离，沿竖直渐近线上下走。</li><li>例 8.2.7：数一数实轴上哪几段是绿色的（从右往左第 1、3、5 个点的左边）。</li><li>例 8.2.8：K 就是阻尼比 ζ。K=1 时两极点在 −2 汇合（临界阻尼）。</li><li>把例 8.2.6 中的一个极点拖到右半平面，观察小 K 时系统不稳定、大 K 时被“拉回来”。</li></ol></div>`);
      },
      protractor(el) {
        const ui = card(el, { title: "相角条件“量角器”", hint: "对象 $G(s)=\\frac{1}{(s+1)(s+3)}$。<b>拖动测试点 s</b>（黄色），从每个极点画向量到 s，计算 ∠G(s)。相角等于 ±π 时点变绿，并给出所需的 K。" });
        const ps = new Plot(ui.plotHost(), { height: 360, x: [-6, 1], y: [-4, 4], equal: true, xlabel: "σ", ylabel: "jω" });
        let s = [-2, -1];
        const poles = [[-1, 0], [-3, 0]];
        ui.buttons([["−2 − j", () => { s = [-2, -1]; draw(); }, "ghost"], ["−3 + 2√3 j", () => { s = [-3, 2 * Math.sqrt(3)]; draw(); }, "ghost"]]);
        ps.drag({ down: (x, y) => { s = [x, y]; draw(); }, move: (x, y) => { s = [clamp(x, -5.9, 0.9), clamp(y, -3.9, 3.9)]; draw(); }, hover: () => true });
        function draw() {
          let ang = 0, mag = 1;
          poles.forEach(p => { ang -= Math.atan2(s[1] - p[1], s[0] - p[0]); mag /= Math.hypot(s[0] - p[0], s[1] - p[1]); });
          let a = ang; while (a <= -Math.PI) a += 2 * Math.PI; while (a > Math.PI) a -= 2 * Math.PI;
          const on = Math.abs(Math.abs(a) - Math.PI) < 0.03;
          ps.render = () => {
            ps.clear(); ps.line([-2, -2], [-10, 10], COL.c1 + "55", 4); ps.line([-3, -1], [0, 0], COL.c1 + "55", 4);
            poles.forEach((p, i) => { ps.arrow(p[0], p[1], s[0], s[1], i ? COL.c2 : COL.c5, 1.8, 9); const th = Math.atan2(s[1] - p[1], s[0] - p[0]); ps.text(p[0] + 0.25, p[1] + (th > 0 ? 0.3 : -0.45), "φ=" + f2(th * 180 / Math.PI, 0) + "°", i ? COL.c2 : COL.c5); ps.cross(p[0], p[1], COL.ink, 7, 2.6); });
            ps.dot(s[0], s[1], on ? COL.c3 : "#e8b400", 8, COL.ink);
          };
          ps.render();
          ui.out(`s = ${cfmt(s)}<br>∠G(s) = −Σφ = <b>${f2(a * 180 / Math.PI, 1)}°</b><br>|G(s)| = ${f2(mag, 3)}<br>${on ? `<span class="good">✔ 在根轨迹上，K = 1/|G| = ${f2(1 / mag, 2)}</span>` : `<span class="bad">不在根轨迹上，相角差 ${f2((Math.PI - Math.abs(a)) * 180 / Math.PI, 1)}°</span>`}`);
        }
        draw();
      },
    },
  });

  Course.add({
    id: "v1-comp", vol: 1, ord: 140, no: "8.4–8.5", chName: "第8章 根轨迹法", ref: "卷1 §8.4–8.6",
    title: "补偿器设计：PD、超前、滞后与 PID",
    sub: "加零点、加极点，就能把根轨迹“掰”到想要的地方。超前补偿器提速，滞后补偿器减小稳态误差，PID 集大成。",
    goals: ["会用相角条件算出 PD 控制器的零点位置", "理解超前补偿器（零点更靠近虚轴）为什么能提速", "理解滞后补偿器（极点更靠近虚轴）为什么能减小稳态误差，以及为什么要贴近原点", "理解 PID 的“现在、过去、未来”"],
    minutes: 35,
    html: String.raw`
<p>对象仍是 $G(s)=\frac{1}{(s+1)(s+3)}$。单纯比例控制时，根轨迹在 −2 处分离后沿竖线 $\sigma=-2$ 走——无论 K 多大，闭环极点实部最多到 −2，系统最快也只能按 $e^{-2t}$ 收敛。<b>想更快，就得改变根轨迹本身</b>：给开环串联一个<span class="kw">补偿器</span>，添加零点或极点。</p>
<p><b>主导极点</b>：$K=0.5$ 时闭环极点为 −1.29 和 −2.71，响应 $x=c+C_1e^{-1.29t}+C_2e^{-2.71t}$。离虚轴更近的 −1.29 衰减得更慢，决定了整体收敛速度，称为主导极点。</p>

<h2><span class="sec">8.4.1</span>比例微分（PD）控制：加一个零点</h2>
<p>设计目标：把闭环极点放在 $-3\pm2\sqrt3j$（按 $e^{-3t}$ 收敛）。上一课算过它的相角是 $-\frac76\pi$，差 $\frac16\pi$。补一个零点，让它到目标点的角度恰为 $\frac\pi6$：由 $\tan\frac\pi6=\frac{2\sqrt3}{6}$，零点应在 −3 左边 6 个单位，即 <b>$s_z=-9$</b>。</p>
<div class="eq" data-no="(8.4.5)" data-title="PD 控制器">
$$C(s)=s+9,\qquad u(t)=K\left(\frac{\mathrm{d}e}{\mathrm{d}t}+9e\right)$$
<div class="say">微分项看的是误差的变化趋势，能“预测”并提前反应，所以更快。（本页补充：此时闭环特征式 $s^2+(4+K)s+3+9K$，取 K=2 恰好得到 $-3\pm2\sqrt3j$。）</div>
</div>
<p>但 PD 有两个毛病：需要主动能源（无源元件做不出纯微分），而且会<b>放大高频噪声</b>（第 9 章细说）。</p>

<h2><span class="sec">8.4.2</span>超前补偿器</h2>
<div class="eq" data-no="(8.4.6)" data-title="超前补偿器（零点比极点更靠近虚轴）">
$$C(s)=\frac{s-s_{zc}}{s-s_{pc}},\qquad s_{pc}<s_{zc}<0,\qquad \text{设计条件：}\;\varphi_z-\varphi_p=\frac{\pi}{6}$$
<div class="say">渐近线交点变为 $\sigma_a=-2+\frac{s_{pc}-s_{zc}}{2}<-2$，被“拉”离虚轴——响应变快。（本页补充示例：$C=\frac{s+5}{s+9}$，K=24。）</div>
</div>

<h2><span class="sec">8.4.3</span>滞后补偿器</h2>
<p>用终值定理算加了补偿器后的阶跃稳态误差：</p>
<div class="eq" data-no="(8.4.10)" data-title="滞后补偿器（极点比零点更靠近虚轴）">
$$e_{ss}=\frac{D(0)}{D(0)+K\,N(0)\,\dfrac{s_{zc}}{s_{pc}}},\qquad s_{zc}<s_{pc}<0$$
<div class="say">$s_{zc}/s_{pc}$ 越大，误差越小。极限情况 $s_{pc}=0$：$C=\frac{s-s_{zc}}{s}$——正是 PI 控制器！</div>
</div>
<p>例：K=2 时原系统 $e_{ss}=\frac{1}{1+2/3}=0.6$。取 $s_{zc}/s_{pc}=5$ 可降到 $\frac{3}{13}\approx0.23$。两种选择 $C_1=\frac{s+5}{s+1}$ 和 $C_2=\frac{s+0.5}{s+0.1}$ 效果一样吗？</p>
<div class="sim" data-id="compbench"></div>
<div class="idea">$C_2$ 的零极点都紧贴原点，对原根轨迹几乎没有影响：瞬态和原来一样，然后由一个很慢的极点（≈ −0.31）把输出“慢慢拉高”到 0.77。$C_1$ 却彻底改变了根轨迹形状，带来明显的超调振荡。<b>实践中滞后补偿器的零极点要选得靠近原点</b>。</div>

<h2><span class="sec">8.5</span>PID 控制器</h2>
<div class="story" data-t="浴缸里的 PID">洗澡调水温时，旋钮到水温之间隔着阀门、冷热水混合、出水温度、和浴缸里原有的水混合……没人会去建这个模型。但人人都会调：现在太冷就往热拧（当前误差）；冷了好一阵了就多拧一点（过去误差的累积）；水温正在快速变热就提前往回收一点（未来的趋势）。PID 正是源自对水手掌舵的类似观察。</div>
<div class="eq" data-no="(8.5.1)–(8.5.2)" data-title="PID 控制器">
$$u(t)=\underbrace{K_Pe(t)}_{\text{现在}}+\underbrace{K_I\int e\,\mathrm{d}t}_{\text{过去}}+\underbrace{K_D\frac{\mathrm{d}e}{\mathrm{d}t}}_{\text{未来}},\qquad C(s)=K_P+\frac{K_I}{s}+K_Ds$$
<div class="say">即使完全不知道系统内部结构，大部分情况下 PID 也能调出满意的效果——这是它统治工业界的原因。</div>
</div>
<div class="sim" data-id="pid"></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>补偿器通过增加零点/极点改造根轨迹；用相角条件计算零极点位置。</li>
<li>PD/超前补偿（零点更靠近虚轴）提高响应速度；滞后补偿（极点更靠近虚轴）减小稳态误差，零极点应靠近原点以免破坏瞬态。</li>
<li>PI 是极点在原点的滞后补偿；PID = 现在 + 过去 + 未来。</li>
</ol></div>`,
    sims: {
      compbench(el) {
        const ui = card(el, { title: "补偿器设计工作台", hint: "对象 $G=\\frac{1}{(s+1)(s+3)}$，红色菱形是目标点 $-3\\pm2\\sqrt3j$。选择补偿器类型并调参数，看根轨迹和阶跃响应。" });
        const [a, b] = ui.row("", "");
        const ps = new Plot(a, { height: 330, x: [-12, 2], y: [-7, 7], equal: true, xlabel: "σ", title: "根轨迹" });
        const pr = new Plot(b, { height: 330, x: [0, 10], y: [0, 1.4], xlabel: "t", title: "单位阶跃响应" });
        let type = "pd";
        const presets = { none: { K: 2 }, pd: { K: 2, z: -9 }, lead: { K: 24, z: -5, p: -9 }, lag1: { K: 2, z: -5, p: -1 }, lag2: { K: 2, z: -0.5, p: -0.1 } };
        const sK = ui.slider({ label: "增益 K", min: 0.1, max: 60, value: 2, log: true, onInput: draw });
        const sz = ui.slider({ label: "补偿器零点 $s_{zc}$", min: -12, max: -0.05, step: 0.05, value: -9, onInput: draw });
        const sp = ui.slider({ label: "补偿器极点 $s_{pc}$", min: -12, max: -0.01, step: 0.01, value: -9, onInput: draw });
        ui.seg({ label: "补偿器", options: [["none", "无 (P)"], ["pd", "PD"], ["lead", "超前"], ["lag1", "滞后 C₁"], ["lag2", "滞后 C₂"]], value: type, onChange: v => { type = v; const q = presets[v]; sK.set(q.K); if (q.z) sz.set(q.z); if (q.p) sp.set(q.p); draw(); } });
        const target = [-3, 2 * Math.sqrt(3)];
        function comp() {
          if (type === "none") return { N: [1], D: [1] };
          if (type === "pd") return { N: [1, -sz.get()], D: [1] };
          return { N: [1, -sz.get()], D: [1, -sp.get()] };
        }
        function draw() {
          const { N: Nc, D: Dc } = comp(), K = sK.get(), G0 = [1, 4, 3];
          const N = Nc, D = P.mul(Dc, G0);
          const pts = locus(N, D, 3e3), cl = P.roots(P.add(D, P.scale(N, K)));
          const num = P.scale(N, K), den = P.add(D, num), r = CE.step(num, den, 10, 0.01), r0 = CE.step([2], [1, 4, 5], 10, 0.01);
          const ess = 1 - r.y[r.y.length - 1];
          ps.render = () => {
            ps.clear(); pts.forEach(q => ps.dot(q[0], q[1], COL.c1 + "88", 1.5));
            ps.cross(-1, 0, COL.ink, 6, 2.4); ps.cross(-3, 0, COL.ink, 6, 2.4);
            if (type !== "none") { ps.circle(sz.get(), 0, COL.c5, 7, 2.6); if (type !== "pd") ps.cross(sp.get(), 0, COL.c5, 7, 2.6); }
            [1, -1].forEach(sg => { const X = ps.fx(target[0]), Y = ps.fy(sg * target[1]), c = ps.ctx; c.strokeStyle = COL.c4; c.lineWidth = 2; c.beginPath(); c.moveTo(X, Y - 8); c.lineTo(X + 8, Y); c.lineTo(X, Y + 8); c.lineTo(X - 8, Y); c.closePath(); c.stroke(); });
            cl.forEach(q => ps.dot(q[0], q[1], COL.c4, 6));
          };
          pr.render = () => { pr.clear(); pr.hline(1, COL.c3, [6, 4]); pr.line(r0.t, r0.y, COL.muted, 1.4, [5, 4]); pr.line(r.t, r.y, COL.c1, 2.6); pr.legend([{ label: "当前设计", color: COL.c1 }, { label: "原系统 K=2", color: COL.muted, dash: [5, 4] }], "tr"); };
          ps.render(); pr.render();
          // 目标点相角
          let ang = -Math.atan2(target[1], target[0] + 1) - Math.atan2(target[1], target[0] + 3);
          if (type !== "none") ang += Math.atan2(target[1], target[0] - sz.get());
          if (type === "lead" || type.startsWith("lag")) ang -= Math.atan2(target[1], target[0] - sp.get());
          const sa = type === "pd" ? null : type === "none" ? -2 : (-1 - 3 + sp.get() - sz.get()) / 2;
          ui.out(`闭环极点：${cl.map(q => cfmt(q)).join("，")}<br>目标点相角 ∠CG = <b>${f2(ang * 180 / Math.PI, 1)}°</b>（需 −180°）<br>${sa !== null ? `渐近线 σa = ${f2(sa)}<br>` : ""}稳态误差 e_ss ≈ <b>${f2(ess, 3)}</b>`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>PD：零点拖到 −9 附近，根轨迹变成以 −9 为圆心的圆，正好穿过目标点（K=2）。</li><li>超前：保持 φz−φp=30°（如 −5/−9），渐近线左移到 −4。</li><li>滞后 C₁ 与 C₂：两者 e_ss 都降到 ≈0.23，但瞬态天差地别。</li><li>把滞后补偿器的极点拖到 −0.01（接近 0）：几乎变成 PI，误差趋于 0。</li></ol></div>`);
      },
      pid(el) {
        const ui = card(el, { title: "PID 调参：现在、过去、未来", hint: "对象：三阶系统 $G(s)=\\frac{1}{(s+1)^3}$（模拟“旋钮→水温”这种中间环节很多的过程）。下图分别显示 P、I、D 三项各自的贡献。" });
        const pr = new Plot(ui.plotHost(), { height: 230, x: [0, 20], y: [0, 1.6], xlabel: "t", ylabel: "水温（归一化）" });
        const pu = new Plot(ui.plotHost(), { height: 200, x: [0, 20], y: [-2, 4], xlabel: "t", ylabel: "u 的组成" });
        const sP = ui.slider({ label: "$K_P$（现在）", min: 0, max: 6, step: 0.05, value: 1.2, onInput: draw });
        const sI = ui.slider({ label: "$K_I$（过去）", min: 0, max: 2, step: 0.01, value: 0.4, onInput: draw });
        const sD = ui.slider({ label: "$K_D$（未来）", min: 0, max: 4, step: 0.05, value: 0.8, onInput: draw });
        const sN = ui.slider({ label: "测量噪声", min: 0, max: 0.05, step: 0.002, value: 0, fmt: v => v.toFixed(3), onInput: draw });
        function draw() {
          const Kp = sP.get(), Ki = sI.get(), Kd = sD.get(), nz = sN.get(), dt = 0.01, rg = CE.rng(7);
          let x = [0, 0, 0], I = 0, ePrev = 1; const ts = [], ys = [], uP = [], uI = [], uD = [];
          for (let t = 0; t <= 20; t += dt) {
            const y = x[0] + nz * rg.n(), e = 1 - y; I += e * dt; const de = (e - ePrev) / dt; ePrev = e;
            const p = Kp * e, i = Ki * I, d = t === 0 ? 0 : Kd * de, u = p + i + d;
            ts.push(t); ys.push(x[0]); uP.push(p); uI.push(i); uD.push(d);
            // (s+1)^3 链式一阶环节
            x = [x[0] + dt * (x[1] - x[0]), x[1] + dt * (x[2] - x[1]), x[2] + dt * (u - x[2])];
          }
          const mx = Math.max(...ys.map(Math.abs)); pr.setRange(null, [0, Math.max(1.5, Math.min(mx * 1.1, 4))]);
          pr.render = () => { pr.clear(); pr.hline(1, COL.c3, [6, 4]); pr.line(ts, ys, COL.c1, 2.6); };
          const cl = v => clamp(v, -6, 8);
          pu.render = () => { pu.clear(); pu.line(ts, uP.map(cl), COL.c2, 1.8); pu.line(ts, uI.map(cl), COL.c3, 1.8); pu.line(ts, uD.map(cl), COL.c5, 1.4); pu.legend([{ label: "P 现在", color: COL.c2 }, { label: "I 过去", color: COL.c3 }, { label: "D 未来", color: COL.c5 }], "tr"); };
          pr.render(); pu.render();
          const info = CE.stepInfo(ts, ys, 1);
          ui.out(`超调 ${f2(info.os, 1)}%<br>稳定时间 ${info.ts > 19.9 ? "&gt;20" : f2(info.ts, 1)}<br>终值 ${f2(ys[ys.length - 1], 3)}`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>只开 P（K_I=K_D=0）：有稳态误差；P 太大会振荡。</li><li>加 I：误差消除，但太大会超调振荡。</li><li>加 D：阻尼振荡、减小超调。然后加一点测量噪声——D 项（紫线）剧烈抖动：微分放大噪声。</li></ol></div>`);
      },
    },
  });
})();
