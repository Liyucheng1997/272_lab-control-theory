/* 卷1 第10章 + 附录：能控性 / 极点配置 / LQR / 轨迹追踪 / 观测器 / 线性化 / 傅里叶 */
(function () {
  const { M, P, Plot, Scene, card, linspace, logspace, f2, COL, clamp, rk4, cfmt, mtex } = CE;
  const GD = 10; // g/d
  const A0 = [[0, 1], [GD, 0]], B0 = [[0], [1]], C0 = [[1, 0]];

  // 画“指尖上的平衡”：手指 + 杆 + 球，xi 为手指位置（相机跟随）
  function drawFinger(sc, phi, xi, opt = {}) {
    const c = sc.ctx, W = sc.W, H = sc.H; sc.clear();
    const base = H - 30, L = Math.min(H - 60, 150), cx = W / 2, scale = 60;
    // 地面刻度（随手指位置滚动）
    c.strokeStyle = COL.axis; c.lineWidth = 1;
    for (let k = -20; k <= 20; k++) { const x = cx + (k - (xi % 1)) * scale; if (x < 0 || x > W) continue; c.beginPath(); c.moveTo(x, base + 8); c.lineTo(x, base + 14); c.stroke(); }
    c.beginPath(); c.moveTo(0, base + 8); c.lineTo(W, base + 8); c.stroke();
    if (opt.target !== undefined) { c.strokeStyle = COL.c3; c.setLineDash([5, 4]); c.beginPath(); c.moveTo(cx, base); c.lineTo(cx + L * Math.sin(opt.target), base - L * Math.cos(opt.target)); c.stroke(); c.setLineDash([]); }
    // 手指
    c.fillStyle = "#e9b48a"; c.beginPath(); c.roundRect(cx - 9, base, 18, 40, 8); c.fill();
    // 杆与球
    const bx = cx + L * Math.sin(phi), by = base - L * Math.cos(phi);
    c.strokeStyle = COL.ink2; c.lineWidth = 3; c.beginPath(); c.moveTo(cx, base); c.lineTo(bx, by); c.stroke();
    c.fillStyle = Math.abs(phi) > 1.2 ? COL.c4 : COL.c1; c.beginPath(); c.arc(bx, by, 13, 0, 7); c.fill();
    c.fillStyle = COL.muted; c.font = "12px Consolas"; c.textAlign = "left";
    c.fillText(`φ = ${(phi * 180 / Math.PI).toFixed(1)}°   手指位置 ξ = ${xi.toFixed(2)}`, 8, 16);
    if (opt.u !== undefined) { c.fillStyle = COL.c2; const ax = -opt.u * 4; c.fillRect(cx, base + 22, clamp(ax, -W / 2, W / 2), 6); c.fillText("手指加速度", cx + 6, base + 40 > H ? H - 2 : base + 22); }
  }

  /* ════════════════ 10.1–10.2 指尖平衡与能控性 ════════════════ */
  Course.add({
    id: "v1-ctrb", vol: 1, ord: 170, no: "10.1–10.2", chName: "第10章 状态空间控制器与观测器", ref: "卷1 §10.1–10.2",
    title: "指尖上的平衡与能控性",
    sub: "用手指顶住一根竖杆——一个天生不稳定的系统。PID 为什么难调？状态空间方法要先回答一个问题：它“能控”吗？",
    goals: ["建立指尖平衡（简化倒立摆）的模型 $\\ddot\\phi-\\frac gd\\phi=u$", "理解为什么只用输出反馈（P/PID）很难稳定它", "掌握能控性的含义与判据 $\\mathrm{rank}[B\\;AB\\cdots A^{n-1}B]=n$"],
    minutes: 25,
    html: String.raw`
<div class="story" data-t="引子：指尖上的平衡">用手指顶住一根长 d 的轻杆，杆顶有个质量为 m 的小球。你需要左右移动手指，让杆的角度 φ 跟随一个目标：0–10 秒向右倾 π/20（9°），10–20 秒竖直，20–30 秒向左倾 π/20，之后再竖直。这本质上就是著名的<b>倒立摆</b>问题。</div>
<p>受力分析：手指沿杆给小球一个力 F。竖直方向小角度下 $F\approx mg$；水平方向 $x_1=\xi+d\sin\phi$，再用 $\sin\phi\approx\phi$ 线性化（附录 A），得到：</p>
<div class="eq" data-no="(10.1.7)–(10.1.9)" data-title="指尖平衡的模型">
$$\ddot\phi-\frac gd\phi=u,\qquad u=-\frac1d\ddot\xi,\qquad G(s)=\frac{1}{s^2-\frac gd}$$
<dl><dt>$\phi$</dt><dd>杆与竖直方向的夹角（顺时针为正）</dd><dt>$\xi$</dt><dd>手指的水平位置；$u$ 正比于手指的加速度</dd><dt>$g/d$</dt><dd>本课取 10</dd></dl>
<div class="say">G(s) 有一个正极点 $+\sqrt{g/d}$：不加控制，杆一定会倒。</div>
</div>
<h3>PID 为什么难？</h3>
<p>纯比例控制的根轨迹：两个极点 $\pm\sqrt{g/d}$ 沿实轴相向而行，在原点汇合后沿虚轴上下发散——<b>无论 K 多大都只能得到等幅振荡</b>，永远稳定不了。要稳定就得加零点（微分）把根轨迹拉到左边，还要加积分消误差，于是变成 PID。但由于存在正极点，需要足够大的增益才能把它拉回来，而增益大又带来振荡、超调和过大的控制量。根本原因在于：PID 只用了输出误差，<b>没有用到全部状态信息</b>。</p>
<div class="sim" data-id="nocontrol"></div>

<h3>状态空间模型</h3>
<div class="eq" data-no="(10.1.13)">
$$\bm z=\begin{bmatrix}\phi\\\dot\phi\end{bmatrix},\qquad \dot{\bm z}=\begin{bmatrix}0&1\\\frac gd&0\end{bmatrix}\bm z+\begin{bmatrix}0\\1\end{bmatrix}u,\qquad y=\begin{bmatrix}1&0\end{bmatrix}\bm z$$
<div class="say">特征值 $\pm\sqrt{g/d}$：一正一负，平衡点是鞍点——偏一点就回不来。</div>
</div>

<h2><span class="sec">10.2</span>能控性</h2>
<p>设计控制器前，先要问：通过输入，能不能把系统的<b>每一个</b>状态都送到任意想要的值？比如两辆用弹簧连着的小车，你只能推第一辆——能让两辆车同时停在指定位置吗？</p>
<div class="eq" data-no="(10.2.3)" data-title="能控性判据">
$$C_o=\begin{bmatrix}B & AB & A^2B & \cdots & A^{n-1}B\end{bmatrix},\qquad \text{系统能控}\iff \mathrm{rank}(C_o)=n$$
<dl><dt>定义</dt><dd>存在输入，能在有限时间内把状态从任意 $\bm z(t_0)$ 转移到任意 $\bm z(t_1)$</dd><dt>直觉</dt><dd>$B$：输入直接推动的方向；$AB$：经过系统“传递”一步后能影响的方向……这些方向要能张满整个状态空间</dd></dl>
</div>
<ul>
<li>例 10.2.1：$A=\begin{bmatrix}2&0\\1&1\end{bmatrix},B=\begin{bmatrix}1\\1\end{bmatrix}$，$C_o=\begin{bmatrix}1&2\\1&2\end{bmatrix}$，秩 1 → <b>不能控</b>。</li>
<li>例 10.2.2（单辆小车）：$C_o=\begin{bmatrix}0&1\\1&0\end{bmatrix}$，秩 2 → 能控。</li>
<li>例 10.2.3（两车 + 弹簧，$m_1=m_2=1$，$k=100$）：4×4 的 $C_o$ 秩为 4 → 能控。推第一辆车的力会通过弹簧传给第二辆。</li>
<li>指尖平衡：$C_o=\begin{bmatrix}0&1\\1&0\end{bmatrix}$，秩 2 → 能控。</li>
</ul>
<div class="sim" data-id="carts"></div>
<div class="warn">两个容易误解的地方：① 能控只保证“从起点能到终点”，<b>不保证能走你指定的路径</b>——比如小车不可能“一边向右开、一边位置向左减小”。② 能控也<b>不代表能稳定停在任意状态</b>——比如速度不为零的状态就不可能“停”在那里（10.3.3 节会看到）。</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>指尖平衡：$G=\frac{1}{s^2-g/d}$ 有正极点；只用输出的比例控制无法稳定它，PID 难调。</li>
<li>能控性判据：$\mathrm{rank}[B\;AB\cdots A^{n-1}B]=n$。</li>
<li>能控 ≠ 能沿任意轨迹走，也 ≠ 能停在任意状态。</li>
</ol></div>`,
    sims: {
      nocontrol(el) {
        const ui = card(el, { title: "试着只用比例控制扶住杆", hint: "控制律 $u=-k\\,\\phi$（只看角度）。右边是根轨迹：增益多大，闭环极点都逃不出右半平面或虚轴。点“踢一脚”给杆一个扰动。" });
        const [a, b] = ui.row("", "");
        const sc = new Scene(a, { height: 260 });
        const ps = new Plot(b, { height: 260, x: [-6, 6], y: [-8, 8], equal: true, xlabel: "σ", title: "u = −kφ 的闭环极点" });
        const sk = ui.slider({ label: "比例增益 k", min: 0, max: 40, step: 0.5, value: 15 });
        let z = [0.05, 0], xi = 0, xid = 0;
        ui.buttons([["👆 踢一脚", () => { z[1] += 0.3; }], ["↺ 扶正重来", () => { z = [0.05, 0]; xi = 0; xid = 0; }, "ghost"]]);
        CE.loop(ui.root, dt => {
          const k = sk.get();
          for (let i = 0; i < 10; i++) { const h = dt / 10, u = -k * z[0]; if (Math.abs(z[0]) < 1.4) { z = [z[0] + h * z[1], z[1] + h * (GD * Math.sin(z[0]) + u)]; xid += -u * h; xi += xid * h; } }
          drawFinger(sc, z[0], xi);
          const rt = P.roots([1, 0, k - GD]);
          ps.render = () => { ps.clear(); for (let kk = 0; kk <= 70; kk += 0.25) P.roots([1, 0, kk - GD]).forEach(r => ps.dot(r[0], r[1], COL.c1 + "77", 1.6)); rt.forEach(r => ps.cross(r[0], r[1], COL.c4, 8, 3)); };
          ps.render();
          ui.out(k < GD ? `k &lt; g/d：极点 ${cfmt(rt[0])}, ${cfmt(rt[1])}<br><span class="bad">有正极点，杆倒下</span>` : `k &gt; g/d：极点 ±${f2(Math.sqrt(k - GD))}j<br><span class="bad">纯虚数：永远来回摆，停不下来</span>`);
        });
      },
      carts(el) {
        const ui = card(el, { title: "能控性实验室：只推左边那辆车", hint: "用按钮（或键盘 ← →）推左车。目标：让两辆车都停在绿色旗子处。太难了？试试“自动驾驶”（一个状态反馈控制器）。把弹簧断开，看看能控矩阵的秩怎样变化。" });
        const sc = new Scene(ui.plotHost(), { height: 150 });
        const pt = new Plot(ui.plotHost(), { height: 180, x: [0, 15], y: [-1, 3], xlabel: "t", ylabel: "位置" });
        const sk = ui.slider({ label: "弹簧劲度 k（书中为 100）", min: 0, max: 100, step: 1, value: 20, fmt: v => v.toFixed(0), onInput: () => { auto = false; info(); } });
        let z = [0, 0, -0.8, 0], t = 0, push = 0, auto = false, H = [], Kfb = null;
        const goal = [1.5, 0, 0.7, 0];
        const rowB = ui.buttons([["◀ 推", () => { }], ["推 ▶", () => { }], ["🤖 自动驾驶", () => { auto = !auto; Kfb = auto ? gain() : null; info(); }], ["↺", () => { z = [0, 0, -0.8, 0]; t = 0; H = []; }, "ghost"]]);
        const [bl, br] = rowB.querySelectorAll("button");
        [[bl, -1], [br, 1]].forEach(([btn, s]) => { btn.onpointerdown = () => push = s; btn.onpointerup = btn.onpointerleave = () => push = 0; });
        const kd = e => { if (!ui.root.isConnected) { window.removeEventListener("keydown", kd); return; } if (e.key === "ArrowLeft") push = -1; if (e.key === "ArrowRight") push = 1; };
        const ku = e => { if (e.key === "ArrowLeft" || e.key === "ArrowRight") push = 0; };
        window.addEventListener("keydown", kd); window.addEventListener("keyup", ku);
        const sys = () => { const k = sk.get(); return { A: [[0, 1, 0, 0], [-k, 0, k, 0], [0, 0, 0, 1], [k, 0, -k, 0]], B: [[0], [1], [0], [0]] }; };
        function gain() { const { A, B } = sys(); if (M.rank(M.ctrb(A, B)) < 4) return null; return M.lqr(A, B, M.eye(4).map((r, i) => r.map(v => v * (i % 2 ? 1 : 30))), [[0.3]]).K[0]; }
        function info() {
          const { A, B } = sys(), Co = M.ctrb(A, B), rk = M.rank(Co);
          ui.out(`$C_o=${mtex(Co, 0)}$<br>rank = <b class="big">${rk}</b> ${rk === 4 ? '<span class="good">能控</span>' : '<span class="bad">不能控！右车不受影响</span>'}${auto && !Kfb ? "<br>不能控，自动驾驶无能为力" : ""}`);
        }
        info();
        CE.loop(ui.root, dt => {
          const { A, B } = sys();
          for (let i = 0; i < 10; i++) {
            const h = dt / 10; let u = push * 6;
            if (auto && Kfb) u = -M.dot(Kfb, M.vsub(z, goal)) + (sk.get() ? 0 : 0);
            u = clamp(u, -40, 40);
            z = rk4(s => M.vadd(M.mulv(A, s), M.flat(B).map(v => v * u)), z, 0, h);
          }
          t += dt; H.push([t, z[0], z[2]]); if (t > 15) { H = H.map(h => [h[0] - dt, h[1], h[2]]); H.shift(); t -= dt; }
          const c = sc.ctx, W = sc.W; sc.clear();
          const x2p = x => W * 0.35 + x * W * 0.22, Y = 90;
          c.strokeStyle = COL.axis; c.beginPath(); c.moveTo(0, Y + 22); c.lineTo(W, Y + 22); c.stroke();
          [[goal[0], COL.c3], [goal[2], COL.c3]].forEach(([g, col]) => { const X = x2p(g); c.strokeStyle = col; c.beginPath(); c.moveTo(X, Y + 22); c.lineTo(X, Y - 40); c.stroke(); c.fillStyle = col; c.beginPath(); c.moveTo(X, Y - 40); c.lineTo(X + 14, Y - 35); c.lineTo(X, Y - 30); c.fill(); });
          const X1 = x2p(z[0]), X2 = x2p(z[2]);
          if (sk.get() > 0) { c.strokeStyle = COL.c5; c.lineWidth = 2; c.beginPath(); const a = X2 + 22, b = X1 - 22; c.moveTo(a, Y); for (let i = 1; i < 12; i++) c.lineTo(a + (b - a) * i / 12, Y + (i % 2 ? -7 : 7)); c.lineTo(b, Y); c.stroke(); }
          [[X1, COL.c1, "m₁"], [X2, COL.c2, "m₂"]].forEach(([X, col, n]) => { c.fillStyle = col; c.beginPath(); c.roundRect(X - 22, Y - 14, 44, 30, 6); c.fill(); c.fillStyle = "#fff"; c.font = "600 13px serif"; c.textAlign = "center"; c.fillText(n, X, Y + 6); });
          if (push || auto) { c.fillStyle = COL.c4; c.font = "20px sans-serif"; c.fillText(auto ? "🤖" : push > 0 ? "→" : "←", X1 + (push < 0 ? 34 : -34), Y + 6); }
          pt.render = () => { pt.clear(); pt.hline(goal[0], COL.c1 + "66", [5, 4]); pt.hline(goal[2], COL.c2 + "66", [5, 4]); pt.line(H.map(h => h[0]), H.map(h => h[1]), COL.c1, 2); pt.line(H.map(h => h[0]), H.map(h => h[2]), COL.c2, 2); pt.legend([{ label: "x₁", color: COL.c1 }, { label: "x₂", color: COL.c2 }]); };
          pt.render();
        });
      },
    },
  });

  /* ════════════════ 10.3 状态反馈：极点配置 / LQR / 轨迹追踪 ════════════════ */
  Course.add({
    id: "v1-place", vol: 1, ord: 180, no: "10.3", chName: "第10章 状态空间控制器与观测器", ref: "卷1 §10.3",
    title: "状态反馈：极点配置、LQR 初探与轨迹追踪",
    sub: "把所有状态都反馈回来：$u=-Kz$。闭环矩阵变成 $A-BK$，它的特征值想放哪儿就放哪儿。",
    goals: ["会用系数对比法设计 $K$ 把极点放到指定位置", "理解 LQR 的代价函数 $J=\\int z^TQz+u^TRu\\,dt$ 与权重的权衡", "会用“前馈 + 误差反馈”$u=Fz_d+K_ee$ 做轨迹追踪"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">10.3.1</span>极点配置</h2>
<p>只反馈角度 $u=-ky$ 时，闭环矩阵 $\begin{bmatrix}0&1\\\frac gd-k&0\end{bmatrix}$ 的特征值是 $\pm\sqrt{g/d-k}$——要么一正一负，要么纯虚数，永远稳定不了。把<b>角速度</b>也反馈回来：</p>
<div class="eq" data-no="(10.3.7)–(10.3.9)" data-title="全状态反馈">
$$u=-K\bm z=-\begin{bmatrix}k_1 & k_2\end{bmatrix}\bm z,\qquad A_{cl}=A-BK=\begin{bmatrix}0&1\\\frac gd-k_1&-k_2\end{bmatrix},\qquad \lambda^2+k_2\lambda-\left(\tfrac gd-k_1\right)=0$$
<div class="say">本质上还是比例控制，只是位置和速度都有自己的增益。$k_2$ 提供了“阻尼”，正是它让稳定成为可能。</div>
</div>
<p>想要闭环极点为 $-1,-1$，即特征多项式 $\lambda^2+2\lambda+1$。系数对比：$k_2=2$，$k_1=1+\frac gd$。取 $g/d=10$ 得 $K=[11,\;2]$，原来的鞍点变成了稳定节点。因为“状态矩阵特征值 = 传递函数极点”，这种方法叫<span class="kw">极点配置</span>。一般地，只要系统<b>能控</b>，极点就可以任意配置。</p>
<div class="sim" data-id="place"></div>

<h2><span class="sec">10.3.2</span>LQR 初探：极点放哪里才好？</h2>
<p>上面的 −1 是随手选的。极点越靠左，收敛越快，但手指的加速度（控制量）越大、越费力。怎样权衡？引入<span class="kw">代价函数</span>：</p>
<div class="eq" data-no="(10.3.15)" data-title="线性二次型调节器 LQR">
$$J=\int_0^{\infty}\left(\bm z^TQ\bm z+\bm u^TR\bm u\right)\mathrm{d}t,\qquad Q=\mathrm{diag}(q_1,\dots,q_n),\;R=\mathrm{diag}(r_1,\dots)$$
<dl><dt>$\bm z^TQ\bm z$</dt><dd>状态偏离平衡点的“惩罚”；$q_i$ 越大，越在意第 i 个状态尽快归零</dd><dt>$\bm u^TR\bm u$</dt><dd>控制量的“能耗惩罚”；$r$ 越大越省力</dd></dl>
<div class="say">选 K 使 J 最小。怎么求？卷 2 第 4 章会详细推导（Riccati 方程），这里先用软件算出结果、体会权重的作用。</div>
</div>
<table class="tbl"><tr><th>书中测试</th><th>Q</th><th>R</th><th>LQR 得到的 K</th><th>特点</th></tr>
<tr><td>组 1</td><td>diag(100, 1)</td><td>1</td><td>[24.14, 7.02]</td><td>最在意角度：回正最快，但控制量最大</td></tr>
<tr><td>组 2</td><td>diag(1, 100)</td><td>1</td><td>[20.05, 11.84]</td><td>最在意角速度：动作缓慢平稳，回正最慢</td></tr>
<tr><td>组 3</td><td>diag(1, 1)</td><td>100</td><td>[20.00, 6.33]</td><td>最在意能耗：控制量最小</td></tr></table>
<div class="sim" data-id="lqr"></div>
<div class="idea" data-t="权衡（trade-off）">控制里的“最优”不像百米赛跑那样只有一个标准：从性能看组 1 最好，从节能看组 3 最好。机器学习、神经网络本质上也是“建代价函数、设权重、求极值”。根据场景设定权重，是优化工程师最重要的工作之一。</div>

<h2><span class="sec">10.3.3</span>轨迹追踪：把平衡点搬到目标处</h2>
<p>前面的控制器都让状态回到 0。现在要让角度跟随目标 $\bm z_d=[\phi_d,\;0]^T$。定义误差 $\bm e=\bm z_d-\bm z$，目标是让 $\bm e=0$ 成为稳定平衡点。设 $\bm z_d$ 为常数，误差动态为：</p>
<div class="eq" data-no="(10.3.22)">
$$\dot e_1 = e_2 - z_{2d},\qquad \dot e_2=\frac gd e_1-\frac gd z_{1d}-u$$
<div class="say">第一个方程里根本没有 u！所以 $e_1$ 要平衡，必须 $e_2=z_{2d}$；而 $e_2$ 想为 0，只能 $z_{2d}=0$。也就是说：<b>角速度唯一能“停住”的值是 0</b>——能控，但不能稳定在任意状态。</div>
</div>
<div class="eq" data-no="(10.3.23)–(10.3.24)" data-title="前馈 + 误差反馈">
$$u=\underbrace{F\bm z_d}_{\text{前馈}}+\underbrace{K_e\bm e}_{\text{反馈}},\quad F=\begin{bmatrix}-\frac gd & 0\end{bmatrix}\;\Rightarrow\;\dot{\bm e}=\begin{bmatrix}0&1\\\frac gd-k_{e1}&-k_{e2}\end{bmatrix}\bm e$$
<dl><dt>前馈 $F\bm z_d$</dt><dd>抵消重力在目标角度处造成的“倾倒力矩”，把平衡点挪到 $\phi_d$</dd><dt>反馈 $K_e\bm e$</dt><dd>让新的平衡点稳定（极点配置或 LQR 设计 $K_e$）</dd></dl>
</div>
<p>书中取 $K_e=[25,\;7]$（闭环极点 $-3.5\pm1.66j$），并与 PID（$K_P=50,K_I=8,K_D=10$）对比：</p>
<div class="sim" data-id="track"></div>
<p>状态反馈利用了全部状态，可调参数更多、更灵活；PID 只用输出误差。但如果角速度测不到怎么办？下一课：观测器。</p>

<div class="summary"><h3>本课要点</h3><ol>
<li>全状态反馈 $u=-Kz$，闭环 $A-BK$；系统能控时可任意配置极点。</li>
<li>LQR：最小化 $J=\int z^TQz+u^TRu\,dt$；Q 管性能，R 管能耗，最优是一种权衡。</li>
<li>轨迹追踪：$u=Fz_d+K_ee$；前馈搬动平衡点，反馈保证稳定。速度只能稳定在 0。</li>
</ol></div>`,
    sims: {
      place(el) {
        const ui = card(el, { title: "极点配置沙盒", hint: "<b>拖动右图中的两个期望极点</b>（保持共轭或都在实轴上），程序自动算出 $K=[\\lambda_1\\lambda_2+g/d,\\;-(\\lambda_1+\\lambda_2)]$ 并用它扶杆。点“踢一脚”测试。" });
        const sc = new Scene(ui.plotHost(), { height: 230 });
        const [a, b] = ui.row("", "");
        const ph = new Plot(a, { height: 240, x: [-0.4, 0.4], y: [-1.5, 1.5], xlabel: "φ", ylabel: "φ̇", title: "相平面" });
        const ps = new Plot(b, { height: 240, x: [-8, 2], y: [-5, 5], equal: true, xlabel: "σ", title: "闭环极点（拖动）" });
        let pl = [[-1, 0], [-1, 0]], z = [0.15, 0], xi = 0, xid = 0, H = [], uNow = 0, mode = "re";
        ui.seg({ label: "极点类型", options: [["re", "两个实极点"], ["cx", "共轭复极点"]], value: "re", onChange: v => { mode = v; pl = v === "re" ? [[-1, 0], [-1, 0]] : [[-2, 2], [-2, -2]]; info(); } });
        ui.buttons([["👆 踢一脚", () => { z[1] += 0.6; H = []; }], ["↺ 从 φ=0.15 释放", () => { z = [0.15, 0]; xi = 0; xid = 0; H = []; }, "ghost"]]);
        const K = () => { const s = pl[0][0] + pl[1][0], p = mode === "re" ? pl[0][0] * pl[1][0] : pl[0][0] ** 2 + pl[0][1] ** 2; return [p + GD, -s]; };
        let drag = -1;
        ps.drag({
          down: (x, y) => { drag = -1; pl.forEach((q, i) => { if (ps.pxDist(q[0], q[1], x, y) < 16) drag = i; }); return drag >= 0 ? undefined : false; },
          move: (x, y) => { x = clamp(x, -7.8, 1.5); if (mode === "re") pl[drag] = [x, 0]; else { const w = Math.max(0.1, Math.abs(y)); pl = [[x, w], [x, -w]]; } info(); },
          hover: (x, y) => pl.some(q => ps.pxDist(q[0], q[1], x, y) < 16),
        });
        function info() { const k = K(); ui.out(`K = [<b>${f2(k[0])}</b>, <b>${f2(k[1])}</b>]<br>闭环极点：${cfmt(pl[0])}，${cfmt(pl[1])}<br>当前手指加速度指令 u = ${f2(uNow)}`); }
        info();
        CE.loop(ui.root, dt => {
          const k = K();
          for (let i = 0; i < 10; i++) { const h = dt / 10; uNow = -(k[0] * z[0] + k[1] * z[1]); if (Math.abs(z[0]) < 1.4) { z = [z[0] + h * z[1], z[1] + h * (GD * z[0] + uNow)]; xid += -uNow * h; xi += xid * h; } }
          xid *= 0.995; H.push(z.slice()); if (H.length > 600) H.shift();
          drawFinger(sc, z[0], xi, { u: uNow });
          ph.render = () => { ph.clear(); ph.line(H.map(p => p[0]), H.map(p => p[1]), COL.c1, 2); ph.dot(z[0], z[1], COL.c4, 5); };
          ps.render = () => { ps.clear(); ps.rect(0, -5, 2, 5, COL.c4 + "12"); ps.cross(-Math.sqrt(GD), 0, COL.muted, 6, 2); ps.cross(Math.sqrt(GD), 0, COL.muted, 6, 2); ps.text(Math.sqrt(GD), 0.6, "开环", COL.muted, "center"); pl.forEach(q => ps.cross(q[0], q[1], COL.c4, 9, 3)); };
          ph.render(); ps.render();
          if (Math.random() < 0.1) info();
        });
      },
      lqr(el) {
        const ui = card(el, { title: "LQR 权重调音台", hint: "初始状态 $z(0)=[\\pi/20,\\;0]$。调 $q_1,q_2,r$（对数滑块），实时求解 Riccati 方程得到 K。灰色曲线是上一次的结果。" });
        const [a, b, c3] = ui.row("", "", ""); a.parentElement.classList.add("r3");
        const p1 = new Plot(a, { height: 200, x: [0, 4], y: [-0.05, 0.17], xlabel: "t", title: "角度 z₁" });
        const p2 = new Plot(b, { height: 200, x: [0, 4], y: [-0.6, 0.1], xlabel: "t", title: "角速度 z₂" });
        const p3 = new Plot(c3, { height: 200, x: [0, 4], y: [-4, 0.5], xlabel: "t", title: "控制量 u" });
        const sq1 = ui.slider({ label: "$q_1$（角度权重）", min: 0.01, max: 1000, value: 100, log: true, onInput: draw });
        const sq2 = ui.slider({ label: "$q_2$（角速度权重）", min: 0.01, max: 1000, value: 1, log: true, onInput: draw });
        const sr = ui.slider({ label: "$r$（控制量权重）", min: 0.01, max: 1000, value: 1, log: true, onInput: draw });
        const set = (q1, q2, r) => { sq1.set(q1); sq2.set(q2); sr.set(r); draw(); };
        ui.buttons([["组1", () => set(100, 1, 1), "ghost"], ["组2", () => set(1, 100, 1), "ghost"], ["组3", () => set(1, 1, 100), "ghost"]]);
        let prev = null;
        function sim(K) { const dt = 0.005; let z = [Math.PI / 20, 0]; const ts = [], a = [], b = [], u = []; let J = [0, 0, 0]; for (let t = 0; t <= 4; t += dt) { const uu = -(K[0] * z[0] + K[1] * z[1]); ts.push(t); a.push(z[0]); b.push(z[1]); u.push(uu); J[0] += sq1.get() * z[0] ** 2 * dt; J[1] += sq2.get() * z[1] ** 2 * dt; J[2] += sr.get() * uu * uu * dt; z = rk4(s => [s[1], GD * s[0] + uu], z, t, dt); } return { ts, a, b, u, J }; }
        function draw() {
          const { K } = M.lqr(A0, B0, [[sq1.get(), 0], [0, sq2.get()]], [[sr.get()]]);
          const k = K[0], r = sim(k);
          [p1, p2, p3].forEach((p, i) => { const key = ["a", "b", "u"][i]; p.autoY([r[key], prev ? prev[key] : [], [0]], 0.08); p.render = () => { p.clear(); if (prev) p.line(prev.ts, prev[key], COL.muted, 1.4, [4, 3]); p.line(r.ts, r[key], [COL.c1, COL.c2, COL.c4][i], 2.4); }; p.render(); });
          const ev = P.roots([1, k[1], k[0] - GD]);
          ui.out(`K = [<b>${f2(k[0], 4)}</b>, <b>${f2(k[1], 4)}</b>]<br>闭环极点：${ev.map(e => cfmt(e)).join("，")}<br>J 的组成：∫q₁z₁² = ${f2(r.J[0], 3)}，∫q₂z₂² = ${f2(r.J[1], 3)}，∫ru² = ${f2(r.J[2], 3)}<br>u 的峰值：${f2(Math.min(...r.u), 2)}`);
          clearTimeout(draw._t); draw._t = setTimeout(() => { prev = r; }, 600);
        }
        draw();
      },
      track(el) {
        const ui = card(el, { title: "方波追踪擂台：状态反馈 + 前馈 vs PID", hint: "目标角度为书中的方波（±π/20）。可以关掉前馈，或修改两种控制器的参数。" });
        const pt = new Plot(ui.plotHost(), { height: 280, x: [0, 40], y: [-0.3, 0.3], xlabel: "t / s", ylabel: "φ / rad" });
        const pu = new Plot(ui.plotHost(), { height: 160, x: [0, 40], y: [-15, 15], xlabel: "t", ylabel: "u" });
        const s1 = ui.slider({ label: "$k_{e1}$", min: 11, max: 60, step: 0.5, value: 25, onInput: draw });
        const s2 = ui.slider({ label: "$k_{e2}$", min: 0.5, max: 20, step: 0.1, value: 7, onInput: draw });
        const cF = ui.check({ label: "启用前馈 F = [−g/d, 0]", value: true, onChange: draw });
        const sP = ui.slider({ label: "PID：$K_P$", min: 10, max: 120, step: 1, value: 50, fmt: v => v.toFixed(0), onInput: draw });
        const sI = ui.slider({ label: "PID：$K_I$", min: 0, max: 30, step: 0.5, value: 8, onInput: draw });
        const sD = ui.slider({ label: "PID：$K_D$", min: 1, max: 30, step: 0.5, value: 10, onInput: draw });
        const ref = t => t < 10 ? Math.PI / 20 : t < 20 ? 0 : t < 30 ? -Math.PI / 20 : 0;
        function draw() {
          const dt = 0.002, ts = [], a = [], b = [], ua = [], ub = [];
          let z = [0, 0], w = [0, 0], I = 0, ep = ref(0);
          for (let t = 0; t <= 40; t += dt) {
            const r = ref(t);
            const u1 = (cF.get() ? -GD * r : 0) + s1.get() * (r - z[0]) + s2.get() * (0 - z[1]);
            const e = r - w[0]; I += e * dt; const de = (e - ep) / dt; ep = e;
            const u2 = sP.get() * e + sI.get() * I + sD.get() * (Math.abs(de) > 50 ? 0 : de);
            ts.push(t); a.push(z[0]); b.push(w[0]); ua.push(u1); ub.push(u2);
            z = [z[0] + dt * z[1], z[1] + dt * (GD * z[0] + u1)]; w = [w[0] + dt * w[1], w[1] + dt * (GD * w[0] + u2)];
          }
          const mx = Math.max(0.25, ...a.map(Math.abs), ...b.map(Math.abs).filter(v => v < 3)); pt.setRange(null, [-mx * 1.1, mx * 1.1]);
          pt.render = () => { pt.clear(); pt.line(ts, ts.map(ref), COL.c3, 2, [6, 4]); pt.line(ts, b, COL.c2, 2); pt.line(ts, a, COL.c1, 2.4); pt.legend([{ label: "目标 φd", color: COL.c3, dash: [6, 4] }, { label: "状态反馈+前馈", color: COL.c1 }, { label: "PID", color: COL.c2 }], "tr"); };
          pu.render = () => { pu.clear(); pu.line(ts, ub.map(v => clamp(v, -15, 15)), COL.c2, 1.4); pu.line(ts, ua.map(v => clamp(v, -15, 15)), COL.c1, 1.6); };
          pt.render(); pu.render();
          const ev = P.roots([1, s2.get(), s1.get() - GD]), err = ts.reduce((s, t, i) => s + Math.abs(a[i] - ref(t)), 0) * dt, err2 = ts.reduce((s, t, i) => s + Math.abs(b[i] - ref(t)), 0) * dt;
          ui.out(`误差系统极点：${ev.map(e => cfmt(e)).join("，")}<br>累计绝对误差：状态反馈 <b>${f2(err, 3)}</b>，PID <b>${f2(err2, 3)}</b>${cF.get() ? "" : '<br><span class="bad">关掉前馈：平衡点不在目标处，有稳态误差</span>'}`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 10.4 观测器 ════════════════ */
  Course.add({
    id: "v1-obs", vol: 1, ord: 190, no: "10.4", chName: "第10章 状态空间控制器与观测器", ref: "卷1 §10.4–10.5",
    title: "观测器：看不见的状态，算出来",
    sub: "传感器装不到的地方怎么办？在后台同步运行一个“影子系统”，用测得的输出不断纠正它——这就是龙伯格观测器。",
    goals: ["理解能观测性与判据 $\\mathrm{rank}[C;CA;\\cdots]=n$", "推导观测器误差动态 $\\dot{\\tilde z}=(A-LC)\\tilde z$ 并设计 L", "理解分离原理，以及“观测器要比控制器快 2–5 倍”"],
    minutes: 30,
    html: String.raw`
<p>很多状态没法直接测：电动汽车电池内芯的温度只能通过表面温度推断；隧道里没有 GPS，车速要靠其他信息估计；指尖平衡里，角速度可能也没有传感器。这时就需要<span class="kw">观测器</span>。</p>
<h2><span class="sec">10.4.1</span>能观测性</h2>
<div class="eq" data-no="(10.4.2)" data-title="能观测性判据">
$$O=\begin{bmatrix}C\\CA\\\vdots\\CA^{n-1}\end{bmatrix},\qquad \text{能观测}\iff\mathrm{rank}(O)=n$$
<div class="say">能观测：根据一段时间内的输入和输出，能唯一确定初始状态。它与能控性是一对“对偶”的概念。</div>
</div>

<h2><span class="sec">10.4.2</span>龙伯格观测器</h2>
<p>最朴素的想法：把模型复制一份在电脑里跑，$\dot{\hat{\bm z}}=A\hat{\bm z}+B\bm u$。如果初值猜对了，估计就永远对；可一旦初值猜错，估计就会一直错下去（对不稳定系统甚至越差越远）。改进：用“实际输出 − 估计输出”来不断纠正：</p>
<div class="eq" data-no="(10.4.4)–(10.4.5)" data-title="龙伯格（Luenberger）观测器">
$$\dot{\hat{\bm z}}=A\hat{\bm z}+B\bm u+L\,(\bm y-\hat{\bm y})=(A-LC)\hat{\bm z}+(B-LD)\bm u+L\bm y$$
<dl><dt>$\hat{\bm z}$</dt><dd>状态估计值</dd><dt>$\bm y-\hat{\bm y}$</dt><dd>测量值与预测值之差：模型“猜错了多少”</dd><dt>$L$</dt><dd>观测器增益：纠正的力度</dd></dl>
</div>
<p>用真实系统减去观测器，得到估计误差 $\tilde{\bm z}=\bm z-\hat{\bm z}$ 的动态：</p>
<div class="eq" data-no="(10.4.8)">
$$\dot{\tilde{\bm z}}=(A-LC)\,\tilde{\bm z}$$
<div class="say">只要选 L 让 $A-LC$ 的特征值都在左半平面，估计误差就会衰减到 0——不管初值猜得多离谱。这和极点配置一模一样，只是把 $A-BK$ 换成了 $A-LC$。</div>
</div>
<p><b>例</b>：弹簧–质量–阻尼（$m=1,k=1,b=0.5$），只能测位移。令 $A-LC$ 的特征值为 $-1,-1$：特征式 $\lambda^2+(0.5+l_1)\lambda+0.5l_1+1+l_2$ 对比 $\lambda^2+2\lambda+1$，得 $l_1=1.5,\;l_2=-0.75$。真实初值 $[1,1]$，观测器从 $[0,0]$ 开始猜：</p>
<div class="sim" data-id="observer"></div>

<h2><span class="sec">10.4.4</span>观测器 + 控制器：分离原理</h2>
<p>回到指尖平衡，假设角速度测不到。用估计值做反馈 $u=-K\hat{\bm z}$，代入得：</p>
<div class="eq" data-no="(10.4.20)" data-title="分离原理">
$$\frac{\mathrm{d}}{\mathrm{d}t}\begin{bmatrix}\tilde{\bm z}\\\bm z\end{bmatrix}=\begin{bmatrix}A-LC & 0\\ BK & A-BK\end{bmatrix}\begin{bmatrix}\tilde{\bm z}\\\bm z\end{bmatrix}$$
<div class="say">块三角矩阵的特征值 = $A-LC$ 的特征值 ∪ $A-BK$ 的特征值。所以 <b>K 和 L 可以分开独立设计</b>。但要注意：观测器的收敛速度一般要比控制器快 2–5 倍，控制器才是在“相对准确”的估计上工作。</div>
</div>
<div class="sim" data-id="separation"></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>能观测判据：$\mathrm{rank}[C;CA;\cdots;CA^{n-1}]=n$。</li>
<li>观测器 $\dot{\hat z}=A\hat z+Bu+L(y-\hat y)$；误差 $\dot{\tilde z}=(A-LC)\tilde z$，设计 L 使其稳定。</li>
<li>分离原理：K 与 L 可独立设计；观测器要比控制器快 2–5 倍。</li>
</ol></div>`,
    sims: {
      observer(el) {
        const ui = card(el, { title: "观测器“追影子”", hint: "实线：真实状态；虚线：观测器的估计。调观测器极点，或者把 L 设为 0（纯靠模型猜），比较估计误差。" });
        const [a, b] = ui.row("", "");
        const p1 = new Plot(a, { height: 220, x: [0, 15], y: [-1.5, 1.5], xlabel: "t", title: "位移 z₁（可测）" });
        const p2 = new Plot(b, { height: 220, x: [0, 15], y: [-1.5, 1.5], xlabel: "t", title: "速度 z₂（不可测，只能估计）" });
        const pe = new Plot(ui.plotHost(), { height: 170, x: [0, 15], y: [-4, 0.5], xlabel: "t", title: "估计误差 log₁₀‖z − ẑ‖" });
        const A = [[0, 1], [-1, -0.5]], C = [[1, 0]];
        const sp = ui.slider({ label: "观测器极点（重根）λ", min: -6, max: -0.1, step: 0.05, value: -1, onInput: draw });
        const sn = ui.slider({ label: "测量噪声", min: 0, max: 0.1, step: 0.005, value: 0, fmt: v => v.toFixed(3), onInput: draw });
        const cZ = ui.check({ label: "L = 0（开环“猜”）", value: false, onChange: draw });
        function draw() {
          const lam = sp.get(); // λ²−2λ·lam... (λ−lam)² = λ² − 2lam λ + lam²
          let l1 = -2 * lam - 0.5, l2 = lam * lam - 1 - 0.5 * l1; if (cZ.get()) { l1 = 0; l2 = 0; }
          const dt = 0.005, rg = CE.rng(5); let z = [1, 1], zh = [0, 0]; const ts = [], Z = [], Zh = [], E = [];
          for (let t = 0; t <= 15; t += dt) {
            const y = z[0] + sn.get() * rg.n(), innov = y - zh[0];
            ts.push(t); Z.push(z.slice()); Zh.push(zh.slice()); E.push(Math.log10(Math.max(1e-6, Math.hypot(z[0] - zh[0], z[1] - zh[1]))));
            z = rk4(s => M.mulv(A, s), z, t, dt);
            zh = [zh[0] + dt * (zh[1] + l1 * innov), zh[1] + dt * (-zh[0] - 0.5 * zh[1] + l2 * innov)];
          }
          p1.render = () => { p1.clear(); p1.line(ts, Z.map(v => v[0]), COL.c1, 2.4); p1.line(ts, Zh.map(v => v[0]), COL.c4, 2.2, [6, 4]); p1.legend([{ label: "真实", color: COL.c1 }, { label: "估计", color: COL.c4, dash: [6, 4] }], "tr"); };
          p2.render = () => { p2.clear(); p2.line(ts, Z.map(v => v[1]), COL.c1, 2.4); p2.line(ts, Zh.map(v => v[1]), COL.c4, 2.2, [6, 4]); };
          pe.render = () => { pe.clear(); pe.line(ts, E, COL.c5, 2); };
          p1.render(); p2.render(); pe.render();
          ui.out(`L = [<b>${f2(l1)}</b>, <b>${f2(l2)}</b>]ᵀ<br>$A-LC=${mtex([[-l1, 1], [-1 - l2, -0.5]])}$<br>${cZ.get() ? "开环观测器：误差只能靠系统自身的阻尼慢慢衰减（极点 −0.25±0.97j）" : "书中取 λ=−1：l₁=1.5，l₂=−0.75"}`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b>：观测器极点越靠左，收敛越快——但加一点测量噪声后，快的观测器估计变得很“抖”。快速与抗噪之间要权衡，这正是卷 2 卡尔曼滤波要解决的问题。</div>`);
      },
      separation(el) {
        const ui = card(el, { title: "只测角度也能扶杆：观测器 + 控制器", hint: "控制器极点固定为 −2 ± 1j。调“观测器速度倍数”：观测器极点 = 倍数 × 控制器极点。观测器初值故意猜错（以为杆是竖直静止的）。" });
        const sc = new Scene(ui.plotHost(), { height: 210 });
        const [a, b] = ui.row("", "");
        const pt = new Plot(a, { height: 220, x: [0, 8], y: [-0.4, 0.4], xlabel: "t", title: "角度：真实 vs 估计" });
        const ps = new Plot(b, { height: 220, x: [-16, 1], y: [-6, 6], equal: true, xlabel: "σ", title: "增广系统的 4 个特征值" });
        const sm = ui.slider({ label: "观测器速度倍数", min: 0.3, max: 6, step: 0.1, value: 3, onInput: restart });
        ui.buttons([["↺ 重来", restart, "ghost"]]);
        const pc = [[-2, 1], [-2, -1]];
        let z, zh, xi, xid, H, t, K, L;
        function restart() {
          const m = sm.get();
          K = M.acker(A0, B0, pc)[0];
          L = M.flat(M.T(M.acker(M.T(A0), M.T(C0), pc.map(p => [p[0] * m, p[1] * m]))));
          z = [0.2, 0.3]; zh = [0, 0]; xi = 0; xid = 0; H = []; t = 0;
          const Aug = M.vcat(M.hcat(M.sub(A0, M.mul(M.col(L), C0)), M.zeros(2, 2)), M.hcat(M.mul(B0, [K]), M.sub(A0, M.mul(B0, [K]))));
          const ev = M.eig(Aug);
          ps.render = () => { ps.clear(); ps.rect(0, -6, 1, 6, COL.c4 + "12"); ev.forEach(e => ps.cross(e[0], e[1], COL.c4, 7, 2.6)); pc.forEach(p => ps.circle(p[0], p[1], COL.c1, 11, 1.6)); ps.text(-15.5, 5.2, "蓝圈：控制器极点", COL.c1); };
          ps.render();
          ui.out(`K = [${f2(K[0])}, ${f2(K[1])}]<br>L = [${f2(L[0])}, ${f2(L[1])}]ᵀ<br>${m < 1 ? '<span class="bad">观测器比控制器还慢：初期基于错误估计在控制</span>' : m <= 5 ? '<span class="good">观测器快 2–5 倍：接近全状态反馈</span>' : "非常快的观测器：对噪声会很敏感"}`);
        }
        restart();
        CE.loop(ui.root, dt => {
          if (t > 8) return;
          for (let i = 0; i < 10; i++) {
            const h = dt / 10, u = -(K[0] * zh[0] + K[1] * zh[1]), y = z[0], inn = y - zh[0];
            if (Math.abs(z[0]) < 1.4) { z = [z[0] + h * z[1], z[1] + h * (GD * z[0] + u)]; xid += -u * h; xi += xid * h; }
            zh = [zh[0] + h * (zh[1] + L[0] * inn), zh[1] + h * (GD * zh[0] + u + L[1] * inn)];
            t += h;
          }
          xid *= 0.995; H.push([t, z[0], zh[0]]);
          drawFinger(sc, z[0], xi);
          pt.render = () => { pt.clear(); pt.line(H.map(h => h[0]), H.map(h => h[1]), COL.c1, 2.4); pt.line(H.map(h => h[0]), H.map(h => h[2]), COL.c4, 2, [6, 4]); pt.legend([{ label: "真实 φ", color: COL.c1 }, { label: "估计 φ̂", color: COL.c4, dash: [6, 4] }], "tr"); };
          pt.render(); ps.render();
        });
      },
    },
  });

  /* ════════════════ 附录 A：线性化 ════════════════ */
  Course.add({
    id: "v1-appA", vol: 1, ord: 200, no: "附录A", chName: "附录", ref: "卷1 附录A",
    title: "非线性系统的线性化",
    sub: "在平衡点附近用一条切线代替曲线：泰勒展开、雅可比矩阵，以及线性化在什么时候会失效。",
    goals: ["会用一阶泰勒展开线性化一个函数", "会求状态方程在平衡点处的雅可比矩阵", "理解线性化模型离平衡点越远误差越大"],
    minutes: 15,
    html: String.raw`
<p>本书几乎所有方法都针对线性系统，而真实世界大多是非线性的。桥梁是<span class="kw">线性化</span>：在工作点附近，用切线近似曲线。</p>
<div class="eq" data-no="(A.1)–(A.2)" data-title="一阶泰勒展开">
$$f(x)\approx f(x_0)+f'(x_0)\,(x-x_0)$$
<div class="say">忽略二次及更高次项。所以 x 离 $x_0$ 越远，误差越大。</div>
</div>
<p>例：$\sin x$ 在 0 附近 ≈ $x$。$x=\pi/6$ 时误差约 4%，$x=\pi/4$ 时约 11%。指尖平衡里的 $\sin\phi\approx\phi$、$\cos\phi\approx1$ 正是这样来的。</p>
<div class="sim" data-id="sinlin"></div>
<p>对状态方程 $\dot{\bm z}=\bm f(\bm z)$，在平衡点 $\bm z_f$（$\bm f(\bm z_f)=0$）附近令 $\bm z=\bm z_f+\bm z_\delta$：</p>
<div class="eq" data-no="(A.9)" data-title="雅可比矩阵">
$$\dot{\bm z}_\delta=\left.\begin{bmatrix}\frac{\partial f_1}{\partial z_1}&\cdots&\frac{\partial f_1}{\partial z_n}\\\vdots&\ddots&\vdots\\\frac{\partial f_n}{\partial z_1}&\cdots&\frac{\partial f_n}{\partial z_n}\end{bmatrix}\right|_{\bm z=\bm z_f}\bm z_\delta$$
</div>
<p>书中例子：$\ddot x+\dot x+\frac1x=1$，平衡点 $x_f=1$。令 $z_1=x,z_2=\dot x$：$f_1=z_2$，$f_2=1-\frac{1}{z_1}-z_2$，雅可比 $\begin{bmatrix}0&1\\\frac{1}{z_1^2}&-1\end{bmatrix}_{z_f}=\begin{bmatrix}0&1\\1&-1\end{bmatrix}$，即 $\ddot x_\delta+\dot x_\delta-x_\delta=0$（特征值一正一负：该平衡点是鞍点）。</p>
<div class="sim" data-id="pend"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>线性化 = 在平衡点处做一阶泰勒展开，多维情形用雅可比矩阵。</li>
<li>线性模型只在平衡点附近可信，偏离越大误差越大。</li>
</ol></div>`,
    sims: {
      sinlin(el) {
        const ui = card(el, { title: "sin x ≈ x 有多准？", hint: "拖动 x，比较真实值与线性化值。" });
        const pl = new Plot(ui.plotHost(), { height: 240, x: [-Math.PI, Math.PI], y: [-2, 2], xlabel: "x" });
        const sx = ui.slider({ label: "x (rad)", min: 0, max: 1.6, step: 0.01, value: Math.PI / 6, onInput: draw });
        function draw() { const x = sx.get(); pl.render = () => { pl.clear(); pl.fn(Math.sin, COL.c1, 2.6); pl.fn(t => t, COL.c2, 2, [6, 4]); pl.vline(x, COL.muted, [3, 3]); pl.dot(x, Math.sin(x), COL.c1, 5); pl.dot(x, x, COL.c2, 5); pl.legend([{ label: "sin x", color: COL.c1 }, { label: "线性化 x", color: COL.c2, dash: [6, 4] }]); }; pl.render(); ui.out(`x = ${f2(x)}（${f2(x * 180 / Math.PI, 1)}°）<br>sin x = ${f2(Math.sin(x), 4)}<br>相对误差 <b class="big">${f2(Math.abs(x - Math.sin(x)) / Math.max(1e-9, Math.sin(x)) * 100, 1)}%</b>`); }
        draw();
      },
      pend(el) {
        const ui = card(el, { title: "单摆：非线性模型 vs 线性化模型", hint: "单摆 $\\ddot\\theta=-\\frac gL\\sin\\theta$（g/L=10）及其在下方平衡点处的线性化 $\\ddot\\theta=-\\frac gL\\theta$。初始角度越大，两者差别越明显。" });
        const pl = new Plot(ui.plotHost(), { height: 250, x: [0, 10], y: [-3.2, 3.2], xlabel: "t", ylabel: "θ" });
        const s = ui.slider({ label: "初始角度 θ₀ (°)", min: 5, max: 175, step: 1, value: 30, fmt: v => v.toFixed(0) + "°", onInput: draw });
        function draw() {
          const th0 = s.get() * Math.PI / 180, dt = 0.005; let a = [th0, 0], b = [th0, 0]; const ts = [], A = [], B = [];
          for (let t = 0; t <= 10; t += dt) { ts.push(t); A.push(a[0]); B.push(b[0]); a = rk4(x => [x[1], -10 * Math.sin(x[0])], a, t, dt); b = rk4(x => [x[1], -10 * x[0]], b, t, dt); }
          pl.render = () => { pl.clear(); pl.line(ts, A, COL.c1, 2.4); pl.line(ts, B, COL.c2, 2, [6, 4]); pl.legend([{ label: "非线性", color: COL.c1 }, { label: "线性化", color: COL.c2, dash: [6, 4] }], "tr"); };
          pl.render(); ui.out(`小角度时两条曲线几乎重合；大角度时非线性单摆的周期明显变长。`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 附录 B：傅里叶 ════════════════ */
  Course.add({
    id: "v1-appB", vol: 1, ord: 210, no: "附录B", chName: "附录", ref: "卷1 附录B",
    title: "傅里叶级数与变换：一切信号都是正弦的合唱",
    sub: "利用三角函数的正交性，把周期信号拆成一系列正弦和余弦；周期趋于无穷，就得到傅里叶变换——它是拉普拉斯变换在 s = jω 上的特例。",
    goals: ["理解三角函数系的正交性", "知道傅里叶系数怎么算，以及 $a_0/2$ 的来历", "理解频谱，以及傅里叶变换与拉普拉斯变换的关系"],
    minutes: 20,
    html: String.raw`
<h2><span class="sec">B.1</span>正交性：函数也能“垂直”</h2>
<p>两个向量垂直 ⇔ 内积为 0：$\bm a\cdot\bm b=\sum a_ib_i=0$。把求和推广成积分，就有了函数的“内积”：$\int_a^b f(x)g(x)\,\mathrm{d}x$。三角函数系 $\{1,\sin x,\cos x,\sin2x,\cos2x,\dots\}$ 中任意两个不同的函数在 $[-\pi,\pi]$ 上都正交，而 $\int_{-\pi}^{\pi}\cos^2nx\,\mathrm{d}x=\pi$。</p>
<h2><span class="sec">B.2</span>傅里叶级数</h2>
<div class="eq" data-no="(B.2.12)" data-title="周期为 2π 的傅里叶级数">
$$f(x)=\frac{a_0}{2}+\sum_{n=1}^{\infty}\left(a_n\cos nx+b_n\sin nx\right),\qquad a_n=\frac1\pi\int_{-\pi}^{\pi}f(x)\cos nx\,\mathrm{d}x,\quad b_n=\frac1\pi\int_{-\pi}^{\pi}f(x)\sin nx\,\mathrm{d}x$$
<div class="say">求某个系数，就把 f 乘上对应的三角函数再积分——正交性让其他项全部消失。</div>
</div>
<details class="derive"><summary>那个“不和谐”的 $a_0/2$ 是怎么来的？</summary><div class="body">
<div class="step">先写成 $f=a_0^*+\sum(a_n^*\cos nx+b_n\sin nx)$，两边积分得 $a_0^*=\frac{1}{2\pi}\int f\,\mathrm{d}x$。</div>
<div class="step">乘 $\cos nx$ 再积分得 $a_n^*=\frac1\pi\int f\cos nx\,\mathrm{d}x$——把 n=0 代进这个公式，会比 $a_0^*$ 多出一个因子 2。</div>
<div class="step">为了用同一个公式算所有 $a_n$，定义 $a_0=2a_0^*$，于是常数项写成 $a_0/2$。</div>
</div></details>
<p>周期为 $T$ 时，基频 $\omega_T=2\pi/T$，用 $\cos n\omega_Tt,\sin n\omega_Tt$ 展开，系数前面变成 $\frac2T$。用欧拉公式还可写成更简洁的复数形式：</p>
<div class="eq" data-no="(B.4.6)–(B.4.9)">
$$f_T(t)=\sum_{n=-\infty}^{\infty}c_ne^{jn\omega_Tt},\qquad c_n=\frac1T\int_{-T/2}^{T/2}f_T(t)e^{-jn\omega_Tt}\,\mathrm{d}t$$
<div class="say">$|c_n|$ 随频率 $n\omega_T$ 的分布就是信号的<b>频谱</b>：每个频率上有多“强”。</div>
</div>
<div class="sim" data-id="fourier"></div>
<h2><span class="sec">B.5</span>傅里叶变换</h2>
<p>让周期 $T\to\infty$：谱线间隔 $\Delta\omega=2\pi/T\to0$，离散的谱线连成连续曲线，求和变成积分：</p>
<div class="eq" data-no="(B.5.7)" data-title="傅里叶变换">
$$F(\omega)=\int_{-\infty}^{\infty}f(t)e^{-j\omega t}\,\mathrm{d}t,\qquad f(t)=\frac{1}{2\pi}\int_{-\infty}^{\infty}F(\omega)e^{j\omega t}\,\mathrm{d}\omega$$
<div class="say">对比拉普拉斯变换 $\int_0^\infty f(t)e^{-st}\mathrm{d}t$：傅里叶变换就是 $s=j\omega$ 的特例。实信号满足 $\overline{F(\omega)}=F(-\omega)$（共轭对称）——第 9 章推导频率响应时用到的正是这一性质。</div>
</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>三角函数系两两正交，所以系数可以“一个一个地”积分求出。</li>
<li>周期信号 = 直流 + 基频 + 各次谐波；$|c_n|$ 构成频谱。</li>
<li>傅里叶变换是周期 → ∞ 的极限，也是拉普拉斯变换在虚轴上的取值。</li>
</ol></div>`,
    sims: {
      fourier(el) {
        const ui = card(el, { title: "傅里叶叠加器", hint: "选择一个目标波形，增加谐波数 N，看部分和怎样逼近它；下方是对应的频谱 |cₙ|。" });
        const pl = new Plot(ui.plotHost(), { height: 250, x: [-2 * Math.PI, 2 * Math.PI], y: [-1.6, 1.6], xlabel: "x" });
        const pspec = new Plot(ui.plotHost(), { height: 160, x: [-0.5, 20.5], y: [0, 1.4], xlabel: "n（谐波次数）", title: "频谱：每个频率成分的振幅" });
        let kind = "square";
        const W = {
          square: { n: "方波", f: x => Math.sign(Math.sin(x)) || 1, c: n => n % 2 ? [0, 4 / (Math.PI * n)] : [0, 0] },
          saw: { n: "锯齿波", f: x => { let y = ((x + Math.PI) % (2 * Math.PI) + 2 * Math.PI) % (2 * Math.PI) - Math.PI; return y / Math.PI; }, c: n => [0, 2 * (n % 2 ? 1 : -1) / (Math.PI * n)] },
          tri: { n: "三角波", f: x => { const y = ((x % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI); return y < Math.PI ? 1 - 2 * y / Math.PI : -1 + 2 * (y - Math.PI) / Math.PI; }, c: n => n % 2 ? [8 / (Math.PI * Math.PI * n * n), 0] : [0, 0] },
          book: { n: "书中图 B.2.1", f: x => Math.sin(x) + 0.5 * Math.cos(2 * x) + 0.2 * Math.sin(3 * x) + 0.3 * Math.cos(4 * x), c: n => n === 1 ? [0, 1] : n === 2 ? [0.5, 0] : n === 3 ? [0, 0.2] : n === 4 ? [0.3, 0] : [0, 0] },
        };
        ui.seg({ label: "目标波形", options: Object.entries(W).map(([k, v]) => [k, v.n]), value: kind, onChange: v => { kind = v; draw(); } });
        const sN = ui.slider({ label: "谐波数 N", min: 1, max: 40, step: 1, value: 3, fmt: v => v.toFixed(0), onInput: draw });
        function draw() {
          const w = W[kind], N = sN.get(), xs = linspace(-2 * Math.PI, 2 * Math.PI, 800);
          const part = x => { let s = 0; for (let n = 1; n <= N; n++) { const [a, b] = w.c(n); s += a * Math.cos(n * x) + b * Math.sin(n * x); } return s; };
          pl.render = () => { pl.clear(); pl.line(xs, xs.map(w.f), COL.muted, 1.6, [5, 4]); for (let n = 1; n <= Math.min(N, 6); n++) { const [a, b] = w.c(n); if (a || b) pl.line(xs, xs.map(x => a * Math.cos(n * x) + b * Math.sin(n * x)), `hsl(${n * 50},50%,65%)`, 1); } pl.line(xs, xs.map(part), COL.c1, 2.6); };
          pspec.render = () => { pspec.clear(); for (let n = 1; n <= 20; n++) { const [a, b] = w.c(n), m = Math.hypot(a, b); if (m > 1e-9) { pspec.line([n, n], [0, m], n <= N ? COL.c1 : COL.axis, 5); } } };
          pl.render(); pspec.render();
          ui.out(`使用前 <b>${N}</b> 次谐波<br>${kind === "square" ? "方波只含奇次谐波，振幅 4/(πn)。注意跳变处的“过冲”始终约 9%（吉布斯现象）。" : kind === "book" ? "这就是书中由 4 个三角函数叠加出的波形。" : ""}`);
        }
        draw();
      },
    },
  });
})();
