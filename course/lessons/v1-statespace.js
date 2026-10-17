/* 卷1 第3章：状态空间方程 / 特征值与解耦 / 相平面 / 爱情故事 */
(function () {
  const { M, P, Plot, Scene, card, linspace, f2, COL, clamp, rk4, mtex, cfmt } = CE;

  /* ── 通用：二维相图绘制器 ── */
  function phasePlot(host, opt = {}) {
    const pl = new Plot(host, Object.assign({ height: 360, x: [-3, 3], y: [-3, 3], equal: true, xlabel: "z₁", ylabel: "z₂" }, opt));
    const st = { A: [[0, 1], [-1, 0]], f: null, trajs: [], parts: [], showField: true, showEig: true, showParts: true };
    const field = z => st.f ? st.f(z) : M.mulv(st.A, z);
    function traj(z0, T = 12, dt = 0.02) {
      const fw = [z0.slice()]; let z = z0.slice();
      for (let t = 0; t < T; t += dt) { z = rk4(field, z, 0, dt); if (Math.abs(z[0]) > 50 || Math.abs(z[1]) > 50) break; fw.push(z); }
      return fw;
    }
    function seed() { const [a, b] = pl.opt.x, [c, d] = pl.opt.y; return [a + Math.random() * (b - a), c + Math.random() * (d - c)]; }
    st.reset = (seeds) => {
      st.trajs = (seeds || defaultSeeds()).map(s => traj(s));
      st.parts = Array.from({ length: 70 }, () => ({ z: seed(), age: Math.random() * 120 }));
    };
    function defaultSeeds() { const out = [], R = (pl.opt.y[1] - pl.opt.y[0]) * 0.45; for (let k = 0; k < 12; k++) { const a = k / 12 * 2 * Math.PI + 0.13; out.push([R * Math.cos(a), R * Math.sin(a)]); } out.push([0.3, 0.05], [-0.3, -0.05]); return out; }
    st.add = z0 => { st.trajs.push(traj(z0)); if (st.trajs.length > 30) st.trajs.shift(); };
    st.tick = dt => {
      st.parts.forEach(p => {
        p.z = rk4(field, p.z, 0, dt * 0.5); p.age++;
        const [a, b] = pl.opt.x, [c, d] = pl.opt.y, out = p.z[0] < a || p.z[0] > b || p.z[1] < c || p.z[1] > d;
        if (out || p.age > 160 || !isFinite(p.z[0]) || (Math.hypot(...p.z) < 0.03)) { p.z = seed(); p.age = 0; }
      });
    };
    pl.render = () => {
      pl.clear();
      const [a, b] = pl.opt.x, [c, d] = pl.opt.y, n = 15;
      if (st.showField) {
        for (let i = 0; i <= n; i++) for (let j = 0; j <= n; j++) {
          const z = [a + (b - a) * i / n, c + (d - c) * j / n], v = field(z), L = Math.hypot(v[0], v[1]);
          if (L < 1e-9) continue; const s = (b - a) / n * 0.38;
          pl.arrow(z[0], z[1], z[0] + v[0] / L * s, z[1] + v[1] / L * s, COL.axis, 1.1, 5);
        }
      }
      if (st.showEig && !st.f) {
        const e = M.eig2(st.A);
        if (e.real) e.vecs.forEach((v, k) => { const L = 20; pl.line([-v[0] * L, v[0] * L], [-v[1] * L, v[1] * L], k ? COL.c5 : COL.c2, 1.6, [7, 5]); });
      }
      st.trajs.forEach(tr => pl.line(tr.map(p => p[0]), tr.map(p => p[1]), COL.c1 + "cc", 1.6));
      st.trajs.forEach(tr => pl.dot(tr[0][0], tr[0][1], COL.c1, 3));
      if (st.showParts) st.parts.forEach(p => pl.dot(p.z[0], p.z[1], COL.c4 + "bb", 2.2));
      pl.dot(0, 0, COL.ink, 4);
    };
    pl.drag({ down: (x, y) => { st.add([x, y]); pl.render(); return false; } });
    pl.cv.style.cursor = "crosshair";
    st.pl = pl;
    return st;
  }
  function classify(A) {
    const e = M.eig2(A), tr = e.tr, det = e.det;
    if (Math.abs(det) < 1e-9) return { name: "退化（有零特征值）", cls: "", e };
    if (!e.real) {
      if (Math.abs(tr) < 1e-6) return { name: "中心点（等幅循环）", cls: "", e };
      return tr < 0 ? { name: "稳定焦点（螺旋收敛）", cls: "good", e } : { name: "不稳定焦点（螺旋发散）", cls: "bad", e };
    }
    if (det < 0) return { name: "鞍点（不稳定）", cls: "bad", e };
    return tr < 0 ? { name: "稳定节点", cls: "good", e } : { name: "不稳定节点", cls: "bad", e };
  }

  /* ════════════════ 3.1 状态空间方程 ════════════════ */
  Course.add({
    id: "v1-ss", vol: 1, ord: 50, no: "3.1", chName: "第3章 动态系统建模——状态空间方程", ref: "卷1 §3.1",
    title: "状态空间方程：把高阶方程拆成一阶方程组",
    sub: "现代控制理论的语言。它把“位置、速度……”这些状态写成一个向量，用矩阵描述它们怎样互相影响。",
    goals: ["会把弹簧–质量–阻尼系统写成 $\\dot z=Az+Bu,\\;y=Cz+Du$", "知道 A、B、C、D 四个矩阵各管什么", "理解“A 的特征值 = 传递函数的极点”"],
    minutes: 20,
    html: String.raw`
<p>传递函数只描述“输入到输出”这一条通路，系统内部发生了什么它不管；而且多输入多输出（MIMO）时很难用。<span class="kw">状态空间方程</span>换了一个思路：选出一组能完整描述系统“此刻处境”的量——<b>状态变量</b>——然后写出它们的变化率。MPC、卡尔曼滤波、LQR……卷 2 的所有主角都建立在这种描述之上。</p>

<h2><span class="sec">3.1.1</span>以弹簧–质量–阻尼系统为例</h2>
<p>方程 $m\ddot x + b\dot x + kx = u$ 是二阶的。诀窍是选两个状态：<b>位置</b> $z_1 = x$ 和<b>速度</b> $z_2 = \dot x$。于是</p>
<div class="eq" data-no="(3.1.3)–(3.1.5)">
$$\dot z_1 = z_2,\qquad \dot z_2 = \ddot x = \frac{1}{m}\left(u - b z_2 - k z_1\right)$$
<div class="say">一个二阶方程，变成了两个一阶方程：位置的变化率是速度；速度的变化率由受力决定。</div>
</div>
<p>写成矩阵形式：</p>
<div class="eq" data-no="(3.1.6)–(3.1.7)">
$$\frac{\mathrm{d}}{\mathrm{d}t}\begin{bmatrix} z_1\\ z_2\end{bmatrix} = \underbrace{\begin{bmatrix}0 & 1\\ -\frac{k}{m} & -\frac{b}{m}\end{bmatrix}}_{A}\begin{bmatrix} z_1\\ z_2\end{bmatrix} + \underbrace{\begin{bmatrix}0\\ \frac1m\end{bmatrix}}_{B} u,\qquad y = \underbrace{\begin{bmatrix}1 & 0\end{bmatrix}}_{C}\begin{bmatrix} z_1\\ z_2\end{bmatrix} + \underbrace{0}_{D}\cdot u$$
</div>
<p>一般地，任何 LTI 系统都可以写成：</p>
<div class="eq" data-no="(3.1.8)" data-title="状态空间方程">
$$\dot{\bm z}(t) = A\,\bm z(t) + B\,\bm u(t),\qquad \bm y(t) = C\,\bm z(t) + D\,\bm u(t)$$
<dl>
<dt>$\bm z$ ($n\times1$)</dt><dd>状态向量：描述系统此刻处境的最少一组量</dd>
<dt>$\bm u$ ($p\times1$)</dt><dd>输入向量；$\bm y$ ($m\times1$) 输出向量</dd>
<dt>$A$ ($n\times n$)</dt><dd><b>状态矩阵</b>：状态之间怎样互相影响（系统的“性格”全在这里）</dd>
<dt>$B$ ($n\times p$)</dt><dd><b>输入矩阵</b>：输入推动了哪些状态</dd>
<dt>$C$ ($m\times n$)</dt><dd><b>输出矩阵</b>：我们能看到（测到）状态的哪些组合</dd>
<dt>$D$ ($m\times p$)</dt><dd><b>直接传递矩阵</b>：输入不经过状态直接到输出的部分，常为 0</dd>
</dl></div>
<p>书中还给了一个两回路共享电阻 $R_3$ 的电路：两个电压源输入、两个输出，用传递函数很别扭，用状态空间只是把 $A,B,C,D$ 变成 $2\times2$ 矩阵而已——这就是状态空间处理 MIMO 系统的优势。</p>

<h2><span class="sec">3.1.2</span>与传递函数的关系</h2>
<p>对状态方程做拉普拉斯变换（零初值）：$sZ = AZ + BU$，所以 $Z=(sI-A)^{-1}BU$，再代入输出方程：</p>
<div class="eq" data-no="(3.1.18)" data-title="从状态空间到传递函数">
$$G(s) = C\,(sI-A)^{-1}B + D = \frac{C\,(sI-A)^{*}B}{\left|sI-A\right|} + D$$
<dl><dt>$(sI-A)^{*}$</dt><dd>伴随矩阵</dd><dt>$|sI-A|$</dt><dd>行列式——它就是传递函数的分母</dd></dl>
<div class="say">令 $|sI-A|=0$ 解出的 $s$，既是传递函数的<b>极点</b>，又是矩阵 A 的<b>特征值</b>。两种描述方法在这里会师。</div>
</div>
<div class="sim" data-id="ssview"></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>状态空间方程用一组一阶微分方程描述系统：$\dot z=Az+Bu$，$y=Cz+Du$。</li>
<li>二阶机械系统的自然状态选择是“位置 + 速度”。</li>
<li>$G(s)=C(sI-A)^{-1}B+D$；A 的特征值就是传递函数的极点，分析 A 的特征值就能判断系统表现。</li>
</ol></div>`,
    sims: {
      ssview(el) {
        const ui = card(el, { title: "同一个系统的三种样子：动画、相平面、特征值", hint: "调 m、b、k，看 A 矩阵怎样变化、特征值（= 极点）跑到哪里、相平面上的轨迹（位置 z₁ vs 速度 z₂）是什么形状。" });
        const sc = new Scene(ui.plotHost(), { height: 90 });
        const [a, b] = ui.row("", "");
        const ph = new Plot(a, { height: 280, x: [-1.6, 1.6], y: [-2.5, 2.5], xlabel: "z₁ 位置", ylabel: "z₂ 速度", title: "相平面" });
        const ps = new Plot(b, { height: 280, x: [-4, 1], y: [-3, 3], equal: true, xlabel: "σ", title: "A 的特征值（= G(s) 的极点）" });
        const sm = ui.slider({ label: "质量 m", min: 0.3, max: 3, step: 0.1, value: 1, onInput: restart });
        const sb = ui.slider({ label: "阻尼 b", min: 0, max: 5, step: 0.05, value: 0.6, onInput: restart });
        const sk = ui.slider({ label: "弹性系数 k", min: 0.3, max: 6, step: 0.1, value: 2, onInput: restart });
        const mat = document.createElement("div"); mat.className = "sim-note"; ui.controls.appendChild(mat);
        let z, H;
        function restart() { z = [1.2, 0]; H = [[...z]]; }
        restart();
        ph.drag({ down: (x, y) => { z = [x, y]; H = [[x, y]]; return false; } }); ph.cv.style.cursor = "crosshair";
        CE.loop(ui.root, dt => {
          const m = sm.get(), bb = sb.get(), k = sk.get(), A = [[0, 1], [-k / m, -bb / m]];
          for (let i = 0; i < 8; i++) z = rk4(s => M.mulv(A, s), z, 0, dt / 8);
          H.push([...z]); if (H.length > 900) H.shift();
          const c = sc.ctx, W = sc.W; sc.clear();
          const cx = W / 2 + z[0] * W * 0.2;
          c.fillStyle = COL.axis; c.fillRect(20, 15, 8, 60);
          c.strokeStyle = COL.c1; c.lineWidth = 2; c.beginPath(); c.moveTo(28, 45); for (let i = 1; i <= 14; i++) c.lineTo(28 + (cx - 50 - 28) * i / 15, 45 + (i % 2 ? -8 : 8)); c.lineTo(cx - 25, 45); c.stroke();
          c.fillStyle = COL.c2; c.beginPath(); c.roundRect(cx - 25, 22, 50, 46, 6); c.fill();
          c.strokeStyle = COL.c4; c.lineWidth = 3; c.beginPath(); c.moveTo(cx, 45); c.lineTo(cx + z[1] * 30, 45); c.stroke();
          ph.render = () => { ph.clear(); ph.line(H.map(p => p[0]), H.map(p => p[1]), COL.c1, 2); ph.dot(z[0], z[1], COL.c4, 5); };
          const e = M.eig2(A);
          ps.render = () => { ps.clear(); ps.rect(0, -3, 1, 3, COL.c4 + "12"); e.vals.forEach(v => ps.cross(v[0], v[1], COL.c4, 8, 3)); };
          ph.render(); ps.render();
        });
        function info() {
          const m = sm.get(), bb = sb.get(), k = sk.get(), A = [[0, 1], [-k / m, -bb / m]], e = M.eig2(A), cl = classify(A);
          mat.innerHTML = `$A=${mtex(A)}$`; CE.renderMath(mat);
          ui.out(`$G(s)=\\frac{1}{${f2(m, 1)}s^2+${f2(bb)}s+${f2(k, 1)}}$<br>特征值：<b>${cfmt(e.vals[0])}</b>，<b>${cfmt(e.vals[1])}</b><br>类型：<span class="${cl.cls}">${cl.name}</span><br><span class="sim-note">在相平面上点击可以设定新的初始状态。</span>`);
        }
        [sm, sb, sk].forEach(s => s.el.querySelector("input").addEventListener("input", info)); info();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>b=0：特征值落在虚轴上，相轨迹是闭合椭圆（永远振荡）。</li><li>慢慢增大 b：特征值向左移动，相轨迹变成向内的螺旋。</li><li>继续增大 b，直到两个特征值在实轴上“相遇”后分开：相轨迹不再绕圈，直接滑向原点——不振荡了。</li></ol></div>`);
      },
    },
  });

  /* ════════════════ 3.2 特征值与解耦 ════════════════ */
  Course.add({
    id: "v1-eig", vol: 1, ord: 60, no: "3.2", chName: "第3章 动态系统建模——状态空间方程", ref: "卷1 §3.2",
    title: "特征值与特征向量：找到系统的“天然方向”",
    sub: "矩阵会把大多数向量转个方向，但有几个特殊方向只被拉长或压短——沿着这些方向，耦合的系统就被拆开了。",
    goals: ["直观理解 $Av=\\lambda v$", "会手算 2×2 矩阵的特征值和特征向量", "理解用 $P^{-1}AP=D$ 把耦合系统解耦"],
    minutes: 25,
    html: String.raw`
<p>状态方程的解析解是 $\bm z(t)=e^{A(t-t_0)}\bm z(t_0)+\int e^{A(t-\tau)}B\bm u\,\mathrm{d}\tau$，里面有矩阵指数，直接算很难。书中的策略是：先借助<span class="kw">特征值与特征向量</span>，把耦合的系统拆开，再用图形化的相平面去理解它。</p>

<h2><span class="sec">3.2.1</span>特征向量：变换前后仍在同一条直线上</h2>
<div class="eq" data-no="(3.2.2)" data-title="特征值与特征向量的定义">
$$A\,\bm v = \lambda\,\bm v$$
<dl><dt>$\bm v$</dt><dd>特征向量：被 A 作用后方向不变（或正好反向）的非零向量</dd><dt>$\lambda$</dt><dd>特征值：沿这个方向被拉伸的倍数（负数表示反向）</dd></dl>
</div>
<p>书中的例子 $A=\begin{bmatrix}1&1\\4&-2\end{bmatrix}$：取 $\bm v_a=[1,2]^T$，得 $A\bm v_a=[3,0]^T$，方向变了；取 $\bm v_b=[1,1]^T$，得 $A\bm v_b=[2,2]^T=2\bm v_b$，方向不变、长度翻倍——所以 $\bm v_b$ 是特征向量，$\lambda=2$。</p>
<div class="sim" data-id="eighunt"></div>
<p><b>怎么求</b>：$A\bm v=\lambda\bm v \Rightarrow (A-\lambda I)\bm v = 0$。要有非零解，矩阵 $A-\lambda I$ 必须“压扁空间”，即行列式为零：</p>
<div class="eq" data-no="(3.2.7)–(3.2.9)" data-title="特征方程">
$$|A-\lambda I| = \begin{vmatrix}1-\lambda & 1\\ 4 & -2-\lambda\end{vmatrix} = \lambda^2+\lambda-6 = 0 \;\Rightarrow\; \lambda_1 = 2,\;\lambda_2=-3$$
</div>
<p>再把每个 $\lambda$ 代回 $(A-\lambda I)\bm v=0$：$\lambda_1=2$ 时得 $v_{11}=v_{12}$，取 $\bm v_1=[1,1]^T$；$\lambda_2=-3$ 时得 $4v_{21}+v_{22}=0$，取 $\bm v_2=[0.5,-2]^T$。特征向量其实是一整条直线，取其中任意一个非零向量都行。</p>

<h2><span class="sec">3.2.2</span>应用：把耦合的方程组解耦</h2>
<p>考虑 $\dot z_1 = z_1+z_2,\;\dot z_2 = 4z_1-2z_2$。$z_1$ 的变化依赖 $z_2$，$z_2$ 的变化又依赖 $z_1$——它们<span class="kw">耦合</span>在一起，单独看哪一个都看不清。办法是换一套坐标，让新坐标轴正好沿着特征向量：</p>
<div class="eq" data-no="(3.2.16)–(3.2.24)" data-title="对角化解耦">
$$P=[\bm v_1,\;\bm v_2],\quad P^{-1}AP = D=\begin{bmatrix}\lambda_1&0\\0&\lambda_2\end{bmatrix},\quad \bm z = P\bar{\bm z}\;\Rightarrow\; \dot{\bar z}_1=\lambda_1\bar z_1,\;\; \dot{\bar z}_2=\lambda_2\bar z_2$$
<dl><dt>$P$</dt><dd>过渡矩阵：以特征向量为列</dd><dt>$\bar{\bm z}$</dt><dd>新坐标（沿特征向量方向的分量）</dd><dt>$D$</dt><dd>特征值组成的对角阵</dd></dl>
<div class="say">在新坐标里，每个分量只和自己有关：$\bar z_1=C_1e^{2t}$，$\bar z_2=C_2e^{-3t}$。再变回去：$\bm z = C_1e^{2t}\bm v_1 + C_2e^{-3t}\bm v_2$。</div>
</div>
<details class="derive"><summary>为什么 $P^{-1}AP=D$？</summary><div class="body">
<div class="step">$AP = [A\bm v_1,\;A\bm v_2] = [\lambda_1\bm v_1,\;\lambda_2\bm v_2] = [\bm v_1,\bm v_2]\begin{bmatrix}\lambda_1&0\\0&\lambda_2\end{bmatrix} = PD$。</div>
<div class="step">两边左乘 $P^{-1}$：$P^{-1}AP=D$。</div>
<div class="step">把 $\bm z=P\bar{\bm z}$ 代入 $\dot{\bm z}=A\bm z$：$P\dot{\bar{\bm z}}=AP\bar{\bm z}$，左乘 $P^{-1}$ 得 $\dot{\bar{\bm z}}=D\bar{\bm z}$。</div>
</div></details>
<div class="idea">解耦后答案一目了然：特征值 2 和 −3 出现在指数上。只要初值在 $\bm v_1$ 方向有一点分量（$C_1\neq0$），$e^{2t}$ 就会无限增长——系统发散。<b>特征值决定了系统沿各个“天然方向”是增长还是衰减。</b></div>
<div class="sim" data-id="decouple"></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>特征向量是被矩阵作用后方向不变的向量，特征值是拉伸倍数；由 $|A-\lambda I|=0$ 求出。</li>
<li>以特征向量为列组成 P，可把 A 对角化：$P^{-1}AP=D$，从而把耦合系统拆成互不相干的一阶系统。</li>
<li>解 $\bm z(t)=\sum C_ie^{\lambda_i t}\bm v_i$：特征值决定各方向的增长/衰减。</li>
</ol></div>`,
    sims: {
      eighunt(el) {
        const ui = card(el, { title: "寻找“不变方向”", hint: "拖动蓝色向量 v，橙色是 Av。当两者共线时会高亮——你就找到了一个特征向量。打开“扫描”看 v 转一圈时 Av 怎么走。" });
        const pl = new Plot(ui.plotHost(), { height: 360, x: [-4, 4], y: [-4, 4], equal: true, xlabel: "x", ylabel: "y" });
        const presets = { book: [[1, 1], [4, -2]], sym: [[2, 1], [1, 2]], rot: [[0, -1], [1, 0]], shear: [[1, 1], [0, 1]] };
        let A = presets.book, v = [1, 2], sweep = false, ang = 0, found = [];
        ui.seg({ label: "矩阵 A", options: [["book", "书中例子"], ["sym", "对称 [[2,1],[1,2]]"], ["rot", "旋转 90°"], ["shear", "剪切"]], value: "book", onChange: k => { A = presets[k]; found = []; draw(); } });
        ui.check({ label: "自动扫描（v 绕单位圆转动）", value: false, onChange: x => sweep = x });
        pl.drag({ down: (x, y) => { v = [x, y]; draw(); }, move: (x, y) => { v = [clamp(x, -4, 4), clamp(y, -4, 4)]; draw(); }, hover: () => true });
        function draw() {
          const Av = M.mulv(A, v), cr = v[0] * Av[1] - v[1] * Av[0], nv = Math.hypot(...v), na = Math.hypot(...Av);
          const col = nv > 0.05 && na > 0.05 && Math.abs(cr) / (nv * na) < 0.04;
          if (col && sweep) { const key = (Math.atan2(v[1], v[0]) * 18 / Math.PI).toFixed(0); if (!found.includes(key)) found.push(key); }
          pl.render = () => {
            pl.clear();
            // 单位圆在 A 下的像：椭圆
            const ts = linspace(0, 2 * Math.PI, 120), im = ts.map(t => M.mulv(A, [Math.cos(t), Math.sin(t)]));
            pl.line(ts.map(Math.cos), ts.map(Math.sin), COL.axis, 1.2, [4, 4]); pl.line(im.map(p => p[0]), im.map(p => p[1]), COL.c2 + "66", 1.2, [4, 4]);
            const e = M.eig2(A); if (e.real) e.vecs.forEach(q => pl.line([-q[0] * 9, q[0] * 9], [-q[1] * 9, q[1] * 9], col ? COL.c3 + "aa" : COL.grid, 1.4));
            pl.arrow(0, 0, Av[0], Av[1], col ? COL.c3 : COL.c2, 3, 11);
            pl.arrow(0, 0, v[0], v[1], col ? COL.c3 : COL.c1, 3, 11);
            pl.text(v[0] * 1.08, v[1] * 1.08 + 0.15, "v", COL.c1, "center", 14, true);
            pl.text(Av[0] * 1.06, Av[1] * 1.06 - 0.35, "Av", COL.c2, "center", 14, true);
          };
          pl.render();
          const lam = col ? (M.dot(Av, v) / M.dot(v, v)) : null;
          const e = M.eig2(A);
          ui.out(`v = [${f2(v[0])}, ${f2(v[1])}]<br>Av = [${f2(Av[0])}, ${f2(Av[1])}]<br>${col ? `<span class="good">✔ 共线！λ = ${f2(lam)}</span>` : "还不共线，继续拖动"}<br><span class="sim-note">该矩阵的特征值：${cfmt(e.vals[0])}，${cfmt(e.vals[1])}${e.real ? "" : "（复数：没有实的不变方向，所有向量都被转动）"}</span>`);
        }
        CE.loop(ui.root, dt => { if (!sweep) return; ang += dt * 0.6; const r = 2; v = [r * Math.cos(ang), r * Math.sin(ang)]; draw(); });
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>书中矩阵：把 v 拖到 [1,1] 方向（λ=2，Av 是 v 的两倍）和 [0.5,−2] 方向（λ=−3，Av 反向且变长）。</li><li>选“旋转 90°”：无论 v 指向哪里，Av 都转了 90°——找不到实特征向量，特征值是 ±j。复特征值意味着“转圈”，这正是振荡的来源。</li></ol></div>`);
      },
      decouple(el) {
        const ui = card(el, { title: "戴上“解耦眼镜”：坐标变换动画", hint: "以例 3.3.2 $A=\\begin{bmatrix}-3&4\\\\-2&3\\end{bmatrix}$（λ=−1, 1）为例。滑块从 0 到 1，坐标系从“特征坐标 z̄”平滑变到原坐标 z。" });
        const pl = new Plot(ui.plotHost(), { height: 360, x: [-3, 3], y: [-3, 3], equal: true });
        const Pm = [[2, 1], [1, 1]];
        const s = ui.slider({ label: "变换进度 s（0 = 解耦坐标，1 = 原坐标）", min: 0, max: 1, step: 0.01, value: 0, onInput: draw });
        ui.buttons([["▶ 播放", () => { let k = 0; const id = setInterval(() => { k += 0.02; s.set(Math.min(1, k)); draw(); if (k >= 1) clearInterval(id); }, 30); }]]);
        function draw() {
          const t = s.get(), T = [[1 + (Pm[0][0] - 1) * t, Pm[0][1] * t], [Pm[1][0] * t, 1 + (Pm[1][1] - 1) * t]];
          const map = p => M.mulv(T, p);
          pl.render = () => {
            pl.clear();
            const ax1 = map([4, 0]), ax2 = map([0, 4]);
            pl.line([-ax1[0], ax1[0]], [-ax1[1], ax1[1]], COL.c2, 2, [7, 5]); pl.line([-ax2[0], ax2[0]], [-ax2[1], ax2[1]], COL.c5, 2, [7, 5]);
            pl.text(ax1[0] * 0.62, ax1[1] * 0.62 + 0.2, t > 0.5 ? "v₁=[2,1]（收敛 λ=−1）" : "z̄₁ 轴（收敛）", COL.c2);
            pl.text(ax2[0] * 0.6 + 0.1, ax2[1] * 0.6, t > 0.5 ? "v₂=[1,1]（发散 λ=1）" : "z̄₂ 轴（发散）", COL.c5);
            // 在 z̄ 坐标中的鞍点轨迹：z̄ = [C1 e^{-t}, C2 e^{t}]
            for (const c1 of [-2, -1, 1, 2]) for (const c2 of [-0.06, 0.06, -0.3, 0.3]) {
              const ts = linspace(-0.2, 4, 80), pts = ts.map(tt => map([c1 * Math.exp(-tt), c2 * Math.exp(tt)]));
              pl.line(pts.map(p => p[0]), pts.map(p => p[1]), COL.c1 + "bb", 1.5);
              const k = 40; pl.arrow(pts[k - 2][0], pts[k - 2][1], pts[k][0], pts[k][1], COL.c1, 1.2, 7);
            }
          };
          pl.render();
          ui.out(t < 0.02 ? "解耦坐标：$\\dot{\\bar z}_1=-\\bar z_1$，$\\dot{\\bar z}_2=\\bar z_2$<br>标准鞍点：沿横轴进、沿纵轴出" : t > 0.98 ? "原坐标：$\\bm z = C_1e^{-t}\\begin{bmatrix}2\\\\1\\end{bmatrix}+C_2e^{t}\\begin{bmatrix}1\\\\1\\end{bmatrix}$<br>轨迹被“拉长、压扁、旋转”，但仍然是鞍点" : "坐标轴正在旋转到特征向量方向……");
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>关键观察</b>：线性变换会扭曲轨迹的形状，但<b>不会改变平衡点的性质</b>——鞍点还是鞍点。所以分析对角化后的简单系统，就等于分析了原系统。</div>`);
        CE.renderMath(ui.root);
      },
    },
  });

  /* ════════════════ 3.3 相平面与相轨迹 ════════════════ */
  const PRESETS = {
    e3214: { A: [[1, 1], [4, -2]], name: "例 3.2.14 解耦例" },
    e332: { A: [[-3, 4], [-2, 3]], name: "例 3.3.2 鞍点" },
    e333: { A: [[3, -2], [1, 0]], name: "例 3.3.3 不稳定节点" },
    e334: { A: [[-3, 2], [-1, 0]], name: "例 3.3.4 稳定节点" },
    e335: { A: [[0, 4], [-1, 0]], name: "例 3.3.5 中心点" },
    e336: { A: [[1, -2], [2, 1]], name: "例 3.3.6 不稳定焦点" },
    e337: { A: [[-4, 2], [-5, -2]], name: "例 3.3.7 稳定焦点" },
  };
  Course.add({
    id: "v1-phase", vol: 1, ord: 70, no: "3.3", chName: "第3章 动态系统建模——状态空间方程", ref: "卷1 §3.3",
    title: "相平面与相轨迹：一张图看懂系统往哪走",
    sub: "不解方程，只画箭头：每一点上状态往哪儿跑、跑多快。平衡点是稳定还是不稳定，一眼可见。",
    goals: ["会用一维相线判断平衡点稳定性", "认识二维平衡点的六种类型：节点、鞍点、焦点、中心", "掌握“实部定稳定、虚部定振荡”"],
    minutes: 30,
    html: String.raw`
<p><span class="kw">相轨迹</span>是状态变量随时间在“状态坐标系”（相平面）里画出的轨迹。它不需要解出方程，只需要知道每一点的变化率，就能把系统的行为画出来，而且对非线性系统同样适用。</p>

<h2><span class="sec">3.3.1</span>一维相线</h2>
<p>例：$\dot z = z^2-1$。以 $z$ 为横轴、$\dot z$ 为纵轴画出抛物线，它和横轴交于 $z=-1$ 和 $z=1$——这两点 $\dot z=0$，叫<span class="kw">平衡点</span>。在横轴上看 $\dot z$ 的正负：$\dot z>0$ 就往右走，$\dot z<0$ 就往左走。</p>
<ul><li>$z=-1$ 两侧的箭头都指向它：<b>稳定平衡点</b>（局部稳定：初值必须小于 1 才能回来）。</li>
<li>$z=1$ 两侧的箭头都背离它：<b>不稳定平衡点</b>。</li></ul>
<p>同样的方法可以分析非线性的 <b>Logistic 人口模型</b> $\dot P = rP(1-P/K)$：$P=0$ 不稳定（只要有一点生命就会生生不息），$P=K$（环境承载力）稳定；增长速度在 $P=K/2$ 时最快，所以人口曲线是 S 形。</p>
<div class="sim" data-id="phase1d"></div>
<div class="tip" data-t="书中的人文注脚">从少量人口增长到 K 需要几代人，而从远超 K 回落到 K 却很快——繁衍需要几十年，灾难却可能在短时间内发生。模型虽简单，却提醒人们敬畏环境承载力。</div>

<h2><span class="sec">3.3.2</span>二维相平面：从对角矩阵开始</h2>
<p>先看已经解耦的情况 $\dot z_1=a_{11}z_1,\;\dot z_2=a_{22}z_2$。分析套路：先看坐标轴上的点怎么动，再推广到整个平面。</p>
<table class="tbl"><tr><th>情形</th><th>坐标轴上的运动</th><th>平衡点类型</th></tr>
<tr><td>$a_{11}>0,\;a_{22}>0$</td><td>两个方向都远离原点</td><td>不稳定节点</td></tr>
<tr><td>$a_{11}>0,\;a_{22}<0$</td><td>沿 $z_2$ 进来，沿 $z_1$ 出去</td><td>鞍点（不稳定）</td></tr>
<tr><td>$a_{11}<0,\;a_{22}<0$</td><td>两个方向都靠近原点</td><td>稳定节点</td></tr></table>
<p>两个系数大小不同时轨迹是曲线：靠近原点时它贴着“慢”的方向（$|a|$ 小的那个）走。</p>

<h2><span class="sec">3.3.3</span>一般形式：先解耦，再变回来</h2>
<p>对一般的 $A$，用上一课的方法 $\bm z=P\bar{\bm z}$ 解耦，在特征坐标里套用对角情形的结论，再映射回原坐标——轨迹会沿特征向量方向被拉伸、旋转，但平衡点性质不变。复特征值时没有实特征向量，解里出现 $e^{\sigma t}(\cos\omega t,\sin\omega t)$：轨迹绕原点<b>转圈</b>，$\sigma$ 决定圈越转越大还是越转越小。</p>
<div class="sim" data-id="phase2d"></div>
<div class="eq" data-title="表 3.3.1　特征值 $\lambda_{1,2}=\sigma\pm j\omega$ 与平衡点类型">
<table class="tbl" style="margin:0"><tr><th>特征值</th><th>条件</th><th>类型</th></tr>
<tr><td rowspan="3">实数（ω=0）</td><td>都为负</td><td>稳定节点</td></tr><tr><td>一正一负</td><td>鞍点</td></tr><tr><td>都为正</td><td>不稳定节点</td></tr>
<tr><td rowspan="3">复数（ω≠0）</td><td>纯虚数 σ=0</td><td>中心点</td></tr><tr><td>σ&gt;0</td><td>不稳定焦点</td></tr><tr><td>σ&lt;0</td><td>稳定焦点</td></tr></table>
<div class="say">实部决定稳定性，虚部决定是否振荡。（小技巧：$\lambda_1\lambda_2=\det A$，$\lambda_1+\lambda_2=\mathrm{tr}A$。）</div>
</div>

<div class="quiz"><div class="q">$A=\begin{bmatrix}-1 & 3\\ -3 & -1\end{bmatrix}$ 的平衡点是什么类型？</div>
<div class="opts"><button class="opt">鞍点</button><button class="opt" data-ok>稳定焦点</button><button class="opt">中心点</button></div>
<div class="why">$\mathrm{tr}A=-2$，$\det A=1+9=10$，特征方程 $\lambda^2+2\lambda+10=0$，$\lambda=-1\pm3j$：实部负 → 稳定；虚部非零 → 螺旋。可以在上面的实验室里输入验证。</div></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>一维：在横轴上看 $\dot z$ 的符号，箭头汇聚为稳定、发散为不稳定。</li>
<li>二维线性系统的平衡点由 A 的特征值决定：节点、鞍点、焦点、中心。</li>
<li>实部 → 稳定性；虚部 → 振荡。线性变换不改变平衡点性质。</li>
</ol></div>`,
    sims: {
      phase1d(el) {
        const ui = card(el, { title: "一维相线探索器", hint: "上图：$\\dot z$ 关于 $z$ 的曲线；<b>点击横轴附近</b>放一个小球，它会按 $\\dot z$ 的符号移动。下图：对应的 $z(t)$。" });
        const p1 = new Plot(ui.plotHost(), { height: 230, x: [-2.5, 2.5], y: [-1.5, 3], xlabel: "z", ylabel: "ż" });
        const p2 = new Plot(ui.plotHost(), { height: 200, x: [0, 8], y: [-2.5, 2.5], xlabel: "t", ylabel: "z(t)" });
        let mode = "quad", balls = [];
        const sr = ui.slider({ label: "Logistic：增长率 r", min: 0.2, max: 2, step: 0.05, value: 1, onInput: setup });
        const sK = ui.slider({ label: "Logistic：承载力 K", min: 50, max: 150, step: 5, value: 100, fmt: v => v.toFixed(0), onInput: setup });
        ui.seg({ label: "方程", options: [["quad", "ż = z² − 1"], ["log", "Logistic 人口"]], value: mode, onChange: v => { mode = v; setup(); } });
        ui.buttons([["清空", () => { balls = []; }, "ghost"]]);
        const f = z => mode === "quad" ? z * z - 1 : sr.get() * z * (1 - z / sK.get());
        function setup() {
          balls = [];
          if (mode === "quad") { p1.setRange([-2.5, 2.5], [-1.5, 3]); p2.setRange([0, 8], [-2.5, 2.5]); }
          else { const K = sK.get(), r = sr.get(); p1.setRange([-10, 1.8 * K], [-0.6 * r * K, 0.4 * r * K]); p2.setRange([0, 14], [0, 1.8 * K]); }
        }
        p1.drag({ down: x => { balls.push({ z: x, t: 0, H: [[0, x]] }); if (balls.length > 8) balls.shift(); return false; } }); p1.cv.style.cursor = "crosshair";
        setup();
        CE.loop(ui.root, dt => {
          balls.forEach(b => { if (b.t < p2.opt.x[1]) { for (let i = 0; i < 10; i++) { b.z += f(b.z) * dt / 10; } b.t += dt; if (Math.abs(b.z) < 1e4) b.H.push([b.t, b.z]); } });
          p1.render = () => {
            p1.clear(); const [a, b] = p1.opt.x; p1.fn(f, COL.c1, 2.4);
            const eq = mode === "quad" ? [[-1, true], [1, false]] : [[0, false], [sK.get(), true]];
            for (let i = 0; i < 14; i++) { const z = a + (b - a) * (i + 0.5) / 14, v = f(z); if (Math.abs(v) < 1e-6) continue; const L = (b - a) / 30; p1.arrow(z - Math.sign(v) * L / 2, 0, z + Math.sign(v) * L / 2, 0, v > 0 ? COL.c3 : COL.c4, 2, 7); }
            eq.forEach(([z, s]) => { if (s) p1.dot(z, 0, COL.ink, 6); else { p1.dot(z, 0, COL.plotbg, 6, COL.ink); } });
            if (mode === "log") { p1.vline(sK.get() / 2, COL.muted, [3, 3]); p1.text(sK.get() / 2, p1.opt.y[1] * 0.9, "K/2：增长最快", COL.muted, "center"); }
            balls.forEach(b => p1.dot(b.z, 0, COL.c2, 6));
          };
          p2.render = () => { p2.clear(); balls.forEach(b => p2.line(b.H.map(p => p[0]), b.H.map(p => p[1]), COL.c2, 2)); if (mode === "log") p2.hline(sK.get(), COL.c3, [6, 4]); else { p2.hline(-1, COL.c3, [6, 4]); p2.hline(1, COL.c4, [6, 4]); } };
          p1.render(); p2.render();
        });
        ui.out("实心点 ● = 稳定平衡点<br>空心点 ○ = 不稳定平衡点<br>横轴上的绿色箭头：ż&gt;0 向右；红色：ż&lt;0 向左");
      },
      phase2d(el) {
        const ui = card(el, { title: "2×2 相图实验室", hint: "红色小点随矢量场流动；蓝线是相轨迹；虚线是实特征向量。<b>点击相平面</b>可添加一条从该点出发的轨迹。" });
        const [a, b] = ui.row("", "");
        const ph = phasePlot(a, { height: 360 });
        const td = new Plot(b, { height: 360, x: [-7, 7], y: [-10, 26], xlabel: "tr A", ylabel: "det A", title: "迹–行列式平面（拖动黄点）" });
        let A = M.clone(PRESETS.e337.A);
        const mk = (i, j, lab) => ui.slider({ label: lab, min: -5, max: 5, step: 0.1, value: A[i][j], onInput: v => { A[i][j] = v; update(true); } });
        const s = [mk(0, 0, "$a_{11}$"), mk(0, 1, "$a_{12}$"), mk(1, 0, "$a_{21}$"), mk(1, 1, "$a_{22}$")];
        const sel = document.createElement("select"); sel.className = "seg"; sel.style.cssText = "padding:6px;border-radius:8px;border:1px solid var(--line);background:var(--panel);color:var(--ink);font:inherit;font-size:.85rem";
        sel.innerHTML = Object.entries(PRESETS).map(([k, v]) => `<option value="${k}" ${k === "e337" ? "selected" : ""}>${v.name}</option>`).join("");
        sel.onchange = () => { A = M.clone(PRESETS[sel.value].A); s[0].set(A[0][0]); s[1].set(A[0][1]); s[2].set(A[1][0]); s[3].set(A[1][1]); update(true); };
        ui.controls.prepend(sel);
        function update(reseed) {
          ph.A = A; if (reseed) ph.reset();
          const cl = classify(A), e = cl.e;
          td.render = () => {
            td.clear();
            const xs = linspace(-7, 7, 200); td.line(xs, xs.map(t => t * t / 4), COL.c5, 1.6);
            td.text(-6.8, 22, "抛物线 det = tr²/4 以上：复特征值（焦点）", COL.c5); td.text(-6.8, -8, "det < 0：鞍点", COL.c4);
            td.text(-5.5, 5, "稳定", COL.c3, "left", 12, true); td.text(4.2, 5, "不稳定", COL.c4, "left", 12, true);
            td.vline(0, COL.muted, [4, 3]); td.text(0.15, 17, "tr=0：中心", COL.muted);
            td.dot(e.tr, e.det, "#e8b400", 7, COL.ink);
          };
          td.render();
          ui.out(`$A=${mtex(A, 1)}$<br>特征值：<b>${cfmt(e.vals[0])}</b>，<b>${cfmt(e.vals[1])}</b><br>类型：<span class="${cl.cls}" style="font-size:1.05rem">${cl.name}</span>`);
        }
        td.drag({
          down: () => { }, move: (tr, det) => {
            // 保持 a12 = 1，用 tr/det 反解一个伴随型矩阵：[[0,1],[-det,tr]]
            A = [[0, 1], [-det, tr]]; s[0].set(0); s[1].set(1); s[2].set(-det); s[3].set(tr); update(true);
          }, hover: () => true,
        });
        update(true);
        CE.loop(ui.root, dt => { ph.tick(dt); ph.pl.render(); });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>依次选择书中 7 个例题，对照表 3.3.1。</li><li>在右图把黄点拖过抛物线：从“节点”变成“焦点”，轨迹开始转圈。</li><li>拖过横轴（det&lt;0）：变成鞍点；拖过纵轴 tr=0：稳定 ↔ 不稳定。</li></ol></div>`);
      },
    },
  });

  /* ════════════════ 3.4 爱情故事 ════════════════ */
  Course.add({
    id: "v1-love", vol: 1, ord: 80, no: "3.4", chName: "第3章 动态系统建模——状态空间方程", ref: "卷1 §3.4",
    title: "相平面案例：雨飞与梦寒的爱情故事",
    sub: "用两个微分方程描述两个人的感情，相轨迹会告诉你他们的结局——以及“第一印象”有多重要。",
    goals: ["练习用相平面分析一个“非物理”的系统", "体会特征值、特征向量在故事里的含义", "学会“描述系统 → 数学分析 → 结果讨论”三步法"],
    minutes: 20,
    html: String.raw`
<div class="story">这个模型最早由康奈尔大学的 Steven Strogatz 引入教学，原型是罗密欧与朱丽叶。书中换成了两位中国朋友：男生<b>雨飞</b>、女生<b>梦寒</b>。用 $Y(t)$ 表示雨飞对梦寒的感情，$M(t)$ 表示梦寒对雨飞的感情：大于 0 是爱，小于 0 是怨恨，等于 0 是无感。状态向量 $\bm z=[Y,M]^T$。</div>
<p>书中顺便示范了一种写论文、做报告都适用的<b>动态系统分析三步法</b>：① 用语言描述系统；② 用数学工具分析；③ 对结果进行讨论。</p>

<h2>情况 1：耿直的他 + 若即若离的她</h2>
<div class="eq" data-title="模型">
$$\dot Y = aM,\qquad \dot M = -bY,\qquad A=\begin{bmatrix}0&a\\-b&0\end{bmatrix},\quad a,b>0$$
</div>
<p><b>① 描述</b>：雨飞投桃报李——她对他好（$M>0$），他的好感就增加；她讨厌他，他也冷下来。梦寒则欲迎还拒——他越热情（$Y>0$），她越疏远；他一冷淡，她反而想起他的好。</p>
<p><b>② 分析</b>：特征值 $\lambda=\pm j\sqrt{ab}$，纯虚数 → <b>中心点</b>，相轨迹是顺时针的椭圆。</p>
<p><b>③ 讨论</b>：永远循环。四分之一的时间相爱（第一象限），四分之一的时间互相看不顺眼（第三象限），一半时间“一半火焰一半海水”。作者认为这其实是很多情侣的常态——即使在冷战时，也充满了对未来的期待。</p>

<h2>情况 2：两个性格相同的人</h2>
<div class="eq" data-title="模型">
$$\dot Y = -aY + bM,\qquad \dot M = bY - aM,\qquad \lambda_1=-a+b,\;\lambda_2=-a-b,\quad \bm v_1=\begin{bmatrix}1\\1\end{bmatrix},\;\bm v_2=\begin{bmatrix}1\\-1\end{bmatrix}$$
<dl><dt>$a$</dt><dd>“自我意识/理性/克制”：对自己的感情有所保留</dd><dt>$b$</dt><dd>“对对方的回应/投入”：感受到爱就更爱，感受到恨就受挫</dd></dl>
</div>
<ul>
<li><b>$a>b$（都很理性）</b>：两个特征值都为负 → <b>稳定节点</b>。无论曾经爱得多深、恨得多深，最终都回到 $Y=M=0$——相忘于江湖，成为路人。</li>
<li><b>$a<b$（都很感性）</b>：$\lambda_1>0,\lambda_2<0$ → <b>鞍点</b>。沿 $[1,1]$ 方向发散：要么双双坠入爱河（第一象限），要么由爱生恨（第三象限）。结局由<b>初始状态</b>决定——分界线就是稳定特征方向所在的直线 $Y+M=0$。</li>
</ul>
<div class="sim" data-id="love"></div>
<div class="idea">鞍点情形里，初始点只要在分界线 $Y+M=0$ 上下稍微挪一点，结局就从“幸福生活”翻转成“成为仇人”。数学上，这是因为初值在不稳定特征向量 $\bm v_1$ 上的投影 $(Y_0+M_0)/2$ 的正负决定了 $e^{\lambda_1 t}$ 把它放大到哪一边。生活里，这叫“珍惜第一印象”。</div>
<div class="tip">书的结尾很温柔：不要试图用数学解决爱情——生活有太多无法建模的不确定性，正是这些不确定性让世界浪漫。少一些套路，多一点真诚。（两位朋友最后幸福地生活在了一起。）</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>纯虚特征值 → 中心点 → 无限循环；两负实特征值 → 稳定节点 → 归于平静；一正一负 → 鞍点 → 走向极端，结果取决于初值。</li>
<li>鞍点的分界线是稳定特征向量所在的直线。</li>
<li>动态系统分析三步法：描述系统 → 数学分析 → 结果与讨论。</li>
</ol></div>`,
    sims: {
      love(el) {
        const ui = card(el, { title: "爱情模拟器", hint: "<b>点击相平面</b>设定“第一印象” $(Y_0,M_0)$，看两个人的感情如何演化。上方小剧场随 Y、M 变化。" });
        const sc = new Scene(ui.plotHost(), { height: 110 });
        const [a, b] = ui.row("", "");
        const ph = phasePlot(a, { height: 320, x: [-3, 3], y: [-3, 3], xlabel: "Y 雨飞", ylabel: "M 梦寒" });
        ph.showParts = false;
        const pt = new Plot(b, { height: 320, x: [0, 10], y: [-3, 3], xlabel: "t", title: "感情随时间" });
        let mode = 1, z = [1.5, -1], H = [];
        ui.seg({ label: "情况", options: [[1, "情况1 耿直+若即若离"], [2, "情况2 性格相同"]], value: 1, onChange: v => { mode = v; setA(); restart(z0); } });
        const sa = ui.slider({ label: "a（情况1：他的回应 / 情况2：理性克制）", min: 0.1, max: 2, step: 0.05, value: 0.8, onInput: () => { setA(); restart(z0); } });
        const sb = ui.slider({ label: "b（情况1：她的反向 / 情况2：感情投入）", min: 0.1, max: 2, step: 0.05, value: 1.2, onInput: () => { setA(); restart(z0); } });
        let z0 = [1.5, -1];
        function setA() { const A_ = sa.get(), B_ = sb.get(); ph.A = mode === 1 ? [[0, A_], [-B_, 0]] : [[-A_, B_], [B_, -A_]]; ph.reset([]); }
        function restart(p) { z0 = p.slice(); z = p.slice(); H = [[0, ...z]]; ph.trajs = []; ph.add(p); }
        ph.pl.drag({ down: (x, y) => { restart([x, y]); return false; } });
        setA(); restart(z0);
        CE.loop(ui.root, dt => {
          const t = H[H.length - 1][0];
          if (t < 10 && Math.abs(z[0]) < 3.5 && Math.abs(z[1]) < 3.5) { for (let i = 0; i < 8; i++) z = rk4(s => M.mulv(ph.A, s), z, 0, dt / 8); H.push([t + dt, ...z]); }
          ph.pl.render(); ph.pl.dot(z[0], z[1], COL.c4, 7);
          if (mode === 2) ph.pl.line([-4, 4], [4, -4], COL.c3, 1.6, [3, 3]);
          pt.render = () => { pt.clear(); pt.line(H.map(h => h[0]), H.map(h => h[1]), COL.c1, 2.4); pt.line(H.map(h => h[0]), H.map(h => h[2]), COL.c4, 2.4); pt.legend([{ label: "Y 雨飞", color: COL.c1 }, { label: "M 梦寒", color: COL.c4 }]); };
          pt.render();
          // 小剧场
          const c = sc.ctx, W = sc.W; sc.clear();
          const Y = z[0], Mv = z[1], gap = 70 + (Math.max(0, -Y) + Math.max(0, -Mv)) * 30 - (Math.max(0, Y) + Math.max(0, Mv)) * 10;
          const cx = W / 2, x1 = cx - gap / 2, x2 = cx + gap / 2;
          const face = (x, v, col, lab) => {
            c.fillStyle = col; c.beginPath(); c.arc(x, 40, 18, 0, 7); c.fill();
            c.strokeStyle = "#fff"; c.lineWidth = 2; c.beginPath();
            if (v > 0.3) c.arc(x, 42, 8, 0.15 * Math.PI, 0.85 * Math.PI); else if (v < -0.3) c.arc(x, 54, 8, 1.15 * Math.PI, 1.85 * Math.PI); else { c.moveTo(x - 7, 47); c.lineTo(x + 7, 47); }
            c.stroke(); c.fillStyle = "#fff"; c.beginPath(); c.arc(x - 6, 35, 2.2, 0, 7); c.arc(x + 6, 35, 2.2, 0, 7); c.fill();
            c.fillStyle = COL.ink2; c.font = "12px 'Microsoft YaHei'"; c.textAlign = "center"; c.fillText(lab, x, 78);
            if (v > 0.5) { c.fillStyle = COL.c4; c.font = "16px serif"; c.fillText("♥", x + (x < cx ? 26 : -26), 22); }
          };
          face(x1, Y, COL.c1, "雨飞 " + f2(Y, 1)); face(x2, Mv, COL.c4, "梦寒 " + f2(Mv, 1));
          let msg = "";
          if (mode === 2) { const A_ = sa.get(), B_ = sb.get(); msg = A_ > B_ ? "a > b：终将成为路人" : (z0[0] + z0[1] > 0 ? "a < b 且 Y₀+M₀>0：走向爱情" : "a < b 且 Y₀+M₀<0：由爱生恨"); }
          else msg = Y > 0 && Mv > 0 ? "相爱中 ❤" : Y < 0 && Mv < 0 ? "互相看不顺眼" : "一半火焰，一半海水";
          c.fillStyle = COL.ink; c.font = "600 13px 'Microsoft YaHei'"; c.textAlign = "center"; c.fillText(msg, cx, 100);
        });
        ui.out("情况 2 中的绿色虚线是分界线 $Y+M=0$。<br>试着把初始点放在它的上方或下方一点点。");
      },
    },
  });
})();
