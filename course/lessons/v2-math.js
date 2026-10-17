/* 卷2 第1–2章：导读 / 状态方程的解 / 离散化 / 矩阵求导 / 线性回归 */
(function () {
  const { M, P, Plot, Scene, card, linspace, logspace, f2, COL, clamp, rk4, cfmt, mtex } = CE;

  Course.add({
    id: "v2-intro", vol: 2, ord: 10, no: "1 · 2.1", chName: "第1–2章 绪论与数学基础", ref: "卷2 第1章、§2.1",
    title: "卷 2 导读与状态方程的解",
    sub: "卷 1 教会我们让系统“稳定”；卷 2 要让它“又快又省”，并在看不清时也能估计出状态。先打好数学地基：矩阵指数。",
    goals: ["理解卷 2 要解决的两个新问题：最优控制与状态估计", "理解矩阵指数 $e^{At}$ 的定义与计算 $e^{At}=Pe^{Dt}P^{-1}$", "读懂状态方程的完整解：零输入响应 + 卷积"],
    minutes: 20,
    html: String.raw`
<div class="story" data-t="引子：万里穿针">2023 年 5 月 10 日，天舟六号货运飞船发射升空。这位携带 7.4 吨物资的“快递小哥”从 200 km 的近地轨道一路爬升到 393 km 的空间站轨道，在距空间站 5 km、400 m、200 m、19 m 四个停泊点逐步修正位置和姿态，最终与天宫空间站精准对接。整个过程依赖两件事：<b>精准的控制</b>，以及<b>对自身状态的精确估计</b>——这正是卷 2 的两个主题。</div>
<p>卷 1 结束时，我们已经能让一辆小车稳定地从起点移动到终点。但现代工程的要求更高：</p>
<ul><li>不仅要到，还要<b>尽快</b>到、<b>尽量省能量</b>——这就是<span class="kw">最优控制</span>（第 3–5 章：动态规划、LQR、MPC）。</li>
<li>模型有误差，传感器有噪声，看不清真实状态——需要<span class="kw">状态估计</span>（第 6 章：卡尔曼滤波器）。</li></ul>
<p>卷 2 统一使用状态空间方程。同一辆小车 $m\ddot x=u$，传递函数是 $\frac{1}{ms^2}$，只能看到输入和位置；状态空间形式 $\frac{\mathrm d}{\mathrm dt}\begin{bmatrix}x\\\dot x\end{bmatrix}=\begin{bmatrix}0&1\\0&0\end{bmatrix}\begin{bmatrix}x\\\dot x\end{bmatrix}+\begin{bmatrix}0\\1/m\end{bmatrix}u$ 则把速度也明确写了出来，便于对每个内部状态做优化和估计。</p>

<h2><span class="sec">2.1</span>状态空间方程的解</h2>
<p>一维时，$\dot x=ax$ 的解是 $x(t)=e^{at}x(0)$。多维时我们猜：$\bm x(t)=e^{At}\bm x(0)$。这个 $e^{At}$ 是什么？用泰勒级数定义：</p>
<div class="eq" data-no="(2.1.7)" data-title="矩阵指数">
$$e^{At}=I+At+\frac{1}{2!}(At)^2+\frac{1}{3!}(At)^3+\cdots,\qquad \frac{\mathrm d}{\mathrm dt}e^{At}=Ae^{At}$$
<div class="say">和标量指数完全一样的级数，只是把数换成了矩阵。逐项求导可以验证它满足 $\frac{\mathrm d}{\mathrm dt}e^{At}=Ae^{At}$。</div>
</div>
<p>怎么算？如果 A 是对角阵 $D=\mathrm{diag}(\lambda_i)$，那么 $e^{Dt}=\mathrm{diag}(e^{\lambda_it})$——对角元各自取指数。一般矩阵先对角化 $A=PDP^{-1}$，由 $A^n=PD^nP^{-1}$ 代入级数：</p>
<div class="eq" data-no="(2.1.15)">
$$e^{At}=P\,e^{Dt}\,P^{-1}=P\begin{bmatrix}e^{\lambda_1t}&&\\&\ddots&\\&&e^{\lambda_nt}\end{bmatrix}P^{-1}$$
</div>
<p>有输入时，两边乘 $e^{-At}$ 凑成 $\frac{\mathrm d}{\mathrm dt}(e^{-At}\bm x)=e^{-At}B\bm u$，积分得到完整的解：</p>
<div class="eq" data-no="(2.1.18)" data-title="LTI 状态方程的解">
$$\bm x(t)=\underbrace{e^{A(t-t_0)}\bm x(t_0)}_{\text{零输入响应}}+\underbrace{\int_{t_0}^{t}e^{A(t-\tau)}B\,\bm u(\tau)\,\mathrm d\tau}_{\text{输入的卷积}}$$
<dl><dt>$e^{A(t-t_0)}$</dt><dd><b>状态转移矩阵</b>（有些书记作 $\Phi$）：状态如何随时间自行演化</dd><dt>积分项</dt><dd>输入 u 与系统的卷积——卷 1 第 2 章的老朋友</dd></dl>
<div class="say">A 的特征值实部全小于 0 时，状态转移矩阵 → 0，初始状态的影响消失：这就是稳定性。</div>
</div>
<div class="sim" data-id="expm"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>卷 2 的目标：在稳定之上追求“最优”，并在噪声中估计状态。</li>
<li>$e^{At}=I+At+\frac{(At)^2}{2!}+\cdots=Pe^{Dt}P^{-1}$。</li>
<li>$\bm x(t)=e^{A(t-t_0)}\bm x(t_0)+\int e^{A(t-\tau)}B\bm u\,\mathrm d\tau$。</li>
</ol></div>`,
    sims: {
      expm(el) {
        const ui = card(el, { title: "矩阵指数：截断级数 vs 精确值", hint: "左：零输入轨迹 $\\bm x(t)=e^{At}\\bm x_0$（点击设置初值）。右：只取泰勒级数前 N 项时，$e^{At}$ 的近似误差随 t 的变化。" });
        const [a, b] = ui.row("", "");
        const ph = new Plot(a, { height: 280, x: [-2, 2], y: [-2, 2], equal: true, xlabel: "x₁", ylabel: "x₂", title: "相平面" });
        const pe = new Plot(b, { height: 280, x: [0, 4], y: [-8, 2], xlabel: "t", title: "log₁₀ 误差 ‖截断 − 精确‖" });
        const pres = { ex: [[0, 1], [-2, -3]], cart: [[0, 1], [0, 0]], osc: [[-0.3, 2], [-2, -0.3]], un: [[0.4, 1], [-1, 0.2]] };
        let A = pres.ex, x0 = [1.5, 1];
        ui.seg({ label: "矩阵 A", options: [["ex", "例 2.2.1 [[0,1],[−2,−3]]"], ["cart", "小车 [[0,1],[0,0]]"], ["osc", "衰减振荡"], ["un", "不稳定"]], value: "ex", onChange: k => { A = pres[k]; draw(); } });
        const sN = ui.slider({ label: "级数项数 N", min: 1, max: 20, step: 1, value: 4, fmt: v => v.toFixed(0), onInput: draw });
        ph.drag({ down: (x, y) => { x0 = [x, y]; draw(); return false; } }); ph.cv.style.cursor = "crosshair";
        const trunc = (A, t, N) => { let E = M.eye(2), T = M.eye(2); for (let k = 1; k < N; k++) { T = M.scale(M.mul(T, A), t / k); E = M.add(E, T); } return E; };
        function draw() {
          const ts = linspace(0, 4, 200), N = sN.get();
          const ex = ts.map(t => M.mulv(M.expm(M.scale(A, t)), x0)), tr = ts.map(t => M.mulv(trunc(A, t, N), x0));
          const err = ts.map(t => Math.log10(Math.max(1e-12, M.norm(M.sub(M.expm(M.scale(A, t)), trunc(A, t, N)).flat()))));
          ph.render = () => { ph.clear(); ph.line(tr.map(p => p[0]), tr.map(p => p[1]), COL.c2, 1.8, [5, 4]); ph.line(ex.map(p => p[0]), ex.map(p => p[1]), COL.c1, 2.6); ph.dot(x0[0], x0[1], COL.c4, 5); ph.legend([{ label: "精确 e^{At}x₀", color: COL.c1 }, { label: `前 ${N} 项`, color: COL.c2, dash: [5, 4] }]); };
          pe.render = () => { pe.clear(); pe.line(ts, err.map(e => clamp(e, -8, 2)), COL.c4, 2.2); };
          ph.render(); pe.render();
          const E1 = M.expm(A);
          ui.out(`特征值：${M.eig2(A).vals.map(v => cfmt(v)).join("，")}<br>$e^{A\\cdot1}=${mtex(E1, 3)}$<br><span class="sim-note">项数越多，近似在越长的时间内都准确；t 越大需要的项数越多。</span>`);
        }
        draw();
      },
    },
  });

  Course.add({
    id: "v2-disc", vol: 2, ord: 20, no: "2.2", chName: "第1–2章 绪论与数学基础", ref: "卷2 §2.2",
    title: "连续系统的离散化：计算机眼中的世界",
    sub: "计算机没有“连续时间”，只有一个个采样时刻。怎样选采样周期？怎样把 $\\dot x=Ax+Bu$ 变成 $x_{[k+1]}=Fx_{[k]}+Gu_{[k]}$？",
    goals: ["理解采样、ADC/DAC、零阶保持", "知道奈奎斯特采样定理与混叠", "会推导 ZOH 离散化：$F=e^{AT_s}$，$G=A^{-1}(F-I)B$"],
    minutes: 25,
    html: String.raw`
<p>现代控制器几乎都是数字的：被控对象（由物理定律决定）是连续的，计算机里的控制器却是离散的。传感器信号经 <b>ADC</b>（模拟→数字）按固定间隔 $T_s$（<span class="kw">采样周期</span>）读入计算机；控制器算出 $u_{[k]}$ 后，经 <b>DAC</b>（数字→模拟）送回执行器。</p>
<h2><span class="sec">2.2.1</span>采样周期怎么选？</h2>
<ul><li><b>太长</b>：丢失信息，甚至产生<span class="kw">混叠</span>——高频信号被采成一个“假的”低频信号。<b>奈奎斯特–香农采样定理</b>：采样频率至少是信号最高频率的 2 倍；工程上常取 5～10 倍。</li>
<li><b>太短</b>：数据量剧增；而且控制算法必须在一个采样周期内算完——算法需要 50 ms，那么小于 50 ms 的采样周期就没有意义。</li></ul>
<div class="tip" data-t="书中的类比：管理体重">一个月称一次体重，中间的变化全被忽略，计划跟不上实际；每 10 分钟称一次，体重根本来不及变，测到的全是重复信息，而且计划还没想好就又要称了。采样频率要兼顾“被控对象变化有多快”和“计划（计算）需要多久”。</div>
<div class="sim" data-id="alias"></div>
<p><b>零阶保持（ZOH）</b>：控制量 $u_{[k]}$ 在一个采样周期内保持不变，$u(t)=u_{[k]},\;kT_s\le t<(k+1)T_s$，形成阶梯状的信号——就像风扇转速在两次调节之间保持不变，而不是归零。</p>

<h2><span class="sec">2.2.2</span>状态空间方程的离散化</h2>
<p>从上一课的解出发，只看一个采样区间 $[t_k,t_{k+1}]$，并利用零阶保持把 u 提出积分：</p>
<div class="eq" data-no="(2.2.5)" data-title="零阶保持离散化">
$$\bm x_{[k+1]}=F\,\bm x_{[k]}+G\,\bm u_{[k]},\qquad F=e^{AT_s},\qquad G=\int_0^{T_s}e^{A\tau}\mathrm d\tau\,B=A^{-1}(F-I)B$$
<dl><dt>$F$</dt><dd>一个采样周期内，状态自己演化成什么样</dd><dt>$G$</dt><dd>一个采样周期内恒定的输入，对状态的累积影响</dd></dl>
<div class="say">这个离散化在采样时刻是<b>精确</b>的（不是近似）。$G=A^{-1}(F-I)B$ 要求 A 可逆；小车这类 A 不可逆的系统可直接计算积分（软件里用 c2d）。</div>
</div>
<p><b>例 2.2.1</b>：$A=\begin{bmatrix}0&1\\-2&-3\end{bmatrix}$，$B=\begin{bmatrix}0\\1\end{bmatrix}$，$C=[1\;0]$。$T_s=0.2$ s 时 $F=\begin{bmatrix}0.9671&0.1484\\-0.2968&0.5219\end{bmatrix}$，$G=\begin{bmatrix}0.01643\\0.1484\end{bmatrix}$；$T_s=1$ s 时 $F=\begin{bmatrix}0.6004&0.2325\\-0.4651&-0.0972\end{bmatrix}$，$G=\begin{bmatrix}0.1998\\0.2325\end{bmatrix}$。</p>
<div class="sim" data-id="c2d"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>采样频率 ≥ 2 倍信号最高频率（工程 5–10 倍），且要留够计算时间。</li>
<li>零阶保持：控制量在采样周期内保持不变。</li>
<li>$\bm x_{[k+1]}=F\bm x_{[k]}+G\bm u_{[k]}$，$F=e^{AT_s}$，$G=A^{-1}(F-I)B$——卷 2 所有离散算法的模型基础。</li>
</ol></div>`,
    sims: {
      alias(el) {
        const ui = card(el, { title: "混叠：采样太慢会“看错”", hint: "灰线：真实信号 sin(2πft)；橙点：采样值；蓝线：用采样点连起来“重建”的信号。" });
        const pl = new Plot(ui.plotHost(), { height: 230, x: [0, 4], y: [-1.4, 1.4], xlabel: "t / s" });
        const sf = ui.slider({ label: "信号频率 f (Hz)", min: 0.2, max: 5, step: 0.05, value: 1.8, onInput: draw });
        const ss = ui.slider({ label: "采样频率 fₛ (Hz)", min: 1, max: 30, step: 0.1, value: 2.2, onInput: draw });
        function draw() {
          const f = sf.get(), fs = ss.get(), ts = linspace(0, 4, 800), Ts = 1 / fs, sp = [];
          for (let t = 0; t <= 4 + 1e-9; t += Ts) sp.push([t, Math.sin(2 * Math.PI * f * t)]);
          pl.render = () => { pl.clear(); pl.line(ts, ts.map(t => Math.sin(2 * Math.PI * f * t)), COL.axis, 1.4); pl.line(sp.map(p => p[0]), sp.map(p => p[1]), COL.c1, 2.4); sp.forEach(p => pl.dot(p[0], p[1], COL.c2, 4.5)); };
          pl.render();
          const r = fs / f;
          ui.out(`fₛ / f = <b class="big">${f2(r, 2)}</b><br>${r < 2 ? '<span class="bad">低于 2 倍：混叠！重建出的是一个假的低频信号</span>' : r < 5 ? "满足奈奎斯特下限，但重建很粗糙" : '<span class="good">5 倍以上：重建良好</span>'}`);
        }
        draw();
      },
      c2d(el) {
        const ui = card(el, { title: "离散化对比器（例 2.2.1）", hint: "单位阶跃输入。蓝线：连续系统；橙色阶梯：ZOH 精确离散化；紫色阶梯：前向欧拉近似 $F\\approx I+AT_s$。" });
        const pl = new Plot(ui.plotHost(), { height: 280, x: [0, 6], y: [-0.1, 0.7], xlabel: "t / s", ylabel: "y" });
        const sT = ui.slider({ label: "采样周期 Tₛ (s)", min: 0.05, max: 1.5, step: 0.01, value: 0.2, onInput: draw });
        const cE = ui.check({ label: "显示前向欧拉离散化", value: true, onChange: draw });
        ui.buttons([["Tₛ = 0.2", () => { sT.set(0.2); draw(); }, "ghost"], ["Tₛ = 1", () => { sT.set(1); draw(); }, "ghost"]]);
        const A = [[0, 1], [-2, -3]], B = [[0], [1]];
        function stairs(T, F, G) { let x = [0, 0]; const xs = [], ys = []; for (let k = 0; k * T <= 6; k++) { xs.push(k * T, (k + 1) * T); ys.push(x[0], x[0]); x = M.vadd(M.mulv(F, x), M.flat(G)); if (Math.abs(x[0]) > 10) break; } return { xs, ys }; }
        function draw() {
          const T = sT.get(), { A: F, B: G } = M.c2d(A, B, T), Fe = M.add(M.eye(2), M.scale(A, T)), Ge = M.scale(B, T);
          const ts = linspace(0, 6, 300), yc = ts.map(t => 0.5 - Math.exp(-t) + 0.5 * Math.exp(-2 * t));
          const z = stairs(T, F, G), e = stairs(T, Fe, Ge);
          pl.render = () => { pl.clear(); pl.hline(0.5, COL.c3, [6, 4]); if (cE.get()) pl.line(e.xs, e.ys, COL.c5, 1.8); pl.line(z.xs, z.ys, COL.c2, 2.2); pl.line(ts, yc, COL.c1, 2.4); pl.legend([{ label: "连续", color: COL.c1 }, { label: "ZOH", color: COL.c2 }, ...(cE.get() ? [{ label: "欧拉", color: COL.c5 }] : [])], "tr"); };
          pl.render();
          const eig = M.eig2(Fe).vals.map(v => Math.hypot(v[0], v[1]));
          ui.out(`$F=e^{AT_s}=${mtex(F, 4)}$<br>$G=${mtex(G, 4)}$<br>${cE.get() && Math.max(...eig) > 1 ? '<span class="bad">欧拉离散化的特征值模 &gt; 1：离散系统发散了！</span>' : "ZOH 的台阶总是落在连续曲线上"}`);
        }
        draw();
      },
    },
  });

  // 线性回归数据（表 2.4.1）
  const Z = [183, 175, 187, 185, 176, 176, 185, 191, 195, 185, 174, 180, 178, 170, 184];
  const X = [75, 71, 83, 74, 73, 67, 79, 73, 88, 80, 81, 78, 73, 68, 71];
  const NAMES = ["孙兴慜", "萨拉赫", "罗纳尔多", "本泽马", "儒尼奥尔", "阿斯帕斯", "莱万多夫斯基", "希克", "哈兰德", "因莫比莱", "劳塔罗", "西蒙尼", "姆巴佩", "本耶德尔", "特雷尔"];
  Course.add({
    id: "v2-matdiff", vol: 2, ord: 30, no: "2.3–2.4", chName: "第1–2章 绪论与数学基础", ref: "卷2 §2.3–2.5",
    title: "矩阵求导与线性回归：最优化的工具箱",
    sub: "最优控制的本质是“求极值”，控制量又往往是向量。学会对向量、矩阵求导，并用它做一次线性回归——解析解和梯度下降各来一遍。",
    goals: ["理解分母布局下标量对向量的导数（梯度）", "记住 5 个常用公式：$\\frac{\\partial u^Tf}{\\partial u}=f$，$\\frac{\\partial Au}{\\partial u}=A^T$，$\\frac{\\partial u^TAu}{\\partial u}=2Au$ 等", "会用正规方程和梯度下降求线性回归"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">2.3.1</span>标量对向量求导</h2>
<p>求 $f(u)=u^2-2u+1$ 的最小值：令导数为 0，得 $u=1$。两个变量 $f(u_1,u_2)=u_1^2+u_2^2+2u_1$ 时，要令两个偏导数同时为 0。把 $\bm u=[u_1,u_2]^T$ 看成一个整体，就可以写成“标量对向量求导”：</p>
<div class="eq" data-no="(2.3.3)" data-title="分母布局（本书默认）">
$$\frac{\partial f(\bm u)}{\partial \bm u}\triangleq\begin{bmatrix}\frac{\partial f}{\partial u_1}\\\vdots\\\frac{\partial f}{\partial u_n}\end{bmatrix}\quad(n\times1)$$
<div class="say">结果的行数和分母 $\bm u$ 一样，所以叫“分母布局”——它就是梯度，指向函数增长最快的方向。另一种“分子布局”结果是它的转置；两种都可以，但必须前后一致。</div>
</div>
<p>向量 $\bm f$（$m\times1$）对向量 $\bm u$（$n\times1$）求导，分母布局下是 $n\times m$ 的矩阵，第 i 行第 j 列是 $\partial f_j/\partial u_i$——它是雅可比矩阵的转置。</p>
<div class="sim" data-id="grad"></div>

<h2><span class="sec">2.3.3</span>常用公式（分母布局）</h2>
<table class="tbl"><tr><th>编号</th><th>公式</th><th>标量类比</th></tr>
<tr><td>2.3.1</td><td>$\dfrac{\partial(\bm u^T\bm f)}{\partial\bm u}=\bm f$</td><td>$\frac{d(uf)}{du}=f$</td></tr>
<tr><td>2.3.2</td><td>$\dfrac{\partial(A\bm u)}{\partial\bm u}=A^T$</td><td>$\frac{d(au)}{du}=a$</td></tr>
<tr><td>2.3.3</td><td>$\dfrac{\partial(\bm u^TA\bm u)}{\partial\bm u}=A\bm u+A^T\bm u\;\overset{A\text{对称}}{=}\;2A\bm u$</td><td>$\frac{d(au^2)}{du}=2au$</td></tr>
<tr><td>2.3.4</td><td>$\dfrac{\partial^2(\bm u^TA\bm u)}{\partial\bm u^2}=A+A^T$</td><td>$\frac{d^2(au^2)}{du^2}=2a$</td></tr>
<tr><td>链式法则</td><td>$\dfrac{\partial J}{\partial\bm u}=\dfrac{\partial\bm y}{\partial\bm u}\,\dfrac{\partial J}{\partial\bm y}$（内层在左！）</td><td>$\frac{dJ}{du}=\frac{dJ}{dy}\frac{dy}{du}$</td></tr>
<tr><td>2.3.5</td><td>$J=\bm y^TB\bm y,\;\bm y=A\bm u\;\Rightarrow\;\dfrac{\partial J}{\partial\bm u}=2A^TB\bm y$</td><td></td></tr>
<tr><td>2.3.6</td><td>$\dfrac{\partial\,\mathrm{Tr}(KA)}{\partial K}=A^T$</td><td>（卡尔曼滤波要用）</td></tr>
<tr><td>2.3.7</td><td>$\dfrac{\partial\,\mathrm{Tr}(KAK^T)}{\partial K}=2KA$（A 对称）</td><td></td></tr></table>
<div class="warn" data-t="链式法则的陷阱">照搬标量写法 $\frac{\partial J}{\partial\bm y}\frac{\partial\bm y}{\partial\bm u}$，维度是 $(m\times1)(n\times m)$，根本乘不了！回到元素层面推导会发现，分母布局下正确的顺序是 $\frac{\partial\bm y}{\partial\bm u}\frac{\partial J}{\partial\bm y}$。<b>矩阵乘法不满足交换律，永远先检查维度。</b></div>
<details class="derive"><summary>例 2.3.2：二次型代价对控制量求导（后面推导 LQR 的核心套路）</summary><div class="body">
<div class="step">离散系统 $\bm x_{[k+1]}=A\bm x_{[k]}+B\bm u_{[k]}$，代价 $J=\bm x_{[k+1]}^TS\bm x_{[k+1]}$。</div>
<div class="step">链式法则：$\dfrac{\partial J}{\partial\bm u_{[k]}}=\dfrac{\partial(A\bm x_{[k]}+B\bm u_{[k]})}{\partial\bm u_{[k]}}\cdot 2S\bm x_{[k+1]}=2B^TS\bm x_{[k+1]}$。</div>
<div class="step">书中数值 $A=\begin{bmatrix}1&0&1\\0&1&0\\1&0&1\end{bmatrix}$，$B=\begin{bmatrix}0&1\\1&0\\0&1\end{bmatrix}$，$S=I$，得 $\begin{bmatrix}2(x_2+u_1)\\4(x_1+x_3+u_2)\end{bmatrix}$，与逐项展开求导完全一致。令其为 0 即得“一步最优”的控制量。</div>
<div class="step">后面常写成 $J=\frac12\bm x^TS\bm x$，系数 ½ 正好消掉求导出来的 2。</div>
</div></details>

<h2><span class="sec">2.4</span>应用：线性回归</h2>
<div class="story" data-t="数据：2021—2022 赛季五大联赛射手榜前三">15 名顶级前锋的身高与体重（孙兴慜 183 cm/75 kg、哈兰德 195/88、姆巴佩 178/73……）。能不能用体重 x 预测身高 z？</div>
<p>模型 $\hat z_i=y_1+y_2x_i$，用<b>最小二乘</b>：让误差平方和最小。写成矩阵：$\bm z$ 是身高列向量，$\bm x=[\bm 1,\;\text{体重}]$（$15\times2$），$\bm y=[y_1,y_2]^T$：</p>
<div class="eq" data-no="(2.4.6)–(2.4.8)" data-title="正规方程">
$$f(\bm y)=(\bm z-\bm x\bm y)^T(\bm z-\bm x\bm y)=\bm z^T\bm z-2\bm z^T\bm x\bm y+\bm y^T\bm x^T\bm x\bm y,\qquad \frac{\partial f}{\partial\bm y}=-2\bm x^T\bm z+2\bm x^T\bm x\bm y=0\;\Rightarrow\;\bm y^*=(\bm x^T\bm x)^{-1}\bm x^T\bm z$$
<div class="say">代入数据：$y_1=127.6$，$y_2=0.71$，即 身高 ≈ 127.6 + 0.71 × 体重。</div>
</div>
<p>当变量很多、$(\bm x^T\bm x)^{-1}$ 难求甚至不存在时，用<span class="kw">梯度下降</span>迭代：沿负梯度方向一步步走下山，$\bm y_{i+1}=\bm y_i-\bm\alpha\nabla$。这也是现代机器学习的基础。书中初值 $[120,1]$，学习率 $\bm\alpha=\mathrm{diag}(0.001,\,0.00001)$，迭代 13061 次停在 $[124.9,\,0.75]$。</p>
<div class="sim" data-id="regress"></div>
<div class="idea" data-t="为什么两个学习率差 100 倍？">截距（≈127）和斜率（≈0.7）的量级差了约 100 倍，代价函数的等高线变成极其狭长的椭圆。用同一个学习率，要么斜率方向步子太大而振荡，要么截距方向几乎不动。工程上常把变量<b>归一化</b>来解决——在上面的实验里勾选试试。</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>分母布局：标量对 n 维向量求导得 n×1 列向量（梯度）。</li>
<li>$\partial(\bm u^TA\bm u)/\partial\bm u=2A\bm u$（A 对称）；链式法则内层导数在左。</li>
<li>最小二乘解 $\bm y^*=(\bm x^T\bm x)^{-1}\bm x^T\bm z$；高维用梯度下降，学习率与归一化很关键。</li>
</ol></div>`,
    sims: {
      grad(el) {
        const ui = card(el, { title: "梯度：垂直于等高线、指向上坡", hint: "函数 $f(\\bm u)=\\bm u^TA\\bm u+\\bm b^T\\bm u$ 的等高线。<b>拖动黄点</b>，红色箭头是梯度 $\\frac{\\partial f}{\\partial\\bm u}=(A+A^T)\\bm u+\\bm b$。" });
        const pl = new Plot(ui.plotHost(), { height: 340, x: [-3, 3], y: [-3, 3], equal: true, xlabel: "u₁", ylabel: "u₂" });
        let u = [1.2, 1.5];
        const s12 = ui.slider({ label: "$a_{12}$（非对称部分）", min: -1.5, max: 1.5, step: 0.05, value: 0, onInput: draw });
        const s22 = ui.slider({ label: "$a_{22}$", min: 0.2, max: 3, step: 0.05, value: 1, onInput: draw });
        pl.drag({ down: (x, y) => { u = [x, y]; draw(); }, move: (x, y) => { u = [clamp(x, -2.9, 2.9), clamp(y, -2.9, 2.9)]; draw(); }, hover: () => true });
        function draw() {
          const A = [[1, s12.get()], [0, s22.get()]], b = [2, 0];
          const f = (x, y) => x * x * A[0][0] + x * y * (A[0][1] + A[1][0]) + y * y * A[1][1] + b[0] * x + b[1] * y;
          const g = M.vadd(M.mulv(M.add(A, M.T(A)), u), b);
          const H = M.add(A, M.T(A)), ust = M.vscale(M.solve(H, b), -1);
          pl.render = () => {
            pl.clear(); const c = pl.ctx;
            const lv = [-0.8, -0.5, 0, 0.5, 1.5, 3, 5, 8, 12, 17];
            // 粗略等高线：扫描网格，每点只算一次 f
            for (let x = -3; x <= 3; x += 0.03) for (let y = -3; y <= 3; y += 0.03) {
              const v = f(x, y);
              for (let i = 0; i < lv.length; i++) if (Math.abs(v - lv[i]) < 0.06 + 0.012 * Math.abs(lv[i])) { c.fillStyle = `hsla(${210 - i * 14},55%,55%,.6)`; c.fillRect(pl.fx(x), pl.fy(y), 1.5, 1.5); break; }
            }
            pl.dot(ust[0], ust[1], COL.c3, 6); pl.text(ust[0] + 0.1, ust[1] - 0.3, "极小点", COL.c3);
            const L = Math.hypot(...g), s = Math.min(1, 1.2 / Math.max(L, 1e-9));
            pl.arrow(u[0], u[1], u[0] + g[0] * s, u[1] + g[1] * s, COL.c4, 2.6, 10);
            pl.arrow(u[0], u[1], u[0] - g[0] * s * 0.6, u[1] - g[1] * s * 0.6, COL.c3, 1.6, 8);
            pl.dot(u[0], u[1], "#e8b400", 7, COL.ink);
          };
          pl.render();
          const ge = [2 * A[0][0] * u[0] + (A[0][1] + A[1][0]) * u[1] + b[0], (A[0][1] + A[1][0]) * u[0] + 2 * A[1][1] * u[1] + b[1]];
          ui.out(`$A=${mtex(A, 2)}$<br>逐元素求导：[${f2(ge[0])}, ${f2(ge[1])}]<br>公式 $(A+A^T)\\bm u+\\bm b$：[${f2(g[0])}, ${f2(g[1])}]<br><span class="sim-note">红箭头 = 梯度（上坡）；绿箭头 = 负梯度（梯度下降走的方向）。令梯度为 0 解出极小点 $\\bm u^*=-(A+A^T)^{-1}\\bm b$。</span>`);
        }
        draw();
      },
      regress(el) {
        const ui = card(el, { title: "线性回归：解析解 vs 梯度下降", hint: "左：15 名球员的散点、解析解（绿）与梯度下降当前直线（蓝）。右：代价函数等高线上的下降路径（纵轴为斜率 y₂，横轴为截距 y₁）。默认设置可复现书中的 13061 次迭代。" });
        const [a, b] = ui.row("", "");
        const ps = new Plot(a, { height: 300, x: [64, 91], y: [166, 198], xlabel: "体重 kg", ylabel: "身高 cm" });
        const pc = new Plot(b, { height: 300, x: [110, 150], y: [0.4, 1.1], xlabel: "y₁ 截距", ylabel: "y₂ 斜率" });
        const sa1 = ui.slider({ label: "学习率 α₁（截距）", min: 1e-5, max: 0.01, value: 0.001, log: true, fmt: v => v.toExponential(1), onInput: reset });
        const sa2 = ui.slider({ label: "学习率 α₂（斜率）", min: 1e-7, max: 1e-3, value: 1e-5, log: true, fmt: v => v.toExponential(1), onInput: reset });
        const cN = ui.check({ label: "体重归一化（减均值、除标准差）", value: false, onChange: reset });
        let speed = 200, run = false;
        ui.buttons([["▶ 运行 / 暂停", () => run = !run], ["⏭ 直接跑完", () => { run = false; while (!st.done && st.it < 200000) stepOnce(); draw(); }, "ghost"], ["↺", reset, "ghost"]]);
        const n = X.length, mx = X.reduce((s, v) => s + v, 0) / n, sd = Math.sqrt(X.reduce((s, v) => s + (v - mx) ** 2, 0) / n);
        // 解析解
        const XtX = [[n, X.reduce((s, v) => s + v, 0)], [X.reduce((s, v) => s + v, 0), X.reduce((s, v) => s + v * v, 0)]], XtZ = [Z.reduce((s, v) => s + v, 0), X.reduce((s, v, i) => s + v * Z[i], 0)];
        const ystar = M.solve(XtX, XtZ);
        let st;
        const cost = (y, xs) => xs.reduce((s, x, i) => s + (Z[i] - y[0] - y[1] * x) ** 2, 0);
        function reset() { const norm = cN.get(); st = { y: norm ? [120, 1] : [120, 1], it: 0, done: false, path: [], norm, f: null }; st.xs = norm ? X.map(v => (v - mx) / sd) : X; st.f = cost(st.y, st.xs); st.path.push(toOrig(st.y)); draw(); }
        function toOrig(y) { return st && st.norm ? [y[0] - y[1] * mx / sd, y[1] / sd] : y.slice(); }
        function stepOnce() {
          if (st.done) return;
          const xs = st.xs, r = xs.map((x, i) => st.y[0] + st.y[1] * x - Z[i]);
          const g = [r.reduce((s, v) => s + v, 0), r.reduce((s, v, i) => s + v * xs[i], 0)]; // xᵀ(xy − z)，与书中代码一致
          const a1 = sa1.get(), a2 = st.norm ? sa1.get() : sa2.get();
          st.y = [st.y[0] - a1 * g[0], st.y[1] - a2 * g[1]]; st.it++;
          const fn = cost(st.y, xs);
          if (!isFinite(fn) || Math.abs(fn - st.f) < 1e-4 || st.it > 200000) st.done = true;
          st.f = fn; if (st.it % 20 === 0 || st.done) st.path.push(toOrig(st.y));
        }
        // 等高线背景
        const grid = []; for (let i = 0; i <= 80; i++) for (let j = 0; j <= 80; j++) { const y1 = 110 + 40 * i / 80, y2 = 0.4 + 0.7 * j / 80; grid.push([y1, y2, Math.log10(cost([y1, y2], X))]); }
        function draw() {
          const yo = toOrig(st.y);
          ps.render = () => { ps.clear(); ps.fn(x => ystar[0] + ystar[1] * x, COL.c3, 2, [6, 4]); ps.fn(x => yo[0] + yo[1] * x, COL.c1, 2.6); X.forEach((x, i) => ps.dot(x, Z[i], COL.c2, 4.5, COL.plotbg)); ps.text(88, Z[8] + 1.2, NAMES[8], COL.muted, "right"); };
          pc.render = () => {
            pc.clear(); const c = pc.ctx; const lo = Math.min(...grid.map(g => g[2])), hi = lo + 2.2;
            grid.forEach(([a, b, v]) => { const t = clamp((v - lo) / (hi - lo), 0, 1); c.fillStyle = `hsla(${200 - 160 * t},60%,${55 + 15 * t}%,.35)`; c.fillRect(pc.fx(a) - 3, pc.fy(b) - 3, 7, 7); });
            pc.line(st.path.map(p => p[0]), st.path.map(p => p[1]), COL.c4, 2); pc.dot(ystar[0], ystar[1], COL.c3, 6); pc.dot(yo[0], yo[1], COL.c1, 5);
          };
          ps.render(); pc.render();
          ui.out(`解析解：y₁ = <b>${f2(ystar[0], 2)}</b>，y₂ = <b>${f2(ystar[1], 4)}</b><br>梯度下降：第 <b>${st.it}</b> 次迭代<br>y₁ = ${f2(yo[0], 2)}，y₂ = ${f2(yo[1], 4)}<br>代价 f = ${f2(cost(yo, X), 2)}${st.done ? '<br><span class="good">已满足停止条件 |Δf| &lt; 10⁻⁴</span>' : ""}`);
        }
        reset();
        CE.loop(ui.root, () => { if (!run || st.done) return; for (let k = 0; k < speed; k++) stepOnce(); draw(); });
      },
    },
  });
})();
