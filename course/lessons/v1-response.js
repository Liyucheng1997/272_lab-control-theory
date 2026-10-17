/* 卷1 第4–6章：一阶系统 / 二阶系统 / 稳定性 */
(function () {
  const { M, P, Plot, Scene, card, linspace, f2, COL, clamp, rk4, cfmt } = CE;

  // 二阶标准系统阶跃响应（数值）
  function step2(z, wn, T, n = 600, x0 = [0, 0]) {
    const dt = T / n, ts = [0], ys = [x0[0]]; let x = x0.slice();
    const f = s => [s[1], -2 * z * wn * s[1] - wn * wn * s[0] + wn * wn];
    for (let i = 1; i <= n; i++) { x = rk4(f, x, 0, dt); ts.push(i * dt); ys.push(x[0]); }
    return { t: ts, y: ys };
  }
  // 书中定义的指标：Tr = 第一次到达 1（达不到则 90%），Mp，Ts(2%)
  function metrics(t, y) {
    let Tr = null; for (let i = 1; i < y.length; i++) if (y[i] >= 1) { Tr = t[i]; break; }
    let Tr90 = null; for (let i = 1; i < y.length; i++) if (y[i] >= 0.9) { Tr90 = t[i]; break; }
    let pk = 0, Tp = null; y.forEach((v, i) => { if (v > pk) { pk = v; Tp = t[i]; } });
    let Ts = 0; for (let i = y.length - 1; i >= 0; i--) if (Math.abs(y[i] - 1) > 0.02) { Ts = t[Math.min(i + 1, t.length - 1)]; break; }
    return { Tr: Tr ?? Tr90, reached: Tr !== null, Mp: Math.max(0, pk - 1) * 100, Tp: pk > 1.0005 ? Tp : null, Ts };
  }

  /* ════════════════ 第4章 一阶系统 ════════════════ */
  Course.add({
    id: "v1-first", vol: 1, ord: 90, no: "4", chName: "第4章 一阶系统的时域响应", ref: "卷1 第4章",
    title: "一阶系统：密室里的案发时间",
    sub: "牛顿冷却定律是最典型的一阶系统。学会冲激响应、阶跃响应和时间常数，你就能当一回侦探。",
    goals: ["写出一阶系统的标准形式 $\\dot x+ax=au$ 和传递函数 $\\frac{a}{s+a}$", "掌握冲激响应 $ae^{-at}$、阶跃响应 $1-e^{-at}$ 与时间常数 τ", "用“冲激 + 阶跃”叠加推算出案发时间"],
    minutes: 25,
    html: String.raw`
<div class="story" data-t="引子：案发时间是几点？">你是一名侦探，来到一间密室凶案现场。下午 <b>2:30</b> 进门，用温度枪测得尸体温度 <b>26℃</b>。房间空调把室温恒定在 <b>20℃</b>。取证一小时后，<b>3:30</b> 再测，降到了 <b>25℃</b>。已知死者生前体温是标准的 <b>37℃</b>。请问：案发时间是几点？</div>
<p>理想化的密室里，可以用<b>牛顿冷却定律</b>：物体降温的速度与它和环境的温差成正比。</p>
<div class="eq" data-no="(4.1.1)–(4.1.2)" data-title="牛顿冷却定律 = 一阶系统">
$$\frac{\mathrm{d}T}{\mathrm{d}t} = -K\,(T - C)\quad\Longleftrightarrow\quad \dot x + Kx = Ku$$
<dl><dt>$x=T$</dt><dd>物体温度（输出）</dd><dt>$u=C$</dt><dd>环境温度（输入）</dd><dt>$K$</dt><dd>冷却系数：越大降温越快</dd></dl>
</div>
<p>保留初始条件做拉普拉斯变换：$sX - x(0) + KX = KU$，即 $(s+K)X = K\left(\frac{x(0)}{K} + U\right)$。于是可以把<b>初始温度当成一个额外的输入</b>：</p>
<div class="eq" data-no="(4.1.6)">
$$X(s) = \frac{K}{s+K}\Big(\underbrace{\tfrac{1}{K}x(0)}_{U_1:\ \text{冲激输入}} + \underbrace{U(s)}_{U_2:\ \text{阶跃输入}}\Big)$$
<div class="say">由叠加原理，输出 = 冲激响应 + 阶跃响应。所以要破案，先得弄懂这两种典型响应。</div>
</div>

<h2><span class="sec">4.2</span>一阶系统的两种典型响应</h2>
<p>典型一阶系统 $\dot x + ax = au$，传递函数 $G(s)=\dfrac{a}{s+a}$，唯一的极点 $s_p=-a$。两种典型输入：单位冲激 $\delta(t)$（$\mathcal L=1$）和单位阶跃（高度为 1 的“台阶”，$\mathcal L = 1/s$）。</p>
<div class="eq" data-no="(4.2.9) / (4.2.15)" data-title="冲激响应与阶跃响应">
$$\text{冲激：}\;X=\frac{a}{s+a}\Rightarrow x(t)=a e^{-at}\qquad\qquad \text{阶跃：}\;X=\frac{a}{s(s+a)}=\frac1s-\frac{1}{s+a}\Rightarrow x(t)=1-e^{-at}$$
</div>
<ul><li><b>对初始条件的响应就是冲激响应</b>：零输入、$x(0)=x_0$ 时 $x(t)=x_0e^{-at}$，形状与冲激响应一样。</li>
<li>收敛还是发散只看极点：$a>0$（极点在左边）收敛，$a<0$ 发散。</li></ul>
<div class="eq" data-title="两个性能指标">
$$\tau = \frac1a\;(\text{到达终值 }63\%),\qquad T_s = 4\tau = \frac4a\;(\text{到达终值 }98\%)$$
<dl><dt>$\tau$</dt><dd><b>时间常数</b>：反映响应速度；$a$ 越大（极点越靠左），τ 越小，反应越快</dd><dt>$T_s$</dt><dd><b>调节时间</b>：误差进入 2% 以内，工程上认为“到了”</dd></dl>
</div>
<div class="sim" data-id="first"></div>
<div class="tip" data-t="工程小实验：测出冷却系数">一杯 0℃ 的水放进 20℃ 的房间，测得 <b>5 小时</b>后水温 19.6℃（终值的 98%）。于是 $T_s=4\tau=5$，$\tau=1.25$ 小时，冷却系数 $K=1/\tau=0.8$。</div>
<p><b>相轨迹视角</b>：阶跃输入下 $\dot z = a - az$，这是一条与横轴交于 $z_f=1$ 的直线；$a>0$ 时两侧的箭头都指向 1——稳定平衡点。传递函数擅长给出<b>解析式</b>，相轨迹擅长<b>直观</b>判断，两者结论一致。</p>

<h2><span class="sec">4.3</span>案发时间揭秘</h2>
<p>冲激部分：$x_1 = 37e^{-Kt}$；阶跃部分（20 倍单位阶跃）：$x_2=20(1-e^{-Kt})$。叠加：</p>
<div class="eq" data-no="(4.3.1)–(4.3.3)">
$$x(t) = 20 + 17e^{-Kt},\qquad \begin{cases}26 = 20+17e^{-Kt_1}\\ 25 = 20+17e^{-K(t_1+1)}\end{cases}\Rightarrow\; K=\ln\tfrac{6}{5}\approx 0.182,\;\; t_1\approx 5.72\text{ h}$$
<div class="say">两式相除消去 $t_1$：$e^{-K}=5/6$。下午 2:30 时，死者已经去世 5.72 小时（约 5 小时 43 分），案发时间约为<b>上午 8:47</b>。</div>
</div>
<div class="sim" data-id="murder"></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>一阶系统：$\dot x+ax=au$，$G(s)=\frac{a}{s+a}$，极点 $-a$。</li>
<li>冲激响应 $ae^{-at}$（= 对初始条件的响应）；阶跃响应 $1-e^{-at}$。</li>
<li>时间常数 $\tau=1/a$（63%），调节时间 $T_s=4\tau$（98%）；极点越靠左越快。</li>
<li>初始条件可以视为 $t=0$ 的冲激输入，再用叠加原理分解问题。</li>
</ol></div>`,
    sims: {
      first(el) {
        const ui = card(el, { title: "一阶系统三联动：极点、时间响应、相轨迹", hint: "拖动 a，观察极点位置、响应快慢、相轨迹直线的斜率如何一起变化。" });
        const pr = new Plot(ui.plotHost(), { height: 230, x: [0, 8], y: [-0.2, 1.4], xlabel: "t", ylabel: "x(t)" });
        const [a1, a2] = ui.row("", "");
        const ps = new Plot(a1, { height: 200, x: [-4, 2], y: [-1.2, 1.2], xlabel: "σ", title: "极点 s = −a" });
        const pp = new Plot(a2, { height: 200, x: [-0.5, 2], y: [-2, 2], xlabel: "z", ylabel: "ż", title: "相轨迹" });
        let inp = "step";
        ui.seg({ label: "输入", options: [["step", "单位阶跃"], ["imp", "单位冲激"], ["init", "零输入, x₀=1"]], value: inp, onChange: v => { inp = v; draw(); } });
        const sa = ui.slider({ label: "参数 a", min: -0.6, max: 3, step: 0.05, value: 1, onInput: draw });
        function draw() {
          const a = sa.get(), ts = linspace(0, 8, 400);
          const f = inp === "step" ? t => 1 - Math.exp(-a * t) : inp === "imp" ? t => a * Math.exp(-a * t) : t => Math.exp(-a * t);
          const ys = ts.map(f); pr.autoY([ys.filter(v => Math.abs(v) < 4), [0, 1.1]], 0.08);
          pr.render = () => {
            pr.clear(); pr.line(ts, ys, COL.c1, 2.6);
            if (inp === "step" && a > 0.05) {
              const tau = 1 / a; pr.hline(1, COL.c3, [6, 4]);
              if (tau < 8) { pr.vline(tau, COL.c5, [3, 3]); pr.dot(tau, 1 - Math.exp(-1), COL.c5, 5); pr.text(tau + 0.08, 0.5, "τ: 63%", COL.c5); }
              if (4 * tau < 8) { pr.vline(4 * tau, COL.c2, [3, 3]); pr.dot(4 * tau, 0.9817, COL.c2, 5); pr.text(4 * tau + 0.08, 0.8, "Ts: 98%", COL.c2); }
            }
          };
          ps.render = () => { ps.clear(); ps.rect(0, -1.2, 2, 1.2, COL.c4 + "14"); ps.cross(-a, 0, COL.c4, 9, 3); };
          const zf = inp === "step" ? 1 : 0, rhs = z => inp === "step" ? a - a * z : -a * z;
          pp.render = () => {
            pp.clear(); pp.fn(rhs, COL.c1, 2.4);
            [-0.3, 0.4, 1.4, 1.8].forEach(z => { const v = rhs(z); if (Math.abs(v) > 1e-3) pp.arrow(z, 0, z + Math.sign(v) * 0.18, 0, v > 0 ? COL.c3 : COL.c4, 2.2, 7); });
            Math.abs(a) > 1e-3 && (a > 0 ? pp.dot(zf, 0, COL.ink, 6) : pp.dot(zf, 0, COL.plotbg, 6, COL.ink));
          };
          pr.render(); ps.render(); pp.render();
          ui.out(a > 0 ? `极点 $s_p=${f2(-a)}$（左半平面）<br><span class="good">收敛</span><br>时间常数 τ = <b>${f2(1 / a)}</b><br>调节时间 Ts = <b>${f2(4 / a)}</b>` : a < 0 ? `极点 $s_p=${f2(-a)}$（右半平面）<br><span class="bad">发散！</span>` : "a = 0：极点在原点，积分器");
        }
        draw();
      },
      murder(el) {
        const ui = card(el, { title: "密室破案计算器", hint: "调整两次测温读数和室温，计算器会解出冷却系数 K 并倒推案发时间。看看测量误差对结论的影响有多大。" });
        const p1 = new Plot(ui.plotHost(), { height: 230, x: [0, 9], y: [0, 40], xlabel: "死后经过的时间 / h", title: "叠加：冲激响应 x₁ + 阶跃响应 x₂" });
        const p2 = new Plot(ui.plotHost(), { height: 230, x: [6, 16.5], y: [18, 38], xlabel: "钟点", ylabel: "体温 ℃", title: "推算出的体温曲线" });
        const sC = ui.slider({ label: "室温 C (℃)", min: 10, max: 24, step: 0.5, value: 20, onInput: draw });
        const s1 = ui.slider({ label: "14:30 读数 (℃)", min: 22, max: 34, step: 0.1, value: 26, onInput: draw });
        const s2 = ui.slider({ label: "15:30 读数 (℃)", min: 21, max: 33.5, step: 0.1, value: 25, onInput: draw });
        const hm = h => { h = ((h % 24) + 24) % 24; const H = Math.floor(h), m = Math.round((h - H) * 60); return `${m === 60 ? H + 1 : H}:${String(m === 60 ? 0 : m).padStart(2, "0")}`; };
        function draw() {
          const C = sC.get(), T1 = s1.get(), T2 = s2.get(), x0 = 37;
          if (!(T2 < T1 && T2 > C && T1 > C)) { ui.out('<span class="bad">需要 室温 &lt; 第二次读数 &lt; 第一次读数</span>'); return; }
          const K = Math.log((T1 - C) / (T2 - C)), t1 = Math.log((x0 - C) / (T1 - C)) / K, death = 14.5 - t1;
          const ts = linspace(0, 9, 300);
          p1.render = () => {
            p1.clear();
            p1.line(ts, ts.map(t => x0 * Math.exp(-K * t)), COL.c2, 2, [6, 4]); p1.line(ts, ts.map(t => C * (1 - Math.exp(-K * t))), COL.c5, 2, [6, 4]);
            p1.line(ts, ts.map(t => C + (x0 - C) * Math.exp(-K * t)), COL.c1, 2.6);
            p1.legend([{ label: "x₁=37e^(−Kt)", color: COL.c2, dash: [6, 4] }, { label: `x₂=${C}(1−e^(−Kt))`, color: COL.c5, dash: [6, 4] }, { label: "x = x₁+x₂", color: COL.c1 }], "tr");
          };
          const hs = linspace(Math.max(0, death), 16.5, 300);
          p2.setRange([Math.min(6, Math.floor(death) - 0.5), 16.5], null);
          p2.render = () => {
            p2.clear(); p2.hline(C, COL.c3, [6, 4]);
            p2.line(hs, hs.map(h => C + (x0 - C) * Math.exp(-K * (h - death))), COL.c1, 2.6);
            p2.dot(death, 37, COL.c4, 6); p2.text(death + 0.15, 36.8, "案发 " + hm(death), COL.c4, "left", 12, true);
            p2.dot(14.5, T1, COL.c2, 5); p2.dot(15.5, T2, COL.c2, 5); p2.text(14.5, T1 + 1.2, "14:30 " + T1 + "℃", COL.c2, "center");
          };
          p1.render(); p2.render();
          ui.out(`$K=\\ln\\frac{${f2(T1 - C, 1)}}{${f2(T2 - C, 1)}}$ = <b>${f2(K, 3)}</b> /h<br>$t_1$ = <b>${f2(t1, 2)} h</b><br>案发时间 ≈ <span class="big">${hm(death)}</span>`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b>：把 15:30 的读数从 25.0 改成 24.8 或 25.2（只差 0.2℃），案发时间会差出一个多小时——两次读数太接近时，推算对测量误差非常敏感。真实法医推断会用多次测量来降低误差。</div>`);
      },
    },
  });

  /* ════════════════ 第5章 二阶系统 ════════════════ */
  Course.add({
    id: "v1-second", vol: 1, ord: 100, no: "5", chName: "第5章 二阶系统的时域响应", ref: "卷1 第5章",
    title: "二阶系统：阻尼比与固有频率",
    sub: "两个参数 ζ 和 ωₙ 决定了几乎所有运动控制系统的样子。学会读上升时间、超调量和稳定时间。",
    goals: ["掌握标准形式 $\\frac{\\omega_n^2}{s^2+2\\zeta\\omega_ns+\\omega_n^2}$ 及 ζ、ωₙ 的物理意义", "区分欠阻尼、临界阻尼、过阻尼、无阻尼", "会计算 $T_r$、$T_p$、$M_p$、$T_s$，理解它们与极点位置的几何关系"],
    minutes: 35,
    html: String.raw`
<p>自动驾驶、无人机、电机、机械臂……绝大多数和力学有关的控制问题都可以近似为二阶系统，很多高阶系统也能近似成二阶。</p>
<h2><span class="sec">5.1</span>一般形式</h2>
<p>从弹簧–质量–阻尼 $m\ddot x+b\dot x+kx=f$ 出发，定义两个参数：</p>
<div class="eq" data-no="(5.1.2)–(5.1.3)" data-title="二阶系统标准形式">
$$\omega_n=\sqrt{\frac km},\quad \zeta=\frac{b}{2\sqrt{km}}\qquad\Longrightarrow\qquad \ddot x+2\zeta\omega_n\dot x+\omega_n^2x=\omega_n^2u,\qquad G(s)=\frac{\omega_n^2}{s^2+2\zeta\omega_n s+\omega_n^2}$$
<dl><dt>$\omega_n$</dt><dd><b>固有频率</b>：弹簧硬、质量小 → ωₙ 大 → 反应快</dd><dt>$\zeta$</dt><dd><b>阻尼比</b>：阻尼相对于弹簧与质量的强弱，决定振不振、超调多少</dd><dt>$u=f/k$</dt><dd>“单位化”后的输入，使稳态输出等于 1</dd></dl>
</div>
<p>状态空间形式（状态 = 位置、速度）：$A=\begin{bmatrix}0&1\\-\omega_n^2&-2\zeta\omega_n\end{bmatrix}$，$B=\begin{bmatrix}0\\\omega_n^2\end{bmatrix}$。</p>

<h2><span class="sec">5.2</span>极点（特征值）与四种阻尼</h2>
<div class="eq" data-no="(5.2.3)">
$$\lambda_{1,2}=-\zeta\omega_n\pm\omega_n\sqrt{\zeta^2-1}$$
</div>
<table class="tbl"><tr><th>阻尼比</th><th>极点</th><th>相轨迹</th><th>名称</th></tr>
<tr><td>$\zeta>1$</td><td>两个不等负实数</td><td>稳定节点</td><td>过阻尼：不振荡，较慢</td></tr>
<tr><td>$\zeta=1$</td><td>两个相等负实数</td><td>稳定节点</td><td>临界阻尼：不振荡里最快的</td></tr>
<tr><td>$0<\zeta<1$</td><td>共轭复数，实部 $-\zeta\omega_n$</td><td>稳定焦点</td><td>欠阻尼：振荡衰减（最常见的弹簧）</td></tr>
<tr><td>$\zeta=0$</td><td>纯虚数 $\pm j\omega_n$</td><td>中心点</td><td>无阻尼：永远振荡</td></tr></table>

<h2><span class="sec">5.3</span>单位阶跃响应</h2>
<p>欠阻尼时用部分分式（三个极点：输入的 0 和系统的一对复极点）逐步推导，可得：</p>
<div class="eq" data-no="(5.3.12)" data-title="欠阻尼二阶系统的单位阶跃响应">
$$x(t)=1-\underbrace{e^{-\zeta\omega_n t}\sqrt{\tfrac{1}{1-\zeta^2}}}_{\text{衰减包络}}\;\underbrace{\sin(\omega_d t+\varphi)}_{\text{振荡}},\qquad \omega_d=\omega_n\sqrt{1-\zeta^2},\quad \varphi=\arccos\zeta$$
<dl><dt>$1$</dt><dd>来自输入的极点 0：稳态值</dd><dt>$\zeta\omega_n$</dt><dd>极点<b>实部</b>的大小：决定衰减（收敛）速度</dd><dt>$\omega_d$</dt><dd><b>阻尼固有频率</b> = 极点<b>虚部</b>：决定振荡快慢</dd></dl>
</div>
<details class="derive"><summary>部分分式推导要点</summary><div class="body">
<div class="step">$X(s)=\dfrac{\omega_n^2}{s(s^2+2\zeta\omega_ns+\omega_n^2)}=\dfrac{A}{s}+\dfrac{B}{s-s_{p2}}+\dfrac{C}{s-s_{p3}}$，$s_{p2,3}=-\zeta\omega_n\pm j\omega_d$。</div>
<div class="step">令 $s=0$ 得 $A=1$；令 $s=s_{p2}$ 得 $B=-\frac12\left(1-\frac{\zeta}{\sqrt{1-\zeta^2}}j\right)$，$C$ 是 $B$ 的共轭。</div>
<div class="step">逆变换并提出 $e^{-\zeta\omega_nt}$，用欧拉公式 $\frac12(e^{j\omega_dt}+e^{-j\omega_dt})=\cos\omega_dt$ 等，得 $x=1-e^{-\zeta\omega_nt}\left[\cos\omega_dt+\frac{\zeta}{\sqrt{1-\zeta^2}}\sin\omega_dt\right]$。</div>
<div class="step">辅助角公式合并为 $\sqrt{\frac{1}{1-\zeta^2}}\sin(\omega_dt+\varphi)$，其中 $\cos\varphi=\zeta$。</div>
<div class="step">其他情形：$\zeta=0$：$x=1-\cos\omega_nt$；$\zeta=1$：$x=1-e^{-\omega_nt}(1+\omega_nt)$。</div>
</div></details>
<div class="sim" data-id="second"></div>

<h2><span class="sec">5.4</span>性能指标</h2>
<p>书中对三个指标的定义（注意本书的上升时间是<b>第一次到达稳态值 100%</b> 的时间；过阻尼到不了 100% 时取 90%）：</p>
<div class="eq" data-no="(5.4.3)–(5.4.9)" data-title="欠阻尼二阶系统性能指标">
$$T_r=\frac{\pi-\varphi}{\omega_d},\qquad T_p=\frac{\pi}{\omega_d},\qquad M_p=e^{-\frac{\zeta\pi}{\sqrt{1-\zeta^2}}}\times100\%,\qquad T_s\approx\frac{4}{\zeta\omega_n}$$
<dl><dt>$T_r$ 上升时间</dt><dd>反应速度；ωₙ 越大越快，ζ 越大越慢</dd><dt>$T_p$ 峰值时间</dt><dd>到达最高点的时间，由极点虚部 $\omega_d$ 决定</dd><dt>$M_p$ 最大超调</dt><dd><b>只和 ζ 有关</b>：“矫枉过正”的程度</dd><dt>$T_s$ 稳定时间</dt><dd>进入 ±2% 误差带的时间，由极点实部 $\zeta\omega_n$ 决定</dd></dl>
</div>
<details class="derive"><summary>指标是怎么推出来的</summary><div class="body">
<div class="step"><b>上升时间</b>：令 $x(T_r)=1$，得 $\sin(\omega_dT_r+\varphi)=0$，第一次取 $\omega_dT_r+\varphi=\pi$。</div>
<div class="step"><b>峰值时间</b>：令 $\dot x(T_p)=0$，化简得 $\tan(\omega_dT_p+\varphi)=\tan\varphi$，即 $\omega_dT_p=\pi$。</div>
<div class="step"><b>超调</b>：代回 $x(T_p)=1+e^{-\zeta\omega_n\pi/\omega_d}\frac{\sin\varphi}{\sqrt{1-\zeta^2}}=1+e^{-\zeta\pi/\sqrt{1-\zeta^2}}$。</div>
<div class="step"><b>稳定时间</b>：取振荡项最大值 1，令包络 $e^{-\zeta\omega_nT_s}\approx0.02$，$-\ln0.02\approx3.9\approx4$。</div>
</div></details>
<div class="idea" data-t="s 平面上的几何">极点 $-\zeta\omega_n\pm j\omega_d$ 到原点的距离是 $\omega_n$，与负实轴的夹角 φ 满足 $\cos\varphi=\zeta$。于是：<b>同一条射线上超调相同</b>；<b>同一条竖线上稳定时间相同</b>；<b>同一条横线上峰值时间相同</b>。设计控制器，就是在 s 平面上为极点划一块“合格区域”。</div>
<div class="sim" data-id="splane2"></div>

<h3>案例：你是无人机定高项目的经理</h3>
<div class="story" data-t="项目经理的选择">工程团队交来三种高度控制算法。实际工程中，大多数运动控制闭环都近似为一个二阶系统——相当于把无人机挂在一根“看不见的弹簧阻尼器”上，设计控制器就是设计这根弹簧的 ζ 和 ωₙ。轨迹 1 上升最快但超调最大；轨迹 2 居中；轨迹 3 最慢但没有超调、最早稳定。你选哪个？</div>
<div class="sim" data-id="drone3"></div>
<p>书中的结论很有启发：如果是<b>编队灯光秀</b>，无人机彼此靠得很近，超调意味着撞机，应选轨迹 3；如果是在能见度很差的区域<b>勘探避障</b>，反应速度最重要，轨迹 1 更合适。就像开车：平时变道慢慢并线（轨迹 3），遇到险情先猛打方向脱离危险（轨迹 1）。<b>没有绝对的对错，只有适合与否</b>——这就是工程中的思辨精神。</p>

<div class="summary"><h3>本课要点</h3><ol>
<li>标准二阶系统由 ωₙ（快慢）和 ζ（振荡程度）完全决定；极点 $-\zeta\omega_n\pm j\omega_n\sqrt{1-\zeta^2}$。</li>
<li>ζ=0 无阻尼、0&lt;ζ&lt;1 欠阻尼、ζ=1 临界阻尼、ζ&gt;1 过阻尼。</li>
<li>$T_r=\frac{\pi-\varphi}{\omega_d}$，$T_p=\frac{\pi}{\omega_d}$，$M_p=e^{-\zeta\pi/\sqrt{1-\zeta^2}}$，$T_s\approx\frac{4}{\zeta\omega_n}$。</li>
<li>指标之间有取舍，评价系统要结合应用场景。</li>
</ol></div>`,
    sims: {
      second(el) {
        const ui = card(el, { title: "二阶系统实验台：ζ 与 ωₙ", hint: "左：阶跃响应（自动标注上升时间 Tr、超调 Mp、稳定时间 Ts，绿色带为 ±2% 误差带）；右：极点在 s 平面上的位置（灰色圆半径 = ωₙ，紫色射线满足 cos φ = ζ）。" });
        const [a, b] = ui.row("", "");
        const pr = new Plot(a, { height: 330, x: [0, 12], y: [0, 2.05], xlabel: "t", ylabel: "x(t)" });
        const ps = new Plot(b, { height: 330, x: [-5, 1], y: [-3.2, 3.2], equal: true, xlabel: "σ", title: "极点" });
        const sz = ui.slider({ label: "阻尼比 ζ", min: 0, max: 2, step: 0.01, value: 0.3, onInput: draw });
        const sw = ui.slider({ label: "固有频率 ωₙ", min: 0.5, max: 3, step: 0.05, value: 1.5, onInput: draw });
        const cAll = ui.check({ label: "叠加书中图 5.3.4 的五条曲线 (ζ=0,0.2,0.5,1,2)", value: false, onChange: draw });
        function draw() {
          const z = sz.get(), wn = sw.get(), r = step2(z, wn, 12, 800), m = metrics(r.t, r.y);
          pr.render = () => {
            pr.clear(); pr.hline(1, COL.c3, [6, 4]); pr.rect(0, 0.98, 12, 1.02, COL.c3 + "20");
            if (cAll.get()) [0, 0.2, 0.5, 1, 2].forEach((zz, i) => { const q = step2(zz, wn, 12, 600); pr.line(q.t, q.y, `hsl(${i * 60 + 190},50%,62%)`, 1.3); });
            pr.line(r.t, r.y, COL.c1, 2.8);
            if (m.Tr) { pr.vline(m.Tr, COL.c5, [3, 3]); pr.text(m.Tr + 0.05, 0.12, "Tr", COL.c5, "left", 12, true); }
            if (m.Tp) { pr.dot(m.Tp, 1 + m.Mp / 100, COL.c4, 5); pr.text(m.Tp + 0.1, 1 + m.Mp / 100 + 0.05, "Mp=" + f2(m.Mp, 1) + "%", COL.c4, "left", 12, true); }
            if (m.Ts < 12 && z > 0) { pr.vline(m.Ts, COL.c2, [3, 3]); pr.text(m.Ts + 0.05, 0.3, "Ts", COL.c2, "left", 12, true); }
          };
          const e = M.eig2([[0, 1], [-wn * wn, -2 * z * wn]]);
          ps.render = () => {
            ps.clear(); const th = linspace(Math.PI / 2, 3 * Math.PI / 2, 80);
            ps.line(th.map(t => wn * Math.cos(t)), th.map(t => wn * Math.sin(t)), COL.axis, 1.3, [4, 4]);
            if (z < 1) { const ang = Math.acos(z); ps.line([0, -6 * Math.cos(ang)], [0, 6 * Math.sin(ang)], COL.c5, 1.2, [3, 3]); ps.line([0, -6 * Math.cos(ang)], [0, -6 * Math.sin(ang)], COL.c5, 1.2, [3, 3]); ps.text(-wn * Math.cos(ang) * 1.15, wn * Math.sin(ang) * 1.15 + 0.15, "ζ=cos φ", COL.c5, "center"); }
            e.vals.forEach(v => ps.cross(v[0], v[1], COL.c4, 8, 3));
          };
          pr.render(); ps.render();
          const typ = z === 0 ? "无阻尼" : z < 1 ? "欠阻尼" : z === 1 ? "临界阻尼" : "过阻尼";
          const wd = wn * Math.sqrt(Math.max(0, 1 - z * z)), phi = Math.acos(Math.min(1, z));
          const fm = z > 0 && z < 1 ? `<br><b>公式值</b>：Tr=${f2((Math.PI - phi) / wd)}，Tp=${f2(Math.PI / wd)}，Mp=${f2(100 * Math.exp(-z * Math.PI / Math.sqrt(1 - z * z)), 1)}%，Ts≈${f2(4 / (z * wn))}` : "";
          ui.out(`<b>${typ}</b>　极点 ${cfmt(e.vals[0])}, ${cfmt(e.vals[1])}<br><b>仿真实测</b>：Tr=${f2(m.Tr)}${m.reached ? "" : "(90%)"}，Tp=${m.Tp ? f2(m.Tp) : "—"}，Mp=${f2(m.Mp, 1)}%，Ts=${z > 0 ? f2(m.Ts) : "∞"}${fm}`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>ζ 从 0 慢慢增大到 2：极点沿灰色圆弧从虚轴滑向实轴，在 ζ=1 处“碰头”，之后沿实轴分开。</li><li>固定 ζ=0.3，改变 ωₙ：曲线被横向缩放，超调量 Mp 完全不变。</li><li>对比“公式值”与“实测”：Ts 公式只是近似（忽略了 $\\sqrt{1-\\zeta^2}$ 因子）。</li></ol></div>`);
      },
      splane2(el) {
        const ui = card(el, { title: "在 s 平面上拖动极点", hint: "拖动上半平面的极点（下半平面的共轭极点自动跟随）。沿射线移动：超调不变；竖直移动：只改振荡快慢；水平移动：只改稳定时间。" });
        const [a, b] = ui.row("", "");
        const ps = new Plot(a, { height: 300, x: [-4, 0.5], y: [-3.5, 3.5], equal: true, xlabel: "σ", title: "s 平面" });
        const pr = new Plot(b, { height: 300, x: [0, 8], y: [0, 2], xlabel: "t", title: "阶跃响应" });
        let p = [-0.8, 2], prev = null;
        ps.drag({ down: (x, y) => ps.pxDist(p[0], Math.abs(p[1]), x, Math.abs(y)) < 16 ? undefined : false, move: (x, y) => { prev = p.slice(); p = [clamp(x, -3.9, -0.05), clamp(Math.abs(y), 0.05, 3.4)]; draw(); }, hover: (x, y) => ps.pxDist(p[0], Math.abs(p[1]), x, Math.abs(y)) < 16 });
        const ghost = [];
        ui.buttons([["📌 保留当前曲线做对比", () => { ghost.push(p.slice()); if (ghost.length > 4) ghost.shift(); draw(); }, "ghost"], ["清除", () => { ghost.length = 0; draw(); }, "ghost"]]);
        function resp(q) { const wn = Math.hypot(q[0], q[1]), z = -q[0] / wn; return step2(z, wn, 8, 500); }
        function draw() {
          const wn = Math.hypot(...p), z = -p[0] / wn, r = resp(p), m = metrics(r.t, r.y);
          ps.render = () => {
            ps.clear();
            const ang = Math.atan2(p[1], -p[0]); ps.line([0, -6 * Math.cos(ang)], [0, 6 * Math.sin(ang)], COL.c5, 1.3, [4, 3]); ps.line([0, -6 * Math.cos(ang)], [0, -6 * Math.sin(ang)], COL.c5, 1.3, [4, 3]);
            ps.vline(p[0], COL.c2, [4, 3]); ps.hline(p[1], COL.c3, [4, 3]); ps.hline(-p[1], COL.c3, [4, 3]);
            ghost.forEach(g => { ps.cross(g[0], g[1], COL.muted, 6, 2); ps.cross(g[0], -g[1], COL.muted, 6, 2); });
            ps.cross(p[0], p[1], COL.c4, 9, 3); ps.cross(p[0], -p[1], COL.c4, 9, 3);
            ps.text(-3.8, 3.1, "紫：等超调射线", COL.c5); ps.text(-3.8, 2.7, "橙：等 Ts 竖线", COL.c2); ps.text(-3.8, 2.3, "绿：等 Tp 横线", COL.c3);
          };
          pr.render = () => { pr.clear(); pr.hline(1, COL.c3, [6, 4]); ghost.forEach(g => { const q = resp(g); pr.line(q.t, q.y, COL.muted, 1.4); }); pr.line(r.t, r.y, COL.c1, 2.8); };
          ps.render(); pr.render();
          ui.out(`极点 ${f2(p[0])} ± ${f2(p[1])}j<br>ωₙ = ${f2(wn)}，ζ = ${f2(z)}<br>Mp = <b>${f2(m.Mp, 1)}%</b>　Tp = <b>${m.Tp ? f2(m.Tp) : "—"}</b>　Ts = <b>${f2(m.Ts)}</b>`);
        }
        draw();
      },
      drone3(el) {
        const ui = card(el, { title: "项目经理挑战：给三种算法打分", hint: "三架无人机同时从 0 升到 10 m 的参考高度。下方雷达图按书中规则打分（每项 1–3 分）。你还可以调第 4 架“自定义”无人机。" });
        const sc = new Scene(ui.plotHost(), { height: 220 });
        const [a, b] = ui.row("", "");
        const pr = new Plot(a, { height: 260, x: [0, 6], y: [0, 15], xlabel: "t / s", ylabel: "高度 m" });
        const rad = new Scene(b, { height: 260 });
        const D = [{ n: "轨迹1", z: 0.3, w: 3, c: "c4" }, { n: "轨迹2", z: 0.6, w: 2, c: "c2" }, { n: "轨迹3", z: 1.0, w: 2.2, c: "c3" }, { n: "自定义", z: 0.45, w: 2.5, c: "c5" }];
        const sz = ui.slider({ label: "自定义 ζ", min: 0.1, max: 1.6, step: 0.01, value: 0.45, onInput: v => { D[3].z = v; recompute(); } });
        const sw = ui.slider({ label: "自定义 ωₙ", min: 0.8, max: 4, step: 0.05, value: 2.5, onInput: v => { D[3].w = v; recompute(); } });
        let t0 = 0;
        ui.buttons([["▶ 重新起飞", () => { t0 = 0; }]]);
        function recompute() { D.forEach(d => { d.r = step2(d.z, d.w, 6, 600); d.m = metrics(d.r.t, d.r.y); }); score(); }
        function score() {
          const rank = (key, better) => { const arr = D.slice(0, 3).map((d, i) => [d.m[key], i]).sort((x, y) => better(x[0], y[0])); arr.forEach(([, i], k) => D[i]["s_" + key] = 3 - k); };
          rank("Tr", (x, y) => x - y); rank("Mp", (x, y) => x - y); rank("Ts", (x, y) => x - y);
          ui.out(`<table class="tbl" style="margin:0;font-size:.8rem"><tr><th></th><th>Tr</th><th>Mp</th><th>Ts</th></tr>${D.map(d => `<tr><td style="color:var(--${d.c})"><b>${d.n}</b></td><td>${f2(d.m.Tr)}</td><td>${f2(d.m.Mp, 0)}%</td><td>${f2(d.m.Ts)}</td></tr>`).join("")}</table>`);
        }
        recompute();
        CE.loop(ui.root, dt => {
          t0 = Math.min(6, t0 + dt);
          const k = Math.round(t0 / 6 * 600);
          const c = sc.ctx, W = sc.W, H = sc.H; sc.clear();
          const y2p = h => H - 20 - h / 15 * (H - 40);
          c.strokeStyle = COL.c3; c.setLineDash([6, 4]); c.beginPath(); c.moveTo(0, y2p(10)); c.lineTo(W, y2p(10)); c.stroke(); c.setLineDash([]);
          c.fillStyle = COL.c3; c.font = "12px 'Microsoft YaHei'"; c.fillText("参考高度 10 m", 6, y2p(10) - 5);
          D.forEach((d, i) => {
            const x = W * (i + 0.6) / 4.3, h = 10 * d.r.y[k], y = y2p(h), col = COL[d.c];
            c.strokeStyle = col; c.lineWidth = 3; c.beginPath(); c.moveTo(x - 22, y); c.lineTo(x + 22, y); c.stroke();
            [-22, 22].forEach(dx => { c.beginPath(); c.ellipse(x + dx, y - 3, 10, 3, 0, 0, 7); c.stroke(); });
            c.fillStyle = col; c.fillRect(x - 6, y - 4, 12, 9);
            c.fillStyle = COL.ink2; c.textAlign = "center"; c.fillText(d.n, x, H - 4); c.textAlign = "left";
          });
          pr.render = () => { pr.clear(); pr.hline(10, COL.c3, [6, 4]); D.forEach(d => pr.line(d.r.t.slice(0, k + 1), d.r.y.slice(0, k + 1).map(v => 10 * v), COL[d.c], 2.4)); pr.legend(D.map(d => ({ label: d.n, color: COL[d.c] }))); };
          pr.render();
          // 雷达图
          const r = rad, cc = r.ctx; r.clear(); const cx = r.W / 2, cy = r.H / 2 + 12, R = Math.min(r.W, r.H) * 0.36;
          const axes = [["Mp", -Math.PI / 2], ["Ts", Math.PI / 6], ["Tr", 5 * Math.PI / 6]];
          cc.strokeStyle = COL.grid; cc.lineWidth = 1;
          [1, 2, 3].forEach(s => { cc.beginPath(); axes.forEach(([, an], i) => { const X = cx + R * s / 3 * Math.cos(an), Y = cy + R * s / 3 * Math.sin(an); i ? cc.lineTo(X, Y) : cc.moveTo(X, Y); }); cc.closePath(); cc.stroke(); });
          cc.fillStyle = COL.ink2; cc.font = "600 12px 'Microsoft YaHei'"; cc.textAlign = "center";
          axes.forEach(([n, an]) => { cc.beginPath(); cc.moveTo(cx, cy); cc.lineTo(cx + R * Math.cos(an), cy + R * Math.sin(an)); cc.stroke(); cc.fillText(n, cx + (R + 16) * Math.cos(an), cy + (R + 16) * Math.sin(an) + 4); });
          D.slice(0, 3).forEach(d => {
            cc.strokeStyle = COL[d.c]; cc.fillStyle = COL[d.c] + "22"; cc.lineWidth = 2.2; cc.beginPath();
            axes.forEach(([n, an], i) => { const s = d["s_" + n], X = cx + R * s / 3 * Math.cos(an), Y = cy + R * s / 3 * Math.sin(an); i ? cc.lineTo(X, Y) : cc.moveTo(X, Y); });
            cc.closePath(); cc.fill(); cc.stroke();
          });
          cc.fillStyle = COL.muted; cc.font = "11px 'Microsoft YaHei'"; cc.fillText("性能雷达图（越大越好，轨迹1–3）", cx, 16);
        });
      },
    },
  });

  /* ════════════════ 第6章 稳定性 ════════════════ */
  Course.add({
    id: "v1-stab", vol: 1, ord: 110, no: "6", chName: "第6章 稳定性", ref: "卷1 第6章",
    title: "稳定性：控制系统的第一要求",
    sub: "系统不稳定，其他一切性能都无从谈起。弄清楚李雅普诺夫稳定、渐近稳定、BIBO 稳定，以及怎么用极点判断。",
    goals: ["区分临界稳定、渐近稳定、不稳定", "理解李雅普诺夫意义下的稳定（ε–δ 定义）和全局/局部稳定", "用极点（或 A 的特征值）判断稳定性，理解 BIBO 稳定与共振"],
    minutes: 25,
    html: String.raw`
<h2><span class="sec">6.1</span>直观理解：轨道上的小球</h2>
<p>一条起伏的轨道上有三个位置：<b>A</b> 是光滑的谷底，<b>B</b> 是峰顶，<b>C</b> 是有摩擦的谷底。小球静止放在这三处都不会动——它们都是<span class="kw">平衡点</span>。但轻轻推一下：</p>
<ul><li>A：小球来回摆动、幅度不增不减——<b>临界稳定</b>（有界，但不回到原位）。</li>
<li>B：一去不复返——<b>不稳定</b>。</li>
<li>C：摩擦消耗能量，摆动越来越小，最终停回 C——<b>（渐近）稳定</b>。</li></ul>
<div class="sim" data-id="ball"></div>

<h3>严格定义（李雅普诺夫，1892）</h3>
<p>对无输入系统 $\dot{\bm z}=\bm f(\bm z)$，平衡点 $\bm z_f$ 满足 $\bm f(\bm z_f)=0$（不妨平移到原点）。</p>
<div class="eq" data-no="(6.1.3)" data-title="李雅普诺夫意义下的稳定">
$$\forall\varepsilon>0,\;\exists\,\delta>0:\quad \|\bm z(t_0)\|<\delta\;\Rightarrow\;\|\bm z(t)\|<\varepsilon,\;\forall t\ge t_0$$
<div class="say">无论你划多小的“警戒圈” ε，总能找到一个“起点圈” δ：只要从 δ 内出发，就永远不会跑出 ε。这是一个“有界”的概念。</div>
</div>
<div class="eq" data-no="(6.1.4)" data-title="渐近稳定">
$$\text{李雅普诺夫稳定，且}\quad \|\bm z(t_0)\|<\delta\;\Rightarrow\;\lim_{t\to\infty}\|\bm z(t)\|=0$$
<div class="say">不但不跑远，最终还要回到平衡点——比李雅普诺夫稳定更严格。</div>
</div>
<p>δ 可以理解为平衡点能承受的最大扰动（收敛域）。如果 δ 可以任意大，叫<b>全局稳定</b>；否则是<b>局部稳定</b>——比如把 A 处的小球推得太猛、越过了 B，它就再也回不到 A 附近了。</p>
<p>注意稳定性是<b>针对平衡点</b>说的。一个单摆有两个平衡点：正下方渐近稳定，正上方不稳定。要让正上方也稳定，就得加控制——比如在两边放吹风机，按偏离角度调风力，这就构成了一个反馈控制系统。</p>

<h2><span class="sec">6.2</span>用传递函数判断稳定性</h2>
<p>看单位冲激响应（就是 $G(s)$ 本身）。把 $G(s)$ 分解成部分分式，逆变换后：</p>
<div class="eq" data-no="(6.2.3)">
$$x(t)=\sum_i A_i e^{s_{pi}t}+\sum_k B_k e^{\sigma_k t}\sin(\omega_k t+\varphi_k)$$
<div class="say">每个实极点贡献一个指数，每对复极点贡献一个“指数 × 正弦”。只要有一个极点的实部大于 0，x(t) 就会发散。</div>
</div>
<div class="sim" data-id="poles"></div>
<div class="eq" data-title="稳定性判据（经典控制中“稳定”即渐近稳定）">
$$\text{系统稳定}\iff \text{(闭环)传递函数的所有极点都在左半平面：}\;\mathrm{Re}(s_p)<0$$
</div>
<p><b>BIBO 稳定</b>（有界输入有界输出）：满足上式时，任何有界输入都只产生有界输出。极点在虚轴上的“临界稳定”系统<b>不是</b> BIBO 稳定的：用和它固有频率相同的正弦去推，会发生<span class="kw">共振</span>，输出越来越大。</p>
<div class="sim" data-id="bibo"></div>
<p>对控制系统而言，<b>控制器设计的首要目标</b>就是让闭环传递函数 $\frac{CG}{1+CG}$ 的极点全部位于左半平面。在稳定的前提下，才谈得上调速度、减误差。</p>

<h2><span class="sec">6.3</span>用状态空间方程判断稳定性</h2>
<p>对 $\dot{\bm z}=A\bm z$，用特征向量解耦后每个分量是 $C_ie^{\lambda_it}$，所以：</p>
<div class="eq" data-title="状态空间稳定性判据">
$$\text{A 的特征值实部都}\le 0\Rightarrow\text{李雅普诺夫稳定};\qquad \text{都}<0\Rightarrow\text{渐近稳定}$$
</div>
<p>引入状态反馈 $\bm u=-K\bm z$，闭环变成 $\dot{\bm z}=(A-BK)\bm z$。只要选对 K，让 $A-BK$ 的特征值全在左半平面，系统就稳定了——这是第 10 章“极点配置”的核心思想。</p>

<div class="summary"><h3>本课要点</h3><ol>
<li>李雅普诺夫稳定 = 始终有界；渐近稳定 = 有界且收敛到平衡点。</li>
<li>传递函数：所有极点实部 &lt; 0 ⇔ 稳定（且 BIBO 稳定）；虚轴上的极点 = 临界稳定，可能共振。</li>
<li>状态空间：A 的特征值实部决定稳定性；状态反馈把 A 变成 A − BK。</li>
</ol></div>`,
    sims: {
      ball(el) {
        const ui = card(el, { title: "轨道上的小球", hint: "<b>拖动小球</b>到任意位置后松手。左半边光滑（A 附近），右半边有摩擦（C 附近）。下图记录小球位置 x(t)。" });
        const sc = new Scene(ui.plotHost(), { height: 220 });
        const pt = new Plot(ui.plotHost(), { height: 180, x: [0, 12], y: [-4, 4], xlabel: "t", ylabel: "x" });
        const g = 6, yf = x => Math.cos(Math.PI * x / 2) * 0.9, dyf = x => -0.9 * Math.PI / 2 * Math.sin(Math.PI * x / 2);
        const fr = x => x > 0.3 ? 0.9 : 0;
        let s = [-1.4, 0], t = 0, H = [], dragging = false;
        const X2P = x => sc.W / 2 + x / 4.2 * (sc.W / 2 - 10), Y2P = y => 120 - y * 70;
        sc.cv.style.touchAction = "none"; sc.cv.style.cursor = "grab";
        const toX = e => { const r = sc.cv.getBoundingClientRect(); return clamp((e.clientX - r.left - sc.W / 2) / (sc.W / 2 - 10) * 4.2, -3.9, 3.9); };
        sc.cv.addEventListener("pointerdown", e => { dragging = true; sc.cv.setPointerCapture(e.pointerId); s = [toX(e), 0]; });
        sc.cv.addEventListener("pointermove", e => { if (dragging) s = [toX(e), 0]; });
        sc.cv.addEventListener("pointerup", () => { dragging = false; t = 0; H = []; });
        CE.loop(ui.root, dt => {
          if (!dragging) {
            for (let i = 0; i < 10; i++) { const h = dt / 10, d = dyf(s[0]); const acc = -g * d / (1 + d * d) - fr(s[0]) * s[1]; s = [s[0] + h * s[1], s[1] + h * acc]; if (s[0] < -4) { s[0] = -4; s[1] = Math.abs(s[1]) * 0.5; } if (s[0] > 4) { s[0] = 4; s[1] = -Math.abs(s[1]) * 0.5; } }
            t += dt; if (t < 12) H.push([t, s[0]]);
          }
          const c = sc.ctx; sc.clear();
          c.lineWidth = 4; c.beginPath(); for (let x = -4.2; x <= 0.3; x += 0.05) { x === -4.2 ? c.moveTo(X2P(x), Y2P(yf(x))) : c.lineTo(X2P(x), Y2P(yf(x))); } c.strokeStyle = COL.c1; c.stroke();
          c.beginPath(); for (let x = 0.3; x <= 4.2; x += 0.05) { x === 0.3 ? c.moveTo(X2P(x), Y2P(yf(x))) : c.lineTo(X2P(x), Y2P(yf(x))); } c.strokeStyle = COL.c2; c.stroke();
          for (let x = 0.4; x < 4.2; x += 0.18) { c.strokeStyle = COL.c2 + "88"; c.lineWidth = 1; c.beginPath(); c.moveTo(X2P(x), Y2P(yf(x)) + 3); c.lineTo(X2P(x) - 5, Y2P(yf(x)) + 10); c.stroke(); }
          [["A", -2, "光滑"], ["B", 0, ""], ["C", 2, "有摩擦"]].forEach(([n, x, d]) => { c.fillStyle = COL.ink2; c.font = "600 14px 'Microsoft YaHei'"; c.textAlign = "center"; c.fillText(n, X2P(x), Y2P(yf(x)) + (n === "B" ? -16 : 28)); if (d) { c.font = "11px 'Microsoft YaHei'"; c.fillText(d, X2P(x), Y2P(yf(x)) + 42); } });
          const bx = X2P(s[0]), by = Y2P(yf(s[0])) - 11; c.fillStyle = COL.c4; c.beginPath(); c.arc(bx, by, 10, 0, 7); c.fill();
          pt.render = () => { pt.clear(); [-2, 0, 2].forEach(v => pt.hline(v, COL.axis, [3, 3])); pt.line(H.map(h => h[0]), H.map(h => h[1]), COL.c4, 2.2); pt.text(11.9, 2.15, "C", COL.muted, "right"); pt.text(11.9, 0.15, "B", COL.muted, "right"); pt.text(11.9, -1.85, "A", COL.muted, "right"); };
          pt.render();
        });
        ui.out("把小球放在 A 附近：等幅来回摆（临界稳定）。<br>放在 C 附近：越摆越小，停在 C（渐近稳定）。<br>放在 B 附近：滚走（不稳定）。<br>放在 A 左侧很远处：会越过 B 滚到 C——A 只是<b>局部</b>稳定。");
      },
      poles(el) {
        const ui = card(el, { title: "极点放置板", hint: "添加实极点或共轭极点对并<b>拖动它们</b>，右边实时显示冲激响应（每个极点的模态相加）。背景：绿色左半平面稳定，红色右半平面不稳定。" });
        const [a, b] = ui.row("", "");
        const ps = new Plot(a, { height: 300, x: [-4, 2], y: [-4, 4], equal: true, xlabel: "σ", title: "s 平面" });
        const pr = new Plot(b, { height: 300, x: [0, 10], y: [-2, 2], xlabel: "t", title: "冲激响应" });
        let poles = [{ r: -1, i: 0 }, { r: -0.4, i: 2 }], drag = -1;
        ui.buttons([["+ 实极点", () => { poles.push({ r: -2, i: 0 }); draw(); }, "ghost"], ["+ 共轭极点对", () => { poles.push({ r: -1, i: 1.5 }); draw(); }, "ghost"], ["清空", () => { poles = []; draw(); }, "ghost"]]);
        ps.drag({
          down: (x, y) => { drag = -1; let bd = 16; poles.forEach((p, k) => { [p.i, -p.i].forEach(im => { const d = ps.pxDist(p.r, im, x, y); if (d < bd) { bd = d; drag = k; } }); }); return drag >= 0 ? undefined : false; },
          move: (x, y) => { const p = poles[drag]; p.r = clamp(x, -3.9, 1.5); if (p.i !== 0) p.i = clamp(Math.abs(y), 0.2, 3.8); draw(); },
          hover: (x, y) => poles.some(p => ps.pxDist(p.r, p.i, x, y) < 16 || ps.pxDist(p.r, -p.i, x, y) < 16),
        });
        function draw() {
          const ts = linspace(0, 10, 500), ys = ts.map(t => poles.reduce((s, p) => s + (p.i === 0 ? Math.exp(p.r * t) : Math.exp(p.r * t) * Math.sin(p.i * t)), 0));
          const lim = Math.max(1, ...ys.slice(0, 350).map(Math.abs).filter(isFinite)); pr.setRange(null, [-Math.min(lim, 30) * 1.1, Math.min(lim, 30) * 1.1]);
          ps.render = () => { ps.clear(); ps.rect(-4, -4, 0, 4, COL.c3 + "14"); ps.rect(0, -4, 2, 4, COL.c4 + "14"); poles.forEach(p => { ps.cross(p.r, p.i, COL.c4, 8, 3); if (p.i) ps.cross(p.r, -p.i, COL.c4, 8, 3); }); };
          pr.render = () => { pr.clear(); pr.line(ts, ys, COL.c1, 2.4); };
          ps.render(); pr.render();
          const mx = poles.length ? Math.max(...poles.map(p => p.r)) : -1;
          ui.out(mx < -1e-3 ? '<span class="good">全部在左半平面：渐近稳定</span>' : mx > 1e-3 ? '<span class="bad">有极点在右半平面：不稳定</span>' : "最右的极点在虚轴上：临界稳定");
        }
        draw();
      },
      bibo(el) {
        const ui = card(el, { title: "共振实验：临界稳定 ≠ BIBO 稳定", hint: "系统 $G(s)=\\frac{\\omega_n^2}{s^2+2\\zeta\\omega_n s+\\omega_n^2}$，ωₙ = 2。输入幅值为 1 的正弦 sin(ωt)。" });
        const pr = new Plot(ui.plotHost(), { height: 260, x: [0, 40], y: [-10, 10], xlabel: "t", ylabel: "输出" });
        const sw = ui.slider({ label: "输入频率 ω", min: 0.5, max: 4, step: 0.01, value: 1.6, onInput: draw });
        const sz = ui.slider({ label: "阻尼比 ζ", min: 0, max: 0.5, step: 0.005, value: 0, fmt: v => v.toFixed(3), onInput: draw });
        function draw() {
          const w = sw.get(), z = sz.get(), wn = 2, r = CE.lsim([wn * wn], [1, 2 * z * wn, wn * wn], t => Math.sin(w * t), 40, 0.01);
          const mx = Math.max(...r.y.map(Math.abs)); pr.setRange(null, [-Math.max(2, mx * 1.1), Math.max(2, mx * 1.1)]);
          pr.render = () => { pr.clear(); pr.line(r.t, r.t.map(t => Math.sin(w * t)), COL.c2, 1.2, [4, 3]); pr.line(r.t, r.y, COL.c1, 2); pr.legend([{ label: "输入 sin ωt", color: COL.c2, dash: [4, 3] }, { label: "输出", color: COL.c1 }], "tr"); };
          pr.render();
          ui.out(`最大输出幅值：<b class="big">${f2(mx, 1)}</b><br>${z === 0 && Math.abs(w - 2) < 0.03 ? '<span class="bad">共振！输出随时间线性增长，无界</span>' : z === 0 ? "ζ=0：极点在虚轴上。把 ω 拖到 2 试试" : "有阻尼：输出总是有界的"}`);
        }
        draw();
      },
    },
  });
})();
