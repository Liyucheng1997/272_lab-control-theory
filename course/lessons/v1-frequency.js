/* 卷1 第9章：频率响应 / 伯德图 / 调音台 */
(function () {
  const { M, P, Plot, Scene, card, linspace, logspace, f2, COL, clamp, rk4, cfmt } = CE;
  const db = v => 20 * Math.log10(Math.max(v, 1e-12));
  const unwrap = arr => { const out = arr.slice(); for (let i = 1; i < out.length; i++) { while (out[i] - out[i - 1] > 180) out[i] -= 360; while (out[i] - out[i - 1] < -180) out[i] += 360; } return out; };

  const SYS = {
    int: { n: "积分器 1/s", tf: () => ({ num: [1], den: [1, 0] }) },
    first: { n: "一阶 a/(s+a)", tf: p => ({ num: [p.a], den: [1, p.a] }) },
    second: { n: "二阶", tf: p => ({ num: [p.wn * p.wn], den: [1, 2 * p.z * p.wn, p.wn * p.wn] }) },
  };

  Course.add({
    id: "v1-freq", vol: 1, ord: 150, no: "9.1–9.4", chName: "第9章 频率响应与分析", ref: "卷1 §9.1–9.4",
    title: "频率响应：正弦进，正弦出",
    sub: "给稳定的 LTI 系统输入正弦波，输出还是同频率的正弦波——只是振幅乘了 |G(jω)|、相位挪了 ∠G(jω)。",
    goals: ["理解频率响应的核心结论及其推导思路", "知道一阶系统是低通滤波器，截止频率处 0.707、−45°", "理解二阶系统的共振与阻尼比的关系"],
    minutes: 30,
    html: String.raw`
<div class="story" data-t="引子：百万调音师">打开音乐播放器里的“均衡器”：10 个频段，从 32 Hz 到 16 kHz，每个都能在 ±12 dB 之间调节。鼓声集中在 32～125 Hz，人声在 500～2000 Hz，一些乐器和清晰的对白在更高的频率上。现代电影配乐的调音师、录音师都是工程师——均衡器背后的数学，正是控制理论里的<b>频率响应</b>。</div>

<h2><span class="sec">9.2</span>核心结论</h2>
<p>输入一般的正弦信号 $u(t)=A\sin\omega_it+B\cos\omega_it=M_i\sin(\omega_it+\varphi_i)$（$M_i=\sqrt{A^2+B^2}$）。只要系统<b>稳定</b>，过渡过程消失后：</p>
<div class="eq" data-no="(9.2.18)" data-title="LTI 系统的频率响应（本章最重要的结论）">
$$x_{ss}(t)=M_o\sin(\omega_it+\varphi_o),\qquad \frac{M_o}{M_i}=\left|G(j\omega_i)\right|,\qquad \varphi_o=\varphi_i+\angle G(j\omega_i)$$
<dl><dt>$G(j\omega_i)$</dt><dd>把 $s=j\omega_i$ 代入传递函数得到的复数</dd><dt>$|G(j\omega_i)|$</dt><dd><b>振幅响应</b>：输出振幅是输入的多少倍</dd><dt>$\angle G(j\omega_i)$</dt><dd><b>相位响应</b>：输出比输入超前（正）或滞后（负）多少</dd></dl>
<div class="say">频率不变；振幅乘 |G|；相位加 ∠G。</div>
</div>
<details class="derive"><summary>推导思路（部分分式 + 共轭对称）</summary><div class="body">
<div class="step">$U(s)=\frac{A\omega_i+Bs}{(s+j\omega_i)(s-j\omega_i)}$，于是 $X=UG=\frac{K_1}{s+j\omega_i}+\frac{K_2}{s-j\omega_i}+\sum\frac{C_k}{s-s_{pk}}$。</div>
<div class="step">系统稳定 ⇒ 所有 $C_ke^{s_{pk}t}\to0$，稳态只剩 $K_1e^{-j\omega_it}+K_2e^{j\omega_it}$。</div>
<div class="step">留数法：$K_1=\frac{B+Aj}{2}G(-j\omega_i)$，$K_2=\frac{B-Aj}{2}G(j\omega_i)$。</div>
<div class="step">实信号的 $G(-j\omega)$ 与 $G(j\omega)$ 共轭：模相同、角相反。代入 $G=|G|e^{j\angle G}$ 并用欧拉公式合并，得 $x_{ss}=|G|M_i\sin(\omega_it+\varphi_i+\angle G)$。</div>
</div></details>
<p><b>例：积分器</b> $G=\frac1s$，$G(j\omega)=\frac{1}{j\omega}=-\frac{j}{\omega}$：振幅 $1/\omega$、相位 $-90°$。所以 $\sin2t$ 积分后是 $0.5\sin(2t-\frac\pi2)$——直接积分 $\int\sin2t\,\mathrm{d}t=-0.5\cos2t$ 验证一致。频率越高输出越小：积分器是<b>低通滤波器</b>。</p>
<div class="sim" data-id="scope"></div>

<h2><span class="sec">9.3</span>一阶系统：低通滤波器</h2>
<div class="eq" data-no="(9.3.2)" data-title="一阶系统 $G(s)=\frac{a}{s+a}$">
$$|G(j\omega)|=\frac{1}{\sqrt{(\omega/a)^2+1}},\qquad \angle G(j\omega)=-\arctan\frac{\omega}{a}$$
<dl><dt>$\omega=0$</dt><dd>$|G|=1$，相位 0：慢变化原样通过</dd><dt>$\omega=a$</dt><dd>$|G|=0.707$，相位 $-45°$：<b>截止频率</b></dd><dt>$\omega\to\infty$</dt><dd>$|G|\to0$，相位 $\to-90°$：快变化被滤掉</dd></dl>
</div>
<p>噪声通常是高频、小振幅的信号，所以一阶低通滤波器最常见的用途就是<b>降噪</b>（信号处理中称一阶巴特沃思滤波器）。反过来，$\frac{s}{s+a}$ 是高通滤波器。</p>
<div class="sim" data-id="lowpass"></div>
<div class="idea" data-t="为什么一阶系统都是低通？——“容器”">RC 电路里的电容、密室里的大量空气、人体的脂肪……一阶系统里总有一个“容器”。短时间暴饮暴食不会立刻体现在体重上，空调快速开关不会立刻体现在体温上——容器提供缓冲，带来延迟，抵消了输入的快速变化。</div>

<h2><span class="sec">9.4</span>二阶系统：共振</h2>
<div class="eq" data-no="(9.4.3)" data-title="二阶系统，归一化频率 $\Omega=\omega/\omega_n$">
$$|G(j\omega)|=\frac{1}{\sqrt{(1-\Omega^2)^2+4\zeta^2\Omega^2}},\qquad \angle G(j\omega)=-\arctan\frac{2\zeta\Omega}{1-\Omega^2}$$
<div class="say">在 $\Omega=1$ 处 $|G|=\frac{1}{2\zeta}$：ζ&lt;0.5 时被放大。当 $\zeta<0.707$ 时幅频曲线有峰值，出现在<b>共振频率</b> $\omega_R=\omega_n\sqrt{1-2\zeta^2}$。</div>
</div>
<div class="sim" data-id="resonance"></div>
<p>共振无处不在：小提琴的琴箱靠共振放大琴弦的声音；塔科马海峡大桥建成四个月就被并不剧烈的风吹垮；坐公交时某些时刻车身剧烈抖动，是因为发动机振动频率接近车体的固有频率。</p>
<div class="warn" data-t="作者的提醒">书中用“找到对方的固有频率才能产生共鸣”来类比沟通，随即又提醒：这些类比只是帮助理解，同一个原理换个说法也能论证相反的观点（比如把“容器”说成“过去的经验成了包袱”）。选择性地借用科学原理来论证观点，是“披着科学外衣的伪科学”——要保持思辨。</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>稳定 LTI 系统对正弦输入的稳态响应：同频，振幅 ×|G(jω)|，相位 +∠G(jω)。</li>
<li>一阶系统是低通滤波器：截止频率 ω=a 处 0.707、−45°；可用于降噪。</li>
<li>二阶系统：阻尼越大响应越迟缓；ζ&lt;0.707 时在 $\omega_R=\omega_n\sqrt{1-2\zeta^2}$ 附近共振。</li>
</ol></div>`,
    sims: {
      scope(el) {
        const ui = card(el, { title: "示波器：正弦进，正弦出", hint: "上：输入（橙）和输出（蓝，含过渡过程）。左下：复数 G(jω) 在复平面上的向量——长度是振幅比、角度是相移。右下：伯德图上当前频率的位置。" });
        const pt = new Plot(ui.plotHost(), { height: 220, x: [0, 20], y: [-1.6, 1.6], xlabel: "t" });
        const [a, b] = ui.row("", "");
        const pc = new Plot(a, { height: 230, x: [-1.4, 1.4], y: [-1.4, 1.4], equal: true, xlabel: "Re", title: "G(jω)" });
        const pb = new Plot(b, { height: 230, x: [0.05, 50], y: [-40, 15], logx: true, xlabel: "ω (rad/s)", title: "幅频 |G| (dB)" });
        let sys = "first";
        ui.seg({ label: "系统", options: Object.entries(SYS).map(([k, v]) => [k, v.n]), value: sys, onChange: v => { sys = v; draw(); } });
        const sw = ui.slider({ label: "输入频率 ω", min: 0.1, max: 20, value: 2, log: true, onInput: draw });
        const sa = ui.slider({ label: "一阶：a", min: 0.2, max: 10, value: 1, log: true, onInput: draw });
        const sz = ui.slider({ label: "二阶：ζ", min: 0.05, max: 2, step: 0.01, value: 0.2, onInput: draw });
        const sn = ui.slider({ label: "二阶：ωₙ", min: 0.5, max: 8, value: 2, log: true, onInput: draw });
        ui.buttons([["例：sin 2t 过积分器", () => { sys = "int"; sw.set(2); draw(); }, "ghost"]]);
        function draw() {
          const w = sw.get(), pr = { a: sa.get(), z: sz.get(), wn: sn.get() }, { num, den } = SYS[sys].tf(pr);
          const g = P.freq(num, den, w), Tend = Math.max(20, 6 * 2 * Math.PI / w);
          const r = CE.lsim(num, den, t => Math.sin(w * t), Math.min(Tend, 40), Math.min(0.01, 0.05 / w));
          const mx = Math.max(1.2, Math.min(5, g.mag * 1.25));
          pt.setRange([0, Math.min(Tend, 40)], [-mx, mx]);
          pt.render = () => { pt.clear(); pt.line(r.t, r.t.map(t => Math.sin(w * t)), COL.c2, 1.8); pt.line(r.t, r.y, COL.c1, 2.4); pt.line(r.t, r.t.map(t => g.mag * Math.sin(w * t + g.ph * Math.PI / 180)), COL.c1 + "55", 5); pt.legend([{ label: "输入 sin ωt", color: COL.c2 }, { label: "输出 x(t)", color: COL.c1 }, { label: "理论稳态 |G|sin(ωt+∠G)", color: COL.c1 + "55" }], "tr"); };
          const R = Math.max(1.2, g.mag * 1.15); pc.setRange([-R, R], [-R, R]);
          pc.render = () => { pc.clear(); const th = linspace(0, 2 * Math.PI, 100); pc.line(th.map(Math.cos), th.map(Math.sin), COL.axis, 1, [3, 3]); pc.arrow(0, 0, g.re, g.im, COL.c5, 2.6, 10); pc.text(g.re / 2 + 0.05, g.im / 2 + 0.12, "|G|=" + f2(g.mag), COL.c5); };
          const ws = logspace(-1.3, 1.7, 300), mags = ws.map(x => db(P.freq(num, den, x).mag));
          pb.render = () => { pb.clear(); pb.line(ws, mags, COL.c1, 2.2); pb.dot(w, db(g.mag), COL.c4, 6); pb.vline(w, COL.c4, [3, 3]); };
          pt.render(); pc.render(); pb.render();
          ui.out(`ω = <b>${f2(w)}</b> rad/s<br>|G(jω)| = <b>${f2(g.mag, 3)}</b>（${f2(db(g.mag), 1)} dB）<br>∠G(jω) = <b>${f2(g.ph, 1)}°</b><br>时间滞后 = ${f2(-g.ph * Math.PI / 180 / w, 2)} s`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>等过渡过程消失后，蓝线与半透明粗线（理论稳态）完全重合。</li><li>一阶系统：把 ω 拖到 a，|G|=0.707，相位 −45°。</li><li>二阶系统 ζ=0.05：把 ω 拖到 ωₙ 附近，输出振幅变得巨大——共振。</li></ol></div>`);
      },
      lowpass(el) {
        const ui = card(el, { title: "一阶低通滤波器降噪", hint: "有效信号 sin t（ω=1）混入高频噪声。调节截止频率 a，找到“噪声滤得干净”和“信号不失真”之间的平衡。" });
        const pt = new Plot(ui.plotHost(), { height: 250, x: [0, 15], y: [-2, 2], xlabel: "t" });
        const pb = new Plot(ui.plotHost(), { height: 170, x: [0.05, 200], y: [-40, 5], logx: true, xlabel: "ω", title: "滤波器幅频 (dB)：绿 = 信号频率，红 = 噪声频率" });
        const sa = ui.slider({ label: "截止频率 a", min: 0.1, max: 100, value: 1, log: true, onInput: draw });
        const sn = ui.slider({ label: "噪声幅度", min: 0, max: 1, step: 0.05, value: 0.6, onInput: draw });
        const sf = ui.slider({ label: "噪声频率 (rad/s)", min: 5, max: 150, step: 1, value: 50, fmt: v => v.toFixed(0), onInput: draw });
        let hp = false;
        ui.check({ label: "改成高通 s/(s+a)", value: false, onChange: v => { hp = v; draw(); } });
        function draw() {
          const a = sa.get(), A = sn.get(), wf = sf.get(), rg = CE.rng(3);
          const ph = Array.from({ length: 6 }, () => rg.u() * 6.28);
          const u = t => Math.sin(t) + A * (Math.sin(wf * t + ph[0]) * 0.6 + Math.sin(1.37 * wf * t + ph[1]) * 0.4);
          const num = hp ? [1, 0] : [a], den = [1, a];
          const r = CE.lsim(num, den, u, 15, 0.002);
          pt.render = () => { pt.clear(); pt.line(r.t, r.t.map(u), COL.c2 + "88", 1); pt.line(r.t, r.t.map(Math.sin), COL.c3, 1.4, [5, 4]); pt.line(r.t, r.y, COL.c1, 2.6); pt.legend([{ label: "带噪输入", color: COL.c2 }, { label: "真实信号 sin t", color: COL.c3, dash: [5, 4] }, { label: "滤波输出", color: COL.c1 }], "tr"); };
          const ws = logspace(-1.3, 2.3, 300);
          pb.render = () => { pb.clear(); pb.line(ws, ws.map(x => db(P.freq(num, den, x).mag)), COL.c1, 2.2); pb.vline(1, COL.c3, [4, 3], 2); pb.vline(wf, COL.c4, [4, 3], 2); };
          pt.render(); pb.render();
          const g1 = P.freq(num, den, 1), gn = P.freq(num, den, wf);
          ui.out(`信号 (ω=1)：保留 <b>${f2(g1.mag * 100, 0)}%</b>，相位 ${f2(g1.ph, 0)}°<br>噪声 (ω=${wf})：保留 <b>${f2(gn.mag * 100, 1)}%</b>`);
        }
        draw();
      },
      resonance(el) {
        const ui = card(el, { title: "二阶系统共振探索器", hint: "灰线为书中图 9.4.2 的六种阻尼比 (ζ=0,0.2,0.5,0.7,1,2)。点“扫频”，让输入频率从低到高慢慢扫过，看输出振幅怎样在共振频率附近鼓起来。" });
        const [a, b] = ui.row("", "");
        const pm = new Plot(a, { height: 260, x: [0, 3], y: [0, 5], xlabel: "Ω = ω/ωₙ", title: "振幅 |G(jω)|（线性坐标）" });
        const pp = new Plot(b, { height: 260, x: [0, 3], y: [-180, 0], xlabel: "Ω", title: "相位 (°)", yticks: [0, -45, -90, -135, -180] });
        const pt = new Plot(ui.plotHost(), { height: 200, x: [0, 60], y: [-5, 5], xlabel: "t", title: "扫频实验：输入频率随时间线性增大 (0 → 3ωₙ)" });
        const sz = ui.slider({ label: "阻尼比 ζ", min: 0.02, max: 2, step: 0.01, value: 0.15, onInput: draw });
        ui.buttons([["▶ 扫频", sweep]]);
        const mag = (z, O) => 1 / Math.sqrt((1 - O * O) ** 2 + 4 * z * z * O * O), pha = (z, O) => Math.atan2(-2 * z * O, 1 - O * O) * 180 / Math.PI;
        let sw = null;
        function sweep() {
          const z = sz.get(), wn = 1, T = 60, dt = 0.01; let x = [0, 0]; const ts = [], ys = [], us = [];
          for (let t = 0; t <= T; t += dt) { const ph = 3 / T * t * t / 2; /* 瞬时频率 3t/T */ const uu = Math.sin(ph); x = rk4(s => [s[1], -2 * z * wn * s[1] - wn * wn * s[0] + wn * wn * uu], x, t, dt); ts.push(t); ys.push(x[0]); us.push(uu); }
          sw = { ts, ys, us }; const m = Math.max(1.5, ...ys.map(Math.abs)) * 1.1; pt.setRange(null, [-m, m]); pt.render();
        }
        pt.render = () => { pt.clear(); if (!sw) { pt.text(30, 0, "点击 ▶ 扫频", COL.muted, "center", 13); return; } pt.line(sw.ts, sw.us, COL.c2 + "99", 1); pt.line(sw.ts, sw.ys, COL.c1, 1.6); pt.vline(20, COL.c4, [4, 3]); pt.text(20.4, pt.opt.y[1] * 0.85, "Ω = 1", COL.c4); };
        function draw() {
          const z = sz.get(), Os = linspace(0.001, 3, 400);
          const pk = z < 0.707 ? mag(z, Math.sqrt(1 - 2 * z * z)) : 1;
          pm.setRange(null, [0, Math.max(2, Math.min(pk * 1.15, 12))]);
          pm.render = () => {
            pm.clear(); [0.2, 0.5, 0.7, 1, 2].forEach(zz => pm.line(Os, Os.map(O => mag(zz, O)), COL.axis, 1.2)); pm.line(Os, Os.map(O => mag(0.0001, O)), COL.axis, 1, [3, 3]);
            pm.line(Os, Os.map(O => mag(z, O)), COL.c1, 2.8); pm.vline(1, COL.muted, [3, 3]);
            if (z < 0.707) { const OR = Math.sqrt(1 - 2 * z * z); pm.dot(OR, pk, COL.c4, 5); pm.text(OR + 0.05, Math.min(pk, pm.opt.y[1] * 0.95), "ω_R, 峰值 " + f2(pk, 2), COL.c4); }
          };
          pp.render = () => { pp.clear(); pp.line(Os, Os.map(O => pha(z, O)), COL.c1, 2.6); pp.hline(-90, COL.muted, [3, 3]); };
          pm.render(); pp.render(); pt.render();
          ui.out(`Ω=1 处 |G| = 1/(2ζ) = <b>${f2(1 / (2 * z), 2)}</b><br>${z < 0.707 ? `共振频率 ω_R = ωₙ√(1−2ζ²) = <b>${f2(Math.sqrt(1 - 2 * z * z), 3)}ωₙ</b><br>峰值 = <b>${f2(pk, 2)}</b>` : "ζ ≥ 0.707：没有共振峰"}`);
        }
        draw();
      },
    },
  });

  /* ═════ 9.5 伯德图 ═════ */
  const BLOCKS = {
    K: { n: "增益 K", tf: p => ({ num: [p.v], den: [1] }), def: 0.1, range: [0.01, 100] },
    int: { n: "积分器 1/s", tf: () => ({ num: [1], den: [1, 0] }) },
    lp: { n: "一阶低通 1/(s/a+1)", tf: p => ({ num: [1], den: [1 / p.v, 1] }), def: 10, range: [0.01, 1000] },
    pd: { n: "PD  s/a+1", tf: p => ({ num: [1 / p.v, 1], den: [1] }), def: 1, range: [0.01, 1000] },
  };
  Course.add({
    id: "v1-bode", vol: 1, ord: 160, no: "9.5", chName: "第9章 频率响应与分析", ref: "卷1 §9.5–9.6",
    title: "伯德图：用积木搭出滤波器和调音台",
    sub: "对数坐标让“乘法变加法”：串联系统的伯德图就是各部分伯德图的叠加。拿着这把钥匙，可以反向设计一个均衡器。",
    goals: ["理解分贝 dB 与对数频率轴", "记住典型环节的伯德图：K、1/s、一阶低通、PD", "理解超前/滞后补偿器名字的由来", "会把目标幅频曲线拆成若干一阶环节，设计一个调音台"],
    minutes: 30,
    html: String.raw`
<p><span class="kw">伯德图</span>（以荷兰裔美国工程师 Hendrik Bode 命名）是最常用的频率响应图：上半是<b>幅频图</b>（纵轴 dB），下半是<b>相频图</b>（纵轴度），横轴都是按 10 倍分度的对数频率。</p>
<h2><span class="sec">9.5.1</span>分贝与“乘法变加法”</h2>
<div class="eq" data-no="(9.5.1)–(9.5.2)" data-title="分贝">
$$L_{dB}=10\log_{10}\frac{P_m}{P_r}=10\log_{10}\left(\frac{M_o}{M_i}\right)^2=20\log_{10}|G(j\omega)|$$
<div class="say">功率与振幅的平方成正比，所以振幅用 20 log。振幅 ×10 ↔ +20 dB；×0.1 ↔ −20 dB；×0.707 ↔ −3 dB；×2 ↔ +6 dB。日常交谈约 60 dB、闹市约 80 dB，相差 20 dB 就是功率差 100 倍。</div>
</div>
<div class="eq" data-no="(9.5.6)" data-title="伯德图最重要的性质">
$$20\log|G_1G_2|=20\log|G_1|+20\log|G_2|,\qquad \angle(G_1G_2)=\angle G_1+\angle G_2$$
<div class="say">串联系统的伯德图 = 各子系统伯德图相加。掌握几个典型环节，就能“搭积木”分析和设计复杂系统。</div>
</div>

<h2><span class="sec">9.5.2</span>典型环节</h2>
<table class="tbl"><tr><th>环节</th><th>幅频</th><th>相频</th><th>备注</th></tr>
<tr><td>比例 $K$</td><td>水平线 $20\log K$</td><td>0°</td><td>只缩放，不延迟</td></tr>
<tr><td>积分器 $\frac1s$</td><td>−20 dB/dec，过 (1, 0 dB)</td><td>−90°</td><td>低频放大 → 需外部能源</td></tr>
<tr><td>一阶 $\frac{1}{s/a+1}$</td><td>0 dB，在 a 处转折为 −20 dB/dec</td><td>0° → −90°，a 处 −45°</td><td>a 处 −3 dB</td></tr>
<tr><td>PD $\frac sa+1$</td><td>0 dB，在 a 处转折为 +20 dB/dec</td><td>0° → +90°</td><td>放大高频噪声</td></tr></table>
<p>例 9.5.5：<b>超前补偿器</b> $\frac{s+1}{s+10}=\frac{1}{10}\cdot(s+1)\cdot\frac{1}{s/10+1}$，由 增益 + PD(转折 1) + 低通(转折 10) 叠加：幅频从 −20 dB 升到 0 dB，相位是一个<b>正</b>的鼓包（约 +55°，在 $\sqrt{10}$ 处）——相位超前，故名“超前”。它不像 PD 那样无限放大高频。例 9.5.6：<b>滞后补偿器</b> $\frac{0.1s+1}{s+1}$ 恰好相反，相位是负的凹谷。</p>
<div class="sim" data-id="blocks"></div>
<div class="idea" data-t="控制器视角">积分器和滞后补偿器相位滞后，关注“过去”的积累，帮助<b>消除稳态误差</b>；PD 和超前补偿器相位超前，会“提前预测”，帮助<b>提高响应速度</b>。第 7、8 章的结论在频域里得到了印证。</div>

<h2><span class="sec">9.5.3</span>调音台的设计</h2>
<p>把上面的思路反过来用：先把均衡器的目标形状画成幅频图（低音 +10 dB、中音 +5 dB、高音 0 dB 的三台阶），再拆成几个一阶环节——一个总增益、64 Hz 处的极点、125 Hz 处的零点、1 kHz 处的极点、2 kHz 处的零点——相乘就是调音台的传递函数。（具体增益数值是本页补充的一种实现。）</p>
<div class="eq" data-title="一种可实现的均衡器">
$$G_{EQ}(s)=K\cdot\frac{s/\omega_3+1}{s/\omega_2+1}\cdot\frac{s/\omega_5+1}{s/\omega_4+1},\quad \omega_k=2\pi f_k$$
</div>
<div class="sim" data-id="mixer"></div>
<div class="tip">一阶环节的斜率只有 ±20 dB/dec，真实的调音台需要更陡的过渡，会使用更高阶的滤波器——思路完全相同。</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>伯德图 = 对数频率轴上的幅频（dB）与相频（度）曲线；$20\log|G|$。</li>
<li>串联系统的伯德图是各部分的叠加：幅值 dB 相加，相位相加。</li>
<li>一阶环节在转折频率处 ±3 dB、±45°，之后 ±20 dB/dec。</li>
<li>超前补偿器相位为正，滞后补偿器相位为负；用积木叠加的思路可以设计滤波器。</li>
</ol></div>`,
    sims: {
      blocks(el) {
        const ui = card(el, { title: "伯德图积木叠加器", hint: "细线是每块积木各自的伯德图，粗线是它们的总和。拖动滑块改转折频率/增益，或载入书中的超前、滞后补偿器。" });
        const pm = new Plot(ui.plotHost(), { height: 230, x: [0.01, 1000], y: [-60, 40], logx: true, xlabel: "ω (rad/s)", title: "幅频 (dB)" });
        const pp = new Plot(ui.plotHost(), { height: 190, x: [0.01, 1000], y: [-100, 100], logx: true, xlabel: "ω", title: "相频 (°)", yticks: [-90, -45, 0, 45, 90] });
        let list = [];
        const box = document.createElement("div"); box.className = "sim-controls"; ui.controls.appendChild(box);
        const asy = ui.check({ label: "显示渐近线（虚线）", value: false, onChange: draw });
        const addSel = document.createElement("div"); addSel.className = "btn-row"; ui.controls.appendChild(addSel);
        Object.entries(BLOCKS).forEach(([k, b]) => { const bt = document.createElement("button"); bt.className = "btn ghost"; bt.textContent = "+ " + b.n.split(" ")[0]; bt.onclick = () => { list.push({ k, v: b.def }); rebuild(); }; addSel.appendChild(bt); });
        ui.buttons([["超前 (s+1)/(s+10)", () => { list = [{ k: "K", v: 0.1 }, { k: "pd", v: 1 }, { k: "lp", v: 10 }]; rebuild(); }], ["滞后 (0.1s+1)/(s+1)", () => { list = [{ k: "lp", v: 1 }, { k: "pd", v: 10 }]; rebuild(); }], ["清空", () => { list = []; rebuild(); }, "ghost"]]);
        const cols = ["c2", "c3", "c5", "c6", "c4"];
        function rebuild() {
          box.innerHTML = "";
          list.forEach((it, i) => {
            const b = BLOCKS[it.k]; const row = document.createElement("div");
            if (b.range) { CE.slider(row, { label: `<span style="color:var(--${cols[i % 5]})">■</span> ${b.n}${it.k === "K" ? "" : "：a"}`, min: b.range[0], max: b.range[1], value: it.v, log: true, onInput: v => { it.v = v; draw(); } }); }
            else row.innerHTML = `<div class="ctl-top"><label><span style="color:var(--${cols[i % 5]})">■</span> ${b.n}</label></div>`;
            const x = document.createElement("button"); x.className = "btn ghost"; x.textContent = "移除"; x.style.cssText = "padding:1px 8px;font-size:.75rem;margin-top:2px"; x.onclick = () => { list.splice(i, 1); rebuild(); };
            row.appendChild(x); box.appendChild(row);
          });
          draw();
        }
        function draw() {
          const ws = logspace(-2, 3, 400);
          const parts = list.map(it => { const { num, den } = BLOCKS[it.k].tf(it); const fr = ws.map(w => P.freq(num, den, w)); return { m: fr.map(g => db(g.mag)), p: unwrap(fr.map(g => g.ph)), it }; });
          const tm = ws.map((_, i) => parts.reduce((s, q) => s + q.m[i], 0)), tp = ws.map((_, i) => parts.reduce((s, q) => s + q.p[i], 0));
          const asym = it => w => it.k === "K" ? [db(it.v), 0] : it.k === "int" ? [-db(w), -90] : it.k === "lp" ? [w < it.v ? 0 : -db(w / it.v), w < it.v ? 0 : -90] : [w < it.v ? 0 : db(w / it.v), w < it.v ? 0 : 90];
          pm.render = () => { pm.clear(); parts.forEach((q, i) => pm.line(ws, q.m, COL[cols[i % 5]] + "aa", 1.3)); if (asy.get() && list.length) pm.line(ws, ws.map(w => list.reduce((s, it) => s + asym(it)(w)[0], 0)), COL.ink2, 1.4, [6, 4]); pm.line(ws, tm, COL.c1, 3); };
          pp.render = () => { pp.clear(); parts.forEach((q, i) => pp.line(ws, q.p, COL[cols[i % 5]] + "aa", 1.3)); if (asy.get() && list.length) pp.line(ws, ws.map(w => list.reduce((s, it) => s + asym(it)(w)[1], 0)), COL.ink2, 1.4, [6, 4]); pp.line(ws, tp, COL.c1, 3); };
          pm.render(); pp.render();
          const pk = tp.reduce((b, v, i) => Math.abs(v) > Math.abs(b[0]) ? [v, ws[i]] : b, [0, 1]);
          ui.out(list.length ? `最大相移 <b>${f2(pk[0], 1)}°</b>，出现在 ω = ${f2(pk[1], 2)}` : "添加一些积木吧");
        }
        list = [{ k: "K", v: 0.1 }, { k: "pd", v: 1 }, { k: "lp", v: 10 }]; rebuild();
      },
      mixer(el) {
        const ui = card(el, { title: "动手做调音台（可以真的听！）", hint: "上方 10 个点是目标均衡器设定（可拖动），蓝线是你用 5 块积木搭出的幅频曲线。点“播放”听一段合成音乐（鼓 + 人声 + 镲片），切换“均衡器 开/关”对比效果。" });
        const pm = new Plot(ui.plotHost(), { height: 300, x: [20, 20000], y: [-15, 15], logx: true, xlabel: "频率 f (Hz)", ylabel: "增益 dB" });
        const bands = [32, 64, 125, 250, 500, 1000, 2000, 4000, 8000, 16000];
        const target = [10, 10, 5, 5, 5, 5, 0, 0, 0, 0];
        const sK = ui.slider({ label: "总增益 K (dB)", min: -12, max: 15, step: 0.1, value: 11.8, onInput: draw });
        const s2 = ui.slider({ label: "极点 f₂ (Hz)", min: 20, max: 500, value: 64, log: true, fmt: v => v.toFixed(0), onInput: draw });
        const s3 = ui.slider({ label: "零点 f₃ (Hz)", min: 30, max: 1000, value: 125, log: true, fmt: v => v.toFixed(0), onInput: draw });
        const s4 = ui.slider({ label: "极点 f₄ (Hz)", min: 200, max: 8000, value: 1000, log: true, fmt: v => v.toFixed(0), onInput: draw });
        const s5 = ui.slider({ label: "零点 f₅ (Hz)", min: 300, max: 16000, value: 2000, log: true, fmt: v => v.toFixed(0), onInput: draw });
        let eqOn = true, ctx = null, src = null, filt = null, gainNode = null;
        const btnRow = ui.buttons([["▶ 播放", play], ["■ 停止", stop, "ghost"]]);
        ui.check({ label: "均衡器 开", value: true, onChange: v => { eqOn = v; if (src) restartAudio(); } });
        let dragI = -1;
        pm.drag({ down: (x, y) => { dragI = -1; bands.forEach((f, i) => { if (pm.pxDist(f, target[i], x, y) < 14) dragI = i; }); return dragI >= 0 ? undefined : false; }, move: (x, y) => { target[dragI] = clamp(Math.round(y), -12, 12); draw(); }, hover: (x, y) => bands.some((f, i) => pm.pxDist(f, target[i], x, y) < 14) });
        function tf() {
          const K = Math.pow(10, sK.get() / 20), w = [s2, s3, s4, s5].map(s => 2 * Math.PI * s.get());
          return { num: P.scale(P.mul([1 / w[1], 1], [1 / w[3], 1]), K), den: P.mul([1 / w[0], 1], [1 / w[2], 1]) };
        }
        function draw() {
          const { num, den } = tf(), fs = logspace(Math.log10(20), Math.log10(20000), 300);
          const resp = fs.map(f => db(P.freq(num, den, 2 * Math.PI * f).mag));
          pm.render = () => {
            pm.clear(); bands.forEach((f, i) => { pm.vline(f, COL.grid, null, 1); });
            [s2, s4].forEach(s => pm.vline(s.get(), COL.c4 + "88", [4, 3])); [s3, s5].forEach(s => pm.vline(s.get(), COL.c5 + "88", [4, 3]));
            pm.line(fs, resp, COL.c1, 3);
            bands.forEach((f, i) => { pm.rect(f / 1.08, -12, f * 1.08, target[i], COL.c2 + "30"); pm.dot(f, target[i], COL.c2, 7, COL.plotbg); });
          };
          pm.render();
          const err = bands.map((f, i) => Math.abs(db(P.freq(num, den, 2 * Math.PI * f).mag) - target[i]));
          ui.out(`与目标的平均偏差：<b class="big">${f2(err.reduce((a, b) => a + b, 0) / err.length, 1)} dB</b><br><span class="sim-note">红色虚线 = 极点（向下转折），紫色虚线 = 零点（向上转折）。</span>`);
          if (src) restartAudio();
        }
        // ─── 音频 ───
        function makeBuffer(ac) {
          const sr = ac.sampleRate, T = 4, buf = ac.createBuffer(1, sr * T, sr), d = buf.getChannelData(0), rg = CE.rng(11);
          for (let i = 0; i < d.length; i++) {
            const t = i / sr, beat = t % 0.5, half = t % 0.25;
            const kick = Math.sin(2 * Math.PI * (55 + 60 * Math.exp(-beat * 30)) * beat) * Math.exp(-beat * 9) * 0.5;
            const note = [220, 247, 262, 294][Math.floor(t / 1) % 4];
            const vib = 1 + 0.004 * Math.sin(2 * Math.PI * 5 * t);
            let voice = 0; for (let h = 1; h <= 6; h++) voice += Math.sin(2 * Math.PI * note * h * vib * t) / (h * 1.3);
            voice *= 0.09 * (0.6 + 0.4 * Math.sin(Math.PI * (t % 1)));
            const hat = (rg.u() * 2 - 1) * Math.exp(-half * 60) * 0.12;
            d[i] = (kick + voice + hat) * 0.5;
          }
          // 简单高通处理镲片噪声：一阶差分
          return buf;
        }
        function bilinear(num, den, fs) {
          // 二阶 s 域 → z 域（num, den 均为 [a2,a1,a0]）
          const pad = c => [0, 0, 0].concat(c).slice(-3), n = pad(num), d = pad(den), k = 2 * fs;
          const B = [n[0] * k * k + n[1] * k + n[2], -2 * n[0] * k * k + 2 * n[2], n[0] * k * k - n[1] * k + n[2]];
          const A = [d[0] * k * k + d[1] * k + d[2], -2 * d[0] * k * k + 2 * d[2], d[0] * k * k - d[1] * k + d[2]];
          return { b: B.map(v => v / A[0]), a: A.map(v => v / A[0]) };
        }
        function restartAudio() { const t = ctx.currentTime; stop(); play(t); }
        function play() {
          try {
            ctx = ctx || new (window.AudioContext || window.webkitAudioContext)();
            if (ctx.state === "suspended") ctx.resume();
            if (!ctx._buf) ctx._buf = makeBuffer(ctx);
            stop();
            src = ctx.createBufferSource(); src.buffer = ctx._buf; src.loop = true;
            gainNode = ctx.createGain(); gainNode.gain.value = 0.5;
            if (eqOn) { const { num, den } = tf(), { b, a } = bilinear(num, den, ctx.sampleRate); filt = ctx.createIIRFilter(b, a); src.connect(filt); filt.connect(gainNode); gainNode.gain.value = 0.5 / Math.pow(10, Math.max(0, sK.get()) / 20) * 1.6; }
            else src.connect(gainNode);
            gainNode.connect(ctx.destination); src.start();
          } catch (e) { ui.out("浏览器不支持音频播放：" + e.message); }
        }
        function stop() { try { if (src) { src.stop(); src.disconnect(); } } catch (e) { } src = null; }
        CE.loop(ui.root, () => { if (!ui.root.isConnected && src) stop(); });
        window.addEventListener("hashchange", stop, { once: true });
        draw();
      },
    },
  });
})();
