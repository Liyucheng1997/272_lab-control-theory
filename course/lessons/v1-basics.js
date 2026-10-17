/* 卷1 第1–2章：绪论 / 卷积 / 微分方程建模 / 拉普拉斯变换与传递函数 */
(function () {
  const { M, P, Plot, Scene, card, linspace, f2, COL, clamp, rk4 } = CE;

  /* ════════════════ 1.1 什么是控制 ════════════════ */
  Course.add({
    id: "v1-intro", vol: 1, ord: 10, no: "1", chName: "第1章 绪论", ref: "卷1 第1章",
    title: "什么是动态系统与控制",
    sub: "从一辆推着走的小车出发，弄清楚“输入、输出、开环、闭环、反馈”这几个最基本的词。",
    goals: ["说清楚什么是动态系统、状态变量、输入与输出", "判断一个系统是不是线性时不变（LTI）", "理解开环控制为什么怕扰动、闭环控制为什么能自我纠正"],
    minutes: 15,
    html: String.raw`
<p>我们每天都被自动控制包围：空调把室温稳在 26℃，冰箱让饮料冰凉却不结冰，扫地机器人自己绕开桌脚，手机根据环境光调屏幕亮度，工厂里的机械臂以人类做不到的精度重复装配。这些东西背后是同一门学问——<span class="kw">控制理论</span>：用数学描述“会随时间变化的东西”，再设计一个“大脑”让它按我们的意愿变化。</p>

<h2><span class="sec">1.1</span>动态系统：状态随时间变化的系统</h2>
<p>控制理论研究的对象叫<span class="kw">动态系统</span>（Dynamic System）。它的特点是：描述它的量（<span class="kw">状态变量</span>）是时间的函数。最简单的例子是一辆放在光滑地面上的小车，我们用力 $u(t)$ 推它，它的位置 $x(t)$ 就会变。</p>
<div class="eq" data-no="(1.1.1)" data-title="牛顿第二定律：小车的数学模型">
$$u(t) = m\,\frac{\mathrm{d}^2 x(t)}{\mathrm{d}t^2}$$
<dl><dt>$u(t)$</dt><dd>推力，我们能施加的量，叫<b>输入</b>（input）</dd>
<dt>$x(t)$</dt><dd>小车位置，我们关心的量，叫<b>输出</b>（output）</dd>
<dt>$m$</dt><dd>小车质量</dd>
<dt>$\frac{\mathrm{d}^2x}{\mathrm{d}t^2}$</dt><dd>加速度（位置对时间求两次导数）</dd></dl>
<div class="say">推力先改变加速度，加速度累积成速度，速度累积成位置。所以输入要“经过一段时间”才体现在输出上——这就是“动态”。</div>
</div>

<h3>线性时不变（LTI）系统</h3>
<p>整套书（以及这套课程）研究的几乎都是 <span class="kw">LTI 系统</span>，它要同时满足两条性质：</p>
<ol>
<li><b>线性（叠加原理）</b>：如果输入 $u_1$ 产生输出 $x_1$、输入 $u_2$ 产生 $x_2$，那么输入 $a u_1 + b u_2$ 产生的输出就是 $a x_1 + b x_2$。也就是“力加倍，效果加倍；两个力一起推，效果相加”。</li>
<li><b>时不变</b>：今天推一下和明天推一下，效果一样，只是晚了一天。数学上表现为方程里的系数都是常数，不随时间变化。</li>
</ol>
<div class="tip">严格说，世界上没有真正“时不变”的系统（人不能两次踏进同一条河流）。但在我们关心的那段时间里参数变化很慢，就可以当成常数处理；非线性系统也常常可以在工作点附近“线性化”（见附录）。</div>

<div class="quiz"><div class="q">方程 $a\ddot x + b\dot x + \sin x = u$ 描述的是什么系统？</div>
<div class="opts"><button class="opt">线性时不变</button><button class="opt">线性时变</button><button class="opt" data-ok>非线性时不变</button></div>
<div class="why">系数 $a,b$ 都是常数，所以“时不变”；但 $\sin x_1 + \sin x_2 \ne \sin(x_1+x_2)$，不满足叠加原理，所以是“非线性”。若把 $\sin x$ 换成 $c(t)\,x$（系数随时间变），就成了线性时变系统。</div></div>

<h2><span class="sec">1.2</span>控制系统：开环与闭环</h2>
<p>知道了输入怎样影响输出，就可以反过来设计输入，让输出达到我们想要的目标。这个目标叫<span class="kw">参考值</span> $r(t)$（reference），负责计算输入的“大脑”叫<span class="kw">控制器</span>（controller）。控制器 + 被控的动态系统，合起来叫<span class="kw">控制系统</span>。</p>
<div data-fig="openloop"></div>
<div class="figcap">开环控制：控制器只看目标，不看结果</div>
<p><b>开环控制</b>：如果模型完全准确，我们可以根据 $r(t)$ 用式 (1.1.1) 反算出需要的力 $u(t)=m\ddot r(t)$，照着推就行。但只要模型有一点不准——比如半路有个东西掉进小车，质量变了——算好的力就不对了，而控制器对此一无所知。</p>
<div data-fig="closedloop"></div>
<div class="figcap">闭环控制：测量输出，与目标比较得到误差，控制器根据误差调整输入</div>
<p><b>闭环控制</b>：用传感器测量实际输出，<span class="kw">反馈</span>（feedback）回来与目标相减得到<span class="kw">误差</span> $e(t)=r(t)-x(t)$，控制器根据误差决定怎么推。偏左了就往右推，偏右了就往左推——它不需要知道扰动是什么，只要看到误差就能纠正。</p>

<div class="sim" data-id="openclosed"></div>

<div class="idea">开环像蒙着眼睛按事先背好的步数走路：路况一变就走偏；闭环像睁着眼睛走路：随时看、随时改。代价是闭环需要传感器，而且如果“改得太猛”还可能来回晃——后面很多章节都在研究“怎样改得又快又稳”。</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>动态系统：状态随时间变化；我们能施加的叫输入 $u$，关心的叫输出 $x$。</li>
<li>LTI 系统 = 满足叠加原理 + 系数不随时间变化。本课程研究的对象几乎都是 LTI 系统。</li>
<li>开环控制结构简单，但对模型误差和扰动无能为力；闭环控制通过反馈误差，自动补偿扰动与建模误差。</li>
</ol></div>`,
    sims: {
      openclosed(el) {
        const ui = card(el, { title: "开环 vs 闭环：半路掉进一个箱子", hint: "两辆质量 1 kg 的小车都要在 3 秒内平稳地从 0 移到 1 m。上面一辆按“事先算好的力”推（开环），下面一辆还会根据误差修正（闭环）。" });
        const sc = new Scene(ui.plotHost(), { height: 200 });
        const pl = new Plot(ui.plotHost(), { height: 220, x: [0, 8], y: [-0.2, 1.6], xlabel: "时间 t / s", ylabel: "位置 x / m" });
        const sDm = ui.slider({ label: "掉进箱子的质量 Δm (kg)", min: 0, max: 3, step: 0.1, value: 1.5, onInput: reset });
        const sW = ui.slider({ label: "持续逆风阻力 (N)", min: 0, max: 0.6, step: 0.05, value: 0, onInput: reset });
        const sK = ui.slider({ label: "闭环反馈增益 $K_p$", min: 2, max: 80, step: 1, value: 25, onInput: reset });
        ui.buttons([["▶ 重新开始", reset]]);
        const T = 3, m0 = 1, tDrop = 1.2;
        const ref = t => { if (t >= T) return [1, 0, 0]; const s = t / T; return [10 * s ** 3 - 15 * s ** 4 + 6 * s ** 5, (30 * s * s - 60 * s ** 3 + 30 * s ** 4) / T, (60 * s - 180 * s * s + 120 * s ** 3) / T / T]; };
        let st;
        function reset() { st = { t: 0, o: [0, 0], c: [0, 0], dropped: false, H: { t: [], o: [], c: [], r: [] } }; }
        reset();
        function stepSim(dt) {
          const dm = sDm.get(), w = sW.get(), Kp = sK.get(), Kd = 2 * Math.sqrt(Kp);
          const h = 0.002;
          for (let k = 0; k < dt / h && st.t < 8; k++) {
            const t = st.t, [r, rd, rdd] = ref(t);
            if (!st.dropped && t >= tDrop) { st.dropped = true; if (dm > 0) { st.o[1] *= m0 / (m0 + dm); st.c[1] *= m0 / (m0 + dm); } }
            const m = m0 + (st.dropped ? dm : 0);
            const uo = m0 * rdd;
            const uc = m0 * rdd + Kp * (r - st.c[0]) + Kd * (rd - st.c[1]);
            st.o = [st.o[0] + h * st.o[1], st.o[1] + h * (uo - w) / m];
            st.c = [st.c[0] + h * st.c[1], st.c[1] + h * (uc - w) / m];
            st.t += h;
          }
          const H = st.H; H.t.push(st.t); H.o.push(st.o[0]); H.c.push(st.c[0]); H.r.push(ref(st.t)[0]);
        }
        function drawCart(c, X, Y, w, dm, show, color) {
          c.fillStyle = color; c.beginPath(); c.roundRect(X - w / 2, Y - 26, w, 22, 5); c.fill();
          c.fillStyle = COL.ink; [X - w / 3, X + w / 3].forEach(x => { c.beginPath(); c.arc(x, Y - 3, 5, 0, 7); c.fill(); });
          if (show && dm > 0) { const s = 10 + dm * 6; c.fillStyle = "#a0723c"; c.fillRect(X - s / 2, Y - 26 - s, s, s); c.strokeStyle = "#6b4a22"; c.strokeRect(X - s / 2, Y - 26 - s, s, s); }
        }
        sc.render = () => {
          const c = sc.ctx, W = sc.W; sc.clear();
          const x2p = x => 50 + (x + 0.2) / 1.8 * (W - 90);
          [["开环", st.o[0], 80, COL.c2], ["闭环", st.c[0], 175, COL.c1]].forEach(([lab, x, Y, col]) => {
            c.strokeStyle = COL.axis; c.lineWidth = 2; c.beginPath(); c.moveTo(30, Y); c.lineTo(W - 20, Y); c.stroke();
            c.fillStyle = COL.ink2; c.font = "600 13px 'Microsoft YaHei',sans-serif"; c.textAlign = "left"; c.fillText(lab, 6, Y - 30);
            // 目标旗
            const fx = x2p(1); c.strokeStyle = COL.c3; c.lineWidth = 2; c.beginPath(); c.moveTo(fx, Y); c.lineTo(fx, Y - 60); c.stroke();
            c.fillStyle = COL.c3; c.beginPath(); c.moveTo(fx, Y - 60); c.lineTo(fx + 18, Y - 54); c.lineTo(fx, Y - 48); c.fill();
            const X = x2p(x), dm = sDm.get();
            // 掉落中的箱子
            if (!st.dropped && st.t > tDrop - 0.4 && dm > 0) { const s = 10 + dm * 6, yy = Y - 26 - s - (tDrop - st.t) / 0.4 * 50; c.fillStyle = "#a0723c"; c.fillRect(X - s / 2, yy, s, s); }
            drawCart(c, X, Y, 54, dm, st.dropped, col);
          });
          c.fillStyle = COL.muted; c.font = "12px Consolas"; c.textAlign = "right"; c.fillText("t = " + st.t.toFixed(2) + " s", W - 10, 18);
        };
        pl.render = () => {
          pl.clear(); const H = st.H;
          pl.vline(tDrop, COL.muted, [3, 4]); pl.text(tDrop + 0.08, -0.12, "箱子掉入", COL.muted);
          pl.line(H.t, H.r, COL.c3, 1.6, [6, 4]); pl.line(H.t, H.o, COL.c2, 2.2); pl.line(H.t, H.c, COL.c1, 2.2);
          pl.legend([{ label: "目标 r(t)", color: COL.c3, dash: [6, 4] }, { label: "开环", color: COL.c2 }, { label: "闭环", color: COL.c1 }]);
        };
        CE.loop(ui.root, dt => {
          if (st.t < 8) stepSim(dt);
          sc.render(); pl.render();
          const eo = 1 - st.o[0], ec = 1 - st.c[0];
          ui.out(`t = <b>${st.t.toFixed(1)} s</b><br>开环误差：<span class="${Math.abs(eo) > 0.05 ? "bad" : "good"}">${f2(eo, 3)} m</span><br>闭环误差：<span class="${Math.abs(ec) > 0.05 ? "bad" : "good"}">${f2(ec, 3)} m</span>`);
        });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>把 Δm 和逆风都调到 0：两辆车都准确到达——模型准确时开环也能完美工作。</li><li>Δm = 1.5 kg：开环车被“拖慢”后永远到不了终点；闭环车短暂落后后追了回来。</li><li>加一点逆风：开环车甚至会倒退；闭环车停在离终点很近（但不完全为零）的地方——这个小误差留到第 7 章用“积分控制”消除。</li><li>把 $K_p$ 调小：闭环纠错变慢、误差变大；调大：纠得更快更准。</li></ol></div>`);
        CE.renderMath(ui.root);
      },
    },
  });

  /* ════════════════ 2.1.1 卷积 ════════════════ */
  // 欠阻尼弹簧的冲激响应：m=1,b=0.6,k=4
  const wn = 2, zeta = 0.15, wd = wn * Math.sqrt(1 - zeta * zeta);
  const hImp = t => t < 0 ? 0 : Math.exp(-zeta * wn * t) * Math.sin(wd * t) / wd;
  const inputs = {
    rise: { name: "平滑上升", f: t => t < 0 ? 0 : t < 5 ? (1 - Math.cos(Math.PI * t / 5)) / 2 : 1 },
    pulse: { name: "短促一推", f: t => (t >= 0 && t < 0.4) ? 2.5 : 0 },
    two: { name: "推两下", f: t => (t >= 0 && t < 0.4) || (t >= 3 && t < 3.4) ? 2.5 : 0 },
    sine: { name: "正弦", f: t => t < 0 ? 0 : 0.5 * Math.sin(1.1 * t) },
    tri: { name: "三角形", f: t => t < 0 ? 0 : t < 2 ? t / 2 : t < 4 ? (4 - t) / 2 : 0 },
  };
  Course.add({
    id: "v1-conv", vol: 1, ord: 20, no: "2.1", chName: "第2章 动态系统建模——传递函数", ref: "卷1 §2.1.1",
    title: "卷积：涟漪是怎样叠加的",
    sub: "LTI 系统的输出 = 过去所有输入各自引起的“涟漪”之和。理解了这句话，就理解了卷积。",
    goals: ["理解单位冲激函数 δ(t) 和冲激响应 h(t)", "把卷积公式拆成“切块 → 放大平移 → 求和”三步", "明白为什么说冲激响应包含了 LTI 系统的全部信息"],
    minutes: 20,
    html: String.raw`
<div class="story">往平静的水面扔一颗石子，会荡起一圈圈涟漪；涟漪还没消失，又扔下第二颗，两次的涟漪就叠在一起。扔石子是<b>输入</b>，涟漪是<b>输出</b>；此刻的水面，是之前每一颗石子留下的影响之和。把这种“过去的输入在现在的叠加”写成数学，就是<b>卷积</b>。</div>

<h2><span class="sec">2.1.1</span>从一个弹簧说起</h2>
<p>想象一根竖直的弹簧，顶端用手向上推一下（力很大、时间很短），松手后弹簧会上下振荡、逐渐停下。弹簧力与形变成正比、任何时候推效果都一样，所以它是 LTI 系统。</p>
<p>问题来了：如果手上的力不是“推一下”，而是一条连续变化的曲线 $u(t)$，输出 $x(t)$ 是什么样子？书中的思路非常巧妙：</p>
<ol>
<li><b>切块</b>：把连续的输入切成很多宽为 $\Delta T$ 的窄矩形，第 $i$ 块的高度是 $u(i\Delta T)$。</li>
<li><b>每块单独算</b>：每一块都是“短促的一推”，因为 LTI，它们引起的响应<b>形状完全相同</b>，只是幅度不同（按面积放大）、起始时间不同（延迟 $i\Delta T$）。</li>
<li><b>求和</b>：所有块的响应加起来就是总输出。</li>
</ol>

<h3>单位冲激函数 δ(t)：“面积为 1 的瞬间一推”</h3>
<p>为了描述“短促的一推”，引入一个宽度为 $\Delta T$、高度为 $1/\Delta T$ 的方波 $\delta_\Delta(t)$——不管多窄，面积永远是 1。把它越压越窄，极限就是<span class="kw">单位冲激函数</span>：</p>
<div class="eq" data-no="(2.1.2)" data-title="单位冲激函数（狄拉克函数）">
$$\delta(t)=0\;(t\neq 0),\qquad \int_{-\infty}^{\infty}\delta(t)\,\mathrm{d}t = 1,\qquad \delta(t)=\lim_{\Delta T\to 0}\delta_\Delta(t)$$
<dl><dt>$\delta(t)$</dt><dd>宽度为 0、面积为 1 的“尖刺”，是纯数学的理想化</dd><dt>$h(t)$</dt><dd>系统对 $\delta(t)$ 的响应，叫<b>冲激响应</b></dd></dl>
</div>
<div class="sim" data-id="squash"></div>

<h3>把三步写成公式</h3>
<p>第 $i$ 块的面积是 $u(i\Delta T)\Delta T$，可以看成“$u(i\Delta T)\Delta T$ 倍的单位方波，延迟 $i\Delta T$”，所以它引起的响应是 $u(i\Delta T)\Delta T\, h_\Delta(t-i\Delta T)$。全部加起来：</p>
<div class="eq" data-no="(2.1.9)" data-title="卷积的离散形式">
$$x(t)=\sum_{i=0}^{n} \underbrace{u(i\Delta T)\,\Delta T}_{\text{第 i 块的“面积”}}\;\underbrace{h_\Delta(t-i\Delta T)}_{\text{延迟 }i\Delta T\text{ 的标准响应}}$$
</div>
<p>让 $\Delta T\to 0$：求和变积分，$i\Delta T$ 变成连续的过去时刻 $\tau$，$\Delta T$ 变成 $\mathrm{d}\tau$：</p>
<div class="eq" data-no="(2.1.10)" data-title="卷积积分">
$$x(t)=\int_0^t u(\tau)\,h(t-\tau)\,\mathrm{d}\tau \;=\; u(t) * h(t)$$
<dl><dt>$\tau$</dt><dd>过去的某个时刻（“石子落水的时刻”）</dd>
<dt>$u(\tau)\,\mathrm{d}\tau$</dt><dd>那一瞬间输入的“量”</dd>
<dt>$h(t-\tau)$</dt><dd>那一瞬间的输入，经过 $t-\tau$ 时间后还剩多少影响</dd>
<dt>$*$</dt><dd>卷积运算符</dd></dl>
<div class="say">现在的输出 = 对所有过去时刻 τ，把“那时的输入”乘以“它至今残留的影响”，再全部加起来。</div>
</div>

<div class="sim" data-id="conv"></div>

<div class="idea">卷积公式里的 $h(t-\tau)$ 把冲激响应“翻转并平移”到了当前时刻 $t$：越久远的输入，对应 $h$ 曲线越靠后的部分（影响已经衰减）；刚刚发生的输入，对应 $h$ 曲线的开头。上面交互图里紫色的曲线就是这个“翻转后的记忆权重”。</div>

<p><b>核心结论</b>：只要知道系统的冲激响应 $h(t)$，任何输入下的输出都能用卷积算出来。所以说<b>冲激响应包含了 LTI 系统的全部特性</b>，它是系统的“指纹”。</p>
<div class="tip" data-t="趣味实验：录下一个房间的“指纹”">在空旷的礼堂里扎破一个气球（短促、能量大的冲激），用麦克风录下回声，就得到了这个礼堂的冲激响应。把任何一段干声和它做卷积，听起来就像在礼堂里演奏——音频软件里的“卷积混响”正是这么做的。</div>

<div class="summary"><h3>本课要点</h3><ol>
<li>单位冲激 δ(t)：宽度趋于 0、面积为 1；系统对它的响应叫冲激响应 h(t)。</li>
<li>LTI 系统的输出与输入是卷积关系：$x(t)=\int_0^t u(\tau)h(t-\tau)\mathrm{d}\tau$。</li>
<li>卷积 = 切块、每块按面积放大并延迟、求和；冲激响应完整描述了 LTI 系统。</li>
</ol></div>`,
    sims: {
      squash(el) {
        const ui = card(el, { title: "把方框压扁：面积不变，响应趋于同一条曲线", hint: "拖动滑块让方波变窄变高（面积始终为 1）。右边的响应很快就和真正的冲激响应 h(t)（虚线）重合。" });
        const [a, b] = ui.row("", "");
        const p1 = new Plot(a, { height: 220, x: [-0.2, 3], y: [0, 6], xlabel: "t", title: "输入 δΔ(t)" });
        const p2 = new Plot(b, { height: 220, x: [0, 12], y: [-0.45, 0.6], xlabel: "t", title: "输出 hΔ(t)" });
        const s = ui.slider({ label: "方框宽度 ΔT", min: 0.2, max: 2.5, step: 0.05, value: 2, onInput: draw });
        function hDelta(t, D) { // (1/D)∫_0^{min(t,D)} h(t-τ)dτ
          const n = 60, up = Math.min(t, D); if (up <= 0) return 0; let s = 0; for (let i = 0; i < n; i++) { const tau = (i + 0.5) / n * up; s += hImp(t - tau); } return s * up / n / D;
        }
        function draw() {
          const D = s.get();
          p1.render = () => { p1.clear(); p1.rect(0, 0, D, 1 / D, COL.c2 + "55", COL.c2); p1.text(D / 2, 1 / D + 0.25, "面积 = 1", COL.ink2, "center"); p1.text(D + 0.1, 1 / D, "高 1/ΔT = " + f2(1 / D), COL.muted); };
          p2.render = () => { p2.clear(); const ts = linspace(0, 12, 300); p2.line(ts, ts.map(hImp), COL.muted, 1.6, [6, 4]); p2.line(ts, ts.map(t => hDelta(t, D)), COL.c1, 2.4); p2.legend([{ label: "hΔ(t)", color: COL.c1 }, { label: "h(t) 极限", color: COL.muted, dash: [6, 4] }], "tr"); };
          p1.render(); p2.render();
        }
        draw();
      },
      conv(el) {
        const ui = card(el, { title: "亲手做一次卷积：切块 → 叠加", hint: "上图：输入被切成矩形块；紫线是当前时刻 t* 的“记忆权重” h(t*−τ)。下图：每块引起的响应（细线）与它们的和（粗线）。<b>拖动下图中的竖线</b>改变 t*。" });
        const p1 = new Plot(ui.plotHost(), { height: 200, x: [0, 14], y: [-0.6, 1.4], xlabel: "τ（过去的时刻）", title: "输入 u(τ) 与翻转平移后的冲激响应" });
        const p2 = new Plot(ui.plotHost(), { height: 240, x: [0, 14], y: [-1, 1.5], xlabel: "t", title: "输出 x(t)" });
        let kind = "rise", tStar = 6;
        ui.seg({ label: "输入波形", options: Object.entries(inputs).map(([k, v]) => [k, v.name]), value: kind, onChange: v => { kind = v; draw(); } });
        const sD = ui.slider({ label: "切块宽度 ΔT", min: 0.1, max: 1.5, step: 0.05, value: 0.8, onInput: draw });
        const cTh = ui.check({ label: "显示每一块的响应（细线）", value: true, onChange: draw });
        const cEx = ui.check({ label: "显示精确卷积（ΔT→0）", value: true, onChange: draw });
        const tg = linspace(0, 14, 701), dt = tg[1] - tg[0];
        function exact(u) { return tg.map(t => { let s = 0; for (let tau = 0; tau <= t; tau += dt) s += u(tau) * hImp(t - tau) * dt; return s; }); }
        let cache = {};
        function draw() {
          const u = inputs[kind].f, D = sD.get();
          if (!cache[kind]) cache[kind] = exact(u);
          const ex = cache[kind];
          const nb = Math.ceil(14 / D), blocks = [];
          for (let i = 0; i < nb; i++) { const ui_ = u(i * D); if (Math.abs(ui_) > 1e-9) blocks.push([i * D, ui_]); }
          const resp = blocks.map(([t0, h0]) => tg.map(t => h0 * D * hImp(t - t0)));
          const sum = tg.map((_, j) => resp.reduce((s, r) => s + r[j], 0));
          let ymx = Math.max(...ex.map(Math.abs), ...sum.map(Math.abs), 0.3) * 1.25;
          p2.setRange(null, [-ymx, ymx]);
          const umx = Math.max(...tg.map(t => Math.abs(u(t))), 0.5);
          p1.setRange(null, [Math.min(-0.3, -umx * 0.6), umx * 1.35]);
          p1.render = () => {
            p1.clear();
            blocks.forEach(([t0, h0]) => p1.rect(t0, Math.min(0, h0), t0 + D, Math.max(0, h0), COL.c2 + "40", COL.c2));
            p1.line(tg, tg.map(u), COL.c2, 1.6);
            // 翻转平移的冲激响应：h(t*-τ)，缩放显示
            const sc = umx / 0.5;
            const taus = tg.filter(t => t <= tStar);
            p1.line(taus, taus.map(t => hImp(tStar - t) * sc), COL.c5, 2.2);
            p1.vline(tStar, COL.c5, [4, 3]);
            p1.legend([{ label: "u(τ)", color: COL.c2 }, { label: "h(t*−τ)（缩放）", color: COL.c5 }], "tr");
          };
          p2.render = () => {
            p2.clear();
            if (cTh.get()) resp.forEach((r, k) => p2.line(tg, r, `hsl(${200 + k * 137 % 120},60%,60%)`, 1, null));
            if (cEx.get()) p2.line(tg, ex, COL.muted, 1.8, [6, 4]);
            p2.line(tg, sum, COL.c1, 2.6);
            p2.vline(tStar, COL.c5, [4, 3], 1.6);
            const j = Math.round(tStar / dt); p2.dot(tStar, sum[j], COL.c5, 5);
            p2.legend([{ label: "Σ 各块响应", color: COL.c1 }, ...(cEx.get() ? [{ label: "精确卷积", color: COL.muted, dash: [6, 4] }] : [])], "tr");
          };
          p1.render(); p2.render();
          const j = Math.round(tStar / dt);
          ui.out(`块数：<b>${blocks.length}</b><br>t* = <b>${f2(tStar)}</b> 时<br>近似输出 x = <b class="big">${f2(sum[j], 3)}</b><br>精确值 = ${f2(ex[j], 3)}<br><span class="sim-note">ΔT 越小，粗线越贴近虚线——这就是 (2.1.9)→(2.1.10) 的极限过程。</span>`);
        }
        p2.drag({ down: (x) => { tStar = clamp(x, 0, 14); draw(); }, move: x => { tStar = clamp(x, 0, 14); draw(); }, hover: () => true });
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>选“短促一推”：输出就是一条冲激响应——系统的“指纹”。</li><li>选“推两下”：第二次推力的响应叠加在第一次的余波上，就像两颗石子的涟漪。</li><li>选“平滑上升”，把 ΔT 从 1.5 拖到 0.1，看近似结果怎样逼近精确卷积。</li><li>拖动竖线 t*：紫色曲线告诉你，越久远的输入（左边）权重越小——系统在慢慢“遗忘”。</li></ol></div>`);
      },
    },
  });

  /* ════════════════ 2.1.2 微分方程建模 ════════════════ */
  Course.add({
    id: "v1-ode", vol: 1, ord: 30, no: "2.1", chName: "第2章 动态系统建模——传递函数", ref: "卷1 §2.1.2",
    title: "用微分方程描述系统：体重、电路与弹簧",
    sub: "建模的套路：分析物理规律 → 列方程 → 消去中间变量 → 整理成标准形式。",
    goals: ["会用“流入−流出”的思路建立体重的一阶模型", "会对弹簧–质量–阻尼系统做受力分析并写出二阶方程", "体会到不同物理系统可以有相同形式的方程"],
    minutes: 20,
    html: String.raw`
<p>上一课说卷积描述了输入与输出的关系，但直接算卷积很麻烦。实际上，<b>微分方程</b>就能直接刻画这种关系。建立微分方程的通用流程是：分析物理特性 → 列出方程组 → 消去中间变量 → 写成标准形式。下面看三个经典例子。</p>

<h2><span class="sec">例 2.1.1</span>体重变化：一个“流入−流出”的系统</h2>
<div class="story">体重的变化取决于热量的进与出：吃饭带来热量摄入 $E_i$，呼吸、日常活动和运动消耗热量 $E_e$。经验上，大约 <b>7000 千卡 ≈ 1 千克脂肪</b>：多吃 7000 千卡长一公斤，多烧 7000 千卡减一公斤。</div>
<div class="eq" data-no="(2.1.11)" data-title="体重的变化率 = 净摄入 / 7000">
$$\frac{\mathrm{d}m}{\mathrm{d}t}=\frac{E_i-E_e}{7000},\qquad E_e = E_a + \alpha P,\qquad P = 10m + 6.25h - 5a + S$$
<dl><dt>$m$</dt><dd>体重（kg），$t$ 以“天”为单位</dd>
<dt>$E_i,\;E_a$</dt><dd>每日摄入热量、每日额外运动消耗（kCal）</dd>
<dt>$\alpha$</dt><dd>劳动强度系数：轻体力 1.3、中体力 1.5、重体力约 1.9</dd>
<dt>$P$</dt><dd>基础代谢率 BMR（Mifflin–St Jeor 公式，误差约 5%）</dd>
<dt>$h,\;a,\;S$</dt><dd>身高（cm）、年龄、性别修正（男 $S=5$，女 $S=-161$）</dd></dl>
</div>
<p>把 $E_e$ 和 $P$ 代进去、把含 $m$ 的项移到左边，就得到一个<b>一阶线性微分方程</b>：</p>
<div class="eq" data-no="(2.1.14)">
$$7000\,\frac{\mathrm{d}m}{\mathrm{d}t} + 10\alpha\, m = E_i - E_a - \alpha\,(6.25h - 5a + S)$$
<div class="say">体重越大，基础代谢越高，消耗越多——体重自己会“刹车”。所以体重不会无限地涨或跌，而是慢慢趋近一个平衡值。</div>
</div>
<div class="sim" data-id="weight"></div>
<div class="tip">从方程可以直接读出两个关键量（本页补充推导）：令 $\dot m = 0$ 得稳态体重 $m_\infty = \dfrac{E_i-E_a-\alpha(6.25h-5a+S)}{10\alpha}$；时间常数 $\tau = \dfrac{7000}{10\alpha}=\dfrac{700}{\alpha}$ 天，轻体力者约 538 天！这解释了为什么减肥是一场持久战。第 7 章会用控制的方法再来“管理”体重。</div>

<h2><span class="sec">例 2.1.2</span>RLC 电路</h2>
<p>电源 $e(t)$ 与电感 $L$、电容 $C$、电阻 $R$ 串联，电流为 $i(t)$。<b>基尔霍夫电压定律</b>说：沿闭合回路，电压降之和等于电源电动势。三个元件的电压分别是 $L\frac{\mathrm{d}i}{\mathrm{d}t}$、$\frac{1}{C}\int i\,\mathrm{d}t$、$iR$，相加等于 $e(t)$，两边求导消去积分：</p>
<div class="eq" data-no="(2.1.18)">
$$L\frac{\mathrm{d}^2 i}{\mathrm{d}t^2} + R\frac{\mathrm{d}i}{\mathrm{d}t} + \frac{1}{C}\,i = \frac{\mathrm{d}e}{\mathrm{d}t}$$
</div>

<h2><span class="sec">例 2.1.3</span>弹簧–质量–阻尼系统</h2>
<p>这是整个控制理论里出镜率最高的模型：绝大多数振动问题都能简化成它，把它竖起来就是汽车的“四分之一悬挂模型”。质量块 $m$ 受外力 $f(t)$，弹簧把它往回拉（胡克定律 $f_k=-kx$），阻尼器阻碍它运动（$f_b=-b\dot x$，就像手在水里划，越快阻力越大）。牛顿第二定律 $m\ddot x = f + f_k + f_b$ 整理得：</p>
<div class="eq" data-no="(2.1.21)" data-title="弹簧–质量–阻尼系统">
$$m\,\ddot x(t) + b\,\dot x(t) + k\,x(t) = f(t)$$
<dl><dt>$m$</dt><dd>质量：惯性，越大越“推不动”</dd><dt>$b$</dt><dd>阻尼系数：耗散能量，让振动停下来</dd><dt>$k$</dt><dd>弹性系数：越大越“硬”，振得越快</dd><dt>$f(t)$</dt><dd>外力（输入）；$x(t)$ 位移（输出），向右为正</dd></dl>
</div>
<div class="sim" data-id="msd"></div>
<p>对比 RLC 电路的方程，你会发现它们的形式一模一样：电感像质量（阻碍变化），电阻像阻尼（耗散能量），电容像弹簧（储存能量）。<b>不同的物理系统，数学上可以是同一个东西</b>——这正是控制理论能“一招通吃”的原因。</p>
<table class="tbl"><tr><th>机械系统</th><th>电路系统</th><th>作用</th></tr>
<tr><td>质量 $m$</td><td>电感 $L$</td><td>惯性，阻碍速度（电流）变化</td></tr>
<tr><td>阻尼 $b$</td><td>电阻 $R$</td><td>耗散能量</td></tr>
<tr><td>弹性系数 $k$</td><td>电容倒数 $1/C$</td><td>储存能量</td></tr>
<tr><td>位移 $x$</td><td>电荷 $q=\int i\,\mathrm{d}t$</td><td>状态</td></tr></table>

<div class="summary"><h3>本课要点</h3><ol>
<li>建模流程：物理规律 → 方程组 → 消去中间变量 → 标准形式。</li>
<li>体重模型是一阶系统，会指数地趋近一个由饮食、运动、身体条件共同决定的稳态值。</li>
<li>弹簧–质量–阻尼 $m\ddot x+b\dot x+kx=f$ 是二阶系统的代表；RLC 电路与它形式相同。</li>
</ol></div>`,
    sims: {
      weight(el) {
        const ui = card(el, { title: "体重模拟器", hint: "调整饮食与运动，看两年内体重怎样变化。注意曲线的形状：先快后慢，趋近一条水平线。" });
        const pl = new Plot(ui.plotHost(), { height: 280, x: [0, 730], y: [50, 100], xlabel: "天", ylabel: "体重 kg" });
        let sex = 5, alpha = 1.3;
        ui.seg({ label: "性别", options: [[5, "男 S=5"], [-161, "女 S=−161"]], value: 5, onChange: v => { sex = v; draw(); } });
        ui.seg({ label: "劳动强度 α", options: [[1.3, "轻 1.3"], [1.5, "中 1.5"], [1.9, "重 1.9"]], value: 1.3, onChange: v => { alpha = v; draw(); } });
        const sEi = ui.slider({ label: "每日摄入 $E_i$ (kCal)", min: 1200, max: 3800, step: 50, value: 2400, fmt: v => v.toFixed(0), onInput: draw });
        const sEa = ui.slider({ label: "每日运动 $E_a$ (kCal)", min: 0, max: 1000, step: 25, value: 200, fmt: v => v.toFixed(0), onInput: draw });
        const sM0 = ui.slider({ label: "初始体重 (kg)", min: 45, max: 120, step: 1, value: 85, fmt: v => v.toFixed(0), onInput: draw });
        const sH = ui.slider({ label: "身高 (cm)", min: 150, max: 200, step: 1, value: 175, fmt: v => v.toFixed(0), onInput: draw });
        const sA = ui.slider({ label: "年龄", min: 18, max: 70, step: 1, value: 30, fmt: v => v.toFixed(0), onInput: draw });
        function draw() {
          const Ei = sEi.get(), Ea = sEa.get(), h = sH.get(), a = sA.get(), m0 = sM0.get();
          const rhs = Ei - Ea - alpha * (6.25 * h - 5 * a + sex), minf = rhs / (10 * alpha), tau = 700 / alpha;
          const ts = linspace(0, 730, 300), ms = ts.map(t => minf + (m0 - minf) * Math.exp(-t / tau));
          pl.setRange(null, [Math.min(m0, minf) - 6, Math.max(m0, minf) + 6]);
          pl.render = () => { pl.clear(); pl.hline(minf, COL.c3, [6, 4], 1.6); pl.text(10, minf + (m0 > minf ? -1.6 : 1.0), "稳态 " + f2(minf, 1) + " kg", COL.c3); pl.line(ts, ms, COL.c1, 2.6); pl.vline(tau, COL.muted, [3, 4]); pl.text(tau + 6, Math.max(m0, minf) + 4, "τ = " + tau.toFixed(0) + " 天", COL.muted); };
          pl.render();
          const P0 = 10 * m0 + 6.25 * h - 5 * a + sex;
          ui.out(`当前 BMR：<b>${P0.toFixed(0)}</b> kCal/天<br>当前每日净摄入：<b>${(Ei - Ea - alpha * P0).toFixed(0)}</b> kCal<br>一年后：<b class="big">${f2(minf + (m0 - minf) * Math.exp(-365 / tau), 1)} kg</b><br>稳态：${f2(minf, 1)} kg　时间常数：${tau.toFixed(0)} 天`);
        }
        draw();
      },
      msd(el) {
        const ui = card(el, { title: "弹簧–质量–阻尼系统：调参数，看它怎么动", hint: "点“踢一脚”给质量块一个冲击；或选择持续外力。观察 m、b、k 各自对运动的影响。" });
        const sc = new Scene(ui.plotHost(), { height: 130 });
        const pl = new Plot(ui.plotHost(), { height: 220, x: [0, 15], y: [-1.5, 1.5], xlabel: "t / s", ylabel: "x(t)" });
        const sm = ui.slider({ label: "质量 m", min: 0.2, max: 5, step: 0.1, value: 1 });
        const sb = ui.slider({ label: "阻尼 b", min: 0, max: 6, step: 0.05, value: 0.5 });
        const sk = ui.slider({ label: "弹性系数 k", min: 0.2, max: 10, step: 0.1, value: 2 });
        let force = 0;
        ui.seg({ label: "外力 f(t)", options: [[0, "无"], [1, "恒力 f=1"], [2, "正弦 sin(t)"]], value: 0, onChange: v => { force = v; } });
        let x = [1, 0], t = 0, H = { t: [], x: [] };
        const restart = (x0 = [1, 0]) => { x = x0; t = 0; H = { t: [], x: [] }; };
        ui.buttons([["👟 踢一脚", () => { x = [x[0], x[1] + 2.5]; }], ["↺ 从 x=1 释放", () => restart([1, 0])]]);
        const fF = tt => force === 1 ? 1 : force === 2 ? Math.sin(tt) : 0;
        CE.loop(ui.root, dt => {
          const m = sm.get(), b = sb.get(), k = sk.get();
          for (let i = 0; i < 10; i++) { x = rk4((s, tt) => [s[1], (fF(tt) - b * s[1] - k * s[0]) / m], x, t, dt / 10); t += dt / 10; }
          H.t.push(t); H.x.push(x[0]);
          if (t > 15) { H.t = H.t.map(v => v - dt); H.t.shift(); H.x.shift(); t -= dt; }
          // 场景
          const c = sc.ctx, W = sc.W; sc.clear();
          const wallX = 30, cx = W / 2 + x[0] * W * 0.18, Y = 70, bw = 60;
          c.fillStyle = COL.axis; c.fillRect(wallX - 12, 20, 12, 100); c.fillRect(wallX - 12, Y + 32, W - 30, 4);
          // 弹簧（上）
          c.strokeStyle = COL.c1; c.lineWidth = 2; c.beginPath(); const y1 = Y - 12, L = cx - bw / 2 - wallX;
          c.moveTo(wallX, y1); for (let i = 1; i <= 12; i++) c.lineTo(wallX + L * i / 13, y1 + (i % 2 ? -9 : 9)); c.lineTo(cx - bw / 2, y1); c.stroke();
          // 阻尼（下）
          const y2 = Y + 14, mid = wallX + L * 0.55; c.strokeStyle = COL.c5; c.beginPath(); c.moveTo(wallX, y2); c.lineTo(mid - 10, y2); c.moveTo(mid - 10, y2 - 9); c.lineTo(mid + 14, y2 - 9); c.moveTo(mid - 10, y2 + 9); c.lineTo(mid + 14, y2 + 9); c.moveTo(mid - 10, y2 - 9); c.lineTo(mid - 10, y2 + 9); c.moveTo(mid, y2 - 6); c.lineTo(mid, y2 + 6); c.moveTo(mid, y2); c.lineTo(cx - bw / 2, y2); c.stroke();
          c.fillStyle = COL.c2; c.beginPath(); c.roundRect(cx - bw / 2, Y - 30, bw, 60, 6); c.fill();
          c.fillStyle = "#fff"; c.font = "600 15px serif"; c.textAlign = "center"; c.fillText("m", cx, Y + 5);
          c.strokeStyle = COL.muted; c.setLineDash([3, 3]); c.beginPath(); c.moveTo(W / 2, 10); c.lineTo(W / 2, Y + 32); c.stroke(); c.setLineDash([]);
          if (force) { const f = fF(t); c.strokeStyle = COL.c4; c.fillStyle = COL.c4; c.lineWidth = 3; const ax = cx + bw / 2, L2 = f * 40; c.beginPath(); c.moveTo(ax, Y); c.lineTo(ax + L2, Y); c.stroke(); }
          pl.render = () => { pl.clear(); pl.line(H.t, H.x, COL.c1, 2.4); };
          pl.render();
          const wn = Math.sqrt(k / m), z = b / (2 * Math.sqrt(k * m));
          ui.out(`方程：${f2(m, 1)}ẍ + ${f2(b, 2)}ẋ + ${f2(k, 1)}x = f(t)<br>固有频率 $\\omega_n=\\sqrt{k/m}$ = <b>${f2(wn)}</b> rad/s<br>阻尼比 $\\zeta=\\frac{b}{2\\sqrt{km}}$ = <b>${f2(z)}</b> ${z < 1 ? "（欠阻尼，会振荡）" : z === 1 ? "（临界阻尼）" : "（过阻尼，不振荡）"}`);
        });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>把 b 调成 0：振动永不停止（没有能量耗散）。</li><li>增大 k：振得更快；增大 m：振得更慢、更“笨重”。</li><li>逐渐增大 b，找到“刚好不振荡”的位置——这就是第 5 章的“临界阻尼”。</li></ol></div>`);
      },
    },
  });

  /* ════════════════ 2.2–2.4 拉普拉斯变换与传递函数 ════════════════ */
  Course.add({
    id: "v1-laplace", vol: 1, ord: 40, no: "2.2–2.4", chName: "第2章 动态系统建模——传递函数", ref: "卷1 §2.2–2.5",
    title: "拉普拉斯变换与传递函数：把微积分变成代数",
    sub: "拉普拉斯变换让“求导”变成“乘 s”、让“卷积”变成“乘法”。极点的位置一眼告诉你系统会怎么动。",
    goals: ["知道拉普拉斯变换的定义和最常用的几条性质", "会用部分分式做逆变换，理解“分母的根出现在指数上”", "理解传递函数、极点、开环与闭环传递函数"],
    pre: "复数 $s=\\sigma+j\\omega$、欧拉公式 $e^{j\\theta}=\\cos\\theta+j\\sin\\theta$",
    minutes: 30,
    html: String.raw`
<p>上一课的 RL 电路方程 $u = L\dot x + Rx$，如果直接解微分方程或者算卷积，都很繁琐，系统一复杂更是如此。经典控制理论的法宝是<span class="kw">拉普拉斯变换</span>：把时间域里的函数 $f(t)$ 搬到复数域里变成 $F(s)$，在那里<b>微分方程变成代数方程，卷积变成乘法</b>。</p>

<h2><span class="sec">2.2.1</span>定义</h2>
<div class="eq" data-no="(2.2.4)" data-title="拉普拉斯变换">
$$\mathcal{L}[f(t)] = F(s) = \int_0^{\infty} f(t)\,e^{-st}\,\mathrm{d}t,\qquad s=\sigma+j\omega$$
<dl><dt>$f(t)$</dt><dd>时间域里的函数（从 $t=0$ 开始，之前的事控制工程不关心）</dd>
<dt>$s$</dt><dd>复数“频率”：实部 $\sigma$ 管衰减/增长，虚部 $\omega$ 管振荡</dd>
<dt>$F(s)$</dt><dd>变换后的函数，叫“象函数”</dd></dl>
<div class="say">用一族 $e^{-st}$ 去“探测”f(t)。当 $\sigma=0$ 时它就是傅里叶变换，所以傅里叶变换是拉普拉斯变换的特例。</div>
</div>
<p>最重要的一个例子：</p>
<div class="eq" data-no="(2.2.5)">
$$\mathcal{L}[e^{-at}] = \int_0^\infty e^{-(s+a)t}\mathrm{d}t = \left.-\frac{1}{s+a}e^{-(s+a)t}\right|_0^\infty = \frac{1}{s+a}$$
</div>
<p>由此和欧拉公式可以推出整张常用变换表。最关键的几条性质：</p>
<table class="tbl">
<tr><th>性质 / 函数</th><th>时间域 $f(t)$</th><th>s 域 $F(s)$</th><th>意义</th></tr>
<tr><td>线性</td><td>$af+bg$</td><td>$aF+bG$</td><td>满足叠加原理</td></tr>
<tr><td>微分</td><td>$\dot f(t)$</td><td>$sF(s)-f(0)$</td><td><b>求导 = 乘 s</b>（再减初值）</td></tr>
<tr><td>卷积</td><td>$f*g$</td><td>$F(s)G(s)$</td><td><b>卷积 = 相乘</b></td></tr>
<tr><td>冲激</td><td>$\delta(t)$</td><td>$1$</td><td></td></tr>
<tr><td>阶跃</td><td>$1$</td><td>$1/s$</td><td></td></tr>
<tr><td>指数</td><td>$e^{-at}$</td><td>$\dfrac{1}{s+a}$</td><td>极点在 $-a$</td></tr>
<tr><td>正弦 / 余弦</td><td>$\sin at$ / $\cos at$</td><td>$\dfrac{a}{s^2+a^2}$ / $\dfrac{s}{s^2+a^2}$</td><td>极点在 $\pm ja$</td></tr>
<tr><td>衰减振荡</td><td>$e^{at}\sin bt$ / $e^{at}\cos bt$</td><td>$\dfrac{b}{(s-a)^2+b^2}$ / $\dfrac{s-a}{(s-a)^2+b^2}$</td><td>极点在 $a\pm jb$</td></tr>
</table>
<details class="derive"><summary>为什么“求导 = 乘 s”、“卷积 = 相乘”？</summary><div class="body">
<div class="step"><b>微分性质</b>（分部积分）：$\displaystyle\int_0^\infty \dot f e^{-st}\mathrm{d}t = \left. f e^{-st}\right|_0^\infty + s\int_0^\infty f e^{-st}\mathrm{d}t = sF(s)-f(0)$。</div>
<div class="step"><b>卷积定理</b>：$\displaystyle\mathcal{L}[f*g]=\int_0^\infty\!\!\int_0^t f(\tau)g(t-\tau)\mathrm{d}\tau\, e^{-st}\mathrm{d}t$。交换积分次序（区域 $0\le\tau\le t<\infty$）得 $\displaystyle\int_0^\infty\!\!\int_\tau^\infty f(\tau)g(t-\tau)e^{-st}\mathrm{d}t\,\mathrm{d}\tau$。</div>
<div class="step">内层令 $u=t-\tau$，$e^{-st}=e^{-s\tau}e^{-su}$，于是两层积分分离：$\displaystyle\int_0^\infty f(\tau)e^{-s\tau}\mathrm{d}\tau\int_0^\infty g(u)e^{-su}\mathrm{d}u = F(s)G(s)$。</div>
<div class="step"><b>正弦</b>：$\sin at=\frac{e^{jat}-e^{-jat}}{2j}$，用线性性质和 (2.2.5)：$\frac{1}{2j}\left(\frac{1}{s-ja}-\frac{1}{s+ja}\right)=\frac{a}{s^2+a^2}$。</div>
</div></details>
<div class="warn" data-t="收敛域">积分要收敛才有意义。比如 $\mathcal{L}[e^{-at}]$ 要求 $e^{-(a+\sigma)t}\to 0$，即 $\sigma>-a$，这叫<b>收敛域</b>。$e^{-j\omega t}$ 只是在单位圆上转圈（模长恒为 1），只带来振荡，不影响收敛。</div>

<h2><span class="sec">2.2.3</span>逆变换：分母的根出现在指数上</h2>
<p>把 $F(s)$ 拆成简单分式之和（<span class="kw">部分分式</span>），每一项查表变回时间域。书中两个例子揭示了一个极其重要的规律：</p>
<div class="eq" data-no="例 2.2.6">
$$F(s)=\frac{-s+5}{s^2+5s+4}=\frac{-3}{s+4}+\frac{2}{s+1}\;\Longrightarrow\; f(t)=-3e^{-4t}+2e^{-t}$$
</div>
<div class="eq" data-no="例 2.2.7">
$$F(s)=\frac{4s+8}{s^2+2s+5}=\frac{2+j}{s+1+2j}+\frac{2-j}{s+1-2j}\;\Longrightarrow\; f(t)=e^{-t}\left(2\sin 2t+4\cos 2t\right)$$
</div>
<p>第一个例子分母的根是 $-4,-1$，它们正好出现在 $e^{-4t},e^{-t}$ 的指数上；第二个的根是 $-1\pm 2j$：实部 $-1$ 变成衰减包络 $e^{-t}$，虚部 $2$ 变成振荡频率 $\sin 2t,\cos 2t$。这些根叫<span class="kw">极点</span>。</p>
<div class="sim" data-id="splane"></div>
<div class="idea">极点的位置就是系统的“性格”：<b>实部</b>决定增长还是衰减（左半平面衰减、右半平面发散），<b>虚部</b>决定振不振、振多快。整个经典控制，很大程度上就是“把极点放到好位置”的学问。</div>

<h2><span class="sec">2.3</span>传递函数</h2>
<div class="eq" data-no="(2.3.1)" data-title="传递函数（零初始条件）">
$$G(s)=\frac{X(s)}{U(s)}$$
<dl><dt>$U(s),X(s)$</dt><dd>输入、输出的拉普拉斯变换</dd><dt>$G(s)$</dt><dd>传递函数：系统把输入“变成”输出的方式</dd></dl>
<div class="say">因为 $\mathcal{L}[\delta]=1$，冲激输入时 $X(s)=G(s)$——传递函数就是冲激响应的拉普拉斯变换。时域的 $x=u*h$ 在 s 域里变成了 $X=UG$。</div>
</div>
<p><b>回到 RL 电路</b>：对 $u = L\dot x + Rx$ 两边做变换（零初值）得 $U(s)=(Ls+R)X(s)$，于是 $G(s)=\dfrac{1}{Ls+R}$。输入恒定电压 $u=C$（$U=C/s$）时：</p>
<div class="eq" data-no="(2.3.6)–(2.3.9)">
$$X(s)=\frac{C}{s}\cdot\frac{1}{Ls+R}=\frac{C}{R}\left(\frac1s-\frac{1}{s+R/L}\right)\;\Longrightarrow\; x(t)=\frac{C}{R}\left(1-e^{-\frac{R}{L}t}\right)$$
<dl><dt>极点 $0$</dt><dd>由输入 $1/s$ 带来，对应稳态常数项 $C/R$（欧姆定律）</dd><dt>极点 $-R/L$</dt><dd>系统自身的极点（令 $G$ 的分母 $Ls+R=0$，即<b>特征方程</b>），决定过渡过程有多快</dd></dl>
</div>
<div class="sim" data-id="rl"></div>

<h3>控制系统的传递函数</h3>
<p>开环时，控制器 $C(s)$ 与对象 $G(s)$ 串联：$X=R\,C\,G$。闭环时误差 $E=R-X$，$X=E\,C\,G=(R-X)CG$，整理得：</p>
<div class="eq" data-no="(2.3.14)" data-title="闭环传递函数">
$$G_{cl}(s)=\frac{X(s)}{R(s)}=\frac{C(s)G(s)}{1+C(s)G(s)}$$
<div class="say">“前向通道”除以“1 + 回路增益”。整个闭环可以看成一个新的系统，它的极点由 $1+C(s)G(s)=0$ 决定——设计控制器 C(s)，本质上就是在挪动这些极点。</div>
</div>
<div class="sim" data-id="loop"></div>

<h2><span class="sec">2.4</span>非零初始状态</h2>
<p>传递函数要求零初值，那初值不为零怎么办？以 $\dot x+ax=u$ 为例，变换后 $sX-x(0)+aX=U$，即 $(s+a)X = U + x(0)$。可见初始状态相当于<b>多了一个输入</b> $x(0)$，它的逆变换是 $x(0)\delta(t)$——在 $t=0$ 的一瞬间给系统注入了一份“能量”。它不改变传递函数，也不改变极点，所以<b>不影响稳定性分析</b>。</p>

<div class="quiz"><div class="q">某系统传递函数的极点是 $-2\pm 5j$，它的冲激响应大概是什么样？</div>
<div class="opts"><button class="opt">单调衰减到 0，不振荡</button><button class="opt" data-ok>以 5 rad/s 振荡，振幅按 $e^{-2t}$ 衰减</button><button class="opt">以 2 rad/s 振荡，振幅越来越大</button></div>
<div class="why">实部 $-2<0$ → 包络 $e^{-2t}$ 衰减；虚部 $\pm5$ → 振荡角频率 5 rad/s。</div></div>

<div class="summary"><h3>本课要点</h3><ol>
<li>拉普拉斯变换是线性变换；求导 ↔ 乘 s，卷积 ↔ 相乘，于是微分方程变成代数方程。</li>
<li>部分分式逆变换揭示：分母的根（极点）出现在时间函数的指数上；实部管衰减/发散，虚部管振荡。</li>
<li>传递函数 = 零初值下输出与输入变换之比 = 冲激响应的变换；闭环传递函数 $G_{cl}=\frac{CG}{1+CG}$。</li>
<li>非零初值相当于 $t=0$ 时刻注入的一份能量，不改变极点。</li>
</ol></div>`,
    sims: {
      splane(el) {
        const ui = card(el, { title: "拖动极点，看时间函数怎么变", hint: "左图是复平面（s 平面），<b>拖动红色 ×</b>。右图是对应的时间函数（$F(s)=\\frac{1}{(s-p_1)(s-p_2)}$ 的逆变换；选例题时为书中的 $F(s)$）。绿色阴影是收敛域 $\\sigma>\\max\\mathrm{Re}(p)$。" });
        const [a, b] = ui.row("", "");
        const ps = new Plot(a, { height: 300, x: [-5, 2], y: [-4, 4], equal: true, xlabel: "σ（实部）", ylabel: "jω", title: "s 平面" });
        const pt = new Plot(b, { height: 300, x: [0, 8], y: [-1, 1], xlabel: "t", title: "时间函数 f(t)" });
        let mode = "cx", p = [[-1, 2], [-1, -2]], preset = "none";
        const seg = ui.seg({ label: "极点类型", options: [["cx", "一对共轭复极点"], ["re", "两个实极点"]], value: mode, onChange: v => { mode = v; preset = "none"; p = v === "cx" ? [[-1, 2], [-1, -2]] : [[-1, 0], [-3, 0]]; draw(); } });
        ui.buttons([["例 2.2.6", () => { mode = "re"; seg.set("re"); preset = "e226"; p = [[-1, 0], [-4, 0]]; draw(); }, "ghost"], ["例 2.2.7", () => { mode = "cx"; seg.set("cx"); preset = "e227"; p = [[-1, 2], [-1, -2]]; draw(); }, "ghost"]]);
        let drag = -1;
        ps.drag({
          down: (x, y) => { let best = -1, bd = 14; p.forEach((q, i) => { const d = ps.pxDist(q[0], q[1], x, y); if (d < bd) { bd = d; best = i; } }); drag = best; return best >= 0 ? undefined : false; },
          move: (x, y) => {
            preset = "none"; x = clamp(x, -4.8, 1.5); y = clamp(y, -3.8, 3.8);
            if (mode === "cx") { const s = drag === 0 ? 1 : -1; const w = Math.max(0.15, Math.abs(y)); p = [[x, w], [x, -w]]; }
            else { p[drag] = [x, 0]; }
            draw();
          },
          hover: (x, y) => p.some(q => ps.pxDist(q[0], q[1], x, y) < 14),
        });
        function ft() {
          if (preset === "e226") return t => -3 * Math.exp(-4 * t) + 2 * Math.exp(-t);
          if (preset === "e227") return t => Math.exp(-t) * (2 * Math.sin(2 * t) + 4 * Math.cos(2 * t));
          if (mode === "cx") { const [s, w] = p[0]; return t => Math.exp(s * t) * Math.sin(w * t) / w; }
          const [p1, p2] = [p[0][0], p[1][0]];
          if (Math.abs(p1 - p2) < 1e-3) return t => t * Math.exp(p1 * t);
          return t => (Math.exp(p1 * t) - Math.exp(p2 * t)) / (p1 - p2);
        }
        function draw() {
          const f = ft(), ts = linspace(0, 8, 500), ys = ts.map(f);
          const lim = Math.max(...ys.slice(0, 400).map(v => Math.abs(v)).filter(isFinite), 0.05);
          pt.setRange(null, [-Math.min(lim, 50) * 1.15, Math.min(lim, 50) * 1.15]);
          const sMax = Math.max(p[0][0], p[1][0]);
          ps.render = () => {
            ps.clear();
            ps.rect(sMax, ps.opt.y[0], ps.opt.x[1], ps.opt.y[1], COL.c3 + "18");
            ps.rect(0, ps.opt.y[0], ps.opt.x[1], ps.opt.y[1], COL.c4 + "10");
            ps.text(0.15, ps.opt.y[1] - 0.4, "右半平面：发散", COL.c4);
            ps.text(ps.opt.x[0] + 0.15, ps.opt.y[1] - 0.4, "左半平面：衰减", COL.c3);
            p.forEach(q => ps.cross(q[0], q[1], COL.c4, 8, 3));
          };
          pt.render = () => {
            pt.clear();
            if (mode === "cx" && preset === "none") { const s = p[0][0], w = p[0][1]; pt.line(ts, ts.map(t => Math.exp(s * t) / w), COL.muted, 1.2, [5, 4]); pt.line(ts, ts.map(t => -Math.exp(s * t) / w), COL.muted, 1.2, [5, 4]); }
            pt.line(ts, ys, COL.c1, 2.6);
          };
          ps.render(); pt.render();
          let form;
          if (preset === "e226") form = "$F(s)=\\frac{-s+5}{s^2+5s+4}$<br>$f(t)=-3e^{-4t}+2e^{-t}$";
          else if (preset === "e227") form = "$F(s)=\\frac{4s+8}{s^2+2s+5}$<br>$f(t)=e^{-t}(2\\sin2t+4\\cos2t)$";
          else if (mode === "cx") form = `$f(t)=\\frac{1}{${f2(p[0][1])}}e^{${f2(p[0][0])}t}\\sin(${f2(p[0][1])}t)$`;
          else form = `$f(t)=\\frac{e^{${f2(p[0][0])}t}-e^{${f2(p[1][0])}t}}{${f2(p[0][0] - p[1][0])}}$`;
          const st = sMax < -0.01 ? '<span class="good">衰减到 0（稳定）</span>' : sMax > 0.01 ? '<span class="bad">发散到无穷（不稳定）</span>' : "既不衰减也不发散（临界）";
          ui.out(`极点：<b>${CE.cfmt(p[0])}</b>，<b>${CE.cfmt(p[1])}</b><br>${form}<br>${st}${Math.abs(p[0][1]) > 0.01 ? `<br>振荡角频率 ω = ${f2(Math.abs(p[0][1]))} rad/s` : "<br>不振荡"}`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>把复极点往左拖：衰减变快；往右拖过虚轴：振幅越来越大。</li><li>把复极点上下拖：离实轴越远，振荡越快。</li><li>正好放在虚轴上：等幅振荡，永不停息（第 6 章的“临界稳定”）。</li></ol></div>`);
      },
      rl(el) {
        const ui = card(el, { title: "RL 电路的阶跃响应", hint: "合上开关（恒定电压 C），电流从 0 指数上升到 C/R。时间常数 τ = L/R 时到达终值的 63%。" });
        const [a, b] = ui.row("", "");
        const pt = new Plot(a, { height: 250, x: [0, 5], y: [0, 2.4], xlabel: "t / s", ylabel: "电流 x(t)" });
        const ps = new Plot(b, { height: 250, x: [-6, 1], y: [-2, 2], xlabel: "σ", title: "极点", equal: true });
        const sL = ui.slider({ label: "电感 L", min: 0.1, max: 3, step: 0.05, value: 1, onInput: draw });
        const sR = ui.slider({ label: "电阻 R", min: 0.5, max: 5, step: 0.05, value: 1, onInput: draw });
        const sC = ui.slider({ label: "电压 C", min: 0.5, max: 3, step: 0.1, value: 2, onInput: draw });
        const sX = ui.slider({ label: "初始电流 x(0)", min: 0, max: 2, step: 0.1, value: 0, onInput: draw });
        function draw() {
          const L = sL.get(), R = sR.get(), C = sC.get(), x0 = sX.get(), tau = L / R;
          const ts = linspace(0, 5, 300);
          pt.setRange(null, [0, Math.max(C / R, x0) * 1.25 + 0.1]);
          pt.render = () => {
            pt.clear(); pt.hline(C / R, COL.c3, [6, 4]); pt.text(4.95, C / R + 0.05 * pt.opt.y[1], "C/R = " + f2(C / R), COL.c3, "right");
            if (tau < 5) { pt.vline(tau, COL.muted, [3, 3]); pt.dot(tau, x0 + (C / R - x0) * (1 - Math.exp(-1)), COL.c5, 5); pt.text(tau + 0.08, 0.06 * pt.opt.y[1], "τ = L/R", COL.c5); }
            if (x0 > 0) { pt.line(ts, ts.map(t => C / R * (1 - Math.exp(-t / tau))), COL.muted, 1.4, [4, 4]); pt.line(ts, ts.map(t => x0 * Math.exp(-t / tau)), COL.c2, 1.4, [4, 4]); }
            pt.line(ts, ts.map(t => C / R * (1 - Math.exp(-t / tau)) + x0 * Math.exp(-t / tau)), COL.c1, 2.6);
            if (x0 > 0) pt.legend([{ label: "总响应", color: COL.c1 }, { label: "零状态", color: COL.muted, dash: [4, 4] }, { label: "x(0)h(t)", color: COL.c2, dash: [4, 4] }], "tr");
          };
          ps.render = () => { ps.clear(); ps.circle(0, 0, COL.c2, 7); ps.text(0.15, 0.4, "0（输入）", COL.c2); ps.cross(-R / L, 0, COL.c4, 8, 3); ps.text(-R / L, -0.6, "−R/L", COL.c4, "center"); };
          pt.render(); ps.render();
          ui.out(`$G(s)=\\frac{1}{${f2(L)}s+${f2(R)}}$<br>系统极点 $-R/L$ = <b>${f2(-R / L)}</b><br>时间常数 τ = <b>${f2(tau)} s</b><br>稳态电流 = <b>${f2(C / R)}</b>`);
        }
        draw();
      },
      loop(el) {
        const ui = card(el, { title: "闭环把极点“挪走了”", hint: "对象 $G(s)=\\frac{1}{s+1}$，控制器是一个放大倍数 $C(s)=K$。闭环传递函数 $G_{cl}=\\frac{K}{s+1+K}$，极点从 −1 移到 −(1+K)。" });
        const [a, b] = ui.row("", "");
        const pt = new Plot(a, { height: 240, x: [0, 5], y: [0, 1.2], xlabel: "t", title: "单位阶跃响应（目标 r=1）" });
        const ps = new Plot(b, { height: 240, x: [-12, 1], y: [-2, 2], xlabel: "σ", title: "闭环极点" });
        const sK = ui.slider({ label: "控制器增益 K", min: 0, max: 10, step: 0.1, value: 2, onInput: draw });
        function draw() {
          const K = sK.get(), ts = linspace(0, 5, 300), pole = -(1 + K), yss = K / (1 + K);
          pt.render = () => { pt.clear(); pt.hline(1, COL.c3, [6, 4]); pt.line(ts, ts.map(t => K > 0 ? yss * (1 - Math.exp(pole * t)) : 0), COL.c1, 2.6); pt.line(ts, ts.map(t => 1 - Math.exp(-t)), COL.muted, 1.4, [4, 4]); pt.legend([{ label: "闭环输出", color: COL.c1 }, { label: "开环 (K=1 直接给)", color: COL.muted, dash: [4, 4] }], "tr"); };
          ps.render = () => { ps.clear(); ps.cross(-1, 0, COL.muted, 6, 2); ps.text(-1, 0.5, "开环 −1", COL.muted, "center"); ps.cross(pole, 0, COL.c4, 8, 3); ps.text(pole, -0.7, f2(pole), COL.c4, "center"); };
          pt.render(); ps.render();
          ui.out(`$G_{cl}(s)=\\frac{${f2(K, 1)}}{s+${f2(1 + K, 1)}}$<br>闭环极点：<b>${f2(pole)}</b>（越靠左越快）<br>稳态输出：<b>${f2(yss, 3)}</b>，误差 ${f2(1 - yss, 3)}`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>观察</b>：K 越大，极点越靠左、响应越快、稳态误差 $\\frac{1}{1+K}$ 越小——但永远不为零。怎样彻底消除误差？第 7 章揭晓。</div>`);
        CE.renderMath(ui.root);
      },
    },
  });
})();
