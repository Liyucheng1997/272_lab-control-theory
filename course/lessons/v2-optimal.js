/* 卷2 第3–4章：最优控制概念 / 动态规划 / HJB / LQR / 轨迹追踪 */
(function () {
  const { M, P, Plot, Scene, card, linspace, logspace, f2, COL, clamp, rk4, cfmt, mtex } = CE;
  const frac = v => { if (!isFinite(v)) return "∞"; for (const d of [1, 2, 3, 4, 5, 6, 12, 15, 30]) { const n = Math.round(v * d); if (Math.abs(n / d - v) < 1e-9) return d === 1 ? String(n) : `${n}/${d}`; } return f2(v, 3); };

  /* ════════════════ 第3章 最优控制的基本概念 ════════════════ */
  Course.add({
    id: "v2-concept", vol: 2, ord: 40, no: "3", chName: "第3章 最优控制的基本概念", ref: "卷2 第3章",
    title: "最优控制：把“想要什么”写成代价函数",
    sub: "停车、省电、沿轨迹走、绕开障碍——每一种需求都可以写成代价函数里的一项。权重矩阵 S、Q、R 就是设计者的“旋钮”。",
    goals: ["理解性能指标（代价函数）的组成：末端代价 + 运行代价 + 控制代价", "理解权重矩阵 S、Q、R 的作用和硬约束/软约束的区别", "认识常见的最优控制问题类型"],
    minutes: 30,
    html: String.raw`
<p>经典控制关心稳定性、超调、稳态误差，手段是“配置极点”。但对高精度或多输入多输出系统，仅仅调特征值远远不够。<span class="kw">最优控制</span>要在可行的条件下，寻找<b>最好</b>的控制策略。</p>
<h2><span class="sec">3.1</span>引子：独轮车模型</h2>
<p>独轮车模型（unicycle model）是车辆的极简模型，广泛用于自动驾驶和 AGV 研究：小车在平面内能前进后退、能原地转向。</p>
<div class="eq" data-no="(3.1.2)" data-title="独轮车模型（非线性、多输入多状态）">
$$\frac{\mathrm d}{\mathrm dt}\begin{bmatrix}p_x\\p_y\\v\\\theta\end{bmatrix}=\begin{bmatrix}v\cos\theta\\v\sin\theta\\0\\0\end{bmatrix}+\begin{bmatrix}0\\0\\\alpha\\\omega\end{bmatrix}$$
<dl><dt>$p_x,p_y$</dt><dd>位置</dd><dt>$v,\theta$</dt><dd>速度、航向角</dd><dt>$\alpha,\omega$</dt><dd>两个输入：加速度、转向角速度</dd></dl>
</div>
<p>书中用四个越来越复杂的场景，一步步“搭”出代价函数：</p>
<table class="tbl"><tr><th>场景</th><th>需求</th><th>代价函数中新增的项</th></tr>
<tr><td>① 停车</td><td>N 步后停到车位、车头朝前</td><td>末端代价 $\|\bm x_{[N]}-\bm x_d\|_S^2=(\bm x_{[N]}-\bm x_d)^TS(\bm x_{[N]}-\bm x_d)$</td></tr>
<tr><td>② 省能量</td><td>在①的基础上少用油门和方向</td><td>控制代价 $\sum\|\bm u_{[k]}\|_R^2$（软约束）</td></tr>
<tr><td>③ 轨迹追踪</td><td>沿规划好的路线走</td><td>运行代价 $\sum\|\bm x_{[k]}-\bm x_{d[k]}\|_Q^2$</td></tr>
<tr><td>④ 避障</td><td>绕开障碍物</td><td>硬约束“不能进入障碍物”，或软约束：距离倒数惩罚 $D^{-1}$</td></tr></table>
<div class="eq" data-no="(3.2.7)" data-title="综合型性能指标（默认形式）">
$$J=\underbrace{\|\bm x_{[N]}-\bm x_{d[N]}\|_S^2}_{\text{末端代价}}+\sum_{k=0}^{N-1}\Big(\underbrace{\|\bm x_{[k]}-\bm x_{d[k]}\|_Q^2}_{\text{运行代价}}+\underbrace{\|\bm u_{[k]}\|_R^2}_{\text{控制代价}}\Big),\qquad \|\bm v\|_W^2\triangleq\bm v^TW\bm v$$
<dl><dt>$S,Q$</dt><dd>对称<b>半正定</b>：哪些状态更重要（比如停车时“停下来”比位置差几厘米更要紧）</dd><dt>$R$</dt><dd>对称<b>正定</b>：控制量的“价格”（舵机比驱动电机省电得多，所以 $r_1>r_2$）</dd><dt>目标</dt><dd>找控制序列 $U^*=[\bm u_{[0]},\dots,\bm u_{[N-1]}]=\arg\min J$</dd></dl>
<div class="say">取平方：正负误差不会抵消、处处可导、便于计算。调权重时要注意各量的单位与数量级，必要时先归一化。</div>
</div>
<div class="sim" data-id="unicycle"></div>
<div class="tip" data-t="硬约束与软约束">硬约束是必须严格满足的物理/安全限制（速度上下限、发动机最大加速度）；软约束写进代价函数，可以“商量”（能少用就少用）。调权重能把软约束变“硬”：如果令 $S=0$，节能成为唯一目标——小车会停在原地一动不动。</div>

<h2><span class="sec">3.2</span>最优控制问题的组成</h2>
<ol><li><b>系统模型</b>：一般是（离散）状态空间方程 $\bm x_{[k+1]}=f(\bm x_{[k]},\bm u_{[k]})$。</li>
<li><b>目标</b>：固定点、轨迹，或一个区域。</li>
<li><b>性能指标</b>：标量函数 $J=h(\bm x_{[N]})+\sum g(\bm x_{[k]},\bm u_{[k]})$，h 末端代价、g 运行代价。</li>
<li><b>约束</b>：状态的容许轨迹 $\bm x\in X$、控制的容许域 $\bm u\in\Omega$。约束太严可能根本无解。</li></ol>
<table class="tbl"><tr><th>问题类型</th><th>离散型性能指标</th></tr>
<tr><td>最短时间</td><td>$J=\sum 1=N$</td></tr><tr><td>末端控制</td><td>$J=\|\bm x_{[N]}-\bm x_d\|_S^2$</td></tr>
<tr><td>最小控制量</td><td>$J=\sum\|\bm u_{[k]}\|_R^2$（绝对最小能耗 = 一动不动，往往没有实际意义）</td></tr>
<tr><td>轨迹追踪</td><td>$J=\sum\|\bm x_{[k]}-\bm x_{d[k]}\|_Q^2$</td></tr>
<tr><td>综合</td><td>上式 (3.2.7)</td></tr>
<tr><td>调节问题</td><td>目标 $\bm x_d=0$：$J=\|\bm x_{[N]}\|_S^2+\sum(\|\bm x_{[k]}\|_Q^2+\|\bm u_{[k]}\|_R^2)$</td></tr></table>
<p>求解时要注意：最优解不一定存在（障碍物正好堵在路上）、也可能不唯一（只关心终点时）；尽量把指标和约束设计成<b>凸</b>的，才好找到全局最优。</p>

<h2><span class="sec">3.3</span>两个贯穿全书的案例</h2>
<p><b>平衡车</b>：简化为倒立摆，$\dot{\bm x}=\begin{bmatrix}0&1\\g/d&0\end{bmatrix}\bm x+\begin{bmatrix}0\\1\end{bmatrix}u$，目标是回到竖直（调节问题），$J=\int(\|\bm x\|_Q^2+\|u\|_R^2)\mathrm dt$。</p>
<p><b>无人机高度控制</b>：$m\ddot h=f-mg$。令 $u=f/m$，方程里多了一个常数 $-g$，不是标准形式。技巧：把 $x_3=-g$ 当作第三个状态（$\dot x_3=0$），得到<span class="kw">增广</span>系统：</p>
<div class="eq" data-no="(3.3.14)–(3.3.15)" data-title="无人机高度：增广状态空间">
$$\frac{\mathrm d}{\mathrm dt}\begin{bmatrix}h\\\dot h\\-g\end{bmatrix}=\begin{bmatrix}0&1&0\\0&0&1\\0&0&0\end{bmatrix}\begin{bmatrix}h\\\dot h\\-g\end{bmatrix}+\begin{bmatrix}0\\1\\0\end{bmatrix}u\;\xrightarrow{T_s=0.1}\;\bm x_{[k+1]}=\begin{bmatrix}1&0.1&0.005\\0&1&0.1\\0&0&1\end{bmatrix}\bm x_{[k]}+\begin{bmatrix}0.005\\0.1\\0\end{bmatrix}u_{[k]}$$
<div class="say">一个小矛盾：悬停时需要 $u=g$，而控制代价 $\|u\|_R^2$ 却总想把 u 往 0 推。怎么办？第 4 章 4.5 节揭晓。</div>
</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>最优控制 = 模型 + 目标 + 性能指标 + 约束。</li>
<li>综合指标 = 末端代价（S）+ 运行代价（Q）+ 控制代价（R）；权重表达“在意什么”。</li>
<li>硬约束必须满足；软约束写进代价、可以权衡。</li>
</ol></div>`,
    sims: {
      unicycle(el) {
        const ui = card(el, { title: "代价函数调音台：独轮车", hint: "选择场景，调权重，点“优化”——程序用数值优化（梯度下降）寻找 40 步的控制序列使 J 最小。蓝线是小车轨迹，橙色是目标/参考路线，灰圈是障碍物。" });
        const [a, b] = ui.row("", "");
        const pl = new Plot(a, { height: 340, x: [-1, 11], y: [-3, 7], equal: true, xlabel: "px", ylabel: "py" });
        const pc = new Plot(b, { height: 340, x: [0, 4], y: [0, 1], title: "代价组成" });
        let scen = 1;
        ui.seg({ label: "场景", options: [[1, "①停车"], [2, "②+省能量"], [3, "③+追踪"], [4, "④+避障"]], value: 1, onChange: v => { scen = v; restart(); } });
        const sS = ui.slider({ label: "末端权重 S（位置/速度/角度）", min: 0, max: 50, step: 0.5, value: 20, onInput: restart });
        const sV = ui.slider({ label: "末端速度权重 s₃（停得住）", min: 0, max: 50, step: 0.5, value: 20, onInput: restart });
        const sR = ui.slider({ label: "控制权重 R", min: 0.001, max: 2, value: 0.05, log: true, onInput: restart });
        const sQ = ui.slider({ label: "追踪权重 Q", min: 0, max: 5, step: 0.05, value: 1, onInput: restart });
        const sP = ui.slider({ label: "避障权重 P", min: 0, max: 5, step: 0.05, value: 1.5, onInput: restart });
        ui.buttons([["▶ 优化 / 暂停", () => run = !run], ["↺ 重置", restart, "ghost"]]);
        const N = 40, Ts = 0.25, x0 = [0, 0, 0, 0], xd = [9, 4, 0, 0];
        const ref = k => { const s = k / N; return [9 * s, 4 * Math.sin(Math.PI * s / 2) * s + 2 * Math.sin(Math.PI * s) * 1.2]; };
        const obs = [[4.2, 3.4, 0.9], [6.8, 1.6, 0.7]];
        let U, run = false, it = 0, m, v, best;
        function sim(U) { let x = x0.slice(); const tr = [x]; for (let k = 0; k < N; k++) { const al = clamp(U[2 * k], -2, 2), om = clamp(U[2 * k + 1], -1.5, 1.5); x = [x[0] + Ts * x[2] * Math.cos(x[3]), x[1] + Ts * x[2] * Math.sin(x[3]), x[2] + Ts * al, x[3] + Ts * om]; tr.push(x); } return tr; }
        function cost(U, parts) {
          const tr = sim(U), xN = tr[N], S = sS.get(), sv = sV.get();
          const term = S * ((xN[0] - xd[0]) ** 2 + (xN[1] - xd[1]) ** 2 + (xN[3] - xd[3]) ** 2) + sv * (xN[2]) ** 2;
          let ctrl = 0, trk = 0, ob = 0;
          for (let k = 0; k < N; k++) {
            if (scen >= 2) ctrl += sR.get() * (U[2 * k] ** 2 + 0.3 * U[2 * k + 1] ** 2);
            if (scen >= 3) { const r = ref(k); trk += sQ.get() * ((tr[k][0] - r[0]) ** 2 + (tr[k][1] - r[1]) ** 2) * 0.2; }
            if (scen >= 4) obs.forEach(([ox, oy, rr]) => { const d = Math.hypot(tr[k][0] - ox, tr[k][1] - oy) - rr; ob += sP.get() * (d < 0.05 ? 20 * (0.05 - d) + 1 / 0.05 : 1 / d) * 0.1; });
          }
          if (parts) return { term, ctrl, trk, ob };
          return term + ctrl + trk + ob;
        }
        function restart() { U = new Array(2 * N).fill(0); m = U.map(() => 0); v = U.map(() => 0); it = 0; draw(); }
        function iterate(n) {
          for (let r = 0; r < n; r++) {
            const f0 = cost(U), g = U.map((_, i) => { const h = 1e-4; U[i] += h; const f1 = cost(U); U[i] -= h; return (f1 - f0) / h; });
            it++; const lr = 0.05, b1 = 0.9, b2 = 0.999;
            for (let i = 0; i < U.length; i++) { m[i] = b1 * m[i] + (1 - b1) * g[i]; v[i] = b2 * v[i] + (1 - b2) * g[i] * g[i]; const mh = m[i] / (1 - b1 ** it), vh = v[i] / (1 - b2 ** it); U[i] = clamp(U[i] - lr * mh / (Math.sqrt(vh) + 1e-8), -3, 3); }
          }
        }
        function draw() {
          const tr = sim(U), pr = cost(U, true);
          pl.render = () => {
            pl.clear();
            if (scen >= 3) { const r = linspace(0, N, 80).map(k => ref(k)); pl.line(r.map(p => p[0]), r.map(p => p[1]), COL.c2, 2, [6, 4]); }
            if (scen >= 4) obs.forEach(([ox, oy, rr]) => { const th = linspace(0, 6.3, 40); pl.area(th.map(t => ox + rr * Math.cos(t)), th.map(t => oy + rr * Math.sin(t)), th.map(() => oy), COL.muted + "55"); pl.line(th.map(t => ox + rr * Math.cos(t)), th.map(t => oy + rr * Math.sin(t)), COL.muted, 1.5); });
            pl.rect(xd[0] - 0.6, xd[1] - 0.35, xd[0] + 0.6, xd[1] + 0.35, null, COL.c3);
            pl.line(tr.map(p => p[0]), tr.map(p => p[1]), COL.c1, 2.6);
            tr.forEach((p, k) => { if (k % 4 === 0) { pl.dot(p[0], p[1], COL.c1, 2.5); pl.arrow(p[0], p[1], p[0] + 0.4 * Math.cos(p[3]), p[1] + 0.4 * Math.sin(p[3]), COL.c1, 1.2, 5); } });
            const e = tr[N]; pl.dot(e[0], e[1], COL.c4, 5);
          };
          const vals = [pr.term, pr.ctrl, pr.trk, pr.ob], mx = Math.max(1e-3, ...vals);
          pc.setRange([0, 4], [0, mx * 1.15]);
          pc.render = () => { pc.clear(); ["末端", "控制", "追踪", "避障"].forEach((n, i) => { pc.rect(i + 0.15, 0, i + 0.85, vals[i], [COL.c4, COL.c5, COL.c2, COL.muted][i] + "aa"); pc.text(i + 0.5, -mx * 0.04, n, COL.ink2, "center"); pc.text(i + 0.5, vals[i] + mx * 0.02, f2(vals[i], 2), COL.ink2, "center"); }); };
          pl.render(); pc.render();
          const e = tr[N];
          ui.out(`迭代 ${it}　J = <b class="big">${f2(pr.term + pr.ctrl + pr.trk + pr.ob, 3)}</b><br>终点：(${f2(e[0], 1)}, ${f2(e[1], 1)})，末速度 ${f2(e[2], 2)}，航向 ${f2(e[3] * 57.3, 0)}°`);
        }
        restart();
        CE.loop(ui.root, () => { if (run) { iterate(6); draw(); if (it > 3000) run = false; } });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>场景①：只有末端代价，小车可能“急刹急转”；把 s₃ 调成 0——到了终点也停不下来。</li><li>场景②：加大 R，轨迹更“懒”、更平滑；把 S 调成 0，小车干脆不动了。</li><li>场景③④：加大 Q 让小车更贴近橙色路线；加大 P 让它离障碍物更远。</li></ol></div>`);
      },
    },
  });

  /* ════════════════ 4.1–4.2 动态规划（数值） ════════════════ */
  function dpGrid({ H = 10, nh = 6, vmax = 3, nv = 4, amin = -3, amax = 2 }) {
    const dh = H / (nh - 1), vs = linspace(0, vmax, nv), J = [], U = [], NX = [];
    for (let s = 0; s < nh; s++) { J.push(new Array(nv).fill(Infinity)); U.push(new Array(nv).fill(NaN)); NX.push(new Array(nv).fill(-1)); }
    J[nh - 1][0] = 0;
    for (let s = nh - 2; s >= 0; s--) {
      for (let i = 0; i < nv; i++) {
        if (s === 0 && i !== 0) continue;
        for (let j = 0; j < nv; j++) {
          if (s === nh - 2 && j !== 0) continue;
          const v1 = vs[i], v2 = vs[j], dt = (v1 + v2) > 0 ? 2 * dh / (v1 + v2) : Infinity, u = (v2 * v2 - v1 * v1) / (2 * dh);
          if (u < amin - 1e-12 || u > amax + 1e-12) continue;
          const c = dt + J[s + 1][j];
          if (c < J[s][i]) { J[s][i] = c; U[s][i] = u; NX[s][i] = j; }
        }
      }
    }
    return { J, U, NX, vs, dh, nh, nv };
  }
  Course.add({
    id: "v2-dp", vol: 2, ord: 50, no: "4.1–4.2", chName: "第4章 动态规划与 LQR", ref: "卷2 §4.1–4.2",
    title: "动态规划：倒着算，正着走",
    sub: "贝尔曼最优性原理：不管之前怎么走到这里，从这里出发的剩余路线必须是最优的。用它解决无人机的最短时间爬升问题。",
    goals: ["理解贝尔曼最优性原理", "会用逆向分级的方法手算一个网格动态规划", "理解动态规划相比穷举法在计算量上的巨大优势"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">4.1</span>贝尔曼最优性原理</h2>
<p>动态规划由 Richard Bellman 在 20 世纪 50 年代提出（programming 是“规划”而非“编程”）。核心原理：</p>
<div class="idea" data-t="贝尔曼最优性原理">一个最优策略具有这样的性质：无论初始状态和初始决策是什么，余下的决策对于由第一个决策所产生的状态而言，必须构成最优策略。<br>——换句话说：<b>最优路线的任何“后半段”，本身也是最优路线</b>。动态规划是面向未来的：不管过去走得好不好，从现在起的决策一定要最优。</div>
<div class="sim" data-id="graph"></div>

<h2><span class="sec">4.2</span>无人机最短时间爬升</h2>
<div class="story" data-t="问题">无人机（m=1 kg，g=10）要在<b>最短时间</b>内从地面静止升到 10 m 并悬停。约束：净加速度 $-3\le u_a\le2$ m/s²（减速能力比加速强），速度 $0\le v\le3$ m/s（不许下降）。</div>
<p>“时间”不是状态变量，很难直接写成控制量的函数。办法是离散化：高度按 2 m 分成 6 个“级”（0, 2, …, 10 m），速度取 0, 1, 2, 3 m/s 四个值。相邻两级之间加速度恒定，于是：</p>
<div class="eq" data-no="(4.2.5)" data-title="两级之间">
$$\Delta t=\frac{2\Delta x_1}{v_k+v_{k+1}},\qquad u_a=\frac{v_{k+1}^2-v_k^2}{2\Delta x_1}\qquad(\Delta x_1=2\text{ m})$$
<div class="say">平均速度 = 两端速度的平均；加速度 = 速度平方之差 / (2 × 距离)，中学物理就够了。</div>
</div>
<p><b>穷举法</b>：每一级 4 种速度，总共 $4^4=256$ 条路径。如果把高度分成 100 级、速度 500 档，就是 $1.9\times10^{267}$ 条——宇宙毁灭也算不完。</p>
<p><b>逆向分级</b>：从终点往回算。每个节点只需要记住两样东西：<b>剩余代价</b> $J^*$（从这里到终点的最短时间）和对应的<b>最优控制量</b> $u^*$。计算某个节点时，只需比较“走到下一级某节点的时间 + 那个节点的剩余代价”，取最小（违反约束的边记为 ∞）：</p>
<div class="eq" data-title="动态规划的核心递推">
$$J^*(\text{节点})=\min_{\text{下一级节点}}\Big[\Delta t(\text{节点}\to\text{下一节点})+J^*(\text{下一节点})\Big]$$
<div class="say">总计算量只有 $(6-2)\times4^2+2\times4=72$ 条边，随级数线性增长。100 级 × 500 档也只要 $2.4\times10^7$ 条。</div>
</div>
<div class="sim" data-id="dpgrid"></div>
<p>算完的表格就是一张“查表控制器”：所有计算都可以<b>离线</b>完成，实时控制时只需查当前状态对应的 $u^*$。书中得到的最优路线：$[0,0]\to[2,2]\to[4,3]\to[6,3]\to[8,3]\to[10,0]$，用时 $82/15\approx5.47$ s。如果在 2 m 处突然失速（状态变成 [2,0]），查表即可得到新的最优路线，并在两级之后与原最优路线重合——贝尔曼原理的生动体现。</p>
<div class="sim" data-id="dpfine"></div>
<div class="tip" data-t="编程技巧：矩阵代替循环">书中建议把两级之间的 16 条边写成矩阵：$v_{dx}$（当前速度，按列变化）和 $v_{dy}$（下一级速度，按行变化），一次算出 $\Delta t$、$u_a$ 矩阵；把违约元素置 ∞；加上下一级的剩余代价后按列取最小值——五步完成一级。</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>贝尔曼原理：最优策略的剩余部分对剩余问题仍然最优。</li>
<li>逆向分级：从终点往回，每个节点只保留最优剩余代价与最优控制。</li>
<li>计算量从指数级降到线性；结果可以离线做成查找表。</li>
<li>最短时间的典型解：最大加速 → 匀速 → 最大减速（bang–coast–bang）。</li>
</ol></div>`,
    sims: {
      graph(el) {
        const ui = card(el, { title: "最短路径与“殊途同归”", hint: "图 4.1.1：从 A 到 C，中间方块是障碍物。数字是剩余代价（到 C 的最短距离），由终点往回逐个算出。<b>点击任意节点</b>，查看从它出发的最优路径。" });
        const pl = new Plot(ui.plotHost(), { height: 320, x: [-0.4, 4.4], y: [-0.2, 3.3], equal: true, grid: false });
        const nodes = { A: [0, 1], B: [1, 1], D: [2, 1.6], E: [3, 1.6], F: [2, 0.4], G: [3, 0.4], C: [4, 1], H: [2.5, 2.9] };
        const edges = [["A", "B", 1], ["B", "D", 2], ["B", "F", 2], ["D", "E", 1], ["F", "G", 1], ["D", "F", 1], ["E", "G", 1], ["E", "C", 2], ["G", "C", 2], ["B", "H", 5], ["H", "E", 2], ["H", "C", 5]];
        const adj = {}; Object.keys(nodes).forEach(n => adj[n] = []); edges.forEach(([a, b, w]) => { adj[a].push([b, w]); adj[b].push([a, w]); });
        const dist = { C: 0 }, nxt = {}; const Q = Object.keys(nodes);
        Object.keys(nodes).forEach(n => { if (n !== "C") dist[n] = Infinity; });
        while (Q.length) { Q.sort((a, b) => dist[a] - dist[b]); const u = Q.shift(); adj[u].forEach(([v, w]) => { if (dist[u] + w < dist[v]) { dist[v] = dist[u] + w; nxt[v] = u; } }); }
        let sel = "A", reveal = 0; const order = Object.keys(nodes).sort((a, b) => dist[a] - dist[b]);
        ui.buttons([["▶ 逐个揭示剩余代价", () => { reveal = Math.min(order.length, reveal + 1); draw(); }], ["全部显示", () => { reveal = order.length; draw(); }, "ghost"]]);
        pl.drag({ down: (x, y) => { for (const [n, p] of Object.entries(nodes)) if (pl.pxDist(p[0], p[1], x, y) < 18) { sel = n; draw(); } return false; } });
        function path(s) { const p = [s]; while (p[p.length - 1] !== "C") p.push(nxt[p[p.length - 1]]); return p; }
        function draw() {
          const pth = path(sel), on = new Set(pth.slice(1).map((n, i) => pth[i] + n));
          pl.render = () => {
            pl.clear(); pl.rect(2, 0.4, 3, 1.6, COL.muted + "40", COL.muted); pl.text(2.5, 0.95, "障碍物", COL.muted, "center");
            edges.forEach(([a, b, w]) => { const A = nodes[a], B = nodes[b], hl = on.has(a + b) || on.has(b + a); pl.line([A[0], B[0]], [A[1], B[1]], hl ? COL.c4 : COL.axis, hl ? 4 : 1.6); pl.text((A[0] + B[0]) / 2 + 0.06, (A[1] + B[1]) / 2 + 0.08, w, COL.ink2); });
            Object.entries(nodes).forEach(([n, p]) => { const shown = order.indexOf(n) < reveal; pl.dot(p[0], p[1], n === sel ? COL.c4 : COL.c1, 13); pl.text(p[0], p[1] - 0.06, n, "#fff", "center", 13, true); if (shown) pl.text(p[0], p[1] + 0.3, "J*=" + dist[n], COL.c3, "center", 12, true); });
          };
          pl.render();
          ui.out(`从 <b>${sel}</b> 出发的最优路径：<br><b class="big">${pth.join("→")}</b><br>长度 ${dist[sel]}<br><span class="sim-note">注意：从 D、F、H 出发，最后一段都汇入 E–C 或 G–C。点 A 有两条等长最优路径（经 D–E 或 F–G）。</span>`);
        }
        draw();
      },
      dpgrid(el) {
        const ui = card(el, { title: "无人机 6×4 网格动态规划：一步一步填表", hint: "每个圆圈是一个状态（高度, 速度）。点“下一级”从终点往回计算；<b>点击任一节点</b>查看它的 4 条候选边（Δt 与 u）以及比较过程。" });
        const pl = new Plot(ui.plotHost(), { height: 420, x: [-0.6, 3.9], y: [-1, 11.3], xlabel: "速度 x₂ (m/s)", ylabel: "高度 x₁ (m)", xticks: [0, 1, 2, 3], yticks: [0, 2, 4, 6, 8, 10] });
        const G = dpGrid({}), nh = 6, nv = 4;
        let shown = nh - 1, sel = null, dist = false;
        ui.buttons([["⏮ 下一级（逆向）", () => { shown = Math.max(0, shown - 1); draw(); }], ["全部算完", () => { shown = 0; draw(); }, "ghost"], ["↺", () => { shown = nh - 1; sel = null; dist = false; draw(); }, "ghost"]]);
        ui.buttons([["⚡ 在 [2,0] 失速", () => { shown = 0; dist = true; draw(); }, "ghost"]]);
        pl.drag({ down: (x, y) => { for (let s = 0; s < nh; s++) for (let i = 0; i < nv; i++) if (pl.pxDist(i, 2 * s, x, y) < 16) { sel = [s, i]; draw(); } return false; } });
        const valid = (s, i) => !(s === 0 && i > 0) && !(s === nh - 1 && i > 0);
        function route(s0, i0) { const r = [[s0, i0]]; let s = s0, i = i0; while (s < nh - 1 && G.NX[s][i] >= 0) { i = G.NX[s][i]; s++; r.push([s, i]); } return r; }
        function draw() {
          pl.render = () => {
            pl.clear();
            if (shown === 0) { const r = route(0, 0); for (let k = 0; k < r.length - 1; k++) { pl.arrow(r[k][1], 2 * r[k][0], r[k + 1][1], 2 * r[k + 1][0], COL.c4, 3.2, 10); pl.text((r[k][1] + r[k + 1][1]) / 2 + 0.08, (2 * r[k][0] + 2 * r[k + 1][0]) / 2, "u=" + frac(G.U[r[k][0]][r[k][1]]), COL.c4); } }
            if (dist) { const r = route(1, 0); for (let k = 0; k < r.length - 1; k++) pl.arrow(r[k][1], 2 * r[k][0], r[k + 1][1], 2 * r[k + 1][0], COL.c5, 2.4, 9); }
            if (sel && valid(...sel) && sel[0] < nh - 1 && sel[0] >= shown) {
              const [s, i] = sel; for (let j = 0; j < nv; j++) { if (!valid(s + 1, j)) continue; const v1 = G.vs[i], v2 = G.vs[j], dt = v1 + v2 > 0 ? 4 / (v1 + v2) : Infinity, u = (v2 * v2 - v1 * v1) / 4, ok = u >= -3 && u <= 2, best = G.NX[s][i] === j; pl.line([i, j], [2 * s, 2 * s + 2], best ? COL.c3 : ok ? COL.c2 + "aa" : COL.c4 + "66", best ? 3.5 : 1.6, ok ? null : [4, 3]); }
            }
            for (let s = 0; s < nh; s++) for (let i = 0; i < nv; i++) {
              if (!valid(s, i)) continue;
              const done = s >= shown, isSel = sel && sel[0] === s && sel[1] === i;
              pl.dot(i, 2 * s, isSel ? "#e8b400" : done ? COL.c1 : COL.axis, 14);
              if (done && s < nh - 1) { pl.text(i, 2 * s + 0.75, "J*=" + frac(G.J[s][i]), COL.ink, "center", 11, true); pl.text(i, 2 * s - 0.95, "u*=" + frac(G.U[s][i]), COL.c5, "center", 10.5); }
            }
            pl.text(0, 10.9, "终点", COL.c3, "center", 12, true); pl.text(0.25, -0.75, "起点", COL.c3, "left", 12, true);
          };
          pl.render();
          let info = shown === nh - 1 ? "从终点 [10, 0] 开始。点击“下一级”。" : `已完成第 ${shown + 1} 级（高度 ${2 * shown} m）及以上的计算。`;
          if (sel && valid(...sel) && sel[0] < nh - 1) {
            const [s, i] = sel; const rows = [];
            for (let j = 0; j < nv; j++) { if (!valid(s + 1, j)) continue; const v1 = G.vs[i], v2 = G.vs[j], dt = v1 + v2 > 0 ? 4 / (v1 + v2) : Infinity, u = (v2 * v2 - v1 * v1) / 4, ok = u >= -3 && u <= 2; rows.push(`<tr><td>→v=${v2}</td><td>${frac(dt)}</td><td>${frac(u)}</td><td>${ok ? frac(dt + G.J[s + 1][j]) : '<span class="bad">违约</span>'}</td></tr>`); }
            info += `<br><b>节点 [${2 * s}, ${G.vs[i]}]</b> 的候选：<table class="tbl" style="margin:4px 0;font-size:.78rem"><tr><th></th><th>Δt</th><th>u</th><th>Δt+J*</th></tr>${rows.join("")}</table>`;
          }
          if (shown === 0) info += `<br>最优总时间 <b class="big">${frac(G.J[0][0])} ≈ ${f2(G.J[0][0], 2)} s</b>`;
          if (dist) info += `<br><span style="color:var(--c5)">紫线：从 [2,0] 出发的新最优路线，两级后并入原路线</span>`;
          ui.out(info);
        }
        draw();
      },
      dpfine(el) {
        const ui = card(el, { title: "精细网格：真正的最短时间轨迹", hint: "把高度和速度分得更细，用同一个动态规划算法求解（在你的浏览器里实时计算）。观察“最大加速—匀速—最大减速”的结构。" });
        const [a, b] = ui.row("", "");
        const ph = new Plot(a, { height: 250, x: [0, 3.2], y: [0, 10.5], xlabel: "速度 m/s", ylabel: "高度 m", title: "相平面" });
        const pt = new Plot(b, { height: 250, x: [0, 7], y: [-3.5, 3.5], xlabel: "t / s", title: "速度 v(t) 与加速度 u(t)" });
        const snh = ui.slider({ label: "高度分级数", min: 6, max: 150, step: 1, value: 100, fmt: v => v.toFixed(0), onInput: run });
        const snv = ui.slider({ label: "速度分档数", min: 4, max: 300, step: 1, value: 150, fmt: v => v.toFixed(0), onInput: run });
        const sam = ui.slider({ label: "最大加速度", min: 0.5, max: 4, step: 0.1, value: 2, onInput: run });
        const sad = ui.slider({ label: "最大减速度", min: 0.5, max: 5, step: 0.1, value: 3, onInput: run });
        function run() {
          const t0 = performance.now();
          const G = dpGrid({ nh: snh.get(), nv: snv.get(), amin: -sad.get(), amax: sam.get() });
          const ms = performance.now() - t0;
          let s = 0, i = 0, t = 0; const hs = [0], vs = [0], ts = [0], us = [];
          while (s < G.nh - 1 && G.NX[s][i] >= 0) { const j = G.NX[s][i], v1 = G.vs[i], v2 = G.vs[j]; const dt = 2 * G.dh / (v1 + v2); us.push([t, G.U[s][i]], [t + dt, G.U[s][i]]); t += dt; s++; i = j; hs.push(s * G.dh); vs.push(v2); ts.push(t); }
          ph.render = () => { ph.clear(); ph.line(vs, hs, COL.c1, 2.6); };
          pt.setRange([0, Math.max(5, t * 1.05)], null);
          pt.render = () => { pt.clear(); pt.line(ts, vs, COL.c1, 2.4); pt.line(us.map(p => p[0]), us.map(p => p[1]), COL.c4, 1.8); pt.legend([{ label: "速度 v", color: COL.c1 }, { label: "加速度 u", color: COL.c4 }], "tr"); };
          ph.render(); pt.render();
          const a = sam.get(), d = sad.get(), vm = 3, ha = vm * vm / (2 * a), hd = vm * vm / (2 * d);
          const th = ha + hd <= 10 ? vm / a + vm / d + (10 - ha - hd) / vm : (() => { const vp = Math.sqrt(2 * 10 * a * d / (a + d)); return vp / a + vp / d; })();
          ui.out(`DP 最短时间：<b class="big">${f2(G.J[0][0], 3)} s</b><br>理论连续最优：${f2(th, 3)} s<br>计算的边数：${((G.nh - 2) * G.nv * G.nv + 2 * G.nv).toLocaleString()}，用时 ${ms.toFixed(0)} ms`);
        }
        run();
      },
    },
  });

  /* ════════════════ 4.3 递归关系与 HJB ════════════════ */
  Course.add({
    id: "v2-hjb", vol: 2, ord: 60, no: "4.3", chName: "第4章 动态规划与 LQR", ref: "卷2 §4.3",
    title: "动态规划的解析形式：贝尔曼方程与 HJB 方程",
    sub: "对线性系统，“查表”可以换成公式。离散时得到递归的贝尔曼方程，连续时得到一个偏微分方程——HJB 方程。",
    goals: ["写出离散贝尔曼递归方程，并用它手算例 4.3.1", "理解连续系统的 HJB 方程和哈密顿项", "发现“剩余代价是二次型、最优控制是线性反馈”这一规律"],
    minutes: 30,
    html: String.raw`
<h2><span class="sec">4.3.1</span>离散系统的递归关系</h2>
<p>系统 $\bm x_{[k+1]}=f(\bm x_{[k]},\bm u_{[k]})$，代价 $J=h(\bm x_{[N]})+\sum_{k=0}^{N-1}g(\bm x_{[k]},\bm u_{[k]})$。从最后一级往回：$J^*_{N\to N}=h(\bm x_{[N]})$；往前一步时，把 $\bm x_{[N]}=f(\bm x_{[N-1]},\bm u_{[N-1]})$ 代进去，就只剩 $\bm u_{[N-1]}$ 一个未知量，求导令其为零即可。一般地：</p>
<div class="eq" data-no="(4.3.12)" data-title="贝尔曼动态规划方程">
$$J^*_{N-k\to N}(\bm x_{[N-k]})=\min_{\bm u_{[N-k]}}\Big(\underbrace{J^*_{N-(k-1)\to N}\big(f(\bm x_{[N-k]},\bm u_{[N-k]})\big)}_{\text{落脚点的最优剩余代价}}+\underbrace{g(\bm x_{[N-k]},\bm u_{[N-k]})}_{\text{这一步的代价}}\Big)$$
<div class="say">从 x 出发的最优总代价 = 在所有可选动作中，取“这一步代价 + 落脚点的最优剩余代价”的最小值。</div>
</div>
<p><b>例 4.3.1</b>：$x_{[k+1]}=x_{[k]}+u_{[k]}$，$x_{[0]}=1$，$J=\frac12x_{[2]}^2+\frac12\sum_{k=0}^{1}(x_{[k]}^2+u_{[k]}^2)$，N=2。</p>
<details class="derive" open><summary>手算过程</summary><div class="body">
<div class="step">$k=1$：$J_{1\to2}=\frac12(x_{[1]}+u_{[1]})^2+\frac12(x_{[1]}^2+u_{[1]}^2)$，求导：$x_{[1]}+2u_{[1]}=0\Rightarrow u^*_{[1]}=-\frac12x_{[1]}$，代回得 $J^*_{1\to2}=\frac34x_{[1]}^2$。</div>
<div class="step">$k=0$：由贝尔曼原理，$J_{0\to2}=\frac34x_{[1]}^2+\frac12(x_{[0]}^2+u_{[0]}^2)=\frac34(x_{[0]}+u_{[0]})^2+\frac12(x_{[0]}^2+u_{[0]}^2)$，求导得 $u^*_{[0]}=-\frac35x_{[0]}$，$J^*_{0\to2}=\frac45x_{[0]}^2$。</div>
<div class="step"><b>求解逆向，施加正向</b>：$x_{[0]}=1\to u=-\frac35\to x_{[1]}=\frac25\to u=-\frac15\to x_{[2]}=\frac15$。验证：$J=0.02+0.78=0.8=\frac45$ ✓。</div>
</div></details>
<div class="sim" data-id="rec1d"></div>
<div class="idea">两个规律：① 剩余代价总是 $\frac12P x^2$ 的二次型；② 最优控制总是线性反馈 $u^*=-Fx$。而且 N 越大，离终点很远时 F 趋于一个常数（这里是黄金分割 0.618！）。下一课的 LQR 正是把这一规律推广到矩阵。</div>

<h2><span class="sec">4.3.3</span>连续系统：HJB 方程</h2>
<p>连续系统不能直接把状态方程代进去消元。办法是把区间 $[t,t_f]$ 拆成 $[t,t+\Delta t]$ 和 $[t+\Delta t,t_f]$，对剩余代价做泰勒展开，再令 $\Delta t\to0$：</p>
<div class="eq" data-no="(4.3.39)–(4.3.43)" data-title="哈密顿–雅可比–贝尔曼（HJB）方程">
$$0=J_t^*(\bm x,t)+\min_{\bm u}\Big\{\underbrace{J_x^*(\bm x,t)^Tf(\bm x,\bm u,t)+g(\bm x,\bm u,t)}_{\text{哈密顿项 }\mathcal H}\Big\},\qquad J^*(\bm x_{(t_f)},t_f)=h(\bm x_{(t_f)},t_f)$$
<dl><dt>$J_t^*$</dt><dd>最优代价对时间的偏导</dd><dt>$J_x^*$</dt><dd>最优代价对状态的偏导（分母布局的列向量）</dd><dt>$\mathcal H$</dt><dd>哈密顿项；最优控制使它最小</dd></dl>
<div class="say">用法：先对 u 求 $\mathcal H$ 的极小，得到 $u^*$ 关于 $J_x^*$ 的表达式，再代回解这个偏微分方程。</div>
</div>
<p><b>例 4.3.2</b>：$\dot x=x+u$（开环不稳定），$h=\frac12x_{(t_f)}^2$，$g=\frac12u^2$。$\mathcal H=J_x^*(x+u)+\frac12u^2$，令 $\partial\mathcal H/\partial u=0$ 得 $u^*=-J_x^*$。猜 $J^*=\frac12k(t)x^2$，代入 HJB 得 $\dot k=k^2-2k$，边界 $k(t_f)=1$，解得：</p>
<div class="eq" data-no="(4.3.53)">
$$k(t)=\frac{2}{1+e^{-2(t_f-t)}},\qquad u^*=-k(t)\,x,\qquad \lim_{t_f\to\infty}k=2$$
</div>
<div class="sim" data-id="hjb"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>离散贝尔曼方程：$J^*_{N-k}=\min_u(J^*_{N-k+1}(f(x,u))+g(x,u))$，逆向求解、正向施加。</li>
<li>连续 HJB：$0=J_t^*+\min_u\mathcal H$，$\mathcal H=J_x^{*T}f+g$。</li>
<li>线性系统 + 二次型代价：剩余代价是二次型，最优控制是线性状态反馈，远离终点时增益趋于常数。</li>
</ol></div>`,
    sims: {
      rec1d(el) {
        const ui = card(el, { title: "一维递归：倒着算 P 和 F，正着跑", hint: "$x_{[k+1]}=x_{[k]}+u_{[k]}$。上图：逆向求出的增益 $F$（u=−Fx）；下图：正向运行的状态与控制量。" });
        const pf = new Plot(ui.plotHost(), { height: 200, x: [-0.5, 10.5], y: [0, 1], xlabel: "距离终点的步数", title: "增益 F（逆向递推）" });
        const px = new Plot(ui.plotHost(), { height: 220, x: [-0.5, 10.5], y: [-0.7, 1.1], xlabel: "k", title: "状态 x（蓝）与控制 u（红）" });
        const sN = ui.slider({ label: "末端时刻 N", min: 1, max: 10, step: 1, value: 2, fmt: v => v.toFixed(0), onInput: draw });
        const sQ = ui.slider({ label: "Q", min: 0, max: 5, step: 0.1, value: 1, onInput: draw });
        const sR = ui.slider({ label: "R", min: 0.1, max: 5, step: 0.1, value: 1, onInput: draw });
        const sS = ui.slider({ label: "S", min: 0, max: 5, step: 0.1, value: 1, onInput: draw });
        function draw() {
          const N = sN.get(), Q = sQ.get(), R = sR.get(), S = sS.get();
          let Pk = S; const F = [], Ps = [S];
          for (let k = 1; k <= N; k++) { const f = Pk / (Pk + R); F.push(f); Pk = (1 - f) ** 2 * Pk + f * f * R + Q; Ps.push(Pk); }
          let x = 1, J = 0; const xs = [1], us = [];
          for (let k = 0; k < N; k++) { const u = -F[N - 1 - k] * x; us.push(u); J += 0.5 * (Q * x * x + R * u * u); x = x + u; xs.push(x); }
          J += 0.5 * S * x * x;
          pf.setRange([-0.5, N + 0.5], [0, 1]); px.setRange([-0.5, N + 0.5], null);
          pf.render = () => { pf.clear(); F.forEach((f, i) => { pf.rect(i + 1 - 0.3, 0, i + 1 + 0.3, f, COL.c5 + "aa"); pf.text(i + 1, f + 0.04, frac(f), COL.ink2, "center", 10.5); }); pf.hline(0.618, COL.c3, [5, 4]); };
          px.render = () => { px.clear(); xs.forEach((v, k) => { px.line([k, k], [0, v], COL.c1, 5); }); us.forEach((v, k) => px.line([k + 0.2, k + 0.2], [0, v], COL.c4, 5)); };
          pf.render(); px.render();
          ui.out(`P 序列：${Ps.map(frac).join(" → ")}<br>F 序列（从终点往回）：${F.map(frac).join(", ")}<br>最优代价 J* = ½P·x₀² = <b>${frac(0.5 * Ps[N])}</b>，逐步累加 J = ${f2(J, 4)}`);
        }
        draw();
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b>：把 N 调到 10（Q=R=S=1）：F 依次为 1/2, 3/5, 8/13, 21/34……斐波那契数之比，迅速收敛到 0.618。</div>`);
      },
      hjb(el) {
        const ui = card(el, { title: "HJB 的解：增益 k(t) 从终点往回“饱和”", hint: "左：$k(t)=\\frac{2}{1+e^{-2(t_f-t)}}$。右：闭环 $\\dot x=(1-k(t))x$ 的状态，以及不加控制时的发散曲线。" });
        const [a, b] = ui.row("", "");
        const pk = new Plot(a, { height: 240, x: [0, 5], y: [0, 2.2], xlabel: "t", title: "k(t)" });
        const px = new Plot(b, { height: 240, x: [0, 5], y: [-0.2, 2], xlabel: "t", title: "x(t)，x₀ = 1" });
        const st = ui.slider({ label: "末端时间 t_f", min: 0.5, max: 5, step: 0.05, value: 3, onInput: draw });
        function draw() {
          const tf = st.get(), ts = linspace(0, tf, 300), k = t => 2 / (1 + Math.exp(-2 * (tf - t)));
          let x = 1; const xs = [1]; for (let i = 1; i < ts.length; i++) { const dt = ts[i] - ts[i - 1]; x = x + dt * (1 - k(ts[i - 1])) * x; xs.push(x); }
          pk.render = () => { pk.clear(); pk.hline(2, COL.c3, [5, 4]); pk.hline(1, COL.muted, [3, 3]); pk.line(ts, ts.map(k), COL.c5, 2.6); pk.vline(tf, COL.muted, [3, 3]); };
          px.render = () => { px.clear(); px.line(ts, ts.map(t => Math.exp(t)), COL.c4, 1.6, [5, 4]); px.line(ts, xs, COL.c1, 2.6); px.legend([{ label: "最优控制", color: COL.c1 }, { label: "不控制 eᵗ", color: COL.c4, dash: [5, 4] }], "tr"); };
          pk.render(); px.render();
          ui.out(`k(0) = ${f2(k(0), 4)}，k(t_f) = 1<br>t_f 越长，前段 k 越接近 2：u* = −2x，闭环 ẋ = −x 稳定`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 4.4 LQR ════════════════ */
  Course.add({
    id: "v2-lqr", vol: 2, ord: 70, no: "4.4", chName: "第4章 动态规划与 LQR", ref: "卷2 §4.4",
    title: "线性二次型调节器 LQR",
    sub: "线性系统 + 二次型代价 + 目标是回到零——最经典、应用最广的最优控制器。离散时是 Riccati 递推，连续时是 Riccati 方程。",
    goals: ["推导离散 LQR 的递推公式 $F=(B^TPB+R)^{-1}B^TPA$", "理解连续 LQR 的 Riccati 方程与代数 Riccati 方程", "会用 Q、R 权衡性能与能耗（平衡车案例）"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">4.4.1</span>离散型 LQR</h2>
<div class="eq" data-no="(4.4.1)–(4.4.2)" data-title="问题">
$$\bm x_{[k+1]}=A\bm x_{[k]}+B\bm u_{[k]},\qquad J=\tfrac12\bm x_{[N]}^TS\bm x_{[N]}+\tfrac12\sum_{k=0}^{N-1}\left(\bm x_{[k]}^TQ\bm x_{[k]}+\bm u_{[k]}^TR\bm u_{[k]}\right)$$
<div class="say">S、Q 对称半正定，R 对称正定；系数 ½ 只是为了求导时消掉 2。</div>
</div>
<p>完全沿用上一课的逆向分级：最后一步 $J^*=\frac12\bm x^TP_{[0]}\bm x$，$P_{[0]}=S$。往前一步，代入状态方程，用第 2 章的链式法则对 $\bm u$ 求导：</p>
<details class="derive"><summary>推导一步（之后每一步都一样）</summary><div class="body">
<div class="step">$J_{N-1\to N}=\frac12(A\bm x+B\bm u)^TP_{[0]}(A\bm x+B\bm u)+\frac12(\bm x^TQ\bm x+\bm u^TR\bm u)$。</div>
<div class="step">$\dfrac{\partial J}{\partial\bm u}=B^TP_{[0]}(A\bm x+B\bm u)+R\bm u=0\Rightarrow\bm u^*=-(B^TP_{[0]}B+R)^{-1}B^TP_{[0]}A\,\bm x$。</div>
<div class="step">二阶导 $B^TP_{[0]}B+R$ 正定（因为 R 正定）→ 确为极小。这就是 R 必须正定的原因。</div>
<div class="step">代回得 $J^*_{N-1\to N}=\frac12\bm x^TP_{[1]}\bm x$，形式与上一步完全相同——归纳得到一般递推。</div>
</div></details>
<div class="eq" data-no="(4.4.23)–(4.4.24)" data-title="离散 LQR（Riccati 递推）">
$$F_{[N-k]}=\left(B^TP_{[k-1]}B+R\right)^{-1}B^TP_{[k-1]}A,\qquad P_{[k]}=(A-BF)^TP_{[k-1]}(A-BF)+F^TRF+Q,\qquad P_{[0]}=S$$
<div class="say">$\bm u^*_{[N-k]}=-F_{[N-k]}\bm x_{[N-k]}$。求解（逆向）可以离线完成，运行（正向）时只做一次矩阵乘法。系统能控且 N→∞ 时 F 收敛到常数——这是 LQR 最常用的形式（书中的 [F1] 模块）。</div>
</div>
<p>用它重算例 4.3.1：$F_{[1]}=\frac12$，$P_{[1]}=\frac32$，$F_{[0]}=\frac35$——和动态规划完全一样。N=8 时 F 依次为 0.5, 0.6, 0.61538, 0.61765, 0.61798, 0.61803……</p>
<div class="sim" data-id="riccati"></div>

<h2><span class="sec">4.4.3</span>连续型 LQR</h2>
<p>对 $\dot{\bm x}=A\bm x+B\bm u$、$J=\frac12\bm x_{(t_f)}^TS\bm x_{(t_f)}+\int\frac12(\bm x^TQ\bm x+\bm u^TR\bm u)\mathrm dt$，用 HJB：哈密顿项对 u 求极小得 $\bm u^*=-R^{-1}B^TJ_x^*$；猜 $J^*=\frac12\bm x^TP(t)\bm x$，再用“对称化技巧” $\bm x^TPA\bm x=\frac12\bm x^T(PA+A^TP)\bm x$，得到：</p>
<div class="eq" data-no="(4.4.40)–(4.4.41)" data-title="连续 LQR">
$$-\dot P=PA+A^TP-PBR^{-1}B^TP+Q,\;\;P(t_f)=S\qquad\xrightarrow{\;\dot P=0\;}\qquad 0=PA+A^TP-PBR^{-1}B^TP+Q\;(\text{代数 Riccati 方程})$$
<div class="say">$\bm u^*=-K\bm x$，$K=R^{-1}B^TP$。一维例子 $\dot x=x+u$，$Q=0$，$R=1$：$0=2P-P^2\Rightarrow P=2$，$u^*=-2x$，与 HJB 结果一致。</div>
</div>

<h2><span class="sec">4.4.5</span>平衡车（倒立摆）</h2>
<p>$A=\begin{bmatrix}0&1\\10&0\end{bmatrix}$（g/d=10），$B=\begin{bmatrix}0\\1\end{bmatrix}$，初始偏离 π/20。书中三组权重：</p>
<table class="tbl"><tr><th></th><th>Q</th><th>R</th><th>K</th><th>∫x₁²</th><th>∫x₂²</th><th>∫u²</th></tr>
<tr><td>测试一</td><td>diag(100,1)</td><td>1</td><td>[24.14, 7.02]</td><td>0.0080</td><td>0.0248</td><td>1.71</td></tr>
<tr><td>测试二</td><td>diag(1,100)</td><td>1</td><td>[20.05, 11.84]</td><td>0.0157</td><td>0.0105</td><td>1.92</td></tr>
<tr><td>测试三</td><td>diag(1,1)</td><td>1</td><td>[20.05, 6.41]</td><td>0.0099</td><td>0.0193</td><td>1.61</td></tr></table>
<p>规律：<b>哪个状态的权重大，哪个状态的积分就最小</b>；Q 相对 R 越小，能耗越低。（注：书中测试一的 K₂ 印作 7.203，按 Riccati 方程计算应为 7.020。）</p>
<div class="sim" data-id="cart"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>离散 LQR：$F=(B^TPB+R)^{-1}B^TPA$，$P\leftarrow(A-BF)^TP(A-BF)+F^TRF+Q$，从 $P=S$ 逆向递推。</li>
<li>连续 LQR：Riccati 方程，稳态为代数 Riccati 方程，$K=R^{-1}B^TP$。</li>
<li>两者都是全状态线性负反馈；权重矩阵是设计者表达“在意什么”的旋钮。</li>
</ol></div>`,
    sims: {
      riccati(el) {
        const ui = card(el, { title: "Riccati 递推：F 怎样收敛", hint: "一维系统 $x_{[k+1]}=a\\,x_{[k]}+b\\,u_{[k]}$。逐次迭代，观察 F 和闭环极点 a−bF。" });
        const pf = new Plot(ui.plotHost(), { height: 230, x: [0, 20], y: [0, 1.2], xlabel: "迭代次数（距终点步数）", title: "F 的收敛" });
        const sa = ui.slider({ label: "a（a>1 开环不稳定）", min: 0.5, max: 1.5, step: 0.01, value: 1, onInput: draw });
        const sb = ui.slider({ label: "b", min: 0.2, max: 2, step: 0.05, value: 1, onInput: draw });
        const sq = ui.slider({ label: "Q", min: 0, max: 10, step: 0.1, value: 1, onInput: draw });
        const sr = ui.slider({ label: "R", min: 0.05, max: 10, value: 1, log: true, onInput: draw });
        function draw() {
          const a = sa.get(), b = sb.get(), Q = sq.get(), R = sr.get(); let Pk = 1; const F = [];
          for (let k = 0; k < 20; k++) { const f = b * Pk * a / (b * b * Pk + R); F.push(f); Pk = (a - b * f) ** 2 * Pk + f * f * R + Q; }
          const mx = Math.max(1, ...F) * 1.15; pf.setRange(null, [0, mx]);
          pf.render = () => { pf.clear(); pf.line(F.map((_, i) => i + 1), F, COL.c5, 2.2); F.forEach((f, i) => pf.dot(i + 1, f, COL.c5, 3.5)); };
          pf.render();
          const Fi = F[F.length - 1];
          ui.out(`前几项：${F.slice(0, 6).map(v => f2(v, 5)).join(", ")}<br>收敛值 F = <b class="big">${f2(Fi, 5)}</b><br>闭环极点 a − bF = <b>${f2(a - b * Fi, 4)}</b>（|·|&lt;1 即稳定）<br><span class="sim-note">R 越大，F 越小——控制器越“懒”。</span>`);
        }
        draw();
      },
      cart(el) {
        const ui = card(el, { title: "平衡车 Q/R 调参实验台", hint: "动画：小车加速度 = 控制量。三条曲线是角度、角速度、控制量。下方柱状图是三项积分（与书中表 4.4.2 对比）。" });
        const sc = new Scene(ui.plotHost(), { height: 170 });
        const [a, b, c3] = ui.row("", "", ""); a.parentElement.classList.add("r3");
        const p1 = new Plot(a, { height: 180, x: [0, 4], y: [-0.05, 0.17], xlabel: "t", title: "角度 x₁" });
        const p2 = new Plot(b, { height: 180, x: [0, 4], y: [-0.3, 0.05], xlabel: "t", title: "角速度 x₂" });
        const p3 = new Plot(c3, { height: 180, x: [0, 4], y: [-4, 0.5], xlabel: "t", title: "u" });
        const sq1 = ui.slider({ label: "$q_1$", min: 0.1, max: 1000, value: 100, log: true, onInput: draw });
        const sq2 = ui.slider({ label: "$q_2$", min: 0.1, max: 1000, value: 1, log: true, onInput: draw });
        const sr = ui.slider({ label: "$R$", min: 0.01, max: 100, value: 1, log: true, onInput: draw });
        const set = (a, b, c) => { sq1.set(a); sq2.set(b); sr.set(c); draw(); };
        ui.buttons([["测试一", () => set(100, 1, 1), "ghost"], ["测试二", () => set(1, 100, 1), "ghost"], ["测试三", () => set(1, 1, 1), "ghost"]]);
        const A = [[0, 1], [10, 0]], B = [[0], [1]];
        let res, anim = 0;
        function draw() {
          const { K } = M.lqr(A, B, [[sq1.get(), 0], [0, sq2.get()]], [[sr.get()]]), k = K[0];
          const dt = 0.002; let x = [Math.PI / 20, 0], xi = 0, xv = 0; const ts = [], X1 = [], X2 = [], U = [], XI = []; let I = [0, 0, 0];
          for (let t = 0; t <= 4; t += dt) { const u = -(k[0] * x[0] + k[1] * x[1]); ts.push(t); X1.push(x[0]); X2.push(x[1]); U.push(u); XI.push(xi); I[0] += x[0] ** 2 * dt; I[1] += x[1] ** 2 * dt; I[2] += u * u * dt; xv += -u * dt; xi += xv * dt; x = rk4(s => [s[1], 10 * s[0] + u], x, t, dt); }
          res = { ts, X1, X2, U, XI }; anim = 0;
          [[p1, X1, COL.c1], [p2, X2, COL.c2], [p3, U, COL.c4]].forEach(([p, ys, c]) => { p.autoY([ys, [0]], 0.08); p.render = () => { p.clear(); p.line(ts, ys, c, 2.4); }; p.render(); });
          ui.out(`K = [<b>${f2(k[0], 4)}</b>, <b>${f2(k[1], 4)}</b>]<br>∫x₁² = <b>${f2(I[0], 4)}</b>　∫x₂² = <b>${f2(I[1], 4)}</b>　∫u² = <b>${f2(I[2], 2)}</b>`);
        }
        draw();
        CE.loop(ui.root, dt => {
          anim = (anim + dt * 0.8) % 5; const i = Math.min(res.ts.length - 1, Math.floor(anim / 0.002));
          const c = sc.ctx, W = sc.W, H = sc.H; sc.clear();
          const phi = res.X1[i], xi = res.XI[i], base = H - 30, cx = W / 2 + clamp(xi, -3, 3) * 50, L = 110;
          c.strokeStyle = COL.axis; c.beginPath(); c.moveTo(0, base + 12); c.lineTo(W, base + 12); c.stroke();
          c.fillStyle = COL.c2; c.beginPath(); c.roundRect(cx - 28, base - 10, 56, 20, 5); c.fill();
          c.fillStyle = COL.ink; [cx - 16, cx + 16].forEach(x => { c.beginPath(); c.arc(x, base + 10, 5, 0, 7); c.fill(); });
          c.strokeStyle = COL.ink2; c.lineWidth = 3; c.beginPath(); c.moveTo(cx, base - 10); c.lineTo(cx + L * Math.sin(phi), base - 10 - L * Math.cos(phi)); c.stroke();
          c.fillStyle = COL.c1; c.beginPath(); c.arc(cx + L * Math.sin(phi), base - 10 - L * Math.cos(phi), 11, 0, 7); c.fill();
          c.fillStyle = COL.muted; c.font = "12px Consolas"; c.fillText("t = " + res.ts[i].toFixed(2) + " s", 8, 16);
        });
      },
    },
  });

  /* ════════════════ 4.5–4.6 轨迹追踪 ════════════════ */
  const MSD = M.c2d([[0, 1], [-1, -0.5]], [[0], [1]], 0.1);
  // 三种追踪方法：返回仿真结果
  function trackSim(method, R, steps = 100, opts = {}) {
    const A = MSD.A, B = MSD.B, n = 2, I2 = M.eye(2), Q = I2, S = I2;
    const AD = opts.AD || I2;
    let Aa, Ba, Ca, xa, ud = 0;
    if (method === 1) { Aa = M.blkdiag(A, AD); Ba = M.vcat(B, M.zeros(2, 1)); Ca = M.hcat(I2, M.scale(I2, -1)); }
    else if (method === 2) { Aa = M.vcat(M.hcat(A, M.sub(I2, A)), M.hcat(M.zeros(2, 2), I2)); Ba = M.vcat(B, M.zeros(2, 1)); Ca = M.hcat(I2, M.scale(I2, -1)); const rhs = M.mulv(M.sub(I2, A), opts.xd || [1, 0]); ud = M.dot(M.flat(B), rhs) / M.dot(M.flat(B), M.flat(B)); }
    else { Aa = M.vcat(M.hcat(A, M.zeros(2, 2), B), M.hcat(M.zeros(2, 2), AD, M.zeros(2, 1)), [[0, 0, 0, 0, 1]]); Ba = M.vcat(B, M.zeros(2, 1), [[1]]); Ca = M.hcat(I2, M.scale(I2, -1), M.zeros(2, 1)); }
    const Qa = M.mul(M.mul(M.T(Ca), Q), Ca), Sa = M.mul(M.mul(M.T(Ca), S), Ca);
    const { K } = M.dlqr(Aa, Ba, Qa, [[R]], Sa);
    let x = [0, 0], xd = (opts.xd0 || [1, 0]).slice(), uPrev = 0; const out = { k: [], x1: [], x2: [], u: [], d1: [], d2: [], K: K[0], ud };
    for (let k = 0; k <= steps; k++) {
      if (opts.vref) xd[1] = opts.vref(k);
      const xaN = method === 3 ? [...x, ...xd, uPrev] : [...x, ...xd];
      let u = -M.dot(K[0], xaN);
      if (method === 2) u += ud; if (method === 3) u += uPrev;
      if (method === 0) u = opts.uOpen ?? 1;
      if (opts.sat) u = clamp(u, opts.sat[0], opts.sat[1]);
      out.k.push(k); out.x1.push(x[0]); out.x2.push(x[1]); out.u.push(u); out.d1.push(xd[0]); out.d2.push(xd[1]);
      x = M.vadd(M.mulv(A, x), M.flat(B).map(v => v * u)); xd = M.mulv(AD, xd); uPrev = u;
    }
    return out;
  }
  Course.add({
    id: "v2-track", vol: 2, ord: 80, no: "4.5–4.6", chName: "第4章 动态规划与 LQR", ref: "卷2 §4.5–4.7",
    title: "轨迹追踪：把“追目标”变成“误差归零”",
    sub: "LQR 只会把状态调到 0。想追一个非零目标或一条轨迹，就要改造状态和代价。三种方法，各有优缺点。",
    goals: ["理解三种追踪方法的增广矩阵构造", "理解方法 1 为什么有稳态误差、方法 2 怎样消除它、方法 3 为什么更平滑", "看到 LQR 处理约束的局限，引出 MPC"],
    minutes: 35,
    html: String.raw`
<p>实验台：弹簧–质量–阻尼系统（$m=1,k=1,b=0.5$，$T_s=0.1$ s 离散化），任务是从 $[0,0]$ 移到 $\bm x_d=[1,0]^T$——停在 1 m 处。思路都是：<b>把参考值也放进状态（增广），把代价改成误差的代价</b>，然后套用标准 LQR 求解器 [F1]，一行都不用改。</p>
<table class="tbl"><tr><th></th><th>方法 1：引入误差</th><th>方法 2：稳态非零参考</th><th>方法 3：输入增量</th></tr>
<tr><td>增广状态</td><td>$[\bm x;\bm x_d]$</td><td>$[\bm x;\bm x_d]$</td><td>$[\bm x;\bm x_d;\bm u_{[k-1]}]$</td></tr>
<tr><td>$A_a$</td><td>$\begin{bmatrix}A&0\\0&A_D\end{bmatrix}$</td><td>$\begin{bmatrix}A&I-A\\0&I\end{bmatrix}$</td><td>$\begin{bmatrix}A&0&B\\0&A_D&0\\0&0&I\end{bmatrix}$</td></tr>
<tr><td>$B_a$</td><td>$\begin{bmatrix}B\\0\end{bmatrix}$</td><td>$\begin{bmatrix}B\\0\end{bmatrix}$</td><td>$\begin{bmatrix}B\\0\\I\end{bmatrix}$</td></tr>
<tr><td>输入代价</td><td>$\bm u^TR\bm u$</td><td>$\delta\bm u^TR\,\delta\bm u$，$\delta\bm u=\bm u-\bm u_d$</td><td>$\Delta\bm u^TR\,\Delta\bm u$，$\Delta\bm u=\bm u_{[k]}-\bm u_{[k-1]}$</td></tr>
<tr><td>特点</td><td>可能有稳态误差</td><td>无稳态误差，需目标为常数</td><td>更平滑，可追踪变化的目标</td></tr></table>
<div class="eq" data-title="三种方法共用的误差与权重">
$$\bm e=\bm x-\bm x_d=C_a\bm x_a,\quad C_a=[I\;\;-I\;\;(0)],\qquad Q_a=C_a^TQC_a,\quad S_a=C_a^TSC_a$$
</div>

<h3>方法 1 的毛病</h3>
<p>要让质量块停在 1 m，就必须<b>持续施加拉力</b>抵抗弹簧。但代价里的 $\bm u^TR\bm u$ 希望输入越小越好——两者矛盾。结果是系统“为了省力”停在半路：R=0.1 时停在约 0.91，R=1 时只到 0.5。</p>
<h3>方法 2：先算出“维持目标所需的输入”</h3>
<div class="eq" data-no="(4.5.10)–(4.5.11)">
$$\bm x_d=A\bm x_d+B\bm u_d\;\Rightarrow\;B\bm u_d=(I-A)\bm x_d$$
<div class="say">本例解得 $u_d=1$：1 N 的力正好抵消弹簧 $k\times1$ m 的弹力。代价改为惩罚“偏离 $u_d$ 的程度” $\delta u$，矛盾就消失了。注意 $u_d$ 要解线性方程，B 不一定可逆。</div>
</div>
<h3>方法 3：惩罚输入的变化</h3>
<p>把上一步的输入也放进状态，代价惩罚 $\Delta u$：输入不会一开始就跳到最大，而是缓缓上升、再回落，状态变化也更平缓、超调更小。而且它能追踪<b>变化</b>的目标，比如匀速运动：$A_D=\begin{bmatrix}1&0.1\\0&1\end{bmatrix}$。</p>
<div class="sim" data-id="three"></div>
<div class="sim" data-id="vtrack"></div>

<h2><span class="sec">4.6</span>无人机高度控制：LQR 的边界</h2>
<p>无人机增广模型（第 3 章）用方法 2：$B u_d=(I-A)\bm x_d$ 解得 $u_d=g=10$——悬停时推力等于重力，<b>第 3 章的“悬停矛盾”就此化解</b>。权重 $Q=S=\mathrm{diag}(1,1,0)$（$x_3=-g$ 是常数，不需要也无法控制）。</p>
<p>无约束时，测试一（R=0.1）开局推力高达约 36（按 Riccati 方程计算）——电机根本做不到。加上饱和函数 $u\in[7,12]$（对应净加速度 −3～2），上升变慢、最高速度降低。但如果还想限制<b>速度</b> ≤ 3 m/s 呢？LQR 推导里根本没有约束，用饱和函数去“剪”状态也不合理——这就是第 5 章 <b>MPC</b> 要解决的问题。</p>
<div class="sim" data-id="drone"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>追踪问题 → 增广状态 + 误差代价 → 标准 LQR。</li>
<li>方法 1 因“省力”与“维持目标”矛盾而有稳态误差；方法 2 用 $\delta u=u-u_d$ 消除它；方法 3 用 $\Delta u$ 让输入更平滑并能追踪变化目标。</li>
<li>LQR 不能处理约束（尤其状态约束），饱和函数只是权宜之计 → MPC。</li>
</ol></div>`,
    sims: {
      three(el) {
        const ui = card(el, { title: "三种追踪方法 PK（弹簧–质量–阻尼）", hint: "目标 $x_d=[1,0]$，S=Q=I，100 步（10 s）。上：位移 x₁；下：输入 u。" });
        const sc = new Scene(ui.plotHost(), { height: 90 });
        const p1 = new Plot(ui.plotHost(), { height: 220, x: [0, 100], y: [0, 1.5], xlabel: "k", ylabel: "x₁" });
        const pu = new Plot(ui.plotHost(), { height: 180, x: [0, 100], y: [-0.5, 3], xlabel: "k", ylabel: "u" });
        const sR = ui.slider({ label: "R", min: 0.01, max: 10, value: 0.1, log: true, onInput: draw });
        const shows = [ui.check({ label: "开环常数输入 u=1", value: false, onChange: draw }), ui.check({ label: "方法 1：引入误差", value: true, onChange: draw }), ui.check({ label: "方法 2：稳态非零参考", value: true, onChange: draw }), ui.check({ label: "方法 3：输入增量", value: true, onChange: draw })];
        const cols = [COL.muted, COL.c4, COL.c1, COL.c3];
        let runs = [], t = 0;
        function draw() {
          const R = sR.get(); runs = [0, 1, 2, 3].map(m => shows[m].get() ? trackSim(m, R, m === 0 ? 100 : 100) : null);
          p1.render = () => { p1.clear(); p1.hline(1, COL.c2, [6, 4]); runs.forEach((r, m) => r && p1.line(r.k, r.x1, cols[m], m ? 2.4 : 1.6, m ? null : [4, 3])); };
          pu.autoY(runs.filter(Boolean).map(r => r.u), 0.1, true);
          pu.render = () => { pu.clear(); runs.forEach((r, m) => r && pu.line(r.k, r.u, cols[m], m ? 2.2 : 1.4)); };
          p1.render(); pu.render();
          ui.out(runs.map((r, m) => r ? `<span style="color:${cols[m]}">●</span> ${["开环", "方法1", "方法2", "方法3"][m]}：终点 x₁ = <b>${f2(r.x1[r.x1.length - 1], 3)}</b>，u 峰值 ${f2(Math.max(...r.u), 2)}` : "").filter(Boolean).join("<br>"));
        }
        draw();
        CE.loop(ui.root, dt => {
          t = (t + dt * 12) % 110; const k = Math.min(100, Math.floor(t)); const c = sc.ctx, W = sc.W; sc.clear();
          const x2p = x => 40 + x * (W - 120) / 1.4;
          c.fillStyle = COL.axis; c.fillRect(20, 10, 8, 70); c.strokeStyle = COL.c3; c.setLineDash([4, 3]); c.beginPath(); c.moveTo(x2p(1), 5); c.lineTo(x2p(1), 85); c.stroke(); c.setLineDash([]);
          runs.forEach((r, m) => { if (!r || m === 0) return; const y = 12 + (m - 1) * 25, X = x2p(r.x1[k]); c.strokeStyle = cols[m]; c.lineWidth = 1.5; c.beginPath(); c.moveTo(28, y + 9); for (let i = 1; i < 10; i++) c.lineTo(28 + (X - 28) * i / 10, y + 9 + (i % 2 ? -5 : 5)); c.lineTo(X, y + 9); c.stroke(); c.fillStyle = cols[m]; c.fillRect(X, y, 22, 18); });
        });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>把 R 从 0.01 调到 10：方法 1 的终点越来越偏离 1（R=1 时约 0.5）。</li><li>R 很大时，方法 2 退化成“开环常数输入”——慢、有振荡（勾选开环对比）。</li><li>比较方法 2 与方法 3 的输入曲线：方法 3 从 0 平滑爬升，没有开局的跳变。</li></ol></div>`);
      },
      vtrack(el) {
        const ui = card(el, { title: "方法 3 追踪匀速运动：来回走", hint: "目标速度 ±0.2 m/s，每 5 s 换一次方向（目标位置是三角波），共 20 s。$A_D=\\begin{bmatrix}1&0.1\\\\0&1\\end{bmatrix}$。" });
        const p1 = new Plot(ui.plotHost(), { height: 200, x: [0, 200], y: [-0.1, 1.2], xlabel: "k", ylabel: "x₁" });
        const p2 = new Plot(ui.plotHost(), { height: 160, x: [0, 200], y: [-0.4, 0.4], xlabel: "k", ylabel: "x₂" });
        const pu = new Plot(ui.plotHost(), { height: 150, x: [0, 200], y: [-0.5, 1.6], xlabel: "k", ylabel: "u" });
        const sR = ui.slider({ label: "R", min: 0.01, max: 10, value: 0.1, log: true, onInput: draw });
        const sV = ui.slider({ label: "目标速度幅值 (m/s)", min: 0.05, max: 0.5, step: 0.01, value: 0.2, onInput: draw });
        const sP = ui.slider({ label: "换向周期 (步)", min: 20, max: 100, step: 5, value: 50, fmt: v => v.toFixed(0), onInput: draw });
        function draw() {
          const v = sV.get(), T = sP.get();
          const r = trackSim(3, sR.get(), 200, { AD: [[1, 0.1], [0, 1]], xd0: [0, v], vref: k => (Math.floor(k / T) % 2 ? -v : v) });
          p1.autoY([r.x1, r.d1], 0.1);
          p1.render = () => { p1.clear(); p1.line(r.k, r.d1, COL.c2, 2, [6, 4]); p1.line(r.k, r.x1, COL.c1, 2.4); };
          p2.setRange(null, [-v * 1.8, v * 1.8]);
          p2.render = () => { p2.clear(); p2.line(r.k, r.d2, COL.c2, 2, [6, 4]); p2.line(r.k, r.x2, COL.c1, 2.2); };
          pu.autoY([r.u], 0.1, true); pu.render = () => { pu.clear(); pu.line(r.k, r.u, COL.c4, 2); };
          p1.render(); p2.render(); pu.render();
          ui.out(`大多数时间 u &gt; 0：质量块在弹簧平衡点一侧运动，需要持续抵抗弹簧力（u ≈ kx₁ + bx₂）。<br>R 小 → 换向时输入尖峰更剧烈但跟得更紧。`);
        }
        draw();
      },
      drone(el) {
        const ui = card(el, { title: "无人机悬停到 10 m：LQR + 饱和", hint: "方法 2，$u_d=g=10$，$Q=S=\\mathrm{diag}(1,1,0)$。阴影是推力的可行范围 [7, 12]。" });
        const [a, b, c3] = ui.row("", "", ""); a.parentElement.classList.add("r3");
        const p1 = new Plot(a, { height: 210, x: [0, 100], y: [0, 12], xlabel: "k", title: "高度 x₁" });
        const p2 = new Plot(b, { height: 210, x: [0, 100], y: [-1, 7], xlabel: "k", title: "速度 x₂" });
        const pu = new Plot(c3, { height: 210, x: [0, 100], y: [5, 28], xlabel: "k", title: "推力 u" });
        const sR = ui.slider({ label: "R", min: 0.01, max: 10, value: 0.1, log: true, onInput: draw });
        const cS = ui.check({ label: "加饱和 u ∈ [7, 12]", value: true, onChange: draw });
        function sim(sat) {
          const A = [[1, 0.1, 0.005], [0, 1, 0.1], [0, 0, 1]], B = [[0.005], [0.1], [0]], I3 = M.eye(3);
          const Aa = M.vcat(M.hcat(A, M.sub(I3, A)), M.hcat(M.zeros(3, 3), I3)), Ba = M.vcat(B, M.zeros(3, 1)), Ca = M.hcat(I3, M.scale(I3, -1));
          const W = [[1, 0, 0], [0, 1, 0], [0, 0, 0]], Qa = M.mul(M.mul(M.T(Ca), W), Ca);
          const { K } = M.dlqr(Aa, Ba, Qa, [[sR.get()]], Qa);
          let x = [0, 0, -10]; const xd = [10, 0, -10], o = { k: [], h: [], v: [], u: [] };
          for (let k = 0; k <= 100; k++) { let u = 10 - M.dot(K[0], [...x, ...xd]); if (sat) u = clamp(u, 7, 12); o.k.push(k); o.h.push(x[0]); o.v.push(x[1]); o.u.push(u); x = M.vadd(M.mulv(A, x), M.flat(B).map(v => v * u)); }
          return o;
        }
        function draw() {
          const free = sim(false), s = cS.get() ? sim(true) : null;
          p1.render = () => { p1.clear(); p1.hline(10, COL.c2, [6, 4]); p1.line(free.k, free.h, COL.muted, 1.6, [4, 3]); if (s) p1.line(s.k, s.h, COL.c1, 2.4); };
          p2.render = () => { p2.clear(); p2.hline(3, COL.c4, [3, 3]); p2.text(98, 3.3, "期望的速度上限 3", COL.c4, "right"); p2.line(free.k, free.v, COL.muted, 1.6, [4, 3]); if (s) p2.line(s.k, s.v, COL.c1, 2.4); };
          pu.setRange(null, [5, Math.max(14, Math.max(...free.u) * 1.05)]);
          pu.render = () => { pu.clear(); pu.rect(0, 7, 100, 12, COL.c3 + "22"); pu.line(free.k, free.u, COL.muted, 1.6, [4, 3]); if (s) pu.line(s.k, s.u, COL.c4, 2.4); };
          p1.render(); p2.render(); pu.render();
          ui.out(`无约束（虚线）：开局推力 <b>${f2(free.u[0], 1)}</b>，最高速度 ${f2(Math.max(...free.v), 2)} m/s${s ? `<br>加饱和（实线）：最高速度 <b>${f2(Math.max(...s.v), 2)}</b> m/s` : ""}<br><span class="bad">速度仍超过 3 m/s：LQR 无法直接处理状态约束 → MPC</span>`);
        }
        draw();
      },
    },
  });
})();
