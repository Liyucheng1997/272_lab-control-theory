/* 卷2 第5章：模型预测控制 MPC */
(function () {
  const { M, P, Plot, Scene, card, linspace, logspace, f2, COL, clamp, rk4, cfmt, mtex, mpcMatrices, mpcConstraints, qpPrep } = CE;
  const INF = Infinity;
  const droneAB = Ts => ({ A: [[1, Ts, Ts * Ts / 2], [0, 1, Ts], [0, 0, 1]], B: [[Ts * Ts / 2], [Ts], [0]] });
  // 稳态非零参考增广（[F2]）
  function augSS(A, B, xd, Q, S) {
    const n = A.length, I = M.eye(n);
    const Aa = M.vcat(M.hcat(A, M.sub(I, A)), M.hcat(M.zeros(n, n), I)), Ba = M.vcat(B, M.zeros(n, B[0].length));
    const Ca = M.hcat(I, M.scale(I, -1)), rhs = M.mulv(M.sub(I, A), xd), bb = M.flat(B);
    const ud = M.dot(bb, rhs) / M.dot(bb, bb);
    return { Aa, Ba, Qa: M.mul(M.mul(M.T(Ca), Q), Ca), Sa: M.mul(M.mul(M.T(Ca), S), Ca), ud };
  }

  /* ════════════════ 5.1–5.2 滚动优化与二次规划 ════════════════ */
  Course.add({
    id: "v2-mpc1", vol: 2, ord: 90, no: "5.1–5.2", chName: "第5章 模型预测控制", ref: "卷2 §5.1–5.2",
    title: "MPC 的基本思想：向前看几步，只走一步",
    sub: "每个时刻都用模型预测未来、规划一整串最优控制，但只执行第一步；下一时刻重新测量、重新规划。求解的核心是二次规划。",
    goals: ["理解预测区间 Np、控制区间 Nc 与滚动优化", "理解为什么“只执行第一步”带来了反馈", "掌握二次规划的标准形式，以及无约束、等式约束、不等式约束三种情况"],
    minutes: 30,
    html: String.raw`
<p><span class="kw">模型预测控制（MPC）</span>一句话：用系统<b>当前的状态</b>和约束，<b>预测</b>未来一段时间内的状态和输入，解出一组最优控制序列，<b>只执行第一个</b>；下一时刻用新测到的状态重新算一遍，如此循环。它相比 LQR 最大的卖点是：<b>约束可以直接放进优化问题</b>，所以在工业界应用极广。</p>
<div data-fig="mpcLoop"></div>
<h2><span class="sec">5.1</span>滚动优化</h2>
<ul><li><b>预测区间 $N_p$</b>：每次向前预测多少步。</li>
<li><b>控制区间 $N_c$</b>：其中有多少个控制量可以自由变化。$N_c<N_p$ 时，超出控制区间的控制量保持不变（节省计算）。本书默认 $N_p=N_c$。</li>
<li><b>滚动优化（receding horizon）</b>：预测窗口随时间一格一格往前“滚”。每次只用计划的第一步，其余全部丢弃——因为模型有误差、有扰动，旧计划很快就过时了。</li>
<li><b>反馈从哪来？</b>每次优化都以<b>实际测得</b>的当前状态为起点——这就是反馈。</li></ul>
<div class="sim" data-id="receding"></div>
<div class="idea" data-t="下棋的智慧">高手下棋会算好几步，但每次只落一子；对手走完后，再根据新局面重新计算。MPC 就是这样：远见 + 随机应变。</div>

<h2><span class="sec">5.2</span>二次规划（QP）</h2>
<p>线性系统 + 二次型代价的 MPC，每一步都可以化成一个<b>二次规划</b>问题。先学会 QP，后面的工作就是把 MPC 问题“翻译”成 QP 的标准形式。</p>
<div class="eq" data-no="(5.2.1)" data-title="二次规划标准形式">
$$\min_{\bm u}\;J=\tfrac12\bm u^TH\bm u+\bm u^T\bm f\qquad \text{s.t.}\;\;M\bm u\le\bm b,\quad M_{eq}\bm u=\bm b_{eq},\quad \mathrm{LB}\le\bm u\le\mathrm{UB}$$
<dl><dt>$H$</dt><dd>对称正定：代价函数是一个“碗”</dd><dt>$\bm f$</dt><dd>线性项：决定碗底偏向哪边</dd><dt>约束</dt><dd>不等式、等式、上下界</dd></dl>
</div>
<table class="tbl"><tr><th>情况</th><th>解法</th><th>书中例题</th></tr>
<tr><td>无约束</td><td>$\frac{\partial J}{\partial\bm u}=H\bm u+\bm f=0\Rightarrow\bm u^*=-H^{-1}\bm f$</td><td>$J=\frac12(u_1^2+u_2^2)+u_1+u_2$：$\bm u^*=(-1,-1)$，$J=-1$</td></tr>
<tr><td>等式约束</td><td>拉格朗日乘数：$\begin{bmatrix}H&M_{eq}^T\\M_{eq}&0\end{bmatrix}\begin{bmatrix}\bm u\\\bm\lambda\end{bmatrix}=\begin{bmatrix}-\bm f\\\bm b_{eq}\end{bmatrix}$</td><td>加上 $u_1-u_2=1$：$\bm u^*=(-0.5,-1.5)$，$\lambda=-0.5$</td></tr>
<tr><td>不等式约束</td><td>数值方法（内点法、梯度法……）或软件 <code>quadprog</code></td><td>$u_1+u_2\le2$，$-u_1+u_2\le1$，$0\le u_1\le1$，$0\le u_2\le2$：$\bm u^*=(0,0)$</td></tr></table>
<div class="sim" data-id="qpbowl"></div>
<div class="idea" data-t="等式约束的几何意义">最优点处，代价函数的等高线与约束直线<b>相切</b>——代价的梯度 $H\bm u+\bm f$ 与约束的梯度 $M_{eq}^T$ 共线，λ 就是两者的比例系数。</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>MPC：预测 $N_p$ 步、解最优序列、只执行第一步、下一时刻重来（滚动优化）。</li>
<li>反馈来自每次用实测状态作为预测起点。</li>
<li>QP：$\min\frac12\bm u^TH\bm u+\bm u^T\bm f$；无约束解 $-H^{-1}\bm f$，等式约束用拉格朗日乘数，不等式约束用数值算法。</li>
</ol></div>`,
    sims: {
      receding(el) {
        const ui = card(el, { title: "滚动窗口：计划总在变，执行的只有第一步", hint: "小车（双积分器，$T_s=0.2$ s）要停到位置 1。半透明区域是预测窗口；虚线是预测状态，柱子是计划的控制量，<b>只有第一根（实心）会被执行</b>。" });
        const p1 = new Plot(ui.plotHost(), { height: 230, x: [0, 40], y: [-0.3, 1.5], xlabel: "k", ylabel: "位置" });
        const pu = new Plot(ui.plotHost(), { height: 170, x: [0, 40], y: [-1.6, 1.6], xlabel: "k", ylabel: "控制量 u" });
        const sN = ui.slider({ label: "预测区间 Np", min: 1, max: 25, step: 1, value: 10, fmt: v => v.toFixed(0), onInput: rebuild });
        const sC = ui.slider({ label: "控制区间 Nc", min: 1, max: 25, step: 1, value: 10, fmt: v => v.toFixed(0), onInput: rebuild });
        const cU = ui.check({ label: "输入约束 |u| ≤ 1", value: true, onChange: rebuild });
        const cD = ui.check({ label: "随机扰动（模型误差）", value: false });
        let auto = false;
        ui.buttons([["⏭ 下一步", () => stepOnce()], ["▶ 自动 / 暂停", () => auto = !auto], ["↺", reset, "ghost"]]);
        const Ts = 0.2, A = [[1, Ts], [0, 1]], B = [[Ts * Ts / 2], [Ts]], Q = [[10, 0], [0, 1]], R = [[0.1]], xd = [1, 0];
        let sol, mm, cons, x, H, plan, k, rg = CE.rng(4), Tm;
        function rebuild() {
          const Np = sN.get(), Nc = Math.min(sC.get(), Np); if (sC.get() > Np) sC.set(Np);
          mm = mpcMatrices(A, B, Q, R, Q, Np);
          // 控制区间：U = T·Uc，超出 Nc 的输入等于最后一个
          Tm = M.zeros(Np, Nc); for (let i = 0; i < Np; i++) Tm[i][Math.min(i, Nc - 1)] = 1;
          const Hc = M.sym(M.mul(M.mul(M.T(Tm), mm.H), Tm)), Fc = M.mul(M.T(Tm), mm.F);
          let Ain = []; if (cU.get()) for (let i = 0; i < Nc; i++) { const r1 = new Array(Nc).fill(0), r2 = new Array(Nc).fill(0); r1[i] = 1; r2[i] = -1; Ain.push(r1, r2); }
          mm.Hc = Hc; mm.Fc = Fc; sol = qpPrep(Hc, Ain); mm.nA = Ain.length;
          plan = solve();
          draw();
        }
        function solve() {
          const e = M.vsub(x, xd), f = M.mulv(mm.Fc, e), b = new Array(mm.nA).fill(1);
          const Uc = sol(f, b).x, U = M.mulv(Tm, Uc);
          // 预测
          const pred = [x.slice()]; let z = e.slice(); U.forEach(u => { z = M.vadd(M.mulv(A, z), M.flat(B).map(v => v * u)); pred.push(M.vadd(z, xd)); });
          return { U, pred };
        }
        function reset() { x = [0, 0]; H = [{ k: 0, x: 0, u: null }]; k = 0; rebuild(); }
        function stepOnce() {
          if (k >= 40) return;
          const u = plan.U[0]; x = M.vadd(M.mulv(A, x), M.flat(B).map(v => v * u));
          if (cD.get()) x = [x[0] + 0.04 * rg.n(), x[1] + 0.15 * rg.n()];
          H[H.length - 1].u = u; k++; H.push({ k, x: x[0], u: null });
          plan = solve(); draw();
        }
        function draw() {
          const Np = sN.get();
          p1.render = () => {
            p1.clear(); p1.rect(k, -0.3, k + Np, 1.5, COL.c5 + "14"); p1.hline(1, COL.c3, [6, 4]);
            p1.line(plan.pred.map((_, i) => k + i), plan.pred.map(p => p[0]), COL.c5, 2, [5, 4]); plan.pred.forEach((p, i) => p1.dot(k + i, p[0], COL.c5, 3));
            p1.line(H.map(h => h.k), H.map(h => h.x), COL.c1, 2.6); p1.dot(k, x[0], COL.c4, 6);
          };
          pu.render = () => {
            pu.clear(); pu.rect(k, -1.6, k + Np, 1.6, COL.c5 + "14"); if (cU.get()) { pu.hline(1, COL.c4, [3, 3]); pu.hline(-1, COL.c4, [3, 3]); }
            plan.U.forEach((u, i) => pu.rect(k + i + 0.1, Math.min(0, u), k + i + 0.9, Math.max(0, u), i === 0 ? COL.c2 : COL.c2 + "40", COL.c2));
            H.forEach(h => { if (h.u !== null) pu.rect(h.k + 0.1, Math.min(0, h.u), h.k + 0.9, Math.max(0, h.u), COL.c1 + "aa"); });
          };
          p1.render(); pu.render();
          ui.out(`k = <b>${k}</b><br>本步执行 u = <b>${f2(plan.U[0], 3)}</b><br>计划序列：${plan.U.slice(0, 6).map(v => f2(v, 2)).join(", ")}${plan.U.length > 6 ? "…" : ""}<br><span class="sim-note">每走一步，紫色窗口右移一格，计划被重新计算。</span>`);
        }
        reset();
        let acc = 0; CE.loop(ui.root, dt => { if (!auto) return; acc += dt; if (acc > 0.25) { acc = 0; stepOnce(); } });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>Np=1：只看一步，“目光短浅”，走得又慢又犹豫；Np=10 以上就很好了。</li><li>Nc 调小：计划中后面的控制量被迫保持不变，自由度减少。</li><li>打开扰动：实际状态偏离预测，但下一步的新计划马上修正——这就是滚动优化的反馈作用。</li></ol></div>`);
      },
      qpbowl(el) {
        const ui = card(el, { title: "QP 的碗与可行域", hint: "$J=\\frac12(u_1^2+u_2^2)+u_1+u_2$ 的等高线（同心圆，碗底在 (−1,−1)）。切换约束类型，观察最优点（红点）跑到哪里。" });
        const pl = new Plot(ui.plotHost(), { height: 380, x: [-2.5, 2.5], y: [-2.5, 2.5], equal: true, xlabel: "u₁", ylabel: "u₂" });
        let mode = "free";
        ui.seg({ label: "约束", options: [["free", "无约束"], ["eq", "等式 u₁−u₂=c"], ["ineq", "不等式（例 5.2.3）"]], value: mode, onChange: v => { mode = v; draw(); } });
        const sc = ui.slider({ label: "等式约束中的 c", min: -2, max: 2, step: 0.05, value: 1, onInput: draw });
        const cons = [["u₁+u₂ ≤ 2", [1, 1], 2], ["−u₁+u₂ ≤ 1", [-1, 1], 1], ["u₁ ≥ 0", [-1, 0], 0], ["u₁ ≤ 1", [1, 0], 1], ["u₂ ≥ 0", [0, -1], 0], ["u₂ ≤ 2", [0, 1], 2]];
        const on = cons.map((c, i) => ui.check({ label: c[0], value: true, onChange: draw }));
        const H = [[1, 0], [0, 1]], f = [1, 1], J = u => 0.5 * (u[0] ** 2 + u[1] ** 2) + u[0] + u[1];
        function draw() {
          let u, extra = "";
          if (mode === "free") u = [-1, -1];
          else if (mode === "eq") { const c = sc.get(); const K = [[1, 0, 1], [0, 1, -1], [1, -1, 0]], sol = M.solve(K, [-1, -1, c]); u = [sol[0], sol[1]]; extra = `λ = ${f2(sol[2], 3)}`; }
          else { const act = cons.filter((_, i) => on[i].get()); u = CE.qp(H, f, act.map(c => c[1]), act.map(c => c[2]), 2000).x; }
          pl.render = () => {
            pl.clear();
            [0.25, 0.5, 1, 1.5, 2.2, 3, 4].forEach(r => { const th = linspace(0, 6.3, 80); pl.line(th.map(t => -1 + r * Math.cos(t)), th.map(t => -1 + r * Math.sin(t)), COL.c1 + "55", 1.2); });
            pl.dot(-1, -1, COL.c1, 4);
            if (mode === "ineq") {
              // 可行域：网格填充
              const c = pl.ctx; c.fillStyle = COL.c3 + "30";
              for (let x = -2.5; x <= 2.5; x += 0.04) for (let y = -2.5; y <= 2.5; y += 0.04) if (cons.every((cc, i) => !on[i].get() || cc[1][0] * x + cc[1][1] * y <= cc[2] + 1e-9)) c.fillRect(pl.fx(x), pl.fy(y + 0.04), Math.ceil(pl.fx(x + 0.04) - pl.fx(x)) + 1, Math.ceil(pl.fy(y) - pl.fy(y + 0.04)) + 1);
              cons.forEach((cc, i) => { if (!on[i].get()) return; const [a, b] = cc[1], r = cc[2]; if (Math.abs(b) > 1e-9) pl.fn(x => (r - a * x) / b, COL.c3, 1.4); else pl.vline(r / a, COL.c3, null, 1.4); });
            }
            if (mode === "eq") { const c = sc.get(); pl.fn(x => x - c, COL.c2, 2); const g = [u[0] + 1, u[1] + 1]; pl.arrow(u[0], u[1], u[0] + g[0] * 0.8, u[1] + g[1] * 0.8, COL.c4, 2, 8); pl.arrow(u[0], u[1], u[0] + 0.5, u[1] - 0.5, COL.c2, 2, 8); }
            pl.dot(u[0], u[1], COL.c4, 7, COL.plotbg);
          };
          pl.render();
          ui.out(`最优点 u* = (<b>${f2(u[0], 3)}</b>, <b>${f2(u[1], 3)}</b>)<br>J* = <b>${f2(J(u), 3)}</b>${extra ? "<br>" + extra : ""}${mode === "eq" ? '<br><span class="sim-note">红箭头：代价梯度；橙箭头：约束法向。二者共线。</span>' : ""}${mode === "ineq" ? '<br><span class="sim-note">绿色区域是可行域。无约束碗底 (−1,−1) 在可行域外，最优点落在边界/角点上。</span>' : ""}`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 5.3 MPC 推导 ════════════════ */
  Course.add({
    id: "v2-mpc2", vol: 2, ord: 100, no: "5.3", chName: "第5章 模型预测控制", ref: "卷2 §5.3",
    title: "MPC 推导：把未来“堆”成一个向量",
    sub: "把预测区间内的所有状态和输入摞成长向量，动态方程就变成一个矩阵等式，代价函数变成标准 QP。无约束时，它和 LQR 殊途同归。",
    goals: ["会构造 $\\Phi,\\Gamma,\\Omega,\\Psi$ 并推出 $H=\\Gamma^T\\Omega\\Gamma+\\Psi$、$F=\\Gamma^T\\Omega\\Phi$", "理解无约束 MPC 等价于一个线性反馈 $u=-k_{mpc}x$，且 Np 足够大时等于 LQR", "理解为什么必须在线滚动计算（反馈特性）"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">5.3.1</span>把预测写成矩阵</h2>
<p>系统 $\bm x_{[k+1]}=A\bm x_{[k]}+B\bm u_{[k]}$。从当前状态 $\bm x_{[k|k]}$ 出发一步步往前推：$\bm x_{[k+2|k]}=A^2\bm x_{[k|k]}+AB\bm u_{[k|k]}+B\bm u_{[k+1|k]}$……每个预测状态只取决于当前状态和未来的输入序列。把它们堆起来：</p>
<div class="eq" data-no="(5.3.4)–(5.3.5)" data-title="预测方程">
$$\underbrace{\begin{bmatrix}\bm x_{[k+1|k]}\\\bm x_{[k+2|k]}\\\vdots\\\bm x_{[k+N_p|k]}\end{bmatrix}}_{X_{[k]}}=\underbrace{\begin{bmatrix}A\\A^2\\\vdots\\A^{N_p}\end{bmatrix}}_{\Phi}\bm x_{[k|k]}+\underbrace{\begin{bmatrix}B&0&\cdots&0\\AB&B&\cdots&0\\\vdots&&\ddots&\vdots\\A^{N_p-1}B&A^{N_p-2}B&\cdots&B\end{bmatrix}}_{\Gamma}\underbrace{\begin{bmatrix}\bm u_{[k|k]}\\\bm u_{[k+1|k]}\\\vdots\\\bm u_{[k+N_p-1|k]}\end{bmatrix}}_{U_{[k]}}$$
<dl><dt>$\bm x_{[k+i|k]}$</dt><dd>在 k 时刻<b>预测</b>的 k+i 时刻状态</dd><dt>$\Gamma$</dt><dd>下三角块矩阵：第 i 行第 j 列的块是 $A^{i-j}B$（j≤i）</dd></dl>
<div class="say">注意 X 从 k+1 开始，U 从 k 开始。LQR 用逆向递推逐个求控制量；MPC 把整串控制量当成一个向量<b>一次性</b>求出来。</div>
</div>

<h2><span class="sec">5.3.2</span>把代价写成 QP</h2>
<p>代价 $J=\frac12\bm x_{[k+N_p]}^TS\bm x_{[k+N_p]}+\frac12\sum(\bm x^TQ\bm x+\bm u^TR\bm u)$。把当前状态那一项（与 U 无关）拿出来，其余写成块对角矩阵：</p>
<div class="eq" data-no="(5.3.8)–(5.3.14)" data-title="MPC 的 QP 形式">
$$\Omega=\mathrm{blkdiag}(Q,\dots,Q,\,S),\quad \Psi=\mathrm{blkdiag}(R,\dots,R)\quad\Longrightarrow\quad J=U^T\underbrace{\Gamma^T\Omega\Phi}_{F}\,\bm x_{[k|k]}+\tfrac12U^T\underbrace{(\Gamma^T\Omega\Gamma+\Psi)}_{H}\,U$$
<div class="say">与 QP 标准形式对照：决策变量是 U，线性项 $\bm f=F\bm x_{[k|k]}$——<b>当前状态出现在线性项里，这就是反馈</b>。Ω 的最后一块是 S，因为 X 的最后一项是末端状态。</div>
</div>
<details class="derive"><summary>代入展开</summary><div class="body">
<div class="step">$J=\frac12\bm x^TQ\bm x+\frac12X^T\Omega X+\frac12U^T\Psi U$，代入 $X=\Phi\bm x+\Gamma U$。</div>
<div class="step">展开得两个交叉项 $\frac12\bm x^T\Phi^T\Omega\Gamma U$ 与 $\frac12U^T\Gamma^T\Omega\Phi\bm x$，它们互为转置的标量，相等，合并为 $U^T\Gamma^T\Omega\Phi\bm x$。</div>
<div class="step">$\frac12\bm x^T(Q+\Phi^T\Omega\Phi)\bm x$ 与 U 无关，优化时忽略。</div>
</div></details>
<div class="sim" data-id="matrices"></div>

<h2><span class="sec">5.3.3</span>无约束解析解</h2>
<div class="eq" data-no="(5.3.16)–(5.3.18)">
$$U^*=-H^{-1}F\bm x_{[k|k]},\qquad \bm u^*_{[k|k]}=-\underbrace{[I\;0\;\cdots\;0]H^{-1}F}_{k_{mpc}}\bm x_{[k|k]}$$
<div class="say">无约束 MPC 本质上就是一个线性状态反馈，形式和 LQR 一样。</div>
</div>
<p>一维例子（A=B=Q=R=S=1，x=1）：$N_p=2$ 时 $H=\begin{bmatrix}3&1\\1&2\end{bmatrix}$，$F=\begin{bmatrix}2\\1\end{bmatrix}$，$U^*=-[0.6,\;0.2]^T$——第一项 −0.6 与 LQR（N=2）完全一致。$N_p$ 增大，第一项收敛到 −0.61803：<b>LQR 与 MPC 求解的是同一个最优化问题</b>，只是 LQR 从末端逆向递推、得到反馈矩阵；MPC 正向预测、用 QP 求一整串控制量。</p>
<div class="sim" data-id="vslqr"></div>
<div class="warn" data-t="为什么不把第一次算出的序列直接用到底？">如果模型完美，离线算好的 $U^*_{[0]}$ 依次执行，和每步在线重算的结果几乎一样。但只要中途出现扰动，离线方案就会“按原计划”一错到底——它本质上是<b>开环</b>控制；在线滚动计算每次用实测状态重新规划，是<b>闭环</b>控制。（LQR 虽然离线算，但算出的是反馈矩阵，每次乘以当前状态，所以也是闭环。）</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>$X=\Phi\bm x+\Gamma U$；$H=\Gamma^T\Omega\Gamma+\Psi$，$F=\Gamma^T\Omega\Phi$；$J=U^TF\bm x+\frac12U^THU$。</li>
<li>无约束：$U^*=-H^{-1}F\bm x$，只用第一项；Np 足够大时与 LQR 相同。</li>
<li>在线滚动计算 = 闭环，能抵抗扰动。</li>
</ol></div>`,
    sims: {
      matrices(el) {
        const ui = card(el, { title: "矩阵搭建器：Φ、Γ、H、F", hint: "选择系统和 Np，查看各矩阵的数值与维度。（例 5.3.1：二输入二状态，Np=2，可与书中数值核对。）", layout: "stack" });
        const box = document.createElement("div"); box.style.cssText = "overflow-x:auto;font-size:.9rem"; ui.stage.appendChild(box);
        let ex = "e531";
        const sys = { e531: { A: [[1, 0.1], [0, -2]], B: [[0, 0.2], [-0.1, 0.5]], Q: M.eye(2), R: M.scale(M.eye(2), 0.1), S: M.eye(2), n: "例 5.3.1" }, one: { A: [[1]], B: [[1]], Q: [[1]], R: [[1]], S: [[1]], n: "一维例（A=B=Q=R=S=1）" } };
        ui.seg({ label: "系统", options: [["e531", "例 5.3.1"], ["one", "一维例"]], value: ex, onChange: v => { ex = v; draw(); } });
        const sN = ui.slider({ label: "Np", min: 1, max: 4, step: 1, value: 2, fmt: v => v.toFixed(0), onInput: draw });
        function draw() {
          const s = sys[ex], mm = mpcMatrices(s.A, s.B, s.Q, s.R, s.S, sN.get());
          const blk = (name, X) => `<div style="display:inline-block;margin:6px 14px 6px 0;vertical-align:top">$${name}=${mtex(X, 4)}$<div class="sim-note" style="text-align:center">${X.length}×${X[0].length}</div></div>`;
          box.innerHTML = blk("\\Phi", mm.Phi) + blk("\\Gamma", mm.Gam) + "<br>" + blk("H", mm.H) + blk("F", mm.F);
          CE.renderMath(box);
          ui.out(`${s.n}：n = ${s.A.length}，p = ${s.B[0].length}，Np = ${sN.get()}<br>U 的维度 = p·Np = ${s.B[0].length * sN.get()}；矩阵随 Np 线性（H 平方）膨胀`);
        }
        draw();
      },
      vslqr(el) {
        const ui = card(el, { title: "MPC → LQR，以及开环 vs 闭环", hint: "左：一维系统不同 Np 下的最优序列 U*（第一项与 LQR 的 −0.61803 对比，复现表 5.3.1）。右：Np=5，在 k=2 时刻给状态一个扰动，比较“在线滚动”与“离线一次算完”。" });
        const [a, b] = ui.row("", "");
        const pb = new Plot(a, { height: 260, x: [0.3, 10.7], y: [-0.7, 0.05], xlabel: "序列中的第 i 项", title: "U* 各元素" });
        const px = new Plot(b, { height: 260, x: [0, 8], y: [-0.1, 1.05], xlabel: "k", title: "状态 x" });
        const sN = ui.slider({ label: "Np（左图）", min: 1, max: 10, step: 1, value: 2, fmt: v => v.toFixed(0), onInput: draw });
        const sD = ui.slider({ label: "扰动大小（右图）", min: 0, max: 0.5, step: 0.01, value: 0.2, onInput: draw });
        function U(Np, x) { const mm = mpcMatrices([[1]], [[1]], [[1]], [[1]], [[1]], Np); return M.solve(mm.H, M.mulv(mm.F, [x])).map(v => -v); }
        function draw() {
          const Np = sN.get(), u = U(Np, 1);
          pb.setRange([0.3, Math.max(Np, 3) + 0.7], null);
          pb.render = () => { pb.clear(); pb.hline(-0.61803, COL.c3, [6, 4]); pb.text(Math.max(Np, 3) + 0.6, -0.65, "LQR: −0.61803", COL.c3, "right"); u.forEach((v, i) => { pb.rect(i + 1 - 0.3, v, i + 1 + 0.3, 0, i === 0 ? COL.c2 : COL.c2 + "55"); pb.text(i + 1, v - 0.035, f2(v, 4), COL.ink2, "center", 10); }); };
          // 在线 vs 离线
          const d = sD.get(), off = U(5, 1); let xo = 1, xn = 1; const XO = [1], XN = [1];
          for (let k = 0; k < 8; k++) { if (k === 2) { xo += d; xn += d; XO[2] = xo; XN[2] = xn; } const un = U(5, xn)[0], uo = k < 5 ? off[k] : 0; xn = xn + un; xo = xo + uo; XN.push(xn); XO.push(xo); }
          px.render = () => { px.clear(); px.line(XO.map((_, i) => i), XO, COL.c4, 2.2, [6, 4]); px.line(XN.map((_, i) => i), XN, COL.c1, 2.6); XN.forEach((v, i) => px.dot(i, v, COL.c1, 3.5)); px.vline(2, COL.muted, [3, 3]); px.text(2.1, 0.95, "扰动", COL.muted); px.legend([{ label: "在线（闭环）", color: COL.c1 }, { label: "离线（开环）", color: COL.c4, dash: [6, 4] }], "tr"); };
          pb.render(); px.render();
          ui.out(`Np = ${Np}：U* = −[${u.map(v => f2(-v, 5)).join(", ")}]<br>第一项 ${f2(u[0], 5)}　（LQR：−0.61803）<br>右图：离线方案终点 x = <b>${f2(XO[XO.length - 1], 3)}</b>，在线方案 <b>${f2(XN[XN.length - 1], 3)}</b>`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 5.4–5.5 追踪与约束 ════════════════ */
  Course.add({
    id: "v2-mpc3", vol: 2, ord: 110, no: "5.4–5.5", chName: "第5章 模型预测控制", ref: "卷2 §5.4–5.5",
    title: "含约束的 MPC：硬约束直接写进优化",
    sub: "MPC 的真正优势：把输入、状态的上下限变成 $MU\\le\\bar\\beta+b\\,x$ 这样的线性约束，交给 QP 求解器。",
    goals: ["会把轨迹追踪问题交给 MPC（复用第 4 章的增广方法）", "会把上下限约束转化为 QP 标准形式", "理解硬约束与软约束的区别与配合"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">5.4</span>轨迹追踪</h2>
<p>与第 4 章完全相同：用稳态非零参考（[F2]）或输入增量（[F3]）构造增广系统，再把增广矩阵交给 MPC 的 [F4]/[F5]。在<b>无约束</b>且 Np 足够大（如 20）时，MPC 与 LQR 的表现几乎一样——MPC 还没有展现出它的优势。优势在约束。</p>

<h2><span class="sec">5.5</span>约束转化为标准形式</h2>
<p>一般线性约束 $\mathcal M\bm x_{[k+i]}+\mathcal F\bm u_{[k+i]}\le\bm\beta$（i=0…Np−1）及末端 $\mathcal M_{N_p}\bm x_{[k+N_p]}\le\bm\beta_{N_p}$。把预测区间内的约束全部摞起来，再代入 $X=\Phi\bm x+\Gamma U$：</p>
<div class="eq" data-no="(5.5.7)–(5.5.10)" data-title="约束的 QP 形式">
$$\overline{\mathcal M}\bm x_{[k|k]}+\overline{\overline{\mathcal M}}X+\overline{\mathcal F}U\le\bar{\bm\beta}\quad\Longrightarrow\quad \underbrace{(\overline{\overline{\mathcal M}}\Gamma+\overline{\mathcal F})}_{M}\,U\;\le\;\bar{\bm\beta}+\underbrace{\left(-\overline{\mathcal M}-\overline{\overline{\mathcal M}}\Phi\right)}_{b}\bm x_{[k|k]}$$
<div class="say">右边含有当前状态——约束本身也随状态变化，同样带有反馈作用。最常见的上下限 $\bm u_{low}\le\bm u\le\bm u_{high}$、$\bm x_{low}\le\bm x\le\bm x_{high}$ 就是 $\begin{bmatrix}-I\\I\end{bmatrix}\bm u\le\begin{bmatrix}-\bm u_{low}\\\bm u_{high}\end{bmatrix}$ 这类形式；±∞ 表示不设限，编程时这些行要剔除。</div>
</div>
<p>完整的含约束 MPC 每一步就是：$\min_U\;U^TF\bm x+\frac12U^THU$，s.t. $MU\le\bar{\bm\beta}+b\bm x$。随着 Np 增大，约束矩阵会变得非常大，需要慎重选择预测区间。</p>
<p><b>例 5.5.1</b>：例 5.3.1 的二维系统，Np=2，约束 $x_2\le0$、$u_2\ge-3$，初始 $\bm x=[1,-1]^T$。无约束时第一步 $\bm u=[0.53,-3.56]$，下一步 $x_2=0.17>0$，两个约束都被违反；加约束后第一步变为 $[5,-3]$——$u_2$ 恰好卡在下限，而 $u_1$ 为了配合<b>大幅增加</b>：输入之间、状态之间是耦合联动的。约束太严时甚至可能找不到可行解。</p>
<div class="sim" data-id="ex551"></div>

<h3>硬约束与软约束</h3>
<p>弹簧–质量–阻尼追踪（R=1）时有约 10% 的超调。原因是 R 大，算法让 $u$ 尽量贴近 $u_d$——R 相当于对 δu 的<b>软约束</b>。如果绝对不允许超调，就加一条<b>硬约束</b> $x_1\le1$。两者冲突时，系统在满足硬约束的前提下尽量满足软约束：此时 u 会较长时间与 $u_d$ 保持距离，以换取“零超调”。</p>
<div class="sim" data-id="msdcons"></div>
<div class="tip">复杂问题拆成简单问题，再用写好的模块组合（[F2] 增广 + [F4] 代价矩阵 + [F6] 约束矩阵 + [F7] 含约束 QP），是非常有效的工程实践。</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>追踪问题：增广后交给 MPC；无约束时与 LQR 等价。</li>
<li>约束化为 $MU\le\bar{\bm\beta}+b\bm x$，与代价一起构成含约束 QP。</li>
<li>硬约束必须满足、软约束（权重）可以权衡；输入与状态之间互相耦合。</li>
</ol></div>`,
    sims: {
      ex551(el) {
        const ui = card(el, { title: "约束开关实验室（例 5.5.1）", hint: "二维系统，Np=2，$\\bm x_0=[1,-1]$。实线：含约束 MPC；虚线：无约束。红色虚线是约束边界。" });
        const [a, b] = ui.row("", "");
        const px = new Plot(a, { height: 250, x: [0, 10], y: [-1.2, 1.2], xlabel: "k", title: "状态 x₁（蓝）、x₂（橙）" });
        const pu = new Plot(b, { height: 250, x: [0, 10], y: [-4, 6], xlabel: "k", title: "输入 u₁（蓝）、u₂（橙）" });
        const sX = ui.slider({ label: "x₂ 上限", min: -0.5, max: 1, step: 0.05, value: 0, onInput: draw });
        const sU = ui.slider({ label: "u₂ 下限", min: -5, max: -1, step: 0.1, value: -3, onInput: draw });
        const A = [[1, 0.1], [0, -2]], B = [[0, 0.2], [-0.1, 0.5]];
        function sim(constrained) {
          const mm = mpcMatrices(A, B, M.eye(2), M.scale(M.eye(2), 0.1), M.eye(2), 2);
          const cs = constrained ? mpcConstraints(mm, [-INF, sU.get()], [INF, INF], [-INF, -INF], [INF, sX.get()]) : { Mc: [], rhs: () => [] };
          const sol = qpPrep(mm.H, cs.Mc);
          let x = [1, -1]; const o = { x1: [1], x2: [-1], u1: [], u2: [], ok: true };
          for (let k = 0; k < 10; k++) { const r = sol(M.mulv(mm.F, x), cs.rhs(x), { iters: 3000 }); const u = r.x.slice(0, 2); o.u1.push(u[0]); o.u2.push(u[1]); x = M.vadd(M.mulv(A, x), M.mulv(B, u)); o.x1.push(x[0]); o.x2.push(x[1]); }
          return o;
        }
        function draw() {
          const c = sim(true), f = sim(false), ks = linspace(0, 10, 11), ku = linspace(0, 9, 10);
          const stair = (ys) => { const xs = [], yy = []; ys.forEach((v, i) => { xs.push(i, i + 1); yy.push(v, v); }); return [xs, yy]; };
          px.render = () => { px.clear(); px.hline(sX.get(), COL.c4, [5, 4], 1.6); px.line(ks, f.x1, COL.c1, 1.4, [5, 4]); px.line(ks, f.x2, COL.c2, 1.4, [5, 4]); px.line(ks, c.x1, COL.c1, 2.4); px.line(ks, c.x2, COL.c2, 2.4); };
          pu.autoY([c.u1, c.u2, f.u1, f.u2, [sU.get()]], 0.1);
          pu.render = () => { pu.clear(); pu.hline(sU.get(), COL.c4, [5, 4], 1.6); let s = stair(f.u1); pu.line(s[0], s[1], COL.c1, 1.4, [5, 4]); s = stair(f.u2); pu.line(s[0], s[1], COL.c2, 1.4, [5, 4]); s = stair(c.u1); pu.line(s[0], s[1], COL.c1, 2.4); s = stair(c.u2); pu.line(s[0], s[1], COL.c2, 2.4); };
          px.render(); pu.render();
          ui.out(`无约束首步 u = [${f2(f.u1[0])}, ${f2(f.u2[0])}]，x₂(1) = ${f2(f.x2[1])}<br>含约束首步 u = [<b>${f2(c.u1[0])}</b>, <b>${f2(c.u2[0])}</b>]，x₂(1) = <b>${f2(c.x2[1])}</b><br><span class="sim-note">把约束收紧，u₁ 会被迫做出更大的“配合”。</span>`);
        }
        draw();
      },
      msdcons(el) {
        const ui = card(el, { title: "零超调：硬约束 x₁ ≤ 1 vs 调小 R", hint: "弹簧–质量–阻尼追踪 $x_d=[1,0]$，MPC（稳态非零参考），Np=20。" });
        const [a, b] = ui.row("", "");
        const px = new Plot(a, { height: 250, x: [0, 100], y: [0, 1.2], xlabel: "k", title: "位移 x₁" });
        const pu = new Plot(b, { height: 250, x: [0, 100], y: [0, 3], xlabel: "k", title: "输入 u（u_d = 1）" });
        const sR = ui.slider({ label: "R（软约束）", min: 0.01, max: 10, value: 1, log: true, onInput: draw });
        const cH = ui.check({ label: "硬约束 0 ≤ x₁ ≤ 1", value: true, onChange: draw });
        const D = CE.M.c2d([[0, 1], [-1, -0.5]], [[0], [1]], 0.1);
        function sim(hard) {
          const g = augSS(D.A, D.B, [1, 0], M.eye(2), M.eye(2)), mm = mpcMatrices(g.Aa, g.Ba, g.Qa, [[sR.get()]], g.Sa, 20);
          const cs = hard ? mpcConstraints(mm, [-INF], [INF], [0, -INF, -INF, -INF], [1, INF, INF, INF]) : { Mc: [], rhs: () => [] };
          const sol = qpPrep(mm.H, cs.Mc); let x = [0, 0]; const o = { x: [], u: [] };
          for (let k = 0; k <= 100; k++) { const xa = [...x, 1, 0], du = sol(M.mulv(mm.F, xa), cs.rhs(xa)).x[0], u = du + g.ud; o.x.push(x[0]); o.u.push(u); x = M.vadd(M.mulv(D.A, x), M.flat(D.B).map(v => v * u)); }
          return o;
        }
        function draw() {
          const f = sim(false), c = cH.get() ? sim(true) : null, ks = linspace(0, 100, 101);
          px.render = () => { px.clear(); px.hline(1, COL.c4, [5, 4]); px.line(ks, f.x, COL.muted, 1.6, [5, 4]); if (c) px.line(ks, c.x, COL.c1, 2.6); };
          pu.autoY([f.u, c ? c.u : [], [0, 1]], 0.1);
          pu.render = () => { pu.clear(); pu.hline(1, COL.c3, [5, 4]); pu.line(ks, f.u, COL.muted, 1.6, [5, 4]); if (c) pu.line(ks, c.u, COL.c2, 2.4); };
          px.render(); pu.render();
          ui.out(`无约束（虚线）最大位移 <b>${f2(Math.max(...f.x), 3)}</b>${c ? `<br>硬约束（实线）最大位移 <b>${f2(Math.max(...c.x), 3)}</b>` : ""}<br><span class="sim-note">另一条路：不加硬约束，把 R 调小（如 0.1），超调同样消失——但控制量更激进。</span>`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 5.6–5.7 无人机案例 ════════════════ */
  Course.add({
    id: "v2-mpc4", vol: 2, ord: 120, no: "5.6–5.7", chName: "第5章 模型预测控制", ref: "卷2 §5.6–5.8",
    title: "MPC 实战：无人机爬升，以及预测区间与采样时间",
    sub: "第 4 章 LQR 没解决的问题：同时约束推力和速度。再看看预测区间和采样时间怎样在“性能”与“算得过来”之间权衡。",
    goals: ["用含约束 MPC 实现带状态约束的无人机爬升", "理解预测区间 Np 太短的后果", "理解采样时间、预测时长与计算时间的权衡"],
    minutes: 25,
    html: String.raw`
<p>无人机模型同第 4 章（增广状态，$u_d=g=10$），约束与动态规划那一节一致：净加速度 $\delta u\in[-3,2]$（即推力 $u\in[7,12]$）、高度 $x_1\in[0,10]$、<b>速度 $x_2\in[0,3]$</b>。权重 $Q=S=\mathrm{diag}(1,1,0)$。</p>
<div class="sim" data-id="uav"></div>
<p><b>结果分析</b>（Np=20，R=0.1）：推力先顶在上限 12，速度一到 3 m/s 就转为匀速（推力回到 10），接近目标时减速——与第 4 章动态规划得到的“最大加速—匀速—最大减速”结构非常接近。与 LQR + 饱和相比，推力处于上限的时间<b>更短</b>，因为速度约束不允许一直加速。与动态规划的差别在于：MPC 把能耗也放进了代价，并且到达后还要保持稳定（滚动优化），所以会缓慢地把加速度降到 0。</p>
<h2><span class="sec">5.6.2</span>预测区间的影响</h2>
<ul><li><b>性能</b>：Np 越长，看得越远，越能提前为约束做准备、抵抗扰动。</li>
<li><b>计算量</b>：Np 越长，矩阵越大，求解越慢，可能来不及在一个采样周期内算完。</li>
<li><b>可行性</b>：Np 太短，可能根本找不到满足约束的解。</li></ul>
<p>Np=20、$T_s=0.1$ 相当于看 2 s 的未来；Np=5 只看 0.5 s——无人机爬升缓慢，100 步后还到不了 10 m。</p>
<h2><span class="sec">5.6.3</span>采样时间</h2>
<p>系统约 8 s 稳定。若想一次“看完”整个过程（10 s），$T_s=0.1$ 时需 Np=100——作者的笔记本上首步计算要 3 秒多，远超 0.1 s 的采样周期：控制量还没算出来，下一时刻已经到了。改用 $T_s=0.5$、Np=20，同样覆盖 10 s，计算时间降到 0.3 s 以内，可以实时运行——代价是控制更新间隔变长、精度下降。</p>
<div class="idea" data-t="人无远虑，必有近忧">只看眼前会错过长远趋势；但“深虑”过远又会拖慢当下的行动。在远见与眼下之间取得平衡——MPC 的参数选择正是如此。</div>
<h2><span class="sec">5.7</span>MPC 的发展方向</h2>
<ul><li><b>非线性 MPC（NMPC）</b>：用非线性模型和非线性规划，适用范围更广，但计算更难。</li>
<li><b>大规模 / 分布式 MPC</b>：并行计算、高效优化算法。</li>
<li><b>数据驱动 MPC</b>：模型难建时，用实时数据建模与控制。</li>
<li><b>多目标 MPC</b>：同时兼顾跟踪、能耗、安全、鲁棒性。</li></ul>
<div class="summary"><h3>本课要点</h3><ol>
<li>含约束 MPC 能同时处理输入约束和状态约束，结果接近动态规划的最优解。</li>
<li>Np 太短：性能差、甚至不可行；太长：算不过来。</li>
<li>预测时长 = Np × Ts；增大 Ts 可以降低计算量，但牺牲控制精度。</li>
</ol></div>`,
    sims: {
      uav(el) {
        const ui = card(el, { title: "无人机爬升：含约束 MPC", hint: "实线：含约束 MPC；灰色虚线：第 4 章的 LQR + 饱和。右下柱状图：每一步的求解时间（在你的浏览器里实测），红线是采样周期 Ts。" });
        const sc = new Scene(ui.plotHost(), { height: 110 });
        const [a, b] = ui.row("", "");
        const ph = new Plot(a, { height: 200, x: [0, 10], y: [0, 11], xlabel: "t / s", title: "高度 x₁" });
        const pv = new Plot(b, { height: 200, x: [0, 10], y: [-0.5, 6], xlabel: "t / s", title: "速度 x₂（上限 3）" });
        const [c3, d4] = ui.row("", "");
        const pu = new Plot(c3, { height: 200, x: [0, 10], y: [6, 14], xlabel: "t / s", title: "推力 u ∈ [7, 12]" });
        const pt = new Plot(d4, { height: 200, x: [0, 10], y: [0, 0.05], xlabel: "t / s", title: "每步求解耗时 (s)" });
        const sN = ui.slider({ label: "预测区间 Np", min: 3, max: 60, step: 1, value: 20, fmt: v => v.toFixed(0), onInput: start });
        let Ts = 0.1;
        ui.seg({ label: "采样时间 Ts", options: [[0.1, "0.1 s"], [0.2, "0.2 s"], [0.5, "0.5 s"]], value: 0.1, onChange: v => { Ts = v; start(); } });
        const sR = ui.slider({ label: "R", min: 0.01, max: 10, value: 0.1, log: true, onInput: start });
        const cV = ui.check({ label: "速度约束 0 ≤ x₂ ≤ 3", value: true, onChange: start });
        let st = null;
        function lqrRun() {
          const { A, B } = droneAB(Ts), g = augSS(A, B, [10, 0, -10], [[1, 0, 0], [0, 1, 0], [0, 0, 0]], [[1, 0, 0], [0, 1, 0], [0, 0, 0]]);
          const { K } = M.dlqr(g.Aa, g.Ba, g.Qa, [[sR.get()]], g.Sa); let x = [0, 0, -10]; const o = { t: [], h: [], v: [], u: [] };
          for (let k = 0; k * Ts <= 10; k++) { const u = clamp(10 - M.dot(K[0], [...x, 10, 0, -10]), 7, 12); o.t.push(k * Ts); o.h.push(x[0]); o.v.push(x[1]); o.u.push(u); x = M.vadd(M.mulv(A, x), M.flat(B).map(v => v * u)); }
          return o;
        }
        function start() {
          const { A, B } = droneAB(Ts), W = [[1, 0, 0], [0, 1, 0], [0, 0, 0]], g = augSS(A, B, [10, 0, -10], W, W);
          const mm = mpcMatrices(g.Aa, g.Ba, g.Qa, [[sR.get()]], g.Sa, sN.get());
          const cs = mpcConstraints(mm, [-3], [2], [0, cV.get() ? 0 : -INF, -INF, -INF, -INF, -INF], [10, cV.get() ? 3 : INF, INF, INF, INF, INF]);
          st = { A, B, g, mm, cs, sol: qpPrep(mm.H, cs.Mc), x: [0, 0, -10], k: 0, o: { t: [], h: [], v: [], u: [], ct: [] }, lqr: lqrRun(), rows: cs.Mc.length };
        }
        function stepSim() {
          const s = st; if (s.k * Ts > 10 + 1e-9) return false;
          const xa = [...s.x, 10, 0, -10], t0 = performance.now(), r = s.sol(M.mulv(s.mm.F, xa), s.cs.rhs(xa), { iters: 600 }), ct = (performance.now() - t0) / 1000;
          const u = clamp(r.x[0] + s.g.ud, 7, 12);
          s.o.t.push(s.k * Ts); s.o.h.push(s.x[0]); s.o.v.push(s.x[1]); s.o.u.push(u); s.o.ct.push(ct);
          s.x = M.vadd(M.mulv(s.A, s.x), M.flat(s.B).map(v => v * u)); s.k++; return true;
        }
        start();
        CE.loop(ui.root, () => {
          const t0 = performance.now(); while (performance.now() - t0 < 12 && stepSim());
          const s = st, o = s.o, L = s.lqr;
          ph.render = () => { ph.clear(); ph.hline(10, COL.c3, [5, 4]); ph.line(L.t, L.h, COL.muted, 1.4, [5, 4]); ph.line(o.t, o.h, COL.c1, 2.6); };
          pv.render = () => { pv.clear(); pv.hline(3, COL.c4, [5, 4]); pv.line(L.t, L.v, COL.muted, 1.4, [5, 4]); pv.line(o.t, o.v, COL.c1, 2.6); };
          pu.render = () => { pu.clear(); pu.rect(0, 7, 10, 12, COL.c3 + "18"); pu.line(L.t, L.u, COL.muted, 1.4, [5, 4]); pu.line(o.t, o.u, COL.c2, 2.4); };
          const mx = Math.max(Ts * 1.2, ...o.ct) * 1.1; pt.setRange(null, [0, mx]);
          pt.render = () => { pt.clear(); o.ct.forEach((c, i) => pt.rect(o.t[i], 0, o.t[i] + Ts * 0.8, c, c > Ts ? COL.c4 : COL.c5 + "aa")); pt.hline(Ts, COL.c4, [5, 4], 1.6); };
          ph.render(); pv.render(); pu.render(); pt.render();
          const c = sc.ctx, W = sc.W, H = sc.H; sc.clear(); const h = o.h[o.h.length - 1] || 0, x = W * 0.5, y2p = hh => H - 12 - hh / 11 * (H - 40), y = y2p(h);
          c.strokeStyle = COL.c3; c.setLineDash([5, 4]); c.beginPath(); c.moveTo(0, y2p(10)); c.lineTo(W, y2p(10)); c.stroke(); c.setLineDash([]);
          c.strokeStyle = COL.c1; c.lineWidth = 3; c.beginPath(); c.moveTo(x - 30, y); c.lineTo(x + 30, y); c.stroke(); [-30, 30].forEach(d => { c.beginPath(); c.ellipse(x + d, y - 4, 14, 3, 0, 0, 7); c.stroke(); }); c.fillStyle = COL.c1; c.fillRect(x - 8, y - 5, 16, 11);
          c.fillStyle = COL.ink2; c.font = "12px Consolas"; c.fillText(`h=${f2(h, 2)} m  v=${f2(o.v[o.v.length - 1] || 0, 2)} m/s  u=${f2(o.u[o.u.length - 1] || 10, 2)}`, 8, 16);
          const late = o.ct.filter(v => v > Ts).length;
          ui.out(`预测时长 Np·Ts = <b>${f2(sN.get() * Ts, 1)} s</b><br>约束行数：${s.rows}，决策变量：${sN.get()}<br>最高速度：<b>${f2(Math.max(...o.v, 0), 2)}</b> m/s<br>10 s 后高度：<b>${f2(o.h[o.h.length - 1], 2)}</b> m<br>最长单步耗时：${f2(Math.max(...o.ct, 0) * 1000, 1)} ms ${late ? `<span class="bad">（${late} 步超过 Ts！）</span>` : '<span class="good">（都在 Ts 内）</span>'}`);
        });
        ui.root.insertAdjacentHTML("beforeend", `<div class="tasks"><b>试一试</b><ol><li>Np=20：推力 12 → 速度到 3 后匀速 → 减速，速度从不超过 3（LQR+饱和会超过）。</li><li>Np=5：只看 0.5 s 的未来，爬升明显变慢。</li><li>Ts=0.5、Np=20：同样看 10 s，但计算量小得多；控制量变成粗台阶。</li></ol></div>`);
      },
    },
  });
})();
