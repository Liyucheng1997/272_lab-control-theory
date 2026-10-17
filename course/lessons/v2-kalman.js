/* 卷2 第6章：卡尔曼滤波器 / 数据融合 / 协方差 / KF+MPC / EKF */
(function () {
  const { M, P, Plot, Scene, card, linspace, f2, COL, clamp, rk4, cfmt, mtex, mpcMatrices, mpcConstraints, qpPrep } = CE;
  const gauss = (x, mu, s) => Math.exp(-((x - mu) ** 2) / (2 * s * s)) / Math.sqrt(2 * Math.PI * s * s);
  const PLAYERS = [[183, 75, 29, "孙兴慜"], [175, 71, 29, "萨拉赫"], [187, 83, 36, "罗纳尔多"], [185, 74, 35, "本泽马"], [176, 73, 22, "儒尼奥尔"], [176, 67, 35, "阿斯帕斯"], [185, 79, 34, "莱万"], [191, 73, 26, "希克"], [195, 88, 22, "哈兰德"], [185, 80, 32, "因莫比莱"], [174, 81, 25, "劳塔罗"], [180, 78, 27, "西蒙尼"], [178, 73, 24, "姆巴佩"], [170, 68, 32, "本耶德尔"], [184, 71, 25, "特雷尔"]];
  // 2 维线性 KF 一步
  function kfStep(A, Bu, Hm, Qc, Rc, xh, Pk, z) {
    const xm = M.vadd(M.mulv(A, xh), Bu), Pm = M.add(M.mul(M.mul(A, Pk), M.T(A)), Qc);
    const S = M.add(M.mul(M.mul(Hm, Pm), M.T(Hm)), Rc), K = M.mul(M.mul(Pm, M.T(Hm)), M.inv(S));
    const xp = M.vadd(xm, M.mulv(K, M.vsub(z, M.mulv(Hm, xm)))), Pp = M.mul(M.sub(M.eye(A.length), M.mul(K, Hm)), Pm);
    return { xm, xp, Pp: M.sym(Pp), K, Pm };
  }

  /* ════════════════ 6.1–6.2 数据融合与协方差 ════════════════ */
  Course.add({
    id: "v2-kf1", vol: 2, ord: 130, no: "6.1–6.2", chName: "第6章 卡尔曼滤波器", ref: "卷2 §6.1–6.2",
    title: "数据融合：两个不太准的传感器，怎样合出一个更准的",
    sub: "卡尔曼滤波器的灵魂是“按可信度加权融合”。先从取平均、递归算法、正态分布和协方差讲起。",
    goals: ["把求平均改写成递归形式 $\\hat z_k=\\hat z_{k-1}+K_k(z_k-\\hat z_{k-1})$", "推导两个传感器的最优融合系数 $K^*=\\frac{\\sigma_1^2}{\\sigma_1^2+\\sigma_2^2}$", "理解协方差与协方差矩阵"],
    minutes: 30,
    html: String.raw`
<p><span class="kw">卡尔曼滤波器</span>是一种最优的、递归的数字处理算法，它同时是<b>观测器</b>（估计看不见的状态）和<b>滤波器</b>（给含噪信号去噪）。之所以需要它，是因为现实中有三类不确定性：模型不完美、扰动无法建模、传感器有误差。它的核心手段就是<b>数据融合</b>。</p>
<h2><span class="sec">6.1</span>递归算法</h2>
<div class="story" data-t="无人机的气压计">无人机用气压计测离地高度，每次读数都不一样：50、50.1、……、50.6。最朴素的办法是取平均——误差有正有负，平均后相互抵消。但每来一个新数据就要把所有历史数据重新加一遍，越来越慢，还要全部存下来。</div>
<div class="eq" data-no="(6.1.2)–(6.1.3)" data-title="把平均写成递归">
$$\hat z_k=\frac1k\sum_{i=1}^kz_i=\frac{k-1}{k}\hat z_{k-1}+\frac1kz_k=\hat z_{k-1}+\underbrace{\frac1k}_{K_k}\left(z_k-\hat z_{k-1}\right)$$
<div class="say">新估计 = 旧估计 + 增益 ×（新测量 − 旧估计）。每步只需两次加法一次乘法，只存上一次的估计。k 越大，K 越小：信心已足，新数据只做微调。这就是卡尔曼滤波更新公式的雏形。</div>
</div>
<div class="sim" data-id="recavg"></div>

<h2><span class="sec">6.2</span>概率论初步与最优融合</h2>
<p>期望 $E(X)$ 是“平均值”，方差 $\mathrm{Var}(X)=E\big((X-E(X))^2\big)=E(X^2)-E(X)^2$ 描述离散程度，标准差 $\sigma=\sqrt{\mathrm{Var}}$。注意方差不是线性的：$\mathrm{Var}(aX)=a^2\mathrm{Var}(X)$；独立时 $\mathrm{Var}(X+Y)=\mathrm{Var}(X)+\mathrm{Var}(Y)$。大多数随机误差服从<b>正态分布</b> $N(\mu,\sigma^2)$：$\mu\pm\sigma$ 内的概率为 68.27%，$\mu\pm3\sigma$ 内为 99.73%。</p>
<div class="story" data-t="两个传感器">无人机再加一个超声波传感器。气压计读 $z_1=50$ m，误差标准差 $\sigma_1=1$ m（真值 99.73% 在 47～53 m）；超声波读 $z_2=52$ m，$\sigma_2=0.5$ m（真值 99.73% 在 50.5～53.5 m）。怎么融合最好？</div>
<p>设融合结果 $\hat z=z_1+K(z_2-z_1)$。估计误差 $\hat e=(1-K)e_1+Ke_2$，均值为 0（无偏），方差为 $(1-K)^2\sigma_1^2+K^2\sigma_2^2$（两传感器独立）。对 K 求导令其为零：</p>
<div class="eq" data-no="(6.2.13)–(6.2.14)" data-title="最优融合">
$$K^*=\frac{\sigma_1^2}{\sigma_1^2+\sigma_2^2}=\frac{1}{1+0.25}=0.8,\qquad \hat z=50+0.8\times(52-50)=51.6,\qquad \sigma_{\hat e}=0.45$$
<div class="say">谁更准就更相信谁；融合后的标准差 0.45 比两个传感器都小！等价地：$\frac{1}{\sigma_{\hat e}^2}=\frac{1}{\sigma_1^2}+\frac{1}{\sigma_2^2}$——“精度”可以相加。</div>
</div>
<div class="sim" data-id="fusion"></div>

<h3>协方差与协方差矩阵</h3>
<p>两个量相互关联时，还要考虑它们的<b>协方差</b>：$\mathrm{Cov}(X_1,X_2)=E\big((X_1-E X_1)(X_2-E X_2)\big)$。正：同涨同跌；负：此消彼长；0：没有线性关系（独立时一定为 0）。方差就是变量与自身的协方差。把所有两两关系排成矩阵：</p>
<div class="eq" data-no="(6.2.25)" data-title="协方差矩阵">
$$C(\bm X)=E\left[(\bm X-\bm\mu)(\bm X-\bm\mu)^T\right]=\begin{bmatrix}\sigma_{1}^2&\sigma_{12}&\cdots\\\sigma_{21}&\sigma_2^2&\cdots\\\vdots&&\ddots\end{bmatrix},\qquad C=C^T$$
<div class="say">对角线是各自的方差，非对角线是两两协方差。书中用 15 名射手的身高、体重、年龄算出：$C=\begin{bmatrix}46.69&24.04&-1.7\\24.04&33.69&-4.13\\-1.7&-4.13&23.70\end{bmatrix}$——身高与体重强正相关，年龄与二者几乎无关。</div>
</div>
<div class="sim" data-id="cov"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>递归估计：$\hat z_k=\hat z_{k-1}+K_k(z_k-\hat z_{k-1})$，省时省内存。</li>
<li>两个独立测量的最优融合：$K^*=\sigma_1^2/(\sigma_1^2+\sigma_2^2)$，融合后比任何一个都准。</li>
<li>协方差描述变量之间的联动；协方差矩阵对称，对角线是方差。</li>
</ol></div>`,
    sims: {
      recavg(el) {
        const ui = card(el, { title: "递归平均：增益 K = 1/k 越来越小", hint: "真实高度 50.3 m，气压计每次测量都带随机误差。蓝线是递归估计，下图是增益 K。" });
        const pz = new Plot(ui.plotHost(), { height: 220, x: [0, 60], y: [47, 54], xlabel: "测量次数 k", ylabel: "高度 m" });
        const pk = new Plot(ui.plotHost(), { height: 140, x: [0, 60], y: [0, 1.05], xlabel: "k", ylabel: "K" });
        const sS = ui.slider({ label: "测量噪声 σ", min: 0.1, max: 2, step: 0.05, value: 0.8, onInput: draw });
        let seed = 1; ui.buttons([["🎲 重新测量", () => { seed++; draw(); }]]);
        const cF = ui.check({ label: "对比固定增益 K = 0.2（会“遗忘”旧数据）", value: false, onChange: draw });
        function draw() {
          const rg = CE.rng(seed), s = sS.get(), z = Array.from({ length: 60 }, () => 50.3 + s * rg.n());
          let h = 0, hf = z[0]; const H = [], HF = [], K = [];
          z.forEach((zz, i) => { const k = i + 1; h = h + (zz - h) / k; H.push(h); K.push(1 / k); hf = i ? hf + 0.2 * (zz - hf) : zz; HF.push(hf); });
          const ks = z.map((_, i) => i + 1);
          pz.render = () => { pz.clear(); pz.hline(50.3, COL.c3, [6, 4]); z.forEach((v, i) => pz.dot(i + 1, v, COL.c2 + "aa", 3)); if (cF.get()) pz.line(ks, HF, COL.c5, 1.8); pz.line(ks, H, COL.c1, 2.6); };
          pk.render = () => { pk.clear(); pk.line(ks, K, COL.c4, 2.2); if (cF.get()) pk.hline(0.2, COL.c5, [4, 3]); };
          pz.render(); pk.render();
          ui.out(`60 次后的估计：<b class="big">${f2(H[59], 3)}</b> m（真值 50.3）<br>${cF.get() ? `固定 K=0.2：${f2(HF[59], 3)} m——一直在抖动` : "递归公式与一次性求平均完全等价"}`);
        }
        draw();
      },
      fusion(el) {
        const ui = card(el, { title: "两个传感器的最优融合", hint: "橙：气压计；紫：超声波；蓝：最优融合结果。右图：融合误差方差随 K 的变化，最低点就是 K*。" });
        const [a, b] = ui.row("", "");
        const pg = new Plot(a, { height: 270, x: [45, 56], y: [0, 1], xlabel: "高度 m", title: "概率密度" });
        const pv = new Plot(b, { height: 270, x: [0, 1], y: [0, 1.2], xlabel: "K", title: "融合误差方差 σ²(K)" });
        const s1 = ui.slider({ label: "气压计读数 z₁", min: 46, max: 55, step: 0.1, value: 50, onInput: draw });
        const d1 = ui.slider({ label: "气压计误差 σ₁", min: 0.2, max: 2.5, step: 0.05, value: 1, onInput: draw });
        const s2 = ui.slider({ label: "超声波读数 z₂", min: 46, max: 55, step: 0.1, value: 52, onInput: draw });
        const d2 = ui.slider({ label: "超声波误差 σ₂", min: 0.2, max: 2.5, step: 0.05, value: 0.5, onInput: draw });
        let Kman = null;
        pv.drag({ down: x => { Kman = clamp(x, 0, 1); draw(); }, move: x => { Kman = clamp(x, 0, 1); draw(); }, up: () => { }, hover: () => true });
        function draw() {
          const z1 = s1.get(), z2 = s2.get(), a1 = d1.get(), a2 = d2.get(), Ks = a1 * a1 / (a1 * a1 + a2 * a2), K = Kman ?? Ks;
          const zh = z1 + K * (z2 - z1), sh = Math.sqrt((1 - K) ** 2 * a1 * a1 + K * K * a2 * a2), zs = z1 + Ks * (z2 - z1), ss = Math.sqrt((1 - Ks) ** 2 * a1 * a1 + Ks * Ks * a2 * a2);
          const xs = linspace(45, 56, 400), ymax = Math.max(gauss(0, 0, Math.min(a1, a2, sh)), 0.1) * 1.15;
          pg.setRange(null, [0, ymax]);
          pg.render = () => { pg.clear(); pg.line(xs, xs.map(x => gauss(x, z1, a1)), COL.c2, 2); pg.line(xs, xs.map(x => gauss(x, z2, a2)), COL.c5, 2); pg.line(xs, xs.map(x => gauss(x, zh, sh)), COL.c1, 3); pg.vline(zh, COL.c1, [4, 3]); pg.legend([{ label: "z₁", color: COL.c2 }, { label: "z₂", color: COL.c5 }, { label: "融合", color: COL.c1 }]); };
          const Kg = linspace(0, 1, 200), V = Kg.map(k => (1 - k) ** 2 * a1 * a1 + k * k * a2 * a2);
          pv.setRange(null, [0, Math.max(...V) * 1.1]);
          pv.render = () => { pv.clear(); pv.line(Kg, V, COL.c4, 2.4); pv.dot(Ks, ss * ss, COL.c3, 6); pv.text(Ks, ss * ss + pv.opt.y[1] * 0.05, "K*", COL.c3, "center", 12, true); if (Kman !== null) pv.dot(Kman, sh * sh, COL.c1, 6, COL.plotbg); };
          pg.render(); pv.render();
          ui.out(`K* = σ₁²/(σ₁²+σ₂²) = <b>${f2(Ks, 3)}</b><br>最优融合 ẑ = <b class="big">${f2(zs, 2)}</b> m，σ = <b>${f2(ss, 3)}</b><br>${Kman !== null ? `你选的 K = ${f2(Kman, 2)}：σ = ${f2(sh, 3)}（${sh > ss + 1e-6 ? "比最优差" : "就是最优"}）` : '<span class="sim-note">在右图中拖动，试试别的 K。</span>'}`);
        }
        draw();
      },
      cov(el) {
        const ui = card(el, { title: "协方差：点云的“形状”", hint: "左：15 名射手的数据（选择两个变量）及其 2σ 协方差椭圆。右：自己调协方差，生成随机点云。" });
        const [a, b] = ui.row("", "");
        const pd = new Plot(a, { height: 300, x: [165, 200], y: [60, 95], title: "射手数据" });
        const pr = new Plot(b, { height: 300, x: [-4, 4], y: [-4, 4], equal: true, title: "随机向量 N(0, C)" });
        let pair = [0, 1];
        ui.seg({ label: "变量对", options: [["01", "身高–体重"], ["02", "身高–年龄"], ["12", "体重–年龄"]], value: "01", onChange: v => { pair = [+v[0], +v[1]]; draw(); } });
        const sv1 = ui.slider({ label: "右图：σ₁²", min: 0.2, max: 3, step: 0.05, value: 1.5, onInput: draw });
        const sv2 = ui.slider({ label: "右图：σ₂²", min: 0.2, max: 3, step: 0.05, value: 1, onInput: draw });
        const sc = ui.slider({ label: "右图：协方差 σ₁₂", min: -1.6, max: 1.6, step: 0.05, value: 0.8, onInput: draw });
        const ellipse = (mu, C, k = 2) => { const e = M.eig2(C), th = linspace(0, 2 * Math.PI, 100); const [l1, l2] = e.vals.map(v => Math.sqrt(Math.max(0, v[0]))), v = e.vecs[0]; return th.map(t => { const u = k * l1 * Math.cos(t), w = k * l2 * Math.sin(t); return [mu[0] + u * v[0] - w * v[1], mu[1] + u * v[1] + w * v[0]]; }); };
        const names = ["身高", "体重", "年龄"], ranges = [[165, 200], [60, 95], [18, 40]];
        function draw() {
          const [i, j] = pair, X = PLAYERS.map(p => p[i]), Y = PLAYERS.map(p => p[j]), n = X.length;
          const mx = X.reduce((s, v) => s + v, 0) / n, my = Y.reduce((s, v) => s + v, 0) / n;
          const cxx = X.reduce((s, v) => s + (v - mx) ** 2, 0) / (n - 1), cyy = Y.reduce((s, v) => s + (v - my) ** 2, 0) / (n - 1), cxy = X.reduce((s, v, k) => s + (v - mx) * (Y[k] - my), 0) / (n - 1);
          pd.setRange(ranges[i], ranges[j]); pd.opt.xlabel = names[i]; pd.opt.ylabel = names[j];
          const el1 = ellipse([mx, my], [[cxx, cxy], [cxy, cyy]]);
          pd.render = () => { pd.clear(); pd.line(el1.map(p => p[0]), el1.map(p => p[1]), COL.c5, 2); X.forEach((x, k) => pd.dot(x, Y[k], COL.c2, 4.5, COL.plotbg)); pd.dot(mx, my, COL.c4, 4); };
          const C = [[sv1.get(), sc.get()], [sc.get(), sv2.get()]], det = C[0][0] * C[1][1] - C[0][1] ** 2;
          const rg = CE.rng(9), pts = [];
          if (det > 0) { const L11 = Math.sqrt(C[0][0]), L21 = C[0][1] / L11, L22 = Math.sqrt(C[1][1] - L21 * L21); for (let k = 0; k < 400; k++) { const a = rg.n(), b = rg.n(); pts.push([L11 * a, L21 * a + L22 * b]); } }
          pr.render = () => { pr.clear(); pts.forEach(p => pr.dot(p[0], p[1], COL.c1 + "77", 2.2)); if (det > 0) { const e = ellipse([0, 0], C); pr.line(e.map(p => p[0]), e.map(p => p[1]), COL.c4, 2); } };
          pd.render(); pr.render();
          ui.out(`数据：σ²(${names[i]}) = ${f2(cxx)}，σ²(${names[j]}) = ${f2(cyy)}<br>协方差 = <b>${f2(cxy)}</b>，相关系数 ${f2(cxy / Math.sqrt(cxx * cyy))}<br>右图 $C=${mtex(C, 2)}$ ${det <= 0 ? '<br><span class="bad">不是正定矩阵（|σ₁₂| 太大），不是合法的协方差</span>' : ""}`);
        }
        draw();
      },
    },
  });

  /* ════════════════ 6.3–6.4 线性卡尔曼滤波 ════════════════ */
  // 统一的噪声数据（仿照书中只生成一组噪声，所有测试复用）
  function droneNoise(seed = 2023, N = 100) { const rg = CE.rng(seed); return { w: Array.from({ length: N + 1 }, () => [Math.sqrt(0.05) * rg.n(), Math.sqrt(0.05) * rg.n()]), v: Array.from({ length: N + 1 }, () => [rg.n(), rg.n()]) }; }
  Course.add({
    id: "v2-kf2", vol: 2, ord: 140, no: "6.3–6.4", chName: "第6章 卡尔曼滤波器", ref: "卷2 §6.3–6.4.2",
    title: "线性卡尔曼滤波器：预测 + 校正",
    sub: "模型算出一个估计，传感器测出一个估计，按各自的不确定性加权融合。最优增益由“最小化估计误差协方差的迹”推出。",
    goals: ["理解过程噪声 Q、测量噪声 R、先验/后验估计、误差协方差 P", "跟随推导得到卡尔曼增益 $K=P^-H^T(HP^-H^T+R)^{-1}$", "掌握五个公式的“预测—校正”循环，并会调 Q、R"],
    minutes: 40,
    html: String.raw`
<h2><span class="sec">6.3.1</span>研究模型</h2>
<div class="eq" data-no="(6.3.1)–(6.3.2)" data-title="带噪声的系统与测量">
$$\bm x_{[k]}=A\bm x_{[k-1]}+B\bm u_{[k-1]}+\bm w_{[k-1]},\quad \bm w\sim N(\bm 0,Q_c);\qquad \bm z_{[k]}=H_m\bm x_{[k]}+\bm v_{[k]},\quad \bm v\sim N(\bm 0,R_c)$$
<dl><dt>$\bm w$, $Q_c$</dt><dd><b>过程噪声</b>及其协方差：模型没描述到的东西（上下坡、阵风……）——“模型有多不可信”</dd>
<dt>$\bm v$, $R_c$</dt><dd><b>测量噪声</b>及其协方差——“传感器有多不可信”</dd>
<dt>$H_m$</dt><dd>测量矩阵：测量值与状态的关系（本书推导取方阵）</dd></dl>
</div>
<p>有 $\bm w$，状态方程算不准；有 $\bm v$，测量反推不准。卡尔曼滤波 = 把两路结果融合：</p>
<ul><li><b>先验估计</b>（模型计算）：$\hat{\bm x}^-_{[k]}=A\hat{\bm x}_{[k-1]}+B\bm u_{[k-1]}$</li>
<li><b>后验估计</b>（融合测量后）：$\hat{\bm x}_{[k]}=\hat{\bm x}^-_{[k]}+K_{[k]}\big(\bm z_{[k]}-H_m\hat{\bm x}^-_{[k]}\big)$，括号里叫“测量残余/新息”。</li></ul>
<p>K→0：完全相信模型；K→$H_m^{-1}$：完全相信测量。最优的 K 是多少？</p>

<h2><span class="sec">6.3.2</span>卡尔曼增益的推导</h2>
<p>定义先验误差 $\bm e^-=\bm x-\hat{\bm x}^-$、后验误差 $\bm e=\bm x-\hat{\bm x}$，以及它们的协方差 $P^-=E(\bm e^-\bm e^{-T})$、$P=E(\bm e\bm e^T)$。$\mathrm{Tr}(P)$ = 各状态估计误差方差之和，越小越好：</p>
<details class="derive" open><summary>推导要点（用到第 2 章的矩阵求导公式 2.3.6、2.3.7）</summary><div class="body">
<div class="step">后验误差：$\bm e=(I-KH_m)\bm e^--K\bm v$。</div>
<div class="step">$\bm v$ 与 $\bm e^-$ 独立、均值为 0，交叉项消失：$P=(I-KH_m)P^-(I-KH_m)^T+KR_cK^T$。</div>
<div class="step">取迹：$\mathrm{Tr}(P)=\mathrm{Tr}(P^-)-2\mathrm{Tr}(KH_mP^-)+\mathrm{Tr}(KH_mP^-H_m^TK^T)+\mathrm{Tr}(KR_cK^T)$。</div>
<div class="step">对 K 求导（标量对矩阵，分母布局）：$-2P^-H_m^T+2K(H_mP^-H_m^T+R_c)=0$。二阶导 $2(H_mP^-H_m^T+R_c)$ 正定，确为最小。</div>
<div class="step">先验误差的传播：$\bm e^-_{[k]}=A\bm e_{[k-1]}+\bm w_{[k-1]}$，所以 $P^-_{[k]}=AP_{[k-1]}A^T+Q_c$。</div>
<div class="step">把最优 K 代回 P 的表达式，化简得 $P=(I-KH_m)P^-$。</div>
</div></details>
<div class="eq" data-title="卡尔曼滤波器的五个公式">
$$\textbf{预测：}\quad\hat{\bm x}^-_{[k]}=A\hat{\bm x}_{[k-1]}+B\bm u_{[k-1]},\qquad P^-_{[k]}=AP_{[k-1]}A^T+Q_c$$
$$\textbf{校正：}\quad K_{[k]}=P^-_{[k]}H_m^T\left(H_mP^-_{[k]}H_m^T+R_c\right)^{-1},\quad \hat{\bm x}_{[k]}=\hat{\bm x}^-_{[k]}+K_{[k]}\left(\bm z_{[k]}-H_m\hat{\bm x}^-_{[k]}\right),\quad P_{[k]}=(I-K_{[k]}H_m)P^-_{[k]}$$
<div class="say">标量情形 $H=1$ 时 $K=\frac{P^-}{P^-+R}$——正是上一课两传感器融合的 $\frac{\sigma_1^2}{\sigma_1^2+\sigma_2^2}$，只是“第一个传感器”换成了模型预测。</div>
</div>
<div data-fig="kalman"></div>
<p>使用时先给定初值 $\hat{\bm x}_{[0]}$（尽量接近真实）和 $P_{[0]}$（可以取大一点），然后不断“预测—测量—校正”。它是递归的：只需保存上一时刻的 $\hat{\bm x}$ 和 P。Q、R 为常数时，K 和 P 会收敛到常数矩阵。</p>

<h2><span class="sec">6.4</span>案例：无人机高度观测</h2>
<p>无人机以 1 m/s 匀速上升（推力等于重力），$T_s=0.1$ s。仿真中真实过程噪声 $Q_{c_a}=\mathrm{diag}(0.05,0.05)$，测量噪声 $R_{c_a}=\mathrm{diag}(1,1)$（高度和速度都能测）。实际中 R 可以从传感器参数得到，Q 却很难知道，只能“猜”。书中三组：</p>
<table class="tbl"><tr><th>组</th><th>滤波器用的 $Q_c$</th><th>书中 MSE</th><th>表现</th></tr>
<tr><td>只用测量</td><td>—</td><td>0.9948</td><td>噪声原样保留</td></tr>
<tr><td>组 1</td><td>diag(1,1) = R</td><td>0.4159</td><td>估计在先验与测量之间“各取一半”，仍然锯齿明显</td></tr>
<tr><td>组 2</td><td>diag(0.05,0.05)（恰好等于真值）</td><td>0.1572</td><td>最好：平滑且紧贴真值</td></tr>
<tr><td>组 3</td><td>diag(0.01,0.01)</td><td>0.1835</td><td>最平滑，但过于相信模型，转弯处跟不上</td></tr></table>
<div class="sim" data-id="kfuav"></div>
<div class="idea" data-t="怎样调 Q 和 R">Q/R 大 → K 大 → 跟着测量走：反应快但噪声多；Q/R 小 → K 小 → 跟着模型走：平滑但对真实变化反应迟钝。<b>实用的调参方法</b>：R 用传感器参数或实验确定；Q 先取大一点，然后不断调小，直到估计效果最好。现实中看不到真值，需要用更高精度的传感器做基准，并多次重复实验。</div>
<div class="summary"><h3>本课要点</h3><ol>
<li>五个公式：预测 $\hat x^-$、$P^-$；校正 K、$\hat x$、P。</li>
<li>K 由最小化 Tr(P) 得到，本质是按不确定性加权融合模型与测量。</li>
<li>Q 表示对模型的不信任，R 表示对传感器的不信任；Q/R 的比值决定平滑与跟踪的权衡。</li>
</ol></div>`,
    sims: {
      kfuav(el) {
        const ui = card(el, { title: "无人机高度：卡尔曼滤波实验台", hint: "灰色虚线是真实高度（仿真中才看得到），橙点是测量值，紫圈是先验估计，蓝线是后验估计。噪声数据固定，方便对比不同 Q。" });
        const ph = new Plot(ui.plotHost(), { height: 280, x: [0, 100], y: [-2, 13], xlabel: "k", ylabel: "高度 x₁" });
        const pk = new Plot(ui.plotHost(), { height: 150, x: [0, 100], y: [0, 1], xlabel: "k", title: "卡尔曼增益 K₁₁ 的收敛" });
        const sQ = ui.slider({ label: "滤波器的 Q（对角元）", min: 0.001, max: 5, value: 0.05, log: true, fmt: v => v.toFixed(3), onInput: draw });
        const sR = ui.slider({ label: "滤波器的 R（对角元）", min: 0.05, max: 10, value: 1, log: true, onInput: draw });
        ui.buttons([["组1 Q=1", () => { sQ.set(1); draw(); }, "ghost"], ["组2 Q=0.05", () => { sQ.set(0.05); draw(); }, "ghost"], ["组3 Q=0.01", () => { sQ.set(0.01); draw(); }, "ghost"]]);
        let seed = 2023; ui.buttons([["🎲 换一组噪声", () => { seed++; draw(); }, "ghost"]]);
        const cP = ui.check({ label: "显示先验估计", value: true, onChange: draw });
        const Ts = 0.1, A = [[1, Ts], [0, 1]];
        function draw() {
          const nz = droneNoise(seed), Q = M.scale(M.eye(2), sQ.get()), R = M.scale(M.eye(2), sR.get());
          let x = [0, 1], xh = [0, 1], Pk = M.eye(2); const T = [], Z = [], XM = [], XP = [], K = [];
          let mse = 0, mseZ = 0;
          for (let k = 1; k <= 100; k++) {
            x = M.vadd(M.mulv(A, x), nz.w[k]); const z = M.vadd(x, nz.v[k]);
            const s = kfStep(A, [0, 0], M.eye(2), Q, R, xh, Pk, z); xh = s.xp; Pk = s.Pp;
            T.push(x[0]); Z.push(z[0]); XM.push(s.xm[0]); XP.push(s.xp[0]); K.push(s.K[0][0]);
            mse += (x[0] - xh[0]) ** 2; mseZ += (x[0] - z[0]) ** 2;
          }
          const ks = T.map((_, i) => i + 1);
          ph.render = () => { ph.clear(); Z.forEach((v, i) => ph.dot(i + 1, v, COL.c2 + "bb", 2.6)); if (cP.get()) XM.forEach((v, i) => ph.circle(i + 1, v, COL.c5 + "99", 2.6, 1)); ph.line(ks, T, COL.ink2, 1.6, [5, 4]); ph.line(ks, XP, COL.c1, 2.6); ph.legend([{ label: "真值", color: COL.ink2, dash: [5, 4] }, { label: "后验估计", color: COL.c1 }], "tr"); };
          pk.render = () => { pk.clear(); pk.line(ks, K, COL.c4, 2.2); };
          ph.render(); pk.render();
          ui.out(`滤波 MSE：<b class="big">${f2(mse / 100, 4)}</b><br>只用测量 MSE：${f2(mseZ / 100, 4)}<br>稳态增益 $K=${mtex(s2(Q, R), 3)}$<br><span class="sim-note">（噪声由程序随机生成，数值与书中不完全相同，但三组的优劣趋势一致。）</span>`);
        }
        function s2(Q, R) { let Pk = M.eye(2), K; for (let i = 0; i < 300; i++) { const s = kfStep(A, [0, 0], M.eye(2), Q, R, [0, 0], Pk, [0, 0]); Pk = s.Pp; K = s.K; } return K; }
        draw();
      },
    },
  });

  /* ════════════════ 6.4.3 KF + MPC，6.5 EKF ════════════════ */
  Course.add({
    id: "v2-kf3", vol: 2, ord: 150, no: "6.4.3–6.5", chName: "第6章 卡尔曼滤波器", ref: "卷2 §6.4.3–6.6",
    title: "卡尔曼滤波 + MPC，以及扩展卡尔曼滤波 EKF",
    sub: "把估计出的状态交给控制器——全书知识的大融合。再把卡尔曼滤波推广到非线性系统：在每一步都线性化。",
    goals: ["理解“估计 + 控制”结构：用卡尔曼滤波的输出作为 MPC 的反馈", "理解采样时间、Q 对控制量抖动的影响", "掌握 EKF：非线性模型做预测、雅可比矩阵算协方差"],
    minutes: 35,
    html: String.raw`
<h2><span class="sec">6.4.3</span>卡尔曼滤波器 + MPC：全书大融合</h2>
<p>第 4、5 章的 LQR、MPC 都假设状态能精确获得。现实中只能拿到带噪声的测量。于是：<b>卡尔曼滤波器实时估计状态，估计值作为 MPC 的反馈</b>。用到的模块几乎是全书的全部：[F2] 稳态非零参考增广、[F4] 代价矩阵、[F6] 约束矩阵、[F7] 含约束 QP、[F8] 卡尔曼滤波器。</p>
<div class="sim" data-id="kfmpc"></div>
<ul><li>以估计值为反馈时，控制量会出现明显的高频抖动——噪声“漏”进了控制器，对执行机构不利。</li>
<li><b>改进 1：减小 Q</b>（更相信模型，进一步平滑估计）：输入更平稳，但估计准确度略有下降（最终高度可能偏离 10 m）。</li>
<li><b>改进 2：增大采样时间</b>（$T_s=0.5$）：输入平滑得多，但响应变慢、精度下降。</li>
<li>估计值有误差时，状态约束往往要放宽，否则 MPC 可能无解。</li></ul>
<div class="tip" data-t="快速响应 vs 抑制噪声">采样越快，控制越及时，但对噪声越敏感；滤波越强，越平滑，但反应越慢。追求快速、精确跟踪的系统倾向短采样、弱滤波；噪声强、要求稳定的系统倾向长采样、强滤波。</div>

<h2><span class="sec">6.5</span>扩展卡尔曼滤波器（EKF）</h2>
<p>线性卡尔曼滤波依赖一个事实：正态分布经过线性变换仍是正态分布。非线性变换会把它“扭曲”，于是不能直接用。EKF 的办法：<b>在当前估计点附近做一阶泰勒展开</b>（附录 A 的线性化），每一步都重新线性化。</p>
<div class="eq" data-no="(6.5.1)–(6.5.7)" data-title="扩展卡尔曼滤波器">
$$\bm x_{[k]}=f(\bm x_{[k-1]},\bm u_{[k-1]},\bm w_{[k-1]}),\quad \bm z_{[k]}=h(\bm x_{[k]},\bm v_{[k]})$$
$$\hat{\bm x}^-_{[k]}=f(\hat{\bm x}_{[k-1]},\bm u_{[k-1]},\bm 0),\qquad P^-_{[k]}=A_{[k]}P_{[k-1]}A_{[k]}^T+W_{[k]}Q_cW_{[k]}^T$$
$$K_{[k]}=P^-H_{m[k]}^T\left(H_{m[k]}P^-H_{m[k]}^T+V_{[k]}R_cV_{[k]}^T\right)^{-1},\quad \hat{\bm x}_{[k]}=\hat{\bm x}^-+K\left(\bm z-h(\hat{\bm x}^-,\bm 0)\right),\quad P=(I-KH_{m[k]})P^-$$
<dl><dt>$A_{[k]}=\frac{\partial f}{\partial\bm x}$, $W_{[k]}=\frac{\partial f}{\partial\bm w}$</dt><dd>在 $\hat{\bm x}_{[k-1]}$ 处的雅可比矩阵（分子布局）</dd><dt>$H_{m[k]}=\frac{\partial h}{\partial\bm x}$, $V_{[k]}=\frac{\partial h}{\partial\bm v}$</dt><dd>在 $\hat{\bm x}^-_{[k]}$ 处的雅可比矩阵</dd></dl>
<div class="say">估计值用<b>非线性模型</b>计算；协方差和增益用<b>线性化模型</b>计算。雅可比矩阵每步都要重算，计算量更大。</div>
</div>
<p><b>案例：单摆</b>。$\theta_{[k]}=\theta_{[k-1]}+\omega_{[k-1]}T_s$，$\omega_{[k]}=\omega_{[k-1]}-\frac gL\sin\theta_{[k-1]}T_s$，测量角度与角速度。唯一的非线性项是 $\sin\theta$，所以</p>
<div class="eq" data-no="(6.5.12)">
$$A_{[k]}=\begin{bmatrix}1&T_s\\-\frac gL\cos\hat\theta_{[k-1]}\,T_s&1\end{bmatrix},\qquad W=H_m=V=I$$
<div class="say">书中参数：g/L=20，$T_s=0.01$，$Q_c=\mathrm{diag}(0.01,0.01)$，$R_c=\mathrm{diag}(0.1,0.1)$，真实初值 $[\pi/4,0]$，估计初值 $[0,0]$，$P_0=I$。</div>
</div>
<div class="sim" data-id="ekf"></div>
<div class="summary"><h3>本课要点</h3><ol>
<li>估计 + 控制：卡尔曼滤波的输出作为控制器反馈；需在估计精度、控制平滑、响应速度之间权衡。</li>
<li>EKF：用非线性模型预测状态，用雅可比矩阵传播协方差和计算增益，每步重新线性化。</li>
<li>恭喜！从传递函数到状态空间，从 PID 到 LQR、MPC，从观测器到卡尔曼滤波——你已经走完了《控制之美》两卷的全部旅程。</li>
</ol></div>`,
    sims: {
      kfmpc(el) {
        const ui = card(el, { title: "估计 + 控制：无人机在噪声中飞到 10 m", hint: "真实系统有过程噪声，测量有噪声；卡尔曼滤波估计高度和速度，含约束 MPC（Np=20，推力 7～12）用估计值做反馈。" });
        const ph = new Plot(ui.plotHost(), { height: 240, x: [0, 10], y: [-1, 13], xlabel: "t / s", ylabel: "高度" });
        const pu = new Plot(ui.plotHost(), { height: 180, x: [0, 10], y: [6, 13], xlabel: "t / s", ylabel: "推力 u" });
        let Ts = 0.1;
        ui.seg({ label: "采样时间 Ts", options: [[0.1, "0.1 s"], [0.5, "0.5 s"]], value: 0.1, onChange: v => { Ts = v; draw(); } });
        const sQ = ui.slider({ label: "滤波器的 Q", min: 0.005, max: 1, value: 0.05, log: true, fmt: v => v.toFixed(3), onInput: draw });
        let seed = 7; ui.buttons([["🎲 换一组噪声", () => { seed++; draw(); }, "ghost"]]);
        const cT = ui.check({ label: "直接用测量值做反馈（不用滤波）", value: false, onChange: draw });
        function draw() {
          const A3 = [[1, Ts, Ts * Ts / 2], [0, 1, Ts], [0, 0, 1]], B3 = [[Ts * Ts / 2], [Ts], [0]], I3 = M.eye(3), W = [[1, 0, 0], [0, 1, 0], [0, 0, 0]];
          const Aa = M.vcat(M.hcat(A3, M.sub(I3, A3)), M.hcat(M.zeros(3, 3), I3)), Ba = M.vcat(B3, M.zeros(3, 1)), Ca = M.hcat(I3, M.scale(I3, -1));
          const Qa = M.mul(M.mul(M.T(Ca), W), Ca), mm = mpcMatrices(Aa, Ba, Qa, [[0.1]], Qa, 20), INF = Infinity;
          const cs = mpcConstraints(mm, [-3], [2], new Array(6).fill(-INF), new Array(6).fill(INF)), sol = qpPrep(mm.H, cs.Mc);
          const A2 = [[1, Ts], [0, 1]], Bq = [Ts * Ts / 2, Ts], rg = CE.rng(seed), Qf = M.scale(M.eye(2), sQ.get()), Rf = M.eye(2);
          let x = [0, 0], xh = [0, 0], Pk = M.eye(2), u = 10; const o = { t: [], h: [], hz: [], he: [], u: [] };
          for (let k = 0; k * Ts <= 10 + 1e-9; k++) {
            const z = [x[0] + rg.n(), x[1] + rg.n()];
            if (k > 0) { const s = kfStep(A2, Bq.map(v => v * (u - 10)), M.eye(2), Qf, Rf, xh, Pk, z); xh = s.xp; Pk = s.Pp; }
            const fb = cT.get() ? z : xh, xa = [fb[0], fb[1], -10, 10, 0, -10];
            u = clamp(sol(M.mulv(mm.F, xa), cs.rhs(xa)).x[0] + 10, 7, 12);
            o.t.push(k * Ts); o.h.push(x[0]); o.hz.push(z[0]); o.he.push(xh[0]); o.u.push(u);
            x = [x[0] + Ts * x[1] + Ts * Ts / 2 * (u - 10) + Math.sqrt(0.05) * rg.n(), x[1] + Ts * (u - 10) + Math.sqrt(0.05) * rg.n()];
          }
          const st = [], su = []; o.t.forEach((t, i) => { st.push(t, t + Ts); su.push(o.u[i], o.u[i]); });
          ph.render = () => { ph.clear(); ph.hline(10, COL.c3, [6, 4]); o.hz.forEach((v, i) => ph.dot(o.t[i], v, COL.c2 + "99", 2.4)); ph.line(o.t, o.h, COL.ink2, 1.6, [5, 4]); if (!cT.get()) ph.line(o.t, o.he, COL.c1, 2.4); ph.legend([{ label: "真值", color: COL.ink2, dash: [5, 4] }, { label: "后验估计", color: COL.c1 }], "tr"); };
          pu.render = () => { pu.clear(); pu.rect(0, 7, 10, 12, COL.c3 + "18"); pu.line(st, su, COL.c4, 1.8); };
          ph.render(); pu.render();
          const tail = o.u.filter((_, i) => o.t[i] > 5), tv = tail.reduce((s, v, i) => i ? s + Math.abs(v - tail[i - 1]) : 0, 0) / Math.max(1, tail.length - 1);
          ui.out(`5 s 后推力平均跳变幅度：<b class="big">${f2(tv, 2)}</b><br>最终真实高度：${f2(o.h[o.h.length - 1], 2)} m<br><span class="sim-note">试试：减小 Q → 推力更平稳；Ts=0.5 → 更平滑但更慢；不用滤波 → 推力剧烈抖动。</span>`);
        }
        draw();
      },
      ekf(el) {
        const ui = card(el, { title: "扩展卡尔曼滤波：单摆", hint: "真实单摆从 π/4 释放（过程噪声使摆动不完全规则）；EKF 的初始估计是 0，很快收敛到真值附近。" });
        const pa = new Plot(ui.plotHost(), { height: 260, x: [0, 200], y: [-2.5, 2.5], xlabel: "k（Ts=0.01 s）", ylabel: "角度 θ" });
        const sc = new Scene(ui.plotHost(), { height: 150 });
        const sR = ui.slider({ label: "滤波器的 R", min: 0.01, max: 2, value: 0.1, log: true, fmt: v => v.toFixed(3), onInput: draw });
        const sQ = ui.slider({ label: "滤波器的 Q", min: 0.0005, max: 0.5, value: 0.01, log: true, fmt: v => v.toFixed(4), onInput: draw });
        const cL = ui.check({ label: "对比：用固定线性化（在 θ=0 处）的普通 KF", value: false, onChange: draw });
        let seed = 3; ui.buttons([["🎲 换一组噪声", () => { seed++; draw(); }, "ghost"]]);
        const Ts = 0.01, gl = 20; let res = null, t0 = 0;
        function draw() {
          const rg = CE.rng(seed); let x = [Math.PI / 4, 0], xh = [0, 0], Pk = M.eye(2), xl = [0, 0], Pl = M.eye(2);
          const Q = M.scale(M.eye(2), sQ.get()), R = M.scale(M.eye(2), sR.get()), o = { k: [], x: [], z: [], e: [], l: [] }; let mse = 0, msel = 0;
          const f = s => [s[0] + s[1] * Ts, s[1] - gl * Math.sin(s[0]) * Ts];
          for (let k = 0; k <= 200; k++) {
            const z = [x[0] + Math.sqrt(0.1) * rg.n(), x[1] + Math.sqrt(0.1) * rg.n()];
            if (k > 0) {
              const Ak = [[1, Ts], [-gl * Math.cos(xh[0]) * Ts, 1]], xm = f(xh), Pm = M.add(M.mul(M.mul(Ak, Pk), M.T(Ak)), Q);
              const K = M.mul(Pm, M.inv(M.add(Pm, R))); xh = M.vadd(xm, M.mulv(K, M.vsub(z, xm))); Pk = M.mul(M.sub(M.eye(2), K), Pm);
              const A0 = [[1, Ts], [-gl * Ts, 1]], s = kfStep(A0, [0, 0], M.eye(2), Q, R, xl, Pl, z); xl = s.xp; Pl = s.Pp;
            }
            o.k.push(k); o.x.push(x[0]); o.z.push(z[0]); o.e.push(xh[0]); o.l.push(xl[0]);
            if (k > 20) { mse += (x[0] - xh[0]) ** 2; msel += (x[0] - xl[0]) ** 2; }
            x = f(x); x = [x[0] + Math.sqrt(0.01) * rg.n() * 0.3, x[1] + Math.sqrt(0.01) * rg.n()];
          }
          res = o;
          pa.render = () => { pa.clear(); o.z.forEach((v, i) => pa.dot(i, v, COL.c2 + "88", 2)); pa.line(o.k, o.x, COL.ink2, 1.6, [5, 4]); if (cL.get()) pa.line(o.k, o.l, COL.c5, 2); pa.line(o.k, o.e, COL.c1, 2.6); pa.legend([{ label: "真值", color: COL.ink2, dash: [5, 4] }, { label: "EKF", color: COL.c1 }, ...(cL.get() ? [{ label: "固定线性化 KF", color: COL.c5 }] : [])], "tr"); };
          pa.render();
          ui.out(`EKF 估计 MSE（k&gt;20）：<b class="big">${f2(mse / 180, 4)}</b>${cL.get() ? `<br>固定线性化 KF：${f2(msel / 180, 4)}（大角度时模型不准）` : ""}`);
        }
        draw();
        CE.loop(ui.root, dt => {
          t0 = (t0 + dt * 60) % 200; const i = Math.floor(t0), c = sc.ctx, W = sc.W; sc.clear();
          const draw1 = (cx, th, col, lab) => { const L = 110, bx = cx + L * Math.sin(th), by = 18 + L * Math.cos(th); c.strokeStyle = col; c.lineWidth = 2.5; c.beginPath(); c.moveTo(cx, 18); c.lineTo(bx, by); c.stroke(); c.fillStyle = col; c.beginPath(); c.arc(bx, by, 10, 0, 7); c.fill(); c.fillStyle = COL.ink2; c.font = "12px 'Microsoft YaHei'"; c.textAlign = "center"; c.fillText(lab, cx, 12); };
          c.fillStyle = COL.axis; c.fillRect(0, 14, W, 4);
          draw1(W * 0.3, res.x[i], COL.ink2, "真实"); draw1(W * 0.7, res.e[i], COL.c1, "EKF 估计");
        });
      },
    },
  });
})();
