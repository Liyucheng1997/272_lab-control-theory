/* 交互仿真集合。每个仿真: (mount) => void，通过 SimEngine 注册生命周期 */
window.Sims = (function () {
  const E = window.SimEngine, num = E.num;

  /* ── 公共工具 ── */
  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0; a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function gauss(rng) {
    let u = 0, v = 0;
    while (u === 0) u = rng(); while (v === 0) v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }
  // 传递函数 -> 能控标准型状态空间
  function tf2ss(numC, denC) {
    const a0 = denC[0];
    const den = denC.map(v => v / a0);
    const n = den.length - 1;
    const nm = new Array(den.length - numC.length).fill(0).concat(numC.map(v => v / a0));
    const D = nm[0];
    const A = num.zeros(n, n);
    for (let i = 0; i < n - 1; i++) A[i][i + 1] = 1;
    for (let j = 0; j < n; j++) A[n - 1][j] = -den[n - j];
    const B = new Array(n).fill(0); B[n - 1] = 1;
    const C = new Array(n);
    for (let j = 0; j < n; j++) C[j] = nm[n - j] - den[n - j] * D;
    return { A, B, C, D, n };
  }
  function stepResponse(ss, T, h, uFn) {
    const f = (x, t) => num.vadd(num.mulv(ss.A, x), num.vscale(ss.B, uFn ? uFn(t) : 1));
    let x = new Array(ss.n).fill(0);
    const ts = [], ys = [];
    for (let t = 0; t <= T; t += h) {
      ts.push(t); ys.push(num.dot(ss.C, x) + ss.D * (uFn ? uFn(t) : 1));
      x = num.rk4(f, x, t, h);
    }
    return { ts, ys };
  }
  // 频率响应幅相 + 裕度
  function margins(numC, denC) {
    const ws = [], mags = [], phs = [];
    let prevPh = null, unwrap = 0;
    for (let lw = -2; lw <= 3; lw += 0.005) {
      const w = Math.pow(10, lw);
      const [re, im] = num.freqResp(numC, denC, w);
      let ph = Math.atan2(im, re) * 180 / Math.PI;
      if (prevPh !== null) {
        while (ph + unwrap - prevPh > 180) unwrap -= 360;
        while (ph + unwrap - prevPh < -180) unwrap += 360;
      }
      ph += unwrap; prevPh = ph;
      ws.push(w); mags.push(Math.sqrt(re * re + im * im)); phs.push(ph);
    }
    let wpc = null, gm = Infinity, wgc = null, pm = null;
    for (let i = 1; i < ws.length; i++) {
      if (phs[i - 1] > -180 && phs[i] <= -180 && wpc === null) {
        wpc = ws[i]; gm = 1 / mags[i];
      }
      if (mags[i - 1] >= 1 && mags[i] < 1 && wgc === null) {
        wgc = ws[i]; pm = 180 + phs[i];
      }
    }
    return { ws, mags, phs, gm, pm, wpc, wgc };
  }

  /* ═══════════ 1. 特征值与相图 ═══════════ */
  function simEigPhase(mount) {
    const ui = E.buildPanel(mount, {
      title: "二维线性系统相图实验",
      desc: "调整矩阵 A 的元素，观察轨迹形态与特征值的关系：实部<0 收敛（稳定焦点/结点），实部>0 发散，纯虚数为中心。"
    });
    const cv = E.addCanvas(ui, 420, 420);
    const cv2 = E.addCanvas(ui, 320, 420);
    const p = new E.Plot(cv, { xmin: -3, xmax: 3, ymin: -3, ymax: 3, xlabel: "x₁", ylabel: "x₂", title: "相平面轨迹" });
    const pe = new E.Plot(cv2, { xmin: -3, xmax: 3, ymin: -3, ymax: 3, xlabel: "Re", ylabel: "Im", title: "特征值位置（s 平面）" });
    const s = {};
    const redraw = () => {
      const A = [[s.a11.get(), s.a12.get()], [s.a21.get(), s.a22.get()]];
      p.clear();
      // 向量场
      for (let i = -2.7; i <= 2.7; i += 0.6) for (let j = -2.7; j <= 2.7; j += 0.6) {
        const d = num.mulv(A, [i, j]); const n0 = Math.hypot(d[0], d[1]) || 1;
        const len = 0.22;
        p.line([i, i + d[0] / n0 * len], [j, j + d[1] / n0 * len], "#ccd6e2", 1);
        p.dot(i + d[0] / n0 * len, j + d[1] / n0 * len, "#aebbcb", 1.2);
      }
      // 多条轨迹
      const colors = ["#5b8def", "#ef5b7b", "#31b5a9", "#e8a838", "#9d6bff", "#7bc74d", "#e06bd8", "#8ab8ff"];
      for (let k = 0; k < 8; k++) {
        const th = k / 8 * 2 * Math.PI;
        let x = [2.4 * Math.cos(th), 2.4 * Math.sin(th)];
        const xs = [x[0]], ys = [x[1]];
        for (let t = 0; t < 6; t += 0.02) {
          x = num.rk4((x) => num.mulv(A, x), x, t, 0.02);
          if (Math.hypot(x[0], x[1]) > 40) break;
          xs.push(x[0]); ys.push(x[1]);
        }
        p.line(xs, ys, colors[k], 1.6);
        p.dot(xs[0], ys[0], colors[k], 3);
      }
      // 特征值
      const tr = A[0][0] + A[1][1], det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
      const disc = tr * tr - 4 * det;
      let eigs, kind;
      if (disc >= 0) {
        const r = Math.sqrt(disc);
        eigs = [[(tr + r) / 2, 0], [(tr - r) / 2, 0]];
      } else {
        const r = Math.sqrt(-disc);
        eigs = [[tr / 2, r / 2], [tr / 2, -r / 2]];
      }
      const maxRe = Math.max(eigs[0][0], eigs[1][0]);
      kind = maxRe < -1e-9 ? "渐近稳定 ✓" : maxRe > 1e-9 ? "不稳定 ✗" : "临界（中心/退化）";
      pe.clear();
      eigs.forEach(ev => pe.cross(ev[0], ev[1], "#ef5b7b", 6));
      ui.readout.textContent =
        `tr(A) = ${tr.toFixed(2)}   det(A) = ${det.toFixed(2)}\n` +
        `特征值: λ₁ = ${eigs[0][0].toFixed(2)}${eigs[0][1] ? (eigs[0][1] > 0 ? "+" : "") + eigs[0][1].toFixed(2) + "j" : ""},  λ₂ = ${eigs[1][0].toFixed(2)}${eigs[1][1] ? (eigs[1][1] > 0 ? "+" : "") + eigs[1][1].toFixed(2) + "j" : ""}\n` +
        `判定: ${kind}`;
    };
    s.a11 = E.addSlider(ui, { label: "a₁₁", min: -3, max: 3, step: 0.1, value: -0.5, onChange: redraw });
    s.a12 = E.addSlider(ui, { label: "a₁₂", min: -3, max: 3, step: 0.1, value: 2, onChange: redraw });
    s.a21 = E.addSlider(ui, { label: "a₂₁", min: -3, max: 3, step: 0.1, value: -2, onChange: redraw });
    s.a22 = E.addSlider(ui, { label: "a₂₂", min: -3, max: 3, step: 0.1, value: -0.5, onChange: redraw });
    redraw();
  }

  /* ═══════════ 2. 二阶系统阶跃响应 ═══════════ */
  function simSecondOrder(mount) {
    const ui = E.buildPanel(mount, {
      title: "标准二阶系统实验",
      desc: "G(s) = ωₙ²/(s²+2ζωₙs+ωₙ²)。拖动 ζ 和 ωₙ，体会阻尼比决定超调、自然频率决定快慢。"
    });
    const cv = E.addCanvas(ui, 480, 320);
    const cv2 = E.addCanvas(ui, 300, 320);
    const s = {};
    const redraw = () => {
      const z = s.zeta.get(), wn = s.wn.get();
      const ss = tf2ss([wn * wn], [1, 2 * z * wn, wn * wn]);
      const { ts, ys } = stepResponse(ss, 12, 0.01);
      const plt = new E.Plot(cv, { xmin: 0, xmax: 12, ymin: 0, ymax: 2, xlabel: "t (s)", ylabel: "y(t)", title: "单位阶跃响应" });
      plt.clear();
      plt.line([0, 12], [1, 1], "#3b4c66", 1, [4, 4]);
      plt.line(ts, ys, "#5b8def", 2.2);
      const pp = new E.Plot(cv2, { xmin: -8, xmax: 2, ymin: -8, ymax: 8, xlabel: "Re", ylabel: "Im", title: "闭环极点" });
      pp.clear();
      const roots = num.polyRoots([1, 2 * z * wn, wn * wn]);
      roots.forEach(r => pp.cross(r[0], r[1], "#ef5b7b", 6));
      const Mp = z < 1 ? Math.exp(-Math.PI * z / Math.sqrt(1 - z * z)) * 100 : 0;
      const tsapprox = z > 0 ? 4 / (z * wn) : Infinity;
      ui.readout.textContent =
        `极点: s = ${(-z * wn).toFixed(2)} ± ${z < 1 ? (wn * Math.sqrt(1 - z * z)).toFixed(2) : 0}j\n` +
        `理论超调量 Mp ≈ ${Mp.toFixed(1)}%   2% 调节时间 ts ≈ ${isFinite(tsapprox) ? tsapprox.toFixed(2) + " s" : "∞"}`;
    };
    s.zeta = E.addSlider(ui, { label: "阻尼比 ζ", min: 0.02, max: 2, step: 0.01, value: 0.3, onChange: redraw });
    s.wn = E.addSlider(ui, { label: "自然频率 ωₙ", min: 0.5, max: 6, step: 0.1, value: 2, onChange: redraw });
    redraw();
  }

  /* ═══════════ 3. PID ═══════════ */
  function simPID(mount) {
    const ui = E.buildPanel(mount, {
      title: "PID 调参实验台",
      desc: "被控对象 G(s)=1/(s²+as+b)。t=0 给单位阶跃设定值，t=10s 加入输入端阶跃扰动 d=-0.5。微分作用于测量值（避免设定值突变引起微分冲击）。"
    });
    const cv = E.addCanvas(ui, 500, 300);
    const cvU = E.addCanvas(ui, 500, 180);
    const s = {};
    const redraw = () => {
      const Kp = s.kp.get(), Ki = s.ki.get(), Kd = s.kd.get();
      const a = s.a.get(), b = s.b.get();
      const N = 20; // 微分滤波系数
      // 状态: [y, ydot, xi(积分), xd(微分滤波器状态)]
      const T = 20, h = 0.002;
      const ts = [], ys = [], us = [];
      let x = [0, 0, 0, 0];
      for (let t = 0; t <= T; t += h) {
        const r = 1, d = t >= 10 ? -0.5 : 0;
        const e = r - x[0];
        // 滤波微分 dy/dt ≈ N(y - xd)
        const dyf = N * (x[0] - x[3]);
        let u = Kp * e + Ki * x[2] - Kd * dyf;
        u = Math.max(-10, Math.min(10, u)); // 执行器饱和
        if (ts.length === 0 || t - ts[ts.length - 1] >= 0.02 - 1e-9) { ts.push(t); ys.push(x[0]); us.push(u); }
        const f = (x) => [
          x[1],
          -a * x[1] - b * x[0] + u + d,
          e,
          N * (x[0] - x[3]),
        ];
        x = num.rk4(f, x, t, h);
      }
      const plt = new E.Plot(cv, { xmin: 0, xmax: 20, ymin: -0.2, ymax: 2, xlabel: "t (s)", ylabel: "y", title: "输出响应（10s 处加扰动）" });
      plt.autoRangeY([ys, [0, 1.6]]);
      plt.clear();
      plt.line([0, 20], [1, 1], "#3b4c66", 1, [4, 4]);
      plt.line(ts, ys, "#5b8def", 2);
      const pu = new E.Plot(cvU, { xmin: 0, xmax: 20, ymin: -2, ymax: 2, xlabel: "t (s)", ylabel: "u", title: "控制量" });
      pu.autoRangeY([us]);
      pu.clear();
      pu.line(ts, us, "#e8a838", 1.6);
      // 性能指标
      let peak = 0, settle = NaN;
      ys.forEach(v => peak = Math.max(peak, v));
      for (let i = ts.length - 1; i >= 0; i--) {
        if (ts[i] < 10 && Math.abs(ys[i] - 1) > 0.02) { settle = ts[i]; break; }
      }
      ui.readout.textContent =
        `超调量 ≈ ${((peak - 1) * 100).toFixed(1)}%   调节时间(2%) ≈ ${isFinite(settle) ? settle.toFixed(2) + " s" : "—"}\n` +
        `尝试: 只用 P（有静差）→ 加 I（消静差但振荡）→ 加 D（阻尼振荡）`;
    };
    s.kp = E.addSlider(ui, { label: "Kp", min: 0, max: 30, step: 0.1, value: 5, onChange: redraw });
    s.ki = E.addSlider(ui, { label: "Ki", min: 0, max: 20, step: 0.1, value: 0, onChange: redraw });
    s.kd = E.addSlider(ui, { label: "Kd", min: 0, max: 10, step: 0.1, value: 0, onChange: redraw });
    s.a = E.addSlider(ui, { label: "对象阻尼 a", min: 0, max: 5, step: 0.1, value: 1, onChange: redraw });
    s.b = E.addSlider(ui, { label: "对象刚度 b", min: 0, max: 5, step: 0.1, value: 1, onChange: redraw });
    redraw();
  }

  /* ═══════════ 4. 极点配置 ═══════════ */
  function simPolePlace(mount) {
    const ui = E.buildPanel(mount, {
      title: "极点配置实验：双积分器",
      desc: "对象 ẍ=u（如无摩擦小车）。选定期望闭环极点（由 ζ、ωₙ 参数化），计算反馈增益 K=[ωₙ², 2ζωₙ]，从 x₀=[1,0] 调节回原点。"
    });
    const cv = E.addCanvas(ui, 300, 300);
    const cv2 = E.addCanvas(ui, 480, 300);
    const s = {};
    const redraw = () => {
      const z = s.zeta.get(), wn = s.wn.get();
      const k1 = wn * wn, k2 = 2 * z * wn;
      const pp = new E.Plot(cv, { xmin: -10, xmax: 2, ymin: -6, ymax: 6, xlabel: "Re", ylabel: "Im", title: "闭环极点（叉）" });
      pp.clear();
      num.polyRoots([1, k2, k1]).forEach(r => pp.cross(r[0], r[1], "#ef5b7b", 6));
      pp.dot(0, 0, "#8b98ab", 3); pp.text(0.15, 0.4, "开环双重极点", "#8b98ab");
      const T = 8, h = 0.01;
      let x = [1, 0];
      const ts = [], x1s = [], x2s = [], us = [];
      for (let t = 0; t <= T; t += h) {
        const u = -k1 * x[0] - k2 * x[1];
        ts.push(t); x1s.push(x[0]); x2s.push(x[1]); us.push(u);
        x = num.rk4((x) => [x[1], -k1 * x[0] - k2 * x[1]], x, t, h);
      }
      const pt = new E.Plot(cv2, { xmin: 0, xmax: 8, ymin: -2, ymax: 1.5, xlabel: "t (s)", title: "状态与控制量" });
      pt.autoRangeY([x1s, x2s, us]);
      pt.clear();
      pt.line(ts, x1s, "#5b8def", 2);
      pt.line(ts, x2s, "#31b5a9", 1.8);
      pt.line(ts, us, "#e8a838", 1.4, [5, 3]);
      pt.legend([{ label: "位置 x₁", color: "#5b8def" }, { label: "速度 x₂", color: "#31b5a9" }, { label: "u", color: "#e8a838", dash: [5, 3] }]);
      ui.readout.textContent = `K = [${k1.toFixed(2)}, ${k2.toFixed(2)}]   期望特征式: s² + ${k2.toFixed(2)}s + ${k1.toFixed(2)}\n极点越靠左收敛越快，但 u 峰值越大（高增益代价）`;
    };
    s.zeta = E.addSlider(ui, { label: "ζ", min: 0.1, max: 1.5, step: 0.05, value: 0.7, onChange: redraw });
    s.wn = E.addSlider(ui, { label: "ωₙ", min: 0.5, max: 8, step: 0.1, value: 2, onChange: redraw });
    redraw();
  }

  /* ═══════════ 5. 状态观测器 ═══════════ */
  function simObserver(mount) {
    const ui = E.buildPanel(mount, {
      title: "Luenberger 观测器实验",
      desc: "双积分器只测位置 y=x₁。观测器从错误初值 [0,0] 出发估计真实状态 [1, -0.5]。观测器极点越快（越靠左），收敛越快，但对测量噪声越敏感。"
    });
    const cv = E.addCanvas(ui, 500, 300);
    const cvE = E.addCanvas(ui, 500, 200);
    const s = {};
    const redraw = () => {
      const wo = s.wo.get(), noise = s.noise.get();
      // 观测器极点配置在 -wo 双重: L = [2wo, wo²]
      const L = [2 * wo, wo * wo];
      const rng = mulberry32(42);
      const T = 6, h = 0.005;
      let x = [1, -0.5], xh = [0, 0];
      const ts = [], x1s = [], x2s = [], xh1s = [], xh2s = [];
      for (let t = 0; t <= T; t += h) {
        const u = Math.sin(1.5 * t); // 激励
        const y = x[0] + noise * gauss(rng);
        if (ts.length === 0 || t - ts[ts.length - 1] >= 0.02 - 1e-9) {
          ts.push(t); x1s.push(x[0]); x2s.push(x[1]); xh1s.push(xh[0]); xh2s.push(xh[1]);
        }
        x = num.rk4((x) => [x[1], u], x, t, h);
        const innov = y - xh[0];
        xh = num.rk4((xh) => [xh[1] + L[0] * innov, u + L[1] * innov], xh, t, h);
      }
      const p1 = new E.Plot(cv, { xmin: 0, xmax: 6, xlabel: "t (s)", title: "真实状态（实线） vs 估计（虚线）" });
      p1.autoRangeY([x1s, x2s, xh1s, xh2s]);
      p1.clear();
      p1.line(ts, x1s, "#5b8def", 2); p1.line(ts, xh1s, "#2a55b0", 1.8, [5, 4]);
      p1.line(ts, x2s, "#ef5b7b", 2); p1.line(ts, xh2s, "#c23a58", 1.8, [5, 4]);
      p1.legend([{ label: "x₁", color: "#5b8def" }, { label: "x̂₁", color: "#2a55b0", dash: [5, 4] }, { label: "x₂", color: "#ef5b7b" }, { label: "x̂₂", color: "#c23a58", dash: [5, 4] }]);
      const errs = ts.map((_, i) => Math.hypot(x1s[i] - xh1s[i], x2s[i] - xh2s[i]));
      const p2 = new E.Plot(cvE, { xmin: 0, xmax: 6, xlabel: "t (s)", ylabel: "‖x-x̂‖", title: "估计误差范数" });
      p2.autoRangeY([errs]);
      p2.clear();
      p2.line(ts, errs, "#31b5a9", 2);
      ui.readout.textContent = `观测器增益 L = [${L[0].toFixed(1)}, ${L[1].toFixed(1)}]（误差极点配置于 s=-${wo}, -${wo}）\n把噪声调大再把极点加快，看看『快收敛 vs 噪声放大』的权衡`;
    };
    s.wo = E.addSlider(ui, { label: "观测器带宽 ω₀", min: 1, max: 20, step: 0.5, value: 4, onChange: redraw });
    s.noise = E.addSlider(ui, { label: "测量噪声 σ", min: 0, max: 0.1, step: 0.005, value: 0.01, fmt: v => (+v).toFixed(3), onChange: redraw });
    redraw();
  }

  /* ═══════════ 6. Lyapunov: 阻尼单摆 ═══════════ */
  function simLyapunov(mount) {
    const ui = E.buildPanel(mount, {
      title: "阻尼单摆的能量视角",
      desc: "θ̈ = -sinθ - cθ̇。背景等高线是能量 V=½θ̇²+(1-cosθ)。有阻尼时 V̇=-cθ̇²≤0：轨迹只能『下山』。注意 V̇ 在 θ̇=0 处为零但轨迹仍收敛——这正是 LaSalle 不变原理的用武之地。"
    });
    const cv = E.addCanvas(ui, 460, 380);
    const cvV = E.addCanvas(ui, 320, 380);
    const s = {};
    let anim = null;
    const start = () => {
      const c = s.c.get();
      const th0 = s.th0.get();
      const p = new E.Plot(cv, { xmin: -7, xmax: 7, ymin: -4, ymax: 4, xlabel: "θ", ylabel: "θ̇", title: "相平面（背景=能量等高线）" });
      p.clear();
      // 能量等高线
      const ctx = p.ctx;
      for (let lev = 0.3; lev <= 4.6; lev += 0.55) {
        ctx.strokeStyle = "rgba(91,141,239,0.22)"; ctx.lineWidth = 1;
        // V = w²/2 + 1-cos(th) = lev  => w = ±sqrt(2(lev-1+cos th))
        ctx.beginPath();
        let st = false;
        for (let th = -7; th <= 7; th += 0.03) {
          const v2 = 2 * (lev - 1 + Math.cos(th));
          if (v2 < 0) { st = false; continue; }
          const w = Math.sqrt(v2);
          const px = p.x2px(th), py = p.y2px(w);
          if (!st) { ctx.moveTo(px, py); st = true; } else ctx.lineTo(px, py);
        }
        ctx.stroke();
        ctx.beginPath(); st = false;
        for (let th = -7; th <= 7; th += 0.03) {
          const v2 = 2 * (lev - 1 + Math.cos(th));
          if (v2 < 0) { st = false; continue; }
          const w = -Math.sqrt(v2);
          const px = p.x2px(th), py = p.y2px(w);
          if (!st) { ctx.moveTo(px, py); st = true; } else ctx.lineTo(px, py);
        }
        ctx.stroke();
      }
      p.dot(0, 0, "#7bc74d", 4);
      [-2 * Math.PI, 2 * Math.PI].forEach(t => p.dot(t, 0, "#7bc74d", 3));
      [-Math.PI, Math.PI].forEach(t => p.cross(t, 0, "#ef5b7b", 4));
      let x = [th0, 0], t = 0;
      const traj = [[x[0], x[1]]];
      const Vs = [], Ts = [];
      const pv = new E.Plot(cvV, { xmin: 0, xmax: 20, ymin: 0, ymax: 3, xlabel: "t (s)", ylabel: "V", title: "能量 V(t) 单调下降" });
      if (anim) anim.alive = false;
      anim = E.loop(() => {
        for (let i = 0; i < 6; i++) {
          x = num.rk4((x) => [x[1], -Math.sin(x[0]) - c * x[1]], x, t, 0.01);
          t += 0.01;
          traj.push([x[0], x[1]]);
        }
        const V = x[1] * x[1] / 2 + 1 - Math.cos(x[0]);
        Ts.push(t); Vs.push(V);
        if (traj.length > 4000) traj.splice(0, traj.length - 4000);
        // 重绘轨迹层（简单方式：整图重绘代价高，改为增量画线）
        const n = traj.length;
        if (n >= 2) {
          const c2 = p.ctx;
          c2.strokeStyle = "#e8a838"; c2.lineWidth = 2;
          c2.beginPath();
          c2.moveTo(p.x2px(traj[n - 2][0]), p.y2px(traj[n - 2][1]));
          c2.lineTo(p.x2px(traj[n - 1][0]), p.y2px(traj[n - 1][1]));
          c2.stroke();
        }
        pv.opt.xmax = Math.max(20, t);
        pv.autoRangeY([Vs, [0]]);
        pv.clear();
        pv.line(Ts, Vs, "#31b5a9", 2);
        ui.readout.textContent = `t = ${t.toFixed(1)} s   V = ${(x[1] * x[1] / 2 + 1 - Math.cos(x[0])).toFixed(4)}   V̇ = ${(-c * x[1] * x[1]).toFixed(4)} ≤ 0`;
        if (t > 40) anim.alive = false;
      });
    };
    s.c = E.addSlider(ui, { label: "阻尼 c", min: 0, max: 1.5, step: 0.05, value: 0.25, onChange: start });
    s.th0 = E.addSlider(ui, { label: "初始角 θ₀", min: -6, max: 6, step: 0.1, value: 2.8, onChange: start });
    E.addButton(ui, "重新开始", start);
    start();
  }

  /* ═══════════ 7. Nyquist ═══════════ */
  function simNyquist(mount) {
    const ui = E.buildPanel(mount, {
      title: "Nyquist 判据实验",
      desc: "L(s) = K/((s+1)(s+2)(s+3))，开环稳定（P=0）。增大 K，Nyquist 曲线膨胀；当它包围 -1 点（N>0）时闭环失稳。临界值 K=60。"
    });
    const cv = E.addCanvas(ui, 420, 380);
    const cvS = E.addCanvas(ui, 380, 380);
    const s = {};
    const redraw = () => {
      const K = s.K.get();
      const nC = [K], dC = num.polyMul(num.polyMul([1, 1], [1, 2]), [1, 3]);
      const res = [], ims = [];
      for (let lw = -2.5; lw <= 3; lw += 0.004) {
        const w = Math.pow(10, lw);
        const [re, im] = num.freqResp(nC, dC, w);
        res.push(re); ims.push(im);
      }
      const resAll = res.concat(res.slice().reverse());
      const imsAll = ims.concat(ims.slice().reverse().map(v => -v));
      const p = new E.Plot(cv, { xmin: -4, xmax: 11, ymin: -7, ymax: 7, xlabel: "Re", ylabel: "Im", title: "Nyquist 图" });
      const spanX = Math.max(3, ...resAll.map(Math.abs)) * 1.15;
      const spanY = Math.max(3, ...imsAll.map(Math.abs)) * 1.15;
      p.opt.xmin = -spanX; p.opt.xmax = spanX; p.opt.ymin = -spanY; p.opt.ymax = spanY;
      p.clear();
      p.line(resAll.slice(0, res.length), imsAll.slice(0, res.length), "#5b8def", 2);
      p.line(resAll.slice(res.length), imsAll.slice(res.length), "#5b8def", 1.4, [4, 4]);
      p.cross(-1, 0, "#ef5b7b", 6);
      p.text(-1, 0.5, "-1", "#ef5b7b", "center");
      // 绕圈数（winding number）
      let wind = 0;
      const pts = resAll.map((r, i) => [r + 1, imsAll[i]]);
      for (let i = 1; i < pts.length; i++) {
        const a1 = Math.atan2(pts[i - 1][1], pts[i - 1][0]);
        const a2 = Math.atan2(pts[i][1], pts[i][0]);
        let d = a2 - a1;
        while (d > Math.PI) d -= 2 * Math.PI;
        while (d < -Math.PI) d += 2 * Math.PI;
        wind += d;
      }
      const N = Math.round(-wind / (2 * Math.PI)); // 顺时针为正
      // 闭环极点
      const clDen = dC.slice(); clDen[clDen.length - 1] += K;
      const roots = num.polyRoots(clDen);
      const ps = new E.Plot(cvS, { xmin: -6, xmax: 3, ymin: -5, ymax: 5, xlabel: "Re", ylabel: "Im", title: "闭环极点（验证）" });
      ps.clear();
      roots.forEach(r => ps.cross(r[0], r[1], r[0] > 0 ? "#ef5b7b" : "#7bc74d", 6));
      const unstable = roots.filter(r => r[0] > 1e-6).length;
      ui.readout.textContent =
        `顺时针包围 -1 圈数 N = ${N}   开环不稳定极点 P = 0\n` +
        `Nyquist 预言闭环不稳定极点 Z = N + P = ${N}   实际求根: ${unstable} 个右半平面极点  ${N === unstable ? "✓ 一致" : ""}\n` +
        `${unstable === 0 ? "闭环稳定 ✓" : "闭环不稳定 ✗"}`;
    };
    s.K = E.addSlider(ui, { label: "增益 K", min: 1, max: 120, step: 1, value: 20, fmt: v => (+v).toFixed(0), onChange: redraw });
    redraw();
  }

  /* ═══════════ 8. 根轨迹 ═══════════ */
  function simRootLocus(mount) {
    const ui = E.buildPanel(mount, {
      title: "根轨迹实验",
      desc: "L(s) = K/((s+1)(s+2)(s+3))。K 从 0 → ∞，闭环极点从开环极点出发沿轨迹移动；两支穿越虚轴（K=60）进入右半平面。"
    });
    const cv = E.addCanvas(ui, 460, 400);
    const cvT = E.addCanvas(ui, 400, 400);
    const s = {};
    const dC = num.polyMul(num.polyMul([1, 1], [1, 2]), [1, 3]);
    // 预计算轨迹
    const locus = [];
    for (let lk = -2; lk <= 2.6; lk += 0.01) {
      const K = Math.pow(10, lk);
      const clDen = dC.slice(); clDen[clDen.length - 1] += K;
      locus.push({ K, roots: num.polyRoots(clDen) });
    }
    const redraw = () => {
      const K = s.K.get();
      const p = new E.Plot(cv, { xmin: -7, xmax: 3, ymin: -6, ymax: 6, xlabel: "Re", ylabel: "Im", title: "根轨迹（红=当前 K 的闭环极点）" });
      p.clear();
      locus.forEach(pt => pt.roots.forEach(r => p.dot(r[0], r[1], "rgba(91,141,239,0.35)", 1.1)));
      [-1, -2, -3].forEach(x0 => p.cross(x0, 0, "#e8a838", 5));
      const clDen = dC.slice(); clDen[clDen.length - 1] += K;
      const roots = num.polyRoots(clDen);
      roots.forEach(r => p.dot(r[0], r[1], "#ef5b7b", 4.5));
      // 时域响应
      const ss = tf2ss([K], clDen);
      const { ts, ys } = stepResponse(ss, 8, 0.005);
      const pt2 = new E.Plot(cvT, { xmin: 0, xmax: 8, xlabel: "t (s)", ylabel: "y", title: "对应闭环阶跃响应" });
      pt2.autoRangeY([ys, [0, 1.2]]);
      pt2.clear();
      pt2.line([0, 8], [K / (K + 6), K / (K + 6)], "#3b4c66", 1, [4, 4]);
      pt2.line(ts, ys, "#31b5a9", 2);
      const unstable = roots.some(r => r[0] > 1e-6);
      ui.readout.textContent = `K = ${K.toFixed(1)}   闭环极点: ${roots.map(r => `${r[0].toFixed(2)}${Math.abs(r[1]) > 1e-6 ? (r[1] > 0 ? "+" : "") + r[1].toFixed(2) + "j" : ""}`).join(",  ")}\n${unstable ? "已越过虚轴，闭环不稳定 ✗" : "闭环稳定 ✓"}（稳态值 = K/(K+6)，K 越大静差越小——但稳定性受限）`;
    };
    s.K = E.addSlider(ui, { label: "增益 K", min: 0.5, max: 150, step: 0.5, value: 10, fmt: v => (+v).toFixed(1), onChange: redraw });
    redraw();
  }

  /* ═══════════ 9. Bode 与裕度 ═══════════ */
  function simBode(mount) {
    const ui = E.buildPanel(mount, {
      title: "Bode 图与稳定裕度实验",
      desc: "L(s) = K/((s+1)(0.5s+1)(0.2s+1))。观察增益裕度 GM 与相位裕度 PM 如何随 K 缩水，以及裕度与超调/稳定的关系。"
    });
    const cvM = E.addCanvas(ui, 440, 240);
    const cvP = E.addCanvas(ui, 440, 240);
    const cvT = E.addCanvas(ui, 440, 240);
    const s = {};
    const redraw = () => {
      const K = s.K.get();
      const nC = [K];
      const dC = num.polyMul(num.polyMul([1, 1], [0.5, 1]), [0.2, 1]);
      const m = margins(nC, dC);
      const lws = m.ws.map(w => Math.log10(w));
      const dbs = m.mags.map(v => 20 * Math.log10(v));
      const pm1 = new E.Plot(cvM, { xmin: -2, xmax: 3, ymin: -60, ymax: 40, xlabel: "log₁₀ω", ylabel: "|L| dB", title: "幅频特性" });
      pm1.clear();
      pm1.line([-2, 3], [0, 0], "#3b4c66", 1, [4, 4]);
      pm1.line(lws, dbs, "#5b8def", 2);
      if (m.wgc) { pm1.dot(Math.log10(m.wgc), 0, "#7bc74d", 4); }
      const pm2 = new E.Plot(cvP, { xmin: -2, xmax: 3, ymin: -280, ymax: 0, xlabel: "log₁₀ω", ylabel: "∠L (°)", title: "相频特性" });
      pm2.clear();
      pm2.line([-2, 3], [-180, -180], "#ef5b7b", 1, [4, 4]);
      pm2.line(lws, m.phs, "#9d6bff", 2);
      if (m.wpc) pm2.dot(Math.log10(m.wpc), -180, "#ef5b7b", 4);
      // 闭环阶跃
      const clNum = nC, clDen = num.T ? null : null;
      const dcl = dC.slice(); // den + num（num 阶数低，加到常数项）
      const padded = new Array(dC.length).fill(0); padded[dC.length - 1] = K;
      const clD = dC.map((v, i) => v + padded[i]);
      const ss = tf2ss([K], clD);
      const { ts, ys } = stepResponse(ss, 10, 0.005);
      const pt = new E.Plot(cvT, { xmin: 0, xmax: 10, xlabel: "t (s)", ylabel: "y", title: "闭环阶跃响应" });
      pt.autoRangeY([ys, [0, 1.3]]);
      pt.clear();
      pt.line([0, 10], [K / (K + 1), K / (K + 1)], "#3b4c66", 1, [4, 4]);
      pt.line(ts, ys, "#31b5a9", 2);
      ui.readout.textContent =
        `增益裕度 GM = ${isFinite(m.gm) ? (20 * Math.log10(m.gm)).toFixed(1) + " dB（还能放大 " + m.gm.toFixed(2) + " 倍）" : "∞"}   ` +
        `相位裕度 PM = ${m.pm !== null ? m.pm.toFixed(1) + "°" : "∞"}\n` +
        `经验法则: PM ≈ 100ζ（%），PM < 30° 时响应明显振荡，PM ≤ 0 失稳（本例 Routh 临界 K = 12.6）`;
    };
    s.K = E.addSlider(ui, { label: "增益 K", min: 0.5, max: 30, step: 0.1, value: 4, onChange: redraw });
    redraw();
  }

  /* ═══════════ 10. 超前校正 ═══════════ */
  function simLeadLag(mount) {
    const ui = E.buildPanel(mount, {
      title: "超前校正实验",
      desc: "对象 G(s)=4/(s(s+1))。超前环节 C(s)=K(Ts+1)/(αTs+1)（α<1）在中频注入正相位，抬高相位裕度、加快响应。对比开/关校正的效果。"
    });
    const cvP = E.addCanvas(ui, 440, 250);
    const cvT = E.addCanvas(ui, 440, 250);
    const s = {};
    const redraw = () => {
      const K = s.K.get(), T = s.T.get(), al = s.alpha.get();
      const gN = [4], gD = [1, 1, 0];
      const cN = [K * T, K], cD = [al * T, 1];
      const LN = num.polyMul(cN, gN), LD = num.polyMul(cD, gD);
      const m0 = margins(gN, gD);
      const m1 = margins(LN, LD);
      const pm2 = new E.Plot(cvP, { xmin: -2, xmax: 2.5, ymin: -200, ymax: -60, xlabel: "log₁₀ω", ylabel: "∠L (°)", title: "相频：超前环节『托起』相位" });
      pm2.clear();
      pm2.line([-2, 2.5], [-180, -180], "#ef5b7b", 1, [4, 4]);
      pm2.line(m0.ws.map(w => Math.log10(w)), m0.phs, "#5d6b80", 1.6, [5, 4]);
      pm2.line(m1.ws.map(w => Math.log10(w)), m1.phs, "#9d6bff", 2);
      pm2.legend([{ label: "未校正", color: "#5d6b80", dash: [5, 4] }, { label: "加超前", color: "#9d6bff" }]);
      // 闭环对比
      const mk = (n0, d0) => {
        const nn = new Array(d0.length - n0.length).fill(0).concat(n0);
        const cl = d0.map((v, i) => v + nn[i]);
        return stepResponse(tf2ss(n0, cl), 10, 0.005);
      };
      const r0 = mk(gN, gD), r1 = mk(LN, LD);
      const pt = new E.Plot(cvT, { xmin: 0, xmax: 10, xlabel: "t (s)", ylabel: "y", title: "闭环阶跃响应对比" });
      pt.autoRangeY([r0.ys, r1.ys, [0, 1.2]]);
      pt.clear();
      pt.line([0, 10], [1, 1], "#3b4c66", 1, [4, 4]);
      pt.line(r0.ts, r0.ys, "#5d6b80", 1.6, [5, 4]);
      pt.line(r1.ts, r1.ys, "#31b5a9", 2.2);
      pt.legend([{ label: "未校正", color: "#5d6b80", dash: [5, 4] }, { label: "加超前", color: "#31b5a9" }]);
      const maxLead = Math.asin((1 - al) / (1 + al)) * 180 / Math.PI;
      ui.readout.textContent =
        `未校正 PM = ${m0.pm !== null ? m0.pm.toFixed(1) : "—"}°  →  校正后 PM = ${m1.pm !== null ? m1.pm.toFixed(1) : "—"}°\n` +
        `该超前环节最大可提供相位 φmax = ${maxLead.toFixed(1)}°（出现在 ω = 1/(T√α) = ${(1 / (T * Math.sqrt(al))).toFixed(2)} rad/s）`;
    };
    s.K = E.addSlider(ui, { label: "K", min: 0.2, max: 4, step: 0.1, value: 1, onChange: redraw });
    s.T = E.addSlider(ui, { label: "T", min: 0.05, max: 2, step: 0.05, value: 0.5, onChange: redraw });
    s.alpha = E.addSlider(ui, { label: "α", min: 0.05, max: 1, step: 0.05, value: 0.1, onChange: redraw });
    redraw();
  }

  /* ═══════════ 11. LQR 倒立摆 ═══════════ */
  function simLQR(mount) {
    const ui = E.buildPanel(mount, {
      title: "LQR 倒立摆（小车摆）实时仿真",
      desc: "非线性小车摆模型，控制器为对上平衡点线性化后求解 Riccati 方程得到的 LQR 增益 u=-Kx。调整权重看『状态收敛速度 vs 控制量大小』的权衡，还可以随时施加干扰。"
    });
    const cv = E.addCanvas(ui, 560, 260);
    const cvT = E.addCanvas(ui, 560, 200);
    const s = {};
    const M = 1.0, mp = 0.2, l = 0.5, g = 9.81;
    let K = null, x = [0, 0, 0.25, 0], t = 0, hist = { t: [], th: [], pos: [], u: [] };
    let anim = null;
    const computeK = () => {
      const A = [[0, 1, 0, 0], [0, 0, -mp * g / M, 0], [0, 0, 0, 1], [0, 0, (M + mp) * g / (M * l), 0]];
      const B = [[0], [1 / M], [0], [-1 / (M * l)]];
      const Q = [[s.qx.get(), 0, 0, 0], [0, 1, 0, 0], [0, 0, s.qth.get(), 0], [0, 0, 0, 1]];
      const R = [[s.r.get()]];
      const res = num.lqrGain(A, B, Q, R);
      K = res.K[0];
      ui.readout.textContent = `K = [${K.map(v => v.toFixed(2)).join(", ")}]  （由求解 CARE 得到）\n点击画布或按钮给摆一脚 👟`;
    };
    const dyn = (x, u) => {
      const [pos, v, th, w] = x;
      const sinT = Math.sin(th), cosT = Math.cos(th);
      const den = M + mp * sinT * sinT;
      // th 从竖直向上量起
      const xdd = (u + mp * sinT * (l * w * w - g * cosT)) / den;
      const thdd = (-u * cosT - mp * l * w * w * sinT * cosT + (M + mp) * g * sinT) / (l * den);
      return [v, xdd, w, thdd];
    };
    const kick = () => { x[3] += (Math.random() > 0.5 ? 1 : -1) * (1.2 + Math.random()); };
    const reset = () => { x = [0, 0, 0.25, 0]; t = 0; hist = { t: [], th: [], pos: [], u: [] }; };
    const start = () => {
      computeK();
      if (anim) anim.alive = false;
      anim = E.loop(() => {
        let u = 0;
        for (let i = 0; i < 8; i++) {
          u = -(K[0] * x[0] + K[1] * x[1] + K[2] * x[2] + K[3] * x[3]);
          u = Math.max(-25, Math.min(25, u));
          x = num.rk4((x) => dyn(x, u), x, t, 0.002);
          t += 0.002;
        }
        if (Math.abs(x[2]) > 1.4) { reset(); return; } // 摔倒重来
        hist.t.push(t); hist.th.push(x[2]); hist.pos.push(x[0]); hist.u.push(u / 10);
        if (hist.t.length > 600) { hist.t.shift(); hist.th.shift(); hist.pos.shift(); hist.u.shift(); }
        // ── 画小车摆 ──
        const c = cv.getContext("2d");
        const W = 560, H = 260;
        c.save();
        c.setTransform(window.devicePixelRatio || 1, 0, 0, window.devicePixelRatio || 1, 0, 0);
        c.clearRect(0, 0, W, H);
        c.fillStyle = "#fcfdfe"; c.fillRect(0, 0, W, H);
        const scale = 90, cx0 = W / 2, groundY = 200;
        c.strokeStyle = "#b9c4d3"; c.beginPath(); c.moveTo(0, groundY); c.lineTo(W, groundY); c.stroke();
        // 刻度
        c.fillStyle = "#8592a5"; c.font = "9px Consolas";
        for (let m = -3; m <= 3; m++) {
          const px = cx0 + m * scale;
          c.fillRect(px, groundY, 1, 5);
          c.textAlign = "center"; c.fillText(m + "m", px, groundY + 16);
        }
        const cartX = cx0 + x[0] * scale;
        // 小车
        c.fillStyle = "#5b8def";
        c.fillRect(cartX - 30, groundY - 24, 60, 24);
        c.fillStyle = "#55617a";
        c.beginPath(); c.arc(cartX - 18, groundY, 7, 0, 7); c.fill();
        c.beginPath(); c.arc(cartX + 18, groundY, 7, 0, 7); c.fill();
        // 摆杆（θ 从竖直向上）
        const tipX = cartX + Math.sin(x[2]) * l * 2 * scale;
        const tipY = groundY - 24 - Math.cos(x[2]) * l * 2 * scale;
        c.strokeStyle = "#e8a838"; c.lineWidth = 5; c.lineCap = "round";
        c.beginPath(); c.moveTo(cartX, groundY - 24); c.lineTo(tipX, tipY); c.stroke();
        c.fillStyle = "#ef5b7b"; c.beginPath(); c.arc(tipX, tipY, 9, 0, 7); c.fill();
        // 控制力箭头
        if (Math.abs(u) > 0.3) {
          c.strokeStyle = "#7bc74d"; c.lineWidth = 3;
          const dir = Math.sign(u), len = Math.min(50, Math.abs(u) * 2.2);
          c.beginPath(); c.moveTo(cartX - dir * 34, groundY - 12); c.lineTo(cartX - dir * (34 + len), groundY - 12); c.stroke();
          c.beginPath(); c.moveTo(cartX - dir * 34, groundY - 12);
          c.lineTo(cartX - dir * 42, groundY - 17); c.lineTo(cartX - dir * 42, groundY - 7); c.closePath();
          c.fillStyle = "#7bc74d"; c.fill();
        }
        c.fillStyle = "#66738a"; c.font = "11px 'Segoe UI'"; c.textAlign = "left";
        c.fillText(`t=${t.toFixed(1)}s  θ=${(x[2] * 180 / Math.PI).toFixed(1)}°  x=${x[0].toFixed(2)}m  u=${u.toFixed(1)}N`, 10, 18);
        c.restore();
        // 曲线
        const pt = new E.Plot(cvT, { xmin: hist.t[0] || 0, xmax: Math.max(6, t), xlabel: "t (s)", title: "θ(rad)·蓝  位置(m)·青  u/10(N)·黄" });
        pt.autoRangeY([hist.th, hist.pos, hist.u, [-0.4, 0.4]]);
        pt.clear();
        pt.line(hist.t, hist.th, "#5b8def", 1.8);
        pt.line(hist.t, hist.pos, "#31b5a9", 1.6);
        pt.line(hist.t, hist.u, "#e8a838", 1.2);
      });
    };
    s.qx = E.addSlider(ui, { label: "位置权重 qₓ", min: 0.1, max: 100, step: 0.1, value: 10, onChange: () => computeK() });
    s.qth = E.addSlider(ui, { label: "角度权重 q_θ", min: 1, max: 500, step: 1, value: 100, fmt: v => (+v).toFixed(0), onChange: () => computeK() });
    s.r = E.addSlider(ui, { label: "控制权重 R", min: 0.01, max: 5, step: 0.01, value: 0.1, onChange: () => computeK() });
    E.addButton(ui, "👟 踢一脚", kick);
    E.addButton(ui, "重置", () => reset(), true);
    cv.addEventListener("click", kick);
    start();
  }

  /* ═══════════ 12. Kalman 滤波 ═══════════ */
  function simKalman(mount) {
    const ui = E.buildPanel(mount, {
      title: "Kalman 滤波：噪声中追踪运动目标",
      desc: "匀速运动模型 + 随机加速度扰动，只测得含噪位置。观察最优估计如何在『相信模型』与『相信测量』之间自动权衡（由 Q/R 之比决定）。"
    });
    const cv = E.addCanvas(ui, 540, 300);
    const cvV = E.addCanvas(ui, 540, 220);
    const s = {};
    const redraw = () => {
      const sr = s.r.get(), sq = s.q.get();
      const dt = 0.1, Nk = 120;
      const A = [[1, dt], [0, 1]], H = [[1, 0]];
      const Q = [[dt * dt * dt / 3 * sq * sq, dt * dt / 2 * sq * sq], [dt * dt / 2 * sq * sq, dt * sq * sq]];
      const R = [[sr * sr]];
      const rng = mulberry32(7);
      // 真值与测量
      let xt = [0, 1];
      const truth = [], meas = [];
      for (let k = 0; k < Nk; k++) {
        truth.push(xt.slice());
        meas.push(xt[0] + sr * gauss(rng));
        const acc = sq * gauss(rng);
        xt = [xt[0] + xt[1] * dt + 0.5 * acc * dt * dt, xt[1] + acc * dt + 0.3 * Math.sin(0.15 * k) * dt];
      }
      // KF
      let xh = [0, 0], P = [[10, 0], [0, 10]];
      const est = [], sig = [], vest = [];
      for (let k = 0; k < Nk; k++) {
        // 预测
        xh = num.mulv(A, xh);
        P = num.add(num.mul(num.mul(A, P), num.T(A)), Q);
        // 更新
        const S = P[0][0] + R[0][0];
        const Kk = [P[0][0] / S, P[1][0] / S];
        const innov = meas[k] - xh[0];
        xh = [xh[0] + Kk[0] * innov, xh[1] + Kk[1] * innov];
        const IKH = [[1 - Kk[0], 0], [-Kk[1], 1]];
        P = num.mul(IKH, P);
        est.push(xh[0]); sig.push(Math.sqrt(Math.max(P[0][0], 0))); vest.push(xh[1]);
      }
      const ts = truth.map((_, k) => k * dt);
      const p = new E.Plot(cv, { xmin: 0, xmax: Nk * dt, xlabel: "t (s)", ylabel: "位置", title: "位置：真值 vs 测量 vs 估计（±2σ 带）" });
      p.autoRangeY([truth.map(x => x[0]), meas, est]);
      p.clear();
      // ±2σ带
      const ctx = p.ctx;
      ctx.fillStyle = "rgba(49,181,169,0.15)";
      ctx.beginPath();
      est.forEach((v, k) => { const px = p.x2px(ts[k]), py = p.y2px(v + 2 * sig[k]); k === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py); });
      for (let k = est.length - 1; k >= 0; k--) ctx.lineTo(p.x2px(ts[k]), p.y2px(est[k] - 2 * sig[k]));
      ctx.closePath(); ctx.fill();
      meas.forEach((m, k) => p.dot(ts[k], m, "rgba(139,152,171,0.55)", 1.8));
      p.line(ts, truth.map(x => x[0]), "#5b8def", 2);
      p.line(ts, est, "#31b5a9", 2);
      p.legend([{ label: "真值", color: "#5b8def" }, { label: "KF 估计", color: "#31b5a9" }, { label: "测量点", color: "#8b98ab" }]);
      const pv = new E.Plot(cvV, { xmin: 0, xmax: Nk * dt, xlabel: "t (s)", ylabel: "速度", title: "速度：从未直接测量，却被估计出来了" });
      pv.autoRangeY([truth.map(x => x[1]), vest]);
      pv.clear();
      pv.line(ts, truth.map(x => x[1]), "#5b8def", 2);
      pv.line(ts, vest, "#31b5a9", 2);
      const rmseM = Math.sqrt(meas.reduce((s2, m, k) => s2 + (m - truth[k][0]) ** 2, 0) / Nk);
      const rmseE = Math.sqrt(est.reduce((s2, e, k) => s2 + (e - truth[k][0]) ** 2, 0) / Nk);
      ui.readout.textContent = `位置 RMSE — 原始测量: ${rmseM.toFixed(3)}   KF 估计: ${rmseE.toFixed(3)}   （提升 ${(100 * (1 - rmseE / rmseM)).toFixed(0)}%）\nR 大→更信模型（平滑但迟钝）；Q 大→更信测量（灵敏但吵）`;
    };
    s.r = E.addSlider(ui, { label: "测量噪声 σᵣ", min: 0.05, max: 2, step: 0.05, value: 0.5, onChange: redraw });
    s.q = E.addSlider(ui, { label: "过程噪声 σ_q", min: 0.05, max: 2, step: 0.05, value: 0.3, onChange: redraw });
    redraw();
  }

  /* ═══════════ 13. MPC ═══════════ */
  function simMPC(mount) {
    const ui = E.buildPanel(mount, {
      title: "MPC vs 饱和 LQR：约束下的较量",
      desc: "双积分器从 x₀=[8,0] 回原点，输入受限 |u|≤u_max。MPC 在优化中『知道』约束，提前减速；LQR 事后裁剪（饱和），可能超调。灰色细线是 MPC 每步的预测轨迹。"
    });
    const cv = E.addCanvas(ui, 540, 280);
    const cvU = E.addCanvas(ui, 540, 180);
    const s = {};
    const redraw = () => {
      const umax = s.umax.get(), N = Math.round(s.N.get()), qp = s.q.get();
      const dt = 0.15;
      const Ad = [[1, dt], [0, 1]], Bd = [[0.5 * dt * dt], [dt]];
      const Qd = [[qp, 0], [0, 1]], Rd = [[1]];
      // 终端权重 = DARE 近似解（用 CARE 连续近似即可教学演示，這裡用迭代 DARE）
      let Pt = num.clone(Qd);
      for (let i = 0; i < 500; i++) {
        const BtP = num.mul(num.T(Bd), Pt);
        const Sinv = 1 / (num.mul(BtP, Bd)[0][0] + Rd[0][0]);
        const Kt = num.scale(num.mul(BtP, Ad), Sinv);
        const Acl = num.sub(Ad, num.mul(Bd, Kt));
        const Pn = num.add(num.add(num.mul(num.mul(num.T(Acl), Pt), Acl),
          num.mul(num.mul(num.T(Kt), Rd), Kt)), Qd);
        Pt = num.scale(num.add(Pn, num.T(Pn)), 0.5);
      }
      const Klqr = (() => {
        const BtP = num.mul(num.T(Bd), Pt);
        const Sinv = 1 / (num.mul(BtP, Bd)[0][0] + Rd[0][0]);
        return num.scale(num.mul(BtP, Ad), Sinv)[0];
      })();
      // 预测矩阵
      const Phi = [], Gam = num.zeros(2 * N, N);
      let Apow = num.eye(2);
      const Apows = [num.eye(2)];
      for (let i = 1; i <= N; i++) { Apow = num.mul(Ad, Apow); Apows.push(Apow); }
      for (let i = 0; i < N; i++) {
        Phi.push(Apows[i + 1]);
        for (let j = 0; j <= i; j++) {
          const blk = num.mul(Apows[i - j], Bd);
          Gam[2 * i][j] = blk[0][0]; Gam[2 * i + 1][j] = blk[1][0];
        }
      }
      // 权重块对角
      const Qbar = num.zeros(2 * N, 2 * N);
      for (let i = 0; i < N; i++) {
        const Qi = (i === N - 1) ? Pt : Qd;
        Qbar[2 * i][2 * i] = Qi[0][0]; Qbar[2 * i][2 * i + 1] = Qi[0][1];
        Qbar[2 * i + 1][2 * i] = Qi[1][0]; Qbar[2 * i + 1][2 * i + 1] = Qi[1][1];
      }
      const GtQ = num.mul(num.T(Gam), Qbar);
      const Hmat = num.add(num.mul(GtQ, Gam), num.scale(num.eye(N), Rd[0][0]));
      // Lipschitz 估计（幂迭代）
      let v = new Array(N).fill(1);
      for (let i = 0; i < 30; i++) { v = num.mulv(Hmat, v); const nv = num.norm(v) || 1; v = num.vscale(v, 1 / nv); }
      const Lip = num.norm(num.mulv(Hmat, v));
      const solveMPC = x0 => {
        // min ½U'HU + f'U, f = Γ'Q̄Φx0
        const phiX = [];
        for (let i = 0; i < N; i++) { const xi = num.mulv(Phi[i], x0); phiX.push(xi[0], xi[1]); }
        const f = num.mulv(GtQ, phiX);
        let U = new Array(N).fill(0);
        const step = 1 / Lip;
        for (let it = 0; it < 120; it++) {
          const grad = num.vadd(num.mulv(Hmat, U), f);
          U = U.map((u, i) => Math.max(-umax, Math.min(umax, u - step * grad[i])));
        }
        return U;
      };
      // 仿真两种控制器
      const Tsteps = 90;
      const runs = {};
      // MPC
      let x = [8, 0];
      const mpcX = [], mpcV = [], mpcU = [], preds = [];
      for (let k = 0; k < Tsteps; k++) {
        const U = solveMPC(x);
        if (k % 6 === 0) {
          // 存预测轨迹
          let xp = x.slice(); const pr = [[k * dt, xp[0]]];
          for (let i = 0; i < N; i++) {
            xp = [xp[0] + xp[1] * dt + 0.5 * dt * dt * U[i], xp[1] + dt * U[i]];
            pr.push([(k + i + 1) * dt, xp[0]]);
          }
          preds.push(pr);
        }
        mpcX.push(x[0]); mpcV.push(x[1]); mpcU.push(U[0]);
        x = [x[0] + x[1] * dt + 0.5 * dt * dt * U[0], x[1] + dt * U[0]];
      }
      // 饱和 LQR
      x = [8, 0];
      const lqrX = [], lqrU = [];
      for (let k = 0; k < Tsteps; k++) {
        let u = -(Klqr[0] * x[0] + Klqr[1] * x[1]);
        u = Math.max(-umax, Math.min(umax, u));
        lqrX.push(x[0]); lqrU.push(u);
        x = [x[0] + x[1] * dt + 0.5 * dt * dt * u, x[1] + dt * u];
      }
      const ts = mpcX.map((_, k) => k * dt);
      const p = new E.Plot(cv, { xmin: 0, xmax: Tsteps * dt, xlabel: "t (s)", ylabel: "位置", title: "位置轨迹（灰细线 = MPC 预测）" });
      p.autoRangeY([mpcX, lqrX, [-1, 8.5]]);
      p.clear();
      preds.forEach(pr => p.line(pr.map(q => q[0]), pr.map(q => q[1]), "rgba(139,152,171,0.35)", 1));
      p.line([0, Tsteps * dt], [0, 0], "#3b4c66", 1, [4, 4]);
      p.line(ts, lqrX, "#e8a838", 1.8, [6, 3]);
      p.line(ts, mpcX, "#5b8def", 2.2);
      p.legend([{ label: "MPC", color: "#5b8def" }, { label: "饱和LQR", color: "#e8a838", dash: [6, 3] }]);
      const pu = new E.Plot(cvU, { xmin: 0, xmax: Tsteps * dt, ymin: -umax * 1.4, ymax: umax * 1.4, xlabel: "t (s)", ylabel: "u", title: "控制量（虚线=约束）" });
      pu.clear();
      pu.line([0, Tsteps * dt], [umax, umax], "#ef5b7b", 1, [4, 4]);
      pu.line([0, Tsteps * dt], [-umax, -umax], "#ef5b7b", 1, [4, 4]);
      pu.line(ts, lqrU, "#e8a838", 1.6, [6, 3]);
      pu.line(ts, mpcU, "#5b8def", 1.8);
      const osM = -Math.min(...mpcX), osL = -Math.min(...lqrX);
      ui.readout.textContent = `超调（冲过原点距离）— MPC: ${Math.max(0, osM).toFixed(2)}   饱和LQR: ${Math.max(0, osL).toFixed(2)}\n把 u_max 调小（约束更紧），差距更明显；视野 N 太短时 MPC 也会『近视』`;
    };
    s.umax = E.addSlider(ui, { label: "u_max", min: 0.3, max: 3, step: 0.1, value: 1, onChange: redraw });
    s.N = E.addSlider(ui, { label: "预测时域 N", min: 3, max: 40, step: 1, value: 20, fmt: v => (+v).toFixed(0), onChange: redraw });
    s.q = E.addSlider(ui, { label: "位置权重 q", min: 0.5, max: 20, step: 0.5, value: 4, onChange: redraw });
    redraw();
  }

  /* ═══════════ 14. 滑模控制 ═══════════ */
  function simSliding(mount) {
    const ui = E.buildPanel(mount, {
      title: "滑模控制实验",
      desc: "ẍ = d(t) + u，扰动 d=0.8sin(3t) 未知（只知 |d|≤0.8）。滑模面 s=cx+ẋ。到达滑模面后，动态退化为 ẋ=-cx，对扰动完全免疫。边界层 φ>0 可消除抖振（代价是精度）。"
    });
    const cv = E.addCanvas(ui, 400, 340);
    const cvT = E.addCanvas(ui, 440, 340);
    const s = {};
    const redraw = () => {
      const c0 = s.c.get(), k = s.k.get(), phi = s.phi.get();
      const T = 10, h = 0.0005;
      let x = [2.5, 1.5];
      const px = [], pv = [], ts = [], xs = [], us = [];
      let cnt = 0;
      for (let t = 0; t <= T; t += h) {
        const d = 0.8 * Math.sin(3 * t);
        const sv = c0 * x[0] + x[1];
        const sw = phi > 0.001 ? Math.max(-1, Math.min(1, sv / phi)) : Math.sign(sv);
        const u = -c0 * x[1] - k * sw;
        if (cnt++ % 40 === 0) { px.push(x[0]); pv.push(x[1]); ts.push(t); xs.push(x[0]); us.push(u); }
        // Euler（保留抖振细节）
        const xdd = d + u;
        x = [x[0] + h * x[1], x[1] + h * xdd];
      }
      const p = new E.Plot(cv, { xmin: -3, xmax: 3, ymin: -4, ymax: 3, xlabel: "x", ylabel: "ẋ", title: "相平面：先打到面上，再滑回原点" });
      p.clear();
      // 滑模面
      p.line([-3, 3], [3 * c0, -3 * c0], "#ef5b7b", 1.6, [6, 4]);
      p.text(1.8, -1.8 * c0 + 0.3, "s=0", "#ef5b7b");
      p.line(px, pv, "#5b8def", 1.8);
      p.dot(px[0], pv[0], "#e8a838", 4);
      const pt = new E.Plot(cvT, { xmin: 0, xmax: T, xlabel: "t (s)", title: "x(t)·蓝  u(t)·黄（注意抖振）" });
      pt.autoRangeY([xs, us]);
      pt.clear();
      pt.line(ts, xs, "#5b8def", 2);
      pt.line(ts, us, "#e8a838", 0.9);
      ui.readout.textContent = `k 必须 > |d|max=0.8 才能保证到达滑模面（到达条件 sṡ ≤ -η|s|）\nφ=0 时 u 高频切换（抖振）；φ>0 平滑但只收敛到边界层内`;
    };
    s.c = E.addSlider(ui, { label: "面斜率 c", min: 0.3, max: 4, step: 0.1, value: 1, onChange: redraw });
    s.k = E.addSlider(ui, { label: "切换增益 k", min: 0.2, max: 6, step: 0.1, value: 2, onChange: redraw });
    s.phi = E.addSlider(ui, { label: "边界层 φ", min: 0, max: 0.5, step: 0.01, value: 0, onChange: redraw });
    redraw();
  }

  /* ═══════════ 15. 反步法 ═══════════ */
  function simBackstepping(mount) {
    const ui = E.buildPanel(mount, {
      title: "反步法镇定严格反馈系统",
      desc: "ẋ₁=x₁²+x₂（x₂ 是『虚拟输入』），ẋ₂=u。反步法先设计 x₂ 的期望值 α(x₁)=-x₁²-k₁x₁，再让 x₂ 追上 α。右图验证复合 Lyapunov 函数 V=½z₁²+½z₂² 单调下降。"
    });
    const cv = E.addCanvas(ui, 440, 320);
    const cvV = E.addCanvas(ui, 400, 320);
    const s = {};
    const redraw = () => {
      const k1 = s.k1.get(), k2 = s.k2.get(), x10 = s.x0.get();
      const T = 6, h = 0.001;
      let x = [x10, 0];
      const ts = [], x1s = [], x2s = [], us = [], Vs = [];
      let cnt = 0;
      let diverged = false;
      for (let t = 0; t <= T; t += h) {
        const [x1, x2] = x;
        const alpha = -x1 * x1 - k1 * x1;
        const z1 = x1, z2 = x2 - alpha;
        const dalpha = (-2 * x1 - k1) * (x1 * x1 + x2);
        const u = -z1 - k2 * z2 + dalpha;
        const V = 0.5 * z1 * z1 + 0.5 * z2 * z2;
        if (cnt++ % 20 === 0) { ts.push(t); x1s.push(x1); x2s.push(x2); us.push(u); Vs.push(V); }
        x = num.rk4((x) => [x[0] * x[0] + x[1], -x[0] - k2 * (x[1] + x[0] * x[0] + k1 * x[0]) + (-2 * x[0] - k1) * (x[0] * x[0] + x[1])], x, t, h);
        if (Math.abs(x[0]) > 1e4) { diverged = true; break; }
      }
      const p = new E.Plot(cv, { xmin: 0, xmax: T, xlabel: "t (s)", title: "状态与控制量" });
      p.autoRangeY([x1s, x2s, us]);
      p.clear();
      p.line(ts, x1s, "#5b8def", 2);
      p.line(ts, x2s, "#31b5a9", 1.8);
      p.line(ts, us, "#e8a838", 1.3, [5, 3]);
      p.legend([{ label: "x₁", color: "#5b8def" }, { label: "x₂", color: "#31b5a9" }, { label: "u", color: "#e8a838", dash: [5, 3] }]);
      const pv = new E.Plot(cvV, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "V", title: "复合 Lyapunov 函数 V(t)" });
      pv.autoRangeY([Vs, [0]]);
      pv.clear();
      pv.line(ts, Vs, "#9d6bff", 2.2);
      ui.readout.textContent = diverged ? "数值发散（初值太大时试试增大 k₁,k₂）" :
        `V̇ = -k₁z₁² - k₂z₂² ≤ 0，理论保证全局渐近稳定\n注意 x₁² 是不稳定非线性（开环 x₁ 会有限时间爆炸），反步法把它精确抵消`;
    };
    s.k1 = E.addSlider(ui, { label: "k₁", min: 0.2, max: 8, step: 0.1, value: 2, onChange: redraw });
    s.k2 = E.addSlider(ui, { label: "k₂", min: 0.2, max: 8, step: 0.1, value: 2, onChange: redraw });
    s.x0 = E.addSlider(ui, { label: "初值 x₁(0)", min: -3, max: 3, step: 0.1, value: 1.5, onChange: redraw });
    redraw();
  }

  /* ═══════════ 16. 自适应控制 MRAC ═══════════ */
  function simAdaptive(mount) {
    const ui = E.buildPanel(mount, {
      title: "MRAC：参数未知也能追上参考模型",
      desc: "被控对象 ẋ = a·x + u，a 未知（甚至不稳定）。参考模型 ẋₘ = -2xₘ + 2r。自适应律 θ̇ = -γex 在线调整反馈增益，使 x → xₘ。观察 θ 收敛到理想值 θ* = -(a+2)。"
    });
    const cv = E.addCanvas(ui, 500, 280);
    const cvP = E.addCanvas(ui, 500, 200);
    const s = {};
    const redraw = () => {
      const gam = s.gamma.get(), a = s.a.get();
      const am = 2, bm = 2;
      const T = 30, h = 0.001;
      let xp = 0, xm = 0, th = 0;
      const ts = [], xps = [], xms = [], ths = [];
      let cnt = 0;
      for (let t = 0; t <= T; t += h) {
        const r = Math.sign(Math.sin(0.5 * t)) * 1; // 方波
        const u = th * xp + bm * r;
        const e = xp - xm;
        if (cnt++ % 30 === 0) { ts.push(t); xps.push(xp); xms.push(xm); ths.push(th); }
        // Euler足够
        const dxp = a * xp + u;
        const dxm = -am * xm + bm * r;
        const dth = -gam * e * xp;
        xp += h * dxp; xm += h * dxm; th += h * dth;
        if (!isFinite(xp)) break;
      }
      const p = new E.Plot(cv, { xmin: 0, xmax: T, xlabel: "t (s)", title: "输出跟踪：x（蓝）追 xₘ（红虚线）" });
      p.autoRangeY([xps, xms]);
      p.clear();
      p.line(ts, xms, "#ef5b7b", 1.8, [6, 4]);
      p.line(ts, xps, "#5b8def", 1.8);
      const pp = new E.Plot(cvP, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "θ", title: "自适应参数 θ(t) → θ*" });
      const thStar = -(a + am);
      pp.autoRangeY([ths, [thStar]]);
      pp.clear();
      pp.line([0, T], [thStar, thStar], "#3b4c66", 1.2, [4, 4]);
      pp.line(ts, ths, "#31b5a9", 2);
      pp.text(T * 0.75, thStar, `θ* = ${thStar.toFixed(1)}`, "#8b98ab");
      ui.readout.textContent = `真实 a = ${a.toFixed(1)}${a > 0 ? "（开环不稳定！）" : ""}   理想增益 θ* = -(a+aₘ) = ${thStar.toFixed(1)}\nγ 太小学得慢，γ 太大瞬态振荡——参数收敛需要激励足够丰富（PE 条件）`;
    };
    s.gamma = E.addSlider(ui, { label: "自适应速率 γ", min: 0.1, max: 20, step: 0.1, value: 2, onChange: redraw });
    s.a = E.addSlider(ui, { label: "未知参数 a", min: -2, max: 3, step: 0.1, value: 1, onChange: redraw });
    redraw();
  }

  /* ═══════════ 17. 反馈线性化 ═══════════ */
  function simFL(mount) {
    const ui = E.buildPanel(mount, {
      title: "反馈线性化 vs 线性控制：大角度摆跟踪",
      desc: "单摆 θ̈ = -10sinθ - 0.5θ̇ + u 跟踪 θd = A·sin(t)。反馈线性化用 u = 10sinθ + 0.5θ̇ + v 精确抵消非线性；线性 PD 只在小角度近似有效。加大幅值 A 看差距。"
    });
    const cv = E.addCanvas(ui, 540, 300);
    const s = {};
    const redraw = () => {
      const A0 = s.A.get(), kp = s.kp.get(), kd = 2 * Math.sqrt(kp) * 0.9;
      const T = 12, h = 0.002;
      const run = useFL => {
        let x = [0.5, 0];
        const ts = [], th = [];
        let cnt = 0;
        for (let t = 0; t <= T; t += h) {
          const thd = A0 * Math.sin(t), thdD = A0 * Math.cos(t), thdDD = -A0 * Math.sin(t);
          const e = x[0] - thd, ed = x[1] - thdD;
          const v = thdDD - kd * ed - kp * e;
          const u = useFL ? 10 * Math.sin(x[0]) + 0.5 * x[1] + v : v;
          if (cnt++ % 10 === 0) { ts.push(t); th.push(x[0]); }
          x = num.rk4((x) => [x[1], -10 * Math.sin(x[0]) - 0.5 * x[1] + u], x, t, h);
        }
        return { ts, th };
      };
      const rFL = run(true), rPD = run(false);
      const thd = rFL.ts.map(t => A0 * Math.sin(t));
      const p = new E.Plot(cv, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "θ (rad)", title: "跟踪对比" });
      p.autoRangeY([thd, rFL.th, rPD.th]);
      p.clear();
      p.line(rFL.ts, thd, "#3b4c66", 1.4, [5, 4]);
      p.line(rPD.ts, rPD.th, "#e8a838", 1.8);
      p.line(rFL.ts, rFL.th, "#31b5a9", 2);
      p.legend([{ label: "目标 θd", color: "#3b4c66", dash: [5, 4] }, { label: "仅线性 PD", color: "#e8a838" }, { label: "反馈线性化", color: "#31b5a9" }]);
      const err = (r) => Math.sqrt(r.th.reduce((s2, v, i) => s2 + (v - A0 * Math.sin(r.ts[i])) ** 2, 0) / r.th.length);
      ui.readout.textContent = `跟踪 RMSE — 反馈线性化: ${err(rFL).toFixed(4)}   仅线性 PD: ${err(rPD).toFixed(4)}\nA 小时两者接近（线性化近似成立）；A 大时非线性项 10sinθ 偏离 10θ，PD 出现系统性误差`;
    };
    s.A = E.addSlider(ui, { label: "目标幅值 A", min: 0.2, max: 2.5, step: 0.1, value: 1.5, onChange: redraw });
    s.kp = E.addSlider(ui, { label: "kp", min: 4, max: 100, step: 1, value: 25, fmt: v => (+v).toFixed(0), onChange: redraw });
    redraw();
  }

  /* ═══════════ 18. EKF 机器人定位 ═══════════ */
  function simEKF(mount) {
    const ui = E.buildPanel(mount, {
      title: "EKF 机器人定位：里程计 + 路标",
      desc: "差速轮机器人转圈，里程计有噪声。灰色虚线 = 纯航位推算（误差累积发散）；青色 = EKF 融合路标测距测角后的估计。调大里程计噪声看 EKF 的『纠偏』能力。"
    });
    const cv = E.addCanvas(ui, 460, 380);
    const cvE = E.addCanvas(ui, 400, 380);
    const s = {};
    const LMS = [[3, 3], [-3, 2], [2, -3]];
    const dt = 0.1, STEPS = 620;
    const redraw = () => {
      const sv = s.odo.get(), sr = s.meas.get(), every = Math.round(s.every.get());
      const rng = mulberry32(11);
      const Q = [[0.02 * 0.02, 0, 0], [0, 0.02 * 0.02, 0], [0, 0, 0.01 * 0.01]];
      const Rm = [[sr * sr, 0], [0, 0.08 * 0.08]];
      const fstep = (x, v, w) => [x[0] + v * Math.cos(x[2]) * dt, x[1] + v * Math.sin(x[2]) * dt, x[2] + w * dt];
      let xt = [0, -2, 0], xd = [0, -2, 0], xe = [0, -2, 0];
      let P = [[0.1, 0, 0], [0, 0.1, 0], [0, 0, 0.1]];
      const tPath = [], dPath = [], ePath = [], errD = [], errE = [], ts = [];
      for (let k = 0; k < STEPS; k++) {
        const v = 1.0, w = 0.5;
        const vn = v + sv * gauss(rng), wn = w + 0.6 * sv * gauss(rng);
        xt = fstep(xt, v, w);
        xd = fstep(xd, vn, wn);
        // EKF 预测
        xe = fstep(xe, vn, wn);
        const F = [[1, 0, -vn * Math.sin(xe[2]) * dt], [0, 1, vn * Math.cos(xe[2]) * dt], [0, 0, 1]];
        P = num.add(num.mul(num.mul(F, P), num.T(F)), Q);
        // 观测更新
        if (k % every === 0) {
          const lm = LMS[Math.floor(k / every) % LMS.length];
          const dx = lm[0] - xt[0], dy = lm[1] - xt[1];
          const z = [Math.hypot(dx, dy) + sr * gauss(rng),
                     Math.atan2(dy, dx) - xt[2] + 0.08 * gauss(rng)];
          const dxh = lm[0] - xe[0], dyh = lm[1] - xe[1];
          const q = dxh * dxh + dyh * dyh, sq = Math.sqrt(q);
          const H = [[-dxh / sq, -dyh / sq, 0], [dyh / q, -dxh / q, -1]];
          const S2 = num.add(num.mul(num.mul(H, P), num.T(H)), Rm);
          const K = num.mul(num.mul(P, num.T(H)), num.inv(S2));
          let innov = [z[0] - sq, z[1] - (Math.atan2(dyh, dxh) - xe[2])];
          while (innov[1] > Math.PI) innov[1] -= 2 * Math.PI;
          while (innov[1] < -Math.PI) innov[1] += 2 * Math.PI;
          const corr = num.mulv(K, innov);
          xe = [xe[0] + corr[0], xe[1] + corr[1], xe[2] + corr[2]];
          P = num.mul(num.sub(num.eye(3), num.mul(K, H)), P);
        }
        tPath.push([xt[0], xt[1]]); dPath.push([xd[0], xd[1]]); ePath.push([xe[0], xe[1]]);
        ts.push(k * dt);
        errD.push(Math.hypot(xd[0] - xt[0], xd[1] - xt[1]));
        errE.push(Math.hypot(xe[0] - xt[0], xe[1] - xt[1]));
      }
      const p = new E.Plot(cv, { xmin: -4.5, xmax: 4.5, ymin: -4.5, ymax: 4.5, xlabel: "x (m)", ylabel: "y (m)", title: "轨迹（菱形=路标）" });
      p.clear();
      LMS.forEach(lm => { p.cross(lm[0], lm[1], "#d98e1f", 6); p.dot(lm[0], lm[1], "#d98e1f", 2.5); });
      p.line(dPath.map(q => q[0]), dPath.map(q => q[1]), "#9aa6b6", 1.6, [5, 4]);
      p.line(tPath.map(q => q[0]), tPath.map(q => q[1]), "#5b8def", 2.2);
      p.line(ePath.map(q => q[0]), ePath.map(q => q[1]), "#199a8c", 2);
      p.legend([{ label: "真值", color: "#5b8def" }, { label: "航位推算", color: "#9aa6b6", dash: [5, 4] }, { label: "EKF", color: "#199a8c" }]);
      const pe = new E.Plot(cvE, { xmin: 0, xmax: STEPS * dt, xlabel: "t (s)", ylabel: "位置误差 (m)", title: "误差对比" });
      pe.autoRangeY([errD, errE, [0]]);
      pe.clear();
      pe.line(ts, errD, "#9aa6b6", 1.8, [5, 4]);
      pe.line(ts, errE, "#199a8c", 2);
      ui.readout.textContent =
        `末端误差 — 航位推算: ${errD[errD.length - 1].toFixed(2)} m   EKF: ${errE[errE.length - 1].toFixed(2)} m\n` +
        `航位推算误差随时间累积；EKF 每次看到路标就把漂移拉回来（观测越频繁越准）`;
    };
    s.odo = E.addSlider(ui, { label: "里程计噪声", min: 0.01, max: 0.25, step: 0.01, value: 0.08, onChange: redraw });
    s.meas = E.addSlider(ui, { label: "测距噪声 σᵣ", min: 0.05, max: 1, step: 0.05, value: 0.25, onChange: redraw });
    s.every = E.addSlider(ui, { label: "观测间隔(步)", min: 2, max: 40, step: 1, value: 5, fmt: v => (+v).toFixed(0), onChange: redraw });
    redraw();
  }

  /* ═══════════ 19. 粒子滤波 ═══════════ */
  function simPF(mount) {
    const ui = E.buildPanel(mount, {
      title: "粒子滤波全局定位（实时动画）",
      desc: "机器人位置完全未知（粒子全图均匀撒开），只能测到两个路标的距离。看粒子云如何从『均匀 → 双峰（镜像歧义）→ 单峰收敛』——这种多峰信念是 EKF 的高斯椭圆无法表达的。"
    });
    const cv = E.addCanvas(ui, 480, 420);
    const s = {};
    const LMS = [[2, 2], [-2, -1]];
    const dt = 0.1;
    let anim = null;
    const start = () => {
      const N = Math.round(s.N.get()), sr = s.meas.get();
      const rng = mulberry32(23);
      let xt = [2.5, -2, 0.8];
      // 粒子: [x, y, th]
      let px = new Float64Array(N), py = new Float64Array(N), pth = new Float64Array(N);
      let w = new Float64Array(N).fill(1 / N);
      for (let i = 0; i < N; i++) {
        px[i] = -4.5 + 9 * rng(); py[i] = -4.5 + 9 * rng(); pth[i] = -Math.PI + 2 * Math.PI * rng();
      }
      let t = 0;
      const p = new E.Plot(cv, { xmin: -5, xmax: 5, ymin: -5, ymax: 5, xlabel: "x (m)", ylabel: "y (m)", title: "粒子云（青点）· 真值（蓝）· 估计（红叉）" });
      if (anim) anim.alive = false;
      anim = E.loop(() => {
        for (let rep = 0; rep < 2; rep++) {
          const v = 1.0, om = 0.4;
          // 真值演化
          xt = [xt[0] + v * Math.cos(xt[2]) * dt, xt[1] + v * Math.sin(xt[2]) * dt, xt[2] + om * dt];
          // 粒子传播
          for (let i = 0; i < N; i++) {
            const vn = v + 0.1 * gauss(rng), on = om + 0.06 * gauss(rng);
            pth[i] += on * dt;
            px[i] += vn * Math.cos(pth[i]) * dt;
            py[i] += vn * Math.sin(pth[i]) * dt;
          }
          // 称重: 两个路标的距离
          let wsum = 0;
          for (const lm of LMS) {
            const z = Math.hypot(xt[0] - lm[0], xt[1] - lm[1]) + sr * gauss(rng);
            for (let i = 0; i < N; i++) {
              const d = Math.hypot(px[i] - lm[0], py[i] - lm[1]);
              const e = (z - d) / sr;
              w[i] *= Math.exp(-0.5 * e * e) + 1e-300;
            }
          }
          for (let i = 0; i < N; i++) wsum += w[i];
          let neff = 0;
          for (let i = 0; i < N; i++) { w[i] /= wsum; neff += w[i] * w[i]; }
          neff = 1 / neff;
          // 系统重采样
          if (neff < N / 2) {
            const nx = new Float64Array(N), ny = new Float64Array(N), nth = new Float64Array(N);
            const step = 1 / N; let u0 = rng() * step, cum = w[0], j = 0;
            for (let i = 0; i < N; i++) {
              const target = u0 + i * step;
              while (cum < target && j < N - 1) { j++; cum += w[j]; }
              nx[i] = px[j] + 0.03 * gauss(rng); ny[i] = py[j] + 0.03 * gauss(rng); nth[i] = pth[j] + 0.02 * gauss(rng);
            }
            px = nx; py = ny; pth = nth; w.fill(1 / N);
          }
          t += dt;
          ui._neff = neff;
        }
        // 估计
        let ex = 0, ey = 0;
        for (let i = 0; i < N; i++) { ex += w[i] * px[i]; ey += w[i] * py[i]; }
        // 绘制
        p.clear();
        LMS.forEach(lm => { p.cross(lm[0], lm[1], "#d98e1f", 6); });
        const c2 = p.ctx;
        c2.fillStyle = "rgba(25,154,140,0.35)";
        for (let i = 0; i < N; i++) {
          const X = p.x2px(px[i]), Y = p.y2px(py[i]);
          c2.fillRect(X - 1, Y - 1, 2, 2);
        }
        p.dot(xt[0], xt[1], "#3565c9", 6);
        p.line([xt[0], xt[0] + 0.5 * Math.cos(xt[2])], [xt[1], xt[1] + 0.5 * Math.sin(xt[2])], "#3565c9", 2);
        p.cross(ex, ey, "#d94f6e", 7);
        const err = Math.hypot(ex - xt[0], ey - xt[1]);
        ui.readout.textContent =
          `t = ${t.toFixed(1)} s   N_eff = ${(ui._neff || 0).toFixed(0)} / ${N}   估计误差 = ${err.toFixed(2)} m\n` +
          `前期出现两团粒子（两路标测距的镜像歧义），机器人转弯打破对称后收敛到单峰`;
        if (t > 45) anim.alive = false;
      });
    };
    s.N = E.addSlider(ui, { label: "粒子数 N", min: 200, max: 3000, step: 100, value: 1200, fmt: v => (+v).toFixed(0), onChange: start });
    s.meas = E.addSlider(ui, { label: "测距噪声 σᵣ", min: 0.1, max: 1, step: 0.05, value: 0.3, onChange: start });
    E.addButton(ui, "重新开始", start);
    start();
  }

  /* ═══════════ 20. Super-Twisting ═══════════ */
  function simHOSM(mount) {
    const ui = E.buildPanel(mount, {
      title: "Super-Twisting vs 一阶滑模",
      desc: "滑模动态 ṡ = ρ(t) + u，扰动 ρ = a·sin(3t) 未知。一阶滑模的 u 直接含 sign(s) → 高频抖振；Super-Twisting 把 sign 藏进 u̇，控制量连续且 s 精确归零（二阶滑模）。"
    });
    const cvS = E.addCanvas(ui, 460, 260);
    const cvU = E.addCanvas(ui, 460, 260);
    const s = {};
    const redraw = () => {
      const a = s.dist.get(), kfo = s.kfo.get();
      const L = 3 * a; // |dρ/dt| ≤ 3a
      const k1 = 1.5 * Math.sqrt(L), k2 = 1.1 * L;
      const T = 8, h = 1e-4;
      const run = mode => {
        let sv = 1.5, v = 0;
        const ts = [], ss = [], us = [];
        let tvU = 0, uPrev = null, errAcc = 0, errN = 0;
        for (let i = 0; i < T / h; i++) {
          const t = i * h;
          const rho = a * Math.sin(3 * t);
          let u;
          if (mode === "fo") u = -kfo * Math.sign(sv);
          else {
            u = -k1 * Math.sqrt(Math.abs(sv)) * Math.sign(sv) + v;
            v += h * (-k2 * Math.sign(sv));
          }
          if (uPrev !== null) tvU += Math.abs(u - uPrev);
          uPrev = u;
          sv += h * (rho + u);
          if (t > 4) { errAcc += Math.abs(sv); errN++; }
          if (i % 40 === 0) { ts.push(t); ss.push(sv); us.push(u); }
        }
        return { ts, ss, us, tvU, meanErr: errAcc / errN };
      };
      const fo = run("fo"), st = run("st");
      const pS = new E.Plot(cvS, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "s", title: "滑模变量 s(t)" });
      pS.autoRangeY([fo.ss, st.ss]);
      pS.clear();
      pS.line(fo.ts, fo.ss, "#d98e1f", 1.4);
      pS.line(st.ts, st.ss, "#199a8c", 2);
      pS.legend([{ label: "一阶滑模", color: "#d98e1f" }, { label: "Super-Twisting", color: "#199a8c" }]);
      const pU = new E.Plot(cvU, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "u", title: "控制量 u(t)（注意抖振 vs 连续）" });
      pU.autoRangeY([fo.us, st.us]);
      pU.clear();
      pU.line(fo.ts, fo.us, "#d98e1f", 0.8);
      pU.line(st.ts, st.us, "#199a8c", 2);
      ui.readout.textContent =
        `稳态 |s| 均值 — 一阶: ${fo.meanErr.toExponential(2)}   ST: ${st.meanErr.toExponential(2)}\n` +
        `u 总变差（抖振度量）— 一阶: ${fo.tvU.toFixed(0)}   ST: ${st.tvU.toFixed(1)}\n` +
        `ST 增益自动整定: k₁=1.5√L=${k1.toFixed(2)}, k₂=1.1L=${k2.toFixed(2)}（L=3a 为扰动导数界）`;
    };
    s.dist = E.addSlider(ui, { label: "扰动幅值 a", min: 0.2, max: 2, step: 0.1, value: 0.8, onChange: redraw });
    s.kfo = E.addSlider(ui, { label: "一阶增益 k", min: 0.5, max: 6, step: 0.1, value: 1.5, onChange: redraw });
    redraw();
  }

  /* ═══════════ 21. 时滞与 Smith 预估器 ═══════════ */
  function simDelay(mount) {
    const ui = E.buildPanel(mount, {
      title: "时滞的杀伤力与 Smith 预估器",
      desc: "对象 G(s) = e^{-τs}/(s+1)，同一组 PI 参数：橙色 = 直接闭环（延迟吃掉相位裕度，τ 大即振荡失稳）；青色 = Smith 预估器（把延迟『搬出』回路，恢复无延迟设计）。"
    });
    const cv = E.addCanvas(ui, 500, 300);
    const cvU = E.addCanvas(ui, 500, 180);
    const s = {};
    const redraw = () => {
      const tau = s.tau.get(), Kp = s.kp.get(), Ki = s.ki.get();
      const h = 0.005, T = 25;
      const nd = Math.max(1, Math.round(tau / h));
      const run = useSmith => {
        let y = 0, xi = 0;             // 对象输出, 积分器
        let ym = 0, ymd = 0;           // Smith 内部模型(无延迟/带延迟)
        const buf = new Array(nd).fill(0);      // 对象输入延迟队列
        const bufM = new Array(nd).fill(0);     // 模型延迟队列
        const ts = [], ys = [], us = [];
        let ptr = 0;
        for (let i = 0; i < T / h; i++) {
          const t = i * h;
          const r = 1;
          const fb = useSmith ? (y + ym - ymd) : y;
          const e = r - fb;
          xi += Ki * e * h;
          xi = Math.max(-3, Math.min(3, xi));
          let u = Kp * e + xi;
          u = Math.max(-5, Math.min(5, u));
          // 对象: y' = -y + u(t-τ)
          const uDel = buf[ptr];
          y += h * (-y + uDel);
          buf[ptr] = u;
          // 内部模型(与对象同构)
          const uDelM = bufM[ptr];
          ym += h * (-ym + u);
          ymd += h * (-ymd + uDelM);
          bufM[ptr] = u;
          ptr = (ptr + 1) % nd;
          if (i % 20 === 0) { ts.push(t); ys.push(y); us.push(u); }
          if (Math.abs(y) > 50) { ys[ys.length - 1] = NaN; break; }
        }
        return { ts, ys, us };
      };
      const direct = run(false), smith = run(true);
      const p = new E.Plot(cv, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "y", title: "阶跃响应对比（同一组 PI 参数）" });
      p.autoRangeY([direct.ys.filter(isFinite), smith.ys, [0, 1.6]]);
      p.opt.ymin = Math.max(p.opt.ymin, -2); p.opt.ymax = Math.min(p.opt.ymax, 4);
      p.clear();
      p.line([0, T], [1, 1], "#c9d3e0", 1, [4, 4]);
      p.line(direct.ts, direct.ys, "#d98e1f", 1.8);
      p.line(smith.ts, smith.ys, "#199a8c", 2.2);
      p.legend([{ label: "直接 PI", color: "#d98e1f" }, { label: "PI + Smith", color: "#199a8c" }]);
      const pu = new E.Plot(cvU, { xmin: 0, xmax: T, xlabel: "t (s)", ylabel: "u", title: "控制量" });
      pu.autoRangeY([direct.us, smith.us]);
      pu.clear();
      pu.line(direct.ts, direct.us, "#d98e1f", 1.2);
      pu.line(smith.ts, smith.us, "#199a8c", 1.6);
      // 相位裕度侵蚀估算: 无延迟回路 L=(Kp s+Ki)/(s(s+1)) 的 wgc
      const m = margins([Kp, Ki], [1, 1, 0]);
      const pmLoss = m.wgc ? (tau * m.wgc * 180 / Math.PI) : 0;
      ui.readout.textContent =
        `延迟吃掉的相位 = ω_gc·τ = ${pmLoss.toFixed(0)}°（无延迟设计的 PM = ${m.pm ? m.pm.toFixed(0) : "—"}°）\n` +
        `PM 被吃穿 ⇒ 直接闭环失稳；Smith 预估器用模型差 (ŷ无延迟 − ŷ带延迟) 补偿反馈，特征方程回到无延迟形态`;
    };
    s.tau = E.addSlider(ui, { label: "延迟 τ (s)", min: 0.1, max: 2.5, step: 0.05, value: 1.0, onChange: redraw });
    s.kp = E.addSlider(ui, { label: "Kp", min: 0.2, max: 6, step: 0.1, value: 2, onChange: redraw });
    s.ki = E.addSlider(ui, { label: "Ki", min: 0, max: 4, step: 0.1, value: 1, onChange: redraw });
    redraw();
  }

  /* ═══════════ 注册表 ═══════════ */
  const registry = {
    "sim-eig-phase": simEigPhase,
    "sim-secondorder": simSecondOrder,
    "sim-pid": simPID,
    "sim-poleplace": simPolePlace,
    "sim-observer": simObserver,
    "sim-lyapunov": simLyapunov,
    "sim-nyquist": simNyquist,
    "sim-rootlocus": simRootLocus,
    "sim-bode": simBode,
    "sim-leadlag": simLeadLag,
    "sim-lqr": simLQR,
    "sim-kalman": simKalman,
    "sim-mpc": simMPC,
    "sim-sliding": simSliding,
    "sim-backstepping": simBackstepping,
    "sim-adaptive": simAdaptive,
    "sim-fl": simFL,
    "sim-ekf": simEKF,
    "sim-pf": simPF,
    "sim-hosm": simHOSM,
    "sim-delay": simDelay,
  };

  function mount(id, el) {
    stopAll();
    try { registry[id] && registry[id](el); }
    catch (err) {
      el.innerHTML = `<div class="box note">仿真加载失败: ${err.message}</div>`;
      console.error(err);
    }
  }
  function stopAll() { E.stopAll(); }

  return { registry, mount, stopAll };
})();
