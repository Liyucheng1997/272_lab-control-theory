/* 仿真引擎：数值计算 + 画图 + 控件 UI + 仿真生命周期管理 */
window.SimEngine = (function () {

  /* ═════════ 数值计算 ═════════ */
  const num = {
    zeros: (r, c) => Array.from({ length: r }, () => new Array(c).fill(0)),
    eye: n => { const I = num.zeros(n, n); for (let i = 0; i < n; i++) I[i][i] = 1; return I; },
    clone: A => A.map(row => row.slice()),
    T: A => A[0].map((_, j) => A.map(row => row[j])),
    add: (A, B) => A.map((r, i) => r.map((v, j) => v + B[i][j])),
    sub: (A, B) => A.map((r, i) => r.map((v, j) => v - B[i][j])),
    scale: (A, s) => A.map(r => r.map(v => v * s)),
    mul: (A, B) => {
      const n = A.length, m = B[0].length, p = B.length;
      const C = num.zeros(n, m);
      for (let i = 0; i < n; i++) for (let k = 0; k < p; k++) {
        const a = A[i][k]; if (a === 0) continue;
        for (let j = 0; j < m; j++) C[i][j] += a * B[k][j];
      }
      return C;
    },
    mulv: (A, x) => A.map(row => row.reduce((s, v, j) => s + v * x[j], 0)),
    dot: (a, b) => a.reduce((s, v, i) => s + v * b[i], 0),
    vadd: (a, b) => a.map((v, i) => v + b[i]),
    vsub: (a, b) => a.map((v, i) => v - b[i]),
    vscale: (a, s) => a.map(v => v * s),
    norm: a => Math.sqrt(a.reduce((s, v) => s + v * v, 0)),
    // 高斯消元解 Ax=b（带部分主元）
    solve: (A0, b0) => {
      const n = A0.length;
      const A = A0.map((row, i) => row.concat([b0[i]]));
      for (let col = 0; col < n; col++) {
        let piv = col;
        for (let r = col + 1; r < n; r++) if (Math.abs(A[r][col]) > Math.abs(A[piv][col])) piv = r;
        [A[col], A[piv]] = [A[piv], A[col]];
        if (Math.abs(A[col][col]) < 1e-12) continue;
        for (let r = 0; r < n; r++) {
          if (r === col) continue;
          const f = A[r][col] / A[col][col];
          for (let c = col; c <= n; c++) A[r][c] -= f * A[col][c];
        }
      }
      return A.map((row, i) => row[n] / (Math.abs(row[i]) > 1e-12 ? row[i] : 1e-12));
    },
    inv: A => {
      const n = A.length, I = num.eye(n);
      const cols = [];
      for (let j = 0; j < n; j++) cols.push(num.solve(A, I.map(r => r[j])));
      return num.T(cols);
    },
    // RK4 一步：f(x, t) -> dx
    rk4: (f, x, t, h) => {
      const k1 = f(x, t);
      const k2 = f(num.vadd(x, num.vscale(k1, h / 2)), t + h / 2);
      const k3 = f(num.vadd(x, num.vscale(k2, h / 2)), t + h / 2);
      const k4 = f(num.vadd(x, num.vscale(k3, h)), t + h);
      return num.vadd(x, num.vscale(num.vadd(num.vadd(k1, num.vscale(k2, 2)), num.vadd(num.vscale(k3, 2), k4)), h / 6));
    },
    /* 求解连续代数 Riccati 方程 A'P+PA-PBR⁻¹B'P+Q=0
       方法：对微分 Riccati 方程 Ṗ = A'P+PA-PBR⁻¹B'P+Q 从 P=Q 积分至收敛 */
    care: (A, B, Q, R, opt = {}) => {
      const At = num.T(A), Bt = num.T(B), Rinv = num.inv(R);
      let P = num.clone(Q);
      const h = opt.h || 0.005, maxIter = opt.maxIter || 40000, tol = opt.tol || 1e-9;
      const f = P => {
        const PB = num.mul(P, B);
        return num.add(num.sub(num.add(num.mul(At, P), num.mul(P, A)),
          num.mul(num.mul(PB, Rinv), num.T(PB))), Q);
      };
      for (let i = 0; i < maxIter; i++) {
        const dP = f(P);
        let maxd = 0;
        for (const row of dP) for (const v of row) maxd = Math.max(maxd, Math.abs(v));
        if (maxd < tol) break;
        P = num.add(P, num.scale(dP, h));
        // 对称化，抑制数值漂移
        P = num.scale(num.add(P, num.T(P)), 0.5);
      }
      return P;
    },
    lqrGain: (A, B, Q, R) => {
      const P = num.care(A, B, Q, R);
      return { K: num.mul(num.mul(num.inv(R), num.T(B)), P), P };
    },
    /* Durand–Kerner 求多项式所有根，coeffs 从最高次到常数项 */
    polyRoots: coeffs => {
      let c = coeffs.slice();
      while (c.length > 1 && Math.abs(c[0]) < 1e-14) c.shift();
      const n = c.length - 1;
      if (n < 1) return [];
      const a0 = c[0];
      c = c.map(v => v / a0);
      // 复数运算
      const cmul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
      const csub = (a, b) => [a[0] - b[0], a[1] - b[1]];
      const cdiv = (a, b) => {
        const d = b[0] * b[0] + b[1] * b[1] || 1e-30;
        return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d];
      };
      const peval = z => {
        let r = [0, 0];
        for (const co of c) r = [r[0] * z[0] - r[1] * z[1] + co, r[0] * z[1] + r[1] * z[0]];
        return r;
      };
      let roots = [];
      for (let i = 0; i < n; i++) roots.push([Math.pow(0.4 + 0.9, 1) * Math.cos(2 * Math.PI * i / n + 0.5), 1.3 * Math.sin(2 * Math.PI * i / n + 0.5)]);
      for (let iter = 0; iter < 300; iter++) {
        let maxStep = 0;
        for (let i = 0; i < n; i++) {
          let denom = [1, 0];
          for (let j = 0; j < n; j++) if (j !== i) denom = cmul(denom, csub(roots[i], roots[j]));
          const delta = cdiv(peval(roots[i]), denom);
          roots[i] = csub(roots[i], delta);
          maxStep = Math.max(maxStep, Math.abs(delta[0]) + Math.abs(delta[1]));
        }
        if (maxStep < 1e-12) break;
      }
      return roots; // [[re, im], ...]
    },
    // 传递函数 G(s) = num/den 在 s=jw 处的值
    freqResp: (numC, denC, w) => {
      const evalP = (c, w) => {
        let re = 0, im = 0;
        for (const co of c) {
          const nre = re * 0 - im * w + co * 0 + 0; // placeholder
        }
        // Horner with s = jw
        let r = [0, 0];
        for (const co of c) r = [-r[1] * w + co, r[0] * w];
        return r;
      };
      const nv = evalP(numC, w), dv = evalP(denC, w);
      const d = dv[0] * dv[0] + dv[1] * dv[1] || 1e-30;
      return [(nv[0] * dv[0] + nv[1] * dv[1]) / d, (nv[1] * dv[0] - nv[0] * dv[1]) / d];
    },
    polyMul: (a, b) => {
      const r = new Array(a.length + b.length - 1).fill(0);
      for (let i = 0; i < a.length; i++) for (let j = 0; j < b.length; j++) r[i + j] += a[i] * b[j];
      return r;
    },
  };

  /* ═════════ 画图 ═════════ */
  class Plot {
    constructor(canvas, opt = {}) {
      this.cv = canvas; this.ctx = canvas.getContext("2d");
      this.opt = Object.assign({
        xlabel: "", ylabel: "", title: "",
        xmin: 0, xmax: 10, ymin: -1, ymax: 1,
        autoY: false, grid: true, equal: false,
        padding: { l: 46, r: 12, t: 26, b: 30 },
      }, opt);
      const dpr = window.devicePixelRatio || 1;
      const w = canvas.getAttribute("width"), h = canvas.getAttribute("height");
      canvas.style.width = w + "px"; canvas.style.height = h + "px";
      canvas.width = w * dpr; canvas.height = h * dpr;
      this.ctx.scale(dpr, dpr);
      this.W = +w; this.H = +h;
      this.series = [];
    }
    x2px(x) { const p = this.opt.padding; return p.l + (x - this.opt.xmin) / (this.opt.xmax - this.opt.xmin) * (this.W - p.l - p.r); }
    y2px(y) { const p = this.opt.padding; return this.H - p.b - (y - this.opt.ymin) / (this.opt.ymax - this.opt.ymin) * (this.H - p.t - p.b); }
    clear() {
      const c = this.ctx;
      c.clearRect(0, 0, this.W, this.H);
      const p = this.opt.padding;
      c.fillStyle = "#fcfdfe"; c.fillRect(0, 0, this.W, this.H);
      if (this.opt.title) {
        c.fillStyle = "#6a7686"; c.font = "12px 'Segoe UI', sans-serif"; c.textAlign = "center";
        c.fillText(this.opt.title, this.W / 2, 16);
      }
      if (this.opt.grid) this.drawAxes();
    }
    niceTicks(min, max, n = 5) {
      const span = max - min || 1;
      const step0 = span / n;
      const mag = Math.pow(10, Math.floor(Math.log10(step0)));
      const norm = step0 / mag;
      const step = (norm < 1.5 ? 1 : norm < 3.5 ? 2 : norm < 7.5 ? 5 : 10) * mag;
      const ticks = [];
      for (let v = Math.ceil(min / step) * step; v <= max + 1e-9; v += step) ticks.push(+v.toFixed(10));
      return ticks;
    }
    drawAxes() {
      const c = this.ctx, o = this.opt, p = o.padding;
      c.strokeStyle = "#eaeef5"; c.lineWidth = 1;
      c.fillStyle = "#93a0b2"; c.font = "10px Consolas, monospace";
      this.niceTicks(o.xmin, o.xmax).forEach(v => {
        const px = this.x2px(v);
        c.beginPath(); c.moveTo(px, p.t); c.lineTo(px, this.H - p.b); c.stroke();
        c.textAlign = "center"; c.fillText(fmtTick(v), px, this.H - p.b + 14);
      });
      this.niceTicks(o.ymin, o.ymax).forEach(v => {
        const py = this.y2px(v);
        c.beginPath(); c.moveTo(p.l, py); c.lineTo(this.W - p.r, py); c.stroke();
        c.textAlign = "right"; c.fillText(fmtTick(v), p.l - 5, py + 3);
      });
      // 零轴
      c.strokeStyle = "#c9d3e0";
      if (o.ymin < 0 && o.ymax > 0) { const py = this.y2px(0); c.beginPath(); c.moveTo(p.l, py); c.lineTo(this.W - p.r, py); c.stroke(); }
      if (o.xmin < 0 && o.xmax > 0) { const px = this.x2px(0); c.beginPath(); c.moveTo(px, p.t); c.lineTo(px, this.H - p.b); c.stroke(); }
      // 轴标签
      c.fillStyle = "#6a7686"; c.font = "11px 'Segoe UI', sans-serif";
      if (o.xlabel) { c.textAlign = "right"; c.fillText(o.xlabel, this.W - p.r, this.H - 6); }
      if (o.ylabel) { c.save(); c.translate(12, p.t + 4); c.textAlign = "left"; c.fillText(o.ylabel, 0, 0); c.restore(); }
    }
    line(xs, ys, color, width = 1.8, dash = null) {
      const c = this.ctx;
      c.strokeStyle = color; c.lineWidth = width;
      c.setLineDash(dash || []);
      c.beginPath();
      let started = false;
      for (let i = 0; i < xs.length; i++) {
        const px = this.x2px(xs[i]), py = this.y2px(ys[i]);
        if (!isFinite(px) || !isFinite(py)) { started = false; continue; }
        const pyc = Math.max(-2000, Math.min(3000, py));
        if (!started) { c.moveTo(px, pyc); started = true; } else c.lineTo(px, pyc);
      }
      c.stroke(); c.setLineDash([]);
    }
    dot(x, y, color, r = 4) {
      const c = this.ctx;
      c.fillStyle = color; c.beginPath();
      c.arc(this.x2px(x), this.y2px(y), r, 0, Math.PI * 2); c.fill();
    }
    cross(x, y, color, r = 5) {
      const c = this.ctx, px = this.x2px(x), py = this.y2px(y);
      c.strokeStyle = color; c.lineWidth = 1.8; c.beginPath();
      c.moveTo(px - r, py - r); c.lineTo(px + r, py + r);
      c.moveTo(px - r, py + r); c.lineTo(px + r, py - r); c.stroke();
    }
    text(x, y, str, color = "#6a7686", align = "left") {
      const c = this.ctx;
      c.fillStyle = color; c.font = "11px 'Segoe UI', sans-serif"; c.textAlign = align;
      c.fillText(str, this.x2px(x), this.y2px(y));
    }
    legend(items) { // [{label, color, dash}]
      const c = this.ctx, p = this.opt.padding;
      let x = p.l + 10;
      const y = p.t + 12;
      c.font = "11px 'Segoe UI', sans-serif";
      items.forEach(it => {
        c.strokeStyle = it.color; c.lineWidth = 2.2;
        c.setLineDash(it.dash || []);
        c.beginPath(); c.moveTo(x, y - 3); c.lineTo(x + 18, y - 3); c.stroke(); c.setLineDash([]);
        c.fillStyle = "#4a5568"; c.textAlign = "left";
        c.fillText(it.label, x + 23, y);
        x += 33 + c.measureText(it.label).width + 14;
      });
    }
    autoRangeY(seriesArr, padFrac = 0.12) {
      let mn = Infinity, mx = -Infinity;
      seriesArr.forEach(ys => ys.forEach(v => { if (isFinite(v)) { mn = Math.min(mn, v); mx = Math.max(mx, v); } }));
      if (!isFinite(mn)) { mn = -1; mx = 1; }
      if (mx - mn < 1e-9) { mn -= 1; mx += 1; }
      const pad = (mx - mn) * padFrac;
      this.opt.ymin = mn - pad; this.opt.ymax = mx + pad;
    }
  }
  function fmtTick(v) {
    if (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0)) return v.toExponential(0);
    return +v.toFixed(3) + "";
  }

  /* ═════════ 控件 UI ═════════ */
  function buildPanel(mount, { title, desc }) {
    const panel = document.createElement("div");
    panel.className = "sim-panel";
    panel.innerHTML = `<div class="sim-title">${title}</div><div class="sim-desc">${desc}</div>
      <div class="sim-canvas-wrap"></div><div class="sim-controls"></div><div class="sim-readout"></div>`;
    mount.appendChild(panel);
    return {
      panel,
      canvasWrap: panel.querySelector(".sim-canvas-wrap"),
      controls: panel.querySelector(".sim-controls"),
      readout: panel.querySelector(".sim-readout"),
    };
  }
  function addCanvas(ui, w, h) {
    const cv = document.createElement("canvas");
    cv.setAttribute("width", w); cv.setAttribute("height", h);
    ui.canvasWrap.appendChild(cv);
    return cv;
  }
  function addSlider(ui, { label, min, max, step, value, fmt, onChange }) {
    const wrap = document.createElement("div");
    wrap.className = "sim-ctl";
    const f = fmt || (v => (+v).toFixed(2));
    wrap.innerHTML = `<label>${label}</label><input type="range" min="${min}" max="${max}" step="${step}" value="${value}"><span class="val">${f(value)}</span>`;
    const input = wrap.querySelector("input"), val = wrap.querySelector(".val");
    input.addEventListener("input", () => { val.textContent = f(+input.value); onChange && onChange(+input.value); });
    ui.controls.appendChild(wrap);
    return { get: () => +input.value, set: v => { input.value = v; val.textContent = f(v); } };
  }
  function addButton(ui, label, onClick, secondary = false) {
    const b = document.createElement("button");
    b.className = "sim-btn" + (secondary ? " secondary" : "");
    b.textContent = label;
    b.onclick = onClick;
    ui.controls.appendChild(b);
    return b;
  }

  /* ═════════ 生命周期 ═════════ */
  const running = [];
  function loop(fn) {
    const state = { alive: true, id: 0 };
    const step = t => {
      if (!state.alive) return;
      fn(t);
      state.id = requestAnimationFrame(step);
    };
    state.id = requestAnimationFrame(step);
    running.push(state);
    return state;
  }
  function stopAll() {
    running.forEach(s => { s.alive = false; cancelAnimationFrame(s.id); });
    running.length = 0;
  }

  return { num, Plot, buildPanel, addCanvas, addSlider, addButton, loop, stopAll };
})();
