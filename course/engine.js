/* 控制之美 · 交互课程 —— 引擎
   数值计算(线性代数/复数/多项式/离散化/LQR/QP/随机数) + 自适应画布绘图 + 控件 + 动画生命周期 */
window.CE = (function () {
  "use strict";

  /* ═════════════════ 线性代数 ═════════════════ */
  const M = {
    zeros: (r, c) => Array.from({ length: r }, () => new Array(c).fill(0)),
    eye: n => { const I = M.zeros(n, n); for (let i = 0; i < n; i++) I[i][i] = 1; return I; },
    clone: A => A.map(r => r.slice()),
    T: A => A[0].map((_, j) => A.map(r => r[j])),
    add: (A, B) => A.map((r, i) => r.map((v, j) => v + B[i][j])),
    sub: (A, B) => A.map((r, i) => r.map((v, j) => v - B[i][j])),
    scale: (A, s) => A.map(r => r.map(v => v * s)),
    mul: (A, B) => {
      const n = A.length, m = B[0].length, p = B.length, C = M.zeros(n, m);
      for (let i = 0; i < n; i++) for (let k = 0; k < p; k++) {
        const a = A[i][k]; if (a === 0) continue;
        for (let j = 0; j < m; j++) C[i][j] += a * B[k][j];
      }
      return C;
    },
    mulv: (A, x) => A.map(r => r.reduce((s, v, j) => s + v * x[j], 0)),
    dot: (a, b) => a.reduce((s, v, i) => s + v * b[i], 0),
    vadd: (a, b) => a.map((v, i) => v + b[i]),
    vsub: (a, b) => a.map((v, i) => v - b[i]),
    vscale: (a, s) => a.map(v => v * s),
    norm: a => Math.sqrt(a.reduce((s, v) => s + v * v, 0)),
    col: v => v.map(x => [x]),
    flat: A => A.map(r => r[0]),
    hcat: (...Ms) => Ms[0].map((_, i) => [].concat(...Ms.map(X => X[i]))),
    vcat: (...Ms) => [].concat(...Ms.map(X => M.clone(X))),
    blkdiag: (...Ms) => {
      const n = Ms.reduce((s, X) => s + X.length, 0), m = Ms.reduce((s, X) => s + X[0].length, 0);
      const R = M.zeros(n, m); let r0 = 0, c0 = 0;
      Ms.forEach(X => { X.forEach((row, i) => row.forEach((v, j) => R[r0 + i][c0 + j] = v)); r0 += X.length; c0 += X[0].length; });
      return R;
    },
    trace: A => A.reduce((s, r, i) => s + r[i], 0),
    sym: A => M.scale(M.add(A, M.T(A)), 0.5),
    solve: (A0, b0) => {               // 高斯消元（部分主元）
      const n = A0.length, A = A0.map((r, i) => r.concat([b0[i]]));
      for (let c = 0; c < n; c++) {
        let p = c;
        for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
        [A[c], A[p]] = [A[p], A[c]];
        if (Math.abs(A[c][c]) < 1e-14) continue;
        for (let r = 0; r < n; r++) {
          if (r === c) continue;
          const f = A[r][c] / A[c][c]; if (f === 0) continue;
          for (let k = c; k <= n; k++) A[r][k] -= f * A[c][k];
        }
      }
      return A.map((r, i) => r[n] / (Math.abs(r[i]) > 1e-14 ? r[i] : 1e-14));
    },
    inv: A => { const n = A.length, I = M.eye(n); return M.T(I.map((_, j) => M.solve(A, I.map(r => r[j])))); },
    det: A0 => {
      const n = A0.length, A = M.clone(A0); let d = 1;
      for (let c = 0; c < n; c++) {
        let p = c; for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
        if (Math.abs(A[p][c]) < 1e-14) return 0;
        if (p !== c) { [A[c], A[p]] = [A[p], A[c]]; d = -d; }
        d *= A[c][c];
        for (let r = c + 1; r < n; r++) { const f = A[r][c] / A[c][c]; for (let k = c; k < n; k++) A[r][k] -= f * A[c][k]; }
      }
      return d;
    },
    rank: (A0, tol = 1e-9) => {
      const A = M.clone(A0), n = A.length, m = A[0].length; let rk = 0;
      for (let c = 0; c < m && rk < n; c++) {
        let p = rk; for (let r = rk + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
        if (Math.abs(A[p][c]) < tol) continue;
        [A[rk], A[p]] = [A[p], A[rk]];
        for (let r = rk + 1; r < n; r++) { const f = A[r][c] / A[rk][c]; for (let k = c; k < m; k++) A[r][k] -= f * A[rk][k]; }
        rk++;
      }
      return rk;
    },
    // 矩阵指数：缩放-平方 + Taylor
    expm: A => {
      const n = A.length;
      let nrm = 0; A.forEach(r => nrm = Math.max(nrm, r.reduce((s, v) => s + Math.abs(v), 0)));
      let s = Math.max(0, Math.ceil(Math.log2(nrm + 1e-300)) + 1);
      const As = M.scale(A, 1 / Math.pow(2, s));
      let E = M.eye(n), term = M.eye(n);
      for (let k = 1; k < 20; k++) { term = M.scale(M.mul(term, As), 1 / k); E = M.add(E, term); }
      for (let i = 0; i < s; i++) E = M.mul(E, E);
      return E;
    },
    // 零阶保持离散化：用增广矩阵 expm([[A,B],[0,0]]·Ts)
    c2d: (A, B, Ts) => {
      const n = A.length, m = B[0].length;
      const Z = M.zeros(n + m, n + m);
      for (let i = 0; i < n; i++) { for (let j = 0; j < n; j++) Z[i][j] = A[i][j] * Ts; for (let j = 0; j < m; j++) Z[i][n + j] = B[i][j] * Ts; }
      const E = M.expm(Z);
      return { A: E.slice(0, n).map(r => r.slice(0, n)), B: E.slice(0, n).map(r => r.slice(n)) };
    },
    // 2×2 特征值（返回复数 [re,im]）与实特征向量
    eig2: A => {
      const a = A[0][0], b = A[0][1], c = A[1][0], d = A[1][1];
      const tr = a + d, det = a * d - b * c, disc = tr * tr / 4 - det;
      if (disc >= 0) {
        const s = Math.sqrt(disc), l1 = tr / 2 + s, l2 = tr / 2 - s;
        const vec = l => {
          let v = Math.abs(b) > 1e-12 ? [b, l - a] : Math.abs(c) > 1e-12 ? [l - d, c] : (Math.abs(l - a) < 1e-12 ? [1, 0] : [0, 1]);
          const nn = Math.hypot(v[0], v[1]) || 1; return [v[0] / nn, v[1] / nn];
        };
        return { vals: [[l1, 0], [l2, 0]], vecs: [vec(l1), vec(l2)], real: true, tr, det, disc };
      }
      const s = Math.sqrt(-disc);
      return { vals: [[tr / 2, s], [tr / 2, -s]], vecs: null, real: false, tr, det, disc };
    },
    // 离散 LQR：有限时域逆向 Riccati 递推，返回每步增益 K_k（k=0..N-1）
    dlqrFinite: (A, B, Q, R, S, N) => {
      let P = M.clone(S); const Ks = new Array(N);
      for (let k = N - 1; k >= 0; k--) {
        const BtP = M.mul(M.T(B), P);
        const K = M.mul(M.inv(M.add(R, M.mul(BtP, B))), M.mul(BtP, A));
        P = M.sym(M.add(Q, M.mul(M.mul(M.T(A), P), M.sub(A, M.mul(B, K)))));
        Ks[k] = K;
      }
      return { Ks, P0: P };
    },
    // 离散代数 Riccati（迭代至收敛）
    // 离散 Riccati 迭代（与书中 [F1] 模块一致：从 P=S 出发，直到增益 F 收敛）
    dlqr: (A, B, Q, R, S = null, iters = 5000) => {
      let P = M.clone(S || Q), K = null, it = 0;
      for (; it < iters; it++) {
        const BtP = M.mul(M.T(B), P);
        const Kn = M.mul(M.inv(M.add(R, M.mul(BtP, B))), M.mul(BtP, A));
        P = M.sym(M.add(M.add(M.mul(M.mul(M.T(M.sub(A, M.mul(B, Kn))), P), M.sub(A, M.mul(B, Kn))), M.mul(M.mul(M.T(Kn), R), Kn)), Q));
        let d = Infinity; if (K) { d = 0; Kn.forEach((r, a) => r.forEach((v, b) => d = Math.max(d, Math.abs(v - K[a][b])))); }
        K = Kn; if (d < 1e-10) break;
      }
      return { K, P, it };
    },
    // 连续 LQR：离散化(小步长)后用 DARE 近似，再还原（对教学足够精确）
    lqr: (A, B, Q, R) => {
      // 对微分 Riccati 方程积分到稳态
      const At = M.T(A), Bt = M.T(B), Ri = M.inv(R);
      let P = M.clone(Q); const h = 0.002;
      for (let i = 0; i < 60000; i++) {
        const PB = M.mul(P, B);
        const dP = M.add(M.sub(M.add(M.mul(At, P), M.mul(P, A)), M.mul(M.mul(PB, Ri), M.T(PB))), Q);
        let mx = 0; dP.forEach(r => r.forEach(v => mx = Math.max(mx, Math.abs(v))));
        P = M.sym(M.add(P, M.scale(dP, h)));
        if (mx < 1e-9 || !isFinite(mx)) break;
      }
      return { K: M.mul(M.mul(Ri, Bt), P), P };
    },
    // Ackermann 极点配置（单输入），poles 为实数或共轭复数对 [[re,im],...]
    acker: (A, B, poles) => {
      const n = A.length;
      let coef = [1]; // 期望特征多项式
      for (let i = 0; i < poles.length; i++) {
        const [re, im] = poles[i];
        if (Math.abs(im) > 1e-12) { coef = P.mul(coef, [1, -2 * re, re * re + im * im]); i++; }
        else coef = P.mul(coef, [1, -re]);
      }
      // φ(A) = A^n + a1 A^{n-1} + ... + an I
      let phi = M.zeros(n, n), Ak = M.eye(n);
      for (let k = n; k >= 0; k--) { phi = M.add(phi, M.scale(Ak, coef[k])); Ak = M.mul(Ak, A); }
      // 能控性矩阵
      const cols = []; let v = M.clone(B);
      for (let k = 0; k < n; k++) { cols.push(M.flat(v)); v = M.mul(A, v); }
      const Wc = M.T(cols);
      const en = new Array(n).fill(0); en[n - 1] = 1;
      const row = M.mulv(M.T(M.inv(Wc)), en); // e_n^T Wc^{-1}
      return [M.mulv(M.T(phi), row)];
    },
    ctrb: (A, B) => { const n = A.length; let v = M.clone(B); const blocks = []; for (let k = 0; k < n; k++) { blocks.push(v); v = M.mul(A, v); } return M.hcat(...blocks); },
    obsv: (A, C) => { const n = A.length; let v = M.clone(C); const blocks = []; for (let k = 0; k < n; k++) { blocks.push(v); v = M.mul(v, A); } return M.vcat(...blocks); },
    // 特征多项式（Faddeev–LeVerrier），返回从最高次到常数
    charpoly: A => {
      const n = A.length; let Mk = M.zeros(n, n); const c = [1];
      for (let k = 1; k <= n; k++) {
        Mk = M.add(M.mul(A, Mk), M.scale(M.eye(n), c[k - 1]));
        c.push(-M.trace(M.mul(A, Mk)) / k);
      }
      return c;
    },
    eig: A => P.roots(M.charpoly(A)),
  };

  /* ═════════════════ 复数与多项式 ═════════════════ */
  const C = {
    add: (a, b) => [a[0] + b[0], a[1] + b[1]],
    sub: (a, b) => [a[0] - b[0], a[1] - b[1]],
    mul: (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]],
    div: (a, b) => { const d = b[0] * b[0] + b[1] * b[1] || 1e-300; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; },
    abs: a => Math.hypot(a[0], a[1]),
    arg: a => Math.atan2(a[1], a[0]),
  };
  const P = {
    mul: (a, b) => { const r = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => r[i + j] += x * y)); return r; },
    add: (a, b) => { const n = Math.max(a.length, b.length), r = new Array(n).fill(0); a.forEach((x, i) => r[n - a.length + i] += x); b.forEach((x, i) => r[n - b.length + i] += x); return r; },
    scale: (a, s) => a.map(x => x * s),
    evalC: (c, z) => { let r = [0, 0]; for (const k of c) r = C.add(C.mul(r, z), [k, 0]); return r; },
    eval: (c, x) => c.reduce((r, k) => r * x + k, 0),
    fromRoots: rs => {               // rs: [[re,im],...]，复根需成对
      let c = [1];
      for (let i = 0; i < rs.length; i++) {
        const [re, im] = rs[i];
        if (Math.abs(im) > 1e-12) { c = P.mul(c, [1, -2 * re, re * re + im * im]); i++; } else c = P.mul(c, [1, -re]);
      }
      return c;
    },
    roots: coeffs => {               // Durand–Kerner
      let c = coeffs.slice(); while (c.length > 1 && Math.abs(c[0]) < 1e-14) c.shift();
      const n = c.length - 1; if (n < 1) return [];
      c = c.map(v => v / c[0]);
      if (n === 1) return [[-c[1], 0]];
      if (n === 2) {
        const b = c[1], k = c[2], d = b * b / 4 - k;
        return d >= 0 ? [[-b / 2 + Math.sqrt(d), 0], [-b / 2 - Math.sqrt(d), 0]] : [[-b / 2, Math.sqrt(-d)], [-b / 2, -Math.sqrt(-d)]];
      }
      let R = 1 + Math.max(...c.slice(1).map(Math.abs));
      let z = Array.from({ length: n }, (_, i) => [R * 0.5 * Math.cos(2 * Math.PI * i / n + 0.4), R * 0.5 * Math.sin(2 * Math.PI * i / n + 0.4)]);
      for (let it = 0; it < 500; it++) {
        let mx = 0;
        for (let i = 0; i < n; i++) {
          let den = [1, 0];
          for (let j = 0; j < n; j++) if (j !== i) den = C.mul(den, C.sub(z[i], z[j]));
          const d = C.div(P.evalC(c, z[i]), den);
          z[i] = C.sub(z[i], d); mx = Math.max(mx, Math.abs(d[0]) + Math.abs(d[1]));
        }
        if (mx < 1e-13) break;
      }
      return z.map(r => [r[0], Math.abs(r[1]) < 1e-9 ? 0 : r[1]]).sort((a, b) => b[0] - a[0] || b[1] - a[1]);
    },
    // G(jw)：返回 {mag, phase(deg), re, im}
    freq: (num, den, w) => {
      const v = C.div(P.evalC(num, [0, w]), P.evalC(den, [0, w]));
      return { re: v[0], im: v[1], mag: C.abs(v), ph: C.arg(v) * 180 / Math.PI };
    },
  };

  /* ═════════════════ 传递函数 → 状态空间 & 仿真 ═════════════════ */
  // 可控标准型实现（严格真分式，num 次数 < den 次数）
  function tf2ss(num, den) {
    const a0 = den[0]; den = den.map(v => v / a0); num = num.map(v => v / a0);
    const n = den.length - 1;
    const nn = new Array(n + 1 - num.length).fill(0).concat(num); // 补齐到 n+1
    const D = nn[0];
    const A = M.zeros(n, n);
    for (let i = 0; i < n - 1; i++) A[i][i + 1] = 1;
    for (let j = 0; j < n; j++) A[n - 1][j] = -den[n - j];
    const B = M.zeros(n, 1); B[n - 1][0] = 1;
    const Cm = [new Array(n).fill(0)];
    for (let j = 0; j < n; j++) Cm[0][j] = nn[n - j] - D * den[n - j];
    return { A, B, C: Cm, D };
  }
  // 线性系统时间响应：u(t) 为函数，返回 {t, y}
  function lsim(num, den, uFn, tEnd, dt = 0.005, x0 = null) {
    const { A, B, C: Cm, D } = tf2ss(num, den);
    const n = A.length;
    const { A: Ad, B: Bd } = M.c2d(A, B, dt);
    let x = x0 ? x0.slice() : new Array(n).fill(0);
    const ts = [], ys = [];
    for (let k = 0, t = 0; t <= tEnd + 1e-9; k++, t = k * dt) {
      const u = uFn(t);
      ts.push(t); ys.push(M.dot(Cm[0], x) + D * u);
      const xn = M.mulv(Ad, x); for (let i = 0; i < n; i++) xn[i] += Bd[i][0] * u; x = xn;
    }
    return { t: ts, y: ys };
  }
  const step = (num, den, tEnd, dt) => lsim(num, den, () => 1, tEnd, dt);
  // 闭环：单位负反馈 G = n/d → n/(d+n)
  const feedback = (num, den) => ({ num: num.slice(), den: P.add(den, num) });

  // 阶跃响应指标
  function stepInfo(t, y, yss) {
    if (yss === undefined) yss = y[y.length - 1];
    let peak = -Infinity, tp = 0;
    y.forEach((v, i) => { if (v > peak) { peak = v; tp = t[i]; } });
    let tr = null, t10 = null, t90 = null;
    for (let i = 0; i < y.length; i++) {
      if (t10 === null && y[i] >= 0.1 * yss) t10 = t[i];
      if (t90 === null && y[i] >= 0.9 * yss) { t90 = t[i]; break; }
    }
    if (t10 !== null && t90 !== null) tr = t90 - t10;
    let ts = 0;
    for (let i = y.length - 1; i >= 0; i--) if (Math.abs(y[i] - yss) > 0.02 * Math.abs(yss)) { ts = t[Math.min(i + 1, t.length - 1)]; break; }
    const os = yss !== 0 ? Math.max(0, (peak - yss) / Math.abs(yss) * 100) : 0;
    return { peak, tp, tr, ts, os, yss };
  }

  /* ═════════════════ 二次规划：Hildreth 对偶坐标上升 ═════════════════
     min ½xᵀHx + fᵀx  s.t.  Aineq x ≤ b ；返回 {x, active} */
  function qp(H, f, Ain, bin, iters = 300) {
    const Hi = M.inv(H);
    const x0 = M.vscale(M.mulv(Hi, f), -1);
    if (!Ain || !Ain.length) return { x: x0, active: [] };
    const viol = Ain.map((r, i) => M.dot(r, x0) - bin[i]);
    if (viol.every(v => v <= 1e-9)) return { x: x0, active: [] };
    const Pm = M.mul(M.mul(Ain, Hi), M.T(Ain));
    const d = M.vadd(M.mulv(M.mul(Ain, Hi), f), bin);
    const m = Ain.length; const lam = new Array(m).fill(0);
    for (let it = 0; it < iters; it++) {
      let ch = 0;
      for (let i = 0; i < m; i++) {
        let w = d[i]; for (let j = 0; j < m; j++) if (j !== i) w += Pm[i][j] * lam[j];
        const nl = Math.max(0, -w / Pm[i][i]);
        ch = Math.max(ch, Math.abs(nl - lam[i])); lam[i] = nl;
      }
      if (ch < 1e-10) break;
    }
    const x = M.vsub(x0, M.mulv(M.mul(Hi, M.T(Ain)), lam));
    return { x, active: lam.map((l, i) => l > 1e-8 ? i : -1).filter(i => i >= 0) };
  }

  /* 预处理版 QP：H 与约束矩阵 Ain 固定（MPC 每步只有 f、b 变化），可热启动 λ */
  function qpPrep(H, Ain) {
    const n = H.length, Hi = M.inv(H), m = Ain ? Ain.length : 0;
    const AHi = m ? M.mul(Ain, Hi) : null, Pm = m ? M.mul(AHi, M.T(Ain)) : null;
    const HiAt = m ? M.T(AHi) : null; // n×m
    let lamPrev = new Float64Array(m);
    return function solve(f, b, opt = {}) {
      const x0 = M.vscale(M.mulv(Hi, f), -1);
      if (!m) return { x: x0, active: 0, iters: 0 };
      const d = new Float64Array(m);
      for (let i = 0; i < m; i++) { let s = b[i]; const r = AHi[i]; for (let j = 0; j < n; j++) s += r[j] * f[j]; d[i] = s; }
      // 无约束解可行则直接返回
      let feas = true; for (let i = 0; i < m && feas; i++) { let s = 0; const r = Ain[i]; for (let j = 0; j < n; j++) s += r[j] * x0[j]; if (s > b[i] + 1e-9) feas = false; }
      if (feas) { lamPrev = new Float64Array(m); return { x: x0, active: 0, iters: 0 }; }
      const lam = opt.cold ? new Float64Array(m) : Float64Array.from(lamPrev);
      const iters = opt.iters || 400, tol = opt.tol || 1e-9; let it = 0;
      for (; it < iters; it++) {
        let ch = 0;
        for (let i = 0; i < m; i++) {
          const Pi = Pm[i]; let w = d[i];
          for (let j = 0; j < m; j++) if (j !== i) w += Pi[j] * lam[j];
          const nl = Math.max(0, -w / Pi[i]); ch = Math.max(ch, Math.abs(nl - lam[i])); lam[i] = nl;
        }
        if (ch < tol) break;
      }
      lamPrev = lam;
      const x = x0.slice(); let act = 0;
      for (let i = 0; i < m; i++) if (lam[i] > 1e-10) { act++; for (let j = 0; j < n; j++) x[j] -= HiAt[j][i] * lam[i]; }
      return { x, active: act, iters: it };
    };
  }

  /* ═════════════════ MPC 矩阵构造（书中 [F4] [F6] 模块） ═════════════════ */
  function mpcMatrices(A, B, Q, R, S, Np) {
    const n = A.length, p = B[0].length;
    const Phi = M.zeros(n * Np, n), Gam = M.zeros(n * Np, p * Np);
    let Ap = M.eye(n); const pw = [M.eye(n)];
    for (let i = 1; i <= Np; i++) { Ap = M.mul(Ap, A); pw.push(Ap); }
    for (let i = 0; i < Np; i++) {
      for (let r = 0; r < n; r++) for (let c = 0; c < n; c++) Phi[i * n + r][c] = pw[i + 1][r][c];
      for (let j = 0; j <= i; j++) { const blk = M.mul(pw[i - j], B); for (let r = 0; r < n; r++) for (let c = 0; c < p; c++) Gam[i * n + r][j * p + c] = blk[r][c]; }
    }
    const Om = M.blkdiag(...Array.from({ length: Np }, (_, i) => i === Np - 1 ? S : Q));
    const Psi = M.blkdiag(...Array.from({ length: Np }, () => R));
    const GtO = M.mul(M.T(Gam), Om);
    const H = M.sym(M.add(M.mul(GtO, Gam), Psi)), F = M.mul(GtO, Phi);
    return { Phi, Gam, Om, Psi, H, F, n, p, Np };
  }
  // 上下限约束 → M U ≤ β̄ + b x（去掉 ±∞ 行）
  function mpcConstraints(mm, ulo, uhi, xlo, xhi) {
    const { Phi, Gam, n, p, Np } = mm, rows = [], rhs0 = [], bx = [];
    // 输入约束（每一步）
    for (let i = 0; i < Np; i++) for (let j = 0; j < p; j++) {
      if (isFinite(uhi[j])) { const r = new Array(p * Np).fill(0); r[i * p + j] = 1; rows.push(r); rhs0.push(uhi[j]); bx.push(new Array(n).fill(0)); }
      if (isFinite(ulo[j])) { const r = new Array(p * Np).fill(0); r[i * p + j] = -1; rows.push(r); rhs0.push(-ulo[j]); bx.push(new Array(n).fill(0)); }
    }
    // 状态约束：x_{k+i|k}, i=1..Np，X = Φx + ΓU
    for (let i = 0; i < Np; i++) for (let j = 0; j < n; j++) {
      const g = Gam[i * n + j], ph = Phi[i * n + j];
      if (isFinite(xhi[j])) { rows.push(g.slice()); rhs0.push(xhi[j]); bx.push(ph.map(v => -v)); }
      if (isFinite(xlo[j])) { rows.push(g.map(v => -v)); rhs0.push(-xlo[j]); bx.push(ph.slice()); }
    }
    return { Mc: rows, beta: rhs0, b: bx, rhs: x => rhs0.map((v, i) => v + M.dot(bx[i], x)) };
  }

  /* ═════════════════ 随机数（可复现） ═════════════════ */
  function rng(seed = 1) {
    let a = seed >>> 0;
    const u = () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
    let spare = null;
    const n = () => {
      if (spare !== null) { const s = spare; spare = null; return s; }
      let x, y, r; do { x = 2 * u() - 1; y = 2 * u() - 1; r = x * x + y * y; } while (r >= 1 || r === 0);
      const k = Math.sqrt(-2 * Math.log(r) / r); spare = y * k; return x * k;
    };
    return { u, n };
  }

  /* ═════════════════ 颜色 ═════════════════ */
  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const COL = {};
  function refreshColors() {
    ["c1", "c2", "c3", "c4", "c5", "c6", "ink", "ink2", "grid", "axis", "plotbg", "muted"].forEach(k => COL[k] = css("--" + k) || "#888");
  }
  refreshColors();

  /* ═════════════════ 绘图：自适应宽度的画布 ═════════════════ */
  const plots = new Set();
  class Plot {
    constructor(host, opt = {}) {
      this.opt = Object.assign({
        height: 260, x: [0, 10], y: [-1, 1], xlabel: "", ylabel: "", title: "",
        logx: false, equal: false, grid: true, pad: { l: 48, r: 14, t: 24, b: 32 }, xticks: null, yticks: null, aspect: null,
      }, opt);
      this.host = host;
      this.cv = document.createElement("canvas");
      this.cv.className = "plot";
      host.appendChild(this.cv);
      this.ctx = this.cv.getContext("2d");
      this.render = () => this.clear();
      this.resize();
      plots.add(this);
    }
    resize() {
      const w = Math.max(160, this.host.clientWidth || 600);
      const h = this.opt.aspect ? Math.round(w * this.opt.aspect) : this.opt.height;
      const dpr = window.devicePixelRatio || 1;
      this.W = w; this.H = h;
      this.cv.style.width = w + "px"; this.cv.style.height = h + "px";
      this.cv.width = Math.round(w * dpr); this.cv.height = Math.round(h * dpr);
      this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      this._fitEqual();
    }
    _fitEqual() {
      if (!this.opt.equal) return;
      // 保持 x/y 等比例：以 y 范围为准扩展 x
      const p = this.opt.pad, pw = this.W - p.l - p.r, ph = this.H - p.t - p.b;
      if (!this._base) this._base = { x: this.opt.x.slice(), y: this.opt.y.slice() };
      const bx = this._base.x, by = this._base.y;
      const sx = (bx[1] - bx[0]) / pw, sy = (by[1] - by[0]) / ph, s = Math.max(sx, sy);
      const cx = (bx[0] + bx[1]) / 2, cy = (by[0] + by[1]) / 2;
      this.opt.x = [cx - s * pw / 2, cx + s * pw / 2]; this.opt.y = [cy - s * ph / 2, cy + s * ph / 2];
    }
    setRange(x, y) { if (x) this.opt.x = x.slice(); if (y) this.opt.y = y.slice(); if (this.opt.equal) { this._base = { x: this.opt.x.slice(), y: this.opt.y.slice() }; this._fitEqual(); } }
    fx(x) { const p = this.opt.pad, [a, b] = this.opt.x; if (this.opt.logx) { x = Math.log10(Math.max(x, 1e-300)); return p.l + (x - Math.log10(a)) / (Math.log10(b) - Math.log10(a)) * (this.W - p.l - p.r); } return p.l + (x - a) / (b - a) * (this.W - p.l - p.r); }
    fy(y) { const p = this.opt.pad, [a, b] = this.opt.y; return this.H - p.b - (y - a) / (b - a) * (this.H - p.t - p.b); }
    ix(px) { const p = this.opt.pad, [a, b] = this.opt.x; const r = (px - p.l) / (this.W - p.l - p.r); return this.opt.logx ? Math.pow(10, Math.log10(a) + r * (Math.log10(b) - Math.log10(a))) : a + r * (b - a); }
    iy(py) { const p = this.opt.pad, [a, b] = this.opt.y; return a + (this.H - p.b - py) / (this.H - p.t - p.b) * (b - a); }
    ticks(a, b, n = 6) {
      const span = b - a || 1, st0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(st0))), nr = st0 / mag;
      const st = (nr < 1.5 ? 1 : nr < 3.5 ? 2 : nr < 7.5 ? 5 : 10) * mag, out = [];
      for (let v = Math.ceil(a / st) * st; v <= b + st * 1e-6; v += st) out.push(+v.toFixed(10));
      return out;
    }
    clear() {
      const c = this.ctx, o = this.opt, p = o.pad;
      c.clearRect(0, 0, this.W, this.H);
      c.fillStyle = COL.plotbg; c.fillRect(0, 0, this.W, this.H);
      if (o.title) { c.fillStyle = COL.ink2; c.font = "600 12px var(--font-ui), sans-serif"; c.font = "600 12px 'Segoe UI','Microsoft YaHei',sans-serif"; c.textAlign = "left"; c.fillText(o.title, p.l, 15); }
      if (!o.grid) return;
      c.lineWidth = 1; c.font = "10.5px Consolas, 'Cascadia Mono', monospace";
      const xt = o.xticks || (o.logx ? this._logTicks() : this.ticks(o.x[0], o.x[1], Math.max(3, Math.round(this.W / 90))));
      const yt = o.yticks || this.ticks(o.y[0], o.y[1], Math.max(3, Math.round(this.H / 50)));
      xt.forEach(v => {
        const x = this.fx(v); if (x < p.l - 1 || x > this.W - p.r + 1) return;
        c.strokeStyle = COL.grid; c.beginPath(); c.moveTo(x, p.t); c.lineTo(x, this.H - p.b); c.stroke();
        c.fillStyle = COL.muted; c.textAlign = "center"; c.fillText(fmt(v), x, this.H - p.b + 13);
      });
      yt.forEach(v => {
        const y = this.fy(v); if (y < p.t - 1 || y > this.H - p.b + 1) return;
        c.strokeStyle = COL.grid; c.beginPath(); c.moveTo(p.l, y); c.lineTo(this.W - p.r, y); c.stroke();
        c.fillStyle = COL.muted; c.textAlign = "right"; c.fillText(fmt(v), p.l - 5, y + 3.5);
      });
      c.strokeStyle = COL.axis; c.lineWidth = 1.1;
      if (o.y[0] < 0 && o.y[1] > 0) { const y = this.fy(0); c.beginPath(); c.moveTo(p.l, y); c.lineTo(this.W - p.r, y); c.stroke(); }
      if (!o.logx && o.x[0] < 0 && o.x[1] > 0) { const x = this.fx(0); c.beginPath(); c.moveTo(x, p.t); c.lineTo(x, this.H - p.b); c.stroke(); }
      c.strokeRect(p.l, p.t, this.W - p.l - p.r, this.H - p.t - p.b);
      c.fillStyle = COL.ink2; c.font = "11.5px 'Segoe UI','Microsoft YaHei',sans-serif";
      if (o.xlabel) { c.textAlign = "right"; c.fillText(o.xlabel, this.W - p.r, this.H - 4); }
      if (o.ylabel) { c.textAlign = "right"; c.fillText(o.ylabel, this.W - p.r, p.t - 7); }
    }
    _logTicks() { const [a, b] = this.opt.x, out = []; for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++) out.push(Math.pow(10, e)); return out; }
    clip(fn) { const c = this.ctx, p = this.opt.pad; c.save(); c.beginPath(); c.rect(p.l, p.t, this.W - p.l - p.r, this.H - p.t - p.b); c.clip(); fn(c); c.restore(); }
    line(xs, ys, color = COL.c1, w = 2, dash) {
      this.clip(c => {
        c.strokeStyle = color; c.lineWidth = w; c.setLineDash(dash || []); c.lineJoin = "round";
        c.beginPath(); let on = false;
        for (let i = 0; i < xs.length; i++) {
          const X = this.fx(xs[i]), Y = this.fy(ys[i]);
          if (!isFinite(X) || !isFinite(Y) || Math.abs(Y) > 1e5) { on = false; continue; }
          on ? c.lineTo(X, Y) : (c.moveTo(X, Y), on = true);
        }
        c.stroke(); c.setLineDash([]);
      });
    }
    fn(f, color, w, dash, n = 400) { const [a, b] = this.opt.x, xs = [], ys = []; for (let i = 0; i <= n; i++) { const x = this.opt.logx ? Math.pow(10, Math.log10(a) + i / n * (Math.log10(b) - Math.log10(a))) : a + (b - a) * i / n; xs.push(x); ys.push(f(x)); } this.line(xs, ys, color, w, dash); }
    area(xs, ys1, ys2, color) {
      this.clip(c => {
        c.fillStyle = color; c.beginPath();
        xs.forEach((x, i) => i ? c.lineTo(this.fx(x), this.fy(ys1[i])) : c.moveTo(this.fx(x), this.fy(ys1[i])));
        for (let i = xs.length - 1; i >= 0; i--) c.lineTo(this.fx(xs[i]), this.fy(ys2[i]));
        c.closePath(); c.fill();
      });
    }
    rect(x0, y0, x1, y1, fill, stroke) { this.clip(c => { const X = this.fx(x0), Y = this.fy(y1), W = this.fx(x1) - X, H = this.fy(y0) - Y; if (fill) { c.fillStyle = fill; c.fillRect(X, Y, W, H); } if (stroke) { c.strokeStyle = stroke; c.lineWidth = 1; c.strokeRect(X, Y, W, H); } }); }
    hline(y, color = COL.muted, dash = [5, 4], w = 1.2) { this.line([this.opt.x[0], this.opt.x[1]], [y, y], color, w, dash); }
    vline(x, color = COL.muted, dash = [5, 4], w = 1.2) { this.line([x, x], [this.opt.y[0], this.opt.y[1]], color, w, dash); }
    inside(x, y, m = 2) { const X = this.fx(x), Y = this.fy(y), p = this.opt.pad; return X >= p.l - m && X <= this.W - p.r + m && Y >= p.t - m && Y <= this.H - p.b + m; }
    dot(x, y, color = COL.c1, r = 4.5, ring) {
      if (!this.inside(x, y)) return;
      const c = this.ctx; c.fillStyle = color; c.beginPath(); c.arc(this.fx(x), this.fy(y), r, 0, 7); c.fill();
      if (ring) { c.strokeStyle = ring; c.lineWidth = 2; c.stroke(); }
    }
    cross(x, y, color = COL.c4, r = 6, w = 2.4) {
      if (!this.inside(x, y)) return;
      const c = this.ctx, X = this.fx(x), Y = this.fy(y);
      c.strokeStyle = color; c.lineWidth = w; c.beginPath();
      c.moveTo(X - r, Y - r); c.lineTo(X + r, Y + r); c.moveTo(X - r, Y + r); c.lineTo(X + r, Y - r); c.stroke();
    }
    circle(x, y, color = COL.c1, r = 6, w = 2.2) { if (!this.inside(x, y)) return; const c = this.ctx; c.strokeStyle = color; c.lineWidth = w; c.beginPath(); c.arc(this.fx(x), this.fy(y), r, 0, 7); c.stroke(); }
    text(x, y, s, color = COL.ink2, align = "left", size = 11.5, bold = false) {
      const c = this.ctx; c.fillStyle = color; c.textAlign = align;
      c.font = `${bold ? "600 " : ""}${size}px 'Segoe UI','Microsoft YaHei',sans-serif`;
      c.fillText(s, this.fx(x), this.fy(y));
    }
    arrow(x0, y0, x1, y1, color = COL.ink2, w = 1.8, head = 8) {
      const c = this.ctx, X0 = this.fx(x0), Y0 = this.fy(y0), X1 = this.fx(x1), Y1 = this.fy(y1);
      const a = Math.atan2(Y1 - Y0, X1 - X0), L = Math.hypot(X1 - X0, Y1 - Y0); if (L < 1) return;
      c.strokeStyle = color; c.fillStyle = color; c.lineWidth = w;
      c.beginPath(); c.moveTo(X0, Y0); c.lineTo(X1, Y1); c.stroke();
      const h = Math.min(head, L * 0.6);
      c.beginPath(); c.moveTo(X1, Y1); c.lineTo(X1 - h * Math.cos(a - 0.4), Y1 - h * Math.sin(a - 0.4)); c.lineTo(X1 - h * Math.cos(a + 0.4), Y1 - h * Math.sin(a + 0.4)); c.closePath(); c.fill();
    }
    legend(items, pos = "tl") { // [{label,color,dash}]
      const c = this.ctx, p = this.opt.pad;
      c.font = "11.5px 'Segoe UI','Microsoft YaHei',sans-serif";
      const widths = items.map(it => c.measureText(it.label).width + 30);
      const tot = widths.reduce((a, b) => a + b, 0);
      let x = pos === "tr" ? this.W - p.r - tot - 6 : p.l + 8; const y = p.t + 14;
      c.fillStyle = COL.plotbg; c.globalAlpha = .85; c.fillRect(x - 4, y - 11, tot + 4, 16); c.globalAlpha = 1;
      items.forEach((it, i) => {
        c.strokeStyle = it.color; c.lineWidth = 2.4; c.setLineDash(it.dash || []);
        c.beginPath(); c.moveTo(x, y - 3.5); c.lineTo(x + 18, y - 3.5); c.stroke(); c.setLineDash([]);
        c.fillStyle = COL.ink2; c.textAlign = "left"; c.fillText(it.label, x + 22, y);
        x += widths[i];
      });
    }
    autoY(arrs, pad = 0.1, incl0 = false) {
      let a = Infinity, b = -Infinity;
      arrs.forEach(ys => ys.forEach(v => { if (isFinite(v) && Math.abs(v) < 1e8) { a = Math.min(a, v); b = Math.max(b, v); } }));
      if (incl0) { a = Math.min(a, 0); b = Math.max(b, 0); }
      if (!isFinite(a)) { a = -1; b = 1; } if (b - a < 1e-9) { a -= 1; b += 1; }
      const d = (b - a) * pad; this.setRange(null, [a - d, b + d]);
    }
    // 拖拽交互：onDown(x,y)->返回是否抓取, onMove(x,y), onUp()
    drag({ down, move, up, hover }) {
      const pos = e => { const r = this.cv.getBoundingClientRect(); return [this.ix(e.clientX - r.left), this.iy(e.clientY - r.top), e.clientX - r.left, e.clientY - r.top]; };
      let active = false;
      this.cv.style.touchAction = "none";
      this.cv.addEventListener("pointerdown", e => { const [x, y, X, Y] = pos(e); if (down && down(x, y, X, Y) !== false) { active = true; this.cv.setPointerCapture(e.pointerId); e.preventDefault(); } });
      this.cv.addEventListener("pointermove", e => { const [x, y, X, Y] = pos(e); if (active) move && move(x, y, X, Y); else if (hover) this.cv.style.cursor = hover(x, y, X, Y) ? "grab" : "default"; });
      const end = () => { if (active) { active = false; up && up(); } };
      this.cv.addEventListener("pointerup", end); this.cv.addEventListener("pointercancel", end);
    }
    // 像素距离（判定拖拽是否命中）
    pxDist(x1, y1, x2, y2) { return Math.hypot(this.fx(x1) - this.fx(x2), this.fy(y1) - this.fy(y2)); }
  }
  function fmt(v) {
    if (v === 0) return "0";
    const a = Math.abs(v);
    if (a >= 1e4 || a < 1e-3) return v.toExponential(0).replace("e+", "e");
    return +v.toFixed(4) + "";
  }
  let rT = 0;
  window.addEventListener("resize", () => { clearTimeout(rT); rT = setTimeout(() => plots.forEach(p => { if (p.cv.isConnected) { p.resize(); p.render(); } else plots.delete(p); }), 80); });

  /* ═════════════════ 自由画布（动画场景） ═════════════════ */
  class Scene {
    constructor(host, { height = 240, aspect = null } = {}) {
      this.host = host; this.height = height; this.aspect = aspect;
      this.cv = document.createElement("canvas"); this.cv.className = "plot scene"; host.appendChild(this.cv);
      this.ctx = this.cv.getContext("2d"); this.render = () => { }; this.resize(); plots.add(this);
    }
    resize() {
      const w = Math.max(160, this.host.clientWidth || 600), h = this.aspect ? Math.round(w * this.aspect) : this.height, dpr = window.devicePixelRatio || 1;
      this.W = w; this.H = h; this.cv.style.width = w + "px"; this.cv.style.height = h + "px";
      this.cv.width = Math.round(w * dpr); this.cv.height = Math.round(h * dpr); this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    }
    clear(bg = COL.plotbg) { const c = this.ctx; c.clearRect(0, 0, this.W, this.H); c.fillStyle = bg; c.fillRect(0, 0, this.W, this.H); }
  }

  /* ═════════════════ 控件 ═════════════════ */
  function h(tag, attrs = {}, ...kids) {
    const el = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs || {})) {
      if (k === "class") el.className = v; else if (k === "html") el.innerHTML = v;
      else if (k.startsWith("on")) el.addEventListener(k.slice(2), v); else if (v !== null && v !== undefined) el.setAttribute(k, v);
    }
    kids.flat().forEach(k => k != null && el.appendChild(typeof k === "string" ? document.createTextNode(k) : k));
    return el;
  }
  /* 构造一个仿真卡片：返回 {root, stage, side, controls, readout, add*} */
  function card(mount, { title = "", hint = "", layout = "wide" } = {}) {
    const root = h("div", { class: "sim-card " + layout });
    if (title) root.appendChild(h("div", { class: "sim-head" }, h("span", { class: "sim-badge" }, "交互"), h("span", { class: "sim-title", html: title })));
    if (hint) root.appendChild(h("div", { class: "sim-hint", html: hint }));
    const body = h("div", { class: "sim-body" });
    const stage = h("div", { class: "sim-stage" });
    const side = h("div", { class: "sim-side" });
    const controls = h("div", { class: "sim-controls" });
    const readout = h("div", { class: "sim-readout" });
    side.append(controls, readout);
    body.append(stage, side); root.appendChild(body);
    mount.appendChild(root);
    const ui = {
      root, stage, side, controls, readout,
      plotHost(cls = "") { const d = h("div", { class: "plot-host " + cls }); stage.appendChild(d); return d; },
      row(...cls) { const d = h("div", { class: "plot-row" }); stage.appendChild(d); return cls.map(c => { const x = h("div", { class: "plot-host " + c }); d.appendChild(x); return x; }); },
      slider(o) { return slider(controls, o); },
      buttons(list) { const d = h("div", { class: "btn-row" }); list.forEach(([lab, fn, cls]) => d.appendChild(h("button", { class: "btn " + (cls || ""), onclick: fn }, lab))); controls.appendChild(d); return d; },
      seg(o) { return seg(controls, o); },
      check(o) { return check(controls, o); },
      out(html) { readout.innerHTML = html; renderMath(readout); },
      note(html) { const d = h("div", { class: "sim-note", html }); controls.appendChild(d); renderMath(d); return d; },
    };
    return ui;
  }
  function slider(mount, { label, min, max, step = 0.01, value, fmt: f, onInput, color, log = false }) {
    const fm = f || (v => (+v).toFixed(step >= 1 ? 0 : step >= 0.1 ? 1 : 2));
    const toPos = v => log ? Math.log10(v) : v, fromPos = p => log ? Math.pow(10, p) : p;
    const inp = h("input", { type: "range", min: toPos(min), max: toPos(max), step: log ? "any" : step, value: toPos(value) });
    const val = h("span", { class: "val" }, fm(value));
    const lab = h("label", { html: label });
    const wrap = h("div", { class: "ctl" }, h("div", { class: "ctl-top" }, lab, val), inp);
    if (color) wrap.style.setProperty("--sc", color);
    mount.appendChild(wrap); renderMath(lab);
    const get = () => fromPos(+inp.value);
    inp.addEventListener("input", () => { val.textContent = fm(get()); onInput && onInput(get()); });
    return { get, set(v, fire) { inp.value = toPos(v); val.textContent = fm(v); if (fire) onInput && onInput(v); }, el: wrap };
  }
  function seg(mount, { label, options, value, onChange }) {
    const wrap = h("div", { class: "ctl seg-wrap" });
    if (label) wrap.appendChild(h("div", { class: "ctl-top" }, h("label", { html: label })));
    const row = h("div", { class: "seg" });
    let cur = value ?? options[0][0];
    const btns = options.map(([v, lab]) => { const b = h("button", { class: v === cur ? "on" : "", html: lab, onclick: () => { cur = v; btns.forEach(x => x.classList.toggle("on", x === b)); onChange && onChange(v); } }); row.appendChild(b); return b; });
    wrap.appendChild(row); mount.appendChild(wrap); renderMath(wrap);
    return { get: () => cur, set(v) { cur = v; btns.forEach((b, i) => b.classList.toggle("on", options[i][0] === v)); } };
  }
  function check(mount, { label, value = false, onChange }) {
    const inp = h("input", { type: "checkbox" }); inp.checked = value;
    const wrap = h("label", { class: "ctl chk" }, inp, h("span", { html: label }));
    inp.addEventListener("change", () => onChange && onChange(inp.checked));
    mount.appendChild(wrap); renderMath(wrap);
    return { get: () => inp.checked, set(v) { inp.checked = v; } };
  }

  /* ═════════════════ 数学渲染 ═════════════════ */
  function renderMath(el) {
    if (!window.renderMathInElement || !el) return;
    renderMathInElement(el, {
      delimiters: [{ left: "$$", right: "$$", display: true }, { left: "$", right: "$", display: false }, { left: "\\[", right: "\\]", display: true }, { left: "\\(", right: "\\)", display: false }],
      throwOnError: false, trust: true, strict: false,
      macros: { "\\bm": "\\boldsymbol{#1}", "\\R": "\\mathbb{R}" },
    });
  }
  function tex(s, display = false) {
    try { return katex.renderToString(s, { displayMode: display, throwOnError: false, trust: true, strict: false, macros: { "\\bm": "\\boldsymbol{#1}" } }); } catch (e) { return s; }
  }

  /* ═════════════════ 动画循环（离开视口自动暂停） ═════════════════ */
  const loops = [];
  function loop(el, fn) {
    const st = { alive: true, visible: true, id: 0, last: 0, paused: false };
    const io = new IntersectionObserver(es => es.forEach(e => st.visible = e.isIntersecting), { rootMargin: "100px" });
    io.observe(el);
    const tick = t => {
      if (!st.alive) return;
      const dt = st.last ? Math.min(0.05, (t - st.last) / 1000) : 0.016; st.last = t;
      if (st.visible && !st.paused) fn(dt, t);
      st.id = requestAnimationFrame(tick);
    };
    st.id = requestAnimationFrame(tick);
    st.stop = () => { st.alive = false; cancelAnimationFrame(st.id); io.disconnect(); };
    st.fn = fn; loops.push(st); return st;
  }
  // 调试/测试用：手动推进所有动画循环（后台标签页不触发 requestAnimationFrame 时使用）
  function tickAll(n = 1, dt = 1 / 60) { for (let i = 0; i < n; i++) loops.forEach(s => s.alive && s.fn(dt, performance.now())); }
  function stopAll() { loops.forEach(s => s.stop()); loops.length = 0; plots.clear(); }

  // 数值格式化
  const f2 = (v, d = 2) => (v === null || v === undefined || !isFinite(v)) ? "—" : (+v).toFixed(d);
  const cfmt = (z, d = 2) => Math.abs(z[1]) < 1e-9 ? f2(z[0], d) : `${f2(z[0], d)} ${z[1] >= 0 ? "+" : "−"} ${f2(Math.abs(z[1]), d)}j`;
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  const logspace = (a, b, n) => Array.from({ length: n }, (_, i) => Math.pow(10, a + (b - a) * i / (n - 1)));
  const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
  const rk4 = (f, x, t, dt) => {
    const k1 = f(x, t), k2 = f(M.vadd(x, M.vscale(k1, dt / 2)), t + dt / 2), k3 = f(M.vadd(x, M.vscale(k2, dt / 2)), t + dt / 2), k4 = f(M.vadd(x, M.vscale(k3, dt)), t + dt);
    return x.map((v, i) => v + dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]));
  };
  // 矩阵 → LaTeX
  const mtex = (A, d = 2) => "\\begin{bmatrix}" + A.map(r => r.map(v => f2(v, d)).join("&")).join("\\\\") + "\\end{bmatrix}";

  return { M, C, P, tf2ss, lsim, step, feedback, stepInfo, qp, qpPrep, mpcMatrices, mpcConstraints, rng, COL, refreshColors, Plot, Scene, h, card, slider, seg, check, renderMath, tex, loop, stopAll, tickAll, f2, cfmt, linspace, logspace, clamp, rk4, mtex };
})();
