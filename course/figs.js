/* 方块图 / 示意图（SVG） */
window.Figs = (function () {
  const box = (x, y, w, h, t, cls = "") => `<rect class="dg-box ${cls}" x="${x}" y="${y}" width="${w}" height="${h}" rx="8"/><text x="${x + w / 2}" y="${y + h / 2 + 5}" text-anchor="middle" font-size="15" font-style="italic">${t}</text>`;
  const sum = (x, y, signs = ["+", "−"]) => `<circle class="dg-sum" cx="${x}" cy="${y}" r="13"/><line class="dg-line" x1="${x - 9}" y1="${y - 9}" x2="${x + 9}" y2="${y + 9}"/><line class="dg-line" x1="${x - 9}" y1="${y + 9}" x2="${x + 9}" y2="${y - 9}"/><text x="${x - 24}" y="${y - 6}" class="dg-lbl">${signs[0]}</text><text x="${x + 5}" y="${y + 30}" class="dg-lbl">${signs[1]}</text>`;
  const arr = (pts, cls = "") => {
    const d = pts.map((p, i) => (i ? "L" : "M") + p[0] + " " + p[1]).join(" ");
    const [x1, y1] = pts[pts.length - 2], [x2, y2] = pts[pts.length - 1];
    const a = Math.atan2(y2 - y1, x2 - x1), hx = x2 - 9 * Math.cos(a), hy = y2 - 9 * Math.sin(a);
    const p1 = [hx - 4.5 * Math.sin(a), hy + 4.5 * Math.cos(a)], p2 = [hx + 4.5 * Math.sin(a), hy - 4.5 * Math.cos(a)];
    return `<path class="dg-line ${cls}" d="${d}"/><polygon class="dg-arrow" points="${x2},${y2} ${p1} ${p2}"/>`;
  };
  const lbl = (x, y, t, a = "middle", size = 14) => `<text x="${x}" y="${y}" text-anchor="${a}" font-size="${size}" class="dg-lbl" font-style="italic">${t}</text>`;
  const note = (x, y, t, a = "middle") => `<text x="${x}" y="${y}" text-anchor="${a}" font-size="12.5" fill="var(--muted)" style="fill:var(--muted)">${t}</text>`;
  const svg = (w, h, inner) => `<svg viewBox="0 0 ${w} ${h}" width="${w}">${inner}</svg>`;

  return {
    openloop: () => svg(560, 110,
      arr([[20, 55], [130, 55]]) + lbl(70, 45, "R(s)") + note(70, 80, "参考/目标") +
      box(130, 30, 110, 50, "C(s)") + note(185, 100, "控制器") +
      arr([[240, 55], [330, 55]]) + lbl(285, 45, "U(s)") +
      box(330, 30, 110, 50, "G(s)") + note(385, 100, "被控对象") +
      arr([[440, 55], [540, 55]]) + lbl(490, 45, "X(s)") + note(490, 80, "输出")),
    closedloop: () => svg(640, 190,
      arr([[10, 60], [97, 60]]) + lbl(50, 50, "R(s)") + note(50, 85, "参考值") +
      sum(110, 60) +
      arr([[123, 60], [210, 60]]) + lbl(165, 50, "E(s)") + note(165, 85, "误差") +
      box(210, 35, 110, 50, "C(s)") + note(265, 105, "控制器") +
      arr([[320, 60], [400, 60]]) + lbl(360, 50, "U(s)") +
      box(400, 35, 110, 50, "G(s)") + note(455, 105, "被控对象") +
      arr([[510, 60], [625, 60]]) + lbl(590, 50, "X(s)") +
      arr([[560, 60], [560, 150], [110, 150], [110, 73]]) + note(335, 170, "反馈：测量输出，与目标比较")),
    stateFeedback: () => svg(640, 200,
      arr([[10, 60], [100, 60]]) + lbl(50, 50, "v") +
      sum(113, 60, ["+", "−"]) + arr([[126, 60], [220, 60]]) + lbl(170, 50, "u") +
      box(220, 30, 200, 60, "ż = Az + Bu") + arr([[420, 60], [520, 60]]) + lbl(470, 50, "z") +
      box(520, 35, 90, 50, "C") + arr([[610, 60], [636, 60]]) +
      arr([[470, 60], [470, 150], [340, 150]]) + box(230, 125, 110, 50, "K") + arr([[230, 150], [113, 150], [113, 73]]) +
      note(170, 190, "u = −Kz（全状态反馈）")),
    observer: () => svg(660, 250,
      box(200, 20, 200, 56, "真实系统  ż = Az + Bu") + arr([[400, 48], [500, 48]]) + lbl(450, 40, "y = Cz") +
      arr([[60, 48], [200, 48]]) + lbl(110, 40, "u") +
      arr([[120, 48], [120, 180], [200, 180]]) +
      box(200, 150, 200, 60, "观测器  ẑ̇ = Aẑ + Bu + L(y − Cẑ)") +
      arr([[500, 48], [540, 48], [540, 180], [400, 180]]) + note(590, 120, "用测量 y 修正") +
      arr([[300, 210], [300, 240]]) + lbl(320, 236, "ẑ  (估计的全部状态)", "start")),
    mpcLoop: () => svg(660, 180,
      box(20, 40, 210, 70, "预测模型 + 优化") + note(125, 130, "每步解一个 QP，得到未来 N 步输入") +
      arr([[230, 75], [330, 75]]) + lbl(280, 65, "只用 u(k|k)") +
      box(330, 50, 130, 50, "被控对象") + arr([[460, 75], [560, 75]]) + lbl(510, 65, "x(k+1)") +
      arr([[520, 75], [520, 160], [125, 160], [125, 110]]) + note(320, 176, "下一时刻测量新状态，重新优化（滚动时域）")),
    kalman: () => svg(660, 170,
      box(20, 50, 190, 64, "① 预测（先验）") + note(115, 140, "用模型推算：x̂⁻, P⁻") +
      arr([[210, 82], [300, 82]]) +
      box(300, 50, 190, 64, "② 校正（后验）") + note(395, 140, "用测量修正：K, x̂, P") +
      arr([[490, 82], [560, 82]]) + lbl(605, 87, "x̂ₖ") +
      arr([[395, 50], [395, 20], [115, 20], [115, 50]]) + note(255, 14, "k → k+1 循环")),
  };
})();
