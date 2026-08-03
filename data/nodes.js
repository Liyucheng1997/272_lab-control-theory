/* 知识图谱结构：节点 + 依赖边 + 层级
 * 每个节点: id, title, en, cat(分类→颜色), tier(层级→图谱列位置), tagline
 * 边 [a, b] 表示: a 是 b 的基础(先学 a 再学 b)
 * tier 层层递进：后一层依赖前面各层
 */
window.KG = (function () {
  const categories = {
    math:      { name: "数学与系统基础", color: "#5b8def" },
    stability: { name: "稳定性理论",     color: "#9d6bff" },
    structure: { name: "结构性质与状态反馈", color: "#31b5a9" },
    classical: { name: "经典控制（频域）", color: "#e8a838" },
    optimal:   { name: "最优控制与估计",  color: "#ef5b7b" },
    nonlinear: { name: "非线性控制",     color: "#7bc74d" },
    advanced:  { name: "鲁棒与先进控制",  color: "#e06bd8" },
  };

  const tiers = {
    1: "① 建模基础",
    2: "② 系统分析",
    3: "③ 经典设计",
    4: "④ 状态空间设计",
    5: "⑤ 最优与估计",
    6: "⑥ 非线性控制",
    7: "⑦ 鲁棒与先进",
  };

  const nodes = [
    // ── ① 建模基础 ──
    { id: "state-space",       tier: 1, title: "状态空间模型",        en: "State-Space Model",        cat: "math",      tagline: "用一阶向量微分方程描述一切动态系统" },
    { id: "matrix-exponential",tier: 1, title: "矩阵指数与系统解",     en: "Matrix Exponential",       cat: "math",      tagline: "线性系统的解析解 x(t)=e^{At}x(0)" },
    { id: "transfer-function", tier: 1, title: "传递函数",            en: "Transfer Function",        cat: "math",      tagline: "频域视角：输入输出的拉氏变换之比" },
    { id: "linearization",     tier: 1, title: "线性化",              en: "Linearization",            cat: "math",      tagline: "在平衡点附近用线性系统近似非线性系统" },
    { id: "discretization",    tier: 1, title: "离散化",              en: "Discretization",           cat: "math",      tagline: "把连续系统搬进数字计算机：ZOH 精确离散" },

    // ── ② 系统分析 ──
    { id: "linear-stability",  tier: 2, title: "线性系统稳定性",       en: "Linear Stability",         cat: "stability", tagline: "特征值实部决定一切" },
    { id: "routh-hurwitz",     tier: 2, title: "Routh–Hurwitz 判据",  en: "Routh–Hurwitz Criterion",  cat: "stability", tagline: "不解方程就能判断多项式根的位置" },
    { id: "lyapunov",          tier: 2, title: "Lyapunov 稳定性理论", en: "Lyapunov Stability",       cat: "stability", tagline: "能量下降 ⇒ 稳定，非线性稳定性的基石" },
    { id: "lasalle",           tier: 2, title: "LaSalle 不变原理",    en: "LaSalle Invariance",       cat: "stability", tagline: "V̇≤0 也能推出渐近稳定的利器" },
    { id: "controllability",   tier: 2, title: "能控性",              en: "Controllability",          cat: "structure", tagline: "输入能否把状态开到任意地方" },
    { id: "observability",     tier: 2, title: "能观性",              en: "Observability",            cat: "structure", tagline: "从输出能否反推出内部状态" },

    // ── ③ 经典设计 ──
    { id: "pid",               tier: 3, title: "PID 控制",            en: "PID Control",              cat: "classical", tagline: "工业界 90% 的回路：比例-积分-微分" },
    { id: "bode",              tier: 3, title: "Bode 图与频域指标",    en: "Bode Plot & Margins",      cat: "classical", tagline: "增益裕度、相位裕度：稳定的『安全余量』" },
    { id: "nyquist",           tier: 3, title: "Nyquist 稳定判据",    en: "Nyquist Criterion",        cat: "classical", tagline: "开环曲线绕 -1 的圈数决定闭环稳定" },
    { id: "root-locus",        tier: 3, title: "根轨迹法",            en: "Root Locus",               cat: "classical", tagline: "增益变化时闭环极点画出的轨迹" },
    { id: "lead-lag",          tier: 3, title: "超前–滞后校正",       en: "Lead–Lag Compensation",    cat: "classical", tagline: "用一个零极点对手工整形频率响应" },
    { id: "discrete-control",  tier: 3, title: "离散控制与 z 域",     en: "Discrete-Time Control",    cat: "classical", tagline: "单位圆内稳定、有限拍控制：数字控制器的地盘" },
    { id: "time-delay",        tier: 3, title: "时滞系统与 Smith 预估", en: "Time-Delay & Smith Predictor", cat: "classical", tagline: "延迟吃相位裕度，预估器把它搬出回路" },

    // ── ④ 状态空间设计 ──
    { id: "pole-placement",    tier: 4, title: "极点配置",            en: "Pole Placement",           cat: "structure", tagline: "能控 ⇒ 闭环极点任你摆" },
    { id: "observer",          tier: 4, title: "状态观测器与分离原理",  en: "Observer & Separation",    cat: "structure", tagline: "估计出看不见的状态，再拿去做反馈" },

    // ── ⑤ 最优与估计 ──
    { id: "dp-hjb",            tier: 5, title: "动态规划与 HJB 方程", en: "Dynamic Programming & HJB", cat: "optimal",  tagline: "最优性原理：从终点倒着想" },
    { id: "lqr",               tier: 5, title: "LQR 最优调节器",      en: "Linear Quadratic Regulator", cat: "optimal", tagline: "二次代价下的最优状态反馈 u=-Kx" },
    { id: "riccati",           tier: 5, title: "Riccati 方程",        en: "Riccati Equation",         cat: "optimal",   tagline: "LQR 与 Kalman 滤波共同的心脏" },
    { id: "pmp",               tier: 5, title: "Pontryagin 极小值原理", en: "Pontryagin's Minimum Principle", cat: "optimal", tagline: "带约束最优控制的一阶必要条件" },
    { id: "kalman-filter",     tier: 5, title: "Kalman 滤波器",       en: "Kalman Filter",            cat: "optimal",   tagline: "噪声中的最优状态估计" },
    { id: "lqg",               tier: 5, title: "LQG 与分离定理",      en: "LQG Control",              cat: "optimal",   tagline: "LQR + Kalman = 输出反馈最优控制" },
    { id: "ekf",               tier: 5, title: "扩展 Kalman 滤波 EKF", en: "Extended Kalman Filter",   cat: "optimal",   tagline: "非线性世界的 Kalman：沿估计点线性化" },
    { id: "particle-filter",   tier: 5, title: "粒子滤波",            en: "Particle Filter",          cat: "optimal",   tagline: "用一群粒子表示任意形状的置信分布" },

    // ── ⑥ 非线性控制 ──
    { id: "feedback-linearization", tier: 6, title: "反馈线性化",     en: "Feedback Linearization",   cat: "nonlinear", tagline: "用非线性反馈把系统变成积分器链" },
    { id: "backstepping",      tier: 6, title: "反步法 Backstepping", en: "Backstepping",             cat: "nonlinear", tagline: "逐级递推构造 Lyapunov 函数与控制律" },
    { id: "sliding-mode",      tier: 6, title: "滑模控制",            en: "Sliding Mode Control",     cat: "nonlinear", tagline: "把状态『钉』在滑动面上，对扰动免疫" },
    { id: "hosm",              tier: 6, title: "高阶滑模 Super-Twisting", en: "Higher-Order Sliding Mode", cat: "nonlinear", tagline: "连续控制量也能有限时间滑模——抖振的解药" },
    { id: "adaptive",          tier: 6, title: "自适应控制 (MRAC)",   en: "Adaptive Control",         cat: "nonlinear", tagline: "参数未知？边跑边学" },
    { id: "iss",               tier: 6, title: "输入-状态稳定 ISS",   en: "Input-to-State Stability", cat: "nonlinear", tagline: "有扰动时的稳定性语言：扰动小则偏差小" },

    // ── ⑦ 鲁棒与先进 ──
    { id: "mpc",               tier: 7, title: "模型预测控制 MPC",    en: "Model Predictive Control", cat: "advanced",  tagline: "滚动求解带约束的有限时域最优控制" },
    { id: "small-gain",        tier: 7, title: "小增益定理",          en: "Small-Gain Theorem",       cat: "advanced",  tagline: "回路增益之积 < 1 ⇒ 互联稳定" },
    { id: "uncertainty",       tier: 7, title: "不确定性建模与鲁棒稳定", en: "Uncertainty & Robust Stability", cat: "advanced", tagline: "把『模型不准』本身建成模型：‖WT‖<1 即稳" },
    { id: "lmi",               tier: 7, title: "LMI 与凸优化方法",    en: "Linear Matrix Inequalities", cat: "advanced", tagline: "把控制设计写成凸优化，一键求解" },
    { id: "hinf",              tier: 7, title: "H∞ 鲁棒控制",         en: "H-infinity Control",       cat: "advanced",  tagline: "最坏情况扰动下的最优抑制" },
    { id: "mu",                tier: 7, title: "μ 分析与综合",        en: "Structured Singular Value", cat: "advanced",  tagline: "结构化不确定性的精确鲁棒性度量" },
  ];

  const edges = [
    // ①内部与①→②
    ["state-space", "matrix-exponential"],
    ["state-space", "transfer-function"],
    ["state-space", "linearization"],
    ["matrix-exponential", "discretization"],
    ["matrix-exponential", "linear-stability"],
    ["state-space", "controllability"],
    ["state-space", "observability"],
    ["linear-stability", "routh-hurwitz"],
    ["linear-stability", "lyapunov"],
    ["linearization", "lyapunov"],
    ["lyapunov", "lasalle"],
    // ②→③
    ["transfer-function", "pid"],
    ["transfer-function", "bode"],
    ["transfer-function", "root-locus"],
    ["routh-hurwitz", "root-locus"],
    ["bode", "nyquist"],
    ["bode", "lead-lag"],
    ["root-locus", "lead-lag"],
    ["discretization", "discrete-control"],
    ["linear-stability", "discrete-control"],
    ["bode", "time-delay"],
    ["nyquist", "time-delay"],
    ["pid", "time-delay"],
    // ②→④
    ["controllability", "pole-placement"],
    ["observability", "observer"],
    ["pole-placement", "observer"],
    // →⑤
    ["state-space", "dp-hjb"],
    ["dp-hjb", "lqr"],
    ["dp-hjb", "pmp"],
    ["lyapunov", "lqr"],
    ["controllability", "lqr"],
    ["pole-placement", "lqr"],
    ["lqr", "riccati"],
    ["lqr", "kalman-filter"],
    ["lqr", "lqg"],
    ["observer", "kalman-filter"],
    ["discretization", "kalman-filter"],
    ["kalman-filter", "lqg"],
    ["kalman-filter", "ekf"],
    ["linearization", "ekf"],
    ["ekf", "particle-filter"],
    // →⑥
    ["linearization", "feedback-linearization"],
    ["feedback-linearization", "backstepping"],
    ["lyapunov", "backstepping"],
    ["lyapunov", "sliding-mode"],
    ["sliding-mode", "hosm"],
    ["lyapunov", "hosm"],
    ["lasalle", "adaptive"],
    ["lyapunov", "iss"],
    // →⑦
    ["lqr", "mpc"],
    ["lyapunov", "mpc"],
    ["discrete-control", "mpc"],
    ["nyquist", "small-gain"],
    ["iss", "small-gain"],
    ["small-gain", "uncertainty"],
    ["bode", "uncertainty"],
    ["uncertainty", "hinf"],
    ["riccati", "hinf"],
    ["lyapunov", "lmi"],
    ["linear-stability", "lmi"],
    ["hinf", "mu"],
    ["uncertainty", "mu"],
    ["lmi", "mu"],
  ];

  // 推荐学习顺序：严格按层级递进
  const learningPath = [
    // ① 建模基础
    "state-space", "matrix-exponential", "transfer-function", "linearization", "discretization",
    // ② 系统分析
    "linear-stability", "routh-hurwitz", "lyapunov", "lasalle", "controllability", "observability",
    // ③ 经典设计
    "pid", "bode", "nyquist", "root-locus", "lead-lag", "discrete-control", "time-delay",
    // ④ 状态空间设计
    "pole-placement", "observer",
    // ⑤ 最优与估计
    "dp-hjb", "lqr", "riccati", "pmp", "kalman-filter", "lqg", "ekf", "particle-filter",
    // ⑥ 非线性控制
    "feedback-linearization", "backstepping", "sliding-mode", "hosm", "adaptive", "iss",
    // ⑦ 鲁棒与先进
    "mpc", "small-gain", "uncertainty", "lmi", "hinf", "mu",
  ];

  return { categories, tiers, nodes, edges, learningPath };
})();
