/* 问题层数据：研究方向 → 控制问题 → 解决方案（目录复用） → 知识点
 * 这些节点直接并入主图谱：方向/问题/方案是图谱的前三列。
 * solutions 目录被多个问题复用；problem.solutions[].ref 指向目录 id，note 是该问题下的具体用法。
 * uses 引用 nodes.js 的知识节点 id。
 */
window.KG_PROBLEMS = {

  /* ── 解决方案目录（17 个，跨问题复用） ── */
  solutions: [
    { id: "s-pid",      name: "级联 PID / 频域整定",       uses: ["pid", "bode", "lead-lag", "transfer-function"],
      desc: "单回路或内快外慢的级联 PID，用频域裕度（Bode/超前滞后）整定。工业与飞控的默认起点：结构简单、不依赖精确模型。" },
    { id: "s-lqr",      name: "LQR / 极点配置状态反馈",     uses: ["lqr", "pole-placement", "linearization", "controllability"],
      desc: "对工作点线性化模型做全状态反馈：要么直接配极点，要么用 LQR 让 Q/R 权衡替你选极点。量产控制器的主力形态。" },
    { id: "s-mpc",      name: "MPC 滚动优化",              uses: ["mpc", "discretization", "dp-hjb"],
      desc: "每拍在线求解带约束的有限时域最优控制，只执行第一步。执行器饱和、安全约束、多目标权衡的系统化答案。" },
    { id: "s-smc",      name: "滑模 / Super-Twisting",     uses: ["sliding-mode", "hosm", "lyapunov"],
      desc: "把误差钉在滑动面上，对匹配不确定性完全免疫；Super-Twisting 用连续控制量消除抖振并保持有限时间收敛。" },
    { id: "s-backstepping", name: "反步法级联设计",         uses: ["backstepping", "iss", "lyapunov"],
      desc: "对『位置→速度→姿态→…』严格反馈链逐级设计虚拟控制并递推构造 Lyapunov 函数，整体稳定性一次证完。" },
    { id: "s-fl",       name: "反馈线性化 / 动态逆 / 微分平坦", uses: ["feedback-linearization", "linearization"],
      desc: "用模型把非线性精确抵消成积分器链（计算力矩、INDI、平坦前馈都是变体），再套线性设计。依赖模型质量。" },
    { id: "s-adaptive", name: "自适应控制",                uses: ["adaptive", "lasalle", "lyapunov"],
      desc: "把未知参数当状态在线学习，自适应律由 Lyapunov 函数反推出来——负载变化、参数漂移场景的正解。" },
    { id: "s-ekf",      name: "EKF / ESKF 多传感融合",     uses: ["ekf", "kalman-filter", "observability", "discretization"],
      desc: "IMU 高频递推 + GNSS/视觉/轮速低频修正，沿估计点线性化套 Kalman 公式。移动机器人与无人机导航的标准架构。" },
    { id: "s-pf",       name: "粒子滤波定位",              uses: ["particle-filter", "ekf"],
      desc: "用粒子云表示多峰置信分布：全局定位、绑架恢复、GPS 拒止环境的地形/指纹匹配。" },
    { id: "s-observer", name: "观测器 / 扰动观测器 / 互补滤波", uses: ["observer", "observability", "transfer-function", "bode"],
      desc: "从有限测量重构状态或把扰动当未知输入估计出来直接补偿（DOB）；互补滤波是它的极简频域形态。" },
    { id: "s-lqg",      name: "Kalman + LQG 最优输出反馈",  uses: ["lqg", "kalman-filter", "lqr", "riccati"],
      desc: "随机噪声下的最优组合：Kalman 滤波估计状态，LQR 拿估计值反馈，分离定理保证各设计各的。" },
    { id: "s-hinf",     name: "H∞ / 不确定性鲁棒设计",      uses: ["hinf", "uncertainty", "small-gain", "mu"],
      desc: "把模型误差建成不确定性包络，对最坏情况优化，带鲁棒稳定证书——性能没那么极致但从不翻车。" },
    { id: "s-lpv",      name: "LMI / LPV 增益调度",        uses: ["lmi", "uncertainty", "linear-stability"],
      desc: "参数（车速、动压、负载）大范围变化时，把一族工作点的设计写成 LMI 凸优化统一求解并保证切换稳定。" },
    { id: "s-smith",    name: "Smith 预估 / 时滞补偿",      uses: ["time-delay", "pid", "nyquist"],
      desc: "延迟吃相位裕度是带宽天花板；用内部模型把延迟搬出回路，恢复无延迟设计的性能。" },
    { id: "s-energy",   name: "能量法 + 切换控制",          uses: ["lyapunov", "lasalle", "linearization"],
      desc: "用能量函数直接当设计工具（摇起倒立摆、无源性整形），进入吸引域后切换给线性控制器接管。" },
    { id: "s-pmp",      name: "最优轨迹 / PMP",            uses: ["pmp", "dp-hjb"],
      desc: "时间最优、能量最优的机动极限由极小值原理给出（bang-bang 结构），离线算轨迹在线跟踪。" },
    { id: "s-impedance", name: "阻抗 / 力位混合控制",       uses: ["transfer-function", "lyapunov", "small-gain", "adaptive"],
      desc: "把机器人端点塑形成期望的质量-弹簧-阻尼，与未知环境接触时用无源性/小增益论证交互稳定。" },
  ],

  /* ── 研究方向 ── */
  domains: [
    { id: "d-vehicle", icon: "🚗", name: "车辆动力学",
      desc: "轮胎-悬架-车身构成的强非线性系统。核心矛盾：模型参数（附着系数、载荷、侧偏刚度）时刻在变，而安全容不得半点闪失。",
      problems: ["p-esc", "p-sideslip", "p-abs", "p-suspension"] },
    { id: "d-ad", icon: "🚙", name: "自动驾驶汽车",
      desc: "感知-规划-控制全栈。控制层两大主题：贴着物理极限精确执行轨迹，以及在噪声与延迟中估计自身状态。",
      problems: ["p-lat", "p-long", "p-loc", "p-limit"] },
    { id: "d-robot", icon: "🦾", name: "机器人",
      desc: "机械臂、移动机器人、腿足平台。共同主题：多刚体非线性动力学 + 参数/环境不确定 + 状态部分可测。",
      problems: ["p-arm", "p-slam", "p-balance", "p-force"] },
    { id: "d-uav", icon: "🚁", name: "无人机",
      desc: "欠驱动、姿态-位置天然级联、强风扰敏感、载荷电量导致参数时变。级联与非线性控制的最佳练兵场。",
      problems: ["p-att", "p-pos", "p-nav", "p-wind", "p-agile"] },
  ],

  /* ── 控制问题 ── */
  problems: [
    /* 车辆动力学 */
    { id: "p-esc", domain: "d-vehicle", title: "横向稳定控制 ESC/DYC",
      desc: "高速过弯或急打方向时车辆甩尾/推头。控制横摆角速度与质心侧偏角贴住期望值，通过差动制动或扭矩分配施加横摆力矩。",
      solutions: [
        { ref: "s-lqr", note: "二自由度自行车模型线性化后做横摆角速度跟踪，期望值由方向盘转角和车速生成" },
        { ref: "s-smc", note: "附着系数、侧偏刚度大幅摄动下仍要稳——滑模对匹配不确定性的不变性正对此症" },
        { ref: "s-mpc", note: "四轮扭矩分配受执行器约束与附着圆约束，滚动优化天然处理" },
      ] },
    { id: "p-sideslip", domain: "d-vehicle", title: "质心侧偏角与轮胎力估计",
      desc: "侧偏角是稳定性控制的核心反馈量，但量产车装不起光学传感器——必须从横摆角速度、侧向加速度里『估』出来。",
      solutions: [
        { ref: "s-ekf", note: "非线性轮胎模型（魔术公式）+ IMU 量测的标准 EKF 配方" },
        { ref: "s-observer", note: "利用轮胎侧偏特性单调性设计带收敛证明的非线性观测器" },
      ] },
    { id: "p-abs", domain: "d-vehicle", title: "防抱死 ABS / 驱动防滑 TCS",
      desc: "把轮胎滑移率控在附着峰值附近（λ≈0.15~0.2）。附着-滑移曲线未知且路面突变，毫秒级时间尺度上强非线性。",
      solutions: [
        { ref: "s-smc", note: "滑移率滑模控制：对路面突变鲁棒 + 有限时间收敛，门限逻辑法的现代化替身" },
        { ref: "s-adaptive", note: "附着峰值位置未知——在线估计并跟踪极值点" },
      ] },
    { id: "p-suspension", domain: "d-vehicle", title: "主动/半主动悬架",
      desc: "舒适性（车身加速度小）与操稳性（轮胎动载荷小）天生矛盾，还受悬架行程约束——典型的多目标频域权衡。",
      solutions: [
        { ref: "s-lqg", note: "四分之一车模型上按二次代价折中舒适与操稳，路面激励视为随机噪声" },
        { ref: "s-hinf", note: "『不同频段抑制不同通道』写成加权灵敏度整形，参数摄动下有保证" },
      ] },

    /* 自动驾驶 */
    { id: "p-lat", domain: "d-ad", title: "轨迹跟踪（横向控制）",
      desc: "把规划层的参考轨迹变成方向盘动作，横向误差压到厘米级，从低速泊车到高速巡航全工况稳定。",
      solutions: [
        { ref: "s-lqr", note: "以横向偏差与航向偏差为状态的误差动力学 + LQR，量产主流（Apollo 同款）" },
        { ref: "s-mpc", note: "预知前方曲率作前馈 + 执行器/侧向加速度约束，大曲率弯道优势明显" },
        { ref: "s-lpv", note: "误差动力学随车速大幅变化——LPV 增益调度给出全速域有保证的设计" },
      ] },
    { id: "p-long", domain: "d-ad", title: "纵向控制 ACC / 队列行驶",
      desc: "油门刹车跟踪期望加速度并保持车距。动力总成有延迟迟滞；队列场景还要误差不沿车队放大（串稳定性）。",
      solutions: [
        { ref: "s-pid", note: "外环车距/车速、内环加速度的分层 PID + 前馈，配抗饱和" },
        { ref: "s-smith", note: "驱动系统 0.2~0.5s 纯延迟直接吃相位裕度，预估器恢复带宽" },
        { ref: "s-mpc", note: "舒适性（加加速度）、安全距离硬约束、燃油最优统一进代价函数" },
      ] },
    { id: "p-loc", domain: "d-ad", title: "定位与状态估计",
      desc: "GNSS 丢星、IMU 漂移、轮速打滑——多源异质信息在噪声中融合出厘米级位姿。",
      solutions: [
        { ref: "s-ekf", note: "IMU 递推 + GNSS/轮速/视觉修正的组合导航标准架构" },
        { ref: "s-pf", note: "初始位置未知或 GNSS 长期失效时的多峰假设定位（激光地图匹配）" },
      ] },
    { id: "p-limit", domain: "d-ad", title: "极限工况运动控制",
      desc: "紧急避障/漂移：轮胎进入饱和区、系统强非线性开环不稳定，常规线性控制直接失效——安全的最后一道防线。",
      solutions: [
        { ref: "s-mpc", note: "带非线性轮胎模型的滚动优化，逼近操稳极限" },
        { ref: "s-fl", note: "已知非线性精确抵消，剩余不确定性交给鲁棒外环" },
        { ref: "s-smc", note: "模型失配与路面突变当匹配扰动一锅端" },
      ] },

    /* 机器人 */
    { id: "p-arm", domain: "d-robot", title: "机械臂轨迹跟踪",
      desc: "多关节强耦合非线性动力学 M(q)q̈+C(q,q̇)q̇+g(q)=τ，高速高精度跟踪，负载一换参数就变。",
      solutions: [
        { ref: "s-fl", note: "计算力矩控制：用模型抵消非线性化成解耦双积分器再配 PD——机器人教科书第一课" },
        { ref: "s-adaptive", note: "抓取未知负载时惯性参数在线学习（Slotine–Li 控制器）" },
        { ref: "s-smc", note: "摩擦、齿隙等难建模项当匹配扰动处理" },
      ] },
    { id: "p-slam", domain: "d-robot", title: "移动机器人定位与建图",
      desc: "里程计漂移 + 地图未知/部分已知。定位是一切自主行为的地基。",
      solutions: [
        { ref: "s-ekf", note: "EKF-SLAM：位姿与路标联合估计，机器人状态估计的起点" },
        { ref: "s-pf", note: "蒙特卡洛定位 MCL：全局定位/绑架恢复，ROS AMCL 即此" },
      ] },
    { id: "p-balance", domain: "d-robot", title: "欠驱动与腿足平衡",
      desc: "轮腿、双足、平衡车的共同内核是倒立摆：开环不稳定 + 欠驱动 + 落脚点/扭矩约束。",
      solutions: [
        { ref: "s-lqr", note: "直立点线性化后的标准配方（本站 LQR 页可直接玩）" },
        { ref: "s-mpc", note: "落脚点、摩擦锥约束下的质心滚动优化，四足/双足现代主流" },
        { ref: "s-energy", note: "能量整形把摆荡到顶，再切换给 LQR 接管" },
      ] },
    { id: "p-force", domain: "d-robot", title: "力控制与物理交互",
      desc: "打磨、装配、人机协作：既控位置又控接触力，环境刚度未知，接触瞬间易失稳。",
      solutions: [
        { ref: "s-impedance", note: "把端点塑形成期望『质量-弹簧-阻尼』，交互稳定性用无源性论证" },
        { ref: "s-observer", note: "无力传感器时用扰动观测器估计接触力" },
      ] },

    /* 无人机 */
    { id: "p-att", domain: "d-uav", title: "姿态控制（内环）",
      desc: "无人机一切能力的地基：角速度/姿态环带宽要做到几十 Hz，大机动下欧拉角奇异、惯量参数只有粗略值。",
      solutions: [
        { ref: "s-pid", note: "角速度内环 + 姿态外环的频率分离设计，PX4/Betaflight 主体结构" },
        { ref: "s-lqr", note: "小角度下姿态解耦为三个双积分器，LQR 一次配平三轴" },
        { ref: "s-smc", note: "惯量不准、桨叶老化、阵风力矩全当匹配扰动，Super-Twisting 不烧电机" },
      ] },
    { id: "p-pos", domain: "d-uav", title: "位置与轨迹跟踪（外环）",
      desc: "位置不能直接驱动——必须先倾斜机身把推力矢量指过去：位置→期望姿态→期望角速度→力矩，典型严格反馈链。",
      solutions: [
        { ref: "s-pid", note: "位置环输出期望加速度，几何分解成期望姿态与总推力，工程主流" },
        { ref: "s-backstepping", note: "沿严格反馈链逐级递推虚拟控制，整体 Lyapunov 证明一气呵成" },
        { ref: "s-fl", note: "微分平坦：由平坦输出解析算出全部状态与输入作前馈，激进轨迹标配" },
      ] },
    { id: "p-nav", domain: "d-uav", title: "状态估计与导航",
      desc: "IMU 高频漂移、GPS 低频丢星、气压计受气流干扰——融合出稳定姿态与位置，姿态还长在旋转流形上。",
      solutions: [
        { ref: "s-ekf", note: "ESKF 误差状态公式回避姿态奇异，PX4 EKF2 即此" },
        { ref: "s-observer", note: "互补滤波：陀螺取高频、加计取低频，算力极小的频域融合" },
        { ref: "s-pf", note: "GPS 拒止时的地形匹配/信号指纹多峰定位" },
      ] },
    { id: "p-wind", domain: "d-uav", title: "抗风扰与鲁棒飞行",
      desc: "阵风是头号敌人：桥下巡检、舰面起降、城市峡谷，风扰可达自重 30% 以上且不可预测。",
      solutions: [
        { ref: "s-observer", note: "扰动观测器把风扰当未知输入实时估计并前馈抵消" },
        { ref: "s-smc", note: "Super-Twisting 只需风扰变化率的界，连续控制量 + 有限时间收敛" },
        { ref: "s-hinf", note: "把风谱与参数摄动写成加权不确定性，带证书的鲁棒设计" },
      ] },
    { id: "p-agile", domain: "d-uav", title: "激进机动与能量最优",
      desc: "穿越机竞速、抢险快达、续航极限：贴着推力饱和与电池极限飞，控制与轨迹优化深度耦合。",
      solutions: [
        { ref: "s-mpc", note: "推力/角速度硬约束下的非线性 MPC，竞速当前最强基线" },
        { ref: "s-pmp", note: "时间最优轨迹：推力恒在约束边界切换（bang-bang），给出机动物理极限" },
        { ref: "s-fl", note: "增量动态逆 INDI：用加速度测量替代模型知识，对模型误差近乎免疫" },
      ] },
  ],
};
