/* 内容：鲁棒与先进控制 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ MPC ════════════════ */
window.KG_CONTENT["mpc"] = {
  sim: "sim-mpc",
  sections: [
    {
      icon: "💡", title: "直觉：开车看远方，只踩当下这一脚",
      html: String.raw`
<p>LQR 很优，但有个致命短板：<strong>不会处理约束</strong>——阀门有上下限、温度不能超、罐子不能溢。而约束恰恰是工业利润所在（最优工作点几乎总在约束边界上）。</p>
<p>MPC 的做法像开车：① 看清前方 $N$ 步，解一个带约束的小型最优化，规划整段动作；② <strong>只执行第一步</strong>；③ 下一拍拿到新测量，从头再规划（滚动时域）。反馈由"反复重算"隐式产生，对扰动和模型误差有天然的纠错力。</p>
<p>1970s 从炼油厂起家（DMC），如今从化工到 SpaceX 火箭着陆、四足机器人步态全是它——因为它是<strong>唯一系统化处理多变量 + 硬约束</strong>的主流方法。理论上，它是在数值地逼近《HJB》的解，而《LQR》是它的无约束特例。</p>`
    },
    {
      icon: "📐", title: "标准形式",
      html: String.raw`
<div class="box def"><span class="box-title">线性 MPC 问题（每拍在线求解）</span>
$$\min_{u_0, \dots, u_{N-1}}\ \sum_{k=0}^{N-1}\big(x_k^{\top}Qx_k + u_k^{\top}Ru_k\big) + x_N^{\top}Px_N$$
$$\text{s.t.}\quad x_{k+1} = Ax_k + Bu_k,\quad x_0 = x(t),\quad u_k \in \mathcal{U},\ x_k \in \mathcal{X},\ x_N \in \mathcal{X}_f$$
<p>$P$（终端代价）、$\mathcal{X}_f$（终端集）是稳定性的关键配件。消去状态后这是一个二次规划（QP），毫秒级可解；解得最优序列后仅施加 $u_0^*$，下拍重来。</p></div>`
    },
    {
      icon: "🧾", title: "稳定性定理与证明（值函数当 Lyapunov 函数）",
      html: String.raw`
<p>滚动时域会不会"短视失稳"？（实验中把 $N$ 调小可复现。）以下经典条件杜绝之。</p>
<div class="box thm"><span class="box-title">定理（MPC 渐近稳定性，Mayne et al. 2000 框架）</span>
<p>设终端配件满足：</p>
<ol>
<li>$\mathcal{X}_f$ 内有局部控制律 $\kappa(x)$（典型取无约束 LQR $\kappa = -Kx$）使 $\mathcal{X}_f$ 不变：$x \in \mathcal{X}_f \Rightarrow Ax + B\kappa(x) \in \mathcal{X}_f$，且约束满足 $\kappa(x) \in \mathcal{U},\ \mathcal{X}_f \subseteq \mathcal{X}$；</li>
<li>终端代价是 $\mathcal{X}_f$ 上该局部闭环的 Lyapunov 函数：
$$P \succ 0:\quad (Ax + B\kappa)^{\top}P(Ax + B\kappa) - x^{\top}Px \le -x^{\top}Qx - \kappa^{\top}R\kappa$$
（取 $P$ 为 DARE 解、$\kappa$ 为 LQR 时恰好<strong>取等</strong>。）</li>
</ol>
<p>则对每个初始可行状态，MPC 闭环渐近稳定且始终满足约束（递推可行）。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>记 $V_N(x)$ 为最优值函数，$\{u_0^*, \dots, u_{N-1}^*\}$、$\{x_0^*, \dots, x_N^*\}$ 为当前最优解。闭环执行 $u_0^*$ 后到达 $x^+ = x_1^*$。</p>
<p><strong>第一步：递推可行性。</strong>为 $x^+$ 构造候选序列——把旧解掐头、尾部续上终端控制律：</p>
$$\tilde{u} = \{u_1^*, \dots, u_{N-1}^*,\ \kappa(x_N^*)\}$$
<p>前 $N-1$ 步沿旧轨迹走（约束已验证）；末步由条件 1：$x_N^* \in \mathcal{X}_f$ ⇒ $\kappa(x_N^*) \in \mathcal{U}$ 且新终态 $Ax_N^* + B\kappa(x_N^*) \in \mathcal{X}_f$。故 $\tilde u$ 可行——<strong>可行性一旦成立便永续</strong>。</p>
<p><strong>第二步：值函数下降。</strong>候选序列的代价 = 旧代价 − 掐掉的头 + 换掉的尾：</p>
$$J(\tilde{u};\, x^+) = V_N(x) - \big(x^{\top}Qx + u_0^{*\top}Ru_0^*\big) - x_N^{*\top}Px_N^* + \ell(x_N^*, \kappa) + x_{N+1}^{\top}Px_{N+1}$$
<p>其中 $\ell(x,u) = x^{\top}Qx + u^{\top}Ru$，$x_{N+1} = Ax_N^* + B\kappa(x_N^*)$。由条件 2，最后三项之和 $\le 0$，于是</p>
$$V_N(x^+) \le J(\tilde{u};\, x^+) \le V_N(x) - \ell(x, u_0^*) \le V_N(x) - x^{\top}Qx$$
<p>（第一个不等号：最优值不超过任何可行候选的值。）</p>
<p><strong>第三步：套 Lyapunov。</strong>$V_N \ge 0$，$V_N(0) = 0$，可证在可行域上有正定上下界；每拍严格下降 $\ge x^{\top}Qx$（$Q \succ 0$ 时）。由离散时间 Lyapunov 定理，闭环渐近稳定。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">工程与理论备注</span>
<ul>
<li>证明的精髓一句话：<strong>"最优值函数本身就是 Lyapunov 函数，尾部拼接给出下降证书"</strong>；</li>
<li>实践中常不显式加 $\mathcal{X}_f$，靠 $N$ 够长 + 终端代价 $P$ 取 DARE 解，同样有稳定性（本页实验即此配置）；</li>
<li>扩展：跟踪 MPC、鲁棒 tube-MPC（不确定性下把"管子"当名义轨迹护栏）、非线性 MPC（SQP/RTI 实时迭代）、显式 MPC（离线把 QP 解成分段仿射律）。</li>
</ul></div>`
    },
  ],
  code: {
    title: "从零实现线性 MPC（cvxpy 求解 QP）",
    text: `import numpy as np
import cvxpy as cp
from scipy.linalg import solve_discrete_are

dt = 0.15
A = np.array([[1, dt], [0, 1]])
B = np.array([[dt**2/2], [dt]])
Q = np.diag([4.0, 1.0]); R = np.array([[1.0]])
P = solve_discrete_are(A, B, Q, R)      # 终端代价 = DARE 解
N, umax = 20, 1.0

def mpc_control(x0):
    x = cp.Variable((2, N+1))
    u = cp.Variable((1, N))
    cost, cons = 0, [x[:, 0] == x0]
    for k in range(N):
        cost += cp.quad_form(x[:, k], Q) + cp.quad_form(u[:, k], R)
        cons += [x[:, k+1] == A @ x[:, k] + B @ u[:, k],
                 cp.abs(u[0, k]) <= umax]
    cost += cp.quad_form(x[:, N], P)
    cp.Problem(cp.Minimize(cost), cons).solve()
    return u.value[0, 0]

# 闭环仿真: 从 x=[8,0] 回原点, |u|<=1
x = np.array([8.0, 0.0])
for k in range(80):
    u0 = mpc_control(x)
    x = A @ x + B.flatten() * u0
print("末态:", np.round(x, 4))          # -> [0, 0], 全程不超约束
# pip install cvxpy 后可运行; 或用本页交互实验(投影梯度法内置求解)`
  }
};

/* ════════════════ 小增益定理 ════════════════ */
window.KG_CONTENT["small-gain"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>两个系统首尾相接成环（控制器-对象、机器人-环境、电网-负载），什么时候闭环不会"越吵越凶"？小增益定理给出朴素而深刻的答案：<strong>只要绕环一圈信号总是被衰减（增益乘积 < 1），回声就会衰灭</strong>。</p>
<p>它不需要知道环里是什么——线性、非线性、时变、甚至只知道一个增益上界的黑箱都行。这正是鲁棒控制的哲学：把"不确定的部分"打包成一个只知增益界的块 $\Delta$，只要名义回路对它的增益够小，怎么不确定都稳。《H∞》就是把这个定理变成可计算的设计方法。</p>`
    },
    {
      icon: "📐", title: "系统增益",
      html: String.raw`
<div class="box def"><span class="box-title">定义（$\mathcal{L}_2$ 增益）</span>
<p>信号范数 $\|u\|_2 = \left(\int_0^\infty \|u(t)\|^2dt\right)^{1/2}$。系统（算子）$H$ 的增益</p>
$$\gamma(H) = \sup_{u \ne 0} \frac{\|Hu\|_2}{\|u\|_2}$$
<p>对稳定 LTI 系统可证 $\gamma(G) = \|G\|_{\infty} \triangleq \sup_\omega \bar\sigma\big(G(j\omega)\big)$——频率响应最大奇异值的峰，即 Bode 幅频图的最高点（此即 H∞ 范数名称的来源）。</p></div>`
    },
    {
      icon: "🧾", title: "定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（小增益定理）</span>
<p>反馈互联（$e_1 = r_1 + H_2e_2,\ e_2 = r_2 + H_1e_1$）中两个子系统各自有限增益稳定：$\|H_ie\|_2 \le \gamma_i\|e\|_2 + \beta_i$。若</p>
$$\gamma_1\gamma_2 < 1$$
<p>则闭环有限增益稳定：内部信号 $e_1, e_2$ 的范数被外部输入 $r_1, r_2$ 的范数线性控制。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>对互联方程取范数（先在截断意义下操作以保证量有限，最后取极限；这里写主干）：</p>
$$\|e_1\| \le \|r_1\| + \gamma_2\|e_2\| + \beta_2,\qquad \|e_2\| \le \|r_2\| + \gamma_1\|e_1\| + \beta_1$$
<p>把第二式代入第一式：</p>
$$\|e_1\| \le \|r_1\| + \gamma_2\|r_2\| + \gamma_2\gamma_1\|e_1\| + \gamma_2\beta_1 + \beta_2$$
<p>移项，除以 $(1 - \gamma_1\gamma_2) > 0$：</p>
$$\|e_1\| \le \frac{\|r_1\| + \gamma_2\|r_2\| + \gamma_2\beta_1 + \beta_2}{1 - \gamma_1\gamma_2}$$
<p>$e_2$ 对称同理。所有内部信号有界且随输入线性缩放，即闭环有限增益 $\mathcal{L}_2$ 稳定。<span class="qed">∎</span></p>
<p><em>严格版需用扩展空间 $\mathcal{L}_{2e}$ 与截断算子处理"解是否先验存在于 $\mathcal{L}_2$"的问题，推理结构完全相同。</em></p></div>
<div class="box note"><span class="box-title">与 Nyquist 的关系及保守性</span>
<ul>
<li>LTI 情形下，$\|L\|_\infty < 1$ ⇒ Nyquist 曲线困在单位圆内、碰不到 $-1$——小增益是 Nyquist 判据"最保守但最普适"的特例；</li>
<li>保守性来源：只用增益、扔掉了相位信息。改进方向：无源性定理（用相位/能量）、积分二次约束 IQC（统一框架）；</li>
<li><strong>鲁棒稳定标准型</strong>：不确定性 $\|\Delta\|_\infty \le 1$ 与名义回路 $M$ 成环，鲁棒稳定判据 $\|M\|_\infty < 1$——这就是《H∞》设计目标的直接来源。</li>
</ul></div>`
    },
  ],
  code: {
    title: "数值验证小增益定理的临界性",
    text: `import numpy as np

def l2_gain_lti(num, den):
    """稳定 LTI 的 L2 增益 = Hinf 范数 = sup |G(jw)|"""
    w = np.logspace(-3, 3, 20000)
    s = 1j * w
    return np.abs(np.polyval(num, s) / np.polyval(den, s)).max()

# 两个稳定一阶系统成负反馈环
g1 = l2_gain_lti([2.0], [1, 1])     # 2/(s+1)  增益 2
for k2 in [0.4, 0.49, 0.51, 0.8]:
    g2 = l2_gain_lti([k2], [1, 2])  # k2/(s+2) 增益 k2/2
    prod = g1 * g2
    # 闭环特征式: (s+1)(s+2) + 2*k2
    clroots = np.roots([1, 3, 2 + 2*k2])
    stable = clroots.real.max() < 0
    print(f"gamma1*gamma2 = {prod:.3f}  小增益{'满足 -> 保证稳定' if prod < 1 else '不满足 -> 无结论'}"
          f" | 实际: {'稳定' if stable else '不稳定'}")
# 注意 prod>1 时实际仍可能稳定 —— 小增益只是充分条件(保守)`
  }
};

/* ════════════════ H∞ ════════════════ */
window.KG_CONTENT["hinf"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>LQG 的教训（Doyle 反例）：对"平均情况"最优的设计可能对"最坏情况"极其脆弱。H∞ 控制换一个哲学：<strong>直接对最坏情况优化</strong>。</p>
<p>把扰动想成一个狡猾的对手，它会挑对你伤害最大的频率进攻；$\|T_{zw}\|_\infty$ 恰是"对手能达到的最大能量放大率"（$\mathcal{L}_2$ 增益的峰值）。H∞ 设计 = 压低这个峰 = 修堤坝按最高洪水位修。配合《小增益定理》，压低对不确定块的传函峰值就<strong>买到了鲁棒稳定的证书</strong>。这是 1980s–90s 控制理论的巅峰工程（Zames 提出，DGKF 1989 给出状态空间解），也是 μ 分析、LMI 时代的入口。</p>`
    },
    {
      icon: "📐", title: "标准问题与混合灵敏度",
      html: String.raw`
<div class="box def"><span class="box-title">H∞ 标准问题</span>
<p>广义对象 $P$：外部输入 $w$（扰动/噪声/参考）、控制 $u$、性能输出 $z$、测量 $y$。求控制器 $K$ 使闭环稳定且</p>
$$\|T_{zw}\|_{\infty} = \sup_{\omega}\bar\sigma\big(T_{zw}(j\omega)\big) < \gamma$$
<p>逐步压小 $\gamma$（γ-迭代）逼近最优。</p></div>
<div class="box def"><span class="box-title">混合灵敏度设计（最常用配方）</span>
<p>记灵敏度 $S = (I + GK)^{-1}$（扰动→误差）、补灵敏度 $T = GK(I+GK)^{-1}$（噪声→输出，也是鲁棒性通道），设计目标</p>
$$\left\|\begin{bmatrix} W_1S \\ W_2KS \\ W_3T \end{bmatrix}\right\|_{\infty} < 1$$
<p>权重即"频域说明书"：$W_1$ 大处逼 $S$ 小（低频跟踪好）；$W_3$ 大处逼 $T$ 小（高频抗噪、容忍模型误差）。注意恒等式 $S + T = I$：不可能处处都小——控制设计的根本守恒律（waterbed 效应）。</p></div>`
    },
    {
      icon: "🧾", title: "核心定理",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（有界实引理 Bounded Real Lemma）</span>
<p>稳定系统 $G = C(sI - A)^{-1}B$ 满足 $\|G\|_\infty < \gamma$ $\iff$ 存在 $P \succ 0$：</p>
$$A^{\top}P + PA + \frac{1}{\gamma^2}PBB^{\top}P + C^{\top}C \prec 0$$
<p>（一个 Riccati 不等式 / LMI。）</p></div>
<div class="box proof"><span class="box-title">证明思路（充分性方向，展示机制）</span>
<p>取储能函数 $V = x^{\top}Px$，对任意 $w \in \mathcal{L}_2$ 沿 $\dot x = Ax + Bw,\ z = Cx$ 计算并配方：</p>
$$\dot{V} - \gamma^2\|w\|^2 + \|z\|^2 = x^{\top}(A^{\top}P + PA + C^{\top}C)x + 2x^{\top}PBw - \gamma^2\|w\|^2$$
$$= x^{\top}\Big(A^{\top}P + PA + \tfrac{1}{\gamma^2}PBB^{\top}P + C^{\top}C\Big)x - \gamma^2\Big\|w - \tfrac{1}{\gamma^2}B^{\top}Px\Big\|^2 \le -\epsilon\|x\|^2$$
<p>积分（$x(0)=0$，$V \ge 0$）：$\|z\|_2^2 \le \gamma^2\|w\|_2^2 - V(\infty) - \epsilon\int\|x\|^2 < \gamma^2\|w\|_2^2$，即 $\mathcal{L}_2$ 增益 $< \gamma$。最坏扰动的形状也顺带读出：$w^* = \frac{1}{\gamma^2}B^{\top}Px$（配方项归零处）——鲁棒控制与微分对策（零和博弈）的联系由此可见。必要性经由 Hamiltonian 矩阵谱条件，同《Riccati》一节的机制。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理（DGKF 中央解，陈述）</span>
输出反馈 H∞ 问题有解 $\iff$ 两个 Riccati 方程（控制型 $X_\infty$、滤波型 $Y_\infty$，各含 $\gamma^{-2}$ 项）有稳定化半正定解，且耦合条件 $\rho(X_\infty Y_\infty) < \gamma^2$。控制器结构 = "最坏情况观测器 + 状态反馈"——与 LQG 同构，但每一处都按最坏情况加固；$\gamma \to \infty$ 时退化为 LQG。</div>
<div class="box note"><span class="box-title">版图与延伸</span>
<ul>
<li>结构化不确定性：μ 分析/综合（小增益按块细化，降低保守性）；</li>
<li>LMI 时代：H∞/H2/极点区域约束统一成凸优化（Boyd et al. 1994）；</li>
<li>实践入口：MATLAB <em>mixsyn / hinfsyn</em>；权重函数的选取是设计功力所在。</li>
</ul></div>`
    },
  ],
  code: {
    title: "计算 H∞ 范数（二分法 + Hamiltonian 谱判据）",
    text: `import numpy as np

# ||G||_inf < gamma  <=>  Hamiltonian 无纯虚特征值 (有界实引理谱形式)
def hinf_norm(A, B, C, tol=1e-6):
    lo = 0.0
    hi = 1.0
    # 先扩大上界
    while not no_imag_eig(A, B, C, hi):
        hi *= 2
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if no_imag_eig(A, B, C, mid):
            hi = mid
        else:
            lo = mid
    return hi

def no_imag_eig(A, B, C, gamma):
    H = np.block([[A, (B @ B.T) / gamma**2],
                  [-C.T @ C, -A.T]])
    eigs = np.linalg.eigvals(H)
    return np.min(np.abs(eigs.real)) > 1e-8

# 例: G(s) = 1/(s^2 + 0.2 s + 1)  谐振峰在 w≈1
A = np.array([[0., 1.], [-1., -0.2]])
B = np.array([[0.], [1.]])
C = np.array([[1., 0.]])

gam = hinf_norm(A, B, C)
# 对照: 直接扫频
w = np.logspace(-2, 2, 40000)
G = np.abs(1 / ((1j*w)**2 + 0.2*(1j*w) + 1))
print(f"二分+Hamiltonian: {gam:.5f}")
print(f"扫频峰值        : {G.max():.5f}")
print(f"理论 1/(0.2*sqrt(1-0.01)) ≈ {1/(0.2*np.sqrt(0.99)):.5f}")`
  }
};
