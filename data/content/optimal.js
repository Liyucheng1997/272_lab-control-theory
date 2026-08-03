/* 内容：最优控制与估计 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ 动态规划与 HJB ════════════════ */
window.KG_CONTENT["dp-hjb"] = {
  sections: [
    {
      icon: "💡", title: "直觉：从终点倒着想",
      html: String.raw`
<p>Bellman（1950s）的最优性原理一句话：<strong>最优路径的任何一段尾巴，本身也是（从该点出发的）最优路径</strong>。若不然，把尾巴换成更优的，整条路就更好——矛盾。</p>
<p>于是可以从终点<strong>倒着</strong>递推"从每个状态出发的最优代价"（值函数 $V$）。连续时间下，这个递推取极限就成为一个偏微分方程——HJB 方程。它是最优控制的"万有理论"：LQR 是它的二次特例，MPC 在数值上逼近它，强化学习的 Bellman 方程是它的采样版本。</p>`
    },
    {
      icon: "📐", title: "问题设定与最优性原理",
      html: String.raw`
<div class="box def"><span class="box-title">最优控制问题</span>
$$V(x_0) = \min_{u(\cdot)} \int_0^{\infty} \ell(x(t), u(t))\,dt,\qquad \dot{x} = f(x, u),\ x(0) = x_0$$
$\ell \ge 0$ 为运行代价（stage cost），$V$ 称为<strong>值函数</strong>（最优代价函数）。</div>
<div class="box thm"><span class="box-title">最优性原理</span>
对任意 $\Delta t > 0$：
$$V(x_0) = \min_{u|_{[0,\Delta t]}} \left\{ \int_0^{\Delta t} \ell(x, u)\,dt + V(x(\Delta t)) \right\}$$</div>
<div class="box proof"><span class="box-title">证明</span>
<p>记右端为 $W(x_0)$。<strong>($V \ge W$)</strong>：任取一条完整策略 $u(\cdot)$，其代价 $= \int_0^{\Delta t}\ell + \big[\text{尾段代价}\big] \ge \int_0^{\Delta t}\ell + V(x(\Delta t)) \ge W(x_0)$（尾段代价不小于从 $x(\Delta t)$ 出发的最优值）。对所有策略取下确界得 $V \ge W$。</p>
<p><strong>($V \le W$)</strong>：任取头段 $u|_{[0,\Delta t]}$ 与 $\varepsilon > 0$，再取从 $x(\Delta t)$ 出发的 $\varepsilon$-最优尾段拼接成完整策略，其代价 $\le \int_0^{\Delta t}\ell + V(x(\Delta t)) + \varepsilon$。故 $V \le W + \varepsilon$，$\varepsilon$ 任意。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "HJB 方程与验证定理",
      html: String.raw`
<div class="box thm"><span class="box-title">HJB 方程（推导）</span>
设 $V$ 连续可微。在最优性原理中令 $\Delta t \to 0$：$V(x(\Delta t)) \approx V(x_0) + \nabla V^{\top} f(x_0, u)\,\Delta t$，$\int_0^{\Delta t}\ell \approx \ell(x_0, u)\Delta t$，代入并消去 $V(x_0)$、除以 $\Delta t$：
$$\boxed{\ 0 = \min_{u}\ \big[\ \ell(x, u) + \nabla V(x)^{\top} f(x, u)\ \big]\ }$$
最优控制律即逐点取 argmin：$u^*(x) = \arg\min_u [\ell + \nabla V^{\top}f]$——<strong>状态反馈</strong>形式自动涌现。</div>
<div class="box thm"><span class="box-title">验证定理（Verification Theorem）</span>
设 $W \in C^1$，$W \ge 0$，$W(0) = 0$，满足 HJB 方程，且反馈 $u^*(x)$ 使闭环轨迹 $x^*(t) \to 0$ 且 $W(x^*(t)) \to 0$。则 $W$ 就是值函数，且 $u^*$ 最优。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>任取容许控制 $u(\cdot)$（相应轨迹 $x(t)$，设 $W(x(t))\to 0$，代价有限）。沿轨迹：</p>
$$\frac{d}{dt}W(x(t)) = \nabla W^{\top} f(x, u) \ge -\ell(x, u)$$
<p>（由 HJB：$\ell(x,u) + \nabla W^{\top}f(x,u) \ge \min_v[\cdots] = 0$。）从 $0$ 到 $T$ 积分：</p>
$$W(x(T)) - W(x_0) \ge -\int_0^T \ell\,dt \ \Longrightarrow\ \int_0^T \ell\,dt \ge W(x_0) - W(x(T))$$
<p>令 $T \to \infty$：总代价 $\ge W(x_0)$。即任何策略的代价都不低于 $W(x_0)$。</p>
<p>而沿 $u^*$，HJB 的 min 恰好取到，微分不等式变等式，同样积分得代价 $= W(x_0)$。达到下界即最优。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">维数灾难与出路</span>
HJB 是 $n$ 维状态空间上的非线性 PDE，一般解不动（"维数灾难"，Bellman 本人命名）。三条活路：① 结构特殊解析解——LQR（下一站）；② 在线数值化、只解当前状态附近——MPC；③ 采样近似值函数——强化学习/ADP。另外 $V$ 常不可微，严格理论需"粘性解"（Crandall–Lions, 1983）。</div>`
    },
  ],
  code: {
    title: "值迭代解一维 HJB（离散化验证 LQR 解析解）",
    text: `import numpy as np

# 问题: xdot = u, cost = ∫ x^2 + u^2 dt
# 解析解: V(x) = x^2 (Riccati: p^2 = 1), u* = -x
dx, du, dt = 0.05, 0.05, 0.02
xs = np.arange(-3, 3 + dx, dx)
us = np.arange(-4, 4 + du, du)
V = np.zeros_like(xs)

for it in range(3000):
    Vnew = np.empty_like(V)
    for i, x in enumerate(xs):
        best = np.inf
        for u in us:
            xn = np.clip(x + u * dt, xs[0], xs[-1])
            Vn = np.interp(xn, xs, V)
            best = min(best, (x*x + u*u) * dt + Vn)
        Vnew[i] = best
    if np.max(abs(Vnew - V)) < 1e-9:
        break
    V = Vnew

print("迭代次数:", it)
for x in [0.5, 1.0, 2.0]:
    print(f"x={x}: V数值={np.interp(x, xs, V):.3f}  解析 x^2={x*x:.3f}")`
  }
};

/* ════════════════ LQR ════════════════ */
window.KG_CONTENT["lqr"] = {
  sim: "sim-lqr",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>极点配置让你"任意"选极点，但该选哪里？LQR 换一种问法：给状态偏差和控制消耗各标个价——</p>
$$J = \int_0^{\infty} \big(\underbrace{x^{\top}Qx}_{\text{偏差罚款}} + \underbrace{u^{\top}Ru}_{\text{能耗罚款}}\big)\,dt$$
<p>然后让数学替你找出<strong>总账最便宜</strong>的控制。答案美得出奇：最优控制是线性状态反馈 $u = -Kx$，增益由一个矩阵二次方程（Riccati 方程）给出，且闭环<strong>自动稳定</strong>、自带优良鲁棒性。调 $Q, R$ 比调极点直观得多：想让角度收敛快，加大对应的 $Q$ 权重；嫌控制太猛，加大 $R$——右侧倒立摆实验直接体验。</p>`
    },
    {
      icon: "🧾", title: "完整推导（经由 HJB）",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（LQR 主定理）</span>
<p>设 $\dot x = Ax + Bu$，$Q \succeq 0$, $R \succ 0$，$(A, B)$ 能镇定，$(A, Q^{1/2})$ 能检测。则代数 Riccati 方程 (CARE)</p>
$$A^{\top}P + PA - PBR^{-1}B^{\top}P + Q = 0$$
<p>存在唯一对称半正定的<strong>镇定解</strong> $P$，最优控制为</p>
$$u^* = -Kx,\qquad K = R^{-1}B^{\top}P$$
<p>最优代价 $J^* = x_0^{\top}Px_0$，且闭环 $A - BK$ 是 Hurwitz 的。</p></div>
<div class="box proof"><span class="box-title">第一步：HJB + 二次猜测 ⇒ Riccati 方程与最优增益</span>
<p>HJB：$0 = \min_u [x^{\top}Qx + u^{\top}Ru + \nabla V^{\top}(Ax + Bu)]$。方括号对 $u$ 是严格凸二次型（$R \succ 0$），求驻点：</p>
$$2Ru + B^{\top}\nabla V = 0 \ \Rightarrow\ u^* = -\tfrac{1}{2}R^{-1}B^{\top}\nabla V$$
<p>猜 $V(x) = x^{\top}Px$（$P$ 对称待定），则 $\nabla V = 2Px$，$u^* = -R^{-1}B^{\top}Px$。代回 HJB：</p>
$$0 = x^{\top}Qx + x^{\top}PBR^{-1}\underbrace{RR^{-1}}_{I}B^{\top}Px + 2x^{\top}P(Ax - BR^{-1}B^{\top}Px)$$
<p>用 $2x^{\top}PAx = x^{\top}(A^{\top}P + PA)x$ 对称化并整理：</p>
$$0 = x^{\top}\big[Q + A^{\top}P + PA - PBR^{-1}B^{\top}P\big]x,\quad \forall x$$
<p>对一切 $x$ 成立且括号内对称 ⇒ 括号 $= 0$，即 CARE。✓</p></div>
<div class="box proof"><span class="box-title">第二步：闭环稳定性（Lyapunov 论证）</span>
<p>取 $V(x) = x^{\top}Px$ 作候选 Lyapunov 函数（$P \succeq 0$ 为镇定解）。沿闭环 $\dot{x} = (A - BK)x$：</p>
$$\dot{V} = x^{\top}\big[(A - BK)^{\top}P + P(A - BK)\big]x$$
<p>展开并用 CARE 替换 $A^{\top}P + PA = PBR^{-1}B^{\top}P - Q$：</p>
$$(A-BK)^{\top}P + P(A-BK) = A^{\top}P + PA - K^{\top}B^{\top}P - PBK = -Q - PBR^{-1}B^{\top}P$$
<p>（代入 $K = R^{-1}B^{\top}P$，两个交叉项各贡献 $-PBR^{-1}B^{\top}P$，与 CARE 的正项合并。）于是</p>
$$\dot{V} = -x^{\top}Qx - u^{*\top}Ru^* \le 0$$
<p>由 LaSalle/可检测性论证：$\dot V \equiv 0$ ⇒ $Q^{1/2}x \equiv 0$ 且 $u \equiv 0$，此时闭环轨迹是开环轨迹且输出 $Q^{1/2}x$ 恒零，可检测性迫使 $x \to 0$。故闭环渐近稳定（进而可证指数稳定）。✓</p></div>
<div class="box proof"><span class="box-title">第三步：最优性（完成平方，独立于 HJB 的初等证明）</span>
<p>对任意使 $x(t) \to 0$ 的容许控制 $u$，计算</p>
$$\frac{d}{dt}(x^{\top}Px) = (Ax+Bu)^{\top}Px + x^{\top}P(Ax+Bu) = x^{\top}(A^{\top}P + PA)x + 2u^{\top}B^{\top}Px$$
<p>用 CARE 替换 $A^{\top}P+PA$，并配方：</p>
$$\frac{d}{dt}(x^{\top}Px) = -x^{\top}Qx - u^{\top}Ru + (u + R^{-1}B^{\top}Px)^{\top}R\,(u + R^{-1}B^{\top}Px)$$
<p>（展开右端最后一项即验证恒等式。）两边从 $0$ 积到 $\infty$（$x(\infty)=0$）移项：</p>
$$J(u) = \int_0^\infty (x^{\top}Qx + u^{\top}Ru)\,dt = x_0^{\top}Px_0 + \int_0^{\infty} \|u + Kx\|_R^2\,dt$$
<p>第二项非负，且当且仅当 $u = -Kx$ 时为零。故 $u^*$ 唯一最优，$J^* = x_0^{\top}Px_0$。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "📐", title: "鲁棒性与调参",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（LQR 的保证裕度，单输入）</span>
LQR 回路 $L(s) = K(sI - A)^{-1}B$ 满足 Kalman 不等式 $|1 + L(j\omega)| \ge 1$，因此增益裕度 $\in [\frac12, \infty)$，相位裕度 $\ge 60°$。</div>
<div class="box proof"><span class="box-title">证明思路</span>
<p>由 CARE 可推出返回差恒等式（把 CARE 夹在 $B^{\top}(-j\omega I - A^{\top})^{-1}$ 与 $(j\omega I - A)^{-1}B$ 之间展开）：</p>
$$(1 + L(-j\omega))^{\top} R\,(1 + L(j\omega)) = R + B^{\top}(-j\omega I - A^{\top})^{-1} Q\,(j\omega I - A)^{-1}B \succeq R$$
<p>单输入取 $R = r$：$|1 + L(j\omega)|^2 \ge 1$。Nyquist 曲线与 $-1$ 点距离 $\ge 1$，由几何关系直接读出上述裕度。<span class="qed">∎</span></p>
<p><em>警告：这些裕度在接上观测器（LQG）后不再保证——见《LQG》。</em></p></div>
<div class="box note"><span class="box-title">调参经验（Bryson 法则）</span>
$Q, R$ 取对角，$Q_{ii} = 1/(\text{状态}_i\text{最大可接受值})^2$，$R_{jj} = 1/(\text{输入}_j\text{最大可用值})^2$，再整体缩放 $R$ 调"松紧"。实验里试试：$q_\theta$ 加大 → 摆角恢复更快；$R$ 加大 → 控制力变温柔、恢复变慢。</div>`
    },
  ],
  code: {
    title: "LQR 全流程：解 Riccati、验证稳定性与最优性",
    text: `import numpy as np
from scipy.linalg import solve_continuous_are

# 倒立摆小车 (线性化, theta 从竖直向上)
M, m, l, g = 1.0, 0.2, 0.5, 9.81
A = np.array([[0, 1, 0, 0],
              [0, 0, -m*g/M, 0],
              [0, 0, 0, 1],
              [0, 0, (M+m)*g/(M*l), 0]])
B = np.array([[0], [1/M], [0], [-1/(M*l)]])
Q = np.diag([10.0, 1.0, 100.0, 1.0])
R = np.array([[0.1]])

P = solve_continuous_are(A, B, Q, R)
K = np.linalg.solve(R, B.T @ P)
print("K =", np.round(K, 2))

# 1) 验证 CARE 残差
res = A.T@P + P@A - P@B@np.linalg.solve(R, B.T@P) + Q
print("CARE 残差:", np.abs(res).max())

# 2) 闭环稳定
print("闭环特征值实部:", np.round(np.linalg.eigvals(A - B@K).real, 3))

# 3) 数值验证最优性: J(u*) vs 扰动过的增益
def cost(Kt, x0, T=20.0, h=1e-3):
    x, J = x0.copy(), 0.0
    for _ in range(int(T/h)):
        u = -Kt @ x
        J += (x @ Q @ x + u @ R @ u) * h
        x = x + h * (A @ x + (B @ u.reshape(1)).flatten())
    return J

x0 = np.array([0, 0, 0.2, 0])
J_star = cost(K, x0)
print(f"J(u*) = {J_star:.4f},  理论 x0'Px0 = {x0 @ P @ x0:.4f}")
for eps in [0.9, 1.1]:
    print(f"扰动增益 {eps}K: J = {cost(eps*K, x0):.4f}  (应更大)")`
  }
};

/* ════════════════ Riccati ════════════════ */
window.KG_CONTENT["riccati"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>LQR、Kalman 滤波、H∞、微分对策……最优控制的每个角落都会撞见同一个方程：Riccati 方程（源于 18 世纪 Jacopo Riccati 研究的标量二次 ODE）。它是"二次值函数假设"在线性动力学下的必然产物：值函数二次 ⇒ 其系数矩阵满足二次矩阵方程。</p>
<p>值得单独立节，因为它的<strong>解的结构</strong>本身就是一门学问：解不唯一（二次方程！），但"好解"（镇定解）唯一；求解的现代方法不是迭代猜测，而是一次特征分解——优雅得像变魔术。</p>`
    },
    {
      icon: "📐", title: "三种形态",
      html: String.raw`
<div class="box def"><span class="box-title">微分 / 代数 / 离散 Riccati</span>
<p><strong>DRE</strong>（有限时域 LQR，终端罚 $P(T) = S$，倒着积分）：</p>
$$-\dot{P} = A^{\top}P + PA - PBR^{-1}B^{\top}P + Q$$
<p><strong>CARE</strong>（无限时域稳态）：$A^{\top}P + PA - PBR^{-1}B^{\top}P + Q = 0$</p>
<p><strong>DARE</strong>（离散时间，MPC/Kalman 用）：</p>
$$P = A^{\top}PA - A^{\top}PB(R + B^{\top}PB)^{-1}B^{\top}PA + Q$$
<p>时域无限化时 DRE 的解收敛到 CARE 的镇定解（能镇定+能检测下），这正是页面各仿真里用"积分 DRE 至收敛"求 CARE 的原理。</p></div>`
    },
    {
      icon: "🧾", title: "Hamiltonian 矩阵解法（定理与证明）",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（不变子空间法）</span>
<p>定义 $2n \times 2n$ Hamiltonian 矩阵</p>
$$H = \begin{bmatrix} A & -BR^{-1}B^{\top} \\ -Q & -A^{\top} \end{bmatrix}$$
<ol>
<li>$H$ 的谱关于虚轴对称（$\lambda$ 是特征值则 $-\lambda$ 也是）；</li>
<li>能镇定+能检测时 $H$ 无虚轴特征值，恰有 $n$ 个稳定特征值；取其稳定不变子空间的基 $\begin{bmatrix} X \\ Y \end{bmatrix}$（$X, Y \in \mathbb{R}^{n\times n}$），则 $X$ 可逆且 $P = YX^{-1}$ 是 CARE 的唯一镇定解。</li>
</ol></div>
<div class="box proof"><span class="box-title">证明（要点完整）</span>
<p><strong>(1) 谱对称。</strong>取 $J = \begin{bmatrix} 0 & I \\ -I & 0 \end{bmatrix}$，直接验证 $J^{-1}HJ = -H^{\top}$：</p>
$$J^{-1}HJ = \begin{bmatrix} -A^{\top} & Q \\ BR^{-1}B^{\top} & A \end{bmatrix}^{\top\text{-check}} = -H^{\top}$$
<p>相似 ⇒ $\mathrm{spec}(H) = \mathrm{spec}(-H^{\top}) = -\mathrm{spec}(H)$。✓</p>
<p><strong>(2) 子空间 ⇒ Riccati 解。</strong>稳定不变子空间满足 $H\begin{bmatrix} X \\ Y\end{bmatrix} = \begin{bmatrix} X \\ Y \end{bmatrix}\Lambda$（$\Lambda$ 稳定）。展开两行：</p>
$$AX - BR^{-1}B^{\top}Y = X\Lambda \qquad (i)$$
$$-QX - A^{\top}Y = Y\Lambda \qquad (ii)$$
<p>设 $X$ 可逆（可由能镇定性证明，略），令 $P = YX^{-1}$。用 $(i)$：$\Lambda' \triangleq X\Lambda X^{-1} = A - BR^{-1}B^{\top}P$。将 $(ii)$ 右乘 $X^{-1}$：</p>
$$-Q - A^{\top}P = PX\Lambda X^{-1} = P(A - BR^{-1}B^{\top}P)$$
<p>整理即 CARE。且闭环矩阵 $A - BR^{-1}B^{\top}P = X\Lambda X^{-1}$ 与稳定的 $\Lambda$ 相似 ⇒ $P$ 是镇定解。✓</p>
<p><strong>(2') 反向。</strong>若 $P$ 解 CARE，直接验证 $H\begin{bmatrix} I \\ P \end{bmatrix} = \begin{bmatrix} I \\ P\end{bmatrix}(A - BR^{-1}B^{\top}P)$：第一行是恒等式，第二行恰为 CARE。即每个解对应一个 $n$ 维不变子空间；镇定解对应稳定子空间，而稳定子空间唯一 ⇒ 镇定解唯一。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">数值实践</span>
工业级求解器（MATLAB <em>icare</em> / SciPy <em>solve_continuous_are</em>）用的正是此定理的 Schur 分解稳健版（Laub 1979）。手写迭代（Kleinman 牛顿法、DRE 积分）适合教学与小规模在线计算——本站仿真用的即 DRE 法。</div>`
    },
  ],
  code: {
    title: "三种方法解 CARE 并互相验证",
    text: `import numpy as np
from scipy.linalg import schur, solve_continuous_are

A = np.array([[0., 1.], [2., -1.]])
B = np.array([[0.], [1.]])
Q = np.diag([4.0, 1.0]); R = np.array([[1.0]])

# 法1: scipy
P1 = solve_continuous_are(A, B, Q, R)

# 法2: Hamiltonian + 有序 Schur 分解
H = np.block([[A, -B @ np.linalg.solve(R, B.T)],
              [-Q, -A.T]])
T, Z, sdim = schur(H, sort=lambda x: x.real < 0, output='real')
n = 2
X, Y = Z[:n, :n], Z[n:, :n]
P2 = Y @ np.linalg.inv(X)

# 法3: 积分微分 Riccati 方程至稳态
P3 = Q.copy().astype(float)
BRB = B @ np.linalg.solve(R, B.T)
for _ in range(200000):
    dP = A.T@P3 + P3@A - P3@BRB@P3 + Q
    P3 += 1e-4 * dP
    if np.abs(dP).max() < 1e-12:
        break

for name, P in [("scipy", P1), ("Hamiltonian", P2), ("DRE", P3)]:
    res = A.T@P + P@A - P@BRB@P + Q
    cl = np.linalg.eigvals(A - BRB @ P).real.max()
    print(f"{name:12s} 残差 {np.abs(res).max():.2e}  闭环max Re {cl:.3f}")`
  }
};

/* ════════════════ PMP ════════════════ */
window.KG_CONTENT["pmp"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>HJB 要解全空间的 PDE，太贵；很多问题只想要<strong>一条</strong>最优轨迹。Pontryagin（1956，苏联学派）给出沿单条轨迹的一阶必要条件——变分法的"带约束终极版"。</p>
<p>核心角色是<strong>协态</strong>（costate）$\lambda(t)$：它是"此刻状态多偏一点，总代价会变多少"的影子价格（$\lambda = \nabla_x V$，与 HJB 打通）。PMP 说：最优控制在每个瞬间都<strong>极小化 Hamiltonian</strong>——一个把"当前代价"和"对未来的影响（用影子价格折算）"加在一起的量。它对控制约束（$|u| \le 1$ 等）毫不畏惧，直接给出 bang-bang 这类非光滑答案，这是它相对 HJB/LQR 的独门本领。</p>`
    },
    {
      icon: "📐", title: "定理陈述",
      html: String.raw`
<div class="box thm"><span class="box-title">Pontryagin 极小值原理</span>
<p>问题：$\min_u J = \varphi(x(T)) + \int_0^T \ell(x, u)\,dt$，$\dot{x} = f(x, u)$，$u(t) \in U$（任意约束集）。定义 Hamiltonian</p>
$$H(x, u, \lambda) = \ell(x, u) + \lambda^{\top} f(x, u)$$
<p>若 $(x^*, u^*)$ 最优，则存在协态 $\lambda(t)$ 使得：</p>
<ol>
<li><strong>状态方程</strong>：$\dot{x}^* = \partial H/\partial \lambda = f(x^*, u^*)$；</li>
<li><strong>协态方程</strong>：$\dot{\lambda} = -\partial H/\partial x = -\ell_x - f_x^{\top}\lambda$；</li>
<li><strong>极小值条件</strong>：$u^*(t) = \arg\min_{v \in U} H(x^*(t), v, \lambda(t))$；</li>
<li><strong>横截条件</strong>：$\lambda(T) = \nabla\varphi(x^*(T))$（终态自由时）。</li>
</ol></div>
<div class="box proof"><span class="box-title">证明梗概（光滑情形的变分推导）</span>
<p>构造增广泛函（用 $\lambda$ 作动力学约束的乘子）：</p>
$$\bar{J} = \varphi(x(T)) + \int_0^T \big[H(x, u, \lambda) - \lambda^{\top}\dot{x}\big]\,dt$$
<p>对 $-\lambda^{\top}\dot x$ 分部积分：$-[\lambda^{\top}x]_0^T + \int \dot\lambda^{\top}x\,dt$。做变分 $u \to u + \delta u$（引起 $x \to x + \delta x$，$\delta x(0) = 0$）：</p>
$$\delta \bar{J} = \big(\nabla\varphi(x(T)) - \lambda(T)\big)^{\top}\delta x(T) + \int_0^T \Big[\big(H_x + \dot{\lambda}\big)^{\top}\delta x + H_u^{\top}\delta u\Big]dt$$
<p>选 $\lambda$ 满足协态方程与横截条件，前两团消失，剩 $\delta\bar J = \int H_u^{\top}\delta u\,dt$。最优性要求对一切容许 $\delta u$ 非负 ⇒ 无约束时 $H_u = 0$；有约束集 $U$ 时，针状变分（needle variation：在小区间上把 $u$ 换成任意 $v \in U$）给出更强的逐点条件 $H(x^*, v, \lambda) \ge H(x^*, u^*, \lambda)$，即极小值条件。完整严格证明（含非光滑情形）见 Pontryagin 原著。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🔬", title: "经典应用：双积分器时间最优（bang-bang）",
      html: String.raw`
<div class="box eg"><span class="box-title">问题与完整求解</span>
<p>小车 $\ddot{y} = u,\ |u| \le 1$，最短时间从 $(y_0, v_0)$ 到原点。$\ell = 1$（计时），$H = 1 + \lambda_1 v + \lambda_2 u$。</p>
<p><strong>协态</strong>：$\dot\lambda_1 = -H_y = 0 \Rightarrow \lambda_1 \equiv c_1$；$\dot\lambda_2 = -H_v = -\lambda_1 \Rightarrow \lambda_2(t) = c_2 - c_1 t$（线性！）。</p>
<p><strong>极小值条件</strong>：$u^* = \arg\min_{|u|\le1} \lambda_2 u = -\mathrm{sign}(\lambda_2)$。$\lambda_2$ 线性至多变号一次 ⇒ <strong>最优控制是 $\pm 1$ 之间至多切换一次的 bang-bang</strong>。</p>
<p><strong>几何</strong>：$u = \pm 1$ 下轨迹是抛物线 $y = \pm\frac{v^2}{2} + C$。经过原点的两条半抛物线拼成<strong>切换曲线</strong> $y = -\frac{v|v|}{2}$：初始在曲线上方先全刹（$u=-1$）滑到曲线，再全油（$u=+1$）冲进原点（下方反之）。这就是电梯、硬盘磁头、卫星姿态机动的"最速"方案。</p></div>
<div class="box note"><span class="box-title">PMP vs HJB</span>
PMP：两点边值问题（ODE），一条轨迹，开环解，数值上用打靶法/配点法（现代轨迹优化如火箭着陆的核心）；HJB：PDE，全空间反馈律。$\lambda(t) = \nabla_x V(x^*(t))$ 把两者缝合：协态就是值函数沿最优轨迹的梯度。</div>`
    },
  ],
  code: {
    title: "bang-bang 时间最优控制（切换曲线反馈实现）",
    text: `import numpy as np

def bang_bang(y, v):
    """双积分器时间最优反馈: 切换曲线 y = -v|v|/2"""
    s = y + v * abs(v) / 2.0
    if abs(s) < 1e-9 and abs(v) < 1e-9:
        return 0.0
    if abs(s) < 1e-9:          # 已在切换曲线上: 沿曲线滑向原点
        return -np.sign(v)
    return -np.sign(s)

h = 0.001
x = np.array([3.0, 0.0])
t, switches, u_prev = 0.0, 0, None
while np.linalg.norm(x) > 1e-2 and t < 20:
    u = bang_bang(*x)
    if u_prev is not None and u != u_prev:
        switches += 1
    u_prev = u
    x = x + h * np.array([x[1], u])
    t += h

print(f"到达原点用时 {t:.3f}s, 切换 {switches} 次")
# 理论: 从 (3,0) 出发 T* = 2*sqrt(3) ≈ 3.464s, 切换 1 次
print(f"理论最短时间 = 2*sqrt(3) = {2*np.sqrt(3):.3f}s")`
  }
};

/* ════════════════ Kalman 滤波 ════════════════ */
window.KG_CONTENT["kalman-filter"] = {
  sim: "sim-kalman",
  sections: [
    {
      icon: "💡", title: "直觉：预测 + 修正的贝叶斯循环",
      html: String.raw`
<p>GPS 有噪声、陀螺会漂移、模型不精确——但飞船还是得知道自己在哪。Kalman（1960）把估计问题变成两拍循环：</p>
<ul>
<li><strong>预测</strong>：用模型把上一时刻的估计外推一步（不确定性变大）；</li>
<li><strong>修正</strong>：新测量到来，按<strong>可信度加权</strong>折中预测与测量（不确定性变小）。</li>
</ul>
<p>权重（Kalman 增益）不是拍脑袋：它逐时刻精确计算"模型 vs 测量谁更可信"，使估计误差方差最小。阿波罗导航、GPS、无人机融合、目标跟踪……它可能是二十世纪被实现次数最多的算法。它还是 LQR 的<strong>对偶</strong>：同一个 Riccati 方程，把 $(A, B, Q, R)$ 换成 $(A^{\top}, C^{\top}, W, V)$。</p>`
    },
    {
      icon: "📐", title: "问题设定与算法",
      html: String.raw`
<div class="box def"><span class="box-title">模型（离散线性高斯）</span>
$$x_{k+1} = Ax_k + Bu_k + w_k,\qquad y_k = Cx_k + v_k$$
$w_k \sim \mathcal{N}(0, W)$、$v_k \sim \mathcal{N}(0, V)$ 独立白噪声。记 $\hat{x}_{k|k-1}$（先验）与 $\hat{x}_{k|k}$（后验），误差协方差 $P_{k|k-1}, P_{k|k}$。</div>
<div class="box def"><span class="box-title">Kalman 滤波算法</span>
<p><strong>预测：</strong> $\hat{x}_{k|k-1} = A\hat{x}_{k-1|k-1} + Bu_{k-1}$，$\quad P_{k|k-1} = AP_{k-1|k-1}A^{\top} + W$</p>
<p><strong>修正：</strong></p>
$$K_k = P_{k|k-1}C^{\top}\big(CP_{k|k-1}C^{\top} + V\big)^{-1}$$
$$\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k\underbrace{(y_k - C\hat{x}_{k|k-1})}_{\text{新息}},\qquad P_{k|k} = (I - K_kC)P_{k|k-1}$$</div>`
    },
    {
      icon: "🧾", title: "完整推导（最小方差最优性证明）",
      html: String.raw`
<div class="box proof"><span class="box-title">第一步：预测方程</span>
<p>误差 $e_{k|k-1} = x_k - \hat{x}_{k|k-1} = A e_{k-1|k-1} + w_{k-1}$（$Bu$ 已知，消去）。$e_{k-1|k-1}$ 与 $w_{k-1}$ 独立（$w$ 是未来噪声），协方差直接相加：</p>
$$P_{k|k-1} = \mathbb{E}[ee^{\top}] = A P_{k-1|k-1}A^{\top} + W \ ✓$$</div>
<div class="box proof"><span class="box-title">第二步：修正——在线性估计类中求最优增益</span>
<p>限定线性修正结构 $\hat{x}_{k|k} = \hat{x}_{k|k-1} + K(y_k - C\hat{x}_{k|k-1})$，$K$ 待定。后验误差：</p>
$$e_{k|k} = x_k - \hat{x}_{k|k} = e_{k|k-1} - K(Ce_{k|k-1} + v_k) = (I - KC)e_{k|k-1} - Kv_k$$
<p>$v_k$ 与 $e_{k|k-1}$ 独立，故</p>
$$P(K) = (I - KC)P_{k|k-1}(I - KC)^{\top} + KVK^{\top}$$
<p>（Joseph 形式，对任意 $K$ 成立。）目标：最小化 $\mathrm{tr}\,P(K)$（均方误差）。展开：</p>
$$\mathrm{tr}\,P(K) = \mathrm{tr}\,P_{k|k-1} - 2\,\mathrm{tr}(KCP_{k|k-1}) + \mathrm{tr}\big(K(CP_{k|k-1}C^{\top} + V)K^{\top}\big)$$
<p>这是 $K$ 的凸二次函数（$S \triangleq CP_{k|k-1}C^{\top} + V \succ 0$）。用矩阵微分 $\frac{\partial}{\partial K}\mathrm{tr}(KM) = M^{\top}$、$\frac{\partial}{\partial K}\mathrm{tr}(KSK^{\top}) = 2KS$，令导数为零：</p>
$$-2P_{k|k-1}C^{\top} + 2K S = 0\ \Longrightarrow\ K_k = P_{k|k-1}C^{\top}S^{-1} \ ✓$$
<p>代回 Joseph 式并化简（用 $K_kS = P_{k|k-1}C^{\top}$）：</p>
$$P_{k|k} = P_{k|k-1} - K_kCP_{k|k-1} - P_{k|k-1}C^{\top}K_k^{\top} + K_kSK_k^{\top} = (I - K_kC)P_{k|k-1} \ ✓$$</div>
<div class="box proof"><span class="box-title">第三步：全局最优性（高斯情形）</span>
<p>高斯先验经线性观测后的条件分布仍高斯，其均值恰为上式（可由联合高斯的条件公式 $\mathbb{E}[x|y] = \mu_x + \Sigma_{xy}\Sigma_{yy}^{-1}(y - \mu_y)$ 逐项对照：$\Sigma_{xy} = P_{k|k-1}C^{\top}$，$\Sigma_{yy} = S$）。条件均值是<strong>一切</strong>（含非线性）估计器中均方误差最小者，故高斯噪声下 Kalman 滤波全局最优；非高斯时仍是<strong>线性</strong>估计器中的最优（以上推导未用高斯性）。归纳起点 $k=0$ 由先验 $\hat x_0, P_0$ 给定，逐步递推保持最优。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">性质与家族</span>
<ul>
<li><strong>新息正交性</strong>：最优性等价于新息序列白且与历史正交（投影定理视角）；工程上监测新息是否"白"即可诊断滤波器健康；</li>
<li>$P_{k|k-1}$ 的递推展开即 <strong>DARE</strong>——与 LQR 的 Riccati 对偶；$(A, C)$ 可检测 + $(A, W^{1/2})$ 可镇定 ⇒ $P$ 收敛到稳态（可用定常增益）；</li>
<li>非线性扩展：EKF（线性化）、UKF（sigma 点）、粒子滤波（蒙特卡洛）。</li>
</ul></div>`
    },
  ],
  code: {
    title: "30 行 Kalman 滤波器 + 新息白性自检",
    text: `import numpy as np

rng = np.random.default_rng(1)
dt = 0.1
A = np.array([[1, dt], [0, 1]])          # 匀速模型
C = np.array([[1., 0.]])
W = np.array([[dt**3/3, dt**2/2], [dt**2/2, dt]]) * 0.3**2
V = np.array([[0.5**2]])

# 生成真值与量测
x = np.array([0., 1.]); xs, ys = [], []
for k in range(200):
    xs.append(x.copy())
    ys.append(C @ x + rng.normal(0, 0.5, 1))
    x = A @ x + rng.multivariate_normal([0, 0], W)

# Kalman 滤波
xh = np.zeros(2); P = np.eye(2) * 10
est, innov = [], []
for y in ys:
    xh = A @ xh                          # 预测
    P = A @ P @ A.T + W
    S = C @ P @ C.T + V                  # 修正
    K = P @ C.T @ np.linalg.inv(S)
    nu = y - C @ xh
    xh = xh + (K @ nu).flatten()
    P = (np.eye(2) - K @ C) @ P
    est.append(xh.copy()); innov.append(nu[0] / np.sqrt(S[0,0]))

xs, est = np.array(xs), np.array(est)
print("位置RMSE  测量:", np.std([y[0]-x[0] for y,x in zip(ys,xs)]).round(3),
      " KF:", np.sqrt(np.mean((est[:,0]-xs[:,0])**2)).round(3))
inn = np.array(innov)
print("归一化新息 均值", inn.mean().round(3), "方差", inn.var().round(3),
      "(白性自检: 应约为 0 和 1)")`
  }
};

/* ════════════════ LQG ════════════════ */
window.KG_CONTENT["lqg"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>现实是"双重不完美"：状态测不全 + 到处是噪声。LQG 给出线性-二次-高斯世界的完整答案：<strong>Kalman 滤波估计状态，LQR 拿估计值当真值用</strong>。</p>
<p>这个"当真值用"（确定性等价）居然是严格最优的——这就是随机版分离定理，其成立高度依赖线性+二次+高斯三件套的特殊结构。它也是一个警世故事的主角：分离设计在名义模型下最优，鲁棒性却可能一塌糊涂（Doyle 1978），由此催生了整个鲁棒控制学科。</p>`
    },
    {
      icon: "📐", title: "问题与解",
      html: String.raw`
<div class="box def"><span class="box-title">LQG 问题</span>
$$\dot{x} = Ax + Bu + w,\quad y = Cx + v,\qquad \min_u\ \mathbb{E}\left[\lim_{T\to\infty}\frac{1}{T}\int_0^T (x^{\top}Qx + u^{\top}Ru)\,dt\right]$$
$w, v$ 为强度 $W, V$ 的白噪声。控制器只能用历史输出 $y|_{[0,t]}$。</div>
<div class="box thm"><span class="box-title">定理（LQG 最优控制器 = LQR ∘ Kalman）</span>
最优控制器为
$$\dot{\hat{x}} = A\hat{x} + Bu + L(y - C\hat{x}),\qquad u = -K\hat{x}$$
其中 $K = R^{-1}B^{\top}P$（控制 CARE，用 $Q, R$），$L = \Sigma C^{\top}V^{-1}$（滤波 CARE，用 $W, V$，对偶）。两个 Riccati 方程<strong>互不耦合</strong>。</div>`
    },
    {
      icon: "🧾", title: "分离定理的证明",
      html: String.raw`
<div class="box proof"><span class="box-title">证明（代价分解论证）</span>
<p><strong>第一步：估计与控制解耦的关键恒等式。</strong>Kalman 滤波的误差 $e = x - \hat{x}$ 满足 $\dot{e} = (A - LC)e + (\text{噪声项})$，其统计特性<strong>与 $u$ 无关</strong>（$u$ 在滤波器中被精确抵消，见《观测器》同款计算）。且最优滤波保证 $\hat x$ 与 $e$ 正交（不相关）。</p>
<p><strong>第二步：代价分解。</strong>把 $x = \hat{x} + e$ 代入代价，利用正交性 $\mathbb{E}[\hat{x}^{\top}Qe] = 0$：</p>
$$\mathbb{E}[x^{\top}Qx + u^{\top}Ru] = \underbrace{\mathbb{E}[\hat{x}^{\top}Q\hat{x} + u^{\top}Ru]}_{J_c:\ \text{只依赖}(\hat x, u)} + \underbrace{\mathbb{E}[e^{\top}Qe]}_{J_e:\ \text{与 } u \text{ 无关}}$$
<p>$J_e$ 是估计的"固定税"，任何控制器都得交，最小化它的正是 Kalman 滤波（已证）。</p>
<p><strong>第三步：对 $J_c$ 应用 LQR。</strong>$\hat{x}$ 的动力学为 $\dot{\hat{x}} = A\hat{x} + Bu + L\nu$，其中新息 $\nu = y - C\hat x$ 是白噪声（Kalman 滤波性质）——这是一个<strong>状态完全已知</strong>（$\hat x$ 就是滤波器自己的状态）的线性系统加白噪声激励。随机 LQR 的最优解与确定性 LQR 相同（白噪声只给代价加一个与控制无关的常数项，由完成平方论证可见交叉项期望为零），即 $u = -K\hat{x}$。</p>
<p>两步各自最优且互不影响对方的可达下界，故组合最优。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">警示：Doyle 1978——"LQG 没有保证的裕度"</span>
<p>LQR 单用有 60° 相位裕度、无穷增益裕度；Kalman 滤波对偶地也有。但 Doyle 的著名单页论文用一个二阶反例证明：LQG 合体后的裕度可以<strong>任意小</strong>。直觉：控制回路的鲁棒性经过观测器的"转述"可能被彻底扭曲。这推动了 LTR（回路传输恢复）与 H∞ 的诞生——见《H∞ 鲁棒控制》。</p></div>`
    },
  ],
  code: {
    title: "LQG 全栈仿真（含 Doyle 式鲁棒性压力测试）",
    text: `import numpy as np
from scipy.linalg import solve_continuous_are

A = np.array([[0., 1.], [3., 0.]])       # 不稳定对象
B = np.array([[0.], [1.]])
C = np.array([[1., 0.]])
Q = np.diag([10., 1.]); R = np.array([[1.]])
W = np.diag([0.0, 1.0]); V = np.array([[0.01]])

P = solve_continuous_are(A, B, Q, R)
K = np.linalg.solve(R, B.T @ P)
S = solve_continuous_are(A.T, C.T, W, V)
L = S @ C.T @ np.linalg.inv(V)

rng = np.random.default_rng(0)
h, T = 1e-3, 10.0
def run(gain_error=1.0):
    """gain_error: 执行器实际增益/名义增益, 测试鲁棒性"""
    x = np.array([1.0, 0.0]); xh = np.zeros(2); J = 0.0
    for k in range(int(T/h)):
        y = C @ x + rng.normal(0, np.sqrt(V[0,0]/h), 1)
        u = float(-K @ xh)
        J += float(x @ Q @ x + u * R[0,0] * u) * h
        x  = x  + h*(A@x + B.flatten()*u*gain_error) \
                + np.sqrt(h)*rng.multivariate_normal([0,0], W)
        xh = xh + h*(A@xh + B.flatten()*u + (L@(y - C@xh)).flatten())
        if not np.isfinite(x).all() or abs(x[0]) > 1e6:
            return np.inf
    return J

for ge in [1.0, 1.15, 1.3, 0.85]:
    J = run(ge)
    print(f"执行器增益偏差 x{ge}: J = {J if J==np.inf else round(J,1)}")
# 名义最优, 但少许模型失配就可能显著劣化甚至失稳`
  }
};
