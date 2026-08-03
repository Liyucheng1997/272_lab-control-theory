/* 内容：结构性质与状态反馈 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ 能控性 ════════════════ */
window.KG_CONTENT["controllability"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>方向盘控制不了油量表的读数，油门控制不了收音机音量——输入"够得着"哪些状态，是系统的<strong>结构性质</strong>，与控制器设计得多聪明无关。</p>
<p>能控性回答的问题是：<strong>用输入 $u$ 能否在有限时间内把状态从任意 $x_0$ 开到任意 $x_1$？</strong>如果不能，那些够不着的方向就是"死角"：死角里若藏着不稳定模态，任何控制器都救不了（不可镇定）。这是 Kalman 1960 年提出的概念，是现代控制理论的第一块基石。</p>
<p>为什么判据里出现 $B, AB, A^2B, \dots$？直觉：$B$ 是输入即刻推动的方向；推动后状态被动力学 $A$ "搬运"，一阶时间后到达 $AB$ 方向；再搬运到 $A^2B$……输入通过动力学的"传播链"能到达的方向，就是这些向量张成的空间。</p>`
    },
    {
      icon: "📐", title: "定义与判据",
      html: String.raw`
<div class="box def"><span class="box-title">定义（能控性）</span>
系统 $\dot{x} = Ax + Bu$ 称为<strong>能控</strong>，若对任意 $x_0, x_1 \in \mathbb{R}^n$ 与任意 $T > 0$，存在输入 $u(\cdot)$ 使 $x(0) = x_0,\ x(T) = x_1$。</div>
<div class="box thm"><span class="box-title">定理（三大等价判据）</span>
以下等价：
<ol>
<li>$(A, B)$ 能控；</li>
<li><strong>Kalman 秩判据</strong>：$\mathrm{rank}\ \mathcal{C} = n$，其中 $\mathcal{C} = [B\ \ AB\ \ A^2B\ \cdots\ A^{n-1}B]$；</li>
<li><strong>能控 Gram 矩阵</strong> $W_c(T) = \displaystyle\int_0^T e^{A\tau}BB^{\top}e^{A^{\top}\tau}d\tau \succ 0$ 对某个（等价于所有）$T > 0$；</li>
<li><strong>PBH 判据</strong>：$\mathrm{rank}\,[\,\lambda I - A \ \ B\,] = n$ 对所有 $\lambda \in \mathbb{C}$（只需检查特征值）。</li>
</ol></div>`
    },
    {
      icon: "🧾", title: "证明",
      html: String.raw`
<div class="box proof"><span class="box-title">证明：(3) ⇒ (1)（Gram 矩阵可逆则能控，构造性）</span>
<p>设 $W_c(T) \succ 0$。给定 $x_0, x_1$，直接给出显式控制律：</p>
$$u(t) = B^{\top}e^{A^{\top}(T-t)}\,W_c(T)^{-1}\big(x_1 - e^{AT}x_0\big)$$
<p>代入解公式验证：</p>
$$x(T) = e^{AT}x_0 + \int_0^T e^{A(T-\tau)}B\,u(\tau)\,d\tau = e^{AT}x_0 + \left[\int_0^T e^{A(T-\tau)}BB^{\top}e^{A^{\top}(T-\tau)}d\tau\right]W_c(T)^{-1}(x_1 - e^{AT}x_0)$$
<p>方括号内换元 $s = T-\tau$ 恰为 $W_c(T)$，与 $W_c(T)^{-1}$ 相消，得 $x(T) = x_1$。✓</p></div>
<div class="box proof"><span class="box-title">证明：(1) ⇒ (2) 的逆否，即 rank 亏 ⇒ 不能控</span>
<p>设 $\mathrm{rank}\,\mathcal{C} < n$，则存在 $v \ne 0$ 使 $v^{\top}A^kB = 0,\ k = 0, \dots, n-1$。由 <strong>Cayley–Hamilton 定理</strong>（$A^n$ 是 $I, A, \dots, A^{n-1}$ 的线性组合），归纳得 $v^{\top}A^kB = 0$ 对<strong>一切</strong> $k \ge 0$ 成立。于是</p>
$$v^{\top}e^{At}B = \sum_{k=0}^{\infty}\frac{t^k}{k!}v^{\top}A^kB = 0,\quad \forall t$$
<p>取 $x_0 = 0$，则对任何输入 $v^{\top}x(T) = \int_0^T v^{\top}e^{A(T-\tau)}Bu(\tau)d\tau = 0$：状态永远到不了 $v$ 方向有分量的目标（如 $x_1 = v$），不能控。✓</p></div>
<div class="box proof"><span class="box-title">证明：(2) ⇒ (3) 的逆否，即 Gram 奇异 ⇒ rank 亏</span>
<p>设 $W_c(T)$ 奇异：存在 $v \ne 0$ 使 $v^{\top}W_c(T)v = \int_0^T \|B^{\top}e^{A^{\top}\tau}v\|^2 d\tau = 0$，被积函数连续非负故恒为零：$v^{\top}e^{A\tau}B \equiv 0$。在 $\tau = 0$ 处逐阶求导：$v^{\top}B = 0,\ v^{\top}AB = 0,\ \dots,\ v^{\top}A^{n-1}B = 0$，即 $v^{\top}\mathcal{C} = 0$，rank 亏。✓</p>
<p>综上 (1)⇒(2)⇒(3)⇒(1) 闭环成立。</p></div>
<div class="box proof"><span class="box-title">证明：(2) ⟺ (4)（PBH）</span>
<p><strong>(4) 不成立 ⇒ (2) 不成立：</strong>若某 $\lambda$ 处 rank 亏，则存在行向量 $v^{\top} \ne 0$：$v^{\top}(\lambda I - A) = 0$ 且 $v^{\top}B = 0$，即 $v$ 是 $A^{\top}$ 的特征向量（$v^{\top}A = \lambda v^{\top}$）且 $\perp B$。于是 $v^{\top}A^kB = \lambda^k v^{\top}B = 0$ 对一切 $k$，故 $v^{\top}\mathcal{C} = 0$。</p>
<p><strong>(2) 不成立 ⇒ (4) 不成立：</strong>设 $\mathrm{rank}\,\mathcal{C} = r < n$。取不可控子空间分解（取 $\mathcal{C}$ 列空间的基扩充为全空间基做相似变换）：</p>
$$\tilde{A} = \begin{bmatrix} A_{11} & A_{12} \\ 0 & A_{22} \end{bmatrix},\qquad \tilde{B} = \begin{bmatrix} B_1 \\ 0 \end{bmatrix}$$
<p>（列空间 $\mathrm{Im}\,\mathcal{C}$ 是 $A$-不变的且包含 $\mathrm{Im}\,B$，故变换后有此形状。）取 $A_{22}$ 的任一特征值 $\lambda$ 与左特征向量 $w^{\top}$，则 $[0\ \ w^{\top}]$ 左零化 $[\lambda I - \tilde{A}\ \ \tilde{B}]$，rank 亏。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">能镇定（Stabilizability）</span>
实践中常不需要完全能控，只需<strong>不能控的模态自己是稳定的</strong>：PBH 只对 $\mathrm{Re}\lambda \ge 0$ 的特征值检查满秩即可。这是 LQR 存在稳定解的标准假设。</div>`
    },
  ],
  code: {
    title: "能控性判定 + 用 Gram 矩阵构造最小能量控制",
    text: `import numpy as np
from scipy.linalg import expm

A = np.array([[0., 1.], [0., 0.]])       # 双积分器
B = np.array([[0.], [1.]])
n = A.shape[0]

# Kalman 秩判据
C = np.hstack([np.linalg.matrix_power(A, k) @ B for k in range(n)])
print("rank C =", np.linalg.matrix_rank(C), "/", n)

# Gram 矩阵数值积分
T, N = 2.0, 2000
h = T / N
Wc = np.zeros((n, n))
for k in range(N):
    t = k * h
    E = expm(A * t) @ B
    Wc += E @ E.T * h

# 最小能量控制: 从 [0,0] 到 [1,0]
x0, x1 = np.zeros(2), np.array([1., 0.])
eta = np.linalg.solve(Wc, x1 - expm(A*T) @ x0)

# 前向仿真验证
x = x0.copy()
for k in range(N):
    t = k * h
    u = (B.T @ expm(A.T * (T - t)) @ eta)[0]
    x = x + h * (A @ x + B.flatten() * u)
print("到达状态:", np.round(x, 4), " 目标:", x1)`
  }
};

/* ════════════════ 能观性 ════════════════ */
window.KG_CONTENT["observability"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>传感器有限而状态很多：只测得小车位置，速度呢？只测得电机电流，转子温度呢？能观性回答：<strong>盯着输出 $y(t)$ 看一段时间，能否唯一反推出内部初始状态 $x_0$？</strong></p>
<p>能推出来，就可以造《状态观测器》，把"看不见的状态"估计出来喂给状态反馈；推不出来，说明有些内部运动在输出上完全隐形——若隐形的恰是不稳定模态，系统会在你毫无察觉时坏掉（不可检测）。</p>
<p>最美的事实是：<strong>能观性 = 转置系统的能控性</strong>。一条对偶原理，让能控性的全部定理免费复用。</p>`
    },
    {
      icon: "📐", title: "定义与判据",
      html: String.raw`
<div class="box def"><span class="box-title">定义（能观性）</span>
系统 $\dot{x} = Ax,\ y = Cx$（输入已知故不失一般性设 $u=0$）称为<strong>能观</strong>，若任意 $T>0$ 内的输出轨迹 $y(\cdot)|_{[0,T]}$ 唯一决定初始状态 $x_0$。</div>
<div class="box thm"><span class="box-title">定理（判据）</span>
以下等价：
<ol>
<li>$(A, C)$ 能观；</li>
<li>$\mathrm{rank}\ \mathcal{O} = n$，其中 $\mathcal{O} = \begin{bmatrix} C \\ CA \\ \vdots \\ CA^{n-1} \end{bmatrix}$；</li>
<li>能观 Gram 矩阵 $W_o(T) = \int_0^T e^{A^{\top}\tau}C^{\top}Ce^{A\tau}\,d\tau \succ 0$；</li>
<li>PBH：$\mathrm{rank}\begin{bmatrix} \lambda I - A \\ C \end{bmatrix} = n,\ \forall \lambda$。</li>
</ol></div>`
    },
    {
      icon: "🧾", title: "证明：对偶原理一步到位",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（对偶原理）</span>
$(A, C)$ 能观 $\iff$ $(A^{\top}, C^{\top})$ 能控。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>先独立证明"能观 ⟺ $\mathcal{O}$ 满秩"，再对比矩阵结构。</p>
<p><strong>($\mathcal{O}$ 满秩 ⇒ 能观)</strong>：$y(t) = Ce^{At}x_0$。在 $t = 0$ 处求各阶导数：</p>
$$y(0) = Cx_0,\quad \dot{y}(0) = CAx_0,\ \dots,\ y^{(n-1)}(0) = CA^{n-1}x_0 \ \Longrightarrow\ \begin{bmatrix} y(0) \\ \vdots \\ y^{(n-1)}(0) \end{bmatrix} = \mathcal{O}x_0$$
<p>$\mathcal{O}$ 列满秩则该线性方程组解唯一，$x_0$ 被 $y$ 唯一决定。（若不想用导数：$W_o(T)x_0 = \int_0^T e^{A^\top \tau}C^{\top}y(\tau)d\tau$ 同样给出显式重构公式。）</p>
<p><strong>(rank 亏 ⇒ 不能观)</strong>：存在 $v \ne 0$，$CA^kv = 0\ (k \le n-1)$；由 Cayley–Hamilton 推广到一切 $k$，于是 $Ce^{At}v \equiv 0$：初值 $x_0 = v$ 与 $x_0 = 0$ 产生<strong>完全相同</strong>的输出（恒为零），无法区分。✓</p>
<p><strong>对偶</strong>：注意 $\mathcal{O}(A, C)^{\top} = [C^{\top}\ \ A^{\top}C^{\top}\ \cdots\ (A^{\top})^{n-1}C^{\top}] = \mathcal{C}(A^{\top}, C^{\top})$，转置不改变秩，故两个满秩条件相同。Gram 判据与 PBH 判据同样在转置下互换。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">对偶的深远回报</span>
<ul>
<li>极点配置定理（能控 ⇒ 反馈极点任意配置）对偶出：<strong>能观 ⇒ 观测器极点任意配置</strong>（对 $(A^{\top}, C^{\top})$ 配置后转置回来）；</li>
<li>LQR 与 Kalman 滤波是同一 Riccati 方程的对偶双生子（$A \leftrightarrow A^{\top}$, $B \leftrightarrow C^{\top}$, $Q \leftrightarrow$ 过程噪声, $R \leftrightarrow$ 测量噪声）；</li>
<li>可检测性（Detectability）：不能观模态自身稳定即可——是 Kalman 滤波收敛的标准假设。</li>
</ul></div>`
    },
  ],
  code: {
    title: "能观性判定 + 从输出重构初始状态",
    text: `import numpy as np
from scipy.linalg import expm

A = np.array([[0., 1.], [-2., -0.6]])
C = np.array([[1., 0.]])            # 只测位置
n = 2

O = np.vstack([C @ np.linalg.matrix_power(A, k) for k in range(n)])
print("rank O =", np.linalg.matrix_rank(O))

# 真实初值(待反推)
x0_true = np.array([0.7, -1.2])

# 记录输出 y(t) = C e^{At} x0
T, N = 3.0, 600
h = T / N
ys, ts = [], []
for k in range(N):
    t = k * h
    ys.append((C @ expm(A*t) @ x0_true)[0])
    ts.append(t)

# 用 Gram 矩阵公式重构: Wo x0 = int e^{A't} C' y dt
Wo = np.zeros((n, n)); rhs = np.zeros(n)
for k in range(N):
    E = expm(A * ts[k])
    Wo  += E.T @ C.T @ C @ E * h
    rhs += E.T @ C.T @ np.atleast_1d(ys[k]) * h
x0_hat = np.linalg.solve(Wo, rhs)
print("真实 x0:", x0_true, " 重构 x0:", np.round(x0_hat, 4))`
  }
};

/* ════════════════ 极点配置 ════════════════ */
window.KG_CONTENT["pole-placement"] = {
  sim: "sim-poleplace",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>系统的性格由极点（$A$ 的特征值）决定：快慢、振荡、稳定与否。状态反馈 $u = -Kx$ 把系统矩阵改写为 $A - BK$——相当于<strong>给系统做性格移植</strong>。</p>
<p>极点配置定理是状态空间方法的"成人礼"：只要系统能控，闭环极点可以配到<strong>任何</strong>你想要的位置（共轭对称即可）。想要多快有多快——代价是控制量随极点左移而暴涨（实验里能直观看到）。这也解释了为什么需要 LQR：它用最优化替你在"快"与"省"之间做系统性的权衡。</p>`
    },
    {
      icon: "🧾", title: "定理与完整证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（极点配置 / Pole Placement）</span>
设 $(A, B)$ 能控（单输入情形 $B = b$）。则对任意给定的共轭对称多项式
$$p^*(s) = s^n + \alpha_1 s^{n-1} + \cdots + \alpha_n$$
存在（且单输入时唯一）反馈增益 $K$，使 $\det(sI - A + bK) = p^*(s)$。反之若 $(A,B)$ 不能控，则不能控模态的特征值无法被反馈移动。</div>
<div class="box proof"><span class="box-title">证明（构造性，单输入）</span>
<p><strong>第一步：化为能控标准型。</strong>设 $A$ 的特征多项式为 $\det(sI - A) = s^n + a_1 s^{n-1} + \cdots + a_n$。定义变换阵</p>
$$T = \mathcal{C}\,\mathcal{M},\qquad \mathcal{C} = [b\ \ Ab\ \cdots\ A^{n-1}b],\qquad
\mathcal{M} = \begin{bmatrix} a_{n-1} & a_{n-2} & \cdots & 1 \\ a_{n-2} & \cdots & 1 & \\ \vdots & ⋰ & & \\ 1 & & & \end{bmatrix}$$
<p>能控性保证 $\mathcal{C}$ 可逆，$\mathcal{M}$ 是主反对角线为 1 的三角阵故可逆，从而 $T$ 可逆。可直接验证（逐列核对 $AT$ 与 $Tb$ 的关系，利用 Cayley–Hamilton 处理最后一列）：</p>
$$\tilde{A} = T^{-1}AT = \begin{bmatrix} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1 \\ -a_n & \cdots & -a_2 & -a_1 \end{bmatrix},\qquad \tilde{b} = T^{-1}b = \begin{bmatrix} 0 \\ \vdots \\ 0 \\ 1 \end{bmatrix}$$
<p>即能控标准型：伴随矩阵形式，特征系数全部暴露在最后一行。</p>
<p><strong>第二步：在标准型下配置。</strong>取 $\tilde{K} = [\alpha_n - a_n,\ \ \alpha_{n-1} - a_{n-1},\ \dots,\ \alpha_1 - a_1]$，则</p>
$$\tilde{A} - \tilde{b}\tilde{K} = \begin{bmatrix} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1 \\ -\alpha_n & \cdots & -\alpha_2 & -\alpha_1 \end{bmatrix}$$
<p>——只改动了最后一行（$\tilde b$ 只有末位非零）。伴随矩阵的特征多项式就是其末行系数：$\det(sI - \tilde{A} + \tilde{b}\tilde{K}) = s^n + \alpha_1 s^{n-1} + \cdots + \alpha_n = p^*(s)$。✓</p>
<p><strong>第三步：变换回去。</strong>取 $K = \tilde{K}T^{-1}$，则 $A - bK = T(\tilde{A} - \tilde{b}\tilde{K})T^{-1}$，相似变换不改变特征多项式。</p>
<p><strong>唯一性</strong>（单输入）：$p^*$ 的 $n$ 个系数与 $K$ 的 $n$ 个分量之间的映射经由标准型是一一仿射对应。</p>
<p><strong>必要性</strong>：不能控时按能控性分解 $\tilde A = \begin{bmatrix} A_{11} & A_{12} \\ 0 & A_{22} \end{bmatrix}, \tilde b = \begin{bmatrix} b_1 \\ 0 \end{bmatrix}$，任何反馈 $K = [K_1\ K_2]$ 给出</p>
$$\tilde A - \tilde b K = \begin{bmatrix} A_{11} - b_1K_1 & A_{12} - b_1K_2 \\ 0 & A_{22} \end{bmatrix}$$
<p>块三角结构使 $A_{22}$（不能控模态）的特征值原封不动。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">实用备注</span>
<ul>
<li><strong>Ackermann 公式</strong>把三步合成一行：$K = [0, \dots, 0, 1]\,\mathcal{C}^{-1}\,p^*(A)$；</li>
<li>多输入时 $K$ 不唯一（自由度可用于特征向量配置/鲁棒性，MATLAB 的 <em>place</em>）；</li>
<li>极点不是越快越好：$u = -Kx$ 峰值随极点模长急剧增大（大约按 $\omega^n$ 量级），执行器饱和、噪声放大、未建模高频动态都会反噬——实验中把 $\omega_n$ 拉到 8 感受一下控制量。</li>
</ul></div>`
    },
  ],
  code: {
    title: "Ackermann 公式实现极点配置",
    text: `import numpy as np

def acker(A, b, desired_poles):
    n = A.shape[0]
    C = np.hstack([np.linalg.matrix_power(A, k) @ b for k in range(n)])
    assert np.linalg.matrix_rank(C) == n, "not controllable"
    pstar = np.poly(desired_poles)          # [1, a1, ..., an]
    # p*(A) = A^n + a1 A^{n-1} + ... + an I
    PA = np.zeros_like(A)
    for i, c in enumerate(pstar):
        PA += c * np.linalg.matrix_power(A, n - i)
    e_n = np.zeros(n); e_n[-1] = 1.0
    return e_n @ np.linalg.inv(C) @ PA

# 倒立摆线性化模型 (theta 从竖直向上)
g, l = 9.81, 1.0
A = np.array([[0., 1.], [g/l, 0.]])         # 不稳定: 特征值 +-sqrt(g/l)
b = np.array([[0.], [1.]])

K = acker(A, b, [-3+2j, -3-2j])
Acl = A - b @ K.reshape(1, -1)
print("K =", np.round(K, 3))
print("闭环极点:", np.round(np.linalg.eigvals(Acl), 3))   # 应为 -3±2j`
  }
};

/* ════════════════ 观测器与分离原理 ════════════════ */
window.KG_CONTENT["observer"] = {
  sim: "sim-observer",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>状态反馈 $u = -Kx$ 很美好，但传感器测不到全部状态。Luenberger（1964）的方案：<strong>在计算机里跑一个系统的仿真副本</strong>，并用"真实输出与仿真输出之差"（新息）去不断修正副本：</p>
<div class="box eg"><span class="box-title">观测器 = 仿真 + 纠错</span>
$$\dot{\hat{x}} = \underbrace{A\hat{x} + Bu}_{\text{照抄模型仿真}} + \underbrace{L(y - C\hat{x})}_{\text{用测量误差纠偏}}$$</div>
<p>纠错增益 $L$ 大，收敛快但放大测量噪声；$L$ 小则相反（实验里可直接体验这个权衡）。更深的问题是：拿估计值 $\hat{x}$ 代替真值做反馈 $u = -K\hat{x}$，两个动态过程互相纠缠，整体还稳吗？<strong>分离原理</strong>给出漂亮回答：控制器与观测器可以各设计各的，合起来自动稳。</p>`
    },
    {
      icon: "🧾", title: "观测器误差动态与极点配置",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（观测器收敛）</span>
定义估计误差 $e = x - \hat{x}$，则
$$\dot{e} = (A - LC)\,e$$
且当 $(A, C)$ 能观时，$L$ 可使 $A - LC$ 的特征值任意配置（从而 $e \to 0$ 以任意指定速率）。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>两式相减：</p>
$$\dot{e} = \dot{x} - \dot{\hat{x}} = (Ax + Bu) - (A\hat{x} + Bu + L(Cx - C\hat{x})) = Ae - LCe = (A - LC)e$$
<p>注意 $Bu$ 完全消去——观测器"免疫"于任何已知输入。误差极点配置：$(A, C)$ 能观 $\iff$ $(A^{\top}, C^{\top})$ 能控（对偶原理），由极点配置定理存在 $M$ 使 $A^{\top} - C^{\top}M$ 特征值任意；取 $L = M^{\top}$，而 $A - LC = (A^{\top} - C^{\top}M)^{\top}$，转置不改变特征值。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "分离原理：定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（分离原理 / Separation Principle）</span>
将状态反馈 $u = -K\hat{x}$ 与观测器联立，得 $2n$ 维闭环系统。其特征值恰为
$$\mathrm{eig}(A - BK)\ \cup\ \mathrm{eig}(A - LC)$$
即：<strong>控制器极点与观测器极点互不干扰，可独立设计</strong>。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>取坐标 $(x, e)$（$e = x - \hat{x}$）。控制律 $u = -K\hat{x} = -K(x - e)$，代入被控对象：</p>
$$\dot{x} = Ax - BK(x - e) = (A - BK)x + BKe$$
<p>与误差方程合并：</p>
$$\begin{bmatrix} \dot{x} \\ \dot{e} \end{bmatrix} = \begin{bmatrix} A - BK & BK \\ 0 & A - LC \end{bmatrix}\begin{bmatrix} x \\ e \end{bmatrix}$$
<p>块上三角矩阵的特征多项式等于对角块特征多项式之积：</p>
$$\det\begin{bmatrix} sI - (A-BK) & -BK \\ 0 & sI - (A-LC) \end{bmatrix} = \det(sI - A + BK)\cdot\det(sI - A + LC)$$
<p>（对块三角阵按 Schur 展开即得。）两组极点均在左半平面 ⇒ 整体渐近稳定。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">工程释义与局限</span>
<ul>
<li>经验法则：观测器极点取控制器极点的 2~5 倍快，让"估计先于控制"到位；太快则噪声穿透。</li>
<li>合并观测器与反馈，从 $y$ 到 $u$ 看是一个动态输出反馈控制器——现代输出反馈（LQG、H∞）都是这个结构。</li>
<li><strong>警告</strong>：分离原理只保证名义模型下的极点；LQR 单独的优良鲁棒裕度在接上观测器后可能全部丢失（Doyle 1978 著名反例，"LQG 没有保证的裕度"）。这是《LQG》与《H∞》一节的出发点。</li>
<li>非线性版本（EKF + 反馈）一般<strong>不</strong>成立分离原理，只有局部/特殊结构结果。</li>
</ul></div>`
    },
  ],
  code: {
    title: "观测器 + 状态反馈闭环：验证分离原理",
    text: `import numpy as np

# 对象: 双积分器, 只测位置
A = np.array([[0., 1.], [0., 0.]])
B = np.array([[0.], [1.]])
C = np.array([[1., 0.]])

def acker(A, b, poles):
    n = A.shape[0]
    Cm = np.hstack([np.linalg.matrix_power(A, k) @ b for k in range(n)])
    pstar = np.poly(poles)
    PA = sum(c * np.linalg.matrix_power(A, n - i) for i, c in enumerate(pstar))
    e = np.zeros(n); e[-1] = 1
    return (e @ np.linalg.inv(Cm) @ PA).reshape(1, -1)

K = acker(A, B, [-2+1j, -2-1j])                 # 控制器极点
L = acker(A.T, C.T, [-8., -9.]).T               # 观测器极点(对偶)

# 组合系统特征值 = 两组之并(分离原理)
top = np.hstack([A - B@K, B@K])
bot = np.hstack([np.zeros((2,2)), A - L@C])
Acl = np.vstack([top, bot])
print("闭环特征值:", np.sort_complex(np.linalg.eigvals(Acl)))
print("理论: -2±1j 与 -8, -9 的并集")

# 时域仿真
x = np.array([1., 0.]); xh = np.zeros(2)
h = 0.001
for k in range(8000):
    y = C @ x
    u = float(-K @ xh)
    x  = x  + h * (A @ x + B.flatten()*u)
    xh = xh + h * (A @ xh + B.flatten()*u + (L @ (y - C @ xh)).flatten())
print("末态 x:", np.round(x, 5), " 估计误差:", np.round(x - xh, 6))`
  }
};
