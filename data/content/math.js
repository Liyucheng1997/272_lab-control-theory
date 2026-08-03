/* 内容：数学与系统基础 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ 状态空间模型 ════════════════ */
window.KG_CONTENT["state-space"] = {
  sim: "sim-eig-phase",
  sections: [
    {
      icon: "💡", title: "直觉：为什么是状态空间",
      html: String.raw`
<p>控制理论研究的对象是<strong>随时间演化的系统</strong>：电机、无人机、化学反应釜、经济系统……它们千差万别，但都能用同一句话概括：</p>
<div class="box eg"><span class="box-title">核心思想</span>
系统在未来如何变化，只取决于它<strong>现在的状态</strong>和<strong>现在的输入</strong>，与历史无关。把"现在需要知道的最少信息"打包成一个向量 $x$，就是<strong>状态</strong>。</div>
<p>例如一辆小车：知道当前位置和速度（$x = [p, v]^{\top}$），加上未来的油门（输入 $u$），就能完全预测它的未来——不需要知道它上周去过哪里。这个"位置+速度"就是状态。</p>
<p>状态空间方法把所有系统统一写成<strong>一阶向量微分方程</strong>，于是整套线性代数（特征值、秩、二次型）都成为分析工具。这是"现代控制理论"（1960s 之后）的出发点，区别于只看输入输出关系的"经典控制"（见 <em>传递函数</em>）。</p>`
    },
    {
      icon: "📐", title: "数学定义",
      html: String.raw`
<div class="box def"><span class="box-title">定义（非线性状态空间模型）</span>
$$\dot{x}(t) = f(x(t), u(t)), \qquad y(t) = h(x(t), u(t))$$
其中 $x \in \mathbb{R}^n$ 为状态，$u \in \mathbb{R}^m$ 为输入（控制量），$y \in \mathbb{R}^p$ 为输出（可测量）。$n$ 称为系统的<strong>阶</strong>。</div>
<div class="box def"><span class="box-title">定义（线性时不变系统 LTI）</span>
$$\dot{x} = Ax + Bu, \qquad y = Cx + Du$$
$A \in \mathbb{R}^{n\times n}$ 为系统矩阵（内部动力学），$B \in \mathbb{R}^{n\times m}$ 为输入矩阵，$C \in \mathbb{R}^{p\times n}$ 为输出矩阵，$D$ 为直通项（通常为 0）。</div>
<p><strong>解的存在唯一性</strong>：若 $f$ 关于 $x$ 局部 Lipschitz（$\|f(x_1,u)-f(x_2,u)\| \le L\|x_1-x_2\|$），则初值问题在局部有唯一解（Picard–Lindelöf 定理）。线性系统的 $f=Ax+Bu$ 全局 Lipschitz，故解全局存在唯一——这就是为什么线性理论如此"干净"。</p>`
    },
    {
      icon: "🧾", title: "高阶方程 → 状态空间（构造性证明）",
      html: String.raw`
<div class="box thm"><span class="box-title">命题</span>
任何 $n$ 阶单输入常微分方程
$$y^{(n)} + a_1 y^{(n-1)} + \cdots + a_n y = b\,u$$
都可等价转化为 $n$ 维一阶状态空间模型。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>取状态为输出及其各阶导数（相变量）：$x_1 = y,\ x_2 = \dot{y},\ \dots,\ x_n = y^{(n-1)}$。则</p>
$$\dot{x}_1 = x_2,\quad \dot{x}_2 = x_3,\ \dots,\ \dot{x}_{n-1} = x_n$$
<p>而最后一个方程直接由原方程给出：</p>
$$\dot{x}_n = y^{(n)} = -a_n x_1 - a_{n-1}x_2 - \cdots - a_1 x_n + b u$$
<p>写成矩阵形式即得能控标准型：</p>
$$A = \begin{bmatrix} 0 & 1 & & \\ & \ddots & \ddots & \\ & & 0 & 1 \\ -a_n & -a_{n-1} & \cdots & -a_1 \end{bmatrix},\quad B = \begin{bmatrix} 0 \\ \vdots \\ 0 \\ b \end{bmatrix},\quad C = [1, 0, \dots, 0]$$
<p>反之，沿此定义逐层求导可还原原方程，故两者等价。<span class="qed">∎</span></p></div>
<div class="box eg"><span class="box-title">例：质量-弹簧-阻尼系统</span>
<p>牛顿第二定律：$m\ddot{y} + c\dot{y} + ky = F$。取 $x = [y, \dot{y}]^{\top}$，$u = F$：</p>
$$\dot{x} = \begin{bmatrix} 0 & 1 \\ -k/m & -c/m \end{bmatrix} x + \begin{bmatrix} 0 \\ 1/m \end{bmatrix} u, \qquad y = [1,\ 0]\,x$$
<p>右侧交互实验就是在玩这个 $2\times 2$ 的 $A$ 矩阵：改变元素，轨迹形态（螺旋收敛、发散、周期）完全由 $A$ 的特征值决定——这正是下一站《矩阵指数》和《线性系统稳定性》要证明的。</p></div>`
    },
  ],
  code: {
    title: "从零实现：用 RK4 数值积分仿真质量-弹簧-阻尼系统",
    text: `import numpy as np
import matplotlib.pyplot as plt

# 质量-弹簧-阻尼: m*y'' + c*y' + k*y = F
m, c, k = 1.0, 0.5, 4.0
A = np.array([[0, 1], [-k/m, -c/m]])
B = np.array([0, 1/m])

def f(x, u):
    return A @ x + B * u

def rk4_step(f, x, u, h):
    k1 = f(x, u)
    k2 = f(x + h/2*k1, u)
    k3 = f(x + h/2*k2, u)
    k4 = f(x + h*k3, u)
    return x + h/6*(k1 + 2*k2 + 2*k3 + k4)

h, T = 0.01, 20.0
ts = np.arange(0, T, h)
x = np.array([1.0, 0.0])        # 初始位移 1m
traj = []
for t in ts:
    u = 0.0                      # 自由响应
    traj.append(x.copy())
    x = rk4_step(f, x, u, h)
traj = np.array(traj)

plt.plot(ts, traj[:, 0], label='position y')
plt.plot(ts, traj[:, 1], label='velocity dy/dt')
plt.xlabel('t (s)'); plt.legend(); plt.grid(True)
plt.title('Mass-Spring-Damper free response')
plt.show()`
  }
};

/* ════════════════ 矩阵指数 ════════════════ */
window.KG_CONTENT["matrix-exponential"] = {
  sim: "sim-eig-phase",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>标量方程 $\dot{x} = ax$ 的解是 $x(t) = e^{at}x(0)$，人人都会。矩阵指数就是把这个结论<strong>原封不动搬到向量情形</strong>：$\dot{x} = Ax$ 的解是 $x(t) = e^{At}x(0)$。</p>
<p>所有线性系统分析——稳定性、能控性、频域响应——最终都建立在这个解析解上。它告诉我们：<strong>线性系统的一切行为都编码在 $A$ 的特征结构里</strong>。特征值是增长/衰减率与振荡频率，特征向量是运动的"模态方向"。</p>`
    },
    {
      icon: "📐", title: "定义与良定性",
      html: String.raw`
<div class="box def"><span class="box-title">定义（矩阵指数）</span>
$$e^{At} \triangleq \sum_{k=0}^{\infty} \frac{(At)^k}{k!} = I + At + \frac{A^2t^2}{2!} + \cdots$$</div>
<div class="box proof"><span class="box-title">收敛性证明</span>
<p>取任意诱导范数，利用 $\|A^k\| \le \|A\|^k$：</p>
$$\sum_{k=0}^{\infty} \left\| \frac{A^k t^k}{k!} \right\| \le \sum_{k=0}^{\infty} \frac{(\|A\| |t|)^k}{k!} = e^{\|A\||t|} < \infty$$
<p>级数绝对收敛，且在任意有界区间上一致收敛（Weierstrass 判别法），故 $e^{At}$ 对一切 $t$ 良定义，且可逐项求导。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "核心定理：线性系统的完整解",
      html: String.raw`
<div class="box thm"><span class="box-title">定理 1</span>
$\dfrac{d}{dt}e^{At} = Ae^{At} = e^{At}A$，且 $x(t) = e^{At}x_0$ 是 $\dot{x}=Ax,\ x(0)=x_0$ 的唯一解。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>由一致收敛性，可逐项求导：</p>
$$\frac{d}{dt}e^{At} = \sum_{k=1}^{\infty} \frac{A^k t^{k-1}}{(k-1)!} = A\sum_{j=0}^{\infty} \frac{A^j t^j}{j!} = Ae^{At}$$
<p>（把 $A$ 提在右边同理得 $e^{At}A$。）于是 $x(t)=e^{At}x_0$ 满足方程与初值。唯一性：设 $z(t)$ 是另一解，令 $g(t) = e^{-At}z(t)$，则</p>
$$\dot{g} = -Ae^{-At}z + e^{-At}\dot{z} = -Ae^{-At}z + e^{-At}Az = 0$$
<p>（$e^{-At}$ 与 $A$ 可交换）。故 $g(t) \equiv g(0) = x_0$，即 $z(t) = e^{At}x_0$。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理 2（常数变易公式，全响应）</span>
$\dot{x} = Ax + Bu$ 的解为
$$x(t) = \underbrace{e^{At}x_0}_{\text{零输入响应}} + \underbrace{\int_0^t e^{A(t-\tau)}Bu(\tau)\,d\tau}_{\text{零状态响应（卷积）}}$$</div>
<div class="box proof"><span class="box-title">证明</span>
<p>猜测形式 $x(t) = e^{At}g(t)$（常数变易）。代入方程：</p>
$$\dot{x} = Ae^{At}g + e^{At}\dot{g} \overset{!}{=} Ae^{At}g + Bu \ \Rightarrow\ \dot{g} = e^{-At}Bu(t)$$
<p>积分得 $g(t) = x_0 + \int_0^t e^{-A\tau}Bu(\tau)d\tau$，乘回 $e^{At}$ 并用半群性质 $e^{At}e^{-A\tau} = e^{A(t-\tau)}$（同一矩阵 $A$ 可交换，二项式展开逐项验证）即得。直接求导可验证它确实满足方程，由解的唯一性得证。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">怎么计算 $e^{At}$</span>
<p>① 若 $A = V\Lambda V^{-1}$ 可对角化，则 $A^k = V\Lambda^k V^{-1}$ 代入级数：$e^{At} = V\,\mathrm{diag}(e^{\lambda_i t})\,V^{-1}$——每个特征值贡献一个模态 $e^{\lambda_i t}$；② 不可对角化时用 Jordan 形，出现 $t^k e^{\lambda t}$ 项；③ 数值上用 Padé 近似 + scaling-squaring（scipy.linalg.expm）。</p>
<p>由 ① 立即看出：<strong>$x(t) \to 0$ 对一切初值成立 $\iff$ 所有特征值实部 $< 0$</strong>——这是下一站《线性系统稳定性》的完整定理。</p></div>`
    },
  ],
  code: {
    title: "三种方式计算矩阵指数并对比",
    text: `import numpy as np
from scipy.linalg import expm

A = np.array([[-0.5, 2.0], [-2.0, -0.5]])   # 螺旋收敛
t = 1.0

# 方法1: 截断级数
def expm_series(M, K=60):
    S, term = np.eye(2), np.eye(2)
    for k in range(1, K):
        term = term @ M / k
        S = S + term
    return S

# 方法2: 特征分解 e^{At} = V diag(e^{lambda t}) V^{-1}
lam, V = np.linalg.eig(A * t)
expm_eig = (V @ np.diag(np.exp(lam)) @ np.linalg.inv(V)).real

# 方法3: scipy (Pade + scaling-squaring)
print(expm_series(A * t))
print(expm_eig)
print(expm(A * t))    # 三者应一致

# 特征值决定行为: lambda = -0.5 +- 2j => 衰减(实部-0.5)的振荡(频率2)
print("eigenvalues:", np.linalg.eigvals(A))`
  }
};

/* ════════════════ 传递函数 ════════════════ */
window.KG_CONTENT["transfer-function"] = {
  sim: "sim-secondorder",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>状态空间关心"内部"，传递函数只关心<strong>输入到输出的黑箱关系</strong>。它是经典控制（1930s–50s，电话放大器与火炮伺服时代）的语言，至今仍是工程师沟通的通用语。</p>
<p>核心洞察：LTI 系统对指数信号 $e^{st}$ 的响应仍是 $e^{st}$ 乘一个复数倍率 $G(s)$——指数函数是 LTI 系统的"特征函数"。把 $s$ 取遍复平面，$G(s)$ 就完整刻画了系统；取 $s = j\omega$ 得频率响应，就是 Bode 图。</p>`
    },
    {
      icon: "📐", title: "定义",
      html: String.raw`
<div class="box def"><span class="box-title">定义（传递函数）</span>
零初始条件下，输出与输入的 Laplace 变换之比：
$$G(s) = \frac{Y(s)}{U(s)} = \frac{b_0 s^m + \cdots + b_m}{s^n + a_1 s^{n-1} + \cdots + a_n}$$
分母根称为<strong>极点</strong>，分子根称为<strong>零点</strong>；$n \ge m$（真分式）对应物理可实现系统。</div>
<div class="box thm"><span class="box-title">命题（与状态空间的关系）</span>
对 $\dot{x} = Ax+Bu,\ y = Cx+Du$：$\quad G(s) = C(sI-A)^{-1}B + D$。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>对方程两边做 Laplace 变换（$x(0)=0$）：$sX(s) = AX(s) + BU(s)$，即 $(sI-A)X(s) = BU(s)$。当 $s$ 不是 $A$ 的特征值时 $sI-A$ 可逆：$X(s) = (sI-A)^{-1}BU(s)$，代入输出方程得 $Y(s) = [C(sI-A)^{-1}B + D]\,U(s)$。<span class="qed">∎</span></p>
<p>由 $(sI-A)^{-1} = \dfrac{\mathrm{adj}(sI-A)}{\det(sI-A)}$ 可见：<strong>$G(s)$ 的极点是 $A$ 的特征值</strong>（无零极对消时）。频域与时域在此汇合。</p></div>`
    },
    {
      icon: "🧾", title: "定理：正弦稳态响应（频率响应的含义）",
      html: String.raw`
<div class="box thm"><span class="box-title">定理</span>
设 $G(s)$ 稳定（极点均在左半平面）。输入 $u(t) = \sin(\omega t)$ 时，输出满足
$$y(t) \to |G(j\omega)|\sin\big(\omega t + \angle G(j\omega)\big) \quad (t\to\infty)$$
即稳态输出是同频正弦，幅值放大 $|G(j\omega)|$ 倍、相位偏移 $\angle G(j\omega)$。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>用状态空间形式（$D=0$ 不失一般性）。由常数变易公式，取复输入 $u(t)=e^{j\omega t}$：</p>
$$x(t) = e^{At}x_0 + \int_0^t e^{A(t-\tau)}Be^{j\omega\tau}d\tau$$
<p>猜特解 $x_p(t) = (j\omega I - A)^{-1}Be^{j\omega t}$（$j\omega$ 不是特征值，因为 $A$ 是 Hurwitz 的），直接代入验证：$\dot{x}_p = j\omega x_p$ 而 $Ax_p + Be^{j\omega t} = (A(j\omega I-A)^{-1} + I)Be^{j\omega t} = j\omega(j\omega I-A)^{-1}Be^{j\omega t} = \dot x_p$ ✓。</p>
<p>齐次部分 $x(t)-x_p(t) = e^{At}(x_0 - x_p(0))$，由稳定性 $e^{At}\to 0$（见《线性系统稳定性》），故 $x \to x_p$，从而</p>
$$y(t) \to Cx_p(t) = C(j\omega I-A)^{-1}Be^{j\omega t} = G(j\omega)e^{j\omega t}$$
<p>写 $G(j\omega) = |G|e^{j\varphi}$，则输出为 $|G|e^{j(\omega t+\varphi)}$。取虚部（线性系统实部虚部独立演化）即得正弦情形结论。<span class="qed">∎</span></p></div>
<div class="box eg"><span class="box-title">标准二阶系统——控制工程的"果蝇"</span>
$$G(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2},\qquad \text{极点 } s = -\zeta\omega_n \pm j\omega_n\sqrt{1-\zeta^2}$$
<p>阶跃响应超调量只由阻尼比决定：$M_p = e^{-\pi\zeta/\sqrt{1-\zeta^2}}$；调节时间 $t_s \approx \dfrac{4}{\zeta\omega_n}$。这两个公式是所有"性能指标 ↔ 极点位置"直觉的来源，请在下方实验中亲手验证。</p></div>`
    },
  ],
  code: {
    title: "由传递函数计算频率响应与阶跃响应",
    text: `import numpy as np
import matplotlib.pyplot as plt

# G(s) = wn^2 / (s^2 + 2*zeta*wn*s + wn^2)
zeta, wn = 0.3, 2.0
num = [wn**2]
den = [1, 2*zeta*wn, wn**2]

def freq_resp(num, den, w):
    s = 1j * w
    return np.polyval(num, s) / np.polyval(den, s)

w = np.logspace(-1, 2, 500)
G = freq_resp(num, den, w)

fig, ax = plt.subplots(3, 1, figsize=(7, 9))
ax[0].semilogx(w, 20*np.log10(abs(G)));  ax[0].set_ylabel('|G| dB')
ax[1].semilogx(w, np.degrees(np.angle(G))); ax[1].set_ylabel('phase deg')

# 阶跃响应（转状态空间后积分）
A = np.array([[0, 1], [-den[2], -den[1]]])
B = np.array([0, 1]); C = np.array([num[0], 0])
x = np.zeros(2); h = 0.001
ts, ys = [], []
for k in range(12000):
    ts.append(k*h); ys.append(C @ x)
    dx = A @ x + B * 1.0
    x = x + h * dx
ax[2].plot(ts, ys); ax[2].axhline(1, ls='--', c='gray')
Mp = np.exp(-np.pi*zeta/np.sqrt(1-zeta**2))
print(f"理论超调 {Mp*100:.1f}%, 实测 {(max(ys)-1)*100:.1f}%")
plt.show()`
  }
};

/* ════════════════ 线性化 ════════════════ */
window.KG_CONTENT["linearization"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>真实世界几乎全是非线性的，而我们最强大的工具箱（特征值、传递函数、LQR……）全是线性的。桥梁就是线性化：<strong>在工作点附近，光滑曲面看起来像平面</strong>。</p>
<p>飞机巡航、倒立摆直立、化工稳态运行——都是"在平衡点附近小偏差运行"，线性化模型 + 线性控制器足以胜任。而线性化的合法性本身是一条严格的定理（Lyapunov 间接法），不是想当然。</p>`
    },
    {
      icon: "📐", title: "定义",
      html: String.raw`
<div class="box def"><span class="box-title">定义（平衡点与线性化）</span>
<p>$\bar{x}$ 称为 $\dot{x} = f(x, u)$ 在恒定输入 $\bar{u}$ 下的<strong>平衡点</strong>，若 $f(\bar{x}, \bar{u}) = 0$。令偏差量 $\delta x = x - \bar{x},\ \delta u = u - \bar{u}$，对 $f$ 做 Taylor 展开并舍去高阶项：</p>
$$\dot{\delta x} = A\,\delta x + B\,\delta u, \qquad A = \frac{\partial f}{\partial x}\Big|_{(\bar{x},\bar{u})},\quad B = \frac{\partial f}{\partial u}\Big|_{(\bar{x},\bar{u})}$$
<p>$A, B$ 是 Jacobian 矩阵：$A_{ij} = \partial f_i/\partial x_j$。</p></div>
<div class="box eg"><span class="box-title">例：单摆的两个平衡点</span>
<p>$\ddot{\theta} = -\frac{g}{l}\sin\theta$，状态 $x = [\theta, \dot\theta]^{\top}$。</p>
<p>下垂点 $\theta = 0$：$\sin\theta \approx \theta$，$A = \begin{bmatrix} 0 & 1 \\ -g/l & 0 \end{bmatrix}$，特征值 $\pm j\sqrt{g/l}$（临界，振荡）。</p>
<p>倒立点 $\theta = \pi$：$\sin(\pi + \delta) \approx -\delta$，$A = \begin{bmatrix} 0 & 1 \\ +g/l & 0 \end{bmatrix}$，特征值 $\pm\sqrt{g/l}$——有一个正实部，<strong>不稳定</strong>，与直觉一致。LQR 页面的倒立摆实验正是对这个不稳定平衡点做镇定。</p></div>`
    },
    {
      icon: "🧾", title: "定理：Lyapunov 间接法（线性化的合法性）",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（间接法 / 一次近似法）</span>
<p>设 $\dot{x} = f(x)$，$f(0) = 0$，$f$ 连续可微，$A = \frac{\partial f}{\partial x}(0)$。</p>
<ol>
<li>若 $A$ 的所有特征值实部 $< 0$（Hurwitz），则原<strong>非线性</strong>系统在 $0$ 处局部指数稳定；</li>
<li>若存在特征值实部 $> 0$，则原点不稳定；</li>
<li>若最大实部 $= 0$，线性化<strong>无法判定</strong>（需 Lyapunov 直接法等）。</li>
</ol></div>
<div class="box proof"><span class="box-title">证明（情形 1，完整证明）</span>
<p><strong>第一步：分解。</strong>写 $f(x) = Ax + g(x)$，其中余项 $g(x) = f(x) - Ax$ 满足（由可微性）</p>
$$\lim_{\|x\|\to 0} \frac{\|g(x)\|}{\|x\|} = 0$$
<p>即：对任意 $\varepsilon > 0$，存在 $r > 0$ 使 $\|x\| < r \Rightarrow \|g(x)\| \le \varepsilon\|x\|$。</p>
<p><strong>第二步：为线性部分构造 Lyapunov 函数。</strong>由于 $A$ 是 Hurwitz 的，取</p>
$$P = \int_0^{\infty} e^{A^{\top}t}e^{At}\,dt$$
<p>积分收敛（被积项范数以 $ce^{-2\alpha t}$ 衰减，$\alpha>0$ 为谱横坐标的余量），$P = P^{\top} > 0$（对 $x \ne 0$，$x^{\top}Px = \int_0^\infty \|e^{At}x\|^2 dt > 0$）。验证它解 Lyapunov 方程：</p>
$$A^{\top}P + PA = \int_0^{\infty} \frac{d}{dt}\left(e^{A^{\top}t}e^{At}\right)dt = \left[e^{A^{\top}t}e^{At}\right]_0^{\infty} = -I$$
<p><strong>第三步：用在非线性系统上。</strong>取 $V(x) = x^{\top}Px$，沿非线性轨迹求导：</p>
$$\dot{V} = x^{\top}(A^{\top}P + PA)x + 2x^{\top}Pg(x) = -\|x\|^2 + 2x^{\top}Pg(x)$$
<p>在球 $\|x\| < r$ 内，$|2x^{\top}Pg(x)| \le 2\|P\|\,\varepsilon\|x\|^2$。选 $\varepsilon = \frac{1}{4\|P\|}$（对应某个 $r$），则</p>
$$\dot{V} \le -\|x\|^2 + \tfrac{1}{2}\|x\|^2 = -\tfrac{1}{2}\|x\|^2 \le -\frac{1}{2\lambda_{\max}(P)}V$$
<p>由比较引理 $V(t) \le V(0)e^{-t/(2\lambda_{\max}(P))}$，再由 $\lambda_{\min}(P)\|x\|^2 \le V$ 得</p>
$$\|x(t)\| \le \sqrt{\tfrac{\lambda_{\max}(P)}{\lambda_{\min}(P)}}\; e^{-t/(4\lambda_{\max}(P))}\,\|x(0)\|$$
<p>即局部指数稳定。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">要点</span>
这个证明是整张知识图谱的一个枢纽：它把《矩阵指数》（构造 $P$）、《Lyapunov 稳定性》（能量函数思想）和工程实践（"线性控制器用在真实非线性系统上为什么有效"）连在了一起。注意结论只是<strong>局部</strong>的——吸引域多大，线性化不告诉你。</div>`
    },
  ],
  code: {
    title: "数值线性化：自动求 Jacobian 并验证间接法",
    text: `import numpy as np

def numeric_jacobian(f, x0, eps=1e-6):
    n = len(x0); fx = f(x0)
    J = np.zeros((len(fx), n))
    for j in range(n):
        dx = np.zeros(n); dx[j] = eps
        J[:, j] = (f(x0 + dx) - f(x0 - dx)) / (2 * eps)
    return J

g, l = 9.81, 1.0
def pendulum(x):           # 无阻尼单摆
    return np.array([x[1], -g/l * np.sin(x[0])])

for name, eq in [("下垂点 theta=0", np.array([0.0, 0.0])),
                 ("倒立点 theta=pi", np.array([np.pi, 0.0]))]:
    A = numeric_jacobian(pendulum, eq)
    lam = np.linalg.eigvals(A)
    verdict = ("渐近稳定" if lam.real.max() < -1e-9 else
               "不稳定" if lam.real.max() > 1e-9 else "临界(无法判定)")
    print(f"{name}: eig = {np.round(lam, 3)} -> {verdict}")`
  }
};
