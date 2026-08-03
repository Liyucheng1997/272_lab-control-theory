/* 内容：稳定性理论 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ 线性系统稳定性 ════════════════ */
window.KG_CONTENT["linear-stability"] = {
  sim: "sim-eig-phase",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>稳定性是控制的第一要务：一个不稳定的系统，谈任何性能都是空话。对线性系统，答案出奇地干净——<strong>看特征值就够了</strong>。</p>
<p>由《矩阵指数》可知解由模态 $e^{\lambda_i t}$ 叠加而成。$\mathrm{Re}\,\lambda < 0$ 的模态衰减，$\mathrm{Re}\,\lambda > 0$ 的爆炸，虚轴上的不增不减。于是整个系统的命运由"最右边的特征值"一票决定。</p>`
    },
    {
      icon: "🧾", title: "主定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（特征值判据）</span>
对 $\dot{x} = Ax$：
<ol>
<li>原点渐近稳定（且必为指数稳定）$\iff$ 所有特征值满足 $\mathrm{Re}\,\lambda_i(A) < 0$（称 $A$ 为 <strong>Hurwitz 矩阵</strong>）；</li>
<li>原点（Lyapunov 意义下）稳定 $\iff$ $\mathrm{Re}\,\lambda_i \le 0$ 且实部为零的特征值对应的 Jordan 块均为 $1\times 1$；</li>
<li>其余情形不稳定。</li>
</ol></div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>(1) 充分性。</strong>设 $\alpha = \max_i \mathrm{Re}\,\lambda_i < 0$，取 $\varepsilon > 0$ 使 $\alpha + \varepsilon < 0$。由 Jordan 分解 $A = TJT^{-1}$，$e^{At} = Te^{Jt}T^{-1}$，而 $e^{Jt}$ 的每个元素形如 $\frac{t^k}{k!}e^{\lambda t}$（$k$ 小于 Jordan 块尺寸）。由于多项式增长慢于任何指数：$t^k e^{\lambda t} = t^k e^{-\varepsilon t}\cdot e^{(\lambda+\varepsilon)t}$，而 $t^k e^{-\varepsilon t}$ 有界（设界为 $M_k$），故存在常数 $c$ 使</p>
$$\|e^{At}\| \le c\,e^{(\alpha+\varepsilon)t} \to 0$$
<p>从而 $\|x(t)\| \le c e^{(\alpha+\varepsilon)t}\|x_0\| \to 0$，指数收敛。</p>
<p><strong>(1) 必要性。</strong>若某 $\mathrm{Re}\,\lambda \ge 0$，取对应特征向量 $v$（复情形取实部处理）：$x(t) = e^{\lambda t}v$ 是一条解，$\|x(t)\| = e^{(\mathrm{Re}\lambda) t}\|v\|$ 不收敛到 0。</p>
<p><strong>(2)</strong> 虚轴特征值若 Jordan 块为 $1\times1$，对应模态 $|e^{j\omega t}| = 1$ 有界，解有界但不衰减——稳定不渐近。若块 $\ge 2\times 2$，出现 $t\,e^{j\omega t}$ 项，范数线性增长，不稳定。例：$A = \begin{bmatrix} 0 & 1 \\ 0 & 0 \end{bmatrix}$（双积分器）特征值全为 0，但 $x_1(t) = x_1(0) + t\,x_2(0)$ 无界。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "Lyapunov 方程判据（代数版）",
      html: String.raw`
<p>特征值判据需要求根；下面的判据只需解线性方程组，而且是通往非线性理论的桥。</p>
<div class="box thm"><span class="box-title">定理（Lyapunov 方程）</span>
$A$ 是 Hurwitz 的 $\iff$ 对任意给定的 $Q = Q^{\top} \succ 0$，Lyapunov 方程
$$A^{\top}P + PA = -Q$$
存在唯一解 $P = P^{\top} \succ 0$。</div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>(⇒) 存在性。</strong>令 $P = \displaystyle\int_0^{\infty} e^{A^{\top}t}\,Q\,e^{At}\,dt$。Hurwitz 性保证被积函数指数衰减，积分收敛；显然 $P$ 对称，且 $x^{\top}Px = \int_0^\infty (e^{At}x)^{\top}Q(e^{At}x)\,dt > 0$（$x\neq 0$ 时被积函数恒正），故 $P \succ 0$。代入验证：</p>
$$A^{\top}P + PA = \int_0^{\infty} \frac{d}{dt}\big(e^{A^{\top}t}Qe^{At}\big)\,dt = \lim_{t\to\infty}e^{A^{\top}t}Qe^{At} - Q = -Q$$
<p><strong>(⇒) 唯一性。</strong>设 $P_1, P_2$ 均为解，$\Delta = P_1 - P_2$ 满足 $A^{\top}\Delta + \Delta A = 0$。于是</p>
$$\frac{d}{dt}\big(e^{A^{\top}t}\Delta e^{At}\big) = e^{A^{\top}t}(A^{\top}\Delta + \Delta A)e^{At} = 0$$
<p>该矩阵恒等于其 $t=0$ 的值 $\Delta$；令 $t \to \infty$ 左边趋于 $0$，故 $\Delta = 0$。</p>
<p><strong>(⇐)</strong> 设方程有解 $P \succ 0$。取 $A$ 的任一特征对 $Av = \lambda v$（$v \ne 0$），左乘 $v^*$、右乘 $v$：</p>
$$v^*(A^{\top}P + PA)v = \bar{\lambda}v^*Pv + \lambda v^*Pv = 2\,\mathrm{Re}(\lambda)\, v^*Pv = -v^*Qv < 0$$
<p>而 $v^*Pv > 0$，故 $\mathrm{Re}\,\lambda < 0$。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">这一定理为什么重要</span>
它说明"二次型能量函数"$V = x^{\top}Px$ 对稳定线性系统<strong>总是存在</strong>（converse Lyapunov 定理的线性版），是《线性化》中间接法证明的核心构件，也是 LQR 稳定性证明、鲁棒性分析的标准工具。</div>`
    },
  ],
  code: {
    title: "用 Lyapunov 方程判稳（不求特征值）",
    text: `import numpy as np
from scipy.linalg import solve_lyapunov  # 解 A P + P A^T = Q 形式

def is_hurwitz_via_lyapunov(A):
    n = A.shape[0]
    Q = -np.eye(n)
    # scipy 解的是 A X + X A^H = Q，对应我们的 A^T P + P A = -I 需转置
    P = solve_lyapunov(A.T, Q)
    # 检查 P 对称正定
    sym_ok = np.allclose(P, P.T, atol=1e-9)
    eigs = np.linalg.eigvalsh((P + P.T) / 2)
    return sym_ok and eigs.min() > 0, P

A_stable   = np.array([[-1.0, 2.0], [0.0, -3.0]])
A_unstable = np.array([[0.1, 1.0], [-1.0, -0.5]])

for A in (A_stable, A_unstable):
    ok, P = is_hurwitz_via_lyapunov(A)
    print("A =", A.tolist())
    print("  Lyapunov 判稳:", "Hurwitz" if ok else "不稳定",
          "| 特征值验证:", np.round(np.linalg.eigvals(A), 3))`
  }
};

/* ════════════════ Routh–Hurwitz ════════════════ */
window.KG_CONTENT["routh-hurwitz"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>1868 年 Maxwell 研究蒸汽机调速器时提出问题：不解特征方程（五次以上没有求根公式！），能否判断所有根都在左半平面？Routh（1877）与 Hurwitz（1895）分别给出了答案——只需对系数做加减乘除。</p>
<p>在计算器时代之前，这是判稳的唯一实用手段；今天它仍然珍贵，因为它能<strong>带着符号参数做判断</strong>：比如"增益 $K$ 取多大范围闭环稳定"，数值求根给不出这种答案。</p>`
    },
    {
      icon: "📐", title: "Routh 表与判据",
      html: String.raw`
<div class="box def"><span class="box-title">Routh 表的构造</span>
<p>对特征多项式 $p(s) = a_0 s^n + a_1 s^{n-1} + \cdots + a_n\ (a_0 > 0)$，前两行填入隔项系数：</p>
$$\begin{array}{c|cccc}
s^n & a_0 & a_2 & a_4 & \cdots\\
s^{n-1} & a_1 & a_3 & a_5 & \cdots\\
s^{n-2} & b_1 & b_2 & \cdots\\
\vdots & \vdots\\
s^0 & *
\end{array}
\qquad b_1 = \frac{a_1 a_2 - a_0 a_3}{a_1},\quad b_2 = \frac{a_1 a_4 - a_0 a_5}{a_1},\ \dots$$
<p>每一行由上面两行按同样的"交叉相乘"规则生成。</p></div>
<div class="box thm"><span class="box-title">定理（Routh–Hurwitz 判据）</span>
$p(s)$ 的所有根都在开左半平面 $\iff$ Routh 表<strong>第一列元素全为正</strong>。且第一列的符号变化次数等于右半平面根的个数。</div>`
    },
    {
      icon: "🧾", title: "证明（必要条件 + 低阶完整证明 + 一般情形思路）",
      html: String.raw`
<div class="box thm"><span class="box-title">引理（必要条件：系数全正）</span>
若实多项式 $p(s)$ 的根全部在开左半平面，则所有系数同号（不失一般性全为正）。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>实多项式可分解为实根与共轭复根对的乘积：</p>
$$p(s) = a_0 \prod_i (s + p_i) \prod_k \big(s^2 + 2\sigma_k s + \sigma_k^2 + \omega_k^2\big)$$
<p>根在左半平面 $\Rightarrow p_i > 0,\ \sigma_k > 0$，因此每个因子的系数都是正数。正系数多项式相乘系数仍为正。<span class="qed">∎</span></p>
<p>（注意反过来不成立：$s^3+s^2+s+1$ 系数全正，但有根 $\pm j$ 在虚轴上。）</p></div>
<div class="box proof"><span class="box-title">二阶与三阶的完整证明</span>
<p><strong>二阶</strong> $s^2 + a_1 s + a_2$：根为 $\frac{-a_1 \pm \sqrt{a_1^2 - 4a_2}}{2}$。若 $a_1, a_2 > 0$：判别式为负时实部 $= -a_1/2 < 0$；为正时 $\sqrt{a_1^2-4a_2} < a_1$，两根均负。反之若某根实部 $\ge 0$，由根与系数关系（$a_1 = -\sum$ 根、$a_2 = \prod$ 根）可推出 $a_1 \le 0$ 或 $a_2 \le 0$。故二阶时"系数全正 $\iff$ 稳定"，与 Routh 表第一列 $(1, a_1, a_2)$ 全正一致。✓</p>
<p><strong>三阶</strong> $s^3 + a_1 s^2 + a_2 s + a_3$：Routh 表第一列为 $1,\ a_1,\ \frac{a_1 a_2 - a_3}{a_1},\ a_3$，判据给出附加条件 $a_1 a_2 > a_3$。验证必要性：设根为 $-p,\ -\sigma \pm j\omega$（$p, \sigma > 0$ 为稳定情形），展开得</p>
$$a_1 = p + 2\sigma,\quad a_2 = 2p\sigma + \sigma^2 + \omega^2,\quad a_3 = p(\sigma^2 + \omega^2)$$
$$a_1 a_2 - a_3 = (p+2\sigma)(2p\sigma + \sigma^2+\omega^2) - p(\sigma^2+\omega^2) = 2\sigma\big[(p+\sigma)^2 + \omega^2\big] > 0 \ ✓$$
<p>且当一对复根移动到虚轴（$\sigma = 0$）时恰有 $a_1a_2 - a_3 = 0$：这个量正是"复根穿越虚轴"的检测器。三实根情形同理可验证。充分性由上式反解同样成立。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">一般情形的证明思路（Routh 的原始路线）</span>
<p>完整证明依赖 <strong>Cauchy 辐角指标 / Sturm 链</strong>：把 $p(j\omega)$ 分解为偶部 $E(\omega)$ 与奇部 $O(\omega)$，右半平面根数可表示为有理函数 $E/O$ 在 $\omega: -\infty \to \infty$ 上的 Cauchy 指标；对 $E, O$ 做欧几里得辗转相除生成 Sturm 链，由 Sturm 定理，该指标等于链首项系数序列的符号变化数——而 Routh 表每一行恰好就是辗转相除的中间余式系数，第一列就是首项系数序列。符号变化数 = 右半平面根数，全正即无右根。</p>
<p>现代替代证法：Hermite–Biehler 定理（稳定 $\iff$ $E,O$ 的根实、互相交错）或与 Lyapunov 方程的联系。</p></div>
<div class="box eg"><span class="box-title">经典应用：求临界增益</span>
<p>单位反馈，开环 $L(s) = \dfrac{K}{s(s+1)(s+2)}$，闭环特征式 $s^3 + 3s^2 + 2s + K$。</p>
<p>由三阶判据：$K > 0$ 且 $3 \times 2 > K$，即 $0 < K < 6$ 时稳定，$K = 6$ 时临界振荡（频率 $\omega = \sqrt{2}$）。这类含参结论正是 Routh 判据的独门优势，也与《根轨迹》《Nyquist》的图形结论互相印证。</p></div>`
    },
  ],
  code: {
    title: "Routh 表自动生成与判稳",
    text: `import numpy as np

def routh_table(coeffs):
    """coeffs: [a0, a1, ..., an] 从最高次到常数项"""
    n = len(coeffs) - 1
    cols = (n // 2) + 1
    T = np.zeros((n + 1, cols))
    T[0, :len(coeffs[0::2])] = coeffs[0::2]
    T[1, :len(coeffs[1::2])] = coeffs[1::2]
    for i in range(2, n + 1):
        for j in range(cols - 1):
            a, b = T[i-2, 0], T[i-2, j+1]
            c, d = T[i-1, 0], T[i-1, j+1]
            if abs(c) < 1e-12:
                c = 1e-9          # epsilon 法处理首列零
            T[i, j] = (c * b - a * d) / c
    return T

def n_rhp_roots(coeffs):
    first_col = routh_table(coeffs)[:, 0]
    signs = np.sign(first_col[np.abs(first_col) > 1e-9])
    return int(np.sum(signs[:-1] != signs[1:]))

# 例: s^3 + 3 s^2 + 2 s + K, 理论临界 K=6
for K in [1, 5.9, 6.1, 20]:
    c = [1, 3, 2, K]
    nr = n_rhp_roots(c)
    print(f"K={K:5}: 右半平面根 {nr} 个 ->",
          "稳定" if nr == 0 else "不稳定",
          "| 验证:", np.round(np.roots(c).real.max(), 3))`
  }
};

/* ════════════════ Lyapunov 稳定性 ════════════════ */
window.KG_CONTENT["lyapunov"] = {
  sim: "sim-lyapunov",
  sections: [
    {
      icon: "💡", title: "直觉：不解方程，看能量",
      html: String.raw`
<p>非线性微分方程绝大多数解不出来。Lyapunov（1892）的天才想法是：<strong>不需要解！</strong>只要能找到一个类似"能量"的函数 $V(x)$——正定、且沿轨迹不增——就能断言稳定。</p>
<p>想象一个碗里的小球：碗的高度函数就是 $V$。小球怎么滚我们不必知道细节，只要确认"高度一直在降低"，小球必然滚向碗底。$V$ 的等高线是一圈圈套住原点的闭曲线，$\dot{V} < 0$ 说明轨迹只能不断穿入更小的等高线，永远出不去。</p>
<p>这一手"绕过求解、直接下结论"的思想，是<strong>整个非线性控制的地基</strong>：反步法在递推构造 $V$，滑模在选取 $V = \frac{1}{2}s^2$，自适应控制在 $V$ 里塞进参数误差，MPC 的稳定性证明把值函数当 $V$……它是本图谱后半部所有节点的公共祖先。</p>`
    },
    {
      icon: "📐", title: "稳定性的严格定义",
      html: String.raw`
<p>考虑 $\dot{x} = f(x)$，$f(0) = 0$（原点是平衡点），解记为 $x(t; x_0)$。</p>
<div class="box def"><span class="box-title">定义（Lyapunov 稳定 / 渐近稳定 / 指数稳定）</span>
<ul>
<li><strong>稳定</strong>：$\forall \varepsilon > 0,\ \exists \delta > 0$：$\|x_0\| < \delta \Rightarrow \|x(t;x_0)\| < \varepsilon,\ \forall t \ge 0$。（想多近就能保持多近）</li>
<li><strong>渐近稳定</strong>：稳定，且存在 $\delta_0$：$\|x_0\| < \delta_0 \Rightarrow x(t;x_0) \to 0$。</li>
<li><strong>全局渐近稳定 (GAS)</strong>：渐近稳定且吸引域为全空间。</li>
<li><strong>指数稳定</strong>：$\exists c, \alpha, \delta$：$\|x(t)\| \le c e^{-\alpha t}\|x_0\|$。</li>
</ul></div>
<div class="box note"><span class="box-title">易错点</span>
"所有轨迹都趋于原点"并<strong>不</strong>蕴含稳定——存在著名反例（Vinograd 系统），轨迹先跑得很远再回来，破坏 $\varepsilon$-$\delta$ 条件。所以渐近稳定的定义必须"稳定 + 吸引"两条同时要。</div>
<div class="box def"><span class="box-title">定义（正定函数）</span>
连续可微 $V: D \to \mathbb{R}$（$0 \in D$）称为<strong>正定</strong>，若 $V(0)=0$ 且 $V(x) > 0\ (x \ne 0)$；称<strong>径向无界</strong>，若 $\|x\|\to\infty \Rightarrow V(x) \to \infty$。沿系统轨迹的导数为
$$\dot{V}(x) = \nabla V(x)^{\top} f(x)$$
（不需要解出轨迹，代入向量场即可计算——这是方法可行的关键。）</div>`
    },
    {
      icon: "🧾", title: "Lyapunov 直接法：定理与完整证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（Lyapunov, 1892）</span>
设存在正定的 $V \in C^1(D)$。
<ol>
<li>若 $\dot{V}(x) \le 0$ 于 $D$，则原点<strong>稳定</strong>；</li>
<li>若 $\dot{V}(x) < 0$ 于 $D\setminus\{0\}$，则原点<strong>渐近稳定</strong>；</li>
<li>若进一步 $D = \mathbb{R}^n$ 且 $V$ 径向无界，则<strong>全局</strong>渐近稳定。</li>
</ol></div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>(1) 稳定性。</strong>给定 $\varepsilon > 0$（不妨小到 $\bar{B}_\varepsilon = \{\|x\| \le \varepsilon\} \subset D$）。令</p>
$$m = \min_{\|x\| = \varepsilon} V(x)$$
<p>球面紧、$V$ 连续，最小值可达；$V$ 正定故 $m > 0$。由 $V$ 连续且 $V(0) = 0$，存在 $\delta \in (0, \varepsilon)$ 使 $\|x\| < \delta \Rightarrow V(x) < m$。</p>
<p>现取 $\|x_0\| < \delta$。沿轨迹 $\frac{d}{dt}V(x(t)) = \dot{V} \le 0$，故 $V(x(t)) \le V(x_0) < m$ 对一切 $t$ 成立。假若轨迹某时刻触及球面 $\|x\| = \varepsilon$，则该点处 $V \ge m$，矛盾。由连续性轨迹无法跳过球面，故 $\|x(t)\| < \varepsilon$ 恒成立（这同时保证解不会有限时间逃逸，全局存在）。稳定性得证。</p>
<p><strong>(2) 吸引性。</strong>由 (1)，取定一个 $\varepsilon$ 及相应 $\delta$，考虑 $\|x_0\| < \delta$ 的轨迹，它停留在紧集 $\bar{B}_\varepsilon$ 内。$V(x(t))$ 单调不增且下有界（$\ge 0$），故收敛：$V(x(t)) \downarrow c \ge 0$。要证 $c = 0$。</p>
<p>反设 $c > 0$。由 $V$ 连续、$V(0)=0$，存在 $r > 0$ 使 $\|x\| < r \Rightarrow V(x) < c$；由于 $V(x(t)) \ge c$，轨迹被困在环形紧集 $K = \{x : r \le \|x\| \le \varepsilon\}$ 中。令</p>
$$\gamma = \max_{x \in K} \dot{V}(x)$$
<p>$\dot V$ 连续、$K$ 紧且不含原点，故 $\gamma < 0$。于是</p>
$$V(x(t)) = V(x_0) + \int_0^t \dot{V}(x(\tau))\,d\tau \le V(x_0) + \gamma t \xrightarrow{t\to\infty} -\infty$$
<p>与 $V \ge 0$ 矛盾。故 $c = 0$，再由 $V$ 正定（在紧集上，$V(x)\to 0$ 迫使 $x \to 0$：否则有子列 $\|x(t_k)\| \ge \eta$，紧性取收敛子列极限点 $\bar x \ne 0$，连续性给出 $V(\bar x) = 0$，与正定矛盾），得 $x(t) \to 0$。</p>
<p><strong>(3) 全局性。</strong>任取 $x_0 \in \mathbb{R}^n$。径向无界性保证水平集 $\Omega = \{x: V(x) \le V(x_0)\}$ 有界（否则有 $\|x_k\| \to \infty$ 而 $V(x_k) \le V(x_0)$，矛盾），且闭，故紧。$\dot V \le 0$ 使 $\Omega$ 正不变，轨迹困于紧集，(2) 的论证逐字适用。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理（指数稳定判据）</span>
若存在 $c_1\|x\|^2 \le V(x) \le c_2\|x\|^2$ 且 $\dot{V} \le -c_3\|x\|^2$，则原点指数稳定，衰减率 $\ge c_3/(2c_2)$。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>$\dot{V} \le -c_3\|x\|^2 \le -\frac{c_3}{c_2}V$。由比较引理（Grönwall）：$V(t) \le V(0)e^{-(c_3/c_2)t}$。于是</p>
$$\|x(t)\|^2 \le \frac{V(t)}{c_1} \le \frac{c_2}{c_1}\|x_0\|^2 e^{-(c_3/c_2)t} \ \Rightarrow\ \|x(t)\| \le \sqrt{\tfrac{c_2}{c_1}}\,\|x_0\|\,e^{-\frac{c_3}{2c_2}t} \qquad ∎$$</div>
<div class="box note"><span class="box-title">局限与出路</span>
<ul>
<li>定理只给<strong>充分</strong>条件：找不到 $V$ 不代表不稳定。但 converse 定理保证：渐近稳定系统必存在 Lyapunov 函数（只是没有构造方法）。</li>
<li>找 $V$ 是艺术：物理能量、二次型 $x^{\top}Px$（线性系统必够用，见《线性系统稳定性》）、变量分组（backstepping）是三大常用套路。</li>
<li>很多实际问题中 $\dot{V}$ 只能做到 $\le 0$（半负定），如阻尼摆的能量在 $\dot\theta = 0$ 瞬间不减。这时渐近稳定性的结论由《LaSalle 不变原理》补上——请前往下一站。</li>
</ul></div>`
    },
    {
      icon: "🔬", title: "例：阻尼单摆（与右侧实验对应）",
      html: String.raw`
<div class="box eg"><span class="box-title">完整分析</span>
<p>$\ddot\theta = -\sin\theta - c\dot\theta$（$c>0$），取机械能</p>
$$V(\theta, \dot\theta) = \tfrac{1}{2}\dot\theta^2 + (1 - \cos\theta)$$
<p>$V$ 在 $(0,0)$ 附近正定（$1-\cos\theta > 0$ 当 $0 < |\theta| < 2\pi$）。求导：</p>
$$\dot{V} = \dot\theta\ddot\theta + \sin\theta\,\dot\theta = \dot\theta(-\sin\theta - c\dot\theta) + \sin\theta\,\dot\theta = -c\dot\theta^2 \le 0$$
<p>由定理 (1)：原点稳定。但 $\dot{V}$ 在整条直线 $\dot\theta = 0$ 上为零（不是负定），定理 (2) 用不上——而实验里轨迹明明收敛到了原点。理论与现象的缺口，正是 LaSalle 不变原理要填的。</p></div>`
    },
  ],
  code: {
    title: "数值验证 Lyapunov 条件（沿大量随机轨迹检查 V 单调性）",
    text: `import numpy as np

c = 0.25
def f(x):     # 阻尼单摆
    return np.array([x[1], -np.sin(x[0]) - c * x[1]])

def V(x):
    return 0.5 * x[1]**2 + (1 - np.cos(x[0]))

def Vdot(x):  # 解析: -c * omega^2
    grad = np.array([np.sin(x[0]), x[1]])
    return grad @ f(x)

rng = np.random.default_rng(0)
violations = 0
for _ in range(200):
    x = rng.uniform([-2.5, -2.0], [2.5, 2.0])
    for _ in range(4000):
        if Vdot(x) > 1e-12:
            violations += 1
            break
        # RK4 一步
        h = 0.01
        k1=f(x); k2=f(x+h/2*k1); k3=f(x+h/2*k2); k4=f(x+h*k3)
        x = x + h/6*(k1+2*k2+2*k3+k4)
print("Vdot > 0 的轨迹数:", violations, "/ 200  (理论应为 0)")
print("示例: V 从", round(V(np.array([2.0, 1.0])), 4), "沿轨迹单调下降")`
  }
};

/* ════════════════ LaSalle ════════════════ */
window.KG_CONTENT["lasalle"] = {
  sim: "sim-lyapunov",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>阻尼摆的能量导数 $\dot{V} = -c\dot\theta^2$ 在 $\dot\theta = 0$ 处为零——摆到最高点的瞬间能量不减。Lyapunov 定理因此卡壳。但物理直觉说：摆<strong>不可能停在最高点</strong>，因为那一瞬间重力立刻把它拽回运动，能量继续流失。</p>
<p>LaSalle（1960）把这个直觉变成定理：轨迹最终只能停留在"$\dot{V}=0$ 集合中<strong>能一直待下去</strong>的部分"（最大不变子集）。$\dot\theta = 0$ 的点大多"待不住"——只有平衡点本身待得住。于是收敛结论得救。</p>`
    },
    {
      icon: "📐", title: "预备概念",
      html: String.raw`
<div class="box def"><span class="box-title">定义（不变集）</span>
集合 $M$ 称为<strong>不变集</strong>，若 $x_0 \in M \Rightarrow x(t; x_0) \in M,\ \forall t \in \mathbb{R}$（正不变：仅要求 $t \ge 0$）。</div>
<div class="box def"><span class="box-title">定义（$\omega$-极限集）</span>
$$\omega(x_0) = \{p : \exists t_k \to \infty,\ x(t_k; x_0) \to p\}$$
即轨迹"最终无限次接近"的点的全体。</div>
<div class="box thm"><span class="box-title">引理（$\omega$-极限集的性质）</span>
若正半轨迹停留在紧集 $K$ 内，则 $\omega(x_0)$ 非空、紧、<strong>不变</strong>，且 $x(t) \to \omega(x_0)$（到集合的距离趋于 0）。</div>
<div class="box proof"><span class="box-title">引理证明</span>
<p><strong>非空</strong>：序列 $x(k)\ (k\in\mathbb{N})$ 落在紧集 $K$ 中，必有收敛子列，其极限属于 $\omega(x_0)$。</p>
<p><strong>紧</strong>：$\omega(x_0) \subset K$ 有界；再证闭：设 $p_i \in \omega(x_0),\ p_i \to p$。对每个 $i$ 取 $t_i > i$ 使 $\|x(t_i) - p_i\| < 1/i$，则 $x(t_i) \to p$ 且 $t_i \to \infty$，故 $p \in \omega(x_0)$。</p>
<p><strong>不变</strong>：设 $p \in \omega(x_0)$，$x(t_k) \to p$。对任意固定 $s$，由解对初值的连续依赖性（在紧集上一致），$x(t_k + s) = x(s;\, x(t_k)) \to x(s; p)$，而 $t_k + s \to \infty$，故 $x(s; p) \in \omega(x_0)$。即从 $p$ 出发的整条轨迹都在 $\omega(x_0)$ 内。</p>
<p><strong>趋近</strong>：反设存在 $\varepsilon > 0$ 与 $t_k \to \infty$ 使 $\mathrm{dist}(x(t_k), \omega(x_0)) \ge \varepsilon$。紧性取 $x(t_k)$ 的收敛子列，其极限一方面按定义属于 $\omega(x_0)$，另一方面与 $\omega(x_0)$ 距离 $\ge \varepsilon$，矛盾。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "LaSalle 不变原理：定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（LaSalle 不变原理）</span>
<p>设 $\Omega$ 为紧的正不变集，$V \in C^1$ 满足 $\dot{V}(x) \le 0$ 于 $\Omega$。定义</p>
$$E = \{x \in \Omega : \dot{V}(x) = 0\},\qquad M = E \text{ 中的最大不变集}$$
<p>则从 $\Omega$ 出发的每条轨迹满足 $x(t) \to M$（$t \to \infty$）。</p>
<p><em>注：这里 $V$ 不要求正定——只要求单调，这也是 LaSalle 比 Lyapunov 更灵活之处。</em></p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>取 $x_0 \in \Omega$。轨迹留在紧集 $\Omega$ 内，$V(x(t))$ 单调不增且在紧集上有下界，故 $V(x(t)) \to c$（某常数）。</p>
<p><strong>第一步：$V$ 在 $\omega(x_0)$ 上恒等于 $c$。</strong>任取 $p \in \omega(x_0)$，有 $t_k \to \infty$，$x(t_k) \to p$；由 $V$ 连续，$V(p) = \lim V(x(t_k)) = c$。</p>
<p><strong>第二步：$\omega(x_0) \subset E$。</strong>由引理 $\omega(x_0)$ 不变：从任一 $p \in \omega(x_0)$ 出发的轨迹 $x(s;p)$ 整条留在 $\omega(x_0)$ 内，于是 $V(x(s;p)) \equiv c$ 对一切 $s$ 成立。对 $s$ 求导得 $\dot{V}(p) = 0$，故 $p \in E$。</p>
<p><strong>第三步：收尾。</strong>$\omega(x_0)$ 是 $E$ 内的不变集，故 $\omega(x_0) \subseteq M$（$M$ 是最大者）。由引理的趋近性质 $x(t) \to \omega(x_0) \subseteq M$。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">推论（Barbashin–Krasovskii 定理）</span>
若 $V$ 正定、$\dot{V} \le 0$，且 <strong>$E = \{\dot{V} = 0\}$ 中除原点外不含任何完整轨迹</strong>，则原点渐近稳定。（加上径向无界则全局渐近稳定。）</div>
<div class="box proof"><span class="box-title">推论证明</span>
<p>正定 + $\dot V \le 0$ 给出稳定性与紧的正不变水平集 $\Omega = \{V \le \rho\}$（Lyapunov 定理第 (1) 部分的论证）。条件说明 $E$ 中最大不变集 $M = \{0\}$，LaSalle 给出 $x(t) \to \{0\}$。稳定 + 吸引 = 渐近稳定。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🔬", title: "例：阻尼单摆收尾 + 更多应用",
      html: String.raw`
<div class="box eg"><span class="box-title">阻尼单摆（接《Lyapunov》的例子）</span>
<p>已有 $\dot{V} = -c\dot\theta^2$，故 $E = \{(\theta, \dot\theta): \dot\theta = 0\}$。找 $E$ 中的完整轨迹：轨迹要一直满足 $\dot\theta(t) \equiv 0$，则 $\ddot\theta \equiv 0$，代入动力学 $\ddot\theta = -\sin\theta - c\dot\theta$ 得 $\sin\theta \equiv 0$，即 $\theta \in \{0, \pm\pi, \dots\}$。在考察区域 $\Omega = \{V < 2\}$（不包含最高点 $\theta = \pm\pi$，那里 $V = 2$）内，$E$ 中唯一完整轨迹是原点。由 Barbashin–Krasovskii：<strong>原点渐近稳定</strong>，吸引域至少包含 $\{V < 2\}$。理论与实验闭环。✓</p></div>
<div class="box note"><span class="box-title">LaSalle 在后续理论中的出场</span>
<ul>
<li><strong>自适应控制</strong>：$V = \frac12 e^2 + \frac{1}{2\gamma}\tilde\theta^2$ 的导数只有 $-e^2$ 项（参数误差项不出现），典型的半负定，靠 LaSalle/Barbalat 得跟踪误差收敛；</li>
<li><strong>机器人 PD+重力补偿</strong>：能量型 $V$ 的 $\dot V = -\dot q^{\top}K_d\dot q$，同样套路；</li>
<li>注意 LaSalle 需要<strong>自治系统</strong>（$f$ 不显含 $t$）；时变系统要用 Barbalat 引理替代。</li>
</ul></div>`
    },
  ],
  code: {
    title: "验证：轨迹收敛到 E 中最大不变集（而非 E 的其它部分）",
    text: `import numpy as np

c = 0.25
def f(x):
    return np.array([x[1], -np.sin(x[0]) - c*x[1]])

def rk4(x, h=0.01):
    k1=f(x); k2=f(x+h/2*k1); k3=f(x+h/2*k2); k4=f(x+h*k3)
    return x + h/6*(k1+2*k2+2*k3+k4)

# 轨迹多次经过 E = {omega=0}（摆的最高点), 但都停不住,
# 最终只收敛到 E 中的不变集: 原点
x = np.array([2.8, 0.0])
crossings = []
for k in range(60000):
    x_new = rk4(x)
    if x[1] * x_new[1] < 0:          # omega 变号: 穿过 E
        crossings.append((k*0.01, x[0]))
    x = x_new
print(f"轨迹穿过 E(omega=0) 共 {len(crossings)} 次")
print("前几次穿越时的 theta:", [round(th,3) for _,th in crossings[:5]])
print("最终状态:", np.round(x, 6), " -> 收敛到原点(E 中最大不变集)")`
  }
};
