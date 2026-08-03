/* 内容：离散化、离散控制、时滞、不确定性建模、LMI、μ 分析、ISS */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ 离散化 ════════════════ */
window.KG_CONTENT["discretization"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>今天的控制器几乎全部跑在数字计算机里：每隔采样周期 $T$ 读一次传感器、算一次、输出保持到下一拍（零阶保持 ZOH）。被控世界是连续的，控制器是离散的——两个世界的换算就是离散化。</p>
<p>关键认识：<strong>在 ZOH 假设下，线性系统的离散化是精确的，不是近似！</strong>因为输入在每拍内恒定，连续系统在采样点上的行为可以用《矩阵指数》逐拍解析积分出来。工程里常见的 Euler、Tustin 才是近似——各有各的失真方式，选错了轻则性能打折，重则把稳定系统离散成不稳定的。</p>`
    },
    {
      icon: "📐", title: "三种离散化方法",
      html: String.raw`
<div class="box def"><span class="box-title">定义（采样与零阶保持）</span>
<p>采样周期 $T$，$x_k \triangleq x(kT)$；ZOH：$u(t) = u_k$ 对 $t \in [kT, (k+1)T)$ 恒定。目标：找 $A_d, B_d$ 使 $x_{k+1} = A_dx_k + B_du_k$ 在采样点<strong>精确</strong>成立。</p></div>
<div class="box def"><span class="box-title">常用方法一览</span>
<ul>
<li><strong>ZOH 精确法</strong>：$A_d = e^{AT},\quad B_d = \left(\int_0^T e^{A\tau}d\tau\right)B$（$A$ 可逆时 $= A^{-1}(e^{AT} - I)B$）；</li>
<li><strong>前向 Euler</strong>：$A_d \approx I + AT$（把 $\dot{x}$ 换成差分）——最粗糙，$T$ 大时可能把稳定系统变不稳定；</li>
<li><strong>Tustin（双线性）</strong>：$s \leftarrow \frac{2}{T}\frac{z-1}{z+1}$——保稳定性、保频率特性形状（有频率翘曲），滤波器与控制器离散化的默认选择。</li>
</ul></div>`
    },
    {
      icon: "🧾", title: "定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理 1（ZOH 离散化的精确性）</span>
在 ZOH 输入下，$x_{k+1} = e^{AT}x_k + \left(\int_0^T e^{A\tau}d\tau\right)Bu_k$ 对采样点<strong>精确</strong>成立。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>由《矩阵指数》的常数变易公式，在区间 $[kT, (k+1)T]$ 上（该区间内 $u \equiv u_k$）：</p>
$$x((k+1)T) = e^{AT}x(kT) + \int_{kT}^{(k+1)T} e^{A((k+1)T - \sigma)}B\,u_k\,d\sigma$$
<p>换元 $\tau = (k+1)T - \sigma$（$d\sigma = -d\tau$，积分限翻转）：</p>
$$\int_{kT}^{(k+1)T} e^{A((k+1)T-\sigma)}d\sigma = \int_0^T e^{A\tau}d\tau$$
<p>提出常量 $u_k$ 即得。全程只用了"区间内输入恒定"，无任何近似。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理 2（特征值映射与稳定性保持）</span>
若 $Av = \lambda v$，则 $e^{AT}v = e^{\lambda T}v$。因此 ZOH 离散化把极点按 $z = e^{sT}$ 搬家，且
$$\mathrm{Re}\,\lambda < 0 \iff |e^{\lambda T}| = e^{(\mathrm{Re}\lambda) T} < 1$$
——左半平面精确映到单位圆内，<strong>稳定性无损保持</strong>。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>对级数逐项作用：$e^{AT}v = \sum_k \frac{T^k}{k!}A^kv = \sum_k \frac{(\lambda T)^k}{k!}v = e^{\lambda T}v$。模长：$|e^{\lambda T}| = e^{\mathrm{Re}(\lambda T)}$，与 1 比较即得等价关系。对照：前向 Euler 的映射 $z = 1 + \lambda T$ 把左半平面映成圆心 $(-1, 0)$ 半径 1 的圆——$\lambda$ 很负而 $T$ 不够小时 $|1 + \lambda T| > 1$，<strong>稳定系统被离散成发散</strong>（刚性系统的经典陷阱）。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理 3（Tustin 保稳定）</span>
双线性变换 $z = \frac{1 + sT/2}{1 - sT/2}$ 把开左半平面一一映入单位圆内。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>设 $s = \sigma + j\omega$，$\sigma < 0$。计算模平方之差：</p>
$$|1 + sT/2|^2 - |1 - sT/2|^2 = \big[(1 + \tfrac{\sigma T}{2})^2 - (1 - \tfrac{\sigma T}{2})^2\big] = 2\sigma T < 0$$
<p>（$\omega$ 项在两边相同，相消。）分子模小于分母模 $\Rightarrow |z| < 1$。反向同理，虚轴映到单位圆周（$\sigma = 0$ 时两模相等）。频率翘曲：圆周上 $\omega_d = \frac{2}{T}\tan\frac{\omega T}{2}$，低频近似恒等，接近 Nyquist 频率处被压缩。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">采样多快才够（工程守则）</span>
<ul>
<li>Shannon 采样定理给<strong>底线</strong>：$\omega_s > 2\omega_{max}$ 才不混叠；控制回路要的远不止"能重建信号"；</li>
<li>经验法则：采样频率取<strong>闭环带宽的 20~30 倍</strong>（即每个上升时间内 4~10 拍）；ZOH 本身引入约 $T/2$ 的等效延迟，相位损失 $\omega T/2$——采样太慢等于白吃一个时滞（衔接《时滞系统》）；</li>
<li>抗混叠滤波器要放在采样之前（模拟侧），数字侧补救不了混叠。</li>
</ul></div>`
    },
  ],
  code: {
    title: "三种离散化对比：精确性与稳定性陷阱",
    text: `import numpy as np
from scipy.linalg import expm

A = np.array([[0., 1.], [-100., -12.]])   # 刚性: 极点 -1.1, -10.9... 实际 -50内
B = np.array([[0.], [1.]])
print("连续极点:", np.round(np.linalg.eigvals(A), 2))

for T in [0.01, 0.1, 0.25]:
    # ZOH 精确
    n = A.shape[0]
    M = np.zeros((n+1, n+1)); M[:n, :n] = A*T; M[:n, n:] = B*T
    E = expm(M)
    Ad, Bd = E[:n, :n], E[:n, n:]
    # 前向 Euler
    Ad_e = np.eye(n) + A*T
    # 稳定性对比
    r_zoh = np.abs(np.linalg.eigvals(Ad)).max()
    r_eul = np.abs(np.linalg.eigvals(Ad_e)).max()
    print(f"T={T:5}: |z|max  ZOH={r_zoh:.3f} ({'稳' if r_zoh<1 else '不稳'})"
          f"  Euler={r_eul:.3f} ({'稳' if r_eul<1 else '不稳!'})")
# ZOH 永远保稳; Euler 在 T 大时把稳定系统离散成发散`
  }
};

/* ════════════════ 离散控制 ════════════════ */
window.KG_CONTENT["discrete-control"] = {
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>拿到离散模型后，整套控制理论在 z 域"平行重演"一遍：Laplace 换成 z 变换、虚轴换成单位圆、$e^{st}$ 模态换成 $z^k$ 模态。绝大多数定理逐条平移——但离散世界还有连续世界<strong>没有的超能力</strong>：把所有极点配到 $z = 0$，系统能在<strong>有限拍内精确归零</strong>（deadbeat，有限拍控制）。连续系统做不到有限时间线性镇定，这是采样带来的独特红利。</p>`
    },
    {
      icon: "📐", title: "z 变换与稳定域",
      html: String.raw`
<div class="box def"><span class="box-title">定义（z 变换与离散传函）</span>
$$X(z) = \sum_{k=0}^{\infty} x_kz^{-k},\qquad G(z) = C(zI - A_d)^{-1}B_d + D$$
<p>$z^{-1}$ 就是"延迟一拍"算子——数字控制器 $u_k = -a_1u_{k-1} + b_0e_k + b_1e_{k-1}$ 与其传函一一对应，代码直接照写。</p></div>
<div class="box thm"><span class="box-title">定理（离散稳定性判据）</span>
$x_{k+1} = A_dx_k$ 渐近稳定 $\iff$ 所有特征值 $|\lambda_i(A_d)| < 1$（单位圆内）。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>与连续情形平行：Jordan 分解下解的模态为 $k^m\lambda^k$，而 $|k^m\lambda^k| = k^m|\lambda|^k \to 0 \iff |\lambda| < 1$（多项式增长败给几何衰减；$|\lambda|=1$ 且有非平凡 Jordan 块时 $k|\lambda|^k = k$ 发散）。必要性取对应特征向量方向的解即可。Lyapunov 版本同样平移：$A_d^{\top}PA_d - P = -Q$ 有正定解 $\iff$ Schur 稳定（证明把积分换级数 $P = \sum_k (A_d^{\top})^kQA_d^k$，逐步同构于连续证明）。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "有限拍（Deadbeat）控制：离散独有的魔法",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（Deadbeat 控制）</span>
<p>设 $(A_d, B_d)$ 能控（单输入，$n$ 维）。取 $K$ 把 $A_d - B_dK$ 的特征多项式配置为 $z^n$（全部极点在原点），则对<strong>任意</strong>初值：</p>
$$x_k = 0,\quad \forall k \ge n$$
<p>——恰好 $n$ 拍精确归零，不是渐近趋近，是<strong>精确到达</strong>。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>由《极点配置》定理（离散版逐字成立），能控性保证存在 $K$ 使 $\det(zI - A_{cl}) = z^n$，$A_{cl} \triangleq A_d - B_dK$。由 <strong>Cayley–Hamilton 定理</strong>，矩阵满足自身特征方程：</p>
$$A_{cl}^n = 0$$
<p>即 $A_{cl}$ 幂零。闭环解 $x_k = A_{cl}^kx_0$，故 $k \ge n$ 时 $x_k = A_{cl}^{k-n}A_{cl}^nx_0 = 0$。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">代价与工程取舍</span>
<ul>
<li>连续系统为何做不到？$e^{A_{cl}t}$ 是解析函数，除非恒为零不可能在有限时间变零；离散的 $A_{cl}^k$ 是矩阵幂，幂零即可；</li>
<li>Deadbeat 的账单：$n$ 拍归零要求的控制量随 $T$ 缩小按 $\sim 1/T^n$ 暴涨，且采样点之间可能有<strong>纹波</strong>（intersample ripple）——采样点上是零、拍间在振；工程常配置"温和版"极点（如 $z = 0.2 \sim 0.5$）；</li>
<li>其余设计法平移自连续：离散 LQR（DARE，见《Riccati》）、离散观测器、w 平面上用 Bode/Nyquist；Jury 判据是 Routh 的 z 域对应物；</li>
<li>MPC 天生离散（见《MPC》）——这条边在图谱上就是这么连的。</li>
</ul></div>`
    },
  ],
  code: {
    title: "Deadbeat：双积分器 2 拍精确归零（对比连续 LQR 的渐近收敛）",
    text: `import numpy as np

T = 0.1
Ad = np.array([[1, T], [0, 1]])          # 离散双积分器
Bd = np.array([[T*T/2], [T]])

# 配置特征式为 z^2: 用 Ackermann
Cm = np.hstack([Bd, Ad @ Bd])
pstar = [1, 0, 0]                        # z^2 + 0 z + 0
PA = sum(c * np.linalg.matrix_power(Ad, 2 - i) for i, c in enumerate(pstar))
K = (np.array([0., 1.]) @ np.linalg.inv(Cm) @ PA).reshape(1, -1)
Acl = Ad - Bd @ K
print("K =", np.round(K, 3), " 闭环特征值:", np.linalg.eigvals(Acl))
print("A_cl^2 =", np.round(np.linalg.matrix_power(Acl, 2), 12), "(幂零)")

x = np.array([1.0, -0.5])
for k in range(5):
    u = float(-K @ x)
    print(f"k={k}: x = {np.round(x, 6)}, u = {u:8.2f}")
    x = Acl @ x
# k=2 起精确为零 —— 连续线性系统永远做不到`
  }
};

/* ════════════════ 时滞系统 ════════════════ */
window.KG_CONTENT["time-delay"] = {
  sim: "sim-delay",
  sections: [
    {
      icon: "💡", title: "直觉：为什么延迟这么难",
      html: String.raw`
<p>物料在管道里流、数据在网络里传、算法在芯片里算——延迟无处不在。淋浴调水温人人体验过它的恶意：调完要等几秒才见效，于是矫枉过正、来回振荡。</p>
<p>数学上延迟环节 $e^{-\tau s}$ 不是有理函数：它是<strong>无穷维</strong>系统（初始条件是一整段函数），特征方程 $1 + C(s)G_0(s)e^{-\tau s} = 0$ 有无穷多个根，Routh、根轨迹全部失灵。频域是它的克星：$e^{-j\omega\tau}$ <strong>幅值恒为 1、只贡献相位 $-\omega\tau$</strong>——延迟是纯粹的"相位小偷"，Bode/Nyquist 照常营业。而 Smith 预估器（1957）的想法堪称漂亮：既然模型知道延迟多长，就用模型把"还在路上的效果"预报出来提前反馈，<strong>把延迟从回路里搬出去</strong>。</p>`
    },
    {
      icon: "📐", title: "延迟对稳定性的侵蚀",
      html: String.raw`
<div class="box def"><span class="box-title">延迟环节的频率特性</span>
$$|e^{-j\omega\tau}| = 1,\qquad \angle e^{-j\omega\tau} = -\omega\tau \ \text{(rad)}$$
<p>幅频曲线纹丝不动，相频曲线随 $\omega$ 线性下坠——高频处相位损失无上限。</p></div>
<div class="box thm"><span class="box-title">命题（最大容忍延迟）</span>
设无延迟回路 $L_0$ 稳定，相位裕度 $\mathrm{PM}$（rad）、增益穿越频率 $\omega_{gc}$。则加入延迟后闭环保持稳定当且仅当
$$\tau < \tau_{\max} = \frac{\mathrm{PM}}{\omega_{gc}}\quad(\text{单穿越情形})$$</div>
<div class="box proof"><span class="box-title">证明</span>
<p>延迟不改幅值，故 $\omega_{gc}$ 不变；它只在每个频率减去相位 $\omega\tau$。失稳临界即 Nyquist 曲线触碰 $-1$：$|L| = 1$ 处相位恰为 $-180°$。原本该处相位余量为 $\mathrm{PM}$，被吃掉 $\omega_{gc}\tau$ 后余 $\mathrm{PM} - \omega_{gc}\tau$，归零即临界，得 $\tau_{\max} = \mathrm{PM}/\omega_{gc}$。<span class="qed">∎</span></p>
<p>这解释了一条铁律：<strong>回路带宽与可容忍延迟成反比</strong>。想快（$\omega_{gc}$ 大）就得延迟小；延迟固定，带宽就有天花板 $\omega_{gc} \lesssim \mathrm{PM}/\tau$。实验里把 τ 拉大，直接看 PM 被吃穿的全过程。</p></div>`
    },
    {
      icon: "🧾", title: "Smith 预估器：定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（Smith 预估器）</span>
<p>设对象 $G(s) = G_0(s)e^{-\tau s}$（$G_0$ 稳定、有理），模型精确。控制器采用内部结构</p>
$$u = C\big(r - \underbrace{[\,y + G_0u - G_0e^{-\tau s}u\,]}_{\text{修正反馈 } \tilde y}\big)$$
<p>（$C$ 按<strong>无延迟对象</strong> $G_0$ 设计）。则闭环传函为</p>
$$\frac{Y(s)}{R(s)} = \frac{C(s)G_0(s)}{1 + C(s)G_0(s)}\,e^{-\tau s}$$
<p>——特征方程 $1 + CG_0 = 0$ <strong>不含延迟</strong>：稳定性与无延迟设计完全一致，延迟只是把响应整体平移 $\tau$。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>修正反馈的妙处：$\tilde{y} = y + G_0u - G_0e^{-\tau s}u$。代入 $y = G_0e^{-\tau s}u$（模型精确）：</p>
$$\tilde{y} = G_0e^{-\tau s}u + G_0u - G_0e^{-\tau s}u = G_0u$$
<p>——控制器"看到"的反馈恰是<strong>无延迟对象的输出</strong>！于是内环关系为 $u = C(r - G_0u)$，解出</p>
$$u = \frac{C}{1 + CG_0}\,r\ \Longrightarrow\ y = G_0e^{-\tau s}u = \frac{CG_0}{1 + CG_0}e^{-\tau s}\,r \qquad ∎$$</div>
<div class="box note"><span class="box-title">细则与局限（每条都有工程血泪）</span>
<ul>
<li><strong>模型失配</strong>：$\tau$ 估错 10~20% 性能就明显劣化——预估器把鲁棒性换成了性能（可与《不确定性建模》的语言精确分析）；</li>
<li><strong>不稳定对象禁用</strong>：结构里 $G_0$ 与 $G_0e^{-\tau s}$ 的不稳定模态在外环不可见（内部不稳定），需改用改进型（如带滤波的预估器、有限谱配置）；</li>
<li>扰动抑制的延迟躲不掉：扰动到输出的通道仍含 $e^{-\tau s}$；</li>
<li>需要有理近似时用 <strong>Padé</strong>：$e^{-\tau s} \approx \frac{1 - \tau s/2}{1 + \tau s/2}$（一阶，右半平面零点暴露了延迟的非最小相位本性）；</li>
<li>现代路线：时滞 Lyapunov–Krasovskii 泛函（LMI 判据，衔接《LMI》）、网络化控制。</li>
</ul></div>`
    },
  ],
  code: {
    title: "验证 τ_max = PM/ω_gc + Smith 预估器仿真",
    text: `import numpy as np

# 回路 L0 = 2/(s+1)^2
num, den = [2.0], np.polymul([1, 1], [1, 1])
w = np.logspace(-2, 2, 20000)
L = np.polyval(num, 1j*w) / np.polyval(den, 1j*w)
i = np.argmin(abs(abs(L) - 1))
wgc, pm = w[i], 180 + np.degrees(np.angle(L[i]))
tau_max = np.radians(pm) / wgc
print(f"PM={pm:.1f}deg @ wgc={wgc:.3f} -> 理论 tau_max={tau_max:.3f}s")

def stable_with_delay(tau, T=60.0, h=1e-3):
    """时域仿真判稳: y''+2y'+y = 2u(t-tau), 单位负反馈"""
    nd = max(1, int(tau/h)); buf = [0.0]*nd
    x = np.zeros(2); ptr = 0
    ymax_late = 0.0
    for k in range(int(T/h)):
        u = 1.0 - x[0]                    # r=1, u=e
        ud = buf[ptr]; buf[ptr] = u; ptr = (ptr+1) % nd
        x = x + h*np.array([x[1], -x[0] - 2*x[1] + 2*ud])
        if k*h > T-10: ymax_late = max(ymax_late, abs(x[0]-2/3))
        if abs(x[0]) > 1e3: return False
    return ymax_late < 0.3

for tau in [0.5*tau_max, 0.95*tau_max, 1.1*tau_max]:
    print(f"tau={tau:.3f}: {'稳定' if stable_with_delay(tau) else '失稳'}")`
  }
};

/* ════════════════ 不确定性建模与鲁棒稳定 ════════════════ */
window.KG_CONTENT["uncertainty"] = {
  sections: [
    {
      icon: "💡", title: "直觉：给『模型不准』本身建模",
      html: String.raw`
<p>"所有模型都是错的"——辨识有误差、参数会漂移、高频动态被截断。经典做法用 GM/PM 留余量，凭经验；鲁棒控制的第一步是把这句话变成数学：<strong>不确定性也有模型</strong>。</p>
<p>思路：真实对象不是一个系统，而是<strong>一族</strong>系统 $\mathcal{G} = \{G_0(1 + W\Delta)\}$——名义模型 $G_0$ 加上"误差包络" $W$（低频小、高频大，因为高频动态最难测准）。设计目标从"让 $G_0$ 闭环稳定"升级为"让<strong>族里每一个</strong> $G$ 闭环稳定"（鲁棒稳定）。美妙的是，这个看似无穷多的要求可以压缩成<strong>一条频域不等式</strong> $\|WT\|_\infty < 1$——补灵敏度 $T$ 在每个频率都要低于误差包络的倒数。这是《小增益定理》落地成设计准则的关键一步，也是 H∞ 权重函数 $W_3$ 的真实出处。</p>`
    },
    {
      icon: "📐", title: "不确定性的标准形",
      html: String.raw`
<div class="box def"><span class="box-title">定义（乘性不确定性）</span>
$$\mathcal{G} = \Big\{\,G(s) = G_0(s)\big(1 + W(s)\Delta(s)\big)\ :\ \Delta \text{ 稳定},\ \|\Delta\|_\infty \le 1\,\Big\}$$
<p>$W$ 是频率相关的误差上界：在频率 $\omega$ 处相对误差 $\Big|\frac{G - G_0}{G_0}\Big| \le |W(j\omega)|$。从实验获取：测多条频响曲线，取相对偏差包络拟合成低阶 $W$（典型形状：低频 5~10%，高频超过 100%）。</p>
<p>其它常用形态：加性 $G_0 + W\Delta$、反乘性（适合描述极点不确定）、参数化（实参数区间，归 《μ 分析》管）。</p></div>`
    },
    {
      icon: "🧾", title: "鲁棒稳定定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（乘性不确定下的鲁棒稳定）</span>
<p>设控制器 $K$ 使名义闭环稳定，$T_0 = \frac{G_0K}{1 + G_0K}$ 为名义补灵敏度。则 $K$ 对整族 $\mathcal{G}$ 鲁棒稳定 $\iff$</p>
$$\|WT_0\|_\infty < 1\quad\text{即}\quad |T_0(j\omega)| < \frac{1}{|W(j\omega)|},\ \forall\omega$$</div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>充分性（Nyquist 几何论证）。</strong>摄动回路 $L = G_0K(1 + W\Delta) = L_0 + L_0W\Delta$。逐频率看，$L(j\omega)$ 落在以 $L_0(j\omega)$ 为圆心、半径 $|L_0W|$ 的圆盘内（$|\Delta| \le 1$ 扫出整圆盘）。Nyquist 曲线连续变形时，只要圆盘<strong>始终不碰 $-1$ 点</strong>，绕圈数不会突变（同伦不变性），闭环极点数与名义情形相同 = 0。不碰的条件：</p>
$$|{-1} - L_0(j\omega)| > |L_0(j\omega)W(j\omega)|\ \forall\omega\ \iff\ \frac{|L_0W|}{|1 + L_0|} < 1\ \iff\ |WT_0| < 1\ \forall\omega$$
<p>（用了 $T_0 = L_0/(1 + L_0)$。）</p>
<p><strong>必要性（构造破坏者）。</strong>若某 $\omega_0$ 处 $|W T_0(j\omega_0)| \ge 1$，取常数相位旋转 $\Delta = e^{j\phi}/\|WT_0\|_\infty$ 类的容许摄动（可实现为稳定全通函数），使 $1 + L(j\omega_0) = (1 + L_0)(1 + WT_0\Delta) = 0$——摄动后的 Nyquist 曲线恰好穿过 $-1$，闭环临界/失稳。即族中存在不稳定成员。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">设计含义</span>
<ul>
<li>$|T_0| < 1/|W|$：误差大的频段（高频）必须把回路增益压下来——<strong>带宽的天花板由模型质量决定</strong>，这是"为什么不能无限追求快"的第一性解释；</li>
<li>鲁棒<strong>性能</strong>（最坏摄动下仍满足 $|W_1S| < 1$）的充要条件是 $\big\||W_1S_0| + |WT_0|\big\|_\infty < 1$（SISO），推导用同一套圆盘几何；</li>
<li>多入多出/结构化摄动时，本定理的"满块 $\Delta$"假设过于悲观 → 《μ 分析》。权重怎么选、控制器怎么自动算 → 《H∞》。</li>
</ul></div>`
    },
  ],
  code: {
    title: "鲁棒稳定判据验证：||WT||<1 vs 蒙特卡洛摄动",
    text: `import numpy as np

rng = np.random.default_rng(5)
w = np.logspace(-2, 2, 4000); s = 1j * w

# 名义对象 G0 = 4/(s(s+1)), 控制器 K = 1(单位反馈)
G0 = 4 / (s * (s + 1))
L0 = G0
T0 = L0 / (1 + L0)
# 误差包络 W = (0.05 + 0.4 s)/(0.2 s + 1): 低频5%, 高频200%
W = (0.05 + 0.4*s) / (0.2*s + 1)

margin = np.abs(W * T0).max()
print(f"||W T0||_inf = {margin:.3f} ->",
      "鲁棒稳定 ✓" if margin < 1 else "存在失稳摄动 ✗")

# 蒙特卡洛: 采样随机相位的边界摄动 Delta, 检查 1+L 是否绕 -1
worst = 1e9
for trial in range(200):
    phi = rng.uniform(0, 2*np.pi)
    L = L0 * (1 + W * np.exp(1j*phi))    # |Delta|=1 边界摄动
    worst = min(worst, np.abs(1 + L).min())
print(f"200 个边界摄动中 |1+L| 最小值 = {worst:.3f}",
      "(>0: 未穿越 -1)" if worst > 0.02 else "(≈0: 有摄动逼近失稳!)")`
  }
};

/* ════════════════ LMI ════════════════ */
window.KG_CONTENT["lmi"] = {
  sections: [
    {
      icon: "💡", title: "直觉：把控制设计变成凸优化",
      html: String.raw`
<p>回顾一路走来的判据：Lyapunov 方程、Riccati 方程/不等式、有界实引理……它们有个共同长相——<strong>关于未知矩阵 $P$ 是线性（仿射）的矩阵不等式</strong>。1990s 的大发现（Boyd 等）：这类问题是<strong>凸</strong>的，而凸问题有多项式时间的内点法求解器——于是"设计控制器"变成"写下 LMI，交给求解器"。</p>
<p>LMI 的价值不在解单个老问题（Riccati 有专用解法），而在<strong>拼积木</strong>：H∞ 性能 + 极点区域 + 多胞不确定性 + 输入饱和，每个要求各写一条 LMI，全部叠在一起仍是凸的——传统解析方法对这种"混搭"束手无策。今天读鲁棒/时滞/网络化控制论文，LMI 是默认语言。</p>`
    },
    {
      icon: "📐", title: "定义与基本事实",
      html: String.raw`
<div class="box def"><span class="box-title">定义（线性矩阵不等式）</span>
$$F(x) = F_0 + x_1F_1 + \cdots + x_mF_m \succ 0$$
<p>$F_i$ 为给定对称矩阵，$x \in \mathbb{R}^m$ 为决策变量（矩阵变量按元素摊开即此形式）。"$\succ 0$" = 正定。</p></div>
<div class="box thm"><span class="box-title">命题（可行集是凸集）</span>
$\{x : F(x) \succ 0\}$ 是凸集。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>设 $F(x) \succ 0$，$F(y) \succ 0$，$\theta \in [0,1]$。由 $F$ 仿射：$F(\theta x + (1-\theta)y) = \theta F(x) + (1-\theta)F(y)$。正定矩阵的非负组合正定（对任意 $v \ne 0$：$v^{\top}[\theta F(x) + (1-\theta)F(y)]v = \theta v^{\top}F(x)v + (1-\theta)v^{\top}F(y)v > 0$）。<span class="qed">∎</span></p>
<p>第一个例子就在《线性系统稳定性》里：$A$ Hurwitz $\iff$ 存在 $P \succ 0$ 使 $A^{\top}P + PA \prec 0$——关于 $P$ 线性，是标准 LMI 可行性问题。</p></div>`
    },
    {
      icon: "🧾", title: "两大核心工具：Schur 补与变量代换",
      html: String.raw`
<div class="box thm"><span class="box-title">引理（Schur 补——非线性变凸的万能钥匙）</span>
<p>设 $R = R^{\top} \succ 0$。则</p>
$$\begin{bmatrix} Q & S \\ S^{\top} & R \end{bmatrix} \succ 0 \iff Q - SR^{-1}S^{\top} \succ 0$$</div>
<div class="box proof"><span class="box-title">证明</span>
<p>取合同变换（congruence）矩阵 $M = \begin{bmatrix} I & -SR^{-1} \\ 0 & I \end{bmatrix}$（可逆），直接计算：</p>
$$M\begin{bmatrix} Q & S \\ S^{\top} & R \end{bmatrix}M^{\top} = \begin{bmatrix} Q - SR^{-1}S^{\top} & 0 \\ 0 & R \end{bmatrix}$$
<p>合同变换保持正定性（$X \succ 0 \iff MXM^{\top} \succ 0$，因 $v^{\top}MXM^{\top}v = (M^{\top}v)^{\top}X(M^{\top}v)$ 且 $M^{\top}$ 可逆）。块对角阵正定 $\iff$ 各块正定，$R \succ 0$ 已设，故等价于 $Q - SR^{-1}S^{\top} \succ 0$。<span class="qed">∎</span></p>
<p><strong>用途</strong>：Riccati 不等式 $A^{\top}P + PA + PBR^{-1}B^{\top}P + Q \prec 0$ 里的二次项 $PBR^{-1}B^{\top}P$ 对 $P$ 非线性——用 Schur 补翻写成</p>
$$\begin{bmatrix} -(A^{\top}P + PA + Q) & PB \\ B^{\top}P & R \end{bmatrix} \succ 0$$
<p>——对 $P$ 线性了！有界实引理（《H∞》）同法变 LMI。</p></div>
<div class="box thm"><span class="box-title">定理（状态反馈镇定的 LMI 综合）</span>
<p>存在 $K$ 使 $A + BK$ Hurwitz $\iff$ 存在 $Y = Y^{\top} \succ 0$ 与矩阵 $Z$ 满足 LMI</p>
$$AY + YA^{\top} + BZ + Z^{\top}B^{\top} \prec 0$$
<p>且此时 $K = ZY^{-1}$。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>闭环稳定 $\iff$ 存在 $P \succ 0$：$(A + BK)^{\top}P + P(A + BK) \prec 0$。难点：$K$ 与 $P$ 相乘，双线性、非凸。<strong>两步标准手术</strong>：</p>
<p>① 左右同乘 $Y \triangleq P^{-1}$（合同变换，保号）：</p>
$$Y(A + BK)^{\top} + (A + BK)Y \prec 0$$
<p>② 变量代换 $Z \triangleq KY$：$BKY = BZ$，得 $AY + YA^{\top} + BZ + Z^{\top}B^{\top} \prec 0$——对 $(Y, Z)$ 联合线性。反向：解出 $(Y, Z)$ 后取 $K = ZY^{-1}$、$P = Y^{-1}$，逐步倒推回去。<span class="qed">∎</span></p>
<p>这套"逆变量 + 换元"的手术刀在 LMI 文献里反复出现（H∞ 综合、极点区域配置、增益调度 LPV 全是它的变奏）。</p></div>
<div class="box note"><span class="box-title">生态</span>
求解器：SeDuMi / Mosek / SCS；建模层：MATLAB 的 YALMIP、Python 的 cvxpy。极点区域约束（LMI region：$\alpha$-稳定度、扇形阻尼比）、多胞系统同一 $P$ 的鲁棒 LMI、时滞 Lyapunov–Krasovskii 判据（衔接《时滞系统》）都是标准积木。局限：LMI 给的常是充分条件（保守），且大规模问题维数增长快。</div>`
    },
  ],
  code: {
    title: "cvxpy 解 LMI：一行式镇定综合 K = Z Y^{-1}",
    text: `import numpy as np
import cvxpy as cp

# 不稳定对象
A = np.array([[1.0, 2.0], [0.0, 0.5]])
B = np.array([[0.], [1.]])
n = 2

Y = cp.Variable((n, n), symmetric=True)
Z = cp.Variable((1, n))
eps = 1e-3
constraints = [
    Y >> eps * np.eye(n),
    A @ Y + Y @ A.T + B @ Z + Z.T @ B.T << -eps * np.eye(n),
]
prob = cp.Problem(cp.Minimize(0), constraints)
prob.solve(solver=cp.SCS)
print("LMI 可行:", prob.status)

K = Z.value @ np.linalg.inv(Y.value)
print("K =", np.round(K, 3))
print("闭环特征值:", np.round(np.linalg.eigvals(A + B @ K), 3))
# 换成 Riccati 也能镇定, 但 LMI 的好处是可以继续叠约束:
# 例如加极点区域 Re(lambda) < -1: A Y + Y A' + BZ + Z'B' + 2*1*Y << 0`
  }
};

/* ════════════════ μ 分析 ════════════════ */
window.KG_CONTENT["mu"] = {
  sections: [
    {
      icon: "💡", title: "直觉：小增益太悲观，μ 量体裁衣",
      html: String.raw`
<p>《小增益定理》和《不确定性建模》把所有不确定性打包成一个"满块" $\Delta$——它允许不确定性在任意通道间随意串扰。可现实中的不确定性是<strong>有结构的</strong>：质量偏差只进质量通道、执行器误差只在自己的对角位置、实参数根本不是复数。按满块算，等于假设敌人无所不能——结论安全但<strong>过度保守</strong>：明明能容忍 50% 的参数摄动，满块判据可能只敢承诺 10%。</p>
<p>Doyle（1982）的 μ（结构化奇异值）把"结构"写进定义：只在<strong>指定结构的摄动集合</strong>里找最小的破坏者。μ 越小越鲁棒，μ < 1 恰好是鲁棒稳定的精确刻画——不多一分保守，不少一分安全。代价：μ 精确计算是 NP-hard，实用靠上下界夹逼（D-scaling）。</p>`
    },
    {
      icon: "📐", title: "定义与主定理",
      html: String.raw`
<div class="box def"><span class="box-title">定义（结构化奇异值 μ）</span>
<p>给定复矩阵 $M$ 与块对角结构 $\boldsymbol{\Delta} = \{\mathrm{diag}(\delta_1I_{r_1}, \dots, \Delta_F)\}$（实/复标量块与满块的混合）：</p>
$$\mu_{\boldsymbol{\Delta}}(M) \triangleq \frac{1}{\min\{\bar\sigma(\Delta) : \Delta \in \boldsymbol{\Delta},\ \det(I - M\Delta) = 0\}}$$
<p>（若无任何 $\Delta$ 使 $\det(I - M\Delta) = 0$，定义 $\mu = 0$。）直白说：<strong>1/μ = 能破坏回路良定性的最小结构化摄动的大小</strong>。</p></div>
<div class="box thm"><span class="box-title">定理（鲁棒稳定的精确判据 / 主回路定理）</span>
把系统整理成 $M$–$\Delta$ 标准互联（$M$ 为名义闭环从摄动输入到摄动输出的传函）。则对一切 $\|\Delta\|_\infty < 1$、结构 $\in \boldsymbol{\Delta}$ 的摄动闭环稳定 $\iff$
$$\sup_{\omega}\ \mu_{\boldsymbol{\Delta}}\big(M(j\omega)\big) \le 1$$
（μ 版 Nyquist：逐频率检查"最小破坏摄动"是否大于 1。证明骨架与《不确定性建模》的同伦论证一致，把"圆盘不碰 -1"换成"det(I−MΔ) 不过零"。）</div>`
    },
    {
      icon: "🧾", title: "上下界与 D-scaling：定理与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（基本夹逼）</span>
$$\rho(M) \le \mu_{\boldsymbol{\Delta}}(M) \le \bar\sigma(M)$$
（$\rho$ = 谱半径，对应结构最"瘦"的情形 $\Delta = \delta I$；$\bar\sigma$ 对应满块。）</div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>上界。</strong>若 $\bar\sigma(\Delta) < 1/\bar\sigma(M)$，则 $\bar\sigma(M\Delta) \le \bar\sigma(M)\bar\sigma(\Delta) < 1$，故 $I - M\Delta$ 的最小奇异值 $> 0$，可逆——没有这么小的破坏者。所以最小破坏摄动 $\ge 1/\bar\sigma(M)$，即 $\mu \le \bar\sigma(M)$。</p>
<p><strong>下界。</strong>结构 $\boldsymbol{\Delta}$ 总包含 $\delta I$（所有块取同一标量）。取 $M$ 的模最大特征值 $\lambda$（$|\lambda| = \rho$），令 $\delta = 1/\lambda$：$\det(I - M\delta) = \det(I - M/\lambda) = 0$（$\lambda$ 是特征值）。该摄动大小 $|\delta| = 1/\rho$，故 $\mu \ge \rho(M)$。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理（D-scaling 不变性与可算上界）</span>
<p>设 $\mathcal{D}$ = 与结构 $\boldsymbol{\Delta}$ 可交换的正定对角缩放集（$D\Delta = \Delta D$）。则</p>
$$\mu_{\boldsymbol{\Delta}}(M) = \mu_{\boldsymbol{\Delta}}(DMD^{-1}) \le \inf_{D \in \mathcal{D}}\ \bar\sigma(DMD^{-1})$$
<p>右端对 $D$（取对数坐标）是凸问题——LMI 可解（衔接《LMI》）。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>不变性：$\det(I - DMD^{-1}\Delta) = \det(I - DMD^{-1}\Delta DD^{-1}) = \det\big(D(I - M\Delta)D^{-1}\big) = \det(I - M\Delta)$，其中用了交换性 $D^{-1}\Delta D = \Delta$。破坏者集合完全相同 ⇒ μ 相同。再对缩放后的矩阵套用上界 $\mu \le \bar\sigma$，对所有 $D$ 取最优即得。<span class="qed">∎</span></p>
<p>复块 ≤ 3 个时上界恰好取等；一般情形间隙通常百分之几——工程足够。<strong>μ 综合（D–K 迭代）</strong>：固定 $D$ 解一次 H∞（得 $K$），固定 $K$ 逐频率优化 $D$，交替至收敛。非凸、无全局保证，但为数不多的"结构化鲁棒性能"系统化设计法，MATLAB <em>musyn</em> 即此。</p></div>
<div class="box note"><span class="box-title">定位</span>
μ 是鲁棒控制金字塔的塔尖：《小增益》给哲学，《不确定性建模》给语言，《H∞》给综合引擎，《LMI》给数值工具，μ 把它们拧成"结构化不确定性下不保守的分析与设计"。实参数占主导时另有 Kharitonov 定理等区间方法可作对照。</div>`
    },
  ],
  code: {
    title: "μ 的上下界数值验证（含 D-scaling 改进）",
    text: `import numpy as np
from scipy.optimize import minimize_scalar

# 2x2 例子, 结构 = diag(delta1, delta2) (两个复标量块)
M = np.array([[0.4, 3.0],
              [0.02, 0.5]])

rho = np.abs(np.linalg.eigvals(M)).max()
sig = np.linalg.svd(M, compute_uv=False).max()

# 对角结构的 mu: 此时 mu = max over 相位的 rho(M diag(e^{j t1}, e^{j t2}))
# 数值扫相位求"真值"
best = 0.0
for t1 in np.linspace(0, 2*np.pi, 60):
    for t2 in np.linspace(0, 2*np.pi, 60):
        D = np.diag([np.exp(1j*t1), np.exp(1j*t2)])
        best = max(best, np.abs(np.linalg.eigvals(M @ D)).max())

# D-scaling 上界: min over d>0 of sigma_max(diag(d,1) M diag(1/d,1))
def ub(logd):
    d = np.exp(logd)
    Ms = np.diag([d, 1.0]) @ M @ np.diag([1/d, 1.0])
    return np.linalg.svd(Ms, compute_uv=False).max()
res = minimize_scalar(ub, bounds=(-6, 6), method='bounded')

print(f"下界 rho(M)        = {rho:.4f}")
print(f"mu (相位扫描)      = {best:.4f}")
print(f"D-scaling 上界     = {res.fun:.4f}   (两复块: 应与 mu 相等)")
print(f"未缩放上界 sigma   = {sig:.4f}   (满块假设的保守度一目了然)")`
  }
};

/* ════════════════ ISS ════════════════ */
window.KG_CONTENT["iss"] = {
  sections: [
    {
      icon: "💡", title: "直觉：有扰动世界的稳定性语言",
      html: String.raw`
<p>Lyapunov 理论回答"没有扰动时会不会回到平衡点"；可现实系统永远泡在扰动里——风、噪声、邻居子系统的输出。我们真正想要的是：<strong>扰动有界则状态有界，扰动越小偏差越小，扰动消失则回到原点</strong>。三个愿望打包成一个定义，就是 Sontag（1989）的输入-状态稳定 ISS。</p>
<p>ISS 之所以成为现代非线性控制的"普通话"：它<strong>可组合</strong>。两个各自 ISS 的子系统级联还是 ISS，回路互联满足非线性小增益条件则闭环稳定——大系统的稳定性可以像电路分析一样逐块推。反步法的鲁棒化、多智能体、网络化控制的理论根基都写在这门语言里。</p>`
    },
    {
      icon: "📐", title: "定义",
      html: String.raw`
<div class="box def"><span class="box-title">比较函数（非线性版的"增益刻度"）</span>
<ul>
<li>$\alpha \in \mathcal{K}$：连续、严格增、$\alpha(0) = 0$；$\mathcal{K}_\infty$：还要 $\alpha(r) \to \infty$；</li>
<li>$\beta \in \mathcal{KL}$：对第一个变量属 $\mathcal{K}$，对第二个变量连续递减趋于 0。</li>
</ul></div>
<div class="box def"><span class="box-title">定义（ISS）</span>
系统 $\dot{x} = f(x, u)$ 称为<strong>输入-状态稳定</strong>，若存在 $\beta \in \mathcal{KL}$、$\gamma \in \mathcal{K}$ 使得对一切 $x_0$ 与有界输入：
$$\|x(t)\| \le \beta(\|x_0\|, t) + \gamma\big(\sup_{s \le t}\|u(s)\|\big),\quad \forall t \ge 0$$
<p>读法：<strong>初始条件的影响衰减掉（$\beta$ 项），留下的偏差被扰动大小控制（$\gamma$ 项，"非线性增益"）</strong>。$u \equiv 0$ 时退化为全局渐近稳定——ISS 是 GAS 的带输入推广。</p></div>
<div class="box eg"><span class="box-title">正反两个例子</span>
<p>✓ $\dot{x} = -x + u$：解的估计 $\|x(t)\| \le e^{-t}\|x_0\| + \sup|u|$，ISS，线性增益 $\gamma(r) = r$。</p>
<p>✗ $\dot{x} = -x + xu$：取 $u \equiv 2$ 得 $\dot x = x$——<strong>有界的小扰动引爆发散</strong>。注意该系统 $u=0$ 时全局渐近稳定！GAS 并不承诺任何抗扰性，这正是需要 ISS 概念的原因。</p></div>`
    },
    {
      icon: "🧾", title: "ISS-Lyapunov 定理：判据与证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（ISS-Lyapunov 函数）</span>
<p>设存在 $V \in C^1$、$\alpha_1, \alpha_2 \in \mathcal{K}_\infty$、$\alpha_3, \sigma \in \mathcal{K}$：</p>
$$\alpha_1(\|x\|) \le V(x) \le \alpha_2(\|x\|),\qquad \dot{V} \le -\alpha_3(\|x\|) + \sigma(\|u\|)$$
<p>则系统 ISS，且增益可取 $\gamma = \alpha_1^{-1} \circ \alpha_2 \circ \alpha_3^{-1} \circ 2\sigma$。</p></div>
<div class="box proof"><span class="box-title">证明（核心机制完整）</span>
<p>记 $\bar{u} = \sup_s\|u(s)\|$，定义"污染半径" $r = \alpha_3^{-1}(2\sigma(\bar{u}))$。把状态空间分成两个区域：</p>
<p><strong>区域一：$\|x\| \ge r$（扰动压不住耗散）。</strong>此处 $\alpha_3(\|x\|) \ge 2\sigma(\bar u)$，于是</p>
$$\dot{V} \le -\alpha_3(\|x\|) + \sigma(\bar{u}) \le -\tfrac{1}{2}\alpha_3(\|x\|) \le -\tfrac{1}{2}\alpha_3\big(\alpha_2^{-1}(V)\big) < 0$$
<p>$V$ 以一致速率下降。由比较引理，存在 $\beta \in \mathcal{KL}$ 使区域一内 $V(t) \le \beta_V(V(0), t)$（严格负定微分不等式的标准 KL 估计）。</p>
<p><strong>区域二：水平集陷阱。</strong>定义 $c = \alpha_2(r)$，考察水平集 $\Omega_c = \{V \le c\}$。在其边界上 $V = c \ge \alpha_2(\|x\|)$... 取边界点：$V(x) = c$ 蕴含 $\|x\| \ge \alpha_2^{-1}(c) = r$，故边界上 $\dot{V} < 0$——<strong>$\Omega_c$ 正不变</strong>：轨迹一旦进入永不出去。而区域一的下降保证有限时间进入 $\Omega_c$。</p>
<p><strong>合并。</strong>进入后 $V \le c = \alpha_2(r)$，换算回状态：</p>
$$\|x\| \le \alpha_1^{-1}(V) \le \alpha_1^{-1}\big(\alpha_2(\alpha_3^{-1}(2\sigma(\bar{u})))\big) = \gamma(\bar{u})$$
<p>进入前由区域一的 KL 估计控制。两段拼接（取逐点最大）恰为 ISS 不等式的形状：$\|x(t)\| \le \beta(\|x_0\|, t) + \gamma(\bar{u})$。<span class="qed">∎</span></p>
<p><em>注：反方向也成立（ISS ⇒ 存在 ISS-Lyapunov 函数，Sontag–Wang），故这是完全刻画。</em></p></div>
<div class="box thm"><span class="box-title">推论（级联定理）</span>
若 $\dot{x}_1 = f_1(x_1, x_2)$ 关于输入 $x_2$ 是 ISS，且 $\dot{x}_2 = f_2(x_2)$ GAS，则级联系统整体 GAS。</div>
<div class="box proof"><span class="box-title">证明思路</span>
<p>$x_2(t) \to 0$（GAS）。对 $x_1$ 的 ISS 估计分段用："过去段"的影响装进衰减的 $\beta$ 项，"近期段"的输入界 $\sup_{s \ge t/2}\|x_2(s)\| \to 0$ 使 $\gamma$ 项趋零。两项都归零 ⇒ $x_1 \to 0$；有界性由 ISS 直接给出。<span class="qed">∎</span></p>
<p><strong>为什么重要</strong>：《反步法》每一步"虚拟控制误差 → 下一级输入"的合法性、观测器误差进入控制回路的无害性（非线性分离原理的替身）、非线性小增益定理（两个 ISS 系统回路互联，增益复合 $\gamma_1 \circ \gamma_2 < \mathrm{id}$ 则稳定，衔接《小增益定理》）——全靠 ISS 的可组合性。</p></div>`
    },
  ],
  code: {
    title: "数值验证 ISS 界与『GAS ≠ 抗扰』反例",
    text: `import numpy as np

def simulate(f, x0, u_func, T=20.0, h=1e-3):
    x, xs = x0, []
    for k in range(int(T/h)):
        xs.append(x)
        x = x + h * f(x, u_func(k*h))
        if abs(x) > 1e6: break
    return np.array(xs)

# 系统A: xdot = -x + u  (ISS, 增益=1)
fA = lambda x, u: -x + u
# 系统B: xdot = -x + x*u (GAS 但非 ISS)
fB = lambda x, u: -x + x*u

u_bounded = lambda t: 2.0 * (t > 5)      # t>5 后加常值扰动 u=2

xa = simulate(fA, 3.0, u_bounded)
xb = simulate(fB, 3.0, u_bounded)
print("系统A (ISS)   : 末态 =", round(float(xa[-1]), 3),
      " <= gamma(2)=2  ✓ 偏差被扰动大小限住")
print("系统B (仅GAS) : 末态 =", f"{float(xb[-1]):.2e}",
      " 有界扰动 -> 发散!")

# 验证A的 ISS 不等式 |x(t)| <= e^{-t}|x0| + sup|u|
h = 1e-3
t = np.arange(len(xa)) * h
bound = np.exp(-t)*3.0 + 2.0
print("A 全程满足 ISS 界:", bool(np.all(np.abs(xa) <= bound + 1e-6)))`
  }
};
