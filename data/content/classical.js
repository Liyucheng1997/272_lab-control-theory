/* 内容：经典控制（频域） */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ PID ════════════════ */
window.KG_CONTENT["pid"] = {
  sim: "sim-pid",
  sections: [
    {
      icon: "💡", title: "直觉：过去、现在、未来",
      html: String.raw`
<p>PID 是世界上部署最广的控制律——工业回路的 90% 以上。它对误差 $e = r - y$ 做三件事：</p>
<ul>
<li><strong>P（现在）</strong>：误差多大就使多大劲，$K_p e$。快，但"够不着"零误差——误差为零时它也出力为零，抵不住恒定负载；</li>
<li><strong>I（过去）</strong>：把历史误差累积起来，$K_i \int e$。只要误差不归零它就一直加码，<strong>保证消灭稳态误差</strong>；代价是"记仇"，容易过冲；</li>
<li><strong>D（未来）</strong>：看误差的变化趋势提前刹车，$K_d \dot{e}$。给系统加阻尼、压超调；代价是放大噪声。</li>
</ul>
<div class="box def"><span class="box-title">控制律</span>
$$u(t) = K_p e(t) + K_i \int_0^t e(\tau)\,d\tau + K_d \frac{de}{dt},\qquad C(s) = K_p + \frac{K_i}{s} + K_d s$$</div>`
    },
    {
      icon: "🧾", title: "定理：积分作用消灭稳态误差",
      html: String.raw`
<div class="box thm"><span class="box-title">定理</span>
设被控对象 $G(s)$ 与 PID 控制器构成的单位负反馈闭环稳定，且 $K_i \ne 0$。则对阶跃参考与阶跃负载扰动，<strong>稳态误差为零</strong>。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>误差传函（参考到误差）：$\dfrac{E(s)}{R(s)} = \dfrac{1}{1 + C(s)G(s)}$。阶跃 $R(s) = 1/s$，用终值定理（闭环稳定保证适用条件：$sE(s)$ 的极点均在左半平面）：</p>
$$e_{\infty} = \lim_{s\to 0} sE(s) = \lim_{s\to 0} \frac{1}{1 + C(s)G(s)}$$
<p>因 $C(s) = \frac{K_d s^2 + K_p s + K_i}{s}$ 含因子 $\frac{1}{s}$，故 $\lim_{s\to0}|C(s)| = \infty$（$K_i \neq 0$，且 $G(0) \ne 0$），于是 $e_\infty = 0$。</p>
<p>对输入端阶跃扰动 $d$：$\dfrac{Y(s)}{D(s)} = \dfrac{G}{1 + CG}$，同理 $y_\infty = \lim_{s\to0} \dfrac{G(s)}{1 + C(s)G(s)} = \lim_{s\to0}\dfrac{1}{1/G(s) + C(s)} = 0$。<span class="qed">∎</span></p>
<p>这是<strong>内模原理</strong>的最简单例：想无静差跟踪/抑制某类信号，回路里必须内置该信号的模型——阶跃信号的模型就是积分器 $1/s$。想无差跟踪斜坡？回路里要有 $1/s^2$。</p></div>
<div class="box note"><span class="box-title">工程细节（实验里都能看到）</span>
<ul>
<li><strong>微分不作用于设定值</strong>：设定值阶跃瞬间 $\dot e$ 是冲激，会打坏执行器；实用形式对测量值微分并加低通 $\frac{K_d s}{s/N + 1}$；</li>
<li><strong>积分抗饱和 (anti-windup)</strong>：执行器饱和时积分继续累积会造成大幅超调，需钳位或反算；</li>
<li><strong>为什么 PID 常常够用</strong>：多数过程在工作频段近似一二阶 + 延迟，三个参数恰好塑形增益穿越频率附近的幅相特性（与《Bode》《超前滞后》视角一致）。</li>
</ul></div>`
    },
    {
      icon: "📐", title: "Ziegler–Nichols 整定",
      html: String.raw`
<p>1942 年提出的经验整定法，至今仍是现场调参的起点（激进版，超调偏大，适合作初值再手调）：</p>
<div class="box eg"><span class="box-title">临界比例法</span>
<p>① 关掉 I、D，只用 P；② 逐渐加大 $K_p$ 直至输出等幅振荡，记下临界增益 $K_u$ 和振荡周期 $T_u$；③ 查表：</p>
$$\begin{array}{c|ccc}
 & K_p & K_i & K_d\\\hline
\text{P} & 0.5K_u & - & -\\
\text{PI} & 0.45K_u & 0.54K_u/T_u & -\\
\text{PID} & 0.6K_u & 1.2K_u/T_u & 0.075K_uT_u
\end{array}$$
<p>临界点 $(K_u, T_u)$ 实质上是 Nyquist 曲线穿过 $-1$ 点的信息——ZN 法是在用最少的实验数据做频域整形。</p></div>`
    },
  ],
  code: {
    title: "带微分滤波与抗饱和的离散 PID（可直接嵌入实际项目）",
    text: `import numpy as np

class PID:
    def __init__(self, kp, ki, kd, dt, N=20.0, umin=-np.inf, umax=np.inf):
        self.kp, self.ki, self.kd = kp, ki, kd
        self.dt, self.N = dt, N
        self.umin, self.umax = umin, umax
        self.integ = 0.0
        self.d_state = 0.0      # 微分滤波器状态
        self.y_prev = None

    def step(self, r, y):
        e = r - y
        # 微分作用于测量值 + 一阶低通 (避免设定值踢与噪声放大)
        if self.y_prev is None:
            self.y_prev = y
        dy = (y - self.y_prev) / self.dt
        self.d_state += self.N * self.dt * (dy - self.d_state)
        self.y_prev = y

        u_unsat = self.kp*e + self.ki*self.integ - self.kd*self.d_state
        u = np.clip(u_unsat, self.umin, self.umax)
        # 条件积分抗饱和: 只在不加剧饱和时累积
        if (u == u_unsat) or (np.sign(e) != np.sign(u_unsat)):
            self.integ += e * self.dt
        return u

# 演示: 控制 G(s)=1/(s^2+s+1), 限幅 |u|<=2
pid = PID(kp=6, ki=4, kd=3, dt=0.001, umax=2, umin=-2)
x = np.zeros(2)
log = []
for k in range(15000):
    t = k * 0.001
    y = x[0]
    u = pid.step(1.0, y)
    d = -0.5 if t > 10 else 0.0          # 10s 后加扰动
    x = x + 0.001 * np.array([x[1], -x[0] - x[1] + u + d])
    log.append((t, y, u))
log = np.array(log)
print("稳态输出:", log[-1, 1].round(4), "(应为 1: 积分作用消差)")`
  }
};

/* ════════════════ Bode ════════════════ */
window.KG_CONTENT["bode"] = {
  sim: "sim-bode",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>由《传递函数》的正弦稳态定理，系统对每个频率的正弦有确定的放大倍数与相移。Bode 图（1930s，贝尔实验室）把它们画成两条曲线：幅频 $20\lg|G(j\omega)|$（分贝）与相频 $\angle G(j\omega)$，横轴取对数。</p>
<p>取对数的妙处：串联环节 $G_1G_2$ 的图形是<strong>两图相加</strong>；每个一阶极点贡献 $-20$ dB/dec 斜率和 $-90°$ 相位，零点相反——复杂系统的频率特性可以像搭积木一样手绘出来。控制设计从"配极点"变成"<strong>整形回路增益曲线 $L(j\omega)$</strong>"：低频高增益（跟踪好）、穿越频率适中（带宽）、高频快速衰减（抗噪）。</p>`
    },
    {
      icon: "📐", title: "稳定裕度",
      html: String.raw`
<div class="box def"><span class="box-title">定义（增益裕度 GM / 相位裕度 PM）</span>
<p>设开环 $L(j\omega)$。记<strong>相位穿越频率</strong> $\omega_{pc}$：$\angle L(j\omega_{pc}) = -180°$；<strong>增益穿越频率</strong> $\omega_{gc}$：$|L(j\omega_{gc})| = 1$。</p>
$$\mathrm{GM} = \frac{1}{|L(j\omega_{pc})|},\qquad \mathrm{PM} = 180° + \angle L(j\omega_{gc})$$
<p>GM：增益还能放大多少倍才失稳；PM：相位还能再滞后多少度才失稳（即还能容忍多大延迟：$\tau_{max} = \mathrm{PM}/\omega_{gc}$）。</p></div>
<div class="box proof"><span class="box-title">为何这两个量决定"离失稳多远"（由 Nyquist 判据）</span>
<p>闭环临界失稳 $\iff$ Nyquist 曲线恰好穿过 $-1$ 点，即存在 $\omega$：$|L| = 1$ 且 $\angle L = -180°$ 同时成立。GM 度量在 $\angle L = -180°$ 处距 $|L|=1$ 的余量，PM 度量在 $|L|=1$ 处距 $-180°$ 的余量——正是到 $-1$ 点的两个方向的距离。更精细的单一指标是 $\displaystyle\min_\omega |1 + L(j\omega)|$（灵敏度峰值的倒数，通往《H∞》）。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">经验设计法则</span>
<ul>
<li>PM ≈ 30°~60°，GM ≥ 6 dB 是常规目标；对二阶系统 PM ≈ 100ζ（%）；</li>
<li>闭环带宽 ≈ $\omega_{gc}$，上升时间 $t_r \approx 1.8/\omega_{gc}$（实验中可验证）；</li>
<li><strong>Bode 增益-相位关系</strong>：最小相位系统的相位由幅频斜率唯一决定，$\angle L \approx$ 斜率 $\times 90°$（-20dB/dec ↔ -90°）。所以设计法则是"让幅频曲线以 -20 dB/dec 的斜率穿越 0 dB"——斜率太陡（-40）意味着相位接近 -180°，PM 不足。</li>
</ul></div>`
    },
  ],
  code: {
    title: "手搓 Bode 图与裕度计算",
    text: `import numpy as np
import matplotlib.pyplot as plt

def bode_margins(num, den, w=None):
    w = np.logspace(-2, 3, 4000) if w is None else w
    s = 1j * w
    L = np.polyval(num, s) / np.polyval(den, s)
    mag_db = 20 * np.log10(abs(L))
    ph = np.unwrap(np.angle(L)) * 180 / np.pi
    # 增益穿越: |L|=1
    idx_g = np.where(np.diff(np.sign(mag_db)))[0]
    wgc, pm = (None, None)
    if len(idx_g):
        i = idx_g[0]
        wgc = np.interp(0, [mag_db[i+1], mag_db[i]], [w[i+1], w[i]])
        pm = 180 + np.interp(np.log10(wgc), np.log10(w), ph)
    # 相位穿越: phase=-180
    idx_p = np.where(np.diff(np.sign(ph + 180)))[0]
    wpc, gm_db = (None, None)
    if len(idx_p):
        i = idx_p[0]
        wpc = np.interp(-180, [ph[i+1], ph[i]], [w[i+1], w[i]])
        gm_db = -np.interp(np.log10(wpc), np.log10(w), mag_db)
    return w, mag_db, ph, (wgc, pm), (wpc, gm_db)

# L(s) = 4/(s(s+1)(0.2s+1))
num = [4]
den = np.polymul([1, 0], np.polymul([1, 1], [0.2, 1]))
w, m, p, (wgc, pm), (wpc, gm) = bode_margins(num, den)
print(f"PM = {pm:.1f} deg @ {wgc:.2f} rad/s")
print(f"GM = {gm:.1f} dB  @ {wpc:.2f} rad/s")

fig, ax = plt.subplots(2, 1, sharex=True)
ax[0].semilogx(w, m); ax[0].axhline(0, c='gray', ls='--')
ax[1].semilogx(w, p); ax[1].axhline(-180, c='r', ls='--')
plt.show()`
  }
};

/* ════════════════ Nyquist ════════════════ */
window.KG_CONTENT["nyquist"] = {
  sim: "sim-nyquist",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>1932 年，贝尔实验室的 Nyquist 要回答：反馈放大器什么时候会啸叫（失稳）？他给出了一个惊人的答案：<strong>不需要知道闭环极点在哪，只要看开环频率响应曲线绕不绕 $-1$ 这个点</strong>。</p>
<p>为什么是 $-1$？闭环极点是 $1 + L(s) = 0$ 的根，即 $L(s) = -1$。开环曲线离 $-1$ 有多远，闭环就离失稳有多远。它是唯一能直接处理<strong>时滞</strong>（$e^{-s\tau}$ 不是有理函数，Routh 用不了）并且能从<strong>实验数据</strong>（扫频测量）直接判稳的方法，也是鲁棒控制的几何起点。</p>`
    },
    {
      icon: "🧾", title: "辐角原理与 Nyquist 判据的证明",
      html: String.raw`
<div class="box thm"><span class="box-title">引理（辐角原理，复分析）</span>
设 $F(s)$ 在简单闭合曲线 $\Gamma$ 内亚纯，在 $\Gamma$ 上无零点极点，内部有 $Z$ 个零点、$P$ 个极点（按重数）。当 $s$ 顺时针绕行 $\Gamma$ 一周，$F(s)$ 的像曲线绕原点<strong>顺时针</strong> $N = Z - P$ 圈。</div>
<div class="box proof"><span class="box-title">引理证明</span>
<p>像曲线绕原点的圈数由辐角总变化决定：</p>
$$N = -\frac{1}{2\pi}\Delta_\Gamma \arg F(s) = -\frac{1}{2\pi}\,\mathrm{Im}\oint_\Gamma \frac{F'(s)}{F(s)}ds \cdot \frac{1}{1} = -\frac{1}{2\pi j}\oint_\Gamma \frac{F'(s)}{F(s)}\,ds \quad (\text{顺时针取负})$$
<p>在零点 $s_0$（重数 $m$）附近 $F(s) = (s-s_0)^m g(s)$，$g(s_0)\neq0$：</p>
$$\frac{F'}{F} = \frac{m}{s - s_0} + \frac{g'}{g}$$
<p>即 $F'/F$ 在零点处有留数 $m$ 的单极点；同理在 $F$ 的 $m$ 重极点处留数为 $-m$。由留数定理（逆时针为正，顺时针加负号）：</p>
$$\frac{1}{2\pi j}\oint_{\Gamma,\circlearrowright} \frac{F'}{F}ds = -(Z - P) \ \Rightarrow\ N = Z - P \qquad ∎$$</div>
<div class="box thm"><span class="box-title">定理（Nyquist 稳定判据）</span>
<p>取 $F(s) = 1 + L(s)$，$\Gamma$ 为 Nyquist 围线（整条虚轴 + 右侧无穷大半圆，顺时针）。设开环 $L(s)$ 在右半平面有 $P$ 个极点。则：</p>
<p style="text-align:center"><strong>闭环右半平面极点数 $Z = N + P$</strong>，其中 $N$ 为 $L(j\omega)$ 曲线绕 $-1$ 点的顺时针圈数。</p>
<p>闭环稳定 $\iff N = -P$（须逆时针绕 $P$ 圈）；开环稳定时即 <strong>不包围 $-1$</strong>。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p>关键对应：① $1 + L(s)$ 的<strong>零点 = 闭环极点</strong>（闭环传函 $\frac{L}{1+L}$ 的分母）；② $1 + L(s)$ 的<strong>极点 = 开环 $L$ 的极点</strong>（加 1 不改变极点）。</p>
<p>对 $F = 1+L$ 在 Nyquist 围线上用辐角原理：$F$ 的像绕<strong>原点</strong>的圈数 $N = Z - P$，其中 $Z$、$P$ 分别是右半平面内的闭环极点数与开环极点数。而 $F = 1 + L$ 的像绕原点 $\iff$ $L$ 的像绕 $-1$（整体平移 1）。无穷大半圆上：物理系统 $L$ 真分式，$L(\infty)$ 为常数（严格真时为 0），像缩为一点不贡献圈数。虚轴上的像正是频率响应曲线 $L(j\omega),\ \omega: -\infty \to \infty$。整理即得 $Z = N + P$。<span class="qed">∎</span></p>
<p>（$L$ 在虚轴上有极点时——如积分器——围线需以小半圆绕行，像中对应大圆弧，实验之外的技术细节。）</p></div>
<div class="box note"><span class="box-title">威力所在</span>
<ul>
<li>时滞系统 $L(s)e^{-\tau s}$：曲线随 $\tau$ 旋转卷向 $-1$，直接读出最大容忍延迟——Routh/根轨迹无能为力；</li>
<li>只需扫频实验数据即可判稳，无需模型；</li>
<li>到 $-1$ 的最近距离 = $1/\|S\|_\infty$（灵敏度峰值倒数），是《小增益定理》《H∞》的桥梁。</li>
</ul></div>`
    },
  ],
  code: {
    title: "数值验证 Nyquist 判据（画曲线、数圈数、对照闭环极点）",
    text: `import numpy as np

def nyquist_winding(num, den, K):
    """计算 K*G 的 Nyquist 曲线绕 -1 的圈数(顺时针为正)"""
    w = np.concatenate([-np.logspace(3, -3, 4000), np.logspace(-3, 3, 4000)])
    s = 1j * w
    L = K * np.polyval(num, s) / np.polyval(den, s)
    z = L + 1                      # 以 -1 为原点
    dtheta = np.diff(np.unwrap(np.angle(z)))
    return -round(dtheta.sum() / (2 * np.pi))   # 顺时针为正

num = [1.0]
den = np.polymul(np.polymul([1, 1], [1, 2]), [1, 3])  # (s+1)(s+2)(s+3)

for K in [10, 59, 61, 100]:
    N = nyquist_winding(num, den, K)
    P = 0                                    # 开环稳定
    Z_pred = N + P
    clden = den.copy(); clden[-1] += K
    Z_true = int(np.sum(np.roots(clden).real > 0))
    print(f"K={K:4}: N={N}, 预测 Z={Z_pred}, 实际右半平面闭环极点={Z_true}",
          "稳定" if Z_true == 0 else "不稳定")
# 理论临界 K=60 (Routh: 6*11-6-K*1... 可自行验证)`
  }
};

/* ════════════════ 根轨迹 ════════════════ */
window.KG_CONTENT["root-locus"] = {
  sim: "sim-rootlocus",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>调大增益，系统会怎样？Evans（1948）的根轨迹给出全景图：把闭环极点随 $K: 0 \to \infty$ 的运动轨迹一次画出。工程师看一眼就知道：增益多大开始振荡、多大失稳、往哪里加零极点能把轨迹"掰"到好位置。</p>
<p>它连接了两个世界：极点位置（时域性能）与开环零极点结构（设计自由度）。在计算机求根只要一行代码的今天，根轨迹的价值在于<strong>解释为什么</strong>——比如"为什么加零点能稳定"、"为什么高增益让相对阶 ≥ 3 的系统必然失稳"。</p>`
    },
    {
      icon: "🧾", title: "基本方程与绘制法则（含推导）",
      html: String.raw`
<div class="box def"><span class="box-title">根轨迹方程</span>
闭环特征方程 $1 + KG(s) = 0$，即 $G(s) = -\frac{1}{K}$。$K > 0$ 时分解为两个条件：
$$\textbf{辐角条件: } \angle G(s) = 180° + k\cdot360°\qquad \textbf{幅值条件: } |G(s)| = \frac{1}{K}$$
辐角条件决定轨迹形状（与 $K$ 无关！），幅值条件标定轨迹上每点对应的 $K$。</div>
<div class="box proof"><span class="box-title">主要法则及推导</span>
<p>设 $G(s) = \frac{\prod_{i=1}^{m}(s - z_i)}{\prod_{j=1}^{n}(s - p_j)}$，$n \ge m$。</p>
<p><strong>法则 1（起点与终点）</strong>：$K \to 0$ 时特征方程 $\to \prod(s - p_j) = 0$，轨迹起于开环极点；$K \to \infty$ 时 $m$ 支终于开环零点（$\prod(s-z_i) = 0$），其余 $n - m$ 支沿渐近线趋于无穷。</p>
<p><strong>法则 2（实轴段）</strong>：实轴上一点 $s_0$ 在轨迹上 $\iff$ 其<strong>右侧</strong>的实零极点总数为奇数。推导：每个右侧实零/极点对 $s_0$ 贡献辐角 $180°$，左侧贡献 $0°$，共轭复对贡献相消（角度和为 $360°$）；辐角条件要求总和为 $180°$ 的奇数倍。</p>
<p><strong>法则 3（渐近线）</strong>：$|s|$ 大时 $G(s) \approx s^{m-n}$，辐角条件给 $(m - n)\theta = 180° + k\,360°$：</p>
$$\theta_k = \frac{180° + k\,360°}{n - m},\qquad \text{交点（质心）} \ \sigma = \frac{\sum p_j - \sum z_i}{n - m}$$
<p>质心公式来自对 $1 + KG = 0$ 做 Laurent 展开保留 $s^{n-m-1}$ 项系数（等价于 Vieta 定理：闭环极点之和在 $n - m \ge 2$ 时不随 $K$ 变——极点"重心守恒"，有的往左跑必有的往右跑！）。</p>
<p><strong>法则 4（分离点）</strong>：实轴上两支相遇分离处为重根，满足 $\frac{d}{ds}\left[\frac{1}{G(s)}\right] = 0$。</p>
<p><strong>法则 5（虚轴穿越）</strong>：令 $s = j\omega$ 代入特征方程，实虚部分别为零联立求 $(K_{crit}, \omega_{crit})$——与 Routh 临界行结果一致。<span class="qed">∎</span></p></div>
<div class="box eg"><span class="box-title">例：$G(s) = \frac{1}{(s+1)(s+2)(s+3)}$（右侧实验）</span>
<p>$n - m = 3$：渐近线 $\pm 60°, 180°$，质心 $\sigma = -2$。两支必然折向右半平面。虚轴穿越：$s^3 + 6s^2 + 11s + 6 + K = 0$ 代 $s = j\omega$：实部 $-6\omega^2 + 6 + K = 0$，虚部 $-\omega^3 + 11\omega = 0$ ⇒ $\omega = \sqrt{11},\ K = 60$。这解释了一条普适规律：<strong>相对阶 $\ge 3$ 的系统增益足够大必失稳</strong>。</p></div>`
    },
  ],
  code: {
    title: "根轨迹绘制（纯 numpy）",
    text: `import numpy as np
import matplotlib.pyplot as plt

den = np.polymul(np.polymul([1, 1], [1, 2]), [1, 3])
num = [1.0]

Ks = np.logspace(-2, 3, 800)
all_roots = []
for K in Ks:
    p = den.copy().astype(float)
    p[-1] += K * num[0]
    all_roots.append(np.roots(p))
all_roots = np.array(all_roots)

plt.figure(figsize=(6, 6))
for i in range(all_roots.shape[1]):
    plt.plot(all_roots[:, i].real, all_roots[:, i].imag, '.', ms=1.5)
plt.plot([-1, -2, -3], [0, 0, 0], 'kx', ms=10, label='open-loop poles')
plt.axvline(0, color='gray', lw=0.5)
plt.axhline(0, color='gray', lw=0.5)
# 理论穿越点
plt.plot([0, 0], [np.sqrt(11), -np.sqrt(11)], 'r*', ms=12,
         label='crossing @ K=60, w=sqrt(11)')
plt.legend(); plt.xlabel('Re'); plt.ylabel('Im'); plt.title('Root Locus')
plt.axis('equal'); plt.grid(alpha=0.3); plt.show()`
  }
};

/* ════════════════ 超前滞后 ════════════════ */
window.KG_CONTENT["lead-lag"] = {
  sim: "sim-leadlag",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>PID 之外，经典设计的另一套积木：用一个"零点 + 极点"对在指定频段<strong>局部整形</strong>回路。</p>
<ul>
<li><strong>超前 (lead)</strong>：零点在前、极点在后（$z < p$）。在中频段注入<strong>正相位</strong>——直接买相位裕度，加快响应。本质是"温和版微分器"（PD 的可实现形式）；</li>
<li><strong>滞后 (lag)</strong>：极点在前、零点在后。低频抬增益（减稳态误差），中频尽量不动相位。本质是"温和版积分器"（PI 的近亲）；</li>
<li>组合使用 = 频域版 PID。区别于 PID 的意义在于：它把设计变成 Bode 图上的<strong>几何操作</strong>，每一步的因果清清楚楚。</li>
</ul>`
    },
    {
      icon: "🧾", title: "超前网络的最大相位：公式与证明",
      html: String.raw`
<div class="box def"><span class="box-title">超前环节</span>
$$C(s) = K\,\frac{Ts + 1}{\alpha Ts + 1},\qquad 0 < \alpha < 1$$
零点 $1/T$ 在极点 $1/(\alpha T)$ 左边（频率更低），中间频段相位为正。</div>
<div class="box thm"><span class="box-title">定理（最大超前角及其位置）</span>
$$\sin\varphi_{\max} = \frac{1 - \alpha}{1 + \alpha},\qquad \omega_{\max} = \frac{1}{T\sqrt{\alpha}}$$
且 $\omega_{\max}$ 是零点频率与极点频率的<strong>几何平均</strong>。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>相位 $\varphi(\omega) = \arctan(\omega T) - \arctan(\alpha\omega T)$。令 $x = \omega T$，用差角公式：</p>
$$\tan\varphi = \frac{x - \alpha x}{1 + \alpha x^2} = \frac{(1-\alpha)x}{1 + \alpha x^2}$$
<p>对 $x$ 求驻点：$\dfrac{d}{dx}\dfrac{(1-\alpha)x}{1+\alpha x^2} = \dfrac{(1-\alpha)(1 - \alpha x^2)}{(1+\alpha x^2)^2} = 0 \Rightarrow x = \dfrac{1}{\sqrt{\alpha}}$，即 $\omega_{\max} = \frac{1}{T\sqrt\alpha} = \sqrt{\frac{1}{T}\cdot\frac{1}{\alpha T}}$（几何平均，Bode 对数轴上正中点）。代回：</p>
$$\tan\varphi_{\max} = \frac{(1-\alpha)/\sqrt\alpha}{2} = \frac{1-\alpha}{2\sqrt\alpha}$$
<p>构造直角三角形（对边 $1-\alpha$，邻边 $2\sqrt\alpha$，斜边 $\sqrt{(1-\alpha)^2 + 4\alpha} = 1+\alpha$）得 $\sin\varphi_{\max} = \frac{1-\alpha}{1+\alpha}$。<span class="qed">∎</span></p></div>
<div class="box eg"><span class="box-title">标准设计流程（实验可复现）</span>
<ol>
<li>按稳态误差要求定增益 $K$；</li>
<li>测未校正 PM，算相位缺口 $\varphi_{need}$（加 5°~10° 余量，因为校正会推高 $\omega_{gc}$）；</li>
<li>$\alpha = \frac{1 - \sin\varphi_{need}}{1 + \sin\varphi_{need}}$；</li>
<li>把 $\omega_{\max}$ 放在<strong>新的</strong>穿越频率上（该处校正环节幅值为 $-10\lg\alpha$ dB，找 $|L|$ 等于其相反数的频率），$T = \frac{1}{\omega_{\max}\sqrt\alpha}$。</li>
</ol>
<p>注意 $\varphi_{\max} < 90°$（单级），要 $70°$ 以上就得两级串联——没有免费的相位。</p></div>`
    },
  ],
  code: {
    title: "按流程自动设计超前校正器",
    text: `import numpy as np

def freq_resp(num, den, w):
    s = 1j * w
    return np.polyval(num, s) / np.polyval(den, s)

def margins(num, den):
    w = np.logspace(-3, 3, 8000)
    L = freq_resp(num, den, w)
    mag, ph = abs(L), np.unwrap(np.angle(L)) * 180/np.pi
    i = np.argmin(abs(mag - 1))
    return w[i], 180 + ph[i]

# 对象 G = 4/(s(s+1)), 目标 PM = 50 deg
Gn, Gd = [4.0], np.polymul([1, 0], [1, 1])
wgc0, pm0 = margins(Gn, Gd)
print(f"未校正: PM={pm0:.1f} deg @ {wgc0:.2f} rad/s")

phi_need = np.radians(50 - pm0 + 8)          # 缺口 + 余量
alpha = (1 - np.sin(phi_need)) / (1 + np.sin(phi_need))
# 新穿越频率: |G| = 10*lg(alpha)/2 dB 处
w = np.logspace(-2, 2, 8000)
magG = abs(freq_resp(Gn, Gd, w))
target = 1 / np.sqrt(1/alpha)                # = sqrt(alpha) 倍
w_new = w[np.argmin(abs(magG - 1/np.sqrt(alpha)**-1 * np.sqrt(1/alpha)*0 - np.sqrt(alpha)**-1))]
w_new = w[np.argmin(abs(magG * np.sqrt(1/alpha) - 1))]   # |G|*|C|max = 1
T = 1 / (w_new * np.sqrt(alpha))

Cn, Cd = [T, 1], [alpha*T, 1]
Ln, Ld = np.polymul(Cn, Gn), np.polymul(Cd, Gd)
wgc1, pm1 = margins(Ln, Ld)
print(f"alpha={alpha:.3f}, T={T:.3f}")
print(f"校正后: PM={pm1:.1f} deg @ {wgc1:.2f} rad/s")`
  }
};
