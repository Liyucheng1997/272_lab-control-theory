/* 内容：EKF、粒子滤波、高阶滑模 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ EKF ════════════════ */
window.KG_CONTENT["ekf"] = {
  sim: "sim-ekf",
  sections: [
    {
      icon: "💡", title: "直觉：机器人为什么离不开它",
      html: String.raw`
<p>Kalman 滤波要求模型线性——但机器人几乎处处非线性：轮式机器人的运动学里有 $\cos\theta, \sin\theta$，激光/相机测的是距离和方位角（含 $\mathrm{atan2}$ 和平方根），无人机姿态在旋转群上。直接套 KF？公式都写不出来。</p>
<p>EKF 的策略朴素而有效：<strong>每一步都在"当前最优估计点"处做一次线性化（Jacobian），然后假装自己是线性系统，套用标准 Kalman 公式</strong>。相当于用一连串切平面去逼近弯曲的流形——只要估计误差不大、非线性不剧烈，就相当好用。</p>
<p>从阿波罗登月的导航（EKF 的首秀！）到今天的机器人 SLAM、无人机 GPS/IMU 组合导航、自动驾驶定位，EKF 统治非线性估计半个多世纪。它与《线性化》的思想一脉相承：那里是"控制器设计前先线性化"，这里是"每个滤波周期都重新线性化一次"。</p>`
    },
    {
      icon: "📐", title: "算法定义",
      html: String.raw`
<div class="box def"><span class="box-title">非线性模型</span>
$$x_k = f(x_{k-1}, u_{k-1}) + w_{k-1},\qquad z_k = h(x_k) + v_k$$
$w \sim \mathcal{N}(0, Q)$，$v \sim \mathcal{N}(0, R)$。定义两个 Jacobian（<strong>在估计点处评估</strong>，这是 EKF 的全部秘密）：
$$F_k = \frac{\partial f}{\partial x}\Big|_{\hat{x}_{k-1|k-1}},\qquad H_k = \frac{\partial h}{\partial x}\Big|_{\hat{x}_{k|k-1}}$$</div>
<div class="box def"><span class="box-title">EKF 算法（与 KF 逐行对照）</span>
<p><strong>预测</strong>（均值走非线性真模型，协方差走线性化模型）：</p>
$$\hat{x}_{k|k-1} = f(\hat{x}_{k-1|k-1}, u_{k-1}),\qquad P_{k|k-1} = F_kP_{k-1|k-1}F_k^{\top} + Q$$
<p><strong>更新</strong>（新息用非线性 $h$，增益用 Jacobian）：</p>
$$K_k = P_{k|k-1}H_k^{\top}\big(H_kP_{k|k-1}H_k^{\top} + R\big)^{-1}$$
$$\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k\big(z_k - h(\hat{x}_{k|k-1})\big),\qquad P_{k|k} = (I - K_kH_k)P_{k|k-1}$$</div>`
    },
    {
      icon: "🧾", title: "推导：误差动态的一阶展开",
      html: String.raw`
<div class="box proof"><span class="box-title">推导（为什么协方差这样传播）</span>
<p>记预测误差 $e_{k|k-1} = x_k - \hat{x}_{k|k-1}$。对真实动态在估计点做 Taylor 展开：</p>
$$x_k = f(x_{k-1}) + w = f(\hat{x}_{k-1|k-1}) + F_k\,e_{k-1|k-1} + \underbrace{O(\|e\|^2)}_{\text{丢弃}} + w$$
<p>减去预测均值 $\hat{x}_{k|k-1} = f(\hat{x}_{k-1|k-1})$：</p>
$$e_{k|k-1} \approx F_k\,e_{k-1|k-1} + w_{k-1}$$
<p>——误差（近似）满足一个<strong>线性时变</strong>系统！于是协方差按线性规则传播：$P_{k|k-1} = F_kP F_k^{\top} + Q$。测量同理：$z_k - h(\hat{x}_{k|k-1}) \approx H_k e_{k|k-1} + v_k$。</p>
<p>结论：<strong>EKF = 把标准 Kalman 滤波逐步应用于"沿估计轨迹线性化"的时变系统</strong>。《Kalman 滤波》里证明的最优增益公式逐条照搬——但注意，最优性证明依赖的"误差动态精确线性"在这里只是近似。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">EKF 不是最优的——什么时候会翻车</span>
<ul>
<li><strong>丢掉的 $O(\|e\|^2)$ 项</strong>是全部误差来源：初始误差大、非线性强（转弯急、传感器近距离）时，线性化点远离真值，滤波器会<strong>过度自信</strong>（$P$ 偏小）甚至发散；</li>
<li>Jacobian 在<strong>估计点</strong>而非真值处评估——估计错了连线性化都错，这是正反馈式风险；</li>
<li>改良谱系：<strong>IEKF</strong>（更新内迭代重新线性化）、<strong>UKF</strong>（sigma 点传播，免求导、二阶精度）、<strong>ESKF</strong>（误差状态，姿态估计标配）、<strong>粒子滤波</strong>（放弃高斯假设，见下一站）；</li>
<li>工程守则：谨慎初始化、监控新息（归一化新息平方 NIS 应服从 $\chi^2$ 分布）、必要时人为放大 $Q$。</li>
</ul></div>`
    },
    {
      icon: "🔬", title: "机器人定位算例（与实验对应）",
      html: String.raw`
<div class="box eg"><span class="box-title">差速轮机器人 + 路标测距测角</span>
<p>状态 $x = [p_x, p_y, \theta]^{\top}$，运动模型（里程计输入 $v, \omega$）：</p>
$$f:\quad p_x^+ = p_x + v\cos\theta\,\Delta t,\quad p_y^+ = p_y + v\sin\theta\,\Delta t,\quad \theta^+ = \theta + \omega\Delta t$$
$$F = \begin{bmatrix} 1 & 0 & -v\sin\theta\,\Delta t \\ 0 & 1 & v\cos\theta\,\Delta t \\ 0 & 0 & 1 \end{bmatrix}$$
<p>对位于 $(l_x, l_y)$ 的路标，测量距离与方位：</p>
$$h(x) = \begin{bmatrix} \sqrt{q} \\ \mathrm{atan2}(l_y - p_y,\ l_x - p_x) - \theta \end{bmatrix},\quad q = (l_x - p_x)^2 + (l_y - p_y)^2$$
$$H = \begin{bmatrix} -\frac{l_x - p_x}{\sqrt{q}} & -\frac{l_y - p_y}{\sqrt{q}} & 0 \\ \frac{l_y - p_y}{q} & -\frac{l_x - p_x}{q} & -1 \end{bmatrix}$$
<p>右侧实验：只靠里程计航位推算（dead reckoning）误差随时间累积发散；EKF 每次看到路标就把漂移"拉"回来。调大里程计噪声看差距。</p></div>`
    },
  ],
  code: {
    title: "EKF 机器人定位（差速轮 + 路标测距测角，可直接运行）",
    text: `import numpy as np

dt = 0.1
landmarks = np.array([[4., 4.], [-4., 3.], [3., -4.]])
Q = np.diag([0.02, 0.02, 0.01])**2       # 过程噪声
R = np.diag([0.25, 0.08])**2             # 距离/方位噪声

def f(x, u):
    v, w = u
    return x + np.array([v*np.cos(x[2])*dt, v*np.sin(x[2])*dt, w*dt])

def F_jac(x, u):
    v, _ = u
    return np.array([[1, 0, -v*np.sin(x[2])*dt],
                     [0, 1,  v*np.cos(x[2])*dt],
                     [0, 0, 1]])

def h(x, lm):
    d = lm - x[:2]; q = d @ d
    return np.array([np.sqrt(q), np.arctan2(d[1], d[0]) - x[2]])

def H_jac(x, lm):
    d = lm - x[:2]; q = d @ d; sq = np.sqrt(q)
    return np.array([[-d[0]/sq, -d[1]/sq, 0],
                     [ d[1]/q,  -d[0]/q, -1]])

rng = np.random.default_rng(2)
xt = np.array([0., 0., 0.])              # 真值
xe = np.array([0., 0., 0.]); P = np.eye(3)*0.1
xd = xe.copy()                            # 纯航位推算对照

for k in range(600):
    u = np.array([1.0, 0.5])                       # 画圆
    u_noisy = u + rng.normal(0, [0.05, 0.03])      # 带噪里程计
    xt = f(xt, u)
    xd = f(xd, u_noisy)
    # EKF 预测
    xe = f(xe, u_noisy)
    Fk = F_jac(xe, u_noisy)
    P = Fk @ P @ Fk.T + Q
    # 每 5 步观测一个路标
    if k % 5 == 0:
        lm = landmarks[(k//5) % 3]
        z = h(xt, lm) + rng.normal(0, [0.25, 0.08])
        Hk = H_jac(xe, lm)
        S = Hk @ P @ Hk.T + R
        K = P @ Hk.T @ np.linalg.inv(S)
        innov = z - h(xe, lm)
        innov[1] = (innov[1] + np.pi) % (2*np.pi) - np.pi   # 角度归一化!
        xe = xe + K @ innov
        P = (np.eye(3) - K @ Hk) @ P

print("末端误差  航位推算:", np.linalg.norm(xd[:2]-xt[:2]).round(3),
      "  EKF:", np.linalg.norm(xe[:2]-xt[:2]).round(3))`
  }
};

/* ════════════════ 粒子滤波 ════════════════ */
window.KG_CONTENT["particle-filter"] = {
  sim: "sim-pf",
  sections: [
    {
      icon: "💡", title: "直觉：用一群粒子代表信念",
      html: String.raw`
<p>EKF 有两个放不下的执念：<strong>分布是高斯的、线性化是够用的</strong>。可机器人的真实处境经常打脸——"绑架的机器人"睁眼一看走廊，长得一样的门有三扇，它的位置信念是<strong>三个峰</strong>的分布；高斯（单峰椭圆）根本表达不了。</p>
<p>粒子滤波干脆放弃解析形式：<strong>用 $N$ 个带权样本（粒子）$\{x^{(i)}, w^{(i)}\}$ 直接近似任意形状的后验分布</strong>。每个粒子是一个"平行宇宙的假设"：预测时每个粒子各自按动力学（加噪声）演化；观测到来时，按"该假设对观测的解释力"$p(z|x^{(i)})$ 重新称重；权重太集中时重采样，把算力聚焦到高可能区域。</p>
<p>三步循环——<strong>传播、称重、重采样</strong>——就是全部。它能表示多峰、能用任何非线性模型、实现只要几十行，代价是算力随维数暴涨。机器人全局定位（MCL/AMCL，ROS 导航栈标配）、SLAM（FastSLAM）、目标跟踪是它的主场。</p>`
    },
    {
      icon: "📐", title: "算法定义",
      html: String.raw`
<div class="box def"><span class="box-title">贝叶斯滤波递推（一切滤波器的母公式）</span>
$$\underbrace{p(x_k | z_{1:k})}_{\text{后验}} \propto \underbrace{p(z_k | x_k)}_{\text{观测似然}} \int \underbrace{p(x_k | x_{k-1})}_{\text{运动模型}}\, p(x_{k-1} | z_{1:k-1})\,dx_{k-1}$$
<p>KF 是它在"线性 + 高斯"下的解析解；EKF 是近似解析解；粒子滤波是<strong>蒙特卡洛数值解</strong>。</p></div>
<div class="box def"><span class="box-title">SIR 粒子滤波（Sequential Importance Resampling）</span>
<ol>
<li><strong>传播</strong>：$x_k^{(i)} \sim p(x_k | x_{k-1}^{(i)}, u_k)$（把每个粒子推过带噪声的运动模型）；</li>
<li><strong>称重</strong>：$w_k^{(i)} \propto w_{k-1}^{(i)}\, p(z_k | x_k^{(i)})$，然后归一化 $\sum_i w_k^{(i)} = 1$；</li>
<li><strong>重采样</strong>：当有效粒子数 $N_{\mathrm{eff}} = 1/\sum_i (w_k^{(i)})^2$ 过小（如 $< N/2$）时，按权重有放回抽取 $N$ 个粒子，权重重置为 $1/N$；</li>
<li>估计：$\hat{x}_k = \sum_i w_k^{(i)}x_k^{(i)}$（或取最大权重簇）。</li>
</ol></div>`
    },
    {
      icon: "🧾", title: "为什么有效：重要性采样的数学",
      html: String.raw`
<div class="box thm"><span class="box-title">引理（重要性采样）</span>
<p>设目标分布 $p$ 难以采样，提议分布 $q$ 易采样且 $p > 0 \Rightarrow q > 0$。取 $x^{(i)} \sim q$，权重 $w(x) = p(x)/q(x)$，则对任意函数 $\varphi$：</p>
$$\mathbb{E}_q\big[w(x)\varphi(x)\big] = \mathbb{E}_p[\varphi(x)]$$
<p>即加权样本均值 $\sum_i w^{(i)}\varphi(x^{(i)})/N$ 是 $\mathbb{E}_p[\varphi]$ 的无偏估计，且由大数定律以 $O(1/\sqrt{N})$（<strong>与维数无关</strong>的速率常数除外）收敛。</p></div>
<div class="box proof"><span class="box-title">证明</span>
$$\mathbb{E}_q[w\varphi] = \int \frac{p(x)}{q(x)}\varphi(x)\,q(x)\,dx = \int \varphi(x)\,p(x)\,dx = \mathbb{E}_p[\varphi] \qquad ∎$$</div>
<div class="box proof"><span class="box-title">推导：递推权重公式从哪来</span>
<p>目标是整条轨迹的后验 $p(x_{0:k}|z_{1:k})$。用链式提议分布 $q(x_{0:k}) = q(x_{0:k-1})\,q(x_k|x_{k-1})$ 逐步采样，权重按贝叶斯展开：</p>
$$w_k \propto \frac{p(x_{0:k}|z_{1:k})}{q(x_{0:k})} = \frac{p(z_k|x_k)\,p(x_k|x_{k-1})\,p(x_{0:k-1}|z_{1:k-1})}{q(x_k|x_{k-1})\,q(x_{0:k-1})} = w_{k-1}\cdot\frac{p(z_k|x_k)\,p(x_k|x_{k-1})}{q(x_k|x_{k-1})}$$
<p>（第一步用了贝叶斯法则与马尔可夫性。）SIR 取<strong>提议 = 运动模型</strong> $q(x_k|x_{k-1}) = p(x_k|x_{k-1})$，分子分母相消，只剩</p>
$$w_k \propto w_{k-1}\,p(z_k|x_k)$$
<p>——这就是"按观测似然称重"的严格出处。<span class="qed">∎</span></p></div>
<div class="box proof"><span class="box-title">退化问题与 $N_{\mathrm{eff}}$ 的由来</span>
<p>可以证明（Doucet et al.）序贯重要性采样的权重方差<strong>单调增长</strong>：迟早只剩一个粒子权重 $\approx 1$，其余全是算力浪费。度量退化的常用量是有效样本数：把加权估计的方差与 $N_{\mathrm{eff}}$ 个等权样本的方差对齐，可得近似</p>
$$N_{\mathrm{eff}} \approx \frac{1}{\sum_i (w^{(i)})^2}$$
<p>（两个极端自查：权重均匀时 $=N$；单粒子独大时 $=1$。）重采样以引入少量额外方差为代价，斩断权重方差的累积——这是粒子滤波能长期运行的关键。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">工程要点</span>
<ul>
<li><strong>粒子贫化</strong>：重采样后多样性下降；对策：低方差(系统)重采样、只在 $N_{\mathrm{eff}}$ 低时重采样、加抖动（roughening）；</li>
<li><strong>维数灾难</strong>：所需粒子数随状态维数指数增长，高维用 Rao-Blackwell 化（FastSLAM：粒子只管位姿，地图用解析 EKF）；</li>
<li>观测太准（$R$ 很小）反而危险：似然尖峰导致几乎所有粒子权重为零——实验里把测量噪声调到最小试试。</li>
</ul></div>`
    },
  ],
  code: {
    title: "粒子滤波机器人定位（仅测距 → 多峰后验，40 行核心）",
    text: `import numpy as np

rng = np.random.default_rng(3)
N = 1000
landmark = np.array([0., 0.])
R_std = 0.3
dt = 0.1

def motion(X, u):
    v, w = u
    n = len(X)
    X = X.copy()
    X[:, 2] += (w + rng.normal(0, 0.05, n)) * dt
    vn = v + rng.normal(0, 0.08, n)
    X[:, 0] += vn * np.cos(X[:, 2]) * dt
    X[:, 1] += vn * np.sin(X[:, 2]) * dt
    return X

def likelihood(X, z):
    d = np.linalg.norm(X[:, :2] - landmark, axis=1)
    return np.exp(-0.5 * ((z - d) / R_std)**2) + 1e-300

def systematic_resample(X, w):
    pos = (rng.random() + np.arange(N)) / N
    idx = np.searchsorted(np.cumsum(w), pos)
    return X[np.clip(idx, 0, N-1)]

# 初始: 完全未知位置(全图撒粒子) -> 全局定位
X = np.column_stack([rng.uniform(-6, 6, N), rng.uniform(-6, 6, N),
                     rng.uniform(-np.pi, np.pi, N)])
w = np.ones(N) / N
xt = np.array([2., -1., 0.7])            # 真实位姿

for k in range(300):
    u = np.array([1.0, 0.4])
    xt = motion(xt[None, :], u)[0]        # 真值也演化(带噪)
    X = motion(X, u)
    z = np.linalg.norm(xt[:2] - landmark) + rng.normal(0, R_std)
    w = w * likelihood(X, z)
    w /= w.sum()
    neff = 1.0 / np.sum(w**2)
    if neff < N / 2:
        X = systematic_resample(X, w)
        w = np.ones(N) / N

est = (X * w[:, None]).sum(axis=0)
print("真实位置:", xt[:2].round(2), " 粒子估计:", est[:2].round(2))
print("误差:", np.linalg.norm(est[:2] - xt[:2]).round(3), "m")`
  }
};

/* ════════════════ 高阶滑模 / Super-Twisting ════════════════ */
window.KG_CONTENT["hosm"] = {
  sim: "sim-hosm",
  sections: [
    {
      icon: "💡", title: "直觉：把不连续藏进积分器",
      html: String.raw`
<p>《滑模控制》的结论：鲁棒性极好，但控制量 $u = -k\,\mathrm{sign}(s)$ 高频抖振，电机、液压阀、机翼都受不了；边界层法能磨平抖振，却牺牲了精度（只收敛到 $O(\varphi)$ 邻域）和不变性。<strong>鱼与熊掌能否兼得？</strong></p>
<p>高阶滑模的核心思想：抖振的根源是"不连续项直接出现在执行器指令里"。那就把它<strong>下移一层——藏进控制量的导数里</strong>：$\dot{u}$ 里含 $\mathrm{sign}(s)$，积分之后送给执行器的 $u$ 就是连续的！同时用分数次幂项 $|s|^{1/2}\mathrm{sign}(s)$ 保住有限时间收敛。</p>
<p>这就是 <strong>Super-Twisting 算法</strong>（Levant, 1993）——机器人和电机驱动领域引用最多的滑模变体：连续控制量、有限时间把 $s$ <strong>和 $\dot{s}$ 同时</strong>压到零（二阶滑模）、只需知道扰动导数的界。工程上几乎是"想用滑模又怕抖振"的默认答案，也是滑模微分器（实时求导的鲁棒方法）的内核。</p>`
    },
    {
      icon: "📐", title: "定义",
      html: String.raw`
<div class="box def"><span class="box-title">定义（$r$ 阶滑模）</span>
<p>称系统在滑模变量 $s$ 上建立 <strong>$r$ 阶滑模</strong>，若有限时间内到达并保持</p>
$$s = \dot{s} = \ddot{s} = \cdots = s^{(r-1)} = 0$$
<p>经典滑模是 $r = 1$（只保证 $s = 0$，$\dot s$ 靠切换来回穿越）；Super-Twisting 实现 $r = 2$：$s$ 与 $\dot s$ 同时归零，收敛更"干净"，稳态精度对采样步长是 $O(\tau^2)$（一阶滑模只有 $O(\tau)$）。</p></div>
<div class="box def"><span class="box-title">定义（Super-Twisting 算法）</span>
<p>对相对阶为 1 的滑模动态 $\dot{s} = \rho(t) + u$（$\rho$ 为集总扰动，满足 $|\dot{\rho}| \le L$）：</p>
$$u = -k_1\,|s|^{1/2}\mathrm{sign}(s) + v,\qquad \dot{v} = -k_2\,\mathrm{sign}(s)$$
<p>两个成分各司其职：$|s|^{1/2}$ 项提供"越接近面越温柔、但永不失速"的推力（有限时间的来源）；积分项 $v$ 负责<strong>连续地</strong>顶住扰动 $\rho$（稳态时 $v \to -\rho$，相当于内置了一个扰动观测器）。注意 $u$ 本身连续——sign 只出现在 $\dot{v}$ 里。</p></div>`
    },
    {
      icon: "🧾", title: "定理：有限时间收敛及其证明",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（Super-Twisting 收敛性，Levant / Moreno）</span>
<p>设 $|\dot\rho| \le L$。若增益满足（一组充分条件）</p>
$$k_2 > L,\qquad k_1^2 > 2(k_2 + L)\cdot\frac{k_2 + L}{k_2 - L}\cdot\frac{1}{2}\ \ \big(\text{常用工程整定 } k_1 = 1.5\sqrt{L},\ k_2 = 1.1L\big)$$
<p>则 $(s, \dot{s})$ 在有限时间内收敛到 $(0, 0)$ 并保持——即建立二阶滑模。</p></div>
<div class="box proof"><span class="box-title">证明（Moreno–Osorio 的 Lyapunov 路线，主干完整）</span>
<p><strong>第一步：换坐标。</strong>令 $\zeta = \begin{bmatrix} \zeta_1 \\ \zeta_2 \end{bmatrix} = \begin{bmatrix} |s|^{1/2}\mathrm{sign}(s) \\ v + \rho \end{bmatrix}$（$\zeta_2$ 是"积分器还差多少才顶平扰动"）。注意 $\dot s = \zeta_2 - k_1\zeta_1$，且 $\dot{\zeta}_1 = \frac{\dot{s}}{2|s|^{1/2}}$。于是在 $s \ne 0$ 处：</p>
$$\dot{\zeta} = \frac{1}{|\zeta_1|}\begin{bmatrix} -\frac{k_1}{2} & \frac{1}{2} \\ -k_2 & 0 \end{bmatrix}\zeta + \begin{bmatrix} 0 \\ \dot\rho \end{bmatrix},\qquad A \triangleq \begin{bmatrix} -\frac{k_1}{2} & \frac{1}{2} \\ -k_2 & 0 \end{bmatrix}$$
<p>——除掉标量因子 $1/|\zeta_1| > 0$，这是一个<strong>线性系统</strong>！$A$ 的特征式 $\lambda^2 + \frac{k_1}{2}\lambda + \frac{k_2}{2}$，$k_1, k_2 > 0$ 时 Hurwitz。</p>
<p><strong>第二步：二次型 Lyapunov 函数。</strong>取 $A^{\top}P + PA = -I$ 的解 $P \succ 0$（存在性见《线性系统稳定性》），令 $V(\zeta) = \zeta^{\top}P\zeta$。沿轨迹（先看无扰动 $\dot\rho = 0$ 情形）：</p>
$$\dot{V} = \frac{1}{|\zeta_1|}\zeta^{\top}(A^{\top}P + PA)\zeta = -\frac{\|\zeta\|^2}{|\zeta_1|}$$
<p><strong>第三步：从负定导数到有限时间。</strong>由 $|\zeta_1| \le \|\zeta\| \le \sqrt{V/\lambda_{\min}(P)}$：</p>
$$\dot{V} \le -\frac{\|\zeta\|^2}{\|\zeta\|} = -\|\zeta\| \le -\sqrt{\frac{V}{\lambda_{\max}(P)}} = -\gamma V^{1/2},\quad \gamma = \lambda_{\max}(P)^{-1/2}$$
<p>分离变量积分：$V^{1/2}(t) \le V^{1/2}(0) - \frac{\gamma}{2}t$，在 $T \le \frac{2V^{1/2}(0)}{\gamma}$ 时刻 $V = 0$，即 $\zeta = 0$，亦即 $s = 0$ 且 $v = -\rho$（从而 $\dot s = 0$）。<strong>幂次 $1/2$ 正是有限时间的引擎</strong>——对照普通渐近稳定的 $\dot V \le -cV$ 只能给指数衰减（回忆《滑模控制》到达条件同款论证）。</p>
<p><strong>第四步：带扰动。</strong>$\dot\rho \ne 0$ 时 $\dot V$ 多出交叉项 $2\zeta^{\top}P[0, \dot\rho]^{\top}$，用 $|\dot\rho| \le L$ 上界化后，可验证当 $k_1, k_2$ 满足定理条件时二次型 $-\frac{1}{|\zeta_1|}\zeta^{\top}(I - \Delta_L)\zeta$ 仍负定（$\Delta_L$ 为扰动贡献的对称阵），同样得 $\dot V \le -\gamma' V^{1/2}$。有限时间结论保持。技术细节（$s = 0$ 处的不可微性）由 Filippov 微分包含意义下的非光滑分析补全，结构不变。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">对比与家族</span>
<ul>
<li><strong>vs 边界层</strong>：两者控制量都连续，但边界层是"放弃精度换平滑"（稳态误差 $O(\varphi)$、失去不变性），Super-Twisting 是"精确收敛 + 平滑"兼得，只多要一个 $|\dot\rho|$ 的界；</li>
<li><strong>Twisting 算法</strong>（需要 $\dot s$ 测量）、<strong>任意阶 HOSM</strong>（Levant 的 quasi-continuous 族）、<strong>滑模微分器</strong>：$\dot{\hat{x}}_0 = -\lambda_1|{\hat x}_0 - f|^{1/2}\mathrm{sign}(\cdot) + \hat x_1, \dot{\hat x}_1 = -\lambda_2\mathrm{sign}(\cdot)$——就是 Super-Twisting 换了个用途，鲁棒实时求导的利器；</li>
<li>实验建议：先把边界层滑模的 $\varphi$ 调大消抖振、观察稳态误差；再切 Super-Twisting——$u$ 连续且 $s$ 精确归零。</li>
</ul></div>`
    },
  ],
  code: {
    title: "Super-Twisting vs 一阶滑模：抖振与精度同台对比",
    text: `import numpy as np

L = 2.4                       # |d rho/dt| 的界
k1, k2 = 1.5*np.sqrt(L), 1.1*L
rho = lambda t: 0.8*np.sin(3*t)     # 集总扰动(未知)

def run(mode, T=8.0, h=1e-4):
    s, v = 1.5, 0.0
    S, U = [], []
    for i in range(int(T/h)):
        t = i*h
        if mode == "first-order":
            u = -1.5 * np.sign(s)             # 不连续 -> 抖振
        else:                                  # super-twisting
            u = -k1*np.sqrt(abs(s))*np.sign(s) + v
            v += h * (-k2*np.sign(s))          # sign 藏在导数里
        s += h * (rho(t) + u)
        if t > 4: S.append(abs(s)); U.append(u)
    U = np.array(U)
    tv = np.abs(np.diff(U)).sum()             # 控制量总变差(抖振度量)
    return np.mean(S), tv

for mode in ["first-order", "super-twisting"]:
    err, tv = run(mode)
    print(f"{mode:14s}  稳态|s|均值={err:.2e}   u总变差={tv:8.1f}")
# super-twisting: 误差低几个数量级, u 平滑(总变差小 ~百倍)`
  }
};
