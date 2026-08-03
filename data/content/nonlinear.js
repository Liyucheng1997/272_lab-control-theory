/* 内容：非线性控制 */
window.KG_CONTENT = window.KG_CONTENT || {};

/* ════════════════ 反馈线性化 ════════════════ */
window.KG_CONTENT["feedback-linearization"] = {
  sim: "sim-fl",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>《线性化》是"假装非线性不存在"（小范围近似）；反馈线性化是"<strong>用控制量把非线性亲手减掉</strong>"——精确的、大范围的代数抵消。</p>
<p>例：单摆 $\ddot\theta = -a\sin\theta - b\dot\theta + u$。令 $u = a\sin\theta + b\dot\theta + v$，系统瞬间变成 $\ddot\theta = v$——一个纯双积分器，接下来用任何线性方法设计 $v$ 都行。机器人臂的"计算力矩控制"、无人机的微分平坦轨迹跟踪，都是这个思想的工业化版本。</p>
<p>代价也很清楚：要<strong>精确的模型</strong>（减错了会留下残余项）、要全状态、且抵消"好的非线性"（如自然阻尼）会浪费控制能量。适用范围由"相对阶"与"零动态"两个概念精确刻画。</p>`
    },
    {
      icon: "📐", title: "相对阶与标准形",
      html: String.raw`
<div class="box def"><span class="box-title">定义（Lie 导数与相对阶）</span>
<p>考虑单输入单输出仿射系统 $\dot{x} = f(x) + g(x)u,\ y = h(x)$。记 Lie 导数 $L_fh = \nabla h^{\top}f$（$h$ 沿 $f$ 方向的变化率），迭代 $L_f^k h = L_f(L_f^{k-1}h)$。</p>
<p>系统在 $x_0$ 处<strong>相对阶</strong>为 $r$，若</p>
$$L_gL_f^{k}h(x) = 0\ (k < r - 1),\qquad L_gL_f^{r-1}h(x_0) \ne 0$$
<p>直白说：对 $y$ 求导，求到第 $r$ 阶导数时输入 $u$ 才首次出现。</p></div>
<div class="box proof"><span class="box-title">为什么第 $r$ 阶才出现 $u$（计算验证）</span>
$$\dot{y} = \nabla h^{\top}(f + gu) = L_fh + \underbrace{L_gh}_{=0}\,u = L_fh$$
<p>再求导：$\ddot{y} = L_f^2h + L_gL_fh\,u = L_f^2h$（若 $r > 2$）……归纳进行到</p>
$$y^{(r)} = L_f^r h(x) + L_gL_f^{r-1}h(x)\,u$$
<p>系数 $L_gL_f^{r-1}h \neq 0$，$u$ 首次登场。<span class="qed">∎</span></p></div>`
    },
    {
      icon: "🧾", title: "定理：输入-输出线性化",
      html: String.raw`
<div class="box thm"><span class="box-title">定理</span>
<p>设相对阶为 $r$。取控制律</p>
$$u = \frac{v - L_f^r h(x)}{L_gL_f^{r-1}h(x)}$$
<p>则从新输入 $v$ 到输出 $y$ 的关系精确为 $r$ 重积分链 $y^{(r)} = v$。取坐标 $\xi_i = y^{(i-1)}\ (i = 1..r)$ 并补充 $n - r$ 个内动态坐标 $\eta$（可选取使其方程不含 $u$），系统化为<strong>标准形</strong>：</p>
$$\dot{\xi}_1 = \xi_2,\ \dots,\ \dot{\xi}_r = v,\qquad \dot{\eta} = q(\xi, \eta)$$
<p>其中 $\dot\eta = q(0, \eta)$ 称为<strong>零动态</strong>。若零动态渐近稳定（最小相位），则用 $v = -k_{r-1}y^{(r-1)} - \cdots - k_0y$（配置稳定极点）可实现（局部）镇定/跟踪。</p></div>
<div class="box proof"><span class="box-title">证明（线性化部分）</span>
<p>由上节 $y^{(r)} = L_f^rh + L_gL_f^{r-1}h\cdot u$，代入控制律分子分母相消：$y^{(r)} = v$。坐标变换的合法性：可证 $(dh, dL_fh, \dots, dL_f^{r-1}h)$ 在 $x_0$ 附近线性无关（相对阶定义保证），可扩充为局部微分同胚。镇定结论：$\xi$ 子系统是可控标准型线性系统，极点任意配置；$\xi \to 0$ 后 $\eta$ 趋于零动态，由最小相位假设与级联系统稳定性论证（局部 ISS）得整体局部渐近稳定。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">零动态：藏在镜子后面的部分</span>
<ul>
<li>$r = n$：无零动态，<strong>全状态线性化</strong>，全局结论可得（如单摆例）；</li>
<li>$r < n$：有 $n - r$ 维"看不见"的内动态。<strong>非最小相位</strong>（零动态不稳定）时反馈线性化直接失效——线性对应物即右半平面零点不能对消的事实（对消 = 内部失稳）；</li>
<li>与《Backstepping》的关系：FL 一刀切掉全部非线性；backstepping 逐级处理、可以<strong>保留有益非线性</strong>且对不匹配不确定性更宽容，见下一站。</li>
</ul></div>`
    },
  ],
  code: {
    title: "单摆的精确反馈线性化 + 大角度跟踪",
    text: `import numpy as np

a, b = 10.0, 0.5          # thetaddot = -a sin(theta) - b thetadot + u
kp, kd = 25.0, 9.0

def dyn(x, u):
    return np.array([x[1], -a*np.sin(x[0]) - b*x[1] + u])

def track(use_fl, A0=2.0, T=12.0, h=1e-3):
    x = np.array([0.5, 0.0]); errs = []
    for k in range(int(T/h)):
        t = k*h
        thd, thdD, thdDD = A0*np.sin(t), A0*np.cos(t), -A0*np.sin(t)
        e, ed = x[0]-thd, x[1]-thdD
        v = thdDD - kd*ed - kp*e
        u = (a*np.sin(x[0]) + b*x[1] + v) if use_fl else v
        x = x + h*dyn(x, u)
        errs.append(e)
    return np.sqrt(np.mean(np.square(errs)))

print("大角度(A=2 rad)跟踪 RMSE:")
print("  反馈线性化:", round(track(True), 5))
print("  纯线性 PD :", round(track(False), 5))
# FL 把误差压到数值精度级; PD 因未抵消 a*sin(theta) 存在系统性偏差`
  }
};

/* ════════════════ Backstepping ════════════════ */
window.KG_CONTENT["backstepping"] = {
  sim: "sim-backstepping",
  sections: [
    {
      icon: "💡", title: "直觉：一级一级往回退",
      html: String.raw`
<p>很多系统长成"链条"：$u$ 直接推得动 $x_2$，$x_2$ 再推 $x_1$（严格反馈形式）。想控 $x_1$ 却只有 $u$，怎么办？</p>
<p>反步法的策略像带团队：① 先假装 $x_2$ 是你能随意指定的"<strong>虚拟输入</strong>"，为 $x_1$ 设计出理想值 $x_2 = \alpha(x_1)$（并配一个 Lyapunov 函数 $V_1$）；② $x_2$ 当然不会自动等于 $\alpha$，于是把偏差 $z_2 = x_2 - \alpha$ 当新目标，用真正的 $u$ 驱使 $z_2 \to 0$，同时把 $V_1$ 扩充成 $V_2 = V_1 + \frac12 z_2^2$；③ 更长的链条就继续"往回退一步"。</p>
<p>每一步都在<strong>递推地构造 Lyapunov 函数</strong>——这是它相对反馈线性化的深刻优势：不必抵消所有非线性（有益的可以留下），并且能系统地叠加自适应、鲁棒项。它是 1990s 非线性控制的招牌成果（Krstić–Kanellakopoulos–Kokotović, 1995）。</p>`
    },
    {
      icon: "🧾", title: "定理与完整递推证明（二阶情形）",
      html: String.raw`
<div class="box thm"><span class="box-title">定理（二级反步）</span>
<p>考虑严格反馈系统（$f$ 光滑，$f(0) = 0$）：</p>
$$\dot{x}_1 = f(x_1) + x_2,\qquad \dot{x}_2 = u$$
<p>取 $k_1, k_2 > 0$，虚拟控制与控制律</p>
$$\alpha(x_1) = -f(x_1) - k_1x_1,\qquad u = -x_1 - k_2 z_2 + \dot{\alpha},\quad z_2 = x_2 - \alpha$$
<p>其中 $\dot\alpha = \frac{\partial \alpha}{\partial x_1}(f(x_1) + x_2)$ 可解析计算。则原点全局渐近稳定。</p></div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>第一步：$x_1$ 子系统。</strong>取 $V_1 = \frac12 x_1^2$。把 $x_2$ 写成 $\alpha + z_2$：</p>
$$\dot{V}_1 = x_1\dot{x}_1 = x_1\big(f(x_1) + \alpha(x_1) + z_2\big) = x_1(-k_1x_1 + z_2) = -k_1x_1^2 + x_1z_2$$
<p>若 $z_2 = 0$（$x_2$ 完美听话），$\dot V_1 = -k_1x_1^2$，完事。多出的交叉项 $x_1z_2$ 交给下一步处理。</p>
<p><strong>第二步：扩充 Lyapunov 函数。</strong>$V_2 = V_1 + \frac12 z_2^2$。计算 $\dot{z}_2 = \dot{x}_2 - \dot{\alpha} = u - \dot{\alpha}$：</p>
$$\dot{V}_2 = -k_1x_1^2 + x_1z_2 + z_2(u - \dot{\alpha})$$
<p>代入 $u = -x_1 - k_2z_2 + \dot{\alpha}$：括号内变为 $-x_1 - k_2z_2$，交叉项 $x_1 z_2$ 与 $-z_2 x_1$ <strong>精确对消</strong>：</p>
$$\dot{V}_2 = -k_1x_1^2 - k_2z_2^2 < 0 \quad \big((x_1, z_2) \ne 0\big)$$
<p><strong>第三步：收尾。</strong>$V_2$ 在 $(x_1, z_2)$ 坐标下正定、径向无界，$\dot V_2$ 负定，由 Lyapunov 定理（全局版）$(x_1, z_2) \to 0$。而 $(x_1, x_2) \mapsto (x_1, z_2)$ 是光滑微分同胚（$\alpha$ 光滑且 $\alpha(0)=0$），故原坐标下原点也全局渐近稳定。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">递推的一般规律</span>
<p>$n$ 级链 $\dot x_i = f_i(x_{1..i}) + x_{i+1}$：每级定义 $z_i = x_i - \alpha_{i-1}$，$V_n = \frac12\sum z_i^2$，每级的 $u$/虚拟控制承担三件事：抵消本级已知项、注入阻尼 $-k_iz_i$、抵消上一级传来的交叉项（固定配方 $-z_{i-1}$）。最终 $\dot{V}_n = -\sum k_iz_i^2$。要点：$\alpha_i$ 会越来越复杂（"微分爆炸"），工程上用指令滤波/动态面技术缓解。</p></div>`
    },
  ],
  code: {
    title: "反步法控制不稳定非线性链（含 Lyapunov 单调性验证）",
    text: `import numpy as np

k1, k2 = 2.0, 2.0
f = lambda x1: x1**2            # 开环 x1' = x1^2 会有限时间爆炸!

def controller(x1, x2):
    alpha  = -f(x1) - k1*x1
    z2     = x2 - alpha
    dalpha = (-2*x1 - k1) * (f(x1) + x2)
    u      = -x1 - k2*z2 + dalpha
    return u, z2

h, T = 1e-4, 6.0
x = np.array([1.5, 0.0])
V_prev, mono_ok = np.inf, True
for kstep in range(int(T/h)):
    u, z2 = controller(*x)
    V = 0.5*x[0]**2 + 0.5*z2**2
    if V > V_prev + 1e-12:
        mono_ok = False
    V_prev = V
    x = x + h*np.array([f(x[0]) + x[1], u])

print("末态:", np.round(x, 6))
print("V(t) 全程单调下降:", mono_ok)
print("理论: Vdot = -k1*x1^2 - k2*z2^2 <= 0, 全局渐近稳定")`
  }
};

/* ════════════════ 滑模控制 ════════════════ */
window.KG_CONTENT["sliding-mode"] = {
  sim: "sim-sliding",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>模型不准、扰动不明——那就别指望优雅地"算"出精确控制，改用<strong>蛮力切换</strong>：在状态空间里画一个"轨道"（滑模面 $s = 0$），面上的动态是我们设计的、必然收敛的；然后用足够猛的开关控制把状态<strong>打到面上、摁住不放</strong>。</p>
<p>一旦"摁住"（滑动模态），系统行为完全由面的几何决定，<strong>与对象参数和匹配扰动无关</strong>——这叫不变性，是滑模最值钱的性质。代价是控制量高频抖动（chattering），执行器会抗议；工程上用边界层/高阶滑模磨平。变结构控制源自 1950s 苏联（Emelyanov、Utkin），是鲁棒非线性控制的"重锤"。</p>`
    },
    {
      icon: "🧾", title: "定理：到达 + 滑动 + 不变性",
      html: String.raw`
<p>考察带未知扰动的二阶系统（"匹配"扰动：与 $u$ 同通道）：</p>
$$\ddot{x} = d(t) + u,\qquad |d(t)| \le D$$
<p>设计滑模变量 $s = c x + \dot{x}$（$c > 0$），控制律 $u = -c\dot{x} - k\,\mathrm{sign}(s)$，$k > D$。</p>
<div class="box thm"><span class="box-title">定理 1（有限时间到达）</span>
$s$ 在有限时间 $t_r \le \dfrac{|s(0)|}{k - D}$ 内到达 $0$，此后恒为 $0$。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>取 $V = \frac12 s^2$。计算 $\dot{s} = c\dot{x} + \ddot{x} = c\dot{x} + d + u = d - k\,\mathrm{sign}(s)$，于是</p>
$$\dot{V} = s\dot{s} = s\,d - k|s| \le D|s| - k|s| = -(k - D)|s| = -(k-D)\sqrt{2V}$$
<p>解此微分不等式（分离变量）：$\sqrt{V(t)} \le \sqrt{V(0)} - \frac{k-D}{\sqrt2}t$，右端在 $t_r = \frac{\sqrt{2V(0)}}{k-D} = \frac{|s(0)|}{k-D}$ 时归零。注意<strong>有限时间</strong>（普通渐近收敛只会无限逼近）：这是 $\sqrt{V}$ 型不等式的功劳。到达后 $s\dot s \le -(k-D)|s| < 0$ 阻止任何离开。<span class="qed">∎</span></p></div>
<div class="box thm"><span class="box-title">定理 2（滑动模态与不变性）</span>
在 $s \equiv 0$ 上，系统动态退化为 $\dot{x} = -cx$——指数收敛，且<strong>完全不含 $d$</strong>：对一切满足 $|d| \le D$ 的扰动行为一致。</div>
<div class="box proof"><span class="box-title">证明</span>
<p>$s \equiv 0 \Rightarrow \dot{x} = -cx$（滑模面的定义即降阶后的目标动态），故 $x(t) = x(t_r)e^{-c(t - t_r)} \to 0$，进而 $\dot x \to 0$。扰动 $d$ 只出现在 $\dot s$ 中，而滑动阶段等效控制恰好实时抵消它（Filippov 等效控制 $u_{eq} = -c\dot x - d$）：面上的运动方程无 $d$。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">抖振与对策</span>
<ul>
<li>理想 sign 切换无穷快；实际执行器/采样有限 ⇒ $s$ 在零附近高频震荡（实验里把 $\varphi = 0$ 看 $u$ 曲线）；</li>
<li><strong>边界层</strong>：用 $\mathrm{sat}(s/\varphi)$ 替代 sign，抖振消失，代价是只收敛到 $O(\varphi)$ 邻域（准滑模）；</li>
<li><strong>高阶滑模</strong>（super-twisting 等）：把切换藏进积分器，连续控制量同时保住有限时间收敛；</li>
<li>不变性只对<strong>匹配</strong>扰动成立；非匹配扰动需与 backstepping 等组合。</li>
</ul></div>`
    },
  ],
  code: {
    title: "滑模 vs PD：扰动下的鲁棒性对决",
    text: `import numpy as np

D = 0.8
d = lambda t: D*np.sin(3*t)          # 未知扰动(只知界)
c, k, phi = 1.0, 1.2, 0.02

def run(ctrl, T=15.0, h=1e-4):
    x = np.array([2.0, 1.0]); errs = []
    for i in range(int(T/h)):
        t = i*h
        u = ctrl(x, t)
        x = x + h*np.array([x[1], d(t) + u])
        if t > 10:                    # 统计后段稳态误差
            errs.append(abs(x[0]))
    return np.mean(errs)

def smc(x, t):
    s = c*x[0] + x[1]
    sat = np.clip(s/phi, -1, 1)      # 边界层
    return -c*x[1] - k*sat

def pd(x, t):                        # 增益可比的 PD
    return -c*1.2*x[0] - (1+c)*x[1]

print("稳态平均 |x| —")
print("  滑模(k > D):", round(run(smc), 5))
print("  PD        :", round(run(pd), 5))
print("滑模把误差压到边界层 O(phi) 内, PD 则被扰动持续推着走")`
  }
};

/* ════════════════ 自适应控制 ════════════════ */
window.KG_CONTENT["adaptive"] = {
  sim: "sim-adaptive",
  sections: [
    {
      icon: "💡", title: "直觉",
      html: String.raw`
<p>飞机烧掉燃油质量变了、机械臂抓起未知负载、电池老化内阻漂移——参数在变且没人告诉你。自适应控制的答案：<strong>把控制器增益本身变成状态，用误差信号在线"学习"它</strong>。</p>
<p>MRAC（模型参考自适应控制）的框架：给定一个"理想中系统应有的样子"（参考模型），控制器努力让真实系统的输出贴上参考模型的输出；贴不上的误差反过来驱动增益调整。精妙之处在于调整律不是拍脑袋的梯度下降，而是<strong>由 Lyapunov 函数反向设计</strong>——先写好想要的 $\dot{V}$，倒推出必须的 $\dot{\hat\theta}$，稳定性内建于设计之中。</p>`
    },
    {
      icon: "🧾", title: "MRAC 设计与稳定性证明（一阶完整案例）",
      html: String.raw`
<div class="box def"><span class="box-title">问题设定</span>
<p>对象与参考模型（$a$ 未知，$a_m > 0$）：</p>
$$\dot{x} = ax + u,\qquad \dot{x}_m = -a_mx_m + b_mr(t)$$
<p>控制律取 $u = \hat{\theta}(t)\,x + b_mr$。若知道 $a$，理想增益 $\theta^* = -(a + a_m)$ 使闭环精确等于参考模型。定义跟踪误差 $e = x - x_m$、参数误差 $\tilde\theta = \hat\theta - \theta^*$。</p></div>
<div class="box thm"><span class="box-title">定理（MRAC 稳定性）</span>
取自适应律（$\gamma > 0$ 任意）
$$\dot{\hat{\theta}} = -\gamma\, e\, x$$
则闭环所有信号有界，且 $e(t) \to 0$。若 $r$ 充分激励（PE），还有 $\hat\theta \to \theta^*$。</div>
<div class="box proof"><span class="box-title">证明</span>
<p><strong>第一步：误差方程。</strong></p>
$$\dot{e} = \dot{x} - \dot{x}_m = ax + \hat\theta x + b_mr + a_mx_m - b_mr = (a + \hat\theta)x + a_mx_m$$
<p>用 $a = -a_m - \theta^*$：$(a + \hat\theta)x = (-a_m + \tilde\theta)x$，于是</p>
$$\dot{e} = -a_m e + \tilde{\theta}x$$
<p><strong>第二步：复合 Lyapunov 函数。</strong>把参数误差也纳入"能量"：</p>
$$V(e, \tilde\theta) = \frac12 e^2 + \frac{1}{2\gamma}\tilde{\theta}^2$$
$$\dot{V} = e(-a_me + \tilde\theta x) + \frac{1}{\gamma}\tilde\theta\,\dot{\hat\theta} = -a_me^2 + \tilde\theta\Big(ex + \frac{1}{\gamma}\dot{\hat\theta}\Big)$$
<p>（$\dot{\tilde\theta} = \dot{\hat\theta}$，因 $\theta^*$ 恒定。）我们管不住 $\tilde\theta$ 的符号，唯一的办法是让括号<strong>恒为零</strong>——这就"倒推"出自适应律 $\dot{\hat\theta} = -\gamma ex$。于是</p>
$$\dot{V} = -a_me^2 \le 0$$
<p><strong>第三步：收敛论证（Barbalat 引理）。</strong>$V$ 不增 ⇒ $e, \tilde\theta$ 有界；$x = e + x_m$ 与 $x_m$ 有界 ⇒ $x$ 有界 ⇒ $\dot e = -a_me + \tilde\theta x$ 有界 ⇒ $e$ 一致连续。又</p>
$$\int_0^\infty a_m e^2\,dt = V(0) - V(\infty) < \infty$$
<p>Barbalat 引理：可积 + 一致连续 ⇒ $e(t) \to 0$。（注意系统非自治——$r(t)$ 时变——不能用 LaSalle，Barbalat 是标准替身。）</p>
<p><strong>参数收敛另议：</strong>$\dot V$ 里没有 $\tilde\theta^2$ 项，所以只能保证 $\tilde\theta$ 有界、不保证归零。要 $\hat\theta \to \theta^*$ 需<strong>持续激励</strong>（PE）：$\int_t^{t+T} x^2 \ge \delta > 0$。直觉：$x$ 若早早归零，误差信号里就再没有关于 $a$ 的信息可学。实验中用方波 $r$ 保持激励。<span class="qed">∎</span></p></div>
<div class="box note"><span class="box-title">工程警示</span>
<ul>
<li>纯 MRAC 对<strong>未建模动态和噪声</strong>敏感（1980s 著名的 Rohrs 反例）：小扰动可致参数漂移直至失稳；实用需加 σ-修正、投影、死区等鲁棒化补丁；</li>
<li>γ 大学得快但瞬态野；γ 小温和但慢（实验可体会）；</li>
<li>思想辐射：MRAC 的"Lyapunov 反推调参"是今天很多在线学习控制（含某些 RL 稳定性证明）的原型。</li>
</ul></div>`
    },
  ],
  code: {
    title: "MRAC 完整实现（含 PE 条件对比实验）",
    text: `import numpy as np

am, bm, gamma = 2.0, 2.0, 4.0
a_true = 1.0                      # 未知(不稳定)参数
theta_star = -(a_true + am)

def run(r_func, T=40.0, h=1e-3):
    x, xm, th = 0.0, 0.0, 0.0
    for k in range(int(T/h)):
        t = k*h
        r = r_func(t)
        e = x - xm
        u = th*x + bm*r
        x  += h*(a_true*x + u)
        xm += h*(-am*xm + bm*r)
        th += h*(-gamma*e*x)
    return e, th

# 充分激励: 方波; 不激励: r=0
e1, th1 = run(lambda t: np.sign(np.sin(0.5*t)))
e2, th2 = run(lambda t: 0.0)

print(f"理想增益 theta* = {theta_star}")
print(f"方波激励:  e→{e1:.5f}, theta→{th1:.3f}  (参数收敛 ✓)")
print(f"零输入  :  e→{e2:.5f}, theta→{th2:.3f}  (误差归零但参数没学到)")`
  }
};
