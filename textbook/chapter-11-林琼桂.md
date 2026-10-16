# 第十一章 球函数\*

## 目录

§1 轴对称球函数　222  
　§1.1 轴对称问题的一般解　222  
　§1.2 Legendre 多项式的基本性质　223  
　§1.3 Legendre 多项式的微分表示　224  
　§1.4 *Legendre 多项式的积分表示　225  
　§1.5 Legendre 多项式的正交关系　226  
　§1.6 Legendre 多项式的模　226  
　§1.7 广义 Fourier 级数　228  
　§1.8 Legendre 多项式的母函数　229  
　§1.9 Legendre 多项式的递推关系　231  
　§1.10 应用　232  
§2 一般球函数　235  
　§2.1 非轴对称问题、连带 Legendre 函数　235  
　§2.2 连带 Legendre 函数的微分表示　237  
　§2.3 函数 $P_l^{-m}(x)$　237  
　§2.4 连带 Legendre 函数的正交关系　238  
　§2.5 连带 Legendre 函数的模　239  
　§2.6 广义 Fourier 级数　239  
　§2.7 连带 Legendre 函数的递推关系　240  
　§2.8 球谐函数　240  
　§2.9 *球谐函数的加法公式　241  
补充习题　244

## §1 轴对称球函数

### §1.1 轴对称问题的一般解

在球坐标下对 Laplace 方程

$$
\nabla^2u(\boldsymbol r)=0 \tag{1}
$$

分离变量，即寻找形如

$$
u(\boldsymbol r)=R(r)Y(\theta,\phi)=R(r)H(\theta)\Phi(\phi) \tag{2}
$$

的解。考虑到对 $\phi$ 的周期性边界条件，可得

$$
\Phi(\phi)=\{e^{im\phi},e^{-im\phi}\},\qquad m\in\mathbb N. \tag{3}
$$

令 $\cos\theta=x$，$H(\theta)=P(x)$，考虑到 $\theta=0,\pi$ 处的自然边界条件，$P(x)$ 应满足本征值问题

$$
\frac{d}{dx}\left[(1-x^2)\frac{dP}{dx}\right]+\left(\lambda-\frac{m^2}{1-x^2}\right)P=0, \tag{4a}
$$

$$
P(\pm1)=0\quad(m\ne0)\qquad\text{或}\qquad |P(\pm1)|<\infty\quad(m=0). \tag{4b}
$$

对于轴对称问题，取对称轴为球坐标的极轴，则问题的解与 $\phi$ 无关，亦即只需考虑 $m=0$，这时上述本征值问题简化为

$$
\frac{d}{dx}\left[(1-x^2)\frac{dP}{dx}\right]+\lambda P=0, \tag{5a}
$$

$$
|P(\pm1)|<\infty. \tag{5b}
$$

我们已经求得该本征值问题的本征值和本征函数为

$$
\lambda=l(l+1),\qquad P(x)=P_l(x),\qquad l\in\mathbb N. \tag{6}
$$

这里 $P_l(x)$ 是 $l$ 次 Legendre 多项式。将本征值代回径向方程

$$
r^2R''(r)+2rR'(r)-\lambda R(r)=0, \tag{7}
$$

容易解出

$$
R(r)=\{r^l,1/r^{l+1}\}. \tag{8}
$$

所以，轴对称问题的一般解为

$$
u(r,\theta)=\sum_{l=0}^{\infty}\left(A_lr^l+\frac{B_l}{r^{l+1}}\right)P_l(\cos\theta), \tag{9}
$$

其中 $A_l$ 和 $B_l$ 是任意常数，由边界条件决定。注意 Laplace 方程只有边界条件，没有初始条件，所以边界条件不需要全部是齐次的，否则就只有平庸解了。

应该指出，式 (9) 的第一部分求和，即 $\sum_{l=0}^{\infty}A_lr^lP_l(\cos\theta)$ 在原点的邻域上满足 Laplace 方程；其第二部分求和，或者再加上第一部分中的常数项，即 $A_0+\sum_{l=0}^{\infty}(B_l/r^{l+1})P_l(\cos\theta)$ 在无穷远点的邻域上（以原点为球心的某球面外部）满足 Laplace 方程。一般解 (9) 在原点的去心邻域上满足 Laplace 方程。在全空间满足 Laplace 方程的解则只能是 $u(r,\theta)=A_0$，即为常数。

本节题目中的“轴对称球函数”即指一般解 (9) 中出现的 Legendre 多项式。为了将该解应用于物理问题，我们需要首先研究 Legendre 多项式的有关性质。

**注**　如果读者无法掌握导出一般解 (9) 的全过程，至少也应该对此过程有所了解。无论如何，读者应该通过反复应用熟练掌握一般解。在其它课程里，比如电动力学，将会直接引用这一结果。

### §1.2 Legendre 多项式的基本性质

我们已经在上一章求得 Legendre 多项式的显明表达式 (explicit expression，简称显式)

$$
P_l(x)=\sum_{k=0}^{[l/2]}(-)^k\frac{(2l-2k)!}{2^lk!(l-k)!(l-2k)!}x^{l-2k},\qquad l\in\mathbb N. \tag{10}
$$

它具有下列基本性质。

1. 奇偶性：

$$
P_l(-x)=(-)^lP_l(x). \tag{11}
$$

由显式立即可以得出这一结论。

2. 原点值：

$$
P_{2n+1}(0)=0,\qquad P_{2n}(0)=(-)^n\frac{(2n)!}{2^{2n}(n!)^2},\qquad n\in\mathbb N. \tag{12}
$$

第一式是明显的，因为奇次的 Legendre 多项式中只包含 $x$ 的奇次幂。至于第二式，它应该等于 $P_{2n}(x)$ 中的常数项，即式 (10) 中 $k=n$ 的那一项（注意此时 $l=2n$）。

3. 端点值：

$$
P_l(1)=1,\qquad P_l(-1)=(-)^l. \tag{13}
$$

其中第一式将在下一小节予以证明，第二式则容易由第一式和奇偶性得出。

4. 头三个 Legendre 多项式的具体形式为

$$
P_0(x)=1,\qquad P_1(x)=x,\qquad P_2(x)=\frac12(3x^2-1). \tag{14}
$$

将 $l$ 值代入显式 (10) 或以下给出的微分表示，即可得出各阶 Legendre 多项式的具体形式。

5. 零点：$P_l(x)$ 在区间 $(-1,1)$ 上有 $l$ 个一阶零点。由于 $P_l(x)$ 是 $l$ 次多项式，它总共只有 $l$ 个零点，所以也可以说，$P_l(x)$ 的 $l$ 个零点均为实数、全部分布在区间 $(-1,1)$ 上且没有重零点。对于低阶的 Legendre 多项式，容易验证这一结论。一般情况可参考下一小节的证明。

### §1.3 Legendre 多项式的微分表示

Legendre 多项式有以下微分表示 (differential representation)

$$
P_l(x)=\frac{1}{2^ll!}\frac{d^l}{dx^l}(x^2-1)^l. \tag{15}
$$

这称为 Rodrigues 公式。

**证明**　由二项式定理，

$$
(x^2-1)^l=\sum_{k=0}^{l}C_l^kx^{2l-2k}(-)^k=\sum_{k=0}^{l}(-)^kC_l^kx^{2l-2k},
$$

将上式对 $x$ 求导 $l$ 次，容易看出，只有满足 $2l-2k\ge l$ 的项求导后才不为零，换句话说，求导后只需对 $k\le l/2$ 求和，利用求导公式 $(d^l/dx^l)x^n=[n!/(n-l)!]x^{n-l}$，可得

$$
\frac{d^l}{dx^l}(x^2-1)^l
=\sum_{k=0}^{[l/2]}(-)^kC_l^k\frac{d^l}{dx^l}x^{2l-2k}
=\sum_{k=0}^{[l/2]}(-)^kC_l^k\frac{(2l-2k)!}{(l-2k)!}x^{l-2k}
$$

$$
=\sum_{k=0}^{[l/2]}(-)^k\frac{l!(2l-2k)!}{k!(l-k)!(l-2k)!}x^{l-2k},
$$

所以

$$
\frac{1}{2^ll!}\frac{d^l}{dx^l}(x^2-1)^l
=\sum_{k=0}^{[l/2]}(-)^k\frac{(2l-2k)!}{2^lk!(l-k)!(l-2k)!}x^{l-2k},
$$

与式 (10) 比较，可见这正是 $P_l(x)$。证毕。

利用微分表示，容易证明 $P_l(1)=1$。事实上，由高阶导数的 Leibniz 公式，有

$$
P_l(x)=\frac{1}{2^ll!}\frac{d^l}{dx^l}[(x-1)^l(x+1)^l]
=\frac{1}{2^ll!}\sum_{k=0}^{l}C_l^k\left[\frac{d^k}{dx^k}(x-1)^l\right]\left[\frac{d^{l-k}}{dx^{l-k}}(x+1)^l\right],
$$

显然，$k<l$ 的各项求导后含有因子 $(x-1)^{l-k}$，因而在 $x=1$ 处为零，所以

$$
P_l(1)=\left.\frac{1}{2^ll!}\left[\frac{d^l}{dx^l}(x-1)^l\right](x+1)^l\right|_{x=1}
=\frac{1}{2^ll!}l!2^l=1.
$$

利用微分表示，还可以证明 $P_l(x)$ 在区间 $(-1,1)$ 上有 $l$ 个零点。首先，当 $n\le l$，我们有

$$
\frac{d^n}{dx^n}(x^2-1)^l
=\frac{d^n}{dx^n}[(x-1)^l(x+1)^l]
=\sum_{k=0}^{n}C_n^k\left[\frac{d^k}{dx^k}(x-1)^l\right]\left[\frac{d^{n-k}}{dx^{n-k}}(x+1)^l\right]
$$

$$
=\sum_{k=0}^{n}C_n^k\frac{l!}{(l-k)!}\frac{l!}{(l-n+k)!}(x-1)^{l-k}(x+1)^{l-n+k},
$$

容易看出，若 $n<l$，则上式在 $x=\pm1$ 处均为零。令 $(x^2-1)^l$ 在 $x=\pm1$ 处均为零（设 $l\ge1$），由微积分中的 Rolle 定理，$(d/dx)(x^2-1)^l$ 在区间 $(-1,1)$ 上至少有一个零点。由上面的讨论，$(d/dx)(x^2-1)^l$ 在 $x=\pm1$ 处仍为零（设 $l\ge2$）。所以 $(d^2/dx^2)(x^2-1)^l$ 在区间 $(-1,1)$ 上至少有两个零点。继续同样的推理可知，$(d^l/dx^l)(x^2-1)^l$ 或 $P_l(x)$ 在区间 $(-1,1)$ 上至少有 $l$ 个零点。然而，$P_l(x)$ 是 $l$ 次多项式，它总共只有 $l$ 个零点，所以，它在区间 $(-1,1)$ 上正好有 $l$ 个零点。由于 $P_l(x)$ 是二阶齐次常微分方程的解，所以这些零点都是一阶的，因为若存在二阶以上的零点 $x_0$，则 $P_l(x_0)=P_l'(x_0)=0$，由此容易导出 $P_l(x)\equiv0$，这显然是错误的。

最后，我们用微分表示推出 Legendre 多项式的其它显式。首先，将 $(x^2-1)^l$ 写成以下形式

$$
(x^2-1)^l=[(x-1)(x-1+2)]^l=(x-1)^l\sum_{k=0}^{l}C_l^k2^{l-k}(x-1)^k
=\sum_{k=0}^{l}2^{l-k}C_l^k(x-1)^{l+k},
$$

然后，对 $x$ 求导 $l$ 次，并除以 $2^ll!$，容易推出

$$
P_l(x)=\sum_{k=0}^{l}\frac{(l+k)!}{(k!)^2(l-k)!}\left(\frac{x-1}{2}\right)^k. \tag{16}
$$

由这一表达式，马上可以得到 $P_l(1)=1$，因为 $x=1$ 时，只有常数项（对应于 $k=0$）有贡献。

另一个类似的表达式是

$$
P_l(x)=(-)^l\sum_{k=0}^{l}\frac{(l+k)!}{(k!)^2(l-k)!}\left(-\frac{x+1}{2}\right)^k. \tag{17}
$$

请有兴趣的读者自己导出上式。由上式立得 $P_l(-1)=(-)^l$。

### §1.4 *Legendre 多项式的积分表示

由微分表示和 Cauchy 高阶导数公式立得

$$
P_l(x)=\frac{1}{2^l2\pi i}\int_C\frac{(\zeta^2-1)^l}{(\zeta-x)^{l+1}}\,d\zeta, \tag{18}
$$

其中 $C$ 是包围 $x$ 的任一围线。上式称为 Schläfli 积分。

将围线 $C$ 取为以 $x$ 为圆心，以 $\sqrt{x^2-1}$ 为半径的圆周，则在 $C$ 上有 $\zeta=x+\sqrt{x^2-1}e^{i\phi}$，$d\zeta=i\sqrt{x^2-1}e^{i\phi}d\phi$，代入上述 Schläfli 积分，化简后得到

$$
P_l(x)=\frac{1}{\pi}\int_0^\pi\left(x+\sqrt{x^2-1}\cos\phi\right)^l d\phi. \tag{19}
$$

上式中 $\sqrt{x^2-1}$ 可取任一分支，因为若取另一分支，则括号中第二项改变符号，这时只需作积分变量置换 $\phi'=\pi-\phi$，然后将积分变量 $\phi'$ 重新写成 $\phi$，则所得结果与上式完全一致。上式称为 Laplace 积分。

应该指出，虽然物理上关心的变量范围是 $-1\le x\le1$（因为 $x=\cos\theta$），但 Legendre 多项式对于自变量取一切实数甚至复数都是有定义的，所以无论是微分表示还是积分表示 (integral representation)，都是对自变量的一切取值有效的。

由 Laplace 积分，立得 $P_l(1)=1$，$P_l(-1)=(-)^l$。另外，取 $x=\cos\theta$，则有

$$
P_l(\cos\theta)=\frac{1}{\pi}\int_0^\pi(\cos\theta+i\sin\theta\cos\phi)^l d\phi. \tag{20}
$$

由于 $|\cos\theta+i\sin\theta\cos\phi|=[\cos^2\theta+\sin^2\theta\cos^2\phi]^{1/2}\le1$，由上式和积分的基本性质立得 $|P_l(\cos\theta)|\le1$，而用其它形式则难以证明这一结论。

### §1.5 Legendre 多项式的正交关系

作为 Sturm–Liouville 本征值问题的特例，Legendre 多项式在区间 $[-1,1]$ 上相互正交，注意到权函数 $\rho(x)=1$，正交关系为

$$
\int_{-1}^{1}P_l(x)P_{l'}(x)\,dx=0,\qquad l\ne l'. \tag{21}
$$

以前我们曾经讲过，如果本征值问题的本征函数比较复杂，则正交性难以由直接计算加以验证。上述 Legendre 多项式的正交关系直接验证起来就比以前验证三角函数的正交关系要来得困难。不过，Legendre 多项式不算太复杂，其正交关系还是可以通过直接计算来验证。在这一过程中，我们可以学习怎样计算一些涉及 Legendre 多项式的积分。

首先，我们计算积分 $\int_{-1}^{1}x^kP_l(x)dx$，其中 $k\le l$。利用 Legendre 多项式的微分表示，有

$$
\begin{aligned}
\int_{-1}^{1}x^kP_l(x)dx
&=\frac{1}{2^ll!}\int_{-1}^{1}x^k\frac{d^l}{dx^l}(x^2-1)^l dx\\
&=\frac{1}{2^ll!}\int_{-1}^{1}x^k d\left[\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l\right]\\
&=\left.\frac{1}{2^ll!}x^k\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l\right|_{-1}^{1}
-\frac{k}{2^ll!}\int_{-1}^{1}x^{k-1}\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l dx\\
&=-\frac{k}{2^ll!}\int_{-1}^{1}x^{k-1}\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l dx,
\end{aligned}
$$

其中作了分部积分，积出部分为零，理由与前面论证 $P_l(1)=1$ 时所用类似，亦可参看下一小节。继续作分部积分，共作 $k$ 次，可以得到

$$
\int_{-1}^{1}x^kP_l(x)dx=\frac{(-)^kk!}{2^ll!}\int_{-1}^{1}\frac{d^{l-k}}{dx^{l-k}}(x^2-1)^l dx. \tag{22}
$$

如果 $l>k$，则 $l-k\ge1$，上面的被积函数仍然是一个导数，其原函数是 $(d^{l-k-1}/dx^{l-k-1})(x^2-1)^l$，它在 $x=\pm1$ 处显然为零，于是，我们得到

$$
\int_{-1}^{1}x^kP_l(x)dx=0,\qquad k<l. \tag{23}
$$

注意式 (22) 和 (23) 对 $k=0$ 也成立，而式 (22) 对 $k=l$ 也成立。

有了式 (23)，就容易得到正交关系 (21)。不失一般性，设 $l'<l$，由于 $P_{l'}(x)$ 是 $l'$ 次多项式，它可以写成 $P_{l'}(x)=\sum_{k=0}^{l'}a_kx^k$，其中每一项的幂都满足 $k\le l'<l$，故每一项都与 $P_l(x)$ 正交，从而 $P_{l'}(x)$ 与 $P_l(x)$ 正交。

如果 $k>l$，则式 (23) 中的积分可能不为零，下一小节我们还会再讨论这一问题。不过，如果 $k$、$l$ 具有不同的奇偶性，该积分显然为零，而与 $k$、$l$ 的相对大小无关。

### §1.6 Legendre 多项式的模

Legendre 多项式的模定义为

$$
\|P_l\|=\left\{\int_{-1}^{1}[P_l(x)]^2dx\right\}^{1/2}. \tag{24}
$$

对于任何本征值问题的本征函数族，都需要计算本征函数的模，因为作广义 Fourier 展开时需要用到它。Legendre 多项式也不例外。

**注**　下一节将用 $P_l^m(x)$ 表示连带 Legendre 函数，所以我们应该注意不要把 $[P_l(x)]^2$ 写成 $P_l^2(x)$，否则会与 $m=2$ 时的 $P_l^m(x)$ 相混淆。

利用 Legendre 多项式的微分表示，有

$$
\begin{aligned}
\|P_l\|^2
&=\int_{-1}^{1}[P_l(x)]^2dx
=\frac{1}{2^ll!}\int_{-1}^{1}P_l(x)\frac{d^l}{dx^l}(x^2-1)^l dx\\
&=\frac{1}{2^ll!}\int_{-1}^{1}P_l(x)d\left[\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l\right]\\
&=\left.\frac{1}{2^ll!}P_l(x)\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l\right|_{-1}^{1}
-\frac{1}{2^ll!}\int_{-1}^{1}P_l'(x)\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l dx\\
&=-\frac{1}{2^ll!}\int_{-1}^{1}P_l'(x)\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l dx,
\end{aligned}
$$

其中作了分部积分，而积出部分为零，这是由于

$$
\frac{d^{l-1}}{dx^{l-1}}(x^2-1)^l
=\sum_{k=0}^{l-1}C_{l-1}^k\frac{l!}{(l-k)!}\frac{l!}{(k+1)!}(x-1)^{l-k}(x+1)^{k+1},
$$

注意到 $k$ 的取值范围，易见上式在 $x=\pm1$ 处均为零。继续作分部积分，共作 $l$ 次，可以得到

$$
\|P_l\|^2=\frac{(-)^l}{2^ll!}\int_{-1}^{1}\frac{d^lP_l(x)}{dx^l}(x^2-1)^l dx. \tag{25}
$$

由于 $P_l(x)$ 是 $l$ 次多项式，求导 $l$ 次后为常数，事实上，

$$
\frac{d^lP_l(x)}{dx^l}=\frac{1}{2^ll!}\frac{d^{2l}}{dx^{2l}}(x^2-1)^l
=\frac{1}{2^ll!}\frac{d^{2l}}{dx^{2l}}x^{2l}=\frac{(2l)!}{2^ll!},
$$

其中幂次小于 $2l$ 的各项求导后为零，代入式 (25) 得到

$$
\|P_l\|^2=\frac{(-)^l(2l)!}{2^{2l}(l!)^2}\int_{-1}^{1}(x^2-1)^l dx
=\frac{2(2l)!}{2^{2l}(l!)^2}\int_0^1(1-x^2)^l dx. \tag{26}
$$

最后一个积分可以通过变量置换 $x=\cos\theta$ 化成三角函数积分来计算：

$$
\int_0^1(1-x^2)^l dx=\int_0^{\pi/2}\sin^{2l+1}\theta\,d\theta
=\frac{2^{2l}(l!)^2}{(2l+1)!}. \tag{27}
$$

最后的结果是读者在微积分中做过的，可以通过作分部积分建立递推关系 $I_{2l+1}=[2l/(2l+1)]I_{2l-1}$ 并利用 $I_1=1$ 而得到，其中 $I_{2l+1}$ 即表示上述积分。计算上式的另一种方法是利用 $(1-x^2)^l=(1-x)^l(1+x)^l$，经过 $l$ 次分部积分，最后化成 $(1+x)^{2l}$ 的积分，即可得出同样的结果。将这一结果代入式 (26) 立得

$$
\|P_l\|=\sqrt{\frac{2}{2l+1}}. \tag{28}
$$

最后，我们可以将式 (21) 和 (28) 统一写成

$$
\int_{-1}^{1}P_l(x)P_{l'}(x)dx=\frac{2}{2l+1}\delta_{ll'}. \tag{29}
$$

我们也可以利用上节的结果 (22) 和 (23) 来计算 $\|P_l\|$。由于 $P_l(x)=\sum_{k=0}^{l}a_kx^k$，利用式 (23)，就有

$$
\|P_l\|^2=\int_{-1}^{1}[P_l(x)]^2dx
=\int_{-1}^{1}P_l(x)\sum_{k=0}^{l}a_kx^kdx
=a_l\int_{-1}^{1}x^lP_l(x)dx,
$$

再利用式 (22)，就有

$$
\|P_l\|^2=a_l\frac{(-)^l}{2^l}\int_{-1}^{1}(x^2-1)^l dx
=2\frac{a_l}{2^l}\int_0^1(1-x^2)^l dx.
$$

代入 $a_l=(2l)!/2^l(l!)^2$，易见结果与式 (26) 一致。

### §1.7 广义 Fourier 级数

作为 Sturm–Liouville 本征值问题的特例，Legendre 多项式在区间 $[-1,1]$ 上是完备的，所以，区间 $[-1,1]$ 上的任一函数 $f(x)$，如果具有良好的解析性质，必可用 $\{P_l(x)\}_{l=0}^{\infty}$ 展开为广义 Fourier 级数

$$
f(x)=\sum_{l=0}^{\infty}f_lP_l(x),\qquad -1\le x\le1, \tag{30}
$$

利用式 (29)，容易求得展开系数为

$$
f_l=\frac{2l+1}{2}\int_{-1}^{1}f(x)P_l(x)dx,\qquad l=0,1,2,\cdots. \tag{31}
$$

下面看一个实例，$f(x)=x^n$，其中 $n$ 为非负整数。由式 (23)，如果 $l>n$，则 $f_l=0$，这就是说，将 $x^n$ 展开为 Legendre 多项式的广义 Fourier 级数，其中只包含 $l\le n$ 的项，这是一个非常直观的结论。现在需要计算

$$
f_l=\frac{2l+1}{2}\int_{-1}^{1}x^nP_l(x)dx,\qquad l\le n.
$$

将微分表示代入，作分部积分，过程与推导式 (22) 的过程类似，所不同者是现在 $n\ge l$，所以应该作 $l$ 次分部积分，而得

$$
f_l=\frac{2l+1}{2}\frac{1}{2^ll!}\frac{n!}{(n-l)!}\int_{-1}^{1}x^{n-l}(1-x^2)^l dx.
$$

如果 $l$ 与 $n$ 奇偶性不同，则容易看出 $f_l=0$，这也是非常直观的结论。所以，只有当 $l=n-2k$，$k=0,1,\cdots,[n/2]$ 时，系数才不为零。这些系数是

$$
f_{n-2k}=\frac{2(n-2k)+1}{2^{n-2k}(n-2k)!}\frac{n!}{(2k)!}\int_0^1x^{2k}(1-x^2)^{n-2k}dx. \tag{32}
$$

令 $t=x^2$，上式中的积分可以化为

$$
\begin{aligned}
\int_0^1x^{2k}(1-x^2)^{n-2k}dx
&=\frac12\int_0^1t^{k-1/2}(1-t)^{n-2k}dt\\
&=\frac12B\left(k+\frac12,n-2k+1\right)\\
&=\frac{\Gamma(k+1/2)\Gamma(n-2k+1)}{2\Gamma(n-k+3/2)}\\
&=\frac{2^{2n-4k}(2k)!(n-k)!(n-2k)!}{k!(2n-2k+1)!},
\end{aligned}
$$

其中用到公式 $B(p,q)=\Gamma(p)\Gamma(q)/\Gamma(p+q)$，$\Gamma(k+1)=k!$，$\Gamma(k+1/2)=\sqrt\pi(2k)!/2^{2k}k!$ 等，代回式 (32) 得到

$$
f_{n-2k}=\frac{(2n-4k+1)2^{n-2k}n!(n-k)!}{k!(2n-2k+1)!}. \tag{33}
$$

最后得到

$$
x^n=\sum_{k=0}^{[n/2]}\frac{(2n-4k+1)2^{n-2k}n!(n-k)!}{k!(2n-2k+1)!}P_{n-2k}(x). \tag{34}
$$

由此可见，即使是很简单的函数，其展开式也可以是颇为复杂的。

由上式容易得到

$$
1=P_0(x),\qquad x=P_1(x),\qquad x^2=\frac23P_2(x)+\frac13P_0(x),\qquad x^3=\frac25P_3(x)+\frac35P_1(x).
$$

这些都可以直接加以验证。通过特例的计算，可以帮助我们发现展开式中可能存在的错误。

### §1.8 Legendre 多项式的母函数

考虑三维空间的 Poisson 方程

$$
\nabla^2u(\boldsymbol r)=-4\pi\delta(x)\delta(y)\delta(z-1). \tag{35}
$$

如果认为 $u(\boldsymbol r)$ 描述一个静电场的电势，将上式与静电场的 Poisson 方程 $\nabla^2u(\boldsymbol r)=-\rho(\boldsymbol r)/\epsilon_0$（其中 $\epsilon_0$ 是真空介电常数）比较，可见 $\rho(\boldsymbol r)=4\pi\epsilon_0\delta(x)\delta(y)\delta(z-1)$，这是一个点电荷的电荷密度，点电荷位于单位球面的北极，其位置 $\boldsymbol r_0$ 用直角坐标表示为 $(0,0,1)$，点电荷的电量为 $Q=4\pi\epsilon_0$，根据电磁学，它产生的静电势为 $u(\boldsymbol r)=Q/4\pi\epsilon_0|\boldsymbol r-\boldsymbol r_0|=1/|\boldsymbol r-\boldsymbol r_0|$，用 $\boldsymbol r$ 的球坐标 $(r,\theta,\phi)$ 表示，即

$$
u(\boldsymbol r)=u(r,\theta)=\frac{1}{\sqrt{1-2r\cos\theta+r^2}}. \tag{36}
$$

这就是上述 Poisson 方程的解。这里是借助于电磁学得到的。但可以从数学上直接证明它确实是式 (35) 的解（参看 Green 函数法一章）。当然，读者可能会想到，这不是唯一的解，因为对式 (36) 加上任一调和函数（即 Laplace 方程的解）后也还是解。不过，我们已经知道，在全空间有限的调和函数只有常数，所以式 (35) 的一般解也只是在式 (36) 上加上一个任意常数，如果将电势零点取为无穷远点，则该常数为零，式 (36) 就是式 (35) 的唯一解。

另一方面，考虑单位球区域 $r<1$，在该区域内，式 (35) 成为 Laplace 方程 $\nabla^2u(\boldsymbol r)=0$，由于轴对称性，该区域的解必定可以表示为式 (9)。因为 $r=0$ 处没有电荷，该处的电势应该有限，这就要求所有 $B_l=0$，于是

$$
u(\boldsymbol r)=u(r,\theta)=\sum_{l=0}^{\infty}A_lr^lP_l(\cos\theta),\qquad r<1. \tag{37}
$$

适当选取式 (37) 中的系数 $A_l$，应该可以使得它与式 (36) 一致。或者说，式 (36) 应该可以展开为

$$
\frac{1}{\sqrt{1-2rx+r^2}}=\sum_{l=0}^{\infty}A_lr^lP_l(x),\qquad r<1, \tag{38}
$$

其中已将 $\cos\theta$ 写成 $x$。对比式 (30) 和 (31)，就可以计算展开系数 $A_l$，但是，这样计算是比较困难的。我们注意到上面展开式中还有另一个变量 $r$，所以我们采用一个取巧的办法，令 $x=1$，注意到 $P_l(1)=1$，上式就成为

$$
\frac{1}{1-r}=\sum_{l=0}^{\infty}A_lr^l,\qquad r<1, \tag{39}
$$

由此立得所有 $A_l=1$。这样，我们就得到

$$
\frac{1}{\sqrt{1-2rx+r^2}}=\sum_{l=0}^{\infty}r^lP_l(x),\qquad r<1. \tag{40}
$$

上式左边的函数就称为 Legendre 多项式的母函数或生成函数 (generating function)。如果 $r>1$，我们可以类似地证明

$$
\frac{1}{\sqrt{1-2rx+r^2}}=\sum_{l=0}^{\infty}\frac{1}{r^{l+1}}P_l(x),\qquad r>1. \tag{41}
$$

这也可以由式 (40) 直接得到，而不必重新推导。事实上，当 $r>1$，有 $\sqrt{1-2rx+r^2}=r\sqrt{1-2(1/r)x+(1/r)^2}$，其中 $1/r<1$，如此即可利用式 (40) 而得到式 (41)。

**思考**　如何利用母函数和正交关系计算 Legendre 多项式的模？（提示：将式 (40) 两边平方并对 $x$ 积分。）

上面借助于电磁学找到了 Legendre 多项式的母函数，并导出了它的展开式。下面我们直接将母函数作 Taylor 展开，导出类似的、但适用范围更广的结果。这一方法在数学上也显得更直接。将母函数写作 $\psi(t)=(1-2tx+t^2)^{-1/2}$，其中变量 $t$ 是复数，而参数 $x$ 也可以是复数。容易看出，母函数在复平面上有两个奇点 $x_\pm=x\pm\sqrt{x^2-1}$，而且是支点。用连接两个支点的适当的曲线（见下）作为割线将复平面割破，即可分出单值分支。由于 $x_+x_-=1$，故 $x_\pm\ne0$，所以总存在 $t=0$ 的某邻域，在其中可以分出 $\psi(t)$ 的单值解析分支。我们规定 $\psi(0)=1$，这样就确定了一个单值分支。在 $t=0$ 处，$\psi(t)$ 可以展开为 Taylor 级数

$$
\psi(t)=(1-2tx+t^2)^{-1/2}=\sum_{l=0}^{\infty}a_l(x)t^l, \tag{42}
$$

其中展开系数是参数 $x$ 的函数。当然，为了使得展开式的收敛圆尽可能大，割线不应该穿过圆 $|t|<\min(|x_+|,|x_-|)$。按 Taylor 展开的系数公式，

$$
a_l(x)=\frac{1}{2\pi i}\int_{|t|=\epsilon}\frac{(1-2tx+t^2)^{-1/2}}{t^{l+1}}\,dt, \tag{43}
$$

其中 $\epsilon$ 取为一个很小的数，以便于下面的讨论（原则上只要 $\epsilon<\min(|x_+|,|x_-|)$ 即可）。今作变量置换，令

$$
\zeta=\frac{1-(1-2tx+t^2)^{1/2}}{t}, \tag{44}
$$

容易求出其反变换为 $t=2(\zeta-x)/(\zeta^2-1)$，两边乘以 $\zeta^2-1$ 并微分得 $dt=2(1-\zeta t)d\zeta/(\zeta^2-1)$，其中的因子 $(1-\zeta t)$ 将与 $(1-2tx+t^2)^{-1/2}=(1-\zeta t)^{-1}$ 相消，将这些关系代入式 (43)，整理得

$$
a_l(x)=\frac{1}{2^l}\frac{1}{2\pi i}\int_{C_\epsilon}\frac{(\zeta^2-1)^l}{(\zeta-x)^{l+1}}\,d\zeta. \tag{45}
$$

其中 $C_\epsilon$ 是圆周 $|t|=\epsilon$ 经过映射 (44) 在 $\zeta$ 平面上得到的围线。将式 (44) 在 $t=0$ 附近展开，保留到 $t$ 的一次项，可得 $\zeta-x=(x^2-1)t/2$，由此可见当 $\epsilon$ 很小时，$C_\epsilon$ 是以 $\zeta=x$ 为圆心、半径很小的圆周，而且是逆时针的。根据 Cauchy 高阶导数公式立得

$$
a_l(x)=\frac{1}{2^ll!}\frac{d^l}{dx^l}(x^2-1)^l=P_l(x). \tag{46}
$$

为清晰起见，我们把结果重新写在下面

$$
(1-2tx+t^2)^{-1/2}=\sum_{l=0}^{\infty}P_l(x)t^l,\qquad |t|<|x_\pm|=|x\pm\sqrt{x^2-1}|. \tag{47}
$$

这样得到的结果虽然形式上与式 (40) 相似，但从计算过程可以看出，现在的结果适用范围更广。首先，$x$ 不一定表示 $\cos\theta$，所以 $-1\le x\le1$ 的限制取消了。其次，$x$ 还可以是复数。再次，即使 $x$ 为实数且 $-1\le x\le1$（这时 $|x\pm\sqrt{x^2-1}|=1$ 从而式 (47) 的收敛圆是 $|t|<1$），其适用范围也大于式 (40)，因为 $t$ 可以是复数。

前面推导式 (40) 时，用到了 $P_l(1)=1$。而在上面的推导中，Legendre 多项式是在最后自然出现的，中间没有用到任何有关的性质。所以我们可以从式 (47) 出发导出 Legendre 多项式的一些性质。比如令 $x=1$ 立得 $P_l(1)=1$。又比如将 $x$ 换为 $-x$，同时将 $t$ 换为 $-t$，则左边不变，而右边变成 $\sum_{l=0}^{\infty}(-)^lP_l(-x)t^l$，再与原式比较，立得 $P_l(-x)=(-)^lP_l(x)$。

最后我们看看 $x$ 为实数且 $x>1$ 的情况，令 $x=\cosh\xi$，其中 $\xi>0$。容易得到 $x_\pm=e^{\pm\xi}$，这时式 (47) 的收敛圆是 $|t|<e^{-\xi}$，可见 $x$ 越大，收敛半径就越小。

### §1.9 Legendre 多项式的递推关系

将式 (40) 两边对 $r$ 求导，然后同乘以 $1-2rx+r^2$ 即得

$$
(x-r)\sum_{l=0}^{\infty}r^lP_l(x)=(1-2rx+r^2)\sum_{l=0}^{\infty}lr^{l-1}P_l(x),
$$

比较 $r^l$ 的系数可得

$$
(2l+1)xP_l(x)=(l+1)P_{l+1}(x)+lP_{l-1}(x). \tag{48}
$$

将式 (40) 两边对 $x$ 求导，然后同乘以 $1-2rx+r^2$ 即得

$$
\sum_{l=0}^{\infty}r^{l+1}P_l(x)=(1-2rx+r^2)\sum_{l=0}^{\infty}r^lP_l'(x),
$$

比较 $r^{l+1}$ 的系数可得

$$
P_l(x)=P_{l+1}'(x)-2xP_l'(x)+P_{l-1}'(x). \tag{49}
$$

将式 (48) 求导并与式 (49) 联立消去 $P_l'(x)$ 即得

$$
(2l+1)P_l(x)=P_{l+1}'(x)-P_{l-1}'(x). \tag{50}
$$

以上是几个基本的递推关系，其它的递推关系就不讨论了。

### §1.10 应用

**例 1**　已知半径为 $a$ 的球面上的电势分布为 $f(\theta)$，球内外无电荷，求空间各处的电势。

由于球内外无电荷，故电势在球内外均满足 Laplace 方程，题意已隐含电势零点取在无穷远处，则定解问题为

$$
\nabla^2u=0,\qquad (r<a,r>a) \tag{51a}
$$

$$
u|_{r=a}=f(\theta),\qquad u|_{r=\infty}=0. \tag{51b}
$$

由于定解条件与 $\phi$ 无关，故本问题为轴对称问题，一般解由式 (9) 给出。

首先求解球内的电势，为了计算方便，将一般解写作

$$
u_1(r,\theta)=\sum_{l=0}^{\infty}\left[A_l\left(\frac{r}{a}\right)^l+B_l\left(\frac{a}{r}\right)^{l+1}\right]P_l(\cos\theta), \tag{52}
$$

由于球内没有电荷，故球心即 $r=0$ 处电势应该有限，所以 $B_l=0\ (l\in\mathbb N)$，于是球内的解应为

$$
u_1(r,\theta)=\sum_{l=0}^{\infty}A_l\left(\frac{r}{a}\right)^lP_l(\cos\theta). \tag{53}
$$

代入 $r=a$ 处的边界条件，得到 $\sum_{l=0}^{\infty}A_lP_l(\cos\theta)=f(\theta)$，或 $\sum_{l=0}^{\infty}A_lP_l(x)=g(x)$，其中 $x=\cos\theta$，$g(x)=f(\theta)$。对比式 (30) 和 (31)，可得

$$
A_l=\frac{2l+1}{2}\int_{-1}^{1}g(x)P_l(x)dx
=\frac{2l+1}{2}\int_0^\pi f(\theta)P_l(\cos\theta)\sin\theta\,d\theta. \tag{54}
$$

给定 $f(\theta)$ 的具体形式，即可按上式计算展开系数 $A_l$，将所得系数代回式 (53)，即得球内的电势。

其次求解球外的电势，一般解的形式仍为

$$
u_2(r,\theta)=\sum_{l=0}^{\infty}\left[C_l\left(\frac{r}{a}\right)^l+D_l\left(\frac{a}{r}\right)^{l+1}\right]P_l(\cos\theta), \tag{55}
$$

由于球外没有电荷，故电势应该处处有限，无穷远处亦然。实际上我们已经取无穷远点为电势零点，所以 $C_l=0\ (l=0,1,2,\cdots)$，于是球外的解应为

$$
u_2(r,\theta)=\sum_{l=0}^{\infty}D_l\left(\frac{a}{r}\right)^{l+1}P_l(\cos\theta). \tag{56}
$$

代入 $r=a$ 处的边界条件，得到 $\sum_{l=0}^{\infty}D_lP_l(\cos\theta)=f(\theta)$，所以应有 $D_l=A_l$。将以上算得的系数代回上式，即得球外的电势。

**特例**　$f(\theta)=\cos^2\theta$，或 $g(x)=x^2$。象这样简单的函数，就不需要使用上面计算系数的公式了。注意到 $P_2(x)=(3x^2-1)/2$，立得 $x^2=[2P_2(x)+1]/3$，又因为 $P_0(x)=1$，所以 $x^2=(2/3)P_2(x)+(1/3)P_0(x)$。于是得到系数

$$
A_0=D_0=\frac13,\qquad A_2=D_2=\frac23,\qquad A_l=D_l=0\quad(l\ne0,2). \tag{57}
$$

将这些系数代回式 (53) 和 (56)，就得到

$$
u_1(r,\theta)=\frac13+\frac23\left(\frac{r}{a}\right)^2P_2(\cos\theta), \tag{58a}
$$

$$
u_2(r,\theta)=\frac13\frac{a}{r}+\frac23\left(\frac{a}{r}\right)^3P_2(\cos\theta). \tag{58b}
$$

**注**　对于球外的结果需要特别小心，因为即使系数算对了，写下最后结果时还常常会出错。

**例 2**　在均匀电场 $\boldsymbol E_0$ 中放入半径为 $a$、带电为 $Q$ 的导体球，求放入导体球后的电场分布。

**注**　这是一个电动力学问题，也是一个使用轴对称球函数的典型应用题。对于尚未学过电动力学的读者，要将它翻译成数学问题可能略有困难，但只要学过电磁学，读懂下面的求解过程则是完全没有问题的。本题的目的是让读者了解如何利用本节的数学知识解决物理问题。电动力学中还有不少涉及导体球或介质球的问题与此类似，读者应该牢固掌握本节的内容，以便在后续课程中能熟悉应用。

为了清楚起见，下面将求解过程分为几步。

**第一步　物理分析。**　静电场中的导体达到静电平衡时，体内电场为零，整个导体为等势体，故球内和球面的电势为常数，记为 $u_0$，通过下面的求解可以确定它。关键是求解球外的电势分布。求出电势，电场的分布也就知道了。

**第二步　建立坐标系。**　取球坐标系，以球心为原点，电场 $\boldsymbol E_0$ 的方向为极轴方向。在这样的坐标系下，本问题的电势分布显然具有轴对称性，故球外电势 $u(\boldsymbol r)=u(r,\theta)$。

**第三步　建立定解问题。**　由于球外没有电荷，故电势 $u(r,\theta)$ 满足 Laplace 方程。本题有两个边界，即球面和无穷远处。球面的电势已取为 $u_0$，此外，球面带电 $Q$（球内电场为零，按 Gauss 定律，球内电荷密度 $\rho=\nabla\cdot\boldsymbol E=0$，故电荷全部分布在球面上）。至于无穷远处，原均匀电场应基本不受影响。由此可列出以下定解问题：

$$
\nabla^2u=0,\qquad (r>a) \tag{59a}
$$

$$
u|_{r=a}=u_0, \tag{59b}
$$

$$
\int_{r=a^+}\left(-\frac{\partial u}{\partial r}\right)dS=\frac{Q}{\epsilon_0}, \tag{59c}
$$

$$
u|_{r\to\infty}=-E_0r\cos\theta. \tag{59d}
$$

**注**　① 按式 (59d)，电势零点已取在 $xy$ 平面上的无穷远处（即 $\theta=\pi/2,r\to\infty$），故球面电势 $u_0$ 是确定的。注意我们在计算极限过程中只是略去无穷小量，而不是只取主要项（即不仅略去无穷小量，而且略去有限量和比主要项低阶的无穷大量）。② 边界条件 (59c) 来源于 Gauss 定律，它表示球面带电 $Q$。这一边界条件是电动力学所特有的，初学的读者在写下上述定解问题时，可能存在的困难即在于此。

**第四步　求解定解问题。**　由于轴对称性，一般解由式 (9) 给出。代入无穷远处的边界条件 (59d)，可得

$$
\sum_{l=0}^{\infty}A_lr^lP_l(\cos\theta)=-E_0r\cos\theta=-E_0rP_1(\cos\theta).
$$

比较两边可得

$$
A_1=-E_0,\qquad A_l=0\quad(l\ne1). \tag{60}
$$

于是得到

$$
u(r,\theta)=-E_0r\cos\theta+\sum_{l=0}^{\infty}\frac{B_l}{r^{l+1}}P_l(\cos\theta). \tag{61}
$$

再代入球面的边界条件 (59b)，可得

$$
-E_0aP_1(\cos\theta)+\sum_{l=0}^{\infty}\frac{B_l}{a^{l+1}}P_l(\cos\theta)=u_0=u_0P_0(\cos\theta).
$$

比较两边可得

$$
B_0=u_0a,\qquad B_1=E_0a^3,\qquad B_l=0\quad(l\ne0,1). \tag{62}
$$

于是

$$
u(r,\theta)=-E_0r\cos\theta+\frac{E_0a^3}{r^2}\cos\theta+\frac{u_0a}{r}. \tag{63}
$$

再代入球面的另一边界条件 (59c)，注意到 $\cos\theta$ 在球面上的积分为零，上式中只有最后一项对积分有贡献，由此容易求出 $u_0a=Q/4\pi\epsilon_0$，所以

$$
u(r,\theta)=-E_0r\cos\theta+\frac{E_0a^3}{r^2}\cos\theta+\frac{Q}{4\pi\epsilon_0r}. \tag{64}
$$

**第五步　结果分析。**　以上结果中，第一项显然是原均匀电场的电势，第二项是球面上的感应面电荷产生的电势，它是电偶极子势，第三项是因导体球带电 $Q$ 而产生的电势。容易看出，当 $Q=0$ 时，结果中没有第三项，这就是电场中放入不带电导体球以后的电势分布；当 $a=0$ 时，结果中没有第二项，因而成为均匀场电势与点电荷电势的线性叠加，这也是一个非常直观的结果。当 $Q=0$ 时，球面附近的电场分布为

$$
E_r=-\left.\frac{\partial u}{\partial r}\right|_{r=a}=3E_0\cos\theta,\qquad
E_\theta=-\left.\frac1r\frac{\partial u}{\partial\theta}\right|_{r=a}=0. \tag{65}
$$

可见电场方向与球面垂直，如所期望。在 $\theta=\pi/2$ 处，场强为零；而在 $\theta=0,\pi$ 处，场强是原场强的三倍。所以，当外场较强时，这两处很容易被击穿。此外，当 $Q=0$ 时，球面上的电荷面密度为

$$
\sigma(\theta)=\epsilon_0E_r|_{r=a}=3\epsilon_0E_0\cos\theta. \tag{66}
$$

显然，球面上 $\theta\in[0,\pi/2)$ 部分带正电，$\theta\in(\pi/2,\pi]$ 部分带负电，而总带电量为 $Q=a^2\int_0^{2\pi}d\phi\int_0^\pi d\theta\sin\theta\sigma(\theta)=0$，如所期望。

**注**　① 本题的结果在 $r\to\infty$ 时有奇性。这是由于无穷远处有源（它产生均匀场），它使得边界条件产生奇性，从而使结果出现奇性。② 极角 $\theta$（其中 $\theta\in(0,\pi/2)$）处 $d\theta$ 角度内的元电荷 $\sigma(\theta)2\pi a^2\sin\theta d\theta$ 与极角 $\pi-\theta$ 处的相应元电荷形成元电偶极子，由于两元电荷的距离为 $2a\cos\theta$，故元电偶极子的偶极矩为 $4\pi a^3\sigma(\theta)\sin\theta\cos\theta d\theta=12\pi a^3\epsilon_0E_0\sin\theta\cos^2\theta d\theta$，总电偶极矩为 $p=12\pi a^3\epsilon_0E_0\int_0^{\pi/2}\sin\theta\cos^2\theta d\theta=4\pi a^3\epsilon_0E_0$，其方向为 $z$ 轴正向。根据电磁学，一个理想电偶极子 $\boldsymbol p$（电偶极子 $p=ql$ 当 $l\to0$、$q\to\infty$ 而 $ql$ 保持固定时称为理想电偶极子）产生的电势为 $\boldsymbol p\cdot\boldsymbol r/4\pi\epsilon_0r^3=E_0a^3\cos\theta/r^2$，这与式 (64) 中的第二项一致，说明感应电荷的分布对于球外来说等价于一个理想电偶极子。

**习题 1**　半径为 $a$ 的导体球壳，被一层绝缘体隔成两个半球壳，充电使上半球面电势为 $v_1$，下半球面电势为 $v_2$，计算球内外的电势分布。

**习题 2**　两个同心的球面，内球面半径为 $r_1$，具有电势 $u_0$，外球面半径为 $r_2$，具有电势 $u_1\cos^2\theta$，其中 $u_0$ 和 $u_1$ 均为常数，求区域 $r_1<r<r_2$ 中的电势分布。

## §2 一般球函数

### §2.1 非轴对称问题、连带 Legendre 函数

我们已经知道，在球坐标下对 Laplace 方程分离变量会出现本征值问题 (4)。如果问题的边界条件并不具有轴对称性，则除了上节已经讨论的特殊情况 $m=0$，我们还需要研究 $m\ne0$ 的情况，即连带 (associated) Legendre 方程的本征值问题。

作变换

$$
P(x)=(1-x^2)^{m/2}y(x),\qquad m>0 \tag{67}
$$

代入式 (4a)，化简可得

$$
(1-x^2)y''-2(m+1)xy'+[\lambda-m(m+1)]y=0. \tag{68}
$$

现将 Legendre 方程 (5a) 中的函数 $P(x)$ 写作 $f(x)$，对该方程求导 $m$ 次，利用高阶导数的 Leibniz 公式，有

$$
[(1-x^2)f'']^{(m)}=(1-x^2)f^{(m)''}-2mxf^{(m)'}-m(m-1)f^{(m)},
$$

$$
(2xf')^{(m)}=2xf^{(m)'}+2mf^{(m)},
$$

整理得

$$
(1-x^2)f^{(m)''}-2(m+1)xf^{(m)'}+[\lambda-m(m+1)]f^{(m)}=0. \tag{69}
$$

比较式 (68) 和 (69) 可知

$$
y(x)=f^{(m)}(x). \tag{70}
$$

这就是说，将 Legendre 方程的解求导 $m$ 次，就得到方程 (68) 的解。为了满足边界条件 $P(\pm1)=0$，$y(\pm1)$ 应该有限，所以 $f(\pm1)$ 也应该有限，我们已经知道，这只有当 $\lambda=l(l+1)\ (l\in\mathbb N)$ 时才有可能，相应的解是 $f(x)=P_l(x)$。于是我们求得了连带 Legendre 方程本征值问题的本征值

$$
\lambda_l=l(l+1),\qquad l=m,m+1,\cdots \tag{71}
$$

和相应的本征函数

$$
P_l^m(x)=(1-x^2)^{m/2}P_l^{(m)}(x),\qquad l=m,m+1,\cdots. \tag{72}
$$

这称为连带 Legendre 函数。由于 $l<m$ 时本征函数为零，所以 $l$ 的取值有以上的限制。注意，对于连带 Legendre 方程的本征值问题来说，$m$ 是事先给定的，它来源于 $\Phi(\phi)$ 的本征值问题。

与 $P_l^m(x)$ 线性独立的另一个解记作 $Q_l^m(x)$，它在 $x=\pm1$ 处有奇性，此处不作详细讨论。

**注**　① 由以上分析可以看出，求解连带 Legendre 方程最关键的一步是作变换 (67)。这一变换的引入是有章可循的。可以验证，$x=\pm1$ 都是方程 (4a) 的正则奇点，且在两点的指标都是 $\pm m/2$。变换 (67) 使得方程 (68) 在 $x=\pm1$ 处的指标都成为 $0$ 和 $-m$（直接分析方程 (68) 的指标可得同样的结果），这样就使得方程 (68) 的两个解在奇点 $x=\pm1$ 处不出现分数指标，而且其中有一个是 Taylor 级数。虽然我们是在常点 $x=0$ 附近求解，但对于方程在奇点处的行为的分析具有重要的作用。② 式 (72) 是连带 Legendre 函数的 Ferrer 定义，另有一种 Hobson 定义，右边多一个因子 $(-)^m$。定义不同，则有关的公式（如递推关系）也会有所不同，所以读者在使用不同的参考书时应该留意其所用的定义。③ 有些作者把连带 Legendre 函数称为连带 Legendre 多项式，这一名称并不恰当，因为连带 Legendre 函数中除了多项式因子 $P_l^{(m)}(x)$ 外，还有一个因子 $(1-x^2)^{m/2}$，后者并不是多项式（只当 $m$ 为偶数时才是）。当然，如果令 $x=\cos\theta$，则后者成为 $\sin^m\theta$，也是多项式。但这样一来，$P_l^m(\cos\theta)$ 成了两个宗量（$\cos\theta$ 和 $\sin\theta$）的多项式，这显然不是通常意义下的多项式。所以，还是称为连带 Legendre 函数比较妥当。

上面指出，只有满足 $f(\pm1)$ 从而 $y(\pm1)$ 也有限的解才能满足边界条件 $P(\pm1)=0$，实际上，即使 $y(x)$ 在 $x=\pm1$ 处有奇性，只要奇性小于 $(1-x^2)^{-m/2}$，还是可以满足 $P(\pm1)=0$，这容易由式 (67) 看出。但实际上并不存在其它满足边界条件的解。我们已经知道，Legendre 方程的任何一个级数解，只要不中断为多项式，它就会在 $x=1$ 或 $x=-1$ 中的至少一点有奇性。容易求出，$x=\pm1$ 处的指标都是 $s_1=s_2=0$，所以，在该处有奇性的解，其奇性为 $\ln(1\mp x)$，求导 $m$ 次后，奇性为 $(1\mp x)^{-m}$，相应的解 $P(x)$ 具有奇性 $(1-x)^{-m/2}$ 或 $(1+x)^{-m/2}$，当然不满足 $P(\pm1)=0$ 的边界条件。所以，只有 Legendre 方程的多项式解才对应于满足边界条件的解 $P(x)$。而只有当 $\lambda=l(l+1)\ (l\in\mathbb N)$ 时 Legendre 方程才有多项式解 $P_l(x)$，这就是上面得到的结论。我们已经知道，这时与 $P_l(x)$ 线性独立的另一个解 $Q_l(x)$ 在 $x=\pm1$ 处都有奇性，所以 $Q_l^m(x)$ 在 $x=\pm1$ 处的奇性为 $(1\mp x)^{-m/2}$。

可能有些读者已经看出，$\lambda=l(l+1)$ 而 $l<m$ 的情况略有微妙之处，所以需要作进一步的讨论。按上面的解法，这时有 $f(x)=\{P_l(x),Q_l(x)\}$，$y(x)=\{P_l^{(m)}(x),Q_l^{(m)}(x)\}=\{0,Q_l^{(m)}(x)\}$，而 $P(x)=\{0,Q_l^m(x)\}$，这是看来应该是矛盾的，因为二阶线性独立常微分方程总有两个线性独立解。实际上，上面的求解程序只保证 $(1-x^2)^{m/2}f^{(m)}(x)$ 是方程 (4a) 的解，而不能保证这样得到的解一定是非平庸的。那么，另一个解形式如何，能否满足边界条件呢？我们不妨直接求解方程 (68)，由于 $x=0$ 是常点，故可令 $y(x)=\sum_{k=0}^{\infty}a_kx^k$，代入方程就得到系数的递推关系

$$
a_{k+2}=\frac{(k+m-l)(k+m+l+1)}{(k+2)(k+1)}a_k.
$$

这样显然可以求出两个线性独立解。容易证明，两个级数解的收敛半径都是 1，如所期望。还可以证明，每个级数解在 $x=\pm1$ 处都发散。于是，任何级数解都在 $x=\pm1$ 中的至少一点发散。由于方程 (68) 在 $x=\pm1$ 处的指标都是 $s_1=0,s_2=-m$，所以，有奇性的解，其奇性为 $(1\mp x)^{-m}$，相应的解 $P(x)$ 具有奇性 $(1-x)^{-m/2}$ 或 $(1+x)^{-m/2}$，当然不满足 $P(\pm1)=0$ 的边界条件。由递推关系可以看出，当 $l<m$ 时，级数解不可能中断为多项式。于是得到结论：即使 $\lambda=l(l+1)$，只要 $l<m$，本征值问题 (4) 就没有非平庸解。

### §2.2 连带 Legendre 函数的微分表示

由 Legendre 多项式的微分表示 (15) 立得连带 Legendre 函数的微分表示

$$
P_l^m(x)=\frac{1}{2^ll!}(1-x^2)^{m/2}\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l. \tag{73}
$$

这也称为 Rodrigues 公式。

### §2.3 函数 $P_l^{-m}(x)$

按式 (72)，$P_l^{-m}(x)$ 并没有意义，因为求导 $-m$ 次没有意义。然而，我们可以按微分表示 (73) 定义函数 $P_l^{-m}(x)$，只要 $m\le l$。由于连带 Legendre 方程中只出现 $m^2$，对于 $m\to -m$ 的替换是不变的，所以我们有理由推测，这样定义的 $P_l^{-m}(x)$ 也是连带 Legendre 方程的解，而且是对应于同一本征值的本征函数（因为 $l$ 没变）。然而，根据 Sturm–Liouville 本征值问题的一般结论，对应于一个本征值只有一个线性独立的本征函数（只有周期性边界条件可能例外），所以，可以进一步推测，$P_l^{-m}(x)$ 与 $P_l^m(x)$ 实质上应该相同，即最多只差一个常数因子。事实上，可以证明

$$
P_l^{-m}(x)=(-)^m\frac{(l-m)!}{(l+m)!}P_l^m(x). \tag{74}
$$

对 $m>0$ 证明上式后，对于 $m<0$，利用上式，就有

$$
P_l^{-m}(x)=P_l^{|m|}(x)=(-)^{|m|}\frac{(l+|m|)!}{(l-|m|)!}P_l^{-|m|}(x)
=(-)^m\frac{(l-m)!}{(l+m)!}P_l^m(x),
$$

可见它对 $m<0$ 同样成立。如果对 $P_l^m(x)$ 中 $m$ 的正负不作限制，则 $l\ge|m|$。

下面证明式 (74)。根据微分表示，所需证明的式子可以改写为

$$
\frac{d^{l-m}}{dx^{l-m}}(x^2-1)^l
=\frac{(l-m)!}{(l+m)!}(x^2-1)^m\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l. \tag{75}
$$

先看右边，由高阶导数的 Leibniz 公式，有

$$
\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l
=\frac{d^{l+m}}{dx^{l+m}}[(x-1)^l(x+1)^l]
=\sum_{k=0}^{l+m}C_{l+m}^k\frac{d^k}{dx^k}(x-1)^l\frac{d^{l+m-k}}{dx^{l+m-k}}(x+1)^l
$$

$$
=\sum_{k=m}^{l}C_{l+m}^k\frac{d^k}{dx^k}(x-1)^l\frac{d^{l+m-k}}{dx^{l+m-k}}(x+1)^l,
$$

最后一步是因为 $k$ 必须满足 $k\le l$，$l+m-k\le l$，否则求导结果为零，这导致 $m\le k\le l$。将求和指标 $k$ 变为 $p=k-m$，并将导数算出来，可得

$$
\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l
=(l!)^2(l+m)!\sum_{p=0}^{l-m}
\frac{(x-1)^{l-m-p}(x+1)^p}{p!(m+p)!(l-p)!(l-m-p)!}.
$$

于是

$$
\frac{(l-m)!}{(l+m)!}(x^2-1)^m\frac{d^{l+m}}{dx^{l+m}}(x^2-1)^l
=(l!)^2(l-m)!\sum_{p=0}^{l-m}
\frac{(x-1)^{l-p}(x+1)^{m+p}}{p!(m+p)!(l-p)!(l-m-p)!}.
$$

再看左边，类似地有

$$
\frac{d^{l-m}}{dx^{l-m}}(x^2-1)^l
=\frac{d^{l-m}}{dx^{l-m}}[(x-1)^l(x+1)^l]
=\sum_{p=0}^{l-m}C_{l-m}^p\frac{d^p}{dx^p}(x-1)^l\frac{d^{l-m-p}}{dx^{l-m-p}}(x+1)^l
$$

$$
=(l!)^2(l-m)!\sum_{p=0}^{l-m}
\frac{(x-1)^{l-p}(x+1)^{m+p}}{p!(m+p)!(l-p)!(l-m-p)!}.
$$

可见式 (75) 左右相等。证毕。

### §2.4 连带 Legendre 函数的正交关系

作为 Sturm–Liouville 本征值问题的特例，相同 $m$ 不同 $l$ 的连带 Legendre 函数在区间 $[-1,1]$ 上有正交关系

$$
\int_{-1}^{1}P_l^m(x)P_{l'}^m(x)dx=0,\qquad l\ne l'. \tag{76}
$$

如同 Legendre 多项式的情况那样，上述正交关系也可以通过直接计算加以验证。如果 $m=0$，那么上式就是 Legendre 多项式的正交关系，这是上节已经得到的。如果上式对 $m>0$ 成立，那么根据式 (74)，它对 $m<0$ 显然也成立。所以只需对 $m>0$ 的情况加以验证。

利用式 (74)，我们就有

$$
\int_{-1}^{1}P_l^m(x)P_{l'}^m(x)dx
=(-)^m\frac{(l+m)!}{(l-m)!}\int_{-1}^{1}P_l^{-m}(x)P_{l'}^m(x)dx. \tag{77}
$$

利用微分表示，又有

$$
\int_{-1}^{1}P_l^{-m}(x)P_{l'}^m(x)dx
=\frac{1}{2^{l+l'}l!l'!}\int_{-1}^{1}
\left[\frac{d^{l-m}}{dx^{l-m}}(x^2-1)^l\right]
\left[\frac{d^{l'+m}}{dx^{l'+m}}(x^2-1)^{l'}\right]dx
$$

$$
=\frac{1}{2^{l+l'}l!l'!}\int_{-1}^{1}
\left[\frac{d^{l-m}}{dx^{l-m}}(x^2-1)^l\right]
d\left[\frac{d^{l'+m-1}}{dx^{l'+m-1}}(x^2-1)^{l'}\right],
$$

对上式分部积分，只要 $m>0$，则积出部分为零，于是有

$$
\int_{-1}^{1}P_l^{-m}(x)P_{l'}^m(x)dx
=-\frac{1}{2^{l+l'}l!l'!}\int_{-1}^{1}
\left[\frac{d^{l-m+1}}{dx^{l-m+1}}(x^2-1)^l\right]
\left[\frac{d^{l'+m-1}}{dx^{l'+m-1}}(x^2-1)^{l'}\right]dx
$$

$$
=-\int_{-1}^{1}P_l^{-(m-1)}(x)P_{l'}^{m-1}(x)dx,
$$

其中最后一步又用了微分表示。继续作分部积分，作了共 $m$ 次后得到

$$
\int_{-1}^{1}P_l^{-m}(x)P_{l'}^m(x)dx
=(-)^m\int_{-1}^{1}P_l(x)P_{l'}(x)dx.
$$

将这一结果代入式 (77)，就得到

$$
\int_{-1}^{1}P_l^m(x)P_{l'}^m(x)dx
=\frac{(l+m)!}{(l-m)!}\int_{-1}^{1}P_l(x)P_{l'}(x)dx. \tag{78}
$$

于是，由 Legendre 多项式的正交关系立得连带 Legendre 函数的正交关系。

### §2.5 连带 Legendre 函数的模

连带 Legendre 函数的模为

$$
\|P_l^m\|=\left\{\int_{-1}^{1}[P_l^m(x)]^2dx\right\}^{1/2}
=\sqrt{\frac{(l+m)!}{(l-m)!}\frac{2}{2l+1}}. \tag{79}
$$

其中第一等号是定义，第二等号是计算结果。上式对 $m>0$ 和 $m<0$ 均成立，当 $m=0$，它就退化为式 (28)。

式 (76) 和式 (79) 可以统一写为

$$
\int_{-1}^{1}P_l^m(x)P_{l'}^m(x)dx
=\frac{(l+m)!}{(l-m)!}\frac{2}{2l+1}\delta_{ll'}
=\|P_l^m\|^2\delta_{ll'}, \tag{80a}
$$

或其等价形式

$$
\int_0^\pi P_l^m(\cos\theta)P_{l'}^m(\cos\theta)\sin\theta\,d\theta
=\frac{(l+m)!}{(l-m)!}\frac{2}{2l+1}\delta_{ll'}
=\|P_l^m\|^2\delta_{ll'}. \tag{80b}
$$

当 $m>0$ 时，由式 (78) 和式 (29) 立得式 (80a)。当 $m<0$ 时，利用式 (74) 和式 (80a)，有

$$
\begin{aligned}
\int_{-1}^{1}P_l^m(x)P_{l'}^m(x)dx
&=\int_{-1}^{1}P_l^{-|m|}(x)P_{l'}^{-|m|}(x)dx\\
&=\left[\frac{(l-|m|)!}{(l+|m|)!}\right]^2
\int_{-1}^{1}P_l^{|m|}(x)P_{l'}^{|m|}(x)dx\\
&=\left[\frac{(l-|m|)!}{(l+|m|)!}\right]^2
\frac{(l+|m|)!}{(l-|m|)!}\frac{2}{2l+1}\delta_{ll'}\\
&=\frac{(l-|m|)!}{(l+|m|)!}\frac{2}{2l+1}\delta_{ll'}
=\frac{(l+m)!}{(l-m)!}\frac{2}{2l+1}\delta_{ll'}.
\end{aligned}
$$

可见式 (80a) 对 $m<0$ 同样成立。

### §2.6 广义 Fourier 级数

作为 Sturm–Liouville 本征值问题的特例，连带 Legendre 函数在区间 $[-1,1]$ 上是完备的，所以，区间 $[-1,1]$ 上的任一函数 $f(x)$，如果具有良好的解析性质，必可用 $\{P_l^m(x)\}_{l=m}^{\infty}$ 展开为广义 Fourier 级数

$$
f(x)=\sum_{l=m}^{\infty}f_lP_l^m(x),\qquad -1\le x\le1, \tag{81}
$$

利用式 (80)，容易求得展开系数为

$$
f_l=\frac{1}{\|P_l^m\|^2}\int_{-1}^{1}f(x)P_l^m(x)dx,\qquad l=m,m+1,m+2,\cdots. \tag{82}
$$

为确定起见，上面认为 $m>0$。给定一个函数 $f(x)$，即使是比较简单的初等函数，要计算上面的展开系数都是不容易的。

### §2.7 连带 Legendre 函数的递推关系

连带 Legendre 函数的递推关系都可以从 Legendre 多项式的递推关系导出。例如对式 (48) 求导 $m$ 次，对式 (50) 求导 $m-1$ 次，联立消去 $P_l^{(m-1)}(x)$，可得

$$
(l-m+1)P_{l+1}^{(m)}(x)-(2l+1)xP_l^{(m)}(x)+(l+m)P_{l-1}^{(m)}(x)=0,
$$

两边同乘以 $(1-x^2)^{m/2}$ 即得

$$
(l-m+1)P_{l+1}^m(x)-(2l+1)xP_l^m(x)+(l+m)P_{l-1}^m(x)=0. \tag{83}
$$

这是连带 Legendre 函数的递推关系中最简单的。其它递推关系均含有因子 $(1-x^2)^{1/2}$ 或其整数幂，此处不再给出。

### §2.8 球谐函数

Laplace 方程或 Helmholtz 方程在球坐标系中作分离变量 $u(\boldsymbol r)=R(r)Y(\theta,\phi)$ 后，角向部分 $Y(\theta,\phi)$ 满足球函数方程

$$
\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial Y}{\partial\theta}\right)
+\frac{1}{\sin^2\theta}\frac{\partial^2Y}{\partial\phi^2}+\lambda Y=0. \tag{84a}
$$

由于球坐标系的特点，函数 $Y(\theta,\phi)$ 应该满足以下两个自然边界条件

$$
Y(0,\phi),\quad Y(\pi,\phi)\quad\text{没有奇性}, \tag{84b}
$$

$$
Y(\theta,\phi+2\pi)=Y(\theta,\phi). \tag{84c}
$$

这里没有奇性指的是具有确切的定义 (well defined)，并且是有限的。上述方程和边界条件构成了一个本征值问题，不过是偏微分方程的本征值问题。进一步作分离变量 $Y(\theta,\phi)=H(\theta)\Phi(\phi)$，我们已经求得，$\Phi(\phi)$ 由式 (3) 给出，相应的本征值 $m^2$ 进入 $H(\theta)$ 的方程；对于确定的 $m$，$H(\theta)=P_l^m(\cos\theta)$，其中 $l\ge m$，本征值则由式 (71) 给出。

总结起来，对于一个确定的本征值 $\lambda_l=l(l+1)$，偏微分方程的本征值问题 (84) 有下列本征函数

$$
P_l^m(\cos\theta)\{e^{im\phi},e^{-im\phi}\},\qquad m=0,1,2,\cdots,l, \tag{85a}
$$

或

$$
P_l^{|m|}(\cos\theta)e^{im\phi},\qquad m=0,\pm1,\pm2,\cdots,\pm l, \tag{85b}
$$

由于 $P_l^{-m}(x)$ 与 $P_l^m(x)$ 存在关系 (74)，所以上述本征函数又等价于

$$
S_{lm}(\theta,\phi)=P_l^m(\cos\theta)e^{im\phi},\qquad m=0,\pm1,\pm2,\cdots,\pm l. \tag{85c}
$$

下面我们将采用最后一种写法，这些函数称为球谐函数 (spherical harmonics)。可见对应于一个本征值，其本征函数是 $2l+1$ 个线性独立的球谐函数，故本征值 $\lambda_l$ 的简并度为 $2l+1$。$l$ 称为球谐函数的阶。

物理上更常用的是归一化的球谐函数

$$
Y_{lm}(\theta,\phi)=(-)^m\frac{P_l^m(\cos\theta)}{\|P_l^m\|}\frac{e^{im\phi}}{\sqrt{2\pi}}
=(-)^m\sqrt{\frac{2l+1}{4\pi}\frac{(l-m)!}{(l+m)!}}P_l^m(\cos\theta)e^{im\phi},
$$

$$
l\in\mathbb N;\qquad m=0,\pm1,\pm2,\cdots,\pm l. \tag{86}
$$

其中各函数的归一化因子是一目了然的，至于 $(-)^m$ 则是相位因子，因为归一化以后相位因子还可以任取。相位因子的取法是历史上沿袭下来的。球谐函数有不同的定义，主要就在于相位因子的取法不同。定义不同，有关的公式也会有所不同，所以读者在使用时应该加以留意。

**注**　如果在定义 (86) 中使用 Hobson 的 $P_l^m(\cos\theta)$，则不需要因子 $(-)^m$。但注意对于 $Y_{lm}(\theta,\phi)$ 来说，这是同一个定义，而无关乎上面所述的相位因子的不同取法。

以上球谐函数在单位球面上满足正交归一关系（简称正一关系，orthonormal relation）

$$
\int_0^{2\pi}\int_0^\pi Y_{l'm'}^*(\theta,\phi)Y_{lm}(\theta,\phi)\sin\theta\,d\theta d\phi
=\delta_{ll'}\delta_{mm'}. \tag{87}
$$

上式的证明非常简单，只要将式 (86) 代入，利用 $\int_0^{2\pi}e^{i(m-m')\phi}d\phi=2\pi\delta_{mm'}$，然后就可以将式子中其余的 $m'$ 换为 $m$，再利用式 (80b)，即得结果。

球谐函数在球面上是完备的，这是因为 $\{e^{im\phi}\}_{m=-\infty}^{\infty}$ 在区间 $\phi\in[0,2\pi]$ 上是完备的，而 $\{P_l^m(\cos\theta)\}_{l=|m|}^{\infty}$ 在区间 $\theta\in[0,\pi]$ 上是完备的。所以，球面上的任意函数 $f(\theta,\phi)$，如果具有良好的解析性质，一定可以展开为

$$
f(\theta,\phi)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}f_{lm}Y_{lm}(\theta,\phi), \tag{88}
$$

由正一关系，容易求出展开系数为

$$
f_{lm}=\int_0^{2\pi}\int_0^\pi Y_{lm}^*(\theta,\phi)f(\theta,\phi)\sin\theta\,d\theta d\phi. \tag{89}
$$

本征值问题 (84) 的本征值与 $m$ 无关，所以，对于 Laplace 方程而言，径向方程的解仍由式 (8) 给出，于是我们得到 Laplace 方程在球坐标系下的一般解

$$
u(r,\theta,\phi)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\left(A_{lm}r^l+\frac{B_{lm}}{r^{l+1}}\right)Y_{lm}(\theta,\phi), \tag{90}
$$

其中的系数 $A_{lm}$ 和 $B_{lm}$ 由边界条件决定。对于轴对称问题，上式只包含 $m=0$ 的项，从而退化为式 (9)。

### §2.9 *球谐函数的加法公式

按式 (90)，在原点的邻域内，Laplace 方程的解是

$$
u(\boldsymbol r)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}A_{lm}r^lY_{lm}(\theta,\phi). \tag{91a}
$$

如果对坐标系作旋转，则几何点 $\boldsymbol r$ 的球坐标由 $(r,\theta,\phi)$ 变为 $(r,\tilde\theta,\tilde\phi)$，其中 $\tilde\theta$ 和 $\tilde\phi$ 是 $\theta$ 和 $\phi$ 的函数，函数的形式取决于绕什么轴旋转多大角度。由于 Laplace 算符是旋转不变量，故在新的球坐标系中，Laplace 方程的形式没有变化，只是球坐标 $\theta$ 和 $\phi$ 分别变成了 $\tilde\theta$ 和 $\tilde\phi$，所以，在新的球坐标系中，

$$
u(\boldsymbol r)=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\tilde A_{lm}r^lY_{lm}(\tilde\theta,\tilde\phi). \tag{91b}
$$

**注**　如果从式 (91b) 出发，将函数关系 $\tilde\theta=\tilde\theta(\theta,\phi)$ 和 $\tilde\phi=\tilde\phi(\theta,\phi)$ 代入，然后展开 $\sum_{m=-l}^{l}\tilde A_{lm}Y_{lm}(\tilde\theta,\tilde\phi)=\sum_{l'=0}^{\infty}\sum_{m'=-l'}^{l'}b_{ll'm'}Y_{l'm'}(\theta,\phi)$，则只能得到 $u(\boldsymbol r)=\sum_{l=0}^{\infty}\sum_{l'=0}^{\infty}\sum_{m'=-l'}^{l'}b_{ll'm'}r^lY_{l'm'}(\theta,\phi)$。可见，这样得到的结果比式 (91a) 复杂得多。$u(\boldsymbol r)$ 实际上由式 (91a) 给出，是源于 Laplace 算符的旋转不变性。可见后者在导出以上结果时具有关键的作用。

将 $u(\boldsymbol r)$ 当作 $r$ 的函数，而将 $\theta,\phi$ 或 $\tilde\theta,\tilde\phi$ 当作参数，则式 (91a) 和 (91b) 都是 $u(\boldsymbol r)$ 对 $r$ 的 Taylor 展开式，因此 $r^l$ 的系数应该相等，从而

$$
\sum_{m=-l}^{l}\tilde A_{lm}Y_{lm}(\tilde\theta,\tilde\phi)
=\sum_{m=-l}^{l}A_{lm}Y_{lm}(\theta,\phi). \tag{92}
$$

上式对于任何调和函数 $u(\boldsymbol r)$ 均必须成立。换句话说，对于任何一组 $\{\tilde A_{lm}\}$，总存在一组 $\{A_{lm}\}$ 使得上式可以成立。特别地，取 $u(\boldsymbol r)$ 使上式左边只有一项，就得到

$$
Y_{lm}(\tilde\theta,\tilde\phi)=\sum_{m'=-l}^{l}C_{m'}Y_{lm'}(\theta,\phi). \tag{93}
$$

对于不同的 $l$ 和 $m$，展开系数 $C_{m'}$ 是不同的，即 $C_{m'}$ 实际上还应该有两个指标 $l$ 和 $m$，为了书写简单，上面没有明显写出来。上式表明，$\{Y_{lm}(\tilde\theta,\tilde\phi)\}_{m=-l}^{l}$ 一定可以用 $\{Y_{lm}(\theta,\phi)\}_{m=-l}^{l}$ 线性表出。其逆命题也容易证明。所以我们得到结论：不同坐标系中的同阶球谐函数是等价的。由此马上得到一个推论：不同阶的球谐函数总是线性无关的。事实上，当 $l'\ne l$，$Y_{l'm'}(\theta,\phi)$ 与 $Y_{lm}(\theta,\phi)$ 互相正交，当然也就线性无关。至于 $Y_{l'm'}(\tilde\theta,\tilde\phi)$，按上面的结论，它可以用 $\{Y_{l'm''}(\theta,\phi)\}_{m''=-l'}^{l'}$ 表出，当然也就与 $Y_{lm}(\theta,\phi)$ 互相正交，从而线性无关了。

容易证明

$$
\delta(\cos\theta-\cos\theta')\delta(\phi-\phi')
=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\theta,\phi)Y_{lm}^*(\theta',\phi'). \tag{94}
$$

这只要将左边的函数按式 (88) 展开，然后按式 (89) 算出展开系数即可。现在，将球坐标系的极轴转到 $(\theta',\phi')$ 方向，在新的球坐标系中，$(\theta',\phi')$ 方向和 $(\theta,\phi)$ 方向的极角分别是 $0$ 和 $\gamma$，其中 $\gamma$ 是两个方向的夹角，满足

$$
\cos\gamma=\cos\theta\cos\theta'+\sin\theta\sin\theta'\cos(\phi-\phi'). \tag{95}
$$

在新的球坐标系中，只要 $\gamma=0$，两个方向便会重合，所以

$$
\delta(\cos\theta-\cos\theta')\delta(\phi-\phi')
=\frac{1}{2\pi}\delta(\cos\gamma-1), \tag{96}
$$

其中因子 $1/2\pi$ 是为了保持两边在球面上的积分均为 1。又容易证明

$$
\frac{1}{2\pi}\delta(\cos\gamma-1)=\sum_{l=0}^{\infty}\frac{2l+1}{4\pi}P_l(\cos\gamma). \tag{97}
$$

综合以上各式可得

$$
\sum_{l=0}^{\infty}\frac{2l+1}{4\pi}P_l(\cos\gamma)
=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}Y_{lm}(\theta,\phi)Y_{lm}^*(\theta',\phi'). \tag{98}
$$

注意到 $P_l(\cos\gamma)\propto Y_{l0}(\gamma,\beta)$（其中 $\beta$ 是 $(\theta,\phi)$ 方向在新的球坐标系中的方位角），以及上面不同阶球谐函数线性无关的结论，就得到

$$
P_l(\cos\gamma)=\frac{4\pi}{2l+1}\sum_{m=-l}^{l}Y_{lm}(\theta,\phi)Y_{lm}^*(\theta',\phi'). \tag{99a}
$$

或

$$
P_l(\cos\gamma)=\sum_{m=-l}^{l}\frac{(l-m)!}{(l+m)!}P_l^m(\cos\theta)P_l^m(\cos\theta')e^{im(\phi-\phi')}. \tag{99b}
$$

上式称为球谐函数的加法公式，它把新坐标系中的球谐函数用旧坐标系中的同阶球谐函数表示出来。

考虑空间任意两点 $\boldsymbol r$ 和 $\boldsymbol r'$，其球坐标分别是 $(r,\theta,\phi)$ 和 $(r',\theta',\phi')$，令 $r_>$（$r_<$）表示 $r$ 和 $r'$ 中之较大（小）者，仍以 $\gamma$ 表示 $\boldsymbol r$ 和 $\boldsymbol r'$ 的夹角，则

$$
|\boldsymbol r-\boldsymbol r'|=\sqrt{r^2+r'^2-2rr'\cos\gamma}
=\sqrt{r_>^2+r_<^2-2r_>r_<\cos\gamma}
=r_>\left[1+\left(\frac{r_<}{r_>}\right)^2-2\left(\frac{r_<}{r_>}\right)\cos\gamma\right]^{1/2}.
$$

利用 Legendre 多项式的母函数展开式 (40)，可得

$$
\frac{1}{|\boldsymbol r-\boldsymbol r'|}
=\sum_{l=0}^{\infty}\frac{r_<^l}{r_>^{l+1}}P_l(\cos\gamma),
$$

再利用上面得到的加法公式，就有

$$
\frac{1}{|\boldsymbol r-\boldsymbol r'|}
=\sum_{l=0}^{\infty}\sum_{m=-l}^{l}\frac{4\pi}{2l+1}\frac{r_<^l}{r_>^{l+1}}Y_{lm}(\theta,\phi)Y_{lm}^*(\theta',\phi'). \tag{100}
$$

这是一个有用的公式。

## 补充习题

1. 设有半径为 $a$ 的导热半球，球面保持恒定温度 $u_0$。底面的情况为：(1) 保持零度；(2) 绝热。求解两种情况下半球内的稳定温度分布。

2. 设有细导线构成的圆环，半径为 $a$，环上均匀带电 $4\pi\epsilon_0q$，求圆环周围空间中的电势。（提示：先计算轴上的电势。）

3. 半径为 $a$ 的均匀导热球，在温度为零度的空气中被平行热流所照射，热流强度为常数 $q_0$，球与周围空气按 Newton 冷却定律交换热量，求解球内的稳定温度分布。如果球只吸收热流，而不按 Newton 冷却定律向周围散发热量，则结果如何？（提示：取球心为坐标原点，极轴反平行于热流方向，则边界条件为 $(\partial u/\partial r+hu)|_{r=a}=f(\theta)$，其中 $h=b/k$，$b$ 为 Newton 冷却定律中的热交换系数，而 $f(\theta)=(q_0/k)\cos\theta\ (0\le\theta\le\pi/2)$，$f(\theta)=0\ (\pi/2\le\theta\le\pi)$。）

4. 长为 $l$ 的柔软轻弦（重力可以忽略），一端固定在以角速度 $\omega$ 旋转的竖直轴上，一端自由。由于离心力的作用，弦的平衡位置应该是水平线。当受到扰动时，弦可以相对于该水平线作横振动。设初始位移为 $\varphi(x)$，初始速度为 $\psi(x)$，求解该横振动问题。

5. 已知半径为 $a$ 的球面上的电势为 $u|_{r=a}=2\sin^2\theta(\sin2\phi+1)$，球内外均没有电荷，求空间各处的电势分布。
