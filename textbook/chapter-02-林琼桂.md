# 第二章　复变函数的积分\*

## 目录

§1　复积分的定义和基本性质　27  
§1.1　曲线的方向　27  
§1.2　复积分的定义与存在性　27  
§1.3　复积分的基本性质　29

§2　复积分的计算　29  
§2.1　\*用定义计算　29  
§2.2　曲线积分方法　30  
§2.3　参数方程法　31

§3　Cauchy 积分定理　32  
§3.1　Cauchy 积分定理　32  
§3.2　\*Cauchy 积分定理的证明　33  
§3.3　Cauchy 积分定理的推论　36  
§3.4　复通区域的 Cauchy 积分定理　37

§4　Cauchy 积分公式及其推论　39  
§4.1　Cauchy 积分公式　39  
§4.2　\*无界区域的 Cauchy 积分公式　41  
§4.3　解析函数的高阶导数　42  
§4.4　\*Liouville 定理　43  
§4.5　\*代数基本定理　43

补充习题　44

本章的 Cauchy 积分定理是整个解析函数理论的基础。

## §1　复积分的定义和基本性质

### §1.1　曲线的方向

本章提到的曲线一般都指光滑或逐段光滑的平面曲线。若一段曲线的方程为 $y=f(x)$，则光滑指的是 $f'(x)$ 连续。若其方程为参数方程 $x=x(t)$、$y=y(t)$，则光滑指的是 $x'(t)$ 和 $y'(t)$ 连续。由有限条光滑曲线衔接而成的曲线称为逐段光滑曲线。折线就是最简单的逐段光滑曲线。曲线的方向定义如下。

1. 简单曲线：没有重点的曲线称为简单曲线。图 1 中的 a 是简单曲线，而 b 则不是。简单曲线的方向由起点指向终点。所以规定了起点和终点就确定了它的方向。

2. 围线（contour）：即简单闭曲线。图 1 中的 c 是简单闭曲线，而 d 则不是。

如果沿着围线走，其所包围的区域在左边，则该方向称为正方向，其实就是逆时针方向。

### §1.2　复积分的定义与存在性

一元实变函数的积分定义在 $x$ 轴上的有限区间内（无限区间的广义积分通过极限过程得到），而复变函数的积分（简称复积分）总是定义在曲线上，参看图 2。

**定义（复积分）**　设函数 $f(z)$ 沿曲线 $C$ 有定义，在 $C$ 上从起点 $a$ 到终点 $b$ 取分点

$$
a=z_0,z_1,\cdots,z_{k-1},z_k,\cdots,z_{n-1},z_n=b
$$

将 $C$ 分为 $n$ 个弧段，在 $z_{k-1}$ 至 $z_k$ 的弧段上任取一点 $\zeta_k$，作和数

$$
S_n=\sum_{k=1}^n f(\zeta_k)\Delta z_k,
$$

其中 $\Delta z_k=z_k-z_{k-1}$。若 $n\to\infty$ 且 $\max_{1\le k\le n}|\Delta z_k|\to0$ 时，$S_n\to I$，则称 $f(z)$ 沿 $C$ 可积，并称 $I$ 为 $f(z)$ 沿 $C$ 的积分，记作

$$
I=\int_C f(z)\,dz.
$$

$C$ 称为积分路径，沿相反方向的积分记作 $\int_{C^-}f(z)\,dz$。

**注**　当 $b=a$，即 $C$ 为围线时，如果没有特别说明，则沿围线 $C$ 积分总是指逆时针方向的积分。

以上积分的定义与实变函数积分的定义是很类似的。在复变函数论中，积分理论具有非常基本的地位。事实上，关于复积分的 Cauchy 积分定理是解析函数理论的基础，解析函数的微分性质的证明大都是由这个定理出发的。

定义了复积分，首先遇到的问题是，给定函数和曲线，如何判断积分是否存在。与实变情况类似，我们有下述

![Figure 1](...)

图 1: a 是简单曲线，而 b 不是；c 是简单闭曲线，而 d 不是

![Figure 2](...)

图 2: $f(z)$ 沿曲线 $C$ 的积分

**定理（积分存在）**　设函数 $f(z)=u(x,y)+iv(x,y)$ 沿曲线 $C$ 连续，则 $f(z)$ 沿 $C$ 可积，且

$$
\int_C f(z)\,dz
=\int_C u\,dx-v\,dy+i\int_C v\,dx+u\,dy. \tag{1}
$$

**注**　形式上，只要将 $f(z)=u(x,y)+iv(x,y)$ 和 $dz=dx+i\,dy$ 代入左边，分开实部和虚部，即可得到右边的结果。这可以帮助我们把握这一公式，虽然不能用这种方法来证明它。

**证明**　设 $z_k=x_k+iy_k$，记 $\Delta x_k=x_k-x_{k-1}$，$\Delta y_k=y_k-y_{k-1}$，则 $\Delta z_k=\Delta x_k+i\Delta y_k$。又设 $\zeta_k=\xi_k+i\eta_k$，记 $u_k=u(\xi_k,\eta_k)$，$v_k=v(\xi_k,\eta_k)$，则 $f(\zeta_k)=u_k+iv_k$，于是

$$
\sum_{k=1}^n f(\zeta_k)\Delta z_k
=\sum_{k=1}^n(u_k\Delta x_k-v_k\Delta y_k)
+i\sum_{k=1}^n(v_k\Delta x_k+u_k\Delta y_k),
$$

由于 $f(z)$ 沿 $C$ 连续，故 $u(x,y)$ 和 $v(x,y)$ 沿 $C$ 连续，所以上式右边两项均有极限，于是左边也有极限，即 $f(z)$ 沿 $C$ 可积。上式的极限正是式（1）。证毕。

式（1）同时也给出了复积分的计算方法，它把复积分的计算转化为实变函数的曲线积分的计算，我们将在下一节进一步讨论。

### §1.3　复积分的基本性质

设 $f(z)$ 和 $g(z)$ 沿 $C$ 连续，复积分有下列基本性质。

1.

$$
\int_C[\alpha f(z)+\beta g(z)]\,dz
=\alpha\int_C f(z)\,dz+\beta\int_C g(z)\,dz.
$$

2.

$$
\int_{C_1+C_2}f(z)\,dz
=\int_{C_1}f(z)\,dz+\int_{C_2}f(z)\,dz.
$$

3.

$$
\int_{C^-}f(z)\,dz=-\int_C f(z)\,dz.
$$

4.

$$
\left|\int_C f(z)\,dz\right|
\le\int_C|f(z)|\,|dz|
=\int_C|f(z)|\,ds.
$$

最后一个积分中的 $ds=|dz|$ 就是弧长的微分。这些性质都可以直接由积分的定义来证明。比如，利用三角不等式，有

$$
\left|\sum_{k=1}^n f(\zeta_k)\Delta z_k\right|
\le\sum_{k=1}^n|f(\zeta_k)||\Delta z_k|,
$$

取极限即得性质 4。

## §2　复积分的计算

本节讨论复积分的计算，这有三方面的目的。第一，加深对复积分的了解；第二，通过计算获得一些有用的具体结果；第三，建立 Cauchy 积分定理的特例，从而对这一抽象的定理获得感性认识。

### §2.1　\*用定义计算

对于比较简单的函数，可以直接用定义来计算。下面就计算两个例子。

**例 1**　$C$ 是 $a$ 到 $b$ 的任一简单曲线，求 $\int_C dz$。

**解**　$f(z)=1$，则 $f(\zeta_k)=1$。按定义作和数，有

$$
S_n=\sum_{k=1}^n f(\zeta_k)\Delta z_k
=\sum_{k=1}^n(z_k-z_{k-1})=z_n-z_0=b-a,
$$

故当 $n\to\infty$、$\max_{1\le k\le n}|\Delta z_k|\to0$ 时，$S_n\to b-a$，即

$$
\int_C dz=b-a. \tag{2}
$$

**例 2**　$C$ 是 $a$ 到 $b$ 的任一简单曲线，求 $\int_C z\,dz$。

**解**　$f(z)=z$，则 $f(\zeta_k)=\zeta_k$。按定义作和数，有

$$
S_n=\sum_{k=1}^n f(\zeta_k)\Delta z_k
=\sum_{k=1}^n\zeta_k(z_k-z_{k-1}).
$$

先取 $\zeta_k=z_k$，得

$$
S_n^{(1)}=\sum_{k=1}^n z_k(z_k-z_{k-1}),
$$

再取 $\zeta_k=z_{k-1}$，得

$$
S_n^{(2)}=\sum_{k=1}^n z_{k-1}(z_k-z_{k-1}).
$$

由于 $f(z)=z$ 是连续函数，按上节定理，积分存在，记作 $I$，则当 $n\to\infty$、$\max_{1\le k\le n}|\Delta z_k|\to0$ 时，应有 $S_n^{(1)}\to I$，$S_n^{(2)}\to I$，于是 $S_n^{(1)}+S_n^{(2)}\to2I$。另一方面，

$$
S_n^{(1)}+S_n^{(2)}
=\sum_{k=1}^n(z_k+z_{k-1})(z_k-z_{k-1})
=\sum_{k=1}^n(z_k^2-z_{k-1}^2)
=z_n^2-z_0^2=b^2-a^2,
$$

因此，$S_n^{(1)}+S_n^{(2)}\to b^2-a^2$。比较即得 $I=(b^2-a^2)/2$，即

$$
\int_C z\,dz=\frac12(b^2-a^2). \tag{3}
$$

### §2.2　曲线积分方法

从上一小节看到，直接用定义来计算复积分，即使对于象 $f(z)=z$ 这样简单的函数，都需要一定的技巧。对于较复杂的函数，这种方法显然是不现实的。注意到上节的定理不仅解决了积分的存在性问题，而且式（1）同时也给出了计算复积分的一种方法，即把复积分转化为实变函数的曲线积分来计算。在简单情况下，可以找到曲线积分的原函数，这时就很容易得到结果。为了叙述方便，我们把这种方法称为曲线积分方法。

**例 3**　用曲线积分方法重新计算例 1。

**解**　$f(z)=1$，即 $u(x,y)=1$，$v(x,y)=0$。按式（1），有

$$
\int_C dz
=\int_C dx+i\int_C dy
=x\big|_a^b+iy\big|_a^b
=(x+iy)\big|_a^b
=z\big|_a^b=b-a.
$$

**例 4**　用曲线积分方法重新计算例 2。

**解**　$f(z)=z$，即 $u(x,y)=x$，$v(x,y)=y$。按式（1），有

$$
\begin{aligned}
\int_C z\,dz
&=\int_C x\,dx-y\,dy+i\int_C y\,dx+x\,dy\\
&=\int_C d\left[\frac12(x^2-y^2)\right]+i\int_C d(xy)\\
&=\frac12(x^2-y^2+i2xy)\big|_a^b
=\frac12(x+iy)^2\big|_a^b
=\frac12z^2\big|_a^b
=\frac12(b^2-a^2).
\end{aligned}
$$

由以上两例看到，积分的结果只与起点和终点有关，而与路径无关，这是与被积函数的解析性紧密相关的。事实上，这一方法可以用于任何具体的解析函数，如指数、三角函数等，并得到类似的结果。后面将会看到，对于解析函数，存在着与实变积分类似的 Newton–Leibniz 公式。

如果被积函数不是解析的，比如 $\bar z$，其积分一般都依赖于路径，而不仅仅是起点和终点。这时用上面的曲线积分方法并不方便。实际上，实变函数中的曲线积分在一般情况下还是要化成定积分来计算的。所以，下面就讨论计算复积分的参数方程法。

### §2.3　参数方程法

设曲线 $C$ 的参数方程为 $z=z(t)$（$\alpha\le t\le\beta$），$f(z)$ 沿 $C$ 连续，则

$$
\int_C f(z)\,dz
=\int_\alpha^\beta f[z(t)]z'(t)\,dt. \tag{4}
$$

这一公式将复积分化为实变数 $t$ 的定积分。形式上，只要将 $z=z(t)$ 和 $dz=z'(t)dt$ 代入左边，即可得到右边的结果。当然这不是证明。右边的定积分可能需要分开实部和虚部来计算，但若能直接找到原函数，也可以不必分开。

式（4）证明如下。由式（1）出发，根据微积分中曲线积分的计算公式，有

$$
\int_C f(z)\,dz
=\int_\alpha^\beta[u(t)x'(t)-v(t)y'(t)]\,dt
+i\int_\alpha^\beta[v(t)x'(t)+u(t)y'(t)]\,dt,
$$

其中 $u(t)=u[x(t),y(t)]$，$v(t)=v[x(t),y(t)]$，上式可以改写为

$$
\int_C f(z)\,dz
=\int_\alpha^\beta[u(t)+iv(t)][x'(t)+iy'(t)]\,dt.
$$

这就是式（4）。

利用式（4），马上可以得到例 1 和例 2 的结果。它们是下述更一般结果的特例。

**例 5**　由于 $z^n(t)z'(t)$ 的原函数显然是 $z^{n+1}(t)/(n+1)$，故有

$$
\int_a^b z^n\,dz
=\frac{1}{n+1}(b^{n+1}-a^{n+1}),\qquad n=0,1,2,\cdots. \tag{5}
$$

结果与曲线 $C$ 的细节无关，而只依赖于起点和终点。如果 $C$ 是围线，即 $b=a$，则

$$
\int_C z^n\,dz=0,\qquad n=0,1,2,\cdots,\quad \forall\text{围线 }C. \tag{6}
$$

式（6）可由式（5）经 $b\to a$ 的极限过程得到。也可以将围线 $C$ 分为两段，再由式（5），两段的积分互相抵消，从而得到式（6）。注意式（6）对任意围线成立。

**注**　由式（6）容易推知，对于任意多项式 $P_n(z)$ 和围线 $C$，均有 $\int_C P_n(z)\,dz=0$。由此我们可以推测，对于一般的解析函数 $f(z)$，亦有 $\int_C f(z)\,dz=0$，这基本上就是 Cauchy 积分定理了。

下面再计算一个例子。

**例 6**　计算积分

$$
\int_{|z-a|=\rho}\frac{dz}{(z-a)^n},
$$

其中 $n\in\mathbb Z$，$a\in\mathbb C$。

**解**　当 $n\le0$，由式（6）可知结果为 0，这与下面的结果一致。圆周 $|z-a|=\rho$ 的参数方程为 $z(\theta)=a+\rho e^{i\theta}$（$0\le\theta\le2\pi$），这里我们将参数写作 $\theta$，因为它是角度。按式（4），有

$$
\int_{|z-a|=\rho}\frac{dz}{(z-a)^n}
=\int_0^{2\pi}\frac{i\rho e^{i\theta}\,d\theta}{\rho^n e^{in\theta}}
=\frac{i}{\rho^{n-1}}\int_0^{2\pi}e^{-i(n-1)\theta}\,d\theta.
$$
如果 $n\neq1$，易得原函数为 $-e^{-i(n-1)\theta}/(n-1)\rho^{n-1}$，代入上下限，结果为 0。若 $n=1$，则结果显然为 $2\pi i$。所以

$$
\int_{|z-a|=\rho}\frac{dz}{(z-a)^n}
=\begin{cases}
2\pi i,&(n=1),\\
0,&(n\in\mathbb Z,\ n\neq1),
\end{cases} \tag{7}
$$

这是一个重要的结果。

**注**　① 这一结果与 $a$ 和 $\rho$ 的数值无关。② 后面将证明，这一结果对任何包围 $a$ 点的围线成立。③ $n=1$ 时结果不为 0，这是因为被积函数在积分围线内存在奇点 $z=a$。由此可以进一步推测，Cauchy 积分定理中的 $f(z)$ 应该在围线 $C$ 及其内部解析，而不仅仅是在围线 $C$ 上解析。④ 虽然 $n=2,3,\cdots$ 时，结果也为 0，但这只是一个具体的结果，没有必然性。

**习题 1**　计算 $\int_{C_1}\bar z\,dz$ 和 $\int_{C_2}\bar z\,dz$，其中 $C_1$ 是上半单位圆（$z=1\to z=-1$），$C_2$ 是下半单位圆（$z=1\to z=-1$）。

**习题 2**　设 $f(z)$ 在原点的邻域内连续，试证

$$
\lim_{r\to0}\int_0^{2\pi}f(re^{i\theta})\,d\theta=2\pi f(0).
$$

## §3　Cauchy 积分定理

### §3.1　Cauchy 积分定理

由上节看到，对于简单的解析函数 $P_n(z)$，有 $\int_C P_n(z)\,dz=0$，$C$ 为任意围线。实际上，这是下面定理的一个特例。

**定理（Cauchy 积分定理）**　设函数 $f(z)$ 在单通区域 $D$ 内解析，$C$ 为 $D$ 内的任一围线，则

$$
\int_C f(z)\,dz=0. \tag{8}
$$

这是复变函数论中最基本的定理。

如果假定 $f'(z)$ 连续，则上述定理的证明是很容易的。这时 $u$、$v$ 的一阶偏导数连续，可将 Green 公式应用于式（1）而得

$$
\int_C f(z)\,dz
=-\iint_G\left(\frac{\partial v}{\partial x}+\frac{\partial u}{\partial y}\right)dx\,dy
+i\iint_G\left(\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}\right)dx\,dy,
$$

其中 $G$ 是围线 $C$ 所包围的区域。由 CR 条件，上式右边两个积分都为 0，故得式（8）。但我们知道，$f(z)$ 解析的定义只是 $f'(z)$ 存在，并不一定连续，所以上面证明的只是一种特殊情况。下一小节将介绍严格的证明，供有兴趣的读者参考。

值得注意的是，Cauchy 在 1825 年给出上述定理，那时候 $f(z)$ 解析的定义是 $f'(z)$ 连续，所以上面的证明是严格的。1900 年，Goursat 发表了新的证明，免去了 $f'(z)$ 连续的条件。也就是说，只要 $f'(z)$ 存在，Cauchy 积分定理就成立。此后 $f(z)$ 解析的定义才改为现在的样子。这无疑是一个实质性的进步，但期间经历了七十多年的时间。

Cauchy 积分定理也可以表述为

**定理（Cauchy 积分定理的等价表述）**　$C$ 为复平面上的围线，$D$ 是 $C$ 所包围的单通区域，函数 $f(z)$ 在闭域 $\overline D=D+C$ 上解析，则

$$
\int_C f(z)\,dz=0.
$$

这与上面的形式是等价的。但从这一形式容易看出它和下面强化的形式有什么区别。

**定理（Cauchy 积分定理的强化形式）**　$C$ 为复平面上的围线，$D$ 是 $C$ 所包围的单通区域，函数 $f(z)$ 在 $D$ 内解析，在 $\overline D=D+C$ 上连续，则

$$
\int_C f(z)\,dz=0.
$$

应当指出，即使证明了上面的定理，要得到这一强化形式也并不是轻而易举的。

### §3.2　\*Cauchy 积分定理的证明

先证明下面的

**引理**　设 $f(z)$ 在区域 $D$ 内连续，$\Gamma$ 是 $D$ 内的简单曲线或围线，则 $\forall\varepsilon>0$，总可找到内接于 $\Gamma$ 且完全位于 $D$ 内的折线 $P$，使得

$$
\left|\int_\Gamma f(z)\,dz-\int_P f(z)\,dz\right|<\varepsilon. \tag{9}
$$

**注**　① 这一引理的大意是沿曲线的积分总可以用沿折线的积分来任意逼近。② 这里只要求 $f(z)$ 连续，而不必解析。③ 区域 $D$ 也不必是单通的。

**证明**　设曲线 $\Gamma$ 的长度为 $L$。在区域 $D$ 内取区域 $G$，使 $\overline G\subset D$ 而 $\Gamma$ 完全在 $G$ 内。由于 $f(z)$ 在区域 $D$ 内连续，故在 $\overline G$ 上一致连续。也就是说，$\forall\varepsilon>0$，$\exists\delta>0$，使得对于 $G$ 内的任意两点 $z'$ 和 $z''$，只要 $|z'-z''|<\delta$，就有 $|f(z')-f(z'')|<\varepsilon/2L$。今在 $\Gamma$ 上从起点 $a$ 到终点 $b$ 取分点

$$
a=z_0,z_1,\cdots,z_{k-1},z_k,\cdots,z_{n-1},z_n=b.
$$

依次连接各分点作成折线 $P$，如图 3 所示。我们将 $\Gamma$ 上从 $z_{k-1}$ 到 $z_k$ 的弧段记作 $\Gamma_k$，相应的线段记作 $P_k$，并记 $\Gamma_k$ 的长度为 $L_k$。只要分点取得足够密，就可以使所有 $L_k<\delta$，并使所有 $P_k$ 完全落在 $G$ 内。这样，对任意 $k$，$\Gamma_k$ 上的任意点 $z$ 都满足 $|f(z)-f(z_k)|<\varepsilon/2L$，$P_k$ 上的任意点亦然。由上节例 1，

$$
\int_{\Gamma_k}f(z_k)\,dz
=\int_{P_k}f(z_k)\,dz
=f(z_k)(z_k-z_{k-1}),
$$

所以

$$
\begin{aligned}
\left|\int_\Gamma f(z)\,dz-\int_P f(z)\,dz\right|
&=\left|\sum_{k=1}^n\int_{\Gamma_k}f(z)\,dz-\sum_{k=1}^n\int_{P_k}f(z)\,dz\right|\\
&=\left|\sum_{k=1}^n\int_{\Gamma_k}[f(z)-f(z_k)]\,dz-\sum_{k=1}^n\int_{P_k}[f(z)-f(z_k)]\,dz\right|\\
&\le\sum_{k=1}^n\left|\int_{\Gamma_k}[f(z)-f(z_k)]\,dz\right|
+\sum_{k=1}^n\left|\int_{P_k}[f(z)-f(z_k)]\,dz\right|\\
&\le\sum_{k=1}^n\int_{\Gamma_k}|f(z)-f(z_k)|\,|dz|
+\sum_{k=1}^n\int_{P_k}|f(z)-f(z_k)|\,|dz|.
\end{aligned}
$$

![Figure 3](...)

图 3: 引理的证明

$$
\begin{aligned}
&<\frac{\varepsilon}{2L}\left(\sum_{k=1}^n\int_{\Gamma_k}|dz|+\sum_{k=1}^n\int_{P_k}|dz|\right)\\
&\le\frac{\varepsilon}{2L}\,2L=\varepsilon.
\end{aligned}
$$

其中第三步用了三角不等式，第四步用了积分的性质 4，第五步用了 $|f(z)-f(z_k)|<\varepsilon/2L$（在 $\Gamma_k$ 或 $P_k$ 上均成立），最后一步用了 $\sum_{k=1}^n\int_{\Gamma_k}|dz|=\int_\Gamma|dz|=L$，而 $\sum_{k=1}^n\int_{P_k}|dz|\le L$。证毕。

现在证明 Cauchy 积分定理。$C$ 是单通区域 $D$ 内的任一围线，由于 $f(z)$ 在 $D$ 内解析，故必在 $D$ 内连续，根据上述引理，$\forall\varepsilon>0$，总可找到内接于 $C$ 且完全位于 $D$ 内的闭折线 $P$（即多边形），使

$$
\left|\int_C f(z)\,dz-\int_P f(z)\,dz\right|<\varepsilon.
$$

如果能证明，对于 $D$ 内的任意闭折线 $P$，都有

$$
\int_P f(z)\,dz=0, \tag{10}
$$

则上式变成 $|\int_C f(z)\,dz|<\varepsilon$，但 $\varepsilon$ 是任意的，故 $|\int_C f(z)\,dz|=0$。反过来，如果 Cauchy 积分定理成立，则式（10）也应该成立，因为 $P$ 也是 $D$ 内的围线。因此问题归结为证明式（10）。

由于 $D$ 是单通区域，所以 $P$ 及其内部都完全在 $D$ 内。今适当连接 $P$ 的各顶点作对角线，将 $P$ 所包围的多边形分解为若干三角形，如图 4a 所示，则沿各三角形边界的积分都存在，且其和等于沿 $P$ 的积分，这是因为在和式中，沿各对角线的积分都出现两次，而且方向相反，所以互相抵消。如果能证明沿 $D$ 内任一三角形边界 $T$ 的积分为 0：

$$
\int_T f(z)\,dz=0, \tag{11}
$$

则式（10）成立。反过来，如果式（10）成立，则式（11）也应成立，因为 $T$ 就是最简单的闭折线。因此问题归结为证明式（11）。为此，记

$$
M=\left|\int_T f(z)\,dz\right|,
$$

只要证明 $M=0$ 即可。

由于 $D$ 是单通区域，所以 $T$ 及其内部都完全在 $D$ 内。将 $T$ 及其内部所构成的闭域，即闭三角形，记作 $\Delta$，则 $\Delta\subset D$。今连接三角形各边中点，将 $\Delta$ 分解为四个三角形，如图 4b 所示，将其边界记作 $T_1$、$T_2$、$T_3$、$T_4$，则

$$
\int_T f(z)\,dz=\sum_{i=1}^4\int_{T_i}f(z)\,dz,
$$

而

$$
M=\left|\int_Tf(z)\,dz\right|
\le\sum_{i=1}^4\left|\int_{T_i}f(z)\,dz\right|,
$$

因此，右边各项中必有一项不小于 $M/4$，记该项所在的三角形为 $T^{(1)}$，所包围的闭域为 $\Delta_1$，则

$$
M_1\equiv\left|\int_{T^{(1)}}f(z)\,dz\right|\ge\frac{M}{4}.
$$

![Figure 4](...)

图 4: Cauchy 积分定理的证明

再将 $\Delta_1$ 分解为四个三角形，重复上面的论证，然后继续同样的过程，我们就可以得到一系列的三角形

$$
\Delta,\ \Delta_1,\ \Delta_2,\ \cdots,\ \Delta_n,\ \cdots,
$$

其边界为

$$
T,\ T^{(1)},\ T^{(2)},\ \cdots,\ T^{(n)},\ \cdots,
$$

在这些边界上的积分满足

$$
M_n\equiv\left|\int_{T^{(n)}}f(z)\,dz\right|
\ge\frac{M_{n-1}}{4}
\ge\frac{M}{4^n},\qquad n=1,2,\cdots, \tag{12}
$$

其中 $M_0\equiv M$。

记 $\Delta$ 的周长为 $L$，则 $\Delta_n$ 的周长为 $L_n=L/2^n$，显然，

$$
\Delta\supset\Delta_1\supset\Delta_2\supset\cdots\supset\Delta_n\supset\cdots,
$$

而 $L_n\to0$（$n\to\infty$），故存在唯一一点 $z_0$ 属于所有的 $\Delta_n$（这里用到了一个定理，称为闭集套定理，又称 Cantor 定理，但这一结论是相当直观的），当然 $z_0\in D$。今 $f(z)$ 在 $z_0$ 解析，则存在有限导数 $f'(z_0)$，成立

$$
\lim_{z\to z_0}\left[\frac{f(z)-f(z_0)}{z-z_0}-f'(z_0)\right]=0,
$$

换句话说，$\forall\varepsilon>0$，$\exists\delta>0$，使得当 $|z-z_0|<\delta$ 时，有

$$
\left|\frac{f(z)-f(z_0)}{z-z_0}-f'(z_0)\right|<\varepsilon,
$$

或

$$
|f(z)-f(z_0)-f'(z_0)(z-z_0)|<\varepsilon|z-z_0|.
$$

由于 $L_n\to0$（$n\to\infty$），故当 $n$ 足够大时，有 $L_n<\delta$，这时，对于 $T^{(n)}$ 上的任意点 $z$，都满足 $|z-z_0|<L_n<\delta$，故上面的不等式在 $T^{(n)}$ 上成立。由上节例 1 和例 2 知道，$\int_{T^{(n)}}f(z_0)\,dz=0$，$\int_{T^{(n)}}f'(z_0)(z-z_0)\,dz=0$，因此

$$
\begin{aligned}
M_n
&=\left|\int_{T^{(n)}}f(z)\,dz\right|
=\left|\int_{T^{(n)}}[f(z)-f(z_0)-f'(z_0)(z-z_0)]\,dz\right|\\
&\le\int_{T^{(n)}}|f(z)-f(z_0)-f'(z_0)(z-z_0)|\,|dz|\\
&<\int_{T^{(n)}}\varepsilon|z-z_0|\,|dz|
<\varepsilon L_n\int_{T^{(n)}}|dz|
=\varepsilon L_n^2
=\frac{L^2\varepsilon}{4^n}. \tag{13}
\end{aligned}
$$

比较式（12）和式（13），可得 $M/4^n\le M_n<L^2\varepsilon/4^n$，即 $M<L^2\varepsilon$，但 $\varepsilon$ 可以任意小，故 $M=0$。这样就完成了 Cauchy 积分定理的证明。

### §3.3　Cauchy 积分定理的推论

由 Cauchy 积分定理，马上可以得到以下

**推论**　设函数 $f(z)$ 在单通区域 $D$ 内解析，$z_1,z_2\in D$，则积分 $\int_{z_1}^{z_2}f(z)\,dz$ 与（$D$ 内）由 $z_1$ 到 $z_2$ 的路径无关。

**证明**　任取 $D$ 内由 $z_1$ 到 $z_2$ 的两条路径 $C_1$ 和 $C_2$，有

$$
\int_{C_1}f(z)\,dz-\int_{C_2}f(z)\,dz
=\int_{C_1}f(z)\,dz+\int_{C_2^-}f(z)\,dz
=\int_{C_1+C_2^-}f(z)\,dz=0,
$$

其中前两步都是用积分的基本性质，最后一步用了 Cauchy 积分定理，因为 $C_1+C_2^-$ 构成 $D$ 内的围线。证毕。

现在我们将积分的起点固定在 $z_0$，而让终点 $z$ 变化，这就是变上限积分。由上面的推论，这一变上限积分是 $z$ 的单值函数，类似于实变函数的情况，我们有下述结论。

**定理（变上限积分）**　设函数 $f(z)$ 在单通区域 $D$ 内解析，$z_0\in D$ 是定点，则由

$$
F(z)=\int_{z_0}^z f(\zeta)\,d\zeta
$$

定义的函数在 $D$ 内解析，且 $F'(z)=f(z)$。

**注**　在一元实变函数的相应定理中，只要求被积函数 $f(x)$ 连续，而不必可导。但这里必须要求 $f(z)$ 解析，否则积分与路径有关，就不是上限 $z$ 的单值函数了。另外，所考虑的区域必须是单通的，这一点也很重要。

**证明**　考虑 $z\in D$，任取 $z_0$ 到 $z$ 的路径，记作 $C_1$，又任取 $z_0$ 到 $z+\Delta z$ 的路径，记作 $C_2$，当然 $C_1$ 和 $C_2$ 都应在 $D$ 内。又记 $z$ 到 $z+\Delta z$ 的线段为 $P$，只要 $\Delta z$ 足够小，$P$ 也必在 $D$ 内。今

$$
F(z)=\int_{C_1}f(\zeta)\,d\zeta,\qquad
F(z+\Delta z)=\int_{C_2}f(\zeta)\,d\zeta.
$$

由上面的推论，$\int_{C_2}f(\zeta)\,d\zeta=\int_{C_1}f(\zeta)\,d\zeta+\int_Pf(\zeta)\,d\zeta$，则

$$
\Delta F\equiv F(z+\Delta z)-F(z)=\int_Pf(\zeta)\,d\zeta,
$$

由上节例 1，又有 $\int_P f(z)\,d\zeta=f(z)\Delta z$，故

$$
\frac{\Delta F}{\Delta z}-f(z)
=\frac1{\Delta z}\int_P[f(\zeta)-f(z)]\,d\zeta,
$$

而

$$
\left|\frac{\Delta F}{\Delta z}-f(z)\right|
=\frac1{|\Delta z|}\left|\int_P[f(\zeta)-f(z)]\,d\zeta\right|
\le\frac1{|\Delta z|}\int_P|f(\zeta)-f(z)|\,|d\zeta|.
$$

由于 $f(z)$ 在 $D$ 内解析，故必在 $D$ 内连续，所以，$\forall\varepsilon>0$，$\exists\delta>0$，使得当 $|\zeta-z|<\delta$ 时，有 $|f(\zeta)-f(z)|<\varepsilon$。今取 $|\Delta z|<\delta$，则线段 $P$ 上的各点均满足 $|f(\zeta)-f(z)|<\varepsilon$，于是

$$
\left|\frac{\Delta F}{\Delta z}-f(z)\right|
<\frac{\varepsilon}{|\Delta z|}\int_P|d\zeta|=\varepsilon,
$$

也就是说，

$$
\lim_{\Delta z\to0}\frac{\Delta F}{\Delta z}=f(z),\qquad z\in D.
$$

证毕。

类似于实变积分，下面引入原函数的概念。

**定义（原函数）**　设函数 $f(z)$ 在区域 $D$ 内解析，若有函数 $\Phi(z)$ 满足 $\Phi'(z)=f(z)$，则 $\Phi(z)$ 称为 $f(z)$ 的一个不定积分或原函数。

**注**　① 本来在这一定义中并没有必要要求 $f(z)$ 解析，但后面我们会看到，当 $\Phi(z)$ 解析（因为 $\Phi'(z)=f(z)$ 存在），则 $\Phi'(z)$ 亦必解析，从而 $f(z)=\Phi'(z)$ 也解析。换句话说，如果 $f(z)$ 不解析，它也不可能存在原函数。所以上述定义只针对解析函数就是很自然的了。② 原函数显然不是唯一的，因为如果 $\Phi(z)$ 是 $f(z)$ 的原函数，则 $\Phi(z)+c$ 也是它的原函数，其中 $c$ 是复常数。但任意两个原函数也只能相差一个常数，因为如果 $\Phi(z)$ 和 $\Psi(z)$ 都是 $f(z)$ 的原函数，则 $[\Psi(z)-\Phi(z)]'=0$，由此易证 $\Psi(z)-\Phi(z)=c$。

定义了原函数，就可以给出一个计算积分的公式，但仍然要注意它的条件。

**定理（Newton–Leibniz 公式）**　设函数 $f(z)$ 在单通区域 $D$ 内 textit解析，$\Phi(z)$ 是 $f(z)$ 的任一原函数，则

$$
\int_{z_0}^{z}f(\zeta)\,d\zeta
=\Phi(z)-\Phi(z_0),\qquad z,z_0\in D. \tag{14}
$$

**证明**　记 $F(z)=\int_{z_0}^z f(\zeta)\,d\zeta$，则 $F'(z)=f(z)$，又已知 $\Phi(z)$ 是 $f(z)$ 的原函数，则根据上面的讨论有 $\Phi(z)=F(z)+c$。以 $z_0$ 代入，可得 $\Phi(z_0)=F(z_0)+c=c$，所以 $F(z)=\Phi(z)-c=\Phi(z)-\Phi(z_0)$，此即式（14）。证毕。

**例 1**　利用这一公式，马上可以得到式（5）。

**例 2**　在不包含原点的单通区域内计算积分

$$
\int_a^b\frac1z\,dz,
$$

其中 $b\neq a$。

**解**　$1/z$ 在去掉原点的复平面上是解析的，但这是复通区域，在其中积分可能与路径有关。上节例 6 证实了这一点。这样，如果不指定路径，则积分就没有意义。但在不包含原点的单通区域内，比如右半平面，则积分与路径无关。$1/z$ 的原函数是 $\operatorname{Ln}z$，根据上面的公式（14），

$$
\int_a^b\frac1z\,dz=\operatorname{Ln}b-\operatorname{Ln}a.
$$

在所考虑的单通区域内，$\operatorname{Ln}z$ 可以分出单值分支，所以结果是确定的，而且与所取的分支无关。事实上，取定一个单值分支，其中 $a$ 和 $b$ 的辐角分别为 $\theta_a$ 和 $\theta_b$，则

$$
\int_a^b\frac1z\,dz
=\ln|b|-\ln|a|+i(\theta_b-\theta_a).
$$

若取定另一分支，其中 $a$ 和 $b$ 的辐角分别为 $\theta_a+2k\pi$ 和 $\theta_b+2k\pi$（$k\in\mathbb Z$），易见结果不变。

但是，对于不同的单通区域，结果却可能不同。比如当 $a=i$，$b=-i$，在单通区域 $0<\operatorname{Arg}z<2\pi$ 内（即去掉 $x$ 轴正半轴的复平面），$\theta_b-\theta_a=\pi$，积分结果为 $\pi i$；但在单通区域 $-\pi<\operatorname{Arg}z<\pi$ 内（即去掉 $x$ 轴负半轴的复平面），$\theta_b-\theta_a=-\pi$，积分结果为 $-\pi i$。

以上结果可以这样理解，在单通区域 $0<\operatorname{Arg}z<2\pi$ 内，沿各种可能路径的积分等于沿左半单位圆由 $i$ 到 $-i$ 的积分，记作 $I_1$，而在单通区域 $-\pi<\operatorname{Arg}z<\pi$ 内，沿各种可能路径的积分等于沿右半单位圆由 $i$ 到 $-i$ 的积分，记作 $I_2$，$I_1-I_2$ 等于沿单位圆正向的积分，由上节例 6，它等于 $2\pi i$，所以 $I_1\neq I_2$ 就是很自然的。而且，上面的计算也给出 $I_1-I_2=2\pi i$，与上节例 6 一致，如所期望。

### §3.4　复通区域的 Cauchy 积分定理

上节讨论 Cauchy 积分定理及其相关结论时，我们多次强调单通区域，因为在复通区域内，它一般是不能成立的。不过在复通区域 $D$ 内，只要围线 $C$ 及其内部完全在 $D$

![Figure 5](...)

图 5: 复通区域及其 Cauchy 积分定理的证明

内，它还是成立的。比如 $1/z$ 在去掉原点的复平面上是解析的，这是复通区域，在其中 Cauchy 积分定理不能成立，因为沿单位圆的积分就不为 0，但对于不包围原点的围线，比如右半平面上的围线，积分为 0。另一方面，我们知道，$1/z$ 在以原点为圆心的一切圆周上，其积分结果相同。这是复通区域的 Cauchy 积分定理的一个特例。本小节就是要讨论 Cauchy 积分定理在复通区域内的推广形式。

考虑 $n+1$ 条围线 $C_0,C_1,\cdots,C_n$，其中 $C_1,\cdots,C_n$ 全在 $C_0$ 内部，且互不相交也互不包含，如图 5a。在 $C_0$ 内部而在 $C_1,\cdots,C_n$ 外部的点集构成一个复通区域 $D$，它的边界包括以上各围线，称为复围线 $C$，即

$$
C=C_0+C_1^-+\cdots+C_n^-,
$$

其中 $C_0$ 取正向，而 $C_1,\cdots,C_n$ 取反向，以便沿边界绕行时，其所包围的区域 $D$ 总在左边。

**定理（复通区域的 Cauchy 积分定理）**　设 $D$ 是由复围线 $C=C_0+C_1^-+\cdots+C_n^-$ 围成的复通区域，函数 $f(z)$ 在 $\overline D$ 上解析，则

$$
\int_C f(z)\,dz
=\int_{C_0}f(z)\,dz+\int_{C_1^-}f(z)\,dz+\cdots+\int_{C_n^-}f(z)\,dz=0, \tag{15a}
$$

或

$$
\int_{C_0}f(z)\,dz
=\int_{C_1}f(z)\,dz+\cdots+\int_{C_n}f(z)\,dz. \tag{15b}
$$

**注**　定理的条件可以减弱为 $f(z)$ 在 $D$ 内解析，在 $\overline D$ 上连续。

**证明**　考虑只有两条围线 $C_0$ 和 $C_1$ 的情况。如图 5b，作线段 $L_1$ 和 $L_2$ 连接 $C_0$ 和 $C_1$，这样区域 $D$ 就被分为 $D'$ 和 $D''$（但注意 $D\neq D'\cup D''$），显然 $D'$ 和 $D''$ 都是单通区域，而且 $f(z)$ 在 $\overline{D'}$ 和 $\overline{D''}$ 上解析，根据 Cauchy 积分定理，有

$$
\int_{\partial D'}f(z)\,dz=0,\qquad
\int_{\partial D''}f(z)\,dz=0,
$$

所以

$$
\int_{\partial D'+\partial D''}f(z)\,dz=0,
$$

![Figure 6](...)

图 6: Cauchy 积分公式的证明

但 $\partial D'+\partial D''=C_0+C_1^-+L_1+L_2+L_1^-+L_2^-$，而 $\int_{L_1+L_1^-}f(z)\,dz=0$，$\int_{L_2+L_2^-}f(z)\,dz=0$，所以

$$
\int_C f(z)\,dz
=\int_{C_0}f(z)\,dz+\int_{C_1^-}f(z)\,dz=0,
$$

此即式（15a），由积分的基本性质易得式（15b）。对于有多条围线的情况，可以类似证明。证毕。

以后我们说到 Cauchy 积分定理，就包括复通区域的情况。比较式（15a）与 Cauchy 积分定理的等价表述，可以看出，只要积分路径包括全部的边界，那么复通与单通情况下的 Cauchy 积分定理在形式上并没有什么区别。

**例 3**　由上面的定理和前面的结果（7）可知

$$
\int_C\frac{dz}{(z-a)^n}
=\begin{cases}
2\pi i,&(n=1),\\
0,&(n\in\mathbb Z,\ n\neq1),
\end{cases}
$$

其中 $C$ 是包围 $a$ 点的任一围线。

如果围线 $C$ 不包围 $a$ 点，则积分显然为 0。于是，当 $n\neq1$ 时，无论围线 $C$ 是否包围 $a$ 点，积分都是 0。但是必须注意，对于 $n>1$，围线不能经过 $a$ 点，否则积分不存在。

## §4　Cauchy 积分公式及其推论

### §4.1　Cauchy 积分公式

由 Cauchy 积分定理，可以推出下面的 Cauchy 积分公式，我们把它写成定理。

**定理（Cauchy 积分公式）**　设区域 $D$ 的边界是围线或复围线 $C$，函数 $f(z)$ 在闭域 $\overline D$ 上解析，则

$$
f(z)=\frac{1}{2\pi i}\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta,
\qquad z\in D. \tag{16}
$$

**注**　① Cauchy 积分公式表明，对于解析函数，只要边界上的函数值给定，则区域内的函数值也就完全确定了。这说明函数值在各点的分布是互相牵制、紧密关联的，其实 Cauchy 积分定理也表明了这种关联。而上一章的 CR 条件则表明解析函数的实部和虚部也互相牵制、紧密关联，给定了其中一个，另一个也就确定了（最多可以差一常数项）。可见解析性对于复变函数是一个很强的限制。在实变函数中，没有任何类似的结论，无论要求函数多么光滑，其变化都还是可以相当任意的，区间端点的函数值完全不能决定区间内部的函数值。② 将公式改写为

$$
\int_C\frac{f(z)}{z-a}\,dz=2\pi i f(a),\qquad a\in D, \tag{16}
$$

则可以用来计算某些积分。③ 定理的条件可以减弱为 $f(z)$ 在 $D$ 内解析，在 $\overline D$ 上连续。④ 由 Cauchy 积分定理可以推出 Cauchy 积分公式（见下面的证明），反过来，由 Cauchy 积分公式也可以推出 Cauchy 积分定理，所以两者是等价的。事实上，设 $F(z)$ 是 $\overline D$ 上的任意解析函数，则 $(z-a)F(z)$ 也是 $\overline D$ 上的解析函数，根据 Cauchy 积分公式，就有

$$
\int_C F(z)\,dz
=\int_C\frac{[(z-a)F(z)]}{z-a}\,dz
=2\pi i[(z-a)F(z)]_{z=a}=0,
$$

这就是 Cauchy 积分定理。

**证明**　设 $D$ 的边界是复围线 $C=C_0+C_1^-+\cdots+C_n^-$。$\forall z\in D$，以 $z$ 为中心，$\rho$ 为半径作圆周 $\Gamma_\rho:|\zeta-z|=\rho$，使 $\Gamma_\rho$ 在 $C_0$ 的内部，而在 $C_1,\cdots,C_n$ 的外部。今 $f(\zeta)/(\zeta-z)$ 在 $C+\Gamma_\rho^-$ 所围成的区域（即区域 $D$ 挖去闭圆 $|\zeta-z|\le\rho$ 后所剩余的点集）及其边界上解析。由 Cauchy 积分定理和积分的基本性质，有

$$
0=\int_{C+\Gamma_\rho^-}\frac{f(\zeta)}{\zeta-z}\,d\zeta
=\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta
+\int_{\Gamma_\rho^-}\frac{f(\zeta)}{\zeta-z}\,d\zeta
=\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta
-\int_{\Gamma_\rho}\frac{f(\zeta)}{\zeta-z}\,d\zeta,
$$

所以

$$
\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta
=\int_{\Gamma_\rho}\frac{f(\zeta)}{\zeta-z}\,d\zeta.
$$

根据式（7），有

$$
\int_{\Gamma_\rho}\frac{f(z)}{\zeta-z}\,d\zeta
=f(z)\int_{\Gamma_\rho}\frac1{\zeta-z}\,d\zeta
=2\pi i f(z).
$$

结合两式，可得

$$
\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta-2\pi i f(z)
=\int_{\Gamma_\rho}\frac{f(\zeta)-f(z)}{\zeta-z}\,d\zeta,
$$

于是

$$
\left|\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta-2\pi i f(z)\right|
\le\frac1\rho\int_{\Gamma_\rho}|f(\zeta)-f(z)|\,|d\zeta|.
$$

由于 $f(\zeta)$ 在 $D$ 内解析，故在 $\zeta=z$ 处连续，故 $\forall\varepsilon>0$，$\exists\delta>0$，当 $|\zeta-z|<\delta$ 时，有 $|f(\zeta)-f(z)|<\varepsilon/2\pi$。取 $\rho<\delta$，则上式在 $\Gamma_\rho$ 上成立，于是

$$
\left|\int_C\frac{f(\zeta)}{\zeta-z}\,d\zeta-2\pi i f(z)\right|
<\frac{\varepsilon}{2\pi\rho}\int_{\Gamma_\rho}|d\zeta|=\varepsilon.
$$

这就是说，只要 $\rho$ 足够小，上式就成立。但上式左边实际上与 $\rho$ 无关，所以它必须为 0，如此即得式（16）。证毕。

**例 1**　计算积分

$$
I=\int_C\frac{dz}{z^2-1},
$$

其中围线 $C$ 是：（1）$|z|=1/2$；（2）$|z-1|=1/2$；（3）$|z+1|=1/2$；（4）$|z|=2$。

**解**　（1）由于被积函数在围线 $C$ 及其内部解析，故 $I=0$。

（2）由式（16），

$$
I=\int_{|z-1|=1/2}\frac{dz}{z^2-1}
=\int_{|z-1|=1/2}\frac{1/(z+1)}{z-1}\,dz
=2\pi i\left.\frac1{z+1}\right|_{z=1}
=\pi i.
$$

（3）由式（16），

$$
I=\int_{|z+1|=1/2}\frac{dz}{z^2-1}
=\int_{|z+1|=1/2}\frac{1/(z-1)}{z+1}\,dz
=2\pi i\left.\frac1{z-1}\right|_{z=-1}
=-\pi i.
$$

（4）由复通区域的 Cauchy 积分定理和（2）、（3）的结果，

$$
I=\int_{|z|=2}\frac{dz}{z^2-1}
=\int_{|z-1|=1/2}\frac{dz}{z^2-1}
+\int_{|z+1|=1/2}\frac{dz}{z^2-1}
=\pi i-\pi i=0.
$$

### §4.2　\*无界区域的 Cauchy 积分公式

既然解析函数在围线上的函数值可以决定其内部的函数值，人们自然会问，围线上的函数值是否也可以决定其外部的函数值？答案就在下面的定理中。

**定理（无界区域的 Cauchy 积分公式）**　设函数 $f(z)$ 在围线 $C$ 及其外部的无界区域 $D$ 上解析，且当 $z\to\infty$ 时，$f(z)\Rightarrow0$，则有

$$
f(z)=\frac{1}{2\pi i}\int_{C^-}\frac{f(\zeta)}{\zeta-z}\,d\zeta,
\qquad z\in D, \tag{17}
$$

其中积分取 $C$ 的反方向，实即 $\partial D$ 的正方向。

**注**　① $z\to\infty$ 时 $f(z)\Rightarrow0$（一致趋于 0）的大意是 $f(z)\to0$ 的速度与 $z$ 的辐角无关。精确地说，就是 $\forall\varepsilon>0$，$\exists R_0>0$ 与 $\theta$ 无关，当 $|z|>R_0$ 时，就有 $|f(z)|<\varepsilon$。② 这一定理的条件可以减弱为 $f(z)$ 在 $D$ 内解析，在 $D+C$ 上连续。③ 如果 $f(z)$ 在围线 $C$ 的内部也解析，则式（17）中的被积函数在围线 $C$ 及其内部解析，根据 Cauchy 积分定理，积分为 0。换句话说，在 $D$ 上，$f(z)=0$。由 $f(z)$ 的连续性，在 $C$ 上也有 $f(z)=0$。再根据 Cauchy 积分公式，则在 $C$ 的内部也有 $f(z)=0$。这就是说，在整个 $z$ 平面上，$f(z)\equiv0$。这是什么原因呢？根据假定，$f(z)$ 在整个 $z$ 平面上解析（称为整函数，entire function），又因为 $z\to\infty$ 时，$f(z)\Rightarrow0$，所以它一定是有界的（$\forall\varepsilon>0$，$\exists R_0>0$，当 $|z|>R_0$ 时，$|f(z)|<\varepsilon$，而在闭圆 $|z|\le R_0$ 上，$\exists M_0>0$，使 $|f(z)|\le M_0$，取 $M=M_0+\varepsilon$，则在整个 $z$ 平面上 $|f(z)|<M$，即有界），根据 Liouville 定理，有界整函数必为常数（见后），而条件 $z\to\infty$ 时，$f(z)\Rightarrow0$ 使得该常数为 0。④ 该定理可以推广到多条围线外的情况，只要将式（17）右边看作各围线的积分之和即可。

**证明**　作大圆 $\Gamma_R:|\zeta|=R$，使得围线 $C$ 和点 $z$ 均在 $\Gamma_R$ 内部，根据 Cauchy 积分公式，

$$
f(z)
=\frac1{2\pi i}\int_{C^-}\frac{f(\zeta)}{\zeta-z}\,d\zeta
+\frac1{2\pi i}\int_{\Gamma_R}\frac{f(\zeta)}{\zeta-z}\,d\zeta, \tag{18}
$$

由此

$$
\int_{\Gamma_R}\frac{f(\zeta)}{\zeta-z}\,d\zeta
=2\pi i f(z)-\int_{C^-}\frac{f(\zeta)}{\zeta-z}\,d\zeta. \tag{19}
$$

由定理条件，$\forall\varepsilon>0$，$\exists R_0>0$，当 $R>R_0$ 时，在 $\Gamma_R$ 上可有 $|f(\zeta)|<\varepsilon/4\pi$，今取 $R>\max\{R_0,2|z|\}$，则在 $\Gamma_R$ 上还成立 $|\zeta-z|>R-|z|>R/2$，于是 $1/|\zeta-z|<2/R$，故

$$
\left|\int_{\Gamma_R}\frac{f(\zeta)}{\zeta-z}\,d\zeta\right|
\le\int_{\Gamma_R}\frac{|f(\zeta)|}{|\zeta-z|}\,|d\zeta|
<\frac{\varepsilon}{4\pi R}\,2\int_{\Gamma_R}|d\zeta|
=\varepsilon.
$$

但由式（19）可知，上式左边实际上与 $R$ 无关，所以它必须为 0，代入式（18），即得式（17）。证毕。

### §4.3　解析函数的高阶导数

我们在前面曾提到，解析函数存在各阶导数，即一次可微导致任意次可微，这是复变函数所特有的结论。而且，边界上的函数值不仅确定了所围区域内的函数值，也确定了其中各阶导数的函数值。下面关于高阶导数的定理给出了这一结论的精确表述。

**定理（Cauchy 高阶导数公式）**　设区域 $D$ 以围线或复围线 $C$ 为边界，函数 $f(z)$ 在闭域 $\overline D$ 上解析，则 $f(z)$ 在区域 $D$ 内有各阶导数，且

$$
f^{(n)}(z)=\frac{n!}{2\pi i}\int_C\frac{f(\zeta)}{(\zeta-z)^{n+1}}\,d\zeta,
\qquad z\in D,\quad n=1,2,\cdots. \tag{20}
$$

**注**　① 将 Cauchy 积分公式（16）两边求导，对右边交换求导与积分的次序，立得一阶导数的 Cauchy 公式。继续求导，重复同样的操作，即得 Cauchy 高阶导数公式。这样的做法显然是不严格的，因为求导与积分交换次序的合法性并未得到证明。然而，这一做法能帮助我们熟悉高阶导数公式，并使得我们能够在记得 Cauchy 积分公式的情况下立即将高阶导数公式“推导”出来。② 让我们再看看怎样可以做得严格一些。以 $n=1$ 为例，对 $z$ 和 $z+\Delta z$ 分别用 Cauchy 积分公式，可得

$$
\frac{\Delta f}{\Delta z}
=\frac{f(z+\Delta z)-f(z)}{\Delta z}
=\frac1{2\pi i}\int_C\frac{f(\zeta)}{(\zeta-z)(\zeta-z-\Delta z)}\,d\zeta. \tag{21}
$$

两边取 $\Delta z\to0$ 的极限，对右边交换求极限与积分的次序，立得一阶导数的 Cauchy 公式。类似可得高阶导数公式。必须指出，这样的做法仍然是不严格的，因为求极限与积分交换次序的合法性也未得到证明。不过，我们已经向严格证明的方向迈进了重要的一步。③ 类似于 Cauchy 积分公式，Cauchy 高阶导数公式也可以用来计算某些积分。④ 定理的条件可以减弱为 $f(z)$ 在区域 $D$ 上解析，在闭域 $\overline D$ 上连续。

**证明**　以下只证明 $n=1$ 的情况，用数学归纳法可以证明一般情况。由式（21），$\forall z\in D$，有

$$
\left|\frac{\Delta f}{\Delta z}-\frac1{2\pi i}\int_C\frac{f(\zeta)}{(\zeta-z)^2}\,d\zeta\right|
=\frac{|\Delta z|}{2\pi}\left|\int_C\frac{f(\zeta)}{(\zeta-z)^2(\zeta-z-\Delta z)}\,d\zeta\right|.
$$

我们的思路是证明左边当 $\Delta z\to0$ 时可以任意小，这样它的极限就必为 0。由于右边含有因子 $|\Delta z|$，故主要需证明右边的积分有界。由于 $f(z)$ 在闭域 $\overline D$ 上连续，故在 $C$ 上有界，设 $\max_{\zeta\in C}|f(\zeta)|=M$（如此则 $|f(\zeta)|\le M$，$\forall\zeta\in C$）。又设 $z$ 与边界 $C$ 的距离为 $d$，即 $\min_{\zeta\in C}|\zeta-z|=d$（如此则 $|\zeta-z|\ge d$，$\forall\zeta\in C$）。今暂取 $|\Delta z|<d/2$，则 $|\zeta-z-\Delta z|\ge|\zeta-z|-|\Delta z|>d-d/2=d/2$。于是

$$
\left|\int_C\frac{f(\zeta)}{(\zeta-z)^2(\zeta-z-\Delta z)}\,d\zeta\right|
\le\int_C\frac{|f(\zeta)|}{|\zeta-z|^2|\zeta-z-\Delta z|}\,|d\zeta|
\le\frac{M}{d^3/2}\int_C|d\zeta|
=\frac{2ML}{d^3},
$$

其中 $L$ 是 $C$ 的长度，而

$$
\left|\frac{\Delta f}{\Delta z}-\frac1{2\pi i}\int_C\frac{f(\zeta)}{(\zeta-z)^2}\,d\zeta\right|
\le\frac{ML}{\pi d^3}|\Delta z|,
$$

上式当 $|\Delta z|<d/2$ 时成立，若同时又有 $|\Delta z|<\varepsilon\pi d^3/ML$，则上式左边 $<\varepsilon$。取 $|\Delta z|<\min\{d/2,\varepsilon\pi d^3/ML\}$，则上式左边 $<\varepsilon$，$\forall\varepsilon>0$，于是

$$
f'(z)=\lim_{\Delta z\to0}\frac{\Delta f}{\Delta z}
=\frac1{2\pi i}\int_C\frac{f(\zeta)}{(\zeta-z)^2}\,d\zeta.
$$

证毕。

### §4.4　\*Liouville 定理

Cauchy 积分公式和高阶导数公式有许多有趣的推论，我们只介绍其中的一个，即 Liouville 定理，在下一小节中，我们要用它来证明代数基本定理。我们在前面曾提到整函数（entire function），下面正式给出定义。

**定义（整函数）**　在整个 $z$ 平面上解析的函数称为整函数。

**例**　多项式、指数函数、正弦、余弦等函数都是整函数。

**定理（Liouville）**　有界整函数 $f(z)$ 必为常数。

**注**　这一定理表明，一个解析函数，要么有奇点，要么当 $z\to\infty$ 时，至少在某些方向是无限的，除非是常数。因此，不存在什么处处解析、处处有限而非平庸的“理想”解析函数。这使得我们对于解析函数的函数值分布的特性有了一个大致的图象。

**证明**　由于 $f(z)$ 有界，故 $\exists M>0$，$\forall z\in\mathbb C$，有 $|f(z)|\le M$。今任取固定点 $z$，作圆周 $\Gamma_R:|\zeta-z|=R$。由高阶导数公式

$$
f'(z)=\frac1{2\pi i}\int_{\Gamma_R}\frac{f(\zeta)}{(\zeta-z)^2}\,d\zeta,
$$

故

$$
|f'(z)|
\le\frac1{2\pi}\int_{\Gamma_R}\frac{|f(\zeta)|}{R^2}|d\zeta|
\le\frac1{2\pi}\frac{M}{R^2}\int_{\Gamma_R}|d\zeta|
=\frac1{2\pi}\frac{M}{R^2}2\pi R
=\frac{M}{R}.
$$

今 $\forall\varepsilon>0$，只要 $R>M/\varepsilon$，就有 $|f'(z)|<\varepsilon$，但 $f'(z)$ 实际上与 $R$ 无关，故 $f'(z)=0$，又因为 $z$ 是任意的，故 $f(z)$ 为常数。证毕。

### §4.5　\*代数基本定理

在复变函数理论中，可以很简单地证明代数基本定理，这是复变函数理论的应用之一。

**定理（代数基本定理）**　在 $z$ 平面上，$n$ 次多项式 $P_n(z)=\sum_{k=0}^n a_kz^k$ 至少有一个零点。

**证明**　设 $P_n(z)$ 在 $z$ 平面上没有零点，则 $f(z)=1/P_n(z)$ 在 $z$ 平面上解析，即为整函数。又

$$
\lim_{z\to\infty}f(z)
=\lim_{z\to\infty}\frac1{P_n(z)}
=\lim_{z\to\infty}\frac1{\sum_{k=0}^n a_kz^k}
=\lim_{z\to\infty}\frac1{z^n\sum_{k=0}^n a_kz^{-(n-k)}}
=0,
$$

所以 $f(z)$ 有界。事实上，由上面的极限可知，$\exists R>0$，当 $|z|>R$ 时，有 $|f(z)|<1$，而由 $f(z)$ 的连续性，在 $|z|\le R$ 上，$\exists M>0$，使 $|f(z)|\le M$，所以在 $z$ 平面上，$|f(z)|<M+1$，即有界。由 Liouville 定理，在 $z$ 平面上，$f(z)$ 为常数，从而 $P_n(z)$ 为常数，这显然是不对的。证毕。

**习题**　计算下列积分

1.

$$
\int_{|z|=2}\frac{1}{z^4-1}\,dz.
$$

2.

$$
\int_{|z|=3}\frac{e^z}{z^2(z-2)^2}\,dz.
$$

## 补充习题

计算下列积分。

1.

$$
\int_{|z|=1}\frac{\cosh z}{z^{n+1}}\,dz,\qquad n=1,2,\cdots.
$$

答案：$\pi i[1+(-)^n]/n!$。

2.

$$
\int_{|z|=1}\frac{\sin^2 z}{z^4}\,dz.
$$

答案：$0$。

3.

$$
\int_{|z|=1}\left(z+\frac1z\right)^{2n}\frac{dz}{z},\qquad n=1,2,\cdots.
$$

答案：$2\pi i C_{2n}^{n}$。

4.

$$
\int_{|z|=2}\frac{\sin(e^z)}{z}\,dz.
$$

答案：$2\pi i\sin1$。

5.

$$
\int_{|z|=2}\frac{e^z}{\cosh z}\,dz.
$$

提示：$\pm i\pi/2$ 是奇点，将原积分化为 $z\mp i\pi/2=\delta e^{i\theta}$ 上的积分之和，其中 $\delta$ 很小，然后作变换 $\zeta=e^z$。第五章补充习题 10 可与本题互相印证。答案：$4\pi i$。

6. （1）

$$
\int_{|z|=1}\frac{dz}{z},
$$

（2）

$$
\int_{|z|=1}\frac{|dz|}{z},
$$

（3）

$$
\int_{|z|=1}\frac{dz}{|z|},
$$

（4）

$$
\int_{|z|=1}\left|\frac{dz}{z}\right|.
$$
