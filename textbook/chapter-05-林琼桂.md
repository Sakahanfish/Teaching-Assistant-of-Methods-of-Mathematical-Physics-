# 第五章 留数定理及其应用*

## §1 留数定理

### §1.1 留数的定义

如果函数 $f(z)$ 在点 $a$ 的邻域 $K:|z-a|<R$ 内解析，围线 $C$ 全含于 $K$（包围 $a$ 或不包围 $a$），则根据 Cauchy 积分定理必有

$$
\oint_C f(z)\,\mathrm dz=0.
$$

但如果 $a$ 是 $f(z)$ 的孤立奇点，即 $f(z)$ 只在点 $a$ 的去心邻域 $K\setminus\{a\}:0<|z-a|<R$ 内解析，围线 $C$ 是 $K\setminus\{a\}$ 中包围 $a$ 的围线，则上式不一定成立。故定义留数如下：

**定义（留数）** 如果函数 $f(z)$ 以 $a$ 为孤立奇点，即 $f(z)$ 在 $K\setminus\{a\}:0<|z-a|<R$ 中解析，则积分

$$
\frac{1}{2\pi i}\int_{\Gamma_\rho} f(z)\,\mathrm dz,
\qquad \Gamma_\rho:|z-a|=\rho\quad(0<\rho<R)
\tag{1}
$$

称为 $f(z)$ 在点 $a$ 处的留数（residue），记作 $\operatorname{Res}_{z=a}f(z)$，简记为 $\operatorname{Res}f(a)$，或 $\operatorname{Res}(f,a)$。

显然，只要 $0<\rho<R$，上述积分的数值与 $\rho$ 的大小无关。

设 $f(z)$ 在 $a$ 处的 Laurent 展开式为

$$
f(z)=\sum_{n=-\infty}^{+\infty}c_n(z-a)^n,
$$

则

$$
c_n=\frac{1}{2\pi i}\int_{\Gamma_\rho}\frac{f(z)}{(z-a)^{n+1}}\,\mathrm dz.
$$

令 $n=-1$，得

$$
c_{-1}=\frac{1}{2\pi i}\int_{\Gamma_\rho}f(z)\,\mathrm dz.
$$

与式 (1) 比较，即得

$$
\operatorname{Res}f(a)=c_{-1}.
\tag{2}
$$

由此可知，可去奇点处的留数为 0。

**注** 有些书上直接用式 (2) 作为留数的定义，这与式 (1) 的定义显然是等价的。

### §1.2 Cauchy 留数定理

利用 Cauchy 积分定理，可以推出下面关于围线积分的

**Cauchy 留数定理** 设函数 $f(z)$ 在围线或复围线 $C$ 所围成的区域 $D$ 中有孤立奇点 $a_1,a_2,\cdots,a_n$，此外 $f(z)$ 在 $\overline D$ 上解析，则有

$$
\int_C f(z)\,\mathrm dz
=2\pi i\sum_{k=1}^{n}\operatorname{Res}f(a_k).
\tag{3}
$$

![Figure 1](...)

图 1：区域 $D$ 的边界为复围线 $C$，其中有 $n$ 个孤立奇点

**证明** 作 $n$ 个小圆周 $\Gamma_k:|z-a_k|=\rho_k$（$k=1,2,\cdots,n$），使各 $\Gamma_k$ 及其内部全含于 $D$，但各 $\Gamma_k$ 互不相交也互不包含，如图 1。根据 Cauchy 积分定理和留数的定义，有

$$
\int_C f(z)\,\mathrm dz
=\sum_{k=1}^{n}\int_{\Gamma_k}f(z)\,\mathrm dz
=\sum_{k=1}^{n}2\pi i\operatorname{Res}f(a_k)
=2\pi i\sum_{k=1}^{n}\operatorname{Res}f(a_k).
$$

证毕。

利用 Cauchy 留数定理，只要算出各孤立奇点处的留数，即可得出围线积分，所以关键在于计算留数。

## §2 留数的算法

计算留数的最一般方法是作 Laurent 展开，求出系数 $c_{-1}$，即得展开中心的留数。但是作 Laurent 展开往往太麻烦，所以我们希望有一些现成的公式可以用来计算留数。下面的方法适用于计算极点的留数。

**定理（极点的留数）** 设 $a$ 是 $f(z)$ 的 $n$ 阶极点，则

$$
\operatorname{Res}f(a)
=\frac{1}{(n-1)!}
\left\{
\frac{\mathrm d^{\,n-1}}{\mathrm dz^{n-1}}
\left[(z-a)^n f(z)\right]
\right\}_{z\to a}.
\tag{4}
$$

**证明** 由于 $a$ 是 $f(z)$ 的 $n$ 阶极点，故在 $a$ 的某去心邻域内有

$$
f(z)=\frac{\varphi(z)}{(z-a)^n},
$$

其中 $\varphi(z)$ 在 $a$ 点解析，且 $\varphi(a)\ne0$。于是

$$
\operatorname{Res}f(a)
=\frac{1}{2\pi i}\int_{\Gamma_\rho}f(z)\,\mathrm dz
=\frac{1}{2\pi i}\int_{\Gamma_\rho}\frac{\varphi(z)}{(z-a)^n}\,\mathrm dz
=\frac{1}{(n-1)!}\varphi^{(n-1)}(a),
$$

最后一个等号用了高阶导数公式。由于 $\varphi(z)=(z-a)^n f(z)$，代入上式即得式 (4)。证毕。

由这一定理，马上可以得到下面两个推论，它们是常用的计算公式。

**推论一（单极点的留数，第一公式）** 若 $a$ 是 $f(z)$ 的单极点，则

$$
\operatorname{Res}f(a)
=\left[(z-a)f(z)\right]_{z\to a}.
\tag{5}
$$

**推论二（二阶极点的留数）** 若 $a$ 是 $f(z)$ 的二阶极点，则

$$
\operatorname{Res}f(a)
=\left[(z-a)^2f(z)\right]'_{z\to a}.
\tag{6}
$$

单极点的留数还有另一个计算公式，通常它比式 (5) 更方便，我们也把它写成定理的形式。

**定理（单极点的留数，第二公式）** 设 $a$ 是 $f(z)=\varphi(z)/\psi(z)$ 的单极点，即 $\varphi(z)$ 和 $\psi(z)$ 均在 $a$ 点解析，且 $\varphi(a)\ne0$，$\psi(a)=0$，$\psi'(a)\ne0$，则

$$
\operatorname{Res}f(a)=\frac{\varphi(a)}{\psi'(a)}.
\tag{7}
$$

**证明** 按式 (5)，并利用 $\psi(a)=0$ 等条件，有

$$
\begin{aligned}
\operatorname{Res}f(a)
&=\lim_{z\to a}\left[(z-a)f(z)\right]
=\lim_{z\to a}\left[(z-a)\frac{\varphi(z)}{\psi(z)}\right]\\
&=\lim_{z\to a}\left[(z-a)\frac{\varphi(z)}{\psi(z)-\psi(a)}\right]
=\lim_{z\to a}\frac{\varphi(z)}{[\psi(z)-\psi(a)]/(z-a)}
=\frac{\varphi(a)}{\psi'(a)}.
\end{aligned}
$$

证毕。

以上公式 (5)～(7) 就是计算极点留数的常用公式。当极点的阶数较高时，用公式 (4) 计算留数并不方便，这时不如作 Laurent 展开，求出系数 $c_{-1}$。对于本性奇点，这也是唯一的计算方法。至于可去奇点，前已指出，其留数为 0。

**注** 作为计算留数的手段，作 Laurent 展开时只需求出系数 $c_{-1}$，其它系数可以不必理会，所以只要把包含 $z^{-1}$ 的项找出来就可以了。§4 有具体的例子。

## §3 用留数定理计算围线积分

下面给出几个用留数定理计算围线积分的例子。

**例 1**

$$
I=\int_{|z|=2}\frac{1}{z(z-1)^2}\,\mathrm dz.
$$

**解** 本题的被积函数 $f(z)=1/[z(z-1)^2]$ 在圆周 $|z|=2$ 的内部有一阶极点 $z=0$ 和二阶极点 $z=1$。按上节的留数计算公式，易得

$$
\operatorname{Res}f(0)
=\left.\frac{1}{(z-1)^2}\right|_{z=0}=1,
\qquad
\operatorname{Res}f(1)
=\left.\left(\frac1z\right)'\right|_{z=1}=-1,
$$

故

$$
I=2\pi i[\operatorname{Res}f(0)+\operatorname{Res}f(1)]=0.
$$

**注** 如果没有留数定理，我们可以作两个小圆周分别包围一阶极点 $z=0$ 和二阶极点 $z=1$，根据 Cauchy 积分定理，原积分等于沿两个小圆周的积分之和，后两个积分可以分别用 Cauchy 积分公式和 Cauchy 高阶导数公式来计算。其实这些做法正是证明留数定理时所采用的思路。

**例 2**

$$
I=\int_{|z|=1}\frac{\sin z}{z}\,\mathrm dz.
$$

**解** 本题的被积函数 $f(z)=\sin z/z$ 在圆周 $|z|=1$ 的内部有一个可去奇点 $z=0$，它在 $z=0$ 处的 Laurent 展开式没有负幂项，故 $\operatorname{Res}f(0)=c_{-1}=0$，而 $I=2\pi i\operatorname{Res}f(0)=0$。

**注** 本题当然也可以用 Cauchy 积分公式来计算，被积函数已经具有公式所要求的形式，故 $I=2\pi i\sin0=0$。

**例 3**

$$
I=\int_{|z|=1}e^{1/z}\,\mathrm dz.
$$

**解** 本题的被积函数 $f(z)=e^{1/z}$ 在圆周 $|z|=1$ 的内部有一个本性奇点 $z=0$，它在 $z=0$ 处的 Laurent 展开式为

$$
f(z)=e^{1/z}=1+\frac1z+\cdots+\frac{1}{n!z^n}+\cdots,
$$

故 $\operatorname{Res}f(0)=c_{-1}=1$，而

$$
I=2\pi i\operatorname{Res}f(0)=2\pi i.
$$

**注** 本题不能直接用 Cauchy 积分公式来计算，但可以对被积函数作 Laurent 展开后逐项积分，只有 $1/z$ 项有贡献，结果当然一样。作 Laurent 展开后逐项积分，这也是与留数定理的精神一致的。

以前，我们曾经用 Cauchy 积分公式和 Cauchy 高阶导数公式计算某些围线积分，现在可以直接应用留数定理来计算。作为计算围线积分的技术，留数定理包括了 Cauchy 积分公式和 Cauchy 高阶导数公式作为特殊情况，所以是更一般的方法。当然，在逻辑上是先有 Cauchy 积分公式和 Cauchy 高阶导数公式的。

另外，我们看到，留数定理把计算围线积分的问题转化为计算围线内各孤立奇点的留数的问题。由上节可以看到，计算极点的留数主要涉及微分运算；对于本性奇点，必须作 Laurent 展开来计算其留数。作 Laurent 展开，通常归结为 Taylor 展开，而计算 Taylor 展开式的系数也是微分运算问题。所以可以说，留数定理把积分运算转化成了比较容易的微分运算，因此它为积分的计算提供了一项非常有用的技术。

**习题** 计算下列围线积分。

(1)
$$
\int_C\frac{\mathrm dz}{(z^2+1)(z-1)^2},
$$
其中 $C$ 的方程是 $x^2+y^2-2x-2y=0$；

(2)
$$
\int_{|z|=1}\frac{\cos z}{z^3}\,\mathrm dz;
$$

(3)
$$
\int_{|z|=2}\frac{2z}{1-2\sin^2 z}\,\mathrm dz.
$$

## §4 *无穷远点的留数

当 $\infty$ 点是函数 $f(z)$ 的孤立奇点时，可以定义其留数如下。

**定义（$\infty$ 的留数）** 如果函数 $f(z)$ 以 $\infty$ 为孤立奇点，即 $f(z)$ 在 $\infty$ 的去心邻域 $r<|z|<\infty$ 中解析（$r$ 为正数），则积分

$$
\frac{1}{2\pi i}\int_{\Gamma^-}f(z)\,\mathrm dz,
\qquad \Gamma:|z|=\rho\quad(r<\rho<\infty)
\tag{8}
$$

称为 $f(z)$ 在 $\infty$ 点的留数，记作 $\operatorname{Res}_{z=\infty}f(z)$，或简记为 $\operatorname{Res}f(\infty)$，或 $\operatorname{Res}(f,\infty)$。

**注** 以上积分取 $\Gamma$ 的负方向，这正是包围 $\infty$ 点的正方向。显然，只要 $r<\rho<\infty$，上述积分的数值与 $\rho$ 的大小无关。

设 $f(z)$ 在 $r<|z|<\infty$ 内的 Laurent 展开式为

$$
f(z)=\sum_{n=-\infty}^{+\infty}c_nz^n,
$$

则

$$
c_n=\frac{1}{2\pi i}\int_\Gamma\frac{f(z)}{z^{n+1}}\,\mathrm dz.
$$

令 $n=-1$，得

$$
c_{-1}=\frac{1}{2\pi i}\int_\Gamma f(z)\,\mathrm dz.
$$

与式 (8) 比较，即得

$$
\operatorname{Res}f(\infty)=-c_{-1}.
\tag{9}
$$

如果函数 $f(z)$ 在 $z$ 平面上除有限个孤立奇点之外没有其它有限的奇点，那么 $\infty$ 点当然也是孤立奇点。这时候，它们的留数满足一个简单的关系，这就是下面的

**定理（总留数）** 设函数 $f(z)$ 在扩充 $z$ 平面上只有有限个孤立奇点 $a_1,a_2,\cdots,a_n$ 和 $\infty$，则 $f(z)$ 在各点的留数之和为 0，即

$$
\sum_{k=1}^{n}\operatorname{Res}f(a_k)+\operatorname{Res}f(\infty)=0.
\tag{10}
$$

**证明** 作大圆周 $\Gamma:|z|=\rho$ 将 $a_1,a_2,\cdots,a_n$ 全部包围在内，根据 Cauchy 留数定理，有

$$
\int_\Gamma f(z)\,\mathrm dz
=2\pi i\sum_{k=1}^{n}\operatorname{Res}f(a_k).
$$

另一方面，由于 $f(z)$ 在 $\Gamma$ 及其外部解析，由定义 (8)，有

$$
\int_{\Gamma^-}f(z)\,\mathrm dz=2\pi i\operatorname{Res}f(\infty).
$$

两式相加，注意到 $\int_{\Gamma^-}f(z)\,\mathrm dz=-\int_\Gamma f(z)\,\mathrm dz$，立得式 (10)。证毕。

这一定理的用处可以通过下面的例子显示出来。

**例** 计算积分

$$
I=\int_{|z|=2}\frac{z^{15}}{(z^4+1)^4}\,\mathrm dz.
$$

**解** 本题的被积函数

$$
f(z)=\frac{z^{15}}{(z^4+1)^4}
$$

在圆周 $|z|=2$ 的内部有四个极点：$z_k=e^{i\pi/4+ik\pi/2}$，$k=0,1,2,3$，它们都是四阶极点，而且都对积分有贡献，所以按照式 (3) 和 (4) 来计算显然是很麻烦的。注意到被积函数在 $z$ 平面上没有其它有限的奇点，则可以应用上述定理而得到

$$
I=2\pi i\sum_{k=0}^{3}\operatorname{Res}f(z_k)
=-2\pi i\operatorname{Res}f(\infty).
$$

这样一来，我们就只需要计算 $\infty$ 的留数。在 $\infty$ 的邻域内，即当 $|z|$ 很大时，有

$$
\begin{aligned}
f(z)
&=\frac{z^{15}}{(z^4+1)^4}
=\frac1z\frac1{(1+1/z^4)^4}
=\frac1z\left(1+\frac1{z^4}\right)^{-4}\\
&=\frac1z\left(1-\frac4{z^4}+\cdots\right),
\end{aligned}
$$

可见，$c_{-1}=1$，$\operatorname{Res}f(\infty)=-1$，而 $I=2\pi i$。这样计算的工作量就大大减少了。

## §5 实积分 $\displaystyle\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta$

留数定理除了可以用来计算复变函数本身的围线积分之外，其更重要的应用在于计算某些实变函数的定积分。当然，为了应用留数定理，我们必须把定积分转化为围线积分。这通常有两种方法，其一是通过变量置换将定积分的积分区间映射为复平面上的围线。其二是将定积分的区间看作复平面的实轴上的一段，然后引入辅助曲线（比如半圆）与该段构成围线，如果辅助曲线上的积分可以比较容易算出（比如可以证明其为 0），或可以与原积分联系起来，问题就有可能解决了。本节就采用前一方法，下面两节则采用后一方法。

可以用留数定理来计算的实变函数定积分是五花八门的，其中涉及的技巧也多种多样，不是一朝一夕所能掌握的。初学者能掌握几种固定的类型即可，学有余力再涉猎其它。对于这些固定的类型，可以导出相应的计算公式。了解这些公式的来源固然是必要的，但更重要的是掌握这些公式并能熟练地使用它们。

本节讨论的积分具有标题给出的形式，其中 $R$ 是两个自变量的有理函数。对于这类积分，关键是作变量置换

$$
z=e^{i\theta}.
\tag{11}
$$

这一变换将 $\theta$ 的积分区间 $[0,2\pi]$ 映射为 $z$ 平面上的单位圆周 $|z|=1$，因此，定积分就变成了围线积分。注意到指数函数和三角函数之间的关系，易得

$$
\int_0^{2\pi}R(\cos\theta,\sin\theta)\,\mathrm d\theta
=\int_{|z|=1}
R\!\left(\frac{z+z^{-1}}{2},\frac{z-z^{-1}}{2i}\right)
\frac{\mathrm dz}{iz}.
\tag{12}
$$

右边的被积函数是 $z$ 的有理函数，所以用留数定理不难计算。对于这一类型的积分，无需记住具体的公式，只要记得变换式 (11) 就可以了。另外，如果题目中 $\theta$ 的积分区间是 $[-\pi,\pi]$，同一变换也将其映射为 $z$ 平面上的单位圆周 $|z|=1$，因此没有必要刻意将 $\theta$ 的积分区间先变换为 $[0,2\pi]$。

**例 1** 计算积分

$$
I=\int_0^{2\pi}\frac{\mathrm d\theta}{1+a\cos\theta},
\qquad -1<a<1.
$$

**解** 按式 (11) 作变量置换，得

$$
I=\int_{|z|=1}\frac{1}{1+a(z+z^{-1})/2}\frac{\mathrm dz}{iz}
=\frac2i\int_{|z|=1}\frac{\mathrm dz}{az^2+2z+a}.
$$

上式右边已是围线积分，被积函数为有理分式，其奇点就是分母的根，满足 $az^2+2z+a=0$，具体来说是

$$
z_1=\frac{-1+\sqrt{1-a^2}}{a},
\qquad
z_2=\frac{-1-\sqrt{1-a^2}}{a}.
$$

显然 $|z_2|>1$。由根与系数的关系，有 $z_1z_2=1$，故 $|z_1||z_2|=1$，从而 $|z_1|<1$。于是只有 $z_1$ 在单位圆内，即只有 $z_1$ 的留数对积分有贡献，根据式 (7)，

$$
\operatorname{Res}\!\left(\frac1{az^2+2z+a},z_1\right)
=\frac1{2az_1+2}
=\frac1{2\sqrt{1-a^2}},
$$

故

$$
I=\frac2i\cdot2\pi i\cdot\frac1{2\sqrt{1-a^2}}
=\frac{2\pi}{\sqrt{1-a^2}}.
$$

**例 2** 计算积分

$$
I=\int_0^{2\pi}\frac{\mathrm d\theta}{a+\cos\theta},
\qquad a>1.
$$

**解** 本题不需重新计算。由于 $a>1$，故 $1/a<1$，利用例 1 的结果，有

$$
I=\frac1a\int_0^{2\pi}\frac{\mathrm d\theta}{1+(1/a)\cos\theta}
=\frac1a\frac{2\pi}{\sqrt{1-1/a^2}}
=\frac{2\pi}{\sqrt{a^2-1}}.
$$

**习题** 计算积分

$$
\int_0^{2\pi}\frac{\cos^2 2\phi}{1-2p\cos\phi+p^2}\,\mathrm d\phi,
\qquad 0<p<1.
$$

## §6 实积分 $\displaystyle\int_{-\infty}^{+\infty}f(x)\,\mathrm dx$

本节讨论的积分具有标题给出的形式，其中的 $f(x)$ 通常是有理分式。关于这类积分，可以推导出一个计算公式，我们把它表述为如下的定理。

**定理** 设函数 $f(z)$ 在上半平面上有有限个孤立奇点，此外它在上半平面和实轴上解析，且当 $z\to\infty$ 时，$zf(z)\Rightarrow0$，则

$$
\int_{-\infty}^{+\infty}f(x)\,\mathrm dx
=2\pi i\sum_{\operatorname{Im}a_k>0}\operatorname{Res}f(a_k).
\tag{13}
$$

**注** ① $z\to\infty$ 时 $zf(z)\Rightarrow0$（一致趋于 0）的大意是 $zf(z)\to0$ 的速度与 $z$ 的辐角无关（在指定的辐角范围内，这里是 $[0,\pi]$）。精确地说，就是 $\forall\varepsilon>0$，$\exists R_0>0$ 与 $\theta$ 无关，当 $|z|>R_0$ 时，就有 $|zf(z)|<\varepsilon$。② 求和号下面的 $\operatorname{Im}a_k>0$ 表示求和只对上半平面的各孤立奇点进行。③ 如果 $f(x)$ 是有理分式，则易知分母的多项式次数应该至少比分子的次数大 2。④ 这一公式应该在理解和反复使用的基础上熟练掌握。

![Figure 2](...)

图 2：围线 $C_R$ 由实轴上 $[-R,R]$ 一段与半圆 $\Gamma_R$ 构成

**证明** 如图 2，作半圆 $\Gamma_R:z=Re^{i\theta}$（$0\le\theta\le\pi$），与实轴上的线段 $[-R,R]$ 一起构成围线 $C_R$，取 $R$ 充分大以使 $f(z)$ 在上半平面的孤立奇点全落在 $C_R$ 内，根据 Cauchy 留数定理，有

$$
\int_{C_R}f(z)\,\mathrm dz
=\int_{-R}^{R}f(x)\,\mathrm dx+\int_{\Gamma_R}f(z)\,\mathrm dz
=2\pi i\sum_{\operatorname{Im}a_k>0}\operatorname{Res}f(a_k).
\tag{14}
$$

由于 $zf(z)\Rightarrow0$（$z\to\infty$），故 $\forall\varepsilon>0$，$\exists R_0>0$ 与 $\theta$ 无关，当 $|z|>R_0$ 时，就有 $|zf(z)|<\varepsilon/\pi$，于是当 $R>R_0$ 时，

$$
\begin{aligned}
\left|\int_{\Gamma_R}f(z)\,\mathrm dz\right|
&\le\int_{\Gamma_R}|f(z)|\,|\mathrm dz|
=\int_0^\pi |f(z)|R\,\mathrm d\theta\\
&=\int_0^\pi|zf(z)|\,\mathrm d\theta
<\frac\varepsilon\pi\int_0^\pi\mathrm d\theta
=\varepsilon,
\end{aligned}
$$

也就是

$$
\int_{\Gamma_R}f(z)\,\mathrm dz\xrightarrow[R\to\infty]{}0.
$$

所以，对式 (14) 取 $R\to\infty$ 的极限，即得式 (13)。证毕。

由以上证明可以看出，这样求出的结果是主值 $\lim_{R\to\infty}\int_{-R}^{R}f(x)\,\mathrm dx$，而不是一般值 $\lim_{R_1\to-\infty,R_2\to+\infty}\int_{R_1}^{R_2}f(x)\,\mathrm dx$。如果一个广义积分的一般值存在，则其主值必定存在且等于一般值。但主值存在并不足以保证一般值存在。比如积分 $\int_{-\infty}^{+\infty}x\,\mathrm dx/(x^2+1)$ 的主值为 0，但一般值不存在。用留数定理导出的关于广义积分的公式，包括本节、下节的公式和本书没有介绍的情况，它们给出的都是积分的主值。

**例 1** 计算积分

$$
I=\int_0^{\infty}\frac{\mathrm dx}{x^4+1}.
$$

**解** 取

$$
f(z)=\frac1{z^4+1},
$$

它满足定理的条件，在 $z$ 平面上有四个一阶极点：$z_k=e^{i\pi/4+ik\pi/2}$，$k=0,1,2,3$，其中 $z_0$ 和 $z_1$ 在上半平面。本题用式 (7) 计算留数比用式 (5) 要方便得多。易得

$$
\operatorname{Res}f(z_k)=\frac1{4z_k^3}
=\frac{z_k}{4z_k^4}
=-\frac{z_k}{4},
$$

$$
\begin{aligned}
\operatorname{Res}f(z_0)+\operatorname{Res}f(z_1)
&=-\frac14(z_0+z_1)\\
&=-\frac{e^{i\pi/4}+e^{i3\pi/4}}4
=-\frac{e^{i\pi/4}-e^{-i\pi/4}}4
=-\frac{i}{2}\sin\frac\pi4,
\end{aligned}
$$

而

$$
I=\int_0^\infty f(x)\,\mathrm dx
=\frac12\int_{-\infty}^{+\infty}f(x)\,\mathrm dx
=\frac12\cdot2\pi i\,[\operatorname{Res}f(z_0)+\operatorname{Res}f(z_1)]
=\frac\pi{2\sqrt2}.
$$

利用

$$
\frac1{x^4+1}
=\frac1{2\sqrt2}
\left[
\frac{x+\sqrt2}{x^2+\sqrt2x+1}
-\frac{x-\sqrt2}{x^2-\sqrt2x+1}
\right],
$$

本题也可以用实变积分的方法来做，但比较麻烦。如果分母的次数越高，则实变方法的优越性就越明显。

如果 $f(z)$ 在实轴上有若干一阶极点，此外满足上述定理的条件，则公式 (13) 可推广为

$$
\int_{-\infty}^{+\infty}f(x)\,\mathrm dx
=2\pi i\left[
\sum_{\operatorname{Im}a_k>0}\operatorname{Res}f(a_k)
+\frac12\sum_{\operatorname{Im}b_k=0}\operatorname{Res}f(b_k)
\right].
\tag{15}
$$

方括号中的第二项求和对实轴上的各一阶极点进行。我们不拟给出这一公式的严格证明，但是可以这样来把握它：形象地看，可以认为实轴上的奇点有一半处于上半平面，因此，它们对积分的贡献是上半平面的奇点所贡献的一半。这样便不难掌握这一公式。另外，注意上述积分具有双重的广义性（无穷限积分，无界函数积分），而所得结果是两种意义上的主值。

**例 2** 计算积分

$$
I=\int_{-\infty}^{+\infty}\frac{\mathrm dx}{x(x-1)(x^2+1)}.
$$

**解** 取

$$
f(z)=\frac1{z(z-1)(z^2+1)},
$$

它满足公式 (15) 的条件，在上半平面上有一阶极点 $z=i$，在实轴上有一阶极点 $z=0$ 和 $z=1$。易得

$$
\operatorname{Res}f(i)=\frac{1+i}{4},
\qquad
\operatorname{Res}f(0)=-1,
\qquad
\operatorname{Res}f(1)=\frac12,
$$

故

$$
\begin{aligned}
I
&=\int_{-\infty}^{+\infty}f(x)\,\mathrm dx\\
&=2\pi i\left[\operatorname{Res}f(i)+\frac12\operatorname{Res}f(0)+\frac12\operatorname{Res}f(1)\right]\\
&=2\pi i\left(\frac{1+i}{4}-\frac12+\frac14\right)
=-\frac\pi2.
\end{aligned}
$$

**习题** 计算积分

$$
\int_{-\infty}^{+\infty}\frac{x^{2m}}{x^{2n}+1}\,\mathrm dx,
\qquad m,n\in\mathbb N^+,\quad m<n.
$$

## §7 积分 $\displaystyle\int_{-\infty}^{+\infty}f(x)e^{imx}\,\mathrm dx$

为了计算这一积分，先介绍下述的

**Jordan 引理** 设函数 $f(z)$ 在半圆 $\Gamma_R:z=Re^{i\theta}$（$0\le\theta\le\pi$）上连续，且当 $R\to\infty$ 时，$f(z)\Rightarrow0$，则

$$
\lim_{R\to\infty}\int_{\Gamma_R}f(z)e^{imz}\,\mathrm dz=0,
\qquad (m>0).
\tag{16}
$$

在给出该引理的严格证明之前，先作一点直观的说明。在上一节计算积分类型 $\int_{-\infty}^{+\infty}f(x)\,\mathrm dx$ 时，曾经用到 $\lim_{R\to\infty}\int_{\Gamma_R}f(z)\,\mathrm dz=0$ 的结论，其条件是当 $R\to\infty$ 时，$zf(z)\Rightarrow0$（虽然没有当作引理列出来）。这里我们看到，Jordan 引理所要求的条件较弱，主要原因是被积函数中含有因子

$$
e^{imz}=e^{-mR\sin\theta}e^{imR\cos\theta}.
$$

其中第一因子随着 $R$ 的增大迅速衰减（$\theta=0$ 和 $\theta=\pi$ 两点例外，但对结论没有影响），这弥补了 $f(z)$ 下降较慢的不足；第二因子在半径很大的圆弧上是 $\theta$ 的快速振荡函数，也有助于积分的收敛，但第一因子起主要作用。由上面的分析容易看出 $m>0$ 这个条件是必不可少的。

下面是 Jordan 引理的严格证明。

**证明** 由积分的基本性质有

$$
\begin{aligned}
\left|\int_{\Gamma_R}f(z)e^{imz}\,\mathrm dz\right|
&\le\int_{\Gamma_R}|f(z)|\left|e^{im(R\cos\theta+iR\sin\theta)}\right|R\,\mathrm d\theta\\
&=\int_{\Gamma_R}|f(z)|e^{-mR\sin\theta}R\,\mathrm d\theta.
\end{aligned}
$$

由于 $R\to\infty$ 时，$f(z)\Rightarrow0$，所以 $\forall\varepsilon>0$，$\exists R_0>0$ 与 $\theta$ 无关，当 $R>R_0$ 时，就有 $|f(z)|<m\varepsilon/\pi$，于是由上式得

$$
\left|\int_{\Gamma_R}f(z)e^{imz}\,\mathrm dz\right|
<\frac{m}{\pi}\varepsilon\int_0^\pi e^{-mR\sin\theta}R\,\mathrm d\theta
=\frac{2m}{\pi}\varepsilon\int_0^{\pi/2}e^{-mR\sin\theta}R\,\mathrm d\theta.
$$

**思考** 上式最后一个等号成立的理由是什么？

当 $0\le\theta\le\pi/2$ 时，有不等式

$$
\sin\theta\ge\frac2\pi\theta.
$$

（该式可由函数的图形看出，或者，由于 $\cos\theta$ 在区间 $[0,\pi/2]$ 上单调下降，就有

$$
\frac1\theta\int_0^\theta\cos x\,\mathrm dx
\ge\frac2\pi\int_0^{\pi/2}\cos x\,\mathrm dx,
$$

算出积分即得不等式），代入上式即得

$$
\begin{aligned}
\left|\int_{\Gamma_R}f(z)e^{imz}\,\mathrm dz\right|
&<\frac{2m}{\pi}\varepsilon\int_0^{\pi/2}
\exp\!\left(-\frac{2mR\theta}{\pi}\right)R\,\mathrm d\theta\\
&=\varepsilon\left.
\exp\!\left(-\frac{2mR\theta}{\pi}\right)
\right|_{\pi/2}^{0}
=\varepsilon(1-e^{-mR})<\varepsilon.
\end{aligned}
$$

由此即得式 (16)。证毕。

有了 Jordan 引理，就很容易推导出下面的公式，我们把它表述成定理。

**定理** 设函数 $f(z)$ 在上半平面上有有限个孤立奇点，此外它在上半平面和实轴上解析，且当 $z\to\infty$ 时，$f(z)\Rightarrow0$，则

$$
\int_{-\infty}^{+\infty}f(x)e^{imx}\,\mathrm dx
=2\pi i\sum_{\operatorname{Im}a_k>0}
\operatorname{Res}[f(z)e^{imz},a_k],
\qquad m>0.
\tag{17}
$$

**注** ① 这里求和也只对上半平面的各孤立奇点进行。② 注意这里的条件“当 $z\to\infty$ 时，$f(z)\Rightarrow0$”比上节的相应条件弱。如果 $f(x)$ 是有理分式，则易知分母的多项式次数应该至少比分子的次数大 1。③ 特别注意公式成立的条件 $m>0$。如果 $m<0$，我们可作变量置换 $x'=-x$，将积分化为 $\int_{-\infty}^{+\infty}f(-x)e^{i|m|x}\,\mathrm dx$，然后再用定理（如果 $f(-z)$ 满足定理的条件）。所以 $m>0$ 这一条件并不构成本质的限制。④ 这一公式应该在理解和反复使用的基础上熟练掌握。

**证明** 如图 2，作半圆 $\Gamma_R:z=Re^{i\theta}$（$0\le\theta\le\pi$），与实轴上的线段 $[-R,R]$ 一起构成围线 $C_R$，取 $R$ 充分大以使 $f(z)$ 在上半平面的孤立奇点全落在 $C_R$ 内，根据 Cauchy 留数定理，有

$$
\int_{C_R}f(z)e^{imz}\,\mathrm dz
=\int_{-R}^{R}f(x)e^{imx}\,\mathrm dx
+\int_{\Gamma_R}f(z)e^{imz}\,\mathrm dz
=2\pi i\sum_{\operatorname{Im}a_k>0}
\operatorname{Res}[f(z)e^{imz},a_k].
$$

取 $R\to\infty$ 的极限，由 Jordan 引理，半圆上的积分为 0，即得式 (17)。证毕。

**例 1** 计算积分

$$
I=\int_0^\infty\frac{\cos mx}{x^2+1}\,\mathrm dx,
\qquad m>0.
$$

**解**

$$
I=\int_0^\infty\frac{\cos mx}{x^2+1}\,\mathrm dx
=\frac12\int_{-\infty}^{+\infty}\frac{\cos mx}{x^2+1}\,\mathrm dx
=\frac12\int_{-\infty}^{+\infty}\frac{e^{imx}}{x^2+1}\,\mathrm dx,
$$

最后一个等号成立是因为最后一个积分的虚部明显为 0。取

$$
f(z)=\frac1{z^2+1},
$$

$f(z)$ 满足定理的条件，在上半平面只有一个一阶极点 $z=i$，易得

$$
\operatorname{Res}[f(z)e^{imz},i]
=\left.\frac{e^{imz}}{2z}\right|_{z=i}
=\frac{e^{-m}}{2i},
$$

故

$$
I=\frac12\cdot2\pi i\cdot\operatorname{Res}[f(z)e^{imz},i]
=\frac12\cdot2\pi i\cdot\frac{e^{-m}}{2i}
=\frac\pi2e^{-m}.
$$

**例 2** 计算积分

$$
I=\int_0^\infty\frac{x\sin mx}{x^2+a^2}\,\mathrm dx,
\qquad a>0,\quad m\in\mathbb R,\quad m\ne0.
$$

**解** 先考虑 $m>0$，

$$
I=\int_0^\infty\frac{x\sin mx}{x^2+a^2}\,\mathrm dx
=\frac12\int_{-\infty}^{+\infty}\frac{x\sin mx}{x^2+a^2}\,\mathrm dx
=\frac12\operatorname{Im}\int_{-\infty}^{+\infty}\frac{xe^{imx}}{x^2+a^2}\,\mathrm dx,
$$

取

$$
f(z)=\frac{z}{z^2+a^2},
$$

$f(z)$ 满足定理的条件，在上半平面只有一个一阶极点 $z=ai$，易得

$$
\operatorname{Res}[f(z)e^{imz},ai]
=\left.\frac{ze^{imz}}{2z}\right|_{z=ai}
=\frac{e^{-ma}}2,
$$

故

$$
\int_{-\infty}^{+\infty}\frac{xe^{imx}}{x^2+a^2}\,\mathrm dx
=2\pi i\operatorname{Res}[f(z)e^{imz},ai]
=2\pi i\frac{e^{-ma}}2
=\pi i e^{-ma},
$$

而

$$
I=\frac\pi2e^{-ma}.
$$

当 $m<0$ 时，直接利用以上 $m>0$ 时的结果，可得

$$
I=\int_0^\infty\frac{x\sin mx}{x^2+a^2}\,\mathrm dx
=-\int_0^\infty\frac{x\sin|m|x}{x^2+a^2}\,\mathrm dx
=-\frac\pi2e^{-|m|a}
=-\frac\pi2e^{ma}.
$$

如果 $f(z)$ 在实轴上有若干一阶极点，此外满足上述定理的条件，则公式 (17) 可推广为

$$
\int_{-\infty}^{+\infty}f(x)e^{imx}\,\mathrm dx
=2\pi i\left\{
\sum_{\operatorname{Im}a_k>0}\operatorname{Res}[f(z)e^{imz},a_k]
+\frac12\sum_{\operatorname{Im}b_k=0}\operatorname{Res}[f(z)e^{imz},b_k]
\right\}.
\tag{18}
$$

这一公式与式 (15) 类似，所以无需更多的说明。

**例 3** 计算积分

$$
I=\int_0^{+\infty}\frac{\sin mx}{x}\,\mathrm dx,
\qquad m\in\mathbb R,\quad m\ne0.
$$

**解** 先考虑 $m>0$，

$$
I=\frac12\int_{-\infty}^{+\infty}\frac{\sin mx}{x}\,\mathrm dx
=\frac12\operatorname{Im}\int_{-\infty}^{+\infty}\frac{e^{imx}}x\,\mathrm dx,
$$

取

$$
f(z)=\frac1z,
$$

它满足公式 (18) 的条件，在上半平面没有奇点，在实轴上有一个一阶极点 $z=0$。易得

$$
\operatorname{Res}[f(z)e^{imz},0]
=\left.e^{imz}\right|_{z=0}=1,
$$

故

$$
\int_{-\infty}^{+\infty}\frac{e^{imx}}x\,\mathrm dx
=2\pi i\cdot\frac12\operatorname{Res}[f(z)e^{imz},0]
=\pi i,
$$

而

$$
I=\frac\pi2.
$$

当 $m<0$ 时，直接利用以上 $m>0$ 时的结果，可得

$$
I=\int_0^\infty\frac{\sin mx}{x}\,\mathrm dx
=-\int_0^\infty\frac{\sin|m|x}{x}\,\mathrm dx
=-\frac\pi2.
$$

可见结果只与 $m$ 的符号有关，而与 $m$ 的大小无关。

**习题** 计算下列积分。

(1)
$$
\int_0^\infty\frac{\cos mx}{x^4+1}\,\mathrm dx,
\qquad m>0;
$$

(2)
$$
\int_0^\infty\frac{\cos mx}{(x^2+a^2)^2}\,\mathrm dx,
\qquad a>0,\quad m>0.
$$

## §8 *含一般幂函数与有理分式的无穷积分

涉及多值函数的很多积分也可以利用留数定理来计算。最简单的多值函数是一般幂函数。本节只考虑涉及一般幂函数的积分，具体形式是

$$
\int_0^\infty x^{\alpha-1}f(x)\,\mathrm dx,
$$

其中 $\alpha\in\mathbb C$，$f(z)$ 一般是有理分式，虽然理论上可以更一般。涉及对数函数的积分暂不考虑。

### §8.1 两个引理

为了本节和后面的需要，先介绍两个引理。

**大圆弧引理** 设函数 $f(z)$ 在圆弧 $S_R:z=Re^{i\theta}$（$\theta_1\le\theta\le\theta_2$）上连续，且当 $R\to\infty$ 时，$zf(z)\Rightarrow\lambda$，则

$$
\lim_{R\to\infty}\int_{S_R}f(z)\,\mathrm dz
=i\lambda\Delta\theta,
\qquad \Delta\theta=\theta_2-\theta_1.
\tag{19}
$$

![Figure 3](...)

图 3：大圆弧引理

**注** ① 这里“大圆弧引理”只是为了引用方便而采用的名称，并不是通用或标准的。② 如果 $\lambda=0$，$\theta_1=0$，$\theta_2=\pi$，与之相应的结果是 $f(z)$ 在上半平面大半圆上的积分为零，这是在 §6 用过的，只是该结果当时是隐含在公式的证明中，而没有列为引理。

**证明** 如图 3，在圆弧 $S_R:z=Re^{i\theta}$（$\theta_1\le\theta\le\theta_2$）上，

$$
\int_{S_R}\frac{\lambda}{z}\,\mathrm dz
=\lambda\int_{S_R}\frac1{Re^{i\theta}}Re^{i\theta}i\,\mathrm d\theta
=i\lambda\int_{\theta_1}^{\theta_2}\mathrm d\theta
=i\lambda\Delta\theta.
$$

故有

$$
\begin{aligned}
\left|\int_{S_R}f(z)\,\mathrm dz-i\lambda\Delta\theta\right|
&=\left|\int_{S_R}\frac{zf(z)-\lambda}{z}\,\mathrm dz\right|\\
&\le\int_{\theta_1}^{\theta_2}|zf(z)-\lambda|\,\mathrm d\theta.
\end{aligned}
$$

由于 $zf(z)\Rightarrow\lambda$（$R\to\infty$），故 $\forall\varepsilon>0$，$\exists R_0>0$ 与 $\theta$ 无关，当 $R>R_0$ 时，就有 $|zf(z)-\lambda|<\varepsilon/\Delta\theta$，$\forall\theta\in[\theta_1,\theta_2]$。于是当 $R>R_0$ 时，

$$
\left|\int_{S_R}f(z)\,\mathrm dz-i\lambda\Delta\theta\right|
<\frac\varepsilon{\Delta\theta}\int_{\theta_1}^{\theta_2}\mathrm d\theta
=\varepsilon,
$$

由此即得结论。证毕。

**小圆弧引理** 设函数 $f(z)$ 在圆弧 $S_r:z-a=re^{i\theta}$（$\theta_1\le\theta\le\theta_2$）上连续，且当 $r\to0$ 时，$(z-a)f(z)\Rightarrow\lambda$，则

$$
\lim_{r\to0}\int_{S_r}f(z)\,\mathrm dz
=i\lambda\Delta\theta,
\qquad \Delta\theta=\theta_2-\theta_1.
\tag{20}
$$

**注** ① 这里“小圆弧引理”也只是为了引用方便而采用的名称。② 如果 $a$ 是 $f(z)$ 的单极点，则 $\lambda=\operatorname{Res}f(a)$。③ 该引理的证明与大圆弧引理类似，故从略。

### §8.2 一般积分公式

现在考虑本节开头提出的积分。一般情况下，$f(z)$ 没有奇偶性。积分公式由下述定理给出。

**定理一** 设函数 $f(z)$ 在复平面上除有限个奇点外单值解析；$f(z)$ 在正实轴上只有单极点；当 $z\to0$，$z^\alpha f(z)\Rightarrow0$；当 $z\to\infty$，$z^\alpha f(z)\Rightarrow0$；则有

$$
\begin{aligned}
\int_0^\infty x^{\alpha-1}f(x)\,\mathrm dx
={}&-\pi\frac{e^{-i\pi\alpha}}{\sin\pi\alpha}
\sum_{0<\arg a_k<2\pi}
\operatorname{Res}\!
\left[z^{\alpha-1}f(z),a_k\right]\\
&-\pi\cot\pi\alpha
\sum_{b_k>0}
\operatorname{Res}\!
\left[z^{\alpha-1}f(z),b_k\right].
\end{aligned}
\tag{21}
$$

其中第一项对 $f(z)$ 在正实轴以外的奇点求和，第二项对其在正实轴上的单极点求和。

**注** ① 该公式通常用于计算 $f(z)$ 是有理分式的情况，我们也只对有理分式给出计算实例，但原则上 $f(z)$ 可以是超越函数，只要满足定理的条件即可。具体来说，如果 $g(z)$ 是有理分式并满足定理的条件，则 $f(z)=g(z)e^{1/(z-a)}$、$f(z)=g(z)\sin[1/(z-a)]$（其中 $\operatorname{Im}a\ne0$）等函数亦满足，它们就不是有理分式。这时候公式仍然适用，但是计算 $a$ 点的留数会比较困难。② $f(z)$ 在正实轴以外可以有本性奇点，当然有本性奇点的情形定理不适于有理分式。比如刚刚举例的两个函数。③ $f(z)$ 具有有限个奇点，则必然都是孤立奇点。

**证明** 为叙述方便，暂设 $f(z)$ 在正实轴上只有一个单极点 $b$。取辅助函数为 $F(z)=z^{\alpha-1}f(z)$。又取图 4 所示的围线 $C$，它由以下各段构成：以原点为中心、半径为 $R$ 的大半圆 $\Gamma_R$；以原点为中心、半径为 $r$ 的小半圆 $\Gamma_r$；实轴上的线段 $[-R,R]$，遇到单极点则以小半圆绕过。所有半圆位于上半平面。正实轴上避开 $b$ 的两段记为 $L_1=[r,b-r]\cup[b+r,R]$，下岸相应的两段记为 $L_2$。注意各圆弧均以逆时针方向为默认方向，故

$$
C=L_1+L_2^-+\Gamma_R+\Gamma_r^-+\Gamma_b^-+\Gamma_{\bar b}^-.
$$

取 $R$ 足够大、$r$ 足够小，使 $f(z)$ 在正实轴以外的奇点均在围线 $C$ 内。根据留数定理，可得

$$
\begin{aligned}
\int_C F(z)\,\mathrm dz
={}&\int_{L_1}F(z)\,\mathrm dz-\int_{L_2}F(z)\,\mathrm dz
+\int_{\Gamma_R}F(z)\,\mathrm dz-\int_{\Gamma_r}F(z)\,\mathrm dz\\
&-\int_{\Gamma_b}F(z)\,\mathrm dz-\int_{\Gamma_{\bar b}}F(z)\,\mathrm dz
=2\pi i\sum_{0<\arg a_k<2\pi}\operatorname{Res}[F(z),a_k].
\end{aligned}
\tag{22}
$$

![Figure 4](...)

图 4：§8 定理一的证明

在 $L_1$ 上，可取 $z=x$，则在 $L_2$ 上，有 $z=xe^{i2\pi}$，于是可得

$$
\begin{aligned}
\int_{L_1}F(z)\,\mathrm dz-\int_{L_2}F(z)\,\mathrm dz
&=(1-e^{i2\pi\alpha})
\left[\int_r^{b-r}+\int_{b+r}^{R}x^{\alpha-1}f(x)\,\mathrm dx\right]\\
&\to(1-e^{i2\pi\alpha})\int_0^\infty x^{\alpha-1}f(x)\,\mathrm dx,
\qquad r\to0,\ R\to\infty.
\end{aligned}
$$

由于按已知条件积分存在，故取极限没有问题。根据已知条件和大圆弧引理，有

$$
\lim_{R\to\infty}\int_{\Gamma_R}F(z)\,\mathrm dz=0.
$$

根据已知条件和小圆弧引理，又有

$$
\lim_{r\to0}\int_{\Gamma_r}F(z)\,\mathrm dz=0.
$$

由于 $b$ 是 $f(z)$ 的单极点，故有

$$
\begin{aligned}
\lim_{r\to0}\int_{\Gamma_b}F(z)\,\mathrm dz
&=\pi i\lim_{z\to b}\left[(z-b)z^{\alpha-1}f(z)\right]\\
&=\pi i b^{\alpha-1}\operatorname{Res}[f(z),b]
=\pi i\operatorname{Res}[z^{\alpha-1}f(z),b].
\end{aligned}
$$

类似地

$$
\lim_{r\to0}\int_{\Gamma_{\bar b}}F(z)\,\mathrm dz
=\pi i e^{i2\pi\alpha}
\operatorname{Res}[z^{\alpha-1}f(z),b].
$$

对式 (22) 取 $r\to0$、$R\to\infty$ 的极限，综合以上诸式，并推广到 $f(z)$ 在正实轴上有多个单极点的情形（各极点可以像 $b$ 一样处理），最后即得式 (21)。证毕。

下面用上述公式计算几个实例。

**例 1** 计算积分

$$
I=\int_0^\infty\frac{x^{\alpha-1}}{x+b}\,\mathrm dx,
\qquad b>0,\quad 0<\operatorname{Re}\alpha<1.
$$

**解** 为了使得积分存在，条件 $0<\operatorname{Re}\alpha<1$ 是必须的。取 $f(z)=1/(z+b)$，它在复平面上只有一个单极点 $z=-b$，更确切地说是 $z=be^{i\pi}$。易得

$$
\operatorname{Res}\left[\frac{z^{\alpha-1}}{z+b},be^{i\pi}\right]
=-b^{\alpha-1}e^{i\pi\alpha}.
$$

代入公式即得

$$
\int_0^\infty\frac{x^{\alpha-1}}{x+b}\,\mathrm dx
=\frac\pi{\sin\pi\alpha}b^{\alpha-1},
\qquad b>0,\quad0<\operatorname{Re}\alpha<1.
$$

如果 $b=1$，则结果可以进一步化简。

**例 2** 计算积分

$$
I=\int_0^\infty\frac{x^{\alpha-1}}{b-x}\,\mathrm dx,
\qquad b>0,\quad0<\operatorname{Re}\alpha<1.
$$

**解** 为了使得积分存在，条件 $0<\operatorname{Re}\alpha<1$ 是必须的。取 $f(z)=1/(b-z)$，它在复平面上只有一个单极点 $z=b$，位于正实轴上。易得

$$
\operatorname{Res}\left[\frac{z^{\alpha-1}}{b-z},b\right]
=-b^{\alpha-1}.
$$

代入公式即得

$$
\int_0^\infty\frac{x^{\alpha-1}}{b-x}\,\mathrm dx
=b^{\alpha-1}\pi\cot\pi\alpha,
\qquad b>0,\quad0<\operatorname{Re}\alpha<1.
$$

如果 $b=1$，则结果可以进一步化简。

以上两个例子都是教材书上常见的，多数书上的做法是每个例子单独选取辅助函数和适当的围线进行计算。而我们的做法是给出一般公式，由此可以计算大量的实例，每个实例只需要选取适当的函数 $f(z)$（这通常可以通过观察被积函数直接得到）并计算留数即可。下面举一个复杂一点的实例。

**例 3** 计算积分

$$
I=\int_0^\infty\frac{x^{\alpha-1}}{x^3+b^3}\,\mathrm dx,
\qquad b>0,\quad0<\operatorname{Re}\alpha<3.
$$

**解** 条件 $0<\operatorname{Re}\alpha<3$ 保证积分存在。取 $f(z)=1/(z^3+b^3)$，它在复平面上只有三个单极点

$$
z_0=be^{i\pi/3},\qquad
z_1=be^{i\pi},\qquad
z_2=be^{i5\pi/3}.
$$

它们都不在正实轴上。易得

$$
\operatorname{Res}\left[\frac{z^{\alpha-1}}{z^3+b^3},z_k\right]
=-\frac{z_k^\alpha}{3b^3},
\qquad k=0,1,2.
$$

代入公式，整理即得

$$
\int_0^\infty\frac{x^{\alpha-1}}{x^3+b^3}\,\mathrm dx
=\frac{\pi b^{\alpha-3}}{3\sin(\pi\alpha/3)},
\qquad b>0,\quad0<\operatorname{Re}\alpha<3.
$$

如果 $b=1$，则结果可以进一步化简。

### §8.3 有理分式具有奇偶性的积分公式

现在考虑 $f(z)$ 具有奇偶性的情形。如果 $f(z)$ 是奇函数，则 $zf(z)$ 是偶函数，所以只需适当重新定义 $\alpha$ 即可化为偶函数的情形。所以下面只考虑偶函数。虽然这时候积分公式 (21) 仍然有效，但是可以化简。

由于 $f(-z)=f(z)$，如果 $a$（$0\le\arg a<\pi$）是 $f(z)$ 的奇点，则 $ae^{i\pi}$ 也是它的奇点。后者的留数与 $a$ 点的留数是相关的。事实上，根据定义，有

$$
\operatorname{Res}[z^{\alpha-1}f(z),ae^{i\pi}]
=\frac1{2\pi i}\int_{|z-ae^{i\pi}|=r}z^{\alpha-1}f(z)\,\mathrm dz,
$$

其中取 $r$ 足够小，以便围线内没有 $f(z)$ 的其它奇点。作变量置换 $z=\zeta e^{i\pi}$，然后将积分变量 $\zeta$ 重新记作 $z$，则得

$$
\operatorname{Res}[z^{\alpha-1}f(z),ae^{i\pi}]
=e^{i\pi\alpha}\frac1{2\pi i}\int_{|z-a|=r}z^{\alpha-1}f(z)\,\mathrm dz
=e^{i\pi\alpha}\operatorname{Res}[z^{\alpha-1}f(z),a].
$$

这个结果对于正实轴上的 $a$（即前述公式中的 $b_k$）当然也是成立的。利用这一关系，将成对的奇点的留数都用上半平面和正实轴上的奇点留数表达出，整理即得下述定理。

**定理二** 设函数 $f(z)$ 在复平面上除有限个奇点外单值解析；$f(-z)=f(z)$；$f(z)$ 在正实轴上只有单极点；当 $z\to0$，$z^\alpha f(z)\Rightarrow0$；当 $z\to\infty$，$z^\alpha f(z)\Rightarrow0$；则有

$$
\begin{aligned}
\int_0^\infty x^{\alpha-1}f(x)\,\mathrm dx
={}&-\pi\frac{e^{-i\pi\alpha/2}}{\sin(\pi\alpha/2)}
\sum_{\text{upper}}\operatorname{Res}[z^{\alpha-1}f(z)]\\
&-\pi\cot\frac{\pi\alpha}{2}
\sum_{\text{positive}}\operatorname{Res}[z^{\alpha-1}f(z)].
\end{aligned}
\tag{23}
$$

其中第一项对 $f(z)$ 在上半平面的奇点求和，第二项对其在正实轴上的单极点求和。

**注** ① 应用这个公式时，需要计算留数的奇点比前面的一般公式少一半，所以当 $f(z)$ 奇点较多时，该公式更为方便。② 取辅助函数为 $F(z)=z^{\alpha-1}f(z)$，并取围线 $C$ 由以下各段构成：以原点为中心的大半圆 $\Gamma_R$ 和小半圆 $\Gamma_r$；实轴上的线段 $[-R,R]$，遇到单极点则以小半圆绕过；所有半圆位于上半平面。然后利用留数定理，也可以导出公式 (23)。

**例 4** 计算积分

$$
I=\int_0^\infty\frac{x^{\alpha-1}}{x^2+a^2}\,\mathrm dx,
\qquad |\arg a|<\frac\pi2,\quad0<\operatorname{Re}\alpha<2.
$$

**解** 条件 $0<\operatorname{Re}\alpha<2$ 保证积分存在。取 $f(z)=1/(z^2+a^2)$，它是偶函数，在上半平面只有一个单极点 $ae^{i\pi/2}$，正实轴上没有奇点，所以只需计算一个留数。易得

$$
\operatorname{Res}\left[\frac{z^{\alpha-1}}{z^2+a^2},ae^{i\pi/2}\right]
=-\frac12a^{\alpha-2}e^{i\pi\alpha/2}.
$$

代入公式，即得

$$
\int_0^\infty\frac{x^{\alpha-1}}{x^2+a^2}\,\mathrm dx
=\frac{\pi a^{\alpha-2}}{2\sin(\pi\alpha/2)},
\qquad |\arg a|<\frac\pi2,\quad0<\operatorname{Re}\alpha<2.
$$

如果 $a=1$，则结果可以进一步化简。

本题当然也可以用公式 (21) 计算，此时尚需计算 $ae^{i3\pi/2}$ 的留数，一定代入公式中的第一和式，整理可得同样的结果。

本题被积函数的分母具有 $z^2+c$ 的形式，如果 $|\arg c|<\pi$，则可以令 $c=a^2$，从而分母成为 $z^2+a^2$，其中 $|\arg a|<\pi/2$。如果 $|\arg c|=\pi$，则可以令 $c=-b^2$，从而分母成为 $z^2-b^2$，其中 $b>0$，此时积分的贡献来自于公式 (23) 中的第二和式（也是只有一项）。结果是

$$
\int_0^\infty\frac{x^{\alpha-1}}{x^2-b^2}\,\mathrm dx
=-\frac\pi2b^{\alpha-2}\cot\frac{\pi\alpha}{2},
\qquad b>0,\quad0<\operatorname{Re}\alpha<2.
$$

具体计算留给读者作为练习。两种情况合起来涵盖了所有非零的复数 $c$。

## §9 *含一般幂函数、有理分式与三角函数的无穷积分

根据 §7 的公式和例题，我们看到，如果 $f(x)$ 是偶函数且满足一定的条件，很容易计算形如 $\int_0^\infty f(x)\cos qx\,\mathrm dx$（本节将 §7 的 $m$ 改写为 $q$）。但是如果 $f(x)$ 是奇函数，如何计算同样的积分呢？对于含有正弦函数的积分，也可以提出类似的问题。更一般地，可以考虑如何计算以下类型的积分：

$$
\int_0^\infty x^\lambda f(x)\cos qx\,\mathrm dx,
\qquad
\int_0^\infty x^\lambda f(x)\sin qx\,\mathrm dx.
\tag{24}
$$

其中 $\lambda$ 一般可以是复数，$q>0$，$f(x)$ 一般是有理分式且为偶函数，暂不考虑更具体限制。如果 $f(x)$ 是奇函数，只需适当调整 $\lambda$ 的定义即可归结为 $f(x)$ 是偶函数的情况。

根据一般的经验，我们尝试选取辅助函数 $F(z)=z^\lambda f(z)e^{iqz}$，并选取包含实轴、上半平面的大半圆和小半圆的围线（即下文定义的围线 $C_0$），考虑 $F(z)$ 沿该围线的积分。计算表明，这样只能得到一个既包含 $\cos qx$ 又包含 $\sin qx$ 的积分，不能分别得到两者。

对于式 (24) 的两个积分，如果我们分别选择辅助函数为 $F(z)=z^\lambda f(z)\cos qz$ 和 $F(z)=z^\lambda f(z)\sin qz$，则可以克服上述困难。但这两个计算中，大半圆上的积分难以计算。下面介绍的引理有助于解决后一困难，从而使得上述积分的计算迎刃而解。[^1]

本节假设 $q,b>0$ 和 $\lambda,\rho,\mu,a\in\mathbb C$。

**引理** 设 $\Gamma_R$ 是上半平面的⼤半圆：$z=Re^{i\theta}$，$0\le\theta\le\pi$，则有

$$
\lim_{R\to\infty}\int_{\Gamma_R}\frac{e^{-iqz}}{z^\rho}\,\mathrm dz
=-\frac{2\pi q^{\rho-1}e^{-i\rho\pi/2}}{\Gamma(\rho)},
\qquad \operatorname{Re}\rho>0.
\tag{25}
$$

[^1]: 该引理可参见以下文献：林琼桂，含幂函数、有理分式与三角函数的无穷积分——一个引理及其应用，大学物理 34 (5) (2015) 1-4。这是一个更一般的引理的特殊情况。一般形式的引理涉及 Bessel 函数，有助于建立若干一般的积分公式，后者使得包含 Bessel 函数或相关函数的几类无穷积分的计算大为简化。参见文献：Qiong-Gui Lin, Infinite integrals involving Bessel functions by an improved approach of contour integration and the residue theorem, Ramanujan J. 35 (2014) 443-466.

![Figure 5](...)

图 5：§9 引理的证明

**注** 过去我们在计算含有三角函数的无穷积分时严重依赖于 Jordan 引理，这极大地限制了围线积分方法的应用范围。上述引理可以与 Jordan 引理形成互补，从而改变现状，大大地提高围线积分方法的威力。必须指出，尽管在上半平面当 $z\to\infty$ 时 $e^{-iqz}\to\infty$，但不可先入为主地认为上述积分的极限为无穷大或不存在。事实上，对于 $\rho=n\in\mathbb N^+$ 的特殊情况，利用约当引理和留数定理，就有

$$
\lim_{R\to\infty}\int_{\Gamma_R}\frac{e^{-iqz}}{z^n}\,\mathrm dz
=\lim_{R\to\infty}\int_{|z|=R}\frac{e^{-iqz}}{z^n}\,\mathrm dz
=2\pi i\operatorname{Res}\!\left[\frac{e^{-iqz}}{z^n},0\right].
$$

算出留数，即可验证结果与式 (25) 右边给出的值一致。

**证明** 取围线 $C_0$ 由以下各段构成：以原点为中心、半径为 $R$ 的大半圆 $\Gamma_R$；以原点为中心、半径为 $r$ 的小半圆 $\Gamma_r$；实轴上两段 $L_1=[r,R]$ 和 $L_2=[-R,-r]$，参见图 5。以上各半圆均位于上半平面。注意 $C_0$ 和各圆弧均以逆时针方向为默认方向，故作为 $C_0$ 的构成部分的小圆弧实际上是 $\Gamma_r^-$。默认 $L_1$ 上 $\arg z=0$，而 $L_2$ 上 $\arg z=\pi$。利用柯西积分定理，可得

$$
\int_{\Gamma_R}\frac{e^{-iqz}}{z^\rho}\,\mathrm dz
=\int_{\Gamma_r}\frac{e^{-iqz}}{z^\rho}\,\mathrm dz
-\int_r^R\frac{e^{-iqx}+e^{-i\rho\pi}e^{iqx}}{x^\rho}\,\mathrm dx.
$$

令 $R\to\infty$，由于 $\operatorname{Re}\rho>0$，右边的极限存在而且是 $\rho$ 的解析函数，故有

$$
\lim_{R\to\infty}\int_{\Gamma_R}\frac{e^{-iqz}}{z^\rho}\,\mathrm dz
=\int_{\Gamma_r}\frac{e^{-iqz}}{z^\rho}\,\mathrm dz
-\int_r^\infty\frac{e^{-iqx}+e^{-i\rho\pi}e^{iqx}}{x^\rho}\,\mathrm dx.
\tag{26}
$$

于是右边也是 $\rho$ 的解析函数，这对于最后一步至关重要。如果我们还假定 $\operatorname{Re}\rho<1$，则可以取极限 $r\to0$。此时右边第一项的极限为零，第二项则可利用以下公式来计算：

$$
\int_0^\infty x^{\mu-1}e^{\pm iqx}\,\mathrm dx
=q^{-\mu}\Gamma(\mu)e^{\pm i\mu\pi/2},
\qquad 0<\operatorname{Re}\mu<1.
\tag{27}
$$

上式可以用围线积分方法证明（参看 §10），也可以用微积分中的方法证明，然后利用 $\Gamma$ 函数的互余关系 $\Gamma(\rho)\Gamma(1-\rho)=\pi/\sin\pi\rho$（证明从略），我们就在假设 $0<\operatorname{Re}\rho<1$ 下获得式 (25)。

由于式 (25) 的右边是 $\rho$ 的解析函数（证明从略），而如上所述，只要 $\operatorname{Re}\rho>0$，左边也是 $\rho$ 的解析函数，根据解析函数的唯一性定理，以上额外假设 $\operatorname{Re}\rho<1$ 就可以去掉。于是式 (25) 在条件 $\operatorname{Re}\rho>0$ 下成立。证毕。

这个引理就可以导出以下两个一般的积分公式。对于定理和实例中出现的幂函数，其辐角规定与引理相同。[^2]

[^2]: 菲赫金哥尔茨 Г. М.《微积分学教程》第二卷 [M]，8 版，徐诚翰、冷生明、梁文骐译（中译本 2 版），北京：高等教育出版社，2006。参看该书第 652 页的计算。

![Figure 6](...)

图 6：§9 定理一的证明

**定理一** 设函数 $f(z)$ 在复平面上除有限个奇点外单值解析；$f(-z)=f(z)$；$f(z)$ 在正实轴上只有极点，使其为 $f(z)\cos qz$ 的不超过一阶的极点；当 $z\to0$，$z^{\lambda+1}f(z)\Rightarrow0$；当 $z\to\infty$，$z^\lambda f(z)\Rightarrow0$；对于足够大的 $|z|$，有 Laurent 展开式（由以上假定所推论）

$$
f(z)=\sum_{k=0}^{\infty}\frac{c_k}{z^{2k+2k_0}},
\qquad k_0\in\mathbb Z,\quad 2k_0>\operatorname{Re}\lambda,
\tag{28}
$$

则有

$$
\begin{aligned}
\int_0^\infty x^\lambda f(x)\cos qx\,\mathrm dx
={}&\pi i e^{-i\lambda\pi/2}\sec\frac{\lambda\pi}{2}
\sum_{\text{upper}}\operatorname{Res}[z^\lambda f(z)\cos qz]\\
&+\pi\tan\frac{\lambda\pi}{2}
\sum_{\text{positive}}\operatorname{Res}[z^\lambda f(z)\cos qz]\\
&+(-1)^{k_0}\frac\pi2\sec\frac{\lambda\pi}{2}
\sum_{k=0}^{\infty}
\frac{(-1)^k c_k q^{2k+2k_0-\lambda-1}}
{\Gamma(2k+2k_0-\lambda)}.
\end{aligned}
\tag{29}
$$

其中第一项对 $f(z)$ 在上半平面的奇点求和，第二项对其在正实轴上的极点求和。

**证明** 为叙述方便，暂设 $f(z)$ 在正实轴上只有一个极点 $b$，因其为偶函数，则 $-b$ 亦必为极点。取辅助函数为 $F(z)=z^\lambda f(z)\cos qz$。又取围线 $C$ 由以下各段构成：大半圆 $\Gamma_R$；小半圆 $\Gamma_r$；以 $b$ 为中心、半径为 $r$ 的小半圆 $\Gamma_b$；以 $-b$ 为中心、半径为 $r$ 的小半圆 $\Gamma_{-b}$；实轴上两段 $L_1=[r,b-r]\cup[b+r,R]$ 和 $L_2=[-R,-b-r]\cup[-b+r,-r]$，参见图 6。以上各半圆均位于上半平面。注意 $C$ 和各圆弧均以逆时针方向为默认方向，故作为 $C$ 的构成部分的小圆弧实际上是 $\Gamma_r^-$、$\Gamma_b^-$ 和 $\Gamma_{-b}^-$。默认 $L_1$ 上 $\arg z=0$。取 $R$ 足够大、$r$ 足够小，使 $f(z)$ 在上半平面的奇点均在围线 $C$ 内。根据留数定理，可得

$$
\begin{aligned}
\int_C F(z)\,\mathrm dz
={}&\int_{L_1}F(z)\,\mathrm dz+\int_{L_2}F(z)\,\mathrm dz
+\int_{\Gamma_R}F(z)\,\mathrm dz\\
&-\int_{\Gamma_r}F(z)\,\mathrm dz
-\int_{\Gamma_b}F(z)\,\mathrm dz
-\int_{\Gamma_{-b}}F(z)\,\mathrm dz\\
={}&2\pi i\sum_{\text{upper}}\operatorname{Res}F(z).
\end{aligned}
\tag{30}
$$

根据被积函数的奇偶性，易得

$$
\int_{L_1}F(z)\,\mathrm dz+\int_{L_2}F(z)\,\mathrm dz
\to(1+e^{i\lambda\pi})\int_0^\infty x^\lambda f(x)\cos qx\,\mathrm dx,
\qquad r\to0,\ R\to\infty.
\tag{31}
$$

由于按已知条件积分存在，故取极限没有问题。根据已知条件，易得

$$
\int_{\Gamma_r}F(z)\,\mathrm dz\to0,
\qquad r\to0.
\tag{32}
$$

由于 $b$ 最多只是 $f(z)\cos qz$ 的单极点，故有

$$
\int_{\Gamma_b}F(z)\,\mathrm dz
\to\pi i\operatorname{Res}[z^\lambda f(z)\cos qz,b],
\qquad r\to0.
\tag{33a}
$$

对于沿 $\Gamma_{-b}$ 的积分，我们作变换 $z=\zeta e^{i\pi}$，则可得

$$
\int_{\Gamma_{-b}}F(z)\,\mathrm dz
=-e^{i\lambda\pi}\int_{\Gamma_b}F(z)\,\mathrm dz
\to-\pi i e^{i\lambda\pi}\operatorname{Res}[z^\lambda f(z)\cos qz,b],
\qquad r\to0,
\tag{33b}
$$

其中 $\Gamma_b$ 是 $\Gamma_{-b}$ 对实轴的镜像。对于沿 $\Gamma_R$ 的积分，由已知条件和引理，易得

$$
\lim_{R\to\infty}\int_{\Gamma_R}F(z)\,\mathrm dz
=\frac12\lim_{R\to\infty}\int_{\Gamma_R}z^\lambda f(z)e^{-iqz}\,\mathrm dz.
\tag{34}
$$

由于取 $R$ 足够大，我们就可以用展开式 (28)，然后利用引理，可得

$$
\lim_{R\to\infty}\int_{\Gamma_R}F(z)\,\mathrm dz
=(-1)^{k_0+1}\pi e^{i\lambda\pi/2}
\sum_{k=0}^{\infty}
\frac{(-1)^k c_k q^{2k+2k_0-\lambda-1}}
{\Gamma(2k+2k_0-\lambda)}.
\tag{35}
$$

对式 (30) 取 $r\to0$、$R\to\infty$ 的极限，综合以上诸式，并推广到 $f(z)$ 在正实轴上有多个极点的情形（各极点可以像 $b$ 一样处理），最后即得式 (29)。证毕。

**定理二** 设函数 $f(z)$ 在复平面上除有限个奇点外单值解析；$f(-z)=f(z)$；$f(z)$ 在正实轴上只有极点，使其为 $f(z)\sin qz$ 的不超过一阶的极点；当 $z\to0$，$z^{\lambda+2}f(z)\Rightarrow0$；当 $z\to\infty$，$z^\lambda f(z)\Rightarrow0$；对于足够大的 $|z|$，有 Laurent 展开式 (28)；则有

$$
\begin{aligned}
\int_0^\infty x^\lambda f(x)\sin qx\,\mathrm dx
={}&-\pi e^{-i\lambda\pi/2}\csc\frac{\lambda\pi}{2}
\sum_{\text{upper}}\operatorname{Res}[z^\lambda f(z)\sin qz]\\
&-\pi\cot\frac{\lambda\pi}{2}
\sum_{\text{positive}}\operatorname{Res}[z^\lambda f(z)\sin qz]\\
&-(-1)^{k_0}\frac\pi2\csc\frac{\lambda\pi}{2}
\sum_{k=0}^{\infty}
\frac{(-1)^k c_k q^{2k+2k_0-\lambda-1}}
{\Gamma(2k+2k_0-\lambda)}.
\end{aligned}
\tag{36}
$$

这一定理的证明与定理一非常类似。具体细节从略。注意 $f(z)$ 在 $z\to0$ 时的条件比定理一宽松。

下面对以上公式作一些补充说明。第一，$f(z)$ 在正实轴上的极点一般只能是单极点。但如果在某极点 $b$ 处 $\cos qb=0$（对定理一）或 $\sin qb=0$（对定理二），则 $b$ 允许是 $f(z)$ 的二阶极点。对于正实轴上的每个极点，取积分都是取柯西主值，这是围线积分方法的惯例。第二，定理中的条件 $f(-z)=f(z)$ 对于实际应用不构成限制，如果 $f(-z)=-f(z)$，只需重新定义 $\lambda$ 和 $f(z)$，并满足定理中的其它条件，即可应用定理。对于没有奇偶性的函数，可以分解为偶函数与奇函数之和，并分别应用定理。第三，公式中最后一项是无穷级数，对于实例，多数情况下可以表示为熟知的特殊函数。第四，虽然我们只对有理分式给出计算实例，但原则上 $f(z)$ 可以是超越函数，只要满足定理的条件即可。

下面给出几个计算实例。

**例 1** 利用公式 (29)，不难计算积分

$$
\int_0^\infty\frac{x^\lambda\cos qx}{x^{2m}+a^{2m}}\,\mathrm dx,
\qquad m\in\mathbb N^+,\quad -1<\operatorname{Re}\lambda<2m,
\quad |\arg a|<\frac\pi{2m}.
$$

由于结果中出现广义超几何函数，为便于理解，此处只考虑 $m=1$ 的情况（下同）。此时 $f(z)=1/(z^2+a^2)$，$-1<\operatorname{Re}\lambda<2$，$|\arg a|<\pi/2$。$f(z)$ 在正实轴上没有奇点，在上半平面只有一个单极点 $z=ae^{i\pi/2}$，容易求得该点的留数为

$$
\operatorname{Res}[z^\lambda f(z)\cos qz,ae^{i\pi/2}]
=-\frac{i}{2}a^{\lambda-1}e^{i\lambda\pi/2}\cosh qa.
$$

当 $|z|>|a|$，$f(z)$ 的 Laurent 级数具有式 (28) 的形式，其中 $k_0=1$，$c_k=(-1)^k a^{2k}$。于是可得结果为

$$
\int_0^\infty\frac{x^\lambda\cos qx}{x^2+a^2}\,\mathrm dx
=\frac{\pi a^{\lambda-1}}{2\cos(\lambda\pi/2)}
\left[
\cosh qa-\sum_{k=0}^{\infty}
\frac{(qa)^{2k+1-\lambda}}{\Gamma(2k+2-\lambda)}
\right].
\tag{37}
$$

方括号内第二项与双曲函数的级数非常类似。特别地，当 $\lambda=0$，它就是双曲正弦函数；当 $\lambda=1$，则是双曲余弦函数。

**注** 上式左边是含参数的无穷积分，右边出现无穷级数。从数值计算的角度来看，给定参数，计算右边未必比计算左边来得容易，这与我们把积分“积出来”的初衷似乎不符。对于本题这样的被积函数，所谓把积分“积出来”通常是指把结果表示为参数的单个（或有限个）初等函数。但是，这在客观上不一定总是可行的。尽管如此，这个结果还是很有意义的。第一，它给出了结果对参数依赖关系的显明形式，尽管函数关系略嫌复杂。至少它在两种数学形式之间建立了一个数量关系，这在数学上就是有意义的。第二，上式中的无穷级数实际上可以用广义超几何函数表示，也就是说，结果是一项初等函数和一项（熟知的）特殊函数之和，给定参数，利用现成的数学软件或函数表即可得出数值结果。第三，如果参数满足 $qa\ll1$，则上式中的无穷级数只需取少数几项甚至一两项就可以获得足够精确的结果，可见上式即使在计算上也是颇有意义的。

当 $\lambda=0$，式 (37) 与应用式 (17) 计算的结果一致。当 $\lambda=1$，它成为 $0/0$ 型，取极限并利用 L'Hospital 法则，可得

$$
\int_0^\infty\frac{x\cos qx}{x^2+a^2}\,\mathrm dx
=\sum_{k=0}^{\infty}\frac{\psi(2k+1)}{(2k)!}(qa)^{2k}
-\cosh qa\ln qa.
\tag{38}
$$

其中 $\psi(x)=\Gamma'(x)/\Gamma(x)$。

**例 2** 考虑 $f(z)=1/(z^2-b^2)$ 和式 (29) 的积分，其中 $-1<\operatorname{Re}\lambda<2$，$b>0$。$f(z)$ 在正实轴上有单极点 $z=b$，在上半平面没有奇点。容易求得 $z=b$ 点的留数为

$$
\operatorname{Res}[z^\lambda f(z)\cos qz,b]
=\frac12b^{\lambda-1}\cos qb.
$$

当 $|z|>b$，$f(z)$ 的 Laurent 级数具有式 (28) 的形式，其中 $k_0=1$，$c_k=b^{2k}$。于是可得结果为

$$
\int_0^\infty\frac{x^\lambda\cos qx}{x^2-b^2}\,\mathrm dx
=\frac{\pi b^{\lambda-1}}{2\cos(\lambda\pi/2)}
\left[
\sin\frac{\lambda\pi}{2}\cos qb
-\sum_{k=0}^{\infty}
\frac{(-1)^k(qb)^{2k+1-\lambda}}{\Gamma(2k+2-\lambda)}
\right].
\tag{39}
$$

上式方括号内第二项与三角函数的级数非常类似。特别地，当 $\lambda=0$，它就是正弦函数；当 $\lambda=1$，则是余弦函数。

当 $\lambda=0$，式 (39) 与应用式 (18) 计算的结果一致。当 $\lambda=1$，它成为 $0/0$ 型，取极限并利用 L'Hospital 法则，可得

$$
\int_0^\infty\frac{x\cos qx}{x^2-b^2}\,\mathrm dx
=\sum_{k=0}^{\infty}\frac{(-1)^k\psi(2k+1)}{(2k)!}(qb)^{2k}
-\cos qb\ln qb.
\tag{40}
$$

**例 3** 类似例 1，计算可得

$$
\int_0^\infty\frac{x^\lambda\sin qx}{x^2+a^2}\,\mathrm dx
=\frac{\pi a^{\lambda-1}}{2\sin(\lambda\pi/2)}
\left[
\sum_{k=0}^{\infty}\frac{(qa)^{2k+1-\lambda}}{\Gamma(2k+2-\lambda)}
-\sinh qa
\right],
\quad -2<\operatorname{Re}\lambda<2,\quad |\arg a|<\frac\pi2.
\tag{41}
$$

读者可自行考虑 $\lambda=1$，$\lambda=-1$ 和 $\lambda=0$ 等特殊情况。特别是，当 $\lambda=0$，它成为 $0/0$ 型，取极限并利用 L'Hospital 法则，可得

$$
\int_0^\infty\frac{\sin qx}{x^2+a^2}\,\mathrm dx
=\frac1a\left[
\sum_{k=0}^{\infty}\frac{\psi(2k+2)}{(2k+1)!}(qa)^{2k+1}
-\sinh qa\ln qa
\right].
\tag{42}
$$

**例 4** 类似例 2，计算可得

$$
\int_0^\infty\frac{x^\lambda\sin qx}{x^2-b^2}\,\mathrm dx
=\frac{\pi b^{\lambda-1}}{2\sin(\lambda\pi/2)}
\left[
\sum_{k=0}^{\infty}\frac{(-1)^k(qb)^{2k+1-\lambda}}{\Gamma(2k+2-\lambda)}
-\cos\frac{\lambda\pi}{2}\sin qb
\right],
\quad -2<\operatorname{Re}\lambda<2.
\tag{43}
$$

读者可自行考虑 $\lambda=1$，$\lambda=-1$ 和 $\lambda=0$ 等特殊情况。特别是，当 $\lambda=0$，它成为 $0/0$ 型，取极限并利用 L'Hospital 法则，可得

$$
\int_0^\infty\frac{\sin qx}{x^2-b^2}\,\mathrm dx
=\frac1b\left[
\sum_{k=0}^{\infty}\frac{(-1)^k\psi(2k+2)}{(2k+1)!}(qb)^{2k+1}
-\sin qb\ln qb
\right].
\tag{44}
$$

我们还可以举出 $f(z)$ 在正实轴上有二阶极点的例子，如果在该极点 $b$ 处 $\cos qb=0$ 或 $\sin qb=0$。由于计算上并无困难，故从略。此外，我们还可以考虑 $f(z)$ 为多值函数的例子，不过此时需要单独计算，而无法直接应用一般公式，亦从略。

最后，如果取 $f(z)=1$，则利用公式 (29) 和 (36) 可以分别得出

$$
\int_0^\infty x^{\mu-1}\cos qx\,\mathrm dx
=q^{-\mu}\Gamma(\mu)\cos\frac{\mu\pi}{2},
\qquad 0<\operatorname{Re}\mu<1,
\tag{45a}
$$

$$
\int_0^\infty x^{\mu-1}\sin qx\,\mathrm dx
=q^{-\mu}\Gamma(\mu)\sin\frac{\mu\pi}{2},
\qquad -1<\operatorname{Re}\mu<1.
\tag{45b}
$$

结合两式，即可重新获得式 (27)。但是我们只能说上式与式 (27) 自洽，而不能把它当作式 (27) 的一种推导，否则便是循环论证了。下一节将给出式 (27) 的一种计算方法。

## §10 *特殊积分

前面几节研究的积分都有一个一般函数 $f(x)$ 或 $R(\cos\theta,\sin\theta)$，除了对其所有限制外，函数的形式可以相当一般。对于这些积分，我们导出了一般公式或给出了一般方法。采用相应的公式或方法就可以计算一类（而不是一个）积分。如果一个积分属于上述几种类型之一，我们就可以套用相应的公式或方法来计算，主要的工作就是找出有关的孤立奇点并计算其留数。对于这种积分，我们可以通俗地说其计算是有套路的。前面我们曾经提到，可以用围线积分和留数定理的方法来计算的积分是五花八门的。但是这其中有许多是没有套路可循的，只能对于每一个具体积分选取适当的辅助函数和积分围线，这对于初学者无疑是比较困难的。对付此类积分，主要凭借经验，只能靠学习模仿许多做练习逐步积累，不存在什么立竿见影的捷径。本节只介绍两个例子。例 1 是证明上节使用过但未加证明的式 (27)，并导出作为其特例的 Fresnel 积分公式。例 2 给出一个重要的积分公式，该公式在数理方程部分有重要的应用。

**例 1** 证明式 (27) 并导出 Fresnel 积分公式。

我们只需证明 $q=1$ 的情况，一般情况可以通过简单的变量置换得出。先证明式 (27) 中取正号的情况。取辅助函数为 $F(z)=z^{\mu-1}e^{iz}$。又取图 7 所示的围线 $C$，它位于第一象限，由以下各段构成：四分之一大圆弧 $\Gamma_R$；四分之一小圆弧 $\Gamma_r$；实轴上线段 $L_1=[r,R]$ 和虚轴上相应的线段 $L_2$，后者由 $L_1$ 逆时针转动 $\pi/2$ 得到。根据围线 $C$ 与各弧段的正方向规定可知 $C=\Gamma_R+\Gamma_r^-+L_1+L_2^-$，默认 $L_1$ 上 $\arg z=0$。由于 $F(z)$ 在围线 $C$ 及其内部单值解析，根据 Cauchy 积分定理，取得

$$
\int_C F(z)\,\mathrm dz
=\int_{L_1}F(z)\,\mathrm dz
-\int_{L_2}F(z)\,\mathrm dz
+\int_{\Gamma_R}F(z)\,\mathrm dz
-\int_{\Gamma_r}F(z)\,\mathrm dz
=0.
\tag{46}
$$

![Figure 7](...)

图 7：§10 例 1 的积分围线

根据辐角规定，易得

$$
\int_{L_1}F(z)\,\mathrm dz
=\int_r^R x^{\mu-1}e^{ix}\,\mathrm dx
\to\int_0^\infty x^{\mu-1}e^{ix}\,\mathrm dx,
\qquad r\to0,\ R\to\infty,
$$

$$
\int_{L_2}F(z)\,\mathrm dz
=e^{i\mu\pi/2}\int_r^R x^{\mu-1}e^{-x}\,\mathrm dx
\to e^{i\mu\pi/2}\int_0^\infty x^{\mu-1}e^{-x}\,\mathrm dx
=e^{i\mu\pi/2}\Gamma(\mu),
\qquad r\to0,\ R\to\infty.
$$

由于按已知条件积分存在，故取极限没有问题。根据已知条件和小圆弧引理，易得

$$
\int_{\Gamma_r}F(z)\,\mathrm dz\to0,
\qquad r\to0.
$$

对于沿 $\Gamma_R$ 的积分，由已知条件和证明 Jordan 引理的方法，易得

$$
\int_{\Gamma_R}F(z)\,\mathrm dz\to0,
\qquad R\to\infty.
$$

对式 (46) 取 $r\to0$、$R\to\infty$ 的极限，综合以上诸式，即得式 (27) 中正号情况的结果。取辅助函数 $F(z)=z^{\mu-1}e^{-iz}$ 并取位于第四象限内的类似围线，即可证明式 (27) 中负号情况。将两式相加和相减，即可获得式 (45) 的结果。

**思考** ① 式 (27) 中正号情况的式下取实部和虚部，是否就可以分别得出式 (45a) 和式 (45b)？② 为什么式 (45b) 的参数范围可以推广到 $-1<\operatorname{Re}\mu<1$？

在式 (45) 中令 $q=1$，$\mu=1/2$，利用 $\Gamma(1/2)=\sqrt\pi$，即可得到

$$
\int_0^\infty\frac{\cos x}{\sqrt x}\,\mathrm dx
=\int_0^\infty\frac{\sin x}{\sqrt x}\,\mathrm dx
=\sqrt{\frac\pi2}.
\tag{47}
$$

在上式中令 $t=\sqrt x$，最后再将积分变量 $t$ 改写为 $x$，我们就得到著名的 Fresnel 积分公式：

$$
\int_0^\infty\cos x^2\,\mathrm dx
=\int_0^\infty\sin x^2\,\mathrm dx
=\sqrt{\frac\pi8}.
\tag{48}
$$

计算 Fresnel 积分的另一种方法是取辅助函数 $F(z)=e^{iz^2}$，并选取由以下各弧段构成的围线 $C$：实轴上 $0\le x\le R$ 的线段 $L_1$；第一象限角平分线上的线段 $L_2:y=x$，$0\le x\le R$；八分之一圆弧 $\Gamma_R:z=Re^{i\theta}$，$0\le\theta\le\pi/4$。显然 $F(z)$ 沿围线 $C$ 的积分为零。采用证明 Jordan 引理的同样方法可以证明当 $R\to\infty$ 时圆弧 $\Gamma_R$ 上的积分极限为零。写出另外两段积分的形式并取 $R\to\infty$ 的极限，即可得出结果。读者可以尝试补充计算细节。这是多数教科书采用的方法。前面的方法可以得出更一般的结果，所以具有一定的优越性。比如作变量置换 $t=x^\alpha$ 并利用式 (45)，不难得出

$$
\int_0^\infty\cos x^\alpha\,\mathrm dx
=\frac1\alpha\Gamma\!\left(\frac1\alpha\right)
\cos\frac\pi{2\alpha},
\qquad \alpha>1,
\tag{49a}
$$

$$
\int_0^\infty\sin x^\alpha\,\mathrm dx
=\frac1\alpha\Gamma\!\left(\frac1\alpha\right)
\sin\frac\pi{2\alpha},
\qquad \alpha>1.
\tag{49b}
$$

这可以看作 Fresnel 积分的推广；当 $\alpha=2$，利用 $\Gamma(1/2)=\sqrt\pi$，由上式即可重新获得 Fresnel 积分。但是因为该方法涉及多值函数，对初学者来说难度略大。

**例 2** 试导出下面的积分公式

$$
I=\int_0^\infty e^{-ax^2}\cos bx\,\mathrm dx
=\frac12\sqrt{\frac\pi a}\,e^{-b^2/4a},
\qquad a>0,\quad b\in\mathbb R.
\tag{50}
$$

上式两边都是 $b$ 的偶函数，所以只需考虑 $b>0$。又因为 $b=0$ 的结果是微积分中熟知的，故只需考虑 $b>0$。将 $\cos bx=(e^{ibx}+e^{-ibx})/2$ 代入上式，对于第一项，作变量置换 $x'=-x$，再将 $x'$ 改写为 $x$，与第二项合并，即可将积分改写为

$$
I=\frac12\int_{-\infty}^{\infty}e^{-ax^2-ibx}\,\mathrm dx.
$$

对指数配方，又可得

$$
I=\frac12e^{-b^2/4a}
\int_{-\infty}^{\infty}e^{-a(x+ic)^2}\,\mathrm dx,
\tag{51}
$$

其中 $c=b/2a>0$。这里不能作变量置换 $x'=x+ic$ 消去指数上的常数项，因为该项是虚数，但是我们可以用围线积分的方法消去这个常数项。具体来说，取辅助函数 $F(z)=e^{-az^2}$ 以及图 8 的围线 $ABCDA$，由于 $F(z)$ 在围线及其内部单值解析，根据 Cauchy 积分定理，即得

$$
\int_{AB}F(z)\,\mathrm dz
-\int_{DC}F(z)\,\mathrm dz
+\int_{BC}F(z)\,\mathrm dz
-\int_{AD}F(z)\,\mathrm dz
=0.
\tag{52}
$$

![Figure 8](...)

图 8：§10 例 2 的积分围线

在 $AB$ 段，$z=x$，在 $DC$ 段，$z=x+ic$，故有

$$
\int_{AB}F(z)\,\mathrm dz
=\int_{-R}^{R}e^{-ax^2}\,\mathrm dx
\to\int_{-\infty}^{\infty}e^{-ax^2}\,\mathrm dx
=\sqrt{\frac\pi a},
\qquad R\to\infty,
$$

$$
\int_{DC}F(z)\,\mathrm dz
=\int_{-R}^{R}e^{-a(x+ic)^2}\,\mathrm dx
\to\int_{-\infty}^{\infty}e^{-a(x+ic)^2}\,\mathrm dx,
\qquad R\to\infty.
$$

在 $BC$ 段，$z=R+iy$，故有

$$
\left|\int_{BC}F(z)\,\mathrm dz\right|
=\left|i\int_0^c e^{a(y^2-R^2)}e^{-i2aRy}\,\mathrm dy\right|
\le ce^{a(c^2-R^2)}\to0,
\qquad R\to\infty.
$$

类似可证

$$
\int_{AD}F(z)\,\mathrm dz\to0,
\qquad R\to\infty.
$$

对式 (52) 取 $R\to\infty$ 的极限，综合以上诸式，即得

$$
\int_{-\infty}^{\infty}e^{-a(x+ic)^2}\,\mathrm dx
=\int_{-\infty}^{\infty}e^{-ax^2}\,\mathrm dx
=\sqrt{\frac\pi a}.
$$

上式的效果类似于通过平移消去了指数上的常数项。代入式 (51)，即得式 (50)。这个结果其实可以推广到 $a,b\in\mathbb C$、$\arg a\in(-\pi/2,\pi/2)$ 的情况，证明从略。

## 补充习题

计算下列积分。

1.
$$
\int_0^{2\pi}\frac{\sin^2\phi}{a+b\cos\phi}\,\mathrm d\phi,
\qquad a>b>0.
$$
答案：
$$
2\pi\frac{a-\sqrt{a^2-b^2}}{b^2}.
$$

2.
$$
\int_0^\pi\frac{a}{a^2+\sin^2\phi}\,\mathrm d\phi,
\qquad a>0.
$$
答案：
$$
\frac\pi{\sqrt{a^2+1}}.
$$
参考 §5 例 2。

3.
$$
\int_0^{\pi/2}\frac{1}{1+\cos^2\phi}\,\mathrm d\phi.
$$
答案：
$$
\frac\pi{2\sqrt2}.
$$
参考 §5 例 2。

4.
$$
\int_0^\infty\frac{1}{(x^2+a^2)^2(x^2+b^2)}\,\mathrm dx,
\qquad a>0,\ b>0.
$$
答案：
$$
\frac{\pi}{(a^2-b^2)^2}\,
\frac{2a^3+b^3-3a^2b}{4a^3b}.
$$

5.
$$
\int_0^\infty\frac{1}{x^4+a^4}\,\mathrm dx,
\qquad a>0.
$$
答案：
$$
\frac\pi{2\sqrt2\,a^3}.
$$
参考 §6 例 1。

6.
$$
\int_0^{+\infty}\frac{x\sin x}{x^2+1}\,\mathrm dx.
$$
答案：
$$
\frac\pi{2e}.
$$
参考 §7 例 2。

7.
$$
\int_0^\infty\frac{\cos x}{(x^2+a^2)(x^2+b^2)}\,\mathrm dx,
\qquad a>0,\ b>0.
$$
答案：
$$
\frac\pi{2(b^2-a^2)}
\left(\frac{e^{-a}}a-\frac{e^{-b}}b\right).
$$

8.
$$
\int_0^\infty\frac{\sin^2x}{x^2}\,\mathrm dx.
$$
答案：
$$
\frac\pi2.
$$

9.
$$
\int_0^\infty\frac{\sin mx}{x(x^2+a^2)}\,\mathrm dx,
\qquad m>0,\ a>0.
$$
答案：
$$
\frac\pi{2a^2}(1-e^{-ma}).
$$

10.
$$
\int_{|z|=2}\frac{e^z}{\cosh z}\,\mathrm dz.
$$
答案：
$$
4\pi i.
$$

11.
$$
\int_0^\infty\frac{1}{x^2}(\cos px-\cos qx)\,\mathrm dx,
\qquad p,q>0.
$$
答案：
$$
\frac\pi2(q-p).
$$

12.
$$
\int_0^\infty\frac{1}{x}(\cos px-\cos qx)\,\mathrm dx,
\qquad p,q>0.
$$
答案：
$$
\ln\frac qp.
$$
（本题较难）

13.
$$
\int_0^\infty\frac{x\sin mx}{(x^2+a^2)(x^2+b^2)}\,\mathrm dx,
\qquad m,a,b>0,\quad a\ne b.
$$
