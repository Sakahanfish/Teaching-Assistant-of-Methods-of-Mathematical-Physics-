# 第十章 $\delta$ 函数

本章将介绍一种新的“函数”——$\delta$ 函数。它是由英国物理学家 Dirac 首先引进的，可用于描写物理学中的点量，例如质点、点电荷、脉冲等，在近代物理学中有着广泛的应用。在数学运算中，按照一定的规则，$\delta$ 函数可以当作普通连续函数一样进行运算，如计算微分和积分，甚至应用于求解微分方程。总之，引入 $\delta$ 函数可以为我们处理有关的数学物理问题，带来极大的便利。

但是，$\delta$ 函数是一类“怪异”的函数，超越了经典的数学概念。它的严格数学理论（广义函数）要涉及泛函分析的知识。本章将从物理学的直观出发，引进 $\delta$ 函数的概念，介绍它的最基本的知识及其初步应用。

## §10.1 $\delta$ 函数的引入

作为 $\delta$ 函数的物理背景，先讨论点源，例如点电荷的密度分布函数的数学表示。

为简单起见，先讨论一维情形。如图 10.1 所示，设有总电量为 1 个单位的电荷，均匀分布在区间 $-l/2<x<l/2$ 内，区间外无电荷，则描述此电荷分布的电荷密度函数为

$$
\delta_l(x)=
\begin{cases}
0, & x\le -l/2,\\[4pt]
1/l, & -l/2<x<l/2,\\[4pt]
0, & x\ge l/2.
\end{cases}
\tag{10.1}
$$

![Figure 10.1](...)

图 10.1　单位点电荷的电荷密度

显然，

$$
\int_{-\infty}^{\infty}\delta_l(x)\,dx=1,
$$

即总电量为 1 个单位。对于在 $-l/2<x<l/2$ 中连续的任意函数 $f(x)$，根据中值定理，有

$$
\int_{-\infty}^{\infty}f(x)\delta_l(x)\,dx=f(\theta l),
\qquad
-\frac12\le\theta\le\frac12.
\tag{10.2}
$$

实际上，积分限不一定是 $\pm\infty$。只要 $a<-l/2,\ b>l/2$，就有

$$
\int_a^b f(x)\delta_l(x)\,dx=f(\theta l),
\qquad
-\frac12\le\theta\le\frac12.
\tag{10.3}
$$

作为极限情形，当 (10.1) 式中 $l\to0$ 时，就得到一维单位点电荷的电荷密度分布函数，记为

$$
\delta(x)=\lim_{l\to0}\delta_l(x)
=
\begin{cases}
0, & x\ne0,\\[4pt]
\infty, & x=0.
\end{cases}
\tag{10.4}
$$

而且，对于任意一个在 $x=0$ 点连续的函数 $f(x)$，有

$$
\int_{-\infty}^{\infty}f(x)\delta(x)\,dx=f(0).
\tag{10.5}
$$

这里的积分限也不一定是 $\pm\infty$。只要 $a<0,\ b>0$，就有

$$
\int_a^b f(x)\delta(x)\,dx=f(0).
\tag{10.6}
$$

总电量为 1 个单位这个条件始终未变，即

$$
\int_{-\infty}^{\infty}\delta(x)\,dx=1.
\tag{10.7}
$$

显然，上面关于点电荷密度分布函数 $\delta_l(x)$ 的描述也适用于其他物理量，例如质量的密度分布函数。重复上面的讨论，作为它们的极限情形，总会得到同样的结果①。

(10.5) 式左端的积分应该理解为

$$
\int_{-\infty}^{\infty}f(x)\delta(x)\,dx
=
\lim_{l\to0}
\int_{-\infty}^{\infty}f(x)\delta_l(x)\,dx.
$$

事实上，对于任意一个检验函数 $f(x)$，凡是具有

$$
\lim_{l\to0}
\int_{-\infty}^{\infty}f(x)\delta_l(x)\,dx
=
f(0)
\tag{10.8}
$$

性质的函数序列 $\{\delta_l(x)\}$，或是具有

$$
\lim_{n\to\infty}
\int_{-\infty}^{\infty}f(x)\delta_n(x)\,dx
=
f(0)
\tag{10.9}
$$

性质的函数序列 $\{\delta_n(x)\}$（例如见图 10.2），它们的极限都是 $\delta$ 函数。我们不妨把 $\delta$ 函数理解为满足 (10.9) 式的任意阶可微函数序列的极限。

$\delta$ 函数不像普通的函数那样具有唯一、确定的表达式。因为这样定义的函数，并不是经典意义下的函数：$\delta$ 函数并不给出普通的数值之间的对应关系，它所给出的对应关系

$$
\delta(x)=
\begin{cases}
0, & x\ne0,\\[4pt]
\infty, & x=0
\end{cases}
$$

按照对经典函数的理解是没有意义的，因为它唯一“有意义”的点是它唯一的奇点。它所给出的“函数值”只是在积分运算中才有所体现。从计算的角度来看，引进 $\delta$ 函数的初衷，即在于简化先对函数序列进行微积分计算，后取极限的过程。由于组成函数序列的函数具有足够好的连续性质，所以，在计算过程中可以把 $\delta$ 函数当作任意阶可微的函数处理，甚至可以定义 $\delta$ 函数的导数 $\delta'(x)$：对于任意一个检验函数 $f(x)$，有

$$
\int_{-\infty}^{\infty}f(x)\delta'(x)\,dx
=
\left.f(x)\delta(x)\right|_{-\infty}^{\infty}
-
\int_{-\infty}^{\infty}f'(x)\delta(x)\,dx
=
-f'(0).
\tag{10.10}
$$

① 关于 $\delta$ 函数的严格定义，需要使用广义函数的语言。出现在 (10.5) 和 (10.6) 式中的连续函数 $f(x)$ 称为检验函数。一维检验函数应当满足：(1) 在整个实轴上任意阶导数均存在，(2) 只在有界区间内函数值不为 0。

![Figure 10.2](...)

图 10.2　$\delta$ 函数的逼近序列举例

这里，就把 $\delta$ 函数当作普通的连续函数一样进行分部积分。

从电荷分布的简单图像或许有助于我们认识 $\delta'(x)$。设在 $x=\pm\varepsilon\ (\varepsilon>0)$ 两点分别有一点电荷，电荷量为 $\pm q\ (q>0)$，于是，电荷密度分布函数即为

$$
q\delta(x-\varepsilon)-q\delta(x+\varepsilon)
=
-2q\varepsilon\,
\frac{\delta(x+\varepsilon)-\delta(x-\varepsilon)}{2\varepsilon}.
$$

令 $\varepsilon\to0$ 而保持 $-2q\varepsilon=-p$（即这一对电荷构成电偶极子，其偶极矩大小为 $p$，方向与 $x$ 轴方向相反）为有限值，就得到电荷密度分布函数的极限值为 $-p\delta'(x)$。换句话说，$\delta'(x)$ 就是（位于 $x=0$ 处的）单位电偶极矩的电荷密度分布函数。

所有有关 $\delta$ 函数的等式，都应当从积分意义下去理解。如

$$
x\delta(x)=0
$$

应理解为

$$
\int_{-\infty}^{\infty}f(x)x\delta(x)\,dx=0,
\tag{10.11}
$$

$$
\delta(-x)=\delta(x)
$$

应理解为

$$
\int_{-\infty}^{\infty}f(x)\delta(-x)\,dx
=
\int_{-\infty}^{\infty}f(x)\delta(x)\,dx,
\tag{10.12}
$$

$$
\delta'(-x)=-\delta'(x)
$$

应理解为

$$
\int_{-\infty}^{\infty}f(x)\delta'(-x)\,dx
=
-\int_{-\infty}^{\infty}f(x)\delta'(x)\,dx,
\tag{10.13}
$$

$$
\delta(ax)=\frac1{|a|}\delta(x)
$$

应理解为

$$
\int_{-\infty}^{\infty}f(x)\delta(ax)\,dx
=
\int_{-\infty}^{\infty}
f(x)\left[\frac1{|a|}\delta(x)\right]dx,
\tag{10.14}
$$

$$
g(x)\delta(x)=g(0)\delta(x)
$$

应理解为

$$
\int_{-\infty}^{\infty}f(x)g(x)\delta(x)\,dx
=
\int_{-\infty}^{\infty}f(x)[g(0)\delta(x)]\,dx.
\tag{10.15}
$$

(10.15) 式中的 $g(x)$ 是经典的可微函数。

根据

$$
[g(x)\delta(x)]'
=
[g(0)\delta(x)]'
=
g(0)\delta'(x)
=
g'(x)\delta(x)+g(x)\delta'(x),
$$

可以导出等式

$$
g(x)\delta'(x)
=
g(0)\delta'(x)-g'(x)\delta(x),
$$

也就应该理解为

$$
\int_{-\infty}^{\infty}f(x)g(x)\delta'(x)\,dx
=
\int_{-\infty}^{\infty}
f(x)[g(0)\delta'(x)-g'(x)\delta(x)]\,dx.
\tag{10.16}
$$

由于

$$
\int_{-\infty}^{x}\delta(\xi)\,d\xi=\eta(x),
\tag{10.17}
$$

右端的 $\eta$ 函数就是 (8.4) 式中定义的单位阶跃函数，所以

$$
\delta(x)=\frac{d\eta(x)}{dx}.
\tag{10.18}
$$

这告诉我们，现在也可以对间断函数求导数，在间断点处就会出现 $\delta$ 函数。例如

$$
\frac{d e^{-|x|}}{dx}
=
-e^{-|x|}\operatorname{sgn}x,
\qquad
\frac{d^2e^{-|x|}}{dx^2}
=
e^{-|x|}-2\delta(x).
$$

$\delta$ 函数也可以表示成初等函数的 Fourier 积分。因为

$$
\int_{-\infty}^{\infty}\delta(x)e^{-ikx}\,dx=1,
$$

所以，根据 Fourier 变换的反演公式，有

$$
\delta(x)
=
\frac1{2\pi}
\int_{-\infty}^{\infty}e^{ikx}\,dk.
\tag{10.19}
$$

上式两边取实部，得

$$
\delta(x)
=
\frac1{2\pi}
\int_{-\infty}^{\infty}\cos kx\,dk.
\tag{10.20}
$$

显然 (10.19) 式和 (10.20) 式右端的积分在经典意义下都是不收敛的。这也从另一个角度表现出 $\delta$ 函数不是经典函数。

还可以对 $\delta$ 函数作 Laplace 变换：

$$
\delta(t-t_0)\doteq
\int_0^\infty
\delta(t-t_0)e^{-pt}\,dt
=
e^{-pt_0},
\qquad
t_0>0.
\tag{10.21}
$$

现在把 $\delta$ 函数推广到二维或三维的情形。显然，如果在平面上 $(x_0,y_0)$ 点处有一个单位点电荷，它的密度分布函数就是 $\delta(x-x_0)\delta(y-y_0)$。同样，如果在三维空间 $(x_0,y_0,z_0)$ 处有一个单位点电荷，它的密度分布函数就是 $\delta(x-x_0)\delta(y-y_0)\delta(z-z_0)$。当然，从三维空间来看，所谓一维点电荷应该是三维空间内的面电荷，二维点电荷就是三维空间内的线电荷。

**例 10.1**　证明

$$
\nabla^2\frac1r=-4\pi\delta(\boldsymbol r),
\tag{10.22}
$$

其中

$$
\nabla^2
\equiv
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}
\tag{10.23}
$$

称为 Laplace 算符，$r\equiv|\boldsymbol r|=\sqrt{x^2+y^2+z^2}$，$\delta(\boldsymbol r)=\delta(x)\delta(y)\delta(z)$。

**证**　正像前面指出的，凡是涉及 $\delta$ 函数的等式都应该从积分意义下去理解，这意味着本题就是应该去证明

$$
\iiint_V
\nabla^2\frac1r\,d\boldsymbol r
=
\begin{cases}
0, & \text{当 }\boldsymbol r=0\notin V,\\[4pt]
-4\pi, & \text{当 }\boldsymbol r=0\in V.
\end{cases}
\tag{10.24}
$$

其中 $d\boldsymbol r=dxdydz$ 是三维空间的体积元。当 $r\ne0$ 时，直接微商可得

$$
\frac{\partial^2}{\partial x^2}
\frac1{\sqrt{x^2+y^2+z^2}}
=
\frac{3x^2-(x^2+y^2+z^2)}
{(x^2+y^2+z^2)^{5/2}},
$$

$$
\frac{\partial^2}{\partial y^2}
\frac1{\sqrt{x^2+y^2+z^2}}
=
\frac{3y^2-(x^2+y^2+z^2)}
{(x^2+y^2+z^2)^{5/2}},
$$

$$
\frac{\partial^2}{\partial z^2}
\frac1{\sqrt{x^2+y^2+z^2}}
=
\frac{3z^2-(x^2+y^2+z^2)}
{(x^2+y^2+z^2)^{5/2}}.
$$

三式相加，即得

$$
\nabla^2\frac1r=0,
\qquad
r\ne0.
$$

这样就证得：当积分体积 $V$ 内不包含原点 $\boldsymbol r=0$ 时，积分恒为 0。

当积分体积 $V$ 内包含原点 $\boldsymbol r=0$ 时，由于函数 $1/r$ 在 $r=0$ 点不可导，上面的结果不成立。这时不妨将 $V$ 就取为整个（三维）空间，因而可以得到

$$
\begin{aligned}
\iiint \nabla^2\frac1r\,d^3\boldsymbol r
&=
\lim_{a\to0}
\iiint
\nabla^2
\frac1{\sqrt{x^2+y^2+z^2+a^2}}
\,dxdydz\\
&=
-\lim_{a\to0}
\iiint
\frac{3a^2}{(r^2+a^2)^{5/2}}
r^2\,dr\sin\theta\,d\theta\,d\phi\\
&=
-12\pi
\lim_{a\to0}
\int_0^\infty
\frac{a^2}{(r^2+a^2)^{5/2}}
r^2\,dr.
\end{aligned}
$$

容易看到，此处求极限和计算积分不能随意交换次序。但是作代换 $r=a\tan x$，即可证明上面的积分与 $a$ 无关，且

$$
\begin{aligned}
\iiint
\nabla^2\frac1r\,d^3\boldsymbol r
&=
-12\pi
\int_0^{\pi/2}
\frac{\tan^2x}{(1+\tan^2x)^{3/2}}\,dx\\
&=
-12\pi
\int_0^{\pi/2}
\sin^2x\cos x\,dx\\
&=
-12\pi\left.\frac13\sin^3x\right|_0^{\pi/2}
=
-4\pi.
\end{aligned}
$$

已知位于原点的单位点电荷产生的电势为

$$
u(\boldsymbol r)
=
\frac1{4\pi\varepsilon_0}\frac1r,
$$

位于原点的单位点电荷的电荷密度为 $\rho(\boldsymbol r)=\delta(\boldsymbol r)$，所以本题就是验证了点电荷情况下静电势所满足的方程

$$
\nabla^2u(\boldsymbol r)
=
-\frac1{\varepsilon_0}\rho(\boldsymbol r).
$$

## *§10.2 利用 $\delta$ 函数计算无穷积分

利用 $\delta$ 函数的常用积分表达式

$$
\delta(x)
=
\frac1{2\pi}
\int_{-\infty}^{\infty}e^{ikx}\,dk
\qquad
\text{或}
\qquad
\delta(x)
=
\frac1{2\pi}
\int_{-\infty}^{\infty}\cos kx\,dk,
$$

也可以计算无穷积分。下面通过几个例题来说明计算的一般步骤。

**例 10.2**　计算积分

$$
\int_{-\infty}^{\infty}\frac{\sin x}{x}\,dx.
$$

**解**　考虑辅助积分

$$
F(\lambda)
=
\int_{-\infty}^{\infty}
\frac{\sin\lambda x}{x}\,dx.
$$

如果许可在积分号下求导，则有

$$
F'(\lambda)
=
\int_{-\infty}^{\infty}
\cos\lambda x\,dx
=
2\pi\delta(\lambda).
$$

所以

$$
F(\lambda)
=
2\pi\eta(\lambda)+C,
$$

其中 $C$ 为积分常数，待定。不妨限定 $\lambda>0$，于是

$$
F(\lambda)=2\pi+C,
\qquad
F(-\lambda)=C.
$$

考虑到 $F(\lambda)$ 是 $\lambda$ 的奇函数，

$$
F(-\lambda)=-F(\lambda),
\qquad
F(0)=0,
$$

即可定出 $C=-\pi$。因此

$$
F(\lambda)
=
\begin{cases}
\pi, & \lambda>0,\\
0, & \lambda=0,\\
-\pi, & \lambda<0.
\end{cases}
$$

特别是，当 $\lambda=1$，就有

$$
\int_{-\infty}^{\infty}
\frac{\sin x}{x}\,dx
=
\pi.
$$

**例 10.3**　计算积分

$$
I
=
\int_{-\infty}^{\infty}
\frac{\sin2x}{x^2+x+1}\,dx.
$$

**解**　可以引进辅助积分

$$
F(\lambda)
=
\int_{-\infty}^{\infty}
\frac{e^{i\lambda x}}{x^2+x+1}\,dx,
$$

同样，在积分号下求导，就能得到它所满足的微分方程

$$
-F''(\lambda)-iF'(\lambda)+F(\lambda)
=
2\pi\delta(\lambda).
\tag{10.25}
$$

这是一个特殊的二阶常微分方程：其非齐次项含有 $\delta$ 函数。这种特殊性表现在两方面：一是当 $\lambda\ne0$ 时，$\delta(\lambda)=0$，方程是齐次的，所以

$$
F(\lambda)
=
\begin{cases}
Ae^{\lambda e^{-\pi i/6}}
+
Be^{\lambda e^{-5\pi i/6}},
& \lambda>0,\\[4pt]
Ce^{\lambda e^{-\pi i/6}}
+
De^{\lambda e^{-5\pi i/6}},
& \lambda<0.
\end{cases}
\tag{10.26}
$$

考虑到 $\lambda\to\pm\infty$ 时 $F(\lambda)$ 的有界性（因为积分收敛），$A$ 和 $D$ 必为 0。二是方程 (10.25) 的非齐次项为 $2\pi\delta(\lambda)$，这一点可以推断出 $F(\lambda)$ 在 $\lambda=0$ 点连续，

$$
\lim_{\varepsilon\to+0}
F(\lambda)\Big|_{0-\varepsilon}^{0+\varepsilon}
=
0,
\qquad
\text{由此导出}
\qquad
B=C,
\tag{10.27}
$$

而 $F'(\lambda)$ 不连续（因此 $F''(\lambda)$ 才会出现 $\delta(\lambda)$）。为了得到 $F'(\lambda)$ 在 $\lambda=0$ 点不连续性的定量描述，可以将微分方程 (10.25) 由 $\lambda=0$ 之左到 $\lambda=0$ 之右积分，于是就有

$$
\int_{0-\varepsilon}^{0+\varepsilon}
[F''(\lambda)+iF'(\lambda)-F(\lambda)]\,d\lambda
=
-2\pi
\int_{0-\varepsilon}^{0+\varepsilon}
\delta(\lambda)\,d\lambda
=
-2\pi.
$$

因 $F(\lambda)$ 在 $\lambda=0$ 点连续，故当 $\varepsilon\to+0$ 时，上式左端第二项和第三项的积分均趋于 0，

$$
\lim_{\varepsilon\to+0}
F'(\lambda)\Big|_{0-\varepsilon}^{0+\varepsilon}
=
-2\pi.
\tag{10.28}
$$

将 (10.26) 式代入 (10.28) 式，得

$$
\frac{-i-\sqrt3}{2}B
-
\frac{-i+\sqrt3}{2}C
=
-2\pi,
$$

因此求得

$$
B=C=\frac{2\pi}{\sqrt3}.
$$

所以

$$
F(\lambda)
=
\begin{cases}
\dfrac{2\pi}{\sqrt3}
e^{-\sqrt3\lambda/2}e^{-i\lambda/2},
& \lambda>0,\\[8pt]
\dfrac{2\pi}{\sqrt3}
e^{\sqrt3\lambda/2}e^{-i\lambda/2},
& \lambda<0.
\end{cases}
$$

而所要求的积分即为

$$
I
=
\operatorname{Im}F(2)
=
-\frac{2\pi}{\sqrt3}e^{-\sqrt3}\sin1.
$$

现在不妨再验证一下 $F'(\lambda)$ 在 $\lambda=0$ 点的不连续性。根据上面所得的结果，可以求出

$$
F'(\lambda)
=
\begin{cases}
-\left(1+\dfrac{i}{\sqrt3}\right)
\pi e^{-\sqrt3\lambda/2}e^{-i\lambda/2},
& \lambda>0,\\[8pt]
\left(1-\dfrac{i}{\sqrt3}\right)
\pi e^{\sqrt3\lambda/2}e^{-i\lambda/2},
& \lambda<0,
\end{cases}
$$

或者统一写成

$$
F'(\lambda)
=
\left[
1-\frac{i}{\sqrt3}-2\eta(\lambda)
\right]
\pi e^{-\sqrt3|\lambda|/2}e^{-i\lambda/2}.
$$

特别是，在 $\lambda=0$ 点，$F'(\lambda)$ 的左、右极限为

$$
\lim_{\varepsilon\to+0}
F'(0\pm\varepsilon)
=
-\frac{\pi i}{\sqrt3}\mp\pi.
$$

需要指出，例 10.3 中的 $F(\lambda)$ 连续而 $F'(\lambda)$ 不连续，这一结论具有普遍意义。为此，我们讨论一般的非齐次项为 $\delta$ 函数的二阶线性常微分方程。

任意一个二阶线性齐次常微分方程都可以整理成如下的形式：

$$
\frac{d}{dx}
\left[
p(x)\frac{dy(x)}{dx}
\right]
+
q(x)y(x)
=
0,
\qquad
a<x<b.
\tag{10.29}
$$

我们考虑非奇异的情形，即假设 $p(x),p'(x)$ 和 $q(x)$ 都是定义域 $a<x<b$ 上的实连续函数，$p(x)$ 无零点。因此，方程 (10.29) 的两个线性无关解 $y_1(x)$ 和 $y_2(x)$ 都在定义域 $a<x<b$ 上连续。而对非齐次项为 $\delta$ 函数的二阶线性常微分方程

$$
\frac{d}{dx}
\left[
p(x)\frac{dg(x;t)}{dx}
\right]
+
q(x)g(x;t)
=
\delta(x-t),
\qquad
a<x<b,\ a<t<b,
\tag{10.30}
$$

当 $x\ne t$ 时，方程 (10.30) 就是齐次方程 (10.29)，因此它的解可以用 $y_1(x)$ 和 $y_2(x)$ 叠加而得到：

$$
g(x;t)
=
\begin{cases}
c_1(t)y_1(x)+c_2(t)y_2(x)\equiv g_<,
& x<t,\\[4pt]
d_1(t)y_1(x)+d_2(t)y_2(x)\equiv g_>,
& x>t
\end{cases}
\tag{10.31}
$$

$$
=
g_<+(g_>-g_<)\eta(x-t).
\tag{10.32}
$$

其中 $c_1,c_2,d_1$ 和 $d_2$ 都是与 $x$ 无关的叠加系数。代回方程 (10.30)，得

$$
\begin{aligned}
&p(x)\left[
\frac{d^2g_<}{dx^2}
+
\left(
\frac{d^2g_>}{dx^2}
-
\frac{d^2g_<}{dx^2}
\right)\eta(x-t)
+
2\left(
\frac{dg_>}{dx}
-
\frac{dg_<}{dx}
\right)\delta(x-t)
+
(g_>-g_<)\delta'(x-t)
\right]\\
&\quad
+
p'(x)\left[
\frac{dg_<}{dx}
+
\left(
\frac{dg_>}{dx}
-
\frac{dg_<}{dx}
\right)\eta(x-t)
+
(g_>-g_<)\delta(x-t)
\right]\\
&\quad
+
q(x)\left[
g_<+(g_>-g_<)\eta(x-t)
\right]
=
\delta(x-t).
\end{aligned}
$$

由于 $y_1(x)$ 和 $y_2(x)$ 满足齐次线性方程 (10.29)，因此由它们叠加得到的 $g_<$ 和 $g_>$ 也都是齐次线性方程 (10.29) 的解，故上式可化简为

$$
\left[
p'(x)(g_>-g_<)
+
2p(x)
\left(
\frac{dg_>}{dx}
-
\frac{dg_<}{dx}
\right)
-
1
\right]\delta(x-t)
=
-p(x)(g_>-g_<)\delta'(x-t).
$$

此式应该在积分意义下理解，也就是说对任意的检验函数 $f(x)$，

$$
\begin{aligned}
&\int
\left[
p'(x)(g_>-g_<)
+
2p(x)
\left(
\frac{dg_>}{dx}
-
\frac{dg_<}{dx}
\right)
-
1
\right]
f(x)\delta(x-t)\,dx\\
&\qquad
=
-\int p(x)(g_>-g_<)f(x)\delta'(x-t)\,dx.
\end{aligned}
$$

运用 (10.5) 式和 (10.10) 式计算上述积分，整理，得

$$
\left[
p(t)
\left(
\frac{dg_>}{dx}
-
\frac{dg_<}{dx}
\right)_{x=t}
-
1
\right]f(t)
=
p(t)(g_>-g_<)_{x=t}f'(t).
$$

由于 $p(x)$ 无零点，函数 $f$ 的任意性，以及 $f$ 的函数值和它的一阶微商 $f'$ 的相互独立性，得

$$
(g_>-g_<)_{x=t}=0,
\qquad
\text{即}
\qquad
g(x,t)\Big|_{x=t-}^{x=t+}=0,
\tag{10.33a}
$$

$$
p(t)
\left(
\frac{dg_>}{dx}
-
\frac{dg_<}{dx}
\right)_{x=t}
-1=0,
\qquad
\text{即}
\qquad
\left.
\frac{dg(x,t)}{dx}
\right|_{x=t-}^{x=t+}
=
\frac1{p(t)}.
\tag{10.33b}
$$

所以，对于非齐次项为 $\delta(x-t)$ 的二阶线性常微分方程 (10.33)，其解 $g(x,t)$ 一定在 $x=t$ 点函数值连续而一阶微商不连续，跃度即为 $1/p(t)$。

## *§10.3 常微分方程初值问题的 Green 函数

上一节中将定积分的计算转化为求解一类特殊的常微分方程，非齐次项为 $\delta$ 函数的常微分方程。这类方程在理论上和实用上都具有特殊的重要性，有必要再做进一步的讨论。下面先讨论这类方程在齐次初始条件（当二阶常微分方程的两个定解条件都是关于微分变量 $t$ 同一个点，则称此定解条件为初始条件）下的解，以及解的某些基本性质。然后讨论如何运用 Green 函数方法求解常微分方程的初值问题。

在开始讨论之前，需要再次强调，第一，在传统意义下，这类非齐次项为 $\delta$ 函数的常微分方程是没有意义的。正像 $\delta$ 函数可以理解为连续函数序列（例如 $\{\delta_n(x)\},n=1,2,3,\cdots$）的极限一样，非齐次项为 $\delta$ 函数的微分方程也不妨理解为非齐次项为 $\delta_n(x)$ 的微分方程的极限，非齐次项为 $\delta$ 函数的微分方程的解也应当理解为非齐次项为 $\delta_n(x)$ 的微分方程的解的极限（先解微分方程再取极限）。引进 $\delta$ 函数的初衷就在于可以直接处理这种极限情形的微分方程求解问题，而不必考虑具体的函数序列以及它的极限过程。第二，正因为 $\delta$ 函数不是传统意义下的函数，这使得非齐次项为 $\delta$ 函数的微分方程的解具有独特的连续性质。就二阶常微分方程而言，正如我们在上一节中已经看到的，它的解是连续的，但是解的一阶导数是不连续的。第三，非齐次项为 $\delta$ 函数的常微分方程是一种特殊的非齐次方程：除了在使 $\delta$ 函数的宗量为零的个别点外，方程是齐次的。这使得这种非齐次常微分方程又很容易求解。在特殊情形下甚至可以直接积分而得到方程的通解。

下面从一个最简单的例子，开始关于常微分方程初值问题 Green 函数的讨论。

**例 10.4**　求解常微分方程初值问题

$$
\frac{d^2g}{dt^2}
=
\delta(t-\tau),
\qquad
t>0,\ \tau>0,
\tag{10.34a}
$$

$$
g\big|_{t=0}=0,
\qquad
\left.\frac{dg}{dt}\right|_{t=0}=0.
\tag{10.34b}
$$

**解**　显然 $g=g(t;\tau)$，即不仅依赖于变量 $t$，还依赖于参量 $\tau$。将方程 (10.34a) 直接积分，即得

$$
\frac{dg}{dt}
=
\eta(t-\tau)+\alpha(\tau).
\tag{10.35}
$$

将 (10.35) 式再积分一次，就得到

$$
g(t;\tau)
=
(t-\tau)\eta(t-\tau)
+
\alpha(\tau)t
+
\beta(\tau).
\tag{10.36}
$$

这里用到了 $g(t;\tau)$ 在 $t=\tau$ 点必须连续的要求，否则在 $dg/dt$ 中就会出现 $\delta$ 函数。

将解式 (10.36) 代入初值 (10.34b)，由此可以定出积分常数 $\alpha(\tau)=0,\ \beta(\tau)=0$，因此即可求得

$$
g(t;\tau)
=
(t-\tau)\eta(t-\tau).
\tag{10.37}
$$

![Figure 10.3](...)

图 10.3　例 10.4 的解 $g(t;\tau)$ 及其一阶导数

在常微分方程初值问题 (10.34) 的基础上，就可以求解常微分方程初值问题

$$
\frac{d^2y}{dt^2}=f(t),
\qquad
t>0,
\tag{10.38a}
$$

$$
y(0)=0,
\qquad
y'(0)=0.
\tag{10.38b}
$$

将方程 (10.38a) 的非齐次项 $f(t)$ 分解成

$$
f(t)
=
\int_0^\infty
f(\tau)\delta(t-\tau)\,d\tau,
$$

和 (10.34) 相比较，（至少在形式上）就可以得到解

$$
y(t)
=
\int_0^\infty
g(t;\tau)f(\tau)\,d\tau
=
\int_0^t
(t-\tau)f(\tau)\,d\tau.
\tag{10.39}
$$

再举一个略为复杂一点的例子。

**例 10.5**　求解常微分方程初值问题

$$
\frac{d^2g(t;\tau)}{dt^2}
+
k^2g(t;\tau)
=
\delta(t-\tau),
\qquad
t>0,\ \tau>0,
\tag{10.40a}
$$

$$
g(0;\tau)=0,
\qquad
\left.
\frac{dg(t;\tau)}{dt}
\right|_{t=0}
=
0.
\tag{10.40b}
$$

**解**　此方程和方程 (10.34a) 不同之处在于难以直接积分，故可介绍另一种解法。

注意到，当 $t\ne\tau$ 时，方程的非齐次项为 0，所以

$$
g(t;\tau)
=
\begin{cases}
A(\tau)\sin kt+B(\tau)\cos kt,
& t<\tau,\\[4pt]
C(\tau)\sin kt+D(\tau)\cos kt,
& t>\tau.
\end{cases}
$$

根据初始条件 (10.40b)，可以定出

$$
A(\tau)=0,
\qquad
B(\tau)=0.
$$

这里的 $g(t;\tau)$ 仍应具有 (10.33) 中相同的连续性质，即在 $t=\tau$ 点，

$$
g(t;\tau)\text{ 连续},
\qquad
g(t;\tau)\Big|_{\tau-0}^{\tau+0}=0;
$$

$$
\frac{dg(t;\tau)}{dt}\text{ 不连续},
\qquad
\left.
\frac{dg(t;\tau)}{dt}
\right|_{\tau-0}^{\tau+0}
=
1.
$$

于是

$$
C(\tau)\sin k\tau+D(\tau)\cos k\tau=0,
$$

$$
C(\tau)\cos k\tau-D(\tau)\sin k\tau=\frac1k.
$$

解之即得

$$
C(\tau)=\frac1k\cos k\tau,
\qquad
D(\tau)=-\frac1k\sin k\tau.
$$

这样，就得到

$$
g(t;\tau)
=
\frac1k
\sin k(t-\tau)\eta(t-\tau).
\tag{10.41}
$$

![Figure 10.4](...)

图 10.4　例 10.5 的解 $g(t;\tau)$ 及其一阶导数

在此基础上，也可以求出常微分方程初值问题

$$
\frac{d^2y}{dt^2}
+
k^2y(t)
=
f(t),
\qquad
t>0,
\tag{10.42a}
$$

$$
y(0)=0,
\qquad
y'(0)=0
\tag{10.42b}
$$

的解

$$
y(t)
=
\frac1k
\int_0^t
f(\tau)\sin k(t-\tau)\,d\tau.
\tag{10.43}
$$

这里顺便指出，(10.41) 式给出的当然是方程 (10.40a) 的一个特解，所以方程 (10.40a) 的通解就是

$$
g(t;\tau)
=
\frac1k
\sin k(t-\tau)\eta(t-\tau)
+
C(\tau)\sin kt
+
D(\tau)\cos kt.
\tag{10.44}
$$

从以上两个例子，可以提炼出一些更普遍的结论。为此不妨讨论更一般的二阶线性常微分方程初值问题

$$
\frac{d}{dt}
\left[
p(t)\frac{dg(t;\tau)}{dt}
\right]
+
q(t)g(t;\tau)
=
\delta(t-\tau),
\qquad
t>0,\ \tau>0;
\tag{10.45a}
$$

$$
g(0;\tau)=0,
\qquad
\left.
\frac{dg(t;\tau)}{dt}
\right|_{t=0}
=
0,
\tag{10.45b}
$$

其中 $p(t),p'(t)$ 和 $q(t)$ 都是区间 $[0,\infty)$ 上的实连续函数，$p(t)$ 无零点。

(1) $t<\tau$ 时 $g(t;\tau)\equiv0$。

由于当 $t<\tau$ 时，方程是齐次的，所以有通解

$$
g(t;\tau)
=
c_1(\tau)y_1(t)
+
c_2(\tau)y_2(t),
\qquad
t<\tau,
$$

其中 $y_1(t)$ 和 $y_2(t)$ 是齐次线性方程 (10.29) 的两个线性无关解。代入齐次初始条件，就得到方程组

$$
c_1(\tau)y_1(0)
+
c_2(\tau)y_2(0)
=
0,
$$

$$
c_1(\tau)y_1'(0)
+
c_2(\tau)y_2'(0)
=
0.
\tag{10.46}
$$

但由于 $y_1(t)$ 和 $y_2(t)$ 是线性无关的，

$$
W[y_1(t),y_2(t)]
\equiv
\begin{vmatrix}
y_1(t)&y_2(t)\\
y_1'(t)&y_2'(t)
\end{vmatrix}
\ne0,
$$

因此，作为 $c_1(\tau)$ 和 $c_2(\tau)$ 的线性齐次代数方程组 (10.46)，一定只有零解 $c_1(\tau)=0,\ c_2(\tau)=0$，即

$$
g(t;\tau)\equiv0,
\qquad
t<\tau.
$$

所以，上面的初始条件 (10.45b) 完全等价于

$$
g(t;\tau)\Big|_{t<\tau}=0,
\qquad
\left.
\frac{dg(t;\tau)}{dt}
\right|_{t<\tau}
=
0.
\tag{10.47}
$$

(2) $g(t;\tau)$ 在 $t=\tau$ 点的连续性质。在上一节 (10.33) 式已经看到 $g(t;\tau)$ 在 $t=\tau$ 点连续而 $dg(t;\tau)/dt$ 在 $t=\tau$ 点不连续，这是非齐次项为 $\delta$ 函数的二阶线性常微分方程解的共同特征，是由微分方程 (10.45a) 决定的。

(3) 可以通过 $y_1(t)$ 和 $y_2(t)$ 构造出常微分方程初值问题 (10.45) 的解 $g(t;\tau)$。

当 $t>\tau$ 时，方程 (10.45a) 也还是齐次的，故通解为

$$
g(t;\tau)
=
c_3(\tau)y_1(t)
+
c_4(\tau)y_2(t),
\qquad
t>\tau.
\tag{10.48}
$$

由连续性 (10.33a) 和 (10.33b) 式，有

$$
c_3(\tau)y_1(\tau)+c_4(\tau)y_2(\tau)=0,
$$

$$
c_3(\tau)y_1'(\tau)+c_4(\tau)y_2'(\tau)
=
\frac1{p(\tau)}.
$$

解之即得

$$
c_3(\tau)
=
-\frac1{p(\tau)}
\frac{y_2(\tau)}
{W[y_1(\tau),y_2(\tau)]},
$$

$$
c_4(\tau)
=
\frac1{p(\tau)}
\frac{y_1(\tau)}
{W[y_1(\tau),y_2(\tau)]}.
$$

代入 (10.48)，并和 $t<\tau$ 时的结果 (10.47) 结合起来，即可写出常微分方程初值问题 (10.45) 的解

$$
g(t;\tau)
=
\frac1{p(\tau)}
\frac{
y_1(\tau)y_2(t)-y_2(\tau)y_1(t)
}{
W[y_1(\tau),y_2(\tau)]
}
\eta(t-\tau).
\tag{10.49}
$$

(4) 在此基础上，可以进一步写出常微分方程初值问题

$$
\frac{d}{dt}
\left[
p(t)\frac{dy(t)}{dt}
\right]
+
q(t)y(t)
=
f(t),
\qquad
t>0,
\tag{10.50a}
$$

$$
y(0)=0,
\qquad
y'(0)=0
\tag{10.50b}
$$

的解

$$
y(t)
=
\int_0^t
g(t;\tau)f(\tau)\,d\tau.
\tag{10.51}
$$

不仅如此，还可以构造出更一般的常微分方程初值问题

$$
\frac{d}{dt}
\left[
p(t)\frac{dy(t)}{dt}
\right]
+
q(t)y(t)
=
f(t),
\qquad
t>0,
\tag{10.52a}
$$

$$
y(0)=A,
\qquad
y'(0)=B
\tag{10.52b}
$$

的解。

初值问题 (10.52) 代表了一般的非奇异二阶线性常微分方程初值问题，其中 $p(t),p'(t)$ 和 $q(t)$ 都是区间 $[0,\infty)$ 上的实连续函数，$p(t)$ 无零点。而且，根据物理上的实际需要，还假设 $p(t)=p(-t),\ q(t)=q(-t)$，这是为了保证算符

$$
\frac{d}{dt}
\left[
p(t)\frac{d}{dt}
\right]
+
q(t)
$$

在变换 $t\to-t$（“时间反演”）下保持不变。与初值问题 (10.52) 相应的 Green 函数 $g(t;\tau)$，它满足以 $\delta$ 函数为非齐次项的常微分方程初值问题 (10.45)，上面我们已经解出它的解 (10.49) 式。

为了将常微分方程的初值问题 (10.52) 的解用相应的 Green 函数 $g(t;\tau)$ 表示出来，先要证明常微分方程初值问题 Green 函数 $g(t;\tau)$ 的对称性

$$
g(t;\tau)
=
g(-\tau;-t).
\tag{10.53}
$$

**证**　首先需要注意，关系式 (10.53) 成立的前提是 $g(-\tau;-t)$ 在 $t,\tau>0$ 的条件下也有意义，这意味着必须要将定解问题 (10.45) 的成立条件延拓为 $-\infty<t,\tau<\infty$。其实这只要去掉 $t,\tau>0$ 的限制条件即可，于是，定解问题就变为①

$$
\frac{d}{dt}
\left[
p(t)\frac{dg(t;\tau)}{dt}
\right]
+
q(t)g(t;\tau)
=
\delta(t-\tau),
\qquad
-\infty<t,\tau<\infty,
\tag{10.45a'}
$$

$$
g(t;\tau)\Big|_{t<\tau}=0,
\qquad
\left.
\frac{dg(t;\tau)}{dt}
\right|_{t<\tau}
=
0.
\tag{10.45b'}
$$

从 Green 函数 $g(t;\tau)$ 的具体表达式 (10.49) 当然可以直接验证 (10.53) 式。一个简单的办法是注意在 $p(t)=p(-t)$ 和 $q(t)=q(-t)$ 的条件下，齐次方程 (10.29) 的两个线性无关解可以具有奇偶性，因此可取一个解为偶函数，另一个解为奇函数，例如

$$
y_1(-t)=y_1(t),
\qquad
y_2(-t)=-y_2(t).
$$

现在介绍另一种方法，它不涉及解的对称性，而是直接从微分方程及初始条件出发来证明。

已知 $g(t;\tau)$ 是常微分方程初值问题 (10.45') 的解。再引入 $g(-t;-t')$，它是定解问题

$$
\frac{d}{d(-t)}
\left[
p(-t)
\frac{dg(-t;-t')}{d(-t)}
\right]
+
q(-t)g(-t;-t')
=
\delta(-t+t'),
$$

$$
g(-t;-t')\Big|_{-t<-t'}=0,
\qquad
\left.
\frac{dg(-t;-t')}{dt}
\right|_{-t<-t'}
=
0
$$

即

$$
\frac{d}{dt}
\left[
p(t)\frac{dg(-t;-t')}{dt}
\right]
+
q(t)g(-t;-t')
=
\delta(t-t'),
\tag{10.54a}
$$

$$
g(-t;-t')\Big|_{t>t'}=0,
\qquad
\left.
\frac{dg(-t;-t')}{dt}
\right|_{t>t'}
=
0
\tag{10.54b}
$$

的解。将方程 (10.45a') 乘以 $g(-t;-t')$，同时将方程 (10.54a) 乘以 $g(t;\tau)$，相减，就得到

$$
\begin{aligned}
&g(-t;-t')
\frac{d}{dt}
\left[
p(t)\frac{dg(t;\tau)}{dt}
\right]
-
g(t;\tau)
\frac{d}{dt}
\left[
p(t)\frac{dg(-t;-t')}{dt}
\right]\\
&=
\frac{d}{dt}
\left\{
p(t)
\left[
g(-t;-t')
\frac{dg(t;\tau)}{dt}
-
g(t;\tau)
\frac{dg(-t;-t')}{dt}
\right]
\right\}\\
&=
g(-t;-t')\delta(t-\tau)
-
g(t;\tau)\delta(t-t').
\end{aligned}
$$

在定义域 $(-\infty,\infty)$ 上对 $t$ 积分，注意到初始条件 (10.45b') 和 (10.54b)，就有

$$
g(-\tau;-t')-g(t';\tau)
=
\left\{
p(t)
\left[
g(-t;-t')
\frac{dg(t;\tau)}{dt}
-
g(t;\tau)
\frac{dg(-t;-t')}{dt}
\right]
\right\}_{t=0}^{t=\infty}
=
0.
$$

将 $t'$ 改写成 $t$，即得到 (10.53) 式。$\square$

① 显然，在 $t,\tau>0$ 的条件下，定解问题 (10.45') 就退化为原始的定解问题 (10.45)。

下面将常微分方程初值问题 (10.52) 的解用相应的 Green 函数 $g(t;\tau)$ 表示出来。首先，根据 $g(-\tau;-t)$ 满足的定解问题 (10.54) 以及对称性 (10.53)，就能得到 $g(t;\tau)$ 作为 $\tau$ 的函数所满足的方程和初始条件

$$
\frac{d}{d\tau}
\left[
p(\tau)\frac{dg(t;\tau)}{d\tau}
\right]
+
q(\tau)g(t;\tau)
=
\delta(t-\tau),
\qquad
\tau>0,\ t>0,
\tag{10.55a}
$$

$$
g(t;\tau)\Big|_{\tau>t}=0,
\qquad
\left.
\frac{dg(t;\tau)}{d\tau}
\right|_{\tau>t}
=
0.
\tag{10.55b}
$$

然后将初值问题 (10.52) 中的自变量 $t$ 改写成 $\tau$，有

$$
\frac{d}{d\tau}
\left[
p(\tau)\frac{dy(\tau)}{d\tau}
\right]
+
q(\tau)y(\tau)
=
f(\tau),
\qquad
\tau>0,
\tag{10.56a}
$$

$$
y(0)=A,
\qquad
y'(0)=B.
\tag{10.56b}
$$

再将方程 (10.55a) 乘以 $y(\tau)$，同时将 (10.56a) 乘以 $g(t;\tau)$，相减，就有

$$
\begin{aligned}
&y(\tau)
\frac{d}{d\tau}
\left[
p(\tau)\frac{dg(t;\tau)}{d\tau}
\right]
-
g(t;\tau)
\frac{d}{d\tau}
\left[
p(\tau)\frac{dy(\tau)}{d\tau}
\right]\\
&=
\frac{d}{d\tau}
\left\{
p(\tau)
\left[
y(\tau)\frac{dg(t;\tau)}{d\tau}
-
g(t;\tau)\frac{dy(\tau)}{d\tau}
\right]
\right\}\\
&=
y(\tau)\delta(t-\tau)
-
g(t;\tau)f(\tau).
\end{aligned}
$$

最后，在 $\tau$ 的定义域 $[0,\infty)$ 上对 $\tau$ 积分，就能得到

$$
\begin{aligned}
y(t)
&=
\int_0^\infty
g(t;\tau)f(\tau)\,d\tau
+
\left\{
p(\tau)
\left[
y(\tau)\frac{dg(t;\tau)}{d\tau}
-
g(t;\tau)\frac{dy(\tau)}{d\tau}
\right]
\right\}_{\tau=0}^{\infty}\\
&=
\int_0^t
g(t;\tau)f(\tau)\,d\tau
-
p(0)
\left[
A\frac{dg(t;\tau)}{d\tau}
-
Bg(t;\tau)
\right]_{\tau=0}.
\end{aligned}
\tag{10.57}
$$

上述计算过程表明，如果 $g(t;\tau)$ 不满足齐次初始条件，最后便不能去掉 $\tau\to\infty$ 或 $\tau=t$ 的项，就无法完成用已知量（函数 $f(t)$ 与初值 $y(0),y'(0)$）及相应的 Green 函数 $g(t;\tau)$ 表示 $y(t)$ 的要求。


总结以上的讨论，读者可以看到，一旦我们求出了常微分方程初值问题 (10.45) 的解 $g(t;\tau)$，就一定能够构造出相应的常微分方程初值问题 (10.52) 的解 $y(t)$，不论它的非齐次项 $f(t)$ 是何形式或是有何初值 $y(0),y'(0)$。这样，这种特殊的非齐次（非齐次项为 $\delta$ 函数的）常微分方程初值问题的重要性，就不言自明了。初值问题 (10.50) 和 (10.52) 相应的 Green 函数满足定解问题 (10.45)，常微分方程初值问题 (10.45) 的解 $g(t;\tau)$ 称为常微分方程初值问题 (10.50) 和 (10.52) 相应的 Green 函数。**常微分方程初值问题相应的 Green 函数，满足的方程就是将原方程的非齐次项换成 $\delta$ 函数，而初始条件就是 (10.47) 式。**上面介绍的利用 Green 函数求解常微分方程初值问题 (10.52) 的方法——Green 函数方法，是求解非齐次常微分方程定解问题的重要方法。作为本章的继续，在数学物理方程部分中我们还要讨论求解偏微分方程定解问题的 Green 函数方法。

## *§10.4 常微分方程边值问题的 Green 函数

本节讨论另一种类型常微分方程定解问题的 Green 函数，即常微分方程边值问题的 Green 函数，它是非齐次项为 $\delta$ 函数的常微分方程在齐次边界条件下的解。之所以称为二阶常微分方程的边值问题而非初值问题，是因为两个定解条件给定的是未知函数在不同点（例如 $x=a$ 和 $x=b$）的函数值或微商值，或者它们的线性组合。最后，本节也要介绍求解常微分方程边值问题的 Green 函数方法。

先讨论几个具体的例子。

**例 10.6**　求解常微分方程边值问题

$$
\frac{d^2g(x;\xi)}{dx^2}
=
\delta(x-\xi),
\qquad
a<x,\xi<b,
\tag{10.58a}
$$

$$
g(a;\xi)=0,
\qquad
g(b;\xi)=0.
\tag{10.58b}
$$

**解**　这个问题和初值问题 (10.34) 的微分方程相同，但定解条件不同，现在要求的是微分方程 (10.58a) 在齐次边界条件 (10.58b) 下的解。

前面在例 10.4 中已经得到方程 (10.58a) 的通解 (10.36) 式，将边界条件 (10.58b) 代入，得

$$
a\alpha(\xi)+\beta(\xi)=0,
$$

$$
b-\xi+b\alpha(\xi)+\beta(\xi)=0.
$$

解之即得

$$
\alpha(\xi)
=
-\frac{b-\xi}{b-a},
\qquad
\beta(\xi)
=
\frac{a(b-\xi)}{b-a}.
$$

这就求出了常微分方程边值问题 (10.58) 的解，

$$
g(x;\xi)
=
(x-\xi)\eta(x-\xi)
-
\frac{b-\xi}{b-a}(x-a).
\tag{10.59}
$$

![Figure 10.5](...)

图 10.5　例 10.6 的解 $g(x;\xi)$ 及其一阶导数

和常微分方程初值问题 (10.34) 的解一样，常微分方程边值问题 (10.58) 的解，也具有下列连续性质：

$$
g(x;\xi)\text{ 在 }x=\xi\text{ 点连续},
\qquad
g(x;\xi)\Big|_{\xi-0}^{\xi+0}=0;
\tag{10.60a}
$$

$$
\frac{dg(x;\xi)}{dx}\text{ 在 }x=\xi\text{ 点不连续},
\qquad
\left.
\frac{dg(x;\xi)}{dx}
\right|_{\xi-0}^{\xi+0}
=
1.
\tag{10.60b}
$$

因为这种特性完全是由微分方程 (10.58a) 决定的。

**例 10.7**　求解常微分方程边值问题

$$
\frac{d^2g(x;\xi)}{dx^2}
+
k^2g(x;\xi)
=
\delta(x-\xi),
\qquad
a<x,\xi<b,
\tag{10.61a}
$$

$$
g(a;\xi)=0,
\qquad
g(b;\xi)=0.
\tag{10.61b}
$$

**解**　例 10.5 中已经给出了方程 (10.61a) 的通解 (10.44) 式。代入边界条件 (10.61b)，有

$$
C(\xi)\sin ka+D(\xi)\cos ka=0,
$$

$$
C(\xi)\sin kb+D(\xi)\cos kb
=
-\frac1k\sin k(b-\xi).
$$

解之即得

$$
C(\xi)
=
-\frac1k
\frac{\sin k(b-\xi)}
{\sin k(b-a)}
\cos ka,
$$

$$
D(\xi)
=
\frac1k
\frac{\sin k(b-\xi)}
{\sin k(b-a)}
\sin ka,
$$

所以

$$
g(x;\xi)
=
\frac1k\sin k(x-\xi)\eta(x-\xi)
-
\frac1k
\frac{\sin k(b-\xi)}
{\sin k(b-a)}
\sin k(x-a).
\tag{10.62}
$$

当 $k(b-a)=n\pi\ne0$ 时，此题无解。

不直接写出方程 (10.61a) 的通解，也能求解这个常微分方程的边值问题。注意当 $x\ne\xi$ 时，微分方程 (10.61a) 是齐次的。于是，当 $a<x<\xi$，要求 $g(a;\xi)=0$，故解为

$$
g(x;\xi)
=
C(\xi)\sin k(x-a),
\qquad
a<x<\xi;
$$

又当 $\xi<x<b$ 时，要求 $g(b;\xi)=0$，又有解

$$
g(x;\xi)
=
D(\xi)\sin k(b-x),
\qquad
\xi<x<b.
$$

利用 $g(x;\xi)$ 在 $x=\xi$ 点的连续性质 (10.60)，可以得到

$$
D(\xi)\sin k(b-\xi)
-
C(\xi)\sin k(\xi-a)
=
0,
$$

$$
D(\xi)\cos k(b-\xi)
+
C(\xi)\cos k(\xi-a)
=
-\frac1k.
$$

![Figure 10.6](...)

图 10.6　例 10.7 的解 $g(x;\xi)$ 及其一阶导数

这样就可以求得

$$
C(\xi)
=
-\frac1k
\frac{\sin k(b-\xi)}
{\sin k(b-a)},
\qquad
D(\xi)
=
-\frac1k
\frac{\sin k(\xi-a)}
{\sin k(b-a)}.
$$

最后，也得到 Green 函数

$$
g(x;\xi)
=
\begin{cases}
-\dfrac1k
\dfrac{\sin k(b-\xi)}
{\sin k(b-a)}
\sin k(x-a),
& a<x<\xi,\\[10pt]
-\dfrac1k
\dfrac{\sin k(\xi-a)}
{\sin k(b-a)}
\sin k(b-x),
& \xi<x<b.
\end{cases}
$$

不难把它化成 (10.62) 式的形式。

容易把这种解法推广到更一般的情形，即根据齐次常微分方程 (10.29) 的线性无关解 $y_1(x)$ 和 $y_2(x)$，去构造常微分方程边值问题

$$
\frac{d}{dx}
\left[
p(x)\frac{dg(x;\xi)}{dx}
\right]
+
q(x)g(x;\xi)
=
\delta(x-\xi),
\qquad
a<x,\xi<b,
\tag{10.63a}
$$

$$
g(a;\xi)=0,
\qquad
g(b;\xi)=0
\tag{10.63b}
$$

的解。方法是先在 $a<x<\xi$ 和 $\xi<x<b$ 的两个区间中写出方程 (10.63a) 在各自的齐次边界条件下的解，然后利用在 $x=\xi$ 点的连续性要求 (10.33) 定出叠加系数。这样便可得到结果

$$
\begin{aligned}
g(x;\xi)
={}&
-\frac1{p(\xi)}
\frac{
y_2(b)y_1(\xi)-y_1(b)y_2(\xi)
}{
y_1(b)y_2(a)-y_1(a)y_2(b)
}
\frac{
y_2(a)y_1(x)-y_1(a)y_2(x)
}{
W[y_1(\xi),y_2(\xi)]
}\\
&+
\frac1{p(\xi)}
\frac{
y_1(\xi)y_2(x)-y_2(\xi)y_1(x)
}{
W[y_1(\xi),y_2(\xi)]
}
\eta(x-\xi).
\end{aligned}
\tag{10.64}
$$

请读者补足这个计算。

需要说明，本节只讨论了一种类型的边界条件，即函数在端点的数值为 0。可以完全类似地讨论其他类型的齐次边界条件。普遍地说，**常微分方程边值问题相应的 Green 函数，满足的方程也是将原方程的非齐次项换成 $\delta$ 函数，而边界条件类型与原定解问题相同，但要一律改成为齐次的。**讨论这些常微分方程边值问题的 Green 函数，目的仍然是求解更一般的常微分方程边值问题。例如，对于常微分方程边值问题

$$
\frac{d}{dx}
\left[
p(x)\frac{dy(x)}{dx}
\right]
+
q(x)y(x)
=
f(x),
\qquad
a<x<b,
\tag{10.65a}
$$

$$
y(a)=A,
\qquad
y(b)=B,
\tag{10.65b}
$$

它所对应的 Green 函数 $g(x;\xi)$ 就满足常微分方程边值问题 (10.63)。问题仍然是要求一定能用 (10.63) 的解 $g(x;\xi)$ 构造出 (10.65) 的解。要圆满地解决这个问题，仍然先要证明常微分方程边值问题 Green 函数 $g(x;\xi)$ 的对称性

$$
g(x;\xi)
=
g(\xi;x).
\tag{10.66}
$$

证明的方法①仍然是再引入一个 Green 函数 $g(x;x')$，它是常微分方程边值问题

$$
\frac{d}{dx}
\left[
p(x)\frac{dg(x;x')}{dx}
\right]
+
q(x)g(x;x')
=
\delta(x-x'),
\qquad
a<x,x'<b,
\tag{10.67a}
$$

$$
g(a;x')=0,
\qquad
g(b;x')=0
\tag{10.67b}
$$

的解。将方程 (10.63a) 和 (10.67a) 分别乘以 $g(x;x')$ 和 $g(x;\xi)$，相减，再在区间 $[a,b]$ 上对 $x$ 积分，即得

$$
g(\xi;x')-g(x';\xi)
=
\left\{
p(x)
\left[
g(x;x')
\frac{dg(x;\xi)}{dx}
-
g(x;\xi)
\frac{dg(x;x')}{dx}
\right]
\right\}_{x=a}^{x=b}.
$$

① 当然，从 $g(x;\xi)$ 的具体表达式 (10.64) 就可以直接看出 (10.66) 式成立。

代入边界条件 (10.63b) 和 (10.67b)，得到上式右方为 0，所以，

$$
g(\xi;x')
=
g(x';\xi).
$$

将 $x'$ 换成 $x$，这正好就是 (10.66) 式。$\square$

下面我们把边值问题 (10.65) 的解用相应的 Green 函数 $g(x;\xi)$ 以及已知条件（非齐次项 $f(x)$ 以及边界条件 $y(a)$ 和 $y(b)$）构造出来，步骤与 Green 函数法求解初值问题时的步骤基本相同。重要的不同之处在于边值问题 Green 函数的对称性和初值问题不同，定解条件的形式也不相同。

首先写出相应的 Green 函数 $g(\xi;x)$ 满足的常微分方程边值问题

$$
\frac{d}{d\xi}
\left[
p(\xi)\frac{dg(\xi;x)}{d\xi}
\right]
+
q(\xi)g(\xi;x)
=
\delta(x-\xi),
\qquad
a<\xi,x<b,
$$

$$
g(\xi;x)\big|_{\xi=a}=0,
\qquad
g(\xi;x)\big|_{\xi=b}=0.
$$

然后利用 $g(x;\xi)$ 的对称性 (10.66)，将 $g(\xi;x)$ 换成 $g(x;\xi)$，就能得到 $g(x;\xi)$ 作为 $\xi$ 的函数所满足的常微分方程边值问题

$$
\frac{d}{d\xi}
\left[
p(\xi)\frac{dg(x;\xi)}{d\xi}
\right]
+
q(\xi)g(x;\xi)
=
\delta(x-\xi),
\qquad
a<\xi,x<b,
\tag{10.68a}
$$

$$
g(x;\xi)\big|_{\xi=a}=0,
\qquad
g(x;\xi)\big|_{\xi=b}=0.
\tag{10.68b}
$$

同时将待求解的边值问题 (10.65) 中的自变量改写成 $\xi$，有

$$
\frac{d}{d\xi}
\left[
p(\xi)\frac{dy(\xi)}{d\xi}
\right]
+
q(\xi)y(\xi)
=
f(\xi),
\qquad
\xi>0,
\tag{10.69a}
$$

$$
y(a)=A,
\qquad
y(b)=B.
\tag{10.69b}
$$

最后将方程 (10.69a) 乘以 $g(x;\xi)$，同时将 (10.68a) 乘以 $y(\xi)$，相减，再在区间 $[a,b]$ 上对 $\xi$ 积分，整理即得

$$
y(x)
=
\int_a^b
g(x;\xi)f(\xi)\,d\xi
+
\left\{
p(\xi)
\left[
y(\xi)\frac{dg(x;\xi)}{d\xi}
-
g(x;\xi)\frac{dy(\xi)}{d\xi}
\right]
\right\}_{\xi=a}^{\xi=b}
\tag{10.70}
$$

$$
=
\int_a^b
g(x;\xi)f(\xi)\,d\xi
+
Bp(b)
\left.
\frac{dg(x;\xi)}{d\xi}
\right|_{\xi=b}
-
Ap(a)
\left.
\frac{dg(x;\xi)}{d\xi}
\right|_{\xi=a}.
\tag{10.71}
$$

在最后一步，利用了边界条件 (10.69b) 和 (10.68b)。因此 (10.70) 式适用于方程为 (10.69a) 式的所有不同类型的边值问题，而 (10.71) 式只适用于边值问题 (10.69)。

最后举一个无界区间的例子。

**例 10.8**　用 Green 函数方法求解无界区间内的常微分方程边值问题

$$
\frac{d^2y(x)}{dx^2}
-
k^2y(x)
=
f(x),
\qquad
-\infty<x<\infty,
\tag{10.72a}
$$

$$
y(x)\big|_{x\to\pm\infty}\text{ 有界},
\tag{10.72b}
$$

其中 $k>0$。

**解**　为此，先求出此问题相应的 Green 函数，即边值问题

$$
\frac{d^2g(x;\xi)}{dx^2}
-
k^2g(x;\xi)
=
\delta(x-\xi),
\qquad
-\infty<x,\xi<\infty,
\tag{10.73a}
$$

$$
g(x;\xi)\big|_{x\to\pm\infty}\text{ 有界}
\tag{10.73b}
$$

的解

$$
g(x;\xi)
=
-\frac1{2k}
e^{-k|x-\xi|}.
\tag{10.74}
$$

![Figure 10.7](...)

图 10.7　例 10.8 的解 $g(x;\xi)$ 及其一阶导数

再模仿有界区间的做法，可以求得 $y(x)$ 和 $g(x;\xi)$ 之间的关系

$$
y(x)
=
\int_{-\infty}^{\infty}
g(x;\xi)f(\xi)\,d\xi
+
\left[
y(\xi)\frac{dg(x;\xi)}{d\xi}
-
g(x;\xi)\frac{dy(\xi)}{d\xi}
\right]_{\xi\to-\infty}^{\xi\to+\infty}.
$$

代入 (10.74) 式以及边界条件 (10.72b) 和 (10.73b)，就得到

$$
y(x)
=
-\frac1{2k}
\int_{-\infty}^{\infty}
e^{-k|x-\xi|}f(\xi)\,d\xi.
\tag{10.75}
$$

这个例子告诉我们，如果待求解的定解问题边界上是有界条件，那么相应的 Green 函数满足的边界条件同样也是有界条件。

以上的讨论还告诉我们，虽然求解微分方程定解问题的 Green 函数法步骤稍显烦琐，但是对每一类定解问题（定解问题中方程和定解条件等号左边都完全相同地属于同一类）只需求解一个 Green 函数，然后按照固定的步骤就能解出这一类中的所有定解问题。而且最后得到的解式还便于我们分析解对于非齐次项以及边界条件等的依赖关系。

## 习题

1. 证明 $\delta$ 函数的下列性质：

   (1) $\delta(x)=\delta(-x)$；

   (2) $x\delta(x)=0$；

   (3) $g(x)\delta(x)=g(0)\delta(x)$；

   (4) $x\delta'(x)=-\delta(x)$；

   (5)
   $$
   \delta(ax)=\frac1a\delta(x),
   \qquad
   a>0;
   $$

   (6)
   $$
   g(x)\delta'(x)
   =
   g(0)\delta'(x)-g'(0)\delta(x);
   $$

   (7)
   $$
   \delta(x^2-a^2)
   =
   \frac1{2a}
   [\delta(x-a)+\delta(x+a)],
   \qquad
   a>0.
   $$

2. 求下列常微分方程初值问题的解：

   (1)

   $$
   \left(
   \frac{d^2}{dt^2}-k^2
   \right)
   g(t;\tau)
   =
   \delta(t-\tau),
   \qquad
   t,\tau>0,\ k>0,
   $$

   $$
   g(0;\tau)=0,
   \qquad
   \left.
   \frac{dg(t;\tau)}{dt}
   \right|_{t=0}
   =
   0;
   $$

   (2)

   $$
   \left(
   \frac{d^2}{dt^2}-t^2
   \right)
   g(t;\tau)
   =
   \delta(t-\tau),
   \qquad
   t,\tau>0,
   $$

   $$
   g(0;\tau)=0,
   \qquad
   \left.
   \frac{dg(t;\tau)}{dt}
   \right|_{t=0}
   =
   0;
   $$

   (3)

   $$
   \left[
   (1+t+t^2)\frac{d^2}{dt^2}
   +
   2(1+2t)\frac{d}{dt}
   +
   2
   \right]
   g(t;\tau)
   =
   \delta(t-\tau),
   \qquad
   t,\tau>0,
   $$

   $$
   g(0;\tau)=0,
   \qquad
   \left.
   \frac{dg(t;\tau)}{dt}
   \right|_{t=0}
   =
   0.
   $$

   提示：以上各变系数常微分方程的解见第 9 章习题。

3. 用三种方法（常数变易法、Laplace 变换和 Green 函数方法）求解常微分方程初值问题：

$$
\begin{cases}
\dfrac{d^2y(t)}{dt^2}+k^2y(t)=f(t),
& t>0,\ k>0,\\[6pt]
y(0)=A,
\qquad
\left.\dfrac{dy(t)}{dt}\right|_{t=0}=B.
\end{cases}
$$

4. 用 Green 函数方法求解下列常微分方程初值问题：

   (1)

   $$
   \frac{d^2y(t)}{dt^2}
   -
   k^2y(t)
   =
   f(t),
   \qquad
   t>0,\ k>0,
   $$

   $$
   y(0)=A,
   \qquad
   \left.
   \frac{dy(t)}{dt}
   \right|_{t=0}
   =
   B;
   $$

   (2)

   $$
   \frac{d^2y(t)}{dt^2}
   -
   t^2y(t)
   =
   f(t),
   \qquad
   t>0,
   $$

   $$
   y(0)=A,
   \qquad
   \left.
   \frac{dy(t)}{dt}
   \right|_{t=0}
   =
   B.
   $$

5. 求下列常微分方程边值问题的解：

   (1)

   $$
   \left(
   \frac{d^2}{dx^2}-k^2
   \right)
   g(x;\xi)
   =
   \delta(x-\xi),
   \qquad
   0<x,\xi<1,\ k>0,
   $$

   $$
   g(0;\xi)=0,
   \qquad
   g(1;\xi)=0;
   $$

   (2)

   $$
   \left(
   \frac{d^2}{dx^2}-x^2
   \right)
   g(x;\xi)
   =
   \delta(x-\xi),
   \qquad
   0<x,\xi<1,
   $$

   $$
   g(0;\xi)=0,
   \qquad
   g(1;\xi)=0;
   $$

   (3)

   $$
   \left[
   (1+x+x^2)\frac{d^2}{dx^2}
   +
   2(1+2x)\frac{d}{dx}
   +
   2
   \right]
   g(x;\xi)
   =
   \delta(x-\xi),
   \qquad
   0<x,\xi<l<1,
   $$

   $$
   g(0;\xi)=0,
   \qquad
   g(l;\xi)=0.
   $$

6. 利用 Laplace 变换，求解常微分方程边值问题：

$$
\begin{cases}
\dfrac{d^2y(x)}{dx^2}-k^2y(x)=\delta(x-\xi),
& x>0,\ k>0,\ \xi>0,\\[6pt]
y(0)=0,
\qquad
\displaystyle\lim_{x\to\infty}y(x)=0.
\end{cases}
$$

7. 求下列常微分方程边值问题相应的 Green 函数：

   (1)

   $$
   \begin{cases}
   \dfrac{d^2y(x)}{dx^2}-k^2y(x)=f(x),
   & x>0,\ k>0,\\[6pt]
   \left.\dfrac{dy(x)}{dx}\right|_{x=0}=v_0,
   \qquad
   y(x)\big|_{x\to\infty}\text{ 有界},
   \end{cases}
   $$

   其中 $v_0$ 是已知常数。

   (2)

   $$
   \begin{cases}
   \dfrac{d^2y(x)}{dx^2}+k^2y(x)=f(x),
   & a<x<b,\ k>0,\\[6pt]
   y(a)=y_0,
   \qquad
   \left.\dfrac{dy(x)}{dx}\right|_{x=b}=v_0,
   \end{cases}
   $$

   其中 $y_0$ 和 $v_0$ 是已知常数。

8. 用 Green 函数方法求解下列常微分方程边值问题：

   (1)

   $$
   \frac{d^2y(x)}{dx^2}
   +
   k^2y(x)
   =
   f(x),
   \qquad
   0<x<1,
   $$

   $$
   y(0)=A,
   \qquad
   y(1)=B;
   $$

   (2)

   $$
   \frac{d^2y(x)}{dx^2}
   -
   k^2y(x)
   =
   f(x),
   \qquad
   0<x<1,\ k>0,
   $$

   $$
   y(0)=A,
   \qquad
   y(1)=B;
   $$

   (3)

   $$
   \frac{d^2y(x)}{dx^2}
   -
   x^2y(x)
   =
   f(x),
   \qquad
   0<x<1,
   $$

   $$
   y(0)=A,
   \qquad
   y(1)=B.
   $$

### 间断函数的导数

重新考察一下 (10.18) 式，它其实就是给出了单位阶跃函数 $\eta(x)$ 的导数公式。严格说来，它是在广义函数的意义下成立。因此，可以预料，在广义函数的意义下，对于具有第一类间断点的任意函数 $f(x)$，也可以类似地求得它们的导数 $f'$ 乃至任意阶导数。例如，因为

$$
\operatorname{sgn}x
=
-1+2\eta(x),
$$

我们就有

$$
\frac{d\,\operatorname{sgn}x}{dx}
=
2\delta(x).
$$

再来计算函数 $e^{-|x|}$ 的各阶导数。根据微积分学知识我们可以判断，函数 $e^{-|x|}$ 连续，除 $x=0$ 之外，处处可导，而且在 $x=0$ 点右侧导数与左侧导数分别为 $\mp1$。这样，我们就得到

$$
\frac{d(e^{-|x|})}{dx}
=
-e^{-|x|}\operatorname{sgn}x.
$$

此函数分段连续，有第一类间断点 $x=0$，因此，在引进 $\delta$ 函数后，我们就有

$$
\begin{aligned}
\frac{d^2(e^{-|x|})}{dx^2}
&=
-\frac{d(e^{-|x|})}{dx}\operatorname{sgn}x
-
e^{-|x|}
\frac{d(\operatorname{sgn}x)}{dx}\\
&=
e^{-|x|}(-\operatorname{sgn}x)^2
-
2e^{-|x|}\delta(x)\\
&=
e^{-|x|}-2\delta(x).
\end{aligned}
$$

其中用到了 $f(x)\delta(x)=f(0)\delta(x)$。

如此继续，还能求得 $e^{-|x|}$ 的高阶导数。例如，

$$
\frac{d^3(e^{-|x|})}{dx^3}
=
-e^{-|x|}\operatorname{sgn}x
-
2\delta'(x),
$$

$$
\frac{d^4(e^{-|x|})}{dx^4}
=
e^{-|x|}
-
2\delta(x)
-
2\delta''(x).
$$
