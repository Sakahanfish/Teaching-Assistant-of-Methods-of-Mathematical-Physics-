# 第十八章 积分变换的应用

本章介绍求解偏微分方程定解问题的另一种方法——将积分变换应用于求解偏微分方程定解问题。常用的积分变换有 Laplace 变换和 Fourier 变换两种。

## §18.1 Laplace 变换

Laplace 变换可用于求解含时间的偏微分方程定解问题。对于系数与 $t$ 无关的偏微分方程，变换后自变量的个数比原来减少一个。例如，原来是 $x$ 和 $t$ 两个自变量的偏微分方程定解问题，变换后就是自变量只有一个 $x$ 的常微分方程定解问题。一般说来，后者总比前者更容易求解。当然，这样求得的是原始定解问题解的像函数，还必须反演才能得到原始问题的解。

下面先举一个无界杆的热传导问题的例子。

**例 18.1　求解无界杆的热传导问题**

$$
\frac{\partial u}{\partial t}
-
\kappa\frac{\partial^2u}{\partial x^2}
=
f(x,t),
\qquad
-\infty<x<\infty,\quad t>0;
\tag{18.1a}
$$

$$
\left.u\right|_{t=0}=0,
\qquad
-\infty<x<\infty.
\tag{18.1b}
$$

**解**　正如例 12.9 的脚注中指出的，在这种无界区间的定解问题中，物理中习惯上往往并不明确列出边界条件。实际上，“无界区间”只是一个物理上的抽象，它只是表明在所考察的限度（时间，精度等）内，两端的影响可以忽略。因此，读者应当明确，如果要完整地列出定解问题的话，则还应当有边界条件

$$
\left.u\right|_{x\to\pm\infty}\text{ 有界}.
\tag{18.2}
$$

现在，对变量 $t$ 作 Laplace 变换。令

$$
u(x,t)\doteq U(x,p)
=
\int_0^\infty u(x,t)e^{-pt}\,dt,
\tag{18.3}
$$

于是

$$
\frac{\partial u}{\partial t}
\doteq
pU(x,p),
\tag{18.4}
$$

$$
\frac{\partial^2u}{\partial x^2}
\doteq
\frac{d^2U(x,p)}{dx^2}.
\tag{18.5}
$$

这里，在写出 (18.4) 时，已经利用了初始条件 (18.1b)。再进一步令

$$
f(x,t)\doteq F(x,p),
\tag{18.6}
$$

这样，在经过 Laplace 变换后，定解问题 (18.1) 就变成常微分方程的边值问题

$$
pU(x,p)
-
\kappa\frac{d^2U(x,p)}{dx^2}
=
F(x,p),
$$

$$
\left.U(x,p)\right|_{x\to\pm\infty}\text{ 有界}.
$$

利用例 10.8 的结果 (10.76) 式，可以得到

$$
U(x,p)
=
\frac1{2\sqrt{\kappa p}}
\int_{-\infty}^{\infty}
F(x',p)
\exp
\left\{
-\sqrt{\frac p\kappa}|x-x'|
\right\}
dx'.
\tag{18.7}
$$

再根据 Laplace 变换的反演公式（见例 8.8，(8.40) 式）以及卷积定理（见 §8.3），得到

$$
u(x,t)
=
\frac1{2\sqrt{\kappa\pi}}
\int_{-\infty}^{\infty}
dx'
\int_0^t
\exp
\left\{
-\frac{(x-x')^2}{4\kappa(t-\tau)}
\right\}
\frac{f(x',\tau)}{\sqrt{t-\tau}}
\,d\tau.
\tag{18.8}
$$

从以上的求解过程可以看出，用 Laplace 变换求解偏微分方程定解问题，除了可以减少自变量的数目以外，某些已知函数（例如方程的非齐次项，不一定容易求得它的像函数）的像函数甚至都不必具体求出，在求反演时只需应用卷积定理即可。

再举一个无界弦的波动问题的例子。

**例 18.2　用 Laplace 变换求解无界弦的波动问题**

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
-\infty<x<\infty,\quad t>0;
\tag{18.9a}
$$

$$
\left.u\right|_{t=0}
=
\phi(x),
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
\psi(x),
\qquad
-\infty<x<\infty.
\tag{18.9b}
$$

**解**　设在 Laplace 变换之下，

$$
u(x,t)\doteq U(x,p),
$$

于是，原来的定解问题就化为

$$
p^2U(x,p)
-
a^2\frac{d^2U(x,p)}{dx^2}
=
p\phi(x)+\psi(x).
\tag{18.10}
$$

考虑到 $U(x,p)$ 在 $\pm\infty$ 的行为，即可援引例 10.8 的结果，从而写出此方程的解

$$
U(x,p)
=
\frac1{2a}
\int_{-\infty}^{\infty}
\left[
\phi(x')
+
\frac{\psi(x')}{p}
\right]
\exp
\left\{
-\frac pa|x-x'|
\right\}
dx'.
\tag{18.11}
$$

因为

$$
e^{-\alpha p}
\doteq
\delta(t-\alpha),
\qquad
\frac1p e^{-\alpha p}
\doteq
\eta(t-\alpha),
$$

所以

$$
\begin{aligned}
u(x,t)
&=
\frac1{2a}
\int_{-\infty}^{\infty}
\phi(x')
\delta
\left(
t-\frac{|x-x'|}{a}
\right)
dx'\\
&\quad+
\frac1{2a}
\int_{-\infty}^{\infty}
\psi(x')
\eta
\left(
t-\frac{|x-x'|}{a}
\right)
dx'\\
&=
\frac12
\int_{-\infty}^{\infty}
\phi(x')
\delta(at-|x-x'|)
dx'
+
\frac1{2a}
\int_{-\infty}^{\infty}
\psi(x')
\eta(at-|x-x'|)
dx'.
\end{aligned}
$$

注意到

$$
\delta(at-|x-x'|)
=
\begin{cases}
0,& |x-x'|\ne at,\\
\infty,& |x-x'|=at,
\end{cases}
$$

及

$$
\eta(at-|x-x'|)
=
\begin{cases}
0,& |x-x'|>at,\\
1,& |x-x'|<at,
\end{cases}
$$

就可以求出

$$
\begin{aligned}
u(x,t)
&=
\frac12
\int_{-\infty}^{\infty}
\phi(x')
\delta(at-|x-x'|)
dx'
+
\frac1{2a}
\int_{x-at}^{x+at}
\psi(x')\,dx'\\
&=
\frac12
[
\phi(x-at)+\phi(x+at)
]
+
\frac1{2a}
\int_{x-at}^{x+at}
\psi(x')\,dx'.
\end{aligned}
\tag{18.12}
$$

这也正是第十二章中应用行波法求解得到的结果。

Laplace 变换也可用于求解半无界空间的定解问题。我们看半无界空间的例子。

**例 18.3　求解半无界问题**

$$
\begin{cases}
\displaystyle
\frac{\partial u}{\partial t}
-
\kappa\frac{\partial^2u}{\partial x^2}
=
0,
&
x>0,\ t>0,
\\[8pt]
\left.u\right|_{x=0}
=
u_0,
&
\left.u\right|_{t=0}
=
0.
\end{cases}
\tag{18.13}
$$

**解**　当然，本题默认无穷远条件为

$$
\left.u\right|_{x\to\infty}\text{ 有界}.
\tag{18.14}
$$

对变量 $t$ 作 Laplace 变换，

$$
u(x,t)\doteq U(x,p)
=
\int_0^\infty
u(x,t)e^{-pt}\,dt,
\tag{18.15}
$$

$$
\frac{\partial u}{\partial t}
\doteq
pU(x,p)
-
\left.u\right|_{t=0}.
\tag{18.16}
$$

于是，在经过 Laplace 变换后，定解问题 (18.13) 就转化为常微分方程的边值问题

$$
\begin{cases}
\displaystyle
pU(x,p)
-
\kappa\frac{d^2U(x,p)}{dx^2}
=
0,
\\[8pt]
\left.U(x,p)\right|_{x=0}
=
\dfrac{u_0}{p},
\qquad
\left.U(x,p)\right|_{x\to\infty}\text{ 有界}.
\end{cases}
\tag{18.17}
$$

解得

$$
U(x,p)
=
\frac{u_0}{p}
e^{-\sqrt{p/\kappa}\,x}.
\tag{18.18}
$$

利用例 8.8 的结果 (8.42) 式，反演即得

$$
u(x,t)
=
u_0\operatorname{erfc}
\frac{x}{2\sqrt{\kappa t}}.
\tag{18.19}
$$

Laplace 变换也可以应用于求解有界区间内的偏微分方程定解问题。

**例 18.4**　设有长为 $l$ 的均匀细杆，一端保持温度为 $u_0$，另一端绝热。杆的初温为 0。求杆中温度的分布和变化。

**解**　首先列出杆中温度 $u(x,t)$ 所满足的定解问题

$$
\frac{\partial u}{\partial t}
-
\kappa\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
0<x<l,\quad t>0,
\tag{18.20a}
$$

$$
\left.u\right|_{x=0}
=
u_0,
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
0,
\qquad
t>0,
\tag{18.20b}
$$

$$
\left.u\right|_{t=0}
=
0,
\qquad
0<x<l.
\tag{18.20c}
$$

作 Laplace 变换 $u(x,t)\doteq U(x,p)$，于是，定解问题 (18.20) 变成

$$
\kappa\frac{d^2U(x,p)}{dx^2}
-
pU(x,p)
=
0,
\tag{18.21a}
$$

$$
U(0,p)
=
\frac{u_0}{p},
\qquad
\left.
\frac{dU(x,p)}{dx}
\right|_{x=l}
=
0.
\tag{18.21b}
$$

解之得

$$
U(x,p)
=
\frac{u_0}{p}
\frac{
\cosh\sqrt{p/\kappa}\,(l-x)
}{
\cosh\sqrt{p/\kappa}\,l
}.
\tag{18.22}
$$

代入普遍反演公式

$$
u(x,t)
=
\frac1{2\pi i}
\int_L
U(x,p)e^{pt}\,dp,
\tag{18.23}
$$

就可以求出定解问题 (18.20) 的解。考虑到 $U(x,p)$ 的具体形式，积分路径 $L$ 应该取为在右半平面上的一条平行于虚轴的无穷直线，$\operatorname{Re}p>0$（保证像函数在 $L$ 之右解析）。

现在用留数定理来计算这个积分。被积函数在左半平面内有无穷多个奇点，它们是 $p=0$ 和 $\cosh\sqrt{p/\kappa}\,l$ 的零点，即

$$
\sqrt{\frac{p_n}{\kappa}}\,l
=
\frac{2n+1}{2}\pi i,
\qquad
p_n
=
-
\left(
\frac{2n+1}{2l}\pi
\right)^2
\kappa,
\qquad
n=0,1,2,\cdots,
$$

且全是一阶极点。容易求出被积函数在孤立奇点 $p=0$ 处的留数为 $u_0$，而在其余孤立奇点处留数为

$$
\left.
\frac{2u_0}{\sqrt{p/\kappa}\,l}
\frac{
\cosh\sqrt{p/\kappa}\,(l-x)
}{
\sinh\sqrt{p/\kappa}\,l
}
e^{pt}
\right|_{p=p_n}
=
-\frac{4u_0}{\pi}
\frac1{2n+1}
\sin\frac{2n+1}{2l}\pi x
\exp
\left[
-
\left(
\frac{2n+1}{2l}\pi
\right)^2
\kappa t
\right].
$$

所以，最后就求得

$$
u(x,t)
=
u_0
-
\frac{4u_0}{\pi}
\sum_{n=0}^{\infty}
\frac1{2n+1}
\sin\frac{2n+1}{2l}\pi x
\exp
\left[
-
\left(
\frac{2n+1}{2l}\pi
\right)^2
\kappa t
\right].
\tag{18.24}
$$

不难验证，这个解式和用分离变量法得到的结果完全相同。

从以上的求解过程可以看出，用 Laplace 变换求解偏微分方程定解问题还有一个优点，这就是不必将非齐次的边界条件齐次化，因为这时原有的偏微分方程定解问题的非齐次边界条件将转化为常微分方程的非齐次边界条件，这并不会带来原则性的困难。

对于这个例子，也还可以用别的方法求反演。例如，可将像函数展开为

$$
\begin{aligned}
\frac{
\cosh\sqrt{p/\kappa}\,(l-x)
}{
\cosh\sqrt{p/\kappa}\,l
}
&=
e^{-\sqrt{p/\kappa}\,l}
\left(
1+e^{-2\sqrt{p/\kappa}\,l}
\right)^{-1}
\left[
e^{\sqrt{p/\kappa}\,(l-x)}
+
e^{-\sqrt{p/\kappa}\,(l-x)}
\right]\\
&=
\left[
e^{-\sqrt{p/\kappa}\,x}
+
e^{-\sqrt{p/\kappa}\,(2l-x)}
\right]
\sum_{n=0}^{\infty}
(-)^n
e^{-2n\sqrt{p/\kappa}\,l}\\
&=
\sum_{n=0}^{\infty}
(-)^n
e^{-\sqrt{p/\kappa}\,(2nl+x)}
-
\sum_{n=1}^{\infty}
(-)^n
e^{-\sqrt{p/\kappa}\,(2nl-x)}.
\end{aligned}
\tag{18.22'}
$$

利用 Laplace 变换的反演公式（见例 8.8，(8.42) 式），就可以得到定解问题 (18.20) 的另一种形式的解

$$
u(x,t)
=
u_0
\sum_{n=0}^{\infty}
(-)^n
\operatorname{erfc}
\frac{2nl+x}{2\sqrt{\kappa t}}
-
u_0
\sum_{n=1}^{\infty}
(-)^n
\operatorname{erfc}
\frac{2nl-x}{2\sqrt{\kappa t}}.
\tag{18.25}
$$

定解问题解的存在唯一性，就保证了这两个不同形式的解实际上是相等的，由余误差函数的有关公式也可以将解式 (18.24) 化为 (18.25)。

再讨论一下半无界杆的极限情形。这时在 (18.22) 或 (18.22') 式中令 $l\to\infty$，就有

$$
U(x,p)
=
\frac{u_0}{p}
\exp
\left\{
-\sqrt{\frac p\kappa}\,x
\right\}.
$$

所以，对于半无界杆，就得到

$$
u(x,t)
=
u_0\operatorname{erfc}
\frac{x}{2\sqrt{\kappa t}},
$$

即 (18.19) 式。

应用 Laplace 变换求解偏微分方程定解问题也有一个明显的局限性，这就是通常只能对时间变量 $t$ 作 Laplace 变换。在一般情形下，我们不建议对空间变量作 Laplace 变换，即使是半无界空间的情形。这是因为对于空间变量，定解问题中只出现边界条件，而 Laplace 变换需要用到的是初始条件。

## §18.2 Fourier 变换

Fourier 变换可对空间变量进行。根据空间变量的变化区间，可以选用 Fourier 变换（包括正弦变换和余弦变换）或有限正弦、余弦变换。例如，对于无界区间 $(-\infty,\infty)$ 上的函数 $f(x)$，如果在任意有限区间上只有有限个极大、极小和有限个第一类间断点，且积分 $\int_{-\infty}^{\infty}f(x)\,dx$ 绝对收敛，则它的 Fourier 变换存在，

$$
\mathscr F[f(x)]
\equiv
F(k)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f(x)e^{-ikx}\,dx,
\tag{18.26}
$$

而逆变换（反演）是

$$
\mathscr F^{-1}[F(k)]
\equiv
f(x)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
F(k)e^{ikx}\,dk.
\tag{18.27}
$$

这里的 Fourier 变换和逆变换的形式可能和读者熟悉的形式略有不同。(18.26) 和 (18.27) 的形式更加对称，更多地为物理学所采用。

为了将 Fourier 变换应用于求解偏微分方程定解问题，必然涉及函数的一、二阶导数的 Fourier 变换。设 $f(x)$ 的 Fourier 变换存在，于是

$$
\begin{aligned}
\mathscr F[f'(x)]
&=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f'(x)e^{-ikx}\,dx\\
&=
\left.
\frac1{\sqrt{2\pi}}
f(x)e^{-ikx}
\right|_{-\infty}^{\infty}
+
\frac{ik}{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f(x)e^{-ikx}\,dx.
\end{aligned}
$$

由于积分 $\int_{-\infty}^{\infty}f(x)\,dx$ 绝对收敛，就一定有 $\lim_{x\to\pm\infty}f(x)=0$，所以

$$
\mathscr F[f'(x)]
=
\frac{ik}{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f(x)e^{-ikx}\,dx
=
ik\mathscr F[f(x)].
\tag{18.28}
$$

更进一步，当然就有

$$
\mathscr F[f''(x)]
=
-k^2\mathscr F[f(x)].
\tag{18.29}
$$

下面用 Fourier 变换重新求解例 18.1 和例 18.2。对于例 18.1，即定解问题 (18.1)，假设 $u(x,t)$ 的 Fourier 变换存在，记

$$
U(k,t)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
u(x,t)e^{-ikx}\,dx,
\tag{18.30}
$$

并设

$$
F(k,t)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f(x,t)e^{-ikx}\,dx.
\tag{18.31}
$$

这样，在对变量 $x$ 作 Fourier 变换后，定解问题 (18.1) 就变为

$$
\frac{dU(k,t)}{dt}
+
\kappa k^2U(k,t)
=
F(k,t),
\tag{18.32a}
$$

$$
\left.U(k,t)\right|_{t=0}
=
0.
\tag{18.32b}
$$

用常数变易法求解这个一阶常微分方程的定解问题，就得到

$$
U(k,t)
=
e^{-\kappa k^2t}
\int_0^t
F(k,\tau)e^{\kappa k^2\tau}\,d\tau.
\tag{18.33}
$$

再求反演，有

$$
\begin{aligned}
u(x,t)
&=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
U(k,t)e^{ikx}\,dk\\
&=
\int_0^t
\left[
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
F(k,\tau)e^{-\kappa k^2(t-\tau)}e^{ikx}\,dk
\right]d\tau.
\end{aligned}
$$

利用 (4.29) 式的结果，可以算出

$$
\begin{aligned}
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
e^{-\kappa k^2(t-\tau)}e^{ikx}\,dk
&=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
e^{-\kappa k^2(t-\tau)}\cos kx\,dk\\
&=
\frac1{\sqrt{2\kappa(t-\tau)}}
\exp
\left[
-\frac{x^2}{4\kappa(t-\tau)}
\right].
\end{aligned}
\tag{18.34}
$$

再利用 (18.31) 的逆变换

$$
f(x,t)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
F(k,t)e^{ikx}\,dk
$$

及 Fourier 变换的卷积公式（请读者自己证明）

$$
\mathscr F[f_1(x)]
\mathscr F[f_2(x)]
=
\mathscr F
\left[
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f_1(\xi)f_2(x-\xi)\,d\xi
\right],
\tag{18.35}
$$

最后就能得到

$$
\begin{aligned}
u(x,t)
&=
\int_0^t
\left\{
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\frac{f(\xi,\tau)}
{\sqrt{2\kappa(t-\tau)}}
\exp
\left[
-\frac{(x-\xi)^2}{4\kappa(t-\tau)}
\right]
d\xi
\right\}d\tau\\
&=
\frac1{2\sqrt{\kappa\pi}}
\int_0^t
\left\{
\int_{-\infty}^{\infty}
f(\xi,\tau)
\exp
\left[
-\frac{(x-\xi)^2}{4\kappa(t-\tau)}
\right]
d\xi
\right\}
\frac{d\tau}{\sqrt{t-\tau}}.
\end{aligned}
\tag{18.36}
$$

和上一节中得到的解式 (18.8) 的形式完全一样。

再来解无界弦上的自由振动问题，定解问题见 (18.9) 式。仍设 $u(x,t)$ 的 Fourier 变换

$$
U(k,t)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
u(x,t)e^{-ikx}\,dx
\tag{18.37}
$$

存在，并设

$$
\Phi(k)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\phi(x)e^{-ikx}\,dx,
\tag{18.38}
$$

$$
\Psi(k)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\psi(x)e^{-ikx}\,dx,
\tag{18.39}
$$

于是，作 Fourier 变换后，定解问题 (18.9) 就变为

$$
\frac{d^2U(k,t)}{dt^2}
+
k^2a^2U(k,t)
=
0,
\tag{18.40a}
$$

$$
\left.U(k,t)\right|_{t=0}
=
\Phi(k),
\qquad
\left.
\frac{dU(k,t)}{dt}
\right|_{t=0}
=
\Psi(k).
\tag{18.40b}
$$

这是一个二阶常微分方程的初值问题，解之得

$$
U(k,t)
=
\Phi(k)\cos kat
+
\Psi(k)\frac{\sin kat}{ka}.
\tag{18.41}
$$

根据 Fourier 变换的反演公式，就可以求出

$$
u(x,t)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\left[
\Phi(k)\cos kat
+
\Psi(k)\frac{\sin kat}{ka}
\right]
e^{ikx}\,dk.
$$

注意

$$
\begin{aligned}
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\Phi(k)\cos kat\,e^{ikx}\,dk
&=
\frac1{\sqrt{2\pi}}
\frac12
\int_{-\infty}^{\infty}
\Phi(k)
\left[
e^{ik(x+at)}
+
e^{ik(x-at)}
\right]dk\\
&=
\frac12
[
\phi(x+at)+\phi(x-at)
],
\end{aligned}
$$

$$
\begin{aligned}
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\Psi(k)\frac{\sin kat}{ka}e^{ikx}\,dk
&=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\Psi(k)
\left[
\int_0^t\cos ka\tau\,d\tau
\right]
e^{ikx}\,dk\\
&=
\int_0^t
\left[
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\Psi(k)\cos ka\tau\,e^{ikx}\,dk
\right]d\tau\\
&=
\frac12
\int_0^t
[
\psi(x+a\tau)+\psi(x-a\tau)
]d\tau\\
&=
\frac1{2a}
\int_{x-at}^{x+at}
\psi(\xi)\,d\xi.
\end{aligned}
$$

最后得到的解式和 (18.12) 式完全相同：

$$
u(x,t)
=
\frac12
[
\phi(x+at)+\phi(x-at)
]
+
\frac1{2a}
\int_{x-at}^{x+at}
\psi(\xi)\,d\xi.
\tag{18.42}
$$

求解偏微分方程定解问题时，还可以在对时间变量作 Laplace 变换的同时，对空间变量作 Fourier 变换，即 Laplace 变换和 Fourier 变换可以混合在一起同时使用。

再来求解无界杆的热传导问题，定解问题见 (18.1) 式。对变量 $t$ 作 Laplace 变换，对变量 $x$ 作 Fourier 变换，令

$$
\mathscr F\{u(x,t)\}
=
U(k,t),
\tag{18.43}
$$

$$
\mathscr F\{f(x,t)\}
=
F(k,t),
\tag{18.44}
$$

$$
U(k,t)\doteq\widetilde U(k,p),
\tag{18.45}
$$

$$
F(k,t)\doteq\widetilde F(k,p).
\tag{18.46}
$$

对变量 $x$ 作 Fourier 变换后定解问题 (18.1) 变为

$$
\begin{cases}
\displaystyle
\frac{dU(k,t)}{dt}
+
\kappa k^2U(k,t)
=
F(k,t),
\\[8pt]
\left.U(k,t)\right|_{t=0}
=
0.
\end{cases}
$$

再对变量 $t$ 作 Laplace 变换后定解问题 (18.1) 变成一个代数方程

$$
p\widetilde U(k,p)
+
\kappa k^2\widetilde U(k,p)
=
\widetilde F(k,p).
$$

解得

$$
\widetilde U(k,p)
=
\frac1{p+\kappa k^2}
\widetilde F(k,p).
$$

应用卷积公式，求反演，得

$$
U(k,t)
=
\int_0^t
F(k,\tau)e^{-\kappa k^2(t-\tau)}\,d\tau.
$$

再作 Fourier 逆变换，利用 (18.34) 式的结果，就得到

$$
\begin{aligned}
u(x,t)
&=
\int_0^t
\left\{
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\frac{f(\xi,\tau)}
{\sqrt{2\kappa(t-\tau)}}
\exp
\left[
-\frac{(x-\xi)^2}{4\kappa(t-\tau)}
\right]d\xi
\right\}d\tau\\
&=
\frac1{2\sqrt{\kappa\pi}}
\int_0^t
\left\{
\int_{-\infty}^{\infty}
f(\xi,\tau)
\exp
\left[
-\frac{(x-\xi)^2}{4\kappa(t-\tau)}
\right]d\xi
\right\}
\frac{d\tau}{\sqrt{t-\tau}}.
\end{aligned}
$$

和上一节中得到的解式 (18.8) 的形式完全一样。


## *§18.3 半无界空间的情形

如果 $f(x)$ 是定义在半无界区间 $[0,\infty)$ 上，则可根据 $x=0$ 端边界条件的不同类型，选用正弦变换

$$
F(k)
=
\sqrt{\frac2\pi}
\int_0^\infty
f(x)\sin kx\,dx,
\tag{18.47}
$$

$$
f(x)
=
\sqrt{\frac2\pi}
\int_0^\infty
F(k)\sin kx\,dk,
\tag{18.48}
$$

或余弦变换

$$
F(k)
=
\sqrt{\frac2\pi}
\int_0^\infty
f(x)\cos kx\,dx,
\tag{18.49}
$$

$$
f(x)
=
\sqrt{\frac2\pi}
\int_0^\infty
F(k)\cos kx\,dk.
\tag{18.50}
$$

在正弦变换下，因为

$$
\sqrt{\frac2\pi}
\int_0^\infty
f'(x)\sin kx\,dx
=
-\sqrt{\frac2\pi}\,
k
\int_0^\infty
f(x)\cos kx\,dx,
\tag{18.51}
$$

$$
\sqrt{\frac2\pi}
\int_0^\infty
f''(x)\sin kx\,dx
=
\sqrt{\frac2\pi}\,
k f(0)
-
k^2F(k).
\tag{18.52}
$$

由此可见，对于二阶偏微分方程的定解问题，则只有当定解问题中仅出现未知函数及其二阶偏导数，且在半无界空间的 $x=0$ 端给出的是第一类边界条件时，才可以选用正弦变换。

同样，对于余弦变换，也有

$$
\sqrt{\frac2\pi}
\int_0^\infty
g'(x)\cos kx\,dx
=
-\sqrt{\frac2\pi}\,
g(0)
+
\sqrt{\frac2\pi}\,
k
\int_0^\infty
g(x)\sin kx\,dx,
\tag{18.53}
$$

$$
\sqrt{\frac2\pi}
\int_0^\infty
g''(x)\cos kx\,dx
=
-\sqrt{\frac2\pi}\,
g'(0)
-
k^2G(k).
\tag{18.54}
$$

所以，如果还限于二阶偏微分方程的定解问题，则只有当定解问题中仅出现未知函数及其二阶偏导数，而且在 $x=0$ 端给出的是第二类边界条件时，才可以选用余弦变换。

**例 18.5　求解半无界空间的稳定问题**

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
0,
\qquad
-\infty<x<\infty,\quad y>0;
\tag{18.55a}
$$

$$
\left.u\right|_{y=0}
=
f(x),
\qquad
-\infty<x<\infty.
\tag{18.55b}
$$

**解**　对变量 $y$ 作正弦变换，有

$$
U(x,k)
\equiv
\mathscr F[u(x,y)]
=
\sqrt{\frac2\pi}
\int_0^\infty
u(x,y)\sin ky\,dy,
\tag{18.56}
$$

$$
\mathscr F
\left[
\frac{\partial^2u(x,y)}{\partial y^2}
\right]
=
\sqrt{\frac2\pi}\,
k f(x)
-
k^2U(x,y).
\tag{18.57}
$$

于是，在经过正弦变换后，定解问题 (18.55) 就转化为求无界区间 $-\infty<x<\infty$ 上

$$
\frac{d^2U(x,k)}{dx^2}
-
k^2U(x,k)
=
\sqrt{\frac2\pi}\,
k f(x)
\tag{18.58}
$$

的有界解。再利用例 10.8 的结果，有

$$
U(x,k)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
e^{-k|x-x'|}
f(x')\,dx'.
\tag{18.59}
$$

最后，根据正弦变换的逆变换公式，就可以求得

$$
\begin{aligned}
u(x,y)
&=
\frac1\pi
\int_0^\infty
\left[
\int_{-\infty}^{\infty}
e^{-k|x-x'|}
f(x')\,dx'
\right]
\sin ky\,dk\\
&=
\frac1\pi
\int_{-\infty}^{\infty}
f(x')
\left[
\int_0^\infty
e^{-k|x-x'|}
\sin ky\,dk
\right]dx'\\
&=
\frac y\pi
\int_{-\infty}^{\infty}
\frac{f(x')}
{(x-x')^2+y^2}
\,dx'.
\end{aligned}
\tag{18.60}
$$

这正是上半平面的 Poisson 公式（见 (3.32b) 式）。

## *§18.4 关于积分变换的一般讨论

以上介绍了两类积分变换（Laplace 变换和 Fourier 变换，包括 Fourier 变换的另外两种特殊形式，即正弦变换和余弦变换）在求解偏微分方程定解问题中的应用。可以把这些变换概括写成

$$
F(k)
=
\int_a^b
K(k,x)f(x)\,dx,
\tag{18.61}
$$

它把自变量 $x\in[a,b]$（这里的 $x$ 也可以代表时间，区间也可以是无界或半无界的）的函数 $f(x)$ 变换为复变量 $k$ 的函数 $F(k)$，其中 $K(k,x)$ 称为积分变换的核：

$$
\text{Laplace 变换}
\qquad
K(k,x)
=
e^{-kx},
\qquad
0\le x<\infty;
\tag{18.62a}
$$

$$
\text{Fourier 变换}
\qquad
K(k,x)
=
e^{-ikx},
\qquad
-\infty<x<\infty;
\tag{18.62b}
$$

$$
\text{正弦变换}
\qquad
K(k,x)
=
\sin kx,
\qquad
0\le x<\infty;
\tag{18.62c}
$$

$$
\text{余弦变换}
\qquad
K(k,x)
=
\cos kx,
\qquad
0\le x<\infty.
\tag{18.62d}
$$

其他还有

$$
\text{Hankel 变换}
\qquad
K(k,x)
=
xJ_n(kx),
\qquad
0\le x<\infty;
\tag{18.62e}
$$

$$
\text{Mellin 变换}
\qquad
K(k,x)
=
x^{k-1},
\qquad
0\le x<\infty.
\tag{18.62f}
$$

就一个具体的偏微分方程（为了叙述的方便，不妨仍限于二阶偏微分方程）定解问题而言，到底应当选用哪一种积分变换，要考虑以下几个原则：

原则 1　所涉及的自变量的变化区间和该变换的要求一致。

原则 2　未知函数的该种积分变换存在。

原则 3　要求函数 $f(x)$ 及其导数 $f'(x),f''(x)$ 在该变换下有简单的代数关系。

原则 4　涉及的未知函数及其导数 $f'(x)$ 的特殊值正好由定解问题中的定解条件给出。

应该说，上面的原则 3 还只是适用于常系数的微分方程。如果讨论变系数的偏微分方程定解问题，那么这一条还需要修改。例如，对于偏微分方程

$$
[
\hat L_1(x)
+
\hat L_2(y)
]
u(x,y)
=
f(x,y),
$$

其中 $\hat L_1(x)$ 和 $\hat L_2(y)$ 分别是 $x$ 和 $y$ 的微分算符，假定算符 $\hat L_1(x)$ 中的系数都是 $x$ 的实值函数，再设 $u(x,y)$ 也是实值函数，则在积分变换

$$
\int_a^b
K(k,x)u(x,y)\,dx
=
U(k,y)
\tag{18.63}
$$

之下，方程变为

$$
\int_a^b
K(k,x)f(x,y)\,dx
=
\int_a^b
K(k,x)
[
\hat L_1(x)u(x,y)
]
dx
+
\hat L_2(y)
\int_a^b
K(k,x)u(x,y)\,dx,
$$

即

$$
\int_a^b
[
\hat M_1(x)K(k,x)
]
u(x,y)\,dx
+
\hat L_2(y)U(k,y)
=
F(k,y),
\tag{18.64}
$$

其中 $\hat M_1(x)$ 是算符 $\hat L_1(x)$ 的伴算符。为了保证方程 (18.64) 是关于 $U(k,y)$ 的微分方程，必须有

$$
\hat M_1(x)K(k,x)
=
\lambda K(k,x),
\tag{18.65}
$$

即 $\hat M_1(x)K(k,x)$ 必须与 $K(k,x)$ 成正比。这就限定了所能选择的变换核 $K(k,x)$。例如，在柱坐标系中求解 Laplace 方程或 Poisson 方程时，对于变量 $r$，就只能选用 Hankel 变换或 Mellin 变换。

下面就举一个应用 Hankel 变换的例子。

**例 18.6　应用 Hankel 变换求解带电导体圆盘的静电势。**

**解**　采用柱坐标系。定解问题是

$$
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial u}{\partial r}
\right)
+
\frac{\partial^2u}{\partial z^2}
=
0,
\qquad
0<r<\infty,\quad z>0;
\tag{18.66a}
$$

$$
\left.u\right|_{r=0}\text{ 有界},
\qquad
\left.u\right|_{r\to\infty}
\to0;
\tag{18.66b}
$$

$$
\left.u\right|_{z=0}
=
u_0,
\qquad
r<a;
\tag{18.66c}
$$

$$
\left.
\frac{\partial u}{\partial z}
\right|_{z=0}
=
0,
\qquad
r>a;
\tag{18.66d}
$$

$$
\left.u\right|_{z\to\infty}
\to0.
\tag{18.66e}
$$

作 Hankel 变换，即令

$$
U(p,z)
=
\int_0^\infty
u(r,z)J_0(pr)r\,dr.
\tag{18.67}
$$

容易证明，在边界条件 (18.66b) 之下，

$$
\int_0^\infty
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial u}{\partial r}
\right)
J_0(pr)r\,dr
=
-p^2U(p,z).
\tag{18.68}
$$

所以，方程 (18.66a) 和边界条件 (18.66b) 就变换为

$$
\frac{d^2U(p,z)}{dz^2}
-
p^2U(p,z)
=
0.
\tag{18.69a}
$$

同样，边界条件 (18.66e) 就变换为

$$
\left.U(p,z)\right|_{z\to\infty}
\to0.
\tag{18.69b}
$$

解之即得

$$
U(p,z)
=
A(p)e^{-pz}.
$$

现在的问题是，难以将平面 $z=0$ 上的边界条件 (18.66c) 和 (18.66d) 也代入变换 (18.67)，因为这一组边界条件给出的是 $0\le r<a$ 时的 $\left.u(r,z)\right|_{z=0}$ 值和 $r>a$ 时的 $\left.\partial u/\partial z\right|_{z=0}$ 值。这样，就只得先求反演，而得到了定解问题 (18.66) 的积分形式的解

$$
u(r,z)
=
\int_0^\infty
A(p)e^{-pz}J_0(pr)p\,dp.
\tag{18.70}
$$

然后再设法定出函数 $A(p)$。为此，将 (18.70) 式代入边界条件 (18.66c) 和 (18.66d)，可以得到一对方程

$$
\int_0^\infty
A(p)J_0(pr)p\,dp
=
u_0,
\qquad
0<r<a;
$$

$$
\int_0^\infty
A(p)J_0(pr)p^2\,dp
=
0,
\qquad
r>a.
$$

关于这种从一对（含 Bessel 函数的）积分方程中求函数 $A(p)$ 的问题，参考书目 [13] 中给出了某些特殊情形下的解（见该书第二册中译本，87—89 页），例如，方程组

$$
\int_0^\infty
f(t)J_\nu(xt)\,dt
=
x^{-\nu}M(x),
\qquad
0<x<1,
$$

$$
\int_0^\infty
f(t)J_\nu(xt)t\,dt
=
0,
\qquad
x>1
$$

的解就是

$$
f(t)
=
\sqrt{\frac{2t}{\pi}}
\int_0^1
x^{\nu+1/2}\xi(x)
J_{\nu-1/2}(xt)\,dx,
$$

$$
\xi(x)
=
x^{-2\nu}M(0)
+
x^{1-2\nu}
\int_0^x
(x^2-y^2)^{-1/2}
M'(y)\,dy.
$$

特别是，当 $\nu=0$ 时，有

$$
f(t)
=
\frac2\pi
\int_0^1
\xi(x)\cos xt\,dx,
$$

$$
\xi(x)
=
M(0)
+
x
\int_0^x
(x^2-y^2)^{-1/2}
M'(y)\,dy.
$$

这样，就能定出

$$
A(p)
=
\frac{2u_0}{\pi}
\frac{\sin ap}{p^2}.
$$

代入 (18.70)，并算出积分，就可以得到带电导体圆盘的静电势

$$
\begin{aligned}
u(r,z)
&=
\frac{2u_0}{\pi}
\int_0^\infty
e^{-pz}J_0(pr)
\frac{\sin ap}{p}
\,dp\\
&=
\frac{2u_0}{\pi}
\arcsin
\frac{
2a
}{
\sqrt{z^2+(r+a)^2}
+
\sqrt{z^2+(r-a)^2}
}.
\end{aligned}
\tag{18.71}
$$

以上讨论的都是无界或半无界区间上的积分变换，它们的共同特点是复变量 $k$ 的取值也是连续的。在实用中还有有界区间上的积分变换，由于这时的 $k$ 只能取离散值，所以可以把变换核记为 $K_n(x)$。例如，常见的有界区间上的积分变换有

$$
\text{有限正弦变换}
\qquad
K_n(x)
=
\sin nx,
\qquad
n=1,2,3,\cdots,
\qquad
0\le x\le\pi;
$$

$$
\text{有限余弦变换}
\qquad
K_n(x)
=
\cos nx,
\qquad
n=0,1,2,\cdots,
\qquad
0\le x\le\pi;
$$

$$
\text{Legendre 变换}
\qquad
K_n(x)
=
P_n(x),
\qquad
n=0,1,2,\cdots,
\qquad
-1\le x\le1.
$$

可以看出，这些积分变换的核恰好就是定义在各自区间上的本征函数。

## *§18.5 小波变换简介

“小波分析”是近年来发展起来的一个比较新的理论课题和应用数学方法。小波，作为表示函数的一种新的基，作为时间-频率分析的一种技术，也已经成为一个新的数学学科，而且还处在迅速发展的过程之中。本节只能从积分变换的角度，对小波变换做一个入门性的介绍。

先分析一下传统的 Fourier 变换。大家知道，在 $(-\infty,\infty)$ 上定义了函数 $f(t)$，只要 $f(t)$ 满足一定条件，则它的 Fourier 变换

$$
F(\omega)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f(t)e^{-i\omega t}\,dt
\tag{18.72a}
$$

以及逆变换

$$
f(t)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
F(\omega)e^{i\omega t}\,d\omega
\tag{18.72b}
$$

都存在。通常，变量 $t$ 代表时间，$\omega$ 代表频率。$F(\omega)$ 就给出了信号 $f(t)$ 的频谱。从 Fourier 变换的公式就可以看出，为了研究一个信号的频谱特性，必须提供信号在整个时间范围 $(-\infty,\infty)$ 内的全部变化情况，甚至包括将来的变化。另一方面如果信号在某一时段内发生了变化，那么，整个频谱都会受到影响。作为一个极端的情形，只出现在 $t_0$ 时刻的脉冲信号 $\delta(t-t_0)$，其频谱是 $e^{-i\omega t_0}/\sqrt{2\pi}$，就覆盖了全部频率范围。


为了弥补 Fourier 变换的不足，早在半个世纪前，就有人提出过“加窗” Fourier 变换，通过引进时间局部化的“窗函数” $g_\alpha(t-b)$，使得我们可以选择任意一个时段（以 $t=b$ 时刻为中心的一定宽度范围内）的信号，而且通过观测信号在某一频率附近的频谱就可以获得此信号的足够精确的信息。而后平移窗函数，即改变 $b$，就可以覆盖整个时域。Gauss 型的函数就是一个这样的窗函数。这种特殊形式的加窗 Fourier 变换（窗函数为 Gauss 型函数）称为 Gabor 变换。由于 Gauss 型函数

$$
g_\alpha(t)
=
\frac1{2\sqrt{\pi\alpha}}
\exp
\left\{
-\frac{t^2}{4\alpha}
\right\}
$$

的 Fourier 变换

$$
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\frac1{2\sqrt{\pi\alpha}}
\exp
\left\{
-\frac{t^2}{4\alpha}
\right\}
e^{-i\omega t}\,dt
=
\frac1{\sqrt{2\pi}}
e^{-\alpha\omega^2}
\tag{18.73}
$$

仍然是 Gauss 型的函数，所以，根据 Fourier 变换的卷积公式

$$
\int_{-\infty}^{\infty}
f_1(t)f_2(t)e^{-i\omega t}\,dt
=
\int_{-\infty}^{\infty}
F_1(\xi)F_2(\omega-\xi)\,d\xi,
\tag{18.74a}
$$

$$
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f_1(t)e^{-i\omega t}\,dt
=
F_1(\omega),
\tag{18.74b}
$$

$$
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f_2(t)e^{-i\omega t}\,dt
=
F_2(\omega),
\tag{18.74c}
$$

对于任意函数 $f(t)$，就有

$$
\int_{-\infty}^{\infty}
f(t)g_\alpha(t-b)e^{-i\omega t}\,dt
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
F(\xi)
e^{-i(\omega-\xi)b}
e^{-\alpha(\omega-\xi)^2}
\,d\xi,
\tag{18.75}
$$

其中，$F(\omega)$ 是 $f(t)$ 的 Fourier 变换，

$$
F(\omega)
=
\frac1{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
f(t)e^{-i\omega t}\,dt.
\tag{18.76}
$$

对于这个结果，可以把 (18.75) 式的左端理解为函数 $g_\alpha(t-b)e^{i\omega t}$ 和 $f(t)$ 的内积

$$
\left(
g_\alpha(t-b)e^{i\omega t},
f(t)
\right)
=
\int_{-\infty}^{\infty}
\left[
g_\alpha(t-b)e^{i\omega t}
\right]^*
f(t)\,dt,
$$

而右端则是它们的变换的内积，

$$
\left(
\frac1{\sqrt{2\pi}}
e^{i(\omega-\xi)b}
e^{-\alpha(\omega-\xi)^2},
F(\xi)
\right)
=
\int_{-\infty}^{\infty}
\left[
\frac1{\sqrt{2\pi}}
e^{i(\omega-\xi)b}
e^{-\alpha(\omega-\xi)^2}
\right]^*
F(\xi)\,d\xi.
$$

这个结果不过是更普遍的 Parseval 方程（见 (17.19) 式）

$$
(f_1(t),f_2(t))
=
(F_1(\xi),F_2(\xi))
\tag{18.77}
$$

的一个特例。从物理上看，由于函数 $g_\alpha(t-b)$ 在 $t=b$ 处有一个尖锐的峰，

$$
\lim_{\alpha\to0}
g_\alpha(t-b)
=
\delta(t-b),
$$

所以，在 (18.75) 式中，对于左端积分的贡献主要来自 $t=b$ 附近，而对于右端积分的贡献则主要来自 $\xi=\omega$ 附近。换句话说，信号 $f(t)$ 在 $t=b$ 时刻的信息可以通过在频率 $\xi=\omega$ 附近观测这个信号的频谱而得到。可以把 (18.75) 式的右端改写为

$$
\left(
\frac1{\sqrt{2\alpha}}
e^{-i\omega b}
\right)
\int_{-\infty}^{\infty}
F(\xi)
g_{1/(4\alpha)}(\xi-\omega)
e^{i\xi b}\,d\xi.
$$

这说明，若 $g_\alpha(t)$ 为信号的时间窗函数，则 $g_{1/(4\alpha)}(\xi)$ 为相应频谱的频率窗函数。不妨定义

$$
\bar t_{g_\alpha}
=
(g_\alpha(t),tg_\alpha(t))
=
\int_{-\infty}^{\infty}
t|g_\alpha(t)|^2\,dt
\tag{18.78}
$$

为 $g_\alpha(t)$ 的中心。因 $g_\alpha(t)$ 为偶函数，故 $\bar t_{g_\alpha}=0$。在此基础上，进一步定义 $g_\alpha(t)$ 的半宽度

$$
\Delta_{g_\alpha}
\equiv
\sqrt{
\frac{
\left(
(t-\bar t_{g_\alpha})g_\alpha(t),
(t-\bar t_{g_\alpha})g_\alpha(t)
\right)
}{
(g_\alpha(t),g_\alpha(t))
}
}.
\tag{18.79}
$$

显然有

$$
\Delta_{g_\alpha}
=
\sqrt{\alpha}.
\tag{18.80}
$$

同样，可以求出 $g_{1/(4\alpha)}$ 中心

$$
\bar\omega_{g_{1/(4\alpha)}}
=
0
\tag{18.81}
$$

和半宽度

$$
\Delta_{g_{1/(4\alpha)}}
=
\frac1{2\sqrt{\alpha}}.
\tag{18.82}
$$

所以，时间窗宽度 $2\Delta_{g_\alpha}$ 与频率窗宽度 $2\Delta_{g_{1/(4\alpha)}}$ 的乘积为常数，

$$
(2\Delta_{g_\alpha})
\times
(2\Delta_{g_{1/(4\alpha)}})
=
2.
\tag{18.83}
$$

我们也可以在 $t-\xi$ 平面上以 $(b,\omega)$ 点为中心作矩形

$$
b-\Delta_{g_\alpha}
\le
t
\le
b+\Delta_{g_\alpha},
\qquad
\omega-\Delta_{g_{1/(4\alpha)}}
\le
\xi
\le
\omega+\Delta_{g_{1/(4\alpha)}}
$$

来形象化地表示时间-频率的局部化（见图 18.1）。这个矩形

$$
[
b-\Delta_{g_\alpha},
b+\Delta_{g_\alpha}
]
\times
[
\omega-\Delta_{g_{1/(4\alpha)}},
\omega+\Delta_{g_{1/(4\alpha)}}
]
$$

就称为 Gabor 变换的时间-频率窗。$2\Delta_{g_\alpha}$ 和 $2\Delta_{g_{1/(4\alpha)}}$ 又分别称为时间-频率窗的宽度和高度。时间-频率窗的宽度、高度以及面积都是固定的。

![Figure 18.1](...)

图 18.1　Gabor 变换的时间-频率窗

还可以取别的形式的窗函数 $w(t)$。如果要求 $w(t)$ 及 $tw(t)$ 均平方可积（因此宽度 $\Delta_w$ 为有限值），同样，它的 Fourier 变换 $W(\xi)$ 及 $\xi W(\xi)$ 也平方可积（因而宽度 $\Delta_W$ 也为有限值），这样得到的加窗 Fourier 变换就称为短时 Fourier 变换。前面提到的 Gabor 变换是短时 Fourier 变换的一种。

Gabor 变换或其他短时 Fourier 变换的缺点是时间-频率窗是固定的，不能随频率的高低而适当地调整。因为频率与单位时间内的周期数成正比，所以，理想的情况是：要精确研究高频现象，就应当取窄的时间窗；而要研究低频现象，则不妨取较宽的时间窗。因此，Gabor 变换或其他短时 Fourier 变换不适合于处理频域宽、变化激烈的信号。而积分小波变换则是针对这类问题发展起来的。

要介绍积分小波变换，首先要建立基小波及小波的概念。如果 $h(t)$ 及其 Fourier 变换 $H(\omega)$ 均为平方可积函数，宽度均为有限值，且满足相容性条件

$$
C_h
=
\int_{-\infty}^{\infty}
\frac{|H(\omega)|^2}{|\omega|}
\,d\omega
<
\infty,
\tag{18.84}
$$

则称 $h(t)$ 为基小波。条件 (18.84) 意味着 $H(0)=0$，即

$$
\int_{-\infty}^{\infty}
h(t)\,dt
=
0.
\tag{18.85}
$$

这正是称其为小波的原因。在建立积分小波变换的反演时，需要用到相容性条件。

给出了基小波 $h(t)$ 后，通过平移和伸缩，可以得到一族函数

$$
h_{b,a}(t)
=
\frac1{\sqrt{|a|}}
h
\left(
\frac{t-b}{a}
\right),
\qquad
a\ne0,
\tag{18.86}
$$

称为小波。而信号 $f(t)$ 的积分小波变换则定义为

$$
(\mathcal W_h f)(b,a)
=
\frac1{\sqrt{|a|}}
\int_{-\infty}^{\infty}
f(t)
h^*
\left(
\frac{t-b}{a}
\right)
dt
=
(h_{b,a},f).
\tag{18.87}
$$

设函数 $h(t)$ 的中心和半宽度分别为 $\bar t$ 和 $\Delta_h$，则函数 $h_{b,a}(t)$ 是中心在 $b+a\bar t$ 且半宽度为 $|a|\Delta_h$ 的窗函数。因此，(18.87) 式表示的积分小波变换给出了信号 $f(t)$ 在时间窗

$$
[
b+a\bar t-|a|\Delta_h,
\,
b+a\bar t+|a|\Delta_h
]
$$

内的局部信息。$|a|$ 变小时时间窗变窄；$|a|$ 变大时时间窗变宽。

容易求出 $h_{b,a}(t)$ 的变换为

$$
H_{b,a}(\xi)
=
\frac1{\sqrt{2\pi a}}
\int_{-\infty}^{\infty}
h
\left(
\frac{t-b}{a}
\right)
e^{-i\xi t}\,dt
=
\sqrt{\frac a{2\pi}}
e^{-ib\xi}
H(a\xi).
\tag{18.88}
$$

设 $H(\xi)$ 的中心和半宽度分别为 $\bar\omega$ 和 $\Delta_H$，令

$$
\eta(\xi-\bar\omega)
\equiv
H(\xi),
\tag{18.89}
$$

则 $\eta(\xi)$ 是中心和半宽度分别为 0 和 $\Delta_H$ 的窗函数。根据 Parseval 方程，可以写出

$$
(\mathcal W_h f)(b,a)
=
\sqrt{\frac a{2\pi}}
\int_{-\infty}^{\infty}
\eta^*(a\xi-\bar\omega)
F(\xi)\,d\xi.
\tag{18.90}
$$

由于 $\eta(a\xi-\bar\omega)$ 的半宽度为 $\Delta_H/|a|$，所以，(18.90) 式给出的又是 $H(\xi)$ 在频率窗

$$
\left[
\frac{\bar\omega}{a}
-
\frac1{|a|}\Delta_H,
\,
\frac{\bar\omega}{a}
+
\frac1{|a|}\Delta_H
\right]
$$

内的局部信息。这样的时间-频率窗就是

$$
[
b+a\bar t-|a|\Delta_h,
\,
b+a\bar t+|a|\Delta_h
]
\times
\left[
\frac{\bar\omega}{|a|}
-
\frac1{|a|}\Delta_H,
\,
\frac{\bar\omega}{|a|}
+
\frac1{|a|}\Delta_H
\right],
$$

宽度为 $2|a|\Delta_h$。因此，积分小波变换具有“变焦”特性：在检测高频现象（即 $|a|$ 小）时窗自动变窄；在检测低频现象（即 $|a|$ 大）时窗自动变宽（见图 18.2）。正是因为这种变焦特性，使得积分小波变换成为许多理论研究以及工程技术应用的有力工具。

![Figure 18.2](...)

图 18.2　小波变换的时间-频率窗

再进一步，由 $f(t)$ 的积分小波变换 $(\mathcal W_h f)(b,a)$ 值还可以重构 $f(t)$（即反演），

$$
f(t)
=
\frac1{C_h}
\int_{-\infty}^{\infty}
\left[
\int_{-\infty}^{\infty}
(\mathcal W_h f)(b,a)
h_{b,a}(t)\,db
\right]
\frac{da}{a^2}.
\tag{18.91}
$$

证明从略。从这个结果，就可以理解为什么对基小波要加上相容性条件的限制。

在实际应用中，还有离散形式的积分小波变换。例如，通过 Haar 函数

$$
h(t)
=
\begin{cases}
1,
&
0\le t<1/2,
\\
-1,
&
1/2\le t<1,
\\
0,
&
t<0\text{ 或 }t\ge1
\end{cases}
\tag{18.92}
$$

的二进制伸缩与平移，得到

$$
h_{j,k}(t)
=
2^{-j/2}
h(2^{-j}t-k),
\qquad
j,k\text{ 为任意整数},
\tag{18.93}
$$

就构成一组正交归一基

$$
(h_{j,k},h_{l,m})
=
\delta_{jl}\delta_{km},
\tag{18.94}
$$

任意一个平方可积的函数 $f(t)$ 都可以（在平均收敛意义下）展开为小波级数，

$$
f(t)
=
\sum_{j,k=-\infty}^{\infty}
c_{j,k}h_{j,k}(t),
\tag{18.95}
$$

展开系数 $c_{j,k}$ 为

$$
c_{j,k}
=
(h_{j,k},f)
=
(\mathcal W_h f)
\left(
\frac{k}{2^j},
\frac1{2^j}
\right).
\tag{18.96}
$$

小波分析及其应用是一门新兴的学科。尽管它的早期思想可以追溯到 20 世纪的上半叶，但是，只是在近一、二十年间才得到蓬勃的发展。它已经或将要广泛应用于信号处理、图像处理、地质勘探、语音识别与合成、雷达、CT 成像、机械故障诊断等等科技领域。本节的介绍只是关于小波分析的一点皮毛。作者摘编这个阅读材料，目的只是希望有读者因此而萌发了解与钻研小波分析的兴趣。

## 习题

1. 用 Laplace 变换求解半无界问题：

$$
\frac{\partial u}{\partial t}
-
\kappa\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
x>0,\ t>0,
$$

$$
\left.u\right|_{x=0}
=
u_0,
\qquad
\left.u\right|_{x\to\infty}
\text{ 有界},
\qquad
t>0,
$$

$$
\left.u\right|_{t=0}
=
0,
\qquad
x>0.
$$

2. 设有两个半无界杆，温度分别为 0 和 $u_0$，在 $t=0$ 时将两杆端点相接，求 $t>0$ 时杆中各点的温度分布。

3. 利用 Laplace 变换求解第 13 章习题第 11 题。

4. 求解下列一维半无界弦上的波动问题：

$$
\begin{cases}
\displaystyle
\frac{\partial^2u(x,t)}{\partial t^2}
-
a^2
\frac{\partial^2u(x,t)}{\partial x^2}
=
0,
&
0<x<\infty,\ t>0,
\\[10pt]
\left.u(x,t)\right|_{x=0}
=
0,
\qquad
\left.u(x,t)\right|_{x\to\infty}
\text{ 有界},
&
t>0,
\\[10pt]
\left.u(x,t)\right|_{t=0}
=
Ae^{-\alpha x}\sin\beta x,
\qquad
\left.
\frac{\partial u(x,t)}{\partial t}
\right|_{t=0}
=
0,
&
0<x<\infty.
\end{cases}
$$

其中 $A,a,\alpha$ 和 $\beta$ 均为已知的与 $x$ 和 $t$ 都无关的正数。

5. 用 Fourier 变换方法求解一维无界弦上的受迫振动问题

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
f(x,t),
$$

$$
\left.u\right|_{t=0}
=
\phi(x),
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
\psi(x).
$$

6. 电子光学中常遇到一种简单的静电透镜——等径双筒镜，它的两极是两个无限接近的等径（设为 $a$）同轴长圆筒，其电势分别为 $V_0$ 和 $-V_0$。求筒内的静电势。

提示：先在边界条件

$$
\left.u\right|_{r=a}
=
V_0e^{-k|z|}
\operatorname{sgn}z
$$

下利用 Fourier 变换求解，而后令 $k\to0$。
