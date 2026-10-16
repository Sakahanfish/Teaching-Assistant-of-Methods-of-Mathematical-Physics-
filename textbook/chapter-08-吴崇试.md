# 第八章 Laplace 变换

Laplace 变换是一种常用的积分变换。在数学、物理及工程科学中有广泛的应用。本章介绍 Laplace 变换的定义及其基本性质，以及它的简单应用。

## §8.1 Laplace 变换的定义

Laplace 变换（简称拉氏变换）是一种积分变换，它把定义在正实轴上的函数 $f(t)$ 变换为定义在复平面上的函数 $F(p)$，

$$
F(p)=\int_0^\infty e^{-pt}f(t)\,dt.
\tag{8.1}
$$

这里的 $t>0$ 是非负实数，$p=s+i\sigma$ 是复数。$F(p)$ 称为 $f(t)$ 的 Laplace 换式，简称拉氏换式。$e^{-pt}$ 是 Laplace 变换的核。通常把 Laplace 变换简写为

$$
F(p)=\mathcal L\{f(t)\}
\qquad\text{或}\qquad
F(p)\doteq f(t);
\tag{8.2}
$$

$$
f(t)=\mathcal L^{-1}\{F(p)\}
\qquad\text{或}\qquad
f(t)\doteq F(p).
\tag{8.3}
$$

$f(t)$ 和 $F(p)$ 有时也分别称为 Laplace 变换的原函数和像函数。已知原函数 $f(t)$，求像函数 $F(p)$ 称为对原函数 $f(t)$ 作 Laplace 变换；已知像函数 $F(p)$，求原函数 $f(t)$ 的运算称为反演。

需要说明，在本章中约定：当 $t<0$ 时应该理解为 $f(t)=0$，或者说，应该将 $f(t)$ 理解为 $f(t)\eta(t)$，其中

$$
\eta(t)=
\begin{cases}
1, & t>0,\\
0, & t<0
\end{cases}
\tag{8.4}
$$

是 Heaviside 的单位阶跃函数。相应地，$f(t-\tau)$ 就应该理解为 $f(t-\tau)\eta(t-\tau)$。

**练习 8.1**　证明：

$$
f(t-\tau)\doteq e^{-p\tau}F(p),\qquad \tau>0;
$$

$$
f(at)\doteq\frac1aF\left(\frac pa\right),\qquad a>0;
$$

$$
e^{p_0t}f(t)\doteq F(p-p_0).
$$

**例 8.1**　函数 $f(t)=1$ 的 Laplace 换式为

$$
1\doteq\int_0^\infty e^{-pt}\,dt
=
-\left.\frac1p e^{-pt}\right|_0^\infty
=
\frac1p,
\qquad
\operatorname{Re}p>0.
\tag{8.5}
$$

这里的限制条件 $\operatorname{Re}p>0$ 是为了保证积分收敛，换言之，也就是函数 $\eta(t)$ 的 Laplace 变换存在的条件。

**例 8.2**　函数 $f(t)=e^{\alpha t}$ 的 Laplace 换式为

$$
e^{\alpha t}\doteq
\int_0^\infty e^{-pt}e^{\alpha t}\,dt
=
\frac1{p-\alpha},
\qquad
\operatorname{Re}p>\operatorname{Re}\alpha.
\tag{8.6}
$$

这里的限制条件 $\operatorname{Re}p>\operatorname{Re}\alpha$，同样是为了保证积分收敛，即保证函数 $e^{\alpha t}\eta(t)$ 的 Laplace 变换存在。

$f(t)$ 的 Laplace 变换存在的条件就是积分

$$
\int_0^\infty e^{-pt}f(t)\,dt
$$

收敛的条件。由于 Laplace 变换的核是 $e^{-pt}$，所以对于相当广泛的函数 $f(t)$，其拉氏换式都存在，甚至当 $t\to\infty,\ f(t)\to\infty$ 时（例如例 8.2 中的 $e^{\alpha t},\ \operatorname{Re}\alpha>0$），拉氏换式也可能存在。

在绝大多数实际问题中，$f(t)$ 都能满足

(1) $f(t)$ 在区间 $0\leqslant t<\infty$ 上除第一类间断点外都是连续的，且有连续导数，在任何有限区间内这种间断点的数目是有限的；

(2) $f(t)$ 有有限的增长指数，即存在正数 $M>0,\ t_0>0$ 及 $s'\geqslant0$，使对于任何 $t>t_0$，

$$
|f(t)|<Me^{s't}.
\tag{8.7}
$$

这是 Laplace 变换存在的充分条件。一般物理问题中遇到的函数都能满足这个要求。

当然，如果 $s'$ 存在，则 $s'$ 并不唯一，因为比 $s'$ 大的任何正数显然都符合要求。符合要求的 $s'$ 的下界称为收敛横标，记为 $s_0$。

## §8.2 Laplace 变换的基本性质

Laplace 变换具有下列基本性质：

**性质 1**　Laplace 变换是一个线性变换，即若 $f_1(t)\doteq F_1(p),\ f_2(t)\doteq F_2(p)$，则对任意复常数 $\alpha_1$ 和 $\alpha_2$，

$$
\boxed{
\alpha_1f_1(t)+\alpha_2f_2(t)
\doteq
\alpha_1F_1(p)+\alpha_2F_2(p)
}
\tag{8.8}
$$

这个性质很容易从 Laplace 变换的定义得到，因为它只不过是积分运算线性性质的反映。

根据这个性质，立即得到

$$
\sin\omega t
=
\frac{e^{i\omega t}-e^{-i\omega t}}{2i}
\doteq
\frac1{2i}
\left(
\frac1{p-i\omega}
-
\frac1{p+i\omega}
\right)
=
\frac{\omega}{p^2+\omega^2};
\tag{8.9}
$$

$$
\cos\omega t
=
\frac{e^{i\omega t}+e^{-i\omega t}}2
\doteq
\frac12
\left(
\frac1{p-i\omega}
+
\frac1{p+i\omega}
\right)
=
\frac{p}{p^2+\omega^2}.
\tag{8.10}
$$

**性质 2**　Laplace 换式的解析性。

如果函数 $f(t)$ 满足 Laplace 变换存在的充分条件，则

$$
|e^{-pt}f(t)|
<
Me^{-(s-s_0)t},
\qquad
s=\operatorname{Re}p.
$$

当 $s-s_0\geqslant\delta>0$ 时，

$$
|e^{-pt}f(t)|<Me^{-\delta t}.
$$

而积分

$$
\int_0^\infty Me^{-\delta t}\,dt
$$

收敛，故

$$
\int_0^\infty e^{-pt}f(t)\,dt
$$

在 $\operatorname{Re}p\geqslant s_0+\delta$ 上一致收敛，因而在 $\operatorname{Re}p>s_0$ 的半平面内代表一个解析函数，即 $F(p)$ 在半平面 $\operatorname{Re}p>s_0$ 内解析。

这个性质可以用来确定收敛横标 $s_0$。

**性质 3**　若 $f(t)$ 满足 Laplace 变换存在的充分条件，则

$$
\boxed{
F(p)\to0,
\qquad
\text{当 }\operatorname{Re}p=s\to+\infty.
}
\tag{8.11}
$$

**证**　因为

$$
|F(p)|
\leqslant
\int_0^\infty|e^{-pt}f(t)|\,dt
\leqslant
M\int_0^\infty e^{-(s-s_0)t}\,dt
=
\frac{M}{s-s_0},
\qquad
s>s_0,
$$

故当 $\operatorname{Re}p=s\to+\infty$ 时，$F(p)\to0$。$\square$

此外，由 Riemann-Lebesgue 定理①还可证明，当 $\operatorname{Re}p=s>s_0$ 时，

$$
\boxed{
\lim_{\operatorname{Im}p\to\pm\infty}F(p)=0.
}
$$

**性质 4**　原函数导数的 Laplace 变换。设 $f(t)$ 及 $f'(t)$ 的 Laplace 变换均存在，$f(t)\doteq F(p)$，则因为

$$
\int_0^\infty f'(t)e^{-pt}\,dt
=
\left.f(t)e^{-pt}\right|_0^\infty
+
p\int_0^\infty f(t)e^{-pt}\,dt,
$$

所以

$$
\boxed{
f'(t)\doteq pF(p)-f(0).
}
\tag{8.12}
$$

同样，只要 $f(t),f'(t),\cdots,f^{(n)}(t)$ 的 Laplace 变换均存在，$f(t)\doteq F(p)$，则

$$
f''(t)\doteq
p^2F(p)-pf(0)-f'(0),
\tag{8.13a}
$$

$$
f^{(3)}(t)\doteq
p^3F(p)-p^2f(0)-pf'(0)-f''(0),
\tag{8.13b}
$$

$$
f^{(n)}(t)\doteq
p^nF(p)-p^{n-1}f(0)-p^{n-2}f'(0)-\cdots-pf^{(n-2)}(0)-f^{(n-1)}(0).
\tag{8.13c}
$$

这个性质表明，对原函数 $f(t)$ 的微商运算就转化为对像函数 $F(p)$ 的乘法运算，而且还自动包括了 $f(t)$ 的初值。所以 Laplace 变换是求解常微分方程初值问题的一种重要方法。

**例 8.3**　$LR$ 串联电路（图 8.1），假设在 K 合上前电路中没有电流，求 K 合上后电路中的电流。

**解**　根据 Kirchhoff 定律，可列出微分方程

$$
L\frac{di}{dt}+Ri=E,
\tag{8.14a}
$$

$$
i(0)=0.
\tag{8.14b}
$$

![Figure 8.1](...)

图 8.1　$LR$ 串联电路

设 $i(t)\doteq I(p)$，则

$$
\frac{di}{dt}\doteq pI(p)-i(0)=pI(p).
$$

所以

$$
LpI(p)+RI(p)=\frac Ep.
\tag{8.15}
$$

这样，经过 Laplace 变换，求解原函数 $i(t)$ 的常微分方程的定解问题 (8.14) 就转化为求解像函数 $I(p)$ 的代数方程 (8.15)，易得

$$
I(p)
=
\frac Ep\frac1{Lp+R}
=
\frac ER
\left(
\frac1p-\frac{L}{Lp+R}
\right).
$$

反演，得

$$
i(t)
=
\frac ER
\left[
1-e^{-(R/L)t}
\right].
\tag{8.16}
$$

Laplace 变换不仅仅可以用于求解常微分方程的初值问题，也可以用来求解半无界空间常微分方程的边值问题。

**例 8.4**　求解下列定解问题

$$
y''(t)-y'(t)-2y(t)=0,
\qquad
t>0,
\tag{8.17a}
$$

$$
y(0)=1,
\qquad
y(t)\big|_{t\to\infty}\text{ 有界}.
\tag{8.17b}
$$

**解**　设 $y(t)$ 的 Laplace 变换存在，$y(t)\doteq Y(p)$，相应地，

$$
y'(t)\doteq pY(p)-1,
\qquad
y''(t)\doteq p^2Y(p)-p-y'(0).
$$

在原题中并未给出 $y'(0)$ 值，不妨当作待定常数，下面我们将看到，它可以根据无穷远条件（$t\to\infty$ 时 $y(t)$ 有界）确定。于是，经过 Laplace 变换后，定解问题 (8.17) 变为

$$
[p^2Y(p)-p-y'(0)]
-
[pY(p)-1]
-
2Y(p)
=
0,
$$

解之得

$$
\begin{aligned}
Y(p)
&=
\frac{p-1+y'(0)}{p^2-p-2}\\
&=
\frac13
\left(
\frac1{p-2}+\frac2{p+1}
\right)
+
\frac13y'(0)
\left(
\frac1{p-2}-\frac1{p+1}
\right).
\end{aligned}
$$

求反演，有

$$
y(t)
=
\frac13[2-y'(0)]e^{-t}
+
\frac13[1+y'(0)]e^{2t}.
$$

因 $t\to\infty$ 时 $e^{2t}$ 无界，故由无穷远条件，可以定出 $y'(0)=-1$。所以最后就求得

$$
y(t)=e^{-t}.
$$

**性质 5**　原函数积分的 Laplace 变换。设 $f(t)$ 满足 Laplace 变换存在的充分条件，则

$$
\left|
\int_0^t f(\tau)\,d\tau
\right|
\leqslant
\int_0^t|f(\tau)|\,d\tau
\leqslant
\int_0^tMe^{s_0\tau}\,d\tau
=
\frac{M}{s_0}(e^{s_0t}-1),
$$

所以 $\int_0^t f(\tau)\,d\tau$ 的 Laplace 变换也存在，记

$$
f(t)\doteq F(p),
\qquad
\int_0^t f(\tau)\,d\tau
\doteq
\mathcal L
\left\{
\int_0^t f(\tau)\,d\tau
\right\}.
$$

因为

$$
\frac{d}{dt}
\int_0^t f(\tau)\,d\tau
=
f(t),
$$

根据性质 4，就有

$$
F(p)
=
p\mathcal L
\left\{
\int_0^t f(\tau)\,d\tau
\right\}
-0.
$$

所以

$$
\boxed{
\int_0^t f(\tau)\,d\tau
\doteq
\frac{F(p)}p.
}
\tag{8.18}
$$

**练习 8.2**　证明：

$$
\int_0^\infty f(t,\tau)\,d\tau
\doteq
\int_0^\infty F(p,\tau)\,d\tau,
$$

$$
\frac{\partial f(x,t)}{\partial x}
\doteq
\frac{dF(x,p)}{dx},
$$

$$
\int_t^\infty\frac{f(\tau)}{\tau}\,d\tau
\doteq
\frac1p\int_0^pF(q)\,dq,
$$

假定有关的积分均存在。

**例 8.5**　求 $LC$ 串联电路（见图 8.2）中电流随时间的变化。

**解**　因为

$$
\frac{q(t)}C
=
L\frac{di}{dt},
\tag{8.19a}
$$

$$
q(t)
=
-\int_0^ti(\tau)\,d\tau+q_0,
\tag{8.19b}
$$

所以

$$
L\frac{di}{dt}
+
\frac1C
\int_0^ti(\tau)\,d\tau
=
\frac{q_0}{C}.
$$

![Figure 8.2](...)

图 8.2　$LC$ 并联电路

这是关于未知函数 $i(t)$ 的微分积分方程。设 $i(t)\doteq I(p)$，则有

$$
LpI(p)
+
\frac1C\frac{I(p)}p
=
\frac{q_0}{C}\frac1p.
$$

所以求解微分积分方程的定解问题 (8.19) 也转化为求解代数方程

$$
I(p)
=
\frac{q_0}{LCp^2+1}.
$$

利用 (8.9) 式求反演，即得

$$
i(t)
=
\frac{q_0}{\sqrt{LC}}
\sin\frac{t}{\sqrt{LC}}.
\tag{8.20}
$$

还值得对上面得到的 (8.18) 式做进一步的分析。如果我们将该式改写成

$$
\int_0^\infty f(t)e^{-pt}\,dt
=
p\int_0^\infty g(t)e^{-pt}\,dt,
\qquad
\text{其中 }
g(t)=\int_0^tf(\tau)\,d\tau,
\tag{8.18'}
$$

上面的推导表明，只要 $f(t)$ 满足 Laplace 变换存在的充分条件，则 $g(t)$ 的 Laplace 变换就存在；反过来说，只要 $g(t)$ 的 Laplace 变换存在，而且 $g(t)$ 可导，则 $g'(t)$ 就一定满足 Laplace 变换存在的充分条件。因此我们在讨论 Laplace 变换的性质时，包括下面在讨论 Laplace 变换的反演时，常常假设原函数满足 Laplace 变换存在的充分条件（或者等价地，假设像函数在 $\operatorname{Re}p\geqslant s_1>s_0$ 的半平面上解析）。这一要求实际上可以减弱，在通常情况下，要求原函数的 Laplace 变换存在即可。

① Riemann-Lebesgue 定理：如果积分

$$
\int_a^bf(t)\,dt
$$

存在，则

$$
\lim_{\omega\to\infty}
\int_a^bf(t)\sin\omega t\,dt
=
0,
\qquad
\lim_{\omega\to\infty}
\int_a^bf(t)\cos\omega t\,dt
=
0.
$$

读者可参阅参考书目 [7] 的（中译本）第 449 页。

## §8.3 Laplace 变换的反演

求 Laplace 变换的反演，首先就遇到反演的唯一性问题，即对于任意给定的像函数 $F(p)$，是否可能存在不止一个原函数？例如是否可能存在两个不同的原函数 $f_1(t)$ 和 $f_2(t)$，使得

$$
f_1(t)\doteq F(p),
\qquad
f_2(t)\doteq F(p).
\tag{8.21}
$$

或者说，是否存在函数

$$
g(t)\equiv f_1(t)-f_2(t)\ne0,
$$

使得

$$
g(t)\doteq0.
\tag{8.22}
$$

回答是：如果限定 $g(t)$ 为连续函数，则由 (8.22) 式一定能推出 $g(t)=0$；但若许可 $g(t)$ 不连续，则 $g(t)$ 可以不恒为 0，但是几乎处处为 0（例如在有限个点处不为 0，而在其余点均为 0）。如果我们也将这种函数称为“零函数”，则对应于同一个 $F(p)$，它的原函数可以有无穷多个，彼此间相差一个零函数。但是，如果限定原函数为连续函数，则 Laplace 变换的反演具有唯一性。以下的讨论中，我们将约定原函数均为连续函数。

**像函数导数的反演**　设 $f(t)$ 满足 Laplace 变换存在的充分条件，$f(t)\doteq F(p)$，则 $F(p)$ 在 $\operatorname{Re}p\geqslant s_1>s_0$ 的半平面上解析，因而可以在积分号下求导：

$$
F^{(n)}(p)
=
\frac{d^n}{dp^n}
\int_0^\infty f(t)e^{-pt}\,dt
=
\int_0^\infty(-t)^nf(t)e^{-pt}\,dt.
$$

所以

$$
F^{(n)}(p)\doteq(-t)^nf(t).
\tag{8.23}
$$

根据这个公式，可以得到

$$
\frac1{p^2}
=
-\frac{d}{dp}\frac1p
\doteq t,
\tag{8.24}
$$

$$
\frac1{p^3}
=
\frac12\frac{d^2}{dp^2}\frac1p
\doteq
\frac12t^2.
\tag{8.25}
$$

这样，若 $F(p)$ 是有理函数，则总可以通过部分分式求反演。例如

$$
\frac{1}{p^3(p+\alpha)}
=
\frac1{\alpha p^3}
-\frac1{\alpha^2p^2}
+\frac1{\alpha^3p}
-\frac1{\alpha^3(p+\alpha)}
$$

$$
\doteq
\frac1{2\alpha}t^2
-\frac1{\alpha^2}t
+\frac1{\alpha^3}
-\frac1{\alpha^3}e^{-\alpha t}.
\tag{8.26}
$$

**像函数的积分的反演**　如果

$$
G(p)
=
\int_p^\infty F(q)\,dq
$$

存在①，且当 $t\to0$ 时，$|f(t)/t|$ 有界，则

$$
\boxed{
\int_p^\infty F(q)\,dq
\doteq
\frac{f(t)}t.
}
\tag{8.27}
$$

**证**　设 $G(p)\doteq g(t)$，则因为 $G'(p)=-F(p)$，故有

$$
(-t)g(t)=-f(t),
$$

因此即证得 (8.27) 式。$\square$

利用这个公式，又可以得到许多函数的 Laplace 变换。例如，

$$
\frac{\sin\omega t}{t}
\doteq
\int_p^\infty\frac{\omega}{q^2+\omega^2}\,dq
=
\frac\pi2-\arctan\frac p\omega.
\tag{8.28}
$$

特别是，如果 $p\to0$ 时，(8.27) 式两端的积分均存在，则有

$$
\boxed{
\int_0^\infty F(p)\,dp
=
\int_0^\infty\frac{f(t)}t\,dt.
}
\tag{8.29}
$$

利用这个结果，可以计算 $\int_0^\infty f(t)/t\,dt$ 型的积分。例如

$$
\int_0^\infty\frac{\sin t}{t}\,dt
=
\int_0^\infty\frac1{p^2+1}\,dp
=
\frac\pi2.
\tag{8.30}
$$

这个积分曾经应用留数定理计算过。但是这里的计算更为简便。

同样，还可以计算

$$
\begin{aligned}
\int_0^\infty\frac{\cos at-\cos bt}{t}\,dt
&=
\int_0^\infty
\left(
\frac{p}{p^2+a^2}
-
\frac{p}{p^2+b^2}
\right)dp\\
&=
\left.
\frac12
\ln\frac{p^2+a^2}{p^2+b^2}
\right|_0^\infty\\
&=
\ln b-\ln a,
\qquad
a>0,\ b>0.
\end{aligned}
\tag{8.31}
$$

**像函数在 $\infty$ 点解析的情形**　现在研究一个更特殊的情形，它特别容易求出原函数。如果可以将 $F(p)$ 由半平面 $\operatorname{Re}p>s_0$（单值地）解析延拓到含有 $p=\infty$ 点在内的一定区域内，且 $F(p)$ 在 $p=\infty$ 点解析，这样，函数 $F(p)$ 就可以在 $p=\infty$ 点作 Taylor 级数展开：

$$
F(p)
=
\sum_{n=1}^\infty c_np^{-n}.
\tag{8.32}
$$

级数中不含 $n=0$ 项，是因为 $F(p)$ 作为 Laplace 换式，应当满足 $\operatorname{Re}p\to+\infty$ 时 $F(p)\to0$ 的要求。将级数逐项求反演，就得到

$$
\boxed{
f(t)
=
\sum_{n=0}^\infty
\frac{c_{n+1}}{n!}t^n.
}
\tag{8.33}
$$

考察这种做法的合法性在于证明级数收敛。为此作圆周 $C_R:|p|=R$，在 $C_R$ 外 $F(p)$ 解析，

$$
c_n
=
\frac{1}{2\pi i}
\oint_{C_R}F(p)p^{n-1}\,dp.
$$

因为 $p=\infty$ 是 $F(p)$ 的零点，故当 $|p|>R$ 时，$|F(p)|<M/R$，因之

$$
|c_n|<MR^{n-1}.
$$

由此即得

$$
\left|
\sum_{n=0}^\infty
\frac{c_{n+1}}{n!}t^n
\right|
\leqslant
\sum_{n=0}^\infty
\frac{|c_{n+1}|}{n!}|t|^n
<
M\sum_{n=0}^\infty
\frac1{n!}R^n|t|^n
=
Me^{R|t|},
$$

故级数收敛。这里同时也证明了 $f(t)$ 具有有限的增长指数，因而 Laplace 变换存在。

应用这个方法可以求出函数 $1/\sqrt{p^2+1}$ 的反演。这是一个多值函数，如果规定单值分支

$$
\left.
\frac1{\sqrt{p^2+1}}
\right|_{p\to\infty}
\to
\frac1p,
$$

则有

$$
\begin{aligned}
\frac1{\sqrt{p^2+1}}
&=
\sum_{k=0}^\infty
(-)^k
\frac{(2k)!}{2^{2k}(k!)^2}
\frac1{p^{2k+1}}\\
&\doteq
\sum_{k=0}^\infty
(-)^k
\frac1{2^{2k}(k!)^2}
t^{2k}\\
&=
\sum_{k=0}^\infty
\frac{(-)^k}{k!\,k!}
\left(\frac t2\right)^{2k}.
\end{aligned}
\tag{8.34}
$$

这正是 §5.5 例 5.9 中见到过的 Bessel 函数 $J_0(t)$。

另一个例子是

$$
\frac1p e^{-1/p}
=
\sum_{n=0}^\infty
(-)^n\frac1{n!}\frac1{p^{n+1}}
\doteq
\sum_{n=0}^\infty
\frac{(-)^n}{n!\,n!}t^n
=
J_0(2\sqrt t).
\tag{8.35}
$$

如果 $F(p)$ 可以分解为 $F_1(p)$ 和 $F_2(p)$ 之积，其反演问题就需要用到下面的卷积定理。

**定理 8.1（卷积定理）**　设 $F_1(p)\doteq f_1(t),\ F_2(p)\doteq f_2(t)$，则

$$
\boxed{
F_1(p)F_2(p)
\doteq
\int_0^t
f_1(\tau)f_2(t-\tau)\,d\tau.
}
\tag{8.36}
$$

**证**　因为

$$
F_1(p)
=
\int_0^\infty f_1(t)e^{-pt}\,dt,
\qquad
F_2(p)
=
\int_0^\infty f_2(t)e^{-pt}\,dt,
$$

所以，

$$
\begin{aligned}
F_1(p)F_2(p)
&=
\int_0^\infty f_1(\tau)e^{-p\tau}\,d\tau
\int_0^\infty f_2(\nu)e^{-p\nu}\,d\nu\\
&=
\int_0^\infty f_1(\tau)\,d\tau
\int_0^\infty f_2(\nu)e^{-p(\tau+\nu)}\,d\nu\\
&=
\int_0^\infty f_1(\tau)\,d\tau
\int_\tau^\infty f_2(t-\tau)e^{-pt}\,dt.
\end{aligned}
$$

在 $Ot\tau$ 平面上画出积分区域（见图 8.3），改变积分次序，即得

$$
F_1(p)F_2(p)
=
\int_0^\infty e^{-pt}\,dt
\int_0^t
f_1(\tau)f_2(t-\tau)\,d\tau,
$$

定理得证。$\square$

![Figure 8.3](...)

图 8.3　卷积定理的证明

![Figure 8.4](...)

图 8.4　$LR$ 串联电路

**例 8.6**　在 $LR$ 串联电路（见图 8.4）中加一方形脉冲电压

$$
\mathcal E(t)=
\begin{cases}
E_0, & 0\leqslant t\leqslant T,\\
0, & t>T.
\end{cases}
\tag{8.37}
$$

求电路中的电流 $i(t)$，设 $i(0)=0$。

**解**　列方程

$$
L\frac{di}{dt}+Ri=\mathcal E(t),
$$

$$
i(0)=0.
$$

作 Laplace 变换：设 $i(t)\doteq I(p),\ \mathcal E(t)\doteq E(p)$，则

$$
LpI(p)+RI(p)=E(p),
\qquad
\text{即}
\qquad
I(p)=\frac1{Lp+R}\,E(p).
$$

所以

$$
\begin{aligned}
i(t)
&=
\int_0^t
\mathcal E(\tau)
\frac1L e^{-R(t-\tau)/L}\,d\tau\\
&=
\begin{cases}
\dfrac{E_0}{R}
\left(
1-e^{-Rt/L}
\right),
& 0\leqslant t\leqslant T,\\[8pt]
\dfrac{E_0}{R}
\left(
e^{RT/L}-1
\right)e^{-Rt/L},
& t>T.
\end{cases}
\end{aligned}
$$

由上述计算过程可以看出，应用卷积定理，在整个求解过程中我们无需求出方程非齐次项 $\mathcal E(t)$ 的像函数 $E(p)$ 的具体表达式。

① 这里的积分上限应理解为 $\operatorname{Re}p\to+\infty$，并且积分路径在 $F(p)$ 的解析区域内，因而积分值与路径无关。

## §8.4 普遍反演公式

若函数 $F(p)=F(s+i\sigma)$ 在区域 $\operatorname{Re}p>s_0$ 内满足：(1) $F(p)$ 解析，(2) 当 $|p|\to\infty$ 时 $F(p)$ 趋于 0，(3) 对于所有的 $\operatorname{Re}p=s>s_0$，沿直线 $L:\operatorname{Re}p=s$ 的无穷积分

$$
\int_{s-i\infty}^{s+i\infty}
|F(p)|\,d\sigma
\qquad
(s>s_0)
$$

收敛，则对于 $\operatorname{Re}p=s>s_0$，$F(p)$ 的原函数为

$$
f(t)
=
\frac1{2\pi i}
\int_{s-i\infty}^{s+i\infty}
F(p)e^{pt}\,dp.
\tag{8.38}
$$

证明从略。读者可参阅参考书目 [29] 的 8.2 节；或见本书（第一版）§10.4。也可转换为 Fourier 变换而证明，见参考书目 [1] 的 10.7 节。

**例 8.7**　用普遍反演公式求 Laplace 换式

$$
F(p)=\frac1{(p^2+\omega^2)^2},
\qquad
\omega>0
$$

的原函数。

**解**　由普遍反演公式，此像函数的原函数为

$$
f(t)
=
\frac1{2\pi i}
\int_{s-i\infty}^{s+i\infty}
\frac1{(p^2+\omega^2)^2}
e^{pt}\,dp.
$$

函数 $1/(p^2+\omega^2)^2$ 的孤立奇点都在虚轴上，所以取 $s>0$ 即可。因为 $t<0$ 时一定有 $f(t)=0$，故只需讨论 $t>0$ 的情形。此时可取围道如图 8.5 所示。由于

$$
\lim_{p\to\infty}
\frac1{(p^2+\omega^2)^2}
=
0,
$$

所以，根据推广的 Jordan 引理①，可以断定

$$
\lim_{R\to\infty}
\int_{C_R}
\frac1{(p^2+\omega^2)^2}
e^{pt}\,dp
=
0.
$$

![Figure 8.5](...)

图 8.5　例 8.7 所用的围道

这样，由留数定理，就得到

$$
\begin{aligned}
f(t)
&=
\frac1{2\pi i}
\int_{s-i\infty}^{s+i\infty}
\frac1{(p^2+\omega^2)^2}
e^{pt}\,dp\\
&=
\sum_{\text{全平面}}
\operatorname{res}
\left\{
\frac1{(p^2+\omega^2)^2}
e^{pt}
\right\}\\
&=
\left\{
\left[
\frac{t}{(p+i\omega)^2}
-
\frac{2}{(p+i\omega)^3}
\right]e^{pt}
\right\}_{p=i\omega}\\
&\quad+
\left\{
\left[
\frac{t}{(p-i\omega)^2}
-
\frac{2}{(p-i\omega)^3}
\right]e^{pt}
\right\}_{p=-i\omega}\\
&=
\frac1{2\omega^3}
(\sin\omega t-\omega t\cos\omega t).
\end{aligned}
\tag{8.39}
$$

**例 8.8**　用普遍反演公式求 Laplace 换式

$$
F(p)
=
\frac1{\sqrt p}e^{-\alpha\sqrt p},
\qquad
\alpha>0
$$

的原函数。

**解**　由普遍反演公式，原函数为

$$
\frac1{2\pi i}
\int_{s-i\infty}^{s+i\infty}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp.
$$

其中的积分路径 $L:\operatorname{Re}p=s>0$ 是右半平面上的一条平行于虚轴的无穷直线。考虑到被积函数是多值函数，$p=0$ 和 $p=\infty$ 是分支点，为保证像函数在积分路径 $L$ 右侧解析，从 $p=0$ 沿负实轴到 $p=\infty$ 作割线。因此，在应用留数定理计算这个积分时，应该取积分围道如图 8.6 所示。因为在积分围道内解析，所以

$$
\begin{aligned}
\oint_C
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
={}&
\int_A^B
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
+
\int_{C_R}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp\\
&+
\int_{C_1}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
+
\int_{C_\delta}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp\\
&+
\int_{C_2}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
+
\int_{C_R'}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
=
0.
\end{aligned}
$$

![Figure 8.6](...)

图 8.6　例 8.8 所用的围道

在 $C_1$ 和 $C_2$ 上，$\arg p=\pm\pi$，故可分别令 $p=re^{\pm i\pi}$ 而得到

$$
\int_{C_1}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
=
-i\int_\delta^R
\frac1{\sqrt r}
e^{-i\alpha\sqrt r}e^{-rt}\,dr,
$$

$$
\int_{C_2}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
=
-i\int_\delta^R
\frac1{\sqrt r}
e^{i\alpha\sqrt r}e^{-rt}\,dr.
$$

又，由推广的 Jordan 引理①，可知

$$
\lim_{R\to\infty}
\int_{C_R}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
=
0,
$$

$$
\lim_{R\to\infty}
\int_{C_R'}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
=
0.
$$

根据引理 3.4，也有

$$
\lim_{\delta\to0}
\int_{C_\delta}
\frac1{\sqrt p}
e^{-\alpha\sqrt p}e^{pt}\,dp
=
0.
$$

所以，在取极限 $R\to\infty,\ \delta\to0$ 后，就有

$$
\begin{aligned}
\frac1{\sqrt p}e^{-\alpha\sqrt p}
&\doteq
\frac1{2\pi}
\int_0^\infty
\frac1{\sqrt r}
\left[
e^{i\alpha\sqrt r}
+
e^{-i\alpha\sqrt r}
\right]e^{-rt}\,dr\\
&=
\frac2\pi
\int_0^\infty
e^{-x^2t}\cos\alpha x\,dx,
\end{aligned}
$$

代入 (4.29) 式的结果即得

$$
\frac1{\sqrt p}e^{-\alpha\sqrt p}
\doteq
\frac1{\sqrt{\pi t}}
e^{-\alpha^2/4t}.
\tag{8.40}
$$

① 这里所谓推广的 Jordan 引理，指的是将原始的 Jordan 引理（见 §6.4）作了如下的变化与扩充：

(1) 将 Jordan 引理所讨论的圆弧旋转了 $90^\circ$。

(2) 现在的圆弧是和直线 $L:\operatorname{Re}p=s>0$ 相交，所以要略大于半圆弧。但可以证明（从略），只要圆弧与虚轴的距离（即 $s$）固定（因而当圆弧的半径 $R\to\infty$ 时，虚轴右方圆弧的张角 $\to0$），则引理仍然成立。

在这个结果的基础上就可以证明下列公式：

$$
\begin{aligned}
\frac1{\sqrt{\pi t}}
\int_0^\infty
f(\tau)e^{-\tau^2/4t}\,d\tau
&\doteq
\int_0^\infty
e^{-pt}
\left\{
\frac1{\sqrt{\pi t}}
\int_0^\infty
f(\tau)e^{-\tau^2/4t}\,d\tau
\right\}dt\\
&=
\int_0^\infty
f(\tau)
\left\{
\int_0^\infty
\frac1{\sqrt{\pi t}}
e^{-\tau^2/4t}e^{-pt}\,dt
\right\}d\tau\\
&=
\int_0^\infty
f(\tau)
\frac1{\sqrt p}
e^{-\tau\sqrt p}\,d\tau\\
&=
\frac1{\sqrt p}
F(\sqrt p).
\end{aligned}
\tag{8.41}
$$

**例 8.9**　求

$$
\frac1p e^{-\alpha\sqrt p},
\qquad
\alpha>0
$$

的原函数。

**解**　在 (8.41) 式中代入

$$
F(p)=\frac1p e^{-\alpha p},
$$

并注意

$$
\frac1p e^{-\alpha p}
\doteq
\eta(t-\alpha),
$$

就有

$$
\frac1p e^{-\alpha\sqrt p}
\doteq
\frac1{\sqrt{\pi t}}
\int_0^\infty
\eta(\tau-\alpha)e^{-\tau^2/4t}\,d\tau
=
\frac1{\sqrt{\pi t}}
\int_\alpha^\infty
e^{-\tau^2/4t}\,d\tau.
$$

由此即得

$$
\frac1p e^{-\alpha\sqrt p}
\doteq
\operatorname{erfc}
\frac{\alpha}{2\sqrt t},
\tag{8.42}
$$

其中的 $\operatorname{erfc}x$ 称为余误差函数，定义为

$$
\operatorname{erfc}x
=
\frac2{\sqrt\pi}
\int_x^\infty e^{-\xi^2}\,d\xi.
\tag{8.43}
$$

相关的还有误差函数 $\operatorname{erf}x$，

$$
\operatorname{erf}x
=
1-\operatorname{erfc}x
=
\frac2{\sqrt\pi}
\int_0^x e^{-\xi^2}\,d\xi.
\tag{8.44}
$$

## *§8.5 利用 Laplace 变换计算级数和

Laplace 变换也可以用来计算某些级数 $\sum F(n)$ 之和。基本思路是利用 Laplace 变换

$$
F(p)
=
\int_0^\infty f(t)e^{-pt}\,dt,
$$

将级数的通项表示成积分，而后交换积分和级数求和的次序，

$$
\sum F(n)
=
\sum\int_0^\infty f(t)e^{-nt}\,dt
=
\int_0^\infty
f(t)
\left(
\sum e^{-nt}
\right)dt,
\tag{8.45}
$$

从而将级数求和的问题转化为计算定积分。由于函数 $e^{-nt}$ 的存在，一般情况下常常可以保证交换积分与求和次序是合法的。

在计算中常用到的 Laplace 变换的结果有

$$
\int_0^\infty
e^{\alpha t}e^{-pt}\,dt
=
\frac1{p-\alpha},
$$

$$
\int_0^\infty
t^{\alpha-1}e^{-pt}\,dt
=
\frac{\Gamma(\alpha)}{p^\alpha},
$$

$$
\int_0^\infty
e^{-pt}\sin\omega t\,dt
=
\frac{\omega}{p^2+\omega^2},
$$

$$
\int_0^\infty
e^{-pt}\cos\omega t\,dt
=
\frac{p}{p^2+\omega^2},
$$

$$
\int_0^\infty
e^{-pt}\sinh at\,dt
=
\frac{a}{p^2-a^2},
$$

$$
\int_0^\infty
e^{-pt}\cosh at\,dt
=
\frac{p}{p^2-a^2}.
$$

**例 8.10**　求级数

$$
\sum_{n=1}^\infty\frac1{n^2}
$$

之和。

**解**　此级数和在第 6 章已经计算过（见例 6.14），现在改用 Laplace 变换方法计算。

容易将待求级数化为积分：

$$
\sum_{n=1}^\infty\frac1{n^2}
=
\sum_{n=1}^\infty
\int_0^\infty xe^{-nx}\,dx
=
\int_0^\infty
\frac{xe^{-x}}{1-e^{-x}}\,dx
=
\int_0^\infty
\frac{x}{e^x-1}\,dx.
$$

下面就用留数定理计算这个积分。取积分围道如图 8.7 所示，则

$$
\begin{aligned}
\oint_C\frac{z^2}{e^z-1}\,dz
={}&
\int_0^R\frac{x^2}{e^x-1}\,dx
+
\int_0^{2\pi}
\frac{(R+iy)^2}{e^{R+iy}-1}i\,dy\\
&+
\int_R^\delta
\frac{(x+2\pi i)^2}{e^x-1}\,dx
+
\int_{C_\delta}
\frac{z^2}{e^z-1}\,dz\\
&+
\int_{2\pi-\delta}^{0}
\frac{(iy)^2}{e^{iy}-1}i\,dy
=
0.
\end{aligned}
$$

取极限 $R\to\infty,\ \delta\to0$，因为

$$
\lim_{R\to\infty}
\int_0^{2\pi}
\frac{(R+iy)^2}{e^{R+iy}-1}i\,dy
=
0,
$$

$$
\lim_{\delta\to0}
\int_{C_\delta}
\frac{z^2}{e^z-1}\,dz
=
-\frac{\pi i}{2}(2\pi i)^2
=
2\pi^3i,
$$

所以就得到

$$
\int_0^\infty
\frac{4\pi ix-4\pi^2}{e^x-1}\,dx
-
\int_0^{2\pi}
\frac{y^2}{e^{iy}-1}i\,dy
=
2\pi^3i.
$$

![Figure 8.7](...)

图 8.7　例 8.10 和例 8.11 所用的围道

注意

$$
\int_0^{2\pi}
\frac{y^2}{e^{iy}-1}i\,dy
=
\int_0^{2\pi}
\frac{y^2e^{-iy/2}}
{2\sin(y/2)}\,dy,
$$

所以比较上式的虚部，即得

$$
4\pi
\int_0^\infty
\frac{x}{e^x-1}\,dx
+
\frac12
\int_0^{2\pi}y^2\,dy
=
2\pi^3,
$$

于是就求得

$$
\int_0^\infty
\frac{x}{e^x-1}\,dx
=
-\frac13\pi^2+\frac12\pi^2
=
\frac{\pi^2}{6},
$$

亦即

$$
\sum_{n=1}^\infty\frac1{n^2}
=
\frac{\pi^2}{6}.
\tag{8.46}
$$

**例 8.11**　求级数

$$
\sum_{n=1}^\infty\frac1{n^4}
$$

之和。

**解**　类似于上题，我们有

$$
\sum_{n=1}^\infty\frac1{n^4}
=
\frac16
\int_0^\infty
\frac{x^3}{e^x-1}\,dx.
$$

考虑复变积分

$$
\oint_C\frac{z^4}{e^z-1}\,dz,
$$

积分围道 $C$ 如图 8.7 所示。因为在积分内无奇点，所以

$$
\begin{aligned}
&\int_0^R\frac{x^4}{e^x-1}\,dx
+
\int_0^{2\pi}
\frac{(R+iy)^4}{e^{R+iy}-1}i\,dy
+
\int_R^\delta
\frac{(x+2\pi i)^4}{e^x-1}\,dx\\
&\qquad+
\int_{C_\delta}
\frac{z^4}{e^z-1}\,dz
+
\int_{2\pi-\delta}^{0}
\frac{(iy)^4}{e^{iy}-1}i\,dy
=
0.
\end{aligned}
$$

取极限 $R\to\infty,\ \delta\to0$，因为

$$
\lim_{R\to\infty}
\int_0^{2\pi}
\frac{(R+iy)^4}{e^{R+iy}-1}i\,dy
=
0,
$$

$$
\lim_{\delta\to0}
\int_{C_\delta}
\frac{z^4}{e^z-1}\,dz
=
-\frac{\pi i}{2}(2\pi i)^4
=
-8\pi^5i,
$$

所以就得到

$$
\int_0^\infty
\frac{
4x^3\,2\pi i
+
6x^2(2\pi i)^2
+
4x(2\pi i)^3
+
(2\pi i)^4
}{e^x-1}\,dx
-
\int_0^{2\pi}
\frac{(iy)^4}{e^{iy}-1}i\,dy
=
-8\pi^5i.
$$

注意

$$
\int_0^{2\pi}
\frac{y^4}{e^{iy}-1}i\,dy
=
\int_0^{2\pi}
\frac{y^4e^{-iy/2}}
{2\sin(y/2)}\,dy,
$$

所以比较上式的虚部，即得

$$
8\pi
\int_0^\infty
\frac{x^3}{e^x-1}\,dx
-
32\pi^3
\int_0^\infty
\frac{x}{e^x-1}\,dx
-
\frac12
\int_0^{2\pi}y^4\,dy
=
-8\pi^5,
$$

于是就求得

$$
\begin{aligned}
\int_0^\infty
\frac{x^3}{e^x-1}\,dx
&=
-\pi^4
+
4\pi^2
\int_0^\infty
\frac{x}{e^x-1}\,dx
+
\frac25\pi^4\\
&=
-\pi^4
+
4\pi^2\cdot\frac{\pi^2}{6}
+
\frac25\pi^4
=
\frac1{15}\pi^4,
\end{aligned}
$$

亦即

$$
\sum_{n=1}^\infty\frac1{n^4}
=
\frac{\pi^4}{90}.
\tag{8.47}
$$

**例 8.12**　计算级数

$$
\sum_{n=0}^\infty
\frac1{n^2-a^2}
$$

之和，其中 $a$ 不为整数，且不妨设 $\operatorname{Re}a>0$。

**解**　这个级数和实际上在 §7.5（见 (7.44) 式）也已经遇到过，当时是作为 Γ 函数无穷乘积表示的直接应用而得到的。这里再采用 Laplace 变换的方法讨论。因为

$$
\int_0^\infty
e^{-pt}\sinh at\,dt
=
\frac{a}{p^2-a^2},
\qquad
\operatorname{Re}p>\operatorname{Re}a,
$$

所以在 $\operatorname{Re}a<1$ 的条件下，级数可化为

$$
\begin{aligned}
\sum_{n=0}^\infty
\frac1{n^2-a^2}
&=
-\frac1{a^2}
+
\sum_{n=1}^\infty
\frac1a
\int_0^\infty
e^{-nt}\sinh at\,dt\\
&=
-\frac1{a^2}
+
\frac1a
\int_0^\infty
\frac{\sinh at}{e^t-1}\,dt.
\end{aligned}
$$

可以直接算出上式左端的积分，

$$
\begin{aligned}
\int_0^\infty
\frac{\sinh at}{e^t-1}\,dt
&=
\frac12
\int_0^\infty
\frac{e^{at}}{e^t-1}\,dt
-
\frac12
\int_0^\infty
\frac{e^{-at}}{e^t-1}\,dt\\
&=
\frac12
\int_1^\infty
\frac{x^{a-1}}{x-1}\,dx
+
\frac12
\int_0^1
\frac{x^a}{x-1}\,dx\\
&=
\frac12
\int_0^\infty
\frac{x^{a-1}}{x-1}\,dx
+
\frac12
\int_0^1x^a\,dx\\
&=
-\frac\pi2\cot\pi a
+
\frac1{2a}.
\end{aligned}
$$

最后一步用到了第 6 章习题 9(1) 的结果。由此就求出了级数和

$$
\sum_{n=0}^\infty
\frac1{n^2-a^2}
=
-\frac1{a^2}
+
\frac1a
\left(
-\frac\pi2\cot\pi a+\frac1{2a}
\right)
=
-\frac1{2a^2}
-\frac{\pi}{2a}\cot\pi a.
\tag{8.48}
$$

此式还可以改写成

$$
\sum_{n=1}^\infty
\frac1{n^2-a^2}
=
\frac1{2a^2}
-\frac{\pi}{2a}\cot\pi a.
\tag{8.48'}
$$

以上结果是在 $0<\operatorname{Re}a<1$ 的条件下得到的。但是容易延拓到 $\operatorname{Re}a>0$。另外，当 $a$ 为纯虚数时也成立。这时不妨将 (8.48) 和 (8.48') 式中的 $a$ 改写成 $ia$，而设 $a>0$，

$$
\sum_{n=0}^\infty
\frac1{n^2+a^2}
=
\frac1{2a^2}
+
\frac{\pi}{2a}\coth\pi a
\tag{8.49}
$$

或

$$
\sum_{n=1}^\infty
\frac1{n^2+a^2}
=
-\frac1{2a^2}
+
\frac{\pi}{2a}\coth\pi a.
\tag{8.49'}
$$

作为 (8.48') 或 (8.49') 式的特殊情况，也还可以取 $a=0$，所得结果当然与例 7.5 一致。

将 (8.49') 式积分，还可以得到

$$
\sum_{n=1}^\infty
\ln\left(1+\frac{a^2}{n^2}\right)
=
\sum_{n=1}^\infty
\int_0^a
\frac{2x}{n^2+x^2}\,dx
=
\int_0^a
\left(
-\frac1x+\pi\coth\pi x
\right)dx
=
\ln\frac{\sinh\pi a}{\pi a}.
\tag{8.50}
$$

## 习题

（下列各题中的原函数 $f(t)$，均应理解为 $f(t)\eta(t)$）

1. 求下列函数的 Laplace 换式：

   (1) $t^n,\ n=0,1,2,\cdots$；

   (2) $t^\alpha,\ \operatorname{Re}\alpha>-1$；

   (3) $e^{\lambda t}\sin\omega t,\ \lambda>0,\ \omega>0$；

   (4) $\dfrac{\sin\omega t}{t},\ \omega>0$；

   (5) $\dfrac{1-\cos\omega t}{t^2},\ \omega>0$；

   (6)
   $$
   \int_t^\infty\frac{\cos\tau}{\tau}\,d\tau.
   $$

2. 若 $f(t)$ 为周期函数，周期为 $\alpha$，即 $f(t+\alpha)=f(t)$。如果 $f(t)\eta(t)$ 的 Laplace 变换存在，证明：像函数是

$$
F(p)
=
\frac{1}{1-e^{-\alpha p}}
\int_0^\alpha e^{-pt}f(t)\,dt.
$$

3. 求下列函数的 Laplace 换式：

   (1) $|\sin\omega t|,\ \omega>0$；

   (2)
   $$
   t-a\left[\frac ta\right],
   \qquad
   a>0.
   $$

4. 求下列 Laplace 换式的原函数：

   (1)
   $$
   \frac{a^3}{p(p+a)^3};
   $$

   (2)
   $$
   \frac{\omega}{p(p^2+\omega^2)},
   \qquad
   \omega>0;
   $$

   (3)
   $$
   \frac{4p-1}{(p^2+p)(4p^2-1)};
   $$

   (4)
   $$
   \frac{p^2+\omega^2}{(p^2-\omega^2)^2},
   \qquad
   \omega>0;
   $$

   (5)
   $$
   \frac{e^{-p\tau}}{p^2},
   \qquad
   \tau>0;
   $$

   (6)
   $$
   \frac1p\frac{e^{-\alpha p}}{1-e^{-\alpha p}},
   \qquad
   \alpha>0.
   $$

5. 利用 Laplace 变换求解下列微分方程（组）或积分方程：

   (1) 如图 8.8 所示，已知 $i(0)=0,\ q(0)=0$，求 $i(t)$；

![Figure 8.8](...)

图 8.8　习题 8.5 (1)

   (2) 如图 8.9 所示，已知 $i(0)=0,\ q(0)=0$，求 $i(t)$；

![Figure 8.9](...)

图 8.9　习题 8.5 (2)

   (3)

   $$
   y(t)
   =
   a\sin t
   -
   2\int_0^t
   y(\tau)\cos(t-\tau)\,d\tau;
   $$

   (4)

   $$
   f(t)
   +
   2\int_0^t
   f(\tau)\cos(t-\tau)\,d\tau
   =
   9e^{2t}.
   $$

6. 利用 Laplace 变换计算下列积分：

   (1)

   $$
   \int_0^\infty
   \frac{e^{-ax}-e^{-bx}}{x}
   \cos cx\,dx,
   \qquad
   a>0,\ b>0,\ c>0;
   $$

   (2)

   $$
   \int_0^\infty
   \frac{1-\cos bx}{x^2}\,dx,
   \qquad
   b>0;
   $$

   (3)

   $$
   \int_0^\infty
   \frac{\sin xt}{x(x^2+1)}\,dx.
   $$

7. 用普遍反演公式求下列 Laplace 换式的原函数：

   (1)

   $$
   \frac{p}{p^2-\omega^2},
   \qquad
   \omega>0;
   $$

   (2)

   $$
   \frac{e^{-p\tau}}{p^4+4\omega^4},
   \qquad
   \tau>0,\ \omega>0;
   $$

   (3)

   $$
   \frac1p e^{-\alpha p},
   \qquad
   \alpha>0;
   $$

   (4)

   $$
   \frac1p
   \frac{
   \cosh(l-x)\sqrt p
   }{
   \cosh l\sqrt p
   },
   \qquad
   0<x<l.
   $$

8. 求下列无穷级数之和：

   (1)

   $$
   \sum_{n=0}^\infty
   \frac{(-1)^n}{3n+1};
   $$

   (2)

   $$
   \sum_{n=0}^\infty
   \frac{(-1)^n}{4n+1};
   $$

   (3)

   $$
   \sum_{n=0}^\infty
   \frac{(-1)^n}
   {(3n+1)(3n+2)(3n+3)};
   $$

   (4)

   $$
   \sum_{n=0}^\infty
   \frac1
   {(3n+1)(3n+2)(3n+3)}.
   $$
