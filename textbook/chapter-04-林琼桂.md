# 第四章 解析函数的 Laurent 展开与孤立奇点*

本章研究解析函数的 Taylor 展开式的推广，即 Laurent 展开式。它是研究解析函数的奇点的重要工具。

## §1 解析函数的 Laurent 展开

### §1.1 双边幂级数

考虑两个级数

$$
\sum_{n=0}^{+\infty} c_n (z-a)^n \tag{1}
$$

和

$$
\sum_{n=1}^{+\infty} \frac{c_{-n}}{(z-a)^n}. \tag{2}
$$

级数 (1) 即是上章研究的幂级数，它和级数 (2) 之和

$$
\sum_{n=-\infty}^{+\infty} c_n (z-a)^n \tag{3}
$$

称为双边幂级数。

假设级数 (1) 的收敛半径为 $R$（$0<R\le +\infty$），则它在圆 $|z-a|<R$ 上绝对收敛且内闭一致收敛，并具有解析的和函数，记作 $f_1(z)$。对于级数 (2)，令 $\zeta=1/(z-a)$，则可改写为

$$
\sum_{n=1}^{+\infty} c_{-n}\zeta^n,
$$

假设它的收敛半径为 $1/r$（$0<1/r\le +\infty$），则它在圆 $|\zeta|<1/r$ 上绝对收敛且内闭一致收敛，并具有解析的和函数。换句话说，级数 (2) 当 $|z-a|>r$（$0\le r<+\infty$）时绝对收敛且内闭一致收敛，并具有解析的和函数，记作 $f_2(z)$。

若 $r>R$，则级数 (1) 和 (2) 没有公共收敛区域，因而双边幂级数 (3) 处处发散。

**例 1** 双边幂级数

$$
\sum_{n=0}^{+\infty}(2z)^n+\sum_{n=-\infty}^{-1}z^n
=\sum_{n=0}^{+\infty}(2z)^n+\sum_{n=1}^{+\infty}\frac1{z^n}
$$

的正幂部分在圆 $|z|<1/2$ 内绝对收敛，负幂部分在单位圆外 $|z|>1$ 绝对收敛。所以原双边幂级数处处发散。

若 $r=R$，则级数 (1) 和 (2) 亦没有公共收敛区域。但在圆周 $|z-a|=R=r$ 上可能存在收敛点。

**例 2** 双边幂级数

$$
\sum_{n=1}^{+\infty}\frac{z^n}{n^2}+\sum_{n=-\infty}^{-1}\frac{z^n}{n^2}
=\sum_{n=1}^{+\infty}\frac{z^n}{n^2}+\sum_{n=1}^{+\infty}\frac{1}{n^2z^n}
$$

的正幂部分在单位圆内 $|z|<1$ 绝对收敛，负幂部分在单位圆外 $|z|>1$ 绝对收敛。所以原双边幂级数没有公共收敛区域，但它在单位圆周 $|z|=1$ 上绝对收敛。

若 $r<R$，级数 (1) 和 (2) 有公共的收敛区域，即环域

$$
H:\qquad r<|z-a|<R\qquad (0\le r<R\le +\infty).
$$

这时，根据上章的 Weierstrass 定理，有如下

**定理** 双边幂级数 (3) 具有下列性质：

(1) 在收敛环 $H$ 内绝对收敛于 $f(z)=f_1(z)+f_2(z)$，且在闭环：$r<\rho_1\le |z-a|\le \rho_2<R$ 上一致收敛（即在 $H$ 上内闭一致收敛）；

(2) 和函数 $f(z)$ 在 $H$ 内解析；

(3) $f(z)=\sum_{n=-\infty}^{+\infty}c_n(z-a)^n$ 在 $H$ 内可以逐项求导和逐项积分。

**例 3** 双边幂级数

$$
\sum_{n=0}^{+\infty}(z/2)^n+\sum_{n=-\infty}^{-1}z^n
=\sum_{n=0}^{+\infty}(z/2)^n+\sum_{n=1}^{+\infty}\frac1{z^n}
$$

的正幂部分在圆 $|z|<2$ 内绝对收敛于函数 $2/(2-z)$，负幂部分在单位圆外 $|z|>1$ 绝对收敛于函数 $1/(z-1)$。所以原双边幂级数在环域 $H:1<|z|<2$ 中绝对收敛于函数 $-z/(z-1)(z-2)$。容易看出，该和函数在 $H$ 中解析。

### §1.2 解析函数的 Laurent 展开

由上面的分析，双边幂级数在收敛环内具有解析的和函数，换句话说，它在收敛环内代表一个解析函数。反过来，在环域内解析的函数是否可以展开为双边幂级数呢？下面的定理给出了肯定的答案。

**定理（Laurent）** 设函数 $f(z)$ 在环域

$$
H:r<|z-a|<R\qquad (0\le r<R\le +\infty)
$$

内解析，则在 $H$ 内可以展开为双边幂级数：

$$
f(z)=\sum_{n=-\infty}^{+\infty}c_n(z-a)^n, \tag{4}
$$

其中

$$
c_n=\frac{1}{2\pi i}\int_\Gamma \frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta,\qquad n\in\mathbb Z, \tag{5}
$$

而 $\Gamma$ 是环内包围内圆的任一围线。且展开式是唯一的。

式 (4) 称为 $f(z)$ 在 $a$ 点的 Laurent 展开式，右边称为 Laurent 级数，式 (5) 称为 Laurent 系数。

**注** ① Laurent 定理在形式上与 Taylor 定理非常相似，证明的方法也相似（参见以下详细的证明）。了解两者的联系和区别对于深入理解这一定理是重要的。② 一般来说，即使正幂项的系数也不能表示为高阶导数的形式。这是因为 $f(z)$ 在闭圆 $|z-a|\le r$ 上有奇点，所以在 $\Gamma$ 上的积分不满足应用 Cauchy 高阶导数公式的条件。③ 如果 $f(z)$ 在闭圆 $|z-a|\le r$ 上确实没有奇点，那么它就在圆 $|z-a|<R$ 上解析，这时 Laurent 级数应该退化为 Taylor 级数。事实上，这时有

$$
c_{-n}=\frac{1}{2\pi i}\int_\Gamma \frac{f(\zeta)}{(\zeta-a)^{-n+1}}\,d\zeta
=\frac{1}{2\pi i}\int_\Gamma f(\zeta)(\zeta-a)^{n-1}\,d\zeta=0,
\qquad n\in\mathbb N^+,
$$

其中最后一步用了 Cauchy 积分定理，因为被积函数中的两个因子均在 $\Gamma$ 所包围的闭域上解析。可见 Taylor 级数是 Laurent 级数的特殊情况。④ 一般情况下，展开式中有负幂项，因为 $f(z)$ 在闭圆 $|z-a|\le r$ 上有奇点。但是，$f(z)$ 的奇点不一定在 $a$。所以，不要因为展开式中有 $z-a$ 的负幂项就误以为 $a$ 是 $f(z)$ 的奇点。例如，当 $1<|z|<+\infty$ 时，有

$$
\frac1{z-1}
=\frac1z\frac1{1-1/z}
=\frac1z\sum_{n=0}^{+\infty}\frac1{z^n}
=\sum_{n=1}^{+\infty}\frac1{z^n}.
$$

展开式中全是 $z$ 的负幂项，但 $z=0$ 并不是被展开函数的奇点。这是因为展开式在 $z=0$ 附近并不成立。⑤ 知道了展开式的唯一性，我们就可以用任何方法来求展开式，而不一定要用式 (5) 来计算系数。最常用的方法是利用已知的 Taylor 级数展开式。具体方法参见下面的展开实例。

下面给出 Laurent 定理的

**证明** $\forall z\in H$，在 $H$ 内作圆 $\Gamma_1:|\zeta-a|=\rho_1$ 和 $\Gamma_2:|\zeta-a|=\rho_2$，使得 $r<\rho_1<|z-a|<\rho_2<R$。由于 $f(\zeta)$ 在闭环 $\rho_1\le |\zeta-a|\le\rho_2$ 上解析，故由 Cauchy 积分公式，有

$$
f(z)=\frac1{2\pi i}\int_{\Gamma_2}\frac{f(\zeta)}{\zeta-z}\,d\zeta
+\frac1{2\pi i}\int_{\Gamma_1}\frac{f(\zeta)}{z-\zeta}\,d\zeta,
$$

其中积分皆沿曲线的正方向进行。对于第一项，仿照 Taylor 定理的证明，易得

$$
\frac1{2\pi i}\int_{\Gamma_2}\frac{f(\zeta)}{\zeta-z}\,d\zeta
=\sum_{n=0}^{+\infty}c_n(z-a)^n,
$$

其中

$$
c_n=\frac1{2\pi i}\int_{\Gamma_2}\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta,
\qquad n\in\mathbb N.
$$

对于第二项，由于 $|\zeta-a|<|z-a|$，故有

$$
\begin{aligned}
\frac{f(\zeta)}{z-\zeta}
&=\frac{f(\zeta)}{(z-a)-(\zeta-a)}
=\frac{f(\zeta)}{z-a}\frac1{1-(\zeta-a)/(z-a)}\\
&=\frac{f(\zeta)}{z-a}\sum_{n=1}^{\infty}\left(\frac{\zeta-a}{z-a}\right)^{n-1}\\
&=\sum_{n=1}^{\infty}\frac{f(\zeta)}{(\zeta-a)^{-n+1}}(z-a)^{-n}
=\sum_{n=-\infty}^{-1}\frac{f(\zeta)}{(\zeta-a)^{n+1}}(z-a)^n.
\end{aligned}
$$

上式右边在 $\Gamma_1$ 上一致收敛，故可逐项积分，而得

$$
\frac1{2\pi i}\int_{\Gamma_1}\frac{f(\zeta)}{z-\zeta}\,d\zeta
=\sum_{n=-\infty}^{-1}c_n(z-a)^n,
$$

其中

$$
c_n=\frac1{2\pi i}\int_{\Gamma_1}\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta,
\qquad n=-1,-2,\cdots.
$$

设 $\Gamma$ 为 $H$ 中包围内圆的任一围线，由 Cauchy 积分定理，有

$$
\int_{\Gamma_2}\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta
=\int_{\Gamma_1}\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta
=\int_\Gamma\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta,
\qquad n\in\mathbb Z.
$$

综合以上结果，即得式 (4) 和式 (5)。

下面证明展开式的唯一性。设有另一展开式

$$
f(z)=\sum_{n=-\infty}^{+\infty}c_n'(z-a)^n,
$$

两边乘以 $(z-a)^{-m-1}$，得

$$
\frac{f(z)}{(z-a)^{m+1}}=\sum_{n=-\infty}^{+\infty}c_n'(z-a)^{n-m-1}.
$$

在 $\Gamma$ 上积分，由于右边在 $\Gamma$ 上一致收敛，故可逐项积分，利用第二章例题结果，右边积分只有一项不为零，即

$$
\int_\Gamma\frac{f(z)}{(z-a)^{m+1}}\,dz
=\sum_{n=-\infty}^{+\infty}c_n'\int_\Gamma(z-a)^{n-m-1}\,dz
=2\pi i c_m',
$$

所以

$$
c_m'=\frac1{2\pi i}\int_\Gamma\frac{f(\zeta)}{(\zeta-a)^{m+1}}\,d\zeta=c_m,
\qquad m\in\mathbb Z.
$$

这就证明了唯一性。证毕。

### §1.3 展开实例

**例 4** 在 $a=0$ 处展开 $f(z)=\dfrac1{z(z-1)}$ 为 Laurent 级数。

**解** $f(z)$ 有奇点 $z=0$ 和 $z=1$，故分别在 $H_1:0<|z|<1$ 和 $H_2:1<|z|<+\infty$ 上解析。在 $H_1$ 上，

$$
f(z)=-\frac1z\frac1{1-z}
=-\frac1z\sum_{n=0}^{+\infty}z^n
=-\sum_{n=0}^{+\infty}z^{n-1}
=-\sum_{n=-1}^{+\infty}z^n.
$$

在 $H_2$ 上，

$$
f(z)=\frac1{z^2}\frac1{1-1/z}
=\frac1{z^2}\sum_{n=0}^{+\infty}\frac1{z^n}
=\sum_{n=2}^{+\infty}\frac1{z^n}.
$$

**例 5** 在 $0<|z|<+\infty$ 中可展开

$$
\cosh\left(z+\frac1z\right)=\sum_{n=-\infty}^{+\infty}c_nz^n,
$$

试证 $c_{-n}=c_n$，$c_{2n-1}=0$（$n\in\mathbb N^+$）。

**解** 在原式中令 $z=1/\zeta$，可得

$$
\cosh\left(\zeta+\frac1\zeta\right)
=\sum_{n=-\infty}^{+\infty}c_n\zeta^{-n}
=\sum_{n=-\infty}^{+\infty}c_{-n}\zeta^n,
\qquad 0<|\zeta|<+\infty.
$$

把自变量重新写作 $z$，则有

$$
\cosh\left(z+\frac1z\right)
=\sum_{n=-\infty}^{+\infty}c_{-n}z^n,
\qquad 0<|z|<+\infty.
$$

与原式比较，利用展开式的唯一性，即得

$$
c_{-n}=c_n,\qquad n\in\mathbb N^+.
$$

另一方面，在原式中令 $z=-\zeta$，可得

$$
\cosh\left(\zeta+\frac1\zeta\right)
=\sum_{n=-\infty}^{+\infty}(-)^n c_n\zeta^n,
\qquad 0<|\zeta|<+\infty.
$$

把自变量重新写作 $z$，则有

$$
\cosh\left(z+\frac1z\right)
=\sum_{n=-\infty}^{+\infty}(-)^n c_n z^n,
\qquad 0<|z|<+\infty.
$$

与原式比较，利用展开式的唯一性，即得

$$
c_n=(-)^n c_n,\qquad n\in\mathbb Z.
$$

考虑到 $c_{-n}=c_n$，独立的结论有

$$
c_{2n-1}=0,\qquad n\in\mathbb N^+.
$$

最后，展开式可以写成

$$
\cosh\left(z+\frac1z\right)
=c_0+\sum_{n=1}^{+\infty}c_{2n}\left(z^{2n}+\frac1{z^{2n}}\right),
\qquad 0<|z|<+\infty,
$$

其中 $c_{2n}$ 可用积分表出，这里从略。以上结论也可以用系数的积分公式证明，但比较麻烦。

**例 6** 求函数 $f(t)=e^{(z/2)(t-1/t)}$ 在 $0<|t|<+\infty$ 中的 Laurent 展开式，其中 $z$ 是参数。

**解** 用指数函数的 Taylor 展开式将 $e^{(z/2)t}$ 和 $e^{-(z/2)(1/t)}$ 分别展开为其宗量的 Taylor 级数。因为两个级数均在 $0<|t|<+\infty$ 中绝对收敛，所以乘积级数也绝对收敛，从而可以重新排列各项，得到

$$
\begin{aligned}
f(t)
&=\sum_{k,l=0}^{\infty}\frac{(-)^l}{k!l!}\left(\frac z2\right)^{k+l}t^{k-l}\\
&=\sum_{l=0}^{\infty}\sum_{k=l}^{\infty}\frac{(-)^l}{k!l!}\left(\frac z2\right)^{k+l}t^{k-l}
+\sum_{k=0}^{\infty}\sum_{l=k+1}^{\infty}\frac{(-)^l}{k!l!}\left(\frac z2\right)^{k+l}t^{k-l},
\end{aligned}
$$

其中第一和式是 $k\ge l$ 部分，第二和式是 $k<l$ 部分。在第一和式中令 $k=l+m$，在第二和式中令 $l=k+m$，可以化为

$$
f(t)=\sum_{m=0}^{\infty}\left[\sum_{l=1}^{\infty}\frac{(-)^l}{l!(l+m)!}\left(\frac z2\right)^{2l+m}\right]t^m
+\sum_{m=1}^{\infty}(-)^m\left[\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}\left(\frac z2\right)^{2k+m}\right]t^{-m}.
$$

两个方括号内的求和实际上是一样的。定义 $m$ 阶 Bessel 函数为

$$
J_m(z)=\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}\left(\frac z2\right)^{2k+m},
$$

则上面的展开式可以写成

$$
f(t)=\sum_{m=0}^{\infty}J_m(z)t^m+\sum_{m=1}^{\infty}(-)^mJ_m(z)t^{-m}.
$$

在定义了更一般的 $\nu$（$\nu\in\mathbb C$）阶 Bessel 函数 $J_\nu(z)$（参看第十章）之后，可以证明 $J_{-m}(z)=(-)^mJ_m(z)$，于是上式也可以写成

$$
f(t)=\sum_{m=-\infty}^{\infty}J_m(z)t^m.
$$

可见当函数 $f(t)=e^{(z/2)(t-1/t)}$ 展开为 $t$ 的 Laurent 级数时，系数正好是参数 $z$ 的整数阶 Bessel 函数，所以该函数称为整数阶 Bessel 函数的母函数。从上述展开式出发，可以导出整数阶 Bessel 函数的一些有用的性质。

**习题** 将下列各函数在指定的环域内展开为 Laurent 级数。(1) $1/z^2(z-1)$，$0<|z-1|<1$；(2) $(z-1)(z-2)/(z-3)(z-4)$，$4<|z|<\infty$；(3) $\cot z$，$0<|z|<\pi$（计算三个非零项）。

## §2 解析函数的零点与孤立奇点

### §2.1 解析函数的零点

为了后面讨论的需要，这里简单介绍一下解析函数的零点的概念。

**定义（$m$ 阶零点）** 若函数 $f(z)$ 在 $a$ 点解析，且

$$
f(a)=f'(a)=\cdots=f^{(m-1)}(a)=0,
$$

而 $f^{(m)}(a)\ne0$，则 $a$ 称为 $f(z)$ 的 $m$ 阶零点。一阶零点亦称为单零点。

**注** $f(z)$ 在 $a$ 点解析，即在某圆 $K:|z-a|<R$ 内解析，在该圆内，$f(z)$ 可展开为 Taylor 级数 $f(z)=\sum_{n=0}^{+\infty}c_n(z-a)^n$，其中 $c_n=f^{(n)}(a)/n!$。如果所有的系数均为零，则 $f(z)$ 在 $K$ 内恒为零。若 $f(z)$ 在 $K$ 内不恒为零，则上述定义中的 $m$ 值总是存在的。

**例 1** $f(z)=\sin z$ 的零点为 $z_n=n\pi$（$n\in\mathbb Z$）。由于 $f'(z_n)=\cos z_n=\cos n\pi=(-)^n\ne0$，故所有的 $z_n$ 均为单零点。

**例 2** 显然 $z=0$ 是函数 $f(z)=z-\sin z$ 的零点。由于 $f'(0)=(1-\cos z)|_{z=0}=0$，$f''(0)=\sin z|_{z=0}=0$，而 $f'''(0)=\cos z|_{z=0}=1\ne0$，故 $z=0$ 是 $f(z)$ 三阶零点。

**例 3** $f(z)=(z-a)^m$，$m\in\mathbb N^+$。显然 $a$ 是 $f(z)$ 的 $m$ 阶零点。

**例 4** $f(z)=(z-a)^m\varphi(z)$，其中 $\varphi(z)$ 在 $a$ 点解析且 $\varphi(a)\ne0$。显然 $a$ 是 $f(z)$ 的 $m$ 阶零点。

多项式是最简单的解析函数。如果 $a$ 是多项式 $P_n(z)$ 的 $m$ 重根（$m\le n$），则有 $P_n(z)=(z-a)^mQ_{n-m}(z)$，其中 $Q_{n-m}(z)$ 是 $n-m$ 次多项式，且 $Q_{n-m}(a)\ne0$。容易验证，$P_n(z)$ 在 $a$ 点满足以上条件，所以 $a$ 是 $P_n(z)$ 的 $m$ 阶零点。一般解析函数的 $m$ 阶零点是多项式 $m$ 重根的推广。所以，类似于多项式的情况，有以下

**定理（$m$ 阶零点）** 若函数 $f(z)$ 在 $a$ 点的邻域 $K:|z-a|<R$ 内解析且不恒为零，则 $f(z)$ 以 $a$ 为 $m$ 阶零点的充要条件是 $f(z)=(z-a)^m\varphi(z)$，其中 $\varphi(z)$ 在 $K$ 内解析，且 $\varphi(a)\ne0$。

**证明** 充分性。若 $f(z)=(z-a)^m\varphi(z)$，且 $\varphi(a)\ne0$，则容易验证 $f(a)=f'(a)=\cdots=f^{(m-1)}(a)=0$，而 $f^{(m)}(a)=m!\varphi(a)\ne0$，故 $a$ 是 $f(z)$ 的 $m$ 阶零点。

必要性。$f(z)$ 在 $K$ 内解析，故可展开为 Taylor 级数 $f(z)=\sum_{n=0}^{+\infty}c_n(z-a)^n$，其中 $c_n=f^{(n)}(a)/n!$。由于 $a$ 是 $m$ 阶零点，故 $c_0=c_1=\cdots=c_{m-1}=0$，而 $c_m\ne0$，所以

$$
\begin{aligned}
f(z)
&=\sum_{n=m}^{+\infty}c_n(z-a)^n
=(z-a)^m\sum_{n=m}^{+\infty}c_n(z-a)^{n-m}\\
&=(z-a)^m\sum_{k=0}^{+\infty}c_{k+m}(z-a)^k
\equiv (z-a)^m\varphi(z).
\end{aligned}
$$

显然 $\varphi(z)=\sum_{k=0}^{+\infty}c_{k+m}(z-a)^k$ 在 $K$ 内解析（注意 $\varphi(z)$ 与 $f(z)$ 只相差一个因子 $(z-a)^m$，且没有奇性），且 $\varphi(a)=c_m\ne0$。证毕。

由上面的定理容易推出零点的孤立性，即有下面的

**定理（零点的孤立性）** 若函数 $f(z)$ 在 $a$ 点的邻域 $K:|z-a|<R$ 内解析且不恒为零，$a$ 是 $f(z)$ 的零点，则必存在邻域 $K':|z-a|<r$，在其中 $f(z)$ 只有 $a$ 一个零点。

**证明** 由于 $f(z)$ 在 $K$ 内解析，故可展开为 Taylor 级数 $f(z)=\sum_{n=0}^{+\infty}c_n(z-a)^n$，其中 $c_n=f^{(n)}(a)/n!$。如果所有的系数均为零，则 $f(z)$ 在 $K$ 内恒为零，这与假设矛盾。故必存在 $m\in\mathbb N^+$，使得 $a$ 满足 $m$ 阶零点的条件。按以上定理，必有 $f(z)=(z-a)^m\varphi(z)$，其中 $\varphi(z)$ 在 $K$ 内解析，且 $\varphi(a)\ne0$。由于 $\varphi(z)$ 连续，故存在邻域 $K':|z-a|<r$，在其中 $\varphi(z)\ne0$，从而 $f(z)$ 在其中只有 $a$ 一个零点。证毕。

这一定理有一个重要的

**推论** 若函数 $f(z)$ 在区域 $D$ 内解析，在 $D$ 内有点列 $\{z_n\}_{n=1}^{+\infty}$，满足 $z_n\ne a$，$z_n\to a\in D$，$f(z_n)=0$，则在 $D$ 内 $f(z)\equiv0$。

**证明** 由于 $f(z)$ 连续，易得 $f(a)=\lim_{n\to+\infty}f(z_n)=0$，故 $a$ 是 $f(z)$ 的非孤立零点。由上述定理，在 $a$ 点的某邻域 $K:|z-a|<R$ 内必有 $f(z)\equiv0$，只要 $K\subset D$。如果 $D$ 本身就是 $a$ 点的邻域，或者 $D$ 是复平面，则结论得证。否则，设 $b\ne a$ 是 $D$ 内任一点，只需证明 $f(b)=0$ 即可。今在 $D$ 内作折线 $L$ 连接 $a$ 和 $b$，设 $L$ 与 $C=\partial D$ 的最小距离为 $d$（$d>0$）。在 $L$ 上由 $a$ 至 $b$ 取分点 $a_1,a_2,\cdots,a_m$，使

$$
|a_1-a|=|a_2-a_1|=\cdots=|a_m-a_{m-1}|=R/2,
$$

而 $|b-a_m|\le R/2$，其中 $0<R<d$。考虑一系列圆 $K:|z-a|<R$ 和 $K_i:|z-a_i|<R$（$i=1,2,\cdots,m$），它们都在 $D$ 内。如上所述，在 $K$ 内必有 $f(z)\equiv0$。按分点的取法，$a_1$ 落在 $K$ 内，且是 $f(z)$ 的非孤立零点，故可推知在 $K_1$ 内必有 $f(z)\equiv0$，继续同样的推理，最后可得在 $K_m$ 内必有 $f(z)\equiv0$，而 $b$ 落在 $K_m$ 内，故 $f(b)=0$。证毕。

由这一推论，马上可以推出非常重要的唯一性定理：

**定理（解析函数的唯一性）** 若函数 $f_1(z)$ 和 $f_2(z)$ 均在区域 $D$ 内解析，在 $D$ 内有点列 $\{z_n\}_{n=1}^{+\infty}$，满足 $z_n\ne a$，$z_n\to a\in D$，$f_1(z_n)=f_2(z_n)$，则在 $D$ 内必有 $f_1(z)=f_2(z)$。

**证明** 令 $F(z)=f_1(z)-f_2(z)$，则 $F(z)$ 在区域 $D$ 内解析，且在点列 $\{z_n\}_{n=1}^{+\infty}$ 上有 $F(z_n)=0$，根据以上推论，在 $D$ 内必有 $F(z)\equiv0$，从而 $f_1(z)=f_2(z)$。证毕。

唯一性定理表明，解析函数在区域 $D$ 内一个点列上的取值，就完全决定了它在区域 $D$ 内的形式。这里我们再次看到，解析性对于函数值的分布给予了极大的限制。如果将定理中的点列改为区域 $D$ 内的一个弧段或子区域，定理当然还是成立的。下面是应用唯一性定理的例子。

**例 5** 函数 $\sin2z$ 和 $2\sin z\cos z$ 均在复平面上解析，而在实轴上有 $\sin2x=2\sin x\cos x$，由唯一性定理，在复平面上必有 $\sin2z=2\sin z\cos z$。按照这种方式，可以将实变函数中的许多恒等式推广到复变函数中。

**例 6** 上章介绍了解析开拓的基本概念，并且指出，如果区域 $D$ 内的解析函数 $f(z)$ 可以解析开拓到 $G$（$G\supset D$），则结果是唯一的。现在用唯一性定理容易证明这一结论。事实上，如果 $F_1(z)$ 和 $F_2(z)$ 都是 $f(z)$ 在 $G$ 内的解析开拓，由于在 $G$ 的子区域 $D$ 上必有 $F_1(z)=F_2(z)=f(z)$，根据唯一性定理，在 $G$ 内必有 $F_1(z)=F_2(z)$。

### §2.2 解析函数的孤立奇点及其分类

**定义（孤立奇点）** 如果函数 $f(z)$ 以 $a$ 为奇点，但在 $a$ 的某去心邻域 $K\backslash\{a\}:0<|z-a|<R$ 中解析，则 $a$ 称为 $f(z)$ 的孤立奇点（isolated singularity）。

粗略地说，如果函数 $f(z)$ 以 $a$ 为奇点，但在 $a$ 附近没有别的奇点，则 $a$ 就是 $f(z)$ 的孤立奇点。所以这一定义是非常直观的。

去心邻域 $K\backslash\{a\}:0<|z-a|<R$ 是内半径为零的环域，$f(z)$ 在其中解析，则可以展开为 Laurent 级数：

$$
f(z)=\sum_{n=-\infty}^{+\infty}c_n(z-a)^n.
$$

展开式中的正幂部分，即 $\sum_{n=0}^{+\infty}c_n(z-a)^n$ 称为正则部分，而负幂部分，即 $\sum_{n=1}^{+\infty}c_{-n}/(z-a)^n$ 称为主要部分。由于现在展开式在 $a$ 点附近成立，故其中负幂项的多少可以刻画 $f(z)$ 在 $a$ 点的奇性的大小。Laurent 展开式的主要用途正在于此。根据主要部分的不同情况，我们有如下

**定义（孤立奇点的分类）** (1) 如果主要部分为零，则 $a$ 称为 $f(z)$ 的可去奇点（removable singularity）。(2) 如果主要部分为 $\sum_{k=1}^{m}c_{-k}/(z-a)^k$，其中 $c_{-m}\ne0$，则 $a$ 称为 $f(z)$ 的 $m$ 阶极点（pole）。$m=1$ 时亦称为单极点。(3) 如果主要部分有无穷多项，则 $a$ 称为 $f(z)$ 的本性奇点（essential singularity）。

## §3 各种孤立奇点的判断

### §3.1 可去奇点

先看看下面的例子。

**例 1** 函数 $\sin z/z$ 在 $0<|z|<+\infty$ 内解析，故可以展开为 Laurent 级数：

$$
\frac{\sin z}{z}
=\frac1z\sum_{n=0}^{+\infty}\frac{(-)^n z^{2n+1}}{(2n+1)!}
=\sum_{n=0}^{+\infty}\frac{(-)^n z^{2n}}{(2n+1)!}. \tag{6}
$$

由于展开式的主要部分为零，按定义，$z=0$ 是可去奇点。如果只看右边的幂级数，易知其收敛半径为 $+\infty$，故它在复平面上解析，包括 $z=0$ 点。由于函数

$$
f(z)=
\begin{cases}
\sin z/z, & z\ne0,\\
1, & z=0
\end{cases}
$$

在整个复平面上与式 (6) 右边的幂级数相等，故 $f(z)$ 在复平面上解析，包括 $z=0$ 点。这就是说，只要对原来的函数 $\sin z/z$ 适当补上 $z=0$ 处的定义，即可构造出一个解析函数。所以 $z=0$ 称为可去奇点是非常恰当的。一般情况下，只要令 $f(a)=c_0$ 即可将可去奇点 $a$ 变为解析点。

下面的定理用于判断可去奇点。

**定理（可去奇点）** 函数 $f(z)$ 的孤立奇点 $a$ 为可去奇点的充要条件是

$$
\lim_{z\to a}f(z)=b\ne\infty.
$$

**证明** 必要性。按定义，此时 Laurent 展开式没有负幂项，即

$$
f(z)=\sum_{n=0}^{+\infty}c_n(z-a)^n,
\qquad 0<|z-a|<R.
$$

显然，$\lim_{z\to a}f(z)=c_0\ne\infty$。

充分性。此时 $f(z)$ 在 $a$ 的某去心邻域 $K:0<|z-a|<R_0$ 内有界：$|f(z)|<M$。（事实上，按极限定义，取 $\varepsilon=1$，$\exists R_0$，当 $0<|z-a|<R_0$，有 $|f(z)-b|<1$，则 $|f(z)|<|b|+1$。）按公式 (5)，Laurent 展开式的负幂项系数

$$
c_{-n}=\frac1{2\pi i}\int_{|\zeta-a|=\rho}f(\zeta)(\zeta-a)^{n-1}\,d\zeta,
\qquad n\in\mathbb N^+,
$$

其中 $\rho>0$ 可以任取，但不超过展开环域的外半径。$c_{-n}$ 实际上与 $\rho$ 的大小无关，取 $0<\rho<R_0$，则易得

$$
|c_{-n}|=
\left|\frac1{2\pi i}\int_{|\zeta-a|=\rho}f(\zeta)(\zeta-a)^{n-1}\,d\zeta\right|
<\frac{M}{2\pi}\rho^{n-1}\int_{|\zeta-a|=\rho}|d\zeta|
=M\rho^n,
\qquad n\in\mathbb N^+.
$$

最后，令 $\rho\to0$，则得 $c_{-n}=0$（$n\in\mathbb N^+$）。所以，展开式的主要部分为零，即 $a$ 是可去奇点。证毕。

下面是关于可去奇点的一个更一般的

**定理** 若 $a$ 是 $f(z)$ 的孤立奇点，则下列三命题等价：

(1) $f(z)$ 在 $a$ 的主要部分为零，即 $a$ 为可去奇点；

(2) $\lim_{z\to a}f(z)=b\ne\infty$；

(3) $f(z)$ 在 $a$ 的某去心邻域内有界。

这类定理的证明方式是 $(1)\Rightarrow(2)\Rightarrow(3)\Rightarrow(1)$。本定理的具体证明基本上已包含在上一定理的证明中，故从略。

### §3.2 极点

先看一个简单的例子。

**例 2** 函数 $e^z/z$ 在 $0<|z|<+\infty$ 内解析，故可以展开为 Laurent 级数：

$$
\frac{e^z}{z}=\frac1z+\sum_{n=0}^{+\infty}\frac{z^n}{(n+1)!}. \tag{7}
$$

由于展开式的主要部分有一项，按定义，$z=0$ 是单极点。

判断极点，主要用下面的

**定理（$m$ 阶极点）** 若 $a$ 是 $f(z)$ 的孤立奇点，则下列三命题等价：

(1) $f(z)$ 在 $a$ 的主要部分为 $\sum_{k=1}^{m}c_{-k}/(z-a)^k$（其中 $c_{-m}\ne0$），即 $a$ 为 $m$ 阶极点；

(2) $f(z)$ 在 $a$ 的某去心邻域内可表为 $f(z)=\varphi(z)/(z-a)^m$，其中 $\varphi(z)$ 在 $a$ 解析，且 $\varphi(a)\ne0$；

(3) $g(z)=1/f(z)$ 以 $a$ 为 $m$ 阶零点。

**注** 其中第一条是 $m$ 阶极点的定义，其它两条也很直观。下面的证明供有兴趣的读者参考。

**证明** $(1)\Rightarrow(2)$。在 $a$ 的某去心邻域内，有

$$
f(z)=\sum_{n=-m}^{+\infty}c_n(z-a)^n
=(z-a)^{-m}\sum_{n=-m}^{+\infty}c_n(z-a)^{n+m}
=(z-a)^{-m}\sum_{k=0}^{+\infty}c_{k-m}(z-a)^k
\equiv\frac{\varphi(z)}{(z-a)^m}.
$$

其中 $\varphi(z)=\sum_{k=0}^{+\infty}c_{k-m}(z-a)^k$ 是 Taylor 级数，且收敛半径不为零（否则 $\varphi(z)$ 除 $z=a$ 外处处发散，由上式，$f(z)$ 在 $a$ 的任一去心邻域内发散，与假设矛盾），故在 $a$ 点解析，且 $\varphi(a)=c_{-m}\ne0$。

$(2)\Rightarrow(3)$。由已知条件，得 $g(z)=(z-a)^m/\varphi(z)$。因为 $\varphi(z)$ 在 $a$ 解析，且 $\varphi(a)\ne0$，故在 $a$ 的某邻域内 $\varphi(z)\ne0$，在该邻域内 $1/\varphi(z)$ 亦解析，且显然 $1/\varphi(a)\ne0$。根据 §2 的定理，$a$ 是 $g(z)$ 的 $m$ 阶零点。

$(3)\Rightarrow(1)$。根据 §2 的定理，有 $g(z)=(z-a)^m\psi(z)$。其中 $\psi(z)$ 在 $a$ 解析且 $\psi(a)\ne0$，故在 $a$ 的某邻域内 $\psi(z)\ne0$，在该邻域内 $1/\psi(z)$ 亦解析，故可展开为 Taylor 级数：$1/\psi(z)=\sum_{n=0}^{+\infty}c_{n-m}(z-a)^n$，其中系数的写法是为了下面的方便。于是

$$
f(z)=\frac{(z-a)^{-m}}{\psi(z)}
=\sum_{n=0}^{+\infty}c_{n-m}(z-a)^{n-m}
=\sum_{n=-m}^{+\infty}c_n(z-a)^n,
$$

其中 $c_{-m}=1/\psi(a)\ne0$。所以 $a$ 是 $f(z)$ 的 $m$ 阶极点。证毕。

下面就用以上定理来分析几个实例。

**例 3** 再回头看例 2 的函数 $e^z/z$，由于 $e^z$ 解析，且 $e^0=1\ne0$，根据定理的第二条，$z=0$ 是该函数的单极点。

**例 4** 考虑函数 $\cot z=\cos z/\sin z$，易知 $z=n\pi$（$n\in\mathbb Z$）都是其孤立奇点，而且它们都是 $1/\cot z=\sin z/\cos z$ 的单零点，根据定理的第三条，$z=n\pi$（$n\in\mathbb Z$）都是 $\cot z$ 的单极点。

**例 5** 考虑函数 $f(z)=1/(z-\sin z)$，由 §2 例 2 知 $z=0$ 是 $z-\sin z$ 的三阶零点，根据定理的第三条，$z=0$ 是 $f(z)$ 三阶极点。

由上述定理，容易得到以下

**推论（极点的判定）** 函数 $f(z)$ 的孤立奇点 $a$ 为极点的充要条件是

$$
\lim_{z\to a}f(z)=\infty.
$$

这一推论主要用来判定一个孤立奇点是否极点，若要进一步判定其阶数，就需要用上面的定理。

### §3.3 本性奇点

由前面关于可去奇点和极点的讨论，容易得到以下判定本性奇点的

**定理（本性奇点）** 函数 $f(z)$ 的孤立奇点 $a$ 为本性奇点的充要条件是 $\lim_{z\to a}f(z)$ 不存在，即不为有限值，亦不为 $\infty$。

这说明在本性奇点处极限 $\lim_{z\to a}f(z)$ 与 $z\to a$ 的方式有关。

**例 6** 函数 $e^{1/z}$ 在 $0<|z|<+\infty$ 内解析，故可以展开为 Laurent 级数：

$$
e^{1/z}=\sum_{n=0}^{+\infty}\frac1{n!}\frac1{z^n}.
$$

展开式的主要部分有无穷多项，按定义，$z=0$ 是本性奇点。当 $z=x$，容易知道

$$
\lim_{x\to0^+}e^{1/x}=+\infty,\qquad
\lim_{x\to0^-}e^{1/x}=0,
$$

可见极限 $\lim_{z\to0}e^{1/z}$ 不存在。

实际上，关于本性奇点，有下述

**定理（Weierstrass）** 若 $a$ 为 $f(z)$ 的本性奇点，则对任意复数 $A$（有限或无穷），都有点列 $\{z_n\}_{n=1}^{\infty}$，$\lim_{n\to\infty}z_n=a$，使 $\lim_{n\to\infty}f(z_n)=A$。

**证明** 首先，如果 $A=\infty$，必能找到所要求的点列。否则，必存在 $M>0$ 和 $\delta>0$，当 $0<|z-a|<\delta$，有 $|f(z)|<M$，如此，则 $a$ 为可去奇点，与假设矛盾。

其次，对于任何有限的 $A$，亦必能找到所要求的点列。否则，必有 $\varepsilon>0$ 和 $\delta>0$，当 $0<|z-a|<\delta$，有 $|f(z)-A|\ge\varepsilon$，于是 $1/[f(z)-A]$ 在 $0<|z-a|<\delta$ 内解析，且 $1/|f(z)-A|\le1/\varepsilon$，如此，则 $a$ 为 $1/[f(z)-A]$ 的可去奇点，而有

$$
\lim_{z\to a}\frac1{f(z)-A}=B,
$$

其中 $B$ 为有限的复数。若 $B=0$，则得 $\lim_{z\to a}f(z)=\infty$，从而 $a$ 为极点；若 $B\ne0$，则得 $\lim_{z\to a}f(z)=A+1/B$，从而 $a$ 为可去奇点，均与假设矛盾。证毕。

这一定理是 Weierstrass 在 1876 年给出的，它表明，在本性奇点附近，函数值可以无限接近于任何事先给定的复数，包括 $\infty$。1879 年，Picard 给出了下述更为精密的

**定理（Picard）** 若 $a$ 为 $f(z)$ 的本性奇点，则除去可能的一个例外值 $A_0$，对于任意的有限复数 $A$，都有点列 $\{z_n\}_{n=1}^{\infty}$，$\lim_{n\to\infty}z_n=a$，使 $f(z_n)=A$。

这一定理表明，在本性奇点附近，函数值可以无限次地取得任何事先给定的有限复数，最多有一个例外值（称为 Picard 例外值）。可见，在本性奇点附近，分布着非常丰富的函数值，而且是“取之不尽，用之不竭”的。

**例 7** 由例 6 知道，$z=0$ 是函数 $f(z)=e^{1/z}$ 的本性奇点。对任意的有限复数 $A\ne0$，存在点列

$$
z_n=\frac1{\ln A+i2\pi n}\qquad (n\in\mathbb N^+),
$$

满足 $\lim_{n\to\infty}z_n=0$，$f(z_n)=A$。由于 $e^{1/z}\ne0$，故 $A_0=0$ 是 Picard 例外值。但存在点列 $z_n=-1/n$（其中 $n\in\mathbb N^+$），满足 $\lim_{n\to\infty}z_n=0$，

$$
\lim_{n\to\infty}f(z_n)=\lim_{n\to\infty}e^{-n}=0.
$$

另外，存在点列 $z_n=1/n$（其中 $n\in\mathbb N^+$），满足 $\lim_{n\to\infty}z_n=0$，

$$
\lim_{n\to\infty}f(z_n)=\lim_{n\to\infty}e^n=\infty.
$$

**习题** 指出下列函数在 $z$ 平面上的奇点及其类型：

(1) $f(z)=\dfrac1{\sin(1/z)}$；

(2) $f(z)=\dfrac{e^{1/(z-1)}}{e^z-1}$。

## §4 *无穷远点

由于函数在无穷远点没有定义，所以无穷远点总是一个奇点。我们关心的是，什么情况下，无穷远点是孤立奇点？下面是一个非常直观的

**定义** 如果函数 $f(z)$ 在 $\infty$ 的去心邻域 $R<|z|<+\infty$ 内解析，则称 $\infty$ 为 $f(z)$ 的孤立奇点。

在 $\infty$ 是孤立奇点的情况下，我们可以作变换 $t=1/z$，将 $z=\infty$ 化作 $t=0$ 来研究。若 $f(z)$ 在 $R<|z|<+\infty$ 内解析，则 $f(1/t)$ 在 $0<|t|<1/R$（若 $R=0$，则令 $1/R=+\infty$）内解析，即 $t=0$ 是 $f(1/t)$ 的孤立奇点。于是有以下

**定义** 如果 $t=0$ 是 $f(1/t)$ 的可去奇点、$m$ 阶极点或本性奇点，则相应地称 $z=\infty$ 是 $f(z)$ 的可去奇点、$m$ 阶极点或本性奇点。

由于 $f(1/t)$ 在 $0<|t|<1/R$ 内解析，故可展开为 Laurent 级数：

$$
f(1/t)=\sum_{n=-\infty}^{+\infty}c_nt^n,
$$

展开式的主要部分是 $\sum_{n=1}^{+\infty}c_{-n}/t^n$，换回自变量 $z$，就是说，$f(z)$ 在 $R<|z|<+\infty$ 内的 Laurent 展开式的主要部分是 $\sum_{n=1}^{+\infty}c_{-n}z^n$。由此可见：

如果 $z=\infty$ 是 $f(z)$ 的可去奇点，则 $f(z)$ 在 $R<|z|<+\infty$ 内的 Laurent 展开式没有正幂项；

如果 $z=\infty$ 是 $f(z)$ 的 $m$ 阶极点，则 $f(z)$ 在 $R<|z|<+\infty$ 内的 Laurent 展开式的正幂部分为 $\sum_{k=1}^{m}c_{-k}z^k$（其中 $c_{-m}\ne0$）；

如果 $z=\infty$ 是 $f(z)$ 的本性奇点，则 $f(z)$ 在 $R<|z|<+\infty$ 内的 Laurent 展开式有无穷多正幂项。

这些结果表明以上的定义是恰当的，因为在 $\infty$ 附近，正幂项的多少才体现奇性的大小。

如果 $z=\infty$ 是 $f(z)$ 的可去奇点，我们也说 $f(z)$ 在 $\infty$ 解析。所以，所谓 $f(z)$ 在 $\infty$ 解析，无关乎 $f(z)$ 在 $\infty$ 的导数，事实上，我们从未定义过 $f(z)$ 在 $\infty$ 的导数。

**例 1** 函数 $f(z)=1/(z-1)$ 以 $z=\infty$ 为可去奇点，因为 $f(1/t)=1/(1/t-1)=t/(1-t)$ 以 $t=0$ 为可去奇点。另一方面，在 $\infty$ 的去心邻域内，即当 $1<|z|<+\infty$，有 Laurent 展开式

$$
\frac1{z-1}
=\frac1z\frac1{1-1/z}
=\frac1z\sum_{n=0}^{+\infty}\frac1{z^n}
=\sum_{n=1}^{+\infty}\frac1{z^n},
$$

它没有正幂项，这正是可去奇点的特征。

**例 2** $n$ 次多项式 $P_n(z)=\sum_{k=0}^{n}a_kz^k$（其中 $a_n\ne0$）以 $z=\infty$ 为 $n$ 阶极点，因为 $P_n(1/t)=\sum_{k=0}^{n}a_k/t^k$ 以 $t=0$ 为 $n$ 阶极点。另一方面，多项式 $P_n(z)=\sum_{k=0}^{n}a_kz^k$ 本身已具有 $z$ 的 Laurent 级数的形式（在复平面上，即 $|z|<+\infty$ 成立），它只有正幂项，且最高正幂项为 $a_nz^n$，正是 $n$ 阶极点的特征。

**例 3** 指数函数 $e^z$ 以 $z=\infty$ 为本性奇点，因为 $e^{1/t}$ 以 $t=0$ 为本性奇点。另一方面，在复平面上，即 $|z|<+\infty$，有 $e^z=\sum_{n=0}^{+\infty}z^n/n!$，这是通常的 Taylor 展开式，但它也是 $\infty$ 的去心邻域内的 Laurent 展开式，它有无穷多正幂项，这正是本性奇点的特征。

**例 4** $f(z)=\cot z$ 以 $z=\infty$ 为非孤立奇点，因为 $z_n=n\pi$（其中 $n\in\mathbb Z$）是其单极点，而 $z_n\to\infty$（$n\to\infty$）。（参看上节例 4。）

## 补充习题

1. 将下列各函数在指定的环域内展开为 Laurent 级数。

(a) $1/(z-a)(z-b)$（其中 $b\ne a$），$|b-a|<|z-a|<\infty$；

(b) $1/(z-a)(z-b)$（其中 $b\ne a$），$0<|z-a|<|b-a|$；

(c) $e^{1/(1-z)}$，$1<|z|<\infty$（计算五个非零项）；

(d) $1/(z^2-3z+2)$，$1<|z|<2$；

(e) $1/(z^2-3z+2)$，$2<|z|<\infty$。

2. 研究下列函数在 $z$ 平面上的奇点及其类型。

(a) $z^2/(e^z-1)$；

(b) $\cot\pi z/(z-1)$；

(c) $(\sin z-z)/z^3$；

(d) $(e^{az}-e^{bz})/z^2$（其中 $b\ne a$）；

(e) $\cos(1/z^2)$；

(f) $1/\cos z$。
