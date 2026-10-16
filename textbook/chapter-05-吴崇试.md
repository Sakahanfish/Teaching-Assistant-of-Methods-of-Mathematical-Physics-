# 第五章 解析函数的局域性展开

## §5.1 解析函数的 Taylor 展开

前面我们看到，幂级数在它的收敛圆内代表一个解析函数。现在提一个相反的问题：如何把一个解析函数表示成幂级数？

**定理 5.1（Taylor 展开）** 设函数 $f(z)$ 在以 $a$ 为圆心的闭圆 $C$ 内及 $C$ 上解析，则对于圆内的任何 $z$ 点，$f(z)$ 可用幂级数展开为（或者说，$f(z)$ 可在 $a$ 点展开为幂级数）

$$
f(z)=\sum_{n=0}^{\infty}a_n(z-a)^n,
\tag{5.1}
$$

其中

$$
a_n=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta
=\frac{f^{(n)}(a)}{n!}.
\tag{5.2}
$$

$C$ 取逆时针方向①。

**证** 根据 Cauchy 积分公式，对于圆 $C$ 内任意一点 $z$，有

$$
f(z)=\frac{1}{2\pi i}\oint_C\frac{f(\zeta)}{\zeta-z}\,d\zeta.
\tag{5.3}
$$

因为

$$
\frac{1}{\zeta-z}
=\frac{1}{(\zeta-a)-(z-a)}
=\frac{1}{\zeta-a}\sum_{n=0}^{\infty}
\left(\frac{z-a}{\zeta-a}\right)^n,
\qquad |z-a|<|\zeta-a|,
\tag{5.4}
$$

由 (4.16) 式知，级数

$$
\sum_{n=0}^{\infty}
\left(\frac{z-a}{\zeta-a}\right)^n
$$

在区域

$$
\left|\frac{z-a}{\zeta-a}\right|\le r<1
$$

中一致收敛，因此可以逐项积分，

$$
\begin{aligned}
f(z)
&=\frac{1}{2\pi i}\oint_C
\left[
\sum_{n=0}^{\infty}\frac{(z-a)^n}{(\zeta-a)^{n+1}}
\right]f(\zeta)\,d\zeta\\
&=\sum_{n=0}^{\infty}
\left[
\frac{1}{2\pi i}\oint_C
\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta
\right](z-a)^n\\
&=\sum_{n=0}^{\infty}a_n(z-a)^n,
\qquad
\left|\frac{z-a}{\zeta-a}\right|\le r<1,
\end{aligned}
$$

$$
a_n=\frac{1}{2\pi i}\oint_C
\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta
=\frac{f^{(n)}(a)}{n!}.
$$

对于这个定理，需要做以下说明：

1. 定理的条件可以放宽，只要 $f(z)$ 在圆域 $C$ 内解析即可。这时对于给定的圆内 $z$ 点，总可以以 $a$ 为圆心作一圆 $C'$，使 $z\in C'$。于是即可在闭圆域 $\overline{C'}$ 中应用上面的结论。

① 以后的围道积分，若无特殊说明，均为逆时针方向。

2. 这里 Taylor 展开的形式和实变函数中的 Taylor 公式相同，但是条件不同。在实变函数中，$f(x)$ 的任何阶导数存在，还不足以保证 Taylor 公式存在（或 Taylor 级数收敛）。在复变函数中，函数解析的要求就足以保证 Taylor 级数收敛。

3. **收敛范围** 设 $b$ 是 $f(z)$ 的离 $a$ 点最近的奇点，则 $f(z)$ 在圆 $|z-a|<|b-a|$ 内处处解析，$f(z)$ 可以在圆内展开为 Taylor 级数（或者说，Taylor 级数在圆 $|z-a|<|b-a|$ 内收敛）。这就是说，$f(z)$ 在 $a$ 点展开得到的 Taylor 级数收敛半径不小于 $|b-a|$。另一方面，收敛半径一般也不能大于 $|b-a|$。否则，$b$ 点就包含在收敛圆内，因而幂级数在 $b$ 点解析，与 $f(z)$ 在 $b$ 点不解析的假设矛盾（除非 $b$ 点是可去奇点，见 §5.6）。所以，一般说来，收敛半径 $R=|b-a|$。函数 $f(z)$ 的奇点完全决定了其 Taylor 级数的收敛半径。例如

$$
\frac{1}{1+z^2}
=\sum_{n=0}^{\infty}(-1)^nz^{2n},
\qquad |z|<1.
\tag{5.5}
$$

左端的函数在 $z=\pm i$ 不解析，就决定了右端 Taylor 级数的收敛半径 $R=|\pm i|=1$。但是在实数范围内，

$$
\frac{1}{1+x^2}
=\sum_{n=0}^{\infty}(-1)^nx^{2n},
\qquad -1<x<1,
$$

就无法直观地讨论级数的收敛区间与函数性质之间的联系，因为函数 $1/(1+x^2)$ 在整个实轴上都是任意阶连续可导的。

4. **Taylor 展开的唯一性** 给定一个在圆 $C$ 内解析的函数，则它的 Taylor 展开是唯一的，即展开系数 $a_n$ 是完全确定的。

**证** 假定有两个 Taylor 级数在圆 $C$ 内都收敛到同一个解析函数 $f(z)$，

$$
\begin{aligned}
f(z)
&=a_0+a_1(z-a)+a_2(z-a)^2+\cdots+a_n(z-a)^n+\cdots\\
&=a_0'+a_1'(z-a)+a_2'(z-a)^2+\cdots+a_n'(z-a)^n+\cdots.
\end{aligned}
$$

取极限 $z\to a$，则由于级数在圆 $C$ 内内闭一致收敛，故有

$$
a_0=a_0'.
$$

逐项微商，再取极限 $z\to a$，又得

$$
a_1=a_1'.
$$

如此继续，即可证得

$$
a_n=a_n',
\qquad n=0,1,2,\cdots.
$$

Taylor 展开的唯一性告诉我们：(1) 不论用什么方法，得到的 $f(z)$ 在同一个圆内的 Taylor 展开是唯一的。因此，不一定非得用公式 (5.2) 去求展开系数。(2) 如果在同一点展开的两个 Taylor 级数相等，则可以逐项比较系数。这里要强调，必须是在同一点展开的两个 Taylor 级数相等，才可以逐项比较系数。

## §5.2 Taylor 级数求法举例

这里介绍求 Taylor 级数的一些常见方法。为简单起见，本小节举例都在 $z=0$ 邻域内展开。

对于最基本的几个初等函数，需要利用系数公式求出展开系数。由于公式的形式和实变函数中完全相同，因此，可以把实变函数中的结果原封不动地改写成复数形式：

$$
e^z
=1+z+\frac{z^2}{2!}+\cdots+\frac{z^n}{n!}+\cdots
=\sum_{n=0}^{\infty}\frac{z^n}{n!},
\qquad |z|<\infty,
\tag{5.6}
$$

$$
\frac{1}{1-z}
=\sum_{n=0}^{\infty}z^n,
\qquad |z|<1.
\tag{5.7}
$$

对于其他函数，总是尽量利用这些已知的结果，包括作变量代换，或者它们的线性组合和微商、积分。例如

$$
\frac{1}{1+z^2}
=\sum_{n=0}^{\infty}(-z^2)^n
=\sum_{n=0}^{\infty}(-1)^nz^{2n},
\qquad |z|<1,
\tag{5.8}
$$

$$
\sin z
=\frac{e^{iz}-e^{-iz}}{2i}
=\sum_{n=0}^{\infty}\frac{(-1)^n}{(2n+1)!}z^{2n+1},
\qquad |z|<\infty,
\tag{5.9}
$$

$$
\cos z
=\frac{e^{iz}+e^{-iz}}{2}
=\sum_{n=0}^{\infty}\frac{(-1)^n}{(2n)!}z^{2n},
\qquad |z|<\infty,
\tag{5.10}
$$

$$
\frac{1}{1-3z+2z^2}
=-\frac{1}{1-z}+\frac{2}{1-2z}
=\sum_{n=0}^{\infty}(2^{n+1}-1)z^n,
\qquad |z|<\frac12,
\tag{5.11}
$$

$$
\frac{1}{(1-z)^2}
=\frac{d}{dz}\frac{1}{1-z}
=\frac{d}{dz}\sum_{n=0}^{\infty}z^n
=\sum_{n=0}^{\infty}(n+1)z^n,
\qquad |z|<1.
\tag{5.12}
$$

下面介绍两种新的方法，即级数乘法和待定系数法。

如果一个函数可以表示成两个（或多个）函数的乘积，而每一个因子的 Taylor 展开比较容易求出，则可以采用级数相乘的方法。例如

$$
\begin{aligned}
\frac{1}{1-3z+2z^2}
&=\frac{1}{1-z}\cdot\frac{1}{1-2z}
=\sum_{k=0}^{\infty}z^k\cdot\sum_{l=0}^{\infty}2^lz^l\\
&=\sum_{k=0}^{\infty}\sum_{l=0}^{\infty}2^lz^{k+l}\\
&=\sum_{n=0}^{\infty}\left(\sum_{l=0}^{n}2^l\right)z^n
=\sum_{n=0}^{\infty}(2^{n+1}-1)z^n,
\qquad |z|<\frac12.
\end{aligned}
\tag{5.13}
$$

幂级数在收敛圆内绝对收敛，故乘积在两收敛圆的公共区域内仍绝对收敛。

关于待定系数法，见下面的例子。

**例 5.1** 求 $\tan z$ 在 $z=0$ 的 Taylor 展开。

**解** 由于 $\tan z$ 是奇函数，故其在 $z=0$ 的 Taylor 展开应只有奇次幂，

$$
\tan z=\sum_{k=0}^{\infty}a_{2k+1}z^{2k+1}.
$$

因此，

$$
\sin z=\cos z\cdot\sum_{k=0}^{\infty}a_{2k+1}z^{2k+1},
$$

即

$$
\sum_{n=0}^{\infty}\frac{(-1)^n}{(2n+1)!}z^{2n+1}
=
\sum_{l=0}^{\infty}\frac{(-1)^l}{(2l)!}z^{2l}
\cdot
\sum_{k=0}^{\infty}a_{2k+1}z^{2k+1}
=
\sum_{n=0}^{\infty}
\left(
\sum_{k=0}^{n}
\frac{(-1)^{n-k}}{(2n-2k)!}a_{2k+1}
\right)z^{2n+1}.
$$

根据 Taylor 展开的唯一性，可以将上式中左右两式比较系数，由此即得

$$
\sum_{k=0}^{n}
\frac{(-1)^k}{(2n-2k)!}a_{2k+1}
=\frac{1}{(2n+1)!},
\qquad n=0,1,2,\cdots.
$$

所以

$$
\begin{array}{lll}
n=0, & a_1=1;\\[4pt]
n=1, & \dfrac{1}{2}a_1-a_3=\dfrac{1}{6}, & a_3=\dfrac13;\\[8pt]
n=2, & \dfrac{1}{24}a_1-\dfrac12a_3+a_5=\dfrac{1}{120}, & a_5=\dfrac{2}{15};\\
\vdots & \vdots & \vdots
\end{array}
$$

最后就得到

$$
\tan z=z+\frac13z^3+\frac{2}{15}z^5+\frac{17}{315}z^7+\cdots.
\tag{5.14}
$$

从 $\tan z$ 的奇点可以判断，级数的收敛半径应为 $\pi/2$。

应用待定系数法，能得到系数之间的递推关系，从而逐个求出展开系数，但一般很难求出级数的通项公式（即展开系数 $a_n$ 的解析表达式）。然而，如果我们只需要求出级数中的某一项或某几项系数，待定系数法还不失为可取的方法之一。

**多值函数的 Taylor 展开** 对于多值函数，在适当规定了单值分支后，即可像单值函数那样在解析点邻域内作 Taylor 展开。

**例 5.2** 求多值函数 $(1+z)^\alpha$ 在 $z=0$ 的 Taylor 展开，规定 $z=0$ 时

$$
(1+z)^\alpha\big|_{z=0}=1.
$$

**解** 可直接求出函数 $(1+z)^\alpha$ 在 $z=0$ 点的各阶导数值，

$$
f(0)=1,
$$

$$
f'(0)=\alpha(1+z)^{\alpha-1}\big|_{z=0}=\alpha,
$$

$$
f''(0)=\alpha(\alpha-1)(1+z)^{\alpha-2}\big|_{z=0}
=\alpha(\alpha-1),
$$

$$
\vdots
$$

$$
f^{(n)}(0)
=\alpha(\alpha-1)(\alpha-2)\cdots(\alpha-n+1)
(1+z)^{\alpha-n}\big|_{z=0}
=\alpha(\alpha-1)\cdots(\alpha-n+1),
$$

$$
\vdots
$$

因此

$$
(1+z)^\alpha
=\sum_{n=0}^{\infty}\binom{\alpha}{n}z^n,
\tag{5.15}
$$

其中

$$
\binom{\alpha}{0}=1,
\qquad
\binom{\alpha}{n}
=\frac{\alpha(\alpha-1)\cdots(\alpha-n+1)}{n!}
$$

称为普遍的二项式展开系数。级数的收敛区域，应视割线的作法而定。收敛半径等于 $z=0$ 到割线的最短距离。最大可能的收敛区域是 $|z|<1$。所以，只要许可，我们总是会将割线作在单位圆外。

**例 5.3** 求多值函数 $\ln(1+z)$ 在 $z=0$ 的 Taylor 展开，规定

$$
\ln(1+z)\big|_{z=0}=0.
$$

**解** 在上述规定下，函数 $\ln(1+z)$ 可表示为定积分，因此

$$
\begin{aligned}
\ln(1+z)
&=\ln(1+z)-\ln(1+z)\big|_{z=0}
=\int_0^z\frac{1}{1+\zeta}\,d\zeta\\
&=\int_0^z\sum_{n=0}^{\infty}(-1)^n\zeta^n\,d\zeta
=\sum_{n=0}^{\infty}(-1)^n\int_0^z\zeta^n\,d\zeta\\
&=\sum_{n=1}^{\infty}\frac{(-1)^{n-1}}{n}z^n.
\end{aligned}
\tag{5.16}
$$

收敛区域也要看割线怎么作。收敛半径等于 $z=0$ 到割线的最短距离，最大可能的收敛区域是 $|z|<1$。同样，只要许可，我们总是将割线作在单位圆外。

**在无穷远点的 Taylor 展开** 如果函数 $f(z)$ 在 $z=\infty$ 点解析，则也可以在 $z=\infty$ 点展开成 Taylor 级数。

所谓 $f(z)$ 在 $\infty$ 点展开成 Taylor 级数，完全等价于作变换 $z=1/t$，而将 $f(1/t)$ 在 $t=0$ 点展开成 Taylor 级数。因为 $f(1/t)$ 在 $t=0$ 点解析，故

$$
f\left(\frac1t\right)
=a_0+a_1t+a_2t^2+\cdots+a_nt^n+\cdots,
\qquad |t|<r;
\tag{5.17a}
$$

即

$$
f(z)
=a_0+\frac{a_1}{z}+\frac{a_2}{z^2}
+\cdots+\frac{a_n}{z^n}+\cdots,
\qquad |z|>\frac1r.
\tag{5.17b}
$$

值得注意的是，$f(z)$ 在 $\infty$ 点的 Taylor 级数中只有常数项及负幂项，没有正幂项，而收敛范围为 $|z|>1/r$，也就是说，级数在以 $\infty$ 为圆心的某个圆内收敛。

**例 5.4** 将函数 $1/(1-z^2)$ 在 $\infty$ 点作 Taylor 展开。

**解** 因为在代换 $z=1/t$ 之下，

$$
\frac{1}{1-z^2}=-\frac{t^2}{1-t^2},
$$

$$
-\frac{t^2}{1-t^2}
=-t^2\sum_{n=0}^{\infty}t^{2n}
=-\sum_{n=1}^{\infty}t^{2n},
\qquad |t|<1,
$$

所以，$1/(1-z^2)$ 在 $\infty$ 点的 Taylor 展开就是

$$
\frac{1}{1-z^2}
=-\sum_{n=1}^{\infty}z^{-2n},
\qquad |z|>1.
\tag{5.18}
$$

其实，在作展开时也可不必写出代换 $z=1/t$，而只需记住当 $|z|>1$ 时 $|1/z|<1$，

$$
\frac{1}{1-z^2}
=-\frac{1}{z^2}\frac{1}{1-\frac1{z^2}}
=-\frac{1}{z^2}\sum_{n=0}^{\infty}\frac{1}{z^{2n}}
=-\sum_{n=1}^{\infty}\frac{1}{z^{2n}},
\qquad |z|>1.
$$

再举一个多值函数的例子。

**例 5.5** 将函数 $1/\sqrt{z^2-1}$ 在 $z=\infty$ 作 Taylor 展开，约定直接连接 $z=\pm1$ 作割线，且规定当 $z$ 位于 $z=1$ 之右的实轴上时函数取正值。

**解** 在此单值分支的规定下，显然有

$$
\begin{aligned}
\frac{1}{\sqrt{z^2-1}}
&=\frac1z\left(1-\frac1{z^2}\right)^{-1/2}
=\frac1z\sum_{n=0}^{\infty}
\binom{-1/2}{n}\left(-\frac1{z^2}\right)^n\\
&=\sum_{n=0}^{\infty}
(-1)^n\binom{-1/2}{n}\frac{1}{z^{2n+1}},
\qquad |z|>1.
\end{aligned}
$$

还可以将展开系数改写为

$$
\begin{aligned}
(-1)^n\binom{-1/2}{n}
&=\frac{(-1)^n}{n!}
\left(-\frac12\right)
\left(-\frac32\right)\cdots
\left(\frac12-n\right)\\
&=\frac{1}{2^n}\frac{(2n-1)!!}{n!}
=\frac{1}{2^{2n}}\frac{(2n)!}{n!\,n!}.
\end{aligned}
$$

## §5.3 解析函数的零点孤立性和解析函数的唯一性

如果 $f(z)$ 在 $a$ 点的邻域内解析且不恒为 0，若 $f(a)=0$，则称 $z=a$ 为 $f(z)$ 的零点。

设 $f(z)$ 在 $z=a$ 点的邻域内解析，则 $f(z)$ 可以在 $z=a$ 的邻域内展成 Taylor 级数，

$$
f(z)=\sum_{n=0}^{\infty}a_n(z-a)^n,
\qquad |z-a|<\rho.
\tag{5.19}
$$

故若 $z=a$ 为零点，则必有

$$
a_0=a_1=\cdots=a_{m-1}=0,
\qquad a_m\ne0.
\tag{5.20}
$$

此时，称 $z=a$ 点为 $f(z)$ 的 $m$ 阶零点，相应地，

$$
f(a)=f'(a)=\cdots=f^{(m-1)}(a)=0,
\qquad f^{(m)}(a)\ne0.
\tag{5.21}
$$

可见 $m$ 不可能是零或者负数。因为 $f(z)$ 不恒为零，所以 $m$ 也不可能为 $\infty$。$m$ 更不可能是非整数①，因此零点的阶数 $m$ 一定是正整数。

解析函数零点的一个重要性质是它的孤立性。

**定理 5.2（解析函数的零点孤立性定理）** 若 $z=a$ 是 $f(z)$ 的零点且 $f(z)$ 在 $z=a$ 的邻域内不恒等于零，则一定 $\exists\rho>0$，使得 $f(z)$ 在空心邻域 $0<|z-a|<\rho$ 内无零点。

① 对于 $(z-a)^s$，当 $s\ne$ 整数时，$z=a$ 是分支点，函数在该点并不解析，因而不是这里所讨论的解析函数的零点。

**证** 设 $a$ 为 $f(z)$ 的 $m$ 阶零点，则

$$
f(z)=(z-a)^m\phi(z),
$$

$\phi(z)$ 在 $|z-a|<R$ 内解析，且 $\phi(a)\ne0$。因为 $\phi(z)$ 在 $z=a$ 点连续，即 $\forall\varepsilon>0$，$\exists\rho>0$，使得当 $|z-a|<\rho$ 时，恒有 $|\phi(z)-\phi(a)|<\varepsilon$。不妨取 $\varepsilon=|\phi(a)|/2$，则得

$$
|\phi(z)|>|\phi(a)|-\varepsilon
=\frac12|\phi(a)|>0.
$$

由此即证得 $f(z)$ 在 $0<|z-a|<\rho$ 内无零点。

解析函数的零点孤立性定理的逆否定理是：若 $\exists R>0$，$f(z)$ 在 $|z-a|<R$ 内解析，$\forall\rho>0$，$f(z)$ 在空心邻域 $0<|z-a|<\rho$ 内都有零点，那么 $f(z)$ 在 $z=a$ 的空心邻域 $0<|z-a|<R$ 内恒等于零。也就是说，如果解析函数 $f(z)$ 的零点是非孤立的，则此函数在其解析区域内一定恒为 0。这个结论还可以表述为下面的几个推论。

**推论 5.1** 设 $f(z)$ 在 $G:|z-a|<R$ 内解析。若在 $G$ 内存在 $f(z)$ 的无穷多个互不相等的零点 $\{z_n\}$，且 $\lim\limits_{n\to\infty}z_n=a$，但 $z_n\ne a$，则在 $G$ 内 $f(z)\equiv0$。

**推论 5.2** 设 $f(z)$ 在 $G:|z-a|<R$ 内解析。若在 $G$ 内存在过 $a$ 点的一段弧 $l$ 或含有 $a$ 点的一个子区域 $g$，在 $l$ 上或 $g$ 内 $f(z)\equiv0$，则在整个区域 $G$ 内 $f(z)\equiv0$。

**推论 5.3** 设 $f(z)$ 在 $G$ 内解析。若在 $G$ 内存在一点 $z=a$ 及过 $a$ 点的一段弧 $l$ 或含有 $a$ 点的一个子区域 $g$，在 $l$ 上或 $g$ 内 $f(z)\equiv0$，则在整个区域 $G$ 内 $f(z)\equiv0$。

**推论 5.4** 设 $f_1(z)$ 和 $f_2(z)$ 都在区域 $G$ 内解析，且在 $G$ 内的一段弧或一个子区域内相等，则在 $G$ 内 $f_1(z)\equiv f_2(z)$。

**推论 5.5** 在实轴上成立的恒等式，在 $z$ 复平面上仍然成立，只要这个恒等式两端的函数在 $z$ 复平面上都是解析的。

也还可以把推论 1 改写成解析函数的唯一性定理。

**定理 5.3（解析函数的唯一性定理）** 设在区域 $G$ 内有两个解析函数 $f_1(z)$ 和 $f_2(z)$，且在 $G$ 内存在序列 $\{z_n\}$，$\forall n$，$f_1(z_n)=f_2(z_n)$。若 $\{z_n\}$ 的一个极限点 $z=a$（$\ne z_n$）也落在 $G$ 内，则在 $G$ 内有 $f_1(z)\equiv f_2(z)$。

## §5.4 解析函数的 Laurent 展开

一个函数除了可在解析点的邻域（单连通区域）内作 Taylor 展开外，有时还需要将它在环形区域（多连通区域）展开成幂级数。这时就得到 Laurent 展开。

**定理 5.4（Laurent 展开）** 设函数 $f(z)$ 在以 $b$ 为圆心的环形区域 $R_1\le |z-b|\le R_2$ 中单值解析，则对于环域内的任何 $z$ 点，$f(z)$ 可以用幂级数展开为

$$
f(z)=\sum_{n=-\infty}^{\infty}a_n(z-b)^n,
\qquad R_1<|z-b|<R_2,
\tag{5.22}
$$

其中

$$
a_n=\frac{1}{2\pi i}\oint_C
\frac{f(\zeta)}{(\zeta-b)^{n+1}}\,d\zeta,
\tag{5.23}
$$

$C$ 是环域内绕内圆一周的任意一条闭合曲线（见图 5.1）。

![Figure 5.1](...)

图 5.1 Laurent 展开

**证** 将环域内、外边界分别记为 $C_1$ 和 $C_2$，则根据（多连通区域的）Cauchy 积分公式，有

$$
f(z)
=\frac{1}{2\pi i}\oint_{C_2}\frac{f(\zeta)}{\zeta-z}\,d\zeta
-\frac{1}{2\pi i}\oint_{C_1}\frac{f(\zeta)}{\zeta-z}\,d\zeta.
$$

对于 $C_2$ 上的积分，可以直接引用证明 Taylor 展开定理时得到的结果，

$$
\frac{1}{2\pi i}\oint_{C_2}
\frac{f(\zeta)}{\zeta-z}\,d\zeta
=\sum_{n=0}^{\infty}a_n(z-b)^n,
\qquad |z-b|<R_2,
$$

$$
a_n=\frac{1}{2\pi i}\oint_{C_2}
\frac{f(\zeta)}{(\zeta-b)^{n+1}}\,d\zeta.
$$

对于 $C_1$ 上的积分，

$$
\begin{aligned}
-\frac{1}{2\pi i}\oint_{C_1}\frac{f(\zeta)}{\zeta-z}\,d\zeta
&=\frac{1}{2\pi i}\oint_{C_1}
\frac{f(\zeta)}{(z-b)-(\zeta-b)}\,d\zeta\\
&=\frac{1}{2\pi i}\oint_{C_1}
\frac{f(\zeta)}{z-b}
\sum_{k=0}^{\infty}
\left(\frac{\zeta-b}{z-b}\right)^k\,d\zeta\\
&=\sum_{k=0}^{\infty}(z-b)^{-k-1}
\frac{1}{2\pi i}\oint_{C_1}f(\zeta)(\zeta-b)^k\,d\zeta\\
&=\sum_{n=-1}^{-\infty}a_n(z-b)^n,
\qquad |z-b|>R_1,
\end{aligned}
$$

其中

$$
a_n=\frac{1}{2\pi i}\oint_{C_1}
\frac{f(\zeta)}{(\zeta-b)^{n+1}}\,d\zeta.
$$

把两部分合并起来，就有

$$
f(z)=\sum_{n=-\infty}^{\infty}a_n(z-b)^n,
\qquad R_1<|z-b|<R_2,
$$

$$
a_n=\frac{1}{2\pi i}\oint_C
\frac{f(\zeta)}{(\zeta-b)^{n+1}}\,d\zeta.
$$

这里把系数 $a_n$ 公式中的积分围道统一写成了 $C$，其中 $C$ 是环域内绕内圆一周的任意一条闭合曲线（为什么能这样做？）。

这个结果称为函数 $f(z)$ 在环域 $R_1<|z-b|<R_2$ 内的 Laurent 展开，展开而得到的级数称为 Laurent 级数。

对于上面的结果，也需要做一些补充讨论。

1. 和 Taylor 展开一样，本定理的条件也可放宽为 $f(z)$ 在 $R_1<|z-b|<R_2$ 内单值解析。

2. 对于 Laurent 展开来说，系数（即使是正幂项的系数）

$$
a_n\ne\frac{1}{n!}f^{(n)}(b).
$$

3. $f(z)$ 在内圆 $\overline{C_1}$ 中不解析。一般说来，在 $C_1$ 圆周上有奇点。在圆心 $b$ 点，$f(z)$ 可能解析，也可能不解析。

如果 $b$ 点是 $f(z)$ 在圆 $\overline{C_1}$ 中唯一的奇点，则 $C_1$ 的半径可以无限缩小，收敛范围就变成 $0<|z-b|<R$。这时圆心 $b$ 点称为函数 $f(z)$ 的孤立奇点，就得到 $f(z)$ 在孤立奇点 $b$ 的空心邻域内的 Laurent 展开。同样，外圆 $C_2$ 的半径也可能为 $\infty$。

4. Laurent 展开既有正幂项，又有负幂项。正幂项在外圆 $C_2$ 内（$|z-b|<R_2$）绝对收敛，并且内闭一致收敛，称为 Laurent 级数的正则部分；负幂项在内圆 $C_1$ 外（$|z-b|>R_1$）绝对收敛，在 $C_1$ 外的任一闭区域中一致收敛，称为 Laurent 级数的主要部分。两部分合起来，就构成 Laurent 级数，在环域 $R_1<|z-b|<R_2$ 内绝对收敛，并且内闭一致收敛。当内圆半径 $R_1=0$ 时，Laurent 级数的主要部分就完全反映了 $f(z)$ 在孤立奇点 $z=b$ 点的奇异性。

5. 如果函数 $f(z)$ 在无穷远点不解析，但是在无穷远点的一个空心邻域内单值解析，则可将 $f(z)$ 在 $\infty$ 点的空心邻域内作 Laurent 展开（或者简单地说成在 $\infty$ 点作 Laurent 展开）。

所谓 $f(z)$ 在 $\infty$ 点的空心邻域内单值解析，就意味着函数 $f(1/t)$ 在 $t=0$ 点的空心邻域内单值解析，因而

$$
f\left(\frac1t\right)
=\sum_{n=-\infty}^{\infty}a_nt^n,
\qquad 0<|t|<r,
$$

即

$$
f(z)=\sum_{n=-\infty}^{\infty}a_nz^{-n},
\qquad \frac1r<|z|<\infty.
\tag{5.24}
$$

这里的收敛范围可以理解为是以 $z=\infty$ 点为圆心的一个环域。$f(1/t)$ 在 $t=0$ 点的 Laurent 级数中 $t$ 的正幂项（包括常数项）部分是正则部分；$t$ 的负幂项是主要部分，完全反映了 $f(1/t)$ 在孤立奇点 $t=0$ 点的奇异性。对应地，我们把 $f(z)$ 在 $z=\infty$ 点空心邻域内的 Laurent 级数中，$z$ 的负幂项称为正则部分；$z$ 的正幂项称为主要部分，完全反映了函数 $f(z)$ 在 $\infty$ 点的奇异性。

6. **Laurent 展开的唯一性** 设 $f(z)$ 在环域 $R_1<|z-b|<R_2$ 内有两个 Laurent 级数，

$$
f(z)
=\sum_{n=-\infty}^{\infty}a_n(z-b)^n
=\sum_{n=-\infty}^{\infty}a_n'(z-b)^n.
$$

两端同乘以 $(z-b)^{-k-1}$，其中 $k$ 为任意整数，沿环域内绕内圆一周的任一围道 $C$ 积分（这两个级数在环域内内闭一致收敛，因而可以逐项积分），则由于

$$
\oint_C(z-b)^{n-k-1}\,dz=2\pi i\delta_{nk},
$$

故有 $a_k=a_k'$。因为 $k$ 任意，故有

$$
a_k=a_k',
\qquad k=0,\pm1,\pm2,\cdots.
\tag{5.25}
$$

即证得 Laurent 展开的唯一性。它告诉我们：两个 Laurent 级数在同一环域内处处相等，则对应项系数相等（即可以比较系数）。

## §5.5 Laurent 级数求法举例

求 Laurent 展开，一般不直接利用公式求系数。由于函数在给定环域内 Laurent 展开唯一性的保证，因此，不论用什么方法，包括引用 Taylor 展开中得到的结果，只要最终得到的是在给定环域内收敛到 $f(z)$ 的幂级数，那它就一定是 $f(z)$ 的 Laurent 展开。

**例 5.6** 求

$$
\frac{1}{z(z-1)}
$$

在 $0<|z|<1$ 内和 $|z|>1$ 内的展开式。

**解** $1/[z(z-1)]$ 在 $0<|z|<1$ 内的展开形式一定是 $\sum_{n=-\infty}^{\infty}a_nz^n$，所以

$$
\frac{1}{z(z-1)}
=-\frac1z\frac{1}{1-z}
=-\frac1z\sum_{n=0}^{\infty}z^n
=-\sum_{n=-1}^{\infty}z^n,
\qquad 0<|z|<1.
$$

也可以用部分分式的方法：

$$
\frac{1}{z(z-1)}
=-\frac1z-\frac{1}{1-z}
=-\frac1z-\sum_{n=0}^{\infty}z^n
=-\sum_{n=-1}^{\infty}z^n,
\qquad 0<|z|<1.
$$

$1/[z(z-1)]$ 在 $|z|>1$ 内解析，在 $|z|>1$ 内幂级数展开形式也是 $\sum_{n=-\infty}^{\infty}a_nz^n$，

$$
\frac{1}{z(z-1)}
=\frac{1}{z^2}\frac{1}{1-\frac1z}
=\frac{1}{z^2}\sum_{n=0}^{\infty}\left(\frac1z\right)^n
=\sum_{n=-\infty}^{-2}z^n,
\qquad |z|>1.
$$

此结果也可以看成是 $1/[z(z-1)]$ 在 $\infty$ 点邻域内的 Taylor 展开。

**例 5.7** 用待定系数法求 $\cot z$ 在 $z=0$ 空心邻域内的 Laurent 展开。

**解** 待定系数法只能用于只有有限个负幂项（或有限个正幂项）的情形。

$$
\cot z=\sum_{n=-1}^{\infty}b_{2n+1}z^{2n+1}.
\tag{5.26}
$$

（为什么只有一个负幂项，这个道理将在 §5.6 讨论）。

$$
\cos z=\sin z\sum_{n=-1}^{\infty}b_{2n+1}z^{2n+1},
$$

$$
\sum_{n=0}^{\infty}\frac{(-1)^n}{(2n)!}z^{2n}
=
\sum_{k=0}^{\infty}\frac{(-1)^k}{(2k+1)!}z^{2k+1}
\sum_{l=0}^{\infty}b_{2l-1}z^{2l-1}
=
\sum_{n=0}^{\infty}
\left[
\sum_{l=0}^{n}
\frac{(-1)^{n-l}}{(2n-2l+1)!}b_{2l-1}
\right]z^{2n}.
$$

由此得到递推关系

$$
\sum_{l=0}^{n}
\frac{(-1)^l}{(2n-2l+1)!}b_{2l-1}
=\frac{1}{(2n)!}.
\tag{5.27}
$$

逐次求解，即得

$$
\begin{array}{lll}
n=0, & b_{-1}=1;\\[4pt]
n=1, & \dfrac1{3!}b_{-1}-\dfrac1{1!}b_1=\dfrac1{2!}, & b_1=-\dfrac13;\\[8pt]
n=2, & \dfrac1{5!}b_{-1}-\dfrac1{3!}b_1+\dfrac1{1!}b_3=\dfrac1{4!}, & b_3=-\dfrac1{45};\\[8pt]
n=3, & \dfrac1{7!}b_{-1}-\dfrac1{5!}b_1+\dfrac1{3!}b_3-\dfrac1{1!}b_5=\dfrac1{6!}, & b_5=-\dfrac2{945};\\
\vdots & \vdots & \vdots
\end{array}
$$

所以

$$
\cot z=\frac1z-\frac13z-\frac1{45}z^3-\frac2{945}z^5-\cdots.
\tag{5.28}
$$

根据 $\cot z$ 的奇点分布，可判断此级数的收敛范围为 $0<|z|<\pi$。

本题还可以采用级数除法。当 $|z|>0$ 而又足够小时，

$$
\begin{aligned}
\cot z
&=\frac1{\tan z}
=\left(z+\frac13z^3+\frac2{15}z^5+\frac{17}{315}z^7+\cdots\right)^{-1}\\
&=\frac1z\left(1+\frac13z^2+\frac2{15}z^4+\frac{17}{315}z^6+\cdots\right)^{-1}\\
&=\frac1z\left[
1-\left(\frac13z^2+\frac2{15}z^4+\frac{17}{315}z^6+\cdots\right)
+\left(\frac13z^2+\frac2{15}z^4+\frac{17}{315}z^6+\cdots\right)^2\right.\\
&\qquad\left.
-\left(\frac13z^2+\frac2{15}z^4+\frac{17}{315}z^6+\cdots\right)^3+\cdots
\right]\\
&=\frac1z\left[
1-\frac13z^2
+\left(-\frac2{15}+\frac19\right)z^4
+\left(-\frac{17}{315}+2\times\frac13\times\frac2{15}-\frac1{27}\right)z^6
+\cdots
\right]\\
&=\frac1z-\frac13z-\frac1{45}z^3-\frac2{945}z^5+\cdots.
\end{aligned}
$$

这个计算过程显然只在区域

$$
0<\left|\frac{\tan z}{z}-1\right|<1
$$

（明显比 $0<|z|<\pi$ 小，因为 $\tan z$ 在 $z=\pi/2$ 点不解析）内有效，但最后结果却在更大的区域 $0<|z|<\pi$ 内成立，理由见 §5.7。

**例 5.8** 求函数

$$
f(z)=\ln\frac{z-2}{z-1}
$$

在 $1<|z|<2$ 及 $2<|z|<\infty$ 内的幂级数展开。

**解** 本题中指定的展开区域是环形区域，是多连通区域，所以，如果能作幂级数展开的话，得到的一定是 Laurent 级数。

$f(z)$ 有两个分支点：$z=1$ 和 $z=2$。$f(z)$ 的割线一定要连接这两个分支点。不论割线怎么作，它一定会穿过环域 $1<|z|<2$，换言之，$f(z)$ 在 $z$ 复平面内的环域 $1<|z|<2$ 内一定不解析，故在此环域内不可能作 Laurent 展开。

在环域 $2<|z|<\infty$ 内，如果割线不通过此环域，则规定单值分支后，$f(z)$ 单值解析，此时可作 Laurent 展开。例如，若沿实轴从 $z=1$ 到 $z=2$ 作割线，规定在割线上岸

$$
\arg(z-2)=\arg(z-1)=\pi,
$$

则

$$
f(\infty)
=\left.\ln\frac{z-2}{z-1}\right|_{z=\infty}=0.
$$

于是有

$$
\begin{aligned}
\ln\frac{z-2}{z-1}
&=\int_{\infty}^{z}
\left(\frac1{\zeta-2}-\frac1{\zeta-1}\right)d\zeta\\
&=\int_{\infty}^{z}\frac1{\zeta}
\sum_{n=0}^{\infty}
\left[
\left(\frac2{\zeta}\right)^n
-\left(\frac1{\zeta}\right)^n
\right]d\zeta\\
&=-\frac1z-\frac32\frac1{z^2}
-\frac73\frac1{z^3}
-\cdots
-\frac{2^n-1}{n}\frac1{z^n}
-\cdots.
\end{aligned}
\tag{5.29}
$$

最大可能的收敛范围是 $|z|>2$。这是 $\ln[(z-2)/(z-1)]$ 在 $\infty$ 点的 Taylor 级数。

**思考题** 若对单值分支做其他规定，函数 $\ln[(z-2)/(z-1)]$ 在 $|z|>2$ 内的幂级数展开和 (5.29) 式有何异同？

**例 5.9** 求函数

$$
\exp\left\{\frac z2\left(t-\frac1t\right)\right\}
$$

在环域 $0<|t|<\infty$ 内的 Laurent 展开。

**解** 用级数乘法。因为

$$
e^{zt/2}
=\sum_{k=0}^{\infty}
\left(\frac z2\right)^k\frac{t^k}{k!},
\qquad |t|<\infty,
$$

$$
e^{-z/2t}
=\sum_{l=0}^{\infty}
\left(\frac z2\right)^l
\frac{(-1)^l}{l!}\left(\frac1t\right)^l,
\qquad
\left|\frac1t\right|<\infty
\quad\text{即}\quad |t|>0,
$$

所以

$$
\begin{aligned}
\exp\left\{\frac z2\left(t-\frac1t\right)\right\}
&=
\sum_{k=0}^{\infty}
\left(\frac z2\right)^k\frac{t^k}{k!}
\sum_{l=0}^{\infty}
\left(\frac z2\right)^l\frac{(-1)^l}{l!}\left(\frac1t\right)^l\\
&=
\sum_{k=0}^{\infty}\sum_{l=0}^{\infty}
\frac{(-1)^l}{k!\,l!}
\left(\frac z2\right)^{k+l}t^{k-l}\\
&=
\sum_{n=0}^{\infty}
\left[
\sum_{l=0}^{\infty}
\frac{(-1)^l}{l!(l+n)!}
\left(\frac z2\right)^{2l+n}
\right]t^n\\
&\qquad+
\sum_{n=-1}^{-\infty}
\left[
\sum_{l=-n}^{\infty}
\frac{(-1)^l}{l!(l+n)!}
\left(\frac z2\right)^{2l+n}
\right]t^n\\
&=\sum_{n=-\infty}^{\infty}J_n(z)t^n,
\end{aligned}
\tag{5.30}
$$

其中

$$
J_n(z)=
\begin{cases}
\displaystyle
\sum_{l=0}^{\infty}
\frac{(-1)^l}{l!(l+n)!}
\left(\frac z2\right)^{2l+n},
& n=0,1,2,\cdots,\\[12pt]
\displaystyle
\sum_{l=-n}^{\infty}
\frac{(-1)^l}{l!(l+n)!}
\left(\frac z2\right)^{2l+n},
& n=-1,-2,-3,\cdots.
\end{cases}
\tag{5.31}
$$

称为 $n$ 阶 Bessel 函数。(5.30) 式既可以看成是函数 $\exp[z(t-1/t)/2]$ 在 $z=0$ 点空心邻域内的 Laurent 展开（有无穷多个负幂项），又可以看成是该函数在 $\infty$ 点空心邻域内的 Laurent 展开（有无穷多个正幂项）。

**思考题** 求 $f(z)=\sin z$ 在 $\infty$ 点空心邻域内的 Laurent 展开。

## §5.6 单值函数的孤立奇点

设 $f(z)$ 为单值函数（或多值函数的一个单值分支），如果 $f(z)$ 在 $b$ 点不解析，但是 $\exists r>0$，$f(z)$ 在 $b$ 点的空心邻域 $0<|z-b|<r$ 内处处可导，则称 $b$ 点为 $f(z)$ 的孤立奇点。反之，如果对 $\forall r>0$，在 $\beta$ 点的空心邻域 $0<|z-\beta|<r$ 内都有 $f(z)$ 的奇点，则称 $\beta$ 点为 $f(z)$ 的非孤立奇点。

孤立奇点的例子已经见过很多，这里举一个非孤立奇点的例子。我们知道，对于 $1/\sin(1/z)$，$z=1/(n\pi)$（即 $1/z=n\pi$），$n=\pm1,\pm2,\cdots$ 均为奇点。显然，$z=0$ 点是这些奇点的聚点：在 $z=0$ 的任一空心邻域中，总存在函数 $1/\sin(1/z)$ 的奇点，故 $z=0$ 是函数 $1/\sin(1/z)$ 的非孤立奇点。

$z=b$ 是单值函数 $f(z)$ 的孤立奇点，意味着一定存在 $R>0$，$f(z)$ 在环域 $0<|z-b|<R$ 内可以展开成 Laurent 级数：

$$
f(z)=\sum_{n=-\infty}^{\infty}a_n(z-b)^n,
\qquad 0<|z-b|<R.
\tag{5.32}
$$

这时可能出现三种情况：

(1) 级数展开式 (5.32) 不含负幂项，此时 $b$ 点称为 $f(z)$ 的可去奇点。例如，$z=0$ 就是函数

$$
\frac{\sin z}{z}
=\sum_{n=0}^{\infty}\frac{(-1)^n}{(2n+1)!}z^{2n},
\qquad |z|<\infty
$$

和

$$
\frac1z-\cot z
=\frac13z+\frac1{45}z^3+\frac2{945}z^5+\cdots,
\qquad |z|<\pi
$$

的可去奇点。

(2) 级数展开式 (5.32) 只含有限个负幂项，此时 $b$ 点称为 $f(z)$ 的极点。

(3) 级数展开式 (5.32) 含有无穷多个负幂项，此时 $b$ 点称为 $f(z)$ 的本性奇点。

下面分别讨论函数在三种孤立奇点处的行为。

**可去奇点** 由于在可去奇点处的级数展开不含负幂项，故级数不但在环域内收敛，而且在环域的中心，即可去奇点 $z=b$ 处也是收敛的。也就是说，这时的收敛区域其实是一个圆形区域，圆心位于可去奇点 $z=b$ 处，级数在收敛圆内内闭一致收敛，因而其和函数连续，

$$
\lim_{z\to b}f(z)
=\lim_{z\to b}\sum_{n=0}^{\infty}a_n(z-b)^n
=a_0,
\tag{5.33}
$$

即函数在可去奇点处的极限值是有限的。如果定义一个新的函数

$$
F(z)=
\begin{cases}
f(z), & z\ne b,\\[4pt]
\displaystyle \lim_{z\to b}f(z), & z=b,
\end{cases}
\tag{5.34}
$$

这样函数 $F(z)$ 在 $b$ 点也是解析的。这正是可去奇点这一称谓的由来。

反过来说，如果 $z=b$ 是函数 $f(z)$ 的孤立奇点，而且 $f(z)$ 在 $z=b$ 的邻域内有界，则 $z=b$ 是 $f(z)$ 的可去奇点。

**极点** 函数在极点空心邻域内的 Laurent 展开有有限个负幂项，

$$
\begin{aligned}
f(z)
&=a_{-m}(z-b)^{-m}+a_{-m+1}(z-b)^{-m+1}
+\cdots+a_{-1}(z-b)^{-1}\\
&\qquad +a_0+a_1(z-b)+\cdots\\
&=(z-b)^{-m}
\left[
a_{-m}+a_{-m+1}(z-b)+a_{-m+2}(z-b)^2+\cdots
\right]\\
&=(z-b)^{-m}\phi(z),
\qquad 0<|z-b|<R,
\end{aligned}
\tag{5.35}
$$

其中 $m$ 是正整数，

$$
\phi(z)=a_{-m}+a_{-m+1}(z-b)+a_{-m+2}(z-b)^2+\cdots
$$

在 $z=b$ 点的邻域内是解析的，而且 $\phi(b)=a_{-m}\ne0$。$b$ 点就称为 $f(z)$ 的 $m$ 阶极点。显然，只要 $|z-b|$ 足够小，$|f(z)|$ 可以大于任何正数，

$$
\lim_{z\to b}f(z)=\infty.
\tag{5.36}
$$

所以，函数在极点处的极限值是 $\infty$，即函数在极点附近无界。

反之，如果 $b$ 是 $f(z)$ 的孤立奇点，且 $\lim\limits_{z\to b}f(z)=\infty$，则 $b$ 是 $f(z)$ 的极点。

从 (5.35) 式可以看到，如果 $z=b$ 是 $f(z)$ 的 $m$ 阶极点，则

$$
\frac1{f(z)}
=(z-b)^m\frac1{\phi(z)},
\qquad
\frac1{\phi(z)}
\text{ 在 }z=b\text{ 点解析},
$$

所以 $z=b$ 必定是 $1/f(z)$ 的 $m$ 阶零点。反之亦然。利用这个关系，可以帮助我们寻找极点，尤其是确定极点的阶数。例如，$z=n\pi$ 是 $1/\sin z$ 的一阶极点；$z=2k\pi i,\ k=0,\pm1,\pm2,\cdots$ 是 $1/(e^z-1)$ 的一阶极点；$z=1$ 是 $1/(z-1)^2$ 的二阶极点。

**练习 5.1** 设 $f(z)$ 和 $g(z)$ 均在 $z=z_0$ 点解析，且 $f(z_0)=g(z_0)=0$，证明 l’Hôpital 法则成立，即

$$
\lim_{z\to z_0}\frac{f(z)}{g(z)}
=
\lim_{z\to z_0}\frac{f'(z)}{g'(z)}.
$$

**练习 5.2** 设 $f(z)$ 和 $g(z)$ 均在 $z=z_0$ 点的空心邻域内解析，且 $z=z_0$ 是 $f(z)$ 的 $n$ 阶极点，是 $g(z)$ 的 $m$ 阶极点，$n\le m$，证明 l’Hôpital 法则成立，即

$$
\lim_{z\to z_0}\frac{f(z)}{g(z)}
=
\lim_{z\to z_0}\frac{f'(z)}{g'(z)}.
$$

**本性奇点** 函数在本性奇点邻域内的 Laurent 展开具有无穷多个负幂项。

如果 $z=b$ 是函数 $f(z)$ 的本性奇点，则当 $z\to b$ 时，$f(z)$ 的极限不存在。更准确地说，$z\to b$ 的方式不同，$f(z)$ 可以逼近不同的数值。例如，$z=0$ 是函数

$$
e^{1/z}
=\sum_{n=0}^{\infty}\frac1{n!}\left(\frac1z\right)^n,
\qquad 0<|z|<\infty
$$

的本性奇点。当 $z$ 以不同方式趋于 0 时，就有不同的结果：

- 当 $z$ 沿正实轴趋于 0 时，$e^{1/z}\to\infty$；
- 当 $z$ 沿负实轴趋于 0 时，$e^{1/z}\to0$；
- 当 $z$ 沿虚轴趋于 0 时，$e^{1/z}$ 不趋于一个确定的数。

事实上，可以证明：对于本性奇点 $z=b$ 来说，任意给定一个数 $A$（有限或 $\infty$），总可以找到一个极限为 $b$ 的序列 $\{z_n\}$，当 $z_n\to b$ 时，使得 $f(z_n)\to A$。Picard 定理告诉我们，如果函数 $f(z)$ 只有一个本性奇点而别无其他奇点，则在本性奇点的任意一个小邻域内，$f(z)$ 可以取（并且无穷多次取）任意的有限复数值，至多可能有一个例外（称为 Picard 例外值）。例如 $e^{1/z}$ 在 $z=0$ 的空心邻域内的 Picard 例外值就是 0。

**练习 5.3** 设计几个序列 $\{z_n\}$，使当 $z\to0$ 时，$e^{1/z}$ 分别趋于 $\pm1$ 和 $\pm i$。

**无穷远点** 仍然通过变换 $z=1/t$ 把 $f(z)$ 化成 $f(1/t)$ 来讨论。若 $t=0$ 是 $f(1/t)$ 的可去奇点，则 $z=\infty$ 是 $f(z)$ 的可去奇点（甚至看成是解析点，就实际函数而言，我们往往可以不必苛刻地加以区分）；若 $t=0$ 是 $f(1/t)$ 的极点，则 $z=\infty$ 是 $f(z)$ 的极点；若 $t=0$ 是 $f(1/t)$ 的本性奇点，则 $z=\infty$ 是 $f(z)$ 的本性奇点。例如，$z=\infty$ 是 $1/(1+z)$ 的可去奇点（或解析点）；$z=\infty$ 是 $1+z^2$ 的二阶极点；$z=\infty$ 是 $e^z,\sin z,\cos z,\cdots$ 的本性奇点。

**练习 5.4** 试就 $z=\infty$ 点是 $f(z)$ 的可去奇点、极点、本性奇点三种情形，分别讨论 $f(z)$ 在 $z=\infty$ 点邻域内 Laurent 展开的特点。

## §5.7 解析延拓

先介绍一个例子。

几何级数

$$
\sum_{n=0}^{\infty}z^n=1+z+z^2+\cdots
\tag{5.37}
$$

在单位圆 $g_1:|z|<1$ 内收敛，代表一个解析函数（记为 $f_1(z)$）。事实上，例 4.1 中已经求出

$$
\sum_{n=0}^{\infty}z^n=\frac{1}{1-z},
\qquad |z|<1.
\tag{5.38}
$$

然而，考察函数 $1/(1-z)$ 本身可以发现，它在全复平面（并不限于单位圆内）都有定义，而且在 $z\ne1$ 的全复平面上处处解析。这样，通过幂级数 (5.37) 的求和，我们不单得出了幂级数在收敛圆内所代表的解析函数，而且得到的解析函数本身，还在更大的范围内解析。从这个例子中，可以抽象出解析延拓的概念。

**定义** 设函数 $f_1(z)$ 在区域 $g_1$ 内解析，函数 $f_2(z)$ 在区域 $g_2$ 内解析，$g_1\cap g_2\ne\varnothing$，而在 $g_1$ 与 $g_2$ 的公共区域 $g_1\cap g_2$ 内，$f_1(z)\equiv f_2(z)$，则称 $f_2(z)$ 为 $f_1(z)$ 在 $g_2$ 内的解析延拓；反之，$f_1(z)$ 是 $f_2(z)$ 在 $g_1$ 内的解析延拓。

对照上面的例子，$f_1(z)$ 就是幂级数 (5.37)，它在区域 $g_1$（单位圆）内解析；$f_2(z)$ 就是函数 $1/(1-z)$，它在区域 $g_2$（$z\ne1$ 的全复平面）上处处解析；而且，在 $g_1$ 和 $g_2$ 的公共区域 $g_1\cap g_2$ 内，$f_1(z)\equiv f_2(z)$。所以，从概念上说，函数 $1/(1-z)$ 就是幂级数 (5.37) 在全复平面（$z\ne1$）上的解析延拓。

求出幂级数的和函数，当然是理想的解析延拓方法。但是在许多问题中，往往难以求出和函数的有限表达式。即使如此，也仍然可以作解析延拓。仍以幂级数 (5.37) 为例。因为这个幂级数在单位圆内解析（记为 $f_1(z)$），所以能求出 $f_1(z)$ 在单位圆内任意一点的函数值及各阶导数值。例如，在 $z=i/2$ 点，有

$$
f_1\left(\frac i2\right)
=1+\frac i2+\left(\frac i2\right)^2+\cdots,
$$

$$
f_1'\left(\frac i2\right)
=1+2\cdot\frac i2+3\cdot\left(\frac i2\right)^2+\cdots,
$$

$$
\vdots
$$

因此，$f_1(z)$ 在 $z=i/2$ 点的 Taylor 展开是

$$
\sum_{n=0}^{\infty}
\frac{1}{n!}
f_1^{(n)}\left(\frac i2\right)
\left(z-\frac i2\right)^n.
\tag{5.39}
$$

这个级数当然也在它的收敛圆 $g_2:|z-i/2|<r$ 内收敛，也代表了一个解析函数，记为 $f_2(z)$。显然，在 $g_1$ 与 $g_2$ 的公共区域 $g_1\cap g_2(\ne\varnothing)$ 内，$f_1(z)\equiv f_2(z)$。这样得到的 $f_2(z)$ 就是 $f_1(z)$ 在区域 $g_2$ 内的解析延拓；反之，$f_1(z)$ 也是 $f_2(z)$ 在区域 $g_1$ 内的解析延拓。事实上，我们知道，

$$
f_1^{(n)}\left(\frac i2\right)
=
\frac{n!}{(1-i/2)^{n+1}},
\tag{5.40}
$$

所以，在区域 $g_2:|z-i/2|<\sqrt5/2$ 内，有

$$
f_2(z)
=\sum_{n=0}^{\infty}
\frac{1}{(1-i/2)^{n+1}}
\left(z-\frac i2\right)^n
=\frac{1}{1-z}.
\tag{5.41}
$$

因此，$f_1(z)$ 和 $f_2(z)$ 只不过是同一个函数（即 $1/(1-z)$）在不同区域内的幂级数表达式。这两个幂级数表达式都有各自的有效范围：$g_1:|z|<1$ 和 $g_2:|z-i/2|<\sqrt5/2$。同时 $g_1$ 和 $g_2$ 也有公共区域 $g_1\cap g_2\ne\varnothing$。在公共区域内，$f_1(z)\equiv f_2(z)$。

![Figure 5.2](...)

图 5.2 解析延拓

这样，我们从定义在一定区域 $g_1$ 内的幂级数出发，有可能得到在另一区域 $g_2$ 内的另一幂级数表达式，在两个区域的公共部分 $g_1\cap g_2(\ne\varnothing)$ 内二者相等。重复这个步骤，就有可能超出原来的定义范围，甚至可能扩展到整个 $z$ 复平面。

解析延拓的办法，可以用来扩大函数的定义域和解析范围。例如，用级数或积分定义的函数，本来都可能只适用于复平面上的某一范围，但是，通过解析延拓，就可能定义出在更大范围内的解析函数。在上面的例子中，就得到了在区域 $g_1\cup g_2$ 内解析的函数 $f(z)$，

$$
f(z)=
\begin{cases}
f_1(z), & z\in g_1,\\
f_2(z), & z\in g_2.
\end{cases}
\tag{5.42}
$$

一个实际例子是 $\Gamma$ 函数（见 §7.1）。由常用的积分定义（只在右半平面解析）出发，经过解析延拓，就得到它在全复平面上的定义。

在作解析函数的 Taylor 展开或 Laurent 展开时，也常常用到解析延拓。我们只要在级数的收敛区域（圆域或环域）内一个子区域乃至一段弧上定出展开系数即可。例如，例 5.7 中介绍的级数除法，实际上就需要用到解析延拓。

另一类要用到解析延拓的问题是常微分方程的求解问题。例如对于二阶常微分方程

$$
\frac{d^2w}{dz^2}
+p(z)\frac{dw}{dz}
+q(z)w=0,
\tag{5.43}
$$

通常只能求得在一定范围内的解式，而通过解析延拓，可以从这个解式推算出在其他范围内的表达式。§9.2 中有相关的讨论。

解析延拓是复变函数理论中最重要的概念之一。本节只是浅显地介绍了一下解析延拓的概念，并没有涉及解析延拓的理论问题。例如，解析延拓能否实现，这取决于函数奇点的分布状况：如果在收敛圆 $g_1$ 的边界上“布满了”函数 $f(z)$ 的奇点（例如，见本章习题 10），即在收敛圆的圆周上任意一点，其任意小的邻域内都有 $f(z)$ 的奇点，而且有无穷多个奇点，那么，在 $g_1$ 内重新作 Taylor 展开，其收敛范围绝不可能超出 $g_1$。再例如，解析延拓的结果是否与路径有关，即沿着不同路径延拓（到同一区域）的结果是否相同，或者说，沿着闭合路径一周回到原处时，延拓的结果是否还原。有兴趣的读者可以参看有关的专著。

## *§5.8 Bernoulli 数和 Euler 数

先讨论函数

$$
f(z)=\frac{z}{e^z-1}
$$

在 $z=0$ 点的幂级数展开。由于此函数分母的零点为 $2n\pi i,\ n=0,\pm1,\pm2,\cdots$，且均为一阶零点；同时 $z=0$ 又是函数分子的一阶零点，因此 $z=0$ 点是 $f(z)$ 的可去奇点，

$$
\lim_{z\to0}\frac{z}{e^z-1}=1.
$$

$f(z)$ 可在圆域 $|z|<2\pi$ 内作 Taylor 展开，

$$
\frac{z}{e^z-1}
=\sum_{n=0}^{\infty}\frac{B_n}{n!}z^n,
\tag{5.44}
$$

其中 $B_n$ 称为 Bernoulli 数。下面就来求 $B_n$。因为

$$
\frac{z}{e^z-1}
=
\frac z2\left(\frac{e^z+1}{e^z-1}-1\right)
=
\frac z2\frac{e^{z/2}+e^{-z/2}}{e^{z/2}-e^{-z/2}}
-\frac z2,
$$

第一项为 $z$ 的偶函数，展开时只有 $z$ 的偶次幂，所以

$$
B_{2n+1}=-\frac12\delta_{n0},
\qquad n=0,1,2,\cdots.
$$

再用待定系数法求出 $B_{2n}$。为此，将 (5.44) 式改写成（其中 $[k/2]$ 表示取 $k/2$ 的整数部分）

$$
\frac{e^z-1}{z}
\left[
-\frac z2+\sum_{n=0}^{\infty}\frac{B_{2n}}{(2n)!}z^{2n}
\right]
=
-\frac12\sum_{k=1}^{\infty}\frac{1}{k!}z^k
+
\sum_{k=0}^{\infty}
\left[
\sum_{n=0}^{[k/2]}
\frac{1}{(k-2n+1)!}\frac{B_{2n}}{(2n)!}
\right]z^{2k}
=1,
$$

就可以求得 $B_0=1$ 以及 $B_{2n}$ 的递推关系

$$
\sum_{n=0}^{[k/2]}
\frac{k!}{(k-2n+1)!}\frac{B_{2n}}{(2n)!}
=\frac12.
$$

从而推得

$$
B_2=\frac16,\qquad
B_4=-\frac1{30},\qquad
B_6=\frac1{42},\qquad
B_8=-\frac1{30},
$$

$$
B_{10}=\frac5{66},\qquad
B_{12}=-\frac{691}{2730},\qquad
B_{14}=\frac76,\qquad
B_{16}=-\frac{3617}{510},\qquad\cdots.
$$

不难看出 $B_{4n-2}>0,\ B_{4n}<0,\ n=1,2,3,\cdots$。

应用上面的结果还可以得到许多有用的展开式，例如

$$
\frac z2\cot\frac z2
=
\frac{iz}{2}\frac{e^{iz/2}+e^{-iz/2}}{e^{iz/2}-e^{-iz/2}}
=
\sum_{n=0}^{\infty}
(-1)^n\frac{B_{2n}}{(2n)!}z^{2n},
\qquad |z|<2\pi;
\tag{5.45}
$$

$$
\frac z2\tan\frac z2
=
\frac z2\cot\frac z2-z\cot z
=
\sum_{n=1}^{\infty}
(-1)^{n-1}
\frac{2^{2n}-1}{(2n)!}B_{2n}z^{2n},
\qquad |z|<\pi;
\tag{5.46}
$$

$$
z\csc z
=
\frac z2\cot\frac z2+\frac z2\tan\frac z2
=
\sum_{n=0}^{\infty}
(-1)^{n-1}
\frac{2(2^{2n-1}-1)}{(2n)!}B_{2n}z^{2n},
\qquad |z|<\pi;
\tag{5.47}
$$

$$
\ln\frac{\sin z}{z}
=
\sum_{n=1}^{\infty}
(-1)^n
\frac{2^{2n-1}}{n(2n)!}B_{2n}z^{2n},
\qquad |z|<\pi;
\tag{5.48}
$$

$$
\ln\cos z
=
\sum_{n=1}^{\infty}
(-1)^n
\frac{2^{2n-1}(2^{2n}-1)}{n(2n)!}B_{2n}z^{2n},
\qquad |z|<\frac{\pi}{2};
\tag{5.49}
$$

$$
\ln\frac{\tan z}{z}
=
\sum_{n=1}^{\infty}
(-1)^{n-1}
\frac{2^{2n}(2^{2n-1}-1)}{n(2n)!}B_{2n}z^{2n},
\qquad |z|<\frac{\pi}{2}.
\tag{5.50}
$$

用同样的办法还可以得到

$$
\frac{2e^{z/2}}{e^z+1}
=\operatorname{sech}\frac z2
=\sum_{n=0}^{\infty}\frac{E_n}{n!}\left(\frac z2\right)^n,
\qquad |z|<\pi,
\tag{5.51}
$$

其中 $E_n$ 称为 Euler 数，

$$
E_{2n+1}=0,
\qquad n=0,1,2,\cdots;
$$

$$
E_0=1,\qquad
E_2=-1,\qquad
E_4=5,\qquad
E_6=-61,\qquad\cdots.
$$

它们满足递推关系

$$
\sum_{l=0}^{k}
\frac{(2k)!}{(2l)!(2k-2l)!}E_{2l}
=0,
\qquad k\ge1.
$$

而 (5.51) 式还可以写成

$$
\frac{2e^{z/2}}{e^z+1}
=\operatorname{sech}\frac z2
=
\sum_{n=0}^{\infty}
\frac{E_{2n}}{(2n)!}
\left(\frac z2\right)^{2n},
\qquad |z|<\pi.
\tag{5.51'}
$$

**练习 5.5** 证明：

$$
\sec z
=\sum_{n=0}^{\infty}
(-1)^n\frac{E_{2n}}{(2n)!}z^{2n},
\qquad |z|<\frac{\pi}{2}.
$$

## 习题

1. 将下列函数在指定点展开为 Taylor 级数，并给出其收敛半径：

   (1) $1-z^2$，在 $z=1$ 展开；

   (2) $\sin z$，在 $z=n\pi$ 展开；

   (3) $\dfrac{1}{1+z+z^2}$，在 $z=0$ 展开；

   (4) $\dfrac{\sin z}{1-z}$，在 $z=0$ 展开；

   (5) $\exp\left\{\dfrac{1}{1-z}\right\}$，在 $z=0$ 展开（可只求前四项）。

2. 将下列函数在指定点展开为 Taylor 级数，并给出其收敛半径：

   (1) $\ln z$，在 $z=i$ 展开，规定 $0\le\arg z<2\pi$；

   (2) $\ln z$，在 $z=i$ 展开，规定 $\left.\ln z\right|_{z=i}=-\dfrac32\pi i$；

   (3) $\arctan z$ 的主值，在 $z=0$ 展开；

   (4) $\ln\dfrac{1+z}{1-z}$，在 $z=\infty$ 展开，规定
   $$
   \left.\ln\frac{1+z}{1-z}\right|_{z=\infty}=(2k+1)\pi i.
   $$

3. 求下列无穷级数之和，注意给出相应的收敛区域：

   (1)
   $$
   \sum_{n=0}^{\infty}\frac{1}{2n+1}z^{2n+1};
   $$

   (2)
   $$
   \sum_{n=0}^{\infty}\frac{1}{(2n)!}z^{2n};
   $$

   (3)
   $$
   \sum_{n=0}^{\infty}\sum_{m=0}^{\infty}
   \frac{(n+m)!}{n!\,m!}\left(\frac z2\right)^{n+m};
   $$

   (4)
   $$
   \sum_{n=0}^{\infty}\sum_{m=0}^{\infty}\sum_{p=0}^{\infty}
   \frac{(n+m+p)!}{n!\,m!\,p!}
   \left(\frac z3\right)^{n+m+p}.
   $$

4. 求下列函数的 Laurent 展开：

   (1) $\dfrac{1}{z^2(z-1)}$，在 $z=1$ 附近展开；

   (2) $\dfrac{1}{z^2(z-1)}$，展开区域为 $1<|z|<\infty$；

   (3) $\dfrac{1}{z^2-3z+2}$，展开区域为 $1<|z|<2$；

   (4) $\dfrac{1}{z^2-3z+2}$，展开区域为 $2<|z|<\infty$；

   (5) $\dfrac{(z-1)(z-2)}{(z-3)(z-4)}$，展开区域为 $3<|z|<4$；

   (6) $\dfrac{(z-1)(z-2)}{(z-3)(z-4)}$，展开区域为 $4<|z|<\infty$。

5. 用级数相乘的方法求下列函数（取主值分支）在 $z=0$ 点附近的级数展开：

   (1) $-\ln(1-z)\ln(1+z)$；

   (2) $\ln(1+z^2)\arctan z$。

6. 判断下列函数孤立奇点的性质，如果是极点，确定其阶数：

   (1) $\dfrac{1}{z^2+a^2}$，$a\ne0$；

   (2) $\dfrac{\cos az}{z^2}$；

   (3) $\dfrac{\cos az-\cos bz}{z^2}$，$a^2\ne b^2$；

   (4) $\dfrac{\sin z}{z^2}-\dfrac1z$；

   (5) $\cos\dfrac{1}{\sqrt z}$；

   (6) $\dfrac{\sqrt z}{\sin\sqrt z}$；

   (7) $\dfrac{1}{(z-1)\ln z}$；

   (8) $\displaystyle \int_0^z\frac{\sinh\sqrt\zeta}{\sqrt\zeta}\,d\zeta$。

7. 判断下列函数在 $\infty$ 点的性质：

   (1) $z^2$；

   (2) $\dfrac1z$；

   (3) $\dfrac{\cos z}{z}$；

   (4) $\dfrac{z}{\cos z}$；

   (5) $\dfrac{z^2+1}{e^z}$；

   (6) $\exp\left\{-\dfrac1{z^2}\right\}$；

   (7) $\dfrac{1}{\cosh\sqrt z}$；

   (8) $\sqrt{(z-1)(z-2)}$。

8. 定义在不同区域内的两个级数可以互为解析延拓。作为一个例子，证明

$$
f_1(z)=1+az+a^2z^2+a^3z^3+\cdots
$$

与

$$
f_2(z)
=\frac{1}{1-z}
-\frac{(1-a)z}{(1-z)^2}
+\frac{(1-a)^2z^2}{(1-z)^3}
-\cdots
$$

互为解析延拓。

9. 无穷级数在不同区域内可以收敛到不同的和函数。这两个和函数尽管（在不同区域内）有相同形式的级数表达式，但却不互为解析延拓。作为一个例子，证明级数

$$
\sum_{n=1}^{\infty}
\left(
\frac{1}{1-z^{2n+1}}
-\frac{1}{1-z^{2n}}
\right)
$$

在区域 $|z|<1$ 与 $|z|>1$ 内分别代表两个解析函数，但不互为解析延拓。

10. 已知级数

$$
\sum_{n=0}^{\infty}z^{2^n}
=z+z^2+z^4+z^8+z^{16}+\cdots
$$

在单位圆内收敛，其和函数为 $f(z)$。

   (1) 证明：$z=1$ 点是 $f(z)$ 的奇点。

   (2) 证明：$f(z)=z+f(z^2)$，因此，$z^2=1$ 的根都是 $f(z)$ 的奇点。

   (3) 类似地证明：$z^{2^k}=1$ 的 $2^k$ 个根都是 $f(z)$ 的奇点，其中 $k$ 为任意正整数。

   (4) 由此证明：不可能将 $f(z)$ 延拓到单位圆外。
