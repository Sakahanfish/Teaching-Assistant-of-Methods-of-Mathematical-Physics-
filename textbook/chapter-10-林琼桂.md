# 第十章 二阶线性常微分方程的级数解法和一般本征值问题\*

## 目录

§1 常点邻域的级数解法　196  
§2 Legendre 方程及其本征值问题　197  
　§2.1 Legendre 方程的级数解　197  
　§2.2 Legendre 方程的本征值问题　199  
§3 正则奇点邻域的级数解法　200  
§4 Bessel 方程　204  
　§4.1 Bessel 方程的级数解　204  
　§4.2 半奇数阶 Bessel 方程　206  
　§4.3 整数阶 Bessel 方程　207  
　§4.4 Neumann 函数的一般定义　207  
　§4.5 \*Neumann 函数的常规级数解法　208  
　§4.6 \*Frobenius 级数解法　210  
　§4.7 \*一般 Neumann 函数的整数阶极限　213  
§5 Sturm–Liouville 本征值问题　214  
　§5.1 SL 本征值问题的一般提法　214  
　§5.2 SL 本征值问题的一般结论　216  
补充习题　220

对偏微分方程分离变量后，马上需要解决的就是常微分方程及其本征值问题的求解。本书遇到的都是二阶线性常微分方程，因为它们来源于二阶线性偏微分方程。虽然常微分方程比偏微分方程简单，但也并不存在什么普遍有效的解析求解的法则。我们知道，一阶线性常微分方程的解可以用系数和非齐次项的积分表出，尽管这些积分不一定能积出来（即其原函数不一定是初等函数）。但对于二阶线性常微分方程，并不存在类似的结果。常系数方程和少数特殊类型的方程（比如 Euler 方程）可以用初等函数求解，另一些方程可以用级数解法或积分解法求解。级数解法可以算是比较系统的一种方法，因为对于那些能够用初等函数求解的简单情况，级数解法通常也有效。积分解法是积分变换的推广，本书不作介绍。应该指出，能够用级数解法或积分解法求解的方程是非常有限的。更多的时候人们只能采用数值解法。

## §1 常点邻域的级数解法

二阶线性齐次常微分方程的一般形式是

$$
y''(x)+p(x)y'(x)+q(x)y(x)=0, \tag{1}
$$

对于物理和工程问题中导出的微分方程，$x$ 通常是实数，$p(x)$ 和 $q(x)$ 是已知函数，$y(x)$ 是未知函数，它们的函数值也都是实数。为了应用复变函数理论来研究微分方程的解，可以把 $x$ 看作复数，并仍记作 $x$，相应地，$p(x)$、$q(x)$ 和 $y(x)$ 就成为复变函数，它们在实轴上取相应的实函数值。方程 (1) 可以附加初始条件

$$
y(x_0)=c_0,\qquad y'(x_0)=c_1. \tag{2}
$$

如果不附加初始条件，则通解中含有两个任意常数。

显然，方程 (1) 的解的行为取决于系数的行为。我们假定在复平面的某区域 $D$ 内，$p(x)$ 和 $q(x)$ 除有限个孤立奇点外是单值解析的。级数解法就是在 $D$ 内某点 $x_0$ 的邻域或去心邻域内将 $y(x)$ 展开为幂级数，即 Taylor 级数、Laurent 级数或更一般的幂级数（见后）。展开式的形式取决于 $p(x)$ 和 $q(x)$ 在 $x_0$ 的性质。如果 $p(x)$ 和 $q(x)$ 在 $x_0$ 解析，则 $x_0$ 称为方程的常点 (ordinary point)。如果 $x_0$ 是 $p(x)$ 和（或）$q(x)$ 的极点或本性奇点，它也就称为方程的奇点。

本节研究常点的级数解法，其理论基础是下面的

**定理**　如果 $p(x)$ 和 $q(x)$ 在圆 $|x-x_0|<R$ 内解析，则在该圆内满足方程 (1) 和初始条件 (2) 的解是存在、唯一而且解析的。

定理的大意是，如果系数是解析的，则方程的解也是解析的。这一结论非常直观，但证明起来却并不容易，所以我们不去探究定理的证明，而是把注意力集中在计算方法上。

既然 $p(x)$、$q(x)$ 和 $y(x)$ 都在圆内解析，那么就可以展开为 Taylor 级数：

$$
p(x)=\sum_{k=0}^{\infty}p_k(x-x_0)^k,\qquad
q(x)=\sum_{k=0}^{\infty}q_k(x-x_0)^k,\qquad
y(x)=\sum_{k=0}^{\infty}a_k(x-x_0)^k. \tag{3}
$$

其中的展开系数 $p_k$ 和 $q_k$ 是已知的，而 $a_k$ 是未知的。将这些展开式代入方程 (1)，合并同幂项，将左边整理成一个幂级数，由于右边为零，故所有 $(x-x_0)^k$ 的系数均必须为零，由此可得 $a_k$ 间的一系列代数方程。求解这些代数方程即可用 $a_0$ 和 $a_1$ 表出 $a_2,a_3,\cdots$，从而得到级数解。容易看出，$a_0=c_0$，$a_1=c_1$。如果不给定初始条件，则级数解中含有两个任意常数 $a_0$ 和 $a_1$，所以是方程 (1) 的通解。

下面补充讨论两个有关问题。它们与级数解法无关，也与常点或奇点无关。

首先，如果我们已经求得方程 (1) 的一个解 $y_1(x)$（不管用什么方法），则第二解就可以用积分表出。事实上，令 $y_2(x)=C(x)y_1(x)$，其中 $C(x)$ 是未知函数。代入方程 (1)，容易得到 $y_1C''+(2y_1'+py_1)C'=0$，这是 $C'(x)$ 的一阶线性方程，容易求出 $C'(x)$，再积分一次即得 $C(x)$，最后得到

$$
y_2(x)=y_1(x)\int^x\frac{1}{y_1^2(u)}\exp\left(-\int^u p(v)\,dv\right)du. \tag{4}
$$

这里包含两次不定积分，所以结果中有两个任意常数，因而已是方程 (1) 的通解。如果采用固定下限，则得到的是第二解。后面会用到这个结果。

其次，考虑非齐次方程

$$
y''(x)+p(x)y'(x)+q(x)y(x)=f(x). \tag{5}
$$

如果已经求得相应的齐次方程 (1) 的两个线性无关解 $y_1(x)$ 和 $y_2(x)$（不管用什么方法），则非齐次方程的一个特解 $Y(x)$ 也可以用积分表出。事实上，令 $Y(x)=C_1(x)y_1(x)+C_2(x)y_2(x)$，其中 $C_1(x)$ 和 $C_2(x)$ 是未知函数，满足附加条件

$$
y_1C_1'+y_2C_2'=0, \tag{6a}
$$

代入非齐次方程 (5)，利用附加条件以及 $y_1(x)$ 和 $y_2(x)$ 满足齐次方程的事实，易得

$$
y_1'C_1'+y_2'C_2'=f. \tag{6b}
$$

由于 $y_1(x)$ 和 $y_2(x)$ 线性无关，故 $\Delta\equiv y_1y_2'-y_2y_1'\ne0$（否则可以证明 $y_1(x)\propto y_2(x)$，则 $y_1(x)$ 与 $y_2(x)$ 线性相关，矛盾）。于是可以解得 $C_1'=-fy_2/\Delta$，$C_2'=fy_1/\Delta$，积分即得 $C_1(x)$ 和 $C_2(x)$，最后得到

$$
Y(x)=y_2(x)\int^x\frac{f(u)y_1(u)}{\Delta(u)}\,du-y_1(x)\int^x\frac{f(u)y_2(u)}{\Delta(u)}\,du. \tag{7}
$$

这里包含两个不定积分，所以结果中有两个任意常数，因而已是非齐次方程 (5) 的通解。如果采用固定下限，则得到的是一个特解。

## §2 Legendre 方程及其本征值问题

### §2.1 Legendre 方程的级数解

现在考虑 Legendre 方程

$$
(1-x^2)y''-2xy'+\lambda y=0. \tag{8}
$$

我们无法找到这一方程的简单解法，所以只能考虑级数解。与标准形式 (1) 比较可见

$$
p(x)=-\frac{2x}{1-x^2},\qquad q(x)=\frac{\lambda}{1-x^2}. \tag{9}
$$

显然，$x=0$ 是常点。又容易看出，$p(x)$ 和 $q(x)$ 在复平面上只有两个奇点 $x=\pm1$，所以，它们在圆 $|x|<1$ 内解析。按上节定理，在该圆内方程的解是解析的，故可设

$$
y(x)=\sum_{k=0}^{\infty}a_kx^k. \tag{10}
$$

容易得到下列各式：

$$
xy'(x)=\sum_{k=0}^{\infty}ka_kx^k,\qquad
x^2y''(x)=\sum_{k=0}^{\infty}k(k-1)a_kx^k, \tag{11a}
$$

$$
y''(x)=\sum_{k=0}^{\infty}k(k-1)a_kx^{k-2}
=\sum_{k=2}^{\infty}k(k-1)a_kx^{k-2}
=\sum_{k=0}^{\infty}(k+2)(k+1)a_{k+2}x^k. \tag{11b}
$$

代入方程并整理得

$$
\sum_{k=0}^{\infty}\big[(k+2)(k+1)a_{k+2}-k(k+1)a_k+\lambda a_k\big]x^k=0. \tag{12}
$$

比较两边，即得递推关系

$$
a_{k+2}=\frac{k(k+1)-\lambda}{(k+2)(k+1)}a_k,\qquad k\in\mathbb N. \tag{13}
$$

由此递推关系，所有 $a_{2k}$（$k\in\mathbb N^+$）均可由 $a_0$ 确定，所有 $a_{2k+1}$（$k\in\mathbb N^+$）均可由 $a_1$ 确定，于是得到级数解

$$
y(x)=a_0y_0(x)+a_1y_1(x), \tag{14}
$$

其中

$$
y_0(x)=\sum_{k=0}^{\infty}c_{2k}x^{2k},\qquad
y_1(x)=\sum_{k=0}^{\infty}c_{2k+1}x^{2k+1}, \tag{15}
$$

而 $c_{2k}=a_{2k}/a_0$，$c_{2k+1}=a_{2k+1}/a_1$，它们都只是 $k$ 和 $\lambda$ 的函数，而与 $a_0$、$a_1$ 无关。显然，$y_0(x)$ 和 $y_1(x)$ 是线性独立的，而 $y(x)$ 就是方程 (8) 的通解。

由于

$$
\lim_{k\to\infty}\left|\frac{a_{k+2}}{a_k}\right|
=\lim_{k\to\infty}\left|\frac{k(k+1)-\lambda}{(k+2)(k+1)}\right|=1, \tag{16}
$$

所以两个级数解的收敛半径都是 1，如所期望。但可以证明（参看第三章 §1.5），$y_0(\pm1)=\infty$（这是简化的写法，表示 $y_0(x)$ 在 $x=\pm1$ 两点均发散），$y_1(\pm1)=\infty$。这一结果对于下面确定本征值问题的解非常重要。

令 $\lambda=\nu(\nu+1)$，则递推关系 (13) 可以写作

$$
a_{k+2}=\frac{(k-\nu)(k+\nu+1)}{(k+2)(k+1)}a_k,\qquad k\in\mathbb N. \tag{17}
$$

注意给定 $\lambda$，方程 $\lambda=\nu(\nu+1)$ 有两个解，记任何一个解为 $\nu$，则另一解为 $-\nu-1$。容易看出，上面的递推关系在变换 $\nu\to-\nu-1$ 下不变，以下其它结果亦然。所以取任何一个解代入，结果都是一样的。重复利用递推关系可得

$$
a_{2k}=\frac{2^{2k}}{(2k)!}\left(-\frac{\nu}{2}\right)_k\left(\frac{\nu+1}{2}\right)_k a_0,
\qquad
a_{2k+1}=\frac{2^{2k}}{(2k+1)!}\left(\frac{-\nu+1}{2}\right)_k\left(\frac{\nu+2}{2}\right)_k a_1, \tag{18}
$$

其中引入了记号

$$
(\alpha)_0=1,\qquad (\alpha)_k=\alpha(\alpha+1)\cdots(\alpha+k-1)=\frac{\Gamma(\alpha+k)}{\Gamma(\alpha)},\qquad k\in\mathbb N^+. \tag{19}
$$

于是，式 (15) 中两个解的显式为

$$
y_0(x)=\sum_{k=0}^{\infty}\frac{2^{2k}}{(2k)!}\left(-\frac{\nu}{2}\right)_k\left(\frac{\nu+1}{2}\right)_k x^{2k}, \tag{20a}
$$

$$
y_1(x)=\sum_{k=0}^{\infty}\frac{2^{2k}}{(2k+1)!}\left(\frac{-\nu+1}{2}\right)_k\left(\frac{\nu+2}{2}\right)_k x^{2k+1}. \tag{20b}
$$

下面的讨论并不需要用到这一显式，所以读者能否掌握它都无关紧要。

### §2.2 Legendre 方程的本征值问题

由上章的讨论知道，物理上要求 Legendre 方程的解满足自然边界条件

$$
|y(\pm1)|<\infty. \tag{21}
$$

一般情况下，上面得到的两个解均不满足这一条件，所以，唯一的出路是让它们中断为多项式。由递推关系 (13) 可以看出，只要 $\lambda$ 取值恰当，这是可能的。这样同时也就确定了本征值。

如果 $\lambda=2n(2n+1)$（$n\in\mathbb N$），则由递推关系 (13) 可以看出，$a_{2n+2}=a_{2n+4}=\cdots=0$，从而 $y_0(x)$ 中断为 $2n$ 次多项式，而另一个解 $y_1(x)$ 仍为无穷级数，不满足边界条件 (21)。

如果 $\lambda=(2n+1)(2n+2)$（$n\in\mathbb N$），则由递推关系 (13) 可以看出，$a_{2n+3}=a_{2n+5}=\cdots=0$，从而 $y_1(x)$ 中断为 $2n+1$ 次多项式，而另一个解 $y_0(x)$ 仍为无穷级数，不满足边界条件 (21)。

综合起来，如果 $\lambda=l(l+1)$（$l\in\mathbb N$），两个线性独立解中就有一个中断为 $l$ 次多项式，它当然满足边界条件 (21)，而另一个解仍为无穷级数，不满足边界条件。

适当选取 $a_0$（当 $l=2n$）或 $a_1$（当 $l=2n+1$），使多项式解的最高次幂系数为

$$
a_l=\frac{(2l)!}{2^l(l!)^2},\qquad l\in\mathbb N, \tag{22}
$$

这样得到的解称为 $l$ 次 Legendre 多项式，记作 $P_l(x)$。

将 $\lambda=l(l+1)$ 代入递推关系 (13)，并将它改写为

$$
a_{k-2}=\frac{k(k-1)}{(k-l-2)(k+l-1)}a_k, \tag{23}
$$

反复应用这一递推关系，可以归纳出一般系数

$$
a_{l-2k}=(-)^k\frac{(2l-2k)!}{2^l k!(l-k)!(l-2k)!}, \tag{24}
$$

然后可以用数学归纳法加以证明。显然，$k$ 的取值应该使得 $l-2k\ge0$，故 $k\le l/2$，于是

$$
P_l(x)=\sum_{k=0}^{[l/2]}(-)^k\frac{(2l-2k)!}{2^lk!(l-k)!(l-2k)!}x^{l-2k},\qquad l\in\mathbb N. \tag{25}
$$

这就是 Legendre 多项式的显式，它的各种性质将在下一章详细讨论。

总结起来，Legendre 方程 (8) 在自然边界条件 (21) 下的本征值是 $\lambda=l(l+1)$，相应的本征函数是 $l$ 次 Legendre 多项式 $P_l(x)$。

当 $\lambda=l(l+1)$，与 $P_l(x)$ 线性独立的另一个解可以取为剩下的一个无穷级数解，也可以取为后者与 $P_l(x)$ 的线性组合，这个解有标准的取法，记作 $Q_l(x)$，它在 $x=\pm1$ 处具有 $\ln(1\mp x)$ 的奇性，此处不作详细讨论。

读者可能想到的一个问题是：当 $\lambda\ne l(l+1)$ 时，虽然 $y_0(x)$ 和 $y_1(x)$ 均不满足自然边界条件 (21)，但是否存在其适当的线性组合可以满足呢？假定存在适当的系数 $a_0$ 和 $a_1$（不全为零）使得式 (14) 满足 $|y(\pm1)|<\infty$，那么 $a_0$ 和 $a_1$ 必定都不为零，否则与已知结论 $y_0(\pm1)=\infty$、$y_1(\pm1)=\infty$ 矛盾。注意到 $y_0(x)$ 是偶函数，而 $y_1(x)$ 是奇函数，就容易推出

$$
y_0(x)=\frac{y(x)+y(-x)}{2a_0},\qquad y_1(x)=\frac{y(x)-y(-x)}{2a_1}. \tag{26}
$$

于是得到 $|y_0(\pm1)|<\infty$，$|y_1(\pm1)|<\infty$，与已知结论矛盾。所以，不存在任何无穷级数解满足自然边界条件 (21)。也可以说，不存在任何异于 $l(l+1)$ 的本征值。

**习题**　(1) 对 Legendre 方程作自变量变换 $t=(1-x)/2$，将其变为对 $t$ 的微分方程。  
(2) 试令 $y(t)=\sum_{k=0}^{\infty}a_kt^k$，求出一个级数解，其中只有 $a_0$ 是任意的。  
(3) 当 $\lambda$ 取何值时，该级数解退化为多项式？取 $a_0=1$，写出该多项式。并思考其与正文中所得 Legendre 多项式是否一致，理由何在？

## §3 正则奇点邻域的级数解法

现在再看方程 (1)。我们已经假定在复平面的某区域 $D$ 内，$p(x)$ 和 $q(x)$ 除有限个孤立奇点外是单值解析的。所以，如果 $D$ 内某点 $x_0$ 是 $p(x)$ 和（或）$q(x)$ 的奇点，那就只能是极点或本性奇点（而不会是支点，至于可去奇点则可当作常点，因为此时系数的展开式实际上是 Taylor 级数），所以一般有

$$
p(x)=\sum_{k=-\infty}^{\infty}p_k(x-x_0)^k,\qquad
q(x)=\sum_{k=-\infty}^{\infty}q_k(x-x_0)^k,\qquad 0<|x-x_0|<R. \tag{27}
$$

当然，如果 $x_0$ 是极点，则上面 Laurent 展开式中只有有限个负幂项。可以证明，这时方程 (1) 的两个解具有下列形式：

$$
y_1(x)=\sum_{k=-\infty}^{\infty}a_k(x-x_0)^{k+s_1},\qquad
y_2(x)=\sum_{k=-\infty}^{\infty}b_k(x-x_0)^{k+s_2}+\beta y_1(x)\ln(x-x_0). \tag{28}
$$

上式中 $s_1$ 和 $s_2$ 通常不是整数（最一般情况下可以是复数），所以 $x_0$ 一般来说是解的支点。当 $s_1-s_2\in\mathbb Z$，可能 $\beta\ne0$，即第二解中可能出现对数函数，否则 $\beta=0$。将上面的解式代入方程 (1)，可以得到一系列递推关系，由这些递推关系原则上可以确定 $s_1$、$s_2$ 和系数 $a_k$、$b_k$、$\beta$ 等。但是，由于每个递推关系都涉及无穷多个系数，所以实际计算是困难的甚至是不可能的。

比较简单的一种情况是上面的解式中只有有限个负幂项，这时适当调整 $s_1$ 和 $s_2$，总可以将它们写成

$$
y_1(x)=\sum_{k=0}^{\infty}a_k(x-x_0)^{k+s_1},\qquad a_0\ne0, \tag{29a}
$$

$$
y_2(x)=\sum_{k=0}^{\infty}b_k(x-x_0)^{k+s_2}+\beta y_1(x)\ln(x-x_0),\qquad b_0\ne0\ \text{或}\ \beta\ne0. \tag{29b}
$$

这样的解称为正则解。方程 (1) 是否有正则解，有几个（一个或两个）正则解，显然取决于 $p(x)$ 和 $q(x)$ 在 $x_0$ 的性质。对此，我们有下列

**定理（Fuchs）**　方程 (1) 在 $x_0$ 处有两个正则解的充要条件是：$(x-x_0)p(x)$ 和 $(x-x_0)^2q(x)$ 在 $x_0$ 解析。

上述条件就是说，$p(x)$ 以 $x_0$ 为不高于一阶的极点，$q(x)$ 以 $x_0$ 为不高于二阶的极点，这样的奇点称为方程的正则奇点（regular singular point，有些早期文献也简称为 regular point，但容易引起误解）。所以，Fuchs 定理也可以叙述为：方程 (1) 在 $x_0$ 处有两个正则解的充要条件是：$x_0$ 是方程的正则奇点。

我们不去研究这一定理的证明，但补充以下几点：① $s_1$ 和 $s_2$ 称为正则奇点或正则解的指标，$\mathrm{Re}\,s_1\ge\mathrm{Re}\,s_2$，即第一解表示指标实部较大者，它总不包含对数函数。② 若 $s_1-s_2\ne0,1,2,\cdots$，则 $\beta=0$，即第二解必定不包含对数函数。③ 若 $s_1=s_2$，则 $\beta\ne0$，即第二解必定包含对数函数。④ 若 $s_1-s_2\in\mathbb N^+$，则第二解可能包含对数函数，也可能不包含。

将正则解的形式和 $p(x)$、$q(x)$ 的 Laurent 展开式代入方程，即可得到一系列递推关系，从而确定正则解中的系数和指标。

下面分析一下求解的过程，这可以帮助我们理解上面几点补充结论。为书写方便，下面令 $x_0=0$，这并不失一般性。

因 $x_0=0$ 是正则奇点，故

$$
p(x)=\sum_{k=0}^{\infty}p_kx^{k-1},\qquad q(x)=\sum_{k=0}^{\infty}q_kx^{k-2},\qquad 0<|x|<R. \tag{30}
$$

将方程 (1) 两边同乘以 $x^2$，得

$$
x^2y''(x)+xp(x)\cdot xy'(x)+x^2q(x)\cdot y(x)=0. \tag{31}
$$

采用第一解的形式

$$
y(x)=\sum_{k=0}^{\infty}a_kx^{k+s},\qquad a_0\ne0, \tag{32}
$$

易得

$$
x^2y''(x)=\sum_{k=0}^{\infty}(k+s)(k+s-1)a_kx^{k+s},
$$

$$
xp(x)\cdot xy'(x)=\sum_{l=0}^{\infty}p_lx^l\sum_{n=0}^{\infty}(n+s)a_nx^{n+s}
=\sum_{k=0}^{\infty}\left[\sum_{n=0}^k(n+s)p_{k-n}a_n\right]x^{k+s},
$$

$$
x^2q(x)\cdot y(x)=\sum_{l=0}^{\infty}q_lx^l\sum_{n=0}^{\infty}a_nx^{n+s}
=\sum_{k=0}^{\infty}\left[\sum_{n=0}^k q_{k-n}a_n\right]x^{k+s}.
$$

全部代入上式，整理即得递推关系

$$
(k+s)(k+s-1)a_k+\sum_{n=0}^k\big[(n+s)p_{k-n}+q_{k-n}\big]a_n=0. \tag{33}
$$

令 $k=0$，由于 $a_0\ne0$，故得

$$
s(s-1)+p_0s+q_0=0. \tag{34}
$$

这就是决定指标的方程，它有两个根，记实部较大的根为 $s_1$，较小的为 $s_2$。当 $k>0$，递推关系可以写为

$$
\big[(k+s)(k+s-1)+p_0(k+s)+q_0\big]a_k
=-\sum_{n=0}^{k-1}\big[(n+s)p_{k-n}+q_{k-n}\big]a_n. \tag{35}
$$

以 $s=s_1$ 代入，可由 $a_0$ 推出所有的系数，即得第一解。显然，$a_k$ 与 $a_0$ 成正比，可以写作 $a_k=a_0f(k,s_1)$，其中 $f(k,s_1)$ 是 $k$ 和 $s_1$ 的函数，当然还依赖于 $p_0,q_0,p_1,q_1,\cdots,p_k,q_k$，但与 $a_0$ 无关。因此第一解的形式为

$$
y_1(x)=\sum_{k=0}^{\infty}a_kx^{k+s_1}=a_0\sum_{k=0}^{\infty}f(k,s_1)x^{k+s_1},\qquad a_0\ne0,\quad f(0,s_1)=1. \tag{36}
$$

这是上面的补充结论①。若 $s_1-s_2\ne0,1,2,\cdots$，以 $s=s_2$ 代入，亦可由 $a_0$ 推出所有的系数，将系数 $a_k$ 改写为 $b_k$，即得第二解，其形式为

$$
y_2(x)=\sum_{k=0}^{\infty}b_kx^{k+s_2}=b_0\sum_{k=0}^{\infty}f(k,s_2)x^{k+s_2},\qquad b_0\ne0,\quad f(0,s_2)=1. \tag{37}
$$

这是上面的补充结论②。若 $s_1=s_2$，式 (36) 与 (37) 实质上相同，所以用上面的方法只能求得一个解，第二解需要用其它方法求出，下面将看到，它包含对数函数。这是上面的补充结论③。下面分析补充结论④。

当 $s_1-s_2=m\in\mathbb N^+$，在用递推关系计算第二解的系数时可能会遇到困难。按上面的做法，将第二解的系数写作 $b_k$，则递推关系 (35) 成为

$$
\big[(k+s_2)(k+s_2-1)+p_0(k+s_2)+q_0\big]b_k
=-\sum_{n=0}^{k-1}\big[(n+s_2)p_{k-n}+q_{k-n}\big]b_n, \tag{38}
$$

其中已代入 $s=s_2$。利用上式可以由 $b_0$ 推出 $b_1,b_2,\cdots,b_{m-1}$，但计算 $b_m$ 时，左边的系数成为

$$(m+s_2)(m+s_2-1)+p_0(m+s_2)+q_0=s_1(s_1-1)+p_0s_1+q_0=0,$$

这时需要分开两种情况来讨论。

如果这时式 (38) 右边不为零，则出现矛盾，这说明第二解不可能具有形式 (32)，而需要用其它方法求出，下面将看到，它包含对数函数。

如果这时式 (38) 右边也为零，则 $b_m$ 可以任意。取 $b_m=0$，即可求出第二解。显然，所有系数均与 $b_0$ 成正比，可表示为 $b_k=b_0g(k)$，故第二解的形式为

$$
y_2(x)=b_0\sum_{k=0}^{\infty}g(k)x^{k+s_2},\qquad b_0\ne0,\quad g(0)=1,\quad g(m)=0. \tag{39}
$$

当 $k\le m-1$，$g(k)=f(k,s_2)$，但 $k\ge m$ 后则不成立。

如果不取 $b_m=0$，则自 $b_{m+1}$ 以后，递推关系成为

$$
\begin{aligned}
&\big[(k+s_1)(k+s_1-1)+p_0(k+s_1)+q_0\big]b_{k+m}\\
&\quad=-\sum_{n=0}^{m-1}\big[(n+s_2)p_{k+m-n}+q_{k+m-n}\big]b_n
-\sum_{n=0}^{k-1}\big[(n+s_1)p_{k-n}+q_{k-n}\big]b_{m+n},\qquad k\ge1.
\end{aligned} \tag{40}
$$

显然，所有系数均为 $b_0$ 与 $b_m$ 的线性组合。注意到 $b_m=0$ 时已经有 $b_k=b_0g(k)$，则 $b_{m+k}=b_0g(m+k)+b_m\tilde f(k)$。当 $b_0=0$，上式成为 $b_{m+k}=b_m\tilde f(k)$。但 $b_0=0$ 时式 (40) 变成

$$
\big[(k+s_1)(k+s_1-1)+p_0(k+s_1)+q_0\big]b_{k+m}
=-\sum_{n=0}^{k-1}\big[(n+s_1)p_{k-n}+q_{k-n}\big]b_{m+n},\qquad k\ge1. \tag{41}
$$

即由 $b_m$ 递推 $b_{m+k}$ 的方程与由 $a_0$ 递推 $a_k$ 的方程完全一样，所以 $b_{m+k}=b_mf(k,s_1)$，从而 $\tilde f(k)=f(k,s_1)$。于是

$$
b_{m+k}=b_0g(m+k)+b_mf(k,s_1),\qquad k\ge0. \tag{42}
$$

注意上式对 $k=0$ 亦成立。故第二解为

$$
\begin{aligned}
y_2(x)&=\sum_{k=0}^{\infty}b_kx^{k+s_2}
=b_0\sum_{k=0}^{\infty}g(k)x^{k+s_2}+b_m\sum_{k=0}^{\infty}f(k,s_1)x^{k+m+s_2}\\
&=b_0\sum_{k=0}^{\infty}g(k)x^{k+s_2}+b_m\sum_{k=0}^{\infty}f(k,s_1)x^{k+s_1}.
\end{aligned} \tag{43}
$$

易见其中第二部分求和与第一解 $y_1(x)$ 成正比，所以上式已经是通解。因此，上面取 $b_m=0$ 而得第二解是恰当的。

最后讲一下式 (32) 失效时如何求出第二解。由于已知第一解总具有式 (36) 的形式，故可由式 (4) 求出第二解。下面作积分时，可以适当选取积分常数。由式 (30)，有

$$
\int^u p(v)\,dv=p_0\ln u+\sum_{k=0}^{\infty}r_ku^k,
$$

其中 $r_k=p_k/k$（$k\ne0$），$r_0$ 为积分常数，故

$$
\exp\left(-\int^u p(v)\,dv\right)=u^{-p_0}\sum_{k=0}^{\infty}c_ku^k,\qquad c_0\ne0.
$$

而

$$
y_1^2(u)=u^{2s_1}\sum_{k=0}^{\infty}d_ku^k,\qquad d_0\ne0,
$$

所以

$$
\frac{1}{y_1^2(u)}\exp\left(-\int^u p(v)\,dv\right)
=u^{-2s_1-p_0}\sum_{k=0}^{\infty}e_ku^k
=u^{-m-1}\sum_{k=0}^{\infty}e_ku^k
=\sum_{k=0}^{\infty}e_ku^{k-m-1},\qquad e_0\ne0,
$$

其中利用了 $s_1+s_2=-p_0+1$ 以及 $2s_1+p_0=s_1+s_2+m+p_0=m+1$。下面积分时要分开两种情况。

若 $m=0$，即两个指标相等，积分得

$$
\int^x\frac{1}{y_1^2(u)}\exp\left(-\int^u p(v)\,dv\right)du
=\beta\ln x+\sum_{k=0}^{\infty}f_kx^k,\qquad \beta=e_0\ne0,
$$

其中 $f_k=e_k/k$（$k>0$），$f_0$ 是积分常数，故

$$
y_2(x)=\sum_{k=0}^{\infty}b_kx^{k+s_1}+\beta y_1(x)\ln x,\qquad \beta\ne0, \tag{44}
$$

其中 $b_0$ 可以为零，但因为 $\beta\ne0$，第二解一定包含对数函数。

若 $m>0$，积分得

$$
\int^x\frac{1}{y_1^2(u)}\exp\left(-\int^u p(v)\,dv\right)du
=\beta\ln x+\sum_{k=0}^{\infty}f_kx^{k-m},\qquad f_0=-\frac{e_0}{m}\ne0,
$$

其中 $\beta=e_m$ 可能为零，也可能不为零，$f_k=e_k/(k-m)$（$k\ne m$），$f_m$ 为积分常数。故

$$
y_2(x)=\sum_{k=0}^{\infty}b_kx^{k+s_2}+\beta y_1(x)\ln x,\qquad b_0\ne0. \tag{45}
$$

这一解可能包含对数函数，也可能不包含。

## §4 Bessel 方程

### §4.1 Bessel 方程的级数解

上章已经看到，在柱坐标下对 Laplace 方程分离变量，会遇到 Bessel 方程及其本征值问题。对波动或热传导方程分离变量，也会遇到类似的问题。Bessel 方程也没有简单的解法，所以只能考虑级数解。

数学上，Bessel 方程的一般形式是

$$
y''+\frac1x y'+\left(1-\frac{\nu^2}{x^2}\right)y=0, \tag{46}
$$

其中 $\nu$ 称为 Bessel 方程的阶，它可以是复数。若将 $\nu$ 换为 $-\nu$，上式不变，故可设 $\mathrm{Re}\,\nu\ge0$ 而不失一般性。与标准形式 (1) 比较可见

$$
p(x)=\frac1x,\qquad q(x)=1-\frac{\nu^2}{x^2}. \tag{47}
$$

显然，$x=0$ 是方程的正则奇点。又容易看出，在复平面上方程没有其它奇点。在 $x=0$ 的去心邻域内，可设解的形式为

$$
y(x)=\sum_{k=0}^{\infty}a_kx^{k+s},\qquad a_0\ne0. \tag{48}
$$

将式 (46) 改写为

$$
Ly\equiv\left(x^2\frac{d^2}{dx^2}+x\frac d{dx}+x^2-\nu^2\right)y
=x^2y''+xy'+(x^2-\nu^2)y=0, \tag{49}
$$

其中引入算符 $L$ 只是为了后面书写方便。容易求出

$$
x^2y''+xy'=\sum_{k=0}^{\infty}(k+s)(k+s-1)a_kx^{k+s}+\sum_{k=0}^{\infty}(k+s)a_kx^{k+s}
=\sum_{k=0}^{\infty}(k+s)^2a_kx^{k+s}.
$$

于是

$$
\begin{aligned}
Ly&=\sum_{k=0}^{\infty}\big[(k+s)^2-\nu^2\big]a_kx^{k+s}+\sum_{k=0}^{\infty}a_kx^{k+s+2}\\
&=(s^2-\nu^2)a_0x^s+\big[(s+1)^2-\nu^2\big]a_1x^{s+1}
+\sum_{k=2}^{\infty}\big[(k+s)^2-\nu^2\big]a_kx^{k+s}+\sum_{k=0}^{\infty}a_kx^{k+s+2}\\
&=(s^2-\nu^2)a_0x^s+\big[(s+1)^2-\nu^2\big]a_1x^{s+1}
+\sum_{k=0}^{\infty}\left\{\big[(k+s+2)^2-\nu^2\big]a_{k+2}+a_k\right\}x^{k+s+2}.
\end{aligned} \tag{50}
$$

上式必须为零，故各项系数均为零。由于 $a_0\ne0$，即得 $s^2-\nu^2=0$，这是决定指标的方程，因 $\mathrm{Re}\,\nu\ge0$，故其解为

$$
s_1=\nu,\qquad s_2=-\nu. \tag{51}
$$

于是 $s_1-s_2=2\nu$，暂时假定 $\nu$ 不等于整数或半奇数，则 $s_1-s_2$ 不为整数，按上节的一般理论，两个解均不包含对数函数。$\nu$ 为整数或半奇数的情况将在随后各小节中加以讨论。

首先讨论第一解，对应于 $s_1=\nu$。容易得到 $a_1=0$ 和递推关系

$$
a_{k+2}=-\frac{a_k}{(k+\nu+2)^2-\nu^2}=-\frac{a_k}{(k+2\nu+2)(k+2)}. \tag{52}
$$

由此递推关系和 $a_1=0$ 可以推出

$$
a_{2k+1}=0,\qquad k\in\mathbb N. \tag{53}
$$

反复利用递推关系又可推出

$$
\begin{aligned}
a_{2k}&=-\frac{a_{2k-2}}{(2\nu+2k)\cdot2k}
=(-)^2\frac{a_{2k-4}}{(2\nu+2k)(2\nu+2k-2)\cdot2k(2k-2)}=\cdots\\
&=(-)^k\frac{a_0}{(2\nu+2k)(2\nu+2k-2)\cdots(2\nu+2)\cdot2k(2k-2)\cdots2}\\
&=(-)^k\frac{a_0}{2^{2k}(\nu+k)(\nu+k-1)\cdots(\nu+1)\cdot k!}\\
&=(-)^k\frac{a_0\Gamma(\nu+1)}{2^{2k}k!\Gamma(\nu+k+1)}.
\end{aligned} \tag{54}
$$

最后一步将分子分母同乘以 $\Gamma(\nu+1)$ 并利用了 $\Gamma$ 函数的性质 $\Gamma(z+1)=z\Gamma(z)$。于是得到第一解为

$$
y_1(x)=\sum_{k=0}^{\infty}(-)^k\frac{a_0\Gamma(\nu+1)}{2^{2k}k!\Gamma(k+\nu+1)}x^{2k+\nu}. \tag{55}
$$

取 $a_0=1/2^\nu\Gamma(\nu+1)$，这样得到的解称为 $\nu$ 阶 Bessel 函数，记作 $J_\nu(x)$，即 $y_1(x)=J_\nu(x)$，其形式为

$$
J_\nu(x)=\sum_{k=0}^{\infty}\frac{(-)^k}{k!\Gamma(k+\nu+1)}\left(\frac x2\right)^{2k+\nu}. \tag{56}
$$

应该强调的是，第一解对任何 $\nu$ 值都是适用的。

其次讨论第二解，对应于 $s_2=-\nu$。此时 $x^{s+1}=x^{-\nu+1}$ 项的系数为 $[(s_2+1)^2-\nu^2]a_1=[(-\nu+1)^2-\nu^2]a_1=(-2\nu+1)a_1$，它应该为零。由于现在 $\nu\ne1/2$，故 $a_1=0$ 仍成立。重复第一解的推导过程，将其中的 $\nu$ 换为 $-\nu$，即得第二解为 $y_2(x)=J_{-\nu}(x)$，后者仍由式 (56) 给出，只是将其中的 $\nu$ 换为 $-\nu$，只要 $\nu$ 不为整数，$J_{-\nu}(x)$ 的定义就是恰当的。

总结起来，当 $\nu$ 不等于整数或半奇数时，Bessel 方程的两个线性独立解为

$$
y(x)=\{J_\nu(x),J_{-\nu}(x)\}. \tag{57}
$$

显然，$x=0$ 是 $J_{\pm\nu}(x)$ 的支点，将 $J_\nu(x)$（或 $J_{-\nu}(x)$）中的 $x^\nu$（或 $x^{-\nu}$）因子提出来，剩下的因子是一个幂级数，由递推关系 (52) 容易看出，该幂级数的收敛半径为无穷大，故 $J_{\pm\nu}(x)$ 在沿 $x=0$ 至 $x=\infty$ 适当割破的 $x$ 平面上是单值解析的。

**习题\***（选做）重新考虑 §2 习题中的 Legendre 方程（作为自变量 $t=(1-x)/2$ 的微分方程）。试用 $y(t)=\sum_{k=0}^{\infty}a_kt^{k+s}$ 的形式求级数解，先确定 $s$ 的数值，然后求出级数解。用这种形式能求出两个解吗？

### §4.2 半奇数阶 Bessel 方程

本小节考虑 $\nu=(2l-1)/2$，$l\in\mathbb N^+$ 的情况。此时 $s_1-s_2=2l-1\in\mathbb N^+$，根据上节的一般理论，第二解可能会包含对数函数。但具体求解可以发现，第二解实际上是不包含对数函数的。

以 $\nu=1/2$ 为例。第一解当然就是 $y_1(x)=J_{1/2}(x)$。对第二解，$s_2=-\nu=-1/2$，此时由 $x^{s+1}=x^{-\nu+1}=x^{1/2}$ 项的系数为零得 $0=[(s_2+1)^2-\nu^2]a_1=[(-\nu+1)^2-\nu^2]a_1=(-2\nu+1)a_1=0\cdot a_1$。由此可见，$a_1$ 可以任取。那么取 $a_1=0$ 显然是最方便的。这样就可以象上一小节一样求得第二解为 $y_2(x)=J_{-1/2}(x)$。如果不取 $a_1=0$，那么求得的解将包含两个任意常数 $a_0$ 和 $a_1$，计算可以发现（或参考上节小字部分的一般讨论），包含 $a_1$ 的部分与 $J_{1/2}(x)$ 成正比，所以这个解已经是通解。由此可见，取 $a_1=0$ 而得第二解是恰当的。

当 $l\ge2$，即 $\nu\ge3/2$，第一解仍然是 $y_1(x)=J_\nu(x)$。对第二解，有 $a_1=0$，由递推关系 $(k+2)(k-2l+3)a_{k+2}=-a_k$ 可得 $a_3=a_5=\cdots=a_{2l-3}=0$，然后就是 $0\cdot a_{2l-1}=-a_{2l-3}=0$，于是 $a_{2l-1}$ 可以任取。如前取 $a_{2l-1}=0$，则所有 $a_{2k+1}=0$（$k\in\mathbb N$），这样也就可以象上一小节一样求得第二解为 $y_2(x)=J_{-\nu}(x)$。

总结起来，当 $\nu$ 为半奇数时，Bessel 方程的两个线性独立解仍由式 (57) 给出，只需将相应的 $\nu$ 值代入即可。

以后我们会证明，半奇数阶的 Bessel 函数都是初等函数。这里证明 $J_{\pm1/2}(x)$ 是初等函数。

将 $\nu=\pm1/2$ 代入式 (56) 并利用 $\Gamma(k+1/2)=\sqrt\pi(2k)!/2^{2k}k!$，$\Gamma(k+3/2)=\sqrt\pi(2k+1)!/2^{2k+1}k!$，整理即得

$$
J_{1/2}(x)=\sqrt{\frac{2}{\pi x}}\sum_{k=0}^{\infty}\frac{(-)^k}{(2k+1)!}x^{2k+1},\qquad
J_{-1/2}(x)=\sqrt{\frac{2}{\pi x}}\sum_{k=0}^{\infty}\frac{(-)^k}{(2k)!}x^{2k},
$$

即

$$
J_{1/2}(x)=\sqrt{\frac{2}{\pi x}}\sin x,\qquad J_{-1/2}(x)=\sqrt{\frac{2}{\pi x}}\cos x. \tag{58}
$$

### §4.3 整数阶 Bessel 方程

本小节考虑 $\nu=m\in\mathbb N$ 的情况。此时 $s_1-s_2=2m$ 为自然数或零。当 $m=0$ 时，$s_1=s_2=0$，我们只能求得一个形如 (48) 的解，即 $y_1(x)=J_0(x)$，根据上节的一般理论，第二解必定包含对数函数。当 $m\in\mathbb N^+$ 时，第一解为 $y_1(x)=J_m(x)$。对于第二解，将 $\nu=m$ 和 $s_2=-m$ 代入式 (50) 并令各项系数为零，即得 $a_1=0$ 和递推关系 $(k+2)(k-2m+2)a_{k+2}=-a_k$，由此即可推出所有 $a_{2k+1}=0$（$k\in\mathbb N$），并用 $a_0$ 表出 $a_{2k}$（$k=1,2,\cdots,m-1$），然而当 $k=2m-2$ 时，递推关系给出 $0\cdot a_{2m}=-a_{2m-2}$，但 $a_{2m-2}\propto a_0$ 不为零，所以导致矛盾，于是第二解不可能具有式 (48) 的形式，因而也包含对数函数。令

$$
y(x)=\beta J_m(x)\ln x+\sum_{k=0}^{\infty}a_kx^{k-m}, \tag{59}
$$

代入式 (49)（其中 $\nu=m$），适当选取其中可以任取的两个常数（$m=0$ 时是 $\beta$ 和 $a_0$，$m\in\mathbb N^+$ 时是 $\beta$ 和 $a_{2m}$，参看后面的求解）可以求得第二解 $y_2(x)=N_m(x)$，称为 $m$ 阶 Neumann 函数，其形式为

$$
\begin{aligned}
N_m(x)&=\frac{2}{\pi}J_m(x)\ln\frac x2-\frac1\pi\sum_{k=0}^{m-1}\frac{(m-k-1)!}{k!}\left(\frac x2\right)^{2k-m}\\
&\quad-\frac1\pi\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}[\psi(k+m+1)+\psi(k+1)]\left(\frac x2\right)^{2k+m}.
\end{aligned} \tag{60}
$$

当 $m=0$ 时，规定去掉其中第二项有限和。上式中出现的 $\psi$ 函数定义为 $\Gamma$ 函数的对数微商，即 $\psi(x)=d\ln\Gamma(x)/dx=\Gamma'(x)/\Gamma(x)$，这里我们只需要知道 $\psi(x+1)=\psi(x)+1/x$（来源于 $\Gamma(x+1)=x\Gamma(x)$）和 $\psi(1)=-\gamma$，而 $\gamma=0.577\,215\,664\,901\cdots$ 称为 Euler 常数，定义为

$$
\gamma=\lim_{n\to\infty}\left(\sum_{k=1}^{n}\frac1k-\ln n\right). \tag{61}
$$

在不同的书上，Neumann 函数的形式略有不同，但实质上是一样的。

总结起来，当 $\nu=m$ 时，Bessel 方程的两个线性独立解为

$$
y(x)=\{J_m(x),N_m(x)\}. \tag{62}
$$

物理上最常遇到的就是这种情况，注意 $N_m(x)$ 在 $x=0$ 处有奇性，所以，如果求解区域包括 $x=0$ 点，就应该舍弃 $N_m(x)$。

### §4.4 Neumann 函数的一般定义

对一般的 $\nu$ 值，我们可以定义 $\nu$ 阶 Neumann 函数为

$$
N_\nu(x)=\frac{J_\nu(x)\cos\nu\pi-J_{-\nu}(x)}{\sin\nu\pi}. \tag{63}
$$

当 $\nu\ne m\in\mathbb N$，它显然与 $J_\nu(x)$ 线性独立，故此时 Bessel 方程的两个线性独立解亦可取为 $y(x)=\{J_\nu(x),N_\nu(x)\}$。当 $\nu\to0$ 时，上式成为 $0/0$ 型。而当 $\nu\to m\in\mathbb N^+$ 时，$J_{-\nu}(x)$ 的表达式中，$k\le m-1$ 各项对求和没有贡献，这是因为对这样的 $k$ 值，$\Gamma(k-\nu+1)\to\infty$。因此，$J_{-\nu}(x)$ 中的求和实际上从 $k=m$ 开始，于是

$$
J_{-m}(x)\equiv\lim_{\nu\to m}J_{-\nu}(x)
=\sum_{k=m}^{\infty}\frac{(-)^k}{k!(k-m)!}\left(\frac x2\right)^{2k-m}
=(-)^m\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}\left(\frac x2\right)^{2k+m},
$$

其中最后一步作了变换 $k=k'+m$ 并在变换后将求和指标 $k'$ 重新换为 $k$，易见右边的和式正是 $J_m(x)$，于是

$$
J_{-m}(x)=(-)^mJ_m(x), \tag{64}
$$

可见当 $\nu\to m\in\mathbb N^+$ 时，式 (63) 亦成为 $0/0$ 型。用 L'Hospital 法则求出极限，可以发现，

$$
\lim_{\nu\to m}N_\nu(x)=N_m(x),\qquad m\in\mathbb N, \tag{65}
$$

其中右边由式 (60) 给出。证明可参看后面第七小节。

总结起来，对于任何 $\nu$ 值，Bessel 方程的两个线性独立解总可以取为

$$
y(x)=\{J_\nu(x),N_\nu(x)\}. \tag{66}
$$

容易看出，当 $\nu$ 为半奇数时，$N_\nu(x)$ 与 $J_{-\nu}(x)$ 只相差一个常数因子。

### §4.5 \*Neumann 函数的常规级数解法

首先考虑 $m=0$ 的情况，将 $y(x)=\beta J_0(x)\ln x+\sum_{k=0}^{\infty}a_kx^k$ 代入式 (49)（其中 $\nu=0$）可得

$$
a_1x+\sum_{k=0}^{\infty}\big[(k+2)^2a_{k+2}+a_k\big]x^{k+2}+2\beta xJ_0'(x)+\beta[x^2J_0''(x)+xJ_0'(x)+x^2J_0(x)]\ln x=0,
$$

因为 $J_0(x)$ 是式 (49)（其中 $\nu=0$）的解，故上式中最后一项为零，代入 $J_0(x)$ 的级数表示，上式可化为

$$
a_1x+\sum_{k=1}^{\infty}[(2k+1)^2a_{2k+1}+a_{2k-1}]x^{2k+1}
+\sum_{k=1}^{\infty}\left[(2k)^2a_{2k}+a_{2k-2}+\frac{(-)^k\beta}{2^{2k-2}(k-1)!k!}\right]x^{2k}=0.
$$

上式对 $a_0$ 没有限制，故 $a_0$ 可以任取。由上式第一项和第二项立得

$$
a_{2k+1}=0,\qquad k\in\mathbb N, \tag{67}
$$

而由第三项可得递推关系

$$
a_{2k}=-\frac{a_{2k-2}}{(2k)^2}-\frac{(-)^k\beta}{2^{2k}(k!)^2},\qquad k\in\mathbb N^+. \tag{68}
$$

由此可以看出，$\beta$ 也是任意常数。反复利用这一递推关系可得

$$
a_{2k}=\frac{(-)^ka_0}{2^{2k}(k!)^2}-\frac{(-)^k\beta}{2^{2k}(k!)^2}\sum_{r=1}^{k}\frac1r,\qquad k\in\mathbb N^+,
$$

利用 $\psi$ 函数，上式可以写作

$$
a_{2k}=\frac{(-)^ka_0}{2^{2k}(k!)^2}-\frac{(-)^k\beta}{2^{2k}(k!)^2}[\psi(k+1)-\psi(1)],\qquad k\in\mathbb N. \tag{69}
$$

注意最后的表达式对于 $k=0$ 也成立。于是得到级数解

$$
y(x)=\beta J_0(x)\ln x-\beta\sum_{k=0}^{\infty}\frac{(-)^k}{(k!)^2}[\psi(k+1)-\psi(1)]\left(\frac x2\right)^{2k}
+a_0\sum_{k=0}^{\infty}\frac{(-)^k}{(k!)^2}\left(\frac x2\right)^{2k}. \tag{70}
$$

容易看出，最后一项求和正是第一解 $J_0(x)$，因此上式已经是通解；同时我们看到，如果取 $\beta=0$，则只能得到第一解，所以，第二解一定包含对数函数。我们当然可以取 $a_0=0$ 和 $\beta=1$ 而得到第二解，但习惯上是取

$$
\beta=\frac2\pi,\qquad a_0=-\frac2\pi[\ln2+\psi(1)], \tag{71}
$$

这样得到的第二解就是 $y_2(x)=N_0(x)$，其形式已经在前面给出。这里系数的取法是为了使得结果与式 (63) 的极限一致。

其次考虑 $m\in\mathbb N^+$ 的情况，将式 (59) 代入式 (49)（其中 $\nu=m$）可得

$$
(1-2m)a_1x^{1-m}+\sum_{k=0}^{\infty}[(k+2)(k+2-2m)a_{k+2}+a_k]x^{k+2-m}
+2\beta xJ_m'(x)
+\beta[x^2J_m''(x)+xJ_m'(x)+(x^2-m^2)J_m(x)]\ln x=0,
$$

因为 $J_m(x)$ 是式 (49)（其中 $\nu=m$）的解，故上式中最后一项为零。将上式两边同乘以 $x^m$，代入 $J_m(x)$ 的级数表示，整理得

$$
\begin{aligned}
&(1-2m)a_1x+\sum_{k=1}^{\infty}[(2k+1)(2k+1-2m)a_{2k+1}+a_{2k-1}]x^{2k+1}\\
&\quad+\sum_{k=1}^{\infty}[(2k)(2k-2m)a_{2k}+a_{2k-2}]x^{2k}
+2\beta\sum_{k=0}^{\infty}\frac{(-)^k(2k+m)}{2^{2k+m}k!(k+m)!}x^{2k+2m}=0.
\end{aligned}
$$

将上式第三项拆成两部分，第一部分从 $k=1$ 到 $k=m-1$ 求和，第二部分是从 $k=m$ 开始的无穷级数，对第二部分作指标变换 $k=k'+m$，再将 $k'$ 改写为 $k$，并与最后一项求和合并，得到

$$
\begin{aligned}
&(1-2m)a_1x+\sum_{k=1}^{\infty}[(2k+1)(2k+1-2m)a_{2k+1}+a_{2k-1}]x^{2k+1}\\
&\quad+\sum_{k=1}^{m-1}[(2k)(2k-2m)a_{2k}+a_{2k-2}]x^{2k}\\
&\quad+\sum_{k=0}^{\infty}\left[(2k)(2k+2m)a_{2k+2m}+a_{2k+2m-2}+2\beta\frac{(-)^k(2k+m)}{2^{2k+m}k!(k+m)!}\right]x^{2k+2m}=0. \tag{72}
\end{aligned}
$$

上式对 $a_0$ 不构成限制，因此 $a_0$ 可以任取。由上式第一项和第二项立得

$$
a_{2k+1}=0,\qquad k\in\mathbb N. \tag{73}
$$

由第三项可得递推关系

$$
a_{2k}=\frac{a_{2k-2}}{2k(m-k)},\qquad k=1,2,\cdots,m-1, \tag{74}
$$

反复利用这一递推关系可得

$$
a_{2k}=\frac{(m-k-1)!}{2^{2k}k!(m-1)!}a_0,\qquad k=0,1,\cdots,m-1. \tag{75}
$$

注意最后的表达式对于 $k=0$ 也成立。由式 (72) 的最后一项可得（由 $k=0$ 项系数为零）

$$
\beta=-2^{m-1}(m-1)!a_{2m-2}=-\frac{a_0}{2^{m-1}(m-1)!}. \tag{76}
$$

和递推关系

$$
a_{2k+2m}=-\frac{a_{2k+2m-2}}{2^{2}k(k+m)}-\frac{(-)^k\beta}{2^{2k+m+1}k!(k+m)!}\left(\frac1k+\frac1{k+m}\right),\qquad k\in\mathbb N^+. \tag{77}
$$

由此可以看出，$a_{2m}$ 也是任意常数。反复利用这一递推关系可得

$$
a_{2k+2m}=\frac{(-)^km!a_{2m}}{2^{2k}k!(k+m)!}
-\frac{(-)^k\beta}{2^{2k+m+1}k!(k+m)!}[\psi(k+m+1)+\psi(k+1)-\psi(m+1)-\psi(1)],
\qquad k\in\mathbb N. \tag{78}
$$

注意最后的表达式对于 $k=0$ 也成立。由于式 (76)，任意常数 $a_0$ 也可以用 $\beta$ 表出，故可将式 (75) 改写为

$$
a_{2k}=-\frac{(m-k-1)!}{2^{2k-m+1}k!}\beta,\qquad k=0,1,\cdots,m-1. \tag{79}
$$

于是得到级数解

$$
\begin{aligned}
y(x)&=\beta J_m(x)\ln x-\frac\beta2\sum_{k=0}^{m-1}\frac{(m-k-1)!}{k!}\left(\frac x2\right)^{2k-m}\\
&\quad-\frac\beta2\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}[\psi(k+m+1)+\psi(k+1)-\psi(m+1)-\psi(1)]\left(\frac x2\right)^{2k+m}\\
&\quad+2^mm!a_{2m}\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}\left(\frac x2\right)^{2k+m}.
\end{aligned} \tag{80}
$$

容易看出，最后一项求和正是第一解 $J_m(x)$，因此上式已经是通解；同时我们看到，如果取 $\beta=0$，则只能得到第一解，所以，$m\in\mathbb N^+$ 时的第二解也一定包含对数函数。我们当然可以取 $a_{2m}=0$ 和 $\beta=1$ 而得到第二解，但习惯上是取

$$
\beta=\frac2\pi,\qquad a_{2m}=-\frac{1}{2^mm!\pi}[\psi(m+1)+\psi(1)+2\ln2], \tag{81}
$$

这样得到的第二解就是 $y_2(x)=N_m(x)$，其形式已经在前面给出。这里系数的取法也是为了使得结果与式 (63) 的极限一致。

### §4.6 \*Neumann 函数的 Frobenius 级数解法

本小节介绍另一种求第二解的方法，称为 Frobenius 方法，这种方法也可以用于求解其它方程在正则奇点处的第二解，尤其是第二解包含对数函数的情况。对于第二解不包含对数函数的情况，常规的级数解法是方便的。所以我们只考虑 $\nu=m\in\mathbb N$ 的情况。

首先考虑 $m=0$ 的情况，将 $\nu=0$ 代入式 (50) 可得

$$
Ly=s^2a_0x^s+(s+1)^2a_1x^{s+1}+\sum_{k=0}^{\infty}\big[(k+s+2)^2a_{k+2}+a_k\big]x^{k+s+2}. \tag{82}
$$

暂时不确定 $s$，取 $a_0$ 为与 $s$ 无关的任意常数，$a_1=0$，而其它系数由下列递推关系确定

$$
(k+s+2)^2a_{k+2}=-a_k, \tag{83}
$$

这样得到的级数

$$
y(x,s)=\sum_{k=0}^{\infty}a_{2k}(s)x^{2k+s}
=a_0\sum_{k=0}^{\infty}\frac{(-)^k\Gamma^2(s/2+1)}{2^{2k}\Gamma^2(s/2+k+1)}x^{2k+s} \tag{84}
$$

满足

$$
Ly(x,s)=s^2a_0x^s. \tag{85}
$$

显然 $Ly(x,0)=0$，即 $y(x,0)$ 是解，容易看出 $y(x,0)=a_0J_0(x)$，取 $a_0=1$ 即得 $y_1(x)=J_0(x)$，这是熟知的第一解。将上式两边对 $s$ 求导，可得

$$
L\frac{\partial y(x,s)}{\partial s}=2sa_0x^s+s^2a_0x^s\ln x. \tag{86}
$$

显然，

$$
L\left.\frac{\partial y(x,s)}{\partial s}\right|_{s=0}=0,
$$

所以 $[\partial y(x,s)/\partial s]_{s=0}$ 也是方程的解，此即第二解。易得

$$
\frac{\partial y(x,s)}{\partial s}=\sum_{k=0}^{\infty}a_{2k}(s)x^{2k+s}\ln x
+\sum_{k=0}^{\infty}a_{2k}(s)\frac{\partial\ln a_{2k}(s)}{\partial s}x^{2k+s},
$$

算出第二项中的导数，然后在各项中代入 $s=0$，可得

$$
\left.\frac{\partial y(x,s)}{\partial s}\right|_{s=0}
=a_0J_0(x)\ln x-a_0\sum_{k=0}^{\infty}\frac{(-)^k}{(k!)^2}[\psi(k+1)-\psi(1)]\left(\frac x2\right)^{2k}.
$$

取 $a_0=2/\pi$，得第二解

$$
y_2(x)=\frac2\pi J_0(x)\ln x-\frac2\pi\sum_{k=0}^{\infty}\frac{(-)^k}{(k!)^2}\psi(k+1)\left(\frac x2\right)^{2k}+\frac2\pi\psi(1)J_0(x). \tag{87}
$$

对上式加上 $-(2/\pi)[\psi(1)+\ln2]J_0(x)$ 即得 $N_0(x)$。

其次考虑 $m\in\mathbb N^+$ 的情况，将 $\nu=m$ 代入式 (50) 可得

$$
Ly=(s^2-m^2)a_0x^s+[(s+1)^2-m^2]a_1x^{s+1}+\sum_{k=0}^{\infty}\left\{[(k+s+2)^2-m^2]a_{k+2}+a_k\right\}x^{k+s+2}. \tag{88}
$$

暂时不确定 $s$，取 $a_0$ 为与 $s$ 无关的任意常数，$a_1=0$，而其它系数由下列递推关系确定

$$
[(k+s+2)^2-m^2]a_{k+2}=-a_k, \tag{89}
$$

这样得到的级数

$$
y(x,s)=\sum_{k=0}^{\infty}a_{2k}(s)x^{2k+s}
=a_0\sum_{k=0}^{\infty}\frac{(-)^k\Gamma((s+m)/2+1)\Gamma((s-m)/2+1)}{2^{2k}\Gamma((s+m)/2+k+1)\Gamma((s-m)/2+k+1)}x^{2k+s} \tag{90}
$$

满足

$$
Ly(x,s)=(s-m)(s+m)a_0x^s. \tag{91}
$$

显然 $Ly(x,\pm m)=0$，即 $y(x,m)$ 和 $y(x,-m)$ 都是解。容易看出 $y(x,m)=2^mm!a_0J_m(x)$，取 $a_0=1/2^mm!$ 即得 $y_1(x)=J_m(x)$，这是熟知的第一解。然而，当我们将 $s=-m$ 代入式 (90) 计算 $y(x,-m)$ 的具体形式时，就会遇到困难。当 $k\le m-1$ 时，$a_{2k}(-m)$ 的分子和分母都是无穷大，如果取 $s\to-m$ 的极限，还可以得到有限的结果，但当 $k\ge m$ 时，$a_{2k}(-m)$ 的分母为有限，而分子是无穷大。

解决上述困难的方法是取一个无穷小的 $a_0$，比如

$$
a_0=(s+m)c_0, \tag{92}
$$

其中 $c_0$ 是与 $s$ 无关的常数。这样可以避免 $k\ge m$ 时 $a_{2k}(-m)$ 成为无穷大，但 $k\le m-1$ 时，$a_{2k}(-m)$ 就必然成为零，所以求和实际上是从 $k=m$ 项开始，而得到的解 $y(x,-m)$ 与 $y_1(x)=J_m(x)$ 成正比（见下），所以并不是线性独立的第二解。这与我们已经知道的结果 (64) 基本上是一回事。

按式 (92) 的取法，$a_0$ 已经成为 $s$ 的函数，于是式 (91) 变成

$$
Ly(x,s)=(s-m)(s+m)^2c_0x^s. \tag{93}
$$

将上式两边对 $s$ 求导，可得

$$
L\frac{\partial y(x,s)}{\partial s}=(s+m)^2c_0x^s+2(s-m)(s+m)c_0x^s+(s-m)(s+m)^2c_0x^s\ln x. \tag{94}
$$

显然，

$$
L\left.\frac{\partial y(x,s)}{\partial s}\right|_{s=-m}=0,
$$

所以 $[\partial y(x,s)/\partial s]_{s=-m}$ 也是方程的解，此即第二解。如前

$$
\frac{\partial y(x,s)}{\partial s}=\sum_{k=0}^{\infty}a_{2k}(s)x^{2k+s}\ln x
+\sum_{k=0}^{\infty}a_{2k}(s)\frac{\partial\ln a_{2k}(s)}{\partial s}x^{2k+s}. \tag{95}
$$

按式 (90) 和 (92)，

$$
a_{2k}(s)=c_0\frac{(-)^k(s+m)\Gamma((s+m)/2+1)\Gamma((s-m)/2+1)}{2^{2k}\Gamma((s+m)/2+k+1)\Gamma((s-m)/2+k+1)}, \tag{96}
$$

而

$$
\frac{\partial\ln a_{2k}(s)}{\partial s}=\frac1{s+m}+\frac12\psi\left(\frac{s+m}{2}+1\right)+\frac12\psi\left(\frac{s-m}{2}+1\right)
-\frac12\psi\left(\frac{s+m}{2}+k+1\right)-\frac12\psi\left(\frac{s-m}{2}+k+1\right). \tag{97}
$$

现在需要计算 $s\to-m$ 时式 (95) 中的各项系数。注意到宗量为零或负整数时，$\Gamma$ 函数和 $\psi$ 函数都有奇性，所以在计算时需要引用下面的公式：

$$
\Gamma(1-x)=\frac{\pi}{\Gamma(x)\sin\pi x},\qquad \psi(1-x)=\psi(x)+\pi\cot\pi x, \tag{98}
$$

其中后者来源于前者，前者证明从略。计算可得当 $k\le m-1$ 时，

$$
a_{2k}(s)\to\frac{(s+m)\Gamma(m-k)}{2^{2k}k!\Gamma(m)}c_0\to0, \tag{99a}
$$

$$
\frac{\partial\ln a_{2k}(s)}{\partial s}\to\frac1{s+m}+\frac12[\psi(1)+\psi(m)-\psi(k+1)-\psi(m-k)]\to\infty, \tag{99b}
$$

$$
a_{2k}(s)\frac{\partial\ln a_{2k}(s)}{\partial s}\to\frac{\Gamma(m-k)}{2^{2k}k!\Gamma(m)}c_0, \tag{99c}
$$

而当 $k\ge m$ 时，

$$
a_{2k}(s)\to\frac{(-)^{k-m+1}}{2^{2k-1}k!(k-m)!\Gamma(m)}c_0, \tag{100a}
$$

$$
\frac{\partial\ln a_{2k}(s)}{\partial s}\to\frac12[\psi(1)+\psi(m)-\psi(k+1)-\psi(k-m+1)]. \tag{100b}
$$

于是得到

$$
\begin{aligned}
\left.\frac{\partial y(x,s)}{\partial s}\right|_{s=-m}
&=-\frac{2c_0}{2^m\Gamma(m)}\sum_{k=m}^{\infty}\frac{(-)^{k-m}}{k!(k-m)!}\left(\frac x2\right)^{2k-m}\ln x\\
&\quad+\frac{c_0}{2^m\Gamma(m)}\sum_{k=0}^{m-1}\frac{(m-k-1)!}{k!}\left(\frac x2\right)^{2k-m}\\
&\quad+\frac{c_0}{2^m\Gamma(m)}\sum_{k=m}^{\infty}\frac{(-)^{k-m}}{k!(k-m)!}[\psi(k+1)+\psi(k-m+1)-\psi(1)-\psi(m)]\left(\frac x2\right)^{2k-m}.
\end{aligned} \tag{101}
$$

其中的第一项如果去掉 $\ln x$ 因子就是 $y(x,-m)$，容易看出其中的求和正是 $J_m(x)$，这就是前面提到的结论。取 $c_0=-2^m\Gamma(m)/\pi$，并对第一项和第三项作求和指标变换 $k=k'+m$，再将 $k'$ 重新写成 $k$，即得

$$
\begin{aligned}
y_2(x)&=\frac2\pi J_m(x)\ln x-\frac1\pi\sum_{k=0}^{m-1}\frac{(m-k-1)!}{k!}\left(\frac x2\right)^{2k-m}\\
&\quad-\frac1\pi\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}[\psi(k+m+1)+\psi(k+1)-\psi(1)-\psi(m)]\left(\frac x2\right)^{2k+m}.
\end{aligned} \tag{102}
$$

对上式加上 $-(1/\pi)[\psi(1)+\psi(m)+2\ln2]J_m(x)$ 即得 $N_m(x)$。

由上述求解过程可以看出，用 Frobenius 级数解法求 $N_0(x)$ 比较简单，但求 $N_m(x)$ 则并不比常规解法来得方便。

### §4.7 \*一般 Neumann 函数的整数阶极限

本小节补充证明式 (65)。

由式 (63) 和 L'Hospital 法则可得

$$
\lim_{\nu\to m}N_\nu(x)=\frac1\pi\lim_{\nu\to m}\left[\frac{\partial J_\nu(x)}{\partial\nu}-(-)^m\frac{\partial J_{-\nu}(x)}{\partial\nu}\right]. \tag{103}
$$

由式 (56) 易得

$$
\lim_{\nu\to m}\frac{\partial J_\nu(x)}{\partial\nu}
=J_m(x)\ln\frac x2-\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}\psi(k+m+1)\left(\frac x2\right)^{2k+m}, \tag{104}
$$

而

$$
\frac{\partial J_{-\nu}(x)}{\partial\nu}
=-\sum_{k=0}^{\infty}\frac{(-)^k}{k!\Gamma(k-\nu+1)}\left(\frac x2\right)^{2k-\nu}\ln\frac x2
+\sum_{k=0}^{\infty}\frac{(-)^k}{k!\Gamma(k-\nu+1)}\psi(k-\nu+1)\left(\frac x2\right)^{2k-\nu}.
$$

对上式取 $\nu\to m$ 的极限时，需要注意 $\Gamma$ 函数和 $\psi$ 函数的奇性。首先看第一项，当 $k\le m-1$ 时，$\Gamma(k-\nu+1)\to\infty$，故求和实际从 $k=m$ 项开始，不难求得该项的极限为 $-(-)^mJ_m(x)\ln(x/2)$。再看第二项，当 $k\le m-1$ 时，$\Gamma(k-\nu+1)\to\infty$，同时 $\psi(k-\nu+1)\to\infty$，利用式 (98) 可以求出有限的结果，不过需要注意，当 $m=0$ 时这一部分是不存在的；当 $k\ge m$ 时，各因子均无奇性。结果为

$$
\begin{aligned}
\lim_{\nu\to m}\frac{\partial J_{-\nu}(x)}{\partial\nu}
&=-(-)^mJ_m(x)\ln\frac x2+(-)^m\sum_{k=0}^{m-1}\frac{(m-k-1)!}{k!}\left(\frac x2\right)^{2k-m}\\
&\quad+(-)^m\sum_{k=0}^{\infty}\frac{(-)^k}{k!(k+m)!}\psi(k+1)\left(\frac x2\right)^{2k+m}.
\end{aligned} \tag{105}
$$

将以上结果代入式 (103)，结果即为式 (60) 给出的 $N_m(x)$，当 $m=0$ 时，需去掉其中第二项有限和。

## §5 Sturm–Liouville 本征值问题

### §5.1 Sturm–Liouville 本征值问题的一般提法

二阶线性齐次常微分方程的一般形式由式 (1) 给出。为了下面符号上的方便，这里将它改写为

$$
y''(x)+P(x)y'(x)+Q(x)y(x)=0. \tag{106}
$$

从数理方程分离变量得到的常微分方程一般含有待定常数，记作 $\lambda$，如果将含有 $\lambda$ 的项单独写出来，方程的形式通常是

$$
y''(x)+P(x)y'(x)-\widetilde Q(x)y(x)+\lambda\widetilde\rho(x)y(x)=0. \tag{107}
$$

两边同乘以

$$
k(x)\equiv\exp\left(\int_{x_0}^{x}P(\xi)d\xi\right), \tag{108}
$$

则得

$$
k(x)y''(x)+k(x)P(x)y'(x)-k(x)\widetilde Q(x)y(x)+\lambda k(x)\widetilde\rho(x)y(x)=0.
$$

显然，$k(x)P(x)=k'(x)$，故上式可以写成

$$
\frac d{dx}\left[k(x)\frac{dy(x)}{dx}\right]-q(x)y(x)+\lambda\rho(x)y(x)=0,\qquad (a<x<b) \tag{109}
$$

其中 $q(x)=k(x)\widetilde Q(x)$，$\rho(x)=k(x)\widetilde\rho(x)$，$(a,b)$ 为求解区间。上式称为 Sturm–Liouville 型方程，其中 $\rho(x)$ 称为权函数。以上推导说明一般形式 (106) 与 Sturm–Liouville 形式是等价的。后者对于本节的讨论是方便的。

由于方程 (109) 是由数理方程分离变量得到的，所以在区间端点 $a$、$b$ 通常附有边界条件，满足边界条件的解并不一定存在，除非 $\lambda$ 取某些特定值，这样的 $\lambda$ 值称为本征值，相应的解称为本征函数。Sturm–Liouville 方程与边界条件一起构成的问题称为 Sturm–Liouville 本征值问题。

本征值问题的类型决定于边界条件的类型，主要有以下几种。

1. 第一、二、三类边界条件。比如以下本征值问题：

$$
y''(x)+\lambda y(x)=0,\qquad (a<x<b) \tag{110a}
$$

$$
(\alpha y'-\beta y)|_{x=a}=0,\qquad (\gamma y'+\delta y)|_{x=b}=0, \tag{110b}
$$

其中 $\alpha$、$\beta$、$\gamma$、$\delta\ge0$，但 $\alpha^2+\beta^2\ne0$，$\gamma^2+\delta^2\ne0$。与 Sturm–Liouville 方程比较可知 $k(x)=1$，$q(x)=0$，$\rho(x)=1$。

2. 自然边界条件。比如 Legendre 方程的本征值问题：

$$
\frac d{dx}\left[(1-x^2)\frac{dy}{dx}\right]+\lambda y=0,\qquad (-1<x<1) \tag{111a}
$$

$$
|y(\pm1)|<\infty. \tag{111b}
$$

其中的边界条件即是自然边界条件。与 Sturm–Liouville 方程比较可知 $k(x)=1-x^2$，$q(x)=0$，$\rho(x)=1$。注意 $k(\pm1)=0$，而 $x=\pm1$ 处均有自然边界条件。

又比如 Bessel 方程的本征值问题：

$$
\frac d{dx}\left(x\frac{dy}{dx}\right)-\frac{m^2}{x}y+\lambda xy=0,\qquad (0<x<a) \tag{112a}
$$

$$
y(0)=0\ \text{或}\ |y(0)|<\infty,\qquad \alpha y'(a)+\beta y(a)=0, \tag{112b}
$$

其中 $\alpha$、$\beta\ge0$，但 $\alpha^2+\beta^2\ne0$；$x=0$ 处的边界条件即是自然边界条件。与 Sturm–Liouville 方程比较可知 $k(x)=x$，$q(x)=m^2/x$，$\rho(x)=x$。注意 $k(0)=0$，而 $x=0$ 处有自然边界条件。

另一个例子是球 Bessel 方程的本征值问题：

$$
\frac d{dx}\left(x^2\frac{dy}{dx}\right)-\lambda_l y+\lambda x^2y=0,\qquad (0<x<a) \tag{113a}
$$

$$
y(0)=0\ \text{或}\ |y(0)|<\infty,\qquad \alpha y'(a)+\beta y(a)=0, \tag{113b}
$$

其中 $\alpha$、$\beta\ge0$，但 $\alpha^2+\beta^2\ne0$；$x=0$ 处的边界条件即是自然边界条件。与 Sturm–Liouville 方程比较可知 $k(x)=x^2$，$q(x)=\lambda_l$，$\rho(x)=x^2$。注意 $k(0)=0$，而 $x=0$ 处有自然边界条件。又注意与上章的记号比较，这里 $x$ 代表径向球坐标 $r$，$y(x)$ 代表径向函数 $R(r)$，$\lambda_l$ 相当于上章的 $\lambda$，它由角向方程，即球函数方程的本征值问题决定，而 $\lambda$ 相当于上章的 $k^2$，它才是径向方程的本征值。

一般来说，端点 $a$（或 $b$）处出现自然边界条件的充要条件是 $k(a)=0$（或 $k(b)=0$）。

Sturm–Liouville 方程 (109) 可以改写为

$$
y''(x)+\frac{k'(x)}{k(x)}y'(x)+\frac{\lambda\rho(x)-q(x)}{k(x)}y(x)=0,
$$

与一般形式 (106) 比较可知 $P(x)=k'(x)/k(x)$，$Q(x)=[\lambda\rho(x)-q(x)]/k(x)$。从以上各例看到，在求解区间上，$\rho(x)\ge0$ 且性质良好；$q(x)\ge0$，只在端点可能有奇性，且最多为一阶极点；$k(x)\ge0$，没有奇性，只在端点可能为零，且最多为二阶零点。而且，当端点为 $k(x)$ 的一阶零点时，$q(x)$ 最多以其为一阶极点；当端点为 $k(x)$ 的二阶零点时，$q(x)$ 在该处性质良好。这样若端点为方程的奇点，则必为正则奇点。物理上遇到的情况大致如此。

不失一般性，以左端点 $a$ 为例，现在基本上可以断定，如果 $k(a)=0$，则 $a$ 处必有自然边界条件。反之，如果 $k(a)\ne0$，则 $a$ 处不会有自然边界条件。由于 $a$ 是方程的常点或正则奇点，我们可以按式 (30) 的形式展开 $P(x)=\sum_{k=0}^{\infty}p_k(x-a)^{k-1}$，$Q(x)=\sum_{k=0}^{\infty}q_k(x-a)^{k-2}$，然后按具体情况讨论如下。

首先，若 $a$ 是 $k(x)$ 的一阶零点，则它最多是 $q(x)$ 的一阶极点，于是有 $k(x)=(x-a)\varphi(x)$，$q(x)=\psi(x)/(x-a)$，其中 $\varphi(a)>0$ 而 $\psi(a)\ge0$。易得 $p_0=[(x-a)P(x)]|_{x=a}=[1+(x-a)\varphi'(x)/\varphi(x)]|_{x=a}=1$，$q_0=[(x-a)^2Q(x)]|_{x=a}=[\lambda(x-a)\rho(x)/\varphi(x)-\psi(x)/\varphi(x)]|_{x=a}=-\psi(a)/\varphi(a)\le0$，于是指标方程 (34) 变成 $s^2+q_0=0$，其中 $q_0$ 为实数且 $q_0\le0$，所以两根一正一负或均为零。若两根一正一负，则对应于 $s_2$ 的解含有 $(x-a)^{s_2}$ 项而在 $x=a$ 处发散；若两根均为零，则对应于 $s_2$ 的解含有 $\ln(x-a)$ 项而在 $x=a$ 处发散。对于物理问题，应该排除在 $x=a$ 处有奇性的解，因而 $x=a$ 处有自然边界条件。

其次，若 $a$ 是 $k(x)$ 的二阶零点，则 $k(x)=(x-a)^2\varphi(x)$，其中 $\varphi(a)>0$，而 $q(x)$ 在 $a$ 处性质良好。易得 $p_0=[2+(x-a)\varphi'(x)/\varphi(x)]|_{x=a}=2$。对于 $q_0$，难以得到一般结论，但物理上遇到的通常是 Helmholtz 方程或中心力场中的定态 Schrödinger 方程在球坐标系中分离变量后得到的径向方程（前者对应的径向方程即上面的球 Bessel 方程），这时 $\rho(a)=0$。对于这种满足 $\rho(a)=0$ 的情况，易得 $q_0=[\lambda\rho(a)-q(a)]/\varphi(a)=-q(a)/\varphi(a)\le0$，于是指标方程变成 $s^2+s+q_0=0$，其中 $q_0$ 为实数且 $q_0\le0$，由此易知其两根满足 $s_1\ge0$，$s_2<0$，于是对应于 $s_2$ 的解含有 $(x-a)^{s_2}$ 项而在 $x=a$ 处发散，因而 $x=a$ 处有自然边界条件。

最后，如果 $k(a)\ne0$，则容易推出 $p_0=0$，$q_0=0$，于是指标方程变成 $s(s-1)=0$，其两根为 $s_1=1$，$s_2=0$，所以两个解在 $a$ 处性质良好，不需要排除有奇性的解，因而在 $a$ 处也就不会有自然边界条件。

3. 周期性边界条件。如果 $k(a)=k(b)$，$q(a)=q(b)$，$\rho(a)=\rho(b)$，则可以对 Sturm–Liouville 方程附加周期性边界条件 $y(a)=y(b)$，$y'(a)=y'(b)$。比如下列本征值问题：

$$
y''(x)+\lambda y(x)=0,\qquad (0<x<2\pi) \tag{114a}
$$

$$
y(0)=y(2\pi),\qquad y'(0)=y'(2\pi). \tag{114b}
$$

其中方程与式 (110a) 相同，故 $k(x)=1$，$q(x)=0$，$\rho(x)=1$。

我们以前用的周期性边界条件是 $y(x+2\pi)=y(x)$，由此容易推出式 (114b)。反过来，由后者虽然不能推出前者，但结合方程 (114a)，易得 $y^{(n)}(0)=y^{(n)}(2\pi)$（$n=0,1,2,\cdots$），于是可得

$$
y(x+2\pi)=\sum_{k=0}^{\infty}\frac{y^{(n)}(2\pi)}{n!}x^n
=\sum_{k=0}^{\infty}\frac{y^{(n)}(0)}{n!}x^n=y(x).
$$

当然，这里假定 $y(x)$ 具有良好的解析性质，从而可以展开为 Taylor 级数。

### §5.2 Sturm–Liouville 本征值问题的一般结论

对于物理问题，Sturm–Liouville 方程 (109) 中的系数满足 $k(x)\ge0$，$q(x)\ge0$，$\rho(x)\ge0$（上面所举的例子均满足这些条件）。在这样的前提下，Sturm–Liouville 本征值问题有以下一般结论。

1. 所有本征值都是非负的，即 $\lambda\ge0$。

**注**　有了这一结论，以后求解本征值问题时，只要方程的系数满足上述条件，就可以立即排除 $\lambda<0$ 的可能性。对于象式 (110) 这样的本征值问题，这并不一定能减少多大的工作量，但是对于方程比较复杂的情况，比如 Bessel 方程的本征值问题，由此带来的方便是非常明显的。

**证明**　将式 (109) 两边同乘以 $y^*(x)$，得

$$
y^*(x)\frac d{dx}\left[k(x)\frac{dy(x)}{dx}\right]-q(x)y^*(x)y(x)+\lambda\rho(x)y^*(x)y(x)=0,
$$

对 $x$ 由 $a$ 到 $b$ 积分，可得

$$
\begin{aligned}
\lambda\int_a^b\rho(x)|y(x)|^2dx
&=\int_a^b q(x)|y(x)|^2dx-\int_a^b y^*(x)\frac d{dx}\left[k(x)\frac{dy(x)}{dx}\right]dx\\
&=\int_a^b q(x)|y(x)|^2dx+\int_a^b k(x)|y'(x)|^2dx-k(x)y^*(x)y'(x)\Big|_a^b\\
&\ge k(x)y^*(x)y'(x)\Big|_b^a
=k(a)y^*(a)y'(a)-k(b)y^*(b)y'(b).
\end{aligned}
$$

其中第二步作了分部积分，第三步是因为第二步中的两项积分均非负。因为 $y(x)$ 是本征函数，所以除了可能的若干零点外应不为零，而 $\rho(x)\ge0$ 而且一般只在端点处才可能为零，所以上式左边的积分是正数。现在我们只需证明右边非负就行了。先看右边第一项，如果 $x=a$ 处是第一类边界条件，则 $y(a)=0$，若是第二类边界条件，则 $y'(a)=0$，若是第三类边界条件 $\alpha y'(a)-\beta y(a)=0$，其中 $\alpha$、$\beta>0$，则 $k(a)y^*(a)y'(a)=(\beta/\alpha)k(a)|y(a)|^2>0$，若是自然边界条件，则 $k(a)=0$，因此在以上各种边界条件下，总有 $k(a)y^*(a)y'(a)\ge0$，同理也有 $-k(b)y^*(b)y'(b)\ge0$。若是周期性边界条件，则由于 $y(a)=y(b)$，$y'(a)=y'(b)$，且 $k(a)=k(b)$，故 $k(a)y^*(a)y'(a)-k(b)y^*(b)y'(b)=0$。因此，不论何种边界条件，上式右边总是非负的。证毕。

2. 存在无穷多分立的本征值：$\lambda_1\le\lambda_2\le\cdots\le\lambda_n\le\lambda_{n+1}\le\cdots$，且 $\lim_{n\to\infty}\lambda_n=+\infty$。除了周期性边界条件的情况，本征值都是非简并的，且 $y_{n+1}(x)$ 比 $y_n(x)$ 多一个零点。

这一结论的证明很困难，读者只要直接承认就可以了。应该指出，由于我们考虑的是二阶常微分方程，所以如果本征值有简并，其简并度只能是 2。

下面说明一下为什么除了周期性边界条件的情况，本征值都是非简并的。

假设对于某个本征值 $\lambda$ 存在两个线性独立的本征函数 $y_1(x)$ 和 $y_2(x)$，则两者均满足方程 (109)，且其中的 $\lambda$ 是相同的。对 $y_1(x)$ 的方程乘以 $y_2(x)$，对 $y_2(x)$ 的方程乘以 $y_1(x)$，所得两式相减，可得

$$
y_1(x)\frac d{dx}\left[k(x)\frac{dy_2(x)}{dx}\right]-y_2(x)\frac d{dx}\left[k(x)\frac{dy_1(x)}{dx}\right]=0,
$$

这可以化为

$$
\frac d{dx}\big[k(x)y_1(x)y_2'(x)-k(x)y_2(x)y_1'(x)\big]=0,
$$

于是

$$
k(x)[y_1(x)y_2'(x)-y_2(x)y_1'(x)]=C,
$$

其中 $C$ 是与 $x$ 无关的常数。现在根据一个端点，比如左端点 $a$ 的边界条件来确定常数 $C$。如果 $a$ 处有自然边界条件，则 $k(a)=0$，于是得 $C=0$。如果 $a$ 处是第一、二、三类边界条件，则 $\alpha y_1'(a)-\beta y_1(a)=0$，$\alpha y_2'(a)-\beta y_2(a)=0$，由于 $\alpha$、$\beta$ 不全为零，则 $y_1(a)y_2'(a)-y_2(a)y_1'(a)=0$，于是得 $C=0$。所以，除了周期性边界条件，总有

$$
k(x)[y_1(x)y_2'(x)-y_2(x)y_1'(x)]=0,
$$

但 $k(x)$ 只在端点才可能为零，于是得到

$$
y_1(x)y_2'(x)-y_2(x)y_1'(x)=0,
$$

由此可以证明 $y_1(x)\propto y_2(x)$，即 $y_1(x)$ 与 $y_2(x)$ 线性相关，这与假设矛盾。如果是周期性边界条件，则不能推出以上结果，因此可能存在简并。实际上，我们已经知道周期性边界条件下是有简并的。

3. 对应于不同本征值的本征函数在区间 $[a,b]$ 上带权正交：

$$
\int_a^b y_m^*(x)y_n(x)\rho(x)\,dx=0,\qquad (\lambda_m\ne\lambda_n). \tag{115}
$$

**注**　① 本征函数族的正交性对于后面计算广义 Fourier 级数的系数非常重要。有时候，通过直接计算来验证本征函数族的正交性有一定的困难，所以上述结论给我们带来了很大的方便。② 与三角函数族的正交性相比，这里有两点推广：一是多了权函数 $\rho(x)$，如果 $\rho(x)=1$，就是普通正交；二是考虑了本征函数是复值函数的情况（注意自变量仍是实数，只是函数值可取复数，故并非复变函数），比如本征值问题 (114) 的本征函数族 $\{e^{imx}\}$，其中 $0\le x\le2\pi$，而 $m=0,\pm1,\pm2,\cdots$。③ 对应于同一本征值的两个本征函数（如果有简并）不一定相互正交，但我们总可以取其适当的线性组合（线性组合后的函数仍是对应于同一本征值的本征函数）使得它们相互正交，通过这样的做法（称为 Schmidt 正交化），可以使所有的本征函数相互正交。

**证明**　将 $y_n(x)$ 的方程，即式 (109) 两边同乘以 $y_m^*(x)$，得

$$
y_m^*(x)\frac d{dx}\left[k(x)\frac{dy_n(x)}{dx}\right]-q(x)y_m^*(x)y_n(x)+\lambda_n\rho(x)y_m^*(x)y_n(x)=0,
$$

对上式交换 $m$ 和 $n$，并取复共轭，得到

$$
y_n(x)\frac d{dx}\left[k(x)\frac{dy_m^*(x)}{dx}\right]-q(x)y_m^*(x)y_n(x)+\lambda_m\rho(x)y_m^*(x)y_n(x)=0.
$$

两式相减并从 $a$ 到 $b$ 积分得

$$
\begin{aligned}
(\lambda_m-\lambda_n)\int_a^b\rho(x)y_m^*(x)y_n(x)\,dx
&=\int_a^b\left\{y_m^*(x)[k(x)y_n'(x)]'-y_n(x)[k(x)y_m^*(x)]'\right\}dx\\
&=\int_a^b\left[k(x)y_m^*(x)y_n'(x)-k(x)y_m^{*\prime}(x)y_n(x)\right]'dx\\
&=k(x)[y_m^*(x)y_n'(x)-y_m^{*\prime}(x)y_n(x)]\Big|_a^b.
\end{aligned}
$$

对周期性边界条件，上式右边显然为零。对其它边界条件，以 $a$ 或 $b$ 代入分别为零。以 $x=b$ 为例，若为自然边界条件，则 $k(b)=0$；若为第一、二、三类边界条件，则 $\gamma y_n'(b)+\delta y_n(b)=0$，$\gamma y_m^{*\prime}(b)+\delta y_m^*(b)=0$，其中第二式取了复共轭并利用了 $\gamma$ 和 $\delta$ 都是实数的事实。由于 $\gamma$ 和 $\delta$ 不全为零，故系数行列式为零，即 $y_m^*(b)y_n'(b)-y_m^{*\prime}(b)y_n(b)=0$。因此，总有 $k(b)[y_m^*(b)y_n'(b)-y_m^{*\prime}(b)y_n(b)]=0$。于是上式右边为零，考虑到 $\lambda_m\ne\lambda_n$，即得式 (115)。

4. 本征函数族 $\{y_n(x)\}_{n=1}^{\infty}$ 在区间 $[a,b]$ 上是完备的，这就是说，区间 $[a,b]$ 上的任意函数 $f(x)$，只要解析性质良好且与本征函数族满足相同的边界条件，就一定可以用 $\{y_n(x)\}_{n=1}^{\infty}$ 展开为广义 Fourier 级数：

$$
f(x)=\sum_{n=1}^{\infty}f_ny_n(x), \tag{116}
$$

其中展开系数为

$$
f_n=\frac{\displaystyle\int_a^b y_n^*(x)f(x)\rho(x)\,dx}{\displaystyle\int_a^b y_n^*(x)y_n(x)\rho(x)\,dx},\qquad n\in\mathbb N^+. \tag{117}
$$

**注**　我们已经知道，本征函数族的完备性是分离变量法的理论基础。所以，这一结论的重要性是显而易见的。

这一性质的证明也比较困难，读者只要掌握结论就可以了。不过，只要承认式 (116)，即承认本征函数族的完备性，就很容易推出展开系数。事实上，将式 (116) 中的求和指标 $n$ 换为 $k$，然后两边同乘以 $y_n^*(x)\rho(x)$ 并积分，由于正交性（假设简并的本征函数也已经正交化），右边的积分对求和有贡献的只有 $k=n$ 一项，由此立得式 (117)。这与我们以前的做法是一致的，只是现在多了权函数 $\rho(x)$，并且出现了本征函数的复共轭。

如果 $y_n(x)$ 是本征函数，则 $u_n(x)=Cy_n(x)$（其中 $C\ne0$ 是常数）也是本征函数，仍然对应于同一个本征值。只要取

$$
C=\left[\int_a^b y_n^*(x)y_n(x)\rho(x)\,dx\right]^{-1/2},
$$

就可使

$$
\int_a^b u_n^*(x)u_n(x)\rho(x)\,dx=1. \tag{118}
$$

$u_n(x)$ 称为归一化的本征函数。如果在展开式 (116) 中采用归一化的本征函数 $u_n(x)$，则展开系数 (117) 的分母为 1（当然分子中的 $y_n(x)$ 应代以 $u_n(x)$）。在量子力学中经常需要计算涉及本征函数的积分，所以采用归一化的本征函数能带来很大的方便。在本课程里，由于相关的计算并不多，所以我们没有强调本征函数的归一化。

假设简并的本征函数也已经正交化，则正交性和归一化关系可以统一写成

$$
\int_a^b u_m^*(x)u_n(x)\rho(x)\,dx=\delta_{mn}. \tag{119}
$$

这称为正交归一关系，简称正交系 (orthonormal relation)。

在展开式 (116) 中采用归一化的本征函数 $u_n(x)$，并将展开系数的表达式代回展开式中，可以得到

$$
f(x)=\int_a^b\left[\sum_{n=1}^{\infty}u_n(x)u_n^*(x')\rho(x')\right]f(x')\,dx',
$$

由于 $f(x)$ 是任意函数，故有

$$
\sum_{n=1}^{\infty}u_n(x)u_n^*(x')\rho(x')=\delta(x-x'). \tag{120}
$$

这是完备性关系 (completeness relation) 的数学表示。

**习题**　长为 $l$ 的均匀细杆，侧面绝热，左端保持恒定温度零度，右端按 Newton 冷却定律与温度为零度的介质交换热量，即右端边界条件为 $(u+h\partial u/\partial x)|_{x=l}=0$（其中 $h>0$），已知 $u|_{t=0}=\varphi(x)$。求杆上温度的变化规律。

## 补充习题

1. 归纳出 $a_{l-2k}$ 的表达式 (24)，并用归纳法加以证明。

2. 从头推导 Bessel 方程的第二解 $J_{-\nu}(x)$，设 $\nu$ 不等于整数或半奇数。

3. 在 $x_0=0$ 的邻域上求解 Hermite 方程 $y''-2xy'+(\lambda-1)y=0$。当 $\lambda$ 取何值时可以使级数解退化为多项式？适当选取级数解中的任意常数，可以使多项式的最高次幂项具有形式 $(2x)^n$，这些多项式称为 Hermite 多项式，记作 $H_n(x)$。写出前几个 Hermite 多项式。

4. 在 $x_0=0$ 的邻域上求解 Laguerre 方程 $xy''+(1-x)y'+\lambda y=0$。当 $\lambda$ 取何值时可以使级数解退化为多项式？适当选取级数解中的任意常数，可以使多项式的最高次幂项具有形式 $(-x)^n$，这些多项式称为 Laguerre 多项式，记作 $L_n(x)$。写出前几个 Laguerre 多项式。

**注**　Hermite 多项式和 Laguerre 多项式在量子力学中有重要的应用。

5. 在 $x_0=0$ 的邻域上求解超几何方程 (Gauss 方程) $x(x-1)y''+[(1+\alpha+\beta)x-\gamma]y'+\alpha\beta y=0$。

6. 在 $x_0=0$ 的邻域上求解合流超几何方程 (Kummer 方程) $xy''+(\gamma-x)y'-\alpha y=0$。

**注**　超几何方程与合流超几何方程都是数学物理中的重要方程，其解超几何函数与合流超几何函数都是重要的特殊函数（也称为高级超越函数）。

7. 将上述各题中的微分方程化为 Sturm–Liouville 型。
