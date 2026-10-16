# 第十七章 分离变量法总结

到现在为止，我们已经处理了几种典型的偏微分方程定解问题，介绍了求解这些定解问题的一种有效方法——分离变量法。这种方法当然有一定的适用条件，例如，要求方程和定解条件都是线性的，因此解具有叠加性。在第十三章中，我们曾经结合具体的求解过程，分析了这种解法对于定解问题的要求。特别是，曾经指出（见 §13.1）这种方法是否能够普遍地应用于求解偏微分方程定解问题，在理论上取决于下列几个关键问题：

1. 本征值问题是否一定有解，或者说，在什么条件下，本征值问题一定有解；
2. 定解问题的解是否一定可以按照某一组本征函数展开，或者说，在什么条件下，本征函数组是完备的；
3. 本征函数是否一定具有正交性。

这一章就要从理论上回答这几个问题，从而为分离变量法奠定一个坚实的理论基础。

为此，我们从自伴算符的本征值问题入手，开始本章的讨论。由于我们所关心的这种特定背景，所以这里讨论的算符仅限于微分算符，并且主要是二阶常微分算符。而在讨论自伴算符的本征值问题之前，先要建立内积空间和函数空间的概念。作为必需的基础知识，读者应当已经掌握线性空间的基本概念与计算方法。

## *§17.1 内积空间

设在数域 $\mathbb K$ 上定义了 $n$ 维线性空间 $V$，它的元素（矢量）用 $x,y,\cdots$ 表示。作为（三维）矢量分析中矢量数量积概念的推广，先定义 $n$ 维矢量的内积。

**定义 17.1**　（定义在实数或复数域 $\mathbb K$ 上的）线性空间 $V$ 中矢量 $x$ 和 $y$ 的内积 $(x,y)$ 是一个 $V\times V\mapsto\mathbb K$ 的映射，$\forall x\in V,\forall y\in V$，都存在数域 $\mathbb K$ 内的唯一一个数与之对应，记为 $(x,y)\in\mathbb K$，满足：

1. $(x,y)$ 与 $(y,x)$ 互为复共轭，即 $(x,y)=(y,x)^*$；
2. $(\alpha x+\beta y,z)=\alpha^*(x,z)+\beta^*(y,z)$，其中 $\alpha\in\mathbb K,\beta\in\mathbb K$；
3. 对于任何 $x\in V$，$(x,x)\ge0$；当且仅当 $x=0\in V$ 时，$(x,x)=0\in\mathbb K$。

定义了内积的线性空间称为内积空间。具有内积的实线性空间称为欧氏（Euclid）空间（Euclidean space）；具有内积的复线性空间称为酉空间（unitary space）。

**例 17.1**　若

$$
x=
\begin{pmatrix}
x_1\\
x_2\\
\vdots\\
x_n
\end{pmatrix},
\qquad
y=
\begin{pmatrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{pmatrix}
$$

是实数域上的 $n$ 维列矢量，全体实数域上的 $n$ 维列矢量的集合在矢量加法以及矢量和实数的乘法下构成一个实数域上的线性空间。$P$ 为（给定的）$n\times n$ 实对角正定矩阵①，可定义矢量 $x$ 和 $y$ 的内积为

$$
(x,y)
=
(x_1,x_2,\cdots,x_n)
\begin{pmatrix}
P_{11}&0&\cdots&0\\
0&P_{22}&\cdots&0\\
\vdots&\vdots&&\vdots\\
0&0&\cdots&P_{nn}
\end{pmatrix}
\begin{pmatrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{pmatrix}.
$$

特别是，如果 $P$ 是单位矩阵，则矢量 $x$ 和 $y$ 的内积即为

$$
(x,y)
=
x_1y_1+x_2y_2+\cdots+x_ny_n
=
\sum_{i=1}^{n}x_i y_i.
$$

可见，内积是我们熟悉的矢量数量积的推广。

**例 17.2**　$0\le t\le1$ 区间内的实变量 $t$ 的所有复系数的多项式的集合，在多项式加法以及多项式和复数的乘法下构成一个复数域上的线性空间，记为 $\mathcal H$。$\rho(t)$ 已知，它是定义在 $0\le t\le1$ 区间内的非负函数，且不恒为 0，称为权函数②。设 $x(t)$ 和 $y(t)$ 是线性空间 $\mathcal H$ 中的两个元素（也可称为矢量，在这里是多项式），它们的内积可以定义为

$$
(x,y)
=
\int_0^1x^*(t)y(t)\rho(t)\,dt.
\tag{17.1}
$$

请读者证明对线性空间 $\mathcal H$ 内的任意两个元素，内积 (17.1) 都存在而且满足内积定义中的三条要求。

根据内积定义中的第 1 条要求，可以看出，不论是欧氏空间或酉空间，矢量和它自身的内积总是实数，这样第 3 条要求中的不等式才有意义。因此定义

$$
(x,x)^{1/2}
\equiv
\|x\|
\tag{17.2}
$$

为矢量 $x$ 的模（即矢量 $x$ 的“长度”）。而由内积定义中的第 1 和第 2 条要求，有

$$
(x,\alpha y)
=
\alpha(x,y).
\tag{17.3}
$$

因此

$$
\|\alpha x\|
=
(\alpha x,\alpha x)^{1/2}
=
[\alpha\alpha^*(x,x)]^{1/2}
=
|\alpha|\,\|x\|.
\tag{17.4}
$$

任何一个非零矢量除以它的模就成为“单位长度”的矢量，或称为归一化的矢量，

$$
\left(
\frac{x}{\|x\|},
\frac{x}{\|x\|}
\right)
=
1.
$$

在建立了内积定义后，就可以引入矢量正交的概念：当且仅当 $(x,y)=0$ 时，两矢量 $x,y$ 正交。显然，零矢量和任何矢量都正交。

任何一组线性无关③的非零矢量 $\{y_1,y_2,y_3,\cdots\}$，均可通过标准步骤使之两两正交。例如，作线性组合

$$
x_1=y_1,
\qquad
x_2=y_2+\alpha_{21}x_1,
\qquad
x_3=y_3+\alpha_{31}x_1+\alpha_{32}x_2,
\qquad
\cdots,
$$

令它们两两正交：

$$
(x_1,x_2)
=
(x_1,y_2+\alpha_{21}x_1)
=
(x_1,y_2)+\alpha_{21}(x_1,x_1)
=
0,
$$

$$
(x_1,x_3)
=
(x_1,y_3+\alpha_{31}x_1+\alpha_{32}x_2)
=
(x_1,y_3)+\alpha_{31}(x_1,x_1)
=
0,
$$

$$
(x_2,x_3)
=
(x_2,y_3+\alpha_{31}x_1+\alpha_{32}x_2)
=
(x_2,y_3)+\alpha_{32}(x_2,x_2)
=
0,
$$

$\vdots$

所以

$$
\alpha_{21}
=
-\frac{(x_1,y_2)}{(x_1,x_1)},
\qquad
\alpha_{31}
=
-\frac{(x_1,y_3)}{(x_1,x_1)},
\qquad
\alpha_{32}
=
-\frac{(x_2,y_3)}{(x_2,x_2)}.
$$

更普遍的结果是

$$
\alpha_{jk}
=
-\frac{(x_k,y_j)}{(x_k,x_k)}.
$$

这样的步骤称为 Schmidt 正交化。

**定义 17.2**　若对于所有正整数 $i$ 和 $j$，都有 $(x_i,x_j)=\delta_{ij}$，则称矢量组 $\{x_1,x_2,\cdots\}$ 是正交归一的。

设在 $n$ 维线性空间 $V$ 中有一组正交归一矢量

$$
\{x_1,x_2,\cdots,x_k\},
\qquad
k\le n.
$$

对于任意一个矢量 $x\in V$，可求出 $\alpha_i=(x_i,x)$。显然，应当有

$$
\begin{aligned}
\left\|
x-\sum_{i=1}^{k}\alpha_i x_i
\right\|^2
&\equiv
\left(
x-\sum_{i=1}^{k}\alpha_i x_i,
x-\sum_{i=1}^{k}\alpha_i x_i
\right)\\
&=
(x,x)
-
\sum_{i=1}^{k}\alpha_i^*\alpha_i
-
\sum_{i=1}^{k}\alpha_i\alpha_i^*
+
\sum_{i,j=1}^{k}\alpha_i^*\alpha_j\delta_{ij}\\
&=
(x,x)
-
\sum_{i=1}^{k}\alpha_i^*\alpha_i
\ge0.
\end{aligned}
$$

由此就得到一个重要的不等式——Bessel 不等式：

$$
(x,x)
\ge
\sum_{i=1}^{k}|(x_i,x)|^2,
\qquad
\text{即}
\qquad
\|x\|^2
\ge
\sum_{i=1}^{k}|(x_i,x)|^2.
\tag{17.5}
$$

由 Bessel 不等式可以推出 Schwarz 不等式：若 $x,y$ 是内积空间中的两个矢量，则

$$
|(x,y)|
\le
\|x\|\,\|y\|.
\tag{17.6}
$$

如果 $y$ 是零矢量，$y=0$，上式中的等号成立。如果 $y\ne0$，则可在 Bessel 不等式 (17.5) 中的 $k$ 取为 1，且 $x_1=y/\|y\|$，于是就有

$$
\left|
\left(
\frac{y}{\|y\|},x
\right)
\right|^2
\le
\|x\|^2,
$$

由此即可证得 (17.6) 式。

**定义 17.3**　在有限维线性空间中，如果一组正交归一的矢量（称为一个正交归一矢量集），并不包含在另一个更大的正交归一矢量集之中，则称该正交归一矢量集是完备的。

正交归一的矢量一定线性无关。反之，任何一组线性无关的矢量都可以正交归一化。所以 $n$ 维线性空间中的任何一组 $n$ 个正交归一矢量都可以构成此空间的基，称为正交归一基（线性代数中称为正交标准基）。

① 设 $M$ 是 $n\times n$ 矩阵，若对任何非零实矢量 $z$，都有 $z^TMz>0$，其中 $z^T$ 表示 $z$ 的转置，则称 $M$ 为正定矩阵。

② 准确地说，定义在区间 $a\le t\le b$ 内的权函数 $\rho(t)$，必须是非负可积函数，且 $\int_a^b\rho(t)\,dt>0$。以下出现的权函数均需满足此要求，不再一一指明。

③ 如果当且仅当 $\alpha_1=\alpha_2=\cdots=\alpha_n=0$ 时 $n$ 个矢量 $\{y_1,y_2,y_3,\cdots,y_n\}$ 的线性组合 $\sum_{i=1}^n\alpha_i y_i=0$，则称这 $n$ 个矢量 $\{y_1,y_2,y_3,\cdots,y_n\}$ 线性无关。无穷个矢量的线性无关则定义为，从其中任意挑选有限个矢量，全部都线性无关。

## *§17.2 函数空间

函数空间是一类特殊的线性空间：空间的元素是函数，更确切地说，本课程中讨论的是定义在一定区间（可以是有界区间、无界区间或半无界区间，为确定起见，以下不妨写为 $a\le x\le b$）上的复值平方可积①函数 $f(x)$。定义元素 $f_1$ 和 $f_2$ 的加法 $f_1+f_2$ 就是两函数值相加，

$$
(f_1+f_2)(x)
=
f_1(x)+f_2(x),
\tag{17.7}
$$

元素 $f$ 和复数 $\alpha$ 的数乘 $\alpha f$ 是

$$
(\alpha f)(x)
=
\alpha f(x),
\tag{17.8}
$$

可以证明，两个平方可积函数之和仍是平方可积的，这是因为

$$
|f_1(x)+f_2(x)|^2
+
|f_1(x)-f_2(x)|^2
=
2[
|f_1(x)|^2
+
|f_2(x)|^2
],
$$

所以

$$
|f_1(x)+f_2(x)|^2
\le
2[
|f_1(x)|^2
+
|f_2(x)|^2
].
$$

因此全体定义在 $a\le x\le b$ 上的复值平方可积函数组成的集合，对于加法 (17.7) 和数乘 (17.8) 是封闭的，构成一个复数域上的线性空间，记为 $\mathcal H$。

**定义 17.4**　设 $f_1(x)$ 和 $f_2(x)$ 是函数空间 $\mathcal H$ 中的两个函数，可以定义它们的内积为

$$
(f_1,f_2)
=
\int_a^b
f_1^*(x)f_2(x)\rho(x)\,dx,
\tag{17.9}
$$

其中 $\rho(x)$ 就是函数空间 $\mathcal H$ 中定义函数平方可积概念时给定的权函数。

由于

$$
|f_1(x)|^2+|f_2(x)|^2-2|f_1(x)|\,|f_2(x)|
=
[
|f_1(x)|-|f_2(x)|
]^2
\ge0,
$$

因此

$$
|f_1^*(x)f_2(x)|
=
|f_1(x)|\,|f_2(x)|
\le
\frac12
[
|f_1(x)|^2+|f_2(x)|^2
],
$$

所以，积分

$$
\int_a^b
|f_1^*(x)f_2(x)|\rho(x)\,dx
$$

存在。又因为

$$
\left|
\int_a^b
f_1^*(x)f_2(x)\rho(x)\,dx
\right|
\le
\int_a^b
|f_1^*(x)f_2(x)|\rho(x)\,dx,
$$

所以，只要 $f_1(x)$ 和 $f_2(x)$ 平方可积，则它们的内积也一定存在。

回顾一下内积定义中的三条要求，前两条显然满足。而且，对于空间中的任意函数 $f(x)$，恒有 $(f,f)\ge0$。现在的问题是：在这样的内积定义下，如果 $(f,f)=0$，$f(x)$ 并不见得在整个区间上处处为 0。准确地说，如果 $(f,f)=0$，则 $f(x)$ 可以在测度为零的点集上取非零值。我们把这类（几乎处处为 0 的）函数称为（广义的）零函数。在这种零函数的意义下，(17.9) 式定义的内积也就符合内积定义中的第 3 条要求。因此 (17.9) 式确实可以作为函数空间 $\mathcal H$ 上的内积定义。定义了内积 (17.9) 的函数空间 $\mathcal H$ 就成为一个内积空间。

利用内积，可以定义函数 $f(x)$ 的“长度”

$$
\|f\|
=
(f,f)^{1/2},
\tag{17.10}
$$

称为函数 $f(x)$ 的范数。

有了内积，就可以定义函数的正交性与归一性，以及正交归一函数集合的概念。

若函数 $f(x)$ 和 $g(x)$ 满足

$$
(f,g)
\equiv
\int_a^b
f^*(x)g(x)\rho(x)\,dx
=
0,
\tag{17.11}
$$

则称它们在区间 $[a,b]$ 上以权函数 $\rho(x)$ 正交。若函数 $f(x)$ 和它自身的内积

$$
(f,f)
\equiv
\int_a^b
f^*(x)f(x)\rho(x)\,dx
=
1,
\qquad
\text{亦即}
\qquad
\|f\|
=
1,
\tag{17.12}
$$

则称 $f(x)$ 是归一化的。而若对于函数集合 $\{f_i\}$，对于所有的正整数 $i$ 和 $j$，恒有

$$
(f_i,f_j)
\equiv
\int_a^b
f_i^*(x)f_j(x)\rho(x)\,dx
=
\delta_{ij},
\tag{17.13}
$$

则称此函数集合在 $[a,b]$ 上是正交归一的。

**例 17.3**　函数集合

$$
\left\{
\frac{e^{inx}}{\sqrt{2\pi}},
\quad
n=0,\pm1,\pm2,\cdots
\right\}
$$

在 $[-\pi,\pi]$ 上是正交归一的。证明请读者完成。

有了长度，可以很自然地给出两个函数之间距离的定义：内积空间 $\mathcal H$ 内任意两个函数 $f(x)$ 和 $g(x)$ 之间的距离为 $\|f-g\|$。而在距离概念的基础上，就可以进一步定义内积空间 $\mathcal H$ 内函数的邻域概念，进而定义函数序列的收敛性。

对于由内积空间 $\mathcal H$ 内的函数组成的任意 Cauchy 序列，如果它们的聚点仍保持在空间 $\mathcal H$ 内，则称该空间完备。可以证明，在内积定义 (17.9) 下，具有加法 (17.7) 和数乘 (17.8) 运算、由全体定义在 $a\le x\le b$ 上的复值平方可积函数组成的酉空间 $\mathcal H$ 是完备的。通常，把完备的内积空间称为 Hilbert 空间。可以证明上面定义的酉空间 $\mathcal H$ 还是可分的（空间 $\mathcal H$ 存在可数基）。可分的 Hilbert 空间，在物理学中有广泛的应用。下面的讨论，实际上都是在（可分的）Hilbert 空间 $\mathcal H$ 内进行的。

关于正交归一函数集的完备性概念，总是和 $\mathcal H$ 内任意函数是否可以按该函数集展开相联系的。如果对于函数空间 $\mathcal H$ 中的任意函数 $f(x)$，总可表示成正交归一函数集 $\{f_i,\ i=1,2,3,\cdots\}$① 的线性组合

$$
f(x)
=
\sum_{i=1}^{\infty}
c_i f_i(x),
\tag{17.14}
$$

则称该正交归一函数集 $\{f_i,\ i=1,2,3,\cdots\}$ 是完备的。它们就构成了函数空间 $\mathcal H$ 的一组基。

这里我们看到，第一，在完备的正交归一函数集中，函数的数目（即空间的维数）可以是有限的（空间的维数有限），也可以是无穷（空间的维数为无穷）。我们更关心无穷维的函数空间。第二，若 (17.14) 式对区间 $[a,b]$ 内的每一点 $x$ 都成立，则表明对于区间 $[a,b]$ 内的每一点 $x$，级数 $\sum_{i=1}^\infty c_i f_i(x)$ 都逐点收敛于 $f(x)$。但是，由内积定义 (17.9) 式知，Hilbert 空间 $\mathcal H$ 的零元素是广义零函数，因此 (17.14) 式应理解为左右两端相差一个广义零函数，换句话说，级数 $\sum_{i=1}^\infty c_i f_i(x)$ 应理解为平均收敛于 $f(x)$，即

$$
\lim_{n\to\infty}
\int_a^b
\left|
f(x)
-
\sum_{i=1}^{n}
c_i f_i(x)
\right|^2
\rho(x)\,dx
=
0.
\tag{17.15}
$$

第三，由函数集 $\{f_i,\ i=1,2,3,\cdots\}$ 的正交归一性，可求得展开系数

$$
c_i
=
\int_a^b
f_i^*(x)f(x)\rho(x)\,dx
=
(f_i,f).
\tag{17.16}
$$

第四，可以证明

$$
\begin{aligned}
\int_a^b
\left|
f(x)
-
\sum_{i=1}^{n}
c_i f_i(x)
\right|^2
\rho(x)\,dx
&=
(f,f)
-
\sum_{i=1}^{n}
c_i^*(f_i,f)
-
\sum_{i=1}^{n}
c_i(f,f_i)
+
\sum_{i=1}^{n}
|c_i|^2\\
&=
(f,f)
-
\sum_{i=1}^{n}
|c_i|^2,
\end{aligned}
$$

因此，只要函数集 $\{f_i,\ i=1,2,3,\cdots\}$ 是完备的，那么，由 (17.15) 式，就可以得到函数集 $\{f_i,\ i=1,2,3,\cdots\}$ 的完备性关系——Parseval 方程：

$$
(f,f)
=
\sum_{i=1}^{\infty}|c_i|^2
=
\sum_{i=1}^{\infty}|(f_i,f)|^2,
\qquad
\forall f\in\mathcal H.
\tag{17.17}
$$

① 这里隐含了此函数集为可数集的假设。所谓集合是可数的，指的是它的元素能与正整数或正整数的某个子集一一对应。只有有限个元素的集合一定是可数的。有理数构成可数集，而实数集则否。

将 (17.16) 代入 (17.14)，有

$$
f(x)
=
\sum_{i=1}^{\infty}
\int_a^b
f_i^*(x')f(x')\rho(x')\,dx'\,
f_i(x)
=
\int_a^b
\left[
\sum_{i=1}^{\infty}
f_i(x)f_i^*(x')\rho(x')
\right]
f(x')\,dx'.
$$

此结果对于 Hilbert 空间 $\mathcal H$ 中的任意函数 $f(x)$ 均成立，所以函数集 $\{f_i,\ i=1,2,\cdots\}$ 的完备性还可表示成

$$
\sum_{i=1}^{\infty}
f_i(x)f_i^*(x')\rho(x')
=
\delta(x-x'),
\qquad
a\le x,x'\le b.
\tag{17.18}
$$

在此基础上，又可以得到（更普遍的 Parseval 方程）

$$
(f,g)
=
\sum_{i=1}^{\infty}
(f,f_i)(f_i,g),
\qquad
\forall f\in\mathcal H,\ \forall g\in\mathcal H.
\tag{17.19}
$$

现在把函数集 $\{f_i,\ i=1,2,\cdots\}$ 的条件放宽，假设此函数集是正交归一的，但不一定完备，仍可以试图用这个函数集的线性组合 $\sum_{i=1}^{n}a_i f_i(x)$ 来逼近 $f(x)$。现在的问题是：如何选择组合系数 $a_i$（与 $n$ 无关），可以得到最佳逼近，使误差

$$
\left\|
f(x)
-
\sum_{i=1}^{n}
a_i f_i(x)
\right\|^2
\equiv
\int_a^b
\left|
f(x)
-
\sum_{i=1}^{n}
a_i f_i(x)
\right|^2
\rho(x)\,dx
\tag{17.20}
$$

取极小？仿照 (17.17) 式的证明，可以求得

$$
\begin{aligned}
\int_a^b
\left|
f(x)
-
\sum_{i=1}^{n}
a_i f_i(x)
\right|^2
\rho(x)\,dx
&=
(f,f)
-
\sum_{i=1}^{n}
a_i^*(f_i,f)
-
\sum_{i=1}^{n}
a_i(f,f_i)
+
\sum_{i=1}^{n}|a_i|^2\\
&=
(f,f)
-
\sum_{i=1}^{n}a_i^*c_i
-
\sum_{i=1}^{n}a_i c_i^*
+
\sum_{i=1}^{n}a_i^*a_i\\
&=
(f,f)
+
\sum_{i=1}^{n}|a_i-c_i|^2
-
\sum_{i=1}^{n}c_i^*c_i.
\end{aligned}
$$

因此，当 $a_i=c_i\equiv(f_i,f)$ 时，误差一定取极小值，

$$
(f,f)
-
\sum_{i=1}^{n}|c_i|^2
\ge0,
$$

而且，随着项数 $n$ 的增加，误差越来越小。但无论如何，总有

$$
(f,f)
\ge
\sum_{i=1}^{\infty}|c_i|^2.
\tag{17.21}
$$

这正好是函数空间中的 Bessel 不等式。等号对应于函数集是完备的情形。

**练习 17.1**　证明函数空间中的 Schwarz 不等式

$$
|(f,g)|
\le
\|f\|\,\|g\|.
$$

① 这里的平方可积定义为积分 $\int_a^b|f(x)|^2\rho(x)\,dx$ 存在，其中 $\rho(x)$ 已知，是定义在 $a\le x\le b$ 区间上的权函数。


## §17.3 自伴算符的本征值问题

上一节已经证明全体定义在 $a\le x\le b$ 上的复值平方可积的函数组成的集合，对于加法 (17.7) 和数乘 (17.8) 是封闭的，构成一个复数域上的线性空间，记为 $\mathcal H$。在 $\mathcal H$ 上定义内积 (17.9)，则 $\mathcal H$ 成为一个内积空间。设 $V\subseteq\mathcal H$，如果对于 $V$ 内的每一个函数 $f$，都有唯一一个函数 $g\in\mathcal H$ 与之对应，$g$ 和 $f$ 之间的这种对应关系记为

$$
g=\hat Lf,
\qquad
f\in V,
$$

则称 $\hat L$ 为定义在 $V$ 上的算符，其中 $f$ 称为算符 $\hat L$ 的自变量函数，$g$ 称为算符 $\hat L$ 在 $f$ 点的像函数，$\mathcal H$ 的子集 $V$ 称为算符 $\hat L$ 的定义域，常记为 $D(\hat L)=V$。因此算符 $\hat L$ 是 $\mathcal H\mapsto\mathcal H$ 的映射。需要注意的是，要完全定义一个算符，对应关系和定义域缺一不可。不明确指明定义域时，默认为在 $\mathcal H$ 内可取到的最大子集。

物理上常用的是线性算符：如果对 $\forall\alpha\in\mathbb C,\ \forall\beta\in\mathbb C,\ \forall u\in D(\hat L),\ \forall v\in D(\hat L)$，算符 $\hat L$ 满足

$$
\hat L(\alpha u+\beta v)
=
\alpha\hat Lu+\beta\hat Lv,
$$

则称算符 $\hat L$ 为线性算符。本书中也只讨论线性算符。

在物理学中，尤其量子力学中有重要应用的是自伴算符。为此，我们必须先介绍伴算符。

**定义 17.5**　$D(\hat L)\subseteq\mathcal H$，设 $\hat L$ 为定义在 $D(\hat L)$ 上的线性（微分）算符，$\forall u\in D(\hat L)$，如果对一个 $v\in\mathcal H$，能找到唯一一个 $w\in\mathcal H$ 满足

$$
(v,\hat Lu)
=
(w,u),
\qquad
\text{即}
\qquad
\int_a^b
v^*(\hat Lu)\rho(x)\,dx
=
\int_a^b
w^*u\rho(x)\,dx,
\tag{17.22}
$$

则建立了函数 $v$ 和函数 $w$ 之间的对应关系，记这种对应关系为

$$
w=\hat Mv,
$$

相应地，$\mathcal H$ 内所有 $v$ 的集合就构成了 $\hat M$ 的定义域。这样就完全确定了一个算符 $\hat M$。算符 $\hat M$ 称为原算符 $\hat L$ 的伴算符。$\hat L$ 的伴算符常记为 $\hat L^\dagger$。

从定义可以看出，不是每一个算符都一定有伴算符。即使算符 $\hat L$ 有伴算符，其伴算符 $\hat M$ 的定义域也可能不同于 $D(\hat L)$。

**例 17.4**　若 $\hat L$ 为微分算符

$$
\begin{cases}
\displaystyle
\hat Lu
=
\frac1{\rho(x)}
\frac{du}{dx},
\\[8pt]
D(\hat L)
=
\{u(x):u\in\mathcal H,\ \hat Lu\in\mathcal H,\ \text{且 }u(a)=u(b)\},
\end{cases}
$$

其中 $\rho(x)$ 为权函数，于是，由

$$
\int_a^b
v^*
\frac1{\rho(x)}
\frac{du}{dx}
\rho(x)\,dx
=
\int_a^b
v^*\frac{du}{dx}\,dx
=
v^*u\bigg|_a^b
+
\int_a^b
\left(
-\frac{dv}{dx}
\right)^*
u\,dx
$$

知 $\hat L$ 的伴算符是

$$
\begin{cases}
\displaystyle
\hat L^\dagger v
=
-\frac1{\rho(x)}
\frac{dv}{dx},
\\[8pt]
D(\hat L^\dagger)
=
\{v(x):v\in\mathcal H,\ \hat L^\dagger v\in\mathcal H,\ \text{且 }v(a)=v(b)\}.
\end{cases}
$$

**例 17.5**　若 $\hat L$ 为微分算符

$$
\begin{cases}
\displaystyle
\hat Lu
=
\frac{i}{\rho(x)}
\frac{du}{dx},
\\[8pt]
D(\hat L)
=
\{u(x):u\in\mathcal H,\ \hat Lu\in\mathcal H,\ \text{且 }u(a)=u(b)=0\},
\end{cases}
$$

其中 $\rho(x)$ 为权函数，类似地，由

$$
\int_a^b
v^*
\frac{i}{\rho(x)}
\frac{du}{dx}
\rho(x)\,dx
=
\int_a^b
v^*i\frac{du}{dx}\,dx
=
iv^*u\bigg|_a^b
+
\int_a^b
\left(
i\frac{dv}{dx}
\right)^*
u\,dx
$$

知 $\hat L$ 的伴算符是

$$
\begin{cases}
\displaystyle
\hat L^\dagger v
=
\frac{i}{\rho(x)}
\frac{dv}{dx},
\\[8pt]
D(\hat L^\dagger)
=
\{v(x):v\in\mathcal H,\ \hat L^\dagger v\in\mathcal H,\ \text{且 }v(a)\text{ 有界},\ v(b)\text{ 有界}\}.
\end{cases}
$$

上面两个例子中的伴算符 $\hat L^\dagger$ 都不同于原来的算符 $\hat L$：例 17.4 中伴算符的对应关系和原算符不同，例 17.5 中伴算符和原算符的定义域不同。但也有一些算符，不仅有伴算符，而且伴算符和原来的算符相同，不仅对应关系相同，而且定义域也完全相同，即

$$
\hat L^\dagger
=
\hat L.
\tag{17.23}
$$

**定义 17.6**　若算符 $\hat L$ 的伴算符就是它自身，则称 $\hat L$ 是自伴算符。

**例 17.6**　若 $\hat L$ 为微分算符

$$
\begin{cases}
\displaystyle
\hat Lu
=
\frac1{\rho(x)}
\frac{d^2u}{dx^2},
\\[8pt]
D(\hat L)
=
\{u(x):u\in\mathcal H,\ \hat Lu\in\mathcal H,\ \text{且 }u(a)=u(b),\ u'(a)=u'(b)\},
\end{cases}
$$

其中 $\rho(x)$ 为权函数，于是，由

$$
\int_a^b
v^*
\frac{d^2u}{dx^2}\,dx
=
\left[
v^*u'-(v^*)'u
\right]_a^b
+
\int_a^b
\left(
\frac{d^2v}{dx^2}
\right)^*
u\,dx
$$

可知，

$$
\begin{cases}
\displaystyle
\hat L^\dagger v
=
\frac1{\rho(x)}
\frac{d^2v}{dx^2},
\\[8pt]
D(\hat L^\dagger)
=
\{v(x):v\in\mathcal H,\ \hat L^\dagger v\in\mathcal H,\ \text{且 }v(a)=v(b),\ v'(a)=v'(b)\}.
\end{cases}
$$

所以这个算符是自伴算符。

**例 17.7**　若 $\hat L$ 为微分算符

$$
\begin{cases}
\displaystyle
\hat Lu
=
\frac{i}{\rho(x)}
\frac{du}{dx},
\\[8pt]
D(\hat L)
=
\{u(x):u\in\mathcal H,\ \hat Lu\in\mathcal H,\ \text{且 }u(a)=u(b)\},
\end{cases}
$$

其中 $\rho(x)$ 为权函数，根据

$$
\int_a^b
v^*i\frac{du}{dx}\,dx
=
iv^*u\bigg|_a^b
+
\int_a^b
\left(
i\frac{dv}{dx}
\right)^*
u\,dx,
$$

我们也能看出，$\hat L$ 的伴算符是

$$
\begin{cases}
\displaystyle
\hat L^\dagger v
=
\frac{i}{\rho(x)}
\frac{dv}{dx},
\\[8pt]
D(\hat L^\dagger)
=
\{v(x):v\in\mathcal H,\ \hat L^\dagger v\in\mathcal H,\ \text{且 }v(a)=v(b)\}.
\end{cases}
$$

所以这个算符也是自伴算符。

自伴算符的概念容易与 Hermite 算符混淆。

**定义 17.7**　如果对 $\forall u\in D(\hat L),\ \forall v\in D(\hat L)$，均有

$$
(u,\hat Lv)
=
(\hat Lu,v),
\tag{17.24}
$$

则称算符 $\hat L$ 为 Hermite 算符。

由 Hermite 算符的定义可知，一个算符是不是 Hermite 算符，与其伴算符甚至有无伴算符都没有关系，只要给出了算符的对应关系和定义域，我们就可以判别它是不是 Hermite 算符。自伴算符比 Hermite 算符的要求更强。从定义可以看出，自伴算符一定是 Hermite 算符，但 Hermite 算符不一定是自伴算符。例 17.5 中定义的算符不是自伴算符，但是简单验算一下就知道它是 Hermite 算符。

**例 17.8**　设

$$
\begin{cases}
\displaystyle
\hat Lu
=
\frac{i}{\rho(x)}
\frac{du}{dx},
\\[8pt]
D(\hat L)
=
\{u(x):u\in\mathcal H,\ \hat Lu\in\mathcal H,\ \text{且 }u(b)=\alpha u(a),\ \alpha\in\mathbb C\},
\end{cases}
$$

其中 $\rho(x)$ 为权函数，由

$$
\begin{aligned}
\int_a^b
v^*
\left(
i\frac{du}{dx}
\right)\,dx
&=
iv^*u\bigg|_a^b
-
i\int_a^b
\frac{dv^*}{dx}u\,dx\\
&=
i[
\alpha v^*(b)-v^*(a)
]u(a)
+
\int_a^b
\left(
\frac{i}{\rho(x)}
\frac{dv}{dx}
\right)^*
u\rho(x)\,dx
\end{aligned}
$$

知

$$
\begin{cases}
\displaystyle
\hat L^\dagger v
=
\frac{i}{\rho(x)}
\frac{dv}{dx},
\\[8pt]
D(\hat L^\dagger)
=
\{v(x):v\in\mathcal H,\ \hat L^\dagger v\in\mathcal H,\ \text{且 }\alpha^*v(b)=v(a),\ \alpha\in\mathbb C\}.
\end{cases}
$$

只有边界条件中的 $\alpha$ 满足 $\alpha\alpha^*=1$ 时，算符 $\hat L$ 才是自伴的。例 17.7 就是 $\alpha=1$ 时的特例。

**定义 17.8**　若 $\hat L$ 为自伴算符，则方程

$$
\hat Ly(x)
=
\lambda y(x)
\tag{17.25}
$$

称为自伴算符的本征值问题。

注意，这里本征方程的定解条件已由算符的定义域给定。只有 $\lambda$ 为一些特定值时，本征值问题才有解。称这些特定的 $\lambda$ 为本征值，相应的满足方程 (17.25) 的非零的 $y(x)$ 称为本征值 $\lambda$ 相应的本征函数。

讨论本征值问题，首先要回答解（即本征值与相应的本征函数）是否存在。就自伴算符而言，需要区别正则的与奇异的两种情形：如果算符的定义域是无界或半无界区间，或者区间的端点是方程的奇点，则此本征值问题是奇异的，否则就属于正则的本征值问题。

可以证明，在可分的 Hilbert 空间 $\mathcal H$ 内，正则的自伴算符本征值问题一定有解，而且本征值是离散的（因而构成可数集）。

对于奇异的自伴算符，其本征值问题则不一定有解，即使有解，也可能是连续谱，或者离散谱与连续谱二者兼而有之。由于篇幅限制，本书就不做进一步的讨论了。

在自伴算符本征值问题有解的前提下，求得的本征值与本征函数具有下列性质。

**性质 1　自伴算符的本征值必为实数。**

**证**　因为

$$
\hat Ly
=
\lambda y,
$$

取复共轭

$$
(\hat Ly)^*
=
\lambda^*y^*.
$$

由于 $\hat L$ 是自伴算符，$\hat L^\dagger=\hat L$，所以

$$
\begin{aligned}
0
&=
\int_a^b
[
y^*\hat Ly-(\hat L^\dagger y)^*y
]
\rho(x)\,dx\\
&=
\int_a^b
[
y^*\hat Ly-(\hat Ly)^*y
]
\rho(x)\,dx\\
&=
(\lambda-\lambda^*)
\int_a^b
y^*y\rho(x)\,dx.
\end{aligned}
$$

又因为

$$
\int_a^b
y^*y\rho(x)\,dx
\ne0,
$$

所以

$$
\lambda=\lambda^*,
$$

即证得本征值 $\lambda$ 为实数。□

**性质 2　自伴算符的本征函数具有正交性，即对应不同本征值的本征函数一定正交。**

**证**　设 $\lambda_i$ 和 $\lambda_j$ 是不相等的两个本征值①，对应的本征函数为 $y_i$ 和 $y_j$，

$$
\hat Ly_i
=
\lambda_i y_i,
\qquad
\hat Ly_j
=
\lambda_j y_j.
$$

注意到本征值 $\lambda_i,\lambda_j$ 为实数，于是

$$
\int_a^b
[
y_i^*\hat Ly_j-(\hat Ly_i)^*y_j
]
\rho(x)\,dx
=
(\lambda_j-\lambda_i)
\int_a^b
y_i^*y_j\rho(x)\,dx.
$$

当 $\lambda_i\ne\lambda_j$ 时，一定有

$$
\int_a^b
y_i^*(x)y_j(x)\rho(x)\,dx
=
0.
\tag{17.26}
$$

这样就证明了对应不同本征值的本征函数一定正交。□

① 这种表述相当于认定本征值是可数的，否则需要稍加改写。

因为 $\hat L$ 是线性算符，所以本征函数乘以一个非零常数仍然是本征函数。如果本征函数平方可积，则适当选择这个常数，总可以使得对于任意一个本征值 $\lambda_i$，相应的本征函数都有

$$
\int_a^b
y_i^*(x)y_i(x)\rho(x)\,dx
=
1.
\tag{17.27}
$$

这样得到的就是一个正交归一的函数组。而且，(17.26) 和 (17.27) 两式还可以统一写成

$$
\int_a^b
y_i^*(x)y_j(x)\rho(x)\,dx
=
\delta_{ij}.
\tag{17.28}
$$

如果本征函数不平方可积（这对应于本征值构成连续谱的情形），则应当将本征函数归一化到 $\delta$ 函数。

**性质 3　自伴算符 $\hat L$ 的全体本征函数构成一个完备函数组，即 $\mathcal H$ 内的任意函数 $f(x)$，均可按本征函数 $\{y_n(x)\}$① 展开为绝对而且一致收敛的级数**

$$
f(x)
=
\sum_{n=1}^{\infty}
c_n y_n(x),
\qquad
\forall f\in\mathcal H,
\tag{17.29}
$$

其中

$$
c_n
=
\frac{
\displaystyle
\int_a^b
y_n^*(x)f(x)\rho(x)\,dx
}{
\displaystyle
\int_a^b
y_n^*(x)y_n(x)\rho(x)\,dx
}.
\tag{17.30}
$$

特别是，如果本征函数是归一化的，(17.30) 式中的分母为 1，展开形式更加简单。

同样，正交归一的本征函数组的完备性也还可以表示成

$$
\sum_{n=1}^{\infty}
y_n(x)y_n^*(x')\rho(x')
=
\delta(x-x'),
\qquad
a\le x,x'\le b.
\tag{17.31}
$$

它就是 (17.18) 式，也可以直接由 (17.29) 和 (17.30) 式推出。

由上面的性质 2 和性质 3 可以看到，只要将本征函数归一化，则本征函数的全体就构成了 $\mathcal H$ 内一个完备的正交归一函数集。§17.2 中有关完备的正交归一函数集的讨论均可适用。这里我们默认对应一个本征值只有一个本征函数。如果对应于一个本征值有不止一个（线性无关的）本征函数，这些对应同一个本征值的不同本征函数可能并不彼此正交。这种情形将在 §17.5 讨论。但即使如此，总可以采用 Schmidt 的正交化步骤（见 §17.1，读者可以方便地把它移植到函数空间中来）使之正交化，因而仍然可以得到一个完备的正交归一函数集。

由于 Hilbert 空间 $\mathcal H$ 中的零元素是广义零函数，因此对于任意定义在 $[a,b]$ 中平方可积的函数，(17.29) 及 (17.31) 式都是在平均收敛

$$
\lim_{N\to\infty}
\int_a^b
\left|
f(x)
-
\sum_{n=1}^{N}
c_n y_n(x)
\right|^2
\rho(x)\,dx
=
0
\tag{17.32}
$$

的意义下成立。

① 这同样假设本征值是离散的。对于连续谱的情形，需将级数改为积分。

最后，需要指出，算符的自伴性也不是该算符的本征值问题有解的必要条件。因为，既然自伴算符 $\hat L$ 本征值一定存在，全是实数，而 $i\hat L$ 一定不是自伴算符，但后者的本征值也一定存在，就是前者的本征值乘以 $i$，它们全部都是纯虚数。

## §17.4 Sturm-Liouville 型方程的本征值问题

在前面几章中，我们讨论过几个常微分方程的本征值问题。涉及的微分方程有

$$
X''+\lambda X=0;
$$

$$
\frac d{dx}
\left[
(1-x^2)\frac{dy}{dx}
\right]
+
\left(
\lambda-\frac{m^2}{1-x^2}
\right)y
=
0;
$$

$$
\frac1r
\frac d{dr}
\left(
r\frac{dR}{dr}
\right)
+
\left(
\lambda-\frac{m^2}{r^2}
\right)R
=
0;
$$

$$
\frac1{r^2}
\frac d{dr}
\left(
r^2\frac{dR}{dr}
\right)
+
\left[
\lambda-\frac{l(l+1)}{r^2}
\right]R
=
0.
$$

它们可以归纳为 Sturm-Liouville 型（简称 S-L 型）方程的一般形式：

$$
\frac d{dx}
\left[
p(x)\frac{dy}{dx}
\right]
+
[
\lambda\rho(x)-q(x)
]y
=
0.
\tag{17.33}
$$

不妨把 S-L 型方程中的函数 $p(x),q(x)$ 和 $\rho(x)$ 均限制为实值函数，（在给定的区间 $[a,b]$ 中）都满足必要的连续性要求。其中 $p(x)\ge0,\ q(x)\ge0$，且 $p(x)$ 的零点只能位于区间的端点。$\rho(x)$ 称为权函数。当 $\rho(x)=$ 常数时，即可取为 1。不恒为常数的权函数，可能来源于正交曲面坐标系的使用（这时可以从 Laplace 算符的具体表达式中追寻到权函数的来历。从根本上说，它反映了坐标长度单位是该变量的函数，不妨称之为来源于空间的几何描述的不均匀性），也可能来源于问题所涉及的物理性质的不均匀性（例如，密度分布的不均匀）。就我们所关心的物理问题而言，都满足 $\rho(x)\ge0$，且

$$
\int_a^b
\rho(x)\,dx
>
0.
$$

S-L 型方程附加上适当的（齐次）边界条件，就构成 S-L 型方程的本征值问题。只有当待定常数 $\lambda$ 取某些特定值（称为本征值）时，才有既满足 S-L 型方程、又满足（齐次）边界条件的非零解（即本征函数）。

**练习 17.2**　将二阶常微分方程

$$
a(x)y''+b(x)y'+[c(x)-\lambda d(x)]y=0
$$

化为 S-L 方程的形式。

为了书写的紧凑，可以引进记号

$$
\hat L
\equiv
\frac1{\rho(x)}
\left\{
-\frac d{dx}
\left[
p(x)\frac d{dx}
\right]
+
q(x)
\right\}.
\tag{17.34}
$$

这样，S-L 型方程 (17.33) 就可以改写成

$$
\hat Ly(x)
=
\lambda y(x).
\tag{17.35}
$$


在采用内积定义

$$
(y_1,y_2)
\equiv
\int_a^b
y_1^*(x)y_2(x)\rho(x)\,dx
\tag{17.36}
$$

的前提条件下，我们分析一下，需要加上什么边界条件，算符

$$
\begin{cases}
\displaystyle
\hat Ly
=
\frac1{\rho(x)}
\left\{
-\frac d{dx}
\left[
p(x)\frac d{dx}
\right]
+
q(x)
\right\}y,
\\[10pt]
D(\hat L)
=
\{y(x):y\in\mathcal H,\ \hat Ly\in\mathcal H,\ \text{且满足适当的边界条件}\}
\end{cases}
\tag{17.37}
$$

是自伴的。或者说，对于 $D(\hat L)$ 中的任意两个函数 $y_i(x),y_j(x)$，我们可以按照自伴算符的定义，考察它们应当满足的边界条件。

$$
\begin{aligned}
(y_i,\hat Ly_j)
&=
\int_a^b
y_i^*(x)
\frac1{\rho(x)}
\left\{
-\frac d{dx}
\left[
p(x)\frac d{dx}
\right]
+
q(x)
\right\}
y_j(x)\rho(x)\,dx\\
&=
-y_i^*(x)p(x)
\frac{dy_j(x)}{dx}
\bigg|_{x=a}^{x=b}
+
\int_a^b
\frac{dy_i^*(x)}{dx}
p(x)
\frac{dy_j(x)}{dx}
\,dx
+
\int_a^b
y_i^*(x)q(x)y_j(x)\,dx\\
&=
-p(x)y_i^*(x)
\frac{dy_j(x)}{dx}
\bigg|_{x=a}^{x=b}
+
p(x)y_j(x)
\frac{dy_i^*(x)}{dx}
\bigg|_{x=a}^{x=b}\\
&\quad
+
\int_a^b
y_j(x)
\left\{
-\frac d{dx}
\left[
p(x)\frac{dy_i^*(x)}{dx}
\right]
+
q(x)y_i^*(x)
\right\}dx\\
&=
p(x)
\left[
y_j(x)\frac{dy_i^*(x)}{dx}
-
y_i^*(x)\frac{dy_j(x)}{dx}
\right]_{x=a}^{x=b}\\
&\quad
+
\int_a^b
\left(
\frac1{\rho(x)}
\left\{
-\frac d{dx}
\left[
p(x)\frac d{dx}
\right]
+
q(x)
\right\}
y_i(x)
\right)^*
y_j(x)\rho(x)\,dx,
\end{aligned}
$$

其中用到了 $p(x),q(x),\rho(x)$ 都是实值函数。因此，只要 $y_i(x),y_j(x)$ 满足同样的边界条件，使得

$$
p(x)
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)
\bigg|_{x=a}^{x=b}
=
0
\tag{17.38}
$$

成立，则算符 (17.37) 就是自伴算符。

下面对本书中出现过的几种类型的边界条件逐一加以检验，我们将会看到，这几种类型的边界条件都符合上述要求，都能保证 (17.38) 式成立，因此，自伴算符本征值问题的一系列结论都能移植到 S-L 型方程 (17.33) 的本征值问题中。

1. 在每个端点处，$y_j\dfrac{dy_i^*}{dx}-y_i^*\dfrac{dy_j}{dx}$ 均为 0。

不妨以端点 $x=a$ 为例。若 $p(a)\ne0$，$y(x)$ 在该点满足第一、二、三类边界条件

$$
y(a)\cos\alpha-y'(a)\sin\alpha=0,
$$

其中 $0\le\alpha\le\pi/2$。因此，根据

$$
y_i^*(a)\cos\alpha-y_i^{*\,\prime}(a)\sin\alpha=0,
$$

$$
y_j(a)\cos\alpha-y_j'(a)\sin\alpha=0,
$$

由于 $\sin\alpha,\cos\alpha$ 不可能同时为 0，我们就能导出

$$
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)_{x=a}
=
0.
$$

$x=b$ 点可做类似的讨论。

2. 在端点处，$y_j\dfrac{dy_i^*}{dx}-y_i^*\dfrac{dy_j}{dx}\ne0$，但 $p(x)$ 之值为 0。

仍不妨设端点为 $x=a$。此时 $p(a)=0$，$x=a$ 是方程 (17.33) 的奇点。如果进一步限制该端点 $x=a$ 是方程 (17.33) 的正则奇点，且指标 $\rho_1\ge0,\rho_2<0$①。这时，只要采用有界条件作为边界条件，本征函数 $y_i,y_j$ 即对应于指标 $\rho_1$。因此

$$
p(x)
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)_{x=a}
=
0.
$$

这种情形，显然构成奇异的本征值问题。

3. 单独在一个端点，$p(x)\left(y_j\dfrac{dy_i^*}{dx}-y_i^*\dfrac{dy_j}{dx}\right)\ne0$。

这发生在 $p(a)=p(b)\ne0$ 时。若 $y(x)$ 在端点满足周期条件

$$
y(a)=y(b),
\qquad
y'(a)=y'(b),
$$

我们立刻就能看到，尽管

$$
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)_{x=a}
$$

和

$$
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)_{x=b}
$$

单独可能均不为 0，但是一定有

$$
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)_{x=a}
=
\left(
y_j\frac{dy_i^*}{dx}
-
y_i^*\frac{dy_j}{dx}
\right)_{x=b},
$$

因而 (17.38) 式成立。

除了上述几种类型的边界条件外，物理学中还会遇到其他类型的边界条件。例如，对于区间无界（或半无界）的情形，这时需要附加上无穷远条件，即要求 $x\to\pm\infty$（或 $+\infty$）时 $y(x)$ 具有适当的渐近行为。量子力学课程中会遇到这种类型的边界条件。也会遇到区间端点为非正则奇点的情形。这些情形也都是构成奇异的本征值问题。正如上一节指出的，这些本征值问题是否有解，需要针对具体问题做具体讨论。

① 这时除了 $p(x),q(x)$ 和 $\rho(x)$ 应满足本节开始处提到的要求外，还应要求 $p(x)$ 和 $\rho(x)$ 在 $x=a$ 点解析，且满足

$$
p(a)=0,
\qquad
p'(a)\ne0,
\qquad
(x-a)q(x)\text{ 在 }x=a\text{ 点解析}
$$

或

$$
p(a)=0,
\qquad
p'(a)=0,
\qquad
p''(a)\ne0,
\qquad
q(x)\text{ 在 }x=a\text{ 点解析},
\qquad
\rho(a)=0.
$$

这在我们讨论过的实际问题中都是能够满足的。直接计算即可导出，方程在正则奇点 $x=a$ 处的指标方程为

$$
\rho(\rho-1)+\rho-b_0=0
$$

或

$$
\rho(\rho-1)+2\rho-b_0=0,
$$

其中

$$
b_0
=
\lim_{x\to a}
(x-a)^2
\frac{\lambda\rho(x)-q(x)}{p(x)}
\ge0.
$$

**例 17.9**　求解本征值问题

$$
\begin{cases}
\displaystyle
\frac1r\frac d{dr}
\left(
r\frac{dR(r)}{dr}
\right)
+
\frac{\lambda}{r^2}R(r)
=
0,
\\[10pt]
R(a)=0,
\qquad
R(b)=0,
\end{cases}
$$

其中 $b>a>0$。

**解**　此微分方程属于 Euler 型，在第十四章及第十五章中均已遇到过。求解此方程，可作代换 $\zeta=\ln r,\ Z(\zeta)=R(r)$，于是本征值问题即化为

$$
\begin{cases}
Z''(\zeta)+\lambda Z(\zeta)=0,
\\
Z(\ln a)=0,
\qquad
Z(\ln b)=0.
\end{cases}
$$

此本征值问题也已多次见过，其解为

本征值

$$
\lambda_n
=
\left(
\frac{n\pi}{\ln b-\ln a}
\right)^2,
\qquad
n=1,2,3,\cdots,
$$

本征函数

$$
Z_n(\zeta)
=
\sin
\left(
\frac{\zeta-\ln a}{\ln b-\ln a}
n\pi
\right).
$$

回到原来的变量 $r$，本征值不变，而本征函数则为

$$
R_n(r)
=
\sin
\left(
\frac{\ln r-\ln a}{\ln b-\ln a}
n\pi
\right).
$$

容易判断，此本征函数的正交权函数为 $1/r$。

在求解第十四章第 6 题时，可以出现此本征值问题。

**例 17.10**　设有本征值问题

$$
\begin{cases}
X^{(4)}+\lambda X=0,
\\
X(0)=0,
\qquad
X(l)=0,
\\
X''(0)=0,
\qquad
X''(l)=0,
\end{cases}
$$

试证明：本征值

$$
\lambda
=
-
\frac{
\displaystyle
\int_0^l
|X''(x)|^2\,dx
}{
\displaystyle
\int_0^l
|X(x)|^2\,dx
}
<
0.
$$

**证**　将 $X(x)$ 的微分方程乘以 $X^*(x)$，积分，并代入 $X(x),X^*(x)$ 所满足的边界条件，就有

$$
\begin{aligned}
\lambda
\int_0^l
|X(x)|^2\,dx
&=
-
\int_0^l
X^{(4)}(x)X^*(x)\,dx\\
&=
\int_0^l
X^{(3)}(x)X^{*\,\prime}(x)\,dx
-
X^{(3)}(x)X^*(x)\bigg|_{x=0}^{x=l}\\
&=
-
\int_0^l
|X''(x)|^2\,dx
+
X''(x)X^{*\,\prime}(x)\bigg|_{x=0}^{x=l}\\
&=
-
\int_0^l
|X''(x)|^2\,dx.
\end{aligned}
$$

因此即证得

$$
\lambda
=
-
\frac{
\displaystyle
\int_0^l
|X''(x)|^2\,dx
}{
\displaystyle
\int_0^l
|X(x)|^2\,dx
}
<
0.
$$

实际上，我们写出微分方程 $X^{(4)}+\lambda X=0$ 的通解

$$
X(x)
=
A\sin\sqrt[4]{-\lambda}\,x
+
B\cos\sqrt[4]{-\lambda}\,x
+
C\sinh\sqrt[4]{-\lambda}\,x
+
D\cosh\sqrt[4]{-\lambda}\,x,
$$

代入边界条件，可以得到 $B=C=D=0$，取 $A=1$，即可求得

本征值

$$
\lambda_n
=
-
\left(
\frac{n\pi}{l}
\right)^4,
\qquad
n=1,2,3,\cdots,
$$

本征函数

$$
X_n(x)
=
\sin\frac{n\pi}{l}x.
$$

在这个例题中，讨论的是 4 阶常微分方程的本征值问题，应该说超出了 Sturm-Liouville 型方程的形式，但是我们可以证明，涉及的 4 阶微分算符是自伴的，对应不同本征值的本征函数也是正交的，甚至从本征函数的具体形式，也可以看出它们的确构成完备的函数组。讨论这个例题，一方面是因为在实际问题中可能会遇到，例如弹性力学中会出现这种类型的 4 阶常微分方程本征值问题。另一方面，也想提请读者注意，所谓自伴算符的本征值问题，绝不仅仅限于 Sturm-Liouville 型方程本征值问题这一种形式。我们不但会遇到 4 阶的自伴算符，前面在例 17.7 中，也遇到过一阶自伴算符。

## §17.5 Sturm-Liouville 型方程本征值问题的简并现象

在第 14 章中，我们曾经遇到过对应一个本征值有不止一个（线性无关的）本征函数的情形。这种现象称为简并或退化。由于 S-L 型方程是二阶线性常微分方程，所以，对应一个本征值最多只能有两个（线性无关的）本征函数。现在的问题是，到底在什么条件下，S-L 型方程的本征值问题是简并的，在什么条件下是非简并的，下面的两个定理可以给予明确的回答。

**定理 17.1**　如果 S-L 型方程本征值问题的本征函数是复值函数，且其实部和虚部线性无关，则此本征值问题是二重简并的。

**证**　根据定理所设，本征函数 $y(x)$ 是复函数，设其实部和虚部分别为 $f(x)$ 和 $g(x)$，

$$
y(x)
=
f(x)+ig(x);
$$

则 S-L 型方程可以写成

$$
\hat L(f+ig)
=
\lambda(f+ig).
$$

由于算符 $\hat L$ 为实算符，当 $\hat L$ 自伴时本征值 $\lambda$ 亦为实数，故将上式分别比较实部和虚部，就得到

$$
\hat Lf
=
\lambda f,
\qquad
\hat Lg
=
\lambda g.
$$

这说明 $f(x)$ 和 $g(x)$ 都是对应于同一个本征值 $\lambda$ 的本征函数。

还必须证明 $f(x)$ 和 $g(x)$ 也满足原本征值问题的边界条件。这只要注意边界条件也是线性齐次的，其中的系数也是实数，于是在边界条件中也分别比较实部和虚部即可。□

**定理 17.2**　设 $y_1(x)$ 和 $y_2(x)$ 都是 S-L 型方程本征值问题 (17.37) 的两个线性无关的实值本征函数，并且在 $x=a$ 和 $x=b$ 点都单独满足边界条件 (17.38)，则 $y_1(x)$ 和 $y_2(x)$ 不可能对应于同一个本征值 $\lambda$。

**证**　用反证法。设 $y_1(x)$ 和 $y_2(x)$ 对应于同一个本征值 $\lambda$，

$$
\hat Ly_1
=
\lambda y_1,
\qquad
\hat Ly_2
=
\lambda y_2,
$$

因此

$$
(y_1,\hat Ly_2)
-
(y_2,\hat Ly_1)
=
0.
$$

注意 $y_1(x)$ 和 $y_2(x)$ 都是实值函数，$y_1^*(x)=y_1(x),\ y_2^*(x)=y_2(x)$，将 (17.34) 式代入，即得

$$
\frac d{dx}
\left[
p(x)
\left(
y_1\frac{dy_2}{dx}
-
y_2\frac{dy_1}{dx}
\right)
\right]
=
0.
$$

于是

$$
p(x)
\left(
y_1\frac{dy_2}{dx}
-
y_2\frac{dy_1}{dx}
\right)
=
\text{常数 }C.
$$

而根据定理给出的已知条件，

$$
p(x)
\left(
y_1\frac{dy_2}{dx}
-
y_2\frac{dy_1}{dx}
\right)_{x=a}
=
0,
$$

$$
p(x)
\left(
y_1\frac{dy_2}{dx}
-
y_2\frac{dy_1}{dx}
\right)_{x=b}
=
0,
$$

所以

$$
p(x)
\left(
y_1\frac{dy_2}{dx}
-
y_2\frac{dy_1}{dx}
\right)
\equiv
0.
$$

但因为 $p(x)\ne0$，故有

$$
y_1\frac{dy_2}{dx}
-
y_2\frac{dy_1}{dx}
=
0,
$$

即

$$
W[y_1(x),y_2(x)]
\equiv
\begin{vmatrix}
y_1(x)&y_2(x)\\
y_1'(x)&y_2'(x)
\end{vmatrix}
=
0.
$$

$y_1(x)$ 和 $y_2(x)$ 线性相关，与已知条件矛盾。故 $y_1(x)$ 和 $y_2(x)$ 不可能对应于同一个本征值。□

这个定理说明，在一、二、三类（齐次）边界条件或有界条件下，S-L 型方程本征值问题不可能是简并的。就本书所讨论过的几种类型的边界条件而言，只有在周期条件之下，本征函数在区间的每一个端点并不单独满足 (17.38)，才有可能发生简并现象。

最后强调一下，这里讨论的是常微分方程的本征值问题。如果是偏微分方程的本征值问题，一般说来，即使在一、二、三类（齐次）边界条件或（和）有界条件下，也会出现简并现象。由于偏微分方程的本征值问题，原则上总可以通过分离变量而化为若干个常微分方程的本征值问题，因此本书不拟做更多的讨论。

## §17.6 从 Sturm-Liouville 型方程的本征值问题看分离变量法

仍以弦的横振动问题为例。对于两端固定弦的自由振动，定解问题是

$$
\frac{\partial^2u}{\partial t^2}
-
a^2
\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0;
\tag{17.39a}
$$

$$
\left.u\right|_{x=0}
=
0,
\qquad
\left.u\right|_{x=l}
=
0,
\qquad
t\ge0;
\tag{17.39b}
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
0\le x\le l.
\tag{17.39c}
$$

根据 §17.3 和 §17.4 的讨论可知，如果存在一个自伴的 S-L 型算符

$$
\hat LX
=
\frac1{\rho(x)}
\left\{
-\frac d{dx}
\left[
p(x)\frac d{dx}
\right]
+
q(x)
\right\}X,
\tag{17.40a}
$$

$$
D(\hat L)
=
\{X(x):X\in\mathcal H,\ \hat LX\in\mathcal H,\ X(0)=0,\ X(l)=0\},
\tag{17.40b}
$$

则相应的本征值问题

$$
\hat LX
=
\lambda X
\tag{17.41}
$$

是自伴的。并且为方便起见，假设本征函数均已正交归一化：

$$
(X_n,X_m)
=
\int_0^l
X_n^*(x)X_m(x)\rho(x)\,dx
=
\delta_{nm}.
\tag{17.42}
$$

由于定解问题 (17.39) 的边界条件 (17.39b) 的形式和边界条件 (17.40b) 相同，因此可以将定解问题 (17.39) 的解 $u(x,t)$ 按照本征函数组 $\{X_n(x),n=1,2,3,\cdots\}$ 展开：

$$
u(x,t)
=
\sum_{n=1}^{\infty}
T_n(t)X_n(x).
\tag{17.43}
$$

代入方程 (17.39a)，有

$$
\sum_{m=1}^{\infty}
T_m''(t)X_m(x)
-
a^2
\sum_{m=1}^{\infty}
T_m(t)X_m''(x)
=
0.
\tag{17.44}
$$

将上式两端同乘以 $X_n^*(x)\rho(x)$，并在区间 $[0,l]$ 上积分，就得到关于 $T_n(t)$ 的常微分方程组

$$
T_n''(t)
-
a^2
\sum_{m=1}^{\infty}
(X_n,X_m'')T_m(t)
=
0,
\qquad
n=1,2,3,\cdots.
\tag{17.45}
$$

同样，再将初始条件 (17.39c) 也按这一组本征函数展开，还可以得到

$$
T_n(0)
=
(X_n,\phi),
\qquad
T_n'(0)
=
(X_n,\psi).
\tag{17.46}
$$

如果能够由 (17.45) 和 (17.46) 求出 $T_n(t)$，也就求出了定解问题 (17.39) 的解 $u(x,t)$。这里，本征函数组的完备性起了决定性的作用。为了保证级数 $\sum_{n=1}^{\infty}T_n(t)X_n(x)$ 收敛（至少是平均收敛）到解 $u(x,t)$，求和必须遍及全部本征函数，绝不可以无理由地摈弃其中的任何一部分。

弄清楚了本征函数组在分离变量法中的决定性作用后，就不难求解由非齐次方程和齐次边界条件构成的定解问题。把定解问题 (17.39) 中的方程 (17.39a) 改为

$$
\frac{\partial^2u}{\partial t^2}
-
a^2
\frac{\partial^2u}{\partial x^2}
=
f(x,t),
\tag{17.47}
$$

那么，现在看来，求解过程并没有太大的差异，不同之处只在于要将方程的非齐次项 $f(x,t)$ 也按本征函数展开，于是，齐次的常微分方程组 (17.45) 变成了非齐次的方程组

$$
T_n''(t)
-
a^2
\sum_{m=1}^{\infty}
(X_n,X_m'')T_m(t)
=
(X_n,f),
\qquad
n=1,2,3,\cdots.
\tag{17.48}
$$

当然要求 $f(x,t)$ 作为 $x$ 的函数，和 $\{X_n(x)\}$ 同属一个函数空间。

同样不难理解，如果定解问题的边界条件是非齐次的，就必须先将边界条件齐次化。

到现在为止，我们着重分析了齐次边界条件 (17.40b) 在分离变量法中的决定性作用。对于本征函数所满足的微分方程，也就是自伴算符 $\hat L$ 的具体形式 (17.40a) 并没有任何限制。在满足齐次边界条件 (17.40b) 的前提下，选择不同的本征函数组 $\{X_n(x),n=1,2,3,\cdots\}$，得到的关于 $T_n(t)$ 的常微分方程组的形式当然绝不相同，因而求得的 $T_n(t)$ 也不相同。但是，定解问题的解的存在唯一性，保证了最后求得的是同一个 $u(x,t)$。

但是，并不是任何常微分方程组都是容易求解的，何况现在遇到的是无穷维的常微分方程组。在实际的求解过程中，我们就应当恰当地选择本征函数组 $\{X_n(x),n=1,2,3,\cdots\}$，使得 $T_n(t)$ 的求解问题尽可能地简单。最简单的情形当然就是 $T_n(t)$ 满足的是常微分方程，而不是常微分方程组。在 (17.48) 式中，就应当要求本征函数满足

$$
(X_n,X_m'')
=
-\lambda_n\delta_{nm},
\tag{17.49}
$$

即本征函数满足常微分方程

$$
-X_n''(x)
=
\lambda_nX_n(x).
\tag{17.49'}
$$

这正是分离变量法中得到的微分方程。相应地，方程组 (17.45) 和 (17.48) 就变成常微分方程

$$
T_n''(t)
+
a^2\lambda_nT_n(t)
=
0
\tag{17.50}
$$

和

$$
T_n''(t)
+
a^2\lambda_nT_n(t)
=
(X_n,f).
\tag{17.51}
$$

所以，分离变量法实际上提供了一个选择本征函数组的理想方案。如果说，本征函数的完备性是在理论上保证了一定可以将定解问题的解按该本征函数组展开（这是有条件的，定解问题和本征函数要满足相同的齐次边界条件），而选用相应齐次问题的本征函数则保证了可以方便地求出展开系数（实际上是函数），保证了这种解法在实用上的可行性。

然而从解式的最终形式 (17.43) 来看，尽管级数中每一项都是分离变量的形式，但是叠加后就不再是分离变量的了。所以，与其说是用分离变量法求解偏微分方程定解问题，还不如说是用本征函数展开法求解偏微分方程定解问题来得更确切。基于这种认识，我们对于各种类型的定解问题（方程齐次或非齐次，边界条件齐次或非齐次）的求解就有了一个统一的更深入的理解，并且对于求解偏微分方程定解问题也就获得了更大的自由，这表现为拓宽了对于某些定解问题的求解思路。例如，对于单位球内的稳定问题

$$
\nabla^2u
=
f,
\qquad
x^2+y^2+z^2<1,
\tag{17.52a}
$$

$$
\left.u\right|_{x^2+y^2+z^2=1}
=
0,
\tag{17.52b}
$$

采用球坐标系求解，按照过去的做法，应当将 $u(r,\theta,\phi)$ 按相应齐次问题的本征函数展开，

$$
u(r,\theta,\phi)
=
\sum_{l=0}^{\infty}
\sum_{m=-l}^{l}
R_{lm}(r)Y_l^m(\theta,\phi),
\tag{17.53}
$$

代入方程 (17.52a)，导出 $R_{lm}(r)$ 满足的非齐次方程，结合边界条件

$$
R_{lm}(0)\text{ 有界},
\qquad
R_{lm}(1)
=
0
\tag{17.54}
$$

求出 $R_{lm}(r)$。但是，按照本节中前面的分析，如果能找到一组本征函数，只要它也满足此定解问题的齐次边界条件，那么，就可以将 $u(r,\theta,\phi)$ 按这一组本征函数展开。具体说来，可以先求解本征值问题

$$
-\nabla^2w
=
\lambda w,
\qquad
x^2+y^2+z^2<1,
\tag{17.55a}
$$

$$
\left.w\right|_{x^2+y^2+z^2=1}
=
0,
\tag{17.55b}
$$

得到本征值①和本征函数

$$
\lambda_{nl}
=
k_{nl}^2,
\qquad
n=1,2,3,\cdots,\quad
l=0,1,2,\cdots,
\tag{17.56}
$$

$$
w_{nlm}(r,\theta,\phi)
=
j_l(k_{nl}r)
Y_l^m(\theta,\phi),
\qquad
m=0,\pm1,\cdots,\pm l,
\tag{17.57}
$$

$k_{nl}$ 是 $l$ 阶球 Bessel 函数 $j_l(x)$ 的第 $n$ 个正零点。然后将 $u(r,\theta,\phi)$ 按 $w_{nlm}(r,\theta,\phi)$ 展开，

$$
u(r,\theta,\phi)
=
\sum_{n=1}^{\infty}
\sum_{l=0}^{\infty}
\sum_{m=-l}^{l}
c_{nlm}
j_l(k_{nl}r)
Y_l^m(\theta,\phi),
\tag{17.58}
$$

代入方程 (17.52a)，就得到

$$
-k_{nl}^2c_{nlm}
\int_0^1
j_l^2(k_{nl}r)r^2\,dr
=
\int_0^1
j_l(k_{nl}r)r^2\,dr
\iint
Y_l^{m*}(\theta,\phi)
f(r,\theta,\phi)
\sin\theta\,d\theta\,d\phi.
$$

根据球 Bessel 函数的定义 (16.99) 及 (16.46) 式的结果，可以求得

$$
\begin{aligned}
\int_0^1
j_l^2(k_{nl}r)r^2\,dr
&=
\frac{\pi}{2k_{nl}}
\int_0^1
J_{l+1/2}^2(k_{nl}r)r\,dr\\
&=
\frac{\pi}{4k_{nl}}
[J_{l+1/2}'(k_{nl})]^2
=
\frac12
[j_l'(k_{nl})]^2,
\end{aligned}
$$

所以

$$
c_{nlm}
=
-\frac{2}
{k_{nl}^2[j_l'(k_{nl})]^2}
\int_0^1
j_l(k_{nl}r)r^2\,dr
\iint
Y_l^{m*}(\theta,\phi)
f(r,\theta,\phi)
\sin\theta\,d\theta\,d\phi.
\tag{17.59}
$$

前面在求解本征函数时，已经用到了全部边界条件（包括在转换到球坐标系时出现的周期条件和有界条件），这样，就自动保证了解 $u(r,\theta,\phi)$ 也满足这些边界条件。

① 注意，这里的本征值与 $m=0,\pm1,\cdots,\pm l$ 无关。换句话说，此本征值问题是对 $m$ 简并的，简并度为 $2l+1$。

这种解法的优点是除了要找到合适的本征函数外，无需再去求解常微分方程。这是以多作了一重级数展开为代价的。这相当于将 $R_{lm}(r)$ 也按球 Bessel 函数 $j_l(k_{nl}r)$ 展开而得的结果。这种做法也有局限性，它只适用于不含时间的稳定问题，并且还要求相应的本征值问题有解，甚至还要求 0 不是本征值。

以上的做法，可以方便地推广到其他几何形状的三维区域，或二维平面上的适当区域。

## 习题

1. 将下列方程化为 Sturm-Liouville 型方程的标准形式：

(1)

$$
x\frac{d^2y}{dx^2}
+
2\frac{dy}{dx}
+
(x+\lambda)y
=
0;
$$

(2)

$$
x(1-x)\frac{d^2y}{dx^2}
+
(a-bx)\frac{dy}{dx}
+
\lambda y
=
0;
$$

(3)

$$
x\frac{d^2y}{dx^2}
+
(1-x)\frac{dy}{dx}
+
\lambda y
=
0;
$$

(4)

$$
\frac{d^2y}{dx^2}
-
2x\frac{dy}{dx}
+
2\lambda y
=
0.
$$

2. 求解本征值问题：

$$
\frac1r\frac d{dr}
\left(
r\frac{dR}{dr}
\right)
+
\frac{\lambda}{r^2}R
=
0,
$$

$$
R(a)=0,
\qquad
R'(b)=0,
$$

其中 $b>a>0$。

3. 设有本征值问题

$$
\frac d{dx}
\left[
p(x)\frac{dy}{dx}
\right]
+
[
\lambda\rho(x)-q(x)
]y
=
0,
$$

$$
y(b)
=
\alpha_{11}y(a)
+
\alpha_{12}y'(a),
$$

$$
y'(b)
=
\alpha_{21}y(a)
+
\alpha_{22}y'(a),
$$

其中 $p(a)=p(b)$。试证明，当

$$
\begin{vmatrix}
\alpha_{11}&\alpha_{12}\\
\alpha_{21}&\alpha_{22}
\end{vmatrix}
=
1
$$

时，对应不同本征值的本征函数正交。

4. 求解本征值问题

$$
\begin{cases}
\displaystyle
i\frac{dy}{dx}
=
\lambda y,
&
a<x<b,
\\[8pt]
y(b)
=
e^{i\theta}y(a),
\end{cases}
$$

其中 $0\le\theta\le2\pi$。

5. 问 $\alpha,\beta$ 和 $\gamma$ 满足什么条件时，本征值问题

$$
\begin{cases}
\displaystyle
\frac{d^2u}{dx^2}
+
\cot x\frac{du}{dx}
+
\lambda u
=
0,
&
\theta<x<\pi-\theta,
\\[8pt]
u(\pi-\theta)\cos\gamma
=
u(\theta)\cos\alpha
-
u'(\theta)\sin\alpha,
\\[8pt]
u'(\pi-\theta)\sin\gamma
=
u(\theta)\cos\beta
-
u'(\theta)\sin\beta
\end{cases}
$$

是自伴的？其中 $\alpha,\beta,\gamma$ 和 $\theta$ 均为已知常数，$0<\theta<\pi/2,\ 0\le\alpha,\beta,\gamma\le\pi/2$。

6. 设本征值问题

$$
\nabla^2\Phi
+
\lambda\Phi
=
0,
$$

$$
\left.\Phi\right|_{\Sigma}
=
0
$$

的解（本征函数）为 $\Phi_k$，对应的本征值为 $\lambda_k$，这里的 $k$ 是本征值的编号。试证明：当 $\lambda=0$ 不是本征值时，Poisson 方程的第一类边值问题

$$
\nabla^2u
=
-f,
$$

$$
\left.u\right|_{\Sigma}
=
0
$$

的解为

$$
u
=
\sum_k
\frac{A_k}{\lambda_k}
\Phi_k,
$$

$A_k$ 是非齐次项 $f$ 按 $\{\Phi_k\}$ 展开的系数，

$$
f
=
\sum_k
A_k\Phi_k.
$$

这里假设 $\Phi_k$ 已归一化。

7. 用第 6 题的方法求解矩形区域 $0\le x\le a,\ 0\le y\le b$ 内 Poisson 方程的定解问题

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
-f(x,y),
$$

$$
\left.u\right|_{x=0}
=
0,
\qquad
\left.u\right|_{x=a}
=
0,
$$

$$
\left.u\right|_{y=0}
=
0,
\qquad
\left.u\right|_{y=b}
=
0.
$$

8. 非自伴算符的本征值问题举例。

以长为 $l$ 的均匀圆杆作微小扭转振动。在振动过程中，杆的各横截面仍保持为平面而绕杆轴扭转，轴向上不发生位移。杆的一端固定，另一端连接在圆盘上，则偏转角 $\theta$ 所满足的方程和边界条件为

$$
\frac{\partial^2\theta}{\partial t^2}
-
a^2
\frac{\partial^2\theta}{\partial x^2}
=
0,
$$

$$
\left.\theta\right|_{x=0}
=
0,
\qquad
\left.
\frac{\partial^2\theta}{\partial t^2}
\right|_{x=l}
=
-c^2
\left.
\frac{\partial\theta}{\partial x}
\right|_{x=l},
$$

已知 $a$ 和 $c$ 均为正数。

(1) 求相应的本征值 $\lambda_n$ 及本征函数 $X_n(x)$；

(2) 计算积分

$$
\int_0^l
X_n(x)X_m(x)\,dx;
$$

(3) 计算积分

$$
\int_0^l
X_n'(x)X_m'(x)\,dx.
$$
