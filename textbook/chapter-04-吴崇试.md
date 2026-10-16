# 第四章 无穷级数

无穷级数，特别是幂级数，是解析函数的重要表达形式之一。事实上，除了代数函数①，许多初等函数和特殊函数都是用幂级数定义的。

在复变函数理论中，无穷级数的许多基本概念，和高等数学中的实数级数完全相似。对于这部分内容，我们将不加证明地叙述一下有关结论。请读者在注意表述的相似性的同时，更要关注内涵上可能存在的差异。

## §4.1 复数级数

给定复数级数②

$$
u_0+u_1+u_2+\cdots+u_n+\cdots=\sum_{n=0}^{\infty}u_n,\tag{4.1}$$

如果它的部分和

$$S_n=u_0+u_1+u_2+\cdots+u_n\tag{4.2}$$

所构成的序列 $\{S_n\}$ 收敛，则称级数 $\sum u_n$ 收敛，而序列 $\{S_n\}$ 的极限 $S=\lim\limits_{n\to\infty}S_n$，称为级数 $\sum u_n$ 的和；否则，级数 $\sum u_n$ 是发散的。级数的收敛性，完全等价于其部分和序列的收敛性。因此，根据序列收敛的充要条件，可以写出无穷级数收敛的 Cauchy 充要条件：$\forall\varepsilon>0$，$\exists$ 正整数 $n$，使对于任意正整数 $p$，有

$$\left|u_{n+1}+u_{n+2}+\cdots+u_{n+p}\right|<\varepsilon.\tag{4.3}$$

特别是，令 $p=1$，就得到级数收敛的必要条件

$$\lim_{n\to\infty}u_n=0.\tag{4.4}$$

一个收敛序列的子序列一定收敛。因此只要不改变求和次序，可将收敛级数并项，也就是说可以给收敛级数任意添加括号。例如

$$u_1+u_2+u_3+u_4+\cdots=(u_1+u_2)+(u_3+u_4)+\cdots.\tag{4.5}$$

但是不能随意去掉收敛级数中的括号。例如

$$(1-1)+(1-1)+(1-1)+\cdots\ne 1-1+1-1+1-1+\cdots.$$

① 代数函数可用代数方程的根定义，可参考 https://en.wikipedia.org/wiki/Algebraic_function。

② 令 $u_n$ 的实部和虚部分别为 $\alpha_n$ 与 $\beta_n$，则

$$\sum u_n=\sum\alpha_n+i\sum\beta_n.$$

因此，一个复数级数 $\sum u_n$ 完全等价于两个实数级数 $\sum\alpha_n$ 和 $\sum\beta_n$，反之亦然。

如果级数 $\sum_{n=0}^{\infty}|u_n|$ 收敛，则称级数 $\sum_{n=0}^{\infty}u_n$ **绝对收敛**。因为

$$|u_{n+1}+u_{n+2}+\cdots+u_{n+p}|\le |u_{n+1}|+|u_{n+2}|+\cdots+|u_{n+p}|,$$

所以绝对收敛的级数一定收敛。反之，收敛级数不一定绝对收敛。

由于 $\sum_{n=0}^{\infty}|u_n|$ 是实数级数，而且是正项级数，所以高等数学中的任何一种正项级数的收敛判别法都可以用来判别复数级数是否绝对收敛。下面列出最常用的几个判别法。

**比较判别法** 若 $\exists N\in\mathbb N$，对 $\forall n>N$，都有 $|u_n|<v_n$，而 $\sum_{n=0}^{\infty}v_n$ 收敛，则 $\sum_{n=0}^{\infty}|u_n|$ 收敛，即 $\sum_{n=0}^{\infty}u_n$ 绝对收敛。若 $|u_n|>v_n>0$，而 $\sum_{n=0}^{\infty}v_n$ 发散，则 $\sum_{n=0}^{\infty}|u_n|$ 发散。

**比值判别法** 若存在与 $n$ 无关的常数 $\rho$，则当 $|u_{n+1}/u_n|<\rho<1$ 时，级数 $\sum_{n=0}^{\infty}u_n$ 绝对收敛；当 $|u_{n+1}/u_n|>\rho>1$ 时，级数 $\sum_{n=0}^{\infty}u_n$ 发散。

**d'Alembert 判别法** 若

$$\overline{\lim}_{n\to\infty}\left|\frac{u_{n+1}}{u_n}\right|<1,$$

则级数 $\sum_{n=0}^{\infty}u_n$ 绝对收敛；若

$$\underline{\lim}_{n\to\infty}\left|\frac{u_{n+1}}{u_n}\right|>1,$$

则 $\sum_{n=0}^{\infty}u_n$ 发散。

若 $\lim\limits_{n\to\infty}|u_{n+1}/u_n|=1$，则 $\sum_{n=0}^{\infty}u_n$ 的绝对收敛性需要利用下面的 Gauss 判别法进一步检验。

**Gauss 判别法** 设级数 $\sum_{n=0}^{\infty}u_n$ 邻项的比值可以写成

$$\frac{u_n}{u_{n+1}}=1+\frac{\mu}{n}+O(n^{-\lambda}),$$

其中 $\mu=a+ib$，$\lambda>1$，符号 $O$ 表示数量级①。若 $a>1$，则级数 $\sum_{n=0}^{\infty}u_n$ 绝对收敛；若 $a\le 1$，则 $\sum_{n=0}^{\infty}|u_n|$ 发散。

**Cauchy 判别法** 若

$$\overline{\lim}_{n\to\infty}|u_n|^{1/n}<1,$$

则级数 $\sum_{n=0}^{\infty}u_n$ 绝对收敛；若

$$\underline{\lim}_{n\to\infty}|u_n|^{1/n}>1,$$

则级数 $\sum_{n=0}^{\infty}u_n$ 发散。

绝对收敛级数具有下列性质：

(1) 改换次序。例如，

$$u_0+u_1+u_2+u_3+u_4+\cdots=u_0+u_1+u_2+u_4+u_3+u_6+u_8+u_5+\cdots.\tag{4.6}$$

(2) 特别是，可以把绝对收敛级数拆成几个子级数，每个子级数仍绝对收敛。例如，

$$\sum_{n=0}^{\infty}u_n=\sum_{n=0}^{\infty}u_{2n}+\sum_{n=0}^{\infty}u_{2n+1}.\tag{4.7}$$

① $\phi(z)=O\{\psi(z)\}$ 的含义是，对于 $z$ 的一定区域，存在常数 $A>0$，使 $|\phi|\le A|\psi|$。

(3) 两个绝对收敛级数之积仍然绝对收敛，

$$\sum_k u_k\cdot\sum_l v_l=\sum_{k,l}u_kv_l.\tag{4.8}$$

这里的乘积是一个二重级数

$$
\begin{aligned}
&u_0v_0+u_0v_1+u_0v_2+u_0v_3+\cdots\\
+{}&u_1v_0+u_1v_1+u_1v_2+u_1v_3+\cdots\\
+{}&u_2v_0+u_2v_1+u_2v_2+u_2v_3+\cdots\\
+{}&\cdots.
\end{aligned}
$$

绝对收敛性意味着可按任意顺序求和，其值不变。例如可按 $k+l=n$ 的大小顺序排列，

$$\sum_{k=0}^{\infty}u_k\cdot\sum_{l=0}^{\infty}v_l=\sum_{n=0}^{\infty}w_n,\qquad w_n=\sum_{k=0}^{n}u_kv_{n-k}.\tag{4.9}$$

而且如果限于这种求和次序，则乘法的条件还可以放宽成：$\sum u_k,\sum v_l$ 都收敛，且其中之一绝对收敛；或 $\sum u_k,\sum v_l$ 和 $\sum w_n$ 都收敛。

**练习 4.1** 设 $\sum_{n=1}^{\infty}a_n$ 与 $\sum_{n=1}^{\infty}b_n$ 皆为正项级数，试举反例，说明下列各种说法不正确：

(1) 若 $\lim\limits_{n\to\infty}na_n=0$，则 $\sum_{n=1}^{\infty}a_n$ 收敛；

(2) 若 $a_{2n}<a_{2n+1}$，则 $\sum_{n=1}^{\infty}a_n$ 发散；

(3) 若 $\lim\limits_{n\to\infty}\dfrac{a_{2n+1}}{a_{2n}}=\infty$，则 $\sum_{n=1}^{\infty}a_n$ 发散；

(4) 若 $\sum_{n=1}^{\infty}a_n$ 与 $\sum_{n=1}^{\infty}b_n$ 发散，则 $\sum_{n=1}^{\infty}\sqrt{a_nb_n}$ 发散。

## §4.2 二重级数

所谓二重级数，指的是排列成

$$
\begin{aligned}
&a_{11}+a_{12}+a_{13}+a_{14}+\cdots+a_{1n}+\cdots\\
+{}&a_{21}+a_{22}+a_{23}+a_{24}+\cdots+a_{2n}+\cdots\\
+{}&\cdots\\
+{}&a_{m1}+a_{m2}+a_{m3}+a_{m4}+\cdots+a_{mn}+\cdots\\
+{}&\cdots
\end{aligned}\tag{4.10}
$$

的方阵，这个方阵的右端和下端都是无限的。方阵的每一项用 $a_{kl}$ 表示，其中的第一个指标 $k$ 表示行，第二个指标 $l$ 表示列。

现在求出方阵的前 $m$ 行 $n$ 列共 $m\times n$ 项之和

$$S_{mn}=\sum_{\substack{1\le k\le m\\1\le l\le n}}a_{kl},\tag{4.11}$$

根据 $S_{mn}$ 就能构造出这个二重级数的部分和序列。容易理解，如果部分和序列收敛，

$$\lim_{\substack{m\to\infty\\n\to\infty}}S_{mn}=S,\tag{4.12}$$

则称此二重级数收敛。$S$ 就是这个二重级数之和：

$$S=\sum_{k,l=1}^{\infty}a_{kl}.\tag{4.13}$$

这时，$\forall\varepsilon>0$，总可以找到正整数 $N$，当 $m,n>N$ 时，恒有

$$|S_{mn}-S|<\varepsilon.\tag{4.14}$$

**例 4.1** 二重级数

$$
\begin{aligned}
&1+1+1+1+\cdots\\
+{}&1-1-1-1-\cdots\\
+{}&1-1+0+0+\cdots\\
+{}&1-1+0+0+\cdots\\
+{}&\cdots
\end{aligned}\tag{4.15}
$$

的 $S_{mn}=2,\quad m,n>1$，所以，这个二重级数之和为 $S=2$。

除了这种求和方式之外，当然还可以考虑其他求和方式。例如，考虑到上面的部分和

$$S_{mn}=\sum_{\substack{1\le k\le m\\1\le l\le n}}a_{kl}=\sum_{k=1}^{m}\left(\sum_{l=1}^{n}a_{kl}\right)=\sum_{l=1}^{n}\left(\sum_{k=1}^{m}a_{kl}\right),$$

还可以有累次求和：

$$\sum_{k=1}^{\infty}\left(\sum_{l=1}^{\infty}a_{kl}\right)=\lim_{m\to\infty}\left(\lim_{n\to\infty}S_{mn}\right),$$

或

$$\sum_{l=1}^{\infty}\left(\sum_{k=1}^{\infty}a_{kl}\right)=\lim_{n\to\infty}\left(\lim_{m\to\infty}S_{mn}\right).$$

前者是先按列求和，再将各列之和相加（可称为逐列求和）；后者的求和次序则相反（可称为逐行求和）。需要注意，即使二重级数收敛，某些行或列的和也不一定存在，因此累次求和的和也不一定存在。例如，上面的二重级数 (4.15) 式就是如此，它的第一列和第二列的列级数都不收敛，第一行和第二行的行级数也不收敛。

而且，如果逐行和逐列求和的和都存在，这两个和数也不一定相等（即和数与求和次序有关）。例如，二重级数

$$a_{kl}=\frac1{k+1}\left(\frac{k}{k+1}\right)^l-\frac1{k+2}\left(\frac{k+1}{k+2}\right)^l$$

的部分和是

$$S_{mn}=\frac12\left[1-\left(\frac12\right)^n\right]-\frac{m+1}{m+2}\left[1-\left(\frac{m+1}{m+2}\right)^n\right],$$

所以

$$\lim_{m\to\infty}\left(\lim_{n\to\infty}S_{mn}\right)=-\frac12,\qquad
\lim_{n\to\infty}\left(\lim_{m\to\infty}S_{mn}\right)=\frac12.$$

即使逐列求和与逐行求和的和数相等，二重级数也不一定收敛。例如，Kelvin 在讨论两个带电球之间的相互作用力时，就曾经得到通项为 $a_{kl}=(-1)^{k+l}kl/(k+l)^2$ 的二重级数

$$\sum_{k=1}^{\infty}\left(\sum_{l=1}^{\infty}a_{kl}\right)=\sum_{l=1}^{\infty}\left(\sum_{k=1}^{\infty}a_{kl}\right)=\frac16\left(\ln2-\frac14\right),$$

但部分和序列是在

$$\frac16\left(\ln2-\frac58\right)\quad\text{和}\quad\frac16\left(\ln2+\frac18\right)$$

之间振荡，所以 $\sum_{k,l=1}^{\infty}a_{kl}$ 并不收敛。

前面讨论级数乘法时，还涉及另一种特殊的求和方式，见 (4.9) 式。从方阵 (4.10) 来看，这相当于按（次）对角线求和。这个和数，也不一定等于逐列或逐行求和的和数。例如，将二重级数 (4.15) 按对角线求和，得到的和数为 4。

二重级数的和是否依赖于求和方式，原则上与级数是否绝对收敛有关。如果二重级数绝对收敛，则级数各项的先后次序可以重新排列，因而不同求和方式得到相同的和数。

## §4.3 函数级数

设 $u_k(z)(k=1,2,\cdots)$ 在区域 $G$ 内有定义。若对于 $G$ 内一点 $z_0$，级数 $\sum_{k=1}^{\infty}u_k(z_0)$ 收敛，则称级数 $\sum_{k=1}^{\infty}u_k(z)$ 在 $z_0$ 点收敛。反之，若 $\sum_{k=1}^{\infty}u_k(z_0)$ 发散，则称级数 $\sum_{k=1}^{\infty}u_k(z)$ 在 $z_0$ 点发散。

如果级数 $\sum_{k=1}^{\infty}u_k(z)$ 在区域 $G$ 内每一点都收敛，则称级数在 $G$ 内（逐点）收敛。其和函数 $S(z)$ 是 $G$ 内的单值函数。

**例 4.2** 考察几何级数 $\sum_{n=0}^{\infty}z^n=1+z+z^2+\cdots+z^n+\cdots$ 的收敛性。

**解** 几何级数的部分和为

$$S_n=1+z+z^2+\cdots+z^{n-1}=\frac{1-z^n}{1-z}.$$

由于当 $|z|<1$ 时，$z^n\to0$，此时部分和序列 $\{S_n\}$ 收敛，

$$\lim_{n\to\infty}S_n=\frac1{1-z}.$$

而当 $|z|\ge1$ 时，$|z^n|\ge1$，不满足级数收敛的必要条件 (4.4) 式，因此几何级数发散。

综上所述，知

$$\sum_{n=0}^{\infty}z^n=1+z+z^2+\cdots+z^n+\cdots=
\begin{cases}
\dfrac1{1-z},&|z|<1,\\
\text{发散},&|z|\ge1.
\end{cases}\tag{4.16}$$

需要强调的是，所有关于收敛函数级数的等式，一定都是有条件的，即只在函数级数收敛的范围内成立。因此关于函数级数的等式，一定要注明等式成立的条件。

若 $\forall\varepsilon>0$，$\exists$ 与 $z$ 无关的 $N(\varepsilon)$，使 $\forall n>N(\varepsilon)$，$\forall z\in G$，都有

$$\left|S(z)-\sum_{k=1}^{n}u_k(z)\right|<\varepsilon,$$

则称级数 $\sum_{k=1}^{\infty}u_k(z)$ 在 $G$ 内一致收敛。

显然，一致收敛的概念总是和一定的区域联系在一起的。级数的一致收敛性是它在一定区域内的性质。

判断级数是否一致收敛，除直接运用定义外，常用 Weierstrass 的 $M$ 判别法：若 $\exists N>0$，$\forall k>N$，$\forall z\in G$，$|u_k(z)|<a_k$，且 $a_k$ 与 $z$ 无关，而 $\sum_{k=1}^{\infty}a_k$ 收敛，则级数 $\sum_{k=1}^{\infty}u_k(z)$ 在 $G$ 内绝对而且一致收敛。

**练习 4.2** 设 $x$ 为实数，证明级数

$$\sum_{n=1}^{\infty}\frac{x^2}{(1+x^2)^n}$$

绝对收敛，但不一致收敛；而级数

$$\sum_{n=1}^{\infty}\frac{(-)^n}{n+x^2}$$

一致收敛，但不绝对收敛。

一致收敛的级数具有下列重要性质。

**1. 连续性** 如果 $u_k(z)$ 在 $G$ 内连续，级数 $\sum_{k=1}^{\infty}u_k(z)$ 在 $G$ 内一致收敛，则其和函数

$$S(z)=\sum_{k=1}^{\infty}u_k(z)$$

也在 $G$ 内连续。

这个性质说明，如果级数的每一项都是连续函数，则一致收敛级数可以逐项求极限：

$$\lim_{z\to z_0}\sum_{k=1}^{\infty}u_k(z)=\sum_{k=1}^{\infty}\left[\lim_{z\to z_0}u_k(z)\right].\tag{4.17}$$

**2. 逐项求积分** 设 $C$ 是区域 $G$ 内的一条分段光滑曲线，如果 $u_k(z)(k=1,2,\cdots)$ 是 $C$ 上的连续函数，则对于 $C$ 上一致收敛的级数 $\sum_{k=1}^{\infty}u_k(z)$ 可以逐项求积分：

$$\int_C\sum_{k=1}^{\infty}u_k(z)\,dz=\sum_{k=1}^{\infty}\int_Cu_k(z)\,dz.\tag{4.18}$$

**3. 逐项求导数（Weierstrass 定理）** 设 $u_k(z)(k=1,2,\cdots)$ 在 $\overline G$ 中单值解析，$\sum_{k=1}^{\infty}u_k(z)$ 在 $\overline G$ 中一致收敛，则此级数之和 $f(z)$ 是 $G$ 内的解析函数，$f(z)$ 的各阶导数可以由 $\sum_{k=1}^{\infty}u_k(z)$ 逐项求导数得到，

$$f^{(p)}(z)=\sum_{k=1}^{\infty}u_k^{(p)}(z),\tag{4.19}$$

求导数后的级数在 $G$ 内内闭一致收敛①。

① 级数在 $G$ 内内闭一致收敛，意即 $\forall\overline{G'}\subset G$，该级数在闭区域 $\overline{G'}$ 中都一致收敛。

## §4.4 幂级数

幂级数是通项为幂函数的函数项级数，

$$\sum_{n=0}^{\infty}c_n(z-a)^n=c_0+c_1(z-a)+c_2(z-a)^2+\cdots+c_n(z-a)^n+\cdots.\tag{4.20}$$

这是一种特殊形式的函数项级数，也是最基本、最常用的一种函数项级数。

**定理 4.1（Abel 第一定理）** 如果级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在某点 $z_0$ 收敛，则在以 $a$ 点为圆心、$|z_0-a|$ 为半径的圆内绝对收敛，而在圆 $|z-a|\le r\ (r<|z_0-a|)$ 中一致收敛。

**证** 因为 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在 $z_0$ 收敛，故一定满足必要条件

$$\lim_{n\to\infty}c_n(z_0-a)^n=0.$$

因此对 $\forall q>0$，存在正整数 $N$，使 $\forall n>N$，$|c_n(z_0-a)^n|<q$。所以，

$$|c_n(z-a)^n|=|c_n(z_0-a)^n|\left|\frac{z-a}{z_0-a}\right|^n<q\left|\frac{z-a}{z_0-a}\right|^n,\qquad n>N.$$

因当 $\left|\dfrac{z-a}{z_0-a}\right|<1$（即 $|z-a|<|z_0-a|$）时，$\sum_{n=0}^{\infty}\left|\dfrac{z-a}{z_0-a}\right|^n$ 收敛，故

$$\sum_{n=0}^{\infty}c_n(z-a)^n\ \text{在圆 }|z-a|<|z_0-a|\text{ 内绝对收敛}.\tag{4.21}$$

而当 $|z-a|\le r<|z_0-a|$ 时，

$$|c_n(z-a)^n|\le q\frac{r^n}{|z_0-a|^n},\qquad n>N,$$

由于常数项级数 $\sum_{n=0}^{\infty}\dfrac{r^n}{|z_0-a|^n}$ 收敛，故

$$\sum_{n=0}^{\infty}c_n(z-a)^n\ \text{在圆 }|z-a|\le r\ (<|z_0-a|)\text{ 中一致收敛}.\tag{4.22}$$

Abel 第一定理简称为 Abel 定理。

**推论** 若级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在某点 $z_1$ 发散，则在圆 $|z-a|>|z_1-a|$ 外处处发散。

**证** 用反证法。假设级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在圆 $|z-a|>|z_1-a|$ 外某一点 $z_2$ 收敛，则按 Abel 定理，级数必然在圆 $|z-a|<|z_2-a|$ 内收敛，由于 $|z_1-a|<|z_2-a|$，因此级数在 $z_1$ 点收敛，与原设矛盾。故级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在圆 $|z-a|>|z_1-a|$ 外处处发散。□

Abel 定理及其推论意味着，在幂级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 的收敛点与发散点之间存在一个分界线，而且这个分界线一定是圆周。圆内区域称为幂级数的收敛圆。收敛圆的圆心一定是 $z=a$ 点，收敛圆的半径称为收敛半径。例如 (4.16) 式的几何级数就是一个幂级数，其收敛圆圆心是 $z=0$ 点，收敛圆的半径是 1，几何级数在单位圆 $|z|<1$ 内绝对收敛，在单位圆内的任意一个闭区域中一致收敛，当然在任意一个比单位圆小的圆 $|z|\le r\ (r<1)$ 中一致收敛。

作为特殊情况，收敛半径可以是 0，即收敛圆退化为一个点。除 $z=a$ 点外，幂级数在全复平面处处发散；也可以是 $\infty$，收敛圆就是全复平面。幂级数在全复平面收敛，但在 $\infty$ 点肯定发散，除非此幂级数只有常数项一项。

求幂级数的收敛半径的办法，常用的有两个：

(1) Cauchy-Hadamard 公式。根据 Cauchy 判别法，

当

$$\overline{\lim}_{n\to\infty}|c_n(z-a)^n|^{1/n}<1,$$

即

$$|z-a|<\frac1{\overline{\lim}_{n\to\infty}|c_n|^{1/n}},$$

时，级数绝对收敛；

当

$$\overline{\lim}_{n\to\infty}|c_n(z-a)^n|^{1/n}>1,$$

即

$$|z-a|>\frac1{\overline{\lim}_{n\to\infty}|c_n|^{1/n}},$$

时，级数发散。

因此，幂级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 的收敛半径是

$$R=\frac1{\overline{\lim}_{n\to\infty}|c_n|^{1/n}}=\underline{\lim}_{n\to\infty}\left|\frac1{c_n}\right|^{1/n}.\tag{4.23}$$

(2) 根据 d'Alembert 判别法，如果

$$\lim_{n\to\infty}\left|\frac{c_{n+1}(z-a)^{n+1}}{c_n(z-a)^n}\right|=|z-a|\lim_{n\to\infty}\left|\frac{c_{n+1}}{c_n}\right|$$

存在，则

当

$$\lim_{n\to\infty}\left|\frac{c_{n+1}(z-a)^{n+1}}{c_n(z-a)^n}\right|<1,$$

即 $|z-a|<\lim\limits_{n\to\infty}\left|\dfrac{c_n}{c_{n+1}}\right|$ 时，级数绝对收敛；

当

$$\lim_{n\to\infty}\left|\frac{c_{n+1}(z-a)^{n+1}}{c_n(z-a)^n}\right|>1,$$

即 $|z-a|>\lim\limits_{n\to\infty}\left|\dfrac{c_n}{c_{n+1}}\right|$ 时，级数发散。

因此，幂级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 的收敛半径是

$$R=\lim_{n\to\infty}\left|\frac{c_n}{c_{n+1}}\right|.\tag{4.24}$$

这两个求收敛半径的公式各有优缺点。Cauchy-Hadamard 公式是普遍成立的，d'Alembert 公式则是有条件的（要求极限 $\lim\limits_{n\to\infty}|c_n/c_{n+1}|$ 存在）。但当后者能适用时，往往计算更简单。

**练习 4.3** 已知 $\sum_{n=1}^{\infty}a_nz^n$ 和 $\sum_{n=1}^{\infty}b_nz^n$ 的收敛半径分别为 $R_1$ 和 $R_2$，求下列幂级数的收敛半径：

(1) $\displaystyle\sum_{n=1}^{\infty}(a_n-b_n)z^n$；

(2) $\displaystyle\sum_{n=1}^{\infty}a_nb_nz^n$；

(3) $\displaystyle\sum_{n=1}^{\infty}\frac1{a_n}z^n,\quad a_n\ne0$；

(4) $\displaystyle\sum_{n=1}^{\infty}\frac{b_n}{a_n}z^n,\quad a_n\ne0$。

由于幂级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 的每一项在复平面 $\mathbb C$ 内都是 $z$ 的解析函数，Abel 定理又告诉我们，幂级数在其收敛圆内内闭一致收敛，因此，根据 §4.3 的 Weierstrass 定理，幂级数在收敛圆内代表了一个解析函数，可以对幂级数逐项积分或逐项求导数，而收敛半径不变：

$$\int_{z_0}^{z}\sum_{n=0}^{\infty}c_n(z-a)^n\,dz
=\sum_{n=0}^{\infty}c_n\int_{z_0}^{z}(z-a)^n\,dz
=\sum_{n=0}^{\infty}\frac{c_n}{n+1}\left[(z-a)^{n+1}-(z_0-a)^{n+1}\right],\tag{4.25}$$

$$\frac{d}{dz}\left[\sum_{n=0}^{\infty}c_n(z-a)^n\right]
=\sum_{n=0}^{\infty}c_n\frac{d(z-a)^n}{dz}
=\sum_{n=0}^{\infty}c_{n+1}(n+1)(z-a)^n.\tag{4.26}$$

幂级数在收敛圆内一定收敛，在收敛圆外一定发散。在收敛圆的圆周上，级数可能在所有点都收敛，可能在所有点都发散，也可能在一部分点收敛，在另一部分点发散。例如：

$$1+z+\cdots+z^n+\cdots\qquad\text{在 }|z|=1\text{ 上处处发散};$$

$$\frac z1+\frac{z^2}{2}+\cdots+\frac{z^n}{n}+\cdots\qquad\text{在 }|z|=1\text{ 上除 }z=1\text{ 外均收敛，而在 }z=1\text{ 点发散};$$

$$\frac{z^2}{1\cdot2}+\frac{z^3}{2\cdot3}+\cdots+\frac{z^n}{(n-1)n}+\cdots\qquad\text{在 }|z|=1\text{ 上处处收敛}.$$

但不论哪种情况，幂级数的收敛圆的圆周上总肯定有和函数的奇点。特别需要说明，即使在和函数的奇点处，幂级数仍然可能收敛。读者可以求出上面三个级数的和函数，验证它们在 $z=1$ 点不解析。

当幂级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在收敛圆周上某点 $z_0$ 收敛时，其和与级数在收敛圆内的和函数之间的关系，有下面的定理。

![Figure 4.1](...)

图 4.1 Abel 第二定理

**定理 4.2（Abel 第二定理）** 若幂级数 $\sum_{n=0}^{\infty}c_n(z-a)^n$ 在收敛圆内收敛到 $f(z)$，且在收敛圆周上某点 $z_0$ 也收敛，和为 $S(z_0)$，则当 $z$ 由收敛圆内趋于 $z_0$ 时，只要保持在以 $z_0$ 为顶点、张角为 $2\phi<\pi$ 的范围内（见图 4.1），$f(z)$ 就一定趋于 $S(z_0)$。

需要明确一下（当收敛半径为有限值时）“和函数的奇点（或解析性）”这种说法的含义。毫无疑问，幂级数在收敛圆内收敛，在收敛圆外发散，因而幂级数只在收敛圆内有定义，并且代表了一个解析函数。另一方面，就和函数而言，尽管幂级数是在收敛圆内才收敛到和函数，但是，我们从来不会认为这个函数只是局限于收敛圆内才有定义（少数函数除外，例如，见第五章习题 10）。事实上，幂级数只不过是这个函数在收敛圆内的一种表达形式，在其他区域内可以有其他表达形式。我们总会自觉或不自觉地在整个定义域上来考察它的解析性。我们也正是在这个前提下，才能谈论函数在收敛圆周上乃至收敛圆外的奇点。

## §4.5 含参量的反常积分的解析性

§4.3 中有关函数级数解析性的结论，也可以用来讨论含参量的反常积分的解析性。

**定理 4.3** 设

(1) $f(t,z)$ 是 $t$ 和 $z$ 的连续函数，$t>a,\ z\in\overline G$，

(2) 对于任何 $t\ge a$，$f(t,z)$ 是 $\overline G$ 中的单值解析函数，

(3) 积分 $\displaystyle\int_a^{\infty}f(t,z)\,dt$ 在 $\overline G$ 中一致收敛，即 $\forall\varepsilon>0$，$\exists T(\varepsilon)$，当 $T_2>T_1>T(\varepsilon)$ 时，对 $\forall z\in\overline G$ 都有

$$\left|\int_{T_1}^{T_2}f(t,z)\,dt\right|<\varepsilon,$$

则 $F(z)=\displaystyle\int_a^{\infty}f(t,z)\,dt$ 在 $G$ 内是解析的，且

$$F'(z)=\int_a^{\infty}\frac{\partial f(t,z)}{\partial z}\,dt.\tag{4.27}$$

**证** 任取一个单调无界序列 $\{a_n\}$，

$$a_0=a<a_1<a_2<a_3<\cdots<a_n<a_{n+1}<\cdots,\qquad \lim_{n\to\infty}a_n=\infty.$$

令

$$u_n(z)=\int_{a_n}^{a_{n+1}}f(t,z)\,dt,$$

则根据定理 3.3（§3.6，含参量定积分的解析性）可知，$u_n(z)$ 是 $G$ 内的单值解析函数。又因为

$$F(z)=\sum_{n=0}^{\infty}u_n(z)=\int_a^{\infty}f(t,z)\,dt$$

在 $\overline G$ 中一致收敛，故根据 Weierstrass 定理，知 $F(z)$ 在 $G$ 内解析，且

$$F'(z)=\sum_{n=0}^{\infty}u_n'(z)=\int_a^{\infty}\frac{\partial f(t,z)}{\partial z}\,dt.$$

□

对于含参量的瑕积分也可以类似地处理。

在应用这个定理时，需要判断无穷积分（或瑕积分）是否一致收敛。常用的判别法是：**如果存在函数 $\phi(t)$，使得 $\exists T>a$，对 $\forall t>T$，$\forall z\in\overline G$，都有 $|f(t,z)|<\phi(t)$，而且 $\int_a^{\infty}\phi(t)dt$ 收敛，则 $\int_a^{\infty}f(t,z)dt$ 在 $\overline G$ 中绝对而且一致收敛。**

作为含参量的无穷积分的一个例子，下面讨论积分

$$F(z)=\int_0^{\infty}e^{-t^2}\cos2zt\,dt.\tag{4.28}$$

这个积分中的被积函数显然满足定理的前两个条件，而且因为对于复数 $z=x+iy$，有

$$|\cos2zt|=\sqrt{\cosh^22yt-\cos^22xt}\le\cosh2|yt|\le e^{2|yt|}.$$

所以，对于 $z$ 复平面内任意一个有界闭区域中的点，都有 $|\operatorname{Im}z|<y_0$，于是

$$|e^{-t^2}\cos2zt|<e^{-t^2+2y_0t},\qquad t>0.$$

而积分 $\int_0^{\infty}e^{-t^2+2y_0t}dt$ 收敛，所以含参量的无穷积分 (4.28) 一致收敛，因此，这个积分作为 $z$ 的函数，在 $z$ 复平面内的任意一个有界区域内解析。更进一步，就有

$$
\begin{aligned}
F'(z)&=-\int_0^{\infty}e^{-t^2}2t\sin2zt\,dt\\
&=e^{-t^2}\sin2zt\bigg|_0^{\infty}-2z\int_0^{\infty}e^{-t^2}\cos2zt\,dt=-2zF(z).
\end{aligned}
$$

解这个微分方程，就可以得到 $F(z)=Ce^{-z^2}$，其中常数 $C$ 是

$$C=F(0)=\int_0^{\infty}e^{-t^2}dt=\frac12\sqrt\pi,$$

这样，最后就得到

$$\int_0^{\infty}e^{-t^2}\cos2zt\,dt=\frac12\sqrt\pi\,e^{-z^2}.\tag{4.29}$$

## \*§4.6 发散级数与渐近级数

迄今为止，我们的注意力都集中在级数的收敛性上。讨论如何判断级数收敛，收敛级数具有哪些性质，甚至还具体讨论如何求级数的和。但是，绝不要认为发散级数毫无用处。在历史上，许多著名的数学家都研究或使用过发散级数，例如，Stirling 在他的《微分法》（Methodus Differentialis）一书中就曾经给出了

$$\ln m!=\left(m+\frac12\right)\ln\left(m+\frac12\right)-\left(m+\frac12\right)+\frac12\ln(2\pi)+\sum_{n=1}^{\infty}\frac{B_{2n}(1/2)}{(2n-1)(2n)}\left(m+\frac12\right)^{-(2n-1)},\tag{4.30}$$

其中的 $B_n(t)$ 是 Bernoulli 多项式

$$\frac{ze^{tz}}{e^z-1}=\sum_{n=0}^{\infty}\frac{B_n(t)}{n!}z^n.\tag{4.31}$$

Stirling 用级数 (4.30) 的前几项就能计算 $\lg1000!$ 到 10 位小数。但这个级数是发散的，不过它的头几项非常迅速地减小，不多几项就可以给出足够好的近似。Euler 也曾讨论过下面几个级数：

$$1-1+1-1+1-1+1-1+\cdots;\tag{4.32}$$

$$1-2+3-4+5-6+\cdots;\tag{4.33}$$

$$1-2+2^2-2^3+2^4-2^5+2^6+\cdots;\tag{4.34}$$

$$1-2!+3!-4!+5!-6!+\cdots.\tag{4.35}$$

求出的“和”分别是

$$\frac12,\qquad\frac14,\qquad\frac13,\qquad0.4036\cdots.$$

至于物理学家更不排斥发散级数。事实上，物理学中广泛使用的渐近展开常常就是发散级数。问题是，从数学上说，必须在收敛级数的和的概念的基础上，建立起发散级数的“和”的概念。

不妨再仔细分析一下级数 (4.32)，它的部分和序列

$$S_N=\sum_{n=0}^{N}(-1)^n=\frac12[1-(-1)^N]$$

并不收敛，不能用部分和序列的极限来定义它的和。但是可以造一个幂级数

$$\sum_{n=0}^{\infty}(-)^nz^n,\tag{4.36}$$

收敛区域是 $|z|<1$。上面的级数 (4.32) 就可以看成是将幂级数 (4.36) 不合法地取极限 $z\to1$ 的结果。事实上，这个发散级数在 Fourier 的 *Theorie Analytique de la Chaleur*（中译本：《热的解析理论》，桂质亮译，武汉出版社，1993）一书中就曾出现过。他在计算 Fourier 展开时，得到（见该书中译本第 170—171 页）

$$
\begin{aligned}
\frac{\pi}{2}\frac{\sinh x}{\sinh\pi}
={}&\left(\sin x-\frac12\sin2x+\frac13\sin3x-\cdots\right)
-\left(\sin x-\frac1{2^3}\sin2x+\frac1{3^3}\sin3x-\cdots\right)\\
&+\left(\sin x-\frac1{2^5}\sin2x+\frac1{3^5}\sin3x-\cdots\right)
-\left(\sin x-\frac1{2^7}\sin2x+\frac1{3^7}\sin3x-\cdots\right)+\cdots\\
={}&\frac{\sin x}{1+\frac11}-\frac{\sin2x}{2+\frac12}+\frac{\sin3x}{3+\frac13}+\cdots.
\end{aligned}
$$

在导出这个结果时，Fourier 用到了和式

$$\frac1n-\frac1{n^3}+\frac1{n^5}-\cdots=\frac{n}{1+n^2},$$

但是这个公式在 $n=1$ 时并不成立，特别是，左端的级数发散。因此，上述展开中 $\sin x$ 项的展开系数 $a_1$ 不可能表示成这样的发散级数，这个发散级数更不可能收敛到后面的和。实际上，正确的结果是

$$a_1=\frac1{\sinh\pi}\int_0^{\pi}\sinh x\sin x\,dx=\frac12,$$

恰恰又和前面的不合法的计算得到的结果一致。这里要注意，级数

$$1+0-1+1+0-1+1+0-1+\cdots\tag{4.37}$$

不同于级数 (4.32)。级数 (4.37) 的“和”是 $2/3$。

可以类似地讨论另外几个发散级数，只不过由级数 (4.34) 和 (4.35) 造出的幂级数

$$\sum_{n=0}^{\infty}(-)^n2^nz^n\qquad\text{和}\qquad\sum_{n=0}^{\infty}(-)^n(n+1)!z^n$$

的收敛半径分别是 $1/2$ 和 0，而级数 (4.34) 和 (4.35) 的“和”则是这两个幂级数的和函数在收敛圆外的 $z=1$ 点的值。

以上讨论的是几个典型的发散级数。在一种特殊的限制条件下，这些发散级数也可以用于合法的计算。上面的求“和”规则只不过是在这种背景下更复杂的运算过程的缩写。

非收敛级数主要出现在渐近级数中。

在介绍渐近级数（或称渐近展开）的概念之前，还要先引进记号 $O$ 和 $o$。设 $f(z)$ 和 $\phi(z)$ 在 $z_0$ 点的邻域内有定义，且 $\phi(z)\ne0$，若 $z\to z_0$ 时，$f(z)/\phi(z)$ 有界，则记为

$$f(z)=O(\phi(z)),\qquad\text{当 }z\to z_0;$$

若 $\lim\limits_{z\to z_0}f(z)/\phi(z)=0$，则记为

$$f(z)=o(\phi(z)),\qquad\text{当 }z\to z_0.$$

**渐近级数** 若当 $z\to z_0$ 时，对于每一个 $m$ 值，都有

$$f(z)-\sum_{n=0}^{m}a_n\phi_n(z)=o(\phi_m(z)),\tag{4.38}$$

则称 $\sum_{n=0}^{\infty}a_n\phi_n(z)$ 是函数 $f(z)$ 相对于 $\{\phi_n(z)\}$ 的渐近级数，记为

$$f(z)\sim\sum_{n=0}^{\infty}a_n\phi_n(z).\tag{4.39}$$

渐近级数的定义说明，$z$ 越接近 $z_0$，有限和 $\sum_{n=0}^{m}a_n\phi_n(z)$（称为 $z\to z_0$ 时 $f(z)$ 的渐近近似）越逼近于 $f(z)$。它区别于通常的级数展开，例如 Taylor 展开（见第 5 章）

$$f(z)=\sum_{n=0}^{N}a_nz^n+\cdots,$$

后者是 $z$ 点固定，级数的项数越多越准确，

$$\lim_{N\to\infty}\left[f(z)-\sum_{n=0}^{N}u_n(z)\right]=0.$$

特别是，在渐近级数的定义中，并未要求级数 $\sum_{n=0}^{\infty}a_n\phi_n(z)$ 收敛。渐近展开级数可以（而且常常）不是收敛级数。因此，对于一定的 $z$，并不能通过多取项数（即增大 $N$）来无限制地改善近似程度。例如，对于指数积分

$$\operatorname{Ei}(-x)=-\int_x^{\infty}\frac{e^{-t}}t\,dt,\qquad x>0,\tag{4.40}$$

用分部积分的方法可以得到

$$
\begin{aligned}
-\operatorname{Ei}(-x)
&=\frac{e^{-x}}x-\int_x^{\infty}\frac{e^{-t}}{t^2}\,dt=\cdots\\
&=\frac{e^{-x}}x\left[1-\frac1x+\frac{2!}{x^2}-\frac{3!}{x^3}+\cdots+\frac{(-)^{n-1}n!}{x^n}\right]
+(-)^{n+1}(n+1)!\int_x^{\infty}\frac{e^{-t}}{t^{n+2}}\,dt.
\end{aligned}\tag{4.41}
$$

容易证明，它的余项

$$\int_x^{\infty}\frac{e^{-t}}{t^{n+2}}\,dt<\frac1{x^{n+2}}\int_x^{\infty}e^{-t}\,dt=\frac{e^{-x}}{x^{n+2}},$$

因此，(4.41) 式的确给出了 $-\operatorname{Ei}(-x)$ 在 $x\to\infty$ 时的渐近展开。然而，级数相邻两项之比

$$\lim_{n\to\infty}\left|\frac{u_{n+1}}{u_n}\right|=\lim_{n\to\infty}\frac nx=\infty,$$

由 Cauchy 判别法可知级数发散，且对于给定的 $x$，$-\operatorname{Ei}(-x)$ 的渐近近似中取 $N\approx x$ 项可以得到最佳逼近：在此项之前，级数各项绝对值递降，而此后各项的绝对值反而递增。

渐近展开不同于 Taylor 展开或 Laurent 展开（见第五章），还在于渐近级数通常都有一定的辐角限制，即渐近展开只在一定的辐角范围内成立①。同一个函数在不同的辐角范围内，渐近展开的形式可以不同；即使在两个不同区域的公共区域内，两个渐近展开也可以有明显不同的结果。

① 这时关于 $O$ 和 $o$ 的定义以及渐近序列的概念都应作相应的修改。

在 $\arg z$ 的一定范围内，渐近展开（如果存在）是唯一的，系数由

$$a_m=\lim_{z\to z_0}\frac1{\phi_m(z)}\left[f(z)-\sum_{n=0}^{m-1}a_n\phi_n(z)\right]\tag{4.42}$$

决定。但不同的函数在同一个区域内可以有相同形式的渐近展开。例如，当 $|\arg z|<\alpha<\pi/2$，而 $|z|\to\infty$ 时，如果有

$$J(z)\sim\sum_{n=0}^{N}a_n\phi_n(z)+\cdots,$$

则对于函数 $J(z)+e^{-z}$，同样有

$$J(z)+e^{-z}\sim\sum_{n=0}^{N}a_n\phi_n(z)+\cdots.$$

如果要讨论扇形区域 $\alpha<\arg z<\beta$ 中，$z\to\infty$ 时的渐近展开，最简单的渐近序列是 $\{\phi(z)z^{-n}\}$，其中 $\phi(z)$ 在 $z\to\infty$ 时的行为已知（它就决定了扇形区域的辐角范围），

$$f(z)\sim\phi(z)\sum_{n=0}^{\infty}a_nz^{-n},\tag{4.43}$$

而将

$$\frac{f(z)}{\phi(z)}\sim\sum_{n=0}^{\infty}a_nz^{-n}\tag{4.44}$$

称为渐近幂级数。

渐近幂级数和收敛幂级数具有相似的运算性质。例如，当 $z\to\infty$ 时，

$$f(z)\sim\sum_{n=0}^{\infty}a_nz^{-n},\qquad g(z)\sim\sum_{n=0}^{\infty}b_nz^{-n},\qquad \alpha<\arg z<\beta,$$

则根据渐近级数的定义，可以证明：

(1) $Af(z)\sim A\displaystyle\sum_{n=0}^{\infty}a_nz^{-n}$，$A$ 为常数；

(2) $f(z)+g(z)\sim\displaystyle\sum_{n=0}^{\infty}(a_n+b_n)z^{-n}$；

(3) $f(z)\cdot g(z)\sim\displaystyle\sum_{n=0}^{\infty}\left(\sum_{k=0}^{n}a_kb_{n-k}\right)z^{-n}$；

(4)

$$\frac1{f(z)}\sim\frac1{a_0}+\sum_{n=1}^{\infty}d_nz^{-n},$$

其中

$$d_1=\lim_{z\to\infty}z\left[\frac1{f(z)}-\frac1{a_0}\right]=-\frac{a_1}{a_0^2},$$

$$d_2=\lim_{z\to\infty}z^2\left[\frac1{f(z)}-\frac1{a_0}-\frac{d_1}{z}\right]=\frac{a_1^2-a_0a_2}{a_0^3},\ \cdots.$$

此外，若 $|z|>R$ 时 $f(z)$ 连续，则当 $|z|>R$ 时，

$$F(z)=\int_z^{\infty}[f(t)-a_0-a_1t^{-1}]\,dt$$

连续；若在区域 $G:|z|>R,\alpha<\arg z<\beta$ 内 $f(z)$ 解析，且在 $G$ 所含的任一闭扇形中，当 $z\to\infty$ 时，对 $\arg z$ 一致地有

$$f(z)\sim a_0+a_1z^{-1}+a_2z^{-2}+a_3z^{-3}+\cdots,$$

则在 $G$ 所含的任一闭扇形中，当 $z\to\infty$ 时，对 $\arg z$ 一致地有

$$f'(z)\sim-a_1z^{-2}-2a_2z^{-3}-3a_3z^{-4}+\cdots.$$

有关求渐近展开的具体方法，请参阅参考书目 [13] 和 [33]，[34]，[35]。

## 习题

1. 判断下列级数的收敛性与绝对收敛性：

(1) $\displaystyle\sum_{n=2}^{\infty}\frac{i^n}{\ln n}$；

(2) $\displaystyle\sum_{n=1}^{\infty}\frac{i^n}{n}$。

2. 证明级数

$$\sum_{n=1}^{\infty}\frac{z^{n-1}}{(1-z^n)(1-z^{n+1})},\qquad |z|\ne1$$

收敛，并求其和。

3. 求下列二重级数之和：

(1) $\displaystyle\sum_{n=2}^{\infty}\sum_{m=2}^{\infty}\frac1{(2m)^n}$；

(2) $\displaystyle\sum_{n=1}^{\infty}\sum_{m=1}^{\infty}\frac1{(4m-1)^{2n+1}}$；

(3) $\displaystyle\sum_{n=1}^{\infty}\sum_{m=1}^{\infty}\frac1{(4m-1)^{2n}}$；

(4) $\displaystyle\sum_{n=1}^{\infty}\sum_{m=1}^{\infty}\frac1{(4m-2)^{2n}}$。

提示：上述级数均绝对收敛，因此它们的和与求和次序无关。不妨先对 $n$ 求和。

4. 试确定下列级数的收敛区域：

(1) $\displaystyle\sum_{n=1}^{\infty}z^{n!}$；

(2) $\displaystyle\sum_{n=1}^{\infty}\left(\frac{z}{1+z}\right)^n$；

(3) $\displaystyle\sum_{n=1}^{\infty}(-)^n(z^2+2z+2)^n$；

(4) $\displaystyle\sum_{n=1}^{\infty}2^n\sin\frac{z}{3^n}$。

5. 证明级数

$$\sum_{n=0}^{\infty}\left(\frac{z^{n+1}}{n+1}-\frac{2z^{2n+3}}{2n+3}\right)$$

的和函数在 $z=1$ 点不连续。

6. 证明：

$$\ln(1-z)=-z-\frac{z^2}{2}-\frac{z^3}{3}-\frac{z^4}{4}-\cdots,\qquad |z|<1,$$

并由此导出

$$r\cos\theta-r^2\frac{\cos2\theta}{2}+r^3\frac{\cos3\theta}{3}-\cdots=\frac12\ln(1+2r\cos\theta+r^2),$$

$$r\sin\theta-r^2\frac{\sin2\theta}{2}+r^3\frac{\sin3\theta}{3}-\cdots=\arctan\frac{r\sin\theta}{1+r\cos\theta},$$

其中 $-1<r<1$。

7. 求下列级数之和：

(1)

$$\cos\theta+\frac{\cos2\theta}{2}+\frac{\cos3\theta}{3}+\frac{\cos4\theta}{4}+\cdots,\qquad 0<\theta<2\pi,$$

$$\sin\theta+\frac{\sin2\theta}{2}+\frac{\sin3\theta}{3}+\frac{\sin4\theta}{4}+\cdots,\qquad 0<\theta<2\pi;$$

(2)

$$\cos\theta+\frac{\cos3\theta}{3}+\frac{\cos5\theta}{5}+\frac{\cos7\theta}{7}+\cdots,\qquad 0<\theta<\pi,$$

$$\sin\theta+\frac{\sin3\theta}{3}+\frac{\sin5\theta}{5}+\frac{\sin7\theta}{7}+\cdots,\qquad 0\le\theta\le\pi;$$

(3)

$$\sin\theta-\frac{\sin3\theta}{3^2}+\frac{\sin5\theta}{5^2}-\frac{\sin7\theta}{7^2}+\cdots,\qquad -\frac\pi2\le\theta\le\frac\pi2;$$

(4)

$$\cos\theta-\frac{\cos5\theta}{5}+\frac{\cos7\theta}{7}-\frac{\cos11\theta}{11}+\cdots,\qquad -\frac\pi3<\theta<\frac\pi3.$$

提示：利用上题结果以及 Abel 第二定理。

8. 试求下列幂级数的收敛半径：

(1) $\displaystyle\sum_{n=1}^{\infty}\frac1{n^n}z^n$；

(2) $\displaystyle\sum_{n=1}^{\infty}\frac1{2^nn^n}z^n$；

(3) $\displaystyle\sum_{n=1}^{\infty}\frac{n!}{n^n}z^n$；

(4) $\displaystyle\sum_{n=1}^{\infty}\frac{(-)^n}{2^{2n}(n!)^2}z^n$；

(5) $\displaystyle\sum_{n=1}^{\infty}n^{\ln n}z^n$；

(6) $\displaystyle\sum_{n=1}^{\infty}\frac1{2^{2n}}z^{2n}$；

(7) $\displaystyle\sum_{n=1}^{\infty}\frac{\ln n^n}{n!}z^n$；

(8) $\displaystyle\sum_{n=1}^{\infty}\left(1-\frac1n\right)^nz^n$。

### 复数级数的判别法（补充）

1. Dirichlet 判别法：若级数 $\sum a_n$ 有界，而级数 $\sum(v_n-v_{n+1})$ 绝对收敛，且 $\lim v_n=0$，则级数 $\sum a_nv_n$ 收敛。

2. 如果

$$\frac{a_n}{a_{n+1}}=1+\frac{\mu}{n}+O\left(\frac1{n^{\lambda}}\right),$$

其中 $\mu=\alpha+i\beta,\ \lambda>1$，则

(1) 若 $\alpha>1$，则 $\sum a_n$ 收敛；

(2) 若 $\alpha=1$，且 $\beta\ne0$，则 $\sum a_n$ 振荡；

(3) 若 $0<\alpha<1$，则 $\sum a_n$ 发散，尽管仍有 $a_n\to0$；

(4) 若 $\alpha\le0$，则 $a_n$ 不趋于 0，因而 $\sum a_n$ 根本不可能收敛。

此结论首先由 Weierstrass 给出，它的特殊情形就是 Gauss 判别法。
