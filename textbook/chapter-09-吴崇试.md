# 第九章　二阶线性常微分方程的幂级数解法

幂级数解法是求解变系数线性常微分方程的一种实用方法。这种解法的理论称为常微分方程解析理论，其基础仍然是解析函数论。求解步骤类似于前面幂级数展开中介绍的待定系数法：首先通过对二阶线性常微分方程的分析，判断出方程的解在待求区域的级数形式（其中的叠加系数是待定常数），然后通过迭代的方法求出常微分方程幂级数解中的叠加系数。

## §9.1　二阶线性常微分方程的常点和奇点

在数学物理问题中，经常会出现一些二阶线性齐次常微分方程，它们的标准形式是

$$
\frac{d^2w}{dz^2}+p(z)\frac{dw}{dz}+q(z)w=0,
\tag{9.1}
$$

$p(z)$ 和 $q(z)$ 称为方程的系数。显然，方程的解是完全由方程的系数决定的。特别是，我们将看到，方程解的解析性是完全由方程系数的解析性决定的。

用级数解法解常微分方程时，得到的解总是某一指定点 $z_0$ 的邻域内收敛的无穷级数。方程系数 $p(z),q(z)$ 在 $z_0$ 点的解析性就决定了级数解在 $z_0$ 点的解析性，或者说，就决定了级数解的形式，例如，是 Taylor 级数还是 Laurent 级数。

如果方程的系数 $p(z)$ 和 $q(z)$ 都在 $z_0$ 点解析，则 $z_0$ 点称为方程的**常点**。

如果方程的系数 $p(z)$ 与 $q(z)$ 中至少有一个在 $z_0$ 点不解析，则 $z_0$ 点称为方程的**奇点**。

**例 9.1　超几何方程**

$$
z(1-z)\frac{d^2w}{dz^2}
+\bigl[\gamma-(1+\alpha+\beta)z\bigr]\frac{dw}{dz}
-\alpha\beta w=0
\tag{9.2}
$$

的系数是

$$
p(z)=\frac{\gamma-(1+\alpha+\beta)z}{z(1-z)},
\qquad
q(z)=-\frac{\alpha\beta}{z(1-z)}.
$$

在 $z$ 复平面 $\mathbb C$ 内，$p(z)$ 和 $q(z)$ 有两个孤立奇点：$z=0$ 和 $z=1$。所以，除了 $z=0$ 和 $z=1$ 是超几何方程的奇点外，复平面内其他点都是超几何方程的常点。

**例 9.2　Legendre 方程**

$$
(1-z^2)\frac{d^2w}{dz^2}
-2z\frac{dw}{dz}
+l(l+1)w=0,
\tag{9.3}
$$

在 $z$ 复平面内的奇点为 $z=\pm1$。

**练习 9.1**　已知二阶线性齐次常微分方程 (9.1) 的两个线性无关解是 $w_1(z)$ 和 $w_2(z)$，试证：

$$
p(z)=-\frac{\Delta_1(z)}{\Delta(z)},
\qquad
q(z)=\frac{\Delta_2(z)}{\Delta(z)},
$$

其中

$$
\Delta(z)=
\begin{vmatrix}
w_1(z)&w_1'(z)\\
w_2(z)&w_2'(z)
\end{vmatrix},
\qquad
\Delta_1(z)=
\begin{vmatrix}
w_1(z)&w_1''(z)\\
w_2(z)&w_2''(z)
\end{vmatrix},
$$

$$
\Delta_2(z)=
\begin{vmatrix}
w_1'(z)&w_1''(z)\\
w_2'(z)&w_2''(z)
\end{vmatrix}.
$$

要判断无穷远点 $z=\infty$ 是不是方程 (9.1) 的奇点，需将方程作自变量替换 $z=1/t$。因为

$$
\frac{dw}{dz}=-t^2\frac{dw}{dt},
\qquad
\frac{d^2w}{dz^2}
=
-t^2\frac{d}{dt}\left(-t^2\frac{dw}{dt}\right)
=
t^4\frac{d^2w}{dt^2}
+2t^3\frac{dw}{dt},
$$

所以关于 $z$ 的微分方程 (9.1) 变为关于 $t$ 的微分方程

$$
\frac{d^2w}{dt^2}
+
\left[
\frac2t-\frac1{t^2}p\left(\frac1t\right)
\right]\frac{dw}{dt}
+
\frac1{t^4}q\left(\frac1t\right)w
=
0.
\tag{9.4}
$$

如果 $t=0$ 是方程 (9.4) 的常点（奇点），则称无穷远点 $z=\infty$ 为方程 (9.1) 的常点（奇点）。$t=0$（即 $z=\infty$）为方程 (9.4) 常点的条件是

$$
p\left(\frac1t\right)
=
2t+a_2t^2+a_3t^3+\cdots,
\qquad
q\left(\frac1t\right)
=
b_4t^4+b_5t^5+\cdots,
\qquad
|t|<r,
\tag{9.5}
$$

也就是说当

$$
p(z)=\frac2z+\frac{a_2}{z^2}+\frac{a_3}{z^3}+\cdots,
\qquad
q(z)=\frac{b_4}{z^4}+\frac{b_5}{z^5}+\cdots,
\qquad
|z|>\frac1r
\tag{9.6}
$$

时，无穷远点 $z=\infty$ 是方程 (9.1) 的常点。由此可见，无穷远点是超几何方程 (9.2) 和 Legendre 方程 (9.3) 的奇点。

## §9.2　方程常点邻域内的解

我们不加证明地介绍下面的定理①。

**定理 9.1**　如果 $p(z)$ 和 $q(z)$ 在圆 $|z-z_0|<R$ 内单值解析，则在此圆内二阶线性齐次常微分方程的初值问题

$$
\frac{d^2w}{dz^2}
+p(z)\frac{dw}{dz}
+q(z)w=0,
\tag{9.7a}
$$

$$
w(z_0)=c_0,
\qquad
w'(z_0)=c_1
\quad
（其中\ c_0,c_1\ 为任意常数，）
\tag{9.7b}
$$

的解存在而且唯一，并且解 $w(z)$ 在这个圆内单值解析。

根据这个定理，方程的解 $w(z)$ 在常点 $z_0$ 的邻域 $|z-z_0|<R$ 内可以展开为 Taylor 级数

$$
w(z)
=
\sum_{k=0}^{\infty}
c_k(z-z_0)^k.
\tag{9.8}
$$

显然，$(z-z_0)^0$ 与 $(z-z_0)^1$ 的系数 $c_0$ 与 $c_1$ 正好和初值条件 (9.7b) 一致。定理 9.1 说明，系数 $c_k\ (k=2,3,\cdots)$ 一定均可用 $c_0,c_1$ 表示，而且表达式唯一。我们通过一个实例来说明具体的求解过程。

① 有关定理 9.1 及定理 9.2、9.3 的证明，可见参考书目 [12] 的 2.2、2.3 和 2.4 节。

**例 9.3　求 Legendre 方程**

$$
(1-z^2)\frac{d^2w}{dz^2}
-2z\frac{dw}{dz}
+l(l+1)w=0
\tag{9.9}
$$

在 $z=0$ 点邻域内的解，其中 $l$ 为已知参数。

**解**　显然 $z=0$ 是 Legendre 方程的常点。因此，可令方程在 $z=0$ 点邻域内的解为

$$
w=\sum_{k=0}^{\infty}c_kz^k,
\qquad
|z|<1.
\tag{9.10}
$$

代入方程 (9.9)，整理合并，就得到

$$
\sum_{k=0}^{\infty}
\left\{
(k+2)(k+1)c_{k+2}
-
[k(k+1)-l(l+1)]c_k
\right\}z^k
=
0.
$$

根据 Taylor 展开的唯一性，有

$$
(k+2)(k+1)c_{k+2}
-
[k(k+1)-l(l+1)]c_k
=
0.
$$

这样就得到了系数之间的递推关系

$$
c_{k+2}
=
\frac{k(k+1)-l(l+1)}{(k+2)(k+1)}c_k
=
\frac{(k-l)(k+l+1)}{(k+2)(k+1)}c_k.
\tag{9.11}
$$

反复利用递推关系，就可以求得系数

$$
\begin{aligned}
c_{2n}
&=
\frac{(2n-l-2)(2n+l-1)}
{2n(2n-1)}c_{2n-2}\\
&=
\frac{(2n-l-2)(2n-l-4)(2n+l-1)(2n+l-3)}
{2n(2n-1)(2n-2)(2n-3)}c_{2n-4}\\
&=\cdots\\
&=
\frac{c_0}{(2n)!}
(2n-l-2)(2n-l-4)\cdots(-l)
(2n+l-1)(2n+l-3)\cdots(l+1),
\end{aligned}
\tag{9.12}
$$

$$
\begin{aligned}
c_{2n+1}
&=
\frac{(2n-l-1)(2n+l)}
{(2n+1)(2n)}c_{2n-1}\\
&=
\frac{(2n-l-1)(2n-l-3)(2n+l)(2n+l-2)}
{(2n+1)(2n)(2n-1)(2n-2)}c_{2n-3}\\
&=\cdots\\
&=
\frac{c_1}{(2n+1)!}
(2n-l-1)(2n-l-3)\cdots(-l+1)
(2n+l)(2n+l-2)\cdots(l+2).
\end{aligned}
\tag{9.13}
$$

显然

$$
\lambda(\lambda+1)(\lambda+2)\cdots(\lambda+n-1)
=
\frac{\Gamma(\lambda+n)}{\Gamma(\lambda)},
$$

其中 $\Gamma(\lambda)$ 是 $\Gamma$ 函数（见 §7.2），则 $c_{2n}$ 和 $c_{2n+1}$ 可以写成

$$
c_{2n}
=
\frac{2^{2n}}{(2n)!}
\frac{
\Gamma\left(n-\dfrac l2\right)
}{
\Gamma\left(-\dfrac l2\right)
}
\frac{
\Gamma\left(\dfrac{l+1}{2}+n\right)
}{
\Gamma\left(\dfrac{l+1}{2}\right)
}
c_0,
\tag{9.14}
$$

$$
c_{2n+1}
=
\frac{2^{2n}}{(2n+1)!}
\frac{
\Gamma\left(n-\dfrac{l-1}{2}\right)
}{
\Gamma\left(-\dfrac{l-1}{2}\right)
}
\frac{
\Gamma\left(\dfrac l2+1+n\right)
}{
\Gamma\left(\dfrac l2+1\right)
}
c_1.
\tag{9.15}
$$

所以，Legendre 方程 (9.9) 在 $|z|<1$ 内的解就是

$$
w(z)=c_0w_1(z)+c_1w_2(z),
\qquad
|z|<1,
\tag{9.16}
$$

其中

$$
w_1(z)
=
\sum_{n=0}^{\infty}
\frac{2^{2n}}{(2n)!}
\frac{
\Gamma\left(n-\dfrac l2\right)
}{
\Gamma\left(-\dfrac l2\right)
}
\frac{
\Gamma\left(\dfrac{l+1}{2}+n\right)
}{
\Gamma\left(\dfrac{l+1}{2}\right)
}
z^{2n},
\tag{9.17}
$$

$$
w_2(z)
=
\sum_{n=0}^{\infty}
\frac{2^{2n}}{(2n+1)!}
\frac{
\Gamma\left(n-\dfrac{l-1}{2}\right)
}{
\Gamma\left(-\dfrac{l-1}{2}\right)
}
\frac{
\Gamma\left(\dfrac l2+1+n\right)
}{
\Gamma\left(\dfrac l2+1\right)
}
z^{2n+1}.
\tag{9.18}
$$

正如定理 9.1 所说，任意给定一组初条件 $c_0$ 和 $c_1$，就一定可以求出方程的一个特解。特别是，分别取 $c_0=1,c_1=0$ 和 $c_0=0,c_1=1$，得到的特解就是上面的 $w_1(z)$ 和 $w_2(z)$。这两个特解显然是线性无关的，也称为 Legendre 方程 (9.9) 的两个独立解。因此，解式 (9.16) 也就是 Legendre 方程 (9.9) 的通解，只要把常数 $c_0$ 和 $c_1$ 看成任意叠加常数即可。

上面求得的特解中，$w_1(z)$ 和 $w_2(z)$ 分别只含有 $z$ 的偶次幂和奇次幂，即 $w_1(z)$ 是 $z$ 的偶函数，$w_2(z)$ 是 $z$ 的奇函数。从求解的过程来看，这是由于递推关系中只出现系数 $c_{k+2}$ 和 $c_k$，而与 $c_{k+1}$ 无关，因此 $c_{2n}$ 完全由 $c_0$ 决定，$c_{2n+1}$ 完全由 $c_1$ 决定。从根本上来说，方程的解的对称性（这里指的是奇偶性）应该是方程的对称性的反映。事实上，在 Legendre 方程中令 $z\to-z$，就有

$$
[1-(-z)^2]
\frac{d^2w(-z)}{d(-z)^2}
-2(-z)\frac{dw(-z)}{d(-z)}
+l(l+1)w(-z)
=
0,
$$

即仍为

$$
(1-z^2)\frac{d^2w(-z)}{dz^2}
-2z\frac{dw(-z)}{dz}
+l(l+1)w(-z)
=
0.
$$

这说明，在变换 $z\to-z$ 之下，Legendre 方程的形式不变。所以，如果 $w(z)$ 是 Legendre 方程的解，$w(-z)$ 也一定是 Legendre 方程的解。因此 $w(z)\pm w(-z)$ 也是 Legendre 方程的解。显然，$w(z)+w(-z)$ 是 $z$ 的偶函数，$w(z)-w(-z)$ 是 $z$ 的奇函数。

通过这个实例，可以看出在常点邻域内求幂级数解的一般步骤：

(1) 将方程常点邻域内的解展开为 Taylor 级数，同时将方程的系数也在同一区域内展成幂级数，一并代入微分方程；

(2) 整理方程，比较系数，得到系数之间的递推关系；

(3) 反复利用递推关系，求出系数 $c_k$ 的普遍表达式（用 $c_0,c_1$ 表示），最后得到幂级数解。

由于递推关系一定是线性的（因为方程是线性的），所以最后的幂级数解（也正是微分方程的通解）一定可以写成

$$
w(z)=c_0w_1(z)+c_1w_2(z)
$$

的形式。$w_1(z)$ 和 $w_2(z)$ 是方程的两个线性无关解。

需要指出，在系数之间的递推关系中，一般会出现 $c_k,c_{k+1},c_{k+2}$ 三个相邻的系数，因此 $c_k$ 会同时依赖于 $c_0$ 和 $c_1$，最后求得的 $w_1(z)$ 或 $w_2(z)$ 就不会只含有 $z$ 的偶次幂或奇次幂。

应用常微分方程的幂级数解法，可以得到方程在一定区域内的解式。我们也可以根据需要，求出方程在不同区域内的解式。可以证明，方程在两个不同区域内的解式互为解析延拓。因此，也可从方程在某一区域内的解式出发，通过解析延拓，推出方程在其他区域内的解式。

**例 9.4**　设 $w_1$ 是方程

$$
\frac{d^2w}{dz^2}
+p(z)\frac{dw}{dz}
+q(z)w
=
0
\tag{9.19}
$$

的解，在区域 $G_1$ 内解析。若 $\widetilde w_1$ 是 $w_1$ 在区域 $G_2$ 内的解析延拓，不妨假设 $G_1\cap G_2$ 非空，即

$$
w_1\equiv\widetilde w_1,
\qquad
z\in G_1\cap G_2,
\tag{9.20}
$$

试证明：$\widetilde w_1$ 是方程在区域 $G_2$ 内的解。

**证**　设

$$
g(z)
=
\frac{d^2\widetilde w_1}{dz^2}
+p(z)\frac{d\widetilde w_1}{dz}
+q(z)\widetilde w_1,
$$

由题意知 $\widetilde w_1$ 在 $G_2$ 内解析，因此 $g(z)$ 在 $G_2$ 内解析。又因为 $w_1$ 是方程 (9.19) 在区域 $G_1$ 内的解，故在公共区域 $G_1\cap G_2$ 内，仍满足方程

$$
\frac{d^2w_1}{dz^2}
+p(z)\frac{dw_1}{dz}
+q(z)w_1
=
0.
$$

而在此子区域内，$w_1(z)\equiv\widetilde w_1(z)$，故

$$
\frac{d^2\widetilde w_1}{dz^2}
+p(z)\frac{d\widetilde w_1}{dz}
+q(z)\widetilde w_1
=
0,
\qquad
z\in G_1\cap G_2,
$$

即 $g(z)\equiv0,\ z\in G_1\cap G_2$。根据解析函数的唯一性，立即证得

$$
g(z)\equiv0,
\qquad
z\in G_2,
$$

亦即 $\widetilde w_1$ 在 $G_2$ 内满足方程

$$
\frac{d^2\widetilde w_1}{dz^2}
+p(z)\frac{d\widetilde w_1}{dz}
+q(z)\widetilde w_1
=
0.
$$

$\square$

**例 9.5**　设 $w_1$ 和 $w_2$ 都是方程 (9.19) 的两个线性无关解，且均在区域 $G_1$ 内解析。若 $\widetilde w_1$ 和 $\widetilde w_2$ 分别是 $w_1$ 和 $w_2$ 在区域 $G_2$ 内的解析延拓，即在 $z\in G_1\cap G_2$ 中

$$
w_1\equiv\widetilde w_1,
\qquad
w_2\equiv\widetilde w_2,
\tag{9.21}
$$

试证：$\widetilde w_1$ 和 $\widetilde w_2$ 线性无关。

**证**　由例 9.4 知，$\widetilde w_1$ 和 $\widetilde w_2$ 是方程 (9.19) 在 $G_2$ 内的解。因为 $w_1$ 和 $w_2$ 线性无关，

$$
\Delta[w_1,w_2]
\equiv
\begin{vmatrix}
w_1&w_2\\
w_1'&w_2'
\end{vmatrix}
\ne0,
\qquad
z\in G_1.
\tag{9.22}
$$

设

$$
\Delta[\widetilde w_1,\widetilde w_2]
\equiv
\begin{vmatrix}
\widetilde w_1&\widetilde w_2\\
\widetilde w_1'&\widetilde w_2'
\end{vmatrix}
=
g(z),
\qquad
z\in G_2,
\tag{9.23}
$$

$g(z)$ 在 $G_2$ 内解析。由于在 $z\in G_1\cap G_2$ 内时，

$$
w_1\equiv\widetilde w_1,
\qquad
w_2\equiv\widetilde w_2,
\tag{9.24}
$$

故 $g(z)\ne0,\ z\in G_1\cap G_2$。仍然根据解析函数的唯一性，就证得

$$
g(z)\ne0,
\qquad
z\in G_2.
\tag{9.25}
$$

所以，$\widetilde w_1$ 和 $\widetilde w_2$ 在 $G_2$ 内线性无关。$\square$

## §9.3　方程正则奇点邻域内的解

一般说来，方程的奇点（即方程系数的奇点）可能同时也是解的奇点。不但可能是解的极点或本性奇点，还可能是解的分支点。这里只讨论方程的奇点是**极点性奇点**（即方程的奇点是方程系数的极点）的情形。我们再次不加证明地介绍另一个定理：

**定理 9.2**　如果 $z_0$ 是二阶线性齐次常微分方程 (9.1) 的奇点，但是方程的系数 $p(z)$ 和 $q(z)$ 在环形区域 $0<|z-z_0|<R$ 内都解析，则方程 (9.1) 在环形区域 $0<|z-z_0|<R$ 内的线性无关解是

$$
w_1(z)
=
(z-z_0)^{\rho_1}
\sum_{k=-\infty}^{\infty}
c_k(z-z_0)^k,
\tag{9.26a}
$$

$$
w_2(z)
=
gw_1(z)\ln(z-z_0)
+
(z-z_0)^{\rho_2}
\sum_{k=-\infty}^{\infty}
d_k(z-z_0)^k,
\tag{9.26b}
$$

其中 $\rho_1,\rho_2$ 和 $g$ 都是常数。

显然，如果 $\rho_1$ 或 $\rho_2$ 不是整数（而且往往不是整数，否则就可以并入后面的级数中），或 $g\ne0$，方程的解均为多值函数，$z_0$ 为其分支点。

现在的问题是，如果我们把解 (9.26a) 或 (9.26b) 代入方程就会发现，尽管仍然能得到系数之间的递推关系，但却难以通过迭代的办法求出系数的普遍表达式。因为这时的幂级数解中，如果有无穷多个正幂项和负幂项，反复利用递推关系将会永无休止。

例外的情形是如果解式 (9.26) 中的级数只有有限个负幂项，这时总可以调整 $\rho$ 值，使得级数中没有负幂项，

$$
w_1(z)
=
(z-z_0)^{\rho_1}
\sum_{k=0}^{\infty}
c_k(z-z_0)^k,
$$

$$
w_2(z)
=
gw_1(z)\ln(z-z_0)
+
(z-z_0)^{\rho_2}
\sum_{k=0}^{\infty}
d_k(z-z_0)^k.
$$

这种形式的解称为**正则解**。方程在奇点邻域内有两个正则解的条件，见定理 9.3。

**定理 9.3**　二阶线性齐次常微分方程 (9.1) 在它的奇点 $z_0$ 的空心邻域 $0<|z-z_0|<R$ 内有两个正则解

$$
w_1(z)
=
(z-z_0)^{\rho_1}
\sum_{k=0}^{\infty}
c_k(z-z_0)^k,
\qquad
c_0\ne0,
\tag{9.27a}
$$

$$
w_2(z)
=
gw_1(z)\ln(z-z_0)
+
(z-z_0)^{\rho_2}
\sum_{k=0}^{\infty}
d_k(z-z_0)^k,
\qquad
g\text{ 或 }d_0\ne0
\tag{9.27b}
$$

的充分必要条件是 $(z-z_0)p(z)$ 和 $(z-z_0)^2q(z)$ 在 $z_0$ 点解析。

这样的奇点称为**正则奇点**。$\rho_1$ 和 $\rho_2$ 称为在方程正则奇点处（或正则解）的**指标**。

显然，$z=0$ 和 $z=1$ 都是超几何方程

$$
z(1-z)\frac{d^2w}{dz^2}
+
[\gamma-(1+\alpha+\beta)z]\frac{dw}{dz}
-
\alpha\beta w
=
0
$$

的正则奇点，$z=\pm1$ 也都是 Legendre 方程

$$
(1-z^2)\frac{d^2w}{dz^2}
-2z\frac{dw}{dz}
+l(l+1)w
=
0
$$

的正则奇点。

为了判断无穷远点是否为正则奇点，同样要作变量替换 $z=1/t$。如果 $t=0$ 是变换后的方程的正则奇点，即 $t=0$ 点是变换后的方程的奇点，且

$$
t\left[
\frac2t
-
\frac1{t^2}
p\left(\frac1t\right)
\right]
=
2-\frac1t p\left(\frac1t\right)
$$

和

$$
t^2\cdot
\frac1{t^4}
q\left(\frac1t\right)
=
\frac1{t^2}q\left(\frac1t\right)
$$

在 $t=0$ 点解析，亦即 $z=\infty$ 点是变换前方程的奇点，且 $zp(z)$ 和 $z^2q(z)$ 在 $z=\infty$ 点解析，则称 $z=\infty$ 点是变换前方程的正则奇点。

容易判断，$z=\infty$ 也都是超几何方程和 Legendre 方程的正则奇点。

根据定理 9.3，在方程的正则奇点处，我们应当将正则解 $w_1(z)$ 或 $w_2(z)$ 代入方程，比较系数，求出指标和递推关系，进而求出系数的普遍表达式，最终求得方程在正则奇点邻域内的两个正则解。但因为 $w_2(z)$ 的表达式中含有 $w_1(z)$，而且，如果 $g=0$，$w_2(z)$ 的表达式中不含对数函数项，解式的形式又和 $w_1(z)$ 相同，所以实际求解时，总是先将 $w_1(z)$ 形式的解代入方程。如果能够同时求得两个线性无关解，当然任务便告完成。如果这时只能求得一个解（例如 $\rho_1=\rho_2$ 时），那么，就必须再将 $w_2(z)$ 形式的解（这时的 $g$ 一定不为 0）代入方程求解。

下面简述一下求解过程。为了书写简单，不妨假设 $z=0$ 点是二阶线性齐次常微分方程 (9.1) 的正则奇点。于是，在 $z=0$ 点的空心邻域内，可将方程的系数展开为

$$
p(z)
=
\sum_{l=0}^{\infty}a_lz^{l-1},
\qquad
q(z)
=
\sum_{l=0}^{\infty}b_lz^{l-2},
\qquad
0<|z|<R.
$$

设解为

$$
w(z)
=
z^\rho
\sum_{k=0}^{\infty}c_kz^k,
\qquad
0<|z|<R,
$$

代入方程，就有

$$
\sum_{k=0}^{\infty}
c_k(k+\rho)(k+\rho-1)z^{k+\rho-2}
+
\sum_{l=0}^{\infty}
a_lz^{l-1}
\sum_{k=0}^{\infty}
c_k(k+\rho)z^{k+\rho-1}
+
\sum_{l=0}^{\infty}
b_lz^{l-2}
\sum_{k=0}^{\infty}
c_kz^{k+\rho}
=
0.
$$

整理，得

$$
\sum_{k=0}^{\infty}
\left\{
(k+\rho)(k+\rho-1)c_k
+
\sum_{l=0}^{k}
[a_{k-l}(l+\rho)+b_{k-l}]c_l
\right\}
z^k
=
0.
$$

或者单独提出 $k=0$ 项，

$$
c_0f_0(\rho)
+
\sum_{k=1}^{\infty}
\left[
f_0(k+\rho)c_k
+
\sum_{l=1}^{k}
f_l(k-l+\rho)c_{k-l}
\right]z^k
=
0,
$$

其中

$$
f_0(\rho)
=
\rho(\rho-1)+a_0\rho+b_0,
\qquad
f_k(\rho)
=
a_k\rho+b_k
\quad(k\geqslant1).
$$

比较等式两端最低次幂，即 $z^0$ 的系数，可得

$$
c_0f_0(\rho)
=
c_0[\rho(\rho-1)+a_0\rho+b_0]
=
0.
$$

由于 $c_0\ne0$，所以

$$
\boxed{
\rho(\rho-1)+a_0\rho+b_0=0.
}
\tag{9.28}
$$

这就是**指标方程**，注意其中的 $a_0$ 和 $b_0$ 为

$$
\boxed{
a_0=\lim_{z\to0}zp(z),
\qquad
b_0=\lim_{z\to0}z^2q(z).
}
\tag{9.29}
$$

根据指标方程可以求出两个指标，$\rho_1$ 和 $\rho_2$。规定 $\operatorname{Re}\rho_1\geqslant\operatorname{Re}\rho_2$。

再比较 $z^n$ 的系数，便可得到系数之间的递推关系：

$$
f_0(n+\rho)c_n
+
\sum_{l=1}^{n}
f_l(n-l+\rho)c_{n-l}
=
0.
$$

反复利用这个递推关系，就可以得到系数 $c_n$ 的普遍表达式（与 $\rho$ 有关）。分别用 $\rho=\rho_1$ 和 $\rho_2$ 代入，就可以得到解 $w_1(z)$ 和 $w_2(z)$。如果 $\rho_1-\rho_2\ne$ 整数，我们就求出了方程的两个（线性无关的）特解。

但是，需要注意，当 $\rho_1=\rho_2$ 时，这样只能得到同一个解，所以第二解一定含对数函数项。

当 $\rho_1-\rho_2=$ 正整数 $m$ 时，对于第二解的系数 $c_m^{(2)}$，有

$$
f_0(m+\rho_2)c_m^{(2)}
+
\sum_{l=1}^{m}
f_l(m-l+\rho_2)c_{m-l}^{(2)}
=
0.
$$

注意 $m+\rho_2=\rho_1$，所以 $f_0(m+\rho_2)=f_0(\rho_1)=0$，上式就是

$$
0\cdot c_m^{(2)}
+
\sum_{l=1}^{m}
f_l(m-l+\rho_2)c_{m-l}^{(2)}
=
0.
$$

因此

$$
\boxed{
\begin{aligned}
&\text{当 }\sum_{l=1}^{m}
f_l(m-l+\rho_2)c_{m-l}^{(2)}\ne0
\text{ 时，}\quad c_m^{(2)}\text{ 无解；}\\
&\text{当 }\sum_{l=1}^{m}
f_l(m-l+\rho_2)c_{m-l}^{(2)}=0
\text{ 时，}\quad c_m^{(2)}\text{ 任意。}
\end{aligned}
}
$$

对于第一种情形，方程的第二解一定含对数函数项，就需要用解式 (9.27b) 代入方程重新求解。对于第二种情形，当然还能继续求解，只是这时以后的各项系数 $c_n^{(2)}\ (n>m)$ 会同时依赖于 $c_0^{(2)}$ 和 $c_m^{(2)}$。第二解 $w_2(z)$ 便有两项，一项正比于 $c_0^{(2)}$，一项正比于 $c_m^{(2)}$。恰恰由于这时 $c_m^{(2)}$ 可取任意值（零或非零），这意味着与 $c_m^{(2)}$ 成正比的项一定是方程的解，而且就是指标为 $\rho_1$ 的第一解（最多可能差一个常数倍数），因而不妨就取 $c_m^{(2)}=0$。

总结以上讨论，就可以得出下列结论：

$$
\boxed{
\begin{array}{ll}
\text{当 }\rho_1-\rho_2\ne\text{ 整数时，}
&\text{第二解一定不含对数函数项；}\\[2pt]
\text{当 }\rho_1=\rho_2\text{ 时，}
&\text{第二解一定含对数函数项；}\\[2pt]
\text{当 }\rho_1-\rho_2=\text{ 正整数时，}
&\text{第二解可能含对数函数项。}
\end{array}
}
$$

需要指出，根据常微分方程的普遍理论，对于二阶线性齐次常微分方程 (9.1)，如果已经求出了一个解 $w_1(z)$，那么，总可以通过积分

$$
w_2(z)
=
Aw_1(z)
\int^z
\left\{
\frac1{[w_1(z)]^2}
\exp
\left[
-\int^z p(\zeta)\,d\zeta
\right]
\right\}dz
\tag{9.30}
$$

来求出第二解。这是因为 $w_1(z)$ 和 $w_2(z)$ 这两个解都满足方程

$$
\frac{d^2w_1}{dz^2}
+p(z)\frac{dw_1}{dz}
+q(z)w_1
=
0,
\qquad
\frac{d^2w_2}{dz^2}
+p(z)\frac{dw_2}{dz}
+q(z)w_2
=
0.
$$

将第一个方程乘以 $w_2$，第二个方程乘以 $w_1$，相减，便可得到

$$
w_1\frac{d^2w_2}{dz^2}
-
w_2\frac{d^2w_1}{dz^2}
+
p(z)
\left(
w_1\frac{dw_2}{dz}
-
w_2\frac{dw_1}{dz}
\right)
=
0,
$$

即

$$
\frac{d}{dz}
\left(
w_1\frac{dw_2}{dz}
-
w_2\frac{dw_1}{dz}
\right)
+
p(z)
\left(
w_1\frac{dw_2}{dz}
-
w_2\frac{dw_1}{dz}
\right)
=
0.
$$


再对 $z$ 积分，可得

$$
w_1\frac{dw_2}{dz}
-
w_2\frac{dw_1}{dz}
=
A\exp
\left[
-\int^z p(\zeta)\,d\zeta
\right].
\tag{9.31}
$$

两端除以 $w_1^2$，又可以得到

$$
\frac{d}{dz}
\left(
\frac{w_2}{w_1}
\right)
=
\frac{A}{w_1^2}
\exp
\left[
-\int^z p(\zeta)\,d\zeta
\right].
$$

再积分一次，就得到 (9.30) 式。

**例 9.6**　求方程

$$
(1-z^2)\frac{d^2w}{dz^2}
+
2(n-1)z\frac{dw}{dz}
+
2nw
=
0
$$

在 $z=1$ 点的有界解，其中 $n$ 为正整数。①

**解**　$z=1$ 为正则奇点，故可设

$$
w(z)
=
\sum_{k=0}^{\infty}
c_k(z-1)^{k+\rho},
$$

代入方程，有

$$
\begin{aligned}
&\bigl[(z-1)^2+2(z-1)\bigr]
\sum_{k=0}^{\infty}
c_k(k+\rho)(k+\rho-1)(z-1)^{k+\rho-2}\\
&\quad
-2(n-1)\bigl[(z-1)+1\bigr]
\sum_{k=0}^{\infty}
c_k(k+\rho)(z-1)^{k+\rho-1}
-2n
\sum_{k=0}^{\infty}
c_k(z-1)^{k+\rho}
=
0.
\end{aligned}
$$

消去 $(z-1)^\rho$，即可将方程整理为

$$
\sum_{k=0}^{\infty}
c_k(k+\rho+1)(k+\rho-2n)(z-1)^{k+1}
+
2
\sum_{k=0}^{\infty}
c_k(k+\rho)(k+\rho-n)(z-1)^k
=
0.
$$

比较 $(z-1)^0$ 项的系数，可以得到指标方程

$$
\rho(\rho-n)=0,
$$

由此求得指标

$$
\rho_1=n,
\qquad
\rho_2=0.
$$

再比较 $(z-1)^m$ 项的系数，得递推关系

$$
c_{m-1}(m+\rho)(m+\rho-2n-1)
+
2c_m(m+\rho)(m+\rho-n)
=
0,
$$

即

$$
c_m
=
-\frac12
\frac{m+\rho-2n-1}{m+\rho-n}
c_{m-1}.
$$

对于 $\rho_1=n$，

$$
c_m
=
-\frac12
\frac{m-n-1}{m}
c_{m-1}
=
\frac12
\frac{n+1-m}{m}
c_{m-1}.
$$

反复利用此递推关系，即可求得

$$
c_k
=
\frac1{2^k}
\frac{(n+1-k)(n+2-k)\cdots n}{k!}c_0
=
\frac1{2^k}
\frac{n!}{k!(n-k)!}c_0.
$$

显然，当 $k>n$ 时，$c_k=0$。因此，

$$
w_1(z)
=
c_0(z-1)^n
\sum_{k=0}^{n}
\frac{n!}{k!(n-k)!}
\left(
\frac{z-1}{2}
\right)^k.
$$

不妨引进新的常数 $A=2^{-n}c_0$，则有

$$
w_1(z)
=
A(z-1)^n
\sum_{k=0}^{n}
\frac{n!}{k!(n-k)!}
2^{n-k}(z-1)^k
=
A(z^2-1)^n.
$$

对于 $\rho_2=0$，

$$
c_m
=
-\frac12
\frac{m-2n-1}{m-n}
c_{m-1}.
$$

显然 $c_n$ 无意义，这说明 $w_2(z)$ 一定含有对数函数项，因此在 $z=1$ 点无界，非本题所要求。其实，按照 (9.30) 式，可以写出

$$
\begin{aligned}
w_2(x)
&=
Bw_1(x)
\int^x
\frac1{w_1^2(x)}
\exp
\left[
-\int^x
\frac{2(n-1)t}{1-t^2}\,dt
\right]dx\\
&=
B(x^2-1)^n
\int^x
\frac1{(x^2-1)^{2n}}
e^{(n-1)\ln(x^2-1)}\,dx\\
&=
B(x^2-1)^n
\int^x
\frac{dx}{(x^2-1)^{n+1}}.
\end{aligned}
$$

另一方面，如果将方程改写为

$$
(z^2-1)\frac{d^2w}{dz^2}
+
2z\frac{dw}{dz}
-
2n\left(
z\frac{dw}{dz}+w
\right)
=
0,
$$

即

$$
\frac{d}{dz}
\left[
(z^2-1)\frac{dw}{dz}
\right]
-
2n\frac{d(zw)}{dz}
=
0,
$$

我们可以直接积分一次，得

$$
(z^2-1)\frac{dw}{dz}
-
2nzw
=
c_1.
$$

这是非齐次方程，相应齐次方程（即取 $c_1=0$）为

$$
\frac1w\frac{dw}{dz}
=
\frac{2nz}{z^2-1},
$$

再积分一次，即得

$$
\ln w
=
n\ln(z^2-1)+A,
\qquad
\text{即}
\qquad
w=A(z^2-1)^n.
$$

这正是上面求得的第一解。在此基础上，采用常数变易法，即可求出第二解。

① 此方程微商 $n$ 次，就是 Legendre 方程。

## §9.4　Bessel 方程的解

**Bessel 方程**

$$
\frac{d^2w}{dz^2}
+
\frac1z\frac{dw}{dz}
+
\left(
1-\frac{\nu^2}{z^2}
\right)w
=
0
\tag{9.32}
$$

是常见的常微分方程之一，其中 $\nu$ 是常数，$\operatorname{Re}\nu\geqslant0$。容易判断，$z=0$ 是方程的正则奇点，$z=\infty$ 是方程的非正则奇点。

本节讨论 Bessel 方程在正则奇点 $z=0$ 点的空心邻域 $|z|>0$ 内的解。设

$$
w(z)
=
z^\rho
\sum_{k=0}^{\infty}
c_kz^k,
\qquad
c_0\ne0,
\tag{9.33}
$$

代入方程 (9.32)，得

$$
\sum_{k=0}^{\infty}
c_k(k+\rho)(k+\rho-1)z^{k+\rho-2}
+
\sum_{k=0}^{\infty}
c_k(k+\rho)z^{k+\rho-2}
+
\sum_{k=0}^{\infty}
c_kz^{k+\rho}
-
\nu^2
\sum_{k=0}^{\infty}
c_kz^{k+\rho-2}
=
0,
$$

约去 $z^{\rho-2}$，即得

$$
\sum_{k=0}^{\infty}
c_k[(k+\rho)^2-\nu^2]z^k
+
\sum_{k=0}^{\infty}
c_kz^{k+2}
=
0.
$$

根据级数展开的唯一性，即可比较系数。

由最低次幂 $z^0$ 项的系数，且因为 $c_0\ne0$，就得到指标方程

$$
\rho^2-\nu^2=0,
\tag{9.34}
$$

因而求得

$$
\rho_1=\nu,
\qquad
\rho_2=-\nu.
\tag{9.35}
$$

因为 $\operatorname{Re}\nu\geqslant0$，所以满足 $\operatorname{Re}\rho_1\geqslant\operatorname{Re}\rho_2$。

由 $z^1$ 的系数，得

$$
c_1[(\rho+1)^2-\nu^2]=0,
\qquad
\text{即}
\qquad
c_1(2\rho+1)=0,
$$

因此

$$
c_1=0,
\qquad
\text{当 }\rho\ne-\frac12;
\tag{9.36a}
$$

$$
c_1\text{ 任意},
\qquad
\text{当 }\rho=-\frac12.
\tag{9.36b}
$$

以后将看到，即使 $\rho=-1/2$，仍可以取 $c_1=0$。

由 $z^n$ 的系数，得

$$
c_n[(\rho+n)^2-\nu^2]
+
c_{n-2}
=
0,
$$

即

$$
c_nn(2\rho+n)+c_{n-2}=0,
$$

因此，得到递推关系

$$
c_n
=
-\frac1{n(n+2\rho)}c_{n-2}.
\tag{9.37}
$$

反复利用递推关系，就可以求得

$$
\begin{aligned}
c_{2n}
&=
-\frac1{n(n+\rho)}\frac1{2^2}c_{2n-2}\\
&=
(-)^2
\frac1{n(n-1)(n+\rho)(n+\rho-1)}
\frac1{2^4}c_{2n-4}\\
&=\cdots\\
&=
\frac{(-)^n}{n!}
\frac{\Gamma(\rho+1)}{\Gamma(\rho+1+n)}
\frac1{2^{2n}}c_0,
\end{aligned}
\tag{9.38}
$$

$$
c_{2n+1}
=
-\frac1{(n+1/2)(n+\rho+1/2)}
\frac1{2^2}c_{2n-1},
\tag{9.39}
$$

$$
\begin{aligned}
c_{2n+1}
&=
\frac{(-)^2}
{(n+1/2)(n-1/2)(n+\rho+1/2)(n+\rho-1/2)}
\frac1{2^4}c_{2n-3}\\
&=\cdots\\
&=
(-)^n
\frac{\Gamma(3/2)}
{\Gamma(n+3/2)}
\frac{\Gamma(\rho+3/2)}
{\Gamma(\rho+3/2+n)}
\frac1{2^{2n}}c_1
=
0.
\end{aligned}
\tag{9.40}
$$

用 $\rho_1=\nu$ 代入，即得

$$
w_1(z)
=
c_0z^\nu
\sum_{k=0}^{\infty}
\frac{(-)^k\Gamma(\nu+1)}
{k!\Gamma(k+\nu+1)}
\left(
\frac z2
\right)^{2k}.
\tag{9.41}
$$

取

$$
c_0
=
\frac1{2^\nu\Gamma(\nu+1)},
$$

就有解

$$
J_\nu(z)
=
\sum_{k=0}^{\infty}
\frac{(-)^k}{k!\Gamma(k+\nu+1)}
\left(
\frac z2
\right)^{2k+\nu}.
\tag{9.42}
$$

用 $\rho_2=-\nu$ 代入，有

$$
w_2(z)
=
c_0z^{-\nu}
\sum_{k=0}^{\infty}
\frac{(-)^k\Gamma(1-\nu)}
{k!\Gamma(k-\nu+1)}
\left(
\frac z2
\right)^{2k}.
\tag{9.43}
$$

取

$$
c_0
=
\frac{2^\nu}{\Gamma(-\nu+1)},
$$

又得

$$
J_{-\nu}(z)
=
\sum_{k=0}^{\infty}
\frac{(-)^k}{k!\Gamma(k-\nu+1)}
\left(
\frac z2
\right)^{2k-\nu}.
\tag{9.44}
$$

现在补充讨论一下 $\rho=-1/2$ 的情形。前面曾经提到，这时仍然可以取 $c_1=0$。因为如果 $c_1\ne0$，则

$$
c_{2n+1}
=
\frac{(-)^n\Gamma(3/2)}
{\Gamma(3/2+n)\Gamma(n+1)}
\frac1{2^{2n}}c_1.
$$

注意到 $\Gamma(n+1)=n!$，这样在 $w_2(z)$ 中只不过是再增加一项

$$
\begin{aligned}
z^{-1/2}
\sum_{n=0}^{\infty}
c_{2n+1}z^{2n+1}
&=
c_1
\sum_{n=0}^{\infty}
\frac{(-)^n}
{n!\Gamma(n+3/2)}
\left(
\frac z2
\right)^{2n+1/2}
\sqrt{\frac\pi2}\\
&=
c_1\sqrt{\frac\pi2}J_{1/2}(z),
\end{aligned}
$$

即在 $w_2(z)$ 中只不过是再叠加上第一解。

但是，当 $\nu=0$ 时，上面的求解过程只给出了同一个解

$$
J_0(z)
=
\sum_{k=0}^{\infty}
\frac{(-)^k}{k!k!}
\left(
\frac z2
\right)^{2k}.
\tag{9.45}
$$

不仅如此，当 $\nu=n,\ n=1,2,3,\cdots$ 时，上面求出的仍然也只是一个解。这可以从下面三个角度来说明。第一，从递推关系 (9.37) 来看，当 $\nu=n$ 时，显然 $c_{2n}$ 无意义，以后各项系数因而也都失去意义。第二，上面在导出 (9.44) 式时，曾经取 $c_0=2^\nu/\Gamma(1-\nu)$，这在 $\nu=n,\ n=1,2,\cdots$ 是不合法的，因为这意味着 $c_0=0$，违反了幂级数解中首项系数不为 0 的约定。第三，退一步说，即使许可取 $c_0=0$，因为这样得到的幂级数解

$$
J_{-n}(z)
=
\sum_{k=0}^{\infty}
\frac{(-)^k}
{k!\Gamma(k-n+1)}
\left(
\frac z2
\right)^{2k-n}
$$

仍然有意义，但由于 $\Gamma(k-n+1)=\infty,\ k=0,1,\cdots,n-1$，所以

$$
\begin{aligned}
J_{-n}(z)
&=
\sum_{k=n}^{\infty}
\frac{(-)^k}
{k!\Gamma(k-n+1)}
\left(
\frac z2
\right)^{2k-n}\\
&=
\sum_{l=0}^{\infty}
\frac{(-)^{n+l}}
{(n+l)!\Gamma(l+1)}
\left(
\frac z2
\right)^{2(n+l)-n}\\
&=
(-)^n
\sum_{l=0}^{\infty}
\frac{(-)^l}
{l!\Gamma(n+l+1)}
\left(
\frac z2
\right)^{2l+n}\\
&=
(-)^nJ_n(z),
\end{aligned}
\tag{9.46}
$$

和第一解 $J_n(z)$ 线性相关，所以也并没有给出新的解。总之，当 $\nu=0,n,\ n=0,1,2,\cdots$ 时，以上都只是求得了第一解。这说明第二解一定含有对数函数项，即

$$
w_2(z)
=
gJ_n(z)\ln z
+
\sum_{k=0}^{\infty}
d_kz^{k-n},
\qquad
g\ne0.
\tag{9.47}
$$

从原则上说，将 $w_2(z)$ 及 $J_n(z)$ 的级数表达式代入 Bessel 方程

$$
\frac{d^2w}{dz^2}
+
\frac1z\frac{dw}{dz}
+
\left(
1-\frac{n^2}{z^2}
\right)w
=
0,
$$

化简整理得

$$
g
\sum_{k=0}^{\infty}
\frac{(-)^k}{k!(k+n)!}
\frac{2k+n}{2^{2k+n-1}}
z^{2k+2n}
+
\sum_{k=0}^{\infty}
d_kk(k-2n)z^k
+
\sum_{k=0}^{\infty}
d_kz^{k+2}
=
0.
$$

逐一比较 $z$ 的各次幂的系数，定出 $g$ 和 $d_k,\ k=0,1,2,\cdots$，即可求得 $w_2(z)$。不再赘述。

下面介绍求 Bessel 方程第二解的另一种方法。为此，先计算 $J_\nu(z)$ 和 $J_{-\nu}(z)$ 的 Wroński 行列式以分析它们的线性相关性。考虑到 Bessel 方程的系数 $p(z)=1/z$，从 (9.31) 式就可得到

$$
W[J_\nu(z),J_{-\nu}(z)]
\equiv
\begin{vmatrix}
J_\nu(z)&J_{-\nu}(z)\\
J_\nu'(z)&J_{-\nu}'(z)
\end{vmatrix}
=
A\exp
\left[
-\int^z\frac{d\zeta}{\zeta}
\right]
=
\frac Az.
$$

为了定出积分常数 $A$，只需将 $J_\nu(z)$ 和 $J_{-\nu}(z)$ 的幂级数展开式 (9.42) 和 (9.44) 代入，找出

$$
W[J_\nu(z),J_{-\nu}(z)]
\equiv
J_\nu(z)J_{-\nu}'(z)
-
J_{-\nu}(z)J_\nu'(z)
$$

中 $z^{-1}$ 项的系数即可。这只来自各级数中的第一项。因此，

$$
\begin{aligned}
A
&=
\frac1{\Gamma(1+\nu)}
\frac1{2^\nu}
\frac1{\Gamma(1-\nu)}
\frac1{2^{-\nu}}
(-\nu)\\
&\qquad
-
\frac1{\Gamma(1-\nu)}
\frac1{2^{-\nu}}
\frac1{\Gamma(1+\nu)}
\frac1{2^\nu}
\nu\\
&=
-\frac{2\nu}{\Gamma(1+\nu)\Gamma(1-\nu)}
=
-\frac2{\Gamma(\nu)\Gamma(1-\nu)}
=
-\frac2\pi\sin\pi\nu.
\end{aligned}
$$

这样就得到

$$
W[J_\nu(z),J_{-\nu}(z)]
=
-\frac2{\pi z}\sin\pi\nu.
\tag{9.48}
$$

上面的计算中用到了 Γ 函数的性质（见 (7.11) 式）

$$
\Gamma(\nu)\Gamma(1-\nu)
=
\frac{\pi}{\sin\pi\nu}.
$$

(9.48) 式再次表明，当 $\nu=n,\ n=0,1,2,\cdots$ 时，$J_\nu(z)$ 和 $J_{-\nu}(z)$ 线性相关。但是如果将 Bessel 方程的第二解取为 $J_\nu(z)$ 和 $J_{-\nu}(z)$ 的线性组合

$$
w_2(z)
=
c_1J_\nu(z)
+
c_2J_{-\nu}(z),
$$

则

$$
W[J_\nu(z),w_2(z)]
=
c_2W[J_\nu(z),J_{-\nu}(z)]
=
-\frac{2c_2}{\pi z}\sin\pi\nu.
$$

只要选择适当的组合系数 $c_2$，使得 $W[J_\nu(z),w_2(z)]$ 对任何 $\nu$ 均不为 0，这样对任何 $\nu$，$w_2(z)$ 就一定与 $J_\nu(z)$ 线性无关。为此，我们就取 $c_2=-1/\sin\pi\nu$，即取第二解为

$$
w_2(z)
=
\frac{cJ_\nu(z)-J_{-\nu}(z)}
{\sin\pi\nu},
$$

这样便有

$$
W[J_\nu(z),w_2(z)]
=
\frac2{\pi z}.
\tag{9.49}
$$

为了保证这样定义的 $w_2(z)$ 有意义（$\sin n\pi=0$，分母为 0），并注意到 $J_{-n}(z)=(-)^nJ_n(z)$，我们便应当进一步选取系数 $c$，使得 $w_2(z)$ 中的分子在 $\nu=n$ 时也为 0，例如取 $c=\cos\pi\nu$ 即可。这样得到 Bessel 方程的第二解便是

$$
N_\nu(z)
=
\frac{
\cos\pi\nu\,J_\nu(z)-J_{-\nu}(z)
}{
\sin\pi\nu
},
\tag{9.50}
$$

称为 $\nu$ 阶 Neumann 函数。当 $\nu=n,\ n=0,1,2,\cdots$，(9.50) 式为不定式，可按 l’Hôpital 法则求极限，

$$
\begin{aligned}
N_n(z)
&=
\lim_{\nu\to n}N_\nu(z)
=
\lim_{\nu\to n}
\frac{
\cos\pi\nu\,J_\nu(z)-J_{-\nu}(z)
}{
\sin\pi\nu
}\\
&=
\left.
\frac1\pi
\left[
\frac{\partial J_\nu(z)}{\partial\nu}
-
(-)^n
\frac{\partial J_{-\nu}(z)}{\partial\nu}
\right]
\right|_{\nu=n}\\
&=
\frac2\pi J_n(z)\ln\frac z2
-
\frac1\pi
\sum_{k=0}^{n-1}
\frac{(n-k-1)!}{k!}
\left(
\frac z2
\right)^{2k-n}\\
&\quad
-
\frac1\pi
\sum_{k=0}^{\infty}
\frac{(-)^k}
{k!(n+k)!}
[\psi(n+k+1)+\psi(k+1)]
\left(
\frac z2
\right)^{2k+n},
\qquad
|\arg z|<\pi,
\end{aligned}
\tag{9.51}
$$

其中 $\psi(\zeta)$ 是 Γ 函数的对数微商（见 (7.17) 式），

$$
\psi(\zeta)
\equiv
\frac{d\ln\Gamma(\zeta)}{d\zeta}
=
\frac{\Gamma'(\zeta)}{\Gamma(\zeta)},
$$

并且约定，当 $n=0$ 时应当去掉 (9.51) 式右端第二项的有限和。

## *§9.5　方程非正则奇点附近的解

前面讨论了在方程正则奇点邻域内的解。我们知道，在正则奇点的邻域内，二阶线性齐次常微分方程 (9.1) 的两个线性无关解都是正则解。如果 $z=z_0$ 是方程的非正则奇点，则在该点的邻域 $0<|z-z_0|<R$ 内最多只能有一个正则解。如果 $p(z)$ 和 $q(z)$ 可以在 $0<|z-z_0|<R$ 内作 Laurent 展开

$$
p(z)
=
(z-z_0)^{-m}
\sum_{k=0}^{\infty}
a_k(z-z_0)^k,
\qquad
q(z)
=
(z-z_0)^{-n}
\sum_{k=0}^{\infty}
b_k(z-z_0)^k,
$$

将正则解 (9.27a) 以及 $p(z)$ 和 $q(z)$ 的 Laurent 级数代入方程，则当且仅当 $m$ 和 $n$ 满足

$$
\boxed{
m\geqslant2,
\qquad
\text{并且}
\qquad
m\geqslant n-1
}
\tag{9.52}
$$

时，得到的指标方程为一次方程。这时可以模仿正则奇点处的求解步骤求解。不再赘述。

如果无穷远点 $z=\infty$ 是非正则奇点，原则上应当作变量替换 $z=1/t$，然后在 $t=0$ 点的邻域内讨论即可。但是，也可以直接令

$$
w(z)
=
z^\rho
\sum_{k=0}^{\infty}
c_kz^{-k},
$$

同时将 $p(z)$ 和 $q(z)$ 在 $z=\infty$ 点的空心邻域内作 Laurent 展开

$$
p(z)
=
z^m
\sum_{k=0}^{\infty}
a_kz^{-k},
\qquad
q(z)
=
z^n
\sum_{k=0}^{\infty}
b_kz^{-k},
$$

代入微分方程，则得到

$$
\sum_{k=0}^{\infty}
c_k(k-\rho)(k-\rho-1)z^{-k}
+
z^{1+m}
\sum_{l=0}^{\infty}
a_lz^{-l}
\sum_{k=0}^{\infty}
c_k(k-\rho)z^{-k}
+
z^{2+n}
\sum_{l=0}^{\infty}
b_lz^{-l}
\sum_{k=0}^{\infty}
c_kz^{-k}
=
0.
$$

所以，方程只有一个正则解的条件是

$$
\boxed{
m\geqslant0,
\qquad
\text{并且}
\qquad
m\geqslant n+1,
}
\tag{9.52'}
$$

即在 $z=\infty$ 点 $p(z)$ 的阶大于 $q(z)$ 的阶。注意，方程只有一个正则解的条件，对于有限远处和无穷远点，在形式上有所不同。

现在讨论另一种情形，即方程非正则奇点邻域内不存在正则解，或方程有一个正则解，但我们要求这个正则解之外的另一解。这时，$z=z_0$ 应该是解的本性奇点（当然还可能是分支点）。为了描写解在 $z=z_0$ 点的这种奇异性，不妨假设

$$
w(z)
=
e^{Q(z)}v(z),
\tag{9.53}
$$

使得 $z=z_0$ 是 $e^{Q(z)}$ 的本性奇点，而 $v(z)$ 可以写成正则解的形式：

$$
v(z)
=
(z-z_0)^\rho
\sum_{k=0}^{\infty}
c_k(z-z_0)^k,
\qquad
c_0\ne0.
\tag{9.54}
$$

(9.53) 这种形式的解称为**常规解**。把这种解代入方程，则得 $v(z)$ 的方程

$$
\frac{d^2v}{dz^2}
+
p^*(z)\frac{dv}{dz}
+
q^*(z)v
=
0,
\tag{9.55}
$$

其中

$$
p^*(z)
=
p(z)+2Q'(z),
$$

$$
q^*(z)
=
q(z)+p(z)Q'(z)+Q''(z)+[Q'(z)]^2.
$$

再根据方程有正则解的要求，就有可能恰当地选择 $Q(z)$，并进而求出解 $w(z)$。

如果求常规解的尝试失败，还可以求次常规解（详见参考书目 [12] 的 2.11 节）。

下面再来讨论 Bessel 方程 (9.32)。容易看出，在非正则奇点 $z=\infty$ 处，$m=-1,\ n=0$，不满足方程有正则解的要求 (9.52')，因而在 $z=\infty$ 的邻域内，方程最多只能有常规解。作变换 (9.53)，则方程 (9.55) 中的系数为

$$
p^*(z)
=
\frac1z+2Q'(z),
\tag{9.56}
$$

$$
q^*(z)
=
1-\frac{\nu^2}{z^2}
+
\frac{Q'(z)}z
+
Q''(z)
+
[Q'(z)]^2.
\tag{9.57}
$$

为了保证 $z=\infty$ 是 $e^{Q(z)}$ 的本性奇点，可取 $Q(z)$ 为 $z$ 的多项式，这也同时满足了 $p^*(z)$ 的阶 $m\geqslant0$ 的要求；再进一步，为了使得 $p^*(z)$ 的阶大于 $q^*(z)$ 的阶，则必须取

$$
Q(z)=\lambda z,
\qquad
\text{并且}
\qquad
1+\lambda^2=0,
$$

即 $\lambda=\pm i$。这样，$v(z)$ 所满足的方程就是

$$
\frac{d^2v}{dz^2}
+
\left(
\frac1z+2\lambda
\right)\frac{dv}{dz}
+
\left(
\frac{\lambda}{z}
-
\frac{\nu^2}{z^2}
\right)v
=
0.
\tag{9.58}
$$

设

$$
v(z)
=
z^\rho
\sum_{l=0}^{\infty}
c_lz^{-l},
\tag{9.59}
$$

则有

$$
\sum_{l=0}^{\infty}
[(\rho-l)^2-\nu^2]c_lz^{-l}
+
\lambda
\sum_{l=0}^{\infty}
[(2\rho+1)-2l]c_lz^{-l+1}
=
0.
\tag{9.60}
$$

比较 $z^1$ 项的系数，即得

$$
2\rho+1=0,
\qquad
\text{即}
\qquad
\rho=-1/2.
$$

所以，在无穷远点 $z=\infty$ 附近，Bessel 方程解的首项便是 $c_0e^{\pm iz}\sqrt{1/z}$。

为了求出整个解，可将 $\rho=-1/2$ 代入 (9.60) 式，再比较 $z^{-k+1}$ 项的系数，就得到递推关系

$$
c_k
=
-\frac{4\nu^2-(2k-1)^2}{2^2k}
\frac1{2\lambda}
c_{k-1}.
\tag{9.61}
$$

反复利用递推关系，就可以得到系数 $c_k$ 的普遍表达式：

$$
c_k
=
\frac{(\nu,k)}{(-2\lambda)^k}c_0,
\tag{9.62}
$$

其中 $(\nu,0)=1$，

$$
(\nu,k)
=
\frac{
[4\nu^2-(2k-1)^2]
[4\nu^2-(2k-3)^2]
\cdots
[4\nu^2-3^2]
[4\nu^2-1^2]
}{
2^{2k}k!
}.
$$

通常令

$$
c_0
=
\sqrt{\frac2{\pi}}
\exp
\left[
-\lambda
\left(
\frac{\nu\pi}{2}
+
\frac\pi4
\right)
\right],
\tag{9.63}
$$

于是就可以求得 Bessel 方程在无穷远点 $z=\infty$ 邻域内的解。

但是，容易判断，这个幂级数解的收敛半径

$$
R
=
\lim_{k\to\infty}
\left|
\frac{c_{k-1}}{c_k}
\right|
=
|2\lambda|
\times
\lim_{k\to\infty}
\left|
\frac{(\nu,k-1)}{(\nu,k)}
\right|
=
0.
$$

这说明，在一般情况下，这样得到的幂级数解是发散的，除非 $\nu$ 是半奇数，即 $\nu=n+1/2$，这时从 $z^{-n-1}$ 项起，各项的系数均为 0，级数截断为多项式，

$$
w(z)
=
\sqrt{\frac2{\pi z}}
\exp
\left[
\lambda
\left(
z-\frac{n+1}{2}\pi
\right)
\right]
\sum_{k=0}^{n}
\frac{(n+1/2,k)}{(-2\lambda z)^k}.
$$

将 $\lambda=\pm i$ 代入，就得到 Bessel 方程的两个解：

$$
w_1(z)
=
(-i)^{n+1}
\sqrt{\frac2{\pi z}}
e^{iz}
\sum_{k=0}^{n}
\frac{(-)^k(n+1/2,k)}
{(2iz)^k},
\tag{9.64}
$$

$$
w_2(z)
=
i^{n+1}
\sqrt{\frac2{\pi z}}
e^{-iz}
\sum_{k=0}^{n}
\frac{(n+1/2,k)}
{(2iz)^k}.
\tag{9.65}
$$

如果 $\nu$ 不是半奇数，得到的级数实际上是 Bessel 方程的解（事实上，是 Hankel 函数）当 $|z|\to\infty$ 时的渐近展开：

$$
H_\nu^{(1)}(z)
\sim
\sqrt{\frac2{\pi z}}
\exp
\left[
i\left(
z-\frac{\nu\pi}{2}-\frac\pi4
\right)
\right]
\sum_{k=0}^{\infty}
\frac{(-)^k(\nu,k)}
{(2iz)^k},
$$

$$
H_\nu^{(2)}(z)
\sim
\sqrt{\frac2{\pi z}}
\exp
\left[
-i\left(
z-\frac{\nu\pi}{2}-\frac\pi4
\right)
\right]
\sum_{k=0}^{\infty}
\frac{(\nu,k)}
{(2iz)^k},
\qquad
(-2\pi<\arg z<\pi).
\tag{9.66}
$$

将它们重新组合，还可以得到 Bessel 函数和 Neumann 函数在 $|z|\to\infty,\ -\pi<\arg z<\pi$ 时的渐近展开：

$$
\begin{aligned}
J_\nu(z)
&=
\frac{
H_\nu^{(1)}(z)+H_\nu^{(2)}(z)
}{2}\\
&\sim
\sqrt{\frac2{\pi z}}
\left[
\cos
\left(
z-\frac{\nu\pi}{2}-\frac\pi4
\right)
\sum_{k=0}^{\infty}
\frac{(-)^k(\nu,2k)}
{(2z)^{2k}}
\right.\\
&\qquad\qquad\left.
-
\sin
\left(
z-\frac{\nu\pi}{2}-\frac\pi4
\right)
\sum_{k=0}^{\infty}
\frac{(-)^k(\nu,2k+1)}
{(2z)^{2k+1}}
\right],
\end{aligned}
\tag{9.67}
$$

$$
\begin{aligned}
N_\nu(z)
&=
\frac{
H_\nu^{(1)}(z)-H_\nu^{(2)}(z)
}{2i}\\
&\sim
\sqrt{\frac2{\pi z}}
\left[
\sin
\left(
z-\frac{\nu\pi}{2}-\frac\pi4
\right)
\sum_{k=0}^{\infty}
\frac{(-)^k(\nu,2k)}
{(2z)^{2k}}
\right.\\
&\qquad\qquad\left.
+
\cos
\left(
z-\frac{\nu\pi}{2}-\frac\pi4
\right)
\sum_{k=0}^{\infty}
\frac{(-)^k(\nu,2k+1)}
{(2z)^{2k+1}}
\right].
\end{aligned}
\tag{9.68}
$$

## 习题

1. 求二阶线性常微分方程，使其解为：

   (1) $w_1(z)=z,\ w_2(z)=e^z$；

   (2) $w_1(z)=e^{1/z},\ w_2(z)=e^{-2/z}$；

   (3) $w_1(z)=\cos\dfrac az,\ w_2(z)=\sin\dfrac az$；

   (4) $w_1(z)=\dfrac{z^2}{z^2-1},\ w_2(z)=\dfrac{z}{z^2-1}$。

2. 求下列方程在 $z=0$ 邻域内的两个幂级数解：

   (1)
   $$
   w''-z^2w=0;
   $$

   (2)
   $$
   w''-zw=0;
   $$

   (3)
   $$
   (z^2-1)w''+zw'-w=0;
   $$

   (4)
   $$
   (1+z+z^2)w''+2(1+2z)w'+2w=0.
   $$

3. 求下列方程在 $z=0$ 邻域内的两个幂级数解：

   (1)
   $$
   z^2(1-z)w''+z(1-3z)w'-(1+z)w=0;
   $$

   (2)
   $$
   9z^2w''-15zw'+(36z^4+7)w=0;
   $$

   (3)
   $$
   zw''-zw'+w=0;
   $$

   (4)
   $$
   zw''+(z-1)w'+w=0.
   $$

4. 求方程

$$
\frac{d^2u}{dz^2}
+
\frac2z\frac{du}{dz}
+
m^2u
=
0
$$

在 $z=0$ 附近的两个独立解。

5. 求方程

$$
\frac{d^2w}{dz^2}
+
\frac1z\frac{dw}{dz}
-
m^2w
=
0
$$

在 $z=0$ 附近的两个独立解。
