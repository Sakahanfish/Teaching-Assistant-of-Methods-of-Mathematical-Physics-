# 第二十章 变分法初步

变分法是数学物理中一个古老而成熟的方法。早在古希腊时期，有人就已经提出了“等周问题”。17 世纪，许多数学家用不同方法，解决了 Galileo 提出的著名“捷线问题”，而 Euler 和 Lagrange 对于这一类问题的研究，则为变分法奠定了理论基础。

变分法在物理学中的应用集中在两个方面：基本物理规律的一种新的表述语言和具体物理问题的近似方法。本章就围绕这两方面的应用，对变分法做一点初步的介绍。

## §20.1 泛函的概念

泛函是从函数到数的映射，简单地，就是以整个函数为自变量的函数。

设在 $x-y$ 平面上有一簇曲线 $y(x)$，其长度为

$$
L=\int_C ds=\int_{x_0}^{x_1}\sqrt{1+y'^2}\,dx.
$$

显然，$y(x)$ 不同，$L$ 值也不同，即 $L$ 的数值依赖于整个函数 $y(x)$。$L$ 和函数 $y(x)$ 之间的这种依赖关系称为泛函关系。类似的例子还有闭合曲线围成的面积，平面曲线绕固定轴而生成的旋转体体积或表面积，等等。它们也都定义了各自的泛函关系。

设对于（某一函数集合内的）每一个函数 $y(x)$，都有唯一一个数 $J[y]$ 与之对应，则称 $J[y]$ 为 $y(x)$ 的泛函。这里的函数集合，即泛函的定义域，通常要求 $y(x)$ 具有连续的二阶导数，并且满足一定的边界条件。这样的 $y(x)$ 称为可取函数。

这里要特别强调，泛函不同于复合函数，例如 $g=g(f(x))$。对于后者，给定一个 $x$ 值，仍然有一个 $g$ 值与之对应；对于前者，则必须给出某一区间上的整个函数 $y(x)$，才能得到一个泛函值 $J[y]$。（定义在同一区间上的）函数不同，泛函值当然不同。为了强调泛函值 $J[y]$ 与函数 $y(x)$ 之间的依赖关系，常常又把函数 $y(x)$ 称为变量函数。

**例 20.1　捷线问题（brachistochrone problem）.** 如图 20.1 所示，在重力作用下，一质点从 $(x_0,y_0)$ 点沿平面曲线 $y(x)$ 无摩擦地自由下滑到 $(x_1,y_1)$ 点，则所需时间

$$
T
=
\int_{(x_0,y_0)}^{(x_1,y_1)}
\frac{ds}{\sqrt{2g(y_0-y)}}
=
\int_{x_0}^{x_1}
\frac{\sqrt{1+y'^2}}{\sqrt{2g(y_0-y)}}\,dx.
\tag{20.1}
$$

就是 $y(x)$ 的泛函。这里，自然要求变量函数 $y(x)$ 一定通过端点 $(x_0,y_0)$ 和 $(x_1,y_1)$。

![Figure 20.1](...)

图 20.1　捷线问题

**例 20.2　弦的横振动问题.** 设在弦上隔离出足够短的一段弦，则该段弦的动能 $E_k$ 和势能 $E_p$ 分别为

$$
E_k
=
\frac12\rho\Delta x
\left(
\frac{\partial u}{\partial t}
\right)^2,
\qquad
E_p
=
\frac12T\Delta x
\left(
\frac{\partial u}{\partial x}
\right)^2,
$$

其中 $u(x,t)$ 是弦的横向位移，$\rho$ 是弦的线密度，$T$ 是张力。这样，弦的 Hamilton 作用量

$$
S
=
\int_{t_0}^{t_1}dt
\int_{x_0}^{x_1}
\frac12
\left[
\rho
\left(
\frac{\partial u}{\partial t}
\right)^2
-
T
\left(
\frac{\partial u}{\partial x}
\right)^2
\right]
 dx
\tag{20.2}
$$

也是位移 $u(x,t)$ 的泛函。这里的

$$
L
=
\int_{x_0}^{x_1}
\frac12
\left[
\rho
\left(
\frac{\partial u}{\partial t}
\right)^2
-
T
\left(
\frac{\partial u}{\partial x}
\right)^2
\right]
 dx
$$

称为 Lagrange 量，而被积函数

$$
\frac12
\left[
\rho
\left(
\frac{\partial u}{\partial t}
\right)^2
-
T
\left(
\frac{\partial u}{\partial x}
\right)^2
\right]
$$

称为 Lagrange 量密度。

泛函的形式可以多种多样，例如求函数的最大值，也是定义了一个泛函关系。但本书中只限于用积分定义的泛函，例如

$$
J[y]
=
\int_{x_0}^{x_1}
F(x,y,y')\,dx,
\tag{20.3}
$$

其中 $F$ 是它的宗量 $x,y(x)$ 及 $y'(x)$ 的已知函数，且具有连续的二阶偏导数。如果变量函数是二元函数 $u(x,y)$，则泛函为

$$
J[u]
=
\iint_S
F(x,y,u,u_x,u_y)\,dxdy,
\tag{20.4}
$$

其中 $u_x\equiv\partial u/\partial x$，$u_y\equiv\partial u/\partial y$。对于更多个自变量的多元函数，也可以有类似的定义。

## §20.2 泛函的极值

泛函极值的概念和函数极值非常类似。例如，“当变量函数为 $y(x)$ 时，泛函 $J[y]$ 取极小值”的含义就是：对于极值函数 $y(x)$ 及其“附近”的变量函数 $y(x)+\delta y(x)$，恒有

$$
J[y+\delta y]\ge J[y].
\tag{20.5}
$$

所谓 $y(x)+\delta y(x)$ 在 $y(x)$ 的“附近”，是函数空间的语言，其含义是它们的差 $\delta y(x)$（称为函数 $y(x)$ 的变分）满足

(1) $\|\delta y(x)\|<\varepsilon$；

(2) 有时还要求 $\|(\delta y)'(x)\|<\varepsilon$，

其中 $\|\delta y(x)\|$ 就是 (17.10) 式定义的函数 $\delta y(x)$ 的范数。

可以仿照函数极值必要条件的导出办法，导出泛函取极值的必要条件。为此，不妨不失普遍性地假定，所考虑的变量函数均通过固定的两个端点 $y(x_0)=a$，$y(x_1)=b$，即

$$
\delta y(x_0)=0,
\qquad
\delta y(x_1)=0.
\tag{20.6}
$$

现在，考虑泛函的差值

$$
J[y+\delta y]-J[y]
=
\int_{x_0}^{x_1}
\left[
F(x,y+\delta y,y'+(\delta y)')
-
F(x,y,y')
\right]dx.
$$

当函数的变分 $\delta y(x)$ 足够小时，可以将被积函数在极值函数附近作 Taylor 展开，有

$$
\begin{aligned}
J[y+\delta y]-J[y]
&=
\int_{x_0}^{x_1}
\left\{
\left[
\delta y\frac{\partial}{\partial y}
+
(\delta y)'\frac{\partial}{\partial y'}
\right]F
\right.\\
&\qquad\left.
+
\frac1{2!}
\left[
\delta y\frac{\partial}{\partial y}
+
(\delta y)'\frac{\partial}{\partial y'}
\right]^2F
+
\cdots
\right\}dx\\
&=
\delta J[y]
+
\frac1{2!}\delta^2J[y]
+
\cdots,
\end{aligned}
$$

其中

$$
\delta J[y]
\equiv
\int_{x_0}^{x_1}
\left[
\frac{\partial F}{\partial y}\delta y
+
\frac{\partial F}{\partial y'}(\delta y)'
\right]dx,
\tag{20.7}
$$

$$
\begin{aligned}
\delta^2J[y]
&\equiv
\int_{x_0}^{x_1}
\left[
\delta y\frac{\partial}{\partial y}
+
(\delta y)'\frac{\partial}{\partial y'}
\right]^2F\,dx\\
&=
\int_{x_0}^{x_1}
\left[
\frac{\partial^2F}{\partial y^2}(\delta y)^2
+
2\frac{\partial^2F}{\partial y\partial y'}\delta y(\delta y)'
+
\frac{\partial^2F}{\partial y'^2}(\delta y)'^2
\right]dx
\end{aligned}
\tag{20.8}
$$

分别是泛函 $J[y]$ 的一级变分和二级变分。这样就得到：泛函 $J[y]$ 取极小值的必要条件是泛函的一级变分为 0，

$$
\delta J[y]
\equiv
\int_{x_0}^{x_1}
\left[
\frac{\partial F}{\partial y}\delta y
+
\frac{\partial F}{\partial y'}(\delta y)'
\right]dx
=
0.
\tag{20.9}
$$

对于泛函 $J[y]$ 取极大值的情形，也可以类似地讨论，并且也会得到相同形式的必要条件。

将 (20.9) 式的积分中的第二项分部积分，同时考虑到边界条件 (20.6)，就有

$$
\begin{aligned}
\delta J[y]
&=
\left.
\frac{\partial F}{\partial y'}\delta y
\right|_{x_0}^{x_1}
+
\int_{x_0}^{x_1}
\left(
\frac{\partial F}{\partial y}
-
\frac d{dx}\frac{\partial F}{\partial y'}
\right)
\delta y\,dx\\
&=
\int_{x_0}^{x_1}
\left(
\frac{\partial F}{\partial y}
-
\frac d{dx}\frac{\partial F}{\partial y'}
\right)
\delta y\,dx
=
0.
\end{aligned}
$$

由于 $\delta y$ 的任意性，我们又可以得到

$$
\frac{\partial F}{\partial y}
-
\frac d{dx}
\frac{\partial F}{\partial y'}
=
0.
\tag{20.10}
$$

这个方程称为 Euler-Lagrange 方程，它是泛函 $J[y]$ 取极小值的必要条件的微分形式。一般说来，这是一个二阶常微分方程。

在导出方程 (20.10) 时，需要用到下面的引理：

**引理 20.1（变分学基本引理）**　设 $\phi(x)$ 是连续函数，$\eta(x)$ 具有连续的二阶导数，并且满足 $\eta(x)|_{x=x_0}=\eta(x)|_{x=x_1}=0$，若对于任意 $\eta(x)$，

$$
\int_{x_0}^{x_1}\phi(x)\eta(x)\,dx=0
$$

均成立，则必有 $\phi(x)\equiv0$。证明从略。

**例 20.3**　设质点在有势力场中沿路径 $q=q(t)$ 由 $(t_0,q(t_0))$ 点运动到 $(t_1,q(t_1))$ 点，它的 Hamilton 作用量是

$$
S
=
\int_{t_0}^{t_1}
L(t,q,\dot q)\,dt,
\tag{20.11}
$$

其中 $q$ 和 $\dot q$ 是描写质点运动的广义坐标和广义动量，$L=T-V$ 是动能 $T$ 和势能 $V$ 之差，即 Lagrange 量。Hamilton 原理告诉我们，在一切（运动学上允许的）可能路径中，真实发生的（即由力学规律决定的）运动路径使作用量 $S$ 取极值。根据上面的讨论可知，作用量 $S$ 取极值的必要条件的积分形式和微分形式分别是

$$
\delta S
=
\int_{t_0}^{t_1}
\left(
\frac{\partial L}{\partial q}\delta q
+
\frac{\partial L}{\partial\dot q}\delta\dot q
\right)dt
=
0
\tag{20.12}
$$

和

$$
\frac{\partial L}{\partial q}
-
\frac d{dt}
\frac{\partial L}{\partial\dot q}
=
0.
\tag{20.13}
$$

在给定的有势力场中，写出 Lagrange 量 $L$ 的具体形式，代入 (20.13) 式就会发现，它和 Newton 力学的动力学方程完全一样。

现在讨论两种常见的特殊情形。一种是泛函 (20.3) 中的 $F=F(x,y')$ 不显含 $y$，这时的 Euler-Lagrange 方程就是

$$
\frac d{dx}
\frac{\partial F}{\partial y'}
=
0,
$$

所以，立即就得到它的首次积分

$$
\frac{\partial F}{\partial y'}
=
\text{常量 }C.
\tag{20.14}
$$

另一种是泛函 (20.3) 中的 $F=F(y,y')$ 不显含 $x$，容易证明，

$$
\frac d{dx}
\left(
y'\frac{\partial F}{\partial y'}-F
\right)
=
y''\frac{\partial F}{\partial y'}
+
y'\frac d{dx}
\frac{\partial F}{\partial y'}
-
\frac{\partial F}{\partial y}y'
-
\frac{\partial F}{\partial y'}y''
=
-y'
\left(
\frac{\partial F}{\partial y}
-
\frac d{dx}
\frac{\partial F}{\partial y'}
\right),
$$

所以，这时的 Euler-Lagrange 方程也可以有首次积分

$$
y'\frac{\partial F}{\partial y'}-F
=
\text{常量 }C.
\tag{20.15}
$$

把这个结果应用到例 20.3 中，如果 Lagrange 量 $L$ 不显含 $t$，则有

$$
\dot q\frac{\partial L}{\partial\dot q}-L
=
\text{常量 }C,
\tag{20.16}
$$

这就是能量守恒。

下面研究变量函数是二元函数的情形。设有二元函数 $u(x,y)$，$(x,y)\in S$，在此基础上可定义泛函

$$
J[u]
=
\iint_S
F(x,y,u,u_x,u_y)\,dxdy.
\tag{20.17}
$$

仍然约定，$u(x,y)$ 在 $S$ 的边界 $\Gamma$ 上的数值给定，即

$$
\left.u\right|_\Gamma
\text{ 固定}.
\tag{20.18}
$$

首先，当然要计算

$$
\begin{aligned}
J[u+\delta u]-J[u]
&=
\iint_S
F(x,y,u+\delta u,(u+\delta u)_x,(u+\delta u)_y)
\,dxdy\\
&\quad-
\iint_S
F(x,y,u,u_x,u_y)\,dxdy\\
&=
\iint_S
\left[
\delta u\frac{\partial}{\partial u}
+
(\delta u)_x\frac{\partial}{\partial u_x}
+
(\delta u)_y\frac{\partial}{\partial u_y}
\right]F\,dxdy\\
&\quad+
\frac1{2!}
\iint_S
\left[
\delta u\frac{\partial}{\partial u}
+
(\delta u)_x\frac{\partial}{\partial u_x}
+
(\delta u)_y\frac{\partial}{\partial u_y}
\right]^2F\,dxdy
+
\cdots.
\end{aligned}
$$

于是，泛函 $J[u]$ 取极值的必要条件就是泛函的一级变分为 0，

$$
\begin{aligned}
\delta J[u]
&=
\iint_S
\left[
\delta u\frac{\partial F}{\partial u}
+
(\delta u)_x\frac{\partial F}{\partial u_x}
+
(\delta u)_y\frac{\partial F}{\partial u_y}
\right]dxdy\\
&=
\iint_S
\left[
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\left(
\frac{\partial F}{\partial u_x}
\right)
-
\frac{\partial}{\partial y}
\left(
\frac{\partial F}{\partial u_y}
\right)
\right]
\delta u\,dxdy\\
&\quad+
\iint_S
\left[
\frac{\partial}{\partial x}
\left(
\frac{\partial F}{\partial u_x}\delta u
\right)
+
\frac{\partial}{\partial y}
\left(
\frac{\partial F}{\partial u_y}\delta u
\right)
\right]dxdy
=
0.
\end{aligned}
\tag{20.19}
$$

利用公式

$$
\iint_S
\left(
\frac{\partial Q}{\partial x}
-
\frac{\partial P}{\partial y}
\right)dxdy
=
\int_\Gamma
(P\,dx+Q\,dy),
$$

取 $Q=\dfrac{\partial F}{\partial u_x}\delta u$，$P=-\dfrac{\partial F}{\partial u_y}\delta u$，就能将上面的结果化为

$$
\delta J[u]
=
\iint_S
\left(
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\frac{\partial F}{\partial u_x}
-
\frac{\partial}{\partial y}
\frac{\partial F}{\partial u_y}
\right)
\delta u\,dxdy
+
\int_\Gamma
\left(
-\frac{\partial F}{\partial u_y}\,dx
+
\frac{\partial F}{\partial u_x}\,dy
\right)
\delta u.
$$

根据 (20.18) 式，$\left.\delta u\right|_\Gamma=0$，可知上式右端第二项的线积分为 0，所以

$$
\delta J[u]
=
\iint_S
\left(
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\frac{\partial F}{\partial u_x}
-
\frac{\partial}{\partial y}
\frac{\partial F}{\partial u_y}
\right)
\delta u\,dxdy
=
0.
$$

再利用 $\delta u$ 的任意性，就可以导出上面的被积函数一定为 0，即

$$
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\frac{\partial F}{\partial u_x}
-
\frac{\partial}{\partial y}
\frac{\partial F}{\partial u_y}
=
0.
\tag{20.20}
$$

这就是二元函数情形下，泛函

$$
J[u]
=
\iint_S
F(x,y,u,u_x,u_y)\,dxdy
$$

取极值的必要条件的微分形式（Euler-Lagrange 方程）。

把这个结果应用到例 20.2 中弦的横振动问题上，就得到使作用量

$$
S
=
\int_{t_0}^{t_1}dt
\int_{x_0}^{x_1}
\frac12
\left[
\rho
\left(
\frac{\partial u}{\partial t}
\right)^2
-
T
\left(
\frac{\partial u}{\partial x}
\right)^2
\right]dx
$$

取极值的必要条件

$$
\frac{\partial^2u}{\partial t^2}
-
\frac T\rho
\frac{\partial^2u}{\partial x^2}
=
0,
\tag{20.21}
$$

正是我们在第 11 章导出的弦的横振动方程。

**练习 20.1**　在 $n$ 个自变量的情形下，导出泛函

$$
J[u]
=
\int\cdots\int
F(x_1,x_2,\cdots,x_n,u,u_{x_1},u_{x_2},\cdots,u_{x_n})
\,dx_1dx_2\cdots dx_n
$$

取极值的必要条件，包括它的积分形式和微分形式。上式中的积分是在 $n$ 维空间中的一定区域内进行的。

下面以一元函数为例，总结一下变分运算的几条简单法则。

(1) 变量函数 $y$ 本身只有一级变分，

$$
\delta^2y=0.
\tag{20.22}
$$

(2) 变分运算是线性运算，

$$
\delta(\alpha F+\beta G)
=
\alpha\delta F+\beta\delta G,
\tag{20.23}
$$

其中 $\alpha$ 和 $\beta$ 是常数。

(3) 直接计算，就可以得到函数乘积的变分法则：

$$
\delta(FG)
=
(\delta F)G+F(\delta G).
\tag{20.24}
$$

(4) 变分运算和微分或微商运算可交换次序，

$$
\delta\frac{dy}{dx}
=
\frac{d(\delta y)}{dx},
\qquad
\text{即}
\qquad
\delta y'
=
(\delta y)'.
\tag{20.25}
$$

(5) 变分运算和积分（微分的逆运算）也可以交换次序，

$$
\delta
\int_a^b
F\,dx
=
\int_a^b
(\delta F)\,dx.
\tag{20.26}
$$

(6) 复合函数的变分运算，其法则和微分运算完全相同，只要简单地将微分运算中的“$d$”换成“$\delta$”即可。例如，

$$
\delta F(x,y,y')
=
\frac{\partial F}{\partial y}\delta y
+
\frac{\partial F}{\partial y'}\delta y'.
\tag{20.27}
$$

这里注意，引起 $F$ 变化的原因，是函数 $y$ 的变分，而非自变量 $x$。

这些运算法则，可以毫不困难地推广到多元函数的情形。

作为完整的泛函极值问题，在列出泛函取极值的必要条件（即 Euler-Lagrange 方程）后，还需要在给定的定解条件下求解微分方程，才求得极值函数。这里需要注意，Euler-Lagrange 方程只是泛函取极值的必要条件，并非充分条件。在给定的定解条件下，Euler-Lagrange 方程的解到底是不是极值函数，还需要进一步加以甄别。和求函数极值的情形一样，现在也可以有两种方法。一种是直接比较所求得的解及其“附近”的函数的泛函值，根据泛函极值的定义加以判断。另一种方法是计算泛函的二级变分 $\delta^2J$，如果对于所求得的解，泛函的二级变分取正（负）值，则该解即为极值函数，泛函取极小（大）。如果二级变分为 0，则需要继续讨论高级变分。

可是，实际问题往往又特别简单，这就是在给定的边界条件下，Euler-Lagrange 方程只有一个解，同时，从物理或数学内容上又能判断，该泛函的极值一定存在，那么，这时求得的唯一解当然就是所要求的极值函数了。

## §20.3 泛函的条件极值

要求泛函

$$
J[y]
=
\int_{x_0}^{x_1}
F(x,y,y')\,dx
\tag{20.28}
$$

在边界条件

$$
y(x_0)=a,
\qquad
y(x_1)=b
\tag{20.29}
$$

以及约束条件

$$
J_1[y]
\equiv
\int_{x_0}^{x_1}
G(x,y,y')\,dx
=
C
\tag{20.30}
$$

下的极值，可采用处理多元函数条件极值问题的 Lagrange 乘子法，即定义

$$
J_0[y]
=
J[y]-\lambda J_1[y],
\tag{20.31}
$$

其中 Lagrange 乘子 $\lambda$ 为待定常数，仍将 $\delta y$ 看成是独立的，于是泛函 $J_0[y]$ 在边界条件 (20.29) 下取极值的必要条件就是

$$
\left(
\frac{\partial}{\partial y}
-
\frac d{dx}
\frac{\partial}{\partial y'}
\right)
(F-\lambda G)
=
0.
\tag{20.32}
$$

由微分方程 (20.32)、边界条件 (20.29) 以及约束条件 (20.30)，必要时经过甄别，就可以求出 Lagrange 乘子的值 $\lambda=\lambda_0$、极值函数 $y=y(x,\lambda_0)$，及相应泛函 $J_0[y]$ 的条件极值。

**例 20.4**　求泛函

$$
I[y]
=
\int_0^1xy'^2\,dx
\tag{20.33}
$$

在边界条件

$$
y(0)\text{ 有界},
\qquad
y(1)=0
\tag{20.34}
$$

和约束条件

$$
\int_0^1xy^2\,dx
=
1
\tag{20.35}
$$

下的极值曲线。

**解**　采用上面描述的 Lagrange 乘子法，可以得到必要条件

$$
\left(
\frac{\partial}{\partial y}
-
\frac d{dx}
\frac{\partial}{\partial y'}
\right)
(xy'^2-\lambda xy^2)
=
0,
\qquad
\text{即}
\qquad
\frac d{dx}
\left(
x\frac{dy}{dx}
\right)
+
\lambda xy
=
0.
\tag{20.36}
$$

此方程及齐次的边界条件 (20.34) 即构成一个本征值问题，它的本征值

$$
\lambda_i
=
\mu_i^2,
\qquad
\mu_i\text{ 是零阶 Bessel 函数 }J_0(x)\text{ 的第 }i\text{ 个正零点},
\quad
i=1,2,3,\cdots
\tag{20.37}
$$

正好就是 Lagrange 乘子，而极值函数就是相应的本征函数

$$
y_i(x)
=
C_iJ_0(\mu_i x).
$$

常量 $C_i$ 可以由约束条件定出。因为

$$
C_i^2
\int_0^1
xJ_0^2(\mu_i x)\,dx
=
\frac12
C_i^2J_1^2(\mu_i)
=
1,
$$

所以 $C_i=\sqrt2/J_1(\mu_i)$。这样就求出了极值函数

$$
y_i(x)
=
\frac{\sqrt2}{J_1(\mu_i)}
J_0(\mu_i x).
\tag{20.38}
$$

值得注意，这里由于 Lagrange 乘子的引进，在 Euler-Lagrange 方程出现了待定参量，和齐次边界条件组合在一起，就构成本征值问题。而作为本征值问题，它的解（本征值和本征函数）有无穷多个。这无穷多个本征函数都是极值函数，因为 $I[y]$ 的二级变分

$$
\delta^2I[y]
=
2
\int_0^1
x(\delta y')^2\,dx
>
0,
$$

所以这些极值函数均使泛函取极小。而且，这无穷个本征值正好也就是泛函的极值。这是因为，将方程 (20.36) 乘以极值函数 $y(x)$，再积分，就有

$$
\lambda
\int_0^1xy^2\,dx
=
-
\int_0^1
y(xy')'\,dx
=
-
\left.y\cdot xy'\right|_0^1
+
\int_0^1xy'^2\,dx
=
\int_0^1xy'^2\,dx.
$$

根据约束条件 (20.35)，就能得到

$$
\lambda
=
\int_0^1xy'^2\,dx.
\tag{20.39}
$$

还可以从另一个角度来理解泛函的条件极值问题。从第十七章的讨论知道，自伴算符本征函数的全体构成完备函数组。因此，作为泛函 (20.33) 中的可取函数 $y(x)$，一定可以按照本征函数 (20.38) 展开（不妨设展开系数 $c_i$ 为实数）：

$$
y(x)
=
\sum_{i=1}^{\infty}
c_i y_i(x).
\tag{20.40}
$$

容易证明：

$$
\begin{aligned}
\int_0^1xy'^2\,dx
&=
\left.xy'y\right|_0^1
-
\int_0^1(xy')'y\,dx\\
&=
-
\int_0^1
\sum_{i=1}^{\infty}
c_i(xy_i')'
\sum_{j=1}^{\infty}
c_jy_j\,dx\\
&=
\sum_{i,j}
c_ic_j\lambda_i
\int_0^1y_i y_jx\,dx
=
\sum_{i=1}^{\infty}
\lambda_i c_i^2,
\end{aligned}
\tag{20.41}
$$

$$
\int_0^1xy^2\,dx
=
\sum_{i=1}^{\infty}c_i^2
=
1.
\tag{20.42}
$$

因此，例 20.4 中的泛函条件极值问题就完全等价于（无穷维）二次型 (20.41) 在约束条件 (20.42) 下的极值问题。

最后，还要提到，这一类泛函的条件极值问题的原型，可以追溯到“闭合曲线周长一定而面积取极大”的原始几何问题。因此，泛函的条件极值问题，常称为等周问题（isoperimetric problem）。

## §20.4 微分方程定解问题和本征值问题的变分形式

在前两节中读者看到，泛函取极值的必要条件的微分形式（Euler-Lagrange 方程）是常微分方程或偏微分方程，它和变量函数的定解条件结合起来，就构成常微分方程或偏微分方程的定解问题。对于泛函的条件极值问题，其必要条件中出现待定参量（Lagrange 乘子），它和齐次边界条件结合起来，就构成微分方程本征值问题。这一节将研究它的反问题：如何将微分方程的定解问题或本征值问题转化为泛函的极值或条件极值问题，即如何将微分方程的定解问题或本征值问题用变分语言表述。

**例 20.5**　写出常微分方程边值问题

$$
\frac d{dx}
\left[
p(x)\frac{dy}{dx}
\right]
+
q(x)y(x)
=
f(x),
\qquad
x_0<x<x_1,
\tag{20.43a}
$$

$$
y(x_0)=y_0,
\qquad
y(x_1)=y_1
\tag{20.43b}
$$

的泛函形式，即求相应的泛函，它在边界条件 (20.43b) 下取极值的必要条件即为 (20.43a)。

**解**　既然泛函极值必要条件的微分形式就是方程 (20.43a)，那么，这个方程一定来自

$$
\int_{x_0}^{x_1}
\left\{
\frac d{dx}
\left[
p(x)\frac{dy}{dx}
\right]
+
q(x)y(x)
-
f(x)
\right\}
\delta y(x)\,dx
=
0.
$$

现在的问题就是要把上式左端化成泛函的变分。对于上式左端积分中被积函数的第二项和第三项，显然有

$$
\int_{x_0}^{x_1}
q(x)y(x)\delta y(x)\,dx
=
\frac12\delta
\int_{x_0}^{x_1}
q(x)y^2(x)\,dx,
$$

$$
\int_{x_0}^{x_1}
f(x)\delta y(x)\,dx
=
\delta
\int_{x_0}^{x_1}
f(x)y(x)\,dx.
$$

注意已知函数 $q(x)$ 和 $f(x)$ 是与 $y(x)$ 的变分无关的，因此它们在变分计算中都是常量。对于被积函数中的第一项，可以通过分部积分而化为

$$
\int_{x_0}^{x_1}
\frac d{dx}
\left[
p(x)\frac{dy}{dx}
\right]
\delta y(x)\,dx
=
-
\int_{x_0}^{x_1}
p(x)\frac{dy}{dx}
\delta\left(
\frac{dy}{dx}
\right)dx
=
-
\frac12\delta
\int_{x_0}^{x_1}
p(x)
\left(
\frac{dy}{dx}
\right)^2dx,
$$

其中用到了 $\left.\delta y(x)\right|_{x_0}=\left.\delta y(x)\right|_{x_1}=0$。综合上面的结果，就得到

$$
\begin{aligned}
&\int_{x_0}^{x_1}
\left\{
\frac d{dx}
\left[
p(x)\frac{dy}{dx}
\right]
+
q(x)y(x)
-
f(x)
\right\}
\delta y(x)\,dx\\
&=-\delta
\int_{x_0}^{x_1}
\left\{
\frac12
\left[
p(x)
\left(
\frac{dy}{dx}
\right)^2
-
q(x)y^2(x)
\right]
+
f(x)y(x)
\right\}dx
=
0.
\end{aligned}
\tag{20.44}
$$

这就说明，方程 (20.43a) 一定就是泛函

$$
J[y]
=
\int_{x_0}^{x_1}
\left\{
\frac12
\left[
p(x)
\left(
\frac{dy}{dx}
\right)^2
-
q(x)y^2(x)
\right]
+
f(x)y(x)
\right\}dx
\tag{20.45}
$$

取极值的必要条件。读者也可以直接验算。

**例 20.6**　写出偏微分方程定解问题

$$
\nabla^2u(\boldsymbol r)
+
k^2u(\boldsymbol r)
=
-\rho(\boldsymbol r),
\qquad
\boldsymbol r\in V,
\tag{20.46a}
$$

$$
\left.u(\boldsymbol r)\right|_\Sigma
=
f(\Sigma)
\tag{20.46b}
$$

的变分形式。

**解**　可以完全仿照例 20.5 的做法，考虑积分

$$
\iiint_V
\left[
\nabla^2u+k^2u+\rho(\boldsymbol r)
\right]
\delta u\,d\boldsymbol r.
$$

对于被积函数中的后两项，有

$$
\iiint_V
k^2u\delta u\,d\boldsymbol r
=
\frac12\delta
\iiint_V
k^2u^2\,d\boldsymbol r,
\qquad
\iiint_V
\rho(\boldsymbol r)\delta u\,d\boldsymbol r
=
\delta
\iiint_V
\rho(\boldsymbol r)u\,d\boldsymbol r.
$$

对于被积函数中的第一项，则需要应用 Green 第一公式以及边界条件 $\left.\delta u(\boldsymbol r)\right|_\Sigma=0$，有

$$
\iiint_V
\nabla^2u\,\delta u\,d\boldsymbol r
=
\iint_\Sigma
\delta u\,\nabla u\cdot d\boldsymbol\Sigma
-
\iiint_V
\nabla u\cdot\nabla(\delta u)\,d\boldsymbol r
=
-
\frac12\delta
\iiint_V
(\nabla u)^2\,d\boldsymbol r.
$$

因此，原方程就转化为

$$
\delta
\iiint_V
\left\{
\frac12
\left[
(\nabla u)^2-k^2u^2
\right]
-
\rho u
\right\}
d\boldsymbol r
=
0.
$$

这说明，定解问题 (20.46) 就等价于在边界条件 (20.46b) 下求泛函

$$
\iiint_V
\left\{
\frac12
\left[
(\nabla u)^2-k^2u^2
\right]
-
\rho u
\right\}
d\boldsymbol r
\tag{20.47}
$$

的极值问题。

**例 20.7**　写出偏微分方程的本征值问题

$$
\nabla^2u(\boldsymbol r)
+
\lambda u(\boldsymbol r)
=
0,
\qquad
\boldsymbol r\in V,
\tag{20.48a}
$$

$$
\left.u(\boldsymbol r)\right|_\Sigma
=
0
\tag{20.48b}
$$

的变分形式。

**解**　可将本问题看成例 20.6 的特殊情形。因此，本征值问题 (20.48) 就等价于泛函

$$
J[u]
=
\iiint_V
\left\{
[\nabla u(\boldsymbol r)]^2
-
\lambda[u(\boldsymbol r)]^2
\right\}
d\boldsymbol r
\tag{20.49}
$$

在边界条件 (20.48b) 下的极值问题。更进一步，把本征值 $\lambda$ 看成 Lagrange 乘子，那么，这个泛函极值问题又等价于泛函

$$
J[u]
=
\iiint_V
[\nabla u(\boldsymbol r)]^2
\,d\boldsymbol r
\tag{20.50}
$$

在边界条件 (20.48b) 和约束条件（本征函数的归一化条件）

$$
J_1[u]
\equiv
\iiint_V
[u(\boldsymbol r)]^2
\,d\boldsymbol r
=
1
\tag{20.51}
$$

下的条件极值问题。

不难证明，这样得到的泛函的条件极值问题的确和本征值问题 (20.48) 同解。这些本征函数正好就是泛函的极值函数，而本征值正好是泛函的极值。由于泛函 $J[u]$ 的二级变分

$$
\delta^2J[u]
=
2
\iiint_V
[\nabla(\delta u(\boldsymbol r))]^2
\,d\boldsymbol r
\tag{20.52}
$$

恒为正，所以泛函的极值是极小值。这些极小值中的最小者，当然就是本征值问题 (20.48) 的最小本征值。

## *§20.5 变边值问题

在实际问题中还会遇到另一类泛函的极值问题，即极值函数在一端或两端的数值并未指定。可以看下面这两个例子。

**例 20.8**　找出连接一固定点 $A$ 到一铅直线 $L$ 的路径（见图 20.2），使质点在重力作用下以最少的时间由 $A$ 到 $L$。在这个问题中，起点 $A$ 的位置是给定的，

$$
y_A=y(x_A).
$$

但是，终点 $B$ 的位置并不完全确定，它只是限定在 $L$ 上变化，因此，只是 $x_B$ 给定，而 $y_B$ 不定。

![Figure 20.2](...)

图 20.2

![Figure 20.3](...)

图 20.3

**例 20.9**　如图 20.3 所示，求两条不相交曲线之间的最短路径。这里，两个端点的位置都是不完全确定的，需要决定端点的坐标，以及连接这两点的最短路径。

从数学上看，这类问题仍然可以归结为求泛函

$$
J[y]
=
\int_{x_0}^{x_1}
F(x,y,y')\,dx
\tag{20.53}
$$

的极值问题，只不过边界条件需要修改。仿照前面的讨论，可得到这个泛函取极值的必要条件

$$
\begin{aligned}
\delta J[y]
&=
\int_{x_0}^{x_1}
\delta F(x,y,y')\,dx\\
&=
\int_{x_0}^{x_1}
\left(
\frac{\partial F}{\partial y}\delta y
+
\frac{\partial F}{\partial y'}\delta y'
\right)dx\\
&=
\left.
\frac{\partial F}{\partial y'}\delta y
\right|_{x_0}^{x_1}
+
\int_{x_0}^{x_1}
\left(
\frac{\partial F}{\partial y}
-
\frac d{dx}\frac{\partial F}{\partial y'}
\right)
\delta y\,dx
=
0.
\end{aligned}
$$

设起点的位置 $(x_0,y_0)$ 给定，因而 $\left.\delta y\right|_{x_0}=0$；但终点的位置 $(x_1,y_1)$ 仅 $x_1$ 给定，而 $y_1$ 不定，所以 $\left.\delta y\right|_{x_1}$ 也不定。但是，只要极值函数存在（因而终点位置 $(x_1,y_1)$ 完全确定），则此极值函数一定也是同一个泛函在两端固定的边界条件

$$
y(x_0)=y_0,
\qquad
y(x_1)=y_1
$$

下的极值函数，所以一定也还是方程

$$
\frac{\partial F}{\partial y}
-
\frac d{dx}
\frac{\partial F}{\partial y'}
=
0
\tag{20.54}
$$

的解。换句话说，这个 Euler-Lagrange 方程也还是泛函在变边值条件下取极值的必要条件。但是，只有这个条件不够，还必须有

$$
\left.
\frac{\partial F}{\partial y'}
\right|_{x=x_1}
=
0
\tag{20.55}
$$

才能保证泛函的一级变分 $\delta J[y]=0$。

结论：泛函 (20.53) 在一端完全固定（$y(x_0)=y_0$），另一端 $x_1$ 给定，而 $y(x_1)$ 不定的条件下取极值的必要条件是

$$
\frac{\partial F}{\partial y}
-
\frac d{dx}
\frac{\partial F}{\partial y'}
=
0,
\tag{20.56a}
$$

$$
y(x_0)=y_0,
\qquad
\left.
\frac{\partial F}{\partial y'}
\right|_{x=x_1}
=
0.
\tag{20.56b}
$$

**练习 20.2**　在两端均为变边值的条件下，求泛函

$$
J[y]
=
\int_{x_0}^{x_1}
F(x,y,y')\,dx
$$

取极值的必要条件。

再推广到二元函数的情形。在自由边值条件（边界 $\Gamma$ 上的 $u(x,y)$ 值自由）下，泛函

$$
J[u]
=
\iint_S
F
\left(
x,y,u,
\frac{\partial u}{\partial x},
\frac{\partial u}{\partial y}
\right)
dxdy
\tag{20.57}
$$

取极值的必要条件，仍然是此泛函的一级变分为 0，

$$
\begin{aligned}
\delta J[u]
&=
\iint_S
\left(
\frac{\partial F}{\partial u}\delta u
+
\frac{\partial F}{\partial u_x}\delta u_x
+
\frac{\partial F}{\partial u_y}\delta u_y
\right)dxdy\\
&=
\iint_S
\left(
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\frac{\partial F}{\partial u_x}
-
\frac{\partial}{\partial y}
\frac{\partial F}{\partial u_y}
\right)
\delta u\,dxdy\\
&\quad+
\iint_S
\left[
\frac{\partial}{\partial x}
\left(
\frac{\partial F}{\partial u_x}\delta u
\right)
+
\frac{\partial}{\partial y}
\left(
\frac{\partial F}{\partial u_y}\delta u
\right)
\right]dxdy\\
&=
\iint_S
\left(
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\frac{\partial F}{\partial u_x}
-
\frac{\partial}{\partial y}
\frac{\partial F}{\partial u_y}
\right)
\delta u\,dxdy\\
&\quad+
\int_\Gamma
\left(
-
\frac{\partial F}{\partial u_y}\,dx
+
\frac{\partial F}{\partial u_x}\,dy
\right)
\delta u
=
0.
\end{aligned}
$$

所以，必要条件的微分形式是

$$
\frac{\partial F}{\partial u}
-
\frac{\partial}{\partial x}
\frac{\partial F}{\partial u_x}
-
\frac{\partial}{\partial y}
\frac{\partial F}{\partial u_y}
=
0,
\tag{20.58a}
$$

$$
\int_\Gamma
\left(
-
\frac{\partial F}{\partial u_y}\,dx
+
\frac{\partial F}{\partial u_x}\,dy
\right)
\delta u
=
0.
\tag{20.58b}
$$

或者把沿边界的第二型线积分改写成第一型线积分的形式，边界条件 (20.58b) 又可改写成

$$
\left[
\frac{\partial F}{\partial u_x}\cos(n,x)
+
\frac{\partial F}{\partial u_y}\cos(n,y)
\right]_\Gamma
=
0.
\tag{20.58b'}
$$

## §20.6 Rayleigh-Ritz 方法

到现在为止，读者应该已经了解了变分法应用于物理问题的大概轮廓。变分法在物理学中的应用，可以分为两个主要的方面。一种应用是作为基本物理规律的表述语言。可以用 Hamilton 原理或其他类似的语言描述力学系统（质点、质点组等）的运动，可以用费马原理描述光线在介质中的传播，也可以用变分的语言描述电磁场乃至微观粒子的运动，等等。在物理学的这些分支中，支配物质运动的各种特定形式的基本规律，无一例外地都可以表述为各自的泛函极值问题。变分法的这种应用，具有重要的理论意义。它可以使我们用统一的语言描述物质世界的运动，可以协调一致地处理涉及多个物理学分支的综合问题，也可以更方便地从已知的物理领域向新的领域扩展。变分法的第二种应用则体现出它的实用价值：它为求解具体的物理问题提供了一种新的灵活手段。在变分法的基础上，可以建立起实用的近似解法。本节就以常微分方程本征值问题为例介绍这种近似方法。

假设有一个一般的本征值问题

$$
\hat L X
=
\lambda X,
\tag{20.59a}
$$

$$
D(\hat L)
=
\left\{
X(x):
a\le x\le b,
\ \cos\alpha\,X(a)-\sin\alpha\,X'(a)=0,
\ \cos\beta\,X(b)+\sin\beta\,X'(b)=0
\right\}.
\tag{20.59b}
$$

这时存在两种可能：一种是常微分方程 (20.59a) 的解已知，很容易求得；另一种可能是还需要用常微分方程级数解法，才能求出常微分方程的解。但是，无论哪种情况，都还要代入边界条件，定出本征值和本征函数。一般说来，除了少数已经熟悉的函数外，难以指望能得到本征值的准确表达式。即使像

$$
\hat L X
\equiv
-
\frac{d^2X}{dx^2}
=
\lambda X,
$$

$$
D(\hat L)
=
\left\{
X(x):
a\le x\le b,
\ \cos\alpha\,X(a)-\sin\alpha\,X'(a)=0,
\ \cos\beta\,X(b)+\sin\beta\,X'(b)=0
\right\}
$$

这样简单的方程，有熟知的两个线性无关解 $\sin\sqrt\lambda x$ 和 $\cos\sqrt\lambda x$，在一般的第三类边界条件下，也无法写出本征值的显明表达式。对于一般的本征值问题，困难可想而知。变分法就为我们提供了求解本征值的近似方法。

用 Rayleigh-Ritz 方法近似求解本征值问题的基本思路是：首先把本征值问题转化为泛函的条件极值问题，然后在一定的函数空间中求解，把问题又转化为函数的条件极值问题。前提条件是选择的函数空间（对于此本征值问题）是完备的。从实用的角度看，要选择一个“好”的函数空间（实际上是一个函数序列），一方面便于计算，一方面又能够足够快地、足够精确地求得本征值的近似值。这就要求函数序列具有本征函数所要求的主要基本特征，要求我们事先从物理上和数学上对于本征函数的性质做出尽可能准确的判断。

为了便于比较，不妨举一个已知精确解的例子。

**例 20.10**　求本征值问题

$$
\frac1x
\frac d{dx}
\left(
x\frac{dy}{dx}
\right)
+
\lambda y(x)
=
0,
\tag{20.60a}
$$

$$
y(0)\text{ 有界},
\qquad
y(1)=0
\tag{20.60b}
$$

的最小本征值。

**解**　这个本征值问题在 §20.3 的例 20.4 中已经讨论过。当时讨论的是泛函

$$
I[y]
=
\int_0^1xy'^2\,dx
\tag{20.61}
$$

在边界条件 (20.60b) 和约束条件

$$
I_1[y]
\equiv
\int_0^1xy^2\,dx
=
1
\tag{20.62}
$$

下的条件极值问题，它的 Euler-Lagrange 方程就是 (20.60a) 式。

现在就用 Rayleigh-Ritz 方法来近似求解这个泛函的条件极值问题。事先，我们对于本征函数的了解是，它除了必须满足边界条件 (20.60b) 之外，还可以具有奇偶性（为什么？请读者证明）。因此，可用多项式序列

$$
y_n(x)
=
\sum_{k=1}^{n}
\alpha_k(1-x^2)^k,
\qquad
n=1,2,3,\cdots
\tag{20.63}
$$

去逼近本征函数。首先取近似的本征函数 $y_2(x)$，即在 (20.63) 中取前两项，代入泛函 (20.61) 及约束条件 (20.62)，得

$$
I[y_2]
=
\int_0^1xy_2'^2\,dx
=
\alpha_1^2
+
\frac43\alpha_1\alpha_2
+
\frac23\alpha_2^2,
\tag{20.64}
$$

$$
I_1[y_2]
=
\int_0^1xy_2^2\,dx
=
\frac16\alpha_1^2
+
\frac14\alpha_1\alpha_2
+
\frac1{10}\alpha_2^2
=
1.
\tag{20.65}
$$

这可以看成是 $\alpha_1$ 和 $\alpha_2$ 的二元函数的条件极值问题，必要条件是

$$
\frac{\partial(I-\lambda I_1)}{\partial\alpha_1}
=
2\alpha_1
+
\frac43\alpha_2
-
\lambda
\left(
\frac13\alpha_1
+
\frac14\alpha_2
\right)
=
0,
\tag{20.66}
$$

$$
\frac{\partial(I-\lambda I_1)}{\partial\alpha_2}
=
\frac43\alpha_2
+
\frac43\alpha_1
-
\lambda
\left(
\frac15\alpha_2
+
\frac14\alpha_1
\right)
=
0.
\tag{20.67}
$$

这又是关于 $\alpha_1$ 和 $\alpha_2$ 的代数方程组，有非零解的充分必要条件是

$$
\begin{vmatrix}
2-\dfrac\lambda3
&
\dfrac43-\dfrac\lambda4
\\[8pt]
\dfrac43-\dfrac\lambda4
&
\dfrac43-\dfrac\lambda5
\end{vmatrix}
=
0,
\qquad
\text{即}
\qquad
3\lambda^2-128\lambda+640=0.
\tag{20.68}
$$

解之得

$$
\lambda
=
\frac{64}{3}
\pm
\frac83\sqrt{34}.
\tag{20.69}
$$

这两个给出的都是 $\lambda$ 的极小值。在 §20.3 和 §20.4 中已经论证过，最小的极小值就对应于最小的本征值。这里得到的当然只是本征值问题 (20.60) 的最小本征值的近似值

$$
\bar\lambda_1
=
\frac{64}{3}
-
\frac83\sqrt{34}
=
5.7841\cdots,
\tag{20.70}
$$

它和精确值

$$
\lambda_1
=
(2.4048\cdots)^2
=
5.7831\cdots
$$

的相对误差不到 $2\times10^{-4}$。相应地，本征函数的近似解是

$$
\bar y_1(x)
=
\alpha_1(1-x^2)
+
\alpha_2(1-x^2)^2,
\tag{20.71a}
$$

$$
\alpha_1
=
2\sqrt{12-33\sqrt{2/17}}
=
1.6505676\cdots,
\tag{20.71b}
$$

$$
\alpha_2
=
\sqrt{80-230\sqrt{2/17}}
=
1.0538742\cdots.
\tag{20.71c}
$$

为了与精确解

$$
y_1(x)
=
\frac{\sqrt2}{J_1(\mu_1)}
J_0(\mu_1x)
$$

做比较，不妨计算

$$
\begin{aligned}
\Delta
&=
\int_0^1
[y_1(x)-\bar y_1(x)]^2x\,dx\\
&=
2
-
2\int_0^1
y_1(x)\bar y_1(x)x\,dx\\
&=
2
\left\{
1
-
\left[
\alpha_1\frac{4\sqrt2}{\mu_1^3}
+
\alpha_2\frac{8\sqrt2}{\mu_1^3}
\left(
\frac8{\mu_1^2}-1
\right)
\right]
\right\}
=
1.66\times10^{-5}.
\end{aligned}
$$

在多项式逼近 (20.63) 中才只取了两项，本征值和本征函数就能达到这个精度，这的确是令人惊异的。可以想象，如果取的项数更多，得到的精度会更高。

从上面的计算可以看出，在应用 Rayleigh-Ritz 方法时，只能求得最低的几个本征值的近似值，本征值的个数和使用的逼近函数中的参数数目相同。这是应用 Rayleigh-Ritz 方法求解本征值问题的一个特点。在实际应用中，并不会因为 Rayleigh-Ritz 方法只能求得有限个本征值而降低它的实用价值，因为有不少问题只需要求出最小的若干个本征值。

再看一下上面得到的第二个本征值

$$
\bar\lambda_2
=
\frac{64}{3}
+
\frac83\sqrt{34}
=
36.883\cdots,
\tag{20.72}
$$

它和精确值 $\lambda_2=30.471\cdots$ 之间的误差竟超过 $20\%$！为了求得足够精确的第二个本征值，就必须

![Figure 20.4](...)

图 20.4　本征函数 $y_1(x)$ 与近似解 $\bar y_1(x)$ 之差

增加逼近函数中的参数，这必然以计算量的急剧增长为代价。在实用中，更好的办法是求解一个新的泛函条件极值问题，它和原来的泛函条件极值问题的差别只在于排除掉第一个本征值。这只要在原来的泛函条件极值问题中再附加上一个正交条件

$$
\int_0^1
y(x)\bar y_1(x)x\,dx
=
0
\tag{20.73}
$$

即可。这样，在这个新的泛函条件极值问题中，最小的本征值当然就是原来的第二个本征值了。读者不难想到，如果需要求得更高的本征值，应该如何处理。

## 习题

1. 写出使下列泛函取极值的 Euler-Lagrange 方程，并求解：

(1)

$$
\int_{x_0}^{x_1}
\sqrt{1+y^2y'^2}\,dx;
$$

(2)

$$
\int_{x_0}^{x_1}
(y^2+y'^2)\,dx;
$$

(3)

$$
\int_{x_0}^{x_1}
\frac{x}{x+y'}\,dx;
$$

(4)

$$
\int_{x_0}^{x_1}
\sqrt{1+x}\sqrt{1+y'^2}\,dx.
$$

规定极值曲线均通过平面上的已知点 $(x_0,y_0)$ 和 $(x_1,y_1)$。

2. 求锥面 $x^2+y^2=z^2$ 上的“短程线”（准确说，称为测地线，geodesic）。

3. 求圆柱面上的测地线，设圆柱的母线平行于 $z$ 轴。

4. 求泛函

$$
J[y]
=
\iiint_{x^2+y^2+z^2<a^2}
\left[
(\nabla u)^2-k^2u^2
\right]
dxdydz
$$

在边界条件

$$
\left.u\right|_{x^2+y^2+z^2=a^2}
=
z
$$

下的极值函数，其中 $k>0$ 和 $a>0$ 都是已知常数，且满足 $ka\ne\tan ka$。

5. 光在折射率为 $n$ 的介质中的传播速率为 $v=ds/dt=c/n$，$c$ 是真空中的速率，于是光由 $A$ 点 $(x_0,y_0)$ 传播到 $B$ 点 $(x_1,y_1)$ 的时间便是

$$
T
=
\int_{(x_0,y_0)}^{(x_1,y_1)}
\frac{ds}{v}
=
\frac1c
\int_{(x_0,y_0)}^{(x_1,y_1)}
nds.
$$

而光由 $A$ 到 $B$ 的实际路径应当使 $T$ 取极值（Fermat 原理）。试求光在下列介质中传播时的实际轨迹：

(1) $n=k(x+1)$；

(2) $n=k\sqrt y$；

(3) $n=\dfrac{k}{2x+3}$；

(4) $n=\dfrac{k}{y}$；

(5) $n=ke^y$；

(6) $n=k\sqrt{x+y}$；

(7) $n=kr^{-1/2}$；

(8) $n=\dfrac{k}{r}$。

其中 $k$ 均为已知常数，$r^2=x^2+y^2$。

6. 求泛函

$$
J[y]
=
\int_0^1
\left[
\left(
x\frac{dy}{dx}
\right)^2
+
2y^2
\right]dx
$$

在边界条件

$$
y(0)=0,
\qquad
y(1)=0
$$

以及约束条件

$$
J_1[y]
=
\int_0^1
y^2(x)x^2\,dx
=
1
$$

下的最小值。

7. 试写出本征值问题

$$
\nabla^2u+\lambda u=0,
$$

$$
\left(
\alpha u
+
\beta\frac{\partial u}{\partial n}
\right)_\Sigma
=
0
$$

所对应的泛函极值问题，设 $\beta\ne0$。

8. 用 Rayleigh-Ritz 方法求出

$$
y''+\lambda y=0,
$$

$$
y(-1)=0,
\qquad
y(1)=0
$$

的最低的两个本征值的近似值，取试探函数为

(1) $y=c_1(1-x^2)+c_2x(1-x^2)$；

(2) $y=c_1(1-x^2)+c_2x^2(1-x^2)$。
