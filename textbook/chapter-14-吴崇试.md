# 第十四章 正交曲面坐标系

上一章中讨论的问题所涉及的空间几何形状，除了一维的直线（线段）外，仅限于二维平面上的矩形区域，以及三维空间中的长方体。对于这些几何形状，总可以适当地放置直角坐标架，使得所讨论区域的边界面都与坐标面重合，从而实现齐次边界条件的分离变量。

如果我们讨论的问题涉及的空间区域具有其他几何形状，例如圆柱形（包括它的特殊情形，二维平面上的圆形区域）或球形，这时无论怎样放置直角坐标架，总不能使得区域的边界面全部都和坐标面重合。因此，即使边界条件是齐次的，也无法分离变量。

解决这个问题的办法是选用别的坐标系。常用的是一些正交曲面坐标系。这一章我们将介绍正交曲面坐标系以及矢量分析的基本知识。

## §14.1 正交曲面坐标系

为方便讨论，我们常常根据所要讨论的空间区域的特点选择不同的坐标系。例如圆形区域，首选的坐标系当然是平面极坐标系。圆柱形区域，当然应当考虑柱坐标系。对于球形区域，就应该选用球坐标系。作为这些坐标系的概括，利用我们熟悉的直角坐标系 $(x,y,z)$，可以定义曲面坐标系 $(x^1,x^2,x^3)$①，

$$
x^1=\xi(x,y,z),\qquad
x^2=\eta(x,y,z),\qquad
x^3=\zeta(x,y,z).
\tag{14.1}
$$

显然 (14.1) 式就是直角坐标系 $(x,y,z)$ 到曲面坐标系 $(x^1,x^2,x^3)$ 的坐标变换②。曲面坐标系 $(x^1,x^2,x^3)$ 的坐标面是三组曲面

$$
x^1=\text{常数},\qquad
x^2=\text{常数},\qquad
x^3=\text{常数}.
$$

空间任意一点的坐标 $(x^1,x^2,x^3)$，就由过该点的三个坐标面决定。

对于空间的任意一点，如果通过该点的三个坐标面总是互相垂直的，这个坐标系就称为正交曲面坐标系。例如，在直角坐标系中，过空间任意一点 $(x_0,y_0,z_0)$ 的三个坐标面

$$
x=x_0,\qquad
y=y_0,\qquad
z=z_0
$$

就是互相垂直的。

① 这里的 $x^i\ (i=1,2,3)$ 中，上标 $i$ 用来标记空间点的坐标（分量），并不表示方次。

② 为了保证 $x^1,x^2$ 和 $x^3$ 是相互独立的，应当要求 Jacobi 行列式

$$
\frac{\partial(x^1,x^2,x^3)}{\partial(x,y,z)}
\equiv
\begin{vmatrix}
\dfrac{\partial x^1}{\partial x} &
\dfrac{\partial x^1}{\partial y} &
\dfrac{\partial x^1}{\partial z}\\[6pt]
\dfrac{\partial x^2}{\partial x} &
\dfrac{\partial x^2}{\partial y} &
\dfrac{\partial x^2}{\partial z}\\[6pt]
\dfrac{\partial x^3}{\partial x} &
\dfrac{\partial x^3}{\partial y} &
\dfrac{\partial x^3}{\partial z}
\end{vmatrix}
\ne0.
$$

判断一个坐标系是不是正交曲面坐标系，当然可以直接从坐标系的定义（即坐标面的方程）去判断。更常用的办法①是由 (14.1) 式计算出弧长②

$$
\begin{aligned}
ds^2
&=dx^2+dy^2+dz^2\\
&=
\left(
\frac{\partial x}{\partial x^1}dx^1
+
\frac{\partial x}{\partial x^2}dx^2
+
\frac{\partial x}{\partial x^3}dx^3
\right)^2\\
&\quad+
\left(
\frac{\partial y}{\partial x^1}dx^1
+
\frac{\partial y}{\partial x^2}dx^2
+
\frac{\partial y}{\partial x^3}dx^3
\right)^2\\
&\quad+
\left(
\frac{\partial z}{\partial x^1}dx^1
+
\frac{\partial z}{\partial x^2}dx^2
+
\frac{\partial z}{\partial x^3}dx^3
\right)^2\\
&=
\sum_{i,j=1,2,3}g_{ij}\,dx^i dx^j,
\end{aligned}
\tag{14.2}
$$

其中

$$
g_{ij}=g_{ji}
=
\frac{\partial x}{\partial x^i}
\frac{\partial x}{\partial x^j}
+
\frac{\partial y}{\partial x^i}
\frac{\partial y}{\partial x^j}
+
\frac{\partial z}{\partial x^i}
\frac{\partial z}{\partial x^j}.
$$

$g_{ij}$ 构成的矩阵

$$
G=
\begin{pmatrix}
g_{11}&g_{12}&g_{13}\\
g_{21}&g_{22}&g_{23}\\
g_{31}&g_{32}&g_{33}
\end{pmatrix}
$$

称为此空间的度规（metric）。如果度规是对角矩阵，

$$
g_{ij}=g_{ii}\delta_{ij},
\tag{14.3}
$$

则此坐标系为正交曲面坐标系。

**例 14.1**　直角坐标系下，

$$
ds^2=dx^2+dy^2+dz^2,
$$

所以

$$
g_{ij}=\delta_{ij},
\qquad
i,j=1,2,3.
$$

因此直角坐标系是正交曲面坐标系。

**例 14.2**　对于柱坐标系，$x=\rho\cos\phi,\ y=\rho\sin\phi,\ z=z$，

$$
\begin{aligned}
ds^2
&=
(\cos\phi\,d\rho-\rho\sin\phi\,d\phi)^2
+
(\sin\phi\,d\rho+\rho\cos\phi\,d\phi)^2
+
dz^2\\
&=
d\rho^2+\rho^2d\phi^2+dz^2.
\end{aligned}
\tag{14.4}
$$

所以，$g_{11}=1,\ g_{22}=\rho^2,\ g_{33}=1$，当 $i\ne j$ 时 $g_{ij}=0$。因此柱坐标系是正交曲面坐标系。

① 这种讨论方法的一个优点是可以直接推广到高维空间的情形。

② 在微分几何中，更常略去 (14.2) 式中的和号，而直接写成

$$
ds^2=g_{ij}dx^i dx^j.
$$

按照 Einstein 规则，此式应理解为需对所有重复指标（并且一个是上指标，一个是下指标）求和。

**例 14.3**　对于球坐标系，$x=r\sin\theta\cos\phi,\ y=r\sin\theta\sin\phi,\ z=r\cos\theta$，

$$
\begin{aligned}
ds^2
&=
(\sin\theta\cos\phi\,dr+r\cos\theta\cos\phi\,d\theta-r\sin\theta\sin\phi\,d\phi)^2\\
&\quad+
(\sin\theta\sin\phi\,dr+r\cos\theta\sin\phi\,d\theta+r\sin\theta\cos\phi\,d\phi)^2\\
&\quad+
(\cos\theta\,dr-r\sin\theta\,d\theta)^2\\
&=
dr^2+r^2d\theta^2+r^2\sin^2\theta\,d\phi^2.
\end{aligned}
\tag{14.5}
$$

所以，$g_{11}=1,\ g_{22}=r^2,\ g_{33}=r^2\sin^2\theta$，$i\ne j$ 时 $g_{ij}=0$。故球坐标系也是正交曲面坐标系。

## *§14.2 正交曲面坐标系中的 Laplace 算符

这一节，通过外微分法介绍正交曲面坐标系中 Laplace 算符的一般形式。这种方法的优点在于它的协变性，即可以脱离开坐标系的具体定义，而得到最普遍的表达式。作为一个最初步的介绍，我们略去数学上的严格定义，只给出有关运算的规则。

**外微分法则**　这里要介绍外微分算符、楔积运算及 $*$ 算符，以及微分形式的概念。

首先定义外微分算符 $d$。它作用在（$n$ 维空间的标量）函数 $f(x^1,x^2,\cdots,x^n)$ 上，

$$
d:f\mapsto df
=
\sum_{i=1}^n
\frac{\partial f}{\partial x^i}dx^i,
\tag{14.6}
$$

得到的 $df$ 称为一次微分形式（简称一次形式）。

设 $u$ 是三维空间的标量函数，例如三维空间的温度分布。对于直角坐标系，

$$
du
=
\frac{\partial u}{\partial x}dx
+
\frac{\partial u}{\partial y}dy
+
\frac{\partial u}{\partial z}dz.
\tag{14.7}
$$

对于柱坐标系，

$$
du
=
\frac{\partial u}{\partial\rho}d\rho
+
\frac{\partial u}{\partial\phi}d\phi
+
\frac{\partial u}{\partial z}dz.
\tag{14.8}
$$

对于球坐标系，

$$
du
=
\frac{\partial u}{\partial r}dr
+
\frac{\partial u}{\partial\theta}d\theta
+
\frac{\partial u}{\partial\phi}d\phi.
\tag{14.9}
$$

外微分算符 $d$ 作用于一次微分形式 $\alpha=\sum_{i=1}^n\alpha_i dx^i$，其中 $\alpha_i$ 都是标量函数，则得到二次微分形式

$$
d\alpha
=
d\left(
\sum_{i=1}^n\alpha_i dx^i
\right)
=
\sum_{i=1}^n d\alpha_i\wedge dx^i
=
\sum_{i=1}^n
\sum_{j=1}^n
\frac{\partial\alpha_i}{\partial x^j}
dx^j\wedge dx^i.
$$

规定

$$
d(dx^i)=0,
\qquad
i=1,2,\cdots,n.
\tag{14.10}
$$

运算 $\wedge$ 称为楔积。

**运算法则 1**

$$
dx^i\wedge dx^j
=
-dx^j\wedge dx^i,
\qquad
i,j=1,2,\cdots,n.
\tag{14.11}
$$

因此，

$$
dx^i\wedge dx^i=0,
\qquad
i=1,2,\cdots,n.
$$

所以，$n$ 维空间最多只有 $n$ 次微分形式。将外微分算符 $d$ 作用于 $p$ 次 $(p<n)$ 微分形式 $\alpha=\sum_I\alpha_I dx^I$，就得到 $(p+1)$ 次微分形式：

$$
d\alpha
=
d\left(
\sum_I\alpha_I dx^I
\right)
=
\sum_I d\alpha_I\wedge dx^I
=
\sum_I\sum_{i=1}^n
\frac{\partial\alpha_I}{\partial x^i}
dx^i\wedge dx^I,
\tag{14.12}
$$

其中

$$
dx^I
\equiv
dx^{i_1}\wedge dx^{i_2}\wedge\cdots\wedge dx^{i_p},
\qquad
I=(i_1,i_2,\cdots,i_p).
\tag{14.13}
$$

$i_1,i_2,\cdots,i_p$ 都是正整数 $1,2,\cdots,n$。

**运算法则 2**　设 $\alpha$ 为 $p$ 次微分形式，$\beta$ 和 $\gamma$ 为 $q$ 次微分形式，

$$
d(\beta+\gamma)=d\beta+d\gamma,
\tag{14.14}
$$

$$
d(\alpha\wedge\beta)
=
(d\alpha)\wedge\beta
+
(-)^p\alpha\wedge(d\beta),
\tag{14.15}
$$

$$
d(d\alpha)=0.
\tag{14.16}
$$

$*$ 运算是一个线性变换，它把 $p$ 次微分形式变换为相应的 $n-p$ 次微分形式，

$$
*
\left(
dx^{i_1}\wedge dx^{i_2}\wedge\cdots\wedge dx^{i_p}
\right)
=
\frac{\sqrt{\det G}}
{g_{i_1i_1}g_{i_2i_2}\cdots g_{i_pi_p}}
dx^{i_{p+1}}\wedge dx^{i_{p+2}}\wedge\cdots\wedge dx^{i_n},
\tag{14.17}
$$

其中 $(i_1,i_2,\cdots,i_p,i_{p+1},i_{p+2},\cdots,i_n)$ 构成 $(1,2,3,\cdots,n)$ 的偶排列。特别是，对于二维空间，有

$$
*dx^1
=
\frac{\sqrt{\det G}}{g_{11}}dx^2,
\qquad
*dx^2
=
-\frac{\sqrt{\det G}}{g_{22}}dx^1
=
-\frac{g_{11}}{\sqrt{\det G}}dx^1,
\tag{14.18}
$$

而对于三维空间，则有

$$
*dx^i
=
\frac{\sqrt{\det G}}{g_{ii}}dx^I,
\qquad
*dx^I
=
\frac{g_{ii}}{\sqrt{\det G}}dx^i,
\qquad
i=1,2,3,
\tag{14.19}
$$

其中 $(i,I)$ 构成 $(1,2,3)$ 的偶排列，$\det G$ 为度规矩阵 $G$ 的行列式值。

**运算法则 3**　三维空间中，有

$$
*1
=
\sqrt{\det G}\,
dx^1\wedge dx^2\wedge dx^3,
\tag{14.20}
$$

$$
*
\left(
\sqrt{\det G}\,
dx^1\wedge dx^2\wedge dx^3
\right)
=
1.
\tag{14.21}
$$

注意 $\sqrt{\det G}\,dx^1\wedge dx^2\wedge dx^3$ 正好是通常的三维空间的体积元。

**例 14.4**　对于直角坐标系，$\det G=1$，

$$
*du
=
\frac{\partial u}{\partial x}dy\wedge dz
+
\frac{\partial u}{\partial y}dz\wedge dx
+
\frac{\partial u}{\partial z}dx\wedge dy.
\tag{14.22}
$$

**例 14.5**　对于柱坐标系，$\det G=\rho^2$，

$$
*du
=
\rho\frac{\partial u}{\partial\rho}d\phi\wedge dz
+
\frac1\rho
\frac{\partial u}{\partial\phi}dz\wedge d\rho
+
\rho\frac{\partial u}{\partial z}d\rho\wedge d\phi.
\tag{14.23}
$$

**例 14.6**　对于球坐标系，$\det G=r^4\sin^2\theta$，

$$
*du
=
r^2\sin\theta
\frac{\partial u}{\partial r}
d\theta\wedge d\phi
+
\sin\theta
\frac{\partial u}{\partial\theta}
d\phi\wedge dr
+
\frac1{\sin\theta}
\frac{\partial u}{\partial\phi}
dr\wedge d\theta.
\tag{14.24}
$$

**外微分运算与矢量分析**　外微分运算的优点是它的协变性，即它的表述形式与坐标系无关。运用外微分运算可以方便地导出各种正交曲面坐标系中的矢量分析公式，只需注意把握微分形式中的 $\sqrt{g_{ii}}\,dx^i$ 和矢量式中坐标单位矢量 $\boldsymbol e_i$ 之间的对应关系。

外微分 $d$（作用于标量函数）是梯度 $\operatorname{grad}$ 的协变微分形式。准确说，将一次微分形式 (14.6) 式改写为

$$
du
=
\sum_i
\frac1{\sqrt{g_{ii}}}
\frac{\partial u}{\partial x^i}
\sqrt{g_{ii}}\,dx^i,
$$

就正好对应于梯度

$$
\nabla u
\equiv
\operatorname{grad}u
=
\sum_i
\frac1{\sqrt{g_{ii}}}
\frac{\partial u}{\partial x^i}
\boldsymbol e_i.
\tag{14.25}
$$

**例 14.7**　直角坐标系中，$\boldsymbol e_x$ 对应 $dx$，$\boldsymbol e_y$ 对应 $dy$，$\boldsymbol e_z$ 对应 $dz$，梯度 $\nabla u$ 的三个分量是

$$
(\nabla u)_x
=
\frac{\partial u}{\partial x},
\qquad
(\nabla u)_y
=
\frac{\partial u}{\partial y},
\qquad
(\nabla u)_z
=
\frac{\partial u}{\partial z}.
\tag{14.26}
$$

**例 14.8**　柱坐标系中，$\boldsymbol e_\rho$ 对应 $d\rho$，$\boldsymbol e_\phi$ 对应 $\rho d\phi$，$\boldsymbol e_z$ 对应 $dz$，梯度 $\nabla u$ 的三个分量是

$$
(\nabla u)_\rho
=
\frac{\partial u}{\partial\rho},
\qquad
(\nabla u)_\phi
=
\frac1\rho
\frac{\partial u}{\partial\phi},
\qquad
(\nabla u)_z
=
\frac{\partial u}{\partial z}.
\tag{14.27}
$$

**例 14.9**　球坐标系中，$\boldsymbol e_r$ 对应 $dr$，$\boldsymbol e_\theta$ 对应 $r\,d\theta$，$\boldsymbol e_\phi$ 对应 $r\sin\theta\,d\phi$，梯度 $\nabla u$ 的三个分量是

$$
(\nabla u)_r
=
\frac{\partial u}{\partial r},
\qquad
(\nabla u)_\theta
=
\frac1r
\frac{\partial u}{\partial\theta},
\qquad
(\nabla u)_\phi
=
\frac1{r\sin\theta}
\frac{\partial u}{\partial\phi}.
\tag{14.28}
$$

$*d$ 是旋度 $\operatorname{curl}$ 的协变微分形式。

**例 14.10**　在直角坐标系中，矢量

$$
\boldsymbol A
=
A_x(x,y,z)\boldsymbol e_x
+
A_y(x,y,z)\boldsymbol e_y
+
A_z(x,y,z)\boldsymbol e_z
\tag{14.29}
$$

所对应的一次微分形式是

$$
A_x(x,y,z)dx
+
A_y(x,y,z)dy
+
A_z(x,y,z)dz.
\tag{14.30}
$$

用 $*d$ 作用，就得到

$$
*d(A_xdx+A_ydy+A_zdz)
\tag{14.31}
$$

$$
\begin{aligned}
&=
*
\left(
\frac{\partial A_x}{\partial y}dy\wedge dx
+
\frac{\partial A_x}{\partial z}dz\wedge dx
+
\frac{\partial A_y}{\partial x}dx\wedge dy
+
\frac{\partial A_y}{\partial z}dz\wedge dy\right.\\
&\qquad\left.
+
\frac{\partial A_z}{\partial x}dx\wedge dz
+
\frac{\partial A_z}{\partial y}dy\wedge dz
\right)\\
&=
\left(
\frac{\partial A_z}{\partial y}
-
\frac{\partial A_y}{\partial z}
\right)dx
+
\left(
\frac{\partial A_x}{\partial z}
-
\frac{\partial A_z}{\partial x}
\right)dy
+
\left(
\frac{\partial A_y}{\partial x}
-
\frac{\partial A_x}{\partial y}
\right)dz.
\end{aligned}
\tag{14.32}
$$

相应地，就有

$$
\operatorname{curl}\boldsymbol A
\equiv
\nabla\times\boldsymbol A
=
\left(
\frac{\partial A_z}{\partial y}
-
\frac{\partial A_y}{\partial z}
\right)\boldsymbol e_x
+
\left(
\frac{\partial A_x}{\partial z}
-
\frac{\partial A_z}{\partial x}
\right)\boldsymbol e_y
+
\left(
\frac{\partial A_y}{\partial x}
-
\frac{\partial A_x}{\partial y}
\right)\boldsymbol e_z.
\tag{14.33}
$$

**例 14.11**　在柱坐标系中，矢量

$$
\boldsymbol A
=
A_\rho(\rho,\phi,z)\boldsymbol e_\rho
+
A_\phi(\rho,\phi,z)\boldsymbol e_\phi
+
A_z(\rho,\phi,z)\boldsymbol e_z
\tag{14.34}
$$

所对应的一次微分形式是

$$
A_\rho(\rho,\phi,z)d\rho
+
A_\phi(\rho,\phi,z)\rho d\phi
+
A_z(\rho,\phi,z)dz.
\tag{14.35}
$$

作用以 $*d$，就得到

$$
\begin{aligned}
&*d(A_\rho d\rho+A_\phi\rho d\phi+A_zdz)\\
&=
\left[
\frac{\partial A_z}{\partial\phi}
-
\frac{\partial(\rho A_\phi)}{\partial z}
\right]\frac1\rho d\rho
+
\left(
\frac{\partial A_\rho}{\partial z}
-
\frac{\partial A_z}{\partial\rho}
\right)\rho d\phi\\
&\quad+
\left[
\frac{\partial(\rho A_\phi)}{\partial\rho}
-
\frac{\partial A_\rho}{\partial\phi}
\right]\frac1\rho dz.
\end{aligned}
\tag{14.36}
$$

相应地，就有

$$
\operatorname{curl}\boldsymbol A
\equiv
\nabla\times\boldsymbol A
=
\frac1\rho
\left[
\frac{\partial A_z}{\partial\phi}
-
\frac{\partial(\rho A_\phi)}{\partial z}
\right]\boldsymbol e_\rho
+
\left(
\frac{\partial A_\rho}{\partial z}
-
\frac{\partial A_z}{\partial\rho}
\right)\boldsymbol e_\phi
+
\frac1\rho
\left[
\frac{\partial(\rho A_\phi)}{\partial\rho}
-
\frac{\partial A_\rho}{\partial\phi}
\right]\boldsymbol e_z.
\tag{14.37}
$$

**例 14.12**　在球坐标系中，矢量

$$
\boldsymbol A
=
A_r(r,\theta,\phi)\boldsymbol e_r
+
A_\theta(r,\theta,\phi)\boldsymbol e_\theta
+
A_\phi(r,\theta,\phi)\boldsymbol e_\phi,
\tag{14.38}
$$

所对应的一次微分形式是

$$
A_r(r,\theta,\phi)dr
+
A_\theta(r,\theta,\phi)r\,d\theta
+
A_\phi(r,\theta,\phi)r\sin\theta\,d\phi.
\tag{14.39}
$$

作用以 $*d$，就得到

$$
\begin{aligned}
&*d(A_rdr+A_\theta r\,d\theta+A_\phi r\sin\theta\,d\phi)\\
&=
\left[
\frac{\partial(r\sin\theta A_\phi)}{\partial\theta}
-
\frac{\partial(rA_\theta)}{\partial\phi}
\right]
\frac1{r^2\sin\theta}dr\\
&\quad+
\left[
\frac{\partial A_r}{\partial\phi}
-
\frac{\partial(r\sin\theta A_\phi)}{\partial r}
\right]
\frac1{\sin\theta}d\theta\\
&\quad+
\left[
\frac{\partial(rA_\theta)}{\partial r}
-
\frac{\partial A_r}{\partial\theta}
\right]
\sin\theta\,d\phi.
\end{aligned}
\tag{14.40}
$$

相应地，就有

$$
\begin{aligned}
\operatorname{curl}\boldsymbol A
\equiv
\nabla\times\boldsymbol A
&=
\frac1{r\sin\theta}
\left[
\frac{\partial(\sin\theta A_\phi)}{\partial\theta}
-
\frac{\partial A_\theta}{\partial\phi}
\right]\boldsymbol e_r\\
&\quad+
\frac1r
\left[
\frac1{\sin\theta}
\frac{\partial A_r}{\partial\phi}
-
\frac{\partial(rA_\phi)}{\partial r}
\right]\boldsymbol e_\theta\\
&\quad+
\frac1r
\left[
\frac{\partial(rA_\theta)}{\partial r}
-
\frac{\partial A_r}{\partial\theta}
\right]\boldsymbol e_\phi.
\end{aligned}
\tag{14.41}
$$

$*d*$ 是散度 $\operatorname{div}$ 的协变微分形式。

**例 14.13**　在直角坐标系中，利用

$$
*dx=dy\wedge dz,\qquad
*(dy\wedge dz)=dx,\qquad
*(dx\wedge dy\wedge dz)=1,
$$

$$
*dy=dz\wedge dx,\qquad
*(dz\wedge dx)=dy,\qquad
*1=dx\wedge dy\wedge dz,
$$

$$
*dz=dx\wedge dy,\qquad
*(dx\wedge dy)=dz,
$$

得

$$
\begin{aligned}
&*d*(A_xdx+A_ydy+A_zdz)\\
&=
*d(A_xdy\wedge dz+A_ydz\wedge dx+A_zdx\wedge dy)\\
&=
*
\left(
\frac{\partial A_x}{\partial x}dx\wedge dy\wedge dz
+
\frac{\partial A_y}{\partial y}dy\wedge dz\wedge dx
+
\frac{\partial A_z}{\partial z}dz\wedge dx\wedge dy
\right)\\
&=
\frac{\partial A_x}{\partial x}
+
\frac{\partial A_y}{\partial y}
+
\frac{\partial A_z}{\partial z}.
\end{aligned}
\tag{14.42}
$$

因此，

$$
\operatorname{div}\boldsymbol A
\equiv
\nabla\cdot\boldsymbol A
=
\frac{\partial A_x}{\partial x}
+
\frac{\partial A_y}{\partial y}
+
\frac{\partial A_z}{\partial z}.
\tag{14.43}
$$

**例 14.14**　在柱坐标系中，利用

$$
*d\rho=\rho\,d\phi\wedge dz,
\qquad
*(d\phi\wedge dz)=\frac1\rho d\rho,
\qquad
*(d\rho\wedge d\phi\wedge dz)=\frac1\rho,
$$

$$
*d\phi=\frac1\rho dz\wedge d\rho,
\qquad
*(dz\wedge d\rho)=\rho\,d\phi,
$$

$$
*dz=\rho\,d\rho\wedge d\phi,
\qquad
*(d\rho\wedge d\phi)=\frac1\rho dz,
$$

得

$$
*d*(A_\rho d\rho+A_\phi\rho d\phi+A_zdz)
=
\frac1\rho
\frac{\partial(\rho A_\rho)}{\partial\rho}
+
\frac1\rho
\frac{\partial A_\phi}{\partial\phi}
+
\frac{\partial A_z}{\partial z}.
\tag{14.44}
$$

因此，

$$
\operatorname{div}\boldsymbol A
\equiv
\nabla\cdot\boldsymbol A
=
\frac1\rho
\frac{\partial(\rho A_\rho)}{\partial\rho}
+
\frac1\rho
\frac{\partial A_\phi}{\partial\phi}
+
\frac{\partial A_z}{\partial z}.
\tag{14.45}
$$

**例 14.15**　在球坐标系中，利用

$$
*dr=r^2\sin\theta\,d\theta\wedge d\phi,
\qquad
*(d\theta\wedge d\phi)=\frac1{r^2\sin\theta}dr,
\qquad
*(dr\wedge d\theta\wedge d\phi)=\frac1{r^2\sin\theta},
$$

$$
*d\theta=\sin\theta\,d\phi\wedge dr,
\qquad
*(d\phi\wedge dr)=\frac1{\sin\theta}d\theta,
$$

$$
*d\phi=\frac1{\sin\theta}dr\wedge d\theta,
\qquad
*(dr\wedge d\theta)=\sin\theta\,d\phi,
$$

得

$$
*d*(A_rdr+A_\theta r\,d\theta+A_\phi r\sin\theta\,d\phi)
=
\frac1{r^2}
\frac{\partial(r^2A_r)}{\partial r}
+
\frac1{r\sin\theta}
\frac{\partial(\sin\theta A_\theta)}{\partial\theta}
+
\frac1{r\sin\theta}
\frac{\partial A_\phi}{\partial\phi}.
\tag{14.46}
$$

因此，

$$
\operatorname{div}\boldsymbol A
\equiv
\nabla\cdot\boldsymbol A
=
\frac1{r^2}
\frac{\partial(r^2A_r)}{\partial r}
+
\frac1{r\sin\theta}
\frac{\partial(\sin\theta A_\theta)}{\partial\theta}
+
\frac1{r\sin\theta}
\frac{\partial A_\phi}{\partial\phi}.
\tag{14.47}
$$

### 正交曲面坐标系中的 Laplace 算符

$*d*d$ 是作用在标量函数上的 Laplace 算符 $\nabla^2\equiv\nabla\cdot\nabla\equiv\operatorname{div}\operatorname{grad}$ 的协变微分形式。

**例 14.16**　直角坐标系，由 (14.22) 式，可得

$$
\begin{aligned}
d*du
&=
d\left(
\frac{\partial u}{\partial x}dy\wedge dz
+
\frac{\partial u}{\partial y}dz\wedge dx
+
\frac{\partial u}{\partial z}dx\wedge dy
\right)\\
&=
\frac{\partial^2u}{\partial x^2}
dx\wedge dy\wedge dz
+
\frac{\partial^2u}{\partial y^2}
dy\wedge dz\wedge dx
+
\frac{\partial^2u}{\partial z^2}
dz\wedge dx\wedge dy\\
&=
\left(
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
+
\frac{\partial^2u}{\partial z^2}
\right)
dx\wedge dy\wedge dz.
\end{aligned}
$$

因此，

$$
*d*du
=
\left(
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}
\right)u.
$$

这就是说，作用在标量函数上的 Laplace 算符在直角坐标系下的表达式是

$$
\nabla^2
\equiv
\nabla\cdot\nabla
\equiv
\operatorname{div}\operatorname{grad}
=
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}.
\tag{14.48}
$$

**例 14.17**　柱坐标系，由 (14.23) 式，可得

$$
d*du
=
\left[
\frac{\partial}{\partial\rho}
\left(
\rho\frac{\partial u}{\partial\rho}
\right)
+
\frac1\rho\frac{\partial^2u}{\partial\phi^2}
+
\rho\frac{\partial^2u}{\partial z^2}
\right]
d\rho\wedge d\phi\wedge dz,
$$

$$
*d*du
=
\frac1\rho
\frac{\partial}{\partial\rho}
\left(
\rho\frac{\partial u}{\partial\rho}
\right)
+
\frac1{\rho^2}
\frac{\partial^2u}{\partial\phi^2}
+
\frac{\partial^2u}{\partial z^2}.
$$

这就是说，作用在标量函数上的 Laplace 算符在柱坐标系下的表达式是

$$
\nabla^2
\equiv
\nabla\cdot\nabla
\equiv
\operatorname{div}\operatorname{grad}
=
\frac1\rho
\frac{\partial}{\partial\rho}
\left(
\rho\frac{\partial}{\partial\rho}
\right)
+
\frac1{\rho^2}
\frac{\partial^2}{\partial\phi^2}
+
\frac{\partial^2}{\partial z^2}.
\tag{14.49}
$$

**例 14.18**　球坐标系，由 (14.24) 式，可得

$$
\begin{aligned}
d*du
&=
\left[
\sin\theta\frac{\partial}{\partial r}
\left(
r^2\frac{\partial u}{\partial r}
\right)
+
\frac{\partial}{\partial\theta}
\left(
\sin\theta\frac{\partial u}{\partial\theta}
\right)
+
\frac1{\sin\theta}
\frac{\partial^2u}{\partial\phi^2}
\right]
dr\wedge d\theta\wedge d\phi,
\end{aligned}
$$

$$
*d*du
=
\frac1{r^2}
\frac{\partial}{\partial r}
\left(
r^2\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2\sin\theta}
\frac{\partial}{\partial\theta}
\left(
\sin\theta\frac{\partial u}{\partial\theta}
\right)
+
\frac1{r^2\sin^2\theta}
\frac{\partial^2u}{\partial\phi^2}.
$$

所以，作用在标量函数上的 Laplace 算符在球坐标系下的表达式是

$$
\nabla^2
\equiv
\nabla\cdot\nabla
\equiv
\operatorname{div}\operatorname{grad}
=
\frac1{r^2}
\frac{\partial}{\partial r}
\left(
r^2\frac{\partial}{\partial r}
\right)
+
\frac1{r^2\sin\theta}
\frac{\partial}{\partial\theta}
\left(
\sin\theta\frac{\partial}{\partial\theta}
\right)
+
\frac1{r^2\sin^2\theta}
\frac{\partial^2}{\partial\phi^2}.
\tag{14.50}
$$

Laplace 算符还可以作用在矢量函数 $\boldsymbol A=\sum_{i=1}^3A_i\boldsymbol e_i$ 上，其定义是

$$
\nabla^2\boldsymbol A
\equiv
\nabla(\nabla\cdot\boldsymbol A)
-
\nabla\times(\nabla\times\boldsymbol A).
\tag{14.51}
$$

因此作用在矢量函数上的 Laplace 算符 $\nabla^2\equiv\nabla(\nabla\cdot)-\nabla\times(\nabla\times)\equiv\operatorname{grad}\operatorname{div}-\operatorname{curl}\operatorname{curl}$ 的协变微分形式是 $d*d*-*d*d$，将 $(d*d*-*d*d)$ 作用在一次微分形式 $\sum_{i=1}^3A_i\sqrt{g_{ii}}dx^i$ 上。

**例 14.19　$\nabla^2\boldsymbol A$ 的直角坐标分解**

由 (14.42) 式和 (14.32) 式，得

$$
d*d*(A_xdx+A_ydy+A_zdz)
=
d\left(
\frac{\partial A_x}{\partial x}
+
\frac{\partial A_y}{\partial y}
+
\frac{\partial A_z}{\partial z}
\right),
$$

$$
\begin{aligned}
&=
\left(
\frac{\partial^2A_x}{\partial x^2}
+
\frac{\partial^2A_y}{\partial x\partial y}
+
\frac{\partial^2A_z}{\partial x\partial z}
\right)dx\\
&\quad+
\left(
\frac{\partial^2A_x}{\partial x\partial y}
+
\frac{\partial^2A_y}{\partial y^2}
+
\frac{\partial^2A_z}{\partial y\partial z}
\right)dy\\
&\quad+
\left(
\frac{\partial^2A_x}{\partial x\partial z}
+
\frac{\partial^2A_y}{\partial y\partial z}
+
\frac{\partial^2A_z}{\partial z^2}
\right)dz,
\end{aligned}
$$

$$
*d*d(A_xdx+A_ydy+A_zdz)
=
*d
\left[
\left(
\frac{\partial A_z}{\partial y}
-
\frac{\partial A_y}{\partial z}
\right)dx
+
\left(
\frac{\partial A_x}{\partial z}
-
\frac{\partial A_z}{\partial x}
\right)dy
+
\left(
\frac{\partial A_y}{\partial x}
-
\frac{\partial A_x}{\partial y}
\right)dz
\right]
$$

$$
\begin{aligned}
&=
\left[
\frac{\partial}{\partial y}
\left(
\frac{\partial A_y}{\partial x}
-
\frac{\partial A_x}{\partial y}
\right)
-
\frac{\partial}{\partial z}
\left(
\frac{\partial A_x}{\partial z}
-
\frac{\partial A_z}{\partial x}
\right)
\right]dx\\
&\quad+
\left[
\frac{\partial}{\partial z}
\left(
\frac{\partial A_z}{\partial y}
-
\frac{\partial A_y}{\partial z}
\right)
-
\frac{\partial}{\partial x}
\left(
\frac{\partial A_y}{\partial x}
-
\frac{\partial A_x}{\partial y}
\right)
\right]dy\\
&\quad+
\left[
\frac{\partial}{\partial x}
\left(
\frac{\partial A_x}{\partial z}
-
\frac{\partial A_z}{\partial x}
\right)
-
\frac{\partial}{\partial y}
\left(
\frac{\partial A_z}{\partial y}
-
\frac{\partial A_y}{\partial z}
\right)
\right]dz.
\end{aligned}
$$

上面两式相减，得 $\nabla^2\boldsymbol A$ 的直角坐标分量表达式

$$
\nabla^2\boldsymbol A
=
(\nabla^2A)_x\boldsymbol e_x
+
(\nabla^2A)_y\boldsymbol e_y
+
(\nabla^2A)_z\boldsymbol e_z
$$

$$
=
\left[
\left(
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}
\right)
A_x
\right]\boldsymbol e_x
+
\left[
\left(
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}
\right)
A_y
\right]\boldsymbol e_y
\tag{14.52}
$$

$$
\quad+
\left[
\left(
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}
\right)
A_z
\right]\boldsymbol e_z.
\tag{14.53}
$$

**例 14.20　$\nabla^2\boldsymbol A$ 的柱坐标分解**

利用 (14.44) 式和 (14.36) 式就能求出 $\nabla^2\boldsymbol A$ 的柱坐标分量表达式

$$
\nabla^2\boldsymbol A
=
(\nabla^2A)_\rho\boldsymbol e_\rho
+
(\nabla^2A)_\phi\boldsymbol e_\phi
+
(\nabla^2A)_z\boldsymbol e_z
$$

$$
=
\left(
\nabla^2A_\rho
-
\frac1{\rho^2}A_\rho
-
\frac2{\rho^2}
\frac{\partial A_\phi}{\partial\phi}
\right)\boldsymbol e_\rho
+
\left(
\nabla^2A_\phi
-
\frac1{\rho^2}A_\phi
+
\frac2{\rho^2}
\frac{\partial A_\rho}{\partial\phi}
\right)\boldsymbol e_\phi
+
\nabla^2A_z\boldsymbol e_z.
\tag{14.54}
$$

**例 14.21　$\nabla^2\boldsymbol A$ 的球坐标分解**

利用 (14.46) 式和 (14.40) 式可以得到 $\nabla^2\boldsymbol A$ 的球坐标分量表达式

$$
\nabla^2\boldsymbol A
=
(\nabla^2A)_r\boldsymbol e_r
+
(\nabla^2A)_\theta\boldsymbol e_\theta
+
(\nabla^2A)_\phi\boldsymbol e_\phi
$$

$$
\begin{aligned}
&=
\left(
\nabla^2A_r
-
\frac2{r^2}A_r
-
\frac2{r^2\sin\theta}
\frac{\partial(\sin\theta A_\theta)}{\partial\theta}
-
\frac2{r^2\sin\theta}
\frac{\partial A_\phi}{\partial\phi}
\right)\boldsymbol e_r\\
&\quad+
\left(
\nabla^2A_\theta
-
\frac1{r^2\sin^2\theta}A_\theta
+
\frac2{r^2}
\frac{\partial A_r}{\partial\theta}
-
\frac{2\cos\theta}{r^2\sin^2\theta}
\frac{\partial A_\phi}{\partial\phi}
\right)\boldsymbol e_\theta\\
&\quad+
\left(
\nabla^2A_\phi
-
\frac1{r^2\sin^2\theta}A_\phi
+
\frac2{r^2\sin\theta}
\frac{\partial A_r}{\partial\phi}
+
\frac{2\cos\theta}{r^2\sin^2\theta}
\frac{\partial A_\theta}{\partial\phi}
\right)\boldsymbol e_\phi.
\end{aligned}
\tag{14.55}
$$

我们将上面用到的矢量分析与外微分运算的对应关系总结在下表中：

| 矢量分析 | 外微分运算 |
|---|---|
| $n$ 维空间体积元 | $\sqrt{\det G}\,d\xi^1\wedge d\xi^2\wedge d\xi^3\cdots\wedge d\xi^n$ |
| 坐标单位矢量 $\boldsymbol e_i$ | 一次微分形式 $\sqrt{g_{ii}}\,d\xi^i$ |
| 矢量 $\boldsymbol A=\sum_{i=1}^nA_i\boldsymbol e_i$ | 一次微分形式 $\sum_{i=1}^nA_i\sqrt{g_{ii}}\,d\xi^i$ |
| 标量 $u$ 的梯度 $\operatorname{grad}u\equiv\nabla u$ | 外微分 $d$ 作用于标量函数 $u$ |
| 矢量 $\boldsymbol A$ 的旋度 $\operatorname{curl}\boldsymbol A\equiv\nabla\times\boldsymbol A$ | $*d$ 作用于一次微分形式 $\sum_{i=1}^nA_i\sqrt{g_{ii}}\,d\xi^i$ |
| 矢量 $\boldsymbol A$ 的散度 $\operatorname{div}\boldsymbol A\equiv\nabla\cdot\boldsymbol A$ | $*d*$ 作用于一次微分形式 $\sum_{i=1}^nA_i\sqrt{g_{ii}}\,d\xi^i$ |
| 作用于标量函数 $u$ 的 Laplace 算符 $\nabla^2\equiv\nabla\cdot\nabla\equiv\operatorname{div}\operatorname{grad}$ | $*d*d$ 作用于标量 $u$ |
| 作用于矢量函数 $\boldsymbol A$ 的 Laplace 算符 $\nabla^2\equiv\nabla(\nabla\cdot)-\nabla\times(\nabla\times)\equiv\operatorname{grad}\operatorname{div}-\operatorname{curl}\operatorname{curl}$ | $(d*d*-*d*d)$ 作用于一次微分形式 $\sum_{i=1}^nA_i\sqrt{g_{ii}}\,d\xi^i$ |


## *§14.3 Laplace 算符的平移、转动和反射不变性

根据问题的需要，我们可以选择合适的坐标系。选定了坐标系以后，在求解定解问题时，往往还需要考虑两个问题。一是坐标架如何放置，包括坐标原点位置和坐标轴取向的选择，以最大限度地利用问题中的对称性，使求解过程得到充分的简化。二是讨论定解问题的对称性与解的对称性之间的必然联系。

坐标架的不同放置，在数学上，就表现为不同坐标系之间的线性变换。这一节就讨论一些具体的变换，证明 Laplace 算符在这些变换下的不变性。为了确定起见，下面的讨论在直角坐标系中进行。当然，得到的结论在任意正交曲面坐标系中都成立。

首先，关于坐标原点的不同选择，涉及的是平移变换

$$
x'=x-a,\qquad
y'=y-b,\qquad
z'=z-c.
\tag{14.56}
$$

容易看出，Laplace 算符在平移变换下是不变的，即

$$
\frac{\partial^2}{\partial x'^2}
+
\frac{\partial^2}{\partial y'^2}
+
\frac{\partial^2}{\partial z'^2}
\equiv
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}.
\tag{14.57}
$$

关于坐标轴的取向问题，涉及坐标系之间的正交变换。设空间一点在变换前后的坐标分别是 $(x,y,z)$ 和 $(x',y',z')$，它们之间的关系是

$$
\begin{pmatrix}
x\\y\\z
\end{pmatrix}
=
\begin{pmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}
\end{pmatrix}
\begin{pmatrix}
x'\\y'\\z'
\end{pmatrix}.
\tag{14.58}
$$

所谓正交变换，指的是变换矩阵

$$
\boldsymbol A
=
\begin{pmatrix}
a_{11}&a_{12}&a_{13}\\
a_{21}&a_{22}&a_{23}\\
a_{31}&a_{32}&a_{33}
\end{pmatrix}
\tag{14.59}
$$

满足正交关系

$$
\sum_{k=1}^3a_{ik}a_{jk}
=
\sum_{k=1}^3a_{ki}a_{kj}
=
\delta_{ij}.
\tag{14.60}
$$

在正交变换之下，Laplace 算符的形式是不变的。这只要证明变换后的度规矩阵仍为单位矩阵即可。为了书写方便，下面把变换前后的坐标改写为 $(x^1,x^2,x^3)$ 和 $(x'^1,x'^2,x'^3)$。根据 (14.58) 式就有

$$
dx^i
=
\sum_{k=1}^3a_{ik}dx'^k.
\tag{14.61}
$$

容易得到

$$
\begin{aligned}
ds^2
&=
\sum_{ij}\delta_{ij}dx^i dx^j\\
&=
\sum_{ij}\sum_{kl}
\delta_{ij}a_{ik}a_{jl}dx'^kdx'^l\\
&=
\sum_{kl}
\left(
\sum_i a_{ik}a_{il}
\right)
dx'^kdx'^l
=
\sum_{kl}
\delta_{kl}dx'^kdx'^l.
\end{aligned}
\tag{14.62}
$$

这就证明了变换后的度规矩阵仍是单位矩阵，所以变换后的 Laplace 算符仍为

$$
\nabla^2
=
\frac{\partial^2}{\partial x'^2}
+
\frac{\partial^2}{\partial y'^2}
+
\frac{\partial^2}{\partial z'^2}.
\tag{14.63}
$$

正交变换的矩阵 $\boldsymbol A$ 中有 $3\times3$ 个元素。由于正交关系含有 $3\times4/2$ 个限制条件，因此只有 3 个矩阵元是独立的，或者说只有 3 个自由参数。这三个参数不妨取为描写刚体转动的三个 Euler 角。这等价于坐标架绕（过原点的）固定轴的转动。因此 Laplace 算符在绕（过原点的）任意固定轴的转动下也是不变的。

最后，容易看出，在空间反射

$$
x'=-x,\qquad
y'=-y,\qquad
z'=-z
\tag{14.64}
$$

下，Laplace 算符也是不变的。

讨论 Laplace 算符的不变性，一方面是求解定解问题的直接需要（选择坐标系），另一方面，这些不变性也必然会影响到相关定解问题的解的对称性。

## §14.4 圆形区域

现在讨论圆形区域中的稳定问题。定解问题为

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
0,
\qquad
x^2+y^2<a^2,
\tag{14.65a}
$$

$$
\left.u\right|_{x^2+y^2=a^2}
=
f.
\tag{14.65b}
$$

在直角坐标系下方程 (14.65a) 可以分离变量（例如见 §13.2），但边界条件 (14.65b) 显然不能分离变量。由于边界的形状是圆形，很自然地，应该采用平面极坐标系。

在平面极坐标系中，定解问题 (14.65) 似乎可以改写为

$$
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2}
\frac{\partial^2u}{\partial\phi^2}
=
0,
\qquad
0<r<a,\quad
0\le\phi\le2\pi,
\tag{14.66a}
$$

$$
\left.u\right|_{r=a}
=
f(\phi),
\tag{14.66b}
$$

但方程 (14.66a) 和边界条件 (14.66b) 并不构成适定的定解问题，或者说，并不完全等价于原来的定解问题 (14.65)。原因是原来的方程 (14.65a) 在圆内处处成立，包括在坐标原点 $(r=0)$ 以及实轴 $(\phi=0$ 或 $2\pi)$ 也成立。但方程 (14.66a) 在 $r=0$ 点并不成立（因为 $\partial u(r,\phi)/\partial r$ 在 $r=0$ 点无定义）。$r=0$ 点作为（平面极坐标系）自变量 $r$ 的一个端点，需要补充上 $u(r,\phi)$ 在 $r=0$ 点所应当满足的边界条件。考虑到方程 (14.65a) 是齐次的（“无源”），因此，$u(r,\phi)$ 在坐标原点应当有界①，要补充上有界条件

$$
\left.u(r,\phi)\right|_{r=0}
\text{ 有界}.
$$

① 不妨用大家熟悉的静电场的知识来理解这个问题。三维空间中，只有在点电荷和线电荷所在处电势分别以 $1/r$ 和 $\ln r$ 的形式发散，除此之外，电势是处处有界的。二维平面上的点电荷，就是三维空间中的均匀线电荷。

同样，自变量 $\phi$ 的变化范围是 $[0,2\pi]$，方程 (14.66a) 在区间的端点 $\phi=0$ 和 $\phi=2\pi$ 并不成立（因为 $\partial u(r,\phi)/\partial\phi$ 在端点 $\phi=0$ 和 $\phi=2\pi$ 处无定义）。考虑到平面极坐标系的特点，$(r,\phi=0)$ 和 $(r,\phi=2\pi)$ 代表的是平面上的同一点，所以，作为完整的定解问题，还应当补充上周期条件

$$
\left.u(r,\phi)\right|_{\phi=0}
=
\left.u(r,\phi)\right|_{\phi=2\pi}
$$

和

$$
\left.
\frac{\partial u(r,\phi)}{\partial\phi}
\right|_{\phi=0}
=
\left.
\frac{\partial u(r,\phi)}{\partial\phi}
\right|_{\phi=2\pi}.
$$

总结上面的讨论，我们看到，在转换到平面极坐标系后，定解问题 (14.65) 应该变为

$$
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2}
\frac{\partial^2u}{\partial\phi^2}
=
0,
\qquad
0<\phi<2\pi,\quad
0<r<a,
\tag{14.67a}
$$

$$
\left.u(r,\phi)\right|_{\phi=0}
=
\left.u(r,\phi)\right|_{\phi=2\pi},
\qquad
0\le r\le a,
\tag{14.67b}
$$

$$
\left.
\frac{\partial u(r,\phi)}{\partial\phi}
\right|_{\phi=0}
=
\left.
\frac{\partial u(r,\phi)}{\partial\phi}
\right|_{\phi=2\pi},
\qquad
0\le r\le a,
\tag{14.67c}
$$

$$
\left.u(r,\phi)\right|_{r=0}
\text{ 有界},
\qquad
0\le\phi\le2\pi,
\tag{14.67d}
$$

$$
\left.u(r,\phi)\right|_{r=a}
=
f(\phi),
\qquad
0\le\phi\le2\pi.
\tag{14.67e}
$$

现在应用分离变量法求上述定解问题的解。令 $u(r,\phi)=R(r)\Phi(\phi)$，代入方程 (14.67a) 和边界条件 (14.67b)、(14.67c)，分离变量，就可以得到

$$
r\frac{d}{dr}
\left(
r\frac{dR}{dr}
\right)
-
\lambda R
=
0
\tag{14.68}
$$

和

$$
\frac{d^2\Phi}{d\phi^2}
+
\lambda\Phi
=
0,
\tag{14.69a}
$$

$$
\Phi(0)=\Phi(2\pi),
\qquad
\Phi'(0)=\Phi'(2\pi).
\tag{14.69b}
$$

其中 $\lambda$ 是分离变量时引入的常数。这样，由带有待定参数的齐次常微分方程 (14.69a) 和周期边界条件 (14.69b) 又构成了一个本征值问题。不难求出，$\lambda=0$ 是本征值，即本征值问题 (14.69) 有一组解

本征值

$$
\lambda_0=0,
\tag{14.70a}
$$

本征函数

$$
\Phi_0(\phi)=1.
\tag{14.70b}
$$

当 $\lambda\ne0$ 时，方程 (14.69a) 的通解为 $\Phi(\phi)=A\sin\sqrt{\lambda}\,\phi+B\cos\sqrt{\lambda}\,\phi$。代入周期条件 (14.69b)，就得到

$$
B
=
A\sin\sqrt{\lambda}\,2\pi
+
B\cos\sqrt{\lambda}\,2\pi,
$$

$$
A
=
A\cos\sqrt{\lambda}\,2\pi
-
B\sin\sqrt{\lambda}\,2\pi.
$$

这是关于 $A$ 和 $B$ 的线性齐次代数方程组，有非零解的充要条件是

$$
\begin{vmatrix}
\sin\sqrt{\lambda}\,2\pi &
\cos\sqrt{\lambda}\,2\pi-1\\
\cos\sqrt{\lambda}\,2\pi-1 &
-\sin\sqrt{\lambda}\,2\pi
\end{vmatrix}
=
0,
$$

即 $2(\cos\sqrt{\lambda}\,2\pi-1)=0$。这样又可以求得本征值

$$
\lambda_m=m^2,
\qquad
m=1,2,3,\cdots,
\tag{14.71a}
$$

相应的非零解是

$$
A\text{ 任意},
\qquad
B\text{ 任意}.
$$

这就是说，对应于一个本征值 $\lambda_m$，有两个本征函数

$$
\Phi_{m1}(\phi)
=
\sin m\phi,
\tag{14.71b}
$$

$$
\Phi_{m2}(\phi)
=
\cos m\phi.
\tag{14.71c}
$$

当然，也还可以将 $\lambda_0=0$ 的结果 (14.70) 和 (14.71) 合并，统一写成 (14.71) 的形式，而将 (14.71a) 和 (14.71c) 中 $m$ 的取值相应地改为 $0,1,2,3,\cdots$。

再来求常微分方程 (14.68) 的解。注意这个常微分方程是一个特殊的变系数常微分方程（Euler 型微分方程），经过自变量的变换

$$
\frac d{dt}
=
r\frac d{dr},
\qquad
\text{即}
\qquad
t=\ln r
\tag{14.72}
$$

后，就可以变为常系数的常微分方程

$$
\frac{d^2R}{dt^2}
-
\lambda R
=
0.
\tag{14.68'}
$$

所以，当 $\lambda_0=0$ 时，通解为

$$
R_0(r)
=
C_0+D_0t
=
C_0+D_0\ln r;
\tag{14.73}
$$

当 $\lambda_m=m^2,\ m\ne0$ 时，通解为

$$
R_m(r)
=
C_me^{mt}
+
D_me^{-mt}
=
C_mr^m+D_mr^{-m}.
\tag{14.74}
$$

这样，就求得了满足齐次方程 (14.67a) 和周期条件 (14.67b)、(14.67c) 的全部特解：

$$
u_0(r,\phi)
=
C_0+D_0\ln r,
\tag{14.75a}
$$

$$
u_{m1}(r,\phi)
=
(C_{m1}r^m+D_{m1}r^{-m})
\sin m\phi,
\tag{14.75b}
$$

$$
u_{m2}(r,\phi)
=
(C_{m2}r^m+D_{m2}r^{-m})
\cos m\phi.
\tag{14.75c}
$$

把它们叠加起来，就得到方程 (14.67a) 在周期条件 (14.67b)、(14.67c) 之下的一般解：

$$
\begin{aligned}
u(r,\phi)
&=
C_0+D_0\ln r
+
\sum_{m=1}^{\infty}
(C_{m1}r^m+D_{m1}r^{-m})
\sin m\phi\\
&\quad+
\sum_{m=1}^{\infty}
(C_{m2}r^m+D_{m2}r^{-m})
\cos m\phi.
\end{aligned}
\tag{14.76}
$$

考虑到 $r=0$ 处的有界条件 (14.67d)，$\ln r$ 和 $r^{-m}$ 项的系数都必须为 0，

$$
D_0=0,
\qquad
D_{m1}=0,
\qquad
D_{m2}=0,
\tag{14.77}
$$

这是因为 $\ln r$ 和 $r^{-m}$ 在 $r=0$ 处均无界。再代入边界条件 (14.67e)，就得到

$$
\left.u(r,\phi)\right|_{r=a}
=
C_0
+
\sum_{m=1}^{\infty}
a^m
(C_{m1}\sin m\phi+C_{m2}\cos m\phi)
=
f(\phi).
$$

下面的问题便是如何定出叠加系数 $C_0,C_{m1}$ 和 $C_{m2}$。尽管现在也可以根据 Fourier 展开的公式定出系数 $C_0,C_{m1}$ 和 $C_{m2}$，但我们还是采用分离变量法的标准做法，利用本征函数的正交性定叠加系数。

容易证明，对于本征值问题 (14.69)，对应不同本征值的本征函数正交：

$$
\int_0^{2\pi}
\sin m\phi\,d\phi
=
0,
\tag{14.78a}
$$

$$
\int_0^{2\pi}
\cos m\phi\,d\phi
=
0,
\qquad
m\ne0,
\tag{14.78b}
$$

$$
\int_0^{2\pi}
\sin n\phi\sin m\phi\,d\phi
=
0,
\qquad
m\ne n,
\tag{14.78c}
$$

$$
\int_0^{2\pi}
\sin n\phi\cos m\phi\,d\phi
=
0,
\qquad
m\ne n,
\tag{14.78d}
$$

$$
\int_0^{2\pi}
\cos n\phi\cos m\phi\,d\phi
=
0,
\qquad
m\ne n.
\tag{14.78e}
$$

而且，对应于同一个本征值 $\lambda_m=m^2$ 的两个本征函数 $\sin m\phi$ 和 $\cos m\phi$ 也正交：

$$
\int_0^{2\pi}
\sin m\phi\cos m\phi\,d\phi
=
0.
\tag{14.78f}
$$

因此，利用本征函数的正交性以及

$$
\int_0^{2\pi}\sin^2m\phi\,d\phi=\pi,
\qquad
\int_0^{2\pi}\cos^2m\phi\,d\phi
=
\pi(1+\delta_{m0}),
$$

就可求得

$$
C_0
=
\frac1{2\pi}
\int_0^{2\pi}
f(\phi)\,d\phi,
\tag{14.79a}
$$

$$
C_{m1}
=
\frac1{a^m\pi}
\int_0^{2\pi}
f(\phi)\sin m\phi\,d\phi,
\tag{14.79b}
$$

$$
C_{m2}
=
\frac1{a^m\pi}
\int_0^{2\pi}
f(\phi)\cos m\phi\,d\phi.
\tag{14.79c}
$$

这样我们就完全给出了定解问题 (14.67) 的解。

**练习 14.1**　直接从本征值问题 (14.69) 的方程和边界条件出发，证明本征函数的正交性 (14.78a)–(14.78e)。

**练习 14.2**　证明本征函数的正交性 (14.78f)。

现在再对上面求解过程中的某些问题做一些补充讨论。

(1) 关于定解问题 (14.67) 和 (14.65) 的等价性。

在直角坐标系中，$u(x,y)$ 在圆内满足 Laplace 方程，因此可以找到 $u(x,y)$ 的共轭调和函数 $v(x,y)$，使 $w(z)=w(x+iy)=u(x,y)+iv(x,y)$ 在圆内 $x^2+y^2<a^2$ 解析。但是，一旦选择了平面极坐标系，则定解问题转变为 (14.67)。Laplace 方程只在 $0<\phi<2\pi,\ 0<r<a$ 的范围内成立。周期条件的引入，排除了 $w(z)=u(r,\phi)+iv(r,\phi)$ 的多值性，从而将 $w(z)$ 的解析区域扩充为环域 $0<|z|<a$。而有界条件 (14.67d) 又排除了 $z=0$ 是 $w(z)$ 的极点或本性奇点的可能性，将 $w(z)$ 的解析区域进一步扩充为圆域 $x^2+y^2<a^2$（因此 $w(z)$ 及其实部与虚部都在此圆域内处处满足 Laplace 方程）。这同样体现在偏微分方程 (14.67a) 的特解

$$
1,\qquad
\ln r,\qquad
r^m\sin m\phi,\qquad
r^m\cos m\phi,\qquad
r^{-m}\sin m\phi,\qquad
r^{-m}\cos m\phi
$$

上，它们正是解析函数

$$
z^0,\qquad
\ln z,\qquad
z^m,\qquad
z^{-m}
$$

的实部或者虚部。周期条件 (14.67b)、(14.67c) 和有界条件 (14.67d) 的作用就在于排除掉在圆内 $|z|<a$ 并不处处解析的函数 $\ln z$ 和 $z^{-m}$。

(2) 圆内 Laplace 方程边值问题中的边界条件。

从上面的求解过程可以看出，圆内 Laplace 方程边值问题的本征函数，是由 Laplace 方程分离变量后得到的 (14.69a) 以及周期条件 (14.69b) 决定的，而圆周上的边界条件，只是用来定叠加系数。因此，哪怕圆周上的边界条件是非齐次的，也无需齐次化。事实上，对于圆内 Laplace 方程的边值问题，如果遇到的是齐次的第一类边界条件，那么这个定解问题就只有零解。至于圆内 Poisson 方程的边值问题，在将方程齐次化时，选择的齐次化函数自然应当是周期函数，同时也应当在 $r=0$ 处有界，当然在圆周上也应该有界，但是就不必过多担心对圆周上边界条件的影响。

这里还要顺便指出，周期条件在构成圆形区域乃至环形区域的定解问题是必不可少的。但是，如果求解区域不是完整的圆形或环形，例如是半圆，甚至是带有割线的圆形，周期条件就不再成立。在作为区域边界一部分的半径上，就必须给出相应的边界条件。

(3) 本征值问题的简并现象。

在求解本征值问题 (14.69) 时，对应于一个本征值得到两个（线性无关的）本征函数。我们把对应一个本征值有不止一个（线性无关的）本征函数的现象，称为简并（或退化）。如果对应一个本征值有 $n$ 个本征函数，则称本征值问题是 $n$ 重简并的，或者说简并度为 $n$。对于二阶常微分方程的本征值问题，最多只能是二重简并的。我们还将证明（见第十七章），在二阶常微分方程的本征值问题中，如果边界条件是一、二、三类，则对应一个本征值，只能有一个本征函数，或者说，本征值问题一定是非简并的。

(4) 简并情形下本征函数的正交性问题。

对于简并的本征值问题，本征函数的选取并不唯一。对应同一个本征值的本征函数也不一定正交，但是一定可以通过适当的重新组合而使它们正交化。如果将本征值问题 (14.69) 的对应于 $\lambda_m=m^2,\ m=1,2,3,\cdots$ 的本征函数取为

$$
\sin m\phi
\qquad\text{和}\qquad
e^{im\phi},
\tag{14.80}
$$

这两个本征函数线性无关但是不正交，因为

$$
\int_0^{2\pi}
(\sin m\phi)^*e^{im\phi}\,d\phi
=
i\pi
\ne0.
\tag{14.81}
$$

现在的本征函数是复函数，在上面的正交关系中需要将其中的一个本征函数取复共轭。这样做的直接原因就是为了保证本征函数的模方恒为正值。而如果将本征值问题 (14.69) 的对应于 $\lambda_m=m^2,\ m=1,2,3,\cdots$ 的本征函数取为

$$
e^{im\phi}
\qquad\text{和}\qquad
e^{-im\phi},
$$

这两个本征函数是正交的：

$$
\int_0^{2\pi}
(e^{im\phi})^*e^{-im\phi}\,d\phi
=
0.
\tag{14.82}
$$

可以简单地将本征值问题 (14.69) 的全部本征值（包括 $\lambda_0=0$）和本征函数统一写成

$$
\lambda_m=m^2,
\qquad
m=0,\pm1,\pm2,\cdots,
\tag{14.83a}
$$

$$
\Phi_m(\phi)
=
e^{im\phi}.
\tag{14.83b}
$$

这时，对应不同本征值的本征函数当然仍然是正交的：

$$
\int_0^{2\pi}
(e^{im\phi})^*e^{in\phi}\,d\phi
=
0,
\qquad
m,n=0,\pm1,\pm2,\cdots,\quad n\ne m.
\tag{14.84}
$$

**练习 14.3**　证明本征函数组 (14.83b) 的正交性 (14.82)。

(5) 圆内 Laplace 方程第一边值问题的解与 Poisson 公式。

把上面求得的系数 (14.77) 和 (14.79) 代入到 (14.76) 式中，还可以得到

$$
\begin{aligned}
u(r,\phi)
&=
\frac1{2\pi}
\int_0^{2\pi}
f(\phi')\,d\phi'
+
\frac1\pi
\sum_{m=1}^{\infty}
\left(
\frac ra
\right)^m
\sin m\phi
\int_0^{2\pi}
f(\phi')\sin m\phi'\,d\phi'\\
&\quad+
\frac1\pi
\sum_{m=1}^{\infty}
\left(
\frac ra
\right)^m
\cos m\phi
\int_0^{2\pi}
f(\phi')\cos m\phi'\,d\phi'\\
&=
\frac1{2\pi}
\int_0^{2\pi}
f(\phi')
\left[
1+
2\sum_{m=1}^{\infty}
\left(
\frac ra
\right)^m
\cos m(\phi-\phi')
\right]d\phi'.
\end{aligned}
$$

显然，当 $r<a$ 时级数收敛。将余弦函数改写为复指数函数，利用等比级数的求和公式就可以求出级数的和，就得到

$$
u(r,\phi)
=
\frac{a^2-r^2}{2\pi}
\int_0^{2\pi}
\frac{f(\phi')}
{r^2+a^2-2ar\cos(\phi-\phi')}
\,d\phi'.
\tag{14.85}
$$

这个结果称为 Poisson 积分公式，它把 Laplace 方程在圆内第一类边值问题的解表示为边值 $f(\phi)$ 的积分。事实上，由解析函数的 Cauchy 积分公式，也可以推出这个结果（见 §3.7，(3.39) 式），而 $u(r,\phi)$ 正好是解析函数的实部（或虚部）。

**例 14.22**　求解圆内 Poisson 方程的定解问题：

$$
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2}
\frac{\partial^2u}{\partial\phi^2}
=
-r^2\cos2\phi,
\qquad
0<r<a,\quad
0<\phi<2\pi,
$$

$$
\left.u\right|_{\phi=0}
=
\left.u\right|_{\phi=2\pi},
\qquad
\left.
\frac{\partial u}{\partial\phi}
\right|_{\phi=0}
=
\left.
\frac{\partial u}{\partial\phi}
\right|_{\phi=2\pi},
\qquad
0\le r\le a,
$$

$$
\left.u\right|_{r=0}
\text{ 有界},
\qquad
\left.u\right|_{r=a}
=
0,
\qquad
0\le\phi\le2\pi.
$$

**解**　这是非齐次方程的定解问题，可以采用两种方法求解。

**解法一　方程齐次化。**　将方程改写到直角坐标系，有

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
-(x^2-y^2),
$$

因此可取齐次化函数为 $A(x^4-y^4)=Ar^4\cos2\phi$，代入 Poisson 方程，可以定出 $A=-1/12$。令

$$
u(r,\phi)
=
-\frac1{12}r^4\cos2\phi+w(r,\phi),
$$

即可得到 $w(r,\phi)$ 所满足的定解问题

$$
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial w}{\partial r}
\right)
+
\frac1{r^2}
\frac{\partial^2w}{\partial\phi^2}
=
0,
\qquad
0<r<a,\quad
0<\phi<2\pi,
$$

$$
\left.w\right|_{\phi=0}
=
\left.w\right|_{\phi=2\pi},
\qquad
\left.
\frac{\partial w}{\partial\phi}
\right|_{\phi=0}
=
\left.
\frac{\partial w}{\partial\phi}
\right|_{\phi=2\pi},
\qquad
0\le r\le a,
$$

$$
\left.w\right|_{r=0}
\text{ 有界},
\qquad
\left.w\right|_{r=a}
=
\frac1{12}a^4\cos2\phi,
\qquad
0\le\phi\le2\pi.
$$

容易写出满足齐次方程、周期条件以及 $r=0$ 处有界条件的一般解

$$
w(r,\phi)
=
\sum_{m=0}^{\infty}
r^m
(A_m\cos m\phi+B_m\sin m\phi).
$$

代入 $r=a$ 处的边界条件，

$$
w(a,\phi)
=
\sum_{m=0}^{\infty}
a^m
(A_m\cos m\phi+B_m\sin m\phi)
=
\frac1{12}a^4\cos2\phi,
$$

即可定出

$$
A_m
=
\frac1{12}a^2\delta_{m2},
\qquad
B_m=0,
\qquad
m=0,1,2,\cdots.
$$

于是就得到原始 Poisson 方程边值问题的解为

$$
u(r,\phi)
=
\frac1{12}r^2(a^2-r^2)\cos2\phi.
$$

**解法二　按相应齐次问题本征函数展开。**　相应齐次问题的本征函数已在上面解法一中给出，所以就可以令

$$
u(r,\phi)
=
R_0(r)
+
\sum_{m=1}^{\infty}
\left[
R_{m1}(r)\cos m\phi
+
R_{m2}(r)\sin m\phi
\right].
$$

代入方程及关于 $r$ 边界条件，利用本征函数的正交性，就能导出 $R_0(r),R_{m1}(r)$ 及 $R_{m2}(r)$ 满足的常微分方程边值问题：

$$
\begin{cases}
\dfrac1r\dfrac d{dr}\left(r\dfrac{dR_0}{dr}\right)=0,\\
R_0(0)\text{ 有界},\qquad R_0(a)=0;
\end{cases}
$$

$$
\begin{cases}
\dfrac1r\dfrac d{dr}\left(r\dfrac{dR_{m1}}{dr}\right)
-\dfrac{m^2}{r^2}R_{m1}
=
-r^2\delta_{m2},\\
R_{m1}(0)\text{ 有界},\qquad R_{m1}(a)=0;
\end{cases}
$$

$$
\begin{cases}
\dfrac1r\dfrac d{dr}\left(r\dfrac{dR_{m2}}{dr}\right)
-\dfrac{m^2}{r^2}R_{m2}
=
0,\\
R_{m2}(0)\text{ 有界},\qquad R_{m2}(a)=0.
\end{cases}
$$

解之即得

$$
R_0(r)=0,
\qquad
R_{m1}(r)
=
\frac1{12}r^2(a^2-r^2)\delta_{m2},
\qquad
R_{m2}(r)=0.
$$

最后也就求得

$$
u(r,\phi)
=
\frac1{12}r^2(a^2-r^2)\cos2\phi.
$$

与解法一完全相同。

**讨论**　其实本题还可以有第三种解法：直接令 $u(r,\phi)=R(r)\cos2\phi$，代入定解问题，求出 $R(r)$ 即可。请读者完成此计算，并且考虑：这种解法的根据是什么？

**例 14.23**　设有扇形区域 $-\alpha<\phi<\alpha,\ r>0$。两边上 $r<a$ 的线段内温度为 1，其余部分温度为 0。求扇形区域内的稳定温度分布。

**解**　设扇形区域内 $(r,\phi)$ 处的温度为 $u(r,\phi)$，它满足的定解问题为

$$
\frac1r
\frac{\partial}{\partial r}
\left(
r\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2}
\frac{\partial^2u}{\partial\phi^2}
=
0,
\qquad
r>0,\quad
-\alpha<\phi<\alpha,
\tag{14.86a}
$$

$$
\left.u\right|_{r=0}
\text{ 有界},
\qquad
\left.u\right|_{r\to\infty}
\to0,
\qquad
-\alpha\le\phi\le\alpha,
\tag{14.86b}
$$

$$
\left.u\right|_{\phi=\pm\alpha}
=
\eta(a-r),
\qquad
r\ge0,
\tag{14.86c}
$$

其中 $\eta(x)$ 是 Heaviside 单位阶跃函数。和圆形区域的定解问题不同，现在需要用 $\phi=\pm\alpha$ 上齐次边界条件以及相应的微分方程构成本征值问题。显然，无论是 $r>a$ 时，或是 $r<a$ 但将边界条件齐次化（齐次化函数就是 1）后，将 Laplace 方程以及齐次边界条件分离变量，都会得到本征值问题

$$
\Phi''+\lambda\Phi=0,
$$

$$
\Phi(\pm\alpha)=0.
$$

这样就能定出

本征值

$$
\lambda_n
=
\left(
\frac{n\pi}{2\alpha}
\right)^2,
\qquad
n=1,2,3,\cdots,
$$

本征函数

$$
\Phi_n(\phi)
=
\sin\frac{n\pi}{2\alpha}(\phi+\alpha).
$$

相应的径向方程

$$
\frac1r
\frac d{dr}
\left(
r\frac{dR_n}{dr}
\right)
-
\frac{\lambda_n}{r^2}R_n
=
0
$$

的解就可以取为

$$
R_n(r)
=
A_n
\left(
\frac ra
\right)^{n\pi/2\alpha}
+
B_n
\left(
\frac ar
\right)^{n\pi/2\alpha}.
$$

考虑到 $r=0$ 处的有界条件以及 $r\to\infty$ 时的无穷远条件，就能写出一般解

$$
u(r,\phi)
=
\begin{cases}
1+
\displaystyle
\sum_{n=1}^{\infty}
C_n
\left(
\frac ra
\right)^{n\pi/2\alpha}
\sin\frac{n\pi}{2\alpha}(\phi+\alpha),
& r<a,\\[14pt]
\displaystyle
\sum_{n=1}^{\infty}
D_n
\left(
\frac ar
\right)^{n\pi/2\alpha}
\sin\frac{n\pi}{2\alpha}(\phi+\alpha),
& r>a.
\end{cases}
$$

$u(r,\phi)$ 应当满足连接条件

$$
\left.u(r,\phi)\right|_{r=a-0}
=
\left.u(r,\phi)\right|_{r=a+0},
$$

$$
\left.
\frac{\partial u(r,\phi)}{\partial r}
\right|_{r=a-0}
=
\left.
\frac{\partial u(r,\phi)}{\partial r}
\right|_{r=a+0},
$$

即

$$
1+
\sum_{n=1}^{\infty}
C_n
\sin\frac{n\pi}{2\alpha}(\phi+\alpha)
=
\sum_{n=1}^{\infty}
D_n
\sin\frac{n\pi}{2\alpha}(\phi+\alpha),
$$

$$
\sum_{n=1}^{\infty}
nC_n
\sin\frac{n\pi}{2\alpha}(\phi+\alpha)
=
-
\sum_{n=1}^{\infty}
nD_n
\sin\frac{n\pi}{2\alpha}(\phi+\alpha).
$$

根据本征函数的正交性

$$
\int_{-\alpha}^{\alpha}
\sin\frac{n\pi}{2\alpha}(\phi+\alpha)
\sin\frac{m\pi}{2\alpha}(\phi+\alpha)
\,d\phi
=
\alpha\delta_{nm},
$$

就能得到

$$
\frac2{n\pi}
[1-(-1)^n]
+
C_n
=
D_n,
\qquad
C_n=-D_n,
$$

从而定出

$$
C_n
=
-\frac1{n\pi}
[1-(-1)^n],
\qquad
D_n
=
\frac1{n\pi}
[1-(-1)^n].
$$

最后就求得

$$
u(r,\phi)
=
\begin{cases}
\displaystyle
1-
\frac2\pi
\sum_{n=0}^{\infty}
\frac{(-1)^n}{2n+1}
\left(
\frac ra
\right)^{(2n+1)\pi/2\alpha}
\cos\frac{2n+1}{2\alpha}\pi\phi,
& r<a,\\[16pt]
\displaystyle
\frac2\pi
\sum_{n=0}^{\infty}
\frac{(-1)^n}{2n+1}
\left(
\frac ar
\right)^{(2n+1)\pi/2\alpha}
\cos\frac{2n+1}{2\alpha}\pi\phi,
& r>a.
\end{cases}
$$


## §14.5 Helmholtz 方程在柱坐标系下的分离变量

在柱坐标系中，Helmholtz 方程

$$
\nabla^2u+k^2u=0
\tag{14.87}
$$

的具体形式是

$$
\frac1\rho
\frac{\partial}{\partial\rho}
\left(
\rho\frac{\partial u}{\partial\rho}
\right)
+
\frac1{\rho^2}
\frac{\partial^2u}{\partial\phi^2}
+
\frac{\partial^2u}{\partial z^2}
+
k^2u
=
0.
\tag{14.88}
$$

这里 $u$ 是三个自变量的函数。要将它分离变量，应采取逐个分离的办法：先分离一个自变量，然后再将其余两个自变量分离。为此，令 $u(\rho,\phi,z)=v(\rho,\phi)Z(z)$，代入方程，即得

$$
\frac1v
\left[
\frac1\rho
\frac{\partial}{\partial\rho}
\left(
\rho\frac{\partial v}{\partial\rho}
\right)
+
\frac1{\rho^2}
\frac{\partial^2v}{\partial\phi^2}
+
k^2v
\right]
=
-\frac1Z\frac{d^2Z}{dz^2}.
$$

等式的左端是 $\rho$ 和 $\phi$ 的函数，与 $z$ 无关；右端是 $z$ 的函数，与 $\rho$ 及 $\phi$ 均无关。所以它们必须等于既与 $\rho,\phi$ 无关又与 $z$ 无关的常数。把这个常数记为 $\lambda$，就得到

$$
\frac1\rho
\frac{\partial}{\partial\rho}
\left(
\rho\frac{\partial v}{\partial\rho}
\right)
+
\frac1{\rho^2}
\frac{\partial^2v}{\partial\phi^2}
+
(k^2-\lambda)v
=
0,
\tag{14.89}
$$

$$
\frac{d^2Z}{dz^2}
+
\lambda Z
=
0.
\tag{14.90}
$$

这样就完成了自变量 $z$ 的分离。接着再令 $v(\rho,\phi)=R(\rho)\Phi(\phi)$，代入方程 (14.89)，又得到

$$
\frac{\rho^2}{R(\rho)}
\left[
\frac1\rho
\frac d{d\rho}
\left(
\rho\frac{dR}{d\rho}
\right)
+
(k^2-\lambda)R
\right]
=
-\frac1{\Phi(\phi)}
\frac{d^2\Phi}{d\phi^2}.
$$

现在，再次看到，等式的左端只是 $\rho$ 的函数，与 $\phi$ 无关；右端只是 $\phi$ 的函数，与 $\rho$ 无关。所以它们必须等于既与 $\rho$ 无关又与 $\phi$ 无关的常数，记为 $\mu$。于是又得到

$$
\frac1\rho
\frac d{d\rho}
\left(
\rho\frac{dR}{d\rho}
\right)
+
\left(
k^2-\lambda-\frac{\mu}{\rho^2}
\right)R
=
0,
\tag{14.91}
$$

$$
\frac{d^2\Phi}{d\phi^2}
+
\mu\Phi
=
0.
\tag{14.92}
$$

这样，就完成了对于自变量 $\rho$ 和 $\phi$ 的分离变量，也就完成了 Helmholtz 方程的分离变量。

在柱坐标系中将 Helmholtz 方程分离变量时，得到了三个常微分方程，即 (14.90)、(14.91) 和 (14.92)。(14.90) 和 (14.92) 是常系数常微分方程。方程 (14.91) 是变系数常微分方程，当 $k^2-\lambda\ne0$ 时可化为 Bessel 方程，在第九章中曾经讨论 Bessel 方程的求解问题。后面第十六章中还要详细讨论 Bessel 方程的解和解的性质。

## §14.6 Helmholtz 方程在球坐标系下的分离变量

在球坐标系中，Helmholtz 方程的具体形式是

$$
\frac1{r^2}
\frac{\partial}{\partial r}
\left(
r^2\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2\sin\theta}
\frac{\partial}{\partial\theta}
\left(
\sin\theta\frac{\partial u}{\partial\theta}
\right)
+
\frac1{r^2\sin^2\theta}
\frac{\partial^2u}{\partial\phi^2}
+
k^2u
=
0.
\tag{14.93}
$$

令 $u(r,\theta,\phi)=R(r)S(\theta,\phi)$，代入方程，即得

$$
\frac{r^2}{R(r)}
\left[
\frac1{r^2}
\frac d{dr}
\left(
r^2\frac{dR(r)}{dr}
\right)
+
k^2R(r)
\right]
=
-\frac1{S(\theta,\phi)}
\left\{
\frac1{\sin\theta}
\frac{\partial}{\partial\theta}
\left[
\sin\theta
\frac{\partial S(\theta,\phi)}{\partial\theta}
\right]
+
\frac1{\sin^2\theta}
\frac{\partial^2S(\theta,\phi)}{\partial\phi^2}
\right\}.
$$

等式的左端只是 $r$ 的函数，与 $\theta,\phi$ 无关；右端只是 $\theta,\phi$ 的函数，与 $r$ 无关。所以它们必须等于既与 $r$ 无关又与 $\theta,\phi$ 无关的常数。把这个常数记为 $\lambda$，就成功地完成了从 Helmholtz 方程中分离去径向部分（即与 $r$ 有关的部分）的任务，

$$
\frac1{r^2}
\frac d{dr}
\left(
r^2\frac{dR(r)}{dr}
\right)
+
\left(
k^2-\frac{\lambda}{r^2}
\right)
R(r)
=
0,
\tag{14.94}
$$

$$
\frac1{\sin\theta}
\frac{\partial}{\partial\theta}
\left[
\sin\theta
\frac{\partial S(\theta,\phi)}{\partial\theta}
\right]
+
\frac1{\sin^2\theta}
\frac{\partial^2S(\theta,\phi)}{\partial\phi^2}
+
\lambda S(\theta,\phi)
=
0.
\tag{14.95}
$$

继续将方程 (14.95) 分离变量。为此，再令 $S(\theta,\phi)=\Theta(\theta)\Phi(\phi)$，代入 (14.95)，又得到

$$
\frac{\sin^2\theta}{\Theta}
\left[
\frac1{\sin\theta}
\frac d{d\theta}
\left(
\sin\theta\frac{d\Theta(\theta)}{d\theta}
\right)
+
\lambda\Theta
\right]
=
-\frac1\Phi
\frac{d^2\Phi}{d\phi^2}.
$$

这样，我们再次看到，等式的左端只是 $\theta$ 的函数，与 $\phi$ 无关；右端只是 $\phi$ 的函数，与 $\theta$ 无关。所以它们必须等于既与 $\theta$ 无关又与 $\phi$ 无关的常数，记为 $\mu$，于是又完成了将 $\theta$ 部分和 $\phi$ 部分的分离，得到的两个常微分方程是

$$
\frac1{\sin\theta}
\frac d{d\theta}
\left[
\sin\theta
\frac{d\Theta(\theta)}{d\theta}
\right]
+
\left(
\lambda-\frac{\mu}{\sin^2\theta}
\right)
\Theta
=
0,
\tag{14.96}
$$

$$
\frac{d^2\Phi}{d\phi^2}
+
\mu\Phi
=
0.
\tag{14.97}
$$

这样又最终完成了 Helmholtz 方程在球坐标系中的分离变量。

在球坐标系中将 Helmholtz 方程分离变量，当然也还是得到三个常微分方程，即 (14.94)、(14.96) 和 (14.97)。其中的 (14.97) 是常系数常微分方程，不难求解。其余两个方程都是变系数常微分方程。方程 (14.96) 称为连带 Legendre 方程。当 $k\ne0$ 时，方程 (14.94) 可化为球 Bessel 方程。在第十五章和第十六章中也将要分别讨论它们。

这里还要讨论一种常见特殊情形，即 $u=u(r,\theta)$ 与 $\phi$ 无关的情形。这就是说，整个定解问题在绕极轴转动任意角时不变。在这种情形下，Helmholtz 方程的形式就化简为

$$
\frac1{r^2}
\frac{\partial}{\partial r}
\left(
r^2\frac{\partial u}{\partial r}
\right)
+
\frac1{r^2\sin\theta}
\frac{\partial}{\partial\theta}
\left(
\sin\theta\frac{\partial u}{\partial\theta}
\right)
+
k^2u
=
0.
\tag{14.98}
$$

令 $u(r,\theta)=R(r)\Theta(\theta)$，代入方程，即得

$$
\frac{r^2}{R(r)}
\left\{
\frac1{r^2}
\frac d{dr}
\left[
r^2\frac{dR(r)}{dr}
\right]
+
k^2R(r)
\right\}
=
-\frac1{\Theta(\theta)}
\frac1{\sin\theta}
\frac d{d\theta}
\left(
\sin\theta
\frac{d\Theta(\theta)}{d\theta}
\right).
$$

等式的左端只是 $r$ 的函数，与 $\theta$ 无关；右端只是 $\theta$ 的函数，与 $r$ 无关。所以它们必须等于既与 $r$ 无关又与 $\theta$ 无关的常数，记作 $\lambda$，这样就完成了分离变量。得到的两个常微分方程，径向方程和 (14.94) 完全相同，另一个常微分方程

$$
\frac1{\sin\theta}
\frac d{d\theta}
\left[
\sin\theta
\frac{d\Theta(\theta)}{d\theta}
\right]
+
\lambda\Theta(\theta)
=
0
\tag{14.99}
$$

仅与方位角 $\theta$ 有关，称为 Legendre 方程，是连带 Legendre 方程 (14.96) 的特殊情形 $(\mu=0)$。

## *§14.7 矢量波动方程和矢量 Helmholtz 方程

前面讨论的波动方程

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\nabla^2u
=
0
$$

中，未知函数代表的物理量都是标量。但在某些物理问题中，未知函数代表的物理量 $\boldsymbol u(\boldsymbol r,t)$ 是矢量，因而可能出现矢量波动方程

$$
\frac{\partial^2\boldsymbol u}{\partial t^2}
-
a^2\nabla^2\boldsymbol u
=
0.
\tag{14.100}
$$

例如，在均匀各向同性介质中，电磁场的 Maxwell 方程组为

$$
\nabla\cdot\boldsymbol E
=
\frac{\rho}{\varepsilon},
\qquad
\nabla\times\boldsymbol E
=
-\mu
\frac{\partial\boldsymbol H}{\partial t},
$$

$$
\nabla\cdot\boldsymbol H
=
0,
\qquad
\nabla\times\boldsymbol H
=
\varepsilon
\frac{\partial\boldsymbol E}{\partial t}
+
\boldsymbol j.
\tag{14.101}
$$

由此即可导出

$$
\frac{\partial^2\boldsymbol H}{\partial t^2}
=
-\frac1\mu
\nabla\times
\left(
\frac{\partial\boldsymbol E}{\partial t}
\right)
=
-\frac1{\mu\varepsilon}
\left[
\nabla\times(\nabla\times\boldsymbol H)
-
\nabla\times\boldsymbol j
\right],
$$

$$
\frac{\partial^2\boldsymbol E}{\partial t^2}
=
\frac1\varepsilon
\left[
\nabla\times
\left(
\frac{\partial\boldsymbol H}{\partial t}
\right)
-
\frac{\partial\boldsymbol j}{\partial t}
\right]
=
-\frac1{\mu\varepsilon}
\nabla\times(\nabla\times\boldsymbol E)
-
\frac1\varepsilon
\frac{\partial\boldsymbol j}{\partial t}.
$$

这样，在 $\rho=0,\ \boldsymbol j=0$ 的区域内，$\boldsymbol E$ 和 $\boldsymbol H$ 就都满足齐次矢量波动方程

$$
\frac{\partial^2\boldsymbol E}{\partial t^2}
-
c^2\nabla^2\boldsymbol E
=
0,
\qquad
\frac{\partial^2\boldsymbol H}{\partial t^2}
-
c^2\nabla^2\boldsymbol H
=
0,
\tag{14.102}
$$

其中 $c=1/\sqrt{\varepsilon\mu}$。

将矢量波动方程 (14.100) 分离变量，$\boldsymbol u(\boldsymbol r,t)=\boldsymbol A(\boldsymbol r)T(t)$，可以进一步得到矢量 $\boldsymbol A$ 满足的方程

$$
\nabla^2\boldsymbol A
+
k^2\boldsymbol A
=
0,
\tag{14.103}
$$

称为矢量 Helmholtz 方程。由 (14.53) 式可知，在直角坐标系下矢量 Helmholtz 方程 (14.103) 就完全等价于三个独立的标量 Helmholtz 方程

$$
\nabla^2A_x+k^2A_x=0,
\qquad
\nabla^2A_y+k^2A_y=0,
\qquad
\nabla^2A_z+k^2A_z=0.
\tag{14.104}
$$

这种做法可以简称为矢量 Helmholtz 方程的直角坐标分解。当然也可以作柱坐标（包括二维情形下的平面极坐标）分解或球坐标分解。

### (1) 矢量 Helmholtz 方程 (14.103) 的柱坐标分解

利用 (14.54) 式，柱坐标下矢量 Helmholtz 方程 (14.103) 就分解为关于 $A_\rho$ 和 $A_\phi$ 的偏微分方程组

$$
\nabla^2A_\rho
-
\frac1{\rho^2}A_\rho
-
\frac2{\rho^2}
\frac{\partial A_\phi}{\partial\phi}
+
k^2A_\rho
=
0,
\tag{14.105a}
$$

$$
\nabla^2A_\phi
-
\frac1{\rho^2}A_\phi
+
\frac2{\rho^2}
\frac{\partial A_\rho}{\partial\phi}
+
k^2A_\phi
=
0,
\tag{14.105b}
$$

以及关于 $A_z$ 的 Helmholtz 方程

$$
\nabla^2A_z+k^2A_z=0.
\tag{14.105c}
$$

### (2) 矢量 Helmholtz 方程 (14.103) 的球坐标分解

利用 (14.55) 式，球坐标下矢量 Helmholtz 方程 (14.103) 就分解为关于 $A_r,A_\theta$ 和 $A_\phi$ 的偏微分方程组

$$
\nabla^2A_r
-
\frac2{r^2}A_r
-
\frac2{r^2\sin\theta}
\frac{\partial(\sin\theta A_\theta)}{\partial\theta}
-
\frac2{r^2\sin\theta}
\frac{\partial A_\phi}{\partial\phi}
+
k^2A_r
=
0,
\tag{14.106a}
$$

$$
\nabla^2A_\theta
-
\frac1{r^2\sin^2\theta}A_\theta
+
\frac2{r^2}
\frac{\partial A_r}{\partial\theta}
-
\frac{2\cos\theta}{r^2\sin^2\theta}
\frac{\partial A_\phi}{\partial\phi}
+
k^2A_\theta
=
0,
\tag{14.106b}
$$

$$
\nabla^2A_\phi
-
\frac1{r^2\sin^2\theta}A_\phi
+
\frac2{r^2\sin\theta}
\frac{\partial A_r}{\partial\phi}
+
\frac{2\cos\theta}{r^2\sin^2\theta}
\frac{\partial A_\theta}{\partial\phi}
+
k^2A_\phi
=
0.
\tag{14.106c}
$$

从以上的讨论可以看出，矢量 Helmholtz 方程的各种分解式中，直角坐标分解最简单。这时仍可采用适当的坐标系，求解三个数学形式完全相同的标量 Helmholtz 方程，而后叠加出所要求的矢量解，包括矢量解的柱坐标或球坐标分量。求解矢量 Helmholtz 方程的另一种有效办法①，需要少许矢量运算。不难验证，如果标量函数 $w$ 是 Helmholtz 方程

$$
\nabla^2w+k^2w=0
$$

的解，则下面三个矢量

$$
\boldsymbol L=\nabla w,
\qquad
\boldsymbol M=\nabla\times(\boldsymbol a w),
\qquad
\boldsymbol N=\frac1k\nabla\times\boldsymbol M
$$

（其中 $\boldsymbol a$ 是任意的单位常矢量，$|\boldsymbol a|=1$）都满足矢量 Helmholtz 方程

$$
\nabla^2\boldsymbol A+k^2\boldsymbol A=0,
$$

而且它们是线性无关的。这种解法的优点除了只需求解标量 Helmholtz 方程外，还在于解 $\boldsymbol L,\boldsymbol M$ 和 $\boldsymbol N$ 具有明显的物理特性，即

$$
\nabla\times\boldsymbol L=0,
\qquad
\nabla\cdot\boldsymbol M=0,
\qquad
\nabla\cdot\boldsymbol N=0.
$$

这表明 $\boldsymbol L$ 是一个无旋场，而 $\boldsymbol M$ 和 $\boldsymbol N$ 是无源场，当然可以有选择地用来描写电磁场。

① 详细内容，请参阅有关著作，例如 J. A. Stratton, *Electromagnetic Theory*, Ch. 7, McGraw-Hill, New York, 1941; P. M. Morse & H. Feshbach, *Methods of Theoretical Physics*, Ch. 13, McGraw-Hill, New York, 1953.

## 习题

1. 一个半径为 $a$ 的无穷长空心导体圆柱，沿柱轴分成两半，互相绝缘。一半电势为 $V$，另一半电势为 $-V$（见图 14.1）。求柱内的电势分布。

![Figure 14.1](...)

图 14.1

2. 求在环形区域 $a\le r\le b$ 内满足边界条件

$$
\left.u\right|_{r=a}
=
f(\phi),
\qquad
\left.u\right|_{r=b}
=
g(\phi)
$$

的调和函数。

3. 半径为 $a$、表面熏黑的均匀金属圆柱，平放在地上，受到阳光照射，在垂直于光线的单位面积上单位时间内吸收热量为 $M$，同时，柱面按 Newton 冷却定律向外散热。试求柱内的稳定温度分布。取外界温度为 0，并设圆柱为无穷长。

4. 在圆域 $0\le x^2+y^2\le a^2$ 上求解：

(1)

$$
\begin{cases}
\nabla^2u=-4,\\
\left.u\right|_{x^2+y^2=a^2}=0;
\end{cases}
$$

(2)

$$
\begin{cases}
\nabla^2u=-4y,\\
\left.u\right|_{x^2+y^2=a^2}=0;
\end{cases}
$$

(3)

$$
\begin{cases}
\nabla^2u=-4xy,\\
\left.u\right|_{x^2+y^2=a^2}=0;
\end{cases}
$$

(4)

$$
\begin{cases}
\nabla^2u=-4(x+y),\\
\left.u\right|_{x^2+y^2=a^2}=0.
\end{cases}
$$

5. 一个由理想导体做成的无穷长波导管，其截面均匀，如图 14.2 所示。管内为真空。

![Figure 14.2](...)

图 14.2

假定一个平面（即图中的一条直边）的电势为 $V$，其余面上的电势均为 0。试求波导管内的电势分布。

6. 求解球内的定解问题

$$
\frac{\partial u}{\partial t}
-
\frac{\kappa}{r^2}
\frac{\partial}{\partial r}
\left(
r^2\frac{\partial u}{\partial r}
\right)
=
0,
$$

$$
\left.u\right|_{r=0}
\text{ 有界},
\qquad
\left.u\right|_{r=1}
=
A\exp\{-(p\pi)^2\kappa t\},
$$

$$
\left.u\right|_{t=0}
=
0.
$$

提示：

$$
\frac1{r^2}
\frac{\partial}{\partial r}
\left(
r^2\frac{\partial u}{\partial r}
\right)
=
\frac1r
\frac{\partial^2(ru)}{\partial r^2}.
$$
