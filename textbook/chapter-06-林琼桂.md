# 第六章　数学物理方程的导出*

## 目录

- §1 简介　113
  - §1.1 波动方程　113
  - §1.2 热传导方程和扩散方程　114
  - §1.3 稳定场方程　114
- §2 *二阶线性偏微分方程的分类与化简　114
  - §2.1 二阶线性偏微分方程的分类　114
  - §2.2 双曲型方程　116
  - §2.3 椭圆型方程　117
  - §2.4 抛物型方程　117
- §3 波动方程　118
  - §3.1 弦的横振动　118
  - §3.2 杆的纵振动　120
  - §3.3 定解条件　121
  - §3.4 定解问题的适定性　123
- §4 热传导方程和扩散方程　124
  - §4.1 热传导方程　124
  - §4.2 扩散方程　126
  - §4.3 定解条件　127
- §5 稳定场方程　128
  - §5.1 稳定温度分布和稳定浓度分布　128
  - §5.2 静电场方程　129
  - §5.3 定解条件　129
- 补充习题　130

## §1　简介

本书的第一部分内容，即复变函数部分至上章已经结束。读者已经看到，复变函数是微积分的内容从实变量到复变量的扩展。它是进一步学习数学物理方程和特殊函数的基础。（由于本书只研究数学物理方程和特殊函数的基本内容，所以较少使用复变函数的知识。但复变函数对于深入研究数学物理方程和特殊函数是非常重要的。）复变函数本身在物理、力学和工程问题中也具有多方面的应用。

本章开始研究本课程的第二部分，即数学物理方程。数学物理方程研究物理或工程问题中所涉及的各种偏微分方程、积分方程和微分—积分方程。本课程只研究偏微分方程，而且只限于二阶线性偏微分方程。内容包括方程的导出、求解和对解的物理分析。即在一定的条件下将物理问题简化（忽略一些次要因素），应用物理学的基本规律将它翻译成数学问题；求解所得的数学问题（这部分是重点）；之后，分析所得解的物理图象、意义和适用范围等。它所提供的方法在经典物理、近代物理和工程技术中都有广泛的应用。这部分内容对于读者学好后续的四大力学课程以及有关的研究生基础课程，乃至于将来从事研究工作都有重要的作用。因为它不仅提供了具体的数学方法，也培养了思维方式和计算能力。这部分内容有一定的难度，主要是计算比较复杂。但读者只要多思考，多动手计算，就能够逐步适应，渐入佳境。

### §1.1　波动方程

含有未知函数及其偏导数的方程称为偏微分方程（partial differential equation）。如果方程只包含未知函数及其偏导数的一次项，则称为线性偏微分方程。方程中出现的偏导数的最高阶数称为方程的阶。本课程只研究二阶线性偏微分方程。以两个自变量的情况为例，记自变量为 $x、y$，未知函数为 $u(x,y)$，则二阶线性偏微分方程的一般形式为

$$
a_{11}\frac{\partial^2u}{\partial x^2}+2a_{12}\frac{\partial^2u}{\partial x\partial y}+a_{22}\frac{\partial^2u}{\partial y^2}+b_1\frac{\partial u}{\partial x}+b_2\frac{\partial u}{\partial y}+cu+f=0,\tag{1}
$$

其中 $a_{11}$ 等系数和 $f$ 都可以是 $x、y$ 的函数。如果 $f=0$，称为齐次方程，否则称为非齐次方程。上式类似于二次曲线的一般形式。通过坐标的旋转和平移，可以将二次曲线化为标准形式。与此类似，通过适当的自变量变换和函数变换，也可以将二阶线性偏微分方程简化为标准形式。二次曲线有三种不同类型，即双曲线、抛物线和椭圆。与此类似，二阶线性偏微分方程也有三种不同类型，也称为双曲型、抛物型和椭圆型。曲线或方程的类型不因变换而改变，但化为标准形式后，各类型的特征可以一目了然。更重要的是，如果我们能够求解标准形式，则一般形式的解就可以通过适当的变换得到。关于分类和化简的具体细节，读者可参看下节的讨论。本节只介绍物理上常见的三类方程的标准形式，这比数学上所指的标准形式还要简单。

在三维空间，波动方程的标准形式为

$$
\frac{\partial^2u}{\partial t^2}-a^2\nabla^2u=f,\tag{2a}
$$

其中 $u=u(\boldsymbol r,t)$ 是我们所要研究的物理量，比如位移或电场（磁场）的一个分量，$\boldsymbol r$ 是空间位置，$t$ 是时间，$a$ 是常数，表示波动的传播速度，$f=f(\boldsymbol r,t)$ 是已知函数，代表波动的源，$\nabla^2$ 是 Laplace 算符，有些书上记作 $\Delta$，在直角坐标系中，

$$
\nabla^2u=\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}+\frac{\partial^2u}{\partial z^2}.
$$

但读者应当注意，这并不是 Laplace 算符的定义。Laplace 算符的定义是 $\nabla^2u=\nabla\cdot\nabla u$，这具有确定的几何意义，而不依赖于坐标系。在一维空间，波动方程的标准形式为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f,\tag{2b}
$$

这可以描写弦的横振动或弹性细杆的纵振动。在数学上，波动方程属于双曲型。

### §1.2　热传导方程和扩散方程

热传导方程和扩散方程统称为输运方程。在三维空间，它们的标准形式为

$$
\frac{\partial u}{\partial t}-a^2\nabla^2u=f,\tag{3}
$$

其中 $u=u(\boldsymbol r,t)$ 是我们所要研究的物理量，对于热传导方程，它表示温度的分布，对于扩散方程，它表示杂质浓度的分布。$f=f(\boldsymbol r,t)$ 是已知函数，代表热源或杂质源。在数学上，输运方程属于抛物型。

### §1.3　稳定场方程

稳定场方程的标准形式为

$$
\nabla^2u=f,\tag{4}
$$

其中 $u=u(\boldsymbol r)$ 可以代表稳定温度分布、稳定浓度分布或静电场的电势，相应地，$f=f(\boldsymbol r)$ 代表不随时间变化的热源、杂质源或自由电荷分布。上式又称为 Poisson 方程，当 $f=0$ 时又称为 Laplace 方程。在数学上，稳定场方程属于椭圆型。

## §2　*二阶线性偏微分方程的分类与化简

### §2.1　二阶线性偏微分方程的分类

我们只讨论具有两个自变量的二阶线性偏微分方程，其一般形式已经在式 (1) 给出。显然，系数 $a_{11}、a_{12}、a_{22}$ 不能全都为零，否则方程就不是二阶的了。如果 $a_{11}$ 和 $a_{22}$ 同时为零，则方程已经属于标准形式（参看以下的进一步讨论）。因此，不妨假定 $a_{11}\ne0$。在下面的讨论中，记 $x_1=x、x_2=y$，并简记 $\partial u/\partial x_i$ 为 $\partial_i u$，这样，式 (1) 可以写作

$$
\sum_{i,j=1}^2a_{ij}\partial_i\partial_j u+\sum_{i=1}^2b_i\partial_i u+cu+f=0,\tag{5}
$$

其中 $a_{21}=a_{12}$。这是近代物理的文献中常见的记号，它可以使书写得以简化。在下面的讨论中，如果读者一时看不清式子的意义，可以在草稿纸上写出其原来的形式加以比较，慢慢就熟悉了。今作自变量变换

$$
\xi_i=\xi_i(x_1,x_2),\qquad i=1,2.\tag{6}
$$

记变换的 Jacobi 矩阵为

$$
D=\begin{pmatrix}
\partial_1\xi_1&\partial_1\xi_2\\
\partial_2\xi_1&\partial_2\xi_2
\end{pmatrix}.\tag{7}
$$

我们要求 $\det D\ne0$，以便变换是可逆的。将 $\partial u/\partial\xi_i$ 简记为 $\partial'_i u$，利用偏导数的求导法则可得

$$
\partial_k u=\sum_{i=1}^2(\partial_k\xi_i)\partial'_i u
$$

和

$$
\partial_k\partial_lu=\sum_{i=1}^2(\partial_k\partial_l\xi_i)\partial'_i u+\sum_{i,j=1}^2(\partial_k\xi_i)(\partial_l\xi_j)\partial'_i\partial'_j u.
$$

将式 (5) 中的求和指标 $i$ 换为 $k$，$j$ 换为 $l$，然后将上面两式代入，略作整理，就得到关于新变量的偏微分方程

$$
\sum_{i,j=1}^2\alpha_{ij}\partial'_i\partial'_j u+\sum_{i=1}^2\beta_i\partial'_i u+cu+f=0,\tag{8}
$$

它仍然具有原来的形式，只是其中的系数变为

$$
\alpha_{ij}=\sum_{k,l=1}^2a_{kl}(\partial_k\xi_i)(\partial_l\xi_j),\qquad
\beta_i=\sum_{k=1}^2b_k(\partial_k\xi_i)+\sum_{k,l=1}^2a_{kl}(\partial_k\partial_l\xi_i).\tag{9}
$$

上式的第一式可以写作

$$
\alpha_{ij}=\sum_{k,l=1}^2(\partial_k\xi_i)a_{kl}(\partial_l\xi_j)=\sum_{k,l=1}^2D_{ki}a_{kl}D_{lj}.
$$

写成矩阵形式，就是

$$
\alpha=D^\tau aD,\tag{10}
$$

其中 $\tau$ 表示转置。于是有 $\det\alpha=(\det D)^2\det a$，再用分量写出来，就是

$$
\alpha_{12}^2-\alpha_{11}\alpha_{22}=(\det D)^2(a_{12}^2-a_{11}a_{22}).\tag{11}
$$

可见 $a_{12}^2-a_{11}a_{22}$ 的符号（或为零）不随变换而改变，因而可以用来对方程进行分类。如果 $a_{12}^2-a_{11}a_{22}>0$，则称方程为双曲型（hyperbolic）；如果 $a_{12}^2-a_{11}a_{22}<0$，则称方程为椭圆型（elliptic）；如果 $a_{12}^2-a_{11}a_{22}=0$，则称方程为抛物型（parabolic）。一般来说，$a_{11}$ 等系数是 $x、y$ 的函数，所以，在一定的区域内，$a_{12}^2-a_{11}a_{22}$ 不一定能保持固定的符号或恒为零。如果在所考虑的区域的不同部分，方程属于不同的类型，则称方程在该区域为混合型。

接下来的任务就是对方程进行化简。化简的思路是选取适当的变换使得 $\alpha_{11}$ 或 $\alpha_{22}$ 为零，或两者同时为零（如果有可能）。注意到 $\alpha_{11}$ 和 $\alpha_{22}$ 的形式类似，所以关键是考察下列方程的解

$$
\sum_{k,l=1}^2a_{kl}\partial_k\varphi\partial_l\varphi=a_{11}(\partial_x\varphi)^2+2a_{12}(\partial_x\varphi)(\partial_y\varphi)+a_{22}(\partial_y\varphi)^2=0,\tag{12}
$$

其中 $\partial_x\varphi$ 即 $\partial_1\varphi$，$\partial_y\varphi$ 即 $\partial_2\varphi$。如果能找到上式的两个函数无关解 $\varphi_1(x,y)$ 和 $\varphi_2(x,y)$，令 $\xi_1(x,y)=\varphi_1(x,y)$，$\xi_2(x,y)=\varphi_2(x,y)$，则可得 $\alpha_{11}=\alpha_{22}=0$。即使只能找到一个解，也可以使 $\alpha_{11}$ 和 $\alpha_{22}$ 中的一个为零。当然，化简是在方程具有确定类型的区域内进行的。

为了寻找方程 (12) 的解，我们介绍一个

**定理**　如果 $\varphi(x,y)$ 是方程 (12) 的解，则由 $\varphi(x,y)=c$ 确定的隐函数满足下列常微分方程

$$
a_{11}(dy)^2-2a_{12}dxdy+a_{22}(dx)^2=0.\tag{13}
$$

反之亦然。

这个定理的证明并不难。这里就从略了。不过，应该指出，在曲线 $\varphi(x,y)=c$ 上各点，$\partial_x\varphi$ 和 $\partial_y\varphi$ 不能同时为零，以便在曲线上各点，$dy/dx$ 和 $dx/dy$ 当中至少有一个是有意义的。

式 (13) 称为式 (1) 的特征方程，其积分曲线 $\varphi(x,y)=c$ 称为式 (1) 的特征线。由于假定 $a_{11}\ne0$，由式 (13) 可以解出

$$
\frac{dy}{dx}=\frac{a_{12}\pm\sqrt{a_{12}^2-a_{11}a_{22}}}{a_{11}}.\tag{14}
$$

下面对各种类型分别进行讨论。

### §2.2　双曲型方程

此时 $a_{12}^2-a_{11}a_{22}>0$，式 (14) 有两族积分曲线，$\varphi_1(x,y)=c_1$（对应于 $+$ 号）和 $\varphi_2(x,y)=c_2$（对应于 $-$ 号）。容易得到

$$
a_{11}\partial_x\varphi_1=-\left(a_{12}+\sqrt{a_{12}^2-a_{11}a_{22}}\right)\partial_y\varphi_1,
\qquad
a_{11}\partial_x\varphi_2=-\left(a_{12}-\sqrt{a_{12}^2-a_{11}a_{22}}\right)\partial_y\varphi_2.\tag{15}
$$

由此可得

$$
\frac{\partial(\varphi_1,\varphi_2)}{\partial(x,y)}=-2\frac{\sqrt{a_{12}^2-a_{11}a_{22}}}{a_{11}}(\partial_y\varphi_1)(\partial_y\varphi_2).
$$

但 $\partial_y\varphi_1\ne0$，否则由式 (15) 将得到 $\partial_x\varphi_1=0$，这与 $\partial_x\varphi_1$ 和 $\partial_y\varphi_1$ 不同时为零的要求矛盾；同理 $\partial_y\varphi_2\ne0$。于是上式左边的 Jacobi 行列式不为零，故 $\varphi_1(x,y)$ 和 $\varphi_2(x,y)$ 是函数独立的。作变换

$$
\xi_1(x,y)=\varphi_1(x,y),\qquad \xi_2(x,y)=\varphi_2(x,y).\tag{16}
$$

上述讨论保证了变换是可逆的。易知该变换导致 $\alpha_{11}=\alpha_{22}=0$。而由式 (11) 可知 $\alpha_{12}\ne0$。所以现在式 (8) 的形式为

$$
\partial'_1\partial'_2u+\sum_{i=1}^2\widetilde\beta_i\partial'_iu+\widetilde cu+\widetilde f=0,\tag{17}
$$

其中 $\widetilde\beta_i=\beta_i/2\alpha_{12}$，$\widetilde c=c/2\alpha_{12}$，$\widetilde f=f/2\alpha_{12}$。就两个自变量的方程来说，上式已经是双曲型方程的标准形式。我们还可以再作变换

$$
t=\xi_1+\xi_2,\qquad s=\xi_1-\xi_2,\tag{18}
$$

将式 (17) 化为

$$
\partial_t^2u-\partial_s^2u+(\widetilde\beta_1+\widetilde\beta_2)\partial_tu+(\widetilde\beta_1-\widetilde\beta_2)\partial_su+\widetilde cu+\widetilde f=0.\tag{19}
$$

上节所介绍的一维波动方程就具有这样的形式。多个自变量的双曲型方程的标准形式也与上式类似，一个典型的例子就是上节给出的三维波动方程。

### §2.3　椭圆型方程

此时 $a_{12}^2-a_{11}a_{22}<0$，式 (14) 化为

$$
\frac{dy}{dx}=\frac{a_{12}\pm i\sqrt{a_{11}a_{22}-a_{12}^2}}{a_{11}}.\tag{20}
$$

其两组解为 $\varphi(x,y)=c_1$（对应于 $+$ 号）和 $\varphi^*(x,y)=c_2$（对应于 $-$ 号），其中 $\varphi(x,y)$ 是复值函数。为了避免引入复变量，作变换

$$
\xi_1(x,y)=\operatorname{Re}\varphi(x,y),\qquad \xi_2(x,y)=\operatorname{Im}\varphi(x,y).\tag{21}
$$

将 $\partial_x\varphi+(\partial_y\varphi)(dy/dx)=0$ 分开实部和虚部，得到

$$
a_{11}\partial_x\xi_1=-a_{12}\partial_y\xi_1+\sqrt{a_{11}a_{22}-a_{12}^2}\,\partial_y\xi_2,
\qquad
a_{11}\partial_x\xi_2=-a_{12}\partial_y\xi_2-\sqrt{a_{11}a_{22}-a_{12}^2}\,\partial_y\xi_1.\tag{22}
$$

由此可得

$$
\frac{\partial(\xi_1,\xi_2)}{\partial(x,y)}=\frac{\sqrt{a_{11}a_{22}-a_{12}^2}}{a_{11}}\big[(\partial_y\xi_1)^2+(\partial_y\xi_2)^2\big].
$$

左边的 Jacobi 行列式不应为零，否则将导致 $\partial_y\xi_1=\partial_y\xi_2=0$，再由式 (22) 就有 $\partial_x\xi_1=\partial_x\xi_2=0$，但这与 $\partial_x\varphi$ 和 $\partial_y\varphi$ 不同时为零的要求矛盾。因此 $\xi_1(x,y)$ 和 $\xi_2(x,y)$ 是函数独立的。从而变换 (21) 是可逆的。由于 $\varphi(x,y)$ 满足式 (12)，分开该式的实部和虚部，可得 $\alpha_{11}=\alpha_{22}$ 和 $\alpha_{12}=0$。而由式 (11) 可知 $\alpha_{11}\ne0$ 且为实数。所以现在式 (8) 的形式为

$$
\partial_1'^2u+\partial_2'^2u+\sum_{i=1}^2\widetilde\beta_i\partial'_iu+\widetilde cu+\widetilde f=0,\tag{23}
$$

其中 $\widetilde\beta_i=\beta_i/\alpha_{11}$，$\widetilde c=c/\alpha_{11}$，$\widetilde f=f/\alpha_{11}$。这就是椭圆型方程的标准形式。二维稳定场方程就具有这样的形式。多个自变量的椭圆型方程的标准形式也与上式类似。

### §2.4　抛物型方程

此时 $a_{12}^2-a_{11}a_{22}=0$，式 (14) 只有一族积分曲线，记作 $\varphi_2(x,y)=c$。任取与 $\varphi_2(x,y)$ 函数独立的 $\varphi_1(x,y)$，作变换

$$
\xi_1(x,y)=\varphi_1(x,y),\qquad \xi_2(x,y)=\varphi_2(x,y).\tag{24}
$$

易知上述变换导致 $\alpha_{22}=0$。而由式 (11) 可知 $\alpha_{12}=0$。此时 $\alpha_{11}$ 必不为零，否则系数变换公式 (10) 将导致所有 $a_{ij}$ 为零，这与 $a_{11}\ne0$ 的假设矛盾。所以现在式 (8) 的形式为

$$
\partial_1'^2u+\sum_{i=1}^2\widetilde\beta_i\partial'_iu+\widetilde cu+\widetilde f=0,\tag{25}
$$

其中 $\widetilde\beta_i=\beta_i/\alpha_{11}$，$\widetilde c=c/\alpha_{11}$，$\widetilde f=f/\alpha_{11}$。再作变换

$$
u(\xi_1,\xi_2)=v(\xi_1,\xi_2)\exp\left[-\frac12\int^{\xi_1}\widetilde\beta_1(\tau,\xi_2)\,d\tau\right],\tag{26}
$$

则式 (25) 进一步化简为

$$
\partial_1'^2v+\widetilde\beta_2\partial'_2v+dv+g=0,\tag{27}
$$

其中

$$
d=\widetilde c-\frac14\widetilde\beta_1^2-\frac12\partial'_1\widetilde\beta_1-\frac12\widetilde\beta_2\int^{\xi_1}\partial'_2\widetilde\beta_1(\tau,\xi_2)\,d\tau,
\qquad
g=f\exp\left[\frac12\int^{\xi_1}\widetilde\beta_1(\tau,\xi_2)\,d\tau\right].\tag{28}
$$

化简后的标准形式 (27) 不再含有对 $\xi_1$ 的一阶偏导数。上节所介绍的一维热传导方程就具有这样的形式。多个自变量的抛物型方程的标准形式也与式 (27) 类似，一个典型的例子就是上节给出的三维热传导方程。

## §3　波动方程

![Figure 1](...)

图 1：弦上一小段的位移和受力分析

### §3.1　弦的横振动

考虑一根柔软的弦，平衡时沿 $x$ 轴绷紧。当受到外界扰动时，它可以在垂直于 $x$ 轴的方向作横向振动。弦上的张力使得各部分互相牵制，从而振动可以在弦上传播而形成机械波。假定初始激励（即外界扰动）所引起的振动是小振动（其意义将在以下的讨论中逐步明确），并且弦上各点的位移方向平行，即振动发生在平面内，因而位移是一个标量函数（只有一个分量）。在这些前提下，我们可以用 Newton 运动定律来推导位移所满足的波动方程。记弦上 $x$ 点在 $t$ 时刻的位移为 $u(x,t)$。通常在振动平面内将 $x$ 轴绕原点逆时针旋转 $90^\circ$ 作为 $u$ 轴，这与观察的位置有关，所以，位移的正向的规定带有一定的任意性。如果振动发生在竖直平面内，一般取向上为正。波动方程就是 $u(x,t)$ 所满足的微分方程，一般来说是偏微分方程。为清楚起见，下面将推导过程分为几步。

1. 由于弦是柔软的，所以它不抗弯曲，放松时可以取任意形状。在张紧的状态下，不管是处于平衡状态或运动状态，弦上的张力必定沿切线方向。这是由柔软的假定所得出的结论。

2. 考虑运动引起的弦长的变化。平衡时位于区间 $[x,x+\Delta x]$ 的一小段，当运动时，其长度为

$$
\Delta s=\int_x^{x+\Delta x}\sqrt{(dx)^2+(du)^2}
=\int_x^{x+\Delta x}\sqrt{1+\left(\frac{\partial u}{\partial x}\right)^2}\,dx.
$$

由于所考虑的是小振动，故 $|\partial u/\partial x|\ll1$，其二次项可以忽略，故由上式得到

$$
\Delta s\approx\int_x^{x+\Delta x}dx=\Delta x.
$$

可见运动没有使得弦进一步伸长。由于张力取决于弦的伸长量，所以弦上的张力不随时间变化。这是由小振动的假定所得出的结论。

3. 考虑弦上平衡时位于区间 $[x,x+\Delta x]$ 的一小段在 $x$ 方向的运动方程。如图 1 所示，其中 $T_1$ 是 $x$ 处的张力，$T_2$ 是 $x+\Delta x$ 处的张力，图中对位移作了夸大。由于弦上各点都只有 $u$ 方向的位移，故 $x$ 方向的运动方程为

$$
T_2\cos\theta_2-T_1\cos\theta_1=0.
$$

注意这里用了 Newton 第二定律。由于所考虑的是小振动，故 $|\theta_1|\ll1$ 和 $|\theta_2|\ll1$，所以 $\cos\theta_1\approx1$（忽略了二阶小量），$\cos\theta_2\approx1$，于是由上式近似得到

$$
T_2=T_1\equiv T.\tag{29}
$$

这就是说，弦上的张力不随位置而变化。这也是由小振动的假定所得出的结论。综合以上两点可知，弦上的张力是一个常数。

4. 考虑上述小段在 $u$ 方向的运动方程。记弦的线密度（即单位长度的质量）为 $\rho$。由 Newton 第二定律，有

$$
\rho\Delta x\frac{\partial^2\bar u}{\partial t^2}=T\sin\theta_2-T\sin\theta_1,\tag{30}
$$

其中 $\bar u$ 是该小段的平均位移，它依赖于 $x、\Delta x$ 和 $t$。再一次利用小振动的假定，忽略三阶小量，就有

$$
\sin\theta_1\approx\tan\theta_1=\left.\frac{\partial u}{\partial x}\right|_x,
\qquad
\sin\theta_2\approx\tan\theta_2=\left.\frac{\partial u}{\partial x}\right|_{x+\Delta x}.
$$

注意一阶小量千万不能忽略。代入式 (30)，得到

$$
\rho\frac{\partial^2\bar u}{\partial t^2}=T\frac1{\Delta x}\left[\left.\frac{\partial u}{\partial x}\right|_{x+\Delta x}-\left.\frac{\partial u}{\partial x}\right|_x\right].
$$

令 $\Delta x\to0$，则 $\bar u\to u(x,t)$，上式成为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0,\tag{31}
$$

其中 $a=\sqrt{T/\rho}$。这正是 §1 所介绍的一维波动方程的标准形式。由于是自由振动，所以其中没有非齐次项。

5. 如果弦在振动过程中受到外力的作用，即强迫振动，则上面的运动方程需要修正。设 $x$ 处单位长度的受力（即力密度）为 $F(x,t)$，方向为 $u$ 轴的正方向，则式 (30) 应修改为

$$
\rho\Delta x\frac{\partial^2\bar u}{\partial t^2}=T\sin\theta_2-T\sin\theta_1+\bar F\Delta x,\tag{32}
$$

其中 $\bar F$ 是该小段的平均力密度，它依赖于 $x、\Delta x$ 和 $t$。当 $\Delta x\to0$，$\bar F\to F(x,t)$。重复上面的计算，得到外力作用下的波动方程为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f,\tag{33}
$$

其中 $f(x,t)=F(x,t)/\rho$ 是 $x$ 处单位质量的受力。这与 §1 给出的一维波动方程的标准形式 (2b) 完全一致。

上面的推导主要用到经典物理的 Newton 运动定律，所以并没有什么深奥之处。但读者应该特别注意其中对研究对象及其运动过程作了若干简化的假定，所以由上面得到的方程只是实际情况的比较粗糙的近似。如果在某些具体问题中，这些假定难以满足，比如振动的振幅较大，就需要重新考虑方程的推导。这时候问题显然要复杂得多。

![Figure 2](...)

图 2：杆上一小段。下图为其平衡位置，上图为运动时的位置和受力分析。

### §3.2　杆的纵振动

本小节研究弹性杆沿着杆长方向的小振动，即纵振动。初始激励（纵向的拉伸或挤压）或外力的作用（作用于内部或端点）都可以引起纵振动。取 $x$ 轴沿着杆长方向且向右。为简单起见，设杆的横截面面积和物理性质不随 $x$ 变化。记杆上 $x$ 点在 $t$ 时刻的位移为 $u(x,t)$，通常规定其正向与 $x$ 轴的方向一致，就目前的情况来说，即以离开平衡位置向右为正。下面推导 $u(x,t)$ 所满足的偏微分方程。

由于形变，杆上各处存在弹性力，记 $x$ 点的应力（单位面积的弹性力）为 $P(x,t)$，它表示 $x$ 点右边部分对 $x$ 点左边部分的作用力（单位面积），向右为正。考虑平衡时位于区间 $[x,x+\Delta x]$ 的一小段，如图 2 所示，其 $x$ 端受到左边部分的弹性力为 $SP(x,t)$，而 $x+\Delta x$ 端受到右边部分的弹性力为 $SP(x+\Delta x,t)$，其中 $S$ 是杆的横截面面积。设杆的质量密度为 $\rho$。对于自由振动情况，由 Newton 第二定律，该段的运动方程为

$$
\rho S\Delta x\frac{\partial^2\bar u}{\partial t^2}=SP(x+\Delta x,t)-SP(x,t),\tag{34}
$$

其中 $\bar u$ 是该小段的平均位移。两边消去 $S$，除以 $\Delta x$，并令 $\Delta x\to0$，得到

$$
\rho\frac{\partial^2u}{\partial t^2}=\frac{\partial P}{\partial x}.\tag{35}
$$

但上式并不是运动方程。只有将应力与位移的关系代入上式，才能得到运动方程。

实验表明：在弹性限度内，应力正比于应变。这称为 Hooke 定律。应变就是相对伸长。首先需要推导应变的表示。仍考虑 $[x,x+\Delta x]$ 段，当振动时，其绝对伸长为 $u(x+\Delta x,t)-u(x,t)$，而相对伸长为 $[u(x+\Delta x,t)-u(x,t)]/\Delta x$，这是该段的平均应变，当 $\Delta x\to0$，即得 $x$ 点的应变为 $\partial u/\partial x$，所以，Hooke 定律可以表示为

$$
P(x,t)=Y\frac{\partial u(x,t)}{\partial x},\tag{36}
$$

其中 $Y$ 是材料的力学参数，称为 Young 模量，对于性质均匀的材料，它与 $x$ 无关，是一个常数。

将式 (36) 代入式 (35)，即得弹性杆纵振动所满足的波动方程为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0,\tag{37}
$$

其中 $a=\sqrt{Y/\rho}$。这与弦的横振动方程 (31) 在形式上完全一致，尽管两者的物理背景颇为不同。所以，一个偏微分方程并不一定只是描述一个具体的物理过程，而可能描述一系列类似的物理现象。这一点在下节会得到进一步的印证。所以，研究偏微分方程的求解具有较普遍的意义。由于是自由振动，所以上式没有非齐次项。

如果杆在振动过程中受到外力的作用，即强迫振动，则上面的运动方程需要修正。设 $x$ 处单位体积的受力（即力密度）为 $F(x,t)$，方向为 $x$ 轴的正方向，则式 (34) 应修改为

$$
\rho S\Delta x\frac{\partial^2\bar u}{\partial t^2}=SP(x+\Delta x,t)-SP(x,t)+\bar F S\Delta x,\tag{38}
$$

其中 $\bar F$ 是该小段的平均力密度。当 $\Delta x\to0$ 时，$\bar F\to F(x,t)$。重复前面的计算并利用 Hooke 定律，得到外力作用下的波动方程为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f,\tag{39}
$$

其中 $f(x,t)=F(x,t)/\rho$ 是 $x$ 处单位质量的受力。这与外力作用下弦的横振动方程 (33) 也完全一致。

值得指出，如果杆在振动过程中一端受到外力的作用，这个外力并不是方程 (39) 中的非齐次项，而是出现在边界条件中（参看下一小节的讨论）。只有作用于杆内部各处的外力才表现为方程中的非齐次项。一个典型的例子是重力。如果杆处于竖直方向，取 $x$ 轴和位移的正方向向下，则 $f(x,t)=g$，其中 $g$ 是重力加速度。

### §3.3　定解条件

确定一个常微分方程的解需要给出初始条件，初始条件的个数与方程的阶数相同。比如确定二阶常微分方程的解需要两个初始条件，即需要同一点的函数值和一阶导数值。如果自变量表示时间，函数值表示位移，那就是要给出某一时刻（不妨取为初始时刻，即 $t=0$）的位移和速度。在数学上，即使方程的自变量并不表示时间，这些条件也称为初始条件。不过，下面提到的初始条件是就时间变量来说的，所以是物理上的初始条件。

波动方程含有对时间的二阶偏导数，所以确定其解也需要两个初始条件（initial condition）。就本节所研究的一维波动方程来说，初始条件为

$$
u|_{t=0}=\varphi(x),\qquad \left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),\tag{40}
$$

其中 $\varphi(x)$ 和 $\psi(x)$ 是已知函数。从物理上说，就是要给定弦上或杆上各点的初始位移和初始速度。对于三维的波动方程，情况与此类似。

由于我们所研究的对象是有界的连续体，所以除了初始条件，还需要知道边界上的约束情况，才能确定问题的解。这就是边界条件（boundary condition）。常见的边界条件有三类，分别讨论如下：

1. 第一类边界条件，就是给定边界上的 $u$ 值。比如弦的横振动，如果两个端点是固定的，则边界条件为

$$
u|_{x=0}=0,\qquad u|_{x=l}=0,\tag{41}
$$

其中 $x=0$ 和 $x=l$ 分别是弦的两个端点的坐标。上式同样适用于杆的两端固定的情况。它是齐次边界条件。如果边界上的 $u$ 值不为零，则称为非齐次边界条件。

![Figure 3](...)

图 3：第二类边界条件

2. 第二类边界条件，就是给定边界上的 $\partial u/\partial n$ 值，即外法向方向的方向导数。比如弹性杆在纵振动的过程中，其 $x=l$ 端受到已知力 $F(t)$（向右为正）的作用，则该端具有第二类边界条件。推导如下：考虑 $[l-\Delta x,l]$ 段，如图 3 所示，其中下图为该小段的平衡位置，上图为其运动时的位置和受力分析。由 Newton 第二定律，其运动方程为

$$
\rho S\Delta x\frac{\partial^2\bar u}{\partial t^2}=F(t)-SP(l-\Delta x,t),
$$

其中 $\bar u$ 是该小段的平均位移。令 $\Delta x\to0$，并利用 Hooke 定律，得到

$$
\left.\frac{\partial u}{\partial x}\right|_{x=l}=\frac1{YS}F(t).\tag{42}
$$

这是第二类非齐次边界条件。特别地，如果 $F(t)=0$，则边界条件为

$$
\left.\frac{\partial u}{\partial x}\right|_{x=l}=0.\tag{43}
$$

即自由端具有第二类齐次边界条件。这是一个常用的结论，读者应该熟悉掌握。

**思考**　弹性杆在纵振动的过程中，其 $x=0$ 端受到已知力 $F(t)$（向右为正）的作用，则该端的边界条件如何？

**答案**　$\displaystyle \left.\frac{\partial u}{\partial x}\right|_{x=0}=\frac1{YS}F(t)$。

3. 第三类边界条件，就是给定边界上的 $u$ 和 $\partial u/\partial n$ 的线性组合。比如弹性杆在纵振动的过程中，其一端与弹簧连接（弹簧的另一端固定），且当该端处于平衡位置（即没有位移）时，弹簧也处于平衡状态，则该端具有第三类齐次边界条件。以 $x=l$ 端为例，推导如下：仍然考虑 $[l-\Delta x,l]$ 段，如图 4 所示，其中下图为该小段的平衡位置，上图为其运动时的位置和受力分析。由 Newton 第二定律，其运动方程为

$$
\rho S\Delta x\frac{\partial^2\bar u}{\partial t^2}=-ku(l,t)-SP(l-\Delta x,t),
$$

其中 $\bar u$ 是该小段的平均位移。令 $\Delta x\to0$，并利用 Hooke 定律，得到

$$
\left.\left(\frac{\partial u}{\partial x}+\frac{k}{YS}u\right)\right|_{x=l}=0.\tag{44a}
$$

![Figure 4](...)

图 4：第三类边界条件

这是第三类齐次边界条件。对 $x=0$ 端，类似可得

$$
\left.\left(\frac{\partial u}{\partial x}-\frac{k}{YS}u\right)\right|_{x=0}=0.\tag{44b}
$$

注意式 (44a) 与式 (44b) 中的符号差异。这是第三类边界条件的特点。如果杆与弹簧连接的端点处于平衡位置时，弹簧已经有一定的形变，则该端将具有第三类非齐次边界条件。

**思考**　① 如果杆的 $x=l$ 端处于平衡位置时，弹簧已经有伸长量 $d$（若 $d<0$，则为缩短），则该端的边界条件如何？② 对 $x=0$ 端，考虑同样的问题。③ 在式 (44) 中，若 $k\to0$，结果如何？若 $k\to\infty$，结果又如何？其物理意义何在？

**答案**　① $\displaystyle \left.\left(\frac{\partial u}{\partial x}+\frac{k}{YS}u\right)\right|_{x=l}=\frac{kd}{YS}$。　② $\displaystyle \left.\left(\frac{\partial u}{\partial x}-\frac{k}{YS}u\right)\right|_{x=0}=\frac{kd}{YS}$。

以上讨论了三类常见的边界条件。初始条件和边界条件统称为定解条件。定解条件和偏微分方程一起构成定解问题。一个定解问题可以在边界的不同部分具有不同类型的边界条件。比如杆的纵振动，可以一端固定，一端与弹簧连接。

最后举例介绍一下所谓的连接条件（linking condition）。考虑弹性杆由材料性质不同的两段构成（但假设横截面均匀），其连接点的坐标为 $d$。由于连接点两边的材料具有不同的力学参数，故其两边的波动方程具有不同的形式（波速不同）。但两边的振动又是互相联系的，这表现为以下的连接条件：

$$
u|_{x=d-0}=u|_{x=d+0},\qquad
Y_1\left.\frac{\partial u}{\partial x}\right|_{x=d-0}=Y_2\left.\frac{\partial u}{\partial x}\right|_{x=d+0},
$$

其中 $d-0$（$d+0$）表示 $d$ 点的左（右）极限，$Y_1$（$Y_2$）是 $d$ 点左（右）边部分的 Young 模量。上式第一式表示 $d$ 点两侧位移连续，否则杆将会断裂；第二式表示两侧应力相等，可以通过考虑 $[d-\Delta x,d+\Delta x]$ 段的运动方程并令 $\Delta x\to0$ 得到。对于此类涉及多种介质的问题，连接条件也是定解条件的一部分，是确定具体问题的解所必不可少的。对于热传导、静电场等方程的定解问题，如果涉及多种介质，则不同介质的交界处也需要连接条件。与边界条件类似，连接条件的推导也需要一定的物理知识，今后对此不再具体讨论。

### §3.4　定解问题的适定性

定解问题的适定性指的是其解的存在性、唯一性和稳定性。存在性指的是有解。唯一性指的是解是确定的，没有任意性。读者对这两点已经比较熟悉。稳定性指的是，当定解条件有微小的变化时，所引起的解的变化也是微小的。这在实际问题中非常重要，因为由测量给出的定解条件与实际情况总会有一定的误差，如果定解问题不具有稳定性，那么理论计算所得到的解将不能反映客体的运动情况，因而是没有实际意义的。

从物理学的角度来看，如果在导出偏微分方程时对客体和物理过程所做的简化和近似是合理的，定解条件恰当地描写了客观情况，那么这样的定解问题应该具有适定性。但是，从数学上研究各类定解问题的适定性也是有实际意义的。如果数学上证明了一个定解问题是不适定的，这可能说明物理学家在推导偏微分方程时作了不合理的近似，或给出了不恰当的定解条件，从而促使他们对研究结果做出改进或修正。

本书除本章研究各类方程的导出外，以后各章主要研究定解问题的求解，并尽可能对解的物理图像作一些分析和说明，至于定解问题的适定性，今后不再考虑。这里只对边界条件作一点说明，虽然本书研究的几类方程都具有对空间变量的二阶偏导，但边界条件只能有一个（虽然在边界的不同部分可以有不同类型的边界条件），否则定解问题将是不适定的。比如，在边界上同时给出 $u$ 和 $\partial u/\partial n$，通常会导致定解条件自相矛盾，从而使得定解问题的解不存在。这与初始条件的情况是颇为不同的，后者的个数与方程对时间的偏导数的阶数相同。

**习题 1**　弦在阻尼介质中振动，$t$ 时刻 $x$ 处单位长度所受阻力为 $F(x,t)=-R\partial u(x,t)/\partial t$（其中 $R$ 称为阻力系数）。试推导弦在该阻尼介质中的波动方程。

**习题 2**　弹性均匀细杆，$x=0$ 端固定，$x=l$ 端被拉长至 $x=l+d$ 并保持静止（$d$ 不超过弹性限度），$t=0$ 时突然放开 $x=l$ 端，写出杆作纵振动的定解问题。

## §4　热传导方程和扩散方程

### §4.1　热传导方程

在导热介质中，如果温度分布不均匀，热量就会从温度高的地方向温度低的地方流动，这就是热传导现象，读者在普通物理课程中已经有所了解。热量的流动可以用热流强度来描写。热流强度定义为单位时间内垂直流过单位面积的热量，记作 $\boldsymbol q$，其方向即热量流动的方向，一般来说，它是 $\boldsymbol r$ 和 $t$ 的函数。设介质中的温度分布为 $u(\boldsymbol r,t)$。实验表明：

$$
\boldsymbol q=-k\nabla u.\tag{45}
$$

即热流强度由温度分布决定。上式称为热传导定律（也称为 Fourier 定律），它表明热量沿着温度下降最快的方向流动。$k$ 称为热导率，它与介质的材料有关，在非均匀介质中，它可以是 $\boldsymbol r$ 的函数。原则上，$k$ 还与温度有关（如此则热传导方程将成为非线性方程），但如果温度的变化范围不大，可以近似地认为它与温度无关。

下面从能量守恒定律和热传导定律出发，推导热传导方程。

在介质中任取一区域 $V$，其边界面为 $S=\partial V$，如图 5 所示。设介质中有热源，热源强度为 $F(\boldsymbol r,t)$，它表示 $t$ 时刻 $\boldsymbol r$ 处单位时间单位体积内放出的热量。记介质的质量密度为 $\rho$，比热容为 $c$，则区域 $V$ 内单位时间由于温度升高而增加的能量为

$$
\int_V c\rho\frac{\partial u}{\partial t}\,d\boldsymbol r,
$$

其中 $d\boldsymbol r$ 是体积元。这一能量有两个来源，其一是边界面流入的热量

$$
-\int_S\boldsymbol q\cdot d\boldsymbol\sigma,
$$

其中 $d\boldsymbol\sigma$ 是边界面的面积元，其方向为边界面的外法向方向；其二是热源产生的热量

$$
\int_VF\,d\boldsymbol r.
$$

由能量守恒定律，就有

$$
\int_Vc\rho\frac{\partial u}{\partial t}\,d\boldsymbol r=-\int_S\boldsymbol q\cdot d\boldsymbol\sigma+\int_VF\,d\boldsymbol r.\tag{46}
$$

![Figure 5](...)

图 5：介质中的区域 $V$，$\boldsymbol n$ 为其边界面的外法向单位矢量。

利用数学上的 Gauss 定理，右边第一项可以化为体积分 $-\int_V\nabla\cdot\boldsymbol q\,d\boldsymbol r$，再由区域 $V$ 的任意性，就得到

$$
c\rho\frac{\partial u}{\partial t}+\nabla\cdot\boldsymbol q=F.\tag{47}
$$

把热传导定律 (45) 代入上式，就得到热传导方程

$$
c\rho\frac{\partial u}{\partial t}-\nabla\cdot(k\nabla u)=F.\tag{48}
$$

在直角坐标系中，上式的具体形式是

$$
c\rho\frac{\partial u}{\partial t}-\frac{\partial}{\partial x}\left(k\frac{\partial u}{\partial x}\right)-\frac{\partial}{\partial y}\left(k\frac{\partial u}{\partial y}\right)-\frac{\partial}{\partial z}\left(k\frac{\partial u}{\partial z}\right)=F.\tag{49}
$$

但是，式 (48) 非常简练，且不依赖于坐标系的选择，读者应该逐步熟悉这样的表达式。

对于均匀介质，$k$ 是常数，则式 (48) 简化为

$$
\frac{\partial u}{\partial t}-a^2\nabla^2u=f,\tag{50}
$$

其中 $a=\sqrt{k/c\rho}$，$f(\boldsymbol r,t)=F(\boldsymbol r,t)/c\rho$。这就是 §1 中介绍的输运方程的标准形式，如果没有热源，就得到相应的齐次方程。

考虑一均匀导热细杆的热传导问题。设细杆的侧面绝热，由于细杆的横截面面积很小，故各横截面上的温度分布可以很快达到均匀，之后，温度在空间上只依赖于杆长方向的坐标 $x$，热量也只沿着 $x$ 方向流动，于是就得到一个一维的热传导问题。此时，式 (50) 简化为

$$
\frac{\partial u}{\partial t}-a^2\frac{\partial^2u}{\partial x^2}=f,\tag{51}
$$

其中 $u=u(x,t)$，$f=f(x,t)$。若无热源，则得相应的齐次方程。

### §4.2　扩散方程

在介质中引入一种杂质，如果杂质浓度分布不均匀，它就会从浓度高的地方向浓度低的地方流动，这就是扩散现象。值得指出，如果是混合气体，则各成分的扩散应该在温度和总压强均匀的条件下进行。如果因总压强不均匀而产生气流，就不是扩散过程。这是读者在普通物理课程中已经了解的。杂质的扩散可以用扩散流强度来描写。扩散流强度定义为单位时间内垂直流过单位面积的杂质质量，记作 $\boldsymbol q$，其方向即杂质流动的方向，一般来说，它是 $\boldsymbol r$ 和 $t$ 的函数。设介质中的杂质浓度分布为 $u(\boldsymbol r,t)$，它表示 $t$ 时刻 $\boldsymbol r$ 处单位体积内的杂质质量。实验表明：

$$
\boldsymbol q=-D\nabla u.\tag{52}
$$

即扩散流强度由杂质浓度分布决定。上式称为扩散定律（也称为 Fick 定律），它表明杂质沿着浓度下降最快的方向流动。$D$ 称为扩散系数，它与介质的材料和扩散时介质中的温度有关，在非均匀介质中，它可以是 $\boldsymbol r$ 的函数。

下面从物质守恒定律和扩散定律出发，推导扩散方程。这与热传导方程的推导非常类似。

在介质中任取一区域 $V$，其边界面为 $S=\partial V$，仍如图 5 所示。设介质中有杂质源（比如由化学反应产生所考虑的杂质），杂质源强度为 $F(\boldsymbol r,t)$，它表示 $t$ 时刻 $\boldsymbol r$ 处单位时间单位体积产生的杂质质量。考虑区域 $V$ 中单位时间内杂质质量的增加，由物质守恒定律可得

$$
\int_V\frac{\partial u}{\partial t}\,d\boldsymbol r=-\int_S\boldsymbol q\cdot d\boldsymbol\sigma+\int_VF\,d\boldsymbol r.\tag{53}
$$

利用数学上的 Gauss 定理将右边第一项化为体积分，并考虑到区域 $V$ 的任意性，就得到

$$
\frac{\partial u}{\partial t}+\nabla\cdot\boldsymbol q=F.\tag{54}
$$

把扩散定律 (52) 代入上式，就得到扩散方程

$$
\frac{\partial u}{\partial t}-\nabla\cdot(D\nabla u)=F.\tag{55}
$$

对于均匀介质，$D$ 是常数，则式 (55) 简化为

$$
\frac{\partial u}{\partial t}-a^2\nabla^2u=F,\tag{56}
$$

其中 $a=\sqrt D$。这与热传导方程在形式上完全一致。如果没有杂质源，就得到相应的齐次方程。

考虑杂质气体在均匀细管内的扩散问题。由于细管的侧面封闭，且横截面面积很小，故各横截面上的杂质浓度分布可以很快达到均匀，之后，浓度在空间上只依赖于管长方向的坐标 $x$，杂质也只沿 $x$ 方向扩散，于是就得到一个一维的扩散问题。这与上一小节所考虑的一维热传导问题也非常类似。

在微观上，热传导和扩散过程都是通过分子的碰撞完成的。碰撞使得能量在分子间重新分布，这就是热传导过程，也就是能量的输运过程。类似地，碰撞改变了不同物质的分子在空间上的分布，这就是扩散过程，也就是分子数的输运过程。由于这两种过程具有类似的微观机制，所以其所满足的方程具有同样的形式就不足为奇了。

### §4.3　定解条件

输运方程具有对时间的一阶偏导数，所以初始条件只有一个，就是给定初始时刻的温度分布或浓度分布

$$
u|_{t=0}=\varphi(\boldsymbol r),\qquad \boldsymbol r\in\Omega.\tag{57}
$$

其中 $\Omega$ 是研究对象在三维空间所占据的区域，其边界面记作 $\partial\Omega$。

除了初始条件，确定一个具体问题的解也还需要边界条件。常见的边界条件有三类，其数学形式与上节所述类似。但是，同样类型的边界条件，在输运问题中对应于不同的物理状况，所以下面对三类边界条件分别举例加以讨论。上节主要以一维问题为例，本节则主要以三维热传导问题为例。对于扩散问题的边界条件可作类似的讨论。

1. 第一类边界条件。将研究对象置于温度已知的环境中，如果该物质的导热性能良好，则其表面温度可以很快达到与环境一致，故边界条件为

$$
u|_{\partial\Omega}=u_0(\boldsymbol r,t),\qquad \boldsymbol r\in\partial\Omega.\tag{58}
$$

其中 $u_0(\boldsymbol r,t)$ 是已知函数，表示环境温度。特别地，如果 $u_0(\boldsymbol r,t)=0$，就得到第一类齐次边界条件。

2. 第二类边界条件。比如边界面上有已知的热流流入，其强度为 $q(\boldsymbol r,t)$ 且垂直于边界面。从边界面的内侧看，垂直流入的热流强度按热传导定律应为 $k(\partial u/\partial n)|_{\partial\Omega}$，从边界面的外侧看，则为已知量 $q(\boldsymbol r,t)$。由于边界面上不能有热量的积聚（注意边界面是一个几何概念，是没有厚度的），故必有

$$
\left.\frac{\partial u}{\partial n}\right|_{\partial\Omega}=\frac1kq(\boldsymbol r,t),\qquad \boldsymbol r\in\partial\Omega.\tag{59}
$$

特别地，如果 $q(\boldsymbol r,t)=0$，则得

$$
\left.\frac{\partial u}{\partial n}\right|_{\partial\Omega}=0.\tag{60}
$$

即绝热的边界面（对一维问题，就是绝热的端点）具有第二类齐次边界条件。这是一个常用的结论，读者应该熟悉掌握。

3. 第三类边界条件。比如边界面按 Newton 冷却定律与外界交换热量，即从介质表面流出的法向热流强度正比于表面内侧与外侧的温度差。即

$$
q_n|_{\partial\Omega}=b[u|_{\partial\Omega}-u_0(\boldsymbol r,t)],\qquad \boldsymbol r\in\partial\Omega,
$$

其中 $q_n=\boldsymbol q\cdot\boldsymbol n$，$\boldsymbol n$ 是边界面的外法向单位矢量，$b$ 是常数，称为介质的热交换系数。另一方面，按热传导定律有 $q_n|_{\partial\Omega}=-k(\partial u/\partial n)|_{\partial\Omega}$，故有

$$
\left.\left(\frac{\partial u}{\partial n}+hu\right)\right|_{\partial\Omega}=hu_0(\boldsymbol r,t),\qquad \boldsymbol r\in\partial\Omega,\tag{61}
$$

其中 $h=b/k$。特别地，如果 $u_0(\boldsymbol r,t)=0$，则得第三类齐次边界条件

$$
\left.\left(\frac{\partial u}{\partial n}+hu\right)\right|_{\partial\Omega}=0.\tag{62}
$$

对于一维细杆的热传导问题，设细杆的两端坐标为 $x=0$ 和 $x=l$，则上式成为

$$
\left.\left(\frac{\partial u}{\partial x}-hu\right)\right|_{x=0}=0,
\qquad
\left.\left(\frac{\partial u}{\partial x}+hu\right)\right|_{x=l}=0.\tag{63}
$$

注意两式中的符号差异。上节已经看到，波动方程的第三类边界条件中也存在类似的符号差异。

**习题**　混凝土浇灌后逐渐放出水化热，放热速率正比于当时尚储存着的水化热密度 $Q$，即 $dQ/dt=-\beta Q$。试推导浇灌后混凝土内的热传导方程。

## §5　稳定场方程

### §5.1　稳定温度分布和稳定浓度分布

考虑热传导方程 (50)，如果 $f=f(\boldsymbol r)$ 与 $t$ 无关，且边界条件也与 $t$ 无关，则长时间后，温度分布有可能达到稳定状态。这时温度只是 $\boldsymbol r$ 的函数，故 $\partial u/\partial t=0$，从而方程 (50) 化为

$$
\nabla^2u=-\frac{f}{a^2}.\tag{64}
$$

可见稳定温度分布满足 Poisson 方程。如果没有热源，则得相应的齐次方程，即 Laplace 方程。对于稳定浓度分布，情况是完全类似的。

**注**　① $f$ 和边界条件与 $t$ 无关只是达到稳定温度分布的必要条件，而不是充分条件。比如介质表面绝热，而内部有稳定热源，这满足上述条件，但一般来说不可能达到稳定状态，除非热源产生的总热量为零，即 $\int_\Omega f(\boldsymbol r)\,d\boldsymbol r=0$。② 稳定状态不一定是平衡状态。比如均匀导热细杆，侧面绝热，左端保持较高温度 $u_1$，右端保持较低温度 $u_2$，长时间后温度分布可以达到稳定状态，但在这种稳定状态下，显然有热量源源不断地从左向右流动，所以这一稳定状态需要靠外部条件来维持，因而不是平衡状态。

### §5.2　静电场方程

在介质中，静电场的基本方程是

$$
\nabla\cdot\boldsymbol D=\rho,\qquad \nabla\times\boldsymbol E=0,\tag{65}
$$

其中 $\rho$ 是自由电荷密度，$\boldsymbol E$ 是电场强度，$\boldsymbol D$ 是电位移矢量。引入静电势 $u$ 使得

$$
\boldsymbol E=-\nabla u,\tag{66}
$$

则方程 (65) 中的第二式自动满足。为了导出 $u$ 所满足的方程，我们需要 $\boldsymbol D$ 和 $\boldsymbol E$ 之间的关系，这称为本构关系（constitutive relation）。本构关系的形式取决于介质的性质。

一般来说，本构关系可以非常复杂。首先，$\boldsymbol D$ 可以是 $\boldsymbol E$ 的非线性函数。当场强不是很大时，可以忽略非线性效应。这时，$\boldsymbol D$ 和 $\boldsymbol E$ 之间虽然是线性关系，但一般来说是各向异性的，即 $\boldsymbol D$ 的每一个分量依赖于 $\boldsymbol E$ 的各分量：$D_i(\boldsymbol r)=\sum_j\epsilon_{ij}(\boldsymbol r)E_j(\boldsymbol r)$。对于时变场，$\boldsymbol D$ 还可以依赖于 $\boldsymbol E$ 的历史。如果介质是各向同性的，则 $\epsilon_{ij}(\boldsymbol r)=\delta_{ij}\epsilon(\boldsymbol r)$，这时有 $\boldsymbol D(\boldsymbol r)=\epsilon(\boldsymbol r)\boldsymbol E(\boldsymbol r)$。如果介质还是均匀的，则 $\epsilon(\boldsymbol r)$ 与 $\boldsymbol r$ 无关，是一个常数。通常教科书上处理的就是这种情形。尽管这是最简单的情形，但却是非常重要的，因为存在若干可以解析求解的具体情况（取决于区域的形状和边界条件）。一方面，通过对解析结果的分析可以使我们对问题的物理图像有深入的了解。另一方面，简单情况的处理方法对于复杂情况的处理具有借鉴意义，而且其结果也常常可以作为近似求解的出发点。

对于线性、各向同性的均匀介质，本构关系具有下列简单形式

$$
\boldsymbol D(\boldsymbol r)=\epsilon\boldsymbol E(\boldsymbol r),\tag{67}
$$

其中 $\epsilon$ 是介质的介电常数。结合上述诸式，易得

$$
\nabla^2u=-\frac{\rho}{\epsilon}.\tag{68}
$$

可见在所考虑的介质中，静电势满足 Poisson 方程，在没有自由电荷的区域则满足 Laplace 方程。对于真空，只需将 $\epsilon$ 换成真空介电常数 $\epsilon_0$。

### §5.3　定解条件

稳定场方程不含对时间的偏导数，所以不需要初始条件。常见的边界条件有三类，与前类似，不再赘述。

## 补充习题

1. 研究截面非均匀（即随着杆长变化）的弹性细杆的纵振动方程。

   **答案：** 如果材料的性质（此处主要指 Young 模量）也随着杆长变化，则振动方程为

   $$
   \rho S(x)\frac{\partial^2u}{\partial t^2}=\frac{\partial}{\partial x}\left[Y(x)S(x)\frac{\partial u}{\partial x}\right]+S(x)F(x,t),
   $$

   其中 $F(x,t)$ 是杆上单位体积的受力。如果材料性质均匀，则 Young 模量为常数，上式简化为

   $$
   \rho S(x)\frac{\partial^2u}{\partial t^2}=Y\frac{\partial}{\partial x}\left[S(x)\frac{\partial u}{\partial x}\right]+S(x)F(x,t).
   $$

2. 长为 $l$ 的柔软重弦（重力不能忽略），上端固定在以角速度 $\omega$ 旋转的竖直轴上，下端自由。由于重力的作用，弦的平衡位置应该是竖直线。当受到扰动时，弦可以相对于该竖直线作横振动。试推导该横振动方程。[^1]

3. 长为 $l$ 的柔软轻弦（重力可以忽略），一端固定在以角速度 $\omega$ 旋转的竖直轴上，一端自由。由于离心力的作用，弦的平衡位置应该是水平线。当受到扰动时，弦可以相对于该水平线作横振动。试推导该横振动方程。[^2]

[^1]: 常见教科书上关于本题的答案，以及学习指导书上关于本题的解法存在错误。参看以下文献：林琼桂，关于转动弦的波动方程，大学物理 **27** (2) (2008) 7-9.

[^2]: 参看上题所引文献，其中对本问题作了比较深入的讨论。
