# 第十三章 Green 函数法\*

**目录**

- §1 Green 函数的一般概念　270
- §2 稳定场方程的 Green 函数法　271
  - §2.1 *Poisson* 方程的基本解　271
  - §2.2 第二 *Green* 公式　273
  - §2.3 *Poisson* 方程第一边值问题　274
  - §2.4 *Laplace* 方程第一边值问题　276
  - §2.5 \*Green 函数的对称性　277
- §3 特殊区域上 Green 函数的制作　277
  - §3.1 球内 *Laplace* 方程第一边值问题　278
  - §3.2 上半平面 *Laplace* 方程第一边值问题　280
  - §3.3 圆内 *Laplace* 方程第一边值问题　281
- §4 \*波动方程的 Green 函数法　283
- §5 \*热传导方程的 Green 函数法　283
- 补充习题　284

前面几章研究数学物理方程定解问题的求解，主要介绍了分离变量法和 Fourier 变换法。Fourier 变换法的精神也是分离变量，只是形式上不甚明显。在曲线坐标系中对数理方程分离变量常常导出一些特殊的常微分方程，而特殊函数大多是这样一些常微分方程的解，所以特殊函数也是由分离变量法派生出来的。本章研究求解数学物理方程定解问题的 Green 函数法。这是一种新的解法。它的思路是先求出与点源对应的解，即 Green 函数，然后将原来的定解问题的解表示为包含 Green 函数和定解条件的积分。

## §1 Green 函数的一般概念

线性偏微分方程可以统一写成以下形式

$$
Lu=f,\tag{1}
$$

其中 $u$ 是未知函数，即我们所研究的物理量，$f$ 是已知函数，$L$ 是线性微分算符。从数学上说，上式可以涉及任意的有限个自变量。对于常见的物理问题，自变量是三维空间变量 $\boldsymbol r$ 和时间变量 $t$（稳定场方程实际上不含 $t$）。因此，$u=u(\boldsymbol r,t)$，$f=f(\boldsymbol r,t)$，$L=L(\boldsymbol r,t;\nabla,\partial/\partial t)$。具体来说，我们熟悉的几类方程所对应的算符如下：

波动方程：

$$
L=\frac{\partial^2}{\partial t^2}-a^2\nabla^2;\tag{2a}
$$

输运方程：

$$
L=\frac{\partial}{\partial t}-a^2\nabla^2;\tag{2b}
$$

Poisson 方程：

$$
L=-\nabla^2;\tag{2c}
$$

Helmholtz 方程：

$$
L=-(\nabla^2+k^2).\tag{2d}
$$

首先考虑不含时的稳定场方程，即 Poisson 方程（包括 Laplace 方程）和 Helmholtz 方程，这时 $L=L(\boldsymbol r,\nabla)$。由 $\boldsymbol r_0$ 处的点源 $f(\boldsymbol r)=\delta(\boldsymbol r-\boldsymbol r_0)$ 所产生的场称为 Green 函数，记作 $G(\boldsymbol r,\boldsymbol r_0)$，即 $G(\boldsymbol r,\boldsymbol r_0)$ 满足方程

$$
LG(\boldsymbol r,\boldsymbol r_0)=\delta(\boldsymbol r-\boldsymbol r_0),\tag{3}
$$

故 Green 函数又称为点源影响函数。对于有界区域，Green 函数除满足上述方程外，还要满足一定的边界条件（参看下节）。无界空间的 Green 函数又称为基本解。

就无界空间来说，如果能求出 Green 函数，则方程 (1) 的解可以表示为

$$
u(\boldsymbol r)=\int G(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r_0)\,d\boldsymbol r_0.\tag{4}
$$

事实上，注意到 $L$ 只对 $\boldsymbol r$ 作用，而不对 $\boldsymbol r_0$ 作用，就有

$$
Lu(\boldsymbol r)=L\int G(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r_0)\,d\boldsymbol r_0
=\int [LG(\boldsymbol r,\boldsymbol r_0)]f(\boldsymbol r_0)\,d\boldsymbol r_0
=\int \delta(\boldsymbol r-\boldsymbol r_0)f(\boldsymbol r_0)\,d\boldsymbol r_0=f(\boldsymbol r).
$$

式 (4) 的物理意义是：$G(\boldsymbol r,\boldsymbol r_0)$ 是 $\boldsymbol r_0$ 处的单位点源（比如单位点电荷）产生的场，乘以 $f(\boldsymbol r_0)d\boldsymbol r_0$ 就是 $\boldsymbol r_0$ 处大小（比如电量）为 $f(\boldsymbol r_0)d\boldsymbol r_0$ 的点源产生的场，再积分就是把所有源产生的场叠加，如此即得总场 $u(\boldsymbol r)$。

**例**　考虑 Poisson 方程 $\nabla^2u(\boldsymbol r)=-f(\boldsymbol r)$，相应的 Green 函数满足 $\nabla^2G(\boldsymbol r,\boldsymbol r_0)=-\delta(\boldsymbol r-\boldsymbol r_0)$，规定 $r\to\infty$ 时 $G(\boldsymbol r,\boldsymbol r_0)\to0$，可得 Green 函数为（证明见下节）

$$
G(\boldsymbol r,\boldsymbol r_0)=\frac{1}{4\pi|\boldsymbol r-\boldsymbol r_0|}.\tag{5}
$$

于是

$$
u(\boldsymbol r)=\int\frac{f(\boldsymbol r_0)}{4\pi|\boldsymbol r-\boldsymbol r_0|}\,d\boldsymbol r_0.\tag{6}
$$

特别地，如果 $f(\boldsymbol r)=\rho(\boldsymbol r)/\varepsilon_0$，则所考虑的 Poisson 方程为真空中静电势所满足的方程，此时

$$
u(\boldsymbol r)=\int\frac{\rho(\boldsymbol r_0)}{4\pi\varepsilon_0|\boldsymbol r-\boldsymbol r_0|}\,d\boldsymbol r_0.\tag{7}
$$

这正是读者在电磁学中已经熟悉的公式。

对于含时的方程，如波动方程和输运方程，$L=L(\boldsymbol r,t;\nabla,\partial/\partial t)$。这时可以定义 Green 函数 $G(\boldsymbol r,t;\boldsymbol r_0,t_0)$ 满足方程

$$
LG(\boldsymbol r,t;\boldsymbol r_0,t_0)=\delta(\boldsymbol r-\boldsymbol r_0)\delta(t-t_0).\tag{8}
$$

上式右边的 $\delta$ 函数不仅是空间上的点源，也是时间上的瞬时源。含时的 Green 函数通常还要满足一定的初始条件。如果空间是有界的，则还要满足适当的边界条件。

含时的问题通常要复杂得多，即使对于无界空间，方程 (1) 的解也不能简单地表示为

$$
u(\boldsymbol r,t)=\int G(\boldsymbol r,t;\boldsymbol r_0,t_0)f(\boldsymbol r_0,t_0)\,d\boldsymbol r_0\,dt_0.\tag{9}
$$

因为初始条件尚未考虑在内。实际上，如果没有给定适当的初始条件，上式中的 Green 函数本身也没有确切的定义。

本章的重点是稳定场方程的 Green 函数法。波动方程和热传导方程的 Green 函数法也加以介绍，供读者参考。

## §2 稳定场方程的 Green 函数法

稳定场方程不涉及时间变量。Poisson 方程（包括 Laplace 方程）和 Helmholtz 方程都是稳定场方程。我们主要研究 Poisson 方程。对于 Helmholtz 方程，可用类似的方法加以研究。

### §2.1 *Poisson* 方程的基本解

考虑无界空间的 Poisson 方程

$$
\nabla^2u(\boldsymbol r)=-f(\boldsymbol r).\tag{10}
$$

它的 Green 函数，即基本解 $G(\boldsymbol r,\boldsymbol r_0)$ 满足下列方程

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_0)=-\delta(\boldsymbol r-\boldsymbol r_0).\tag{11}
$$

在三维情况下，规定 $r\to\infty$ 时 $G(\boldsymbol r,\boldsymbol r_0)\to0$，则

$$
G(\boldsymbol r,\boldsymbol r_0)=\frac{1}{4\pi|\boldsymbol r-\boldsymbol r_0|}.\tag{12}
$$

**证明**　令 $\boldsymbol R=\boldsymbol r-\boldsymbol r_0$，$R=|\boldsymbol R|$，则 $G(\boldsymbol r,\boldsymbol r_0)=1/4\pi R$。显然，$\nabla=\nabla_{\boldsymbol R}$，后者是对 $\boldsymbol R$ 求导的梯度算符。当 $\boldsymbol r\ne\boldsymbol r_0$，即 $R\ne0$，有

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_0)=\nabla^2_{\boldsymbol R}\left(\frac{1}{4\pi R}\right)
=\frac{1}{R^2}\frac{d}{dR}\left[R^2\frac{d}{dR}\left(\frac{1}{4\pi R}\right)\right]=0.
$$

故式 (11) 成立。但是，当 $\boldsymbol r=\boldsymbol r_0$，即 $R=0$，Green 函数有奇性，而且该点是球坐标系的奇点，上述计算不成立。考虑以 $\boldsymbol r_0$ 为球心，$\varepsilon$ 为半径的球 $K_\varepsilon$，其表面记作 $S_\varepsilon$，则

$$
\int_{K_\varepsilon}\nabla^2G(\boldsymbol r,\boldsymbol r_0)\,d\boldsymbol r
=\int_{K_\varepsilon}\nabla\cdot\nabla\left(\frac{1}{4\pi R}\right)d\boldsymbol R
=\int_{S_\varepsilon}\nabla\left(\frac{1}{4\pi R}\right)\cdot d\boldsymbol\sigma
=-\int_{S_\varepsilon}\frac{1}{4\pi\varepsilon^2}\,d\sigma=-1,
$$

其中已略去了 $\nabla_{\boldsymbol R}$ 的下标 $\boldsymbol R$，第二步利用了 Gauss 定理，$d\sigma$ 是面积元的大小。上式表明式 (11) 在 $\boldsymbol r=\boldsymbol r_0$ 处也成立，因为两边在该点的奇性一致。故式 (11) 处处成立。证毕。

读者可能已经注意到，上面利用 Gauss 定理将球内的体积分化为球面上的面积分，按 Gauss 定理成立的条件，必须要求 $G(\boldsymbol r,\boldsymbol r_0)$ 在球内具有连续的二阶偏导数，但 $G(\boldsymbol r,\boldsymbol r_0)$ 在球心有奇性，并不满足这一条件，所以上面的证明是不严格的。下面给出一个比较严格的论证，供读者参考。

考虑方程

$$
\nabla^2G_\varepsilon(\boldsymbol r,\boldsymbol r_0)=-\delta_\varepsilon(|\boldsymbol r-\boldsymbol r_0|)=-\delta_\varepsilon(R),\tag{13}
$$

其中 $\delta_\varepsilon(R)$ 是连续函数，只当 $R<\varepsilon$ 时不为零，且满足 $\int\delta_\varepsilon(R)d\boldsymbol r=1$。因此当 $\varepsilon\to0$ 时，$\delta_\varepsilon(R)\to\delta(\boldsymbol r-\boldsymbol r_0)$，此时应有 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)\to G(\boldsymbol r,\boldsymbol r_0)$。因此，我们的任务是求出 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)$ 并证明当 $\varepsilon\to0$ 时其极限形式与式 (12) 一致。

由于源具有球对称性，且没有非球对称边界条件的限制，故 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)$ 应该只依赖于 $R$ 和 $\varepsilon$。记 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)=g(R,\varepsilon)$，则式 (13) 化为

$$
\frac{1}{R^2}\frac{d}{dR}\left[R^2\frac{dg(R,\varepsilon)}{dR}\right]=-\delta_\varepsilon(R).
$$

当 $R>\varepsilon$ 时，上式右边为零，故易得 $g(R,\varepsilon)=c/4\pi R+c'$，规定 $r\to\infty$ 时 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)\to0$，则 $c'=0$。于是

$$
g(R,\varepsilon)=
\begin{cases}
\dfrac{c}{4\pi R}, & R>\varepsilon,\\[4pt]
g_0(R,\varepsilon), & R<\varepsilon,
\end{cases}\tag{14}
$$

其中 $g_0(R,\varepsilon)$ 的具体形式取决于 $\delta_\varepsilon(R)$ 的具体形式，它与 $c/4\pi R$ 在 $R=\varepsilon$ 处应该光滑连接。由于我们最后要取 $\varepsilon\to0$ 的极限，所以不需要用到 $g_0(R,\varepsilon)$ 的具体形式。现在需要确定常数 $c$ 的取值。任取 $R_0\ge\varepsilon$，在球 $K_0:R<R_0$ 上对式 (13) 左边积分，利用 Gauss 定理，有

$$
\int_{K_0}\nabla^2G_\varepsilon(\boldsymbol r,\boldsymbol r_0)\,d\boldsymbol r
=\int_{S_0}\nabla G_\varepsilon(\boldsymbol r,\boldsymbol r_0)\cdot d\boldsymbol\sigma
=\int_{S_0}\nabla\left(\frac{c}{4\pi R}\right)\cdot d\boldsymbol\sigma
=-\int_{S_0}\frac{c}{4\pi R_0^2}\,d\sigma=-c,
$$

其中 $S_0$ 是 $K_0$ 的边界面，即球面 $R=R_0$。由于 $\delta_\varepsilon(R)$ 是连续函数，式 (13) 的解 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)$ 应具有连续的二阶偏导数，故上述演算中利用 Gauss 定理是合理的。按假定，式 (13) 右边的积分为 $-1$，故 $c=1$。令 $\varepsilon\to0$，则得到 $\boldsymbol r\ne\boldsymbol r_0$ 时，$G_\varepsilon(\boldsymbol r,\boldsymbol r_0)\to1/4\pi R$，这正是我们期望的结果。

读者可能还会想，尽管对于任何有限的 $\varepsilon$，式 (13) 均成立，但对 $G_\varepsilon(\boldsymbol r,\boldsymbol r_0)$ 和 $\delta_\varepsilon(R)$ 分别取 $\varepsilon\to0$ 的极限后，它是否还继续成立？或者说，Laplace 算符 $\nabla^2$ 与求极限运算是否可以交换次序？其实，正确的理解应该是，$G(\boldsymbol r,\boldsymbol r_0)$，即式 (11) 的解，正是用 $\lim_{\varepsilon\to0}G_\varepsilon(\boldsymbol r,\boldsymbol r_0)$ 来定义的。

基本解的形式与空间的维数有关。在二维空间，我们把位置矢量记作 $\boldsymbol\rho$，以示与三维空间的区别。相应地，式 (11) 写作

$$
\nabla^2G(\boldsymbol\rho,\boldsymbol\rho_0)=-\delta(\boldsymbol\rho-\boldsymbol\rho_0).\tag{15}
$$

其解可取为

$$
G(\boldsymbol\rho,\boldsymbol\rho_0)=-\frac{1}{2\pi}\ln|\boldsymbol\rho-\boldsymbol\rho_0|
=\frac{1}{2\pi}\ln\frac{1}{|\boldsymbol\rho-\boldsymbol\rho_0|},\tag{16}
$$

或再加上一个常数。在二维空间，一般不能取 $\rho\to\infty$ 为 $G(\boldsymbol\rho,\boldsymbol\rho_0)$ 的零点。上式的第一表式一般用于计算，第二表式与式 (12) 形式上类似，可以帮助记忆。证明上式满足式 (15) 的过程与三维情况完全类似，读者应该尝试自己做出证明。下面的证明供参考。

**证明**　令 $\boldsymbol R=\boldsymbol\rho-\boldsymbol\rho_0$，$R=|\boldsymbol R|$，则 $G(\boldsymbol\rho,\boldsymbol\rho_0)=-(1/2\pi)\ln R$。当 $\boldsymbol\rho\ne\boldsymbol\rho_0$，即 $R\ne0$，有

$$
\nabla^2G(\boldsymbol\rho,\boldsymbol\rho_0)
=-\frac{1}{2\pi}\nabla^2_{\boldsymbol R}\ln R
=-\frac{1}{2\pi}\frac{1}{R}\frac{d}{dR}\left(R\frac{d}{dR}\ln R\right)=0.
$$

故式 (15) 成立。但是，当 $\boldsymbol\rho=\boldsymbol\rho_0$，即 $R=0$，Green 函数有奇性，而且该点是极坐标系的奇点，上述计算不成立。考虑以 $\boldsymbol\rho_0$ 为圆心，$\varepsilon$ 为半径的圆 $K_\varepsilon$，其边界（即圆周）记作 $L_\varepsilon$，则

$$
\int_{K_\varepsilon}\nabla^2G(\boldsymbol\rho,\boldsymbol\rho_0)\,d\boldsymbol\rho
=\int_{K_\varepsilon}\nabla\cdot\nabla\left(-\frac{1}{2\pi}\ln R\right)d\boldsymbol R
=\int_{L_\varepsilon}\nabla\left(-\frac{1}{2\pi}\ln R\right)\cdot\hat{\boldsymbol R}\,ds
=-\int_{L_\varepsilon}\frac{1}{2\pi\varepsilon}\,ds=-1,
$$

其中已略去了 $\nabla_{\boldsymbol R}$ 的下标 $\boldsymbol R$，第二步利用了 Gauss 定理（其实就是 Green 公式），$ds$ 是弧元，$\hat{\boldsymbol R}$ 是圆周的外法向单位矢量。上式表明式 (15) 在 $\boldsymbol\rho=\boldsymbol\rho_0$ 处也成立，因为两边在该点的奇性一致。故式 (15) 处处成立。证毕。

当然，读者还可以仿照三维情况做出更严格的讨论。

容易看出，二维基本解的奇性比三维情况为弱。可以证明，一维空间基本解的奇性更弱。

### §2.2 第二 *Green* 公式

设 $u(\boldsymbol r)$ 和 $v(\boldsymbol r)$ 在区域 $D$ 上具有连续的二阶偏导数，在闭域 $\bar D$ 上具有连续的一阶偏导数，$S=\partial D$，$\boldsymbol n$ 是 $S$ 的外法向单位矢量，则

$$
\int_D(u\nabla^2v-v\nabla^2u)\,d\boldsymbol r
=\int_S(u\nabla v-v\nabla u)\cdot d\boldsymbol\sigma
=\int_S\left(u\frac{\partial v}{\partial n}-v\frac{\partial u}{\partial n}\right)d\sigma.\tag{17}
$$

这称为第二 Green 公式。下面将会引用这一公式。

**证明**　由算符 $\nabla$ 的矢量性质和微分性质，易得（也可用直角坐标系的形式直接验证）

$$
\nabla\cdot(u\nabla v)=\nabla u\cdot\nabla v+u\nabla^2v,\qquad
\nabla\cdot(v\nabla u)=\nabla u\cdot\nabla v+v\nabla^2u.
$$

两式相减，可得

$$
u\nabla^2v-v\nabla^2u=\nabla\cdot(u\nabla v-v\nabla u).
$$

在区域 $D$ 上积分，利用 Gauss 定理，并注意到 $\nabla v\cdot d\boldsymbol\sigma=(\partial v/\partial n)d\sigma$ 等关系，即得第二 Green 公式。证毕。

**注**　第二 Green 公式对二维空间也成立。

### §2.3 *Poisson* 方程第一边值问题

考虑有界区域 $D$ 上的 Poisson 方程第一边值问题

$$
\nabla^2u(\boldsymbol r)=-f(\boldsymbol r),\qquad \boldsymbol r\in D,\tag{18a}
$$

$$
u(\boldsymbol r)|_S=\varphi(\boldsymbol r).\tag{18b}
$$

在无界空间，我们曾经用基本解 $G(\boldsymbol r,\boldsymbol r_0)$ 和 $f(\boldsymbol r)$ 表出一般解 $u(\boldsymbol r)$，见式 (4) 和式 (6)。现在由于边界条件 (18b) 的限制，式 (4) 或式 (6) 不再成立。为了得到类似的结果，我们需要修改 Green 函数的定义。

假设 Green 函数仍满足方程 (11)。将式 (11) 乘以 $u(\boldsymbol r)$，而将式 (18a) 乘以 $G(\boldsymbol r,\boldsymbol r_0)$，相减即得

$$
G(\boldsymbol r,\boldsymbol r_0)\nabla^2u(\boldsymbol r)-u(\boldsymbol r)\nabla^2G(\boldsymbol r,\boldsymbol r_0)
=u(\boldsymbol r)\delta(\boldsymbol r-\boldsymbol r_0)-G(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r).\tag{19}
$$

在区域 $D$ 上积分并利用第二 Green 公式，得到

$$
u(\boldsymbol r_0)=\int_DG(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r)\,d\boldsymbol r
+\int_S\left[G(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial n}
-u(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial n}\right]d\sigma.\tag{20}
$$

表面上看，上式似乎给出了问题的解，但它要求已知边界上的 $u$ 值和 $\partial u/\partial n$ 值，而定解问题的边界条件只能给出其中一个或两者的线性组合，所以上式不是最后的解。不难看出，只要在 Green 函数的定义中加上适当的边界条件，即可解决这一困难。

对于 Poisson 方程第一边值问题 (18)，我们定义相应的 Green 函数满足下列定解问题

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_0)=-\delta(\boldsymbol r-\boldsymbol r_0),\qquad \boldsymbol r,\boldsymbol r_0\in D,\tag{21a}
$$

$$
G(\boldsymbol r,\boldsymbol r_0)|_{\boldsymbol r\in S}=0.\tag{21b}
$$

则由式 (20) 可得定解问题 (18) 的解为

$$
u(\boldsymbol r_0)=\int_DG(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r)\,d\boldsymbol r
-\int_S\varphi(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial n}\,d\sigma.\tag{22}
$$

可见只要能求出 Green 函数，则定解问题 (18) 的解即可用积分表出，尽管这一积分不一定能积出来（指用初等函数表出结果）。因此问题归结为求解 Green 函数。由于边界条件 (21b) 的限制，它不再是原来的基本解。容易看出，这是一个特殊的 Poisson 方程定解问题，一般来说，这一定解问题并不一定比定解问题 (18) 容易求解，但对于某些特殊的区域，存在求解 Green 函数的特殊技巧（即下节研究的镜像法），所以 Green 函数法在实用上是有意义的。此外，上述解的积分表达式也有助于研究解的唯一性和稳定性，因而具有理论上的意义。

细心的读者可能已经发现，式 (22) 的第一项与式 (4) 或式 (6) 不同，后者是对源点（即 Green 函数的第二变量）积分，而前者是对场点（即 Green 函数的第一变量）积分，这似乎使得式 (22) 中第一项的物理诠释发生困难。但可以证明，Green 函数对源点和场点具有对称性，即满足 $G(\boldsymbol r,\boldsymbol r_0)=G(\boldsymbol r_0,\boldsymbol r)$，所以这一困难实际上是不存在的。前面已经看到，基本解确实具有这样的对称性。关于由定解问题 (21) 定义的 Green 函数的对称性，读者可以参看第五小节的证明。利用这一对称性，我们可以将式 (22) 写成更易于理解的形式（参看式 (33)），但在实际计算中只要用式 (22) 就可以了。

常常有读者提出这样的问题：为什么要如此定义 Green 函数？如果读者自己动手将上述推导过程做一遍，就可以体会其精神所在。其实，Green 函数是一个辅助工具，它的作用类似于读者以前证明几何问题时所作的辅助线。如何作辅助线取决于如何有助于问题的解决。Green 函数的定义也是这样。

**思考**　对有界区域 $D$ 上的 Poisson 方程第二边值问题

$$
\nabla^2u(\boldsymbol r)=-f(\boldsymbol r),\qquad \boldsymbol r\in D,
$$

$$
\left.\frac{\partial u(\boldsymbol r)}{\partial n}\right|_S=\psi(\boldsymbol r).
$$

仿照式 (21) 定义 Green 函数的定解问题为

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_0)=-\delta(\boldsymbol r-\boldsymbol r_0),\qquad \boldsymbol r,\boldsymbol r_0\in D,
$$

$$
\left.\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial n}\right|_{\boldsymbol r\in S}=0.
$$

然后由式 (20) 可得

$$
u(\boldsymbol r_0)=\int_DG(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r)\,d\boldsymbol r
+\int_SG(\boldsymbol r,\boldsymbol r_0)\psi(\boldsymbol r)\,d\sigma.
$$

这样的做法是否有问题？分别从数学和物理上考虑这一问题。

上面推导式 (20) 时用到了第二 Green 公式，但 Green 函数在区域 $D$ 内的点 $\boldsymbol r_0$ 处有奇性，并不满足应用第二 Green 公式的条件，所以上述推导严格说来并不成立。下面以三维情况为例给出一个较为严格的推导，供读者参考。

为了避免 Green 函数在 $\boldsymbol r_0$ 处的奇性引起的困难，在区域 $D$ 中挖去以 $\boldsymbol r_0$ 为球心、$\varepsilon$ 为半径的小球 $K_\varepsilon$，对式 (19) 积分时，取积分区域为 $D-K_\varepsilon$，则右边第一项没有贡献，左边可以利用第二 Green 公式，于是得到

$$
\begin{aligned}
\int_{S_\varepsilon}\left[G(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}
-u(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial R}\right]d\sigma
={}&\int_{D-K_\varepsilon}G(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r)\,d\boldsymbol r\\
&+\int_S\left[G(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial n}
-u(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial n}\right]d\sigma,
\end{aligned}\tag{23}
$$

其中 $S_\varepsilon$ 是 $K_\varepsilon$ 的边界面（即球面 $|\boldsymbol r-\boldsymbol r_0|=\varepsilon$），作为区域 $D-K_\varepsilon$ 的边界面，其外法向方向与球的半径方向（即 $\boldsymbol R=\boldsymbol r-\boldsymbol r_0$ 方向）相反。令

$$
G(\boldsymbol r,\boldsymbol r_0)=\frac{1}{4\pi|\boldsymbol r-\boldsymbol r_0|}+g(\boldsymbol r,\boldsymbol r_0),\tag{24}
$$

则 $g(\boldsymbol r,\boldsymbol r_0)$ 满足定解问题

$$
\nabla^2g(\boldsymbol r,\boldsymbol r_0)=0,\qquad \boldsymbol r,\boldsymbol r_0\in D,\tag{25a}
$$

$$
g(\boldsymbol r,\boldsymbol r_0)|_{\boldsymbol r\in S}=-\frac{1}{4\pi|\boldsymbol r-\boldsymbol r_0|}.\tag{25b}
$$

既然 $g(\boldsymbol r,\boldsymbol r_0)$ 满足 Laplace 方程，它在 $\boldsymbol r_0$ 处应该没有奇性（从物理上看，如果把 $g(\boldsymbol r,\boldsymbol r_0)$ 看作静电势，则产生该电势的电荷均分布在区域 $D$ 的边界面和边界面以外，故静电势在 $D$ 内部没有奇性）。利用式 (24)，易得

$$
\begin{aligned}
\int_{S_\varepsilon}G(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}\,d\sigma
&=\frac{1}{4\pi\varepsilon}\int_{S_\varepsilon}\frac{\partial u(\boldsymbol r)}{\partial R}\,d\sigma
+\int_{S_\varepsilon}g(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}\,d\sigma\\
&=-\frac{1}{4\pi\varepsilon}\int_{K_\varepsilon}f(\boldsymbol r)\,d\boldsymbol r
+\int_{S_\varepsilon}g(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}\,d\sigma\\
&=-\frac13\varepsilon^2f(\boldsymbol r_1)+\int_{S_\varepsilon}g(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}\,d\sigma,
\end{aligned}
$$

其中第二步利用了 Gauss 定理和 Poisson 方程，第三步利用了积分中值定理，$\boldsymbol r_1$ 是球 $K_\varepsilon$ 内某点。类似有

$$
\begin{aligned}
\int_{S_\varepsilon}u(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial R}\,d\sigma
&=-\frac{1}{4\pi\varepsilon^2}\int_{S_\varepsilon}u(\boldsymbol r)\,d\sigma
+\int_{S_\varepsilon}u(\boldsymbol r)\frac{\partial g(\boldsymbol r,\boldsymbol r_0)}{\partial R}\,d\sigma\\
&=-u(\boldsymbol r^*)+\int_{S_\varepsilon}u(\boldsymbol r)\frac{\partial g(\boldsymbol r,\boldsymbol r_0)}{\partial R}\,d\sigma,
\end{aligned}
$$

其中第二步积分中值定理，$\boldsymbol r^*$ 是球面 $S_\varepsilon$ 上某点。于是

$$
\begin{aligned}
&\int_{S_\varepsilon}\left[G(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}
-u(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial R}\right]d\sigma\\
={}&u(\boldsymbol r^*)-\frac13\varepsilon^2f(\boldsymbol r_1)
+\int_{S_\varepsilon}\left[g(\boldsymbol r,\boldsymbol r_0)\frac{\partial u(\boldsymbol r)}{\partial R}
-u(\boldsymbol r)\frac{\partial g(\boldsymbol r,\boldsymbol r_0)}{\partial R}\right]d\sigma\\
={}&u(\boldsymbol r^*)-\frac13\varepsilon^2f(\boldsymbol r_1)
+\int_{K_\varepsilon}[g(\boldsymbol r,\boldsymbol r_0)\nabla^2u(\boldsymbol r)-u(\boldsymbol r)\nabla^2g(\boldsymbol r,\boldsymbol r_0)]\,d\boldsymbol r\\
={}&u(\boldsymbol r^*)-\frac13\varepsilon^2f(\boldsymbol r_1)-\frac{4\pi}{3}\varepsilon^3f(\boldsymbol r_2)g(\boldsymbol r_2,\boldsymbol r_0),
\end{aligned}\tag{26}
$$

其中第二步利用了第二 Green 公式，第三步利用了 $u(\boldsymbol r)$ 的 Poisson 方程、$g(\boldsymbol r,\boldsymbol r_0)$ 的 Laplace 方程和积分中值定理，$\boldsymbol r_2$ 是球 $K_\varepsilon$ 内某点。将上式代入 (23)，令 $\varepsilon\to0$，则上式中 $\boldsymbol r^*\to\boldsymbol r_0$，而后两项没有贡献，注意到 $G(\boldsymbol r,\boldsymbol r_0)f(\boldsymbol r)$ 的积分一般没有奇性，以及边界条件 (18b) 和 (21b)，就得到前面的结果 (22)。

### §2.4 *Laplace* 方程第一边值问题

考虑有界区域 $D$ 上的 Laplace 方程第一边值问题

$$
\nabla^2u(\boldsymbol r)=0,\qquad \boldsymbol r\in D,\tag{27a}
$$

$$
u(\boldsymbol r)|_S=\varphi(\boldsymbol r).\tag{27b}
$$

由于这只是 Poisson 方程第一边值问题 (18) 的一种特殊情况，所以其 Green 函数的定解问题仍然由式 (21) 定义，其解仍然由式 (22) 给出，将 $f(\boldsymbol r)=0$ 代入其中即得

$$
u(\boldsymbol r_0)=-\int_S\varphi(\boldsymbol r)\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial n}\,d\sigma.\tag{28}
$$

这一结果比较简单，但读者不必记忆它。

需要特别注意的是，尽管 Laplace 方程是无源的，但相应的 Green 函数的定解问题仍然是有源的，否则由于 Green 函数的边界条件是齐次的，将会得到 $G(\boldsymbol r,\boldsymbol r_0)=0$。读者只要记得 Green 函数又称为点源影响函数，就可以避免在这个问题上出错。

### §2.5 \*Green 函数的对称性

下面证明 Poisson 方程第一边值问题的 Green 函数对源点和场点具有对称性，即由定解问题 (21) 定义的 Green 函数满足

$$
G(\boldsymbol r_1,\boldsymbol r_2)=G(\boldsymbol r_2,\boldsymbol r_1).\tag{29}
$$

写下两个定解问题：

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_1)=-\delta(\boldsymbol r-\boldsymbol r_1),\qquad \boldsymbol r,\boldsymbol r_1\in D,\tag{30a}
$$

$$
G(\boldsymbol r,\boldsymbol r_1)|_{\boldsymbol r\in S}=0,\tag{30b}
$$

和

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_2)=-\delta(\boldsymbol r-\boldsymbol r_2),\qquad \boldsymbol r,\boldsymbol r_2\in D,\tag{31a}
$$

$$
G(\boldsymbol r,\boldsymbol r_2)|_{\boldsymbol r\in S}=0.\tag{31b}
$$

由方程 (30a) 和 (31a) 易得

$$
G(\boldsymbol r,\boldsymbol r_1)\nabla^2G(\boldsymbol r,\boldsymbol r_2)
-G(\boldsymbol r,\boldsymbol r_2)\nabla^2G(\boldsymbol r,\boldsymbol r_1)
=\delta(\boldsymbol r-\boldsymbol r_1)G(\boldsymbol r,\boldsymbol r_2)
-\delta(\boldsymbol r-\boldsymbol r_2)G(\boldsymbol r,\boldsymbol r_1),
$$

在区域 $D$ 上积分，利用第二 Green 公式，得到

$$
G(\boldsymbol r_1,\boldsymbol r_2)-G(\boldsymbol r_2,\boldsymbol r_1)
=\int_S\left[G(\boldsymbol r,\boldsymbol r_1)\frac{\partial G(\boldsymbol r,\boldsymbol r_2)}{\partial n}
-G(\boldsymbol r,\boldsymbol r_2)\frac{\partial G(\boldsymbol r,\boldsymbol r_1)}{\partial n}\right]d\sigma.\tag{32}
$$

再由边界条件 (30b) 和 (31b) 易知上式右边为零，于是得到式 (29)。

原则上，上面使用第二 Green 公式也是存在问题的。读者可以参照第三小节末尾的讨论给出一个比较严格的做法。

有了式 (29)，就可以将结果 (22) 改写为

$$
u(\boldsymbol r_0)=\int_DG(\boldsymbol r_0,\boldsymbol r)f(\boldsymbol r)\,d\boldsymbol r
-\int_S\varphi(\boldsymbol r)\frac{\partial G(\boldsymbol r_0,\boldsymbol r)}{\partial n}\,d\sigma.\tag{33}
$$

这样就是对源点积分，物理意义更加清楚。

对 Poisson 方程第二边值问题，以及 Helmholtz 方程的三类边值问题，都可以证明其 Green 函数的对称性。Poisson 方程第二边值问题的 Green 函数不存在。对其它不含时问题，Green 函数的对称性需要具体讨论。

## §3 特殊区域上 Green 函数的制作

由上节的讨论可知，Green 函数法的关键在于求解 Green 函数的定解问题。前已指出，求解 Green 函数的定解问题，通常并不比求解一般的定解问题来得容易，但对于一些特殊区域上的 Poisson 方程（或 Laplace 方程）第一边值问题，存在求解 Green 函数的特殊技巧，即镜像法。本节讨论几个典型的例子。

![Figure 1](...)

图 1：源电荷与像电荷不能位于球心的两侧

### §3.1 球内 *Laplace* 方程第一边值问题

考虑球内 Laplace 方程第一边值问题

$$
\nabla^2u(\boldsymbol r)=0,\qquad r<a,\tag{34a}
$$

$$
u(\boldsymbol r)|_{r=a}=f(\theta,\phi).\tag{34b}
$$

**思考**　如果没有 Green 函数法，如何求解这一问题？

相应的 Green 函数的定解问题为

$$
\nabla^2G(\boldsymbol r,\boldsymbol r_0)=-\delta(\boldsymbol r-\boldsymbol r_0),\qquad r,r_0<a,\tag{35a}
$$

$$
G(\boldsymbol r,\boldsymbol r_0)|_{r=a}=0.\tag{35b}
$$

如果将上式看作静电场问题，则可以表述为：球内某点 $\boldsymbol r_0$ 处有一个点电荷，其电量为 $\varepsilon_0$，而球面上的电势为零，求解球内的电势分布。为了保持球面上的电势为零，则球面上或球外必须有某种电荷分布，以便抵消 $\boldsymbol r_0$ 处的点电荷在球面上产生的电势。现在我们希望在球外某点 $\boldsymbol r'_0$ 处放置一个点电荷，其电量为 $-q\varepsilon_0$，其中参数 $q$ 有待确定，使得两个点电荷在球面上产生的电势为零，由于数学上可以证明有唯一性定理，那么它们在球内产生的电势就是所求的解。这种方法称为镜像法，在电动力学中也称为电像法。$\boldsymbol r'_0$ 处放置的点电荷称为像电荷。对于某些区域，如果一个像电荷不能解决问题，可以考虑用多个像电荷。当然，镜像法是否可行，取决于区域的形状。通常它只对某些非常规则的区域才能奏效。

具体来说，对于定解问题 (35)，如果镜像法能够奏效，那么根据对称性，像电荷的位置必定在通过球心和点 $\boldsymbol r_0$ 的直线上。而且像电荷与源电荷必定在球心的同一侧。因为，如果像电荷与源电荷在球心的两侧，如图 1 所示，那么源电荷的电势在 $A$ 点较大（指绝对值，下同），在 $B$ 点较小，而像电荷的电势在 $A$ 点较小，在 $B$ 点较大，所以两者的电势在 $A$、$B$ 两点不可能抵消，从而球面上电势不能为零。

现在像电荷的位置 $\boldsymbol r'_0$ 在过 $\boldsymbol r_0$ 的半径的延长线上，如图 2 所示，所以它与 $\boldsymbol r_0$ 具有相同的角向坐标。如果 $\boldsymbol r_0$ 的球坐标为 $(r_0,\theta_0,\phi_0)$，则 $\boldsymbol r'_0$ 的球坐标为 $(r'_0,\theta_0,\phi_0)$。在球内任取一点 $\boldsymbol r$，其球坐标为 $(r,\theta,\phi)$，$\boldsymbol r$ 与 $\boldsymbol r_0$ 的夹角记作 $\psi$。记 $\boldsymbol R=\boldsymbol r-\boldsymbol r_0$，$\boldsymbol R'=\boldsymbol r-\boldsymbol r'_0$。源电荷与像电荷在 $\boldsymbol r$ 处产生的电势为

![Figure 2](...)

图 2：球内区域的镜像法

$$
g(\boldsymbol r,\boldsymbol r_0)=\frac{1}{4\pi R}-\frac{q}{4\pi R'}
=\frac{1}{4\pi\sqrt{r^2+r_0^2-2rr_0\cos\psi}}
-\frac{q}{4\pi\sqrt{r^2+r_0'^2-2rr_0'\cos\psi}}.\tag{36}
$$

当然它还依赖于 $r'_0$ 和 $q$ 两个参数。在球面上，$r=a$，则

$$
g(\boldsymbol r,\boldsymbol r_0)|_{r=a}
=\frac{1}{4\pi r_0\sqrt{1+(a/r_0)^2-2(a/r_0)\cos\psi}}
-\frac{q}{4\pi a\sqrt{1+(r'_0/a)^2-2(r'_0/a)\cos\psi}}.\tag{37}
$$

注意到 $a/r_0$ 和 $r'_0/a$ 均大于 1，由上式可见，只要取

$$
r'_0=\frac{a^2}{r_0},\qquad q=\frac{a}{r_0},\tag{38}
$$

则 $g(\boldsymbol r,\boldsymbol r_0)|_{r=a}=0$，可见镜像法对于本问题确实是成功的。将上式所得参数值代入 $g(\boldsymbol r,\boldsymbol r_0)$，就得到所求的 Green 函数为

$$
G(\boldsymbol r,\boldsymbol r_0)=\frac{1}{4\pi R}-\frac{a}{r_0}\frac{1}{4\pi R'},\tag{39a}
$$

其中

$$
R=\sqrt{r^2+r_0^2-2rr_0\cos\psi},\qquad
R'=\sqrt{r^2+\frac{a^4}{r_0^2}-2r\frac{a^2}{r_0}\cos\psi},\tag{39b}
$$

而

$$
\cos\psi=\cos\theta\cos\theta_0+\sin\theta\sin\theta_0\cos(\phi-\phi_0).\tag{39c}
$$

最后一式可以这样计算：$\cos\psi=\hat{\boldsymbol r}\cdot\hat{\boldsymbol r}_0$，其中 $\hat{\boldsymbol r}$ 和 $\hat{\boldsymbol r}_0$ 分别是 $\boldsymbol r$ 和 $\boldsymbol r_0$ 方向的单位矢量，所以它们可用直角坐标分量表示为 $\hat{\boldsymbol r}=(\sin\theta\cos\phi,\sin\theta\sin\phi,\cos\theta)$ 和 $\hat{\boldsymbol r}_0=(\sin\theta_0\cos\phi_0,\sin\theta_0\sin\phi_0,\cos\theta_0)$，由此可得结果。

**注**　现在看看像电荷的位置和大小。由前面的讨论和式 (38) 可知，像电荷与源电荷具有同样的角向坐标，而径向坐标满足 $r'_0r_0=a^2$，这样的两点称为关于球面的对称点，所以位置关系很容易理解和掌握。由式 (38) 又可知，像电荷与源电荷之比为 $-a/r_0$，负号容易理解，因为两个同号的电荷在球面上的总电势不可能为零。另外，$a/r_0$ 是一个大于 1 的无量纲数，这只要看图 2 的 $B$ 点就容易理解，因为源电荷离该点较近，而像电荷离该点较远，所以像电荷应该较大（指绝对值），该点的电势才能为零。这些考虑可以帮助读者记忆以上结果，但读者应该首先对结果的推导有充分的理解，如果能够自己推导一遍就更好了。

有了 Green 函数，就可以写出定解问题 (34) 的解。注意球面的外法向方向就是半径方向，经过一些运算，可以得到

$$
\left.\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial n}\right|_{r=a}
=\left.\frac{\partial G(\boldsymbol r,\boldsymbol r_0)}{\partial r}\right|_{r=a}
=\left.\frac{r_0^2-a^2}{4\pi aR^3}\right|_{r=a}
=\frac{r_0^2-a^2}{4\pi a(a^2+r_0^2-2ar_0\cos\psi)^{3/2}}.\tag{40}
$$

代入式 (28)，注意到 $\varphi(\boldsymbol r)=f(\theta,\phi)$ 和 $d\sigma=a^2\sin\theta\,d\theta\,d\phi$，就得到

$$
u(r_0,\theta_0,\phi_0)=\frac{a(a^2-r_0^2)}{4\pi}
\int\frac{f(\theta,\phi)\sin\theta\,d\theta\,d\phi}
{(a^2+r_0^2-2ar_0\cos\psi)^{3/2}}.\tag{41}
$$

其中 $\theta$ 的积分范围是 $0$ 到 $\pi$，$\phi$ 的积分范围是 $0$ 到 $2\pi$。上式称为球的 Poisson 积分。可以证明，它与用分离变量法得到的结果是一致的。

**注**　按式 (39c)，$\cos\psi$ 中包含变量 $\theta$、$\phi$ 和参数 $\theta_0$、$\phi_0$，上式对 $\theta$、$\phi$ 积分后，最后是 $r_0$、$\theta_0$ 和 $\phi_0$ 的函数。注意后者也是可以任意变化的，所以在求 Green 函数时，不能将 $r_0$ 置于特殊的位置，比如极轴上。

### §3.2 上半平面 *Laplace* 方程第一边值问题

考虑上半平面 Laplace 方程第一边值问题

$$
\nabla^2u(\boldsymbol\rho)=0,\qquad -\infty<x<+\infty,\quad 0<y<+\infty,\tag{42a}
$$

$$
u(\boldsymbol\rho)|_{y=0}=\varphi(x).\tag{42b}
$$

**思考**　如果没有 Green 函数法，如何求解这一问题？

相应的 Green 函数的定解问题为

$$
\nabla^2G(\boldsymbol\rho,\boldsymbol\rho_0)=-\delta(\boldsymbol\rho-\boldsymbol\rho_0),\qquad -\infty<x,x_0<+\infty,\quad 0<y,y_0<+\infty,\tag{43a}
$$

$$
G(\boldsymbol\rho,\boldsymbol\rho_0)|_{y=0}=0.\tag{43b}
$$

如果将上式看作静电场问题，则可以表述为：上半平面上某点 $\boldsymbol\rho_0=(x_0,y_0)$（括号中是其直角坐标，下同）处有一个二维点电荷，其电量为 $\varepsilon_0$（从三维空间看，实际上是垂直于 $xy$ 平面的无穷长线电荷，其线密度为 $\varepsilon_0$），而 $x$ 轴上的电势为零，求解上半平面的电势分布。为了保持 $x$ 轴上的电势为零，显然应该在下半平面的点 $\boldsymbol\rho'_0=(x_0,-y_0)$ 处放置一个二维点电荷，其电量为 $-\varepsilon_0$，如图 3 所示。它们在上半平面上产生的电势就是所求的 Green 函数：

$$
G(\boldsymbol\rho,\boldsymbol\rho_0)
=-\frac{1}{2\pi}\ln|\boldsymbol\rho-\boldsymbol\rho_0|
+\frac{1}{2\pi}\ln|\boldsymbol\rho-\boldsymbol\rho'_0|
=\frac{1}{4\pi}\ln\frac{(x-x_0)^2+(y+y_0)^2}{(x-x_0)^2+(y-y_0)^2}.\tag{44}
$$

这里正好不需要常数项。

有了 Green 函数，就可以写出定解问题 (42) 的解。注意上半平面的外法向方向是 $-y$ 方向，经过简单的运算，可以得到

$$
\left.\frac{\partial G(\boldsymbol\rho,\boldsymbol\rho_0)}{\partial n}\right|_{y=0}
=-\left.\frac{\partial G(\boldsymbol\rho,\boldsymbol\rho_0)}{\partial y}\right|_{y=0}
=-\frac{y_0}{\pi}\frac{1}{(x-x_0)^2+y_0^2}.\tag{45}
$$

![Figure 3](...)

图 3：上半平面的镜像法

代入式 (28)，注意到 $\varphi(\boldsymbol r)$ 现在就是 $\varphi(x)$，而 $d\sigma=dx$，就得到

$$
u(x_0,y_0)=\frac{y_0}{\pi}\int_{-\infty}^{+\infty}\frac{\varphi(x)}{(x_0-x)^2+y_0^2}\,dx.\tag{46}
$$

**思考**　在以上结果中代入 $y_0=0$，似乎得到 $u(x_0,0)=0$，这与边界条件 (42b) 矛盾，事实果真如此吗？如何由上式得到边界条件 (42b)？

**提示**　证明当 $y_0\to0^+$，$y_0/\pi[(x_0-x)^2+y_0^2]\to\delta(x_0-x)$。

### §3.3 圆内 *Laplace* 方程第一边值问题

考虑平面上圆内 Laplace 方程第一边值问题

$$
\nabla^2u(\boldsymbol\rho)=0,\qquad \rho<a,\tag{47a}
$$

$$
u(\boldsymbol\rho)|_{\rho=a}=f(\phi).\tag{47b}
$$

**思考**　如果没有 Green 函数法，如何求解这一问题？

相应的 Green 函数的定解问题为

$$
\nabla^2G(\boldsymbol\rho,\boldsymbol\rho_0)=-\delta(\boldsymbol\rho-\boldsymbol\rho_0),\qquad \rho,\rho_0<a,\tag{48a}
$$

$$
G(\boldsymbol\rho,\boldsymbol\rho_0)|_{\rho=a}=0.\tag{48b}
$$

如果将上式看作静电场问题，则可以表述为：圆内某点 $\boldsymbol\rho_0=(\rho_0,\phi_0)$（括号中是其极坐标，下同）处有一个二维点电荷，其电量为 $\varepsilon_0$（从三维空间看，实际上是垂直于 $xy$ 平面的无穷长线电荷，其线密度为 $\varepsilon_0$），而圆周上的电势为零，求解圆内的电势分布。今在 $\boldsymbol\rho_0$ 点关于圆周的对称点 $\boldsymbol\rho'_0=(a^2/\rho_0^2)\boldsymbol\rho_0=(a^2/\rho_0,\phi_0)$ 处放置一个二维点电荷，其电量为 $-\varepsilon_0$，如图 4 所示。在圆内任取一点 $\boldsymbol\rho=(\rho,\phi)$，$\boldsymbol\rho$ 与 $\boldsymbol\rho_0$ 的夹角为 $\psi=\phi-\phi_0$。记 $\boldsymbol R=\boldsymbol\rho-\boldsymbol\rho_0$，$\boldsymbol R'=\boldsymbol\rho-\boldsymbol\rho'_0$。源电荷与像电荷在 $\boldsymbol\rho$ 处产生的电势为

$$
g(\boldsymbol\rho,\boldsymbol\rho_0)
=-\frac{1}{2\pi}\ln R+\frac{1}{2\pi}\ln R'+c
=\frac{1}{4\pi}\ln\frac{\rho^2+a^4/\rho_0^2-2\rho(a^2/\rho_0)\cos\psi}{\rho^2+\rho_0^2-2\rho\rho_0\cos\psi}+c.\tag{49}
$$

![Figure 4](...)

图 4：圆内区域的镜像法

注意其中有一个常数项，它依赖于电势零点的选择，这是二维问题的特点。在圆周上，$\rho=a$，易得

$$
g(\boldsymbol\rho,\boldsymbol\rho_0)|_{\rho=a}=\frac{1}{2\pi}\ln\frac{a}{\rho_0}+c.\tag{50}
$$

由上式可见，只要取

$$
c=-\frac{1}{2\pi}\ln\frac{a}{\rho_0},\tag{51}
$$

则 $g(\boldsymbol\rho,\boldsymbol\rho_0)|_{\rho=a}=0$。将上式所得 $c$ 值代入 $g(\boldsymbol\rho,\boldsymbol\rho_0)$，就得到所求的 Green 函数为

$$
G(\boldsymbol\rho,\boldsymbol\rho_0)
=\frac{1}{2\pi}\ln\frac{1}{R}-\frac{1}{2\pi}\ln\left(\frac{a}{\rho_0}\frac{1}{R'}\right)
=\frac{1}{4\pi}\ln\frac{\rho_0^2\rho^2/a^2+a^2-2\rho\rho_0\cos\psi}{\rho^2+\rho_0^2-2\rho\rho_0\cos\psi}.\tag{52}
$$

上式第一步的写法形式上与式 (39a) 类似，可以帮助读者把握这一结果。

**注**　在本问题中，我们当然也可以象第一小节处理球内区域那样，只是假定 $\boldsymbol\rho'_0$ 位于过 $\boldsymbol\rho_0$ 的半径的延长线上，而像电荷的大小为 $-q\varepsilon_0$，然后经过计算确定像电荷的具体位置和参数 $q$ 的大小。不过，用镜像法求解 Green 函数，本来就带有猜测和尝试的成分，所以有了一定的经验以后，在安排像电荷的位置和大小时，可以尽量减少不确定因素，只要最后的结果确实满足方程（只要源电荷大小正确，而像电荷均在所考虑的区域以外即可）和边界条件（可以直接验证）就可以了。极端地说，如果有办法，一步写出结果也未尝不可。

有了 Green 函数，就可以写出定解问题 (47) 的解。注意圆周的外法向方向就是半径方向，经过一些运算，可以得到

$$
\left.\frac{\partial G(\boldsymbol\rho,\boldsymbol\rho_0)}{\partial n}\right|_{\rho=a}
=\left.\frac{\partial G(\boldsymbol\rho,\boldsymbol\rho_0)}{\partial\rho}\right|_{\rho=a}
=\frac{\rho_0^2-a^2}{2\pi a[a^2+\rho_0^2-2a\rho_0\cos(\phi-\phi_0)]}.\tag{53}
$$

代入式 (28)，注意到 $\varphi(\boldsymbol r)$ 现在是 $f(\phi)$，而 $d\sigma=a\,d\phi$，就得到

$$
u(\rho_0,\phi_0)=\frac{a^2-\rho_0^2}{2\pi}
\int_0^{2\pi}\frac{f(\phi)\,d\phi}{a^2+\rho_0^2-2a\rho_0\cos(\phi-\phi_0)}.\tag{54}
$$

上式称为圆的 Poisson 积分。可以证明，它与用分离变量法得到的结果是一致的。

**习题**　对三维半无界空间的定解问题

$$
\nabla^2u(\boldsymbol r)=0,\qquad -\infty<x,y<+\infty,\quad 0<z<+\infty,
$$

$$
u(\boldsymbol r)|_{z=0}=\varphi(x,y).
$$

(1) 写出相应的 Green 函数的定解问题。

(2) 求出这一 Green 函数。

(3) 求出 $u(\boldsymbol r)$ 的积分公式，即用 Green 函数和定解条件表示出 $u(\boldsymbol r)$。

## §4 \*波动方程的 Green 函数法

暂略。

## §5 \*热传导方程的 Green 函数法

暂略。

## 补充习题

1. 对于 Poisson 方程第三边值问题

$$
\nabla^2u(\boldsymbol r)=-f(\boldsymbol r),\qquad \boldsymbol r\in D,
$$

$$
\left.\left(\alpha u+\beta\frac{\partial u}{\partial n}\right)\right|_{\boldsymbol r\in\partial D}=\varphi(\boldsymbol r).
$$

定义相应的 Green 函数，并给出解的积分公式。

2. 对三维空间的定解问题

$$
\nabla^2u(\boldsymbol r)=0,\qquad 0<x,y<+\infty,\quad -\infty<z<+\infty,
$$

$$
u(\boldsymbol r)|_{x=0}=\varphi_1(y,z),\qquad u(\boldsymbol r)|_{y=0}=\varphi_2(x,z),
$$

写出相应的 Green 函数的定解问题并求出该 Green 函数。

3. 对三维空间的定解问题

$$
\nabla^2u(\boldsymbol r)=0,\qquad -\infty<x,y<+\infty,\quad 0<z<h,
$$

$$
u(\boldsymbol r)|_{z=0}=\varphi_1(x,y),\qquad u(\boldsymbol r)|_{z=h}=\varphi_2(x,y),
$$

写出相应的 Green 函数的定解问题并求出该 Green 函数。

4. 对于平面上的定解问题

$$
\nabla^2u(\boldsymbol\rho)=0,\qquad 0<\rho<\infty,\quad 0<\phi<\pi/4,
$$

$$
u(\boldsymbol\rho)|_{\phi=0}=f_1(\rho),\qquad u(\boldsymbol\rho)|_{\phi=\pi/4}=f_2(\rho),
$$

写出相应的 Green 函数的定解问题并求出该 Green 函数。

5. 对于平面上的定解问题

$$
\nabla^2u(\boldsymbol\rho)=0,\qquad 0<\rho<\infty,\quad 0<\phi<2\pi/n,\qquad n\in\mathbb N^+,
$$

$$
u(\boldsymbol\rho)|_{\phi=0}=f_1(\rho),\qquad u(\boldsymbol\rho)|_{\phi=2\pi/n}=f_2(\rho),
$$

写出相应的 Green 函数的定解问题并求出该 Green 函数。
