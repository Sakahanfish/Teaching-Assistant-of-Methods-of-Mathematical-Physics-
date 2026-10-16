# 第十九章 Green 函数方法

在第十章中，已经讨论过常微分方程 Green 函数的定义、对称性质及其求法。本章将继续这一话题，讨论偏微分方程定解问题 Green 函数的概念、对称性质、Green 函数的常用求法以及偏微分方程定解问题的 Green 函数解法。作为必需的预备知识，读者应该熟悉与掌握第十章中有关 $\delta$ 函数以及常微分方程 Green 函数的相关知识。

在本章中，经常要用到 Green 第二公式（或简称 Green 公式）

$$
\iiint_V
\left[
u(\boldsymbol r)\nabla^2v(\boldsymbol r)
-
v(\boldsymbol r)\nabla^2u(\boldsymbol r)
\right]
\,d\boldsymbol r
=
\iint_\Sigma
(u\nabla v-v\nabla u)\cdot d\boldsymbol\Sigma,
\tag{19.1a}
$$

其中 $u(\boldsymbol r)\equiv u(x,y,z)$，$v(\boldsymbol r)\equiv v(x,y,z)$，$d\boldsymbol r=dxdydz$，$\Sigma$ 是 $V$ 的边界面，并且规定边界面的外法线方向为正。这个公式是 Green 第一公式

$$
\iiint_V
u(\boldsymbol r)\nabla^2v(\boldsymbol r)
\,d\boldsymbol r
=
\iint_\Sigma
u\nabla v\cdot d\boldsymbol\Sigma
-
\iiint_V
\nabla u\cdot\nabla v
\,d\boldsymbol r
\tag{19.1b}
$$

的直接推论。有关这两个公式的成立条件及证明从略，读者可参阅高等数学中的有关章节。

## §19.1 Green 函数的概念

以静电场为例。设在无界空间中有一定的电荷分布，电荷密度为 $\rho(\boldsymbol r)$。这样，位于 $\boldsymbol r'=(x',y',z')$ 处的体元 $d\boldsymbol r'$ 内的电量即为 $\rho(\boldsymbol r')d\boldsymbol r'$，它在空间 $\boldsymbol r=(x,y,z)$ 点的电势是

$$
\frac1{4\pi\varepsilon_0}
\frac{\rho(\boldsymbol r')}
{|\boldsymbol r-\boldsymbol r'|}
\,d\boldsymbol r'.
$$

根据电势叠加原理，把空间中全部电荷产生的电势叠加起来，就得到在 $\boldsymbol r$ 点的总电势为

$$
\phi(\boldsymbol r)
=
\frac1{4\pi\varepsilon_0}
\iiint
\frac{\rho(\boldsymbol r')}
{|\boldsymbol r-\boldsymbol r'|}
\,d\boldsymbol r'.
\tag{19.2}
$$

这个结果说明，只要知道了单位点电荷在空间的电势分布，通过电荷的分割与叠加，就可以得到任意电荷分布时的电势。其实，这种做法只不过是利用了线性偏微分方程（以及隐含的边界条件）的线性性质。

这个例子的简单之处是无界空间。如果是有界空间，尽管仍然可以把空间内的电荷无限分割，但由于边界的存在，在边界面上会有（单层或偶极层的）感生面电荷分布，也需要将这些面电荷无限分割。这样，在有界空间的情形下，我们面临的问题就是：如何通过（适当边界条件下的）点电荷电势的叠加，给出任意电荷分布和任意边界条件时的电势。这就是说，要用定解问题

$$
\nabla^2G(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\qquad
\boldsymbol r,\boldsymbol r'\in V,
\tag{19.3a}
$$

$$
\text{适当的边界条件}
\tag{19.3b}
$$

的解 $G(\boldsymbol r;\boldsymbol r')$ 叠加出

$$
\nabla^2u(\boldsymbol r)
=
-\frac1{\varepsilon_0}\rho(\boldsymbol r),
\qquad
\boldsymbol r\in V,
\tag{19.4a}
$$

$$
\left.u\right|_\Sigma
=
f(\Sigma)
\tag{19.4b}
$$

的解 $u(\boldsymbol r)$，即把 $u(\boldsymbol r)$ 用 $\rho(\boldsymbol r)$、$f(\Sigma)$ 以及 $G(\boldsymbol r;\boldsymbol r')$ 表示出来。

为此，将方程 (19.3a) 和 (19.4a) 分别乘以 $u(\boldsymbol r)$ 和 $G(\boldsymbol r;\boldsymbol r')$，相减，再在空间 $V$ 内积分，得

$$
\iiint_V
\left[
u(\boldsymbol r)\nabla^2G(\boldsymbol r;\boldsymbol r')
-
G(\boldsymbol r;\boldsymbol r')\nabla^2u(\boldsymbol r)
\right]
d\boldsymbol r
=
-\frac1{\varepsilon_0}
\left[
u(\boldsymbol r')
-
\iiint_V
G(\boldsymbol r;\boldsymbol r')
\rho(\boldsymbol r)
\,d\boldsymbol r
\right].
$$

根据 Green 公式，可以将上式左端的体积分化为面积分。经过移项、整理，就有

$$
u(\boldsymbol r')
=
\iiint_V
G(\boldsymbol r;\boldsymbol r')
\rho(\boldsymbol r)
\,d\boldsymbol r
-
\varepsilon_0
\iint_\Sigma
\left[
u(\boldsymbol r)\nabla G(\boldsymbol r;\boldsymbol r')
-
G(\boldsymbol r;\boldsymbol r')\nabla u(\boldsymbol r)
\right]\cdot d\boldsymbol\Sigma.
$$

在上面的面积分中，第一项 $u(\boldsymbol r)$ 在边界面 $\Sigma$ 上的数值由边界条件 (19.4b) 给出，是已知的；$G(\boldsymbol r;\boldsymbol r')$ 可由定解问题 (19.3) 求出，故其梯度 $\nabla G(\boldsymbol r;\boldsymbol r')$ 在边界面上的数值也可求；第二项中，$\nabla u(\boldsymbol r)$ 在边界面上的数值未知，所以，为了要能够完成把 $u(\boldsymbol r)$ 用 $\rho(\boldsymbol r)$、$f(\Sigma)$ 以及 $G(\boldsymbol r;\boldsymbol r')$ 表示这一要求，即不允许出现 $\left.\nabla u(\boldsymbol r)\right|_\Sigma$ 这个未知函数，就必须对 $G(\boldsymbol r;\boldsymbol r')$ 加上齐次边界条件

$$
\left.G(\boldsymbol r;\boldsymbol r')\right|_\Sigma
=
0.
\tag{19.3b}
$$

最后就得到

$$
u(\boldsymbol r')
=
\iiint_V
G(\boldsymbol r;\boldsymbol r')
\rho(\boldsymbol r)
\,d\boldsymbol r
-
\varepsilon_0
\iint_\Sigma
f(\Sigma)
\left.
\nabla G(\boldsymbol r;\boldsymbol r')
\right|_\Sigma
\cdot d\boldsymbol\Sigma,
\tag{19.5}
$$

或者把 $\boldsymbol r$ 和 $\boldsymbol r'$ 对换一下，

$$
u(\boldsymbol r)
=
\iiint_V
G(\boldsymbol r';\boldsymbol r)
\rho(\boldsymbol r')
\,d\boldsymbol r'
-
\varepsilon_0
\iint_\Sigma
f(\Sigma')
\left.
\nabla'G(\boldsymbol r';\boldsymbol r)
\right|_{\Sigma'}
\cdot d\boldsymbol\Sigma'
\tag{19.5'}
$$

$$
=
\iiint_V
G(\boldsymbol r';\boldsymbol r)
\rho(\boldsymbol r')
\,d\boldsymbol r'
-
\varepsilon_0
\iint_\Sigma
f(\Sigma')
\left.
\frac{\partial G(\boldsymbol r';\boldsymbol r)}
{\partial n'}
\right|_{\Sigma'}
d\Sigma',
\tag{19.5''}
$$

其中的 $\nabla'$ 和 $\partial/\partial n'$ 表示对自变量 $\boldsymbol r'$ 微商。

熟悉 Green 公式的读者会对上面的做法提出质疑：$G(\boldsymbol r;\boldsymbol r')$ 在 $\boldsymbol r=\boldsymbol r'$ 点肯定不连续，根本不能应用 Green 公式。为了克服这一缺陷，可以将方程 (19.3a) 右端的电荷密度函数改为足够好的连续函数，它在一定尺度内明显不为 0，而总电量保持为 1 个单位（因而可以应用 Green 公式），在求出了 $u(\boldsymbol r)$ 的表达式后再令这一尺度趋于 0。这种做法正是重复了引进 $\delta$ 函数的极限过程。第十章曾经指出，引入 $\delta$ 函数的好处恰恰就在于略去这种极限过程，恰恰就在于可以把 $\delta$ 函数当成连续函数来处理。因此，上面的做法是严格的、正确的。另外一种做法是在点电荷所在的 $\boldsymbol r'$ 点的附近挖去一个小体积，在这个新的空间区域中应用 Green 公式（必须注意，现在的边界面除了原来的表面 $\Sigma$ 之外，还有新增加的在 $\boldsymbol r'$ 点处的界面），然后再令这个小体积趋于 0。毫无疑问，这样得到的结果和 (19.5') 或 (19.5'') 式完全一致。

以上通过静电场的实例引入了 Poisson 方程在第一类边界条件下（简称 Poisson 方程的第一边值问题）的 Green 函数。简言之，所谓 Green 函数就是单位点电荷在齐次边界条件下的电势。对于第二类或第三类边界条件，原则上也可以类似地讨论。从数学上说，不含时间（稳定问题）的偏微分方程（Laplace 方程，Poisson 方程，Helmholtz 方程等）在一定边界条件下的 Green 函数就可以定义为一个特殊的定解问题的解：

(1) 方程和原来定解问题的方程一样，只是非齐次项改为 $\delta$ 函数（点源）；

(2) 同种类型的齐次边界条件。

但是在某些特殊情形下，这样定义的 Green 函数可能无解。例如对于 Poisson 方程定解问题 (19.4)，若边界条件 (19.4b) 改为

$$
\left.
\frac{\partial u(\boldsymbol r)}
{\partial n}
\right|_\Sigma
=
f(\Sigma),
\tag{19.4b'}
$$

则按照上面的讨论，Green 函数 $G(\boldsymbol r;\boldsymbol r')$ 在边界面上应当满足齐次的第二类边界条件

$$
\left.
\frac{\partial G(\boldsymbol r;\boldsymbol r')}
{\partial n}
\right|_\Sigma
=
0.
\tag{19.4b''}
$$

在 Green 公式 (19.1a) 中，令 $u(\boldsymbol r)=1$，$v(\boldsymbol r)=G(\boldsymbol r;\boldsymbol r')$，就应该有

$$
\iiint_V
\nabla^2G(\boldsymbol r;\boldsymbol r')
\,d\boldsymbol r
=
\iint_\Sigma
\nabla G(\boldsymbol r;\boldsymbol r')
\cdot d\boldsymbol\Sigma
=
\iint_\Sigma
\frac{\partial G(\boldsymbol r;\boldsymbol r')}
{\partial n}
\,d\Sigma
=
0.
$$

可是，将方程 (19.3a) 积分，又得到

$$
\iiint_V
\nabla^2G(\boldsymbol r;\boldsymbol r')
\,d\boldsymbol r
=
-\frac1{\varepsilon_0},
\qquad
\text{即}
\qquad
\iint_\Sigma
\frac{\partial G(\boldsymbol r;\boldsymbol r')}
{\partial n}
\,d\Sigma
=
-\frac1{\varepsilon_0},
$$

显然和边界条件 (19.4b'') 矛盾。这说明，在齐次边界条件 (19.4b'') 下，方程 (19.4a) 一定无解，换句话说，这样的 Green 函数一定不存在。在这种情形下，需要引进广义的 Green 函数。有兴趣的读者请参阅参考书目 [1] 的第 17.4 节。这里从略。

## §19.2 稳定问题 Green 函数的一般性质

建立了稳定问题的 Green 函数概念之后，就需要讨论它的一般性质：Green 函数在点源附近的行为及 Green 函数的对称性。

不妨仍然用静电场的语言来描述 Poisson 方程第一边值问题的 Green 函数。电磁学的知识告诉我们，在空间 $V$ 中的点电荷，必然要在边界面上产生一定的感生（面）电荷分布。因此，决定 Green 函数的定解问题又可以（在 $V$ 内）等价地写成无界空间中的 Poisson 方程

$$
\nabla^2G(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\left[
\delta(\boldsymbol r-\boldsymbol r')
+
\sigma(\Sigma)\delta_\Sigma
\right],
\tag{19.6}
$$

其中 $\sigma(\Sigma)$ 是边界面 $\Sigma$ 上的感生面电荷密度，$\delta_\Sigma$ 代表某种 $\delta$ 函数（例如边界面为球面时就是 $\delta(r-a)$）。相应地，（定义在 $V$ 内的）Green 函数 $G(\boldsymbol r;\boldsymbol r')$ 就应该是这两部分电荷电势的叠加：单位点电荷 $\delta(\boldsymbol r-\boldsymbol r')$ 的电势 $G_0(\boldsymbol r;\boldsymbol r')$ 和边界面上的感生电荷 $\sigma(\Sigma)$ 的电势 $g(\boldsymbol r;\boldsymbol r')$，

$$
G(\boldsymbol r;\boldsymbol r')
=
G_0(\boldsymbol r;\boldsymbol r')
+
g(\boldsymbol r;\boldsymbol r'),
\tag{19.7}
$$

$$
\nabla^2G_0(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\tag{19.8}
$$

$$
\nabla^2g(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\sigma(\Sigma)\delta_\Sigma.
\tag{19.9}
$$

显然，方程 (19.8) 的解为

$$
G_0(\boldsymbol r;\boldsymbol r')
=
\frac1{4\pi\varepsilon_0}
\frac1{|\boldsymbol r-\boldsymbol r'|},
\tag{19.10}
$$

所以 $G_0(\boldsymbol r;\boldsymbol r')$ 在 $\boldsymbol r=\boldsymbol r'$ 点不连续。对于方程 (19.9)，因为感生电荷 $\sigma(\Sigma)$ 只分布在曲面 $\Sigma$ 上，所以 $g(\boldsymbol r;\boldsymbol r')$ 及其一阶偏导数在曲面 $\Sigma$ 外（特别是，在 $V$ 内）处处连续。

对于第三类边界条件，也有同样的结果，只不过 $g(\boldsymbol r;\boldsymbol r')$ 的具体表达式会有所不同。

对于其他类型的稳定问题，也可以类似地讨论。例如，Helmholtz 方程的 Green 函数就是

$$
\nabla^2\widetilde G(\boldsymbol r;\boldsymbol r')
+
k^2\widetilde G(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\qquad
\boldsymbol r,\boldsymbol r'\in V,
\tag{19.11a}
$$

$$
\left.
\widetilde G(\boldsymbol r;\boldsymbol r')
\right|_\Sigma
=
0
\tag{19.11b}
$$

的解。可以看出，除了 $\boldsymbol r=\boldsymbol r'$ 点外，$\widetilde G(\boldsymbol r;\boldsymbol r')$ 在 $V$ 内是处处连续的。而且，和 Poisson 方程的 Green 函数 $G(\boldsymbol r;\boldsymbol r')$ 所满足的定解问题 (19.3) 相比，齐次边界条件完全相同，不同之处只是在方程 (19.11a) 中多了一项 $k^2\widetilde G(\boldsymbol r;\boldsymbol r')$。因为 Green 函数的奇异性完全是由微分方程决定的，而且乘法运算不可能改变函数的奇异性，所以我们可以断定，$\widetilde G(\boldsymbol r;\boldsymbol r')$ 和 $G(\boldsymbol r;\boldsymbol r')$ 一样，在 $\boldsymbol r=\boldsymbol r'$ 点都是以 $1/|\boldsymbol r-\boldsymbol r'|$ 的形式发散。事实上，从后面 §19.3 的讨论可知，在 $\boldsymbol r=\boldsymbol r'$ 点附近，一定有

$$
\widetilde G(\boldsymbol r;\boldsymbol r')
\sim
\frac1{4\pi\varepsilon_0}
\frac{\cos(k|\boldsymbol r-\boldsymbol r'|)}
{|\boldsymbol r-\boldsymbol r'|}.
$$

在第十章中，我们曾经指出，一维空间中的 Green 函数是处处连续的，而它的一阶导数不连续。上面讨论的三维空间中 Green 函数在点源处的行为显然和一维空间中 Green 函数的行为不同。原因是空间的维数不同，“点源”的奇异程度也不相同：一维空间中的点源实际上是三维空间中的面源。不难预料，二维空间中的 Green 函数也应该表现出不同的行为。

对于二维空间中 Poisson 方程第一边值问题，它的 Green 函数 $G(x,y;x',y')$ 是定解问题

$$
\left(
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
\right)
G(x,y;x',y')
=
-\frac1{\varepsilon_0}
\delta(x-x')
\delta(y-y'),
\qquad
(x,y),(x',y')\in S,
\tag{19.12a}
$$

$$
\left.
G(x,y;x',y')
\right|_C
=
0
\tag{19.12b}
$$

的解，其中 $C$ 是平面区域 $S$ 的边界。模仿上面三维情形的讨论，可以得出，这时的 Green 函数 $G(x,y;x',y')$ 应当是

$$
G(x,y;x',y')
=
-\frac1{2\pi\varepsilon_0}
\ln
\sqrt{
(x-x')^2+(y-y')^2
}
+
g(x,y;x',y'),
\tag{19.13}
$$

其中第一项是单位点电荷在无界空间中的电势（还可以加上一个常数，取决于电势零点的选取），在“点源”（实际上是三维空间中的线源）$\delta(x-x')\delta(y-y')$ 处是对数发散的，第二项 $g(x,y;x',y')$ 是边界上的感生电荷产生的电势，在 $S$ 内处处连续。

再来讨论 Green 函数的对称性。为此考察一下解式 (19.5')，容易发现这个结果在物理上有令人费解之处：在右端的体积分中，$G(\boldsymbol r';\boldsymbol r)$ 代表 $\boldsymbol r$ 处的单位点电荷在 $\boldsymbol r'$ 处的电势，它乘上在观测点 $\boldsymbol r'$ 处的电荷 $\rho(\boldsymbol r')d\boldsymbol r'$，并对观测点积分，却给出 $\boldsymbol r$ 处的电势。对这个问题的回答要涉及 Green 函数的对称性。因为，如果像无界空间的 Green 函数那样，关系式

$$
G(\boldsymbol r';\boldsymbol r)
=
G(\boldsymbol r;\boldsymbol r')
\tag{19.14}
$$

成立的话，那么，(19.5') 式就能改写成

$$
u(\boldsymbol r)
=
\iiint_V
G(\boldsymbol r;\boldsymbol r')
\rho(\boldsymbol r')
\,d\boldsymbol r'
-
\varepsilon_0
\iint_\Sigma
f(\Sigma')
\left.
\nabla'G(\boldsymbol r;\boldsymbol r')
\right|_{\Sigma'}
\cdot d\boldsymbol\Sigma'.
\tag{19.5'''}
$$

体积分就具有明确的物理意义，第二项的面积分就是来自边界面上的感生面电荷的贡献。

下面就来证明 (19.14) 式。和第十章中的做法一样，再引进 $G(\boldsymbol r;\boldsymbol r'')$，它满足定解问题

$$
\nabla^2G(\boldsymbol r;\boldsymbol r'')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r''),
\qquad
\boldsymbol r,\boldsymbol r''\in V,
\tag{19.15a}
$$

$$
\left.
G(\boldsymbol r;\boldsymbol r'')
\right|_\Sigma
=
0.
\tag{19.15b}
$$

将方程 (19.3a) 和 (19.15a) 分别乘以 $G(\boldsymbol r;\boldsymbol r'')$ 和 $G(\boldsymbol r;\boldsymbol r')$，相减，然后在区域 $V$ 内积分，即得

$$
\iiint_V
\left[
G(\boldsymbol r;\boldsymbol r'')
\nabla^2G(\boldsymbol r;\boldsymbol r')
-
G(\boldsymbol r;\boldsymbol r')
\nabla^2G(\boldsymbol r;\boldsymbol r'')
\right]
d\boldsymbol r
=
-\frac1{\varepsilon_0}
\left[
G(\boldsymbol r';\boldsymbol r'')
-
G(\boldsymbol r'';\boldsymbol r')
\right].
$$

根据 Green 公式，将上式左端的体积分化为面积分，就有

$$
G(\boldsymbol r';\boldsymbol r'')
-
G(\boldsymbol r'';\boldsymbol r')
=
-\varepsilon_0
\iint_\Sigma
\left[
G(\boldsymbol r;\boldsymbol r'')
\nabla G(\boldsymbol r;\boldsymbol r')
-
G(\boldsymbol r;\boldsymbol r')
\nabla G(\boldsymbol r;\boldsymbol r'')
\right]\cdot d\boldsymbol\Sigma.
$$

代入边界条件 (19.3b) 和 (19.15b)，立即得出右端的面积分为 0。这样就证明了

$$
G(\boldsymbol r';\boldsymbol r'')
=
G(\boldsymbol r'';\boldsymbol r').
$$

将 $\boldsymbol r''$ 改写为 $\boldsymbol r$，这就是 (19.14) 式。

如果是第三类边界条件，上面的结论仍然正确。请读者自己补上证明。

对于其他类型的稳定问题，它们的 Green 函数是否仍然有对称关系 (19.14) 需要具体讨论。至少从形式上说，这涉及 $G(\boldsymbol r;\boldsymbol r')$ 和 $G(\boldsymbol r';\boldsymbol r)$ 是否都是同一方程的解。

## §19.3 三维无界空间 Helmholtz 方程的 Green 函数

§19.1 已经给出了三维无界空间 Poisson 方程的 Green 函数的具体形式，即 (19.10) 式。现在再来求三维无界空间中 Helmholtz 方程的 Green 函数，即在三维无界空间中求解

$$
\nabla^2G(\boldsymbol r;\boldsymbol r')
+
k^2G(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\qquad
\boldsymbol r,\boldsymbol r'\in V.
\tag{19.16}
$$


无穷远处的边界条件暂缺，后面再讨论。

方程 (19.16) 是一个非齐次方程，因此，可以按照求解非齐次方程的标准做法，或是先求出方程的一个特解，而将方程齐次化，或是将 $G(\boldsymbol r;\boldsymbol r')$ 按相应齐次问题的本征函数展开。这两种做法，特别是第二种做法，没有原则困难，这里不拟做具体的介绍。现在要强调的是，这是一个特殊的非齐次方程：只在 $\boldsymbol r=\boldsymbol r'$ 点，非齐次项才不为 0。同时，由于是在无界空间，可以适当地安置坐标架，使问题得到充分的简化。

为此，首先作坐标平移

$$
\xi=x-x',
\qquad
\eta=y-y',
\qquad
\zeta=z-z',
\tag{19.17}
$$

即将点电荷所在点取为新坐标系的原点。令 $G(\boldsymbol r;\boldsymbol r')=g(\xi,\eta,\zeta)$，则 $g(\xi,\eta,\zeta)$ 满足方程

$$
\nabla_{\xi,\eta,\zeta}^2
g(\xi,\eta,\zeta)
+
k^2g(\xi,\eta,\zeta)
=
-\frac1{\varepsilon_0}
\delta(\xi)\delta(\eta)\delta(\zeta),
\tag{19.18}
$$

其中

$$
\nabla_{\xi,\eta,\zeta}^2
\equiv
\frac{\partial^2}{\partial\xi^2}
+
\frac{\partial^2}{\partial\eta^2}
+
\frac{\partial^2}{\partial\zeta^2}
$$

是以直角坐标 $\xi,\eta,\zeta$ 为自变量的 Laplace 算符。容易看出，方程 (19.18) 是旋转不变的，存在解 $g(\xi,\eta,\zeta)$ 只是 $R=\sqrt{\xi^2+\eta^2+\zeta^2}$ 的函数，$g(\xi,\eta,\zeta)=f(R)$。因此，如果将直角坐标系 $(\xi,\eta,\zeta)$ 转换为球坐标系，则方程 (19.18) 将变为 $R\ne0$ 点处的齐次方程

$$
\frac1{R^2}
\frac d{dR}
\left[
R^2\frac{df(R)}{dR}
\right]
+
k^2f(R)
=
0
\tag{19.19}
$$

（原因是在 $R=0$ 点只存在单侧导数）以及 $R=0$ 点处的边界条件（在 $R=0$ 点处有一单位点电荷）。方程 (19.19) 可化为零阶球 Bessel 方程，它的通解是①

$$
f(R)
=
A(k)\frac{e^{ikR}}R
+
B(k)\frac{e^{-ikR}}R.
\tag{19.20}
$$

现在就根据 $R=0$ 和无穷远处的边界条件定出常数 $A(k)$ 和 $B(k)$。先讨论无穷远处。考虑到 Helmholtz 方程的实际背景，比如说，它是由波动方程经过分离变量（分离去时间部分）得到的。这时，作为一个例子，假设要求得到的解在无穷远处为发散波。取时间因子为 $e^{-i\omega t}$，则应保留 (19.20) 中的第一项（即 $A(k)\ne0$），而弃去第二项（即令 $B(k)=0$）。常数 $A(k)$ 应该由 $R=0$ 处的边界条件决定，更准确地说，由 $R=0$ 处点源的强度决定。注意，这时并不能直接将 (19.20) 式代入方程 (19.18) 而定出 $A(k)$，原因是 $f(R)$ 或 $g(\xi,\eta,\zeta)$ 在 $R=0$ 处的导数并不存在。另一方面，我们已经约定，凡是涉及 $\delta$ 函数的等式都应该从积分意义下去理解。于是，很自然地，应当将方程 (19.18) 在 $R=0$ 附近的小体积内积分，

$$
\iiint
\nabla_{\xi,\eta,\zeta}^2f(R)
\,d\xi\,d\eta\,d\zeta
+
k^2
\iiint
f(R)
\,d\xi\,d\eta\,d\zeta
=
-\frac1{\varepsilon_0}.
\tag{19.21}
$$

① 这两个特解就是第一种和第二种零阶球 Hankel 函数。或者作变换 $f(R)=w(R)/R$，则方程 (19.19) 化为 $w''(R)+k^2w(R)=0$，也容易写出这个结果。

将左端第一项的体积分化为面积分

$$
\iiint
\nabla_{\xi,\eta,\zeta}^2f(R)
\,d\xi\,d\eta\,d\zeta
=
\iint
[\nabla_{\xi,\eta,\zeta}f(R)]
\cdot d\boldsymbol\Sigma,
$$

从而回避掉在 $R=0$ 点的求导问题。取小体积为以 $R=0$ 点为球心、$\rho$ 为半径的球体，则

$$
\begin{aligned}
\iiint
\nabla_{\xi,\eta,\zeta}^2f(R)
\,d\xi\,d\eta\,d\zeta
&=
\left.
\iint
[\nabla_{\xi,\eta,\zeta}f(R)]
\cdot d\boldsymbol\Sigma
\right|_{R=\rho}\\
&=
\iint
\frac{df(\rho)}{d\rho}
\rho^2\sin\theta\,d\theta\,d\phi\\
&=
-4\pi A(k)(1-ik\rho)e^{ik\rho}.
\end{aligned}
$$

第二项的体积分可以直接算出：

$$
\begin{aligned}
\iiint
f(R)\,d\xi\,d\eta\,d\zeta
&=
4\pi A(k)
\int_0^\rho
e^{ikR}R\,dR\\
&=
\frac{4\pi A(k)}{k^2}
\left[
(e^{ik\rho}-1)
-
ik\rho e^{ik\rho}
\right].
\end{aligned}
$$

将这些结果代回到 (19.21) 式，就有

$$
-4\pi A(k)
=
-\frac1{\varepsilon_0},
$$

所以 $A(k)=1/(4\pi\varepsilon_0)$，与 $k$ 无关。这样就求出了三维无界空间 Helmholtz 方程的 Green 函数

$$
g(\xi,\eta,\zeta)
=
f(R)
=
\frac1{4\pi\varepsilon_0}
\frac{e^{ikR}}R
\tag{19.22}
$$

或

$$
G(\boldsymbol r;\boldsymbol r')
=
\frac1{4\pi\varepsilon_0}
\frac{
e^{ik|\boldsymbol r-\boldsymbol r'|}
}{
|\boldsymbol r-\boldsymbol r'|
}.
\tag{19.23}
$$

当 $k=0$ 时，这个结果就回到 Poisson 方程的 Green 函数 (19.9)。

需要说明，这个结果是在无穷远处为发散的球面波，并取时间因子为 $e^{-i\omega t}$ 的条件下得到的。可以设想，如果要求无穷远处为会聚的球面波（且仍取时间因子为 $e^{-i\omega t}$），则 Green 函数是

$$
G(\boldsymbol r;\boldsymbol r')
=
\frac1{4\pi\varepsilon_0}
\frac{
e^{-ik|\boldsymbol r-\boldsymbol r'|}
}{
|\boldsymbol r-\boldsymbol r'|
}.
\tag{19.24}
$$

如果是其他形式的无穷远条件，当然还会得到其他形式的解。

下面介绍另一种解法。考虑到这是无界空间中的定解问题，可采用 Fourier 变换。令

$$
g(\boldsymbol K;\boldsymbol r')
=
\frac1{(2\pi)^{3/2}}
\iiint
G(\boldsymbol r;\boldsymbol r')
e^{-i\boldsymbol K\cdot\boldsymbol r}
\,d\boldsymbol r,
\tag{19.25}
$$

则方程 (19.16) 变为

$$
(-K^2+k^2)
g(\boldsymbol K;\boldsymbol r')
=
-\frac1{(2\pi)^{3/2}}
\frac1{\varepsilon_0}
e^{-i\boldsymbol K\cdot\boldsymbol r'}.
\tag{19.26}
$$

其中 $K^2=\boldsymbol K\cdot\boldsymbol K=|\boldsymbol K|^2$。注意代数方程 (19.26) 来源于偏微分方程 (19.16)，因此，从根本上说，应该在广义函数的意义下解此代数方程，从而得到①

$$
g(\boldsymbol K;\boldsymbol r')
=
\frac1{(2\pi)^{3/2}}
\frac1{\varepsilon_0}
\left[
\frac1{K^2-k^2}
+
C\delta(K^2-k^2)
\right]
e^{-i\boldsymbol K\cdot\boldsymbol r'}.
\tag{19.27}
$$

再求反演，就有

$$
\begin{aligned}
G(\boldsymbol r;\boldsymbol r')
&=
\frac1{(2\pi)^3}
\frac1{\varepsilon_0}
\iiint
\frac{
e^{i\boldsymbol K\cdot(\boldsymbol r-\boldsymbol r')}
}{
K^2-k^2
}
\,d\boldsymbol K\\
&\quad+
\frac1{(2\pi)^3}
\frac C{\varepsilon_0}
\iiint
\delta(K^2-k^2)
e^{i\boldsymbol K\cdot(\boldsymbol r-\boldsymbol r')}
\,d\boldsymbol K.
\end{aligned}
\tag{19.28}
$$

可以看出，$G(\boldsymbol r;\boldsymbol r')$ 只是 $\boldsymbol r-\boldsymbol r'$ 的函数。不妨令 $\boldsymbol R=\boldsymbol r-\boldsymbol r'$，然后改用 $\boldsymbol K$ 空间中的球坐标计算上面的三重积分。对于 (19.28) 式右端的第一项，可以得到

$$
\begin{aligned}
\frac1{(2\pi)^3}
\frac1{\varepsilon_0}
\iiint
\frac{e^{i\boldsymbol K\cdot\boldsymbol R}}
{K^2-k^2}
\,d\boldsymbol K
&=
\frac1{(2\pi)^3}
\frac1{\varepsilon_0}
\iiint
\frac{e^{iKR\cos\theta}}
{K^2-k^2}
K^2\sin\theta
\,dK\,d\theta\,d\phi\\
&=
\frac1{(2\pi)^2}
\frac1{\varepsilon_0}
\int_0^\infty
\frac{K^2}{K^2-k^2}
\,dK
\int_0^\pi
e^{iKR\cos\theta}
\sin\theta\,d\theta\\
&=
\frac1{(2\pi)^2}
\frac1{iR\varepsilon_0}
\int_0^\infty
\frac K{K^2-k^2}
(e^{iKR}-e^{-iKR})
\,dK\\
&=
\frac1{(2\pi)^2}
\frac1{iR\varepsilon_0}
\int_{-\infty}^{\infty}
\frac K{K^2-k^2}
e^{iKR}
\,dK.
\end{aligned}
$$

应用留数定理容易计算出这个定积分。考虑到被积函数在实轴上有两个奇点（均为一阶极点），$K=\pm k$，故可以采用图 19.1 中的围道。这样，当小圆弧的半径趋于 0 时，沿奇点 $K=\pm k$ 处半圆弧积分的极限值为

$$
-\pi i
\lim_{K\to\pm k}
(K\mp k)
\frac K{K^2-k^2}
e^{iKR}
=
-\frac{\pi i}{2}
e^{\pm ikR}.
$$

![Figure 19.1](...)

图 19.1

于是就可以计算得 (19.28) 式右端的第一项

$$
\begin{aligned}
\frac1{(2\pi)^3}
\frac1{\varepsilon_0}
\iiint
\frac{e^{i\boldsymbol K\cdot\boldsymbol R}}
{K^2-k^2}
\,d\boldsymbol K
&=
\frac1{(2\pi)^2}
\frac1{iR\varepsilon_0}
\frac{\pi i}{2}
(e^{ikR}+e^{-ikR})\\
&=
\frac1{4\pi\varepsilon_0}
\frac{\cos(kR)}R\\
&=
\frac1{4\pi\varepsilon_0}
\frac{
\cos(k|\boldsymbol r-\boldsymbol r'|)
}{
|\boldsymbol r-\boldsymbol r'|
}.
\end{aligned}
\tag{19.29}
$$

类似地，对于 (19.28) 式右端的第二项，有

$$
\begin{aligned}
\frac1{(2\pi)^3}
\frac C{\varepsilon_0}
\iiint
e^{i\boldsymbol K\cdot\boldsymbol R}
\delta(K^2-k^2)
\,d\boldsymbol K
&=
\frac1{(2\pi)^3}
\frac C{\varepsilon_0}
\int_0^\infty
\delta(K^2-k^2)
K^2\,dK
\iint
e^{iKR\cos\theta}
\sin\theta\,d\theta\,d\phi\\
&=
\frac1{2\pi^2}
\frac C{\varepsilon_0R}
\int_0^\infty
\delta(K^2-k^2)
\sin(KR)K\,dK.
\end{aligned}
$$

注意在 $0\le K,k<\infty$ 的条件下，

$$
\delta(K^2-k^2)
=
\frac1{2K}
\delta(K-k),
$$

所以

$$
\begin{aligned}
\frac1{(2\pi)^3}
\frac C{\varepsilon_0}
\iiint
e^{i\boldsymbol K\cdot\boldsymbol R}
\delta(K^2-k^2)
\,d\boldsymbol K
&=
\frac1{(2\pi)^2}
\frac C{\varepsilon_0}
\frac{\sin(kR)}R\\
&=
\frac1{(2\pi)^2}
\frac C{\varepsilon_0}
\frac{
\sin(k|\boldsymbol r-\boldsymbol r'|)
}{
|\boldsymbol r-\boldsymbol r'|
}.
\end{aligned}
\tag{19.30}
$$

把两项合并起来，就得到了三维无界空间中 Helmholtz 方程的 Green 函数

$$
G(\boldsymbol r;\boldsymbol r')
=
\frac1{4\pi\varepsilon_0}
\frac1{|\boldsymbol r-\boldsymbol r'|}
\left[
\cos(k|\boldsymbol r-\boldsymbol r'|)
+
\frac C\pi
\sin(k|\boldsymbol r-\boldsymbol r'|)
\right],
\tag{19.31}
$$

常数 $C$ 由无穷远条件决定。对照 (19.23) 式或 (19.24) 式，可以看出：如果取 $C=\pi i$，则给出发散的球面波；如果取 $C=-\pi i$，则给出会聚的球面波。

① 如果限定在连续函数，则齐次代数方程 $xf(x)=0$ 只有零解 $f(x)=0$；但在广义函数的范围内，由于 $x\delta(x)=0$（见第十章，(10.11) 式），所以此方程有“非零解” $f(x)=C\delta(x)$，$C$ 为任意常数。由此可见，(19.26) 式所对应的齐次方程 $(-K^2+k^2)g(\boldsymbol K;\boldsymbol r')=0$ 应当有解 $g(\boldsymbol K;\boldsymbol r')=C\delta(-K^2+k^2)=C\delta(K^2-k^2)$。

## §19.4 圆内 Poisson 方程第一边值问题的 Green 函数

圆内 Poisson 方程第一边值问题 Green 函数的定义是

$$
\nabla_2^2G(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\qquad
|\boldsymbol r|<a,\quad
|\boldsymbol r'|<a,
\tag{19.32a}
$$

$$
\left.
G(\boldsymbol r;\boldsymbol r')
\right|_{r=a}
=
0,
\tag{19.32b}
$$

其中

$$
r^2=x^2+y^2,
\qquad
\nabla_2^2
=
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}.
$$

首先介绍标准的解法。考虑到 (19.32a) 是非齐次方程，所以应该将 Green 函数按相应齐次问题的本征函数展开，为此，采用平面极坐标系，坐标原点放在圆心，

$$
G(\boldsymbol r;\boldsymbol r')
=
R_0(r)
+
\sum_{m=1}^{\infty}
\left[
R_{m1}(r)\cos m\phi
+
R_{m2}(r)\sin m\phi
\right].
\tag{19.33}
$$

同样，将 $\delta$ 函数也按该组本征函数展开：

$$
\begin{aligned}
\delta(\boldsymbol r-\boldsymbol r')
&=
\delta(x-x')\delta(y-y')\\
&=
\frac1{r'}
\delta(r-r')
\delta(\phi-\phi')\\
&=
\frac1{2\pi r'}
\delta(r-r')
\left[
1+
2\sum_{m=1}^{\infty}
(
\cos m\phi\cos m\phi'
+
\sin m\phi\sin m\phi'
)
\right].
\end{aligned}
\tag{19.34}
$$

于是，定解问题 (19.32) 就转化为求解 $R_0(r)$、$R_{m1}(r)$ 和 $R_{m2}(r)$。决定 $R_0(r)$ 的定解问题是

$$
\frac1r
\frac d{dr}
\left[
r\frac{dR_0(r)}{dr}
\right]
=
-\frac1{2\pi\varepsilon_0}
\frac1{r'}
\delta(r-r'),
\tag{19.35a}
$$

$$
R_0(0)\text{ 有界},
\qquad
R_0(a)=0.
\tag{19.35b}
$$

当 $r\ne r'$ 时，方程 (19.35a) 是齐次的，所以，在考虑到边界条件 (19.35b) 后，就有

$$
R_0(r)
=
\begin{cases}
A_0,
&
r<r',
\\[6pt]
B_0\ln\dfrac ra,
&
r>r'.
\end{cases}
$$

再根据 $R_0(r)$ 在 $r=r'$ 点的连续性，即

$$
R_0(r)\text{ 在 }r=r'\text{ 点连续},
\qquad
\left.
\frac{dR_0(r)}{dr}
\right|_{r'-0}^{r'+0}
=
-\frac1{2\pi\varepsilon_0}
\frac1{r'},
$$

定出 $A_0$ 和 $B_0$，于是就求得

$$
R_0(r)
=
\begin{cases}
-\dfrac1{2\pi\varepsilon_0}
\ln\dfrac{r'}a,
&
r<r',
\\[10pt]
-\dfrac1{2\pi\varepsilon_0}
\ln\dfrac ra,
&
r>r'.
\end{cases}
\tag{19.36}
$$

完全模仿上面的做法，可以由

$$
\left[
\frac1r
\frac d{dr}
\left(
r\frac d{dr}
\right)
-
\frac{m^2}{r^2}
\right]
R_{m1}(r)
=
-\frac1{\pi\varepsilon_0r'}
\delta(r-r')
\cos m\phi',
\tag{19.37a}
$$

$$
R_{m1}(0)\text{ 有界},
\qquad
R_{m1}(a)=0
\tag{19.37b}
$$

求出 $R_{m1}(r)$，

$$
R_{m1}(r)
=
\begin{cases}
-\dfrac1{2\pi\varepsilon_0}
\dfrac1m
\left[
\left(
\dfrac{rr'}{a^2}
\right)^m
-
\left(
\dfrac r{r'}
\right)^m
\right]
\cos m\phi',
&
r<r',
\\[12pt]
-\dfrac1{2\pi\varepsilon_0}
\dfrac1m
\left[
\left(
\dfrac{rr'}{a^2}
\right)^m
-
\left(
\dfrac{r'}r
\right)^m
\right]
\cos m\phi',
&
r>r'.
\end{cases}
\tag{19.38}
$$

同样，由

$$
\left[
\frac1r
\frac d{dr}
\left(
r\frac d{dr}
\right)
-
\frac{m^2}{r^2}
\right]
R_{m2}(r)
=
-\frac1{\pi\varepsilon_0r'}
\delta(r-r')
\sin m\phi',
\tag{19.39a}
$$

$$
R_{m2}(0)\text{ 有界},
\qquad
R_{m2}(a)=0
\tag{19.39b}
$$

求出 $R_{m2}(r)$，

$$
R_{m2}(r)
=
\begin{cases}
-\dfrac1{2\pi\varepsilon_0}
\dfrac1m
\left[
\left(
\dfrac{rr'}{a^2}
\right)^m
-
\left(
\dfrac r{r'}
\right)^m
\right]
\sin m\phi',
&
r<r',
\\[12pt]
-\dfrac1{2\pi\varepsilon_0}
\dfrac1m
\left[
\left(
\dfrac{rr'}{a^2}
\right)^m
-
\left(
\dfrac{r'}r
\right)^m
\right]
\sin m\phi',
&
r>r'.
\end{cases}
\tag{19.40}
$$

这样，就求得了圆内 Poisson 方程第一边值问题的 Green 函数：

$$
G(\boldsymbol r;\boldsymbol r')
=
\begin{cases}
-\dfrac1{2\pi\varepsilon_0}
\left\{
\ln\dfrac{r'}a
+
\displaystyle\sum_{m=1}^{\infty}
\frac1m
\left[
\left(
\dfrac{rr'}{a^2}
\right)^m
-
\left(
\dfrac r{r'}
\right)^m
\right]
\cos m(\phi-\phi')
\right\},
&
r<r',
\\[18pt]
-\dfrac1{2\pi\varepsilon_0}
\left\{
\ln\dfrac ra
+
\displaystyle\sum_{m=1}^{\infty}
\frac1m
\left[
\left(
\dfrac{rr'}{a^2}
\right)^m
-
\left(
\dfrac{r'}r
\right)^m
\right]
\cos m(\phi-\phi')
\right\},
&
r>r'.
\end{cases}
\tag{19.41}
$$

上面这种方法，得到的解式会是无穷级数。当然，不排除在某些特殊情形下可以将级数求和。例如，现在得到的解式 (19.41) 就是如此。读者不妨尝试求出它们的和函数。

下面再介绍一种方法，如果能够成功的话，它将直接给出解的有限形式。

大家知道，一旦在接地圆中放上点电荷后，在圆周上必然出现感生电荷。圆内任意一点的电势，就是点电荷的电势和感生电荷的电势的叠加。前者在点电荷所在点是对数发散的，而后者在圆内是处处连续的。如果我们能够方便地求出感生电荷在圆内所产生的电势，当然也就求出了整个圆内 Poisson 方程第一边值问题的 Green 函数。现在要介绍的这种方法（称为**电像法**），其基本思想是力图将边界上的感生电荷用一个等价的点电荷代替。换句话说，就把接地圆内的点电荷的问题等价地转化为无界空间中的两个点电荷（一个是真实的点电荷，另一个是等价的“虚”电荷）的问题。这个“虚”电荷的等价性，就表现在它和圆内的真实的点电荷一起，在圆内能给出和原来问题同样的解。而由于边值问题解的存在唯一性，我们知道，只要这两个点电荷也能产生圆周 $r=a$ 接地（电势为 0）的效果，只要圆内的电荷分布不变，就能保证这样得到的解和原来问题的解在圆内一定是一致的。可以预见，这个等价电荷如果存在的话，它一定位于圆外，否则圆内的电荷分布就和原来的问题不同，就不能保证等价性。或者换一种说法，由于感生电荷的电势在圆内是处处连续的，在圆内的任何等价电荷都不可能产生同样的效果。应用电像法成败的关键，就在于能否成功地求出等价电荷的电量和它所在的空间位置。

根据对称性的考虑，我们还可以进一步断定，如果这个等价电荷存在的话，它还一定位于真实电荷所处的半径的延长线上。如图 19.2 所示，设这个等价电荷的位置为 $\boldsymbol r_1=(x_1,y_1)$，电量为 $e$，于是，它和真实点电荷一起，在圆内的电势就是

![Figure 19.2](...)

图 19.2　电像法

$$
G(\boldsymbol r;\boldsymbol r')
=
-\frac1{2\pi\varepsilon_0}
\left[
\ln|\boldsymbol r-\boldsymbol r'|
+
e\ln|\boldsymbol r-\boldsymbol r_1|
+
C
\right],
\tag{19.42}
$$

其中常数 $C$ 与电势零点的选择有关。现在的问题就是要从要求圆周 $r=a$ 上的电势为 0，

$$
-\frac1{2\pi\varepsilon_0}
\left(
\ln|\boldsymbol r-\boldsymbol r'|
+
e\ln|\boldsymbol r-\boldsymbol r_1|
+
C
\right)_{r=a}
=
0,
\tag{19.43}
$$

求出 $\boldsymbol r_1$、$e$ 和 $C$。注意等式 (19.43) 应该对圆周上的一切点均成立。如果采用平面极坐标来写出 (19.42) 式及 (19.43) 式中各项的具体形式，即令

$$
x=r\cos\phi,
\qquad
x'=r'\cos\phi',
\qquad
x_1=r_1\cos\phi',
$$

$$
y=r\sin\phi,
\qquad
y'=r'\sin\phi',
\qquad
y_1=r_1\sin\phi',
$$

则 (19.43) 式化为

$$
\ln
[
a^2+r'^2-2ar'\cos(\phi-\phi')
]
+
e\ln
[
a^2+r_1^2-2ar_1\cos(\phi-\phi')
]
+
2C
=
0,
\tag{19.44}
$$

它应该对一切 $\phi$ 均成立。注意当 $|t|<1$ 时有展开式

$$
\begin{aligned}
\ln(1+t^2-2t\cos\phi)
&=
\ln(1-te^{i\phi})
+
\ln(1-te^{-i\phi})\\
&=
-2
\sum_{m=1}^{\infty}
\frac1m
t^m\cos m\phi,
\end{aligned}
\tag{19.45}
$$

于是就可以将 (19.44) 式化为

$$
\begin{aligned}
&2\ln a
+
\ln
\left[
1+
\left(
\frac{r'}a
\right)^2
-
2\frac{r'}a
\cos(\phi-\phi')
\right]\\
&\quad+
2e\ln r_1
+
e\ln
\left[
1+
\left(
\frac a{r_1}
\right)^2
-
2\frac a{r_1}
\cos(\phi-\phi')
\right]
+
2C\\
&=
2\ln a
+
2e\ln r_1
+
2C
-
2\sum_{m=1}^{\infty}
\frac1m
\left[
\left(
\frac{r'}a
\right)^m
+
e
\left(
\frac a{r_1}
\right)^m
\right]
\cos m(\phi-\phi')
=
0.
\end{aligned}
$$

这样就得到

$$
\ln a
+
e\ln r_1
+
C
=
0
\tag{19.46}
$$

和

$$
\left(
\frac{r'}a
\right)^m
+
e
\left(
\frac a{r_1}
\right)^m
=
0,
\qquad
m=1,2,3,\cdots.
\tag{19.47}
$$

将 (19.47) 式化成

$$
e
=
-
\left(
\frac{r_1r'}{a^2}
\right)^m,
\qquad
m=1,2,3,\cdots,
$$

就可以得到 $e=-1$ 和

$$
r_1
=
\frac{a^2}{r'},
\qquad
\text{或}
\qquad
\boldsymbol r_1
=
\left(
\frac a{r'}
\right)^2
\boldsymbol r'.
$$

这样，我们的确求出了这个等价电荷，它位于真实电荷所在半径的延长线上，并且满足

$$
r'r_1
=
a^2.
$$

凡处于同一半径方向且满足这个关系的两个点，称为关于圆 $r=a$ 的反演点对，即这一对点互为关于圆 $r=a$ 的反演点。上面的结果说明，等价电荷和真实电荷构成对于圆 $r=a$ 的反演点对，它们的电量相等，而极性相反。将 $e$ 和 $r_1$ 的结果代入 (19.46) 式，又可以求得

$$
C
=
-\ln a+\ln r_1
=
\ln\frac a{r'}.
$$

将 $e$、$\boldsymbol r_1$ 和 $C$ 的结果代回 (19.42) 式，就求得圆内 Poisson 方程第一边值问题的 Green 函数

$$
G(\boldsymbol r;\boldsymbol r')
=
-\frac1{2\pi\varepsilon_0}
\left[
\ln|\boldsymbol r-\boldsymbol r'|
-
\ln
\left|
\boldsymbol r
-
\left(
\frac a{r'}
\right)^2
\boldsymbol r'
\right|
+
\ln\frac a{r'}
\right],
\tag{19.48}
$$

或者在极坐标系中的表达式

$$
\begin{aligned}
G(\boldsymbol r;\boldsymbol r')
=
-\frac1{4\pi\varepsilon_0}
\Bigg\{
&
\ln
\left[
r^2+r'^2-2rr'\cos(\phi-\phi')
\right]\\
&-
\ln
\left[
r^2
+
\left(
\frac{a^2}{r'}
\right)^2
-
2r\frac{a^2}{r'}
\cos(\phi-\phi')
\right]\\
&+
2\ln\frac a{r'}
\Bigg\}.
\end{aligned}
\tag{19.49}
$$

利用展开式 (19.45)，求出 (19.41) 中的级数和，其结果也正是 (19.49) 式。

以上介绍了电像法求解 Green 函数的基本思路。这种方法只适用于某些非常特殊的几何形状（例如球形，半无界空间等等）。即使放宽到同时使用多个像电荷来等效地代替边界上的感生电荷，对空间的几何形状仍然有相当严格的限制。所以说，电像法的优点是可以给出有限形式的解，缺点是适用范围有限，而且多限于二维稳定问题。

求出圆内 Poisson 方程第一边值问题的 Green 函数后，就可以导出一般的定解问题

$$
\nabla_2^2u(\boldsymbol r)
=
-\frac1{\varepsilon_0}
\rho(\boldsymbol r),
\qquad
|\boldsymbol r|<a,
\tag{19.50a}
$$

$$
\left.
u(\boldsymbol r)
\right|_{r=a}
=
f(\phi)
\tag{19.50b}
$$

的解。为此，将方程 (19.50) 的自变量改写成 $\boldsymbol r'$，

$$
\nabla_2'^2u(\boldsymbol r')
=
-\frac1{\varepsilon_0}
\rho(\boldsymbol r'),
\qquad
|\boldsymbol r'|<a,
\tag{19.51a}
$$

$$
\left.
u(\boldsymbol r')
\right|_{r'=a}
=
f(\phi').
\tag{19.51b}
$$

另一方面，还可以写出 $G(\boldsymbol r';\boldsymbol r)$ 所应该满足的定解问题，

$$
\nabla_2'^2G(\boldsymbol r';\boldsymbol r)
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\qquad
|\boldsymbol r|<a,\quad
|\boldsymbol r'|<a,
$$

$$
\left.
G(\boldsymbol r';\boldsymbol r)
\right|_{r'=a}
=
0.
$$

再利用 Green 函数的对称性（它可以看成是 (19.14) 的特殊情形，也能从上面求出的 $G(\boldsymbol r;\boldsymbol r')$ 的具体表达式直接看出），

$$
G(\boldsymbol r;\boldsymbol r')
=
G(\boldsymbol r';\boldsymbol r),
$$

进一步将其改写成

$$
\nabla_2'^2G(\boldsymbol r;\boldsymbol r')
=
-\frac1{\varepsilon_0}
\delta(\boldsymbol r-\boldsymbol r'),
\qquad
|\boldsymbol r|<a,\quad
|\boldsymbol r'|<a,
\tag{19.52a}
$$

$$
\left.
G(\boldsymbol r;\boldsymbol r')
\right|_{r'=a}
=
0.
\tag{19.52b}
$$

将方程 (19.51a) 和 (19.52a) 分别乘以 $G(\boldsymbol r;\boldsymbol r')$ 和 $u(\boldsymbol r')$，然后相减，再在圆内积分，就得到

$$
\iint_{r'<a}
\rho(\boldsymbol r')
G(\boldsymbol r;\boldsymbol r')
\,d\boldsymbol r'
-
u(\boldsymbol r)
=
-\varepsilon_0
\iint_{r'<a}
\left[
G(\boldsymbol r;\boldsymbol r')
\nabla_2'^2u(\boldsymbol r')
-
u(\boldsymbol r')
\nabla_2'^2G(\boldsymbol r;\boldsymbol r')
\right]
d\boldsymbol r'.
$$

把上面的面积分化为沿圆周 $r=a$ 的线积分，并代入边界条件 (19.51b) 和 (19.52b)，就有

$$
\begin{aligned}
u(\boldsymbol r)
&=
\iint_{r'<a}
\rho(\boldsymbol r')
G(\boldsymbol r;\boldsymbol r')
\,d\boldsymbol r'\\
&\quad+
\varepsilon_0
\int_0^{2\pi}
\left[
G(\boldsymbol r;\boldsymbol r')
\frac{\partial u(\boldsymbol r')}{\partial r'}
-
u(\boldsymbol r')
\frac{\partial G(\boldsymbol r;\boldsymbol r')}{\partial r'}
\right]_{r'=a}
a\,d\phi'\\
&=
\iint_{r'<a}
\rho(\boldsymbol r')
G(\boldsymbol r;\boldsymbol r')
\,d\boldsymbol r'
-
\varepsilon_0
\int_0^{2\pi}
f(\phi')
\left.
\frac{\partial G(\boldsymbol r;\boldsymbol r')}{\partial r'}
\right|_{r'=a}
a\,d\phi'.
\end{aligned}
\tag{19.53}
$$

显然，右端的第一项表示圆内电荷分布的贡献；第二项则是来自圆周上的感生电荷产生的电势，感生电荷的分布当然与给定的边界条件（圆周上电势值的分布）有关。为了更清楚地看出圆周上的电荷分布，可以将第二项中的线积分再改写成

$$
\begin{aligned}
&
\int_0^{2\pi}
f(\phi')
\left.
\frac{\partial G(\boldsymbol r;\boldsymbol r')}{\partial r'}
\right|_{r'=a}
a\,d\phi'\\
&=
\int_0^{2\pi}
f(\phi')
\lim_{\Delta r\to0}
\frac1{\Delta r}
\left[
\left.G(\boldsymbol r;\boldsymbol r')\right|_{r'=a}
-
\left.G(\boldsymbol r;\boldsymbol r')\right|_{r'=a-\Delta r}
\right]
a\,d\phi'\\
&=
\iint_{r'<a}
f(\phi')
\lim_{\Delta r\to0}
\frac{G(\boldsymbol r;\boldsymbol r')}{\Delta r}
\left[
\delta(r'-a)
-
\delta(r'-a+\Delta r)
\right]
r'\,dr'\,d\phi'\\
&=
-
\iint_{r'<a}
f(\phi')
G(\boldsymbol r;\boldsymbol r')
\delta'(r'-a)
r'\,dr'\,d\phi'.
\end{aligned}
$$

这样，就可以把 (19.53) 式写成

$$
u(\boldsymbol r)
=
\iint_{r'<a}
\left[
\rho(\boldsymbol r')
+
\varepsilon_0
f(\phi')
\delta'(r'-a)
\right]
G(\boldsymbol r;\boldsymbol r')
\,d\boldsymbol r'.
\tag{19.54}
$$

这个结果表明，圆周上的感生电荷密度就是 $\varepsilon_0f(\phi')\delta'(r'-a)$。函数 $\delta'(r'-a)$ 的出现，说明在圆周上的感生电荷为偶极层。

这里还可以引申出一个重要结论。设想有另外一个定解问题

$$
\nabla_2^2u(\boldsymbol r)
=
-\frac1{\varepsilon_0}
\left[
\rho(\boldsymbol r)
+
\varepsilon_0
f(\phi)
\delta'(r-a)
\right],
\qquad
|\boldsymbol r|<a,
\tag{19.52a'}
$$

$$
\left.
u(\boldsymbol r)
\right|_{r=a}
=
0,
\tag{19.52b'}
$$

显然，它也将会有同样的解 (19.54)。这说明，在引进 $\delta$ 函数及其导数的前提下，非齐次边界条件的定解问题，也可以改写为齐次的边界条件，相应地在方程中增加一项特殊的非齐次项（含 $\delta$ 函数）。上面的分析就提供了这种非齐次项的写法。当然，非齐次边界条件可以转化为方程的特殊形式的非齐次项，丝毫不意味着可以混淆非齐次边界条件（描写边界面上的源的分布）与方程非齐次项（区域内部的源的分布）的区别。即使把非齐次边界条件改写成方程的非齐次项，它描写的仍然是存在于边界面的源。

现在再回到 (19.53) 式，代入 $G(\boldsymbol r;\boldsymbol r')$ 的表达式，就得到

$$
\begin{aligned}
u(\boldsymbol r)
&=
-\frac1{2\pi\varepsilon_0}
\iint_{r'<a}
\rho(\boldsymbol r')
\left[
\ln|\boldsymbol r-\boldsymbol r'|
-
\ln
\left|
\boldsymbol r
-
\left(
\frac a{r'}
\right)^2
\boldsymbol r'
\right|
+
\ln\frac a{r'}
\right]
d\boldsymbol r'\\
&\quad+
\frac{a^2-r^2}{2\pi}
\int_0^{2\pi}
\frac{
f(\phi')
}{
a^2+r^2-2ar\cos(\phi-\phi')
}
\,d\phi'.
\end{aligned}
\tag{19.55}
$$

当 $\rho(\boldsymbol r)=0$，就有圆内 Laplace 方程第一边值问题的 Poisson 公式

$$
u(\boldsymbol r)
=
\frac{a^2-r^2}{2\pi}
\int_0^{2\pi}
\frac{
f(\phi')
}{
a^2+r^2-2ar\cos(\phi-\phi')
}
\,d\phi'.
\tag{19.56}
$$

用复变函数方法也能得到这个结果，见 (3.39) 式。


## *§19.5 波动方程的 Green 函数

下面用 Green 函数方法研究与时间有关的定解问题。为了确定起见，以有界弦的波动问题为例。最一般的定解问题就是

$$
\frac{\partial^2u(x,t)}{\partial t^2}
-
a^2
\frac{\partial^2u(x,t)}{\partial x^2}
=
f(x,t),
\qquad
0<x<l,\quad t>0,
\tag{19.57a}
$$

$$
\left.u(x,t)\right|_{x=0}
=
\mu(t),
\qquad
\left.u(x,t)\right|_{x=l}
=
\nu(t),
\qquad
t>0,
\tag{19.57b}
$$

$$
\left.u(x,t)\right|_{t=0}
=
\phi(x),
\qquad
\left.
\frac{\partial u(x,t)}{\partial t}
\right|_{t=0}
=
\psi(x),
\qquad
0<x<l.
\tag{19.57c}
$$

可以预料，相应的 Green 函数 $G(x,t;x',t')$ 应该是瞬时（仅存在于某时刻）点（仅存在于空间某点）源问题

$$
\left(
\frac{\partial^2}{\partial t^2}
-
a^2
\frac{\partial^2}{\partial x^2}
\right)
G(x,t;x',t')
=
\delta(x-x')\delta(t-t'),
\qquad
0<x,x'<l,\quad t'>0
\tag{19.58a}
$$

在齐次定解条件

$$
\left.G(x,t;x',t')\right|_{x=0}
=
0,
\qquad
\left.G(x,t;x',t')\right|_{x=l}
=
0,
\qquad
t'>0,
\tag{19.58b}
$$

$$
\left.G(x,t;x',t')\right|_{t<t'}
=
0,
\qquad
\left.
\frac{\partial G(x,t;x',t')}{\partial t}
\right|_{t<t'}
=
0,
\qquad
0<x,x'<l
\tag{19.58c}
$$

下的解。这里初始条件的物理意义是很清楚的：在 $t=t'$ 时刻驱动力出现之前，弦一直保持静止。注意在上面的 (19.58a) 和 (19.58b) 式中去掉了 $t>0$ 的限制，换言之，我们已将 $G(x,t;x',t')$ 延拓到了 $-\infty<t<\infty$。

和稳定问题的 Green 函数问题一样，现在需要分别讨论下列三个问题：

**1. Green 函数 $G(x,t;x',t')$ 的对称性**

为此，再列出关于 Green 函数 $G(x,-t;x'',-t'')$ 的定解问题

$$
\left(
\frac{\partial^2}{\partial t^2}
-
a^2
\frac{\partial^2}{\partial x^2}
\right)
G(x,-t;x'',-t'')
=
\delta(x-x'')\delta(t-t''),
\qquad
0<x,x''<l,\quad t''>0,
\tag{19.59a}
$$

$$
\left.G(x,-t;x'',-t'')\right|_{x=0}
=
0,
\qquad
\left.G(x,-t;x'',-t'')\right|_{x=l}
=
0,
\qquad
t''>0,
\tag{19.59b}
$$

$$
\left.G(x,-t;x'',-t'')\right|_{-t<-t''}
=
0,
\qquad
\left.
\frac{\partial G(x,-t;x'',-t'')}{\partial t}
\right|_{-t<-t''}
=
0,
\qquad
0<x,x''<l.
\tag{19.59c}
$$

将方程 (19.58a) 和 (19.59a) 分别乘以 Green 函数 $G(x,-t;x'',-t'')$ 和 $G(x,t;x',t')$，然后相减，再在区间 $[0,l]$ 和 $(-\infty,\infty)$ 上对 $x$ 和 $t$ 积分，有

$$
\begin{aligned}
&G(x',-t';x'',-t'')
-
G(x'',t'';x',t')\\
&=
\int_0^l
dx
\int_{-\infty}^{\infty}
\left[
G(x,-t;x'',-t'')
\frac{\partial^2G(x,t;x',t')}{\partial t^2}
-
G(x,t;x',t')
\frac{\partial^2G(x,-t;x'',-t'')}{\partial t^2}
\right]dt\\
&\quad-
a^2
\int_{-\infty}^{\infty}
dt
\int_0^l
\left[
G(x,-t;x'',-t'')
\frac{\partial^2G(x,t;x',t')}{\partial x^2}
-
G(x,t;x',t')
\frac{\partial^2G(x,-t;x'',-t'')}{\partial x^2}
\right]dx\\
&=
\int_0^l
\left[
G(x,-t;x'',-t'')
\frac{\partial G(x,t;x',t')}{\partial t}
-
G(x,t;x',t')
\frac{\partial G(x,-t;x'',-t'')}{\partial t}
\right]_{-\infty}^{\infty}
dx\\
&\quad-
a^2
\int_{-\infty}^{\infty}
\left[
G(x,-t;x'',-t'')
\frac{\partial G(x,t;x',t')}{\partial x}
-
G(x,t;x',t')
\frac{\partial G(x,-t;x'',-t'')}{\partial x}
\right]_0^l
dt.
\end{aligned}
$$

代入边界条件 (19.58b)、(19.59b) 和初始条件 (19.58c)、(19.59c)，可以看出，右端的积分为 0，这样就导出了 Green 函数在空间上的对称性与时间上的倒易性：

$$
G(x'',t'';x',t')
=
G(x',-t';x'',-t''),
$$

或者将 $x''$ 和 $t''$ 改写成 $x$ 和 $t$，

$$
G(x,t;x',t')
=
G(x',-t';x,-t).
\tag{19.60}
$$

在此式中，将 $t$ 和 $t'$ 易位时出现的负号，正好保证了时间的先后次序不变，否则就会有悖于因果律的要求。

**2. 用 Green 函数及已知条件表示相关定解问题的解**

为了用 Green 函数及已知条件表示定解问题 (19.57) 的解，可先将该定解问题中的自变量改写为 $x'$ 和 $t'$，

$$
\frac{\partial^2u(x',t')}{\partial t'^2}
-
a^2
\frac{\partial^2u(x',t')}{\partial x'^2}
=
f(x',t'),
\qquad
0<x'<l,\quad t'>0,
\tag{19.61a}
$$

$$
\left.u(x',t')\right|_{x'=0}
=
\mu(t'),
\qquad
\left.u(x',t')\right|_{x'=l}
=
\nu(t'),
\qquad
t'>0,
\tag{19.61b}
$$

$$
\left.u(x',t')\right|_{t'=0}
=
\phi(x'),
\qquad
\left.
\frac{\partial u(x',t')}{\partial t'}
\right|_{t'=0}
=
\psi(x'),
\qquad
0<x'<l,
\tag{19.61c}
$$

再写出 $G(x',-t';x,-t)$ 满足的定解问题，并利用 Green 函数的对称性与倒易性关系 (19.60) 而改写为

$$
\left(
\frac{\partial^2}{\partial t'^2}
-
a^2
\frac{\partial^2}{\partial x'^2}
\right)
G(x,t;x',t')
=
\delta(x-x')\delta(t-t'),
\qquad
0<x,x'<l,\quad t,t'>0,
\tag{19.62a}
$$

$$
\left.G(x,t;x',t')\right|_{x'=0}
=
0,
\qquad
\left.G(x,t;x',t')\right|_{x'=l}
=
0,
\qquad
t,t'>0,
\tag{19.62b}
$$

$$
\left.G(x,t;x',t')\right|_{t'>t}
=
0,
\qquad
\left.
\frac{\partial G(x,t;x',t')}{\partial t'}
\right|_{t'>t}
=
0,
\qquad
0<x,x'<l.
\tag{19.62c}
$$

将方程 (19.61a) 和 (19.62a) 分别乘以 $G(x,t;x',t')$ 和 $u(x',t')$，然后相减，再积分，有

$$
\begin{aligned}
&\int_0^l
dx'
\int_0^\infty
G(x,t;x',t')
f(x',t')
\,dt'
-
u(x,t)\\
&=
\int_0^l
dx'
\int_0^\infty
\left[
G(x,t;x',t')
\frac{\partial^2u(x',t')}{\partial t'^2}
-
u(x',t')
\frac{\partial^2G(x,t;x',t')}{\partial t'^2}
\right]dt'\\
&\quad-
a^2
\int_0^\infty
dt'
\int_0^l
\left[
G(x,t;x',t')
\frac{\partial^2u(x',t')}{\partial x'^2}
-
u(x',t')
\frac{\partial^2G(x,t;x',t')}{\partial x'^2}
\right]dx'\\
&=
\int_0^l
\left[
G(x,t;x',t')
\frac{\partial u(x',t')}{\partial t'}
-
u(x',t')
\frac{\partial G(x,t;x',t')}{\partial t'}
\right]_0^\infty
dx'\\
&\quad-
a^2
\int_0^\infty
\left[
G(x,t;x',t')
\frac{\partial u(x',t')}{\partial x'}
-
u(x',t')
\frac{\partial G(x,t;x',t')}{\partial x'}
\right]_0^l
dt'.
\end{aligned}
$$

代入边界条件 (19.61b)、(19.62b) 和初始条件 (19.61c)、(19.62c)，可以将上面的结果化简为

$$
\begin{aligned}
u(x,t)
&=
\int_0^l
dx'
\int_0^t
G(x,t;x',t')
f(x',t')
\,dt'\\
&\quad+
\int_0^l
\left[
G(x,t;x',0)\psi(x')
-
\phi(x')
\left.
\frac{\partial G(x,t;x',t')}{\partial t'}
\right|_{t'=0}
\right]dx'\\
&\quad-
a^2
\int_0^t
\left[
\nu(t')
\left.
\frac{\partial G(x,t;x',t')}{\partial x'}
\right|_{x'=l}
-
\mu(t')
\left.
\frac{\partial G(x,t;x',t')}{\partial x'}
\right|_{x'=0}
\right]dt'.
\end{aligned}
\tag{19.63}
$$

**练习 19.1**　如果直接从 (19.57) 和 (19.60) 式出发，而不将这些方程的自变量换成 $x'$ 和 $t'$，是否能够用 Green 函数以及已知条件 $f(x,t)$、$\mu(t)$、$\nu(t)$ 和 $\phi(x)$、$\psi(x)$ 将定解问题 (19.57) 的解 $u(x,t)$ 表示出来？

**3. 求解 Green 函数**

现在讨论如何由定解问题 (19.58) 求出 Green 函数。因为这是非齐次方程的定解问题，第一种解法仍然是按相应齐次问题的本征函数展开，

$$
G(x,t;x',t')
=
\sum_{n=1}^{\infty}
T_n(t)
\sin\frac{n\pi}{l}x.
\tag{19.64}
$$

同时，将 $\delta$ 函数也按该组本征函数展开，

$$
\delta(x-x')
=
\frac2l
\sum_{n=1}^{\infty}
\sin\frac{n\pi}{l}x'
\sin\frac{n\pi}{l}x.
\tag{19.65}
$$

于是，$T_n(t)$ 就满足常微分方程的初值问题

$$
T_n''(t)
+
\left(
\frac{n\pi a}{l}
\right)^2
T_n(t)
=
\frac2l
\sin\frac{n\pi}{l}x'
\delta(t-t'),
\tag{19.66a}
$$

$$
T_n(t<t')
=
0,
\qquad
T_n'(t<t')
=
0.
\tag{19.66b}
$$

解之即得

$$
T_n(t)
=
\frac2{n\pi a}
\sin\frac{n\pi}{l}x'
\sin
\left[
\frac{n\pi}{l}a(t-t')
\right]
\eta(t-t').
\tag{19.67}
$$

所以，就求得 Green 函数

$$
G(x,t;x',t')
=
\frac2{\pi a}
\sum_{n=1}^{\infty}
\frac1n
\sin\frac{n\pi}{l}x'
\sin\frac{n\pi}{l}x
\sin
\left[
\frac{n\pi}{l}a(t-t')
\right]
\eta(t-t').
\tag{19.68}
$$

第二种方法是将定解问题 (19.58) 作 Laplace 变换。令

$$
g(x,p;x',t')
=
\int_0^\infty
G(x,t;x',t')
e^{-pt}
\,dt,
\tag{19.69}
$$

则 $g(x,p;x',t')$ 满足常微分方程的边值问题

$$
\left[
\frac{d^2}{dx^2}
-
\left(
\frac pa
\right)^2
\right]
g(x,p;x',t')
=
-\frac1{a^2}
e^{-pt'}
\delta(x-x'),
\tag{19.70a}
$$

$$
\left.
g(x,p;x',t')
\right|_{x=0}
=
0,
\qquad
\left.
g(x,p;x',t')
\right|_{x=l}
=
0.
\tag{19.70b}
$$

由此也可求得

$$
g(x,p;x',t')
=
\begin{cases}
\displaystyle
\frac{
\sinh\frac pa(l-x')
\sinh\frac pa x
}{
pa\sinh\frac pa l
}
e^{-pt'},
&
0\le x<x',
\\[16pt]
\displaystyle
\frac{
\sinh\frac pa x'
\sinh\frac pa(l-x)
}{
pa\sinh\frac pa l
}
e^{-pt'},
&
x'<x\le l.
\end{cases}
\tag{19.71}
$$

最后，求反演，

$$
G(x,t;x',t')
=
\frac1{2\pi i}
\int_L
g(x,p;x',t')e^{pt}
\,dp.
$$

应用留数定理算出这个积分，也可以得到和 (19.68) 相同的结果。

再讨论一个三维空间的例子。这时的 Green 函数 $G(\boldsymbol r,t;\boldsymbol r',t')$ 满足定解问题

$$
\left(
\frac{\partial^2}{\partial t^2}
-
a^2\nabla^2
\right)
G(\boldsymbol r,t;\boldsymbol r',t')
=
\delta(\boldsymbol r-\boldsymbol r')
\delta(t-t'),
\qquad
t,t'>0,
\tag{19.72a}
$$

$$
\left.
G(\boldsymbol r,t;\boldsymbol r',t')
\right|_{t<t'}
=
0,
\qquad
\left.
\frac{\partial G(\boldsymbol r,t;\boldsymbol r',t')}
{\partial t}
\right|_{t<t'}
=
0.
\tag{19.72b}
$$

在 Fourier 变换

$$
g(\boldsymbol k,t;\boldsymbol r',t')
=
\frac1{(2\pi)^{3/2}}
\iiint
G(\boldsymbol r,t;\boldsymbol r',t')
e^{-i\boldsymbol k\cdot\boldsymbol r}
\,d\boldsymbol r
\tag{19.73}
$$

之下，定解问题 (19.72) 化为

$$
\left[
\frac{d^2}{dt^2}
+
(ka)^2
\right]
g(\boldsymbol k,t;\boldsymbol r',t')
=
\frac1{(2\pi)^{3/2}}
e^{-i\boldsymbol k\cdot\boldsymbol r'}
\delta(t-t'),
\tag{19.74a}
$$

$$
\left.
g(\boldsymbol k,t;\boldsymbol r',t')
\right|_{t<t'}
=
0,
\qquad
\left.
\frac{dg(\boldsymbol k,t;\boldsymbol r',t')}
{dt}
\right|_{t<t'}
=
0.
\tag{19.74b}
$$

根据例 10.5 的结果（(10.42) 式），可以得到

$$
g(\boldsymbol k,t;\boldsymbol r',t')
=
\frac1{(2\pi)^{3/2}}
\frac{
\sin ka(t-t')
}{
ka
}
e^{-i\boldsymbol k\cdot\boldsymbol r'}
\eta(t-t').
\tag{19.75}
$$

作逆变换，就有

$$
G(\boldsymbol r,t;\boldsymbol r',t')
=
\frac{\eta(t-t')}{(2\pi)^3}
\iiint
\frac{
\sin ka(t-t')
}{
ka
}
e^{-i\boldsymbol k\cdot(\boldsymbol r-\boldsymbol r')}
\,d\boldsymbol k.
$$

完全模仿第 §19.3 中的做法，在 $\boldsymbol k$ 空间的球坐标系中计算上面的积分，就得到

$$
G(\boldsymbol r,t;\boldsymbol r',t')
=
\frac1{2\pi^2a}
\frac{\eta(t-t')}
{|\boldsymbol r-\boldsymbol r'|}
\int_0^\infty
\sin ka(t-t')
\sin k|\boldsymbol r-\boldsymbol r'|
\,dk.
\tag{19.76}
$$

这个积分在通常的意义下是不存在的。出现这种积分的原因，从根本上说，是由于现在定解问题 (19.72) 中的非齐次项也不是通常意义下的函数。为了算出这个积分，可以将第十八章中的 (18.47) 式代入 (18.48) 式，

$$
f(x)
=
\sqrt{\frac2\pi}
\int_0^\infty
\left[
\sqrt{\frac2\pi}
\int_0^\infty
f(x')
\sin kx'
\,dx'
\right]
\sin kx
\,dk
=
\frac2\pi
\int_0^\infty
f(x')
\left[
\int_0^\infty
\sin kx'\sin kx
\,dk
\right]
dx',
$$

这说明

$$
\frac2\pi
\int_0^\infty
\sin kx'\sin kx
\,dk
=
\delta(x-x'),
\qquad
x,x'>0.
$$

把这个结果代入 (19.76) 式中，就求出了 Green 函数

$$
G(\boldsymbol r,t;\boldsymbol r',t')
=
\frac1{4\pi a}
\frac1{|\boldsymbol r-\boldsymbol r'|}
\delta
\left(
|\boldsymbol r-\boldsymbol r'|
-
a(t-t')
\right).
\tag{19.77}
$$

这里去掉了函数 $\eta(t-t')$，因为 $\delta$ 函数已经保证了 $t-t'<0$ 时 $G(\boldsymbol r,t;\boldsymbol r',t')=0$。这个解式的物理意义明确：$t'$ 时刻在 $\boldsymbol r'$ 处发射的信号，$t$ 时刻一定只到达距 $\boldsymbol r'$ 点为 $a(t-t')$ 的球面上。

利用这个 Green 函数，当然就可以得到三维无界空间中波动方程的初值问题

$$
\frac{\partial^2u(\boldsymbol r,t)}{\partial t^2}
-
a^2\nabla^2u(\boldsymbol r,t)
=
f(\boldsymbol r,t),
\qquad
t>0,
\tag{19.78a}
$$

$$
\left.
u(\boldsymbol r,t)
\right|_{t=0}
=
\phi(\boldsymbol r),
\qquad
\left.
\frac{\partial u(\boldsymbol r,t)}{\partial t}
\right|_{t=0}
=
\psi(\boldsymbol r)
\tag{19.78b}
$$

的解

$$
\begin{aligned}
u(\boldsymbol r,t)
&=
\frac1{4\pi a^2}
\iiint_{|\boldsymbol r'-\boldsymbol r|<at}
\frac{
f\left(
\boldsymbol r',
t-
|\boldsymbol r'-\boldsymbol r|/a
\right)
}{
|\boldsymbol r'-\boldsymbol r|
}
\,d\boldsymbol r'\\
&\quad+
\frac1{4\pi a}
\left[
\iint_\Sigma
\frac{
\psi(\boldsymbol r')
}{
|\boldsymbol r'-\boldsymbol r|
}
\,d\Sigma'
+
\frac{\partial}{\partial t}
\iint_\Sigma
\frac{
\phi(\boldsymbol r')
}{
|\boldsymbol r'-\boldsymbol r|
}
\,d\Sigma'
\right],
\end{aligned}
\tag{19.79}
$$

其中 $\Sigma$ 是以 $\boldsymbol r$ 点为球心、$at$ 为半径的球面 $|\boldsymbol r'-\boldsymbol r|=at$。这个结果的证明留给读者完成。

## *§19.6 热传导方程的 Green 函数

对于热传导问题的 Green 函数，可以仿照波动问题的做法。例如，对于三维有界空间中的热传导问题

$$
\frac{\partial u(\boldsymbol r,t)}{\partial t}
-
\kappa\nabla^2u(\boldsymbol r,t)
=
f(\boldsymbol r,t),
\qquad
\boldsymbol r\in V,\quad t>0,
\tag{19.80a}
$$

$$
\left.
u(\boldsymbol r,t)
\right|_\Sigma
=
\mu(\Sigma,t),
\qquad
t>0,
\tag{19.80b}
$$

$$
\left.
u(\boldsymbol r,t)
\right|_{t=0}
=
\phi(\boldsymbol r),
\qquad
\boldsymbol r\in V,
\tag{19.80c}
$$

相应的 Green 函数 $G(\boldsymbol r,t;\boldsymbol r',t')$ 就可以定义为 $t'$ 时刻在 $\boldsymbol r'$ 处有一个（单位强度的）瞬时点热源所产生的温度分布，换句话说，就是定解问题

$$
\left(
\frac{\partial}{\partial t}
-
\kappa\nabla^2
\right)
G(\boldsymbol r,t;\boldsymbol r',t')
=
\delta(\boldsymbol r-\boldsymbol r')
\delta(t-t'),
\qquad
\boldsymbol r,\boldsymbol r'\in V,\quad
t,t'>0,
\tag{19.81a}
$$

$$
\left.
G(\boldsymbol r,t;\boldsymbol r',t')
\right|_\Sigma
=
0,
\qquad
t,t'>0,
\tag{19.81b}
$$

$$
\left.
G(\boldsymbol r,t;\boldsymbol r',t')
\right|_{t=0}
=
0,
\qquad
\boldsymbol r,\boldsymbol r'\in V
\tag{19.81c}
$$

的解。容易理解，初始条件也可以写成 $\left.G(\boldsymbol r,t;\boldsymbol r',t')\right|_{t<t'}=0$。

仿照波动问题的做法，也可以证明 Green 函数同样具有空间对称性和时间倒易性

$$
G(\boldsymbol r,t;\boldsymbol r',t')
=
G(\boldsymbol r',-t';\boldsymbol r,-t).
\tag{19.82}
$$

请读者完成这个证明。

为了应用 Green 函数方法解定解问题 (19.80)，也必须模仿上一节的做法：首先将定解问题 (19.80) 的自变量改写成 $\boldsymbol r'$ 和 $t'$，

$$
\frac{\partial u(\boldsymbol r',t')}{\partial t'}
-
\kappa\nabla'^2u(\boldsymbol r',t')
=
f(\boldsymbol r',t'),
\qquad
\boldsymbol r'\in V,\quad
t'>0,
\tag{19.83a}
$$

$$
\left.
u(\boldsymbol r',t')
\right|_{\Sigma'}
=
\mu(\Sigma',t'),
\qquad
t'>0,
\tag{19.83b}
$$

$$
\left.
u(\boldsymbol r',t')
\right|_{t'=0}
=
\phi(\boldsymbol r'),
\qquad
\boldsymbol r'\in V,
\tag{19.83c}
$$

然后写出 Green 函数 $G(\boldsymbol r',-t';\boldsymbol r,-t)$ 满足的定解问题

$$
\left[
\frac{\partial}{\partial(-t')}
-
\kappa\nabla'^2
\right]
G(\boldsymbol r',-t';\boldsymbol r,-t)
=
\delta(\boldsymbol r-\boldsymbol r')
\delta(t-t'),
\qquad
\boldsymbol r,\boldsymbol r'\in V,\quad
t,t'>0,
$$

$$
\left.
G(\boldsymbol r',-t';\boldsymbol r,-t)
\right|_{\Sigma'}
=
0,
\qquad
t,t'>0,
$$

$$
\left.
G(\boldsymbol r',-t';\boldsymbol r,-t)
\right|_{-t'<-t}
=
0,
\qquad
\boldsymbol r,\boldsymbol r'\in V.
$$

进一步再利用 Green 函数对于空间的对称性和时间的倒易性关系 (19.82)，将定解问题改写成

$$
\left(
-\frac{\partial}{\partial t'}
-
\kappa\nabla'^2
\right)
G(\boldsymbol r,t;\boldsymbol r',t')
=
\delta(\boldsymbol r-\boldsymbol r')
\delta(t-t'),
\qquad
\boldsymbol r,\boldsymbol r'\in V,\quad
t,t'>0,
\tag{19.84a}
$$

$$
\left.
G(\boldsymbol r,t;\boldsymbol r',t')
\right|_{\Sigma'}
=
0,
\qquad
t,t'>0,
\tag{19.84b}
$$

$$
\left.
G(\boldsymbol r,t;\boldsymbol r',t')
\right|_{t'>t}
=
0,
\qquad
\boldsymbol r,\boldsymbol r'\in V.
\tag{19.84c}
$$

将方程 (19.83a) 和 (19.84a) 分别乘以 $G(\boldsymbol r,t;\boldsymbol r',t')$ 和 $u(\boldsymbol r',t')$，相减并积分，就得到

$$
\begin{aligned}
&
\int_0^\infty
dt'
\iiint_V
f(\boldsymbol r',t')
G(\boldsymbol r,t;\boldsymbol r',t')
\,d\boldsymbol r'
-
u(\boldsymbol r,t)\\
&=
\iiint_V
d\boldsymbol r'
\int_0^\infty
\left[
G(\boldsymbol r,t;\boldsymbol r',t')
\frac{\partial u(\boldsymbol r',t')}{\partial t'}
+
u(\boldsymbol r',t')
\frac{\partial G(\boldsymbol r,t;\boldsymbol r',t')}{\partial t'}
\right]dt'\\
&\quad-
\kappa
\int_0^\infty
dt'
\iiint_V
\left[
G(\boldsymbol r,t;\boldsymbol r',t')
\nabla'^2u(\boldsymbol r',t')
-
u(\boldsymbol r',t')
\nabla'^2G(\boldsymbol r,t;\boldsymbol r',t')
\right]
d\boldsymbol r'\\
&=
\iiint_V
\left.
G(\boldsymbol r,t;\boldsymbol r',t')
u(\boldsymbol r',t')
\right|_{t'=0}^{t'=\infty}
d\boldsymbol r'\\
&\quad-
\kappa
\int_0^\infty
dt'
\iint_{\Sigma'}
\left[
G(\boldsymbol r,t;\boldsymbol r',t')
\nabla'u(\boldsymbol r',t')
-
u(\boldsymbol r',t')
\nabla'G(\boldsymbol r,t;\boldsymbol r',t')
\right]\cdot d\boldsymbol\Sigma'.
\end{aligned}
$$

代入边界条件 (19.83b)、(19.84b) 和初始条件 (19.83c)、(19.84c)，最后就得到

$$
\begin{aligned}
u(\boldsymbol r,t)
&=
\int_0^t
dt'
\iiint_V
f(\boldsymbol r',t')
G(\boldsymbol r,t;\boldsymbol r',t')
\,d\boldsymbol r'\\
&\quad+
\iiint_V
\phi(\boldsymbol r')
G(\boldsymbol r,t;\boldsymbol r',0)
\,d\boldsymbol r'\\
&\quad-
\kappa
\int_0^t
dt'
\iint_\Sigma
\mu(\Sigma',t')
\left.
\frac{\partial G(\boldsymbol r,t;\boldsymbol r',t')}
{\partial n'}
\right|_{\Sigma'}
d\Sigma'.
\end{aligned}
\tag{19.85}
$$

**练习 19.2**　如果直接从 (19.80) 和 (19.81) 式出发，而不将这些方程的自变量换成 $\boldsymbol r'$ 和 $t'$，是否能够用 Green 函数以及已知条件 $f(\boldsymbol r,t)$、$\mu(\Sigma,t)$ 和 $\phi(\boldsymbol r)$ 将定解问题 (19.80) 的解 $u(\boldsymbol r,t)$ 表示出来？

现在来求解第十二章中遗留下的一个定解问题，即关于一维无界空间热传导方程的 Green 函数问题

$$
\left(
\frac{\partial}{\partial t}
-
\kappa\frac{\partial^2}{\partial x^2}
\right)
G(x,t;x',t')
=
\delta(x-x')
\delta(t-t'),
\qquad
t'>0,
\tag{19.86a}
$$

$$
\left.
G(x,t;x',t')
\right|_{x\to\pm\infty}
\text{ 有界},
\tag{19.86b}
$$

$$
\left.
G(x,t;x',t')
\right|_{t=0}
=
0.
\tag{19.86c}
$$

作 Laplace 变换，即令

$$
g(p;x,x')
=
\int_0^\infty
G(x,t;x',t')
e^{-pt}
\,dt,
\tag{19.87}
$$

于是，定解问题化为

$$
pg(p;x,x')
-
\kappa
\frac{d^2g(p;x,x')}{dx^2}
=
\delta(x-x')e^{-pt'},
\tag{19.88a}
$$

$$
\left.
g(p;x,x')
\right|_{x\to\pm\infty}
\text{ 有界}.
\tag{19.88b}
$$

援引例 10.8 中的结果 (10.76) 式，就求得

$$
g(p;x,x')
=
\frac1{2\sqrt{\kappa p}}
e^{-pt'}
\exp
\left\{
-\sqrt{\frac p\kappa}
|x-x'|
\right\}.
\tag{19.89}
$$

再利用 (8.40) 式的结果，求出反演，得

$$
G(x,t;x',t')
=
\frac1{
2\sqrt{\kappa\pi(t-t')}
}
\exp
\left\{
-\frac{(x-x')^2}
{4\kappa(t-t')}
\right\}
\eta(t-t').
\tag{19.90}
$$

**练习 19.3**　求解定解问题 (12.54)。

## 习题

1. 氢原子中的电子处于 $1s$ 轨道上时，其电荷密度为

$$
\rho(r,\theta,\phi)
=
\frac{q_0}{\pi a_0^3}
e^{-2r/a_0},
$$

其中 $q_0$ 为电子电荷，$a_0$ 为 Bohr 半径。应用电势叠加原理，求空间任意一点处的静电势。

提示：积分时可能会用到 (15.88) 式。

2. (1) 用电像法求出球内 Laplace 方程第一类边值问题的 Green 函数 $G(\boldsymbol r;\boldsymbol r')$；

(2) 求出边界面（球面 $r=a$）上各点的感生电荷密度 $\sigma(\theta,\phi)$；

(3) 证明像电荷和感生电荷在球内完全等效；

(4) 证明球内 Laplace 方程第一类边值问题

$$
\nabla^2u=0,
$$

$$
\left.u\right|_{r=a}
=
f(\theta,\phi)
$$

的解是

$$
u(r,\theta,\phi)
=
\frac{
a(a^2-r^2)
}{4\pi}
\int_0^{2\pi}
\left[
\int_0^\pi
\frac{
f(\theta',\phi')
}{
(a^2+r^2-2ar\cos\psi)^{3/2}
}
\sin\theta'
\,d\theta'
\right]
d\phi',
$$

其中 $\psi$ 是 $\boldsymbol r(r,\theta,\phi)$ 与 $\boldsymbol r'(r',\theta',\phi')$ 的夹角，

$$
\cos\psi
=
\cos\theta\cos\theta'
+
\sin\theta\sin\theta'
\cos(\phi-\phi').
$$

3. 利用电像法求三维 Poisson 方程球外问题

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
+
\frac{\partial^2u}{\partial z^2}
=
f(x,y,x),
\qquad
x^2+y^2+z^2>a^2,
$$

$$
\left.
u
\right|_{x^2+y^2+z^2=a^2}
=
g
$$

的 Green 函数，其中 $f$ 和 $g$ 是已知函数。

4. 用 Green 函数方法解第十三章习题第 6 题。

5. 一无穷长弦，$t=t_0$ 时在 $x=x_0$ 处受到瞬时的打击，冲量为 $I$。试求解弦的横振动，设初位移和初速度均为 0。

6. 用 Green 函数方法解无界弦的横振动问题，定解问题为

$$
\frac{\partial^2u}{\partial t^2}
-
a^2
\frac{\partial^2u}{\partial x^2}
=
0,
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
\psi(x).
$$

7. 用 §19.5 中求解三维无界空间波动方程 Green 函数的方法，求解三维无界空间热传导方程的 Green 函数。

8. 证明三维无界空间热传导问题 Green 函数 $G(\boldsymbol r,t;\boldsymbol r',t')$ 的空间对称性和时间倒易性

$$
G(\boldsymbol r,t;\boldsymbol r',t')
=
G(\boldsymbol r',-t';\boldsymbol r,-t).
$$

9. 用 Green 函数方法求解三维无界空间热传导问题

$$
\frac{\partial u}{\partial t}
-
\kappa\nabla^2u
=
f(x,y,z,t),
$$

$$
\left.u\right|_{t=0}
=
\phi(x,y,y).
$$
