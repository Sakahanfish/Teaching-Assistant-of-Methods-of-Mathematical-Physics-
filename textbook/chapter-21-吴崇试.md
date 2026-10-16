# 第二十一章 数学物理方程综述

## §21.1 二阶线性偏微分方程的分类

在本课程的数学物理方程部分中，我们总共讨论了三种类型偏微分方程定解问题的解。这三类方程（即波动方程，热传导方程，以及稳定问题中出现的 Laplace 方程或 Poisson 方程，包括分离变量导出的 Helmholtz 方程）描写了不同的物理过程，在数学上，也分属双曲型、抛物型和椭圆型三类（见 §11.4）。它们的解也都表现出各自不同的特点（见 §12.5～§12.8 的讨论）。下面证明：在两个自变量的情形下，二阶线性偏微分方程就只有这三种类型①。

两个自变量 $(x,y)$ 的二阶线性偏微分方程的普遍形式是：

$$
a\frac{\partial^2u}{\partial x^2}
+2b\frac{\partial^2u}{\partial x\partial y}
+c\frac{\partial^2u}{\partial y^2}
+d\frac{\partial u}{\partial x}
+e\frac{\partial u}{\partial y}
+fu+g=0,
\tag{21.1}
$$

其中 $a,b,c,d,e,f$ 和 $g$ 是 $x,y$ 的已知函数。通常假设它们是连续可微的。显然，函数 $a,b,c$ 中至少有一个不恒为 0，否则就不成其为二阶偏微分方程。

首先考虑 $a$ 和（或）$c$ 不恒为 0 的情形。不妨设 $a\ne0$。这时可作变换

$$
\xi=\phi(x,y),
\qquad
\eta=\psi(x,y).
\tag{21.2}
$$

为了保证 $\xi$ 和 $\eta$ 仍然是独立变量，这一组变换必须满足

$$
\frac{\partial(\xi,\eta)}{\partial(x,y)}
\ne0.
\tag{21.3}
$$

在这一组变换下，方程 (21.1) 变为

$$
A\frac{\partial^2u}{\partial\xi^2}
+2B\frac{\partial^2u}{\partial\xi\partial\eta}
+C\frac{\partial^2u}{\partial\eta^2}
+D\frac{\partial u}{\partial\xi}
+E\frac{\partial u}{\partial\eta}
+Fu+G=0,
\tag{21.4}
$$

其中

$$
A
=
a\left(\frac{\partial\phi}{\partial x}\right)^2
+2b\frac{\partial\phi}{\partial x}\frac{\partial\phi}{\partial y}
+c\left(\frac{\partial\phi}{\partial y}\right)^2,
$$

$$
B
=
a\frac{\partial\phi}{\partial x}\frac{\partial\psi}{\partial x}
+b\left(
\frac{\partial\phi}{\partial x}\frac{\partial\psi}{\partial y}
+
\frac{\partial\phi}{\partial y}\frac{\partial\psi}{\partial x}
\right)
+c\frac{\partial\phi}{\partial y}\frac{\partial\psi}{\partial y},
$$

$$
C
=
a\left(\frac{\partial\psi}{\partial x}\right)^2
+2b\frac{\partial\psi}{\partial x}\frac{\partial\psi}{\partial y}
+c\left(\frac{\partial\psi}{\partial y}\right)^2,
$$

$$
D
=
a\frac{\partial^2\phi}{\partial x^2}
+2b\frac{\partial^2\phi}{\partial x\partial y}
+c\frac{\partial^2\phi}{\partial y^2}
+d\frac{\partial\phi}{\partial x}
+e\frac{\partial\phi}{\partial y},
$$

$$
E
=
a\frac{\partial^2\psi}{\partial x^2}
+2b\frac{\partial^2\psi}{\partial x\partial y}
+c\frac{\partial^2\psi}{\partial y^2}
+d\frac{\partial\psi}{\partial x}
+e\frac{\partial\psi}{\partial y},
$$

$$
F=f,
\qquad
G=g.
$$

容易证明

$$
B^2-AC
=
\left(
\frac{\partial\phi}{\partial x}\frac{\partial\psi}{\partial y}
-
\frac{\partial\phi}{\partial y}\frac{\partial\psi}{\partial x}
\right)^2
(b^2-ac)
=
\left|
\frac{\partial(\xi,\eta)}{\partial(x,y)}
\right|^2
(b^2-ac).
\tag{21.5}
$$

为了书写简便起见，令

$$
\Phi
\left(
\xi,\eta,u,
\frac{\partial u}{\partial\xi},
\frac{\partial u}{\partial\eta}
\right)
\equiv
D\frac{\partial u}{\partial\xi}
+
E\frac{\partial u}{\partial\eta}
+
Fu+G,
\tag{21.6}
$$

则方程 (21.4) 变为

$$
A\frac{\partial^2u}{\partial\xi^2}
+
2B\frac{\partial^2u}{\partial\xi\partial\eta}
+
C\frac{\partial^2u}{\partial\eta^2}
+
\Phi
\left(
\xi,\eta,u,
\frac{\partial u}{\partial\xi},
\frac{\partial u}{\partial\eta}
\right)
=0.
\tag{21.4'}
$$

我们希望，通过适当选择变换 (21.2)，使得 $A,B,C$ 中有一个或几个为 0，达到使方程简化的目的。为此，介绍一个定理。

**定理 21.1**　如果 $\phi(x,y)=C$ 是方程

$$
a(dy)^2-2b\,dy\,dx+c(dx)^2=0
\tag{21.7}
$$

的一般积分，则 $\xi=\phi(x,y)$ 是方程

$$
a\left(\frac{\partial\phi}{\partial x}\right)^2
+
2b\frac{\partial\phi}{\partial x}\frac{\partial\phi}{\partial y}
+
c\left(\frac{\partial\phi}{\partial y}\right)^2
=0
\tag{21.8}
$$

的一个特解。

**证**　因为 $\phi(x,y)=C$，故有

$$
\frac{\partial\phi}{\partial x}\,dx
+
\frac{\partial\phi}{\partial y}\,dy
=
0,
\qquad
\text{即}
\qquad
dy
=
-\frac{\partial\phi/\partial x}{\partial\phi/\partial y}\,dx.
$$

不妨设 $\partial\phi/\partial y\ne0$。代入方程 (21.7)，就有

$$
\begin{aligned}
a(dy)^2-2b\,dy\,dx+c(dx)^2
&=
\left[
a\left(\frac{\partial\phi}{\partial x}\right)^2
+
2b\frac{\partial\phi}{\partial x}\frac{\partial\phi}{\partial y}
+
c\left(\frac{\partial\phi}{\partial y}\right)^2
\right]
\left(\frac{\partial\phi}{\partial y}\right)^{-2}
(dx)^2\\
&=0,
\end{aligned}
$$

所以 (21.8) 成立。□

这个定理告诉我们，如果选择变换 $\xi=\phi(x,y)$ 使 $A=0$，或是选择变换 $\eta=\psi(x,y)$ 使 $C=0$，就可以通过求解常微分方程

$$
a\left(\frac{dy}{dx}\right)^2
-
2b\frac{dy}{dx}
+c
=
0,
\qquad
\text{或}
\qquad
\frac{dy}{dx}
=
\frac ba
\pm
\frac1a\sqrt{b^2-ac}
\tag{21.9}
$$

来得到。在一般情况下，这样能得到两个无关解，称为偏微分方程 (21.1) 的特征线。

在具体求解方程 (21.9) 时，又需要区别下列三种情形：

1. $b^2-ac>0$。这时，从方程 (21.9) 可以求得两个实值函数解

$$
\phi(x,y)=C_1
\qquad
\text{及}
\qquad
\psi(x,y)=C_2,
$$

也就是说，偏微分方程 (21.1) 有两条实的特征线。于是，令

$$
\xi=\phi(x,y),
\qquad
\eta=\psi(x,y),
$$

就可使 $A=C=0$。同时，根据 (21.5) 式，还可以断定 $B$ 一定不为 0。所以，方程 (21.4') 就变成

$$
\frac{\partial^2u}{\partial\xi\partial\eta}
+
\Phi_1
\left(
\xi,\eta,u,
\frac{\partial u}{\partial\xi},
\frac{\partial u}{\partial\eta}
\right)
=
0.
\tag{21.10}
$$

或者进一步作变换

$$
\rho=\xi+\eta,
\qquad
\sigma=\xi-\eta,
$$

于是方程 (21.10) 可以化为

$$
\frac{\partial^2u}{\partial\rho^2}
-
\frac{\partial^2u}{\partial\sigma^2}
+
\Phi_2
\left(
\rho,\sigma,u,
\frac{\partial u}{\partial\rho},
\frac{\partial u}{\partial\sigma}
\right)
=
0.
\tag{21.11}
$$

这种类型的方程称为**双曲型方程**。波动方程就属于这种类型。

2. $b^2-ac<0$。这时可以重复上面的讨论，只不过得到的 $\phi(x,y)$ 和 $\psi(x,y)$ 是一对共轭的复函数，或者说，偏微分方程 (21.1) 的两条特征线都不是实的。于是

$$
\xi=\phi(x,y),
\qquad
\eta=\psi(x,y)
$$

是一对共轭的复变量。这样也能够得到以复变量 $\xi$ 和 $\eta$ 为自变量的方程 (21.10)。进一步引进两个新的实变量

$$
\rho=\xi+\eta,
\qquad
\sigma=i(\xi-\eta),
$$

方程 (21.10) 也可以进一步化为

$$
\frac{\partial^2u}{\partial\rho^2}
+
\frac{\partial^2u}{\partial\sigma^2}
+
\Phi_3
\left(
\rho,\sigma,u,
\frac{\partial u}{\partial\rho},
\frac{\partial u}{\partial\sigma}
\right)
=
0.
\tag{21.12}
$$

称为**椭圆型方程**。Laplace 方程、Poisson 方程和 Helmholtz 方程都属于这种类型。

3. $b^2-ac=0$。这时，方程 (21.9) 一定有重根

$$
\frac{dy}{dx}
=
\frac ba,
$$

因而只能求得一个解，例如 $\phi(x,y)=C$。作变换 $\xi=\phi(x,y)$ 就可以使 $A=0$。但是，由 (21.5) 式可以断定，一定有 $B^2-AC=0$，这意味着 $B$ 也一定为 0。所以，我们完全可以任意选取另一个变换，$\eta=\psi(x,y)$，只要它和 $\xi=\phi(x,y)$ 彼此独立，即

$$
\frac{\partial(\xi,\eta)}{\partial(x,y)}
\ne0,
$$

这样，方程 (21.4') 就化为

$$
\frac{\partial^2u}{\partial\eta^2}
+
\Phi_4
\left(
\xi,\eta,u,
\frac{\partial u}{\partial\xi},
\frac{\partial u}{\partial\eta}
\right)
=
0.
\tag{21.13}
$$

这种类型的方程称为**抛物型方程**。热传导方程就属于这种类型。

以上的讨论是在 $a$ 和 $c$ 不恒为 0 的前提下进行的。适当选择变换 (21.2)，总可以使 $A,B,C$ 中有一个（$B$）或两个（$A,C$ 或 $B,C$）为 0。再作进一步的变换，还可以把不为 0 的系数变为 1 或 $-1$。当 $A=C=1,\ B=0$ 时，方程是椭圆型；当 $A=-C=\pm1,\ B=0$ 时，方程为双曲型；当 $A=B=0,\ C=1$ 或 $A=1,\ B=C=0$ 时，方程为抛物型。

如果 $a$ 和 $c$ 恒为 0。那么，一定有 $b\ne0$。这正属于双曲型方程，不必再讨论。

综合以上的讨论，可以得出结论：要判断二阶线性偏微分方程属于何种类型，只要讨论判别式 $b^2-ac$ 即可。如果方程的系数 $a,b,c$ 为常数，偏微分方程一定属于上述三种类型之一。如果 $a,b,c$ 是 $x,y$ 的函数，那么，在 $x-y$ 平面上的一定区域内，一般说来，$b^2-ac$ 并不会保持为恒正、恒负、或恒为 0，因此方程并不能简单地归结为固定的一种类型。换句话说，方程可能在区域的不同部分属于不同类型。这时，不妨先求出 $b^2-ac=0$，即

$$
\frac{dy}{dx}
=
\frac ba
$$

的解。这条曲线称为抛物型曲线，因在此曲线上方程属于抛物型。整个区域就可能被这条曲线分割为两部分，方程分属于椭圆型和双曲型。例如，对于方程

$$
(1-x^2)\frac{\partial^2u}{\partial x^2}
-
2xy\frac{\partial^2u}{\partial x\partial y}
-
(1+y^2)\frac{\partial^2u}{\partial y^2}
-
2x\frac{\partial u}{\partial x}
-
2y\frac{\partial u}{\partial y}
=
0,
$$

容易求出

$$
b^2-ac
=
1-x^2+y^2.
$$

因此，此方程的抛物型曲线就是一对双曲线 $x^2-y^2=1$。在双曲线上，方程属于抛物型。整个 $x-y$ 平面被这两条曲线分割开。在 $1-x^2+y^2>0$ 的部分，方程属于双曲型；在 $1-x^2+y^2<0$ 的部分，方程属于椭圆型。

对于多个自变量的偏微分方程，原则上也可以选择适当的自变量变换，把方程中混合二阶偏导数项的系数变为 0。如果其余的（二阶偏导数项的）系数（事实上，可以化为 1 或 $-1$）全部同号，则方程为椭圆型；如果其中一个与其余的异号，则方程为双曲型；如果有多个与其余的异号，则方程为超双曲型；如果有一个或多个为 0，则方程为抛物型。当然，除非方程的系数为常数，否则，自变量变换的具体选择总还需要具体讨论。

① 对于更多个自变量的情形，问题要复杂一些，但讨论的基本方法是一样的。

## §21.2 线性偏微分方程解法述评

在本书中，介绍了二阶线性偏微分方程定解问题的几种主要解法，关于这些解法的解题思想、应用条件以及理论根据，以前也都分别做过讨论，这里再集中地对它们做一点综合性的评述，以便于读者有一个横向的比较。

**1. 分离变量法。** 这是求解线性偏微分方程定解问题的主要方法。从理论上说，分离变量法的依据是 Sturm-Liouville 型方程的本征值问题。这在第十七章中已做了较系统的阐述，不再重复。从解题步骤上看，除了留待确定叠加系数的部分定解条件外，要求其余的边界条件都必须是齐次的（如果它是非齐次的，则首先必须齐次化）。这样，对于定解问题中微分方程的具体形式就有一定的限制，对于所讨论问题的空间区域形状更有明显的限制。这又涉及正交曲面坐标系的选取（空间区域的边界面必须是正交曲面坐标系的坐标面）。关于这些问题，在 P. M. Morse 和 H. Feshbach 的 *Methods of Theoretical Physics* (McGraw-Hill, New York, 1941) 一书中有详细的讨论，见该书的 5.1 节。

在具体求解时，当然还必须求解相应的常微分方程本征值问题。除了本书中介绍过的几个本征函数外，还可能出现其他的特殊函数，读者可以查阅参考书目 [12, 13]。此外，值得查阅的还有 E. Kamke 的 *Differentialgleichungen, Lösungsmethoden und Lösungen, Band 1, Gewöhnliche Differentialgleichungen*。该书相当全地收录了各种常微分方程的解。

**2. 积分变换方法。** 这种方法的优点是减少方程的微分变量的数目。从原则上说，无论是对于时间变量或是空间变量，无论是无界空间或是有界空间，都可以采用积分变换的方法求解线性偏微分方程的定解问题。但从实际计算看，就需要根据方程和定解条件的类型，选择最合适的积分变换。反演问题，也是关系所拟采用的积分变换是否实际可行的关键问题。如果反演时涉及的积分很简单，甚至有现成的结果（包括工具书）可供引用，采用积分变换的确可以带来极大的便利。但如果涉及的积分比较复杂，也没有现成的结果（包括工具书）可供引用，那么，反演问题也可以成为积分变换的难点。

积分变换方法和分离变量法存在密切的联系。例如，当本征值过渡到连续谱时，分离变量法就变为相应的积分变换方法。

从实用的角度说，如果是有界空间，一般说来，积分变换和分离变量法没有什么差别，而且只要本征值问题能解出，那么分离变量法没有什么特别的难点，故仍不妨采用分离变量法。

积分变换方法不要求定解问题中的边界条件是齐次的，不需要先将非齐次边界条件齐次化，另外积分变换方法还具有分离变量法所没有的优点：它还可用于求解非线性偏微分方程。

**3. Green 函数方法。** 应该说，这种方法具有极大的理论意义。它给出了定解问题的解和方程的非齐次项以及定解条件之间的关系，因而便于讨论方程的非齐次项或定解条件发生变化时，解如何相应地变化。而且，不止如此，在讨论本征值问题的普遍性质时，也离不开 Green 函数，只不过在本书中未做介绍而已。Green 函数方法已经成为理论物理研究中的常用方法之一。

应用 Green 函数方法，最重要的是，要能够求出 Green 函数的具体形式。尽管 Green 函数所满足的是一种特别简单的定解问题，方程的非齐次项为 $\delta$ 函数，定解条件均为齐次，在少数情形下，能够求得 Green 函数的简单表达式，但是，一般说来，可求出 Green 函数的解析解的情形，仍只限于若干种空间区域形状，和分离变量法没有什么差别。

Green 函数方法的另一个优点是便于进行近似计算。例如，对于某一类偏微分方程的定解问题，由于区域形状的限制，不能求出它的 Green 函数的解析表达式。但是，如果必要的话，总还可以求出 Green 函数的足够精确的近似解（例如数值解）。这样，也就可以进一步求出这一类偏微分方程定解问题的近似解。这在工程上还是具有实际意义的。

**4. 变分法。** 这个方法具有理论价值和实用价值。在理论上，它可以把不同类型的偏微分方程定解问题用相同的泛函语言表达出来（当然不同问题中出现的泛函是不同的），或者说，把不同的物理问题用相同的泛函语言表达出来。正是由于这个原因，变分或泛函语言已经成为表述物理规律的常用工具之一。在实用上，变分法又提供了一种近似计算的好办法。有效地利用物理知识，灵活巧妙地选取试探函数，可以使计算大为简化。在第二十章中，我们已经看到过这样的例子。在物理学中，过去或现在，变分法都是常用的一种近似计算方法。例如，在原子和分子光谱的计算中，就广泛地采用了变分法。

**5. Poisson 积分公式。** 对于二维和三维 Laplace 方程的边值问题，也还可以将解表示为特殊的积分公式。对于二维 Laplace 方程，它的解一定是解析函数的实部或虚部，因此，可以采用复变函数的方法求解。例如，圆内或上半平面的第一类边值问题，Laplace 方程的解就可以表示为 Poisson 积分（见 §3.7，也可以从 Green 函数方法得到，见 §19.4）。三维 Laplace 方程第一类边值问题的解，也可以表示为沿边界面的积分。

除了上面提到的这几种方法外，还有下面几种方法。

**6. 保角变换。** 这种方法的理论基础，是解析函数所代表的变换的保角性。这种解法，主要用于二维 Laplace 方程或 Poisson 方程的边值问题，因为在保角变换下，前者的形式不变，后者也只是非齐次项做相应的改变。粗略地说，运用保角变换，可以把“不规则”的边界形状化为规则的边界形状（但是难以在“不规则”和“规则”之间划定一个界限），例如，可以把多边形化为上半平面或单位圆内。再结合上半平面或圆内的 Poisson 公式，就能直接求出二维 Laplace 方程的解。运用保角变换，的确可以解决一些有意义的物理问题或工程问题，例如，有限大小尺寸的平行板电容器的边缘效应问题，空气动力学中的机翼问题，以及某些流体力学问题。又如，应用保角变换方法，可以把偏心圆化为同心圆。有兴趣的读者，可参阅参考书目 [1, 2]。在 H. Kober 所编的 *Dictionary of Conformal Representations* (Dover Publications, Inc., 1957) 一书中，收录了各种主要的保角变换，包括初等函数和椭圆函数所代表的保角变换，也可供参考。

**7. 双曲型方程定解问题的特殊解法，** 例如平均值法，降维法，等等。在理论上说，双曲型方程的解的存在唯一性，可以通过所谓 Cauchy 型边界条件（即要求解在边界上同时满足给定的函数值与法向微商值）得到保证①。相应地，双曲型方程，就可以采用特征线法（或称 Riemann 方法）求解。由于篇幅限制，这些方法都未做介绍。读者也可参阅参考书目 [1]。

① 椭圆型方程就不同。对于椭圆型方程，只要指定未知函数在边界上的函数值或法向微商值，就足以唯一地确定解。同时指定未知函数在边界上的函数值和法向微商值，反而是过分了，会造成问题无解。


## *§21.3 非线性偏微分方程问题

本书中讨论的偏微分方程定解问题，全部都是由线性方程和线性定解条件构成的。这一类问题的解法特别简单，因为可以运用叠加原理。从实际问题看，这是和物理学的发展状况密切相关的。迄今为止，线性近似仍然是物理学中大量采用的最基本的近似。例如，在 Newton 力学中，质点的加速度与外力成正比，比例系数（质量）是常数，与质点运动的速度大小无关。在弹性力学中，在弹性限度内，应力与应变成正比（Hooke 定律），比例系数（弹性系数）是物质常数，与应变的大小无关。又如，在涉及输运过程的分子动理论中，也是着重讨论相对于平衡状态的线性偏离：由温度的分布不均匀而产生热传导现象，热流密度与温度梯度成正比，比例系数（导热率）是物质常数，与温度高低无关；由物质密度的分布不均匀而产生扩散现象，物质流密度（单位时间通过单位面积的质量）与密度梯度成正比，比例系数（扩散系数）是常数，与物质密度的高低无关。在电磁学中，Ohm 定律说的也是电流密度与电场强度成正比，比例系数（电导率）是物质常数，与电场强度的高低无关。这类例子，在物理学中，可以说俯拾皆是。相应地，在描写连续介质或场的运动的数学物理方程中，就出现了波动方程、热传导方程和 Laplace 方程、Poisson 方程、Helmholtz 方程等线性偏微分方程，以及各种类型的线性定解条件。正是由于采用了线性近似，所以得到的方程形式具有普适性。无论是弹性体中发生的纵振动或横振动，或是电磁场随时间、空间的变化与分布，都遵从同样形式的波动方程，介质的性质只体现在波的传播速率上。无论是热传导过程，或是扩散过程，也都遵从同样的热传导方程，不同的过程，以及有关的介质性质，同样也只表现在方程中的常数（扩散率）上。

对上面各种现象的线性描述，当然都只是在一定限度内的近似。随着科学技术的发展，以及人们对于自然规律认识的深化，不可避免地会超出线性近似的限制。研究各种极端条件（例如，高温、高压、高密度等）下的物理过程，研究物理过程随时间的长期演变，或是在空间上的大尺度范围内的变化，都使得非线性效应变得不可忽略。例如，当介质表面的温度和环境温度相差不大时，单位时间内通过单位表面积散出的热量与温差成正比（Newton 冷却定律），但如果介质表面的温度 $T$ 足够高，热辐射的效应不可忽略，以辐射方式散出的热量便与 $T^4$ 成正比（Stefan-Boltzmann 定律）。

下面再讨论一下无穷直线上的波动问题。正如第十二章中指出的，波动方程

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0
\tag{21.14}
$$

的解

$$
u(x,t)
=
f(x-at)+g(x+at)
\tag{21.15}
$$

表示的是在 $\pm x$ 方向上独立传播的行波。如果只关注其中的一个行波，例如，$u(x,t)=f(x-at)$，它满足的一阶偏微分方程

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}
=
0
\tag{21.16}
$$

可以改写成连续性方程

$$
\frac{\partial u}{\partial t}
+
\frac{\partial j}{\partial x}
=
0,
\tag{21.17}
$$

其中的 $j=au$ 表示“流”（粒子流、能量流等）的强度。如果要考虑非线性的影响，下一级的近似便会有 $u^2$ 项：

$$
j
=
au+\frac{\alpha}{2}u^2.
\tag{21.18}
$$

代入连续性方程 (21.17)，波动方程就变为

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}
+
\alpha u\frac{\partial u}{\partial x}
=
0,
\tag{21.19}
$$

方程中就出现了非线性项。如果同时还存在色散（见 §12.6），流的强度变为

$$
j
=
au
+
\beta\frac{\partial^2u}{\partial x^2}
+
\frac{\alpha}{2}u^2,
\tag{21.20}
$$

代入连续性方程 (21.17)，波动方程又变为

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}
+
\beta\frac{\partial^3u}{\partial x^3}
+
\alpha u\frac{\partial u}{\partial x}
=
0.
\tag{21.21}
$$

为了将方程 (21.21) 的形式化简，可以进一步作变换

$$
\tau=At,
\qquad
\xi=A(x-at),
\qquad
v=Bu,
$$

取 $A^2=1/\beta,\ B=-6/\alpha$，就得到标准的 KdV 方程（Korteweg-de Vries，1895 年）

$$
\frac{\partial v}{\partial\tau}
+
\frac{\partial^3v}{\partial\xi^3}
-
6v\frac{\partial v}{\partial\xi}
=
0.
\tag{21.22}
$$

这是典型的非线性偏微分方程之一。它可以描写浅水波的传播。

在非线性偏微分方程中，经常提到的典型方程还有正弦 Gordon 方程

$$
\frac{\partial^2u}{\partial x\partial t}
=
\sin u,
\qquad
\text{或}
\qquad
\frac{\partial^2u}{\partial x^2}
-
\frac{\partial^2u}{\partial t^2}
=
\sin u
\tag{21.23}
$$

和非线性 Schrödinger 方程

$$
i\frac{\partial u}{\partial t}
=
-\frac{\hbar^2}{2m}
\frac{\partial^2u}{\partial x^2}
+
\alpha|u|^2u.
\tag{21.24}
$$

前者最早出现在 19 世纪的几何问题中。

非线性方程的最大特点，就是解不再具有线性叠加性质。例如，即使对于齐次的非线性方程，如果 $u$ 是方程的解，它的常数倍 $Au$ 也不一定是方程的解；如果 $u_1$ 和 $u_2$ 是方程的解，它们的和 $u_1+u_2$ 也不见得是方程的解。因此，求解非线性方程，需要特殊的技巧。下面就简单介绍 KdV 方程 (21.22) 的几个特解。

为了叙述的方便，不妨撇开 KdV 方程的上述背景，而是简单地把 $\xi$ 和 $\tau$ 仍称为空间和时间变量。最容易求的是

$$
v(\xi,\tau)
=
f(\xi-c\tau)
\tag{21.25}
$$

形式的行波解，因为这样可以转化为常微分方程的求解问题。令 $\eta=\xi-c\tau$，于是

$$
\frac{\partial v}{\partial\tau}
=
-c\frac{df}{d\eta},
\qquad
\frac{\partial v}{\partial\xi}
=
\frac{df}{d\eta},
$$

代入 KdV 方程 (21.22)，得

$$
-c\frac{df}{d\eta}
+
\frac{d^3f}{d\eta^3}
-
6f(\eta)\frac{df}{d\eta}
=
0.
$$

积分一次，有

$$
-cf(\eta)
+
\frac{d^2f}{d\eta^2}
-
3[f(\eta)]^2
=
A,
\tag{21.26}
$$

$A$ 为积分常数。两端乘以 $df/d\eta$，再积分，就得到

$$
-\frac c2[f(\eta)]^2
+
\frac12
\left(
\frac{df}{d\eta}
\right)^2
-
[f(\eta)]^3
=
Af(\eta)+B,
\tag{21.27}
$$

$B$ 是第二个积分常数。如果我们加上边界条件

$$
\eta\to\pm\infty
\text{ 时},
\qquad
f(\eta),
\frac{df}{d\eta},
\frac{d^2f}{d\eta^2}
\text{ 均}\to0,
$$

则可定出 $A=B=0$。于是

$$
\left(
\frac{df}{d\eta}
\right)^2
=
[f(\eta)]^2[2f(\eta)+c],
\qquad
\text{即}
\qquad
\pm
\frac{df}{f\sqrt{2f+c}}
=
d\eta.
\tag{21.28}
$$

这里一定有 $2f(\eta)+c\ge0$。作变换 $\sqrt{2f+c}=\sqrt c\,w$，方程就化为

$$
\mp
\frac{2}{\sqrt c}
\frac{dw}{1-w^2}
=
d\eta.
\tag{21.29}
$$

先考虑上式中取负号的情形。解之即得

$$
-\frac1{\sqrt c}
\ln\frac{1+w}{1-w}
=
\eta-\eta_0,
\qquad
\text{即}
\qquad
\frac{\sqrt c+\sqrt{2f+c}}
{\sqrt c-\sqrt{2f+c}}
=
e^{-\sqrt c(\eta-\eta_0)}.
$$

进一步化简，就得到解

$$
f(\xi-c\tau)
=
-\frac c2
\operatorname{sech}^2
\left\{
\frac{\sqrt c}{2}
[
(\xi-\xi_0)-c(\tau-\tau_0)
]
\right\}.
\tag{21.30}
$$

这是一个行波解，在 $\tau>\tau_0$ 的任意一个时刻，仍然保持 $\tau=\tau_0$ 时刻的波形，只不过向右平移了 $c(\tau-\tau_0)$。在非线性方程中，常把这种不受干扰地传播的波称为孤波，或孤[立]子。(21.30) 式的波形只有一个极值，所以称为单孤波或单孤子。图 21.1 中给出了 $f(\xi-c\tau)$ 的图形。

![Figure 21.1](...)

图 21.1　单孤波

值得注意，与线性波动方程不同，KdV 方程的孤波解的传播速度 $c$ 并不是一个固定的常数，而是任意常数，只要大于 0 即可。对于任意一个 $c$ 值，KdV 方程有一个单孤波解。KdV 方程有无穷多个单孤波解。

再讨论 (21.29) 式中取正号的情形。重复上面的步骤，又可以得到

$$
f(\xi-c\tau)
=
\frac c2
\operatorname{csch}^2
\left\{
\frac{\sqrt c}{2}
[
(\xi-\xi_0)-c(\tau-\tau_0)
]
\right\}.
\tag{21.31}
$$

应该说，这只是一个形式解，它在 $(\xi-\xi_0)-c(\tau-\tau_0)=0$ 处具有奇异性。

KdV 方程还可以有双孤波解

$$
v(\xi,\tau)
=
-2
\frac{
k_1^2E_1
+
k_2^2E_2
+
2(k_2-k_1)^2E_1E_2
+
A(k_2^2E_1+k_1^2E_2)E_1E_2
}{
(1+E_1+E_2+AE_1E_2)^2
},
\tag{21.32}
$$

其中

$$
E_1
=
\exp\{k_1\xi-k_1^3\tau+\alpha_1\},
\qquad
E_2
=
\exp\{k_2\xi-k_2^3\tau+\alpha_2\},
\qquad
A
=
\left(
\frac{k_2-k_1}{k_2+k_1}
\right)^2.
$$

求出这种解的方法和步骤，本书不再介绍。对于非线性偏微分方程的求解问题，读者可以阅读陆振球教授所著的《经典和现代数学物理方法》（上海科学技术出版社，2004 年）。在姚端正和梁家宝教授所著的《数学物理方法》（武汉大学出版社，1997 年第二版）中也介绍了求解非线性方程的某些初等解法。

## §21.4 结束语

在结束本书的时候，需要再次指出，本书关于数学物理方程的介绍，实际上只限于线性偏微分方程的部分内容，因此，远远不是完整的。大的说来，数学物理方程，既包括偏微分方程，还包括积分方程，后者在本书中完全没有触及。就偏微分方程来说，既有线性偏微分方程，又有非线性偏微分方程，对后者只是在 §21.3 中做了一点序言式的介绍。而且，就线性偏微分方程来说，也不单有定解问题的求解问题，即由“源”求“场”的问题，还有它的逆问题，即由“场”求“源”的问题。后者同样具有理论和实用上的重要性。遥感、地质勘探、机器故障诊断、人体 CT 成像等等这些问题，无一不是由场求源的逆问题。这类问题，本书中也丝毫未曾涉及。有兴趣的读者，也可参阅上节末提到的《经典和现代数学物理方法》一书。

## 习题

1. 讨论下列方程的类型，并将它们化为标准形式：

(1)

$$
\frac{\partial^2u}{\partial x^2}
+
y\frac{\partial^2u}{\partial y^2}
+
\frac12\frac{\partial u}{\partial y}
=
0;
$$

(2)

$$
(1+x^2)
\frac{\partial^2u}{\partial x^2}
+
(1+y^2)
\frac{\partial^2u}{\partial y^2}
+
x\frac{\partial u}{\partial x}
+
y\frac{\partial u}{\partial y}
=
0;
$$

(3)

$$
\tan^2x
\frac{\partial^2u}{\partial x^2}
-
2y\tan x
\frac{\partial^2u}{\partial x\partial y}
+
y^2
\frac{\partial^2u}{\partial y^2}
+
y^2
\frac{\partial u}{\partial y}
=
0;
$$

(4)

$$
\frac{\partial^2u}{\partial x^2}
-
2\sin x
\frac{\partial^2u}{\partial x\partial y}
-
\cos^2x
\frac{\partial^2u}{\partial y^2}
-
\cos x
\frac{\partial u}{\partial y}
=
0.
$$

2. 有些方程，对未知函数作适当的变换后，可以消去一阶偏导数项。

(1) 证明：在变换

$$
u(x,y)
=
e^{-(ax+by)}v(x,y)
$$

下，方程

$$
\nabla^2u
+
2a\frac{\partial u}{\partial x}
+
2b\frac{\partial u}{\partial y}
=
0
$$

（其中 $a,b$ 为常数）化为 Helmholtz 方程

$$
\nabla^2v
-
(a^2+b^2)v
=
0.
$$

(2) 寻求适当的变换，使方程

$$
\frac{\partial^2u}{\partial x^2}
-
\frac{\partial^2u}{\partial y^2}
+
2a\frac{\partial u}{\partial x}
+
2b\frac{\partial u}{\partial y}
=
0
$$

在变换后不再含有一阶偏导数项。

(3) 设有方程

$$
a\frac{\partial^2u}{\partial x^2}
+
2b\frac{\partial^2u}{\partial x\partial y}
+
c\frac{\partial^2u}{\partial y^2}
+
d\frac{\partial u}{\partial x}
+
e\frac{\partial u}{\partial y}
+
fu
=
\frac{\partial u}{\partial t},
$$

其中 $a,b,c,d,e,f$ 为常数，且 $b^2-ac\ne0$。试证明：在变换

$$
u(x,y,t)
=
e^{\alpha x+\beta y+\gamma t}
v(x,y,t)
$$

下，可使 $v(x,y,t)$ 满足方程

$$
a\frac{\partial^2v}{\partial x^2}
+
2b\frac{\partial^2v}{\partial x\partial y}
+
c\frac{\partial^2v}{\partial y^2}
=
\frac{\partial v}{\partial t}.
$$

3. 求解弦振动方程的 Goursat 问题：

$$
\frac{\partial^2u}{\partial t^2}
-
a^2
\frac{\partial^2u}{\partial x^2}
=
0,
$$

$$
\left.u\right|_{x-at=0}
=
\phi(x),
\qquad
\left.u\right|_{x+at=0}
=
\psi(x),
$$

其中 $\phi(x),\psi(x)$ 满足 $\phi(0)=\psi(0)$。

4. 在波动方程

$$
\frac{\partial^2u}{\partial t^2}
-
a^2
\frac{\partial^2u}{\partial x^2}
=
0
$$

中用 $iy$ 代替 $at$，就能得到 Laplace 方程的“初值”问题

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
0,
$$

$$
\left.u\right|_{y=0}
=
\phi(x),
\qquad
\left.
\frac{\partial u}{\partial y}
\right|_{y=0}
=
\psi(x)
$$

的形式解为

$$
u
=
\frac12
[
\phi(x+iy)
+
\phi(x-iy)
]
+
\frac1{2i}
\int_{x-iy}^{x+iy}
\psi(\xi)\,d\xi.
$$

(1) 令 $\phi(x)=x,\ \psi(x)=e^{-x}$，则可得

$$
u(x,y)
=
x+e^{-x}\sin y.
$$

验证这个表达式处处满足 Laplace 方程，也满足 $y=0$ 时的“初始”条件。

(2) 如果

$$
\phi(x)
=
\frac1{1+x^2},
\qquad
\psi(x)=0,
$$

则形式解变为

$$
u(x,y)
=
\frac{
1+x^2-y^2
}{
(1+x^2-y^2)^2
+
4x^2y^2
}.
$$

试证明：这个函数在 $(0,\pm1)$ 点不连续，因此，至少在这些点上，并不满足 Laplace 方程。这说明在一般情况下，Laplace 方程的“初值”问题无解。
