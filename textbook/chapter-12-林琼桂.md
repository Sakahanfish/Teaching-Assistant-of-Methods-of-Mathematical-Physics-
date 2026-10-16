# 第十二章　柱函数初步*

当物理问题具有圆柱形边界面时，比如电动力学中的圆柱形波导和同轴传输线，就应该采用柱坐标系求解。在柱坐标系中对常见的数理方程分离变量，就会遇到 Bessel 方程：

$$
y''+\frac{1}{x}y'+\left(1-\frac{\nu^2}{x^2}\right)y=0.
\tag{1}
$$

物理上遇到的通常是 $\nu=m\in\mathbb N$ 的情况，数学上 $\nu$ 可以是复数，但可设 $\operatorname{Re}\nu\geqslant0$ 而不失一般性。

我们已经在第十章用级数解法求解了 Bessel 方程，其线性独立解是

$$
y(x)=\{J_\nu(x),N_\nu(x)\}.
\tag{2}
$$

其中

$$
J_\nu(x)=\sum_{k=0}^{\infty}\frac{(-)^k}{k!\Gamma(k+\nu+1)}
\left(\frac{x}{2}\right)^{2k+\nu},
\qquad
N_\nu(x)=\frac{J_\nu(x)\cos\nu\pi-J_{-\nu}(x)}{\sin\nu\pi}.
\tag{3}
$$

是 $\nu$ 阶 Bessel 函数和 $\nu$ 阶 Neumann 函数。有些书上将 Neumann 函数记作 $Y_\nu(x)$。当 $\nu=m\in\mathbb N$，$N_\nu(x)$ 为 $0/0$ 型，可以通过取极限并应用 L' Hospital 法则确定其形式，这与用级数解法直接求第二解（并适当选取任意常数）所得的 $N_m(x)$ 相同，即

$$
\lim_{\nu\to m}N_\nu(x)=N_m(x).
$$

Bessel 函数也称为第一类柱函数，Neumann 函数也称为第二类柱函数，下列 Hankel 函数也称为第三类柱函数：

$$
H_\nu^{(1)}(x)=J_\nu(x)+iN_\nu(x),
\qquad
H_\nu^{(2)}(x)=J_\nu(x)-iN_\nu(x).
\tag{4}
$$

$H_\nu^{(1)}(x)$ 和 $H_\nu^{(2)}(x)$ 也分别称为第一种和第二种 Hankel 函数。Bessel 方程的线性独立解可以取为 $\{J_\nu(x),N_\nu(x),H_\nu^{(1)}(x),H_\nu^{(2)}(x)\}$ 中的任意两个。

柱函数在物理和工程中有广泛的应用。本章介绍其基本知识和初步应用。人们研究 Bessel 函数和相关函数已经超过两百年，积累了大量的知识。本章虽然内容不少，但都是比较初步的知识。

容易知道，函数 $\{J_{-\nu}(x),N_{-\nu}(x),H_{-\nu}^{(1)}(x),H_{-\nu}^{(2)}(x)\}$ 中的任何一个都是方程 (1) 的解，而且它们彼此线性独立，所以方程的解也可以取为其中的任意两个。但是，这一组解和上面那组解之间则不是互相独立的。事实上，容易证明

$$
J_{-\nu}(x)=J_\nu(x)\cos\nu\pi-N_\nu(x)\sin\nu\pi,
$$

所以它是 $J_\nu(x)$ 和 $N_\nu(x)$ 的线性组合。特别地，

$$
J_{-m}(x)=(-)^mJ_m(x),\qquad
J_{-(l+1/2)}(x)=(-)^{l+1}N_{l+1/2}(x),
\qquad m,l\in\mathbb N.
$$

类似地，

$$
N_{-\nu}(x)=J_\nu(x)\sin\nu\pi+N_\nu(x)\cos\nu\pi,
$$

也是 $J_\nu(x)$ 和 $N_\nu(x)$ 的线性组合。特别地，

$$
N_{-(l+1/2)}(x)=(-)^lJ_{l+1/2}(x),
\qquad
N_{-m}(x)=(-)^mN_m(x).
$$

最后，易证

$$
H_{-\nu}^{(1)}(x)=e^{i\nu\pi}H_\nu^{(1)}(x),
\qquad
H_{-\nu}^{(2)}(x)=e^{-i\nu\pi}H_\nu^{(2)}(x),
$$

可见两组 Hankel 函数之间的线性相关性尤其明显。在实践上，当 $\nu\notin\mathbb N$，可能采用 $J_{-\nu}(x)$ 作为方程的解，其他 $-\nu$ 阶的函数一般都不采用。

如果 $\nu\notin\mathbb Z$，则 $J_\nu(x)$ 是多值函数，通常将宗量的辐角限制在主值范围内，使得函数值有确切的定义。但可以扩展函数的定义使得当宗量辐角超出主值范围时函数值仍然有确切的定义。实际上，在 $J_\nu(x)$ 的级数形式中，只有因子 $x^\nu$ 可能是多值函数，所以只需定义好这个因子即可。设 $x=re^{i\theta}$，其中 $\theta$ 不一定是主值，可以是主值加上 $2\pi$ 的整数倍，按第一章，定义 $x^\nu=e^{\nu\ln x}=r^\nu e^{i\nu\theta}$，则 $J_\nu(x)$ 对宗量的任何辐角就都有了确切的定义。于是，其他柱函数也就对于宗量的任何辐角都有了确切的定义。根据上述定义易得 $(e^{i\alpha}x)^\nu=e^{i\nu\alpha}x^\nu$，从而有

$$
J_\nu(e^{\pm i\pi}x)=e^{\pm i\nu\pi}J_\nu(x),
$$

其中 $x$ 不必是正实数，其辐角也不必是主值。

## §1　柱函数的基本性质

### §1.1　原点附近的行为

根据定义，并利用 $\Gamma(\nu)\Gamma(1-\nu)=\pi/\sin\nu\pi$（证明从略），容易得到，在 $x=0$ 附近，

$$
J_0(x)=1-\left(\frac{x}{2}\right)^2\to1,
\qquad
N_0(x)=\frac{2}{\pi}\ln\frac{x}{2}\to\infty.
\tag{5a}
$$

$$
J_\nu(x)=\frac{1}{\Gamma(1+\nu)}
\left(\frac{x}{2}\right)^\nu\to0,
\qquad
J_{-\nu}(x)=\frac{1}{\Gamma(1-\nu)}
\left(\frac{x}{2}\right)^{-\nu}\to\infty,
\qquad \operatorname{Re}\nu>0.
\tag{5b}
$$

$$
N_\nu(x)=-\frac{\Gamma(\nu)}{\pi}
\left(\frac{x}{2}\right)^{-\nu}\to\infty,
\qquad \operatorname{Re}\nu>0.
\tag{5c}
$$

至于 $H_\nu^{(1)}(x)$ 和 $H_\nu^{(2)}(x)$，在 $\operatorname{Re}\nu\geqslant0$ 的条件下，其行为与 $N_\nu(x)$ 一样，只是多了一个系数 $\pm i$，故不必单独列出。初学的时候，只要记住以上各式的幂次（或为对数函数）即可，系数不是很重要。甚至只是记住哪个有限，哪个有奇性，也已经足以应付不少应用的需要，比如物理问题中解的取舍。

上述性质对 $x\in\mathbb C$ 亦成立。

### §1.2　递推关系

以 $Z_\nu(x)$ 表示 $\{J_\nu(x),N_\nu(x),H_\nu^{(1)}(x),H_\nu^{(2)}(x)\}$ 中的任一函数，可以证明以下递推关系：

$$
\frac{d}{dx}[x^\nu Z_\nu(x)]=x^\nu Z_{\nu-1}(x),
\qquad
\frac{d}{dx}[x^{-\nu}Z_\nu(x)]=-x^{-\nu}Z_{\nu+1}(x).
\tag{6}
$$

对于 $Z_\nu(x)=J_\nu(x)$，将 $J_\nu(x)$ 的级数表达式代入左边，求导，作一些简单的计算和整理，即可证明。利用所得结果和 $N_\nu(x)$ 的一般定义，即可证明两式对 $Z_\nu(x)=N_\nu(x)$ 亦成立。至于 $Z_\nu(x)$ 为 $H_\nu^{(1)}(x)$ 或 $H_\nu^{(2)}(x)$ 的情况，则可由前述结果直接导出。令 $\nu=0$，可得一个非常简单的特例：

$$
J_0'(x)=-J_1(x).
\tag{7}
$$

利用上述递推关系，可以导出一些简单的积分公式。

求出式 (6) 中的导数，两式相加减，可以得出递推关系的另一种形式：

$$
Z_{\nu-1}(x)+Z_{\nu+1}(x)=\frac{2\nu}{x}Z_\nu(x),
\qquad
Z_{\nu-1}(x)-Z_{\nu+1}(x)=2Z_\nu'(x).
\tag{8}
$$

满足递推关系 (6) 或 (8) 的函数称为柱函数。显然，$\{J_\nu(x),N_\nu(x),H_\nu^{(1)}(x),H_\nu^{(2)}(x)\}$ 中的任何一个或者它们的线性组合，只要系数与 $\nu$ 无关，都是柱函数。不难证明，柱函数一定满足 Bessel 方程。事实上，将式 (6) 中的导数展开，然后将第一式中的 $\nu$ 换成 $\nu+1$，再结合第二式，即可证明。反过来，Bessel 方程的解不一定是柱函数。比如 $\nu J_\nu(x)$ 是 Bessel 方程的解，但不是柱函数，由于存在系数 $\nu$，容易验证其不满足递推关系。

### §1.3　渐近形式

对于特殊函数，人们常常希望知道，当自变量很大时，其行为像哪个初等函数，这有助于把握其性质并应用于具体问题。对于柱函数，有如下结果：当 $x\to\infty$ 且 $-\pi<\arg x<\pi$，

$$
H_\nu^{(1)}(x)\sim
\sqrt{\frac{2}{\pi x}}\,
e^{i(x-\nu\pi/2-\pi/4)},
\tag{9}
$$

$$
H_\nu^{(2)}(x)\sim
\sqrt{\frac{2}{\pi x}}\,
e^{-i(x-\nu\pi/2-\pi/4)},
\tag{10}
$$

$$
J_\nu(x)\sim
\sqrt{\frac{2}{\pi x}}
\cos\left(x-\frac{\nu\pi}{2}-\frac{\pi}{4}\right),
\tag{11}
$$

$$
N_\nu(x)\sim
\sqrt{\frac{2}{\pi x}}
\sin\left(x-\frac{\nu\pi}{2}-\frac{\pi}{4}\right).
\tag{12}
$$

证明从略。

**注**　① 由上述渐近形式可以看出，四个函数 $\{J_\nu(x),N_\nu(x),H_\nu^{(1)}(x),H_\nu^{(2)}(x)\}$ 之间的关系类似于 $\{\cos x,\sin x,e^{ix},e^{-ix}\}$ 之间的关系。② 由渐近形式大致可以看出，如果 $\nu\in\mathbb R$，$J_\nu(x)$ 和 $N_\nu(x)$ 在正实轴上都有无穷多零点。③ 这几个函数在正实轴上的图像都是衰减振荡的，振幅随着 $x$ 的增大按 $1/\sqrt{x}$ 衰减，所以和柱面波有关。④ 式 (9) 的成立范围可以推广到 $-\pi<\arg x<2\pi$，而式 (10) 的成立范围可以推广到 $-2\pi<\arg x<\pi$。不过，多数情况下我们只需要用到上半平面的结果，甚至只是 $x\to+\infty$ 的结果。

### §1.4　零点

首先指出，柱函数 $Z_\nu(x)$ 的任何零点 $x_0$ 都是简单零点，只有 $x_0=0$ 可能例外。事实上，如果 $x_0\neq0$ 是二阶零点，则 $Z_\nu(x_0)=Z_\nu'(x_0)=0$，利用 Bessel 方程易得 $Z_\nu''(x_0)=0$（如果 $x_0=0$ 则不成立），再对方程求导并代入 $x=x_0$，易得 $x_0$ 处各阶导数均为零，那么利用 Taylor 展开式，可得 $x_0$ 的邻域 $N(x_0,|x_0|)$ 内 $Z_\nu(x)\equiv0$，这显然是错误的，故前述结论成立。

其次，易知 $J_\nu(e^{i\pi}x)=e^{i\nu\pi}J_\nu(x)$，其中 $x>0$（实际上该式对任何 $x\in\mathbb C$ 成立，$x$ 的辐角也不必是主值，但初学的读者只需理解其对 $x>0$ 成立即可）。如果 $x_0>0$ 是 $J_\nu(x)$ 的零点，则 $e^{i\pi}x_0$ 显然也是其零点。所以，对于 $J_\nu(x)$ 在实轴上的零点，只需要关注正零点。

最后，在讨论零点时，我们只考虑 $\nu\in\mathbb R$ 的情况。

下面给出关于柱函数零点的几个结论，并作一些相关的讨论或说明。

(1) $J_\nu(x)$ 有无穷多正零点。

这一结论的证明并不深奥，但略为繁琐，从略。但是这个结论是容易理解的。一方面，从渐近形式 (11) 可以大致推测这一结果，前面也曾经指出过。另一方面，我们知道结论对于 $J_{1/2}(x)=\sqrt{2/\pi x}\sin x$ 是显然的，而由 $J_\nu(x)$ 的级数形式可以看出，不同 $\nu$ 值的级数其实具有类似的行为。

这一结论将会保证相关物理问题在第一类边界条件下具有无穷多本征值。

由上述结论和 Rolle 定理易知，$J_\nu'(x)$ 有无穷多正零点。这保证了相关物理问题在第二类边界条件下具有无穷多本征值。容易证明这些正零点也都是简单零点。

(2) $J_\nu(x)$ 和 $J_{\nu+1}(x)$ 的正零点是两两相间的。如果 $\nu>-1$，则 $J_\nu(x)$ 的最小正零点更靠近原点。

考察递推关系式 (6) 的第一式，将 $\nu$ 换成 $\nu+1$，写成

$$
\frac{d}{dx}[x^{\nu+1}J_{\nu+1}(x)]=x^{\nu+1}J_\nu(x).
\tag{13}
$$

由该式和 Rolle 定理，易知在 $J_{\nu+1}(x)$ 的两个相邻正零点之间，至少有 $J_\nu(x)$ 的一个正零点；而由式 (6) 的第二式和 Rolle 定理可知，在 $J_\nu(x)$ 的两个相邻正零点之间，至少有 $J_{\nu+1}(x)$ 的一个正零点。所以，$J_\nu(x)$ 和 $J_{\nu+1}(x)$ 的正零点是两两相间的。如果 $\nu>-1$，则 $x^{\nu+1}J_{\nu+1}(x)|_{x=0}=0$，于是，由式 (13) 和 Rolle 定理可知，在 $x=0$ 和 $J_{\nu+1}(x)$ 的最小正零点之间，必有 $J_\nu(x)$ 的一个正零点（也只能有一个）。所以，$J_\nu(x)$ 的最小正零点更靠近原点。

**思考 1**　如果 $\nu>-1$，在 $x=0$ 和 $J_{\nu+1}(x)$ 的最小正零点之间，是否可以存在 $J_\nu(x)$ 的多个正零点？

**思考 2**　如果 $\nu<-1$，则 $J_{\nu+1}(x)$ 的最小正零点比 $J_\nu(x)$ 的更靠近原点。如何证明？

(3) $AJ_\nu(x)+BxJ_\nu'(x)$ 有无穷多正零点，其中 $A、B\in\mathbb R$ 且不全为零。

证明从略。如果 $A、B$ 之一为零，则结论已经在前面给出。所以我们主要关注 $A、B$ 均不为零的情况。此时上述结论保证了相关物理问题在第三类边界条件下存在无穷多本征值。如果 $\nu>-1$，可以证明这些正零点都是简单零点，而且当 $AD-BC\neq0$，$AJ_\nu(x)+BxJ_\nu'(x)$ 与 $CJ_\nu(x)+DxJ_\nu'(x)$ 的正零点是两两相间的。初学者只需要了解基本结论即可。

(4) $N_\nu(x)$ 有无穷多正零点，而且与 $J_\nu(x)$ 的正零点两两相间。

证明从略。可以更一般地证明，$AJ_\nu(x)+BN_\nu(x)$ 有无穷多正零点，其中 $A、B\in\mathbb R$ 且不全为零，而且当 $AD-BC\neq0$，它与 $CJ_\nu(x)+DN_\nu(x)$ 的正零点是两两相间的。当然这些正零点都是简单零点，因为 $AJ_\nu(x)+BN_\nu(x)$ 是一个柱函数。

(5) 当 $\nu>-1$，$J_\nu(x)$ 的零点都是实数。

此时 $J_\nu(x)$ 的零点都在实轴上，类似于 $\cos x$ 和 $\sin x$。不过，当 $\nu<-1$，$J_\nu(x)$ 存在一些复零点。

## §2　整数阶柱函数

当 $\nu=m\in\mathbb N$，Bessel 函数的级数可以简化为

$$
J_m(x)=\sum_{k=0}^{\infty}
\frac{(-)^k}{k!(k+m)!}
\left(\frac{x}{2}\right)^{2k+m}.
\tag{14}
$$

Neumann 函数 $N_m(x)$ 的形式比较复杂，已经在第十章给出。前已指出，$J_{-m}(x)=(-)^mJ_m(x)$，$N_{-m}(x)=(-)^mN_m(x)$。所以，负整数阶的情况不需要过多关注。

§1 的所有结论对整数阶柱函数成立，但整数阶柱函数最常用，所以本节加以详细讨论。虽然有些结果只是 §1 相关结果的简单推论，也会重新提到。

### §2.1　对称性

由式 (14) 容易证明

$$
J_m(-x)=(-)^mJ_m(x).
\tag{15}
$$

即 $J_m(x)$ 的奇偶性与 $m$ 一致。

### §2.2　零点

根据 §1 的结论，易得以下结果。

(1) $J_m(x)$ 有无穷多实数零点，且正负成对。$J_m'(x)$ 亦然。

注意 $J_m(x)$ 是单值函数，所以负实轴上与 $x>0$ 对称的点可以直接写成 $-x$，而不必写成 $e^{i\pi}x$。

$J_m(x)$ 的正零点从小到大排列，记作 $\{x_{mn}\}_{n=1}^\infty$，满足

$$
0<x_{m1}<x_{m2}<\cdots.
\tag{16}
$$

$J_m'(x)$ 的正零点从小到大排列，记作 $\{\tilde{x}_{mn}\}_{n=1}^\infty$，满足

$$
0<\tilde{x}_{m1}<\tilde{x}_{m2}<\cdots.
\tag{17}
$$

(2) $J_m(x)$ 和 $J_{m+1}(x)$ 的正零点是两两相间的。$J_m(x)$ 的最小正零点比 $J_{m+1}(x)$ 的更靠近原点，即 $x_{m1}<x_{m+1,1}$。

### §2.3　柱函数的图像

$J_m(x)$ 和 $N_m(x)$ 的图像都是衰减振荡的，但是 $x\to0$ 时，它们具有不同的特征，此时 $N_m(x)\to-\infty$，而 $J_0(0)=1$，对于 $m\in\mathbb N^+$，则 $J_m(0)=0$。

下图分别给出 $J_0(x)$、$J_1(x)$ 的图像和 $N_0(x)$、$N_1(x)$ 的图像。[^1] $J_m(x)$（$m\geqslant2$）的图像与 $J_1(x)$ 类似，$N_m(x)$（$m\geqslant2$）的图像与 $N_1(x)$ 类似。

### §2.4　母函数

我们在第四章曾经证明

$$
\exp\left[\frac{z}{2}\left(t-\frac{1}{t}\right)\right]
=\sum_{m=-\infty}^{\infty}J_m(z)t^m,
\qquad 0<|t|<+\infty.
\tag{18}
$$

这个展开式的意义是：左边的函数在环域 $0<|t|<+\infty$ 内是 $t$ 的解析函数，可以展开为 $t$ 的 Laurent 级数，展开系数与参数 $z$ 有关，计算表明其正好是 Bessel 函数。初学者如果不能掌握计算的细节，则只需理解上述意义即可。

母函数展开式很有用，后面会看到其部分应用。这里举一个简单的例子。将上式两边对参数 $z$ 求导，对左边求导所得结果再次利用展开式并整理可得

$$
\frac12\sum_{m=-\infty}^{\infty}
[J_{m-1}(z)-J_{m+1}(z)]t^m
=
\sum_{m=-\infty}^{\infty}J_m'(z)t^m.
\tag{19}
$$

比较两边 $t^m$ 的系数，就得到 Bessel 函数的递推关系之一，即式 (8) 第二式的特殊情况。

### §2.5　加法公式

利用指数函数的性质，有

$$
\exp\left[\frac{x+y}{2}\left(t-\frac1t\right)\right]
=
\exp\left[\frac{x}{2}\left(t-\frac1t\right)\right]
\exp\left[\frac{y}{2}\left(t-\frac1t\right)\right].
\tag{20}
$$

两边均用母函数公式展开，右边为 $\sum_{k,l=-\infty}^{\infty}J_k(x)J_l(y)t^{k+l}$，将 $k+l=m$ 的各项（有无穷多个）合并到一起，整理可得

$$
\sum_{m=-\infty}^{\infty}J_m(x+y)t^m
=
\sum_{m=-\infty}^{\infty}t^m
\sum_{k=-\infty}^{\infty}J_k(x)J_{m-k}(y).
\tag{21}
$$

比较两边 $t^m$ 的系数，就得到 Bessel 函数的加法公式

$$
J_m(x+y)=\sum_{k=-\infty}^{\infty}J_k(x)J_{m-k}(y).
\tag{22}
$$

柱函数还有更复杂的加法公式，我们就不介绍了。

### §2.6　积分表示

在母函数公式 (18) 中令 $t=e^{i\phi}$，其中 $\phi\in\mathbb R$，可得

$$
\exp(ix\sin\phi)=\sum_{m=-\infty}^{\infty}J_m(x)e^{im\phi}.
\tag{23}
$$

这是 $\phi$ 的复形式的 Fourier 级数，只不过左边的函数含有一个参数 $x$，因而展开系数是 $x$ 的函数。由系数公式，可得

$$
J_m(x)=\frac{1}{2\pi}\int_{-\pi}^{\pi}
e^{ix\sin\phi-im\phi}\,d\phi
=
\frac1\pi\int_0^\pi\cos(x\sin\phi-m\phi)\,d\phi.
\tag{24}
$$

这就是 Bessel 函数的积分表示。它把 Bessel 函数表示成一个含参数积分，被积函数是初等函数，这是特殊函数积分表示的一般特点。

### §2.7　含 Bessel 函数的积分

含 Bessel 函数的积分五花八门，其计算涉及各种各样的技巧。下面仅举一例。试计算

$$
I=\int_0^\infty e^{-ax}J_0(bx)\,dx,
\qquad a,b>0.
\tag{25}
$$

第一种方法是将积分表示

$$
J_0(x)=\frac{1}{2\pi}\int_{-\pi}^{\pi}e^{ix\sin\phi}\,d\phi
\tag{26}
$$

代入上式并交换积分次序，可得

$$
I=\frac{1}{2\pi}\int_{-\pi}^{\pi}d\phi
\int_0^\infty dx\,e^{-(a-ib\sin\phi)x}.
\tag{27}
$$

完成对 $x$ 的积分，得到

$$
I=\frac{1}{2\pi}\int_{-\pi}^{\pi}\frac{d\phi}{a-ib\sin\phi}.
\tag{28}
$$

这个积分可以用围线积分和留数定理的方法来计算，是第五章讲过的标准类型之一，读者可以自行计算，结果是

$$
\int_0^\infty e^{-ax}J_0(bx)\,dx
=
\frac{1}{\sqrt{a^2+b^2}},
\qquad a,b>0.
\tag{29}
$$

如果 $b<0$，由于 $J_0(x)$ 是偶函数，故在被积函数中 $J_0(bx)$ 可代以 $J_0(|b|x)$，从而不必重新计算，此时上式结果仍然成立。如果 $b=0$，结果也显然成立。所以 $b>0$ 的条件实际上可以去掉。如果将积分中的 $J_0(bx)$ 改为 $J_m(bx)$（其中 $m\in\mathbb N^+$），仍可用以上方法计算，只是较上例略为繁琐而已。具体计算留作练习。

第二种方法是先设 $b<a$，将 $J_0(bx)$ 的级数形式代入，逐项积分，可得

$$
I
=
\sum_{k=0}^{\infty}\frac{(-)^k}{(k!)^2}
\left(\frac b2\right)^{2k}
\int_0^\infty x^{2k}e^{-ax}\,dx
=
\frac1a\sum_{k=0}^{\infty}
\frac{(-)^k(2k)!}{2^{2k}(k!)^2}
\left(\frac ba\right)^{2k}.
\tag{30}
$$

不难验证，上式中的系数等于 $C_{-1/2}^k$（但是认识到这一点则不容易），由此立得式 (29) 的结果。最后，可以证明式 (29) 左边在 $\operatorname{Re}a>0$ 是 $a$ 的解析函数（证明从略），右边亦然，根据解析开拓理论，条件 $b<a$ 的限制就可以去掉。这一方法的主要困难在于将逐项积分得到的级数求和得出一个初等函数。如果将积分中的 $J_0(bx)$ 改为 $J_m(bx)$（其中 $m\in\mathbb N^+$），则最后一步的难度会显著增加。

以上两种方法是计算涉及 Bessel 函数的积分时常用的方法。


### §2.8　本征值问题

对 Laplace 方程或 Helmholtz 方程在柱坐标系中分离变量，会遇到下述方程

$$
y''+\frac1x y'+\left(\lambda-\frac{m^2}{x^2}\right)y=0,
\qquad m\in\mathbb N.
\tag{31}
$$

将它改写成 Sturm–Liouville 形式，结合边界条件，典型的本征值问题是：

$$
\frac{d}{dx}\left(x\frac{dy}{dx}\right)
-\frac{m^2}{x}y+\lambda xy=0,
\qquad (0<x<a)
\tag{32a}
$$

$$
y(0)=0\quad(m=1,2,\cdots)
\quad\text{或}\quad
|y(0)|<\infty\quad(m=0),
\qquad
\alpha y'(a)+\beta y(a)=0,
\tag{32b}
$$

其中 $\alpha、\beta\geqslant0$，但 $\alpha^2+\beta^2\neq0$；$x=0$ 处的边界条件是自然边界条件。与 Sturm–Liouville 方程比较可知 $k(x)=x\geqslant0$，$q(x)=m^2/x\geqslant0$，$\rho(x)=x\geqslant0$。根据 Sturm–Liouville 本征值问题的一般结论可知，$\lambda\geqslant0$，故可令 $\lambda=k^2$，其中 $k\geqslant0$。代入方程 (31)，其解为

$$
y(x)=\{J_m(kx),N_m(kx)\}.
\tag{33}
$$

考虑到 $x=0$ 处的自然边界条件，$N_m(kx)$ 显然不符合要求，应该舍弃。故

$$
y(x)=J_m(kx).
\tag{34}
$$

由 $x=a$ 处的边界条件，可得确定本征值的超越方程为

$$
\alpha kJ_m'(ka)+\beta J_m(ka)=0.
\tag{35}
$$

由上节关于零点的结论，上式有无穷多正根。由于左边作为 $k$ 的函数有确定的奇偶性，故实根是正负成对出现的，但是负根既不能给出新的本征值，也不能给出独立的本征函数，所以是多余的。事实上，我们引入 $k$ 替代 $\lambda$ 的时候，已经限制 $k\geqslant0$。对此一问题，今后将不再重复说明。将正根从小到大排列，记作 $\{k_n\}_{n=1}^{\infty}$，满足

$$
0<k_1<k_2<\cdots.
\tag{36}
$$

相应的本征值和本征函数为

$$
\lambda_n=k_n^2,
\qquad
y_n(x)=J_m(k_nx),
\qquad
n=1,2,\cdots.
\tag{37}
$$

**注**　① 在实际问题中，常常需要把对应于不同 $m$ 值的解叠加，此时需要同时处理对应于不同 $m$ 值的一系列而不是一个本征值问题，则方程 (35) 的正根应该记作 $\{k_{mn}\}_{n=1}^{\infty}$，相应的本征值和本征函数记作 $\lambda_{mn}=k_{mn}^2$ 和 $y_{mn}(x)=J_m(k_{mn}x)$。以上讨论中为了记号方便，认为 $m$ 是固定的。② 以上求解中未考虑 $k=0$ 或 $\lambda=0$，容易验证（第十章也讨论过），只当 $\beta=0$（第二类齐次边界条件）且 $m=0$ 时，$\lambda=0$ 才是一个本征值。即对于第二类齐次边界条件，$m=0$ 对应的结果需要作修正，增加一个本征值 $\lambda_{00}=0$ 和相应的本征函数 $y_{00}(x)=1$。注意 $\lambda=0$ 方程需要单独求解，虽然形式上 $J_0(k_{00}x)=1$ 正好也是正确的。

### §2.9　本征函数的正交关系

为了记号方便，我们仍然考虑 $m$ 固定的情况。本征值和本征函数由式 (37) 给出。作为 Sturm–Liouville 本征值问题的特例，本征函数在区间 $[0,a]$ 上满足正交关系：

$$
\int_0^a y_{n'}(x)y_n(x)x\,dx
=
\int_0^a J_m(k_{n'}x)J_m(k_nx)x\,dx
=0,
\qquad n'\neq n.
\tag{38}
$$

注意到权函数 $\rho(x)=x$，所以上式是带权正交。

### §2.10　本征函数的模

对本征函数 $y_n(x)$ 满足的方程 (32a) 两边同乘以 $xy_n'(x)$，可得

$$
xy_n'\frac{d}{dx}(xy_n')
-m^2y_ny_n'
+k_n^2x^2y_ny_n'=0.
\tag{39}
$$

从 $0$ 到 $a$ 积分，第一项和第二项均可直接得出原函数，对第三项作一次分部积分，注意到 $m^2y_n^2(0)=m^2J_m^2(0)=0$ 对所有 $m\in\mathbb N$ 成立，可得

$$
\|y_n(x)\|^2
\equiv
\int_0^a y_n^2(x)x\,dx
=
\frac{a^2}{2k_n^2}y_n'^2(a)
+\frac12\left(a^2-\frac{m^2}{k_n^2}\right)y_n^2(a),
\tag{40}
$$

或

$$
\|J_m(k_nx)\|^2
\equiv
\int_0^a J_m^2(k_nx)x\,dx
=
\frac{a^2}{2}J_m'^2(k_na)
+\frac12\left(a^2-\frac{m^2}{k_n^2}\right)J_m^2(k_na).
\tag{41}
$$

注意 $y_n'(a)$ 中的导数是对自变量 $x$ 求导，而 $J_m'(k_na)$ 中的导数是对宗量 $k_nx$ 求导（最后都代入 $x=a$），所以两式中对应项的系数不一样。

对于第一类齐次边界条件，式 (35) 简化为 $J_m(ka)=0$，故

$$
k_n=\frac{x_{mn}}{a},
\qquad n=1,2,\cdots,
\tag{42}
$$

而

$$
\|J_m(k_nx)\|^2
=
\frac{a^2}{2}J_m'^2(x_{mn})
=
\frac{a^2}{2}J_{m\pm1}^2(x_{mn}),
\qquad n=1,2,\cdots,
\tag{43}
$$

其中第二步用了递推关系和 $J_m(x_{mn})=0$ 的事实。

对于第二类齐次边界条件，式 (35) 简化为 $J_m'(ka)=0$，故

$$
k_n=\frac{\tilde{x}_{mn}}{a},
\qquad n=1,2,\cdots,
\tag{44}
$$

而

$$
\|J_m(k_nx)\|^2
=
\frac{a^2}{2}
\left(1-\frac{m^2}{\tilde{x}_{mn}^2}\right)
J_m^2(\tilde{x}_{mn}),
\qquad n=1,2,\cdots.
\tag{45}
$$

**注**　① 严格来说，以上记号 $k_n$ 应该记作 $k_{mn}$，$y_n(x)=J_m(k_nx)$ 应该记作 $y_{mn}(x)=J_m(k_{mn}x)$。② 对于第二类齐次边界条件且 $m=0$ 的情况，需要增加一个本征值 $\lambda_{00}=0$ 和相应的本征函数 $y_{00}(x)=1$。易得 $\|y_{00}(x)\|^2=a^2/2$。以下不再重复强调。

### §2.11　广义 Fourier 级数

作为 Sturm–Liouville 本征值问题的特例，$\{J_m(k_nx)\}_{n=1}^{\infty}$ 在区间 $[0,a]$ 上是完备的，也就是说，区间 $[0,a]$ 上的任一函数 $f(x)$，如果满足与本征值问题同样的边界条件且足够光滑，必可用 $\{J_m(k_nx)\}_{n=1}^{\infty}$ 展开为广义 Fourier 级数

$$
f(x)=\sum_{n=0}^{\infty}f_nJ_m(k_nx),
\qquad 0\leqslant x\leqslant a,
\tag{46}
$$

利用正交关系，容易求得展开系数为

$$
f_n=
\frac{1}{\|J_m(k_nx)\|^2}
\int_0^a f(x)J_m(k_nx)x\,dx,
\qquad n=1,2,\cdots.
\tag{47}
$$

### §2.12　应用

**例 1**　无穷长均匀导热圆柱体，半径为 $\rho_0$，侧面保持温度为零度，柱内的初始温度分布为已知函数 $f(\rho,\phi)$，求柱内各处的温度变化。

采用柱坐标系，以圆柱的对称轴为 $z$ 轴，任取轴上一点为原点。柱内的温度分布 $u(\boldsymbol r,t)$ 满足定解问题

$$
\frac{\partial u}{\partial t}-a^2\nabla^2u=0,
\qquad \rho<\rho_0,\ t>0,
\tag{48a}
$$

$$
u|_{\rho=\rho_0}=0,
\tag{48b}
$$

$$
u|_{t=0}=f(\rho,\phi).
\tag{48c}
$$

分离变量，令 $u(\boldsymbol r,t)=v(\boldsymbol r)T(t)$，代入上式中的方程和边界条件，可得常微分方程

$$
T'+\lambda a^2T=0
\tag{49}
$$

和本征值问题

$$
\nabla^2v+\lambda v=0,
\qquad \rho<\rho_0,
\tag{50a}
$$

$$
v|_{\rho=\rho_0}=0.
\tag{50b}
$$

由于圆柱无穷长，初始条件和边界条件均与 $z$ 无关，所以 $u$ 和 $v$ 也与 $z$ 无关。进一步分离变量，令 $v(\boldsymbol r)=R(\rho)\Phi(\phi)$，代入式 (50)，考虑到 $\Phi(\phi)$ 满足周期性边界条件，可得

$$
\Phi(\phi)=\{\cos m\phi,\sin m\phi\},
\qquad m=0,1,2,\cdots.
$$

将求解 $\Phi(\phi)$ 所得的本征值代入 $R(\rho)$ 的方程，并考虑到边界条件，可知 $R(\rho)$ 满足本征值问题

$$
R''+\frac1\rho R'+
\left(\lambda-\frac{m^2}{\rho^2}\right)R=0,
\qquad 0<\rho<\rho_0,
$$

$$
R(0)=0\quad(m\neq0)
\quad\text{或}\quad
|R(0)|<\infty\quad(m=0);
\qquad
R(\rho_0)=0.
$$

其中 $\rho=0$ 处的边界条件是自然边界条件，第九章曾经详细讨论过。这已经是前面讨论过的本征值问题，现在是第一类齐次边界条件，所以得到下列结果。对于给定的 $m$，本征值为

$$
\lambda_{mn}=k_{mn}^2,
\qquad
k_{mn}=\frac{x_{mn}}{\rho_0},
\qquad n=1,2,\cdots.
\tag{51}
$$

相应的本征函数为

$$
R_{mn}(\rho)=J_m(k_{mn}\rho),
\qquad n=1,2,\cdots.
$$

将本征值代回式 (49)，易得 $T_{mn}(t)=\exp(-\lambda_{mn}a^2t)$。将各本征解叠加，得到一般解为

$$
\begin{aligned}
u(\rho,\phi,t)
={}&
\sum_{n=1}^{\infty}
A_{0n}J_0(k_{0n}\rho)\exp(-\lambda_{0n}a^2t)\\
&+
\sum_{m=1}^{\infty}\sum_{n=1}^{\infty}
J_m(k_{mn}\rho)
(A_{mn}\cos m\phi+B_{mn}\sin m\phi)
\exp(-\lambda_{mn}a^2t).
\end{aligned}
\tag{52}
$$

代入初始条件，得到

$$
f(\rho,\phi)
=
\sum_{n=1}^{\infty}A_{0n}J_0(k_{0n}\rho)
+
\sum_{m=1}^{\infty}\sum_{n=1}^{\infty}
J_m(k_{mn}\rho)
(A_{mn}\cos m\phi+B_{mn}\sin m\phi).
\tag{53}
$$

这是双重 Fourier 展开式。对 $\phi$ 是普通的 Fourier 展开式，由系数公式，得到

$$
\sum_{n=1}^{\infty}A_{0n}J_0(k_{0n}\rho)
=
\frac{1}{2\pi}\int_0^{2\pi}f(\rho,\phi)\,d\phi;
\tag{54a}
$$

$$
\sum_{n=1}^{\infty}A_{mn}J_m(k_{mn}\rho)
=
\frac1\pi\int_0^{2\pi}f(\rho,\phi)\cos m\phi\,d\phi,
\qquad m=1,2,\cdots;
\tag{54b}
$$

$$
\sum_{n=1}^{\infty}B_{mn}J_m(k_{mn}\rho)
=
\frac1\pi\int_0^{2\pi}f(\rho,\phi)\sin m\phi\,d\phi,
\qquad m=1,2,\cdots.
\tag{54c}
$$

以上各式中，右边是 $\rho$ 的函数，左边是其广义 Fourier 展开式，根据本节的系数公式，可得

$$
A_{0n}
=
\frac{1}{\pi\rho_0^2J_1^2(x_{0n})}
\int_0^{\rho_0}d\rho\,\rho J_0(k_{0n}\rho)
\int_0^{2\pi}d\phi\,f(\rho,\phi),
\qquad n=1,2,\cdots;
\tag{55a}
$$

$$
A_{mn}
=
\frac{2}{\pi\rho_0^2J_{m+1}^2(x_{mn})}
\int_0^{\rho_0}d\rho\,\rho J_m(k_{mn}\rho)
\int_0^{2\pi}d\phi\,f(\rho,\phi)\cos m\phi,
\qquad m,n=1,2,\cdots;
\tag{55b}
$$

$$
B_{mn}
=
\frac{2}{\pi\rho_0^2J_{m+1}^2(x_{mn})}
\int_0^{\rho_0}d\rho\,\rho J_m(k_{mn}\rho)
\int_0^{2\pi}d\phi\,f(\rho,\phi)\sin m\phi,
\qquad m,n=1,2,\cdots.
\tag{55c}
$$

对计算比较熟悉以后，可以由展开式 (53) 直接写出这些系数。如果 $f(\rho,\phi)$ 给定具体形式，那么原则上可以由上式计算各系数，然后代回式 (52)，就得到最后的结果。容易看出，当 $t\to\infty$，$u(\rho,\phi,t)\to0$，即达到侧面保持的温度，也就是环境的温度，与初始温度分布 $f(\rho,\phi)$ 的具体形式无关。

**特例 1**　如果 $f(\rho,\phi)=f(\rho)$ 与 $\phi$ 无关，那么式 (55b) 和 (55c) 的结果显然为零，此时结果简化为

$$
u(\rho,t)=\sum_{n=1}^{\infty}
A_nJ_0(k_{0n}\rho)\exp(-\lambda_{0n}a^2t),
\tag{56}
$$

其中

$$
A_n=
\frac{2}{\rho_0^2J_1^2(x_{0n})}
\int_0^{\rho_0}f(\rho)J_0(k_{0n}\rho)\rho\,d\rho,
\qquad n=1,2,\cdots.
\tag{57}
$$

**特例 2**　如果在特例 1 的条件下，进一步给定 $f(\rho)=u_0$，其中 $u_0$ 是常数，那么利用 $[xJ_1(x)]'=xJ_0(x)$ 不难求得

$$
A_n=\frac{2u_0}{x_{0n}J_1(x_{0n})},
\qquad n=1,2,\cdots.
\tag{58}
$$

此时结果 (56) 简化为

$$
u(\rho,t)
=
2u_0\sum_{n=1}^{\infty}
\frac{J_0(k_{0n}\rho)}{x_{0n}J_1(x_{0n})}
\exp(-\lambda_{0n}a^2t).
\tag{59}
$$

**习题**　高为 $H$ 的均匀导热圆柱体，半径为 $a$，侧面绝热，上底的温度分布为已知函数 $f_2(\rho)$，下底为 $f_1(\rho)$，求柱内的稳定温度分布。计算特殊情况 $f_1(\rho)=0$，$f_2(\rho)=u_0$ 的结果，其中 $u_0$ 是常数。


## §3　虚宗量 Bessel 函数

在柱坐标系中对常见的数理方程分离变量，除了遇到 Bessel 方程 (1)，还会遇到下列方程（在第九章曾经提到）：

$$
y''+\frac1x y'-\left(1+\frac{\nu^2}{x^2}\right)y=0.
\tag{60}
$$

这称为 $\nu$ 阶虚宗量 Bessel 方程，物理上遇到的通常是 $\nu=m\in\mathbb N$ 的情况。上式可以改写为

$$
y''+\frac1x y'
+\left(i^2-\frac{\nu^2}{x^2}\right)y=0.
\tag{61}
$$

如果令 $t=ix$，那么它就是 $t$ 的 $\nu$ 阶 Bessel 方程，所以它的解就是 $Z_\nu(t)$，或

$$
y(x)=Z_\nu(ix),
\tag{62}
$$

即 $\nu$ 阶的虚宗量柱函数。应该指出，柱函数的宗量本来就可以是复数，所以 $Z_\nu(ix)$ 与 $Z_\nu(x)$ 并无本质区别，那么人们为什么要专门定义虚宗量柱函数呢？这可能有两个原因。第一，对于物理问题，$x$ 通常是实数，那么 $ix$ 就是纯虚数，宗量为纯虚数时柱函数的性质与宗量为实数时颇为不同，类似于双曲函数与三角函数的区别。第二，作为物理问题的解，宗量为纯虚数的柱函数与宗量为实数的柱函数是最常见的两种情况，故两者均受到重视。

为了便于应用，还需要对 $Z_\nu(ix)$ 乘以恰当的系数，使得 $x$ 为实数时 $y(x)$ 也是实数。人们定义了两类虚宗量 Bessel 函数，以下分别讨论。虚宗量 Bessel 函数也称为变型 Bessel 函数。

### §3.1　虚宗量 Bessel 函数

$\nu$ 阶虚宗量 Bessel 函数定义为

$$
I_\nu(x)=e^{-i\nu\pi/2}J_\nu(e^{i\pi/2}x).
\tag{63}
$$

它是方程 (60) 的解。易得其级数形式为

$$
I_\nu(x)=
\sum_{k=0}^{\infty}
\frac{1}{k!\Gamma(k+\nu+1)}
\left(\frac{x}{2}\right)^{2k+\nu}.
\tag{64}
$$

当 $x$ 和 $\nu$ 均为实数时，$I_\nu(x)$ 的函数值为实数，所以上述定义是方便的。虚宗量 Bessel 函数也称为变型 Bessel 函数、第一类虚宗量 Bessel 函数或第一类变型 Bessel 函数。当 $\nu=m\in\mathbb N$，上述级数可以简化为

$$
I_m(x)=
\sum_{k=0}^{\infty}
\frac{1}{k!(k+m)!}
\left(\frac{x}{2}\right)^{2k+m}.
\tag{65}
$$

容易看出，当 $x\geqslant0$，$I_m(x)$ 是单调增函数，这与 $J_m(x)$ 的衰减振荡图像截然不同。粗略地说，$J_m(x)$ 类似于三角函数，而 $I_m(x)$ 类似于指数函数。易得

$$
I_0(0)=1;
\qquad
I_m(0)=0,\quad m\in\mathbb N^+.
\tag{66}
$$

$I_0(x)$ 与 $I_1(x)$ 的图像如图所示。[^2] $I_m(x)$（$m\geqslant2$）的图像与 $I_1(x)$ 类似。

类似于 $J_{-m}(x)=(-)^mJ_m(x)$，不难证明

$$
I_{-m}(x)=I_m(x).
\tag{67}
$$

读者可能会在一些书上看到 $I_\nu(x)$ 的定义包含两个式子，即除了式 (63)，还有另外一个式子，

$$
I_\nu(x)=e^{i3\nu\pi/2}J_\nu(e^{-i3\pi/2}x),
$$

分别用于 $x$ 的不同辐角范围，这是为了保证 $J_\nu(e^{i\pi/2}x)$ 或 $J_\nu(e^{-i3\pi/2}x)$ 的宗量辐角不要超出主值范围。我们已经在 §1 中对 $J_\nu(x)$ 给出了适用于宗量任何辐角值的恰当定义，所以 $I_\nu(x)$ 的定义只需由式 (63) 给出，因为容易证明另一个式子是与其等价的。对于多数应用场景，$x$ 是实数，所以初学者只需了解式 (63) 的定义即可。

根据 $J_\nu(x)$ 的递推关系，可以证明 $I_\nu(x)$ 满足以下递推关系：

$$
\frac{d}{dx}[x^\nu I_\nu(x)]=x^\nu I_{\nu-1}(x),
\qquad
\frac{d}{dx}[x^{-\nu}I_\nu(x)]=x^{-\nu}I_{\nu+1}(x).
\tag{68}
$$

或表示为等价的形式：

$$
I_{\nu-1}(x)+I_{\nu+1}(x)=2I_\nu'(x),
\qquad
I_{\nu-1}(x)-I_{\nu+1}(x)=\frac{2\nu}{x}I_\nu(x).
\tag{69}
$$

一个简单的特例是 $I_0'(x)=I_1(x)$。

### §3.2　虚宗量 Hankel 函数

$\nu$ 阶虚宗量 Hankel 函数定义为

$$
K_\nu(x)=
\frac{\pi}{2\sin\nu\pi}
[I_{-\nu}(x)-I_\nu(x)].
\tag{70}
$$

它也是方程 (60) 的解，而且显然与 $I_\nu(x)$ 线性独立。虚宗量 Hankel 函数也称为第二类虚宗量 Bessel 函数、第二类变型 Bessel 函数或 Macdonald 函数。由定义易知

$$
K_{-\nu}(x)=K_\nu(x).
\tag{71}
$$

考虑方程 (60) 的解，可以假设 $\operatorname{Re}\nu\geqslant0$ 而不失一般性。如果 $\nu\notin\mathbb N$，则易知当 $x\to0$ 时，$I_\nu(x)\to0$，而 $K_\nu(x)\to\infty$，所以，如果求解区间包含 $x=0$，则应当舍弃 $K_\nu(x)$。当 $\nu=m\in\mathbb N$，则有 $I_0(x)\to1$，$I_m(x)\to0$（$m\in\mathbb N^+$），而 $K_m(x)\to\infty$（$m\in\mathbb N$），所以结论也是类似的。由于 $\nu=m\in\mathbb N$ 时，$K_\nu(x)$ 成为 $0/0$ 型，所以它的奇性并不直观，但是我们看到其定义类似于 Neumann 函数的定义，由此可以推测其行为也与 Neumann 函数类似，从而不难理解所述的结论。具体情况可以参考以下的讨论。

当 $\nu=m\in\mathbb N$ 时，$K_\nu(x)$ 成为 $0/0$ 型，我们可以通过取极限 $\nu\to m$ 并用 L' Hospital 法则来计算其级数形式，具体计算过程与 Neumann 函数的情况非常类似，细节从略，结果是

$$
\begin{aligned}
K_m(x)
={}&
(-)^{m+1}I_m(x)\ln\frac{x}{2}
+\frac12\sum_{k=0}^{m-1}
(-)^k\frac{(m-k-1)!}{k!}
\left(\frac{x}{2}\right)^{2k-m}\\
&+
\frac{(-)^m}{2}\sum_{k=0}^{\infty}
\frac{\psi(k+m+1)+\psi(k+1)}
{k!(k+m)!}
\left(\frac{x}{2}\right)^{2k+m}.
\end{aligned}
\tag{72}
$$

当 $m=0$ 时，第一项求和是不存在的，此时结果简化为

$$
K_0(x)=
\sum_{k=0}^{\infty}
\frac{\psi(k+1)}{(k!)^2}
\left(\frac{x}{2}\right)^{2k}
-I_0(x)\ln\frac{x}{2}.
\tag{73}
$$

可以证明

$$
K_\nu(x)
=
\frac{i\pi e^{i\nu\pi/2}}{2}
H_\nu^{(1)}(e^{i\pi/2}x).
\tag{74}
$$

这是它被称为虚宗量 Hankel 函数的原因。

根据 $K_\nu(x)$ 的定义和 $I_\nu(x)$ 的递推关系，或根据式 (74) 和 $Z_\nu(x)$ 的递推关系，可以证明 $K_\nu(x)$ 满足以下递推关系：

$$
\frac{d}{dx}[x^\nu K_\nu(x)]=-x^\nu K_{\nu-1}(x),
\qquad
\frac{d}{dx}[x^{-\nu}K_\nu(x)]=-x^{-\nu}K_{\nu+1}(x).
\tag{75}
$$

或表示为等价的形式：

$$
K_{\nu-1}(x)+K_{\nu+1}(x)=-2K_\nu'(x),
\qquad
K_{\nu-1}(x)-K_{\nu+1}(x)=-\frac{2\nu}{x}K_\nu(x).
\tag{76}
$$

一个简单的特例是 $K_0'(x)=-K_1(x)$。

### §3.3　渐近形式

$I_\nu(x)$ 和 $K_\nu(x)$ 的渐近形式可以由 $J_\nu(x)$ 和 $H_\nu^{(1)}(x)$ 的渐近形式导出，结果如下：

$$
I_\nu(x)\sim\frac{e^x}{\sqrt{2\pi x}},
\qquad
x\to\infty,\quad
-\frac{\pi}{2}<\arg x<\frac{\pi}{2}.
\tag{77}
$$

$$
K_\nu(x)\sim\sqrt{\frac{\pi}{2x}}e^{-x},
\qquad
x\to\infty,\quad
-\pi<\arg x<\pi.
\tag{78}
$$

我们需要用到的情况主要是 $x\to+\infty$ 的情况。$K_\nu(x)$ 的定义使得 $x\to+\infty$ 时不出现发散的部分，从而 $K_\nu(x)\to0$，这使得它成为柱外问题的恰当解。

下图分别给出 $K_0(x)$、$K_1(x)$ 的图像。[^3] $K_m(x)$（$m\geqslant2$）的图像与 $K_1(x)$ 类似。

### §3.4　应用

**例 1**　均匀导热圆柱体，半径为 $a$，高为 $h$，侧面有恒定热流 $q(z)$ 流入，上下底保持温度 $u_0$，$u_0$ 为常数，求柱内的稳定温度分布。

采用柱坐标，以圆柱的对称轴为 $z$ 轴，方向由下底指向上底，以下底中心为原点。由于边界条件与 $\phi$ 无关，故柱内的稳定温度分布 $u(\boldsymbol r)$ 也与 $\phi$ 无关，即 $u(\boldsymbol r)=u(\rho,z)$，它满足定解问题

$$
\nabla^2u
=
\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial u}{\partial\rho}\right)
+
\frac{\partial^2u}{\partial z^2}
=0,
\qquad
\rho<a,\ 0<z<h,
\tag{79a}
$$

$$
\left.\frac{\partial u}{\partial\rho}\right|_{\rho=a}
=
\frac1k q(z),
\tag{79b}
$$

$$
u|_{z=0}=u_0,
\qquad
u|_{z=h}=u_0.
\tag{79c}
$$

为了将边界条件齐次化，令 $u(\rho,z)=u_0+v(\rho,z)$，则 $v(\rho,z)$ 满足的定解问题为

$$
\nabla^2v
=
\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial v}{\partial\rho}\right)
+
\frac{\partial^2v}{\partial z^2}
=0,
\qquad
\rho<a,\ 0<z<h,
\tag{80a}
$$

$$
\left.\frac{\partial v}{\partial\rho}\right|_{\rho=a}
=
\frac1k q(z),
\tag{80b}
$$

$$
v|_{z=0}=0,
\qquad
v|_{z=h}=0.
\tag{80c}
$$

分离变量，令 $v(\rho,z)=R(\rho)Z(z)$，可得方程

$$
R''+\frac1\rho R'-\lambda R=0,
\tag{81}
$$

和本征值问题

$$
Z''+\lambda Z=0,
$$

$$
Z(0)=0,\qquad Z(h)=0.
$$

这个本征值问题是我们所熟悉的。本征值和本征函数为

$$
\lambda_n=\left(\frac{n\pi}{h}\right)^2,
\qquad
Z_n(z)=\sin\frac{n\pi}{h}z,
\qquad n=1,2,\cdots.
\tag{82}
$$

将本征值代回式 (81)，易得其解为

$$
R_n(\rho)=
\left\{
I_0\left(\frac{n\pi}{h}\rho\right),
K_0\left(\frac{n\pi}{h}\rho\right)
\right\}.
$$

但是 $K_0(n\pi\rho/h)$ 在 $\rho=0$ 处有奇性，应该舍弃，故

$$
R_n(\rho)=C_nI_0\left(\frac{n\pi}{h}\rho\right).
\tag{83}
$$

一般解可以取为

$$
v(\rho,z)=
\sum_{n=1}^{\infty}
A_n
\frac{I_0(n\pi\rho/h)}
{I_1(n\pi a/h)}
\sin\frac{n\pi}{h}z.
\tag{84}
$$

其中分母引入 $I_1(n\pi a/h)$ 是为了下面的方便。代入边界条件 (80b)，注意到 $I_0'(x)=I_1(x)$，可得

$$
\sum_{n=1}^{\infty}
A_n\frac{n\pi}{h}
\sin\frac{n\pi}{h}z
=
\frac1k q(z).
\tag{85}
$$

由此解出系数为

$$
A_n=
\frac{2}{n\pi k}
\int_0^h q(z)\sin\frac{n\pi}{h}z\,dz,
\qquad n=1,2,\cdots.
\tag{86}
$$

如果给定 $q(z)$ 具体形式，那么原则上可以由上式计算各系数，然后代回式 (84)，就得到 $v(\rho,z)$，最后结果为

$$
u(\rho,z)
=
u_0+
\sum_{n=1}^{\infty}
A_n
\frac{I_0(n\pi\rho/h)}
{I_1(n\pi a/h)}
\sin\frac{n\pi}{h}z.
\tag{87}
$$

**特例**　如果 $q(z)=q_0\sin(\pi z/h)$，其中 $q_0$ 是常数，则

$$
A_1=\frac{q_0h}{\pi k},
\qquad
A_2=A_3=\cdots=0.
\tag{88}
$$

从而结果简化为

$$
u(\rho,z)
=
u_0+
\frac{q_0h}{\pi k}
\frac{I_0(\pi\rho/h)}
{I_1(\pi a/h)}
\sin\frac{\pi}{h}z.
\tag{89}
$$

**习题**　均匀导热圆柱体，半径为 $a$，高为 $h$，上下底的温度保持零度，侧面温度分布为 $f(\phi)$，求柱内的稳定温度分布。并计算特例 $f(\phi)=u_0\cos\phi$，其中 $u_0$ 为常数。已经取定柱坐标系，$z$ 轴为圆柱的对称轴，方向由下底指向上底，原点在下底中心。

## §4　*半奇数阶柱函数

半奇数阶柱函数有一个重要的性质：它们都是初等函数。

我们在第十章已经证明 $J_{\pm1/2}(x)$ 是初等函数：

$$
J_{1/2}(x)=\sqrt{\frac{2}{\pi x}}\sin x,
\qquad
J_{-1/2}(x)=\sqrt{\frac{2}{\pi x}}\cos x.
\tag{90}
$$

由式 (6) 的第二式，有

$$
x^{-(\nu+1)}J_{\nu+1}(x)
=
-\frac1x\frac{d}{dx}[x^{-\nu}J_\nu(x)].
$$

反复应用上式，并略作整理，易得

$$
J_{\nu+l}(x)
=
(-)^l x^{\nu+l}
\left(\frac1x\frac{d}{dx}\right)^l
[x^{-\nu}J_\nu(x)],
\qquad l\in\mathbb N.
\tag{91}
$$

上式对 $\nu$ 没有任何限制。取 $\nu=1/2$，并利用式 (90) 第一式，可得

$$
J_{l+1/2}(x)
=
(-)^l\sqrt{\frac{2}{\pi}}\,x^{l+1/2}
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{\sin x}{x}\right),
\qquad l\in\mathbb N.
\tag{92}
$$

这显然是初等函数。类似地，由式 (6) 的第一式出发并利用式 (90) 第二式，可以导出

$$
J_{-(l+1/2)}(x)
=
\sqrt{\frac{2}{\pi}}\,x^{l+1/2}
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{\cos x}{x}\right),
\qquad l\in\mathbb N.
\tag{93}
$$

这显然也是初等函数。由于
$N_{l+1/2}(x)=(-)^{l+1}J_{-(l+1/2)}(x)$，
$N_{-(l+1/2)}(x)=(-)^lJ_{l+1/2}(x)$，
故半奇数阶 Neumann 函数也都是初等函数，特别地

$$
N_{l+1/2}(x)
=
(-)^{l+1}\sqrt{\frac{2}{\pi}}\,x^{l+1/2}
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{\cos x}{x}\right),
\qquad l\in\mathbb N.
\tag{94}
$$

结合式 (92) 和 (94)，可得

$$
H_{l+1/2}^{(1)}(x)
=
i(-)^{l+1}\sqrt{\frac{2}{\pi}}\,x^{l+1/2}
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{e^{ix}}{x}\right),
\qquad l\in\mathbb N,
\tag{95}
$$

$$
H_{l+1/2}^{(2)}(x)
=
i(-)^l\sqrt{\frac{2}{\pi}}\,x^{l+1/2}
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{e^{-ix}}{x}\right),
\qquad l\in\mathbb N.
\tag{96}
$$

这些当然也都是初等函数。但是，要算出以上各式中的高阶导数，将各函数写成有限形式，并不是一件容易的事，后面将会再讨论这一问题。


## §5　*球 Bessel 函数

在球坐标系中对 Helmholtz 方程分离变量，得到角向的球函数方程和下面的径向方程：

$$
r^2\frac{d^2y}{dr^2}+2r\frac{dy}{dr}
+[k^2r^2-l(l+1)]y=0,
\qquad l\in\mathbb N.
\tag{97}
$$

令 $x=kr$，上式化为

$$
\frac{d^2y}{dx^2}+\frac2x\frac{dy}{dx}
+\left[1-\frac{l(l+1)}{x^2}\right]y=0.
\tag{98}
$$

这是球 Bessel 方程的标准形式，物理书上也常常直接把式 (97) 称为球 Bessel 方程。

令 $y(x)=\tilde y(x)/\sqrt{x}$，式 (98) 化为

$$
\frac{d^2\tilde y}{dx^2}
+\frac1x\frac{d\tilde y}{dx}
+
\left[
1-\frac{(l+1/2)^2}{x^2}
\right]\tilde y=0.
$$

这是 $l+1/2$ 阶 Bessel 方程，其解为 $l+1/2$ 阶柱函数 $Z_{l+1/2}(x)$，故球 Bessel 方程的解为

$$
y(x)=\frac{Z_{l+1/2}(x)}{\sqrt{x}}.
\tag{99}
$$

定义 $l$ 阶的球 Bessel 函数 $j_l(x)$、球 Neumann 函数 $n_l(x)$、球 Hankel 函数 $h_l^{(1)}(x)$ 和 $h_l^{(2)}(x)$ 如下：

$$
j_l(x)=\sqrt{\frac{\pi}{2x}}J_{l+1/2}(x),
\tag{100}
$$

$$
n_l(x)=\sqrt{\frac{\pi}{2x}}N_{l+1/2}(x),
\tag{101}
$$

$$
h_l^{(1)}(x)=\sqrt{\frac{\pi}{2x}}H_{l+1/2}^{(1)}(x),
\tag{102}
$$

$$
h_l^{(2)}(x)=\sqrt{\frac{\pi}{2x}}H_{l+1/2}^{(2)}(x).
\tag{103}
$$

它们都是球 Bessel 方程 (98) 的解。数学上 $l$ 可以是复数，我们主要考虑 $l\in\mathbb N$ 的情况，物理上用到的主要就是这种情况。

### §5.1　整数阶球 Bessel 函数

当 $l\in\mathbb N$，利用上节得到的半奇数阶柱函数的表达式，可以得到

$$
j_l(x)=(-)^lx^l
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{\sin x}{x}\right),
\qquad l\in\mathbb N,
\tag{104}
$$

$$
n_l(x)=(-)^{l+1}x^l
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{\cos x}{x}\right),
\qquad l\in\mathbb N,
\tag{105}
$$

$$
h_l^{(1)}(x)=i(-)^{l+1}x^l
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{e^{ix}}{x}\right),
\qquad l\in\mathbb N,
\tag{106}
$$

$$
h_l^{(2)}(x)=i(-)^lx^l
\left(\frac1x\frac{d}{dx}\right)^l
\left(\frac{e^{-ix}}{x}\right),
\qquad l\in\mathbb N.
\tag{107}
$$

这些函数都是初等函数。特别地

$$
j_0(x)=\frac{\sin x}{x},
\qquad
n_0(x)=-\frac{\cos x}{x}.
\tag{108}
$$

对于一般的 $l$，如果能将以上各式中的函数写成有限形式，那么也就很容易得到半奇数阶柱函数的有限形式。下面的讨论供有兴趣的读者参考。

存在多种方法可以计算以上各函数的有限形式，但有的方法需要较多的预备知识，以下方法是最方便的一种。

先求解球 Hankel 函数 $h_l^{(1)}(x)$。容易看出，求出高阶导数后它具有以下形式：

$$
h_l^{(1)}(x)
=
\frac{e^{ix}}{x}
\left[
(-i)^{l+1}
+\frac{a_1}{x}
+\frac{a_2}{x^2}
+\cdots
+\frac{a_l}{x^l}
\right],
\tag{109}
$$

且 $a_l=-i(2l-1)!!$（可用于验证下面的计算结果）。我们通过求解方程 (98) 来确定其中的未定系数。作变换 $y(x)=e^{ix}w(x)$，方程 (98) 化为

$$
\frac{d^2w}{dx^2}
+
2\left(\frac1x+i\right)\frac{dw}{dx}
+
\left[
\frac{2i}{x}-\frac{l(l+1)}{x^2}
\right]w=0.
\tag{110}
$$

$x=0$ 仍然是正则奇点，但是我们的目的并不是求一般解（实际上 $y(x)$ 的一般解已经在前面求出，没必要再计算 $w(x)$ 的一般解），而是要求一个具有式 (109) 形式的特解，故可设

$$
w(x)=\sum_{k=0}^{l}\frac{a_k}{x^{k+1}}.
\tag{111}
$$

代入式 (110)，经过计算可以发现，只要上式中的系数满足以下递推关系，方程即可成立：

$$
a_{k+1}
=
i\frac{(l-k)(l+k+1)}{2(k+1)}a_k,
\qquad
k=0,1,\cdots,l-1.
\tag{112}
$$

反复应用这一递推关系，可以得出系数通项为

$$
a_k
=
i^k\frac{(l+k)!}{2^k k!(l-k)!}a_0,
\qquad
k=0,1,\cdots,l.
\tag{113}
$$

比较式 (109) 可知应该取 $a_0=(-i)^{l+1}$，代入上式，即得所有系数的显式，并可以验证 $a_l$ 与前面求导所得一致。整理即得

$$
h_l^{(1)}(x)
=
\sum_{k=0}^{l}
\frac{(l+k)!}{2^k k!(l-k)!}
\frac{\exp[ix-i(l+1-k)\pi/2]}{x^{k+1}}.
\tag{114}
$$

显然，当 $x\in\mathbb R$ 时，$h_l^{(2)}(x)$ 是 $h_l^{(1)}(x)$ 的复共轭。故当 $x\in\mathbb R$ 时，有

$$
h_l^{(2)}(x)
=
\sum_{k=0}^{l}
\frac{(l+k)!}{2^k k!(l-k)!}
\frac{\exp[-ix+i(l+1-k)\pi/2]}{x^{k+1}}.
\tag{115}
$$

然而，由定义可知当 $x\neq0$ 时 $h_l^{(2)}(x)$ 是 $x$ 的解析函数，而上式右边亦然，所以，根据解析开拓理论，上式右边即是 $h_l^{(2)}(x)$ 的正确形式。

取式 (114) 式 (115) 的线性组合，即得

$$
j_l(x)
=
\sum_{k=0}^{l}
\frac{(l+k)!}{2^k k!(l-k)!}
\frac{\cos[x-(l+1-k)\pi/2]}{x^{k+1}}.
\tag{116}
$$

$$
n_l(x)
=
\sum_{k=0}^{l}
\frac{(l+k)!}{2^k k!(l-k)!}
\frac{\sin[x-(l+1-k)\pi/2]}{x^{k+1}}.
\tag{117}
$$

### §5.2　球 Bessel 方程的本征值问题

球 Bessel 方程 (97) 一般来源于在球坐标系中对 Helmholtz 方程分离变量，将 $r$ 改记为 $x$，将方程改写成 Sturm–Liouville 形式，结合边界条件，典型的本征值问题是：

$$
\frac{d}{dx}\left(x^2\frac{dy}{dx}\right)
-l(l+1)y+\lambda x^2y=0,
\qquad (0<x<a)
\tag{118a}
$$

$$
y(0)=0\quad(l=1,2,\cdots)
\quad\text{或}\quad
|y(0)|<\infty\quad(l=0),
\qquad
\alpha y'(a)+\beta y(a)=0,
\tag{118b}
$$

其中 $\alpha、\beta\geqslant0$，但 $\alpha^2+\beta^2\neq0$；$x=0$ 处的边界条件是自然边界条件。与 Sturm–Liouville 方程比较可知 $k(x)=x^2\geqslant0$，$q(x)=l(l+1)\geqslant0$，$\rho(x)=x^2\geqslant0$。根据 Sturm–Liouville 本征值问题的一般结论可知，$\lambda\geqslant0$，故可令 $\lambda=k^2$，其中 $k\geqslant0$。代入方程 (118)，其解为

$$
y(x)=\{j_l(kx),n_l(kx)\}.
\tag{119}
$$

考虑到 $x=0$ 处的自然边界条件，$n_l(kx)$ 显然不符合要求，应该舍弃。故

$$
y(x)=j_l(kx).
\tag{120}
$$

由 $x=a$ 处的边界条件，可得确定本征值的超越方程为

$$
\alpha k j_l'(ka)+\beta j_l(ka)=0.
\tag{121}
$$

上式也可以用 $J_{l+1/2}(ka)$ 及其导数表出。根据关于 Bessel 函数零点的结论，上式有无穷多正根。将这些正根从小到大排列，记作 $\{k_n\}_{n=1}^{\infty}$，满足

$$
0<k_1<k_2<\cdots.
\tag{122}
$$

相应的本征值和本征函数为

$$
\lambda_n=k_n^2,
\qquad
y_n(x)=j_l(k_nx),
\qquad
n=1,2,\cdots.
\tag{123}
$$

作为 Sturm–Liouville 本征值问题的特例，上述本征函数在区间 $[0,a]$ 上满足正交关系：

$$
\int_0^a y_{n'}(x)y_n(x)x^2\,dx
=
\int_0^a j_l(k_{n'}x)j_l(k_nx)x^2\,dx
=0,
\qquad n'\neq n.
\tag{124}
$$

注意到权函数 $\rho(x)=x^2$，所以上式是带权正交。

本征函数的模方由下式给出

$$
\begin{aligned}
\|j_l(k_nx)\|^2
&\equiv
\int_0^a j_l^2(k_nx)x^2\,dx\\
&=
\left[
1-\frac{l(l+1)}{(k_na)^2}
\right]j_l^2(k_na)
+
\frac{a^2}{2k_n}j_l'(k_na)j_l(k_na)
+
\frac{a^3}{2}[j_l'(k_na)]^2.
\end{aligned}
\tag{125}
$$

这可以由式 (41) 导出，也可以对式 (118a)［其中 $y(x)=j_l(k_nx)$］两边乘以 $k_nxj_l'(k_nx)+j_l(k_nx)/2$ 并积分而得出。对第一类或第二类边界条件，上式可以大大简化。

**注**　① 与 Bessel 方程的本征值问题类似，在以上讨论中把 $l$ 当作是固定的。在实际问题中，可能需要对记号加以修正。② 对于第二类齐次边界条件且 $l=0$ 的情况，需要增加一个本征值 $\lambda_{00}=0$ 和相应的本征函数 $y_{00}(x)=1$。易得 $\|y_{00}(x)\|^2=a^3/3$。

作为 Sturm–Liouville 本征值问题的特例，$\{j_l(k_nx)\}_{n=1}^{\infty}$ 在区间 $[0,a]$ 上是完备的，也就是说，区间 $[0,a]$ 上的任一函数 $f(x)$，如果满足与本征值问题同样的边界条件且足够光滑，必可用 $\{j_l(k_nx)\}_{n=1}^{\infty}$ 展开为广义 Fourier 级数

$$
f(x)=\sum_{n=0}^{\infty}f_nj_l(k_nx),
\qquad 0\leqslant x\leqslant a,
\tag{126}
$$

展开系数为

$$
f_n=
\frac{1}{\|j_l(k_nx)\|^2}
\int_0^a f(x)j_l(k_nx)x^2\,dx,
\qquad n=1,2,\cdots.
\tag{127}
$$

### §5.3　球坐标系中 Helmholtz 方程的解

至此，我们得到了 Helmholtz 方程

$$
\nabla^2v+k^2v=0
\tag{128}
$$

在球坐标系中的一系列特解

$$
v(\boldsymbol r)
=
\{j_l(kr),n_l(kr)\}Y_{lm}(\theta,\phi),
\qquad
l=0,1,2,\cdots;
\quad
m=0,\pm1,\cdots,\pm l.
\tag{129}
$$

如果求解区域包含坐标原点 $r=0$，则应该去掉 $n_l(kr)$。Helmholtz 方程中的 $k$ 并不是由研究对象的性质决定的物理参数，而往往是对波动方程、热传导方程等分离出时间变量 $t$ 时引入的，是需要结合边界条件来确定的本征值。原本应该将 $k^2$ 写成 $\lambda$，但物理上常见的边界条件［比如本征值问题 (118) 中的边界条件］导致 $\lambda\geqslant0$，故习惯上直接将 $\lambda$ 写成 $k^2$。由边界条件解出 $k$ 值之后，需要将其代回时间函数 $T(t)$ 的方程并解出 $T(t)$，再写出原问题的一般解，最后由初始条件确定一般解中的系数。因此，Helmholtz 方程不同于 Laplace 方程，我们没有必要写出其一般解。

### §5.4　应用

**例 1**　均匀导热球，半径为 $r_0$，球内的初始温度分布为已知函数 $f(r)$，放入冰水中冷却，使球面保持温度为零度，求球内各处的温度变化。

采用球坐标系，以球心为原点，极轴方向可以任取。球内的温度分布 $u(\boldsymbol r,t)$ 满足定解问题

$$
\frac{\partial u}{\partial t}-a^2\nabla^2u=0,
\qquad r<r_0,\ t>0,
\tag{130a}
$$

$$
u|_{r=r_0}=0,
\tag{130b}
$$

$$
u|_{t=0}=f(r).
\tag{130c}
$$

分离变量，令 $u(\boldsymbol r,t)=v(\boldsymbol r)T(t)$，代入上式中的方程和边界条件，可得常微分方程

$$
T'+k^2a^2T=0
\tag{131}
$$

和本征值问题

$$
\nabla^2v+k^2v=0,
\qquad r<r_0,
\tag{132a}
$$

$$
v|_{r=r_0}=0.
\tag{132b}
$$

这是 Helmholtz 方程的本征值问题，$v(\boldsymbol r)$ 的各种特解由式 (129) 给出，由于边界条件和初始条件都具有球对称性，故 $l=0$，$m=0$。又因为求解区域包含坐标原点 $r=0$，故应该去掉 $n_0(kr)$。于是

$$
v(\boldsymbol r)=v(r)=j_0(kr)=\frac{\sin kr}{kr}.
\tag{133}
$$

由边界条件 (132b) 可得本征值（严格来说 $k_n^2$ 才是本征值）和本征函数为

$$
k_n=\frac{n\pi}{r_0},
\qquad
v_n(r)=j_0(k_nr)=\frac{\sin k_nr}{k_nr},
\qquad
n=1,2,\cdots.
\tag{134}
$$

将本征值代回式 (131)，易得 $T_n(t)=\exp(-k_n^2a^2t)$。将各本征解叠加，得到一般解为

$$
u(r,t)
=
\sum_{n=1}^{\infty}
A_nj_0(k_nr)\exp(-k_n^2a^2t)
=
\sum_{n=1}^{\infty}
A_n\frac{\sin k_nr}{k_nr}\exp(-k_n^2a^2t).
\tag{135}
$$

代入初始条件，得到

$$
f(r)
=
\sum_{n=1}^{\infty}A_nj_0(k_nr)
=
\sum_{n=1}^{\infty}
A_n\frac{\sin k_nr}{k_nr}.
\tag{136}
$$

由于这种情况非常简单，我们就不一定要动用式 (125) 来计算展开系数。事实上，将上式改写为

$$
rf(r)=
\sum_{n=1}^{\infty}\frac{A_n}{k_n}\sin k_nr,
\tag{137}
$$

这就变成了我们非常熟悉的展开式，由此容易得出系数为

$$
A_n=
\frac{2k_n}{r_0}
\int_0^{r_0}rf(r)\sin k_nr\,dr.
\tag{138}
$$

如果利用式 (125) 来计算，结果当然也是一样的。如果 $f(r)$ 给定具体形式，原则上可以由上式计算各系数，然后代回式 (135)，就得到最后的结果。容易看出，当 $t\to\infty$ 时，$u(r,t)\to0$，即达到环境的温度，与 $f(r)$ 的具体形式无关。

## §6　*平面波用球面波展开

暂略。

## §7　*可用柱函数求解的微分方程

暂略。

## 补充习题

1. 尝试直接计算式 (24) 右边的积分，证明它确实给出 $J_m(x)$。即欲证明

$$
\frac{1}{2\pi}\int_{-\pi}^{\pi}
e^{ix\sin\phi-im\phi}\,d\phi
=
J_m(x).
$$

（提示：将 $e^{ix\sin\phi}$ 展开为 Taylor 级数，利用 $\sin\phi=(e^{i\phi}-e^{-i\phi})/2$ 以及 $\{e^{ik\phi}\}_{k=-\infty}^{\infty}$ 的正交性。）

2. 试计算积分

$$
I=\int_0^\infty e^{-ax}J_m(bx)\,dx,
\qquad a,b>0,\ m\in\mathbb N.
$$

答案：

$$
\frac{(\sqrt{a^2+b^2}-a)^m}
{b^m\sqrt{a^2+b^2}}.
$$

3. 计算不定积分

$$
\int x^3Z_0(x)\,dx,
$$

用 $Z_0(x)$ 和 $Z_1(x)$ 表出结果。

答案：

$$
(x^3-4x)Z_1(x)+2x^2Z_0(x)+C.
$$

4. 无穷长均匀导热圆柱体，半径为 $\rho_0$，侧面保持温度为零度，柱内的初始温度分布为已知函数

$$
f(\rho)=u_0\left(1-\frac{\rho^2}{\rho_0^2}\right),
$$

其中 $u_0$ 是常数，求柱内各处的温度变化。

答案：

$$
u(\rho,t)
=
8u_0\sum_{n=1}^{\infty}
\frac{J_0(k_{0n}\rho)}
{x_{0n}^3J_1(x_{0n})}
\exp(-\lambda_{0n}a^2t).
$$

式中各符号的意义与正文一致。

5. 试证明 $I_\nu(x)$ 的递推关系式 (68) 和 (69)。

6. 试证明 $K_\nu(x)$ 与 Hankel 函数的关系 (74)。

7. 试证明 $K_\nu(x)$ 的递推关系式 (75) 和 (76)。

8. 均匀导热圆柱体，半径为 $a$，高为 $h$，上下底绝热，侧面温度分布为 $f(\phi,z)$，求柱内的稳定温度分布。并计算特例 $f(\phi,z)=u_0\cos\phi\cos(\pi z/h)$，其中 $u_0$ 为常数。已经取定柱坐标系，$z$ 轴为圆柱的对称轴，方向由下底指向上底，原点在下底中心。

答案：

$$
u(\rho,\phi,z)
=
u_0
\frac{I_1(\pi\rho/h)}
{I_1(\pi a/h)}
\cos\phi\cos\frac{\pi}{h}z.
$$

9. 均匀导热球，半径为 $r_0$，表明绝热，球内的初始温度分布为已知函数 $f(r)$，求球内各处的温度变化。并考虑 (1) $t\to\infty$ 时的结果和意义；(2) $f(r)=u_0$ 的特例，其中 $u_0$ 为常数。

答案：(2)

$$
u(\boldsymbol r,t)=u(r,t)=u_0.
$$

[^1]: 图像暂略。
[^2]: 图像暂略。
[^3]: 图像暂略。
