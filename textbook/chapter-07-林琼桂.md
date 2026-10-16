# 第七章　分离变量法*

## 目录

- §1 一维波动方程　133
  - §1.1 分离变量　133
  - §1.2 本征值问题　134
  - §1.3 本征函数族的正交完备性　135
  - §1.4 一般解　137
  - §1.5 用初始条件确定系数　138
  - §1.6 解的物理图像　139
- §2 一维热传导方程　139
  - §2.1 分离变量与本征值问题　139
  - §2.2 本征函数族的正交完备性　141
  - §2.3 一般解　142
  - §2.4 用初始条件确定系数　142
  - §2.5 解的物理图像　143
- §3 非齐次方程与非齐次边界条件　143
  - §3.1 非齐次边界条件的齐次化　143
  - §3.2 非齐次方程的本征函数展开法　146
  - §3.3 *本征函数族的选择　148
  - §3.4 边界条件齐次化的一般方法　150
- §4 矩形区域上的稳定场方程　151
- §5 平面极坐标系中的稳定场方程　154
  - §5.1 一般解　155
  - §5.2 用边界条件确定系数　157
  - §5.3 均匀电场中的导体圆柱　159
  - §5.4 Poisson 方程的处理　160
- §6 *解的唯一性与稳定性　162
  - §6.1 解的唯一性　162
  - §6.2 解的稳定性　164
- 补充习题　166

本章开始研究偏微分方程的求解。读者应该记得，求解常微分方程的思路通常是先求出其通解，其中含有若干任意常数（任意常数的个数与微分方程的阶数相同），然后用附加条件（通常是同一点的函数值和导数值，这在数学上称为初始条件，不管自变量是否表示时间）来确定这些任意常数。对于偏微分方程，读者可能会想到按照类似的思路，先求出通解，这时通解中会含有任意函数，然后用定解条件来确定这些任意函数。但是，这通常是做不到的。首先，除了极少数简单情况，人们并不知道怎样求出偏微分方程的通解。其次，即使能求出通解，如果定解条件比较复杂，要想确定其中的任意函数也是非常困难的。所以，对于偏微分方程，只能根据方程和定解条件的各种具体情况分别加以研究。本章介绍的分离变量法（method of separation of variables）是一种比较有效的方法，它在有关的物理和工程问题的求解上被广泛采用。它的基本思路是设法将偏微分方程问题转化为若干个常微分方程问题来求解。对于这一方法，需要通过具体问题的求解逐步积累，才能有深入的认识。不过，应该指出，不管什么方法，能够解析求解的问题都是有限的，分离变量法也不例外。首先，能够分离变量的问题非常有限，这不仅与方程有关，也与定解条件有关。其次，即使分离变量成功，所得到的常微分方程也不一定能够解析求解。

## §1　一维波动方程

本节研究一维波动方程的求解，以第一边值问题为例。

### §1.1　分离变量

考虑两端固定的弦在初始激励下的自由振动，即是给定初始位移和初始速度，求以后各时刻的位移 $u(x,t)$。定解问题为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0,
\qquad (0<x<l,\ t>0)\tag{1a}
$$

$$
u\big|_{x=0}=0,\qquad u\big|_{x=l}=0,
\qquad (t\ge 0)\tag{1b}
$$

$$
u\big|_{t=0}=\varphi(x),\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (0\le x\le l)\tag{1c}
$$

其中 $\varphi(x)$ 和 $\psi(x)$ 是已知函数，分别给出弦上各点的初始位移和初始速度。请注意各式的变量范围。

现在尝试寻找下列形式的特解

$$
u(x,t)=X(x)T(t).\tag{2}
$$

代入式 (1a)，易得

$$
X''(x)/X(x)=T''(t)/a^2T(t),
$$

其中 $'$ 号表示对各自的自变量求导。上式左边与 $t$ 无关，右边与 $x$ 无关，两边相等，故应与 $x、t$ 均无关，即为常数，记作 $-\lambda$，于是得到两个常微分方程：

$$
T''(t)+\lambda a^2T(t)=0,\tag{3}
$$

$$
X''(x)+\lambda X(x)=0.\tag{4}
$$

将式 (2) 代入式 (1b) 和式 (1c)，分别得到

$$
X(0)T(t)=X(l)T(t)=0\tag{5}
$$

和

$$
X(x)T(0)=\varphi(x),\qquad X(x)T'(0)=\psi(x).\tag{6}
$$

先看式 (6)，它是不可能成立的。首先，按式 (6)，$\varphi(x)$ 和 $\psi(x)$ 均与 $X(x)$ 成正比，从而两者成正比，但 $\varphi(x)$ 和 $\psi(x)$ 是可以相当任意的，一般不可能成正比。其次，如果 $\varphi(x)$ 或 $\psi(x)$ 与 $X(x)$ 成正比，则必满足式 (4)，但式 (4) 的解只可能是三角函数（若 $\lambda>0$）、指数函数（若 $\lambda<0$）或线性函数（若 $\lambda=0$），而 $\varphi(x)$ 或 $\psi(x)$ 一般不可能正好是这些函数形式。所以我们暂时不要求式 (2) 满足初始条件 (1c)。再看式 (5)，由于 $T(t)$ 不能恒为零，否则式 (2) 是平庸解，故必须要求

$$
X(0)=0,\qquad X(l)=0.\tag{7}
$$

如果式 (3)、式 (4) 和式 (7) 同时得到满足，则形如式 (2) 的解能够满足方程 (1a) 和边界条件 (1b)。至于初始条件 (1c)，我们稍后再作考虑。

### §1.2　本征值问题

先考虑 $X(x)$ 的方程 (4)，它是常系数二阶线性齐次常微分方程，很容易求出它的通解。过去，我们求解这类常微分方程，是用初始条件 $X(0)=c_0$ 和 $X'(0)=c_1$ 来确定通解中的任意常数。但现在的附加条件是边界条件 (7)。在这样的条件下，非平庸解并不一定存在，除非方程 (4) 中的系数 $\lambda$ 取某些特定值。

实际上，方程 (4) 中的系数 $\lambda$ 是在分离变量时引入的，我们只知道它是常数（即与自变量 $x、t$ 无关），具体的取值则尚不清楚，需要在求解过程中确定。所以，我们也就可以选择适当的 $\lambda$ 值，使得方程 (4) 有满足边界条件 (7) 的非平庸解。

为清楚起见，我们把方程 (4) 和边界条件 (7) 重新写在下面。

$$
X''(x)+\lambda X(x)=0,\tag{8a}
$$

$$
X(0)=0,\qquad X(l)=0.\tag{8b}
$$

这是一个典型的本征值问题。它的特点是常微分方程中含有一个未定参数 $\lambda$，并附有边界条件。使问题存在非平庸解的 $\lambda$ 值称为本征值（eigenvalue），相应的非平庸解称为本征函数（eigenfunction）。我们将在以后总结本征值问题的一般提法并阐述其一般结论。不过，求解一些具体的本征值问题，积累一些经验，对于理解和掌握后面的一般结论无疑是有益的。下面讨论本征值问题 (8) 的解。

1. 如果 $\lambda=0$，则式 (8a) 的解为 $X(x)=Cx+D$，其中 $C、D$ 是任意常数。代入式 (8b)，易得 $C=D=0$，从而 $X(x)\equiv0$，是为平庸解。故 $\lambda=0$ 不是本征值。

2. 如果 $\lambda<0$，令 $\lambda=-\mu^2$（其中 $\mu>0$），则式 (8a) 的解为

$$
X(x)=C\cosh\mu x+D\sinh\mu x,
$$

其中 $C、D$ 是任意常数。代入式 (8b)，由 $X(0)=0$，得 $C=0$；再由 $X(l)=0$，得 $D\sinh\mu l=0$，但 $\mu l>0$，$\sinh\mu l\ne0$，故 $D=0$。于是 $X(x)\equiv0$，是为平庸解。故 $\lambda<0$ 亦不是本征值。

3. 如果 $\lambda>0$，令 $\lambda=\mu^2$（其中 $\mu>0$），则式 (8a) 的解为

$$
X(x)=C\cos\mu x+D\sin\mu x,
$$

其中 $C、D$ 是任意常数。代入式 (8b)，由 $X(0)=0$，得 $C=0$；再由 $X(l)=0$，得 $D\sin\mu l=0$。为得到非平庸解，必须

$$
\sin\mu l=0,\tag{9}
$$

如此则 $D$ 可任意。式 (9) 即是决定本征值的方程，容易解出 $\mu_n=n\pi/l$（其中 $n\in\mathbb N^+$，注意 $\mu>0$），于是得到全部的本征值和本征函数为

$$
\lambda_n=\left(\frac{n\pi}{l}\right)^2,
\qquad X_n(x)=\sin\frac{n\pi}{l}x,
\qquad n\in\mathbb N^+.\tag{10}
$$

其中已取 $D=1$，因为我们还要将 $\lambda_n$ 代回式 (3) 求解 $T_n(t)$，后者会包含两个任意常数，我们关心的是 $u_n(x,t)=X_n(x)T_n(t)$，而不是 $X_n(x)$ 或 $T_n(t)$ 本身，故可以将常数 $D$ 并入 $T_n(t)$ 中。这一点在后面可以看得更清楚。

### §1.3　本征函数族的正交完备性

在继续求解之前，我们先介绍本征函数族 $\{X_n(x)=\sin(n\pi x/l)\}_{n=1}^{\infty}$ 的两个重要性质。这对于随后的讨论是必不可少的。

首先，该本征函数族在区间 $[0,l]$ 上是正交的。这是指

$$
(X_m(x),X_n(x))\equiv\int_0^lX_m(x)X_n(x)\,dx
=\int_0^l\sin\frac{m\pi}{l}x\sin\frac{n\pi}{l}x\,dx=0,
\qquad m\ne n.\tag{11}
$$

这里 $(X_m(x),X_n(x))$ 表示两个函数 $X_m(x)$ 和 $X_n(x)$ 的内积，由随后的积分定义。上式第二个等号只是将 $X_m(x)$ 和 $X_n(x)$ 的具体形式代入，最后一个等号容易由直接计算加以验证，其间用到下列积化和差的三角函数公式

$$
\sin\frac{m\pi}{l}x\sin\frac{n\pi}{l}x
=-\frac12\left[\cos\frac{m+n}{l}\pi x-\cos\frac{m-n}{l}\pi x\right].
$$

其次，可以证明，该本征函数族在区间 $[0,l]$ 上是完备的。这指的是，$\forall f(x)$，它在区间 $[0,l]$ 上具有良好的解析性质（这通常指具有分段连续、连续或可微等解析性质，今后我们对此只用上述笼统说法，而不作细致的描述）且与本征函数族具有相同的边界条件，则一定可以展开为

$$
f(x)=\sum_{n=1}^{\infty}a_n\sin\frac{n\pi}{l}x,
\qquad x\in[0,l],\tag{12}
$$

其中的展开系数为

$$
a_n=\frac{2}{l}\int_0^lf(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{13}
$$

**注**　① 在上面的本征函数族中去掉任何一个，剩下的就不完备了。例如，在一般的三维空间中，任何矢量 $\boldsymbol A$ 都可以用 $\boldsymbol i、\boldsymbol j、\boldsymbol k$ 三个单位矢量来展开，所以后者是完备的，但若只有 $\boldsymbol i、\boldsymbol j$ 两个，则是不完备的。② 即使 $f(x)$ 并不与本征函数族具有相同的边界条件，上述展开式基本上还是正确的，只是在端点可能不成立。对于今后遇到的其它本征函数族，也有类似情况。

在承认本征函数族完备，即展开式 (12) 可以成立的前提下，很容易求出展开系数 (13)。事实上，将式 (12) 中的求和指标 $n$ 写作 $k$，两边乘以 $\sin(n\pi x/l)$，并对 $x$ 从 0 到 $l$ 积分，可得

$$
\int_0^lf(x)\sin\frac{n\pi}{l}x\,dx
=\int_0^l\sum_{k=1}^{\infty}a_k\sin\frac{k\pi}{l}x\sin\frac{n\pi}{l}x\,dx
=\sum_{k=1}^{\infty}a_k\int_0^l\sin\frac{k\pi}{l}x\sin\frac{n\pi}{l}x\,dx,
$$

（当然，其中的级数应当一致收敛，才能逐项积分，对此我们不作深究）注意到本征函数族的正交性，上式右边的级数只有一项不为零，即 $k=n$ 的那一项，该项中的积分为

$$
\|X_n(x)\|^2\equiv(X_n(x),X_n(x))
=\int_0^l\sin^2\frac{n\pi}{l}x\,dx=\frac l2,\tag{14}
$$

由此立得式 (13)。

式 (11) 和式 (14) 可以合并写成

$$
(X_m(x),X_n(x))
=\int_0^l\sin\frac{m\pi}{l}x\sin\frac{n\pi}{l}x\,dx
=\frac l2\delta_{mn},
\qquad \forall m,n\in\mathbb N^+,\tag{15}
$$

其中

$$
\delta_{mn}=\begin{cases}
1,&m=n,\\
0,&m\ne n
\end{cases}\tag{16}
$$

称为 Kronecker $\delta$ 符号。读者应该逐渐熟悉这种简洁的写法。

本书所遇到的本征值问题均属于 Sturm–Liouville 型（后面将会详细讨论），这类本征值问题的本征函数族都是正交和完备的。如果本征函数的形式比较复杂，则其正交性就难以由直接计算加以验证，但可以利用其所满足的方程和边界条件加以证明。至于完备性，其证明比较困难，本书将直接引用结论而不作证明。

对于定义在 $[0,l]$ 上的函数 $f(x)$，我们可以把它延拓为 $[-l,l]$ 上的奇函数，然后再以 $2l$ 为周期延拓为 $(-\infty,+\infty)$ 上的周期函数，记作 $\tilde f(x)$，在 $[0,l]$ 上，$\tilde f(x)=f(x)$。在微积分中，我们知道，周期函数 $\tilde f(x)$ 可以展开为下列 Fourier 级数：

$$
\tilde f(x)=b_0+\sum_{n=1}^{\infty}\left(b_n\cos\frac{n\pi}{l}x+a_n\sin\frac{n\pi}{l}x\right),
\qquad x\in(-\infty,+\infty),\tag{17}
$$

其中的展开系数为

$$
b_0=\frac{1}{2l}\int_{-l}^l\tilde f(x)\,dx,
\qquad
b_n=\frac1l\int_{-l}^l\tilde f(x)\cos\frac{n\pi}{l}x\,dx,
\qquad
a_n=\frac1l\int_{-l}^l\tilde f(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{18}
$$

由于 $\tilde f(x)$ 是 $[-l,l]$ 上的奇函数，故 $b_0=0$，$b_n=0$（$n\in\mathbb N^+$），于是上面的展开式简化为

$$
\tilde f(x)=\sum_{n=1}^{\infty}a_n\sin\frac{n\pi}{l}x,
\qquad x\in(-\infty,+\infty).\tag{19}
$$

而非零系数 $a_n$ 的表示成为

$$
a_n=\frac{2}{l}\int_0^l\tilde f(x)\sin\frac{n\pi}{l}x\,dx
=\frac{2}{l}\int_0^lf(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{20}
$$

当 $x\in[0,l]$，式 (19) 就化为式 (12)，而系数公式 (20) 与 (13) 则完全相同。这样，我们就看清楚了展开式 (12) 与 Fourier 级数的关系。

对于一般的本征值问题，其本征函数族不一定是三角函数，但由于其完备性，任意解析性质良好的函数都可以用它们来展开，这称为广义 Fourier 级数，它一般不能由上述三角函数的 Fourier 级数得到。可见，读者比较熟悉的三角函数的 Fourier 级数只是一种特殊情况。所以，在上面的叙述中，我们是由本征函数族的完备性直接写出式 (12)，而不是借助于三角函数的 Fourier 级数来导出它。后者显然不是一种普遍的思路。

### §1.4　一般解

将式 (10) 中的本征值代入式 (3)，可得

$$
T_n''(t)+\left(\frac{n\pi a}{l}\right)^2T_n(t)=0.\tag{21}
$$

它的解是

$$
T_n(t)=A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t,\tag{22}
$$

其中 $A_n$ 和 $B_n$ 是任意常数。于是我们得到一系列满足方程 (1a) 和边界条件 (1b) 的特解

$$
u_n(x,t)=X_n(x)T_n(t)
=\left(A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t\right)\sin\frac{n\pi}{l}x,
\qquad n\in\mathbb N^+.\tag{23}
$$

这称为弦的本征振动。但是，这些解一般并不满足初始条件 (1c)。除非 $\varphi(x)$ 和 $\psi(x)$ 同时是某个 $\sin(n\pi x/l)$ 的倍数。

现在，我们注意到方程 (1a) 和边界条件 (1b) 都是齐次的，所以上述各特解叠加后仍然满足方程和边界条件，于是我们可以写出下面的一般解

$$
u(x,t)=\sum_{n=1}^{\infty}u_n(x,t)
=\sum_{n=1}^{\infty}\left(A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t\right)\sin\frac{n\pi}{l}x.\tag{24}
$$

它代表弦在两端固定时的最一般振动形式。我们希望，通过适当选取其中的任意常数 $A_n$ 和 $B_n$，可以使给定的初始条件 (1c) 也得到满足。这样一来，它就满足定解问题 (1) 的所有条件了。

**注**　① 一般解 (24) 已不具有分离变量的形式。② 一般解不是通解。通解应该只满足偏微分方程，调节其中的任意函数，就可以满足不同的定解条件（包括边界条件和初始条件）。而上面的一般解除满足偏微分方程，还满足给定的边界条件。如果边界条件有所改变，那么它就不成立了，这时我们需要从头开始重新求解。所以，这个一般解的通用性非常有限，它只能满足不同的初始条件（通过选择任意常数 $A_n$ 和 $B_n$，参看下一小节的讨论）。③ 上面提到，各特解叠加后得到的一般解仍然满足方程和边界条件。这里实际上假定了级数可以逐项求导。但我们知道，在实变的函数项级数中，逐项求导的条件是非常苛刻的。这些条件归结为对 $A_n$ 和 $B_n$ 的限制，而后者由 $\varphi(x)$ 和 $\psi(x)$ 决定（参看下一小节），所以最终归结为对 $\varphi(x)$ 和 $\psi(x)$ 的解析性质的限制，具体情况不再深入讨论。今后我们假定所有形式操作（如逐项求导、逐项积分等）都是合法的，而不对其条件加以深究。

### §1.5　用初始条件确定系数

将一般解 (24) 代入初始条件 (1c)，可得

$$
\sum_{n=1}^{\infty}A_n\sin\frac{n\pi}{l}x=\varphi(x),
\qquad
\sum_{n=1}^{\infty}B_n\frac{n\pi a}{l}\sin\frac{n\pi}{l}x=\psi(x).\tag{25}
$$

现在的问题是，是否适当选取 $A_n$ 和 $B_n$ 就能使上面两式成立？由前面的讨论知道，本征函数族 $\{\sin(n\pi x/l)\}_{n=1}^{\infty}$ 是完备的，任何解析性质良好的函数都可以表示成它们的线性组合。所以，只要 $\varphi(x)$ 和 $\psi(x)$ 具有良好的解析性质，上述问题就有了肯定的答案。参照式 (12) 和 (13)，容易得到

$$
A_n=\frac{2}{l}\int_0^l\varphi(x)\sin\frac{n\pi}{l}x\,dx,
\qquad
B_n=\frac{2}{n\pi a}\int_0^l\psi(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{26}
$$

给定 $\varphi(x)$ 和 $\psi(x)$ 的具体形式，就可以计算出以上积分。将这些系数代回一般解 (24)，就得到定解问题 (1) 的最终答案。

由于以上方法的出发点是特解 (2)，它具有分离变量的形式，即 $x$ 的函数与 $t$ 的函数是分离的（或因式化的），而一般解 (24) 是一系列这样的特解的叠加，故该方法称为分离变量法。它是 Fourier 在 18 世纪初最先引入的，所以也称为 Fourier 方法。实际上，式 (24) 就是 $x$ 的正弦 Fourier 级数。此外，式 (24) 中的每一项都是驻波（参见下一小节的讨论），所以该方法又称为驻波法。特解 (2) 的形式正是受到驻波形式的启发而引入的。

**注**　以上是分离变量法的全过程，读者第一次接触时，可能会觉得比较冗长。应该指出，掌握一般解 (24) 的形式并学会将它应用于不同的初始条件是没有太大用处的。读者应该尝试写下定解问题 (1)，然后自己一步一步地将求解过程做出来，遇到困难，实在做不下去时才参考书本。完成之后，再尝试将定解问题中的方程改为热传导方程（这样 $T_n(t)$ 的形式将有所不同），或将边界条件换成其它类型（这样本征函数的形式将有所不同），然后，重复上面的求解过程。下节的内容可以看作这样的一个例子，读者完全可以自己将结果做出来，再作比较。本节的习题也是类似的问题。通过这样的学习方式，读者可以比较全面地掌握这种方法。这种学习方式还可以推广到类似的内容。实际上，这也是做研究的入门功夫。

### §1.6　解的物理图像

由式 (24) 可以看出，弦上的一般振动是一系列本征振动的叠加，每一本征振动可以写成下面的形式：

$$
u_n(x,t)=\left(A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t\right)\sin\frac{n\pi}{l}x
=C_n\sin k_nx\cos(\omega_nt-\delta_n),\tag{27}
$$

其中

$$
k_n=\frac{n\pi}{l},\qquad
\omega_n=k_na,\qquad
C_n=\sqrt{A_n^2+B_n^2},\qquad
\delta_n=\arctan\frac{B_n}{A_n}.\tag{28}
$$

式 (27) 代表一列驻波。对于给定的 $n$，弦上各点以同样的频率振动，初相位也是各点相同，只是在节点（见下）两侧相差 $\pi$；但振幅则是位置的函数。

$\omega_n$ 只决定于系统的性质（弦的长度、线密度和张力）和边界条件，而与初始条件无关，故称为本征频率。而 $C_n$ 和 $\delta_n$ 则与初始条件有关。

$n=1$ 的振动称为基波，$\omega_1=\pi a/l$。由于 $\omega_n=n\omega_1$，故 $n>1$ 的振动称为 $n$ 次谐波。

驻波 (27) 上有一些点振幅为零，这些点永远不动，称为波节或节点（node），它们满足 $k_nx=n\pi x/l=m\pi$，所以这些点的位置是

$$
x=\frac{m}{n}l,\qquad m=0,1,\cdots,n.\tag{29}
$$

可见 $n$ 次谐波包括端点共有 $n+1$ 个节点，基波除端点外没有节点。振幅最大的点称为波腹（antinode）。由于相邻两个节点之间必有一个波腹，故 $n$ 次谐波共有 $n$ 个波腹。

弦上的一般振动是各次谐波的线性叠加。

**习题**　长为 $l$ 的弹性细杆，上端固定在电梯天花板，杆身竖直，下端自由。电梯下降，当速度为 $v_0$ 时突然停止。不考虑重力的作用，求解杆的振动。

## §2　一维热传导方程

本节研究一维热传导方程的求解，以第二边值问题为例。

### §2.1　分离变量与本征值问题

考虑两端绝热的细杆的热传导问题，即是给定初始温度分布，求以后各时刻的温度 $u(x,t)$。定解问题为

$$
\frac{\partial u}{\partial t}-a^2\frac{\partial^2u}{\partial x^2}=0,
\qquad (0<x<l,\ t>0)\tag{30a}
$$

$$
\left.\frac{\partial u}{\partial x}\right|_{x=0}=0,
\qquad
\left.\frac{\partial u}{\partial x}\right|_{x=l}=0,
\qquad (t\ge0)\tag{30b}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad (0\le x\le l)\tag{30c}
$$

其中 $\varphi(x)$ 是已知函数，给出杆上各点的初始温度分布。请注意各式的变量范围。

现在尝试寻找下列形式的特解

$$
u(x,t)=X(x)T(t).\tag{31}
$$

代入式 (30a) 和 (30b)，仿照 §1 的推导，易得方程

$$
T'(t)+\lambda a^2T(t)=0\tag{32}
$$

和本征值问题

$$
X''(x)+\lambda X(x)=0,\tag{33a}
$$

$$
X'(0)=0,\qquad X'(l)=0,\tag{33b}
$$

其中 $\lambda$ 是在分离变量时引入的常数。

下面先讨论本征值问题 (33) 的解。

1. 如果 $\lambda<0$，令 $\lambda=-\mu^2$（其中 $\mu>0$），则式 (33a) 的解为

$$
X(x)=C\cosh\mu x+D\sinh\mu x,
$$

其中 $C、D$ 是任意常数。代入式 (33b)，得 $\mu D=0$，$\mu C\sinh\mu l=0$，但 $\mu>0、l>0$，故 $C=D=0$。于是 $X(x)\equiv0$，是为平庸解。故 $\lambda<0$ 不是本征值。

2. 如果 $\lambda=0$，则式 (33a) 的解为 $X(x)=Cx+D$，其中 $C、D$ 是任意常数。代入式 (33b)，易得 $C=0$，而 $D$ 可任意。于是得到本征值和本征函数

$$
\lambda_0=0,\qquad X_0(x)=1,\tag{34}
$$

其中已取 $D=1$，理由已在 §1 讨论过。注意这一情况与第一边值问题的区别。如果遗漏了这个解，本征函数族将是不完备的。

3. 如果 $\lambda>0$，令 $\lambda=\mu^2$（其中 $\mu>0$），则式 (33a) 的解为

$$
X(x)=C\cos\mu x+D\sin\mu x,
$$

其中 $C、D$ 是任意常数。代入式 (33b)，由 $X'(0)=0$，得 $D=0$；再由 $X'(l)=0$，得 $C\sin\mu l=0$。为得到非平庸解，必须 $\sin\mu l=0$，如此则 $C$ 可任意。上式是决定本征值的方程，与上节所得相同。容易解出 $\mu_n=n\pi/l$（其中 $n\in\mathbb N^+$），于是得到另一部分本征值和本征函数

$$
\lambda_n=\left(\frac{n\pi}{l}\right)^2,
\qquad X_n(x)=\cos\frac{n\pi}{l}x,
\qquad n\in\mathbb N^+,\tag{35}
$$

其中已取 $C=1$。如果在上式中取 $n=0$，则可得到式 (34) 中的结果。所以，我们可以把所有的本征值和本征函数写成

$$
\lambda_n=\left(\frac{n\pi}{l}\right)^2,
\qquad X_n(x)=\cos\frac{n\pi}{l}x,
\qquad n\in\mathbb N.\tag{36}
$$

但注意，在上面的求解中，$n=0$ 与 $n\ne0$ 两种情况需要分开讨论。

### §2.2　本征函数族的正交完备性

与 §1 类似，我们先讨论本征函数族 $\{X_n(x)=\cos(n\pi x/l)\}_{n=0}^{\infty}$ 的正交完备性。

首先，该本征函数族在区间 $[0,l]$ 上是正交的：

$$
(X_m(x),X_n(x))=\int_0^l\cos\frac{m\pi}{l}x\cos\frac{n\pi}{l}x\,dx=0,
\qquad m,n\in\mathbb N^+,\ m\ne n,\tag{37a}
$$

$$
(X_0(x),X_n(x))=\int_0^l\cos\frac{n\pi}{l}x\,dx=0,
\qquad n\in\mathbb N^+.\tag{37b}
$$

以上两式均可由直接计算加以验证，其间用到下列积化和差的三角函数公式

$$
\cos\frac{m\pi}{l}x\cos\frac{n\pi}{l}x
=\frac12\left[\cos\frac{m+n}{l}\pi x+\cos\frac{m-n}{l}\pi x\right].
$$

其次，该本征函数族在区间 $[0,l]$ 上也是完备的。即 $\forall f(x)$ 在区间 $[0,l]$ 上具有良好的解析性质且与本征函数族具有相同的边界条件，则一定可以展开为

$$
f(x)=\sum_{n=0}^{\infty}a_n\cos\frac{n\pi}{l}x
=a_0+\sum_{n=1}^{\infty}a_n\cos\frac{n\pi}{l}x,
\qquad x\in[0,l],\tag{38}
$$

其中的展开系数为

$$
a_0=\frac1l\int_0^lf(x)\,dx,
\qquad
a_n=\frac2l\int_0^lf(x)\cos\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{39}
$$

注意 $a_0$ 与 $a_n$（$n\in\mathbb N^+$）的表示并不统一，即在 $a_n$ 的表示中令 $n=0$ 并不能得到 $a_0$，所以我们在式 (38) 的后一表示中将 $a_0$ 项单独写出来以示强调。

在承认本征函数族完备，即展开式 (38) 可以成立的前提下，很容易求出展开系数 (39)。事实上，将式 (38) 前一表示中的求和指标 $n$ 写作 $k$，两边乘以 $\cos(n\pi x/l)$，并对 $x$ 从 0 到 $l$ 积分，可得

$$
\int_0^lf(x)\cos\frac{n\pi}{l}x\,dx
=\sum_{k=0}^{\infty}a_k\int_0^l\cos\frac{k\pi}{l}x\cos\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N.
$$

注意到本征函数族的正交性，上式右边的级数只有一项不为零，即 $k=n$ 的那一项，该项中的积分为

$$
\|X_n(x)\|^2=(X_n(x),X_n(x))
=\int_0^l\cos^2\frac{n\pi}{l}x\,dx
=\begin{cases}
l,&n=0,\\
l/2,&n\ne0.
\end{cases}\tag{40}
$$

由此立得式 (39)。

对于定义在 $[0,l]$ 上的函数 $f(x)$，我们可以把它延拓为 $[-l,l]$ 上的偶函数，然后再以 $2l$ 为周期延拓为 $(-\infty,+\infty)$ 上的周期函数，记作 $\tilde f(x)$，在 $[0,l]$ 上，$\tilde f(x)=f(x)$。将周期函数 $\tilde f(x)$ 展开为下列 Fourier 级数：

$$
\tilde f(x)=a_0+\sum_{n=1}^{\infty}\left(a_n\cos\frac{n\pi}{l}x+b_n\sin\frac{n\pi}{l}x\right),
\qquad x\in(-\infty,+\infty),\tag{41}
$$

其中的展开系数为

$$
a_0=\frac{1}{2l}\int_{-l}^l\tilde f(x)\,dx,
\qquad
a_n=\frac1l\int_{-l}^l\tilde f(x)\cos\frac{n\pi}{l}x\,dx,
\qquad
b_n=\frac1l\int_{-l}^l\tilde f(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{42}
$$

由于 $\tilde f(x)$ 是 $[-l,l]$ 上的偶函数，故 $b_n=0$（$n\in\mathbb N^+$），于是上面的展开式简化为

$$
\tilde f(x)=a_0+\sum_{n=1}^{\infty}a_n\cos\frac{n\pi}{l}x,
\qquad x\in(-\infty,+\infty).\tag{43}
$$

而非零系数的表示成为

$$
a_0=\frac1l\int_0^l\tilde f(x)\,dx=\frac1l\int_0^lf(x)\,dx,
\qquad
a_n=\frac2l\int_0^l\tilde f(x)\cos\frac{n\pi}{l}x\,dx
=\frac2l\int_0^lf(x)\cos\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{44}
$$

当 $x\in[0,l]$，式 (43) 就化为式 (38)，而系数公式 (44) 与 (39) 则完全相同。这就是展开式 (38) 与 Fourier 级数的关系。

### §2.3　一般解

将式 (36) 中的本征值代入式 (32)，容易解出

$$
T_0(t)=A_0,
\qquad
T_n(t)=A_n\exp(-\lambda_na^2t)
=A_n\exp\left[-\left(\frac{n\pi a}{l}\right)^2t\right],
\qquad n\in\mathbb N^+,\tag{45}
$$

其中 $A_0$ 和 $A_n$ 是任意常数。于是我们得到一系列满足方程 (30a) 和边界条件 (30b) 的特解

$$
u_0(x,t)=X_0(x)T_0(t)=A_0,
\qquad
u_n(x,t)=X_n(x)T_n(t)=A_n\exp(-\lambda_na^2t)\cos\frac{n\pi}{l}x,
\qquad n\in\mathbb N^+.\tag{46}
$$

将这些特解叠加，得到满足方程和边界条件的一般解

$$
u(x,t)=\sum_{n=0}^{\infty}u_n(x,t)
=A_0+\sum_{n=1}^{\infty}A_n\exp(-\lambda_na^2t)\cos\frac{n\pi}{l}x.\tag{47}
$$

它代表细杆两端绝热时的温度分布的最一般形式。适当选取其中的任意常数 $A_0$ 和 $A_n$，可以满足给定的初始条件 (30c)。

### §2.4　用初始条件确定系数

将一般解 (47) 代入初始条件 (30c)，可得

$$
A_0+\sum_{n=1}^{\infty}A_n\cos\frac{n\pi}{l}x=\varphi(x).\tag{48}
$$

由前面的讨论知道，本征函数族 $\{\cos(n\pi x/l)\}_{n=0}^{\infty}$ 是完备的，所以，只要 $\varphi(x)$ 具有良好的解析性质，上式就可以成立。参照式 (38) 和 (39)，容易得到

$$
A_0=\frac1l\int_0^l\varphi(x)\,dx,
\qquad
A_n=\frac2l\int_0^l\varphi(x)\cos\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N.\tag{49}
$$

给定 $\varphi(x)$ 的具体形式，就可以计算出以上积分。将这些系数代回一般解 (47)，就得到定解问题 (30) 的最终答案。

### §2.5　解的物理图像

我们注意到，一般解 (47) 中除了第一项，各项都有指数衰减因子 $\exp(-\lambda_na^2t)$，这是输运问题的解的特点，它体现了过程的不可逆性。而第一项常数项则是第二类齐次边界条件所特有的。

当 $t\to\infty$，容易得到

$$
u(x,t)\to A_0=\frac1l\int_0^l\varphi(x)\,dx.\tag{50}
$$

注意，在式 (47) 中如果不把 $A_0$ 项单独写出来，而且忽视了 $\lambda_0=0$ 这一事实，则很容易得出 $u(x,t)\to0$ 的错误结果。而且，式 (49) 中 $A_0$ 的表示也与 $A_n$（$n\in\mathbb N^+$）不同，所以把 $A_0$ 项单独写出来有一定的好处。

式 (50) 的物理意义是很清楚的，它表明，长时间后，杆上各点的温度将趋于均匀，并且是初始温度的平均值。这是因为，杆的侧面是绝热的（这是一维热传导问题的前提），并且两端也是绝热的（本问题具有第二类齐次边界条件），所以它与外界没有热量交换，于是热量只能在内部流动，由热传导的 Fourier 定律，热量总是从温度高处流向低处，所以最后温度必定趋于均匀，再由能量守恒定律，最终的温度只能是初始温度的平均值。

**习题**　矩形薄板，板面绝热，$x=0、l$ 两边保持恒定温度零度，$y=0、d$ 两边绝热，初始温度分布为 $f(x,y)=u_0x/l$，其中 $u_0$ 为常数，求以后的温度分布，并讨论 $t\to\infty$ 的极限。

## §3　非齐次方程与非齐次边界条件

前两节介绍的分离变量法适用于齐次方程与齐次边界条件，但许多实际问题具有非齐次方程或非齐次边界条件，所以研究这些问题的求解显然是必要的。如果问题具有齐次边界条件，只是方程是非齐次的，我们可以用改进的分离变量法，即本征函数展开法来求解。但如果边界条件是非齐次的，则不管方程齐次与否，首先应该将边界条件齐次化。如果原来是非齐次方程，边界条件齐次化后，方程一般仍然是非齐次的；即使方程原来是齐次的，边界条件齐次化后，通常也会变成非齐次方程。不管是哪种情况，边界条件齐次化以后，都可以用本征函数展开法来求解。至于初始条件齐次与否，则是无关紧要的。当然，非齐次方程问题的求解并非只有本征函数展开法一种，但由于每种方法都不是轻而易举就能够掌握的，所以本书只介绍一种方法。读者应该在确实掌握一种方法之后，学有余力，再涉猎其它。如果了解很多方法，但都是一知半解，不能学以致用，则是不足为训的。

### §3.1　非齐次边界条件的齐次化

考虑弦振动问题，$x=0$ 端固定，$x=l$ 端作已知的简谐振动 $A\sin\omega t$（其中振幅 $A$ 很小），弦的初始位移和初始速度均为零，且不受外力作用，求解弦的振动。定解问题为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0,
\qquad (0<x<l,\ t>0)\tag{51a}
$$

$$
u\big|_{x=0}=0,
\qquad u\big|_{x=l}=A\sin\omega t,
\qquad (t\ge0)\tag{51b}
$$

$$
u\big|_{t=0}=0,
\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=0,
\qquad (0\le x\le l)\tag{51c}
$$

将边界条件齐次化的思路是构造一个尽可能简单的函数 $u_0(x,t)$，它与 $u(x,t)$ 满足相同的边界条件，然后令

$$
u(x,t)=v(x,t)+u_0(x,t),\tag{52}
$$

其中 $v(x,t)$ 是新的未知函数，它显然满足齐次的边界条件，但它满足的方程和初始条件当然就会不同于 $u(x,t)$。对于本问题，可取

$$
u_0(x,t)=\frac{x}{l}A\sin\omega t,\tag{53}
$$

则 $v(x,t)$ 满足定解问题

$$
\frac{\partial^2v}{\partial t^2}-a^2\frac{\partial^2v}{\partial x^2}
=\frac{\omega^2A}{l}x\sin\omega t,
\qquad (0<x<l,\ t>0)\tag{54a}
$$

$$
v\big|_{x=0}=0,
\qquad v\big|_{x=l}=0,
\qquad (t\ge0)\tag{54b}
$$

$$
v\big|_{t=0}=0,
\qquad
\left.\frac{\partial v}{\partial t}\right|_{t=0}=-\frac{\omega A}{l}x,
\qquad (0\le x\le l).\tag{54c}
$$

这样边界条件就齐次化了，但方程变成非齐次的。不过这不要紧，下面我们会介绍非齐次方程问题的解法。

如果能够在将边界条件齐次化的同时，保持方程的齐次性，那无疑是值得尝试的。对于上述问题，我们可以取

$$
u_0(x,t)=A\frac{\sin(\omega x/a)}{\sin(\omega l/a)}\sin\omega t.\tag{55}
$$

容易看出，它不仅满足边界条件 (51b)，而且满足偏微分方程 (51a)，这样一来，$v(x,t)$ 的定解问题就具有齐次的方程和边界条件：

$$
\frac{\partial^2v}{\partial t^2}-a^2\frac{\partial^2v}{\partial x^2}=0,
\qquad (0<x<l,\ t>0)\tag{56a}
$$

$$
v\big|_{x=0}=0,
\qquad v\big|_{x=l}=0,
\qquad (t\ge0)\tag{56b}
$$

$$
v\big|_{t=0}=0,
\qquad
\left.\frac{\partial v}{\partial t}\right|_{t=0}
=-\omega A\frac{\sin(\omega x/a)}{\sin(\omega l/a)},
\qquad (0\le x\le l).\tag{56c}
$$

这是 §1 解过的定解问题的一种特殊情况，它的一般解就是式 (24)，我们把它重新写在下面。

$$
v(x,t)=\sum_{n=1}^{\infty}\left(A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t\right)\sin\frac{n\pi}{l}x.\tag{57}
$$

代入初始条件 (56c)，得

$$
\sum_{n=1}^{\infty}A_n\sin\frac{n\pi}{l}x=0,
\qquad
\sum_{n=1}^{\infty}B_n\frac{n\pi a}{l}\sin\frac{n\pi}{l}x
=-\omega A\frac{\sin(\omega x/a)}{\sin(\omega l/a)}.\tag{58}
$$

容易求出

$$
A_n=0,
\qquad
B_n=\frac{2\omega A}{la}\frac{(-1)^{n+1}}{(\omega/a)^2-(n\pi/l)^2}.\tag{59}
$$

最后得到

$$
u(x,t)=\frac{A}{\sin(\omega l/a)}\sin\omega t\sin\frac{\omega x}{a}
+\frac{2\omega A}{la}\sum_{n=1}^{\infty}
\frac{(-1)^{n+1}}{(\omega/a)^2-(n\pi/l)^2}
\sin\frac{n\pi a}{l}t\sin\frac{n\pi}{l}x.\tag{60}
$$

容易看出，如果 $\omega=n\pi a/l$，即右端点的振动频率等于弦的某一本征振动频率，就会发生共振，因为此时式 (60) 中有两项分母变成零。实际上，这时式 (55) 是不能使用的，因为它的分母也变成零。不过用式 (53) 就没有这个问题。

这里简单介绍一下式 (55) 的求法。当然通过观察和尝试也可以得到它，但毕竟有一定的困难。由式 (51b) 的形式，我们可以尝试下列函数（这也有猜测的成分，但相对来说比较有章可循）

$$
u_0(x,t)=F(x)\sin\omega t.\tag{61}
$$

要求它满足式 (51a) 和 (51b)，就得到

$$
F''(x)+\frac{\omega^2}{a^2}F(x)=0,\tag{62a}
$$

$$
F(0)=0,\qquad F(l)=A.\tag{62b}
$$

式 (62a) 的通解是 $F(x)=C\cos(\omega x/a)+D\sin(\omega x/a)$，代入边界条件 (62b)，易得 $C=0$，$D=A/\sin(\omega l/a)$。如此即得 $F(x)=A\sin(\omega x/a)/\sin(\omega l/a)$，代入式 (61)，即得式 (55)。

给定一个具有非齐次边界条件的一维波动或热传导方程的定解问题，如何选择式 (52) 中的 $u_0(x,t)$，使得在将边界条件齐次化的同时保持方程的齐次性（如果方程原来就是齐次的），或将方程也齐次化（如果方程原来也是非齐次的）？应该指出，这通常是做不到的。所以我们一般只能将边界条件齐次化，然后求解非齐次方程的定解问题。另一方面，虽然在某些情况下可以找到函数 $u_0(x,t)$ 满足以上要求，但如果该函数的形式非常复杂，那么 $v(x,t)$ 的定解问题就会具有复杂的初始条件，这样一来，计算 $v(x,t)$ 就会有很大的工作量，这时我们还不如找一个简单的函数 $u_0(x,t)$，只求将边界条件齐次化。后面我们会介绍如何构造这样的简单函数 $u_0(x,t)$。当然，这样的函数不是唯一的，如果选择的 $u_0(x,t)$ 不同，则最后解的形式也不一样。但是，可以一般地证明解的唯一性。所以，尽管解的形式不一样，但它们实际上是相同的（直接验证这一点往往是不太方便的）。所以，我们应该选择尽可能简单的 $u_0(x,t)$。

再举一个例子。弦振动方程第一边值问题的最一般形式是

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f(x,t),
\qquad (0<x<l,\ t>0)\tag{63a}
$$

$$
u\big|_{x=0}=\mu(t),
\qquad u\big|_{x=l}=\nu(t),
\qquad (t\ge0)\tag{63b}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (0\le x\le l)\tag{63c}
$$

其中 $\mu(t)、\nu(t)、\varphi(x)$ 和 $\psi(x)$ 都是已知函数。为了将边界条件齐次化，取

$$
u_0(x,t)=\mu(t)+\frac{x}{l}[\nu(t)-\mu(t)],\tag{64}
$$

则 $v(x,t)$ 的定解问题为

$$
\frac{\partial^2v}{\partial t^2}-a^2\frac{\partial^2v}{\partial x^2}
=f(x,t)-\mu''(t)-\frac{x}{l}[\nu''(t)-\mu''(t)],\tag{65a}
$$

$$
v\big|_{x=0}=0,
\qquad v\big|_{x=l}=0,\tag{65b}
$$

$$
v\big|_{t=0}=\varphi(x)-\mu(0)-\frac{\nu(0)-\mu(0)}{l}x,
\qquad
\left.\frac{\partial v}{\partial t}\right|_{t=0}
=\psi(x)-\mu'(0)-\frac{\nu'(0)-\mu'(0)}{l}x.\tag{65c}
$$

这里边界条件已经齐次化，至于方程和初始条件则与原来没有本质差异，只是形式上复杂了一些。我们当然也可以选择

$$
u_0(x,t)=\mu(t)+\frac{x^2}{l^2}[\nu(t)-\mu(t)]\tag{66}
$$

来将边界条件齐次化，但这显然不如式 (64) 来得简单。

### §3.2　非齐次方程的本征函数展开法

假设边界条件已经齐次化，则弦振动方程第一边值问题的一般形式是

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f(x,t),
\qquad (0<x<l,\ t>0)\tag{67a}
$$

$$
u\big|_{x=0}=0,
\qquad u\big|_{x=l}=0,
\qquad (t\ge0)\tag{67b}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (0\le x\le l)\tag{67c}
$$

事实上，式 (65) 就具有这样的形式。这也可以看作两端固定的弦的受迫振动问题。现在考虑这一定解问题的求解。

$u(x,t)$ 是 $x$ 和 $t$ 两个变量的函数，如果固定时刻 $t$，那么它就只是 $x$ 的函数，这时它可以展开为某一完备的本征函数族的线性组合，比如前两节中出现的本征函数族 $\{\sin(n\pi x/l)\}_{n=1}^{\infty}$ 或 $\{\cos(n\pi x/l)\}_{n=0}^{\infty}$。当时刻 $t$ 变化时，展开系数也就相应地不同，所以在一般时刻 $t$ 也应该可以作这样的展开，只不过展开系数是 $t$ 的函数，即

$$
u(x,t)=\sum_{n=1}^{\infty}T_n(t)\sin\frac{n\pi}{l}x,\tag{68}
$$

或

$$
u(x,t)=\sum_{n=0}^{\infty}T_n(t)\cos\frac{n\pi}{l}x.\tag{69}
$$

当然还可以展开为其它的形式。我们的思路是，通过选择适当的函数族 $\{T_n(t)\}$，可以使得这样的展开式满足定解问题 (67)。如此，则展开式 (69) 显然是不可取的，因为它不满足边界条件 (67b)。但是满足边界条件的展开式也可能还有其它的形式，我们稍后会进一步讨论选择的依据。现在我们只是指出，应该用相应的齐次方程的定解问题中所解出的本征函数族来展开未知函数 $u(x,t)$。所以，如果相应的齐次方程的定解问题不是我们所熟悉的，就应该先求解它，即应该先令 $f(x,t)=0$ 进行求解（当然不必求一般解和其中的系数，只需求出本征函数族即可），获得本征函数族之后再作类似于式 (68) 的展开。对于目前的问题，展开式 (68) 正好是恰当的选择。它已经满足边界条件 (67b)，现在应该将它代入方程 (67a) 和初始条件 (67c) 以确定函数族 $\{T_n(t)\}$ 所应该满足的条件。为此，把各式右边的函数作相应的展开：

$$
f(x,t)=\sum_{n=1}^{\infty}f_n(t)\sin\frac{n\pi}{l}x,
\qquad
f_n(t)=\frac2l\int_0^lf(x,t)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+;\tag{70a}
$$

$$
\varphi(x)=\sum_{n=1}^{\infty}\varphi_n\sin\frac{n\pi}{l}x,
\qquad
\varphi_n=\frac2l\int_0^l\varphi(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+;\tag{70b}
$$

$$
\psi(x)=\sum_{n=1}^{\infty}\psi_n\sin\frac{n\pi}{l}x,
\qquad
\psi_n=\frac2l\int_0^l\psi(x)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+.\tag{70c}
$$

这里各展开系数 $f_n(t)、\varphi_n$ 和 $\psi_n$ 都是已知的，给定各函数 $f(x,t)、\varphi(x)$ 和 $\psi(x)$ 的具体形式，原则上就可以通过积分将它们计算出来。将这些展开式和式 (68) 一起代入方程 (67a) 和初始条件 (67c)，可得

$$
\sum_{n=1}^{\infty}\left[T_n''(t)+\left(\frac{n\pi a}{l}\right)^2T_n(t)\right]\sin\frac{n\pi}{l}x
=\sum_{n=1}^{\infty}f_n(t)\sin\frac{n\pi}{l}x,\tag{71a}
$$

$$
\sum_{n=1}^{\infty}T_n(0)\sin\frac{n\pi}{l}x
=\sum_{n=1}^{\infty}\varphi_n\sin\frac{n\pi}{l}x,\tag{71b}
$$

$$
\sum_{n=1}^{\infty}T_n'(0)\sin\frac{n\pi}{l}x
=\sum_{n=1}^{\infty}\psi_n\sin\frac{n\pi}{l}x.\tag{71c}
$$

由本征函数的正交性，容易知道，以上各式两边的展开系数应该分别相等（只要将两边同时乘以某个本征函数，比如 $\sin(k\pi x/l)$，然后从 0 到 $l$ 积分，就可以得到结论）。所以有

$$
T_n''(t)+\left(\frac{n\pi a}{l}\right)^2T_n(t)=f_n(t),\tag{72a}
$$

$$
T_n(0)=\varphi_n,
\qquad T_n'(0)=\psi_n.\tag{72b}
$$

这是一系列常微分方程的初值问题，按照常微分方程的解法，可以解得

$$
T_n(t)=\varphi_n\cos\frac{n\pi a}{l}t
+\frac{l\psi_n}{n\pi a}\sin\frac{n\pi a}{l}t
+\frac{l}{n\pi a}\int_0^tf_n(\tau)\sin\left[\frac{n\pi a}{l}(t-\tau)\right]d\tau.\tag{73}
$$

将这一结果代回式 (68)，即得到定解问题 (67) 的解。如果 $f(x,t)=0$，则所有 $f_n(t)=0$，于是上式中最后一项为零，这时上式就退化为 §1 中的式 (22)，如所期望。

以上虽然只讨论了波动方程的非齐次边界条件的齐次化和非齐次方程的求解，但所用的方法完全可以用于处理热传导方程的类似的定解问题。至于稳定场方程，后面再作讨论。

**注**　① 式 (73) 主要用于显示一般结果，并与相应的齐次方程的结果进行比较。但在实际计算上，它并不方便。所以读者不必试图记住这一结果。对于具体的问题，可以用导出式 (73) 的方法直接求解。当 $f_n(t)$ 比较简单的时候，则可以用一些特殊的方法来求解。比如若某一项 $f_n(t)=c$ 是常数，则容易看出方程有一个特解 $c(l/n\pi a)^2$，于是可以写出它的通解为

$$
T_n(t)=A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t+c\left(\frac{l}{n\pi a}\right)^2,
$$

然后再用初始条件 (72b) 确定其中的任意常数 $A_n$ 和 $B_n$。当 $f_n(t)=ct$ 是线性函数时，方程有一个特解为 $ct(l/n\pi a)^2$，所以也容易写出其通解，并用初始条件确定系数。又比如当某一项 $f_n(t)=c\cos\omega t+d\sin\omega t$（其中 $\omega$ 不等于相应的本征频率 $n\pi a/l$）是三角函数时，根据方程 (72a) 的形式，可以假定一个特解为 $C\cos\omega t+D\sin\omega t$，代入其中定出常数 $C$ 和 $D$，从而得到通解为

$$
T_n(t)=A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t
+\frac{c\cos\omega t+d\sin\omega t}{(n\pi a/l)^2-\omega^2},
$$

然后用初始条件确定系数 $A_n$ 和 $B_n$。② 如果读者熟悉 Laplace 变换法，则可以很简捷地求解初值问题 (72)，或其相应的具体形式。Laplace 变换法也可以用于求解偏微分方程的定解问题，本书不作介绍。

### §3.3　*本征函数族的选择

现在再回头看看以上求解的出发点，即展开式 (68)。这一展开式之所以能够成功，是因为其中采用了相应的齐次方程的定解问题中解出的本征函数。如果方程 (67a) 的左边含有 $\partial u/\partial x$ 项，或者其中的系数 $a$ 与 $x$ 有关，那么 (68) 虽然满足边界条件，但当我们将它代入方程而试图重复前面的求解过程时，就会遇到困难。为了更清楚地说明这一问题，并找出解决的途径，我们仍考虑以上弦的受迫振动问题，但假定弦的线密度不是均匀的，这时定解问题 (67) 应该修正为

$$
\rho(x)\frac{\partial^2u}{\partial t^2}-K\frac{\partial^2u}{\partial x^2}=F(x,t),
\qquad (0<x<l,\ t>0)\tag{74a}
$$

$$
u\big|_{x=0}=0,
\qquad u\big|_{x=l}=0,
\qquad (t\ge0)\tag{74b}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (0\le x\le l).\tag{74c}
$$

这里我们把弦中张力写作 $K$，以避免引起符号的混淆。现在我们将式 (68) 代入以上方程，并展开

$$
F(x,t)=\sum_{n=1}^{\infty}F_n(t)\sin\frac{n\pi}{l}x,
\qquad
F_n(t)=\frac2l\int_0^lF(x,t)\sin\frac{n\pi}{l}x\,dx,
\qquad n\in\mathbb N^+,\tag{75}
$$

容易得到

$$
\sum_{n=1}^{\infty}\left[\rho(x)T_n''(t)+K\left(\frac{n\pi}{l}\right)^2T_n(t)\right]\sin\frac{n\pi}{l}x
=\sum_{n=1}^{\infty}F_n(t)\sin\frac{n\pi}{l}x.\tag{76}
$$

由于左边的方括号中含有 $x$ 的函数 $\rho(x)$，我们就无法由上式导出一个类似于式 (72a) 的方程，这就是困难所在。当然，我们可以退而求其次，在上式两边同时乘以本征函数 $\sin(m\pi x/l)$，从 0 到 $l$ 积分，利用本征函数族的正交性，并定义

$$
\rho_{mn}=\frac2l\int_0^l\rho(x)\sin\frac{m\pi}{l}x\sin\frac{n\pi}{l}x\,dx,\tag{77}
$$

就可以导出

$$
\sum_{n=1}^{\infty}\rho_{mn}T_n''(t)+K\left(\frac{m\pi}{l}\right)^2T_m(t)=F_m(t),
\qquad m\in\mathbb N^+.\tag{78}
$$

这是关于函数族 $\{T_m(t)\}_{m=1}^{\infty}$ 的微分方程组，结合初始条件 (72b)，原则上可以决定这些函数，但实际求解无疑是困难的。所以，在目前的情况下，我们不应该使用展开式 (68)，或者说，不应该使用本征函数族 $\{\sin(n\pi x/l)\}_{n=1}^{\infty}$。

现在，我们形式上写下展开式

$$
u(x,t)=\sum_{n=1}^{\infty}T_n(t)X_n(x).\tag{79}
$$

将它代入方程 (74a) 和边界条件 (74b)，可得

$$
\sum_{n=1}^{\infty}\left[\rho(x)T_n''(t)X_n(x)-KT_n(t)X_n''(x)\right]=F(x,t)\tag{80}
$$

和

$$
\sum_{n=1}^{\infty}T_n(t)X_n(0)=\sum_{n=1}^{\infty}T_n(t)X_n(l)=0.\tag{81}
$$

如果 $X_n''(x)\propto\rho(x)X_n(x)$，并且 $X_n(0)=X_n(l)=0$，那么问题就有望解决了。实际上，如果我们在定解问题 (74) 中令 $F(x,t)=0$，然后用分离变量法求解，即令 $u(x,t)=X(x)T(t)$，那么 $X(x)$ 所满足的正是这样的一个本征值问题：

$$
KX''(x)+\lambda\rho(x)X(x)=0,\tag{82a}
$$

$$
X(0)=0,\qquad X(l)=0.\tag{82b}
$$

（相应地，$T(t)$ 则满足方程 $T''(t)+\lambda T(t)=0$，不过以下不需用到该方程。当 $\rho(x)$ 为常数时，这里的方程形式上与 §1 略有不同，但求解的结果则是一样的。）这是一个标准的 Sturm–Liouville 本征值问题（以后会详细讨论）。它与前两节遇到的本征值问题的不同之处是出现了函数 $\rho(x)$，称为权（weight）函数，这是最简单的带权的本征值问题。权函数会改变本征函数族的正交关系的形式（见下）。根据 Sturm–Liouville 本征值问题的一般结论和该本征值问题的具体形式可知

(1) 所有的本征值都是正的：$\lambda_n>0$；

(2) 存在无穷多分立的本征值：$\lambda_1<\lambda_2<\cdots<\lambda_n<\lambda_{n+1}<\cdots$，且 $\lim_{n\to\infty}\lambda_n=+\infty$；

(3) 对应于不同本征值的本征函数在区间 $[0,l]$ 上带权正交：

$$
\int_0^lX_m(x)X_n(x)\rho(x)\,dx=0,
\qquad m\ne n;\tag{83}
$$

(4) 本征函数族 $\{X_n(x)\}_{n=1}^{\infty}$ 在区间 $[0,l]$ 上是完备的。

当 $\rho(x)$ 为常数时，结论 (3) 就退化为前面两节遇到的情况，即普通正交。

结论 (4) 保证了展开式 (79) 是可行的。

将式 (82a) 代入式 (80)，可得

$$
\sum_{n=1}^{\infty}[T_n''(t)+\lambda_nT_n(t)]\rho(x)X_n(x)=F(x,t).\tag{84}
$$

利用结论 (3) 就可以导出 $T_n(t)$ 所满足的方程为

$$
T_n''(t)+\lambda_nT_n(t)=f_n(t),\tag{85}
$$

其中

$$
f_n(t)=\frac{\int_0^lX_n(x)F(x,t)\,dx}{\int_0^lX_n^2(x)\rho(x)\,dx}.\tag{86}
$$

将展开式 (79) 代入初始条件 (74c)，利用结论 (3)，容易导出

$$
T_n(0)=\frac{\int_0^lX_n(x)\varphi(x)\rho(x)\,dx}{\int_0^lX_n^2(x)\rho(x)\,dx},
\qquad
T_n'(0)=\frac{\int_0^lX_n(x)\psi(x)\rho(x)\,dx}{\int_0^lX_n^2(x)\rho(x)\,dx}.\tag{87}
$$

给定 $\rho(x)$ 的具体形式，我们就可以求解本征值问题 (82) 得出 $X_n(x)$，再由方程 (85) 和初始条件 (87) 解出 $T_n(t)$，代入展开式 (79)，就可以得到定解问题 (74) 的最终答案。不过，即使 $\rho(x)$ 的形式并不复杂，求解本征值问题 (82) 也是相当困难的。

由以上讨论可以看出，用本征函数展开法求解非齐次方程的定解问题时，应该用相应的齐次方程的定解问题中所解出的本征函数族来展开未知函数 $u(x,t)$。这正是前面所强调的。

### §3.4　边界条件齐次化的一般方法

一维波动方程或热传导方程定解问题的边界条件的最一般形式是

$$
\left.\left(\alpha\frac{\partial u}{\partial x}-\beta u\right)\right|_{x=0}=\mu(t),
\qquad
\left.\left(\gamma\frac{\partial u}{\partial x}+\delta u\right)\right|_{x=l}=\nu(t),\tag{88}
$$

其中 $\alpha、\beta、\gamma、\delta\ge0$，但 $\alpha^2+\beta^2\ne0$，$\gamma^2+\delta^2\ne0$；$\mu(t)$ 和 $\nu(t)$ 是已知函数。注意两式中的符号差异。上述形式包括了第一、二、三类边界条件的各种情况。比如对 $x=0$ 端，若 $\alpha=0$，则是第一类边界条件，若 $\beta=0$ 则是第二类边界条件，若 $\alpha\beta\ne0$ 则是第三类边界条件。由于条件 $\alpha^2+\beta^2\ne0$ 的限制，两者不会同时为零。

前面已经指出，一般情况下，不可能在将边界条件齐次化的同时保持方程的齐次性（如果方程原来就是齐次的），或将方程也齐次化（如果方程原来也是非齐次的）。所以我们一般只能考虑边界条件的齐次化，然后求解非齐次方程的定解问题。

我们已经知道，将边界条件齐次化的思路是构造一个尽可能简单的函数 $u_0(x,t)$，它与未知函数 $u(x,t)$ 满足相同的边界条件，然后按式 (52)，令 $u(x,t)=v(x,t)+u_0(x,t)$，将原来的定解问题化为新的未知函数 $v(x,t)$ 的定解问题，后者显然满足齐次的边界条件。所以问题的关键是如何构造 $u_0(x,t)$。

按式 (88)，$u_0(x,t)$ 作为 $t$ 的函数，必定会与函数 $\mu(t)$ 和 $\nu(t)$ 有关，但式 (88) 只涉及 $x$ 的两个点，所以我们可以假设下列对 $x$ 为线性的形式

$$
u_0(x,t)=A(t)+xB(t),\tag{89}
$$

然后代入式 (88) 确定 $B(t)$ 和 $A(t)$。如果两端都是第二类边界条件，则由于 $\partial u_0(x,t)/\partial x=B(t)$，上式不能奏效，除非 $\mu(t)/\alpha=\nu(t)/\gamma$。这时可以假设下列对 $x$ 为二次的形式

$$
u_0(x,t)=xA(t)+x^2B(t),\tag{90}
$$

再代入式 (88) 确定 $B(t)$ 和 $A(t)$。除了这种两端都是第二类边界条件的情况，都可以用线性形式 (89)（参看下面的讨论）。

将式 (89) 然后代入式 (88)，可得

$$
\alpha B(t)-\beta A(t)=\mu(t),\tag{91a}
$$

$$
(\gamma+l\delta)B(t)+\delta A(t)=\nu(t).\tag{91b}
$$

这是关于 $B(t)$ 和 $A(t)$ 的线性方程组，其系数行列式为

$$
\begin{vmatrix}
\alpha&-\beta\\
\gamma+l\delta&\delta
\end{vmatrix}
=\alpha\delta+\beta(\gamma+l\delta).\tag{92}
$$

只要该行列式不为零，就可以唯一确定 $B(t)$ 和 $A(t)$。由于上式中出现的所有参数都是非负的，所以只当右边两项都是零时，该行列式才等于零。由第二项为零，得 $\beta=0$（因为 $\gamma$ 和 $\delta$ 不同时为零，故 $\gamma+l\delta>0$）。再由第一项为零，得 $\delta=0$（因为已经有 $\beta=0$，则 $\alpha\ne0$）。可见只当 $\beta=\delta=0$，即两端都是第二类边界条件时，式 (89) 才会失效。这时可以改用式 (90)。

下面举几个例子。

**例 1**　第一类边界条件：$u|_{x=0}=\mu(t)$，$u|_{x=l}=\nu(t)$。

将式 (89) 代入以上边界条件，易得 $A(t)=\mu(t)$，$A(t)+lB(t)=\nu(t)$。由此定出 $B(t)$ 和 $A(t)$，即得

$$
u_0(x,t)=\mu(t)+\frac{x}{l}[\nu(t)-\mu(t)].
$$

此即前面用过的式 (64)。

**例 2**　第二类边界条件：$\partial u/\partial x|_{x=0}=\mu(t)$，$\partial u/\partial x|_{x=l}=\nu(t)$。

前已指出，此时应该用式 (90)。代入以上边界条件，易得 $A(t)=\mu(t)$，$A(t)+2lB(t)=\nu(t)$。由此定出 $B(t)$ 和 $A(t)$，即得

$$
u_0(x,t)=x\mu(t)+\frac{x^2}{2l}[\nu(t)-\mu(t)].
$$

**例 3**　混合边界条件：$u|_{x=0}=\mu(t)$，$\partial u/\partial x|_{x=l}=\nu(t)$。

将式 (89) 代入以上边界条件，易得 $A(t)=\mu(t)$，$B(t)=\nu(t)$，故

$$
u_0(x,t)=\mu(t)+x\nu(t).
$$

**习题**　长为 $l$ 的柱形管，左端开放，右端封闭。管外空气中含有某种杂质气体，浓度为 $u_0e^{-\alpha^2t}$，其中 $\alpha>0$。$\alpha\ne\mu_n a$，而 $\mu_n=(2n-1)\pi/2l$，$n\in\mathbb N^+$。杂质向管内扩散，设初始时管内没有该种杂质，求以后时刻该杂质在管内的浓度分布 $u(x,t)$。

## §4　矩形区域上的稳定场方程

前面几节研究了一维波动或热传导方程的定解问题的求解，它们都涉及两个自变量。平面上的稳定场方程的定解问题也涉及两个自变量，因此与上述问题相似。但稳定场方程的定解问题中只有边界条件，没有初始条件，所以其解法也有一些不同之处。本节求解两个矩形区域上的稳定场问题，以加深对分离变量法的了解。矩形区域是平面上最简单的区域。

本节的方法与前两节类似，只有细节上的差别。故本节内容可以看作两个例题，主要供读者自学。读者也可以把定解问题写下来，尝试自己求解，遇到技术上的困难再参考下面的解法。

考虑以下矩形区域上的定解问题

$$
\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}=0,
\qquad (0<x<l,\ 0<y<d) \tag{93a}
$$

$$
u\big|_{x=0}=0,\qquad
\left.\frac{\partial u}{\partial x}\right|_{x=l}=0,
\qquad (0\le y\le d) \tag{93b}
$$

$$
u\big|_{y=0}=\varphi(x),\qquad
u\big|_{y=d}=\psi(x),
\qquad (0\le x\le l). \tag{93c}
$$

物理上，这可以看作一块矩形薄板的稳定温度分布的问题。矩形薄板一边保持零度，一边绝热，另外两边的温度分布已知，那么整块薄板上的温度分布就是确定的，可以通过求解上述定解问题得出结果。

现在尝试寻找下列形式的特解

$$
u(x,y)=X(x)Y(y). \tag{94}
$$

代入式 (93a) 和 (93b)，仿照 §1 的推导，易得方程

$$
Y''(y)-\lambda Y(y)=0 \tag{95}
$$

和本征值问题

$$
X''(x)+\lambda X(x)=0, \tag{96a}
$$

$$
X(0)=0,\qquad X'(l)=0, \tag{96b}
$$

其中 $\lambda$ 是在分离变量时引入的常数。

下面先讨论本征值问题 (96) 的解。

1. 如果 $\lambda<0$，令 $\lambda=-\mu^2$（其中 $\mu>0$），则式 (96a) 的解为

$$
X(x)=C\cosh\mu x+D\sinh\mu x,
$$

其中 $C、D$ 是任意常数。代入式 (96b)，得 $C=0$，$\mu D\cosh\mu l=0$。但 $\cosh\mu l\ge1$，$\mu>0$，故 $D=0$。于是 $X(x)=0$，是为平庸解。故 $\lambda<0$ 不是本征值。

2. 如果 $\lambda=0$，则式 (96a) 的解为 $X(x)=Cx+D$，其中 $C、D$ 是任意常数。代入式 (96b)，易得 $D=0$，$C=0$。于是 $X(x)=0$，是为平庸解。故 $\lambda=0$ 亦不是本征值。

3. 如果 $\lambda>0$，令 $\lambda=\mu^2$（其中 $\mu>0$），则式 (96a) 的解为

$$
X(x)=C\cos\mu x+D\sin\mu x,
$$

其中 $C、D$ 是任意常数。代入式 (96b)，由 $X(0)=0$，得 $C=0$；再由 $X'(l)=0$，得 $\mu D\cos\mu l=0$。为得到非平庸解，必须 $\cos\mu l=0$，如此则 $D$ 可任意。上式是决定本征值的方程。容易解出 $\mu_n=(2n-1)\pi/2l$（其中 $n\in\mathbb N^+$），于是得到本征值和本征函数

$$
\lambda_n=\mu_n^2=\left(\frac{2n-1}{2l}\pi\right)^2,
\qquad
X_n(x)=\sin\mu_nx=\sin\frac{2n-1}{2l}\pi x,
\qquad n\in\mathbb N^+. \tag{97}
$$

其中已取 $D=1$。

**注**　① 如果将上面的 $\mu_n$ 写成 $\mu_n=(2n+1)\pi/2l$，则应该取 $n\in\mathbb N$。这是一个必须注意的细节。② 求出本征函数之后，最好再代回边界条件，看看是否确实满足。这样可以尽量避免计算错误或推算过程中的抄写错误。比如上面的问题，如果最后将本征函数误写成 $\cos\mu_nx$，那么通过边界条件的验证就可以立即发现错误。③ 为了方便初学的读者，我们总是尽可能给出计算的细节，但读者熟悉求解过程之后，某些计算细节只需在草稿纸上演算，而不必详细写出来。比如上面 1. 和 2. 两步，只需指出 $\lambda\le0$ 时没有非平庸解即可。进一步熟悉之后，甚至可以直接写出下面的一般解 (101)。

以上本征函数族的正交完备性可由 Sturm–Liouville 本征值问题的一般结论得到保证。正交性也可以直接验证，实际上，容易证明

$$
\int_0^l\sin\mu_mx\sin\mu_nx\,dx
=\int_0^l\sin\frac{2m-1}{2l}\pi x\sin\frac{2n-1}{2l}\pi x\,dx
=\frac l2\delta_{mn},
\qquad m,n\in\mathbb N^+. \tag{98}
$$

将式 (97) 中的本征值代入式 (95)，可得

$$
Y_n''(y)-\mu_n^2Y_n(y)=0,
\qquad n\in\mathbb N^+. \tag{99}
$$

它的两个线性独立解可以取为 $\{\exp(\mu_ny),\exp(-\mu_ny)\}$，或 $\{\cosh\mu_ny,\sinh\mu_ny\}$。为了下面计算的方便，我们将两个线性独立解取为 $\{\sinh\mu_n(d-y),\sinh\mu_ny\}$。也就是说，我们将上式的通解写作

$$
Y_n(y)=A_n\sinh\mu_n(d-y)+B_n\sinh\mu_ny,
\qquad n\in\mathbb N^+, \tag{100}
$$

其中 $A_n$ 和 $B_n$ 是任意常数。于是我们可以写下满足方程 (93a) 和边界条件 (93b) 的一般解

$$
u(x,y)=\sum_{n=1}^{\infty}
[A_n\sinh\mu_n(d-y)+B_n\sinh\mu_ny]\sin\mu_nx. \tag{101}
$$

将上式代入另外两个边界条件 (93c)，即得

$$
\sum_{n=1}^{\infty}A_n\sinh\mu_nd\sin\mu_nx=\varphi(x),
\qquad
\sum_{n=1}^{\infty}B_n\sinh\mu_nd\sin\mu_nx=\psi(x). \tag{102}
$$

由于本征函数族的完备性，只要适当选取系数 $A_n$ 和 $B_n$，上面两式就可以成立。由式 (98)，易得

$$
A_n=\frac{2}{l\sinh\mu_nd}\int_0^l\varphi(x)\sin\mu_nx\,dx,
\qquad
B_n=\frac{2}{l\sinh\mu_nd}\int_0^l\psi(x)\sin\mu_nx\,dx. \tag{103}
$$

给定 $\varphi(x)$ 和 $\psi(x)$ 的具体形式，就可以计算以上积分。将这些系数代回一般解 (101)，就得到原定解问题的解。

这里一部分边界条件（齐次）用于构成本征值问题，另一部分边界条件（齐次或非齐次）用于确定一般解中的系数，这是与波动或热传导问题的不同之处，后者用初始条件确定系数。所以，就 Laplace 方程而言，不必要求所有的边界条件都是齐次的，否则就只有平庸解了。

下面再考虑全部边界条件都是非齐次的 Poisson 方程定解问题

$$
\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}=f(x,y),
\qquad (0<x<l,\ 0<y<d) \tag{104a}
$$

$$
u\big|_{x=0}=\mu(y),\qquad u\big|_{x=l}=\nu(y),
\qquad (0\le y\le d) \tag{104b}
$$

$$
\left.\frac{\partial u}{\partial y}\right|_{y=0}=\varphi(x),
\qquad
\left.\frac{\partial u}{\partial y}\right|_{y=d}=\psi(x),
\qquad (0\le x\le l). \tag{104c}
$$

对于这一定解问题，我们当然可以仿照波动或热传导问题的解法，令 $u(x,y)=v(x,y)+u_0(x,y)$，并取

$$
u_0(x,y)=\mu(y)+\frac{x}{l}[\nu(y)-\mu(y)],
$$

这样一来，$v(x,y)$ 在 $x=0$ 和 $x=l$ 处就满足齐次的边界条件，虽然方程仍然是非齐次的，但可以用本征函数展开法来求解。不过，由于稳定场问题的所有定解条件都是边界条件，所以更简单的方法是把上述定解问题化为两个问题来求解：

$$
\frac{\partial^2u_1}{\partial x^2}+\frac{\partial^2u_1}{\partial y^2}=f(x,y),
\qquad (0<x<l,\ 0<y<d) \tag{105a}
$$

$$
u_1\big|_{x=0}=0,\qquad u_1\big|_{x=l}=0,
\qquad (0\le y\le d) \tag{105b}
$$

$$
\left.\frac{\partial u_1}{\partial y}\right|_{y=0}=\varphi(x),
\qquad
\left.\frac{\partial u_1}{\partial y}\right|_{y=d}=\psi(x),
\qquad (0\le x\le l), \tag{105c}
$$

$$
\frac{\partial^2u_2}{\partial x^2}+\frac{\partial^2u_2}{\partial y^2}=0,
\qquad (0<x<l,\ 0<y<d) \tag{106a}
$$

$$
u_2\big|_{x=0}=\mu(y),\qquad u_2\big|_{x=l}=\nu(y),
\qquad (0\le y\le d) \tag{106b}
$$

$$
\left.\frac{\partial u_2}{\partial y}\right|_{y=0}=0,
\qquad
\left.\frac{\partial u_2}{\partial y}\right|_{y=d}=0,
\qquad (0\le x\le l). \tag{106c}
$$

显然 $u(x,y)=u_1(x,y)+u_2(x,y)$ 就是原定解问题的解。$u_1(x,y)$ 的定解问题可以用本征函数展开法求解

$$
u_1(x,y)=\sum_{n=1}^{\infty}Y_n(y)\sin\frac{n\pi}{l}x. \tag{107}
$$

至于 $u_2(x,y)$ 的定解问题，可以直接用分离变量法求出其一般解

$$
u_2(x,y)=A_0(l-x)+B_0x
+\sum_{n=1}^{\infty}\left[
A_n\sinh\frac{n\pi}{d}(l-x)+B_n\sinh\frac{n\pi}{d}x
\right]\cos\frac{n\pi}{d}y, \tag{108}
$$

然后用边界条件 (106b) 确定系数。作为练习，请读者做出完整的求解过程。

## §5　平面极坐标系中的稳定场方程

如前所述，分离变量法是一种比较具有普遍性的方法，它可以用于求解多个自变量的偏微分方程定解问题，比如二维、三维空间的各类方程的定解问题。分离变量法能否成功，主要取决于边界的形状。在三维空间，比较容易处理的是长方体、球形或圆柱区域，它们的边界都是坐标面。比如球形区域，它的边界是球面，正好是球坐标系中 $r$ 坐标为常数的坐标面。根据边界的形状，在求解时应该采用相应的曲线坐标系，比如对球形区域上的定解问题，就应该采用球坐标系，如果采用直角坐标系，则无法对边界条件分离变量。分离变量成功之后，我们就需要求解常微分方程及其本征值问题。前面几节遇到的常微分方程，其解都是初等函数，所以很容易处理。然而，在更多的问题中，所出现的微分方程都不是我们所熟悉的，其解并不是初等函数。即使是最简单的弦振动问题，如果弦的线密度不是常数，那么相应的本征值问题就困难得多了。对于曲线坐标系中的定解问题，情况往往就更复杂。如果一个微分方程在许多问题中出现，而且可以解析求解（比如用级数解法），人们就会对它的解进行深入的研究。所谓特殊函数，常常就是这样一些微分方程的解。后面我们会作系统的介绍。

对于平面区域上的问题，如果边界涉及圆弧（比如圆内、圆外或扇形区域上的定解问题），那么显然应该采用平面极坐标系。在平面极坐标系中求解波动或热传导问题，也会涉及特殊函数。但稳定场方程的求解则只需用到初等函数。这可能是曲线坐标系中最简单的问题了。这一问题涉及两个自变量，这与前面几节类似，但边界条件则完全不同，方程的形式也有新的特点。读者应该注意其中的不同之处。

### §5.1　一般解

在平面极坐标系中，Laplace 方程的形式为

$$
\nabla^2u=\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial u}{\partial\rho}\right)
+\frac1{\rho^2}\frac{\partial^2u}{\partial\phi^2}=0. \tag{109}
$$

其中 Laplace 算符的形式可以由直角坐标系中的形式通过坐标变换 $x=\rho\cos\phi$，$y=\rho\sin\phi$ 得到。这种方法思路简单，但计算较为繁琐，当用于计算球坐标系的 Laplace 算符时尤其如此。以后我们会介绍其它方法。暂时不给定定解条件，先求上式的一般解。下面求解时假定 $\phi$ 的取值不受限制。因此，应该注意下面的一般解不适用于横截面为扇形的柱体。

现在尝试寻找下列形式的特解

$$
u(\rho,\phi)=R(\rho)\Phi(\phi). \tag{110}
$$

代入式 (109)，仿照 §1 的推导，易得两个方程

$$
\rho^2R''(\rho)+\rho R'(\rho)-\lambda R(\rho)=0, \tag{111}
$$

$$
\Phi''(\phi)+\lambda\Phi(\phi)=0, \tag{112}
$$

其中 $\lambda$ 是在分离变量时引入的常数。

我们应当注意，在直角坐标系中，坐标与几何点是一一对应的，但曲线坐标系却往往不是这样，比如平面极坐标系中，坐标 $(\rho,\phi)$ 与 $(\rho,\phi+2\pi)$ 代表同一个几何点。一个物理量在一个确定的几何点应该有确定的取值，而不应该以该点数学描述方式的不同为转移，所以，就上述问题而言，必须要求

$$
u(\rho,\phi+2\pi)=u(\rho,\phi). \tag{113}
$$

将式 (110) 代入，得到 $R(\rho)[\Phi(\phi+2\pi)-\Phi(\phi)]=0$，但 $R(\rho)$ 不恒为零，故

$$
\Phi(\phi+2\pi)=\Phi(\phi). \tag{114}
$$

这是由物理量应该满足的单值性和曲线坐标的特性导出的一个边界条件，称为自然边界条件 (natural boundary conditions)。自然边界条件还有其它形式。上式又称为周期性边界条件，它与方程 (112) 构成本征值问题。下面求解这一本征值问题。

1. 如果 $\lambda<0$，令 $\lambda=-\mu^2$（其中 $\mu>0$），则式 (112) 的解为

$$
\Phi(\phi)=Ce^{\mu\phi}+De^{-\mu\phi},
$$

其中 $C、D$ 是任意常数。代入式 (114)，得

$$
C(e^{\mu2\pi}-1)e^{\mu\phi}+D(e^{-\mu2\pi}-1)e^{-\mu\phi}=0.
$$

由于上式对所有的 $\phi$ 值都必须成立，所以必须有 $C(e^{\mu2\pi}-1)=0$，$D(e^{-\mu2\pi}-1)=0$。但 $\mu>0$，故 $e^{\pm\mu2\pi}-1\ne0$，从而 $C=D=0$。于是 $\Phi(\phi)=0$，是为平庸解。故 $\lambda<0$ 不是本征值。

2. 如果 $\lambda=0$，则式 (112) 的解为 $\Phi(\phi)=C+D\phi$，其中 $C、D$ 是任意常数。代入式 (114)，易得 $D=0$，而 $C$ 可任意。于是得到本征值和本征函数

$$
\lambda_0=0,\qquad \Phi_0(\phi)=1, \tag{115}
$$

其中已取 $C=1$。注意这个解不能遗漏，否则本征函数族将是不完备的。

3. 如果 $\lambda>0$，令 $\lambda=\mu^2$（其中 $\mu>0$），则式 (112) 的解为

$$
\Phi(\phi)=Ce^{i\mu\phi}+De^{-i\mu\phi},
$$

其中 $C、D$ 是任意常数。代入式 (114)，得

$$
C(e^{i\mu2\pi}-1)e^{i\mu\phi}+D(e^{-i\mu2\pi}-1)e^{-i\mu\phi}=0.
$$

由于上式对所有的 $\phi$ 值都必须成立，所以必须有 $C(e^{i\mu2\pi}-1)=0$，$D(e^{-i\mu2\pi}-1)=0$，后面一式亦可以写作 $D(e^{i\mu2\pi}-1)=0$。如果 $e^{i\mu2\pi}-1\ne0$，则 $C=D=0$，从而 $\Phi(\phi)=0$，是为平庸解。为得到非平庸解，必须 $e^{i\mu2\pi}=1$，如此则 $C、D$ 均可任意。上式是决定本征值的方程。容易解出 $\mu=m$（其中 $m\in\mathbb N^+$），于是得到另一部分本征值和本征函数

$$
\lambda_m=m^2,\qquad
\Phi_m(\phi)=\{\cos m\phi,\sin m\phi\},
\qquad m\in\mathbb N^+. \tag{116}
$$

现在对应于一个本征值，有两个独立的本征函数。这种一个本征值对应于多个独立本征函数的现象称为简并 (degeneration)，独立的本征函数的个数称为简并度 (degree of degeneracy)。这里简并度为 2。

**注**　对于上面 1. 和 3. 两种情况，我们也可以将解的形式分别取为双曲函数和三角函数，读者不妨用这种形式重新讨论，结果当然是一样的。但容易发现，在目前的边界条件下，用指数形式更加方便。

如果在式 (116) 令 $m=0$，则可得到式 (115) 中的结果。所以，我们可以把所有的本征值和本征函数统一写成

$$
\lambda_m=m^2,\qquad
\Phi_m(\phi)=\{\cos m\phi,\sin m\phi\},
\qquad m\in\mathbb N. \tag{117}
$$

但注意，在上面的求解中，$m=0$ 与 $m\ne0$ 两种情况需要分开讨论，而且，$m=0$ 时实际上只有一个本征函数，即简并度为 1，而 $m\ne0$ 时简并度为 2。式 (117) 中的本征函数族在区间 $[0,2\pi]$ 上是正交完备的，实际上，它们正是标准 Fourier 展开的基。

将上面解得的本征值代入式 (111)，得到

$$
\rho^2R_m''(\rho)+\rho R_m'(\rho)-m^2R_m(\rho)=0. \tag{118}
$$

这是读者在微积分中学过的 Euler 方程，它虽然不是常系数的，但很容易求解，只要作自变量变换 $\rho=e^t$，就可以将其化为常系数方程 $d^2R_m/dt^2-m^2R_m=0$，由此易得 $R_m=\{e^{mt},e^{-mt}\}$（当 $m\ne0$）和 $R_0=\{1,t\}$。变换回原来的自变量，即

$$
R_m(\rho)=\{\rho^m,\rho^{-m}\}\quad(m\in\mathbb N^+),
\qquad
R_0(\rho)=\{1,\ln\rho\}. \tag{119}
$$

更简单的方法是令 $R(\rho)=\rho^k$，代入式 (118) 可以定出 $k=\pm m$。不过，当 $m=0$ 时，这样只能得到一个解。这时我们可以将式 (118) 写成 $[\rho R_0'(\rho)]'=0$，于是 $R_0'(\rho)=B_0/\rho$，从而 $R_0(\rho)=A_0+B_0\ln\rho$，其中 $A_0、B_0$ 是任意常数。这与上面的结果相同。

结合式 (117) 和式 (119)，即得方程 (109) 的一般解为

$$
u(\rho,\phi)=A_0+B_0\ln\rho
+\sum_{m=1}^{\infty}\rho^m(A_m\cos m\phi+B_m\sin m\phi)
+\sum_{m=1}^{\infty}\rho^{-m}(C_m\cos m\phi+D_m\sin m\phi). \tag{120}
$$

### §5.2　用边界条件确定系数

考虑圆内的 Laplace 方程定解问题

$$
\nabla^2u=\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial u}{\partial\rho}\right)
+\frac1{\rho^2}\frac{\partial^2u}{\partial\phi^2}=0,
\qquad (\rho<a) \tag{121a}
$$

$$
u\big|_{\rho=a}=F(\phi). \tag{121b}
$$

其中 $a$ 是圆的半径。由于物理量应该取有限值，而一般解 (120) 中的 $\ln\rho$ 和 $\rho^{-m}$ 在 $\rho=0$ 处有奇性，应该舍弃，故 $B_0=0$，$C_m=D_m=0$（$m\in\mathbb N^+$），从而

$$
u(\rho,\phi)=A_0+\sum_{m=1}^{\infty}
\left(\frac\rho a\right)^m(A_m\cos m\phi+B_m\sin m\phi), \tag{122}
$$

其中用 $(\rho/a)^m$ 代替 $\rho^m$ 是为了下面计算的方便。代入边界条件 (121b)，可得

$$
A_0+\sum_{m=1}^{\infty}(A_m\cos m\phi+B_m\sin m\phi)=F(\phi). \tag{123}
$$

由于本征函数族的完备性，只要适当选取系数，上式总可得到满足。再由本征函数族的正交性，易得

$$
A_0=\frac1{2\pi}\int_0^{2\pi}F(\phi)\,d\phi,
\qquad
A_m=\frac1\pi\int_0^{2\pi}F(\phi)\cos m\phi\,d\phi,
\qquad
B_m=\frac1\pi\int_0^{2\pi}F(\phi)\sin m\phi\,d\phi. \tag{124}
$$

给定 $F(\phi)$ 的具体形式，即可计算上述积分。将这些系数代回式 (122)，即得定解问题 (121) 的解。

类似地，考虑圆外的 Laplace 方程定解问题

$$
\nabla^2u=\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial u}{\partial\rho}\right)
+\frac1{\rho^2}\frac{\partial^2u}{\partial\phi^2}=0,
\qquad (\rho>a) \tag{125a}
$$

$$
u\big|_{\rho=a}=G(\phi). \tag{125b}
$$

由于一般解 (120) 中的 $\ln\rho$ 和 $\rho^m$ 在 $\rho\to\infty$ 处有奇性，应该舍弃，故 $B_0=0$，$A_m=B_m=0$（$m\in\mathbb N^+$），从而

$$
u(\rho,\phi)=A_0+\sum_{m=1}^{\infty}
\left(\frac\rho a\right)^{-m}(C_m\cos m\phi+D_m\sin m\phi). \tag{126}
$$

代入边界条件 (125b)，可得

$$
A_0+\sum_{m=1}^{\infty}(C_m\cos m\phi+D_m\sin m\phi)=G(\phi). \tag{127}
$$

于是可定出系数为

$$
A_0=\frac1{2\pi}\int_0^{2\pi}G(\phi)\,d\phi,
\qquad
C_m=\frac1\pi\int_0^{2\pi}G(\phi)\cos m\phi\,d\phi,
\qquad
D_m=\frac1\pi\int_0^{2\pi}G(\phi)\sin m\phi\,d\phi. \tag{128}
$$

给定 $G(\phi)$ 的具体形式，即可计算上述积分。将这些系数代回式 (126)，即得定解问题 (125) 的解。

如果所考虑的问题是有界的，且不包含原点，比如环域 $a<\rho<b$ 上的定解问题，则一般解 (120) 中的各项均必须保留。不过，如果问题的边界条件对 $x$ 轴具有反射对称性，则可以推知 $u(\rho,-\phi)=u(\rho,\phi)$，从而一般解 (120) 可以简化为

$$
u(\rho,\phi)=A_0+B_0\ln\rho
+\sum_{m=1}^{\infty}(A_m\rho^m+C_m\rho^{-m})\cos m\phi. \tag{129}
$$

应该指出，式 (122) 在原点的邻域上满足 Laplace 方程，而式 (126) 在无穷远点的邻域上满足 Laplace 方程。一般解 (120) 在原点的去心邻域上满足 Laplace 方程。在全平面上满足 Laplace 方程的解则只能是 $u(\rho,\phi)=A_0$，即为常数。

### §5.3　均匀电场中的导体圆柱

本小节利用前面的结果求解下述物理问题：

无限长圆柱导体，半径为 $a$，不带电，置于均匀外电场 $E_0$ 中，$E_0$ 的方向垂直于柱轴。求空间各处的电势分布。

这是一个综合应用题，通过此类问题的求解，读者可以逐步培养怎样从物理到数学，再由数学到物理的分析解决问题的能力。为了清楚起见，下面将求解过程分为几步。

**第一步　物理分析。** 由于静电场中的导体为等势体，故柱内和柱面的电势为常数，可取为零。关键是求解柱外的电势分布。

**注**　由于存在无限长圆柱，选无穷远为电势零点一般是不合适的。在本题中，不同方向的无穷远点显然具有不同的电势，所以更不能取无穷远点为电势零点。

**第二步　建立坐标系。** 取柱轴为 $z$ 轴，任取柱轴上一点为原点，又取电场 $E_0$ 的方向为 $x$ 轴的方向。由于所有的条件不随 $z$ 变化，故任何垂直于 $z$ 轴的平面上的电势分布都是相同的。不妨选 $xy$ 平面为代表来研究。故本题是二维稳定场问题。鉴于边界的形状，只能采用极坐标系。

**第三步　建立定解问题。** 由于柱外没有电荷，故电势 $u(\rho,\phi)$ 满足 Laplace 方程。本题有两个边界，即柱面（或圆周）和无穷远处。柱面的电势已取为零。至于无穷远处，原均匀电场应基本不受影响。由此可列出以下定解问题：

$$
\nabla^2u=\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial u}{\partial\rho}\right)
+\frac1{\rho^2}\frac{\partial^2u}{\partial\phi^2}=0,
\qquad (\rho>a) \tag{130a}
$$

$$
u\big|_{\rho=a}=0, \tag{130b}
$$

$$
u\big|_{\rho\to\infty}=-E_0\rho\cos\phi+u_0. \tag{130c}
$$

其中 $u_0$ 为常数，需要通过求解来确定。

**注**　我们已取柱面为电势零点，故以上常数 $u_0$ 不能任取。注意我们在计算极限过程时只是略去无穷小量，而不是只取主要项（即不仅略去无穷小量，而且略去有限量和低阶的无穷大量）。

**第四步　求解定解问题。** 由前面求得的一般解，考虑到本题具有对 $x$ 轴的反射对称性，可知本题的一般解为式 (129)。代入无穷远处的边界条件 (130c)，则得

$$
A_0+B_0\ln\rho+\sum_{m=1}^{\infty}A_m\rho^m\cos m\phi
=u_0-E_0\rho\cos\phi.
$$

比较两边可得

$$
A_0=u_0,\qquad B_0=0,\qquad A_1=-E_0,\qquad A_m=0\quad(m>1). \tag{131}
$$

于是得到

$$
u(\rho,\phi)=u_0-E_0\rho\cos\phi
+\sum_{m=1}^{\infty}C_m\rho^{-m}\cos m\phi. \tag{132}
$$

再代入柱面的边界条件 (130b)，可得

$$
u_0-E_0a\cos\phi+\sum_{m=1}^{\infty}C_ma^{-m}\cos m\phi=0.
$$

比较两边可得

$$
u_0=0,\qquad C_1=E_0a^2,\qquad C_m=0\quad(m>1). \tag{133}
$$

最后得到

$$
u(\rho,\phi)=-E_0\rho\cos\phi+\frac{E_0a^2}{\rho}\cos\phi. \tag{134}
$$

**第五步　结果分析。** 以上结果中，第一项显然是原均匀电场的电势，而第二项则是由柱面上的感应面电荷产生的。不难求得柱面上的电荷面密度为

$$
\sigma(\phi)=-\epsilon_0\left.\frac{\partial u}{\partial n}\right|_{\rho=a}
=-\epsilon_0\left.\frac{\partial u}{\partial\rho}\right|_{\rho=a}
=2\epsilon_0E_0\cos\phi. \tag{135}
$$

显然，柱面上 $\phi\in(-\pi/2,\pi/2)$ 部分带正电，$\phi\in(\pi/2,3\pi/2)$ 部分带负电，而每单位长度的总带电量为 $Q=\int_0^{2\pi}\sigma(\phi)a\,d\phi=2a\epsilon_0E_0\int_0^{2\pi}\cos\phi\,d\phi=0$，一切如所期望。

**注**　① 本题的求解区域是圆外 $\rho>a$，即无穷远点的邻域，但电势的求解结果式 (134) 中出现的第一项在 $\rho\to\infty$ 时有奇性。这是由于本题中存在均匀场，后者可以看作由无穷远处的源所产生，所以在求解区域中实际上存在有电荷，也就是说电势在求解区域中并不完全满足 Laplace 方程，所以前面的讨论（即圆外的解不包含 $\rho$ 的正幂项和对数项）对此题并不适用。不过，由于源在无穷远处，所以它是通过边界条件影响结果，而不是作为非齐次项出现在方程中。② 方位角 $\phi$（其中 $\phi\in(-\pi/2,\pi/2)$）处 $d\phi$ 角度内的元电荷 $\sigma(\phi)a\,d\phi$ 与方位角 $\pi-\phi$ 处的相应元电荷形成元电偶极子，由于两元电荷的距离为 $2a\cos\phi$，故元电偶极子的偶极矩为 $\sigma(\phi)a\,d\phi\cdot2a\cos\phi=4a^2\epsilon_0E_0\cos^2\phi\,d\phi$，总电偶极矩为 $p=\int_{-\pi/2}^{\pi/2}4a^2\epsilon_0E_0\cos^2\phi\,d\phi=2\pi a^2\epsilon_0E_0$，其方向为 $x$ 轴正向。根据电磁学，一个二维的理想电偶极子 $p$（电偶极子 $p=ql$ 当 $l\to0$、$q\to\infty$ 而 $ql$ 保持固定时称为理想电偶极子）产生的电势为 $\boldsymbol p\cdot\boldsymbol\rho/2\pi\epsilon_0\rho^2=p\cos\phi/2\pi\epsilon_0\rho$，将以上所得 $p$ 值代入，结果正是式 (134) 中的第二项。③ 柱面附近电场的法向分量为 $E_\rho=2E_0\cos\phi$，在 $\phi=0$ 和 $\phi=\pi$ 两处，电场的大小是原均匀场的两倍，所以这两处很容易被击穿。

**思考**　① 如果圆柱原来带电，每单位长度的带电量为 $Q$，试重新求解该问题。② 如果将圆柱换为球，该问题如何求解？（提示：这需要用到球函数，有兴趣的读者可以先自学相关内容。）

### §5.4　Poisson 方程的处理

本小节简单介绍一下非齐次方程、即 Poisson 方程的处理，以下述圆内的定解问题为例。

$$
\nabla^2u=\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial u}{\partial\rho}\right)
+\frac1{\rho^2}\frac{\partial^2u}{\partial\phi^2}=f(\rho,\phi),
\qquad (\rho<a) \tag{136a}
$$

$$
u\big|_{\rho=a}=F(\phi). \tag{136b}
$$

这一问题可以有两种处理方法，分别介绍如下。

如果 $f(\rho,\phi)$ 的形式比较简单，则可以通过观察和试探寻找一个简单函数 $u_0(\rho,\phi)$，使得 $\nabla^2u_0(\rho,\phi)=f(\rho,\phi)$，然后令 $u(\rho,\phi)=v(\rho,\phi)+u_0(\rho,\phi)$，则 $v(\rho,\phi)$ 满足定解问题

$$
\nabla^2v=\frac1\rho\frac{\partial}{\partial\rho}
\left(\rho\frac{\partial v}{\partial\rho}\right)
+\frac1{\rho^2}\frac{\partial^2v}{\partial\phi^2}=0,
\qquad (\rho<a) \tag{137a}
$$

$$
v\big|_{\rho=a}=G(\phi)=F(\phi)-u_0(a,\phi). \tag{137b}
$$

当然，$u_0(\rho,\phi)$ 不是唯一的，因为 $u_0(\rho,\phi)$ 加上任一调和函数（比如 $\rho^m\cos m\phi$ 或 $\rho^m\sin m\phi$，$m\in\mathbb N$）后仍然满足要求。因此应该选择 $u_0(\rho,\phi)$ 使得 $G(\phi)$ 的形式尽可能简单。

如果 $f(\rho,\phi)$ 的形式比较复杂，则难以找到 $u_0(\rho,\phi)$，这时应该用本征函数展开法，即将未知函数 $u(\rho,\phi)$ 展开为

$$
u(\rho,\phi)=A_0(\rho)+\sum_{m=1}^{\infty}
[A_m(\rho)\cos m\phi+B_m(\rho)\sin m\phi], \tag{138}
$$

其中 $A_0(\rho)、A_m(\rho)、B_m(\rho)$ 是未知函数。以上展开式满足周期性边界条件。将非齐次项和边界条件也作类似的展开：

$$
f(\rho,\phi)=a_0(\rho)+\sum_{m=1}^{\infty}
[a_m(\rho)\cos m\phi+b_m(\rho)\sin m\phi], \tag{139}
$$

$$
F(\phi)=\alpha_0+\sum_{m=1}^{\infty}
(\alpha_m\cos m\phi+\beta_m\sin m\phi), \tag{140}
$$

其中

$$
a_0(\rho)=\frac1{2\pi}\int_0^{2\pi}f(\rho,\phi)\,d\phi,
\qquad
a_m(\rho)=\frac1\pi\int_0^{2\pi}f(\rho,\phi)\cos m\phi\,d\phi,
$$

$$
b_m(\rho)=\frac1\pi\int_0^{2\pi}f(\rho,\phi)\sin m\phi\,d\phi \tag{141}
$$

是已知函数，而

$$
\alpha_0=\frac1{2\pi}\int_0^{2\pi}F(\phi)\,d\phi,
\qquad
\alpha_m=\frac1\pi\int_0^{2\pi}F(\phi)\cos m\phi\,d\phi,
\qquad
\beta_m=\frac1\pi\int_0^{2\pi}F(\phi)\sin m\phi\,d\phi \tag{142}
$$

是已知常数。将这些展开式代回定解问题 (136)，可得未知函数 $A_0(\rho)、A_m(\rho)、B_m(\rho)$ 所满足的方程和边界条件：

$$
(\rho A_0')'=\rho a_0,\qquad
\rho^2A_m''+\rho A_m'-m^2A_m=\rho^2a_m,\qquad
\rho^2B_m''+\rho B_m'-m^2B_m=\rho^2b_m, \tag{143}
$$

$$
A_0(a)=\alpha_0,\qquad A_m(a)=\alpha_m,\qquad B_m(a)=\beta_m. \tag{144}
$$

上述方程中的撇号表示对自变量 $\rho$ 求导。上面的方程均不难求解：因为 $A_0(\rho)$ 的方程可以直接积分；$A_m(\rho)$ 的相应齐次方程为 Euler 方程，容易求得其通解，而对应于非齐次项的特解可以用常数变易法求出；$B_m(\rho)$ 的方程与 $A_m(\rho)$ 类似。但我们注意到，对每个未知函数，边界条件只有一个，不足以确定通解中的两个任意常数。实际上还存在另外一组边界条件，由于在 $\rho=0$ 处 $\phi$ 没有定义，所以展开式 (138) 本身就要求

$$
A_m(0)=0,\qquad B_m(0)=0,\qquad m\in\mathbb N^+, \tag{145a}
$$

否则 $u(\rho,\phi)$ 在 $\rho=0$ 处没有定义。上式是自然边界条件。另外，对展开式 (138) 的第一项虽然没有类似的要求，但它总必须有限，即满足

$$
|A_0(0)|<\infty, \tag{145b}
$$

这也是自然边界条件。注意，上式并不意味着 $A_0(\rho)$ 在 $\rho\ne0$ 处可以不必有限，只不过因为 $\rho=0$ 是坐标系的奇点（在该点 $\phi$ 没有定义），偏微分方程的解在此处很容易出现奇性，所以需要对其有限性特别加以强调。实际上，我们在前面由一般解 (120) 推出圆内的解 (122) 时，所用的理由与上述考虑是类似的。对于具体的问题，按上述方法确定 $A_0(\rho)、A_m(\rho)$ 和 $B_m(\rho)$，再代回展开式 (138)，即得原定解问题的解。此处不再作进一步的讨论。

本小节的方法只需略作修改，即可用于圆外或环域上的 Poisson 方程定解问题。

**习题**　扇形薄板，半径为 $a$，圆心角为 $\beta$，用坐标描述即 $\rho\le a$，$0\le\phi\le\beta$。板面绝热，两条直边保持温度为 0 度，弧边温度为 $f(\phi)=u_0\phi(\beta-\phi)$，其中 $u_0$ 为常数，求板面上的温度分布。

## §6　*解的唯一性与稳定性

本节只就有限区间上的一维波动方程讨论解的唯一性与稳定性。

### §6.1　解的唯一性

考虑一维波动方程在有限区间 $[0,l]$ 上的定解问题，其一般形式为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f(x,t),
\qquad (0<x<l,\ t>0) \tag{146a}
$$

$$
\left.\left(\alpha\frac{\partial u}{\partial x}-\beta u\right)\right|_{x=0}=\mu(t),
\qquad
\left.\left(\gamma\frac{\partial u}{\partial x}+\delta u\right)\right|_{x=l}=\nu(t),
\qquad (t\ge0) \tag{146b}
$$

$$
u\big|_{t=0}=\varphi(x),\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (0\le x\le l). \tag{146c}
$$

其中 $f(x,t)、\mu(t)、\nu(t)、\varphi(x)$ 和 $\psi(x)$ 都是已知函数；常数 $\alpha、\beta、\gamma、\delta\ge0$，但 $\alpha^2+\beta^2\ne0$，$\gamma^2+\delta^2\ne0$。以上形式包括了各种可能的边界条件。

今设定解问题 (146) 存在两个解 $u_1(x,t)$ 和 $u_2(x,t)$，则两者皆满足式 (146a)–(146c)。令 $v(x,t)=u_1(x,t)-u_2(x,t)$，则 $v(x,t)$ 满足定解问题

$$
\frac{\partial^2v}{\partial t^2}-a^2\frac{\partial^2v}{\partial x^2}=0,
\qquad (0<x<l,\ t>0) \tag{147a}
$$

$$
\left.\left(\alpha\frac{\partial v}{\partial x}-\beta v\right)\right|_{x=0}=0,
\qquad
\left.\left(\gamma\frac{\partial v}{\partial x}+\delta v\right)\right|_{x=l}=0,
\qquad (t\ge0) \tag{147b}
$$

$$
v\big|_{t=0}=0,\qquad
\left.\frac{\partial v}{\partial t}\right|_{t=0}=0,
\qquad (0\le x\le l). \tag{147c}
$$

现在我们要证明 $v(x,t)\equiv0$，如此则 $u_1(x,t)=u_2(x,t)$，即解是唯一的。

如果两端均不出现第三类边界条件，则体系的能量积分为

$$
E(t)=\frac12\rho\int_0^l\left[
\left(\frac{\partial v}{\partial t}\right)^2
+a^2\left(\frac{\partial v}{\partial x}\right)^2
\right]dx, \tag{148}
$$

其中 $\rho$ 是体系单位长度的质量。上式第一项是动能，而第二项是势能（参见下面的注释）。就下面讨论的需要来说，$E(t)$ 是否表示体系的能量是无关紧要的，所以上式中的因子 $\rho/2$ 也可以去掉。将上式对时间求导得到

$$
\frac{dE}{dt}=\rho\int_0^l\left[
\frac{\partial v}{\partial t}\frac{\partial^2v}{\partial t^2}
+a^2\frac{\partial v}{\partial x}\frac{\partial}{\partial x}
\left(\frac{\partial v}{\partial t}\right)
\right]dx,
$$

其中已对第二项的第二个因子交换了 $x$ 和 $t$ 的求导次序。对上式的第二项分部积分，得到

$$
\frac{dE}{dt}=\rho\int_0^l\frac{\partial v}{\partial t}
\left[\frac{\partial^2v}{\partial t^2}-a^2\frac{\partial^2v}{\partial x^2}\right]dx
+\left.\rho a^2\frac{\partial v}{\partial x}\frac{\partial v}{\partial t}\right|_0^l.
$$

利用方程 (147a) 和边界条件 (147b)（注意现在每一端都是第一或第二类边界条件），即得 $dE/dt=0$，从而 $E(t)=E(0)$。由初始条件 (147c) 易知 $E(0)=0$，于是 $E(t)=0$，但 $E(t)$ 中的被积函数是非负的，故必有 $\partial v/\partial t=0$，$\partial v/\partial x=0$，于是 $v(x,t)$ 为常数，再由初始条件，易知该常数为零，从而 $v(x,t)=0$。

**注**　① 式 (148) 的第一项表示动能是显然的，记该项为 $U$。至于势能，可以以弦振动为例推导如下。在自由振动情况下，弦只受到张力的作用。考虑张力所做的功。作用在 $x\sim x+\Delta x$ 段的张力在位移 $v$ 方向的分量为 $T[\partial v/\partial x|_{x+\Delta x}-\partial v/\partial x|_x]=T(\partial^2v/\partial x^2)\Delta x$，当运动使位移函数由 $v(x,t)$ 变为 $v(x,t)+\delta v(x,t)$ 时，张力所做的功为

$$
\int_0^lT\frac{\partial^2v}{\partial x^2}\delta v\,dx
=T\left.\delta v\frac{\partial v}{\partial x}\right|_0^l
-T\int_0^l\frac{\partial v}{\partial x}\frac{\partial\delta v}{\partial x}\,dx, \tag{149}
$$

其中作了一次分部积分。由于 $\delta v$ 很小，准确到一阶小量，有

$$
\left[\frac{\partial(v+\delta v)}{\partial x}\right]^2
-\left(\frac{\partial v}{\partial x}\right)^2
=2\frac{\partial v}{\partial x}\frac{\partial\delta v}{\partial x}. \tag{150}
$$

如果每一端都是第一类或第二类边界条件，则式 (149) 中的第一项为零，将式 (150) 代入式 (149)，即得张力所做的功为

$$
\int_0^lT\frac{\partial^2v}{\partial x^2}\delta v\,dx
=\frac12T\int_0^l\left(\frac{\partial v}{\partial x}\right)^2dx
-\frac12T\int_0^l\left[\frac{\partial(v+\delta v)}{\partial x}\right]^2dx. \tag{151}
$$

令

$$
V=\frac12T\int_0^l\left(\frac{\partial v}{\partial x}\right)^2dx
=\frac12\rho a^2\int_0^l\left(\frac{\partial v}{\partial x}\right)^2dx. \tag{152}
$$

则张力所做的功等于 $-\delta V$。可见该体系为保守系，且 $V$ 为其势能。另一方面，张力做功引起动能的增加。事实上，利用运动方程 (147a)，并注意到 $\delta v$ 是运动引起的位移，而有 $\delta v=(\partial v/\partial t)\delta t$，就可以得到

$$
\int_0^lT\frac{\partial^2v}{\partial x^2}\delta v\,dx
=\int_0^l\rho\frac{\partial^2v}{\partial t^2}\delta v\,dx
=\frac12\rho\int_0^l\left[\frac{\partial(v+\delta v)}{\partial t}\right]^2dx
-\frac12\rho\int_0^l\left(\frac{\partial v}{\partial t}\right)^2dx
=\delta U. \tag{153}
$$

所以 $\delta(U+V)=0$。这进一步证实了该体系为保守系，且 $V$ 为其势能。② 如果弦的某一端为第三类边界条件，比如 $x=l$ 端与劲度系数为 $k$ 的弹簧连接，当弦的 $x=l$ 端没有位移时，弹簧也处于平衡状态。分析弦上 $(l-\Delta x,l)$ 小段的运动并令 $\Delta x\to0$，可得 $x=l$ 端的边界条件为 $(T\partial v/\partial x+kv)|_{x=l}=0$。利用该边界条件，可以改写式 (149) 中的第一项，从而得到张力所做的功为

$$
\int_0^lT\frac{\partial^2v}{\partial x^2}\delta v\,dx
=\frac12T\int_0^l\left(\frac{\partial v}{\partial x}\right)^2dx
+\frac12kv_l^2
-\frac12T\int_0^l\left[\frac{\partial(v+\delta v)}{\partial x}\right]^2dx
-\frac12k(v_l+\delta v_l)^2, \tag{154}
$$

其中 $v_l=v|_{x=l}$。可见此时势能 $V$ 中应该加入一项 $kv_l^2/2$，后者显然是弹簧的势能，而能量积分也应修正为

$$
E(t)=\frac12\rho\int_0^l\left[
\left(\frac{\partial v}{\partial t}\right)^2
+a^2\left(\frac{\partial v}{\partial x}\right)^2
\right]dx+\frac12kv_l^2. \tag{155}
$$

一般情况下，如果 $x=l$ 端有第三类边界条件，则能量积分应修正为

$$
E(t)=\frac12\rho\int_0^l\left[
\left(\frac{\partial v}{\partial t}\right)^2
+a^2\left(\frac{\partial v}{\partial x}\right)^2
\right]dx+\frac{\rho a^2\delta}{2\gamma}v_l^2. \tag{156}
$$

当 $x=0$ 端有第三类边界条件时，应加上类似的修正项 $(\rho a^2\beta/2\alpha)v_0^2$，其中 $v_0=v|_{x=0}$。利用修正的能量积分，容易证明在相应的边界条件下也有 $dE/dt=0$ 和 $v(x,t)\equiv0$，即相应的定解问题的解是唯一的。

### §6.2　解的稳定性

我们这里讨论的是解关于初始条件的稳定性。首先明确问题的提法。考虑两个定解问题：

$$
\frac{\partial^2u_i}{\partial t^2}-a^2\frac{\partial^2u_i}{\partial x^2}=f(x,t),
\qquad (0<x<l,\ t>0) \tag{157a}
$$

$$
\left.\left(\alpha\frac{\partial u_i}{\partial x}-\beta u_i\right)\right|_{x=0}=\mu(t),
\qquad
\left.\left(\gamma\frac{\partial u_i}{\partial x}+\delta u_i\right)\right|_{x=l}=\nu(t),
\qquad (t\ge0) \tag{157b}
$$

$$
u_i\big|_{t=0}=\varphi_i(x),\qquad
\left.\frac{\partial u_i}{\partial t}\right|_{t=0}=\psi_i(x),
\qquad (0\le x\le l), \tag{157c}
$$

其中 $i=1,2$。现在要研究的问题是：如果 $\varphi_1(x)-\varphi_2(x)$ 和 $\psi_1(x)-\psi_2(x)$ 都很小，则 $u_1(x,t)-u_2(x,t)$ 是否也很小？令 $u(x,t)=u_1(x,t)-u_2(x,t)$，则 $u(x,t)$ 满足定解问题

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0,
\qquad (0<x<l,\ t>0) \tag{158a}
$$

$$
\left.\left(\alpha\frac{\partial u}{\partial x}-\beta u\right)\right|_{x=0}=0,
\qquad
\left.\left(\gamma\frac{\partial u}{\partial x}+\delta u\right)\right|_{x=l}=0,
\qquad (t\ge0) \tag{158b}
$$

$$
u\big|_{t=0}=\varphi(x),\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (0\le x\le l), \tag{158c}
$$

其中 $\varphi(x)=\varphi_1(x)-\varphi_2(x)$，$\psi(x)=\psi_1(x)-\psi_2(x)$。现在问题归结为：如果 $\varphi(x)$ 和 $\psi(x)$ 都很小，则 $u(x,t)$ 是否也很小？我们定义函数 $u(x,t)$ 的模为

$$
\|u\|=\left[\frac1l\int_0^lu^2(x,t)\,dx\right]^{1/2}, \tag{159}
$$

类似可定义 $\|\varphi\|$ 和 $\|\psi\|$。这里在积分外面引入因子 $1/l$ 是为了使得函数的模与函数本身具有同样的量纲，这样就比通常的定义多出了因子 $1/\sqrt l$。注意 $\|\varphi\|$ 和 $\|\psi\|$ 是常数，而 $\|u\|$ 是时间 $t$ 的函数。上面提到的 $\varphi(x)、\psi(x)$ 和 $u(x,t)$ 很小，是指 $\|\varphi\|、\|\psi\|$ 和 $\|u\|$ 很小，这是从平均的意义上来说的，因为即使 $\|\varphi\|$ 很小，$\varphi(x)$ 仍可能在局部有较大的函数值。

下面只讨论两端都是第一类边界条件的情况。这时一般解由式 (24) 给出，其中的系数由式 (26) 给出。由本征函数族的正交性，容易得到

$$
\|u\|^2=\frac12\sum_{n=1}^{\infty}
\left(A_n\cos\frac{n\pi a}{l}t+B_n\sin\frac{n\pi a}{l}t\right)^2
\le\frac12\sum_{n=1}^{\infty}(A_n^2+B_n^2). \tag{160}
$$

利用 Schwarz 不等式

$$
\left|\int_0^lf(x)g(x)\,dx\right|^2
\le\int_0^lf^2(x)\,dx\int_0^lg^2(x)\,dx, \tag{161}
$$

由式 (26) 可得

$$
B_n^2\le\frac{2l^2}{n^2\pi^2a^2}\|\psi\|^2,
$$

于是

$$
\frac12\sum_{n=1}^{\infty}B_n^2
\le\frac{l^2}{\pi^2a^2}\|\psi\|^2\sum_{n=1}^{\infty}\frac1{n^2}
=\frac{l^2}{6a^2}\|\psi\|^2. \tag{162}
$$

其中用到了熟知的结果 $\sum_{n=1}^{\infty}1/n^2=\pi^2/6$。

注意到定解问题中的初始条件和边界条件必须相容，则应该有 $\varphi(0)=\varphi(l)=0$。将式 (26) 中 $A_n$ 的积分表达式作一次分部积分，利用上述条件，即得

$$
A_n=\frac2{n\pi}\int_0^l\varphi'(x)\cos\frac{n\pi}{l}x\,dx.
$$

利用 Schwarz 不等式，由上式可得

$$
A_n^2\le\frac{2l^2}{n^2\pi^2}\|\varphi'\|^2,
$$

于是

$$
\frac12\sum_{n=1}^{\infty}A_n^2
\le\frac{l^2}{\pi^2}\|\varphi'\|^2\sum_{n=1}^{\infty}\frac1{n^2}
=\frac{l^2}{6}\|\varphi'\|^2. \tag{163}
$$

综合以上结果即得

$$
\|u\|^2\le\frac{l^2}{6}\|\varphi'\|^2+\frac{l^2}{6a^2}\|\psi\|^2. \tag{164}
$$

按上式，$\forall\varepsilon>0$，如果 $\|\varphi'\|$ 和 $\|\psi\|$ 都足够小，则在任何时刻都可有 $\|u\|<\varepsilon$，这就是解的稳定性。这种稳定性是长期的，而不仅仅在有限时间内成立。这里有两个问题需要讨论。第一个问题是，乍一看，可能会出现这样一种情况：$\varphi'(x)$ 和 $\psi(x)$ 都很小，但 $\varphi(x)$ 并不小（比如 $\varphi(x)$ 为常数），按上式 $u(x,t)$ 却可以很小，这显然是不合理的。不过，由于边界条件的限制，这种情况实际上是不会出现的。事实上，由于 $\varphi(0)=0$，就有

$$
\varphi(x)=\int_0^x\varphi'(\xi)\,d\xi. \tag{165}
$$

由于 $\varphi'(x)$ 很小，所以 $\varphi(x)$ 必然也很小。第二个问题是，可能 $\varphi(x)$ 和 $\psi(x)$ 都很小，但 $\varphi'(x)$ 是有限的（比如 $\varphi(x)=l\kappa\sin(x/l\kappa)$，$\varphi'(x)=\cos(x/l\kappa)$，其中 $\kappa$ 是一个无量纲数，当 $\kappa$ 很小时，$\varphi(x)$ 也很小，但 $\varphi'(x)$ 是有限的），这时按上面的结果，只能得到 $\|u\|$ 的上限为有限的结论，而不能断定 $u(x,t)$ 很小，而按直观的判断，$u(x,t)$ 应该可以很小，至少我们尚未找到一个实例，其中 $\varphi(x)$ 和 $\psi(x)$ 都很小，且与边界条件相容，但 $\varphi'(x)$ 和 $u(x,t)$ 是有限的。这可能意味着存在比式 (164) 更强的结果，但目前我们还未能找到它。

如果出现第二类边界条件，上面的讨论只需略作修改。即使有第三类边界条件，也可以证明类似的结论。¹

¹ 参看以下文献：林琼桂，一维波动方程的解的长期稳定性，大学物理 **24** (5) (2005) 18-22。

## 补充习题

1. 弹性均匀细杆，长为 $l$，处于竖直方向，上端固定，下端自由，在重力作用下作纵振动。已知 $t=0$ 时杆静止且没有形变，求 $t>0$ 时的位移函数。设杆的横截面、质量密度、杨氏模量等参数均为已知。

2. 弹性均匀细杆，长为 $l$，两端固定，处于竖直方向。由于重力的作用，杆的各部分可发生相对形变而形成纵振动。已知 $t=0$ 时杆静止且没有形变，求 $t>0$ 时的位移函数。设杆的横截面、质量密度、杨氏模量等参数均为已知。

3. 弹性均匀细杆，长为 $l$，左端固定，右端自由，杆静止且没有形变。从 $t=0$ 开始用恒定拉力 $F$ 作用于其右端，求 $t>0$ 时的位移函数。设杆的横截面、质量密度、杨氏模量等参数均为已知。

4. 均匀导热细杆，长为 $l$，两端绝热，杆上有热源 $f(x,t)=f_0\cos(\pi x/l)$，初始温度分布是 $u|_{t=0}=2u_0\cos^2(\pi x/l)$，其中 $x$ 是杆上一点到一端的距离，$0\le x\le l$，$f_0$ 和 $u_0$ 是常数。求 $t>0$ 以后的温度分布。

5. 在铀块中除了中子的扩散过程以外，还存在增殖过程。设中子浓度为 $u$，则单位时间单位体积内产生的中子量为 $\beta u$，其中 $\beta$ 为常数。对于层状铀块，当其厚度达到某一临界值时，中子量将随着时间不断增长，最后发生爆炸。设扩散系数为 $D=a^2$，求临界厚度。

6. 均匀导热薄板在平面上占据区域 $0<x<a$，$0<y<\infty$，板面绝热，边界的温度分布为 $u|_{x=0}=0$，$u|_{x=a}=0$，$u|_{y=0}=u_0$，$u|_{y=\infty}=0$，求板上的稳定温度分布。

7. 细圆环，半径为 $a$，表面绝热，初始温度分布为 $f(\phi)$（$\phi$ 是极坐标，原点在环心）。求环上的温度变化。

8. 半圆形薄板，半径为 $a$，板面绝热，直边温度保持零度，圆周上温度为常数 $u_0$，求稳定状态下板上的温度分布。

9. 求解圆内 Poisson 方程定解问题 $\nabla^2u=-xy$（$\rho<a$），$u|_{\rho=a}=0$。
