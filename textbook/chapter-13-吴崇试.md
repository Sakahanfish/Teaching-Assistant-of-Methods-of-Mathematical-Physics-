# 第十三章 分离变量法

上一章介绍的行波法虽然得到的解表达式简洁、物理图像清晰，但是能够用行波法求解的定解问题有限。因此我们必须寻找新的求解一般偏微分方程定解问题的方法。

偏微分方程和常微分方程的一个显著不同是：满足偏微分方程的未知函数，自变量一定有两个或者更多个，而常微分方程未知函数的自变量只有一个。我们能否设法把目前还不会求解的偏微分方程转化为已经会求解的常微分方程呢？例如，如果满足偏微分方程的未知函数 $u(x,t)$ 是某种特别的形式，例如

$$
u(x,t)=X(x)T(t),
$$

也就是说是由多个一元函数相乘得到，那么每个一元函数 $X(x)$ 和 $T(t)$ 满足的方程应该都是常微分方程，只要能求出 $X(x)$ 和 $T(t)$，就可能求出 $u(x,t)$。这样我们就把求解多元函数 $u(x,t)$ 的问题转变为求解一元函数 $X(x)$ 和 $T(t)$ 的问题，把求解偏微分方程定解问题转变为求解常微分方程问题。本章介绍的分离变量法就是采用的这种思路，是解偏微分方程定解问题基本方法。在学习本章内容之前，要求读者已经能够熟练求解常微分方程定解问题。

我们在解常微分方程定解问题时，通常总是先求出微分方程的特解，由线性无关的特解叠加出通解，而后用定解条件（例如，初始条件）定出适合于该问题的叠加系数。对于一阶线性偏微分方程的求解问题，基本的方法也是转化为一阶线性常微分方程组的求解问题。对于二阶以及更高阶偏微分方程定解问题，情况有些不同：即使可以先求出偏微分方程的通解，但由于其中含有待定函数，一般也难以直接根据定解条件定出。这样，为了求解偏微分方程的定解问题，就必须把求解步骤加以适当的修改。本章介绍的分离变量法，就是先求出满足方程及一部分定解条件的全部特解，把这全部特解叠加起来，再利用其余的定解条件定出叠加系数，从而求出该定解问题的解。

## §13.1 两端固定弦的自由振动

考虑长为 $l$、两端固定的弦的自由振动，方程及定解条件为

$$
\frac{\partial^2u}{\partial t^2}
-a^2\frac{\partial^2u}{\partial x^2}
=0,
\qquad
0<x<l,\ t>0,
\tag{13.1a}
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.1b}
$$

$$
\left.u\right|_{t=0}=\phi(x),
\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad
0\le x\le l.
\tag{13.1c}
$$

这里的方程和边界条件都是齐次的，而初始条件是非齐次的。

求解过程可分解为下列四步。

第一步，分离变量。

我们希望求得的特解具有分离变量的形式，即①

$$
u(x,t)=X(x)T(t).
\tag{13.2}
$$

将 (13.2) 式代入方程，即得

$$
X(x)T''(t)=a^2X''(x)T(t).
$$

由于 $u(x,t)$ 不恒为 0，等式两端除以 $X(x)T(t)$，就有

$$
\frac1{a^2}\frac{T''(t)}{T(t)}
=
\frac{X''(x)}{X(x)}.
$$

注意在这个等式中，左端只是 $t$ 的函数（即与 $x$ 无关），右端只是 $x$ 的函数（即与 $t$ 无关），左端和右端相等，就必须共同等于一个既与 $x$ 无关、又与 $t$ 无关的常数。令这个常数为 $-\lambda$，上面的结果又可以写成

$$
T''(t)+\lambda a^2T(t)=0,
\tag{13.3a}
$$

$$
X''(x)+\lambda X(x)=0.
\tag{13.3b}
$$

这样我们就得到了两个常微分方程，同时引入了一个待定常数 $\lambda$。同样，将 (13.2) 式代入齐次边界条件，得

$$
X(0)T(t)=0,
\qquad
X(l)T(t)=0.
$$

因为 $T(t)$ 不可能恒为 0（否则 $u(x,t)$ 恒为 0），所以这时必须有

$$
X(0)=0,
\qquad
X(l)=0.
\tag{13.3c}
$$

这样，就完成了分离变量。分离变量之所以能够实现，是因为原来的偏微分方程和边界条件都是齐次的。

第二步，求解本征值问题。

分离变量的结果，是得到了多个含有待定常数的齐次常微分方程和齐次边界条件，齐次常微分方程 (13.3b) 和相应的齐次边界条件 (13.3c) 构成常微分方程的本征值问题。其特点是：常微分方程 (13.3b) 中含有一个待定常数 $\lambda$，而定解条件 (13.3c) 是一对齐次边界条件。由常微分方程的理论可知，常微分方程初值问题的解一定存在并且唯一，而且，如果微分方程和初始条件都是齐次的，就一定只有零解。在齐次边界条件下，齐次常微分方程常常也只有零解。但是，下面将看到，当待定参数 $\lambda$ 取某些特定值时，可以有既满足齐次常微分方程 (13.3b)、又满足齐次边界条件 (13.3c) 的非零解 $X(x)$。$\lambda$ 的这些特定值称为本征值，相应的非零解称为本征函数。

我们需要在复数域来求解本征值问题，即不能先验地认定本征值是实数，本征函数也是实函数。考虑到常微分方程 (13.3b) 解的特点（该方程的特征根是否为重根），应该区别 $\lambda=0$ 与 $\lambda\ne0$ 的两种情形。

① 显然应当假设特解中 $X(x)$ 和 $T(t)$ 不恒等于 0，即 $u(x,t)$ 不恒等于 0。恒等于 0 的解称为平凡（trivial）解。尽管它的确是齐次偏微分方程的解，但并不具有任何实用价值。我们要求的是非零解。

当 $\lambda=0$ 时，方程 (13.3b) 的通解是

$$
X(x)=Ax+B,
$$

代入边界条件 (13.3c)，可以定出 $A=B=0$，所以 $\lambda=0$ 不是本征值。

当 $\lambda\ne0$ 时，常微分方程 (13.3b) 的通解则是

$$
X(x)=A\sin\sqrt{\lambda}\,x+B\cos\sqrt{\lambda}\,x,
$$

代入边界条件 (13.3c)，就有

$$
B=0,
\qquad
A\sin\sqrt{\lambda}\,l=0.
$$

因为 $A\ne0$（否则 $X(x)\equiv0,\ u(x,t)\equiv0$），故必有 $\sqrt{\lambda}\,l=n\pi$，即本征值

$$
\lambda_n
=
\left(
\frac{n\pi}{l}
\right)^2,
\qquad
n=1,2,3,\cdots,
\tag{13.4a}
$$

相应的本征函数就是

$$
X_n(x)=\sin\frac{n\pi}{l}x.
\tag{13.4b}
$$

这里取 $A=1$，因为所要求的只是线性无关解①。也不必考虑 $n$ 为负整数的情形，因为 $n$ 为负整数时给出的本征值与 $-n$（正整数）时相同，相应的本征函数也线性相关。

这样我们求得了所有的本征值 $\lambda_n$ 以及相应的本征函数 $X_n(x)$，有无穷多组，它们可以用脚标 $n$ 标记。

第三步，求特解，并进一步叠加出一般解。

在求解本征值问题后，将本征值 $\lambda_n$ 代入方程 (13.3a)，求出相应的 $T_n(t)$，

$$
T_n(t)
=
C_n\sin\frac{n\pi}{l}at
+
D_n\cos\frac{n\pi}{l}at.
\tag{13.5}
$$

因此，就得到了同时满足偏微分方程 (13.1a) 和边界条件 (13.1b) 的特解

$$
u_n(x,t)
=
\left(
C_n\sin\frac{n\pi}{l}at
+
D_n\cos\frac{n\pi}{l}at
\right)
\sin\frac{n\pi}{l}x,
\qquad
n=1,2,3,\cdots.
\tag{13.6}
$$

这样的特解有无穷多个。每个特解都满足齐次偏微分方程和齐次边界条件。但一般说来，单独任何一个特解不可能也满足给定的初始条件，除非 $\phi(x)$ 和 $\psi(x)$ 恰好同时是同一个 $\sin(n\pi x/l)$ 的倍数。然而由于偏微分方程和边界条件都是线性、齐次的，把全部特解叠加起来，得到

$$
u(x,t)
=
\sum_{n=1}^{\infty}
\left(
C_n\sin\frac{n\pi}{l}at
+
D_n\cos\frac{n\pi}{l}at
\right)
\sin\frac{n\pi}{l}x,
\tag{13.7}
$$

只要级数收敛且可以逐项求二阶偏微商，那么，这样的 $u(x,t)$ 也仍然是线性齐次偏微分方程在齐次边界条件下的解。这种形式的解称为（满足齐次偏微分方程和齐次边界条件的）一般解。它不同于偏微分方程的通解，因为一般解不仅仅满足偏微分方程，而且还要满足齐次边界条件，因此并没有涵盖偏微分方程的全部解式。相应地，从解的形式看，它也不是含有两个任意函数。也不宜称为定解问题的一般解，因为对于定解问题来说，只要是适定的，它就只有唯一、确定的一个解，不再含有任何未定常数。

① 从概念上说，构成本征值问题的常微分方程和边界条件都是齐次的，因此，如果 $X(x)$ 是对应于 $\lambda$ 的本征函数，那么 $cX(x)$ 也一定是本征函数，其中 $c$ 是任意常数。因此，对应于每一个本征值 $\lambda$，就只要求出线性无关的本征函数即可。

现在遗留的问题就是，如何适当选择叠加系数 $C_n$ 和 $D_n$，使解式 (13.7) 满足初始条件

$$
\sum_{n=1}^{\infty}
D_n\sin\frac{n\pi}{l}x
=
\phi(x),
\tag{13.8}
$$

$$
\sum_{n=1}^{\infty}
C_n\frac{n\pi a}{l}
\sin\frac{n\pi}{l}x
=
\psi(x).
\tag{13.9}
$$

或者说，应当如何根据初始条件中的已知函数 $\phi(x)$ 和 $\psi(x)$ 定出叠加系数 $C_n$ 和 $D_n$。

第四步，利用本征函数的正交性定叠加系数。

这里求得的本征函数具有一个重要特性，即本征函数的正交性：对应不同本征值的本征函数正交，

$$
\int_0^l
X_n(x)X_m(x)\,dx
=
0,
\qquad
n\ne m.
\tag{13.10}
$$

设 $X_n$ 是对应本征值 $\lambda_n$ 的本征函数，即

$$
X_n''(x)+\lambda_nX_n(x)=0,
\tag{13.11a}
$$

$$
X_n(0)=0,
\qquad
X_n(l)=0;
\tag{13.11b}
$$

又设 $X_m$ 是对应本征值 $\lambda_m$ 的本征函数，即

$$
X_m''(x)+\lambda_mX_m(x)=0,
\tag{13.12a}
$$

$$
X_m(0)=0,
\qquad
X_m(l)=0.
\tag{13.12b}
$$

将方程 (13.11a) 和 (13.12a) 分别乘以 $X_m(x)$ 和 $X_n(x)$，相减，再积分，得

$$
\int_0^l
\left[
X_n''(x)X_m(x)-X_m''(x)X_n(x)
\right]dx
+
(\lambda_n-\lambda_m)
\int_0^l
X_n(x)X_m(x)\,dx
=
0.
$$

上面第一项积分可以积出，并代入边界条件 (13.11b) 式和 (13.12b) 式，得

$$
(\lambda_n-\lambda_m)
\int_0^l
X_n(x)X_m(x)\,dx
=
\left[
X_n'(x)X_m(x)-X_m'(x)X_n(x)
\right]_{x=0}^{x=l}
=
0.
$$

因此，当 $\lambda_n-\lambda_m\ne0$ 即 $\lambda_n\ne\lambda_m$ 时，一定有

$$
\int_0^l
X_n(x)X_m(x)\,dx
=
0.
$$

即对应不同本征值的本征函数正交。

现在，在 (13.8) 式两端同乘以 $\sin(m\pi x/l)$，并逐项积分，根据本征函数的正交性，就得到

$$
\begin{aligned}
\int_0^l
\phi(x)
\sin\frac{m\pi}{l}x\,dx
&=
\int_0^l
\sum_{n=1}^{\infty}
D_n
\sin\frac{n\pi}{l}x
\sin\frac{m\pi}{l}x\,dx\\
&=
\sum_{n=1}^{\infty}
D_n
\int_0^l
\sin\frac{n\pi}{l}x
\sin\frac{m\pi}{l}x\,dx\\
&=
D_m
\int_0^l
\sin^2\frac{m\pi}{l}x\,dx.
\end{aligned}
$$

上式右端的积分称为本征函数 (13.4b) 的模方①。具体计算可以得到

$$
\|X_n\|^2
\equiv
\int_0^l
|X_n(x)|^2\,dx
=
\int_0^l
\sin^2\frac{n\pi}{l}x\,dx
=
\frac l2.
\tag{13.13}
$$

因此就得到叠加系数

$$
D_n
=
\frac2l
\int_0^l
\phi(x)\sin\frac{n\pi}{l}x\,dx.
\tag{13.14}
$$

同样，可以求得

$$
C_n
=
\frac{2}{n\pi a}
\int_0^l
\psi(x)\sin\frac{n\pi}{l}x\,dx.
\tag{13.15}
$$

这样，根据已知函数 $\phi(x)$ 和 $\psi(x)$ 计算出积分，就求得了整个定解问题的解。

严格说来，上面得到的还只是形式解。对于具体问题，还必须验证：

(1) $u(x,t)$ 是否满足偏微分方程——级数解是否可以逐项求二阶偏微商；

(2) $u(x,t)$ 是否满足边界条件——级数解的和函数是否连续；

(3) 在定叠加系数时，逐项积分是否合法。

这三个问题都涉及级数解的收敛性。由于系数 $C_n$ 和 $D_n$ 是由 $\phi(x)$ 和 $\psi(x)$ 决定的，因而 $\phi(x)$ 和 $\psi(x)$ 的性质就决定了对这三个问题的回答。详细的讨论见参考书目 [1]。假定以后所讨论的问题都能够满足收敛性的要求，不再重复讨论。

从理论上说，分离变量法的成功，要取决于下列几个条件：

(1) 本征值问题有解；

(2) 定解问题的解一定可以按照本征函数展开——本征函数的全体是完备的；

(3) 本征函数一定具有正交性。

以后将在适当时候回答这几个问题。

**练习 13.1**　若定解问题中的偏微分方程和初始条件都是齐次的，仿照上面分离变量的步骤一定可以得到

$$
T''(t)+\lambda a^2T(t)=0,
$$

$$
T(0)=0,
\qquad
T'(0)=0.
$$

是否存在特定的 $\lambda$ 值，使此定解问题有非零解 $T(t)$？读者可从两种角度来回答这个问题：(1) 直接求解；(2) 运用第九章中介绍的常微分方程理论。

从这个问题的回答中，你能得出什么普遍性的结论？

① $\|X_n\|$ 的倒数常称为本征函数的**归一因子**。这是因为

$$
\frac1{\|X_n\|^2}
\int_0^l
X_n^2(x)\,dx
=
1,
$$

即本征函数 $X_n(x)/\|X_n\|$ 的模为 1。另外，(13.10) 和 (13.13) 两式还可以合并写成

$$
\int_0^l
X_n(x)X_m(x)\,dx
=
\frac l2\delta_{nm},
$$

称为本征函数的正交归一性。

**练习 13.2**　证明下列本征值问题本征函数的正交性：

$$
X''(x)+\lambda X(x)=0,
$$

$$
X(0)\cos\alpha-X'(0)\sin\alpha=0,
$$

$$
X(l)\cos\beta+X'(l)\sin\beta=0,
$$

其中 $0\le\alpha\le\pi/2,\ 0\le\beta\le\pi/2$ 都是已知实数。

现在分析一下弦的总能量。在任一时刻 $t$，有界弦的总能量是

$$
E(t)
=
\frac12
\int_0^l
\rho
\left(
\frac{\partial u}{\partial t}
\right)^2dx
+
\frac12
\int_0^l
T
\left(
\frac{\partial u}{\partial x}
\right)^2dx,
\tag{13.16}
$$

其中第一项为动能，第二项为势能。将解式 (13.7) 代入，利用本征函数的正交归一性即得

$$
E(t)
=
\frac{m\pi^2a^2}{4l^2}
\sum_{n=1}^{\infty}
n^2
\left[
|C_n|^2+|D_n|^2
\right].
\tag{13.17}
$$

等式右端显然是常数，与 $t$ 无关，即弦的总能量守恒①。

根据弦的总能量守恒，还可以证明定解问题 (13.1) 的解的唯一性。这是因为，如果此定解问题有两个解 $u_1(x,t)$ 和 $u_2(x,t)$，那么，$v(x,t)\equiv u_1(x,t)-u_2(x,t)$ 就一定满足定解问题

$$
\frac{\partial^2v}{\partial t^2}
-a^2\frac{\partial^2v}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
$$

$$
\left.v\right|_{x=0}=0,
\qquad
\left.v\right|_{x=l}=0,
\qquad
t\ge0,
$$

$$
\left.v\right|_{t=0}=0,
\qquad
\left.\frac{\partial v}{\partial t}\right|_{t=0}=0,
\qquad
0\le x\le l.
$$

只要能够证明 $v(x,t)=0$ 即可。从物理上可以判断，这肯定是正确的。从能量守恒的要求来看，当 $t=0$ 时弦的总能量为 0，因此以后任一时刻的 $E(t)$ 均为 0。这意味着一定有 $\partial v/\partial x=0,\ \partial v/\partial t=0$，即 $v(x,t)$ 为常数。由初始条件或边界条件都能定出此常数为 0。

① 更严格的办法是仿照 §13.6 的作法，直接由 (13.16) 式推出 $dE/dt=0$，而不依赖于具体的求解方法（例如，分离变量法）。这只要将 (13.16) 式对 $t$ 求导，利用方程 (13.1a) 及由边界条件 (13.1b) 导出的 $(\partial u/\partial t)_{x=0}=0$ 和 $(\partial u/\partial t)_{x=l}=0$ 即可证得。

## *§13.2 分离变量法的物理诠释

现在讨论解式 (13.7) 的物理意义。先看特解 (13.6) 式

$$
u_n(x,t)
=
\left(
C_n\sin\frac{n\pi}{l}at
+
D_n\cos\frac{n\pi}{l}at
\right)
\sin\frac{n\pi}{l}x
=
A_n\sin(\omega_nt+\delta_n)\sin k_nx,
$$

其中

$$
\omega_n=\frac{n\pi}{l}a,
\qquad
k_n=\frac{n\pi}{l},
\qquad
A_n\cos\delta_n=C_n,
\qquad
A_n\sin\delta_n=D_n.
$$

因此，$u_n(x,t)$ 代表一个驻波，$A_n\sin k_nx$ 表示弦上各点的振幅分布，$\sin(\omega_nt+\delta_n)$ 表示相位因子。$\omega_n$ 是驻波的角频率，称为两端固定弦的固有频率或本征频率，与初始条件无关；$k_n$ 称为波数，即单位长度上波的周期数；$\delta_n$ 是初相位，由初始条件决定。在 $k_nx=m\pi$，即 $x=m\pi/k_n=(m/n)l,\ m=0,1,2,3,\cdots,n$ 的各点上，振动的振幅恒为 0，称为波节。包括弦的两个端点在内，波节点共有 $n+1$ 个。在 $k_nx=(m+1/2)\pi,\ m=0,1,2,3,\cdots,n-1$ 的各点上，振动振幅的绝对值恒为最大，称为波腹。波腹点共 $n$ 个。相邻两个波节或波腹之间的距离为半波长。整个问题的解则是这些驻波的叠加。正是因为这个原因，这种解法也称为驻波法。

就两端固定的弦来说，固有频率中有一个最小值，即

$$
\omega_1=\frac{\pi}{l}a,
$$

称为基频，其他固有频率 $\omega_n$ 都是基频 $\omega_1$ 的整数倍，

$$
\omega_n=n\omega_1,
\qquad
n=2,3,\cdots,
$$

称为倍频。弦的基频就决定了所发声音的音调，在弦乐器中，当弦的材料一定（即 $\rho$ 一定）时，通过改变弦的绷紧程度（即改变张力 $T$ 的大小），就可以调节基频 $\omega_1$ 的大小。在解式 (13.7) 中，基频和倍频的叠加系数 $\{C_n\}$ 和 $\{D_n\}$ 的相对大小决定了声音的频谱分布，即决定了声音的音色。由 (13.17) 式还看到，和数

$$
\sum_{n=1}^{\infty}
n^2
\left(
|C_n|^2+|D_n|^2
\right)
$$

与弦的总能量成正比，所以就决定了声音的强度。

还可以进一步讨论分离变量法的解和行波解 (12.34) 式的联系。为此，先将初始条件 $\phi(x)$ 和 $\psi(x)$ 作奇延拓

$$
\Phi(x)
=
\begin{cases}
-\phi(-x), & -l\le x\le0,\\
\phi(x), & 0\le x\le l,
\end{cases}
\qquad
\Psi(x)
=
\begin{cases}
-\psi(-x), & -l\le x\le0,\\
\psi(x), & 0\le x\le l.
\end{cases}
$$

然后再延拓为周期为 $2l$ 的周期函数（仍记为 $\Phi(x)$ 和 $\Psi(x)$）。可以看出，这样延拓的结果保证了在端点 $x=l$ 也是奇延拓。将 $\Phi(x)$ 和 $\Psi(x)$ 展开为 Fourier 级数

$$
\Phi(x)
=
\sum_{n=1}^{\infty}
\alpha_n
\sin\frac{n\pi}{l}x,
\qquad
\Psi(x)
=
\sum_{n=1}^{\infty}
\beta_n
\sin\frac{n\pi}{l}x,
$$

其中

$$
\alpha_n
=
\frac1l
\int_{-l}^{l}
\Phi(x)\sin\frac{n\pi}{l}x\,dx
=
\frac2l
\int_0^l
\phi(x)\sin\frac{n\pi}{l}x\,dx,
$$

$$
\beta_n
=
\frac1l
\int_{-l}^{l}
\Psi(x)\sin\frac{n\pi}{l}x\,dx
=
\frac2l
\int_0^l
\psi(x)\sin\frac{n\pi}{l}x\,dx.
$$

与 (13.10)、(13.13) 两式相比较，就可以看出

$$
\alpha_n=D_n,
\qquad
\beta_n=\frac{n\pi a}{l}C_n.
$$

所以

$$
\begin{aligned}
u(x,t)
&=
\sum_{n=1}^{\infty}
\left(
\frac{\beta_nl}{n\pi a}
\sin\frac{n\pi}{l}at
+
\alpha_n
\cos\frac{n\pi}{l}at
\right)
\sin\frac{n\pi}{l}x\\
&=
\frac12
\sum_{n=1}^{\infty}
\alpha_n
\left[
\sin\frac{n\pi}{l}(x-at)
+
\sin\frac{n\pi}{l}(x+at)
\right]\\
&\quad+
\frac12
\sum_{n=1}^{\infty}
\frac{\beta_nl}{n\pi a}
\left[
\cos\frac{n\pi}{l}(x-at)
-
\cos\frac{n\pi}{l}(x+at)
\right]\\
&=
\frac12[\Phi(x-at)+\Phi(x+at)]
+
\frac1{2a}
\int_{x-at}^{x+at}
\Psi(x)\,dx.
\end{aligned}
$$

它和行波解的形式完全一致，只不过这里的 $\Phi(x)$ 和 $\Psi(x)$ 是由初始条件 $\phi(x)$ 和 $\psi(x)$ 按照前面的法则延拓而得的。这样得到的解式 $u(x,t)$ 也只适用于区间 $0\le x\le l$ 中。

从上面的讨论还可以看出波动在两端固定弦上的传播过程。为了简单起见，仍以单纯由初位移引起的波动为例。当 $t>0$ 时，初位移也像在无界弦上分别向左右传播，不同之处是到达端点 $x=0$ 或 $x=l$ 时，必须反射回来，并伴有额外的相位损失 $\pi$（即在端点 $x=0$ 和 $x=l$ 必须作奇延拓，这是由两端固定这样的边界条件决定的）。就弦上任意一点在任意一个时刻的位移而言，它就是初位移在两个端点间多次反复反射而叠加出的结果。对于初速度激发的波动，当然也可以类似地讨论。

**练习 13.3**　如果定解问题中的边界条件改为

$$
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=0}
=
0,
\qquad
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=l}
=
0,
$$

这时应该如何将初始条件作延拓？

**练习 13.4**　如果定解问题中的边界条件改为

$$
\left.u(x,t)\right|_{x=0}
=
0,
\qquad
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=l}
=
0,
$$

这时又该如何将初始条件作延拓？

## §13.3 矩形区域内的稳定问题

分离变量法也能用于求解热传导问题和稳定问题（例如 Laplace 方程）的定解问题。

设有定解问题

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
0,
\qquad
0<x<a,\ 0<y<b,
\tag{13.18a}
$$

$$
\left.u\right|_{x=0}
=
0,
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=a}
=
0,
\qquad
0\le y\le b,
\tag{13.18b}
$$

$$
\left.u\right|_{y=0}
=
f(x),
\qquad
\left.
\frac{\partial u}{\partial y}
\right|_{y=b}
=
0,
\qquad
0\le x\le a.
\tag{13.18c}
$$

仍用分离变量法求解。为此，令

$$
u(x,y)=X(x)Y(y),
\tag{13.19}
$$

代入方程 (13.18a)，分离变量，即得

$$
\frac{X''(x)}{X(x)}
=
-\frac{Y''(y)}{Y(y)}.
$$

在此等式中，左端是 $x$ 的函数（与 $y$ 无关），右端是 $y$ 的函数（与 $x$ 无关），因此，必须共同等于一个既与 $x$ 无关、又与 $y$ 无关的常数。令此常数为 $-\lambda$，上面的结果又可写成

$$
X''(x)+\lambda X(x)=0,
\tag{13.20a}
$$

$$
Y''(y)-\lambda Y(y)=0.
\tag{13.20b}
$$

同样，将解式 (13.19) 代入关于 $x$ 的一对齐次边界条件 (13.18b)，也可以分离变量得

$$
X(0)=0,
\qquad
X'(a)=0.
\tag{13.20c}
$$

这样，齐次常微分方程 (13.20a) 和齐次边界条件 (13.20c) 又构成了一个本征值问题。方程 (13.20a) 的通解是

$$
X(x)
=
\begin{cases}
A_0x+B_0, & \lambda=0,\\[4pt]
A\sin\sqrt{\lambda}\,x+B\cos\sqrt{\lambda}\,x,
& \lambda\ne0.
\end{cases}
$$

代入边界条件 (13.20c)，当 $\lambda=0$ 时，可以定出

$$
A=B=0,
$$

因此 $\lambda=0$ 不是本征值。而若 $\lambda\ne0$，则可得到

$$
B=0,
\qquad
A\ne0,
\qquad
\cos\sqrt{\lambda}\,a=0,
$$

于是就求出了

本征值

$$
\lambda_n
=
\left(
\frac{2n+1}{2a}\pi
\right)^2,
\qquad
n=0,1,2,3,\cdots;
\tag{13.21}
$$

本征函数

$$
X_n(x)
=
\sin\frac{2n+1}{2a}\pi x.
\tag{13.22}
$$

相应地，由方程 (13.20b) 可求出

$$
Y_n(y)
=
C_n\sinh\frac{2n+1}{2a}\pi y
+
D_n\cosh\frac{2n+1}{2a}\pi y.
$$

于是，就得到了既满足 Laplace 方程 (13.18a)，又满足齐次边界条件 (13.18b) 的特解

$$
u_n(x,y)
=
\left(
C_n\sinh\frac{2n+1}{2a}\pi y
+
D_n\cosh\frac{2n+1}{2a}\pi y
\right)
\sin\frac{2n+1}{2a}\pi x.
\tag{13.23}
$$

将这无穷多个特解叠加起来，就得到一般解

$$
u(x,y)
=
\sum_{n=0}^{\infty}
\left(
C_n\sinh\frac{2n+1}{2a}\pi y
+
D_n\cosh\frac{2n+1}{2a}\pi y
\right)
\sin\frac{2n+1}{2a}\pi x.
\tag{13.24}
$$

代入关于 $y$ 的一对边界条件 (13.18c)，

$$
\left.u\right|_{y=0}
=
\sum_{n=0}^{\infty}
D_n
\sin\frac{2n+1}{2a}\pi x
=
f(x),
$$

$$
\left.
\frac{\partial u}{\partial y}
\right|_{y=b}
=
\sum_{n=0}^{\infty}
\frac{2n+1}{2a}\pi
\left(
C_n\cosh\frac{2n+1}{2a}\pi b
+
D_n\sinh\frac{2n+1}{2a}\pi b
\right)
\sin\frac{2n+1}{2a}\pi x
=
0.
$$

再次根据本征函数的正交归一性，

$$
\int_0^a
\sin\frac{2n+1}{2a}\pi x
\sin\frac{2m+1}{2a}\pi x\,dx
=
\frac a2\delta_{nm},
\tag{13.25}
$$

就可以求得

$$
D_n
=
\frac2a
\int_0^a
f(x)
\sin\frac{2n+1}{2a}\pi x\,dx
\tag{13.26}
$$

和

$$
C_n\cosh\frac{2n+1}{2a}\pi b
+
D_n\sinh\frac{2n+1}{2a}\pi b
=
0,
$$

由此得

$$
C_n
=
-D_n\tanh\frac{2n+1}{2a}\pi b.
\tag{13.27}
$$

这样，就求出了矩形区域内 Laplace 方程边值问题 (13.18) 的级数解。知道了 $f(x)$ 的具体形式，就可以求出叠加系数 $C_n$ 和 $D_n$。

在这个问题中，因为与 $t$ 无关，不出现初始条件。用分离变量法求解时，则是采用关于 $x$ 的一对齐次边界条件构成本征值问题，而用另一对边界条件定叠加系数。

**练习 13.5**　将关于 $Y(y)$ 的常微分方程的通解写成别的形式，例如

$$
(1)\quad
Y_n(y)
=
A_n\exp
\left\{
\frac{2n+1}{2a}\pi y
\right\}
+
B_n\exp
\left\{
-\frac{2n+1}{2a}\pi y
\right\},
$$

$$
(2)\quad
Y_n(y)
=
A_n\sinh\frac{2n+1}{2a}\pi y
+
B_n\cosh\frac{2n+1}{2a}\pi(b-y),
$$

再重新求解上述定解问题。试问：(1) 结果有何不同？(2) 采用何种形式的 $Y_n(y)$ 计算更为便捷？

**练习 13.6**　如果定解问题的边界条件全是非齐次的，例如

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
0,
\qquad
0<x<a,\ 0<y<b,
$$

$$
\left.u\right|_{x=0}
=
f(y),
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=a}
=
g(y),
\qquad
0\le y\le b,
$$

$$
\left.u\right|_{y=0}
=
\phi(x),
\qquad
\left.
\frac{\partial u}{\partial y}
\right|_{y=b}
=
\psi(x),
\qquad
0\le x\le a,
$$

这时应该如何利用分离变量法求解？


## §13.4 多于两个自变量的定解问题

当定解问题的自变量超过两个时（当然时间变量最多只能一个），仍可采用分离变量法求解，但有一些技术性的细节需要注意。

作为示例，下面扼要讨论一下矩形介质的热传导问题。假设介质四周绝热，定解问题是

$$
\frac{\partial u}{\partial t}
-
\kappa
\left(
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
\right)
=
0,
\qquad
0<x<a,\ 0<y<b,\ t>0,
\tag{13.28a}
$$

$$
\left.
\frac{\partial u}{\partial x}
\right|_{x=0}
=
0,
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=a}
=
0,
\qquad
0\le y\le b,\ t\ge0,
\tag{13.28b}
$$

$$
\left.
\frac{\partial u}{\partial y}
\right|_{y=0}
=
0,
\qquad
\left.
\frac{\partial u}{\partial y}
\right|_{y=b}
=
0,
\qquad
0\le x\le a,\ t\ge0,
\tag{13.28c}
$$

$$
\left.u\right|_{t=0}
=
\phi(x,y),
\qquad
0\le x\le a,\ 0\le y\le b.
\tag{13.28d}
$$

设 $u(x,y,t)=X(x)Y(y)T(t)$，代入齐次方程 (13.28a)，即得

$$
\frac{X''(x)}{X(x)}
+
\frac{Y''(y)}{Y(y)}
-
\frac1\kappa
\frac{T'(t)}{T(t)}
=
0.
\tag{13.29}
$$

此等式成立，意味着式中各项均应为常数，即

$$
X''(x)+\mu X(x)=0,
\tag{13.30a}
$$

$$
Y''(y)+\nu Y(y)=0,
\tag{13.31a}
$$

$$
T'(t)+\lambda\kappa T(t)=0,
\tag{13.32}
$$

其中 $\mu,\nu,\lambda$ 是分离变量时引进的待定常数，$\mu+\nu=\lambda$。

再将齐次边界条件 (13.28b)、(13.28c) 分离变量，又可以得到

$$
X'(0)=0,
\qquad
X'(a)=0,
\tag{13.30b}
$$

$$
Y'(0)=0,
\qquad
Y'(b)=0.
\tag{13.31b}
$$

(13.30a)、(13.30b) 和 (13.31a)、(13.31b) 就分别构成了两个本征值问题。

现在着手解本征值问题 (13.30)。当 $\mu=0$ 时，方程 (13.30a) 的通解为 $X(x)=Ax+B$，由边界条件 (13.30b)，可以求得 $A=0,\ B$ 任意。因此，对于本征值问题 (13.30)，

$$
\mu_0=0
\tag{13.33}
$$

是本征值，相应的本征函数可以取为

$$
X_0(x)=1.
\tag{13.34}
$$

当 $\mu\ne0$ 时，方程 (13.30a) 的通解为 $X(x)=A\sin\sqrt{\mu}\,x+B\cos\sqrt{\mu}\,x$，代入边界条件 (13.30b)，有 $A=0,\ B\ne0,\ \sin\sqrt{\mu}\,a=0$，因而可以求得 $\sqrt{\mu}\,a=n\pi$，即

$$
\mu_n
=
\left(
\frac{n\pi}{a}
\right)^2,
\qquad
n=1,2,3,\cdots.
\tag{13.35}
$$

相应的本征函数是（取 $B=1$）

$$
X_n(x)
=
\cos\frac{n\pi}{a}x.
\tag{13.36}
$$

把 $\mu=0$ 和 $\mu>0$ 的结果合并起来，就可以统一写成

本征值

$$
\mu_n
=
\left(
\frac{n\pi}{a}
\right)^2,
\qquad
n=0,1,2,3,\cdots,
\tag{13.37a}
$$

本征函数

$$
X_n(x)
=
\cos\frac{n\pi}{a}x.
\tag{13.37b}
$$

同样可以解得本征值问题 (13.31) 的解为

本征值

$$
\nu_m
=
\left(
\frac{m\pi}{b}
\right)^2,
\qquad
m=0,1,2,3,\cdots;
\tag{13.38a}
$$

本征函数

$$
Y_m(x)
=
\cos\frac{m\pi}{b}y.
\tag{13.38b}
$$

注意这里的 $n$ 和 $m$ 互相独立。对于给定的 $n$ 和 $m$，再进一步求出方程 (13.32) 的通解

$$
T_{nm}(t)
=
A_{nm}e^{-\lambda_{nm}\kappa t},
\qquad
n=0,1,2,\cdots,\quad
m=0,1,2,\cdots,
\tag{13.39}
$$

其中

$$
\lambda_{nm}
=
\mu_n+\nu_m
=
\left(
\frac{n\pi}{a}
\right)^2
+
\left(
\frac{m\pi}{b}
\right)^2.
$$

因此，就求得同时满足方程 (13.28a) 和边界条件 (13.28b)、(13.28c) 的特解

$$
u_{nm}(x,y,t)
=
X_n(x)Y_m(y)T_{nm}(t)
=
A_{nm}
\cos\frac{n\pi}{a}x
\cos\frac{m\pi}{b}y
e^{-\lambda_{nm}\kappa t},
\tag{13.40}
$$

并由此进一步叠加得一般解

$$
u(x,y,t)
=
\sum_{n=0}^{\infty}
\sum_{m=0}^{\infty}
A_{nm}
\cos\frac{n\pi}{a}x
\cos\frac{m\pi}{b}y
\exp
\left\{
-
\left[
\left(
\frac{n\pi}{a}
\right)^2
+
\left(
\frac{m\pi}{b}
\right)^2
\right]
\kappa t
\right\}.
\tag{13.41}
$$

代入初始条件 (13.28d)，有

$$
\left.u(x,y,t)\right|_{t=0}
=
\sum_{n=0}^{\infty}
\sum_{m=0}^{\infty}
A_{nm}
\cos\frac{n\pi}{a}x
\cos\frac{m\pi}{b}y
=
\phi(x,y).
\tag{13.42}
$$

下一步应当根据本征函数的正交性定出叠加系数。需要注意，现在既要用到 $\{X_n(x),n=0,1,2,\cdots\}$ 的正交性，又要用到 $\{Y_m(y),m=0,1,2,\cdots\}$ 的正交性，缺一不可。其次，考虑到它们的正交归一性

$$
\int_0^a
X_n(x)X_{n'}(x)\,dx
=
\frac a2
(1+\delta_{n0})
\delta_{nn'},
\tag{13.43}
$$

$$
\int_0^b
Y_m(y)Y_{m'}(y)\,dy
=
\frac b2
(1+\delta_{m0})
\delta_{mm'},
\tag{13.44}
$$

需要留心区分 $n=0$ 与 $n\ne0$ 和 $m=0$ 与 $m\ne0$ 的情形。最后的结果是

$$
A_{nm}
=
\frac4{ab}
\frac1{(1+\delta_{n0})(1+\delta_{m0})}
\int_0^a
\int_0^b
\phi(x,y)
\cos\frac{n\pi}{a}x
\cos\frac{m\pi}{b}y
\,dx\,dy.
\tag{13.45}
$$

## §13.5 两端固定弦的受迫振动

在前几节的讨论中，特别强调了齐次偏微分方程和齐次边界条件在分离变量法中的关键作用：因为方程和边界条件是齐次的，分离变量才得以实现。如果方程和边界条件不是齐次的，分离变量法的基本原则仍能适用，但实际求解步骤需做适当调整。本节就先讨论非齐次方程、齐次边界条件的情形。

为了突出对于方程非齐次项的处理，不妨研究纯粹由外力引起的两端固定弦的受迫振动，弦的初位移和初速度均为 0。这样，定解问题就是

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
f(x,t),
\qquad
0<x<l,\ t>0,
\tag{13.46a}
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.46b}
$$

$$
\left.u\right|_{t=0}=0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l.
\tag{13.46c}
$$

对于这种非齐次方程、齐次边界条件构成的定解问题，可以有两种处理方法。

### 1. 边界条件保持齐次，而将方程齐次化

这种解法的中心思想是把求解非齐次方程、齐次边界条件构成的定解问题转化为求解齐次方程、齐次边界条件构成的定解问题。设

$$
u(x,t)=v(x,t)+w(x,t),
\tag{13.47}
$$

如果我们能够求出满足

$$
\frac{\partial^2v}{\partial t^2}
-
a^2\frac{\partial^2v}{\partial x^2}
=
f(x,t),
\tag{13.48a}
$$

$$
\left.v(x,t)\right|_{x=0}=0,
\qquad
\left.v(x,t)\right|_{x=l}=0
\tag{13.48b}
$$

的函数 $v(x,t)$，换言之，也就是能够求得既满足非齐次方程 (13.46a)，而且还满足齐次边界条件 (13.46b) 的特解 $v(x,t)$，则 $w(x,t)$ 一定满足

$$
\frac{\partial^2w}{\partial t^2}
-
a^2\frac{\partial^2w}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.49a}
$$

$$
\left.w\right|_{x=0}=0,
\qquad
\left.w\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.49b}
$$

$$
\left.w\right|_{t=0}
=
-\left.v\right|_{t=0},
\qquad
\left.
\frac{\partial w}{\partial t}
\right|_{t=0}
=
-
\left.
\frac{\partial v}{\partial t}
\right|_{t=0},
\qquad
0\le x\le l.
\tag{13.49c}
$$

因此，我们可以重复 §13.1 中的步骤，求出满足方程 (13.49a) 和边界条件 (13.49b) 的特解，并进一步叠加出一般解，从而得到

$$
u(x,t)
=
v(x,t)
+
\sum_{n=1}^{\infty}
\left(
C_n\sin\frac{n\pi}{l}at
+
D_n\cos\frac{n\pi}{l}at
\right)
\sin\frac{n\pi}{l}x.
\tag{13.50}
$$

再代入初始条件 (13.49c)，有

$$
\sum_{n=1}^{\infty}
D_n\sin\frac{n\pi}{l}x
=
-\left.v(x,t)\right|_{t=0},
$$

$$
\sum_{n=1}^{\infty}
C_n\frac{n\pi a}{l}
\sin\frac{n\pi}{l}x
=
-
\left.
\frac{\partial v(x,t)}{\partial t}
\right|_{t=0}.
$$

利用本征函数的正交归一性，就可以定出叠加系数

$$
C_n
=
-\frac2{n\pi a}
\int_0^l
\left.
\frac{\partial v(x,t)}{\partial t}
\right|_{t=0}
\sin\frac{n\pi}{l}x\,dx,
\tag{13.51}
$$

$$
D_n
=
-\frac2l
\int_0^l
v(x,0)\sin\frac{n\pi}{l}x\,dx.
\tag{13.52}
$$

这种解法可以简单地称为方程齐次化法。需要特别注意，在将非齐次方程齐次化的同时，必须保持原有的齐次边界条件不变。这种解法的关键就在于求得特解 $v(x,t)$。如果方程的非齐次项 $f(x,t)$ 的形式比较简单，可以尝试采用这种解法。

从求解过程可以看出，初始条件是否齐次不影响求解步骤。

**例 13.1**　求解定解问题

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
f(x),
\qquad
0<x<l,\ t>0,
\tag{13.53a}
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.53b}
$$

$$
\left.u\right|_{t=0}=0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l,
\tag{13.53c}
$$

其中 $f(x)$ 为已知函数。

**解**　这里只给出解题的主要思路。由于方程的非齐次项只是 $x$ 的函数，不妨把齐次化函数也取为只是 $x$ 的函数，即设

$$
u(x,t)=v(x)+w(x,t),
\tag{13.54}
$$

其中 $v(x)$ 可以由常微分方程的边值问题

$$
v''(x)
=
-\frac1{a^2}f(x),
\tag{13.55a}
$$

$$
v(0)=0,
\qquad
v(l)=0
\tag{13.55b}
$$

求出，而 $w(x,t)$ 则满足齐次方程、齐次边界条件型的定解问题

$$
\frac{\partial^2w}{\partial t^2}
-
a^2\frac{\partial^2w}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.56a}
$$

$$
\left.w\right|_{x=0}=0,
\qquad
\left.w\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.56b}
$$

$$
\left.w\right|_{t=0}
=
-v(x),
\qquad
\left.
\frac{\partial w}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l,
\tag{13.56c}
$$

可以用 §13.1 的方法求解。

**例 13.2**　求解定解问题

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
A_0\sin\omega t,
\qquad
0<x<l,\ t>0,
\tag{13.57a}
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.57b}
$$

$$
\left.u\right|_{t=0}=0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l,
\tag{13.57c}
$$

其中 $a,A_0$ 及 $\omega$ 均为已知常数。

**解**　设

$$
u(x,t)=v(x,t)+w(x,t),
\tag{13.58}
$$

其中的齐次化函数 $v(x,t)$ 取为

$$
v(x,t)=f(x)\sin\omega t.
\tag{13.59}
$$

选择 $f(x)$，使得 $v(x,t)$ 满足非齐次方程 (13.57a) 及齐次边界条件 (13.57b)，即要求 $f(x)$ 满足

$$
-\omega^2f(x)-a^2f''(x)=A_0,
\tag{13.60a}
$$

$$
f(0)=0,
\qquad
f(l)=0.
\tag{13.60b}
$$

非齐次常微分方程 (13.60a) 的通解为

$$
f(x)
=
-\frac{A_0}{\omega^2}
+
A\sin\frac{\omega}{a}x
+
B\cos\frac{\omega}{a}x.
$$

利用齐次边界条件 (13.60b) 可以定出

$$
B=\frac{A_0}{\omega^2},
\qquad
A=\frac{A_0}{\omega^2}\tan\frac{\omega l}{2a}.
$$

于是

$$
\begin{aligned}
f(x)
&=
-\frac{A_0}{\omega^2}
\left[
\left(
1-\cos\frac{\omega}{a}x
\right)
-
\tan\frac{\omega l}{2a}
\sin\frac{\omega}{a}x
\right]\\
&=
-\frac{A_0}{\omega^2}
\left[
1-
\frac{
\cos[\omega(x-l/2)/a]
}{
\cos(\omega l/2a)
}
\right].
\end{aligned}
\tag{13.61}
$$

这样就能导出 $w(x,t)$ 所满足的齐次方程、齐次边界条件类型的定解问题

$$
\frac{\partial^2w}{\partial t^2}
-
a^2\frac{\partial^2w}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.62a}
$$

$$
\left.w\right|_{x=0}=0,
\qquad
\left.w\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.62b}
$$

$$
\left.w\right|_{t=0}=0,
\qquad
\left.
\frac{\partial w}{\partial t}
\right|_{t=0}
=
-\omega f(x),
\qquad
0\le x\le l.
\tag{13.62c}
$$

容易求出满足齐次方程 (13.62a) 及齐次边界条件 (13.62b) 的一般解为

$$
w(x,t)
=
\sum_{n=1}^{\infty}
\left(
C_n\sin\frac{n\pi}{l}at
+
D_n\cos\frac{n\pi}{l}at
\right)
\sin\frac{n\pi}{l}x.
\tag{13.63}
$$

利用上面的初始条件 (13.62c) 就可以定出

$$
D_n=0,
\tag{13.64a}
$$

$$
\begin{aligned}
C_n
&=
-\frac{2\omega}{n\pi a}
\int_0^l
f(x)
\sin\frac{n\pi}{l}x\,dx\\
&=
-\frac{2A_0\omega l^3}{\pi^2a}
\frac{1-(-)^n}{n^2}
\frac1{(n\pi a)^2-(\omega l)^2}.
\end{aligned}
\tag{13.64b}
$$

这说明，只有 $n=$ 奇数时，$C_n$ 才不为 0。这样，最后就求出了

$$
w(x,t)
=
-\frac{4A_0\omega l^3}{\pi^2a}
\sum_{n=0}^{\infty}
\frac1{(2n+1)^2}
\frac1{[(2n+1)\pi a]^2-(\omega l)^2}
\sin\frac{2n+1}{l}\pi x
\sin\frac{2n+1}{l}\pi at
\tag{13.65}
$$

和

$$
\begin{aligned}
u(x,t)
={}&
-\frac{A_0}{\omega^2}
\left[
1-
\frac{\cos[\omega(x-l/2)/a]}
{\cos(\omega l/2a)}
\right]
\sin\omega t\\
&-
\frac{4A_0\omega l^3}{\pi^2a}
\sum_{n=0}^{\infty}
\frac1{(2n+1)^2}
\frac1{[(2n+1)\pi a]^2-(\omega l)^2}
\sin\frac{2n+1}{l}\pi x
\sin\frac{2n+1}{l}\pi at.
\end{aligned}
\tag{13.66}
$$

在上面的解题过程中，还忽略了一个特殊情形，即驱动力的角频率 $\omega$ 正好是弦的某个固有频率，例如 $\omega_{2k+1}=(2k+1)\pi a/l$，$k$ 为某个确定的自然数时，弦在驱动力的作用下会发生共振现象。这表现在求解过程中，作为常微分方程边值问题 (13.60) 的解，(13.61) 式失去意义，此后的解题过程均需作相应的修改。事实上，最后的解应该就是 (13.66) 式在 $\omega\to\omega_{2k+1}$ 下的极限值。最简单的作法是将 (13.61) 式中的 $f(x)$ 也按本征函数组 $\{\sin(n\pi x/l),n=1,2,3,\cdots\}$ 展开（其实在 (13.64b) 式中已经完成此项计算），

$$
f(x)
=
-\sum_{n=1}^{\infty}
\frac{n\pi a}{\omega l}
C_n
\sin\frac{n\pi}{l}x,
$$

其中的 $C_n$ 见 (13.64b) 式。整理即得

$$
\begin{aligned}
u(x,t)
={}&
\frac{4A_0l^2}{\pi^2a}
\sum_{n=0}^{\infty}
\frac1{(2n+1)^2}
\frac1{[(2n+1)\pi a]^2-(\omega l)^2}
\sin\frac{2n+1}{l}\pi x\\
&\times
\left[
(2n+1)\pi a\sin\omega t
-
(\omega l)
\sin\frac{2n+1}{l}\pi at
\right].
\end{aligned}
\tag{13.67}
$$

当 $\omega\to\omega_{2k+1}$ 时，此和式中的 $n=k$ 项为不定式，应单独提出，而用 l’Hôpital 法则求出极限值，

$$
\begin{aligned}
u(x,t)
={}&
\frac{4A_0l^2}{\pi^2a}
\sum_{n=0}^{\infty}{}'
\left\{
\frac1{(2n+1)^2}
\frac1{[(2n+1)\pi a]^2-(\omega l)^2}
\sin\frac{2n+1}{l}\pi x\right.\\
&\left.\qquad\qquad\times
\left[
(2n+1)\pi a\sin\omega t
-
(\omega l)
\sin\frac{2n+1}{l}\pi at
\right]
\right\}\\
&-
\frac{2A_0l}{\pi^2a}
\frac1{(2k+1)^2}
\sin\frac{2k+1}{l}\pi x
\left[
t\cos\frac{2k+1}{l}\pi at
-
\frac{l}{(2k+1)\pi a}
\sin\frac{2k+1}{l}\pi at
\right],
\end{aligned}
\tag{13.68}
$$

其中 $\sum'$ 表示和式中不含 $n=k$ 项。

**练习 13.7**　如果定解问题中的方程和初始条件都是非齐次的，但边界条件仍是齐次的，例如

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
f(x,t),
\qquad
0<x<l,\ t>0,
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0,
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
\psi(x),
\qquad
0\le x\le l,
$$

重复上面的办法，利用分离变量法求解。

**练习 13.8**　如果定解问题中的方程和边界条件都是非齐次的，例如

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
f(x,t),
\qquad
0<x<l,\ t>0,
$$

$$
\left.u\right|_{x=0}
=
\mu(t),
\qquad
\left.u\right|_{x=l}
=
\nu(t),
\qquad
t\ge0,
$$

$$
\left.u\right|_{t=0}
=
0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l,
$$

这时是否仍然可以仿照上面的方法求解？

### 2. 按相应齐次问题本征函数展开

如果方程非齐次项 $f(x,t)$ 的形式比较复杂，难以求得非齐次方程的特解，就可以采用按相应齐次问题本征函数展开的解法。其中心思想是寻找一组本征函数 $\{X_n(x),n=1,2,3,\cdots\}$，只要这组本征函数是完备的，就可以将所要求的解 $u(x,t)$ 及方程的非齐次项 $f(x,t)$ 均按这组本征函数展开：

$$
u(x,t)
=
\sum_{n=1}^{\infty}
T_n(t)X_n(x),
\tag{13.69}
$$

$$
f(x,t)
=
\sum_{n=1}^{\infty}
g_n(t)X_n(x)
\tag{13.70}
$$

然后再求出 $T_n(t)$ 即可。由于 $T_n(t)$ 是一元函数，它满足的是常微分方程（组），有可能比求解偏微分方程来得简单。

仍以定解问题 (13.46) 为例。本征函数组 $\{X_n(x)\}$ 必须满足齐次边界条件 (13.46b)，即

$$
X_n(0)=0,
\qquad
X_n(l)=0.
\tag{13.71}
$$

对于 $\{X_n(x)\}$ 所满足的微分方程，原则上没有什么限制。但从实际可行性来看，最佳选择是取 $\{X_n(x)\}$ 为相应齐次定解问题的本征函数，即满足由相应齐次偏微分方程和齐次边界条件

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0
$$

分离变量而得到的本征值问题

$$
X_n''(x)+\lambda_nX_n(x)=0,
$$

$$
X_n(0)=0,
\qquad
X_n(l)=0.
$$

将 $u(x,t)$ 和 $f(x,t)$ 的展开式 (13.69) 和 (13.70) 代入偏微分方程 (13.46a)，并逐项微商，有

$$
\sum_{n=1}^{\infty}
T_n''(t)X_n(x)
-
a^2
\sum_{n=1}^{\infty}
T_n(t)X_n''(x)
=
\sum_{n=1}^{\infty}
g_n(t)X_n(x).
$$

利用 $X_n(x)$ 所满足的常微分方程，又可以化成

$$
\sum_{n=1}^{\infty}
T_n''(t)X_n(x)
+
a^2
\sum_{n=1}^{\infty}
\lambda_nT_n(t)X_n(x)
=
\sum_{n=1}^{\infty}
g_n(t)X_n(x).
$$

再根据本征函数的正交性，就得到 $T_n(t)$ 所满足的常微分方程

$$
T_n''(t)
+
\lambda_na^2T_n(t)
=
g_n(t).
\tag{13.72}
$$

同样将 $u(x,t)$ 的展开式 (13.69) 代入初始条件 (13.46c)，也可得到

$$
\sum_{n=1}^{\infty}
T_n(0)X_n(x)
=
0,
\qquad
\sum_{n=1}^{\infty}
T_n'(0)X_n(x)
=
0.
$$

根据本征函数的正交性，即能导出

$$
T_n(0)=0,
\qquad
T_n'(0)=0.
\tag{13.73}
$$

用解非齐次常微分方程的常数变易法，或者用 Laplace 变换，就可以最后求出非齐次方程 (13.72) 在初始条件 (13.73) 下的解

$$
T_n(t)
=
\frac{l}{n\pi a}
\int_0^t
g_n(\tau)
\sin\frac{n\pi}{l}a(t-\tau)\,d\tau.
\tag{13.74}
$$

**练习 13.9**　为什么本征函数组 $\{X_n(x)\}$ 必须满足定解问题的齐次边界条件 (13.46b) 或 (13.71)？这个限制是否可以取消？

也还可以用这种方法求解例 13.2 中的定解问题。相应齐次问题的本征函数已经在 §13.1 中给出，因此可以将 $u(x,t)$ 及方程的非齐次项 $A_0\sin\omega t$ 也按这一组本征函数展开：

$$
u(x,t)
=
\sum_{n=1}^{\infty}
T_n(t)\sin\frac{n\pi}{l}x,
\tag{13.75}
$$

$$
A_0\sin\omega t
=
\frac{2A_0}{\pi}
\sum_{n=1}^{\infty}
\frac{1-(-1)^n}{n}
\sin\frac{n\pi}{l}x
\sin\omega t.
\tag{13.76}
$$

代入方程 (13.57a) 和初始条件 (13.57c)，就得到

$$
T_n''(t)
+
\left(
\frac{n\pi a}{l}
\right)^2
T_n(t)
=
\frac{2A_0}{\pi}
\frac{1-(-1)^n}{n}
\sin\omega t,
$$

$$
T_n(0)=0,
\qquad
T_n'(0)=0.
$$

解之即得

$$
T_n(t)
=
\frac{2A_0l^2}{n\pi}
\frac{1-(-1)^n}{(n\pi a)^2-(\omega l)^2}
\sin\omega t
-
\frac{2A_0\omega l^3}{n^2\pi^2a}
\frac{1-(-1)^n}{(n\pi a)^2-(\omega l)^2}
\sin\frac{n\pi}{l}at.
\tag{13.77}
$$

代入 (13.75) 式，就求出了例 13.2 的解。不难看出，这正是 (13.67) 式。


## §13.6 非齐次边界条件的齐次化

到目前为止，除了在稳定问题中需要有一部分边界条件用于定叠加系数，因而允许是非齐次的以外，在波动问题或热传导问题中，我们总是要求边界条件是齐次的。上面列举出各种技术性的理由来说明为什么边界条件必须是齐次的，例如说，非齐次边界条件不能分离变量，又说，只有满足齐次方程和齐次边界条件的特解叠加起来才仍能满足齐次方程和齐次边界条件。但最根本的原因涉及本征函数的完备性，我们将在第十七章中讨论这个问题。

现在讨论非齐次边界条件的处理。为了叙述方便，下面仍以波动方程的定解问题为例。为了突出非齐次边界条件的处理，假定方程和初始条件都是齐次的。

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.78a}
$$

$$
\left.u\right|_{x=0}
=
\mu(t),
\qquad
\left.u\right|_{x=l}
=
\nu(t),
\qquad
t\ge0,
\tag{13.78b}
$$

$$
\left.u\right|_{t=0}=0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l.
\tag{13.78c}
$$

这时为了应用分离变量法，别无选择，只有先将非齐次边界条件齐次化，即令

$$
u(x,t)=v(x,t)+w(x,t),
\tag{13.79}
$$

适当选择 $v(x,t)$，使之满足和未知函数 $u(x,t)$ 的非齐次边界条件 (13.78b) 相同的非齐次边界条件

$$
\left.v(x,t)\right|_{x=0}
=
\mu(t),
\qquad
\left.v(x,t)\right|_{x=l}
=
\nu(t).
\tag{13.80}
$$

这样，$w(x,t)$ 当然就一定满足齐次边界条件

$$
\left.w(x,t)\right|_{x=0}=0,
\qquad
\left.w(x,t)\right|_{x=l}=0.
\tag{13.81}
$$

一般说来，$w(x,t)$ 所满足的方程和初始条件都将是非齐次的，

$$
\frac{\partial^2w}{\partial t^2}
-
a^2\frac{\partial^2w}{\partial x^2}
=
-
\left(
\frac{\partial^2v}{\partial t^2}
-
a^2\frac{\partial^2v}{\partial x^2}
\right),
$$

$$
\left.w\right|_{t=0}
=
-\left.v\right|_{t=0},
\qquad
\left.
\frac{\partial w}{\partial t}
\right|_{t=0}
=
-
\left.
\frac{\partial v}{\partial t}
\right|_{t=0}.
$$

采用 §13.5 中的办法，求出 $w(x,t)$，再代回到 (13.79) 式中，就给出了定解问题 (13.78) 的解 $u(x,t)$。

齐次化函数 $v(x,t)$ 有相当大的选择余地，因为仅要求它满足 (13.80) 式。如果把 $t$ 看成参数，这就只要求在 $x-y$ 平面上的曲线 $y=v(x,t)$ 通过给定的两点 $(0,\mu(t))$ 和 $(l,\nu(t))$ 即可。这样的曲线当然有无穷多条。例如，可取直线

$$
v(x,t)=A(t)x+B(t),
$$

由 (13.80) 式即可定出

$$
B(t)=\mu(t),
\qquad
A(t)=\frac1l[\nu(t)-\mu(t)].
$$

也可取抛物线

$$
v(x,t)=A(t)x^2+B(t),
$$

$$
A(t)=\frac1{l^2}[\nu(t)-\mu(t)],
\qquad
B(t)=\mu(t)
$$

或

$$
v(x,t)
=
A(t)(l-x)^2+B(t)x^2,
$$

$$
A(t)=\frac1{l^2}\mu(t),
\qquad
B(t)=\frac1{l^2}\nu(t).
$$

**练习 13.10**　如果将定解问题 (13.78) 中的边界条件改为第二类，即

$$
\left.
\frac{\partial u}{\partial x}
\right|_{x=0}
=
\mu(t),
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
\nu(t).
$$

试问：这时的齐次化函数 $v(x,t)$ 应如何选取？

**练习 13.11**　如果将定解问题中的边界条件是混合型的，例如

$$
\left.u\right|_{x=0}
=
\mu(t),
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
\nu(t).
$$

这时的齐次化函数 $v(x,t)$ 又该如何选取？

**例 13.3**　求解定解问题

$$
\frac{\partial u}{\partial t}
-
\kappa\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.82a}
$$

$$
\left.u\right|_{x=0}
=
A\sin\omega t,
\qquad
\left.u\right|_{x=l}
=
0,
\qquad
t\ge0,
\tag{13.82b}
$$

$$
\left.u\right|_{t=0}
=
0,
\qquad
0\le x\le l.
\tag{13.82c}
$$

**解**　考虑到非齐次边界条件的具体形式，可设齐次化函数 $v(x,t)=A(1-x/l)\sin\omega t$，即令

$$
u(x,t)
=
A\left(
1-\frac xl
\right)
\sin\omega t
+
w(x,t),
\tag{13.83}
$$

则 $w(x,t)$ 满足定解问题

$$
\frac{\partial w}{\partial t}
-
\kappa\frac{\partial^2w}{\partial x^2}
=
-A\omega
\left(
1-\frac xl
\right)
\cos\omega t,
\qquad
0<x<l,\ t>0,
\tag{13.84a}
$$

$$
\left.w\right|_{x=0}=0,
\qquad
\left.w\right|_{x=l}=0,
\qquad
t\ge0,
\tag{13.84b}
$$

$$
\left.w\right|_{t=0}
=
0,
\qquad
0\le x\le l.
\tag{13.84c}
$$

将 $w(x,t)$ 和方程的非齐次项 $1-x/l$ 都按相应齐次问题的本征函数展开，有

$$
w(x,t)
=
\sum_{n=1}^{\infty}
T_n(t)
\sin\frac{n\pi}{l}x,
\qquad
1-\frac xl
=
\sum_{n=1}^{\infty}
\frac2{n\pi}
\sin\frac{n\pi}{l}x.
\tag{13.85}
$$

根据 $T_n(t)$ 应该满足非齐次一阶常微分方程

$$
T_n'(t)
+
\kappa
\left(
\frac{n\pi}{l}
\right)^2
T_n(t)
=
-\frac{2A\omega}{n\pi}
\cos\omega t
$$

和初始条件 $T_n(0)=0$，容易求出

$$
T_n(t)
=
\frac{2A\omega l^2}
{\kappa^2(n\pi)^4+\omega^2l^4}
\frac1{n\pi}
\left[
\kappa(n\pi)^2e^{-(n\pi/l)^2\kappa t}
-
\kappa(n\pi)^2\cos\omega t
-
\omega l^2\sin\omega t
\right].
\tag{13.86}
$$

这样就求得了 $w(x,t)$。代回 (13.83) 式，就得到定解问题的解 $u(x,t)$。

选择不同的齐次化函数 $v(x,t)$，当然导出的 $w(x,t)$ 的定解问题也就不同，于是求出的 $w(x,t)$ 也就不同。但是，由于定解问题解的存在唯一性，就保证了最后给出的 $u(x,t)$ 一定是相同的，尽管表达式的形式可能有所不同。

在某些特殊情形下，可以选择合适的齐次化函数 $v(x,t)$，使得 $w(x,t)$ 所满足的定解问题容易求解。最理想的就是不论原来 $u(x,t)$ 的方程是否齐次，最终 $w(x,t)$ 满足的方程是齐次的①。我们把这种方法称为将方程和边界条件同时齐次化。

**例 13.4**　求解定解问题

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.87a}
$$

$$
\left.u\right|_{x=0}=0,
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
A\sin\omega t,
\qquad
t\ge0,
\tag{13.87b}
$$

$$
\left.u\right|_{t=0}=0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0,
\qquad
0\le x\le l.
\tag{13.87c}
$$

**解**　现在就试图找到齐次化函数，将方程和边界条件同时齐次化。为此，设 $u(x,t)=v(x,t)+w(x,t)$，考虑到非齐次边界条件的具体函数形式，可取齐次化函数 $v(x,t)$ 为

$$
v(x,t)=f(x)\sin\omega t,
\tag{13.88}
$$

且 $f(x)$ 是下列常微分方程边值问题

$$
f''(x)
+
\left(
\frac{\omega}{a}
\right)^2
f(x)
=
0,
$$

$$
f(0)=0,
\qquad
f'(l)=A
$$

的解。容易求得

$$
f(x)
=
\frac{Aa}{\omega}
\frac1{\cos(\omega l/a)}
\sin\frac{\omega}{a}x.
\tag{13.89}
$$

于是，就得到 $w(x,t)$ 所满足的齐次方程、齐次边界条件类型的定解问题

$$
\frac{\partial^2w}{\partial t^2}
-
a^2\frac{\partial^2w}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
\tag{13.90a}
$$

$$
\left.w\right|_{x=0}
=
0,
\qquad
\left.
\frac{\partial w}{\partial x}
\right|_{x=l}
=
0,
\qquad
t\ge0,
\tag{13.90b}
$$

$$
\left.w\right|_{t=0}
=
0,
\qquad
\left.
\frac{\partial w}{\partial t}
\right|_{t=0}
=
-\frac{Aa}{\cos(\omega l/a)}
\sin\frac{\omega}{a}x,
\qquad
0\le x\le l.
\tag{13.90c}
$$

利用 §13.1 的方法解得

$$
w(x,t)
=
\sum_{n=0}^{\infty}
\left(
C_n
\sin\frac{2n+1}{2l}\pi at
+
D_n
\cos\frac{2n+1}{2l}\pi at
\right)
\sin\frac{2n+1}{2l}\pi x,
\tag{13.91}
$$

其中

$$
\begin{aligned}
C_n
&=
-\frac{4A}{\pi\cos(\omega l/a)}
\frac1{2n+1}
\int_0^l
\sin\frac{\omega}{a}x
\sin\frac{2n+1}{2l}\pi x\,dx\\
&=
(-)^n
\frac{16A\omega l^2a}{(2n+1)\pi}
\frac1{(\omega l)^2-[(2n+1)\pi a]^2},
\end{aligned}
\tag{13.92}
$$

$$
D_n=0.
\tag{13.93}
$$

将 $v(x,t)$ 和 $w(x,t)$ 代回到 (13.88) 式中，就最后给出了定解问题 (13.87) 的解 $u(x,t)$。

在这个例子中，之所以能够容易地找到齐次化函数 $v(x,t)$，使方程和边界条件同时齐次化，当然是由边界条件中非齐次项的特殊形式决定的：非齐次项的函数形式是 $\sin\omega t$，而方程正好有 $(A\sin kx+B\cos kx)(C\sin kat+D\cos kat)$ 形式的特解。

现在，对于具有非齐次边界条件的波动问题或热传导问题，就有两种解法可供选择。到底只是简单地将边界条件齐次化，还是力求使方程也同时齐次化，需要视具体条件而定。有时尽管也能使方程和边界条件同时齐次化，但如果齐次化函数 $v(x,t)$ 的形式过于复杂，不易求得，而导出的 $w(x,t)$ 的定解问题又比较复杂，也许还是找一个形式比较简单的函数，只将边界条件齐次化来得方便。

**练习 13.12**　对于最一般的定解问题，例如

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
f(x,t),
\qquad
0<x<l,\ t>0,
$$

$$
\left.u\right|_{x=0}
=
\mu(t),
\qquad
\left.u\right|_{x=l}
=
\nu(t),
\qquad
t\ge0,
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
\psi(x),
\qquad
0\le x\le l,
$$

试用分离变量法求解。

这里值得重新讨论一下前面的例 13.3。我们是否也可以采用将方程和非齐次边界条件同时齐次化的方法求解？问题是，这里出现的是热传导方程，只含有对于时间 $t$ 的一阶偏导数，因此如果简单地将齐次化函数取为 $f(x)\sin\omega t$，不可能使方程仍然保持齐次。这时比较简单的办法是采用“复数解法”，即引进复函数 $\mathcal U(x,t)$，它满足定解问题

$$
\frac{\partial\mathcal U}{\partial t}
-
\kappa\frac{\partial^2\mathcal U}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0;
\tag{13.94a}
$$

$$
\left.\mathcal U\right|_{x=0}
=
Ae^{i\omega t},
\qquad
\left.\mathcal U\right|_{x=l}
=
0,
\qquad
t\ge0;
\tag{13.94b}
$$

$$
\left.\mathcal U\right|_{t=0}
=
0,
\qquad
0\le x\le l.
\tag{13.94c}
$$

取齐次化函数为

$$
\mathcal V(x,t)=F(x)e^{i\omega t},
$$

则 $F(x)$ 满足常微分方程边值问题

$$
F''(x)-i\omega F(x)=0,
$$

$$
F(0)=A,
\qquad
F'(l)=0.
$$

求出 $\mathcal V(x,t)$ 后，再令

$$
\mathcal U(x,t)
=
\mathcal V(x,t)
+
\mathcal W(x,t),
$$

则 $\mathcal W(x,t)$ 满足齐次方程、齐次边界条件的定解问题

$$
\frac{\partial\mathcal W}{\partial t}
-
\kappa\frac{\partial^2\mathcal W}{\partial x^2}
=
0,
\qquad
0<x<l,\ t>0,
$$

$$
\left.\mathcal W\right|_{x=0}
=
0,
\qquad
\left.\mathcal W\right|_{x=l}
=
0,
\qquad
t\ge0,
$$

$$
\left.\mathcal W\right|_{t=0}
=
-\mathcal V(x,0),
\qquad
0\le x\le l.
$$

不难求出 $\mathcal W(x,t)$，因而也就求得了 $\mathcal U(x,t)$。再求虚部，就得到原始定解问题 (13.82) 的解 $u(x,t)$。请读者完成这个计算。

读者可能会想到，对于定解问题 (13.78)，应该也可以将齐次方程 (13.78a) 和齐次初始条件 (13.78c) 分离变量而得到

$$
T''(t)+\lambda a^2T(t)=0,
\tag{13.95a}
$$

$$
T(0)=0,
\qquad
T'(0)=0.
\tag{13.95b}
$$

如果这也构成本征值问题，并且也能够求出本征值和本征函数的话，那么似乎可以因此求得与之相关的一般解，从而利用非齐次边界条件 (13.78b) 定出叠加系数即可。其实，在本章一开始的练习 13.1 中，我们就提出了同样的问题。从实际计算的结果来看，这是不可行的，原因是在齐次初始条件 (13.95b) 下，对于任何 $\lambda$ 值，方程 (13.95a) 总只有零解。容易看出，任何有限的 $t$ 点都是方程 (13.95a) 的常点，因此，根据常微分方程理论，解 $T(t)$ 在全 $t$ 平面上解析，故而可以在 $t=0$ 点展开为 Taylor 级数。齐次的初始条件就决定了级数的全部系数均为 0。正是基于这个原因，我们要再次强调，构成本征值问题，一定是含有待定参数的齐次微分方程和齐次边界条件。边界条件和初始条件的这种差别，在物理意义上，是由于它们分别对应于空间变量和时间变量，在数学上，是由于边界条件给出的是函数（或者导数，或者它们的线性组合）在不同点（例如 $x=0$ 和 $x=l$）的数值，初始条件则是给出函数在同一点（即 $t=0$）的数值与导数值。边界条件和初始条件在分离变量法中起着不同的作用。

对于稳定问题，不含时间变量 $t$，定解条件只是边界条件，这样，即使全部边界条件都是非齐次的，在选择将哪些边界条件齐次化而用于构成本征值问题，其余的边界条件留作定叠加系数时，就有充分的自由。实际上，对于定解问题

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
f(x,y),
\qquad
0<x<a,\ 0<y<b,
$$

$$
\left.u(x,y)\right|_{x=0}
=
\xi(y),
\qquad
\left.u(x,y)\right|_{x=a}
=
\eta(y),
\qquad
0\le y\le b,
$$

$$
\left.u(x,y)\right|_{y=0}
=
\phi(x),
\qquad
\left.u(x,y)\right|_{y=b}
=
\psi(x),
\qquad
0\le x\le a,
$$

还可以更简单地令 $u(x,y)=u_1(x,y)+u_2(x,y)$，使得 $u_1(x,y)$ 和 $u_2(x,y)$ 满足的定解问题中，都只有一组边界条件非齐次的，即

$$
\frac{\partial^2u_1}{\partial x^2}
+
\frac{\partial^2u_1}{\partial y^2}
=
f(x,y),
\qquad
0<x<a,\ 0<y<b,
$$

$$
\left.u_1(x,y)\right|_{x=0}
=
0,
\qquad
\left.u_1(x,y)\right|_{x=a}
=
0,
\qquad
0\le y\le b,
$$

$$
\left.u_1(x,y)\right|_{y=0}
=
\phi(x),
\qquad
\left.u_1(x,y)\right|_{y=b}
=
\psi(x),
\qquad
0\le x\le a
$$

和

$$
\frac{\partial^2u_2}{\partial x^2}
+
\frac{\partial^2u_2}{\partial y^2}
=
0,
\qquad
0<x<a,\ 0<y<b,
$$

$$
\left.u_2(x,y)\right|_{x=0}
=
\xi(y),
\qquad
\left.u_2(x,y)\right|_{x=a}
=
\eta(y),
\qquad
0\le y\le b,
$$

$$
\left.u_2(x,y)\right|_{y=0}
=
0,
\qquad
\left.u_2(x,y)\right|_{y=b}
=
0,
\qquad
0\le x\le a.
$$

读者不难采用前几节中已经讲过的相应方法求出 $u_1(x,y)$ 和 $u_2(x,y)$。

## 习题

1. 求解第十一章习题第 1 题。

2. 长为 $l$、两端固定的均匀弦，初始时，弦被拉开如图 13.1，达到平衡后突然放手。求解此问题。

![Figure 13.1](...)

图 13.1

3. 求解细杆的导热问题。杆长 $l$，两端 $(x=0,l)$ 均保持为零度，初始温度分布为

$$
\left.u\right|_{t=0}
=
b\frac{x(l-x)}{l^2}.
$$

4. 一均匀各向同性的弹性薄膜，$0\le x\le l,\ 0\le y\le l$，四周夹紧，初始位移为 $Axy(l-x)(l-y)$，初始速度为 0。求解膜的横振动。

5. 求解：

$$
\frac{\partial^2u}{\partial x^2}
+
\frac{\partial^2u}{\partial y^2}
=
0,
$$

$$
\left.u\right|_{x=0}
=
u_0,
\qquad
\left.u\right|_{x=a}
=
u_0
\left[
3\left(
\frac yb
\right)^2
-
2\left(
\frac yb
\right)^3
\right],
$$

$$
\left.
\frac{\partial u}{\partial y}
\right|_{y=0}
=
0,
\qquad
\left.
\frac{\partial u}{\partial y}
\right|_{y=b}
=
0,
$$

其中 $b$ 为已知常数。

6. 求解：

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
bx(l-x),
$$

$$
\left.u\right|_{x=0}
=
0,
\qquad
\left.u\right|_{x=l}
=
0,
$$

$$
\left.u\right|_{t=0}
=
0,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
0.
$$

7. 在矩形区域 $0\le x\le a,\ -b/2\le y\le b/2$ 中求解：

(1)

$$
\begin{cases}
\nabla^2u=-2,\\
\left.u\right|_{x=0,a}=0,\\
\left.u\right|_{y=\pm b/2}=0;
\end{cases}
$$

(2)

$$
\begin{cases}
\nabla^2u=-x^2y,\\
\left.u\right|_{x=0,a}=0,\\
\left.u\right|_{y=\pm b/2}=0.
\end{cases}
$$

8. 一细长杆，$x=0$ 端固定，$x=l$ 端受周期力 $A\sin\omega t$ 作用。设初位移和初速度均为 0，求解此杆的纵振动问题。

9. 试求下列定解问题之解：

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
$$

$$
\left.u\right|_{x=0}
=
\cos\frac{\pi}{l}at,
\qquad
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
0,
$$

$$
\left.u\right|_{t=0}
=
\cos\frac{\pi}{l}x,
\qquad
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
\sin\frac{\pi}{2l}x.
$$

10. 求解下列定解问题：

$$
\frac{\partial u}{\partial t}
-
\kappa
\frac{\partial^2u}{\partial x^2}
=
0,
$$

$$
\left.u\right|_{x=0}
=
Ae^{i\omega t},
\qquad
\left.u\right|_{x=l}
=
0,
$$

$$
\left.u\right|_{t=0}
=
0.
$$

11. 求解下列定解问题：

$$
\frac{\partial u}{\partial t}
-
\kappa
\frac{\partial^2u}{\partial x^2}
=
0,
$$

$$
\left.u\right|_{x=0}
=
A\exp\{-\alpha^2\kappa t\},
\qquad
\left.u\right|_{x=l}
=
B\exp\{-\beta^2\kappa t\},
$$

$$
\left.u\right|_{t=0}
=
0.
$$

12. 当层状铀块的厚度超过一定临界值时，中子浓度将随时间而增高，以致引起铀块爆炸。这就是原子弹爆炸的基本过程。试估计层状铀块的临界厚度。中子浓度满足的偏微分方程见第十一章习题第 2 题，假定边界条件为齐次的第一类边界条件。
