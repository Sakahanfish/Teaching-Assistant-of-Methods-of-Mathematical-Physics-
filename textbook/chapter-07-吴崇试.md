# 第七章 Γ 函数

## §7.1 Γ 函数的定义

Γ 函数是最基本的特殊函数。常用的定义是

$$
\Gamma(z)=\int_0^\infty e^{-t}t^{z-1}\,dt,
\qquad
\operatorname{Re}z>0,
\tag{7.1}
$$

称为第二类 Euler 积分，其中的积分变量 $t$ 应该理解为 $\arg t=0$。

首先证明积分在右半平面代表一个解析函数。因为这是一个反常积分，既是瑕积分（在 $t=0$ 端），又是无穷积分，所以要拆成两部分来分别讨论①：

$$
\int_0^\infty e^{-t}t^{z-1}\,dt
=
\int_0^1 e^{-t}t^{z-1}\,dt
+
\int_1^\infty e^{-t}t^{z-1}\,dt.
\tag{7.2}
$$

① 函数

$$
\Gamma(z,\alpha)=\int_0^\alpha e^{-t}t^{z-1}\,dt
$$

和

$$
\gamma(z,\alpha)=\int_\alpha^\infty e^{-t}t^{z-1}\,dt
$$

均称为不完全 Γ 函数。

先看第二部分。显然当 $t\geqslant1$ 时，被积函数 $e^{-t}t^{z-1}$ 是 $t$ 的连续函数，并且作为 $z$ 的函数，在全复平面解析。由定理 4.3 可知，要证明它代表一个解析函数，只需证明积分在某个区域内一致收敛。因为

$$
e^t=\sum_{n=0}^\infty\frac{t^n}{n!},
\qquad
|t|<\infty,
$$

所以 $t>0$ 时，$\forall$ 正整数 $N$，

$$
e^t>\frac{t^N}{N!},
\qquad
e^{-t}<\frac{N!}{t^N}.
$$

对于 $z$ 复平面内任一有界闭区域 $\overline G$，$\exists x_0$，$\forall z\in\overline G$，均有 $\operatorname{Re}z<x_0$（见图 7.1），因此

$$
\left|e^{-t}t^{z-1}\right|
<
N!\,t^{x_0-N-1},
\qquad
t\geqslant1.
$$

![Figure 7.1](...)

图 7.1　Γ 函数的解析区域

这样，只要选择足够大的 $N$（使得 $N>x_0$），积分

$$
\int_1^\infty t^{x_0-N-1}\,dt
$$

就收敛，故

$$
\Gamma(z,1)=\int_1^\infty e^{-t}t^{z-1}\,dt
$$

在 $z$ 复平面内的任一有界闭区域中一致收敛，因此在全 $z$ 复平面解析。

要证明第一部分的瑕积分在 $z$ 右半平面解析，关键也是证明它的一致收敛性。因为

$$
\left|e^{-t}t^{z-1}\right|
=
e^{-t}t^{x-1},
\qquad
x=\operatorname{Re}z,
\qquad
t\in\mathbb R,
$$

所以，对于 $z$ 右半平面的任一闭区域，$\exists\delta>0$，使得对 $z$ 右半平面闭区域内的任意一点 $z$，都有

$$
\operatorname{Re}z=x\geqslant\delta>0,
$$

因此

$$
\left|e^{-t}t^{z-1}\right|
\leqslant t^{\delta-1},
\qquad
0<t<1,
$$

而

$$
\int_0^1 t^{\delta-1}\,dt
$$

收敛，由此即可推知积分

$$
\gamma(z,1)=\int_0^1 e^{-t}t^{z-1}\,dt
$$

在 $z$ 右半平面的任一闭区域中一致收敛，故在 $z$ 右半平面解析。

把两部分合起来，就证得 Γ 函数 (7.1) 在 $z$ 的右半平面解析。$\square$

下面再把 Γ 函数的定义加以扩充。首先，积分路径并不需要限定在实轴上，而可以修改为

$$
\Gamma(z)=\int_L e^{-t}t^{z-1}\,dt,
\qquad
\operatorname{Re}z>0,
\tag{7.3}
$$

积分路径 $L$ 是 $t$ 复平面内从 $t=0$ 出发的半射线，$\arg t=\alpha$ 为常数，$|\alpha|<\pi/2$。取围道 $C$ 如图 7.2，应用留数定理讨论复变积分

$$
\oint_C e^{-t}t^{z-1}\,dt,
$$

就能证明这个结论。请读者补足证明。

![Figure 7.2](...)

图 7.2　Γ 函数定义的扩充

这个结果还可以进一步修改：积分路径 $L$ 可以是 $t$ 复平面内从 $t=0$ 出发的任意分段光滑曲线，只要最后以 $\operatorname{Re}t\to+\infty$ 的方式趋于无穷远点即可。也请读者补足这个证明。

上面介绍的 Γ 函数的定义 (7.1) 式当然只适用于 $\operatorname{Re}z>0$。为了延拓到 $z$ 的全复平面，注意在前面的证明中，积分的第二部分 $\Gamma(z,1)$ 在全复平面解析，因此，只要用适当的方法将积分的第一部分 $\gamma(z,1)$ 延拓到全复平面即可。比较直接的方法是将指数函数作 Taylor 展开：

$$
\begin{aligned}
\int_0^1 e^{-t}t^{z-1}\,dt
&=
\sum_{n=0}^\infty\frac{(-)^n}{n!}
\int_0^1 t^{n+z-1}\,dt\\
&=
\sum_{n=0}^\infty
\frac{(-)^n}{n!}\frac{1}{n+z},
\qquad
\operatorname{Re}z>0.
\end{aligned}
$$

这个结果是在 $\operatorname{Re}z>0$ 的条件下得到的。但等式右端的无穷级数显然在全复平面内（$z\ne0,-1,-2,\cdots$）一致收敛，因此在全复平面解析（$z\ne0,-1,-2,\cdots$）。也就是说，等式右端的级数表达式就是左端积分表达式在全复平面内的解析延拓，

$$
\Gamma(z)
=
\int_1^\infty e^{-t}t^{z-1}\,dt
+
\sum_{n=0}^\infty
\frac{(-)^n}{n!}\frac{1}{n+z},
\qquad
z\ne0,-1,-2,\cdots.
\tag{7.4}
$$

在许多学科中都可能用到 Γ 函数。下面分别举一个物理学和概率统计方面的例子。

**例 7.1**　按照 Maxwell 速率分布率，速率处于 $v$ 与 $v+dv$ 间的粒子数 $dN$ 为

$$
dN
=
4\pi N
\left(
\frac{m}{2\pi kT}
\right)^{3/2}
e^{-mv^2/2kT}v^2\,dv,
$$

其中 $N$ 为总粒子数。求平均值

$$
\langle v^n\rangle
=
\frac1N\int_0^\infty v^n\,dN.
$$

**解**　只需直接计算积分

$$
\langle v^n\rangle
=
\frac1N\int_0^\infty v^n\,dN
=
4\pi
\left(
\frac{m}{2\pi kT}
\right)^{3/2}
\int_0^\infty
e^{-mv^2/2kT}v^{n+2}\,dv
$$

即可。令 $x=mv^2/2kT$，就能算出

$$
\begin{aligned}
\langle v^n\rangle
&=
4\pi
\left(
\frac{m}{2\pi kT}
\right)^{3/2}
\frac12
\left(
\frac{2kT}{m}
\right)^{(n+3)/2}
\int_0^\infty
e^{-x}x^{(n+1)/2}\,dx\\
&=
\frac{2}{\sqrt\pi}
\left(
\frac{2kT}{m}
\right)^{n/2}
\Gamma\left(\frac{n+3}{2}\right).
\end{aligned}
$$

**例 7.2**　$\gamma$ 分布的概率密度为

$$
f(x)=
\begin{cases}
\dfrac{1}{\beta^\alpha\Gamma(\alpha)}
x^{\alpha-1}e^{-x/\beta},
& x>0,\\[6pt]
0,
& x\leqslant0,
\end{cases}
$$

其中 $\alpha,\beta$ 均为已知正数。求 $x$ 的平均值

$$
\langle x\rangle
=
\int_{-\infty}^{\infty}xf(x)\,dx
$$

及标准偏差

$$
\sqrt{\langle x^2\rangle-\langle x\rangle^2}.
$$

**解**　可以直接计算出平均值

$$
\langle x\rangle
=
\frac{1}{\beta^\alpha\Gamma(\alpha)}
\int_0^\infty
x^\alpha e^{-x/\beta}\,dx
=
\frac{\beta^{\alpha+1}\Gamma(\alpha+1)}
{\beta^\alpha\Gamma(\alpha)}
=
\beta\alpha,
$$

$$
\langle x^2\rangle
=
\frac{1}{\beta^\alpha\Gamma(\alpha)}
\int_0^\infty
x^{\alpha+1}e^{-x/\beta}\,dx
=
\frac{\beta^{\alpha+2}\Gamma(\alpha+2)}
{\beta^\alpha\Gamma(\alpha)}
=
\beta^2\alpha(\alpha+1).
$$

由此又能算出标准偏差

$$
\sqrt{\langle x^2\rangle-\langle x\rangle^2}
=
\sqrt{\beta^2\alpha(\alpha+1)-(\beta\alpha)^2}
=
\beta\sqrt\alpha.
$$

## §7.2 Γ 函数的基本性质

**性质 1**

$$
\Gamma(1)=1.
\tag{7.5}
$$

直接在 Γ 函数的定义 (7.1) 式中代入 $z=1$ 即可得到这个结果。

**性质 2**

$$
\Gamma(z+1)=z\Gamma(z).
\tag{7.6}
$$

**证**　根据 Γ 函数的定义，

$$
\begin{aligned}
\Gamma(z+1)
&=
\int_0^\infty e^{-t}t^z\,dt\\
&=
-\left.e^{-t}t^z\right|_0^\infty
+
z\int_0^\infty e^{-t}t^{z-1}\,dt\\
&=
z\int_0^\infty e^{-t}t^{z-1}\,dt
=
z\Gamma(z).
\end{aligned}
$$

$\square$

这个结果可以从两个角度来理解。一方面，尽管在证明中用到了条件 $\operatorname{Re}z>0$，但由于 $\Gamma(z+1)$ 和 $z\Gamma(z)$ 都在全复平面解析（$z=0,-1,-2,\cdots$ 除外），因此，根据解析延拓原理可以断定，这个递推关系在全复平面均成立。另一方面，也可以直接通过递推关系来完成 Γ 函数的解析延拓。这时，因为 $\Gamma(z+1)$ 在半平面 $\operatorname{Re}z>-1$ 内解析，因此就可以把

$$
\Gamma(z)
=
\frac1z\Gamma(z+1),
\qquad
z\ne0
\tag{7.7}
$$

看成是 $\Gamma(z)$ 在区域 $\operatorname{Re}z>-1$ 内的定义，而 $z=0$ 点是 Γ 函数的一阶极点，

$$
\operatorname{res}\Gamma(0)=1.
$$

重复上述步骤，还可以将 Γ 函数解析延拓到区域 $\operatorname{Re}z>-2$，

$$
\Gamma(z)
=
\frac{1}{z(z+1)}\Gamma(z+2),
\qquad
z\ne0,-1.
\tag{7.8}
$$

$z=-1$ 也是 Γ 函数的一阶极点，

$$
\operatorname{res}\Gamma(-1)=-1.
$$

如此继续，就可以将 Γ 函数延拓到全复平面，而 $z=0,-1,-2,\cdots$ 都是它的一阶极点，

$$
\operatorname{res}\Gamma(-n)
=
\frac{(-1)^n}{n!}.
\tag{7.9}
$$

**推论 7.1**　对于正整数 $n$，

$$
\Gamma(n)=(n-1)!.
\tag{7.10}
$$

正是因为这个原因，Γ 函数又称为阶乘函数。

**性质 3　互余宗量关系**

$$
\Gamma(z)\Gamma(1-z)
=
\frac{\pi}{\sin\pi z}.
\tag{7.11}
$$

证明见 §7.4。

**推论 7.2**

$$
\Gamma\left(\frac12\right)=\sqrt\pi.
\tag{7.12}
$$

只需在 (7.11) 式中代入 $z=1/2$，并注意 $\Gamma(1/2)>0$（因被积函数值恒正）即可证得。

**推论 7.3**　Γ 函数在全复平面无零点。

**证**　（反证法）因为 $\pi/\sin\pi z\ne0$，所以

$$
\Gamma(z)\Gamma(1-z)\ne0.
$$

假设在 $z=z_0$ 点有 $\Gamma(z_0)=0$，则必有 $\Gamma(1-z_0)=\infty$。这只能发生在 $1-z_0=-n$（亦即 $z_0=n+1$），$n=0,1,2,\cdots$ 处。但

$$
\Gamma(z_0)
=
\Gamma(n+1)
=
n!,
$$

与所设矛盾。故 $\Gamma(z)$ 在全复平面无零点。$\square$

图 7.3 中给出了 $\Gamma(x)$（$x$ 为实数）的图形。

![Figure 7.3](...)

图 7.3　自变量为实数时的 Γ 函数值

**性质 4　倍乘公式**

$$
\Gamma(2z)
=
2^{2z-1}\pi^{-1/2}
\Gamma(z)
\Gamma\left(z+\frac12\right).
\tag{7.13}
$$

这个公式的证明亦见 §7.4。

**性质 5**　Γ 函数的渐近展开，即 Stirling 公式：当 $|z|\to\infty$，$|\arg z|<\pi$ 时，有

$$
\Gamma(z)
\sim
z^{z-1/2}e^{-z}\sqrt{2\pi}
\left\{
1+\frac{1}{12z}
+\frac{1}{288z^2}
-\frac{139}{51840z^3}
-\frac{571}{2488320z^4}
+\cdots
\right\},
\tag{7.14}
$$

$$
\ln\Gamma(z)
\sim
\left(z-\frac12\right)\ln z-z
+\frac12\ln(2\pi)
+\frac{1}{12z}
-\frac{1}{360z^3}
+\frac{1}{1260z^5}
-\frac{1}{1680z^7}
+\cdots.
\tag{7.15}
$$

在物理中更常用的结果是

$$
\ln n!\sim n\ln n-n.
\tag{7.16}
$$

下面的 §7.6 中，将就实数的情形推导 (7.14) 式。复数的一般情形下的推导，可参阅参考书目 [1]，9.6 节。

## §7.3 ψ 函数

$\psi$ 函数是 Γ 函数的对数微商

$$
\psi(z)
=
\frac{d\ln\Gamma(z)}{dz}
=
\frac{\Gamma'(z)}{\Gamma(z)}.
\tag{7.17}
$$

根据 Γ 函数的性质，可以得出 $\psi(z)$ 的下列性质：

(1) $z=0,-1,-2,\cdots$ 都是 $\psi(z)$ 的一阶极点，留数均为 $-1$；除了这些点以外，$\psi(z)$ 在全复平面解析。

(2)

$$
\psi(z+1)
=
\psi(z)+\frac1z.
\tag{7.18}
$$

$$
\psi(z+n)
=
\psi(z)
+\frac1z
+\frac1{z+1}
+\cdots
+\frac1{z+n-1},
\qquad
n=2,3,4,\cdots.
\tag{7.19}
$$

(3)

$$
\psi(1-z)-\psi(z)
=
\pi\cot\pi z.
\tag{7.20}
$$

(4)

$$
\psi(z)-\psi(-z)
=
-\frac1z-\pi\cot\pi z.
\tag{7.21}
$$

(5)

$$
\psi(2z)
=
\frac12\psi(z)
+\frac12\psi\left(z+\frac12\right)
+\ln2.
\tag{7.22}
$$

(6)

$$
\psi(z)
\sim
\ln z
-\frac{1}{2z}
-\frac{1}{12z^2}
+\frac{1}{120z^4}
-\frac{1}{252z^6}
+\cdots,
\qquad
z\to\infty,\ |\arg z|<\pi.
\tag{7.23}
$$

(7)

$$
\lim_{n\to\infty}
[\psi(z+n)-\ln n]
=
0.
\tag{7.24}
$$

$\psi$ 函数的特殊值有

$$
\psi(1)=-\gamma,
\qquad
\psi'(1)=\frac{\pi^2}{6},
$$

$$
\psi\left(\frac12\right)
=
-\gamma-2\ln2,
\qquad
\psi'\left(\frac12\right)
=
\frac{\pi^2}{2},
$$

$$
\psi\left(-\frac12\right)
=
-\gamma-2\ln2+2,
\qquad
\psi'\left(-\frac12\right)
=
\frac{\pi^2}{2}+4,
$$

$$
\psi\left(\frac14\right)
=
-\gamma-\frac\pi2-3\ln2,
\qquad
\psi\left(\frac34\right)
=
-\gamma+\frac\pi2-3\ln2,
$$

$$
\psi\left(\frac13\right)
=
-\gamma-\frac{\pi}{2\sqrt3}
-\frac32\ln3,
\qquad
\psi\left(\frac23\right)
=
-\gamma+\frac{\pi}{2\sqrt3}
-\frac32\ln3,
$$

其中 $\gamma=-\psi(1)$ 是数学中的一个基本常数，称为 Euler 常数①，

$$
\gamma
=
0.577\,215\,664\,901\,532\,860\,606\,512\,090\,082\,402\,431\,042\,159\,335\,939\,92\cdots.
\tag{7.25}
$$

① 又称 Euler-Mascheroni 常数，是最基本的数学常数之一。至今还不知道它是不是无理数，但猜想它是一个超越数。连分数的分析表明，如果它是有理数的话，则其分母必须大于 $10^{242\,080}$。到 2017 年，其数值已经计算到 $477\,511\,832\,674$ 位小数（见 https://en.wikipedia.org/wiki/Euler_constant）。

**例 7.3**　计算积分

$$
\int_0^\pi
\frac{\sin^2n\theta}{\sin\theta}\,d\theta.
$$

**解**　取半径为 1 的半圆形围道，于是，根据留数定理有

$$
\oint
\frac{z^{2n}-1}{z^2-1}\,dz
=
\int_0^\pi
\frac{e^{i2n\theta}-1}{e^{i2\theta}-1}
e^{i\theta}i\,d\theta
+
\int_{-1}^1
\frac{x^{2n}-1}{x^2-1}\,dx
=
0.
$$

因为

$$
\begin{aligned}
\int_0^\pi
\frac{e^{i2n\theta}-1}{e^{i2\theta}-1}
e^{i\theta}i\,d\theta
&=
\int_0^\pi
\frac{\cos2n\theta-1+i\sin2n\theta}
{2\sin\theta}\,d\theta\\
&=
-\int_0^\pi
\frac{\sin^2n\theta}{\sin\theta}\,d\theta
+
\frac{i}{2}
\int_0^\pi
\frac{\sin2n\theta}{\sin\theta}\,d\theta
\end{aligned}
$$

以及

$$
\begin{aligned}
\int_{-1}^1
\frac{x^{2n}-1}{x^2-1}\,dx
&=
\int_{-1}^1
(x^{2n-2}+x^{2n-4}+\cdots+1)\,dx\\
&=
2\left(
\frac{1}{2n-1}
+\frac{1}{2n-3}
+\cdots
+\frac13
+1
\right)\\
&=
\psi\left(n+\frac12\right)
-
\psi\left(\frac12\right),
\end{aligned}
$$

因此，

$$
-\int_0^\pi
\frac{\sin^2n\theta}{\sin\theta}\,d\theta
+
\frac{i}{2}
\int_0^\pi
\frac{\sin2n\theta}{\sin\theta}\,d\theta
+
\psi\left(n+\frac12\right)
-
\psi\left(\frac12\right)
=
0.
$$

比较实部，即得

$$
\int_0^\pi
\frac{\sin^2n\theta}{\sin\theta}\,d\theta
=
\psi\left(n+\frac12\right)
-
\psi\left(\frac12\right).
\tag{7.26}
$$

与此同时，比较虚部，还能得到

$$
\int_0^\pi
\frac{\sin2n\theta}{\sin\theta}\,d\theta
=
0.
\tag{7.27a}
$$

顺便指出，还可以计算类似于 (7.27a) 的积分

$$
\begin{aligned}
\int_0^\pi
\frac{\sin(2n+1)\theta}{\sin\theta}\,d\theta
&=
\int_0^\pi
\frac{e^{i(2n+1)\theta}-e^{-i(2n+1)\theta}}
{e^{i\theta}-e^{-i\theta}}\,d\theta\\
&=
\int_0^\pi
\left(
e^{i2n\theta}
+e^{i(2n-2)\theta}
+\cdots
+e^{-i(2n-2)\theta}
+e^{-i2n\theta}
\right)d\theta.
\end{aligned}
$$

这 $2n+1$ 项中只有一项对积分有贡献，因此

$$
\int_0^\pi
\frac{\sin(2n+1)\theta}{\sin\theta}\,d\theta
=
\int_0^\pi1\cdot d\theta
=
\pi,
\qquad
n=0,1,2,\cdots.
\tag{7.27b}
$$

利用 $\psi$ 函数，可以方便地求出通项为有理式的无穷级数

$$
\sum_{n=0}^\infty u_n
=
\sum_{n=0}^\infty\frac{p(n)}{d(n)}
\tag{7.28a}
$$

之和，其中 $p(n)$ 和 $d(n)$ 都是 $n$ 的多项式。设 $d(n)$ 是 $n$ 的 $m$ 次多项式，并且全部零点都是一阶零点，

$$
d(n)
=
(n+\alpha_1)(n+\alpha_2)\cdots(n+\alpha_m),
$$

即 $u_n$ 只有一阶极点，则可部分分式为

$$
u_n
=
\frac{p(n)}{d(n)}
=
\sum_{k=1}^m
\frac{a_k}{n+\alpha_k}.
\tag{7.28b}
$$

为了保证级数收敛，必须有

$$
\lim_{n\to\infty}u_n
=
\lim_{n\to\infty}n\cdot u_n
=
0,
\qquad
\text{即}
\qquad
\sum_{k=1}^m a_k=0.
\tag{7.28c}
$$

利用 $\psi$ 函数的递推关系 (7.19)，即可求得

$$
\begin{aligned}
\sum_{n=0}^N u_n
&=
\sum_{k=1}^m
a_k[\psi(\alpha_k+N)-\psi(\alpha_k)]\\
&=
\sum_{k=1}^m
a_k[\psi(\alpha_k+N)-\psi(\alpha_k)]
-
\ln N\sum_{k=1}^m a_k\\
&=
\sum_{k=1}^m
a_k[\psi(\alpha_k+N)-\ln N-\psi(\alpha_k)].
\end{aligned}
$$

取极限 $N\to\infty$，注意到 (7.24) 式，即得

$$
\sum_{n=0}^\infty u_n
=
\lim_{N\to\infty}
\sum_{k=1}^m
a_k[\psi(\alpha_k+N)-\ln N]
-
\sum_{k=1}^m a_k\psi(\alpha_k)
=
-\sum_{k=1}^m a_k\psi(\alpha_k).
\tag{7.29}
$$

**例 7.4**　求无穷级数

$$
\sum_{n=0}^\infty
\frac{1}{(3n+1)(3n+2)(3n+3)}
$$

之和。

**解**　因为

$$
\frac{1}{(3n+1)(3n+2)(3n+3)}
=
\frac16\frac{1}{n+1/3}
-
\frac13\frac{1}{n+2/3}
+
\frac16\frac{1}{n+1},
$$

所以，根据上面给出的求和公式，有

$$
\sum_{n=0}^\infty
\frac{1}{(3n+1)(3n+2)(3n+3)}
=
-\frac16
\left[
\psi\left(\frac13\right)
-
2\psi\left(\frac23\right)
+
\psi(1)
\right].
$$

代入 $\psi$ 函数的特殊值，即得

$$
\sum_{n=0}^\infty
\frac{1}{(3n+1)(3n+2)(3n+3)}
=
\frac14
\left(
\frac{\pi}{\sqrt3}
-
\ln3
\right).
$$

**例 7.5**　求无穷级数

$$
\sum_{n=0}^\infty
\frac{1}{n^2+a^2}
$$

之和，其中 $a>0$。

**解**　因为

$$
\begin{aligned}
\sum_{n=0}^\infty
\frac{1}{n^2+a^2}
&=
\frac{i}{2a}
\sum_{n=0}^\infty
\left(
\frac{1}{n+ia}
-
\frac{1}{n-ia}
\right)\\
&=
-\frac{i}{2a}
[\psi(ia)-\psi(-ia)],
\end{aligned}
$$

利用上面列出的 $\psi$ 函数的性质 4，有

$$
\psi(ia)-\psi(-ia)
=
-\frac{1}{ia}
-\pi\cot i\pi a
=
i\left(
\frac1a+\pi\coth\pi a
\right),
$$

就可以求得

$$
\sum_{n=0}^\infty
\frac{1}{n^2+a^2}
=
\frac{1}{2a^2}
(1+\pi a\coth\pi a).
\tag{7.30}
$$

这个结果也可以用其他方法得到，如 §8.5 的例 8.12。

如果 $u_n$ 还有二阶极点，如

$$
d(n)
=
(n+\alpha_1)(n+\alpha_2)\cdots(n+\alpha_m)
(n+\beta_1)^2(n+\beta_2)^2\cdots(n+\beta_l)^2,
$$

则

$$
u_n
=
\frac{p(n)}{d(n)}
=
\sum_{k=1}^m
\frac{a_k}{n+\alpha_k}
+
\sum_{k=1}^l
\left[
\frac{b_{1k}}{n+\beta_k}
+
\frac{b_{2k}}{(n+\beta_k)^2}
\right].
\tag{7.31a}
$$

相应地，级数收敛的条件是

$$
\sum_{k=1}^m a_k
+
\sum_{k=1}^l b_{1k}
=
0.
\tag{7.31b}
$$

根据 $\psi$ 函数的递推关系，即得

$$
\sum_{n=0}^\infty u_n
=
-\sum_{k=1}^m a_k\psi(\alpha_k)
-
\sum_{k=1}^l
[b_{1k}\psi(\beta_k)-b_{2k}\psi'(\beta_k)].
\tag{7.31c}
$$

**例 7.6**　求无穷级数

$$
\sum_{n=0}^\infty
\frac{1}{(n+1)^2(2n+1)^2}
$$

之和。

**解**　因为

$$
\frac{1}{(n+1)^2(2n+1)^2}
=
\left[
\frac{4}{n+1}
+
\frac{1}{(n+1)^2}
\right]
-
\left[
\frac{4}{n+1/2}
-
\frac{1}{(n+1/2)^2}
\right],
$$

所以

$$
\begin{aligned}
\sum_{n=0}^\infty
\frac{1}{(n+1)^2(2n+1)^2}
&=
-[4\psi(1)-\psi'(1)]
+
\left[
4\psi\left(\frac12\right)
+
\psi'\left(\frac12\right)
\right]\\
&=
\frac{2\pi^2}{3}
-
8\ln2.
\end{aligned}
\tag{7.32}
$$

## §7.4 B 函数

B 函数是由第一类 Euler 积分定义的：

$$
\mathrm B(p,q)
=
\int_0^1
t^{p-1}(1-t)^{q-1}\,dt,
\qquad
\operatorname{Re}p>0,\ \operatorname{Re}q>0,
\tag{7.33}
$$

其中的积分变量 $t$ 应该理解为

$$
\arg t=\arg(1-t)=0.
$$

令 $t=\sin^2\theta$，还可以得到 B 函数的另一个表达式

$$
\mathrm B(p,q)
=
2\int_0^{\pi/2}
\sin^{2p-1}\theta
\cos^{2q-1}\theta\,d\theta.
\tag{7.33'}
$$

B 函数可以用 Γ 函数表示出来，

$$
\mathrm B(p,q)
=
\frac{\Gamma(p)\Gamma(q)}
{\Gamma(p+q)}.
\tag{7.34}
$$

**证**　在 $\operatorname{Re}p>0,\ \operatorname{Re}q>0$ 的条件下，显然有

$$
\Gamma(p)
=
\int_0^\infty e^{-t}t^{p-1}\,dt
=
2\int_0^\infty e^{-x^2}x^{2p-1}\,dx,
$$

$$
\Gamma(q)
=
2\int_0^\infty e^{-y^2}y^{2q-1}\,dy,
$$

于是

$$
\Gamma(p)\Gamma(q)
=
4\int_0^\infty\int_0^\infty
e^{-(x^2+y^2)}
x^{2p-1}y^{2q-1}\,dx\,dy.
$$

令

$$
x=r\sin\theta,
\qquad
y=r\cos\theta,
$$

得

$$
\begin{aligned}
\Gamma(p)\Gamma(q)
&=
4\int_0^\infty\int_0^{\pi/2}
e^{-r^2}
(r\sin\theta)^{2p-1}
(r\cos\theta)^{2q-1}
r\,dr\,d\theta\\
&=
4\int_0^\infty
e^{-r^2}r^{2(p+q)-1}\,dr
\int_0^{\pi/2}
\sin^{2p-1}\theta
\cos^{2q-1}\theta\,d\theta\\
&=
\Gamma(p+q)\mathrm B(p,q).
\end{aligned}
$$

$\square$

利用这个关系式，可把 B 函数解析延拓到 $p$ 或 $q$ 的全复平面。

从 B 函数的定义 (7.33)，或者 (7.34) 式，立即可以看出 $\mathrm B(p,q)$ 对于 $p$ 和 $q$ 是对称的：

$$
\mathrm B(p,q)
=
\mathrm B(q,p).
\tag{7.35}
$$

现在根据 B 函数和 Γ 函数的关系式 (7.34)，补证 Γ 函数的两个性质，即互余宗量关系 (7.11) 和倍乘公式 (7.13)。首先，在 (7.34) 式中令 $p=z,\ q=1-z$，即得

$$
\mathrm B(z,1-z)
=
\frac{\Gamma(z)\Gamma(1-z)}{\Gamma(1)}
=
\Gamma(z)\Gamma(1-z).
$$

另一方面，

$$
\mathrm B(z,1-z)
=
\int_0^1
t^{z-1}(1-t)^{-z}\,dt.
$$

令 $x=t/(1-t)$，上式即可化为

$$
\mathrm B(z,1-z)
=
\int_0^\infty
\frac{x^{z-1}}{1+x}\,dx.
$$

这个积分在 §6.6 中已经计算过（见 (6.32) 式），这样就证得

$$
\Gamma(z)\Gamma(1-z)
=
\mathrm B(z,1-z)
=
\frac{\pi}{\sin\pi z}.
$$

这个证明是在 $0<\operatorname{Re}z<1$ 的条件下得到的。但由于等式的两端在全复平面都解析，因此这个等式在全复平面均成立。$\square$

要证明 (7.13) 式，则可以通过积分

$$
\int_{-1}^1
(1-x^2)^{z-1}\,dx,
\qquad
\operatorname{Re}z>0
$$

的计算得到。令 $x^2=t$，则得

$$
\begin{aligned}
\int_{-1}^1
(1-x^2)^{z-1}\,dx
&=
\int_0^1
(1-t)^{z-1}t^{-1/2}\,dt\\
&=
\mathrm B\left(z,\frac12\right)
=
\frac{\Gamma(z)\Gamma(1/2)}
{\Gamma(z+1/2)}.
\end{aligned}
$$

若作变换

$$
1+x=2t,
\qquad
1-x=2(1-t),
$$

则有另一种形式的结果：

$$
\begin{aligned}
\int_{-1}^1
(1-x^2)^{z-1}\,dx
&=
2^{2z-1}
\int_0^1
t^{z-1}(1-t)^{z-1}\,dt\\
&=
2^{2z-1}\mathrm B(z,z)
=
2^{2z-1}
\frac{\Gamma(z)\Gamma(z)}
{\Gamma(2z)}.
\end{aligned}
$$

于是

$$
\frac{\Gamma(z)\Gamma(1/2)}
{\Gamma(z+1/2)}
=
2^{2z-1}
\frac{\Gamma(z)\Gamma(z)}
{\Gamma(2z)},
$$

此即 (7.13) 式。这里的证明仍然是在 $\operatorname{Re}z>0$ 的条件下进行的。但是，正如前面多次论证过的，结果在全复平面都成立。$\square$

**例 7.7**　应用 B 函数重新计算积分

$$
\int_0^1
\frac{\sqrt[4]{x(1-x)^3}}
{(1+x)^3}\,dx.
$$

**解**　此积分曾经应用留数定理计算过（见例 6.12）。现在选择作特定的自变量变换而化为 B 函数。为此，令

$$
t=\frac{2x}{1+x},
\qquad
1-t=\frac{1-x}{1+x},
\qquad
dt=\frac{2\,dx}{(1+x)^2},
$$

于是直接计算得

$$
\begin{aligned}
\int_0^1
\frac{\sqrt[4]{x(1-x)^3}}
{(1+x)^3}\,dx
&=
\int_0^1
\left(\frac{x}{1+x}\right)^{1/4}
\left(\frac{1-x}{1+x}\right)^{3/4}
\frac{dx}{(1+x)^2}\\
&=
2^{-5/4}
\int_0^1
t^{1/4}(1-t)^{3/4}\,dt\\
&=
\frac{2^{-5/4}}{\Gamma(3)}
\Gamma\left(\frac54\right)
\Gamma\left(\frac74\right)\\
&=
2^{-9/4}\cdot\frac14\cdot\frac34
\cdot\frac{\pi}{\sin(\pi/4)}\\
&=
\frac{3\sqrt[4]{2}}{64}\pi.
\end{aligned}
$$

## *§7.5 Γ 函数的普遍表达式

Γ 函数的定义（第二类 Euler 积分 (7.1) 式）只适用于右半平面。为了弥补这一缺陷，本节不加证明地介绍 Γ 函数的另外几种表达式，包括围道积分表示和无穷乘积表示，它们都在全复平面成立（孤立奇点除外）。有关证明见参考书目 [12] 的第 3 章。

1. Γ 函数的围道积分表示

$$
\Gamma(z)
=
-\frac{1}{2i\sin\pi z}
\int_{-\infty}^{(0+)}
e^{-t}(-t)^{z-1}\,dt,
\qquad
|\arg(-t)|<\pi.
\tag{7.36}
$$

其中的积分围道为：从上半平面挨近正实轴无穷远处出发，左行绕原点正向一周，再右行到下半平面挨近正实轴无穷远处（见图 7.4）。此式在全 $z$ 复平面成立，但 $z=$ 整数除外。

![Figure 7.4](...)

图 7.4　Γ 函数的围道积分表示 (1)

Γ 函数的另一个围道积分表示是

$$
\frac{1}{\Gamma(z)}
=
\frac{1}{2\pi i}
\int_{-\infty}^{(0+)}
e^t t^{-z}\,dt,
\qquad
|\arg t|<\pi.
\tag{7.37}
$$

积分围道从下半平面挨近负实轴无穷远处出发，右行绕原点正向一周，再左行到上半平面挨近负实轴无穷远处（见图 7.5）。此式在全 $z$ 复平面成立，包括 $z=$ 整数。

![Figure 7.5](...)

图 7.5　Γ 函数的围道积分表示 (2)

2. Γ 函数的 Euler 无穷乘积表示

$$
\Gamma(z)
=
\frac1z
\prod_{n=1}^\infty
\left\{
\left(1+\frac zn\right)^{-1}
\left(1+\frac1n\right)^z
\right\}.
\tag{7.38}
$$

此式对任何 $z$ 均成立，但极点 $z=$ 负整数除外。

3. Γ 函数的 Weierstrass 无穷乘积表示

$$
\frac{1}{\Gamma(z)}
=
ze^{\gamma z}
\prod_{n=1}^\infty
\left[
\left(1+\frac zn\right)e^{-z/n}
\right],
\tag{7.39}
$$

其中

$$
\gamma
=
\lim_{n\to\infty}
\left(
\sum_{k=1}^n\frac1k-\ln n
\right)
=
0.577\,215\,664\,901\cdots
\tag{7.40}
$$

就是 Euler 常数（§7.3，(7.25) 式）。这个无穷乘积给出了任何 $z$ 的 $\Gamma(z)$，同时指明了 $\Gamma(z)$ 的孤立奇点为一阶极点 $z=0,-1,-2,\cdots$ 而无零点。

从 Γ 函数的无穷乘积表示可得到一系列有意义的结果。例如

$$
\sin\pi z
=
\frac{\pi}{\Gamma(z)\Gamma(1-z)}
=
\pi z
\prod_{m=1}^\infty
\left(
1-\frac{z^2}{m^2}
\right),
\tag{7.41}
$$

$$
\cos\pi z
=
\frac{\sin2\pi z}{2\sin\pi z}
=
\prod_{m=1}^\infty
\left[
1-\frac{4z^2}{(2m-1)^2}
\right].
\tag{7.42}
$$

将这两式求对数微商，又可以得到

$$
\pi\tan\pi z
=
-8z\sum_{m=1}^\infty
\frac{1}{4z^2-(2m-1)^2},
\qquad
z\ne\pm\frac12,\pm\frac32,\cdots,
\tag{7.43}
$$

$$
\pi\cot\pi z
=
\frac1z
+
2z\sum_{m=1}^\infty
\frac{1}{z^2-m^2},
\qquad
z\ne0,\pm1,\pm2,\cdots,
\tag{7.44}
$$

$$
\pi\csc\pi z
=
\frac{\pi}{2}
\left[
\tan\frac{\pi z}{2}
+
\cot\frac{\pi z}{2}
\right]
=
\frac1z
+
2z\sum_{m=1}^\infty
\frac{(-1)^m}{z^2-m^2},
\qquad
z\ne0,\pm1,\pm2,\cdots,
\tag{7.45}
$$

$$
\pi\sec\pi z
=
\pi\csc\pi\left(\frac12-z\right)
=
4\sum_{m=1}^\infty
\frac{(-)^m(2m-1)}
{4z^2-(2m-1)^2},
\qquad
z\ne\pm\frac12,\pm\frac32,\cdots.
\tag{7.46}
$$

这些展开式称为有理分式展开。这种展开的形式不同于我们过去见过的 Taylor 展开和 Laurent 展开，收敛范围都是全复平面（孤立奇点处除外），也不是某一圆域或环域。在这种展开中，能把函数在它的全部极点的奇异性同时表现无遗。能够作有理分式展开的，仅限于在有限区域内只有极点的单值函数，即亚纯函数。

将 (7.43) 和 (7.44) 两式再求微商，还可以进一步得到

$$
\pi^2\sec^2\pi z
=
4\sum_{m=-\infty}^{\infty}
\frac{1}{[2z-(2m-1)]^2},
\qquad
z\ne\pm\frac12,\pm\frac32,\cdots,
\tag{7.47}
$$

$$
\pi^2\csc^2\pi z
=
\sum_{m=-\infty}^{\infty}
\frac{1}{(z-m)^2},
\qquad
z\ne0,\pm1,\pm2,\cdots.
\tag{7.48}
$$

将 (7.43)～(7.46) 诸式右端在 $z=0$ 点作 Taylor 展开，并交换求和次序，又可以得到

$$
\pi z\tan\pi z
=
2\sum_{n=1}^\infty
\left[
\sum_{m=1}^\infty
\frac{1}{(2m-1)^{2n}}
\right]
(2z)^{2n},
$$

$$
\pi z\cot\pi z
=
1-
2\sum_{n=1}^\infty
\left(
\sum_{m=1}^\infty
\frac{1}{m^{2n}}
\right)
z^{2n},
$$

$$
\pi\sec\pi z
=
4\sum_{n=0}^\infty
\left[
\sum_{m=1}^\infty
\frac{(-)^{m-1}}
{(2m-1)^{2n+1}}
\right]
(2z)^{2n},
$$

$$
\pi z\csc\pi z
=
1+
2\sum_{n=1}^\infty
\left[
\sum_{m=1}^\infty
\frac{(-)^{m-1}}
{m^{2n}}
\right]
z^{2n}.
$$

再与 (5.44)～(5.46) 诸式及 (5.51) 式相比较，就可以得出下列几个级数和：

$$
\sum_{m=1}^\infty
\frac{1}{m^{2n}}
=
\frac{(-)^{n-1}}{2}
\frac{(2\pi)^{2n}}{(2n)!}
B_{2n},
\tag{7.49}
$$

$$
\sum_{m=1}^\infty
\frac{(-1)^{m-1}}{m^{2n}}
=
(-)^{n-1}
\frac{2^{2n-1}-1}{(2n)!}
\pi^{2n}B_{2n},
\tag{7.50}
$$

$$
\sum_{m=1}^\infty
\frac{1}{(2m-1)^{2n}}
=
\frac{(-)^{n-1}}{2}
\frac{2^{2n}-1}{(2n)!}
\pi^{2n}B_{2n},
\tag{7.51}
$$

$$
\sum_{m=1}^\infty
\frac{(-1)^{m-1}}{(2m-1)^{2n+1}}
=
\frac{(-)^n}{2^{2n+2}}
\frac{\pi^{2n+1}}{(2n)!}
E_{2n}.
\tag{7.52}
$$

将 (7.49) 式改写成

$$
B_{2n}
=
(-)^{n-1}
\frac{2(2n)!}{(2\pi)^{2n}}
\sum_{m=1}^\infty
\frac{1}{m^{2n}},
\tag{7.49'}
$$

就可以看出，当 $n\to\infty$ 时 $B_{2n}\to\infty$。例如，

$$
B_{20}=-5.291\times10^2,
\qquad
B_{30}=6.016\times10^8,
\qquad
B_{40}=-1.930\times10^{16},
$$

$$
B_{50}=7.501\times10^{24},
\qquad
\cdots,
\qquad
B_{200}=-3.647\times10^{215}.
$$

## *§7.6 Γ 函数的渐近展开

§7.2 中曾经介绍了 $z\to\infty$ 时 $\Gamma(z)$ 的渐近展开。这一节就 $z$ 为实数 $x$，即

$$
\Gamma(x+1)
=
\int_0^\infty e^{-t}t^x\,dx,
\qquad
x>0
\tag{7.53}
$$

的情形，推导一下这个公式①。分析一下积分 (7.53) 的被积函数，它在 $t=0$ 时为 0，随着 $t$ 的增大而增大，当 $t=x$ 时达到极大，而后又单调下降。由于指数函数的变化特点，被积函数对积分的贡献主要来自 $t=x$ 附近的一个很窄的区间。这时，可以将被积函数写成

$$
e^{-t}t^x
=
e^{-t+x\ln t},
$$

再将函数 $-t+x\ln t$ 在 $t=x$ 点作 Taylor 展开：

$$
-t+x\ln t
=
x\ln x-x
-\frac{(t-x)^2}{2x}
+\cdots.
\tag{7.54}
$$

所以，就得到 Γ 函数的近似表达式

$$
\begin{aligned}
\Gamma(x+1)
&\approx
x^xe^{-x}
\int_0^\infty
e^{-(t-x)^2/2x}\,dt\\
&\approx
x^xe^{-x}
\int_{-\infty}^{\infty}
e^{-\xi^2/2x}\,d\xi\\
&=
x^xe^{-x}\sqrt{2\pi x}.
\end{aligned}
\tag{7.55}
$$

这样，可以预料，$\Gamma(x+1)$ 的渐近展开式应该为

$$
\Gamma(x+1)
=
x^xe^{-x}\sqrt{2\pi x}
\left(
1+\frac Ax+\frac B{x^2}+\frac C{x^3}+\cdots
\right).
\tag{7.56}
$$

问题是如何确定系数 $A,B,C,\cdots$。巧妙的办法是将 (7.56) 式代入 Γ 函数的递推关系

$$
\Gamma(x+1)=x\Gamma(x),
$$

从而得到

$$
1+\frac Ax+\frac B{x^2}+\frac C{x^3}+\cdots
=
e\left(1-\frac1x\right)^{x-1/2}
\left[
1+\frac A{x-1}
+\frac B{(x-1)^2}
+\frac C{(x-1)^3}
+\cdots
\right].
\tag{7.57}
$$

由于

$$
\left(1-\frac1x\right)^{x-1/2}
=
\exp
\left\{
\left(x-\frac12\right)
\ln\left(1-\frac1x\right)
\right\},
$$

当 $x$ 很大时，

$$
\begin{aligned}
\left(x-\frac12\right)
\ln\left(1-\frac1x\right)
&=
\left(x-\frac12\right)
\left(
-\frac1x
-\frac1{2x^2}
-\frac1{3x^3}
-\frac1{4x^4}
-\cdots
\right)\\
&=
-1
-\frac1{12x^2}
-\frac1{12x^3}
-\cdots,
\end{aligned}
$$

所以

$$
\left(1-\frac1x\right)^{x-1/2}
=
\exp
\left\{
-1
-\frac1{12x^2}
-\frac1{12x^3}
-\cdots
\right\}
=
e^{-1}
\left(
1-\frac1{12x^2}
-\frac1{12x^3}
+\cdots
\right).
$$

同时，将 $(x-1)^{-n}$ 也作展开，

$$
\frac1{x-1}
=
\frac1x+\frac1{x^2}+\frac1{x^3}+\cdots,
$$

$$
\frac1{(x-1)^2}
=
\frac1{x^2}+\frac2{x^3}+\cdots,
$$

$$
\frac1{(x-1)^3}
=
\frac1{x^3}+\cdots,
\qquad\cdots.
\tag{7.58}
$$

代入 (7.57) 式，有

$$
\begin{aligned}
1+\frac Ax+\frac B{x^2}+\frac C{x^3}+\cdots
&=
\left(
1-\frac1{12x^2}
-\frac1{12x^3}
+\cdots
\right)
\left(
1+\frac Ax+\frac A{x^2}+\frac A{x^3}
+\frac B{x^2}
+\frac{2B}{x^3}
+\frac C{x^3}
+\cdots
\right)\\
&=
1+\frac Ax
+\left(A+B-\frac1{12}\right)\frac1{x^2}
+\left(
\frac{11}{12}A
+2B+C
-\frac1{12}
\right)\frac1{x^3}
+\cdots.
\end{aligned}
\tag{7.59}
$$

① 本节的推导方法引自 T. C. Bradbury, *Mathematical Methods with Applications to Problems in the Physical Sciences*, John Wiley & Sons, Inc., New York, 1984。基本思想仍是鞍点法，只不过由于限于实数情形，因而比较简单。

比较系数，就求得

$$
A=\frac1{12},
\qquad
B=\frac1{288}.
\tag{7.60}
$$

当然，如果在 (7.57)～(7.57) 式中，写出更高次负幂项，就可以定出系数 $C,D,\cdots$。

## 习题

1. 将下列连乘积用 Γ 函数表示出来：

   (1) $(2n)!!$；

   (2) $(2n-1)!!$；

   (3) $(1+\nu)(2+\nu)(3+\nu)\cdots(n+\nu)$；

   (4) $[n(n+1)-\nu(\nu+1)][(n-1)n-\nu(\nu+1)]\cdots[0-\nu(\nu+1)]$。

2. 计算下列积分：

   (1)

   $$
   \int_0^\infty x^{-\alpha}\sin x\,dx,
   \qquad
   0<\alpha<2;
   $$

   $$
   \int_0^\infty x^{-\alpha}\cos x\,dx,
   \qquad
   0<\alpha<1.
   $$

   (2)

   $$
   \int_0^\infty
   x^{\alpha-1}
   e^{-x\cos\theta}
   \cos(x\sin\theta)\,dx,
   \qquad
   \alpha>0,
   \quad
   -\frac\pi2<\theta<\frac\pi2;
   $$

   $$
   \int_0^\infty
   x^{\alpha-1}
   e^{-x\cos\theta}
   \sin(x\sin\theta)\,dx,
   \qquad
   \alpha>0,
   \quad
   -\frac\pi2<\theta<\frac\pi2.
   $$

3. 设

$$
\psi(z)
=
\frac{d}{dz}\ln\Gamma(z)
=
\frac{\Gamma'(z)}{\Gamma(z)},
$$

证明：

   (1)

   $$
   \psi(z+1)=\frac1z+\psi(z);
   $$

   (2)

   $$
   \psi(z+n)-\psi(z)
   =
   \frac1z+\frac1{z+1}+\cdots+\frac1{z+n-1};
   $$

   (3)

   $$
   \psi(1-z)-\psi(z)
   =
   \pi\cot\pi z;
   $$

   (4)

   $$
   2\psi(2z)
   -
   \psi(z)
   -
   \psi\left(z+\frac12\right)
   =
   2\ln2.
   $$

4. 计算下列积分：

   (1)

   $$
   \int_{-1}^{1}
   (1-x)^p(1+x)^q\,dx,
   \qquad
   \operatorname{Re}p>-1,\ 
   \operatorname{Re}q>-1;
   $$

   (2)

   $$
   \int_0^{\pi/2}\tan^\alpha\theta\,d\theta,
   \qquad
   \int_0^{\pi/2}\cot^\alpha\theta\,d\theta,
   \qquad
   -1<\alpha<1.
   $$

5. 计算积分

$$
\iiint_V
x^{\alpha-1}
y^{\beta-1}
z^{\gamma-1}\,dx\,dy\,dz
\qquad
(\alpha,\beta,\gamma\text{ 均为正数}),
$$

其中积分区域 $V$ 为

   (1) 平面 $x=0,\ y=0,\ z=0$ 及 $x+y+z=1$ 所包围的体积；

   (2) 平面 $x=0,\ y=0,\ z=0$ 及曲面 $x^p+y^q+z^r=1$ 所包围的体积，$p,q,r$ 均为已知正数。

6. 求下列无穷级数之和：

   (1)

   $$
   \sum_{n=1}^{\infty}
   \frac{1}{n(4n^2-1)};
   $$

   (2)

   $$
   \sum_{n=-\infty}^{\infty}
   \frac{1}{(n^2+1)^2}.
   $$
