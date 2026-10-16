# 第二章 解析函数

## §2.1 复变函数的极限和连续

设函数 $f(z)$ 在 $z_0$ 点的空心邻域①内有定义，若存在复数 $A$，$\forall\varepsilon>0$，$\exists\delta(\varepsilon)>0$，使当 $0<|z-z_0|<\delta$ 时，恒有 $|f(z)-A|<\varepsilon$，则称 $z\to z_0$ 时 $f(z)$ 的极限（$=A$）存在，表示为

$$
\lim_{z\to z_0}f(z)=A.
$$

设函数 $f(z)$ 在 $z_0$ 点的邻域②内有定义，且 $f(z)$ 在 $z_0$ 点的极限值等于在该点的函数值 $f(z_0)$，记为 $\lim_{z\to z_0}f(z)=f(z_0)$，即 $\forall\varepsilon>0$，$\exists\delta(\varepsilon,z_0)>0$，使当 $|z-z_0|<\delta$ 时，恒有 $|f(z)-f(z_0)|<\varepsilon$，则称 $f(z)$ 在 $z_0$ 点连续. 由定义可知，函数的连续性是逐点定义的，是函数在复平面内某一点的性质.

复变函数中极限和连续概念的表述，形式上和实变函数中完全相同，但由于涉及的数域不同，因此实际含义并不完全相同.

若函数 $f$ 在区域 $G$ 内每一点都连续，则称 $f$ 为 $G$ 内的连续函数. 连续函数的和、差、积、商（在分母不为零的点）仍为连续函数，连续函数的复合函数也仍为连续函数.

在有界闭区域 $\overline G$ 中连续的函数 $f(z)$ 具有两个重要性质：

1. $|f(z)|$ 在 $\overline G$ 中有界，并达到它的上下界；
2. $f(z)$ 在 $\overline G$ 中一致连续，即 $\forall\varepsilon>0$，$\exists\delta(\varepsilon)>0$，使得对 $\forall z_1\in\overline G$，$\forall z_2\in\overline G$，只要满足 $|z_1-z_2|<\delta$，就有 $|f(z_1)-f(z_2)|<\varepsilon$.

由一致连续的定义可知，函数 $f(z)$ 的一致连续性是和区域相关联的，不能脱离区域谈一个函数的一致连续性. 一致连续性比上面讨论的连续性对函数的约束更强.

## §2.2 可导与可微

设 $w=f(z)$ 是区域 $G$ 内的单值函数，如果在 $G$ 内的某点 $z$，

$$
\lim_{\Delta z\to0}\frac{\Delta w}{\Delta z}
=\lim_{\Delta z\to0}\frac{f(z+\Delta z)-f(z)}{\Delta z}
\tag{2.1}
$$

存在，则称函数 $f(z)$ 在 $z$ 点可导，此极限值即称为 $f(z)$ 在 $z$ 点的导数，记为 $f'(z)$.

导数的定义在形式上和实数中一样，只是把实自变量换成了复自变量，因此高等数学中的各种求导数的公式都可搬用到复变函数中来. 例如

$$
(z^n)'=nz^{n-1},\qquad n=0,1,2,\cdots.
$$

① 所谓 $z_0$ 点的空心邻域，指的是以 $z_0$ 点为圆心的环域 $0<|z-z_0|<\varepsilon$.

② $z_0$ 点的 $\varepsilon$ 邻域（简称邻域），则是以 $z_0$ 点为圆心的圆域 $|z-z_0|<\varepsilon$；有时也称作 $z_0$ 点及其邻域，以示强调.

需要强调，上面所说的 $\lim_{\Delta z\to0}(\Delta w/\Delta z)$ 存在，就意味着 $\Delta z$ 以任意方式趋于 0 时，$\Delta w/\Delta z$ 都趋于同样的有限值. 反过来说，如果当 $\Delta z$ 以不同方式趋于 0，$\Delta w/\Delta z$ 趋于不同的值的话，则 $\lim_{\Delta z\to0}(\Delta w/\Delta z)$ 是不存在的.

特别是，考虑 $\Delta z\to0$ 的不同特殊方式，就可以得到关于函数可导的必要条件. 例如：

(1) $\Delta x\to0,\ \Delta y=0$，

$$
f'(z)=\lim_{\Delta z\to0}\frac{\Delta w}{\Delta z}
=\lim_{\Delta x\to0}\frac{\Delta u+i\Delta v}{\Delta x}
=\frac{\partial u}{\partial x}+i\frac{\partial v}{\partial x};
$$

(2) $\Delta x=0,\ \Delta y\to0$，

$$
f'(z)=\lim_{\Delta z\to0}\frac{\Delta w}{\Delta z}
=\lim_{\Delta y\to0}\frac{\Delta u+i\Delta v}{i\Delta y}
=\frac{\partial v}{\partial y}-i\frac{\partial u}{\partial y}.
$$

因此，函数实部和虚部在 $z$ 点的 4 个偏导数值必须满足 Cauchy-Riemann 条件：

$$
\boxed{
\frac{\partial u}{\partial x}=\frac{\partial v}{\partial y},\qquad
\frac{\partial u}{\partial y}=-\frac{\partial v}{\partial x}
}
\tag{2.2}
$$

Cauchy-Riemann 条件是函数可导的必要条件，但不是充分条件. 它保证了当 $\Delta z$ 以平行于实轴和虚轴这两种特殊方式趋于 0 时 $\Delta w/\Delta z$ 逼近同一值，但并不足以保证当 $\Delta z$ 以任意方式趋于 0 时 $\Delta w/\Delta z$ 逼近同一数值. 可以证明，如果函数 $f(z)=u(x,y)+iv(x,y)$ 的实部 $u(x,y)$ 和虚部 $v(x,y)$ 均在 $(x,y)$ 点可微①，且满足 Cauchy-Riemann 条件，则函数 $f(z)$ 在 $z$ 点可导.

和实数情形一样，函数可导是比函数连续更强的条件. 如果函数 $f(z)$ 在 $z$ 点可导，则在 $z$ 点必连续；但是函数 $f(z)$ 在 $z$ 点连续，并不能推出函数 $f(z)$ 在 $z$ 点可导. 甚至有函数在某区域内处处连续，却处处不可导.

如果函数 $w=f(z)$ 在 $z$ 点函数值的改变量 $\Delta w=f(z+\Delta z)-f(z)$ 可以写成

$$
\Delta w=A(z)\Delta z+\rho(\Delta z),\qquad
\text{其中}\quad \lim_{\Delta z\to0}\frac{\rho(\Delta z)}{\Delta z}=0,
$$

则称 $w=f(z)$ 在 $z$ 点可微，$\Delta w$ 的线性部分 $A(z)\Delta z$ 称为函数 $w$ 在 $z$ 点的微分，记作

$$
dw=A(z)\Delta z.
$$

可以证明，如果函数 $w=f(z)$ 在 $z_0$ 点可导，则一定在该点可微，反之亦然，并且 $A(z_0)=f'(z_0)$. 由于对函数 $f(z)=z$，$f'(z)=1$，从而 $dz=\Delta z$，因此有

$$
dw=A(z_0)dz,
$$

即

$$
dw=f'(z_0)dz\qquad\text{或}\qquad
\left.\frac{dw}{dz}\right|_{z=z_0}=f'(z_0).
\tag{2.3}
$$

① $u(x,y)$ 和 $v(x,y)$ 均可微的充分条件是：四个偏导数 $\dfrac{\partial u}{\partial x}$、$\dfrac{\partial u}{\partial y}$ 和 $\dfrac{\partial v}{\partial x}$、$\dfrac{\partial v}{\partial y}$ 存在且连续.

因此导数也称作微商.

**导数的几何意义**　设 $w=f(z)$ 在 $z_0$ 点可导，则由 (2.3) 式可以看出（见图 2.1），$f'(z_0)$ 的模 $|f'(z_0)|$ 是将微元 $dz$ 映射为 $dw$ 的伸缩率（放大倍数），而 $f'(z_0)$ 的辐角 $\arg f'(z_0)$ 则给出映射的偏转角（$dw$ 与 $dz$ 的辐角差）.

![Figure 2.1](...)

图 2.1　导数的几何意义

$$
|dw|=|f'(z_0)|\,|dz|,\qquad \arg dw=\arg f'(z_0)+\arg dz.
$$

## §2.3 解析函数

在区域 $G$ 内每一点都可导的函数，称为 $G$ 内的解析函数，或者说函数在 $G$ 内解析. 例如，$z^2$ 就是 $\mathbb C$ 内的解析函数.

显然，函数在 $G$ 内解析的必要条件是在 $G$ 内处处满足 Cauchy-Riemann 条件，或者说，在 $G$ 内 Cauchy-Riemann 方程成立. Cauchy-Riemann 方程反映了解析函数的实部与虚部之间的联系. 更准确地说，解析函数的实部和虚部不是相互独立的，知道其中之一，例如实部 $u(x,y)$，就可以唯一地（可相差一个可加实常数）确定其虚部. 这是因为，根据 Cauchy-Riemann 方程，可以求出虚部 $v(x,y)$ 的全微分

$$
dv=\frac{\partial v}{\partial x}dx+\frac{\partial v}{\partial y}dy
=-\frac{\partial u}{\partial y}dx+\frac{\partial u}{\partial x}dy,
$$

因此，通过二维平面内的线积分

$$
\int^{(x,y)}\left(-\frac{\partial u}{\partial y}dx+\frac{\partial u}{\partial x}dy\right),
$$

可以完全确定 $v(x,y)$ 的函数形式，除了一个任意的实常数.

**例 2.1**　已知 $u(x,y)=x^2-y^2$，求 $f(z)$.

**解**　

$$
dv=-\frac{\partial u}{\partial y}dx+\frac{\partial u}{\partial x}dy
=2(ydx+xdy),
$$

所以 $v=2xy+C$，其中 $C$ 为任意实常数.

$$
f(z)=(x^2-y^2)+i(2xy+C)=z^2+iC.
$$

这个问题，还可以有另一种解法，即在 $u(x,y)$ 中直接代入

$$
x=\frac{z+z^*}{2},\qquad y=\frac{z-z^*}{2i},
$$

而后就能将 $u(x,y)$ 化成 $[f(z)+f^*(z)]/2$ 的形式，即

$$
u(x,y)=\left(\frac{z+z^*}{2}\right)^2-
\left(\frac{z-z^*}{2i}\right)^2
=\frac12\left[z^2+(z^2)^*\right].
$$

再经过甄别，弃去不满足 Cauchy-Riemann 方程的函数 $(z^2)^*=(z^*)^2$，同样也能求出 $f(z)=z^2+iC$.

解析函数的实部和虚部之间的这种依赖关系，还可以形象化地表现出来. 如果在 $x-y$ 平面内作一族曲线，$u(x,y)=$ 常数，那么这族曲线的切线的方向矢量便是 $(\partial u/\partial y,-\partial u/\partial x)$. 同样，再作一族曲线，$v(x,y)=$ 常数，它们的切线的方向矢量是 $(\partial v/\partial y,-\partial v/\partial x)$. 由 Cauchy-Riemann 方程，可以求得这两族方向矢量之间的标积

$$
\left(\frac{\partial u}{\partial y},-\frac{\partial u}{\partial x}\right)
\left(\frac{\partial v}{\partial y},-\frac{\partial v}{\partial x}\right)
=\frac{\partial u}{\partial y}\frac{\partial v}{\partial y}
+\frac{\partial u}{\partial x}\frac{\partial v}{\partial x}=0.
\tag{2.4}
$$

![Figure 2.2](...)

图 2.2

这表明，这两族曲线是互相正交的. 图 2.2 中给出了两个例子. 它们分别是函数 $w=z^2$ 和 $w=1/z^2$. 图中的粗黑线表示实部 $u(x,y)=$ 常数，细灰线表示虚部 $v(x,y)=$ 常数.

**练习 2.1**　证明：Cauchy-Riemann 方程等价于

$$
i\frac{\partial f}{\partial x}=\frac{\partial f}{\partial y}.
$$

**练习 2.2**　如果把复变函数 $f=u+iv$ 看成是 $(x,y)$ 的二元函数，即

$$
f(x,y)=u(x,y)+iv(x,y),
$$

再进一步看成是 $z=x+iy$ 和 $z^*=x-iy$ 的二元函数，证明 Cauchy-Riemann 方程等价于

$$
\frac{\partial f}{\partial z^*}=0.
$$

**练习 2.3**　证明：

$$
|f'(z)|^2
=\left(\frac{\partial u}{\partial x}\right)^2+
\left(\frac{\partial v}{\partial x}\right)^2
=\left(\frac{\partial u}{\partial y}\right)^2+
\left(\frac{\partial v}{\partial y}\right)^2
$$

$$
=\left(\frac{\partial u}{\partial x}\right)^2+
\left(\frac{\partial u}{\partial y}\right)^2
=\left(\frac{\partial v}{\partial x}\right)^2+
\left(\frac{\partial v}{\partial y}\right)^2.
$$

我们知道，复变函数 $w=f(z)=u(x,y)+iv(x,y)$ 的实部 $u(x,y)$ 和虚部 $v(x,y)$ 都是二元实函数. 那么是不是任意一个二元函数都可以用来作为解析函数的实部或虚部呢？回答是否定的. §3.5 中将证明，解析函数的实部 $u(x,y)$ 和虚部 $v(x,y)$ 的二阶偏导数一定存在并且连续，因此，根据 Cauchy-Riemann 方程，有

$$
\frac{\partial^2u}{\partial x^2}
=\frac{\partial}{\partial x}\frac{\partial v}{\partial y}
=\frac{\partial^2v}{\partial x\partial y},
\qquad
\frac{\partial^2u}{\partial y^2}
=\frac{\partial}{\partial y}\left(-\frac{\partial v}{\partial x}\right)
=-\frac{\partial^2v}{\partial x\partial y},
$$

$$
\frac{\partial^2v}{\partial x^2}
=\frac{\partial}{\partial x}\left(-\frac{\partial u}{\partial y}\right)
=-\frac{\partial^2u}{\partial x\partial y},
\qquad
\frac{\partial^2v}{\partial y^2}
=\frac{\partial}{\partial y}\frac{\partial u}{\partial x}
=\frac{\partial^2u}{\partial x\partial y}.
$$

这说明，$u(x,y)$ 和 $v(x,y)$ 都必须满足二维 Laplace 方程

$$
\boxed{
\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}=0,
\qquad
\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}=0.
}
\tag{2.5}
$$

即解析函数的实部和虚部都必须是调和函数. 而且，因为一个解析函数的实部和虚部必须受到 Cauchy-Riemann 方程的制约，所以，解析函数的实部和虚部就构成一对共轭调和函数.

函数的解析性，总是和一定的区域联系在一起的. 有时也称函数在 $z_0$ 点解析，这应理解为存在 $z_0$ 点的一个邻域，函数在 $z_0$ 点以及 $z_0$ 的这个邻域内处处可导. 如果要讨论函数 $f(z)$ 在 $z=\infty$ 点是否解析，则需作变换 $z=1/t$，然后讨论函数 $f(1/t)$ 在 $t=0$ 点是否解析.

函数的解析性，对函数是一个高要求，这表现为解析函数具有一系列的重要性质. 讨论解析函数的各种特殊性质，就是复变函数论的中心课题.

**练习 2.4**　证明：

$$
\frac{d}{dz}[f(z)+g(z)]=\frac{df(z)}{dz}+\frac{dg(z)}{dz};
$$

$$
\frac{d}{dz}[f(z)g(z)]=\frac{df(z)}{dz}g(z)+f(z)\frac{dg(z)}{dz};
$$

$$
\frac{d}{dz}\frac{f(z)}{g(z)}
=\frac{f'(z)g(z)-f(z)g'(z)}{g^2(z)},\qquad g(z)\ne0;
$$

$$
\frac{d}{dz}f(g(z))=f'(g(z))g'(z).
$$

**练习 2.5**　举例说明中值定理不适用于解析函数：若函数 $f(z)$ 在 $G$ 中解析，$z_1$ 和 $z_2$ 以及连接两点的线段均在 $G$ 中，在此线段上不一定存在 $z_0$ 点，使得

$$
\frac{f(z_1)-f(z_2)}{z_1-z_2}=f'(z_0).
$$

**练习 2.6**　假设函数 $f(z)$ 在区域 $G$ 内的任何一点都满足 $f'(z)=0$，证明 $f(z)$ 在 $G$ 内为常数.

**练习 2.7**　若函数 $f(z)$ 在区域 $G$ 内解析，且 $\operatorname{Im}f(z)=0$，证明 $f(z)$ 在 $G$ 内为常数.

**练习 2.8**　若函数 $f(z)=u(x,y)+iv(x,y)$ 在区域 $G$ 内解析，且 $au(x,y)+bv(x,y)=c$，$a,b$ 和 $c$ 是不为 0 的实常数，证明 $f(z)$ 必为常数.

如果 $a,b$ 和 $c$ 是不为 0 的复常数，这个结论还成立吗？

## §2.4 初等函数

本节介绍一些基本的解析函数，例如幂函数 $z^n$；指数函数 $e^z$；三角函数 $\sin z,\cos z,\cdots$；双曲函数 $\sinh z,\cosh z,\cdots$；等等. 它们都可以看成是相应实变函数在复数域中的推广. 这里将着重讨论这些函数作为复变函数所特有的那些性质.

### 1. 幂函数 $z^n$

当 $n=0,1,2,\cdots$ 时，$z^n$ 在 $\mathbb C$ 内解析；并且当 $n=1,2,\cdots$ 时，$z^n$ 在 $z=\infty$ 不解析.

当 $n=-1,-2,\cdots$ 时，$z^n$ 在 $z=0$ 不解析，在包括 $\infty$ 点在内的 $\overline{\mathbb C}\setminus\{0\}$ 内处处解析.

在 $z^n$ 的解析区域内

$$
(z^n)'=nz^{n-1}.
$$

由幂函数还可以进一步定义（$n$ 次，$n=0,1,2,\cdots$）多项式（函数）

$$
P_n(z)=a_nz^n+a_{n-1}z^{n-1}+\cdots+a_1z+a_0
$$

和有理函数

$$
R(z)=\frac{P_n(z)}{Q_m(z)},\qquad Q_m\ne0,
$$

其中 $P_n(z)$ 和 $Q_m(z)$ 分别是 $n$ 次和 $m$ 次多项式，$n$ 和 $m$ 都是非负整数.

### 2. 指数函数 $e^z$

$$
e^z=e^{x+iy}=e^x(\cos y+i\sin y).
$$

由实指数函数及纯虚数指数函数的性质容易看出，“指数函数相乘等于指数相加”这个法则，对于复指数函数仍然成立：

$$
e^{z_1}\cdot e^{z_2}
=e^{x_1+iy_1}\cdot e^{x_2+iy_2}
=e^{x_1+x_2}\cdot e^{i(y_1+y_2)}
=e^{(x_1+x_2)+i(y_1+y_2)}
=e^{z_1+z_2}.
$$

$e^z$ 在 $\mathbb C$ 内解析，

$$
(e^z)'=e^z.
$$

但 $e^z$ 在无穷远点无定义，当然也不解析. 例如，当 $z$ 沿正实轴或负实轴趋于 $\infty$ 时，$e^z$ 逼近不同的数值.

复指数函数特有而实指数函数不具备的一个性质是周期性，周期为 $2\pi i$：

$$
e^{z+2\pi i}=e^{x+i(y+2\pi)}
=e^x[\cos(y+2\pi)+i\sin(y+2\pi)]
$$

$$
=e^x(\cos y+i\sin y)=e^{x+iy}=e^z.
$$

**练习 2.9**　如果 $z$ 沿不同辐角方向趋于 $\infty$ 点，试讨论函数 $e^z$ 的变化趋势.

又设常数 $\alpha\ne0$，试设计一个无穷序列 $\{z_n\}$，使 $n\to\infty$ 时，函数序列 $\{e^{z_n}\}$ 趋于 $\alpha$.

### 3. 三角函数 $\sin z,\cos z,\cdots$

复三角函数 $\sin z,\cos z$ 可以用复指数函数定义：

$$
\sin z=\frac{e^{iz}-e^{-iz}}{2i},\qquad
\cos z=\frac{e^{iz}+e^{-iz}}{2}.
\tag{2.6}
$$

由于 $e^{iz}$ 与 $e^{-iz}$ 在 $\mathbb C$ 内解析，所以 $\sin z,\cos z$ 也在 $\mathbb C$ 内解析，

$$
(\sin z)'=\cos z,\qquad (\cos z)'=-\sin z.
$$

$z=\infty$ 是它们的唯一奇点.

和实三角函数一样，$\sin z$ 和 $\cos z$ 也都是周期函数，周期为 $2\pi$.

和实三角函数不同，$\sin z$ 和 $\cos z$ 的模可以大于 1. 例如，

$$
i\sin i=\frac{e^{-1}-e^1}{2}=-1.1752012\cdots,
\qquad
\cos i=\frac{e^{-1}+e^1}{2}=1.5430806\cdots.
$$

其他三角函数，$\tan z,\cot z,\sec z,\csc z$ 可以用 $\sin z$ 和 $\cos z$ 定义，形式和实数情形一样，

$$
\tan z=\frac{\sin z}{\cos z},\qquad
\cot z=\frac{\cos z}{\sin z},\qquad
\sec z=\frac{1}{\cos z},\qquad
\csc z=\frac{1}{\sin z}.
$$

根据这些定义可以证明，实三角函数的各种恒等式对于复三角函数仍然成立.

### 4. 双曲函数 $\sinh z,\cosh z,\cdots$

双曲函数 $\sinh z,\cosh z$ 也是通过复指数函数定义的.

$$
\sinh z=\frac{e^z-e^{-z}}{2},\qquad
\cosh z=\frac{e^z+e^{-z}}{2},\qquad
\tanh z=\frac{\sinh z}{\cosh z},
$$

$$
\coth z=\frac{\cosh z}{\sinh z},\qquad
\operatorname{sech}z=\frac{1}{\cosh z},\qquad
\operatorname{csch}z=\frac{1}{\sinh z}.
\tag{2.7}
$$

由定义可以直接证明，双曲函数和三角函数可以互化，

$$
\sinh z=-i\sin iz,\qquad
\cosh z=\cos iz,\qquad
\tanh z=-i\tan iz.
$$

因此，双曲函数的性质完全可以由三角函数推出. 这里只想特别指出两点：一是周期性，双曲函数 $\sinh z,\cosh z,\operatorname{sech}z$ 和 $\operatorname{csch}z$ 的周期是 $2\pi i$，$\tanh z$ 和 $\coth z$ 的周期是 $\pi i$；二是导数公式

$$
(\sinh z)'=\cosh z,\qquad
(\cosh z)'=\sinh z,\qquad
(\tanh z)'=\operatorname{sech}^2z.
\tag{2.8}
$$

**练习 2.10**　证明下列公式：

$$
\cosh^2z-\sinh^2z=1;\qquad 1-\tanh^2z=\operatorname{sech}^2z;
$$

$$
|\sinh y|\le |\sin(x+iy)|\le\cosh y;\qquad
|\sinh y|\le |\cos(x+iy)|\le\cosh y;
$$

$$
\sinh(z_1\pm z_2)=\sinh z_1\cosh z_2\pm\cosh z_1\sinh z_2;
$$

$$
\cosh(z_1\pm z_2)=\cosh z_1\cosh z_2\pm\sinh z_1\sinh z_2.
$$

## \*§2.5 解析函数的保角性

函数的解析性是研究复变函数的微积分学性质的最基本前提. 这一看似简单的要求，实际上对函数施加了相当强的限制. 这从 Cauchy-Riemann 方程已经可以看出一些端倪. 从解析性出发，可以导出函数的其他许多重要特性. 复变函数理论的主要研究对象，就是解析函数. 以后各章将逐步介绍解析函数的方方面面. 这一节，先侧重于从几何的角度简要讨论解析函数的应用.

复变函数代表了一个变换，或称映射. 在映射 $\zeta=f(z)$ 之下，$z$ 复平面内的一点映射为 $\zeta$ 复平面内的相应一点. 如果函数 $\zeta=f(z)$ 是连续的，$z$ 点邻域内的一点当然也应该映射为相应的 $\zeta$ 点邻域内的一点. 但是，$z$ 复平面内的一个区域，是否也映射为 $\zeta$ 复平面内的一个区域，区域的边界是否仍映射为区域的边界，并不是任何复变函数都能够保证做到的. 可以举一个极端的例子：函数 $\zeta=\operatorname{Re}z$ 就把整个 $z$ 复平面映射为 $\zeta$ 复平面内的实轴，后者甚至不构成一个区域.

设 $\zeta=f(z)$ 在区域 $G$ 内解析，$z_0$ 为 $G$ 内一点. $\zeta=f(z)$ 把 $z_0$ 点映射为相应的 $\zeta_0=f(z_0)$ 点. 在 §2.2 中曾经提到，当 $f'(z_0)\ne0$ 时，$|f'(z_0)|$ 是把 $z_0$ 处的微元 $dz$ 映射为 $\zeta_0$ 处的微元 $d\zeta$ 的伸缩率（放大倍数），而 $\arg f'(z_0)$（这里不妨限制在主值范围内）则给出映射的偏转角，即 $d\zeta$ 与 $dz$ 的辐角差. 可以设想，若在 $z$ 复平面内有两条曲线 $l_1$ 和 $l_2$ 相交于 $z_0$，

![Figure 2.3](...)

图 2.3　解析函数所代表的变换的保角性

那么，在 $\zeta$ 复平面内，必然也有相应的两条曲线 $l_1'$ 和 $l_2'$ 相交于 $\zeta_0$ 点. 由于

$$
f'(z_0)=\lim_{z\to z_0}\frac{f(z)-f(z_0)}{z-z_0}
$$

的数值，包括 $|f'(z_0)|$ 和 $\arg f'(z_0)$，与 $z$ 趋于 $z_0$ 的方式无关，因此，与 $l_1$ 和 $l_2$ 相比，$l_1'$ 和 $l_2'$ 应该放大或缩小了相同倍数，并且偏转了同样大小的角度. 这就是说，在解析函数所代表的变换之下，过同一点的两条曲线，它们的伸缩率相同，并且夹角保持不变（见图 2.3），正是由于这个原因，解析函数所代表的变换（映射）就称为保角变换（保角映射）. 当然，在不同点处，即使 $f'(z)$ 均不为 0，变换都具有保角性，但由于数值不同，因而各处的伸缩率和偏转角不同，所以，区域的几何形状就会发生变化. 我们正是要选择合适的保角变换，把 $z$ 复平面内形状比较复杂的区域变换为 $\zeta$ 复平面内形状比较简单的区域. 例如，把 $z$ 复平面内两个不相交的偏心圆周所围成的区域变换为 $\zeta$ 复平面内同心圆周所围成的区域.

解析函数所代表的变换的保角性，是有条件的：只在 $f'(z)\ne0$ 处才一定有保角性. 在 $f'(z)=0$ 的点，由于 $\arg f'(z)$ 没有确定值，因而变换可能是保角的，也可能是不保角的. 巧妙地利用变换在 $f'(z)=0$ 处的不保角性，更可以把（$z$ 复平面内的）复杂图形变换为（$\zeta$ 复平面内的）简单图形. 例如，可以把多边形变为圆周. Riemann 映射定理告诉我们，如果 $C$ 是 $z$ 平面上有界单连通区域 $G$ 的边界，$C'$ 是 $\zeta$ 平面上单连通区域 $G'$ 的边界，则一定存在一个 $G$ 内的解析函数 $\zeta=f(z)$，使得 $G$ 内的每个点各自变换成 $G'$ 内一个对应点，使得 $C$ 上的每个点各自变换成 $C'$ 的一个对应点，而且这样的变换是双向单值对应的.

将保角变换应用于解决实际问题，不仅涉及区域形状的变换，而且还要涉及数学表述形式（例如微分方程）的变换. 可以想象，只有在既让区域的形状变得很简单，也没有让微分方程的形式变得更复杂的条件下，保角变换才具有真正的实用价值.

就本课程而言，我们特别有兴趣于讨论二维 Laplace 算符

$$
\nabla^2=\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}
\tag{2.9}
$$

在解析函数 $\zeta=\xi+i\eta=f(z)$ 所代表的变换 $(x,y)\mapsto(\xi,\eta)$ 下的变化. 根据偏微商的链式法则，有

$$
\frac{\partial}{\partial x}
=\frac{\partial\xi}{\partial x}\frac{\partial}{\partial\xi}
+\frac{\partial\eta}{\partial x}\frac{\partial}{\partial\eta},
\qquad
\frac{\partial}{\partial y}
=\frac{\partial\xi}{\partial y}\frac{\partial}{\partial\xi}
+\frac{\partial\eta}{\partial y}\frac{\partial}{\partial\eta},
$$

$$
\frac{\partial^2}{\partial x^2}
=\frac{\partial^2\xi}{\partial x^2}\frac{\partial}{\partial\xi}
+\frac{\partial^2\eta}{\partial x^2}\frac{\partial}{\partial\eta}
+\left(\frac{\partial\xi}{\partial x}\right)^2\frac{\partial^2}{\partial\xi^2}
+\left(\frac{\partial\eta}{\partial x}\right)^2\frac{\partial^2}{\partial\eta^2}
+2\frac{\partial\xi}{\partial x}\frac{\partial\eta}{\partial x}\frac{\partial^2}{\partial\xi\partial\eta},
$$

$$
\frac{\partial^2}{\partial y^2}
=\frac{\partial^2\xi}{\partial y^2}\frac{\partial}{\partial\xi}
+\frac{\partial^2\eta}{\partial y^2}\frac{\partial}{\partial\eta}
+\left(\frac{\partial\xi}{\partial y}\right)^2\frac{\partial^2}{\partial\xi^2}
+\left(\frac{\partial\eta}{\partial y}\right)^2\frac{\partial^2}{\partial\eta^2}
+2\frac{\partial\xi}{\partial y}\frac{\partial\eta}{\partial y}\frac{\partial^2}{\partial\xi\partial\eta}.
$$

所以

$$
\begin{aligned}
\nabla^2
={}&\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\\
={}&\left[\left(\frac{\partial\xi}{\partial x}\right)^2+
\left(\frac{\partial\xi}{\partial y}\right)^2\right]\frac{\partial^2}{\partial\xi^2}
+\left[\left(\frac{\partial\eta}{\partial x}\right)^2+
\left(\frac{\partial\eta}{\partial y}\right)^2\right]\frac{\partial^2}{\partial\eta^2}\\
&+\left(\frac{\partial^2\xi}{\partial x^2}+\frac{\partial^2\xi}{\partial y^2}\right)\frac{\partial}{\partial\xi}
+\left(\frac{\partial^2\eta}{\partial x^2}+\frac{\partial^2\eta}{\partial y^2}\right)\frac{\partial}{\partial\eta}\\
&+2\left(\frac{\partial\xi}{\partial x}\frac{\partial\eta}{\partial x}
+\frac{\partial\xi}{\partial y}\frac{\partial\eta}{\partial y}\right)
\frac{\partial^2}{\partial\xi\partial\eta}.
\end{aligned}
$$

再利用练习 2.3 中的结果及 (2.4) 和 (2.5) 式，就能得到

$$
\nabla^2=|f'(z)|^2\left(\frac{\partial^2}{\partial\xi^2}+\frac{\partial^2}{\partial\eta^2}\right).
\tag{2.10}
$$

这个结果表明：在解析函数 $\zeta=f(z)$ 所代表的保角变换之下，二维 Laplace 方程

$$
\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)u(x,y)=0
\tag{2.11}
$$

在 $f'(z)\ne0$ 的点，仍保持为二维 Laplace 方程：

$$
\left(\frac{\partial^2}{\partial\xi^2}+\frac{\partial^2}{\partial\eta^2}\right)
 u(x(\xi,\eta),y(\xi,\eta))=0;
\tag{2.12}
$$

二维 Poisson 方程

$$
\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)u(x,y)=\rho(x,y)
\tag{2.13}
$$

在变化后也仍然是二维 Poisson 方程：

$$
\left(\frac{\partial^2}{\partial\xi^2}+\frac{\partial^2}{\partial\eta^2}\right)
 u(x(\xi,\eta),y(\xi,\eta))
=\frac{1}{|f'(z)|^2}\rho(x(\xi,\eta),y(\xi,\eta)).
\tag{2.14}
$$

**练习 2.11**　证明：在解析函数 $\zeta=\xi+i\eta=f(z)$ 所代表的保角变换之下，面积元的变化公式是

$$
dxdy=|f'(z)|^{-2}d\xi d\eta.
$$

因此，如果把 Poisson 方程设想为平面静电场的方程，非齐次项是面电荷密度，则电荷守恒：

$$
\rho(x,y)dxdy=\rho(\xi,\eta)d\xi d\eta.
$$

在解析函数 $\zeta=f(z)$ 所代表的保角变换之下，二维 Helmholtz 方程

$$
\left(\frac{\partial^2}{\partial x^2}+\frac{\partial^2}{\partial y^2}\right)u(x,y)+k^2u(x,y)=0
\tag{2.15}
$$

在 $f'(z)\ne0$ 的点，也保持为二维 Helmholtz 方程：

$$
\left(\frac{\partial^2}{\partial\xi^2}+\frac{\partial^2}{\partial\eta^2}\right)
 u(x(\xi,\eta),y(\xi,\eta))
+\frac{k^2}{|f'(z)|^2}u(x(\xi,\eta),y(\xi,\eta))=0.
\tag{2.16}
$$

这样，应用保角变换，就可以把形状比较复杂的区域内的 Laplace 方程、Poisson 方程或 Helmholtz 方程的求解问题，转换为形状比较简单的区域（例如圆）内的求解问题.

## §2.6 多值函数

在复变函数中还存在另一类函数，即多值函数：设区域 $G\subseteq\mathbb C$，如果对于复数 $z\in G$，有多个复数 $w$ 与之对应，$w$ 和 $z$ 之间的这种对应关系记为 $f$，则称 $f$ 为定义在 $G$ 上的多值函数. 这部分内容在复变函数理论中占有重要地位，是复变函数中不可或缺的一部分.

从实际计算看，除了 §2.4 介绍过的初等函数（如幂函数、指数函数、三角函数等）外，我们还不可避免地要用到它们的逆运算，即开方、求对数等，它们都是多值函数. 本节只介绍开方和对数这两类基本的多值函数，通过它们阐述多值函数的一些基本概念. 别的多值函数都可以用这两类多值函数表示.

### 1. 根式函数 $\sqrt{z-a}$

以开平方为例. 首先给出平方根的定义：给定一个自变量值 $z$，凡是满足等式

$$
w^2=z-a
\tag{2.17}
$$

的 $w$ 值，就定义为 $z-a$ 的平方根 $\sqrt{z-a}$. 稍后我们会看到，给定一个 $z$ 值，通常会有两个 $w$ 值与之对应.

为了更清楚地看出多值函数的性质，下面仔细分析一下函数

$$
w=\sqrt{z-a}.
\tag{2.18}
$$

如果采用极坐标表达式

$$
w=\rho e^{i\phi},\qquad z-a=re^{i\theta},
$$

代入 (2.17) 式则有 $\rho^2e^{2i\phi}=re^{i\theta}$. 所以 $\rho^2=r$，$2\phi=\theta+2n\pi$，

$$
\rho=\sqrt r,\qquad \phi=\frac{\theta}{2}+n\pi,\qquad n=0,\pm1,\pm2,\cdots.
$$

这里 $\rho$ 和 $r$ 都是非负实数，所以上面的 $\sqrt r$ 表示算术平方根. 这样我们就看到，对于给定的一个 $z$ 值，有两个 $w$ 值与之对应：

$$
w_1(z)=\sqrt r\,e^{i\theta/2}\qquad（相当于\ n=0,\pm2,\cdots），
$$

$$
w_2(z)=\sqrt r\,e^{i(\pi+\theta/2)}=-\sqrt r\,e^{i\theta/2}
\qquad（相当于\ n=\pm1,\pm3,\cdots）.
$$

需要特别强调，$w$ 值的多值性来源于 $z-a$（而非自变量 $z$）辐角的多值性. $w$ 值的多值性表现为 $w$ 辐角的多值性. 我们称引起多值性的 $z-a$ 为宗量①. 为了确定起见，以后就把 $w=\sqrt{z-a}$ 明确表示成

$$
|w|=\sqrt{|z-a|},\qquad \arg w=\frac12\arg(z-a).
\tag{2.19}
$$

为了进一步揭示多值函数 $w=\sqrt{z-a}$ 的性质，不妨在 $z$ 复平面画一个不经过 $a$ 点的简单闭合曲线（即自身不相交的闭合曲线，例如圆周），研究当自变量 $z$ 在 $z$ 复平面从闭合曲线上某一点 $z_0$ 出发，沿此曲线逆时针连续变化一周回到 $z_0$ 点时，$w$ 复平面内 $w$ 值相应的连续变化情况. 我们发现，可能出现两种情形. 一种是 $w$ 值也还原. 例如当 $z$ 复平面内的 $z$ 沿图 2.4 中左上的圆 $C_1$ 一周回到原处时，由于 $a$ 点在圆外，$\arg(z-a)$ 不变，因此 $w$ 复平面内 $w$ 值不变. 另一种情形是 $w$ 值不还原. 例如当 $z$ 复平面内的 $z$ 从 $A$ 点出发沿图 2.4 中左下的 $C_2$ 经 $B,C,D$ 各点回到 $A$ 点时，由于 $a$ 点在闭合围道内，$\arg(z-a)$ 变化 $2\pi$，则由 (2.19) 式知 $\arg w$ 只变化 $\pi$，因此在 $w$ 复平面内由 $A$ 点经 $B,C,D$ 各点到达 $A'$ 点，$w$ 值并不还原.

上面的分析表明，$z=a$ 点在多值函数 $w=\sqrt{z-a}$ 中具有特殊的地位：它是否位于简单闭合路径内就决定了当 $z$ 沿这个路径一周回到原处时，相应的 $w$ 值是否能还原. 我们可以按照如下方法寻找这些特殊点：已知 $w$ 在某一点 $z_0$ 的空心邻域内每一点都有对应值，如果 $\exists r>0$，当 $z$ 绕圆周 $|z-z_0|=r$ 一圈回到原处时，$w$ 值不还原，而且当 $r\to0$ 时，$w$ 值始终不还原. 这种 $z_0$ 点就称为多值函数 $w(z)$ 的分支点②. 例如，$z=a$ 点就是 $w=\sqrt{z-a}$ 的分支点. 而 $z$ 复平面内其他的点，例如 $z=b\ (\ne a)$ 点，就不是 $w=\sqrt{z-a}$ 的分支点. 因为只要 $r<|a-b|\ne0$，$z$ 绕圆周 $|z-b|=r$（此时 $a$ 点在圆外）一周回到原处时，$w=\sqrt{z-a}$ 的值都能还原.

![Figure 2.4](...)

图 2.4　$z$ 沿闭合曲线一周回到原处时，$w=\sqrt{z-a}$ 值的不同变化

需要强调的是，这里所讨论的所有圆周上的点都应该位于多值函数自变量 $z$ 的取值范围以内，也就是说 $z$ 位于圆周上任意一点时，都能找到与之对应的 $w$ 值.

要考察 $z=\infty$ 是不是 $w=f(z)$ 的分支点，有两种方法可供选择. 一种方法是作足够大的简单闭合路径（使得在路径的外部区域内，有限远处均无分支点），考察 $z$ 沿此路径一周，$w$ 值是否不还原，即可判断 $\infty$ 点是否是分支点. 就函数 $\sqrt{z-a}$ 而言，只要这样的路径足够大（因此 $z=a$ 就一定在其内部区域内），$z$ 沿这样的路径一周回到原处，$w$ 值一定不还原，这就说明 $\infty$ 点也是 $\sqrt{z-a}$ 的分支点. 另一种方法是作变换 $z=1/t$，考察 $t=0$ 是不是 $f(1/t)$ 的分支点. 因为

$$
\sqrt{\frac1t-a}=\sqrt{\frac{1-at}{t}},
$$

我们看到，在 $t$ 复平面，$\forall r<1/|a|$，$t$ 绕圆 $|t|=r$ 逆时针一圈回到原处时，$\arg(1-at)$ 不变，但是 $\arg t$ 增加 $2\pi$，因此 $\arg\sqrt{(1/t)-a}$ 减少 $\pi$，$\sqrt{(1/t)-a}$ 的值始终不还原. 因此 $t=0$ 是 $\sqrt{(1/t)-a}$ 的分支点，也就是说 $z=\infty$ 是 $\sqrt{z-a}$ 的分支点.

一个多值函数通常有不止一个分支点. 例如 $w=\sqrt{z-a}$ 就有 $z=a$ 和 $z=\infty$ 两个分支点.

为了能够对 $w=\sqrt{z-a}$ 像单值函数一样研究其连续性、解析性等性质，我们需要将其单值化，也就是说要想办法让对每一个给定的 $z$，只有唯一的 $w$ 与之对应.

比较简单的办法是限制宗量 $z-a$ 的辐角变化范围. 当宗量 $z-a$ 的辐角限制在某个 $2\pi$ 周期内时，自变量 $z$ 取定义域内的任意一个值时，$w=\sqrt{z-a}$ 的辐角都能被唯一地确定，因而 $w$ 值也就唯一地确定. 例如，规定 $0\le\arg(z-a)<2\pi$ 或 $2\pi\le\arg(z-a)<4\pi$，等等.

**例 2.2**　设 $w=\sqrt{z-1}$，规定 $0\le\arg(z-1)<2\pi$，求 $w(2),w(i),w(0),w(-i)$.

**解**　$\arg w=\dfrac12\arg(z-1)$. 因为 $0\le\arg(z-1)<2\pi$，所以

$$
\left.\arg(z-1)\right|_{z=2}=0,\qquad w(2)=1,
$$

$$
\left.\arg(z-1)\right|_{z=i}=\frac34\pi,
\qquad w(i)=\sqrt[4]{2}\,e^{3\pi i/8},
$$

$$
\left.\arg(z-1)\right|_{z=0}=\pi,
\qquad w(0)=e^{\pi i/2}=i,
$$

$$
\left.\arg(z-1)\right|_{z=-i}=\frac54\pi,
\qquad w(-i)=\sqrt[4]{2}\,e^{5\pi i/8}.
$$

显然，在辐角规定 $0\le\arg(z-a)<2\pi$ 下，$w$ 的辐角被限制在 $0\le\arg w<\pi$，即限制在上半平面. 于是 $w=\sqrt{z-a}$ 在 $w$ 上半平面的值与自变量 $z$ 值之间就有一一对应的关系，这样就实现了将多值函数 $w=\sqrt{z-a}$ 单值化的要求. 如果规定 $2\pi\le\arg(z-a)<4\pi$，则 $\pi\le\arg w<2\pi$，$w$ 将限制在下半平面，$w$ 下半平面的值与 $z$ 值又有新的一一对应关系. 在 $4\pi\le\arg(z-a)<6\pi$，$6\pi\le\arg(z-a)<8\pi,\cdots$ 或 $-2\pi\le\arg(z-a)<0$，$-4\pi\le\arg(z-a)<-2\pi,\cdots$ 的规定下，还会不断重复出现这些结果. 宗量辐角变化的各个周期，给出多值函数的各个单值分支. 每个单值分支都是单值函数，整个多值函数就是它的各个单值分支的总和. 在上面的讨论中，多值函数 $w=\sqrt{z-a}$ 有两个单值分支，分别是 $w$ 的上半平面和下半平面：

$$
0\le\arg(z-a)<2\pi\ \text{给出单值分支 I：}\quad 0\le\arg w<\pi,
$$

$$
2\pi\le\arg(z-a)<4\pi\ \text{给出单值分支 II：}\quad \pi\le\arg w<2\pi.
$$

应当指出，宗量辐角变化范围的规定不是唯一的. 例如，也可以规定

$$
-\pi\le\arg(z-a)<\pi\qquad\text{和}\qquad
\pi\le\arg(z-a)<3\pi,
$$

或

$$
-\frac{3\pi}{2}\le\arg(z-a)<\frac{\pi}{2}
\qquad\text{和}\qquad
\frac{\pi}{2}\le\arg(z-a)<\frac{5\pi}{2}.
$$

将多值函数划分为若干个（甚至无穷个）单值分支，其实质就是限制 $z$ 在复平面内的变化方式，例如在上面的例子中，就是限制 $z$ 不得单独绕 $z=a$ 点或 $\infty$ 点转圈. 这种规定可以用几何方法形象化地表现出来（见图 2.5）. 在 $z$ 复平面内平行于实轴从 $z=a$ 点向右作一割线，一直延续到 $\infty$ 点. 如果规定在割线上岸 $\arg(z-a)=0$，就给出单值分支 I；规定在割线上岸 $\arg(z-a)=2\pi$，就给出单值分支 II. 两个单值分支合起来，就得到完整的 $w$ 复平面，也就是整个多值函数 $w$. 在作了割线之后，自变量 $z$ 在 $z$ 复平面内变化的路径就不允许穿越割线. 由于割线连接了多值函数的两个分支点，$z=a$ 和 $z=\infty$，因此，$z$ 不再能够绕单独一个分支点一周（但不禁止同时绕两个或者更多个分支点一周）.

![Figure 2.5](...)

图 2.5　$w=\sqrt{z-a}$ 的两个单值分支

割线的作法也是多种多样，甚至不必是直线（例如本章习题第 12 题）. 只要用割线正确连接了多值函数的分支点，同时适当规定割线一侧相关宗量的辐角值（例如本章习题第 10 题），或者等效地，规定分支点以外某一点处 $w$ 的值（例如本章习题第 11 题），多值函数在复平面内任意一点的值就唯一确定了.

对于更复杂的根式函数，例如 $w=\sqrt[3]{z-a}$ 或 $\sqrt{(z-a)(z-b)(z-c)}$ 等，也可以类似地讨论，只是需要找出全部分支点，并且正确地确定割线的作法. 要注意的是，如果函数有三个以上分支点，不一定需要用一条割线将所有分支点都连接起来. 例如当 $a,b$ 和 $c$ 三个复数两两不等时，函数 $\sqrt{(z-a)(z-b)(z-c)}$ 有四个分支点：$a,b,c$ 和 $\infty$. 可以作两条割线，一条连接 $z=a$ 和 $z=b$ 两个分支点，另一条连接其余两个分支点 $z=c$ 和 $z=\infty$ 即可. 这两条线可以不相交.

一旦涉及多值函数，涉及多值函数的运算，就需要特别留意有关宗量辐角的规定（或者不妨就称之为相关多值函数的定义域）. 有些运算公式，在多值函数的情况下就可能只有条件地成立，或者反过来说，就不会无条件成立. 例如，$\sqrt{(z-a)(z-b)}$ 就不能无条件地写成为 $\sqrt{z-a}$ 和 $\sqrt{z-b}$ 的乘积. 严格说来，这里涉及三个多值函数. 如果它们是互相独立的，可以独立地规定单值分支，那么 $\sqrt{(z-a)(z-b)}$ 和 $\sqrt{z-a}\cdot\sqrt{z-b}$ 可以并不相等.

将多值函数划分为单值分支之后，每个单值分支都是单值函数，因而可以讨论它们的解析性. 在分支点以及割线上的点都是奇点，这是因为分支点以及连接分支点的割线为多个单值分支所共有，而不只属于某一个单值分支，也不存在只属于一个单值分支的邻域. 由于割线作法具有一定的任意性，因此存在这种情况：复平面内分支点以外的点，在一种割线作法下，由于被割线穿过，因而不解析，但是如果换一种作法，使得割线不通过这一点，那么规定单值分支以后，这一点就可以是解析点.

将多值函数划分单值分支的方法，简单易行，但是有一个明显的缺点，就是限制了宗量的辐角变化范围，因而不能讨论比较复杂的问题. 为了克服这个缺点，另一种确定自变量 $z$ 与函数 $w$ 值对应关系的办法是：规定 $w$ 在分支点以外的某一点 $z_0$ 的值，并明确说明由 $z_0$ 出发到达 $z$ 点的连续变化路线. 当 $z$ 沿这曲线连续变化时，$w$ 值也随之连续变化. 请看下面的例题.

**例 2.3**　已知 $f(z)=\sqrt{z+1}$，规定 $f(0)=-1$，求 $f(3)$.

(1) 若 $z$ 沿图 2.6 中的路径 $C_1$ 从原点到达 $z=3$ 点；

(2) 若 $z$ 沿图 2.6 中的路径 $C_2$ 从原点到达 $z=3$ 点.

![Figure 2.6](...)

图 2.6

**解**　按照题意，

$$
f(z)=\sqrt{z+1}=\sqrt{|z+1|}\,e^{i\arg(z+1)/2}.
$$

由已知条件 $f(0)=-1$，知

$$
\left.\arg(z+1)\right|_{z=0}=2(2n+1)\pi,\qquad n\in\mathbb Z.
$$

(1) 若 $z$ 沿图 2.6 中的路径 $C_1$ 从原点到达 $z=3$ 点，则 $\Delta\arg(z+1)=0$. 因此，

$$
\left.\arg(z+1)\right|_{z=3}=2(2n+1)\pi,
\qquad
f(3)=\left.\sqrt{|z+1|}\,e^{i\arg(z+1)/2}\right|_{z=3}=-2.
$$

(2) 若 $z$ 沿图 2.6 中的路径 $C_2$ 从原点到达 $z=3$ 点，则 $\Delta\arg(z+1)=2\pi$. 因此，

$$
\left.\arg(z+1)\right|_{z=3}=4(n+1)\pi,
\qquad
f(3)=\left.\sqrt{|z+1|}\,e^{i\arg(z+1)/2}\right|_{z=3}=2.
$$

采用这种办法，$z$ 的变化路线不受限制，因而就可以从一个单值分支移动到另一个单值分支. 在几何图形上，这相当于将两个割开的 $z$ 复平面粘接起来. 以 $\sqrt{z-a}$ 为例，将图 2.6 中第一个 $z$ 复平面的割线下岸（$\arg(z-a)=2\pi$）和第二个 $z$ 复平面的割线上岸（$\arg(z-a)=2\pi$）合并，同时，将第一个 $z$ 复平面的割线上岸（$\arg(z-a)=0$）和第二个 $z$ 复平面的割线下岸（$\arg(z-a)=4\pi$）合并. 这就构成了二叶 Riemann 面（见图 2.7）. 显然，这个二叶 Riemann 面无法在三维空间中实现，我们只能到更高维空间中去想象. 对于 $w=\sqrt{z-a}$ 来说，二叶 Riemann 面上的点和 $w$ 复平面内的点是一一对应的. Riemann 面是 Riemann 首先研究的，因而后来以他的名字命名的. Riemann 面概念的提出，使得数学家们困扰了很多年的多值函数问题得以圆满解决.

![Figure 2.7](...)

图 2.7　$w=\sqrt{z-a}$ 的 Riemann 面

① 宗量通常不同于自变量. 例如，多值函数 $\sqrt{z-a}$ 的宗量就是 $z-a$，多值函数 $\sqrt[3]{(z-a)(z-b)}$ 的宗量就是 $(z-a)(z-b)$. 当然，也有宗量就是自变量的情形. 例如多值函数 $\sqrt z$ 的宗量就是自变量 $z$.

② 即 branch point，原译作“枝点”，按照《数学名词》的规范，应译为“支点”或“分支点”，前者易与力学中的支点相混淆，故本书中一律采用“分支点”.

### 2. 对数函数 $\ln z$

给定自变量值 $z$，凡满足 $e^w=z$ 的所有 $w$ 值均称为 $w=\ln z$ 的值.

令 $w=u+iv$，$z=re^{i\theta}$，即得 $e^u\cdot e^{iv}=re^{i\theta}$. 所以

$$
u=\ln r=\ln|z|,\qquad v=\theta+2n\pi\qquad(n=0,\pm1,\pm2,\cdots).
$$

这里 $r>0$，因此 $\ln r$ 就是正数的实值对数，有唯一值. 以后就把 $w=\ln z$ 明确表示为

$$
w=\ln z=\ln|z|+i(\theta+2n\pi)=\ln|z|+i\arg z.
\tag{2.20}
$$

对数函数 $w=\ln z$ 也是多值函数，$w$ 值的多值性来源于宗量① $z$ 辐角的多值性，$w$ 值的多值性表现在 $w$ 值虚部的多值性. 对应每一个 $z$ 值，有无穷多个 $w$ 值，它们的实部相同，虚部相差 $2\pi$ 的整数倍. 图 2.8 给出了 $w=\ln z$ 的示意图.

① $\ln z$ 的宗量就是自变量 $z$.

$w=\ln z$ 的分支点是 $z=0$ 和 $\infty$. 作割线连接 0 与 $\infty$，并规定割线一侧的 $\arg z$ 值，即可得到 $w=\ln z$ 的单值分支.

![Figure 2.8](...)

图 2.8　$w=\ln z$：实部和虚部的等值线

![Figure 2.9](...)

图 2.9　$w=\ln z$ 的 Riemann 面

$w=\ln z$ 有无穷多个单值分支. 每个单值分支内，$w=\ln z$ 都是解析函数，而且都有

$$
\frac{d}{dz}(\ln z)=\frac1z.
\tag{2.21}
$$

相应地，$w=\ln z$ 的 Riemann 面是无穷多叶的，见图 2.9. $w=\ln z$ 的 Riemann 面上的点与 $w$ 复平面内的点也是一一对应的.

**练习 2.12**　证明：如果 $f(z)$ 是 $G$ 内的解析函数，且 $f(z)$ 的模或辐角为常数，则 $f(z)$ 必为常数.

### 3. 其他多值函数

例如，还有复反三角函数和指数为任意复数的幂函数

$$
\arcsin z=\frac1i\ln\left(iz+\sqrt{1-z^2}\right),
\tag{2.22}
$$

$$
\arccos z=\frac1i\ln\left(z+\sqrt{z^2-1}\right),
\tag{2.23}
$$

$$
\arctan z=\frac1{2i}\ln\frac{1+iz}{1-iz},
\tag{2.24}
$$

$$
z^\alpha=e^{\alpha\ln z}\qquad(\alpha\ \text{为任意复数}).
\tag{2.25}
$$

它们也都是多值函数，并且是对数函数或对数函数与根式函数的组合，因此它们的多值性可以根据这两种基本的多值函数来讨论. 对于 $\arctan z$ 和 $z^\alpha$ 来说，只涉及对数函数，因而它们的多值性就只需重复上面的讨论. 对于 $\arcsin z$ 和 $\arccos z$，同时涉及两种多值函数，这里不再一一讨论.

最后，还需要提到，对一般的多值函数，通常难以直接得到全局的 Riemann 面，需要用到我们后面将介绍的解析延拓. 而前面讨论的根式函数和对数函数这两种最基本的多值函数，只是两个简单的例子，能够画出全局的 Riemann 面.

## 习题

1. 判断下列函数在何处可导（并求出其导数），在何处解析：

   (1) $|z|$；

   (2) $z^*$；

   (3) $z\operatorname{Re}z$；

   (4) $(x^2+2y)+i(x^2+y^2)$；

   (5) $3x^2+2iy^3$；

   (6) $(x-y)^2+2i(x+y)$.

2. 证明平面极坐标系 $(r,\theta)$ 下的 Cauchy-Riemann 方程：

$$
\frac{\partial u}{\partial r}=\frac1r\frac{\partial v}{\partial\theta},
\qquad
\frac{\partial v}{\partial r}=-\frac1r\frac{\partial u}{\partial\theta},
$$

并由此证明：

$$
f'(z)=\frac rz\left(\frac{\partial u}{\partial r}+i\frac{\partial v}{\partial r}\right)
=\frac1z\left(\frac{\partial v}{\partial\theta}-i\frac{\partial u}{\partial\theta}\right).
$$

其中 $u(r,\theta)$ 和 $v(r,\theta)$ 分别为复变函数的实部和虚部.

3. 设 $z=x+iy$，已知解析函数 $f(z)=u(x,y)+iv(x,y)$ 的实部 $u(x,y)$ 如下，试求出解析函数 $f(z)$：

   (1) $x^2-y^2+x$；

   (2) $\dfrac{x}{x^2+y^2}$；

   (3) $e^y\cos x$；

   (4) $\cos x\cosh y$.

4. 设 $z=x+iy$，已知解析函数 $f(z)=u(x,y)+iv(x,y)$ 的实部或虚部如下，试求 $f'(z)$：

   (1) $u=x+y$；

   (2) $u=\sin x\cosh y$.

5. 若 $f(z)=u(x,y)+iv(x,y)$ 解析，且 $u-v=(x-y)(x^2+4xy+y^2)$，试求 $f(z)$.

6. 解下列方程：

   (1) $\sin z=\dfrac34+\dfrac{i}{4}$；

   (2) $\cos z=4$；

   (3) $\tan z=i$；

   (4) $2\cosh^2z-3\cosh z+1=0$.

7. (1) 求解方程 $(1-z)^n-1=0$；

   (2) 证明：

$$
\sin\frac\pi n\sin\frac{2\pi}n\cdots\sin\frac{(n-1)\pi}n
=\frac{n}{2^{n-1}}.
$$

8. 判断下列哪些是函数，哪些是多值函数：

   (1) $\sqrt{z^2-1}$；

   (2) $z+\sqrt{z-1}$；

   (3) $\sin\sqrt z$；

   (4) $\cos\sqrt z$；

   (5) $\dfrac{\sin\sqrt z}{\sqrt z}$；

   (6) $\dfrac{\cos\sqrt z}{\sqrt z}$；

   (7) $\ln\sin z$；

   (8) $\sin(i\ln z)$.

9. 找出下列多值函数的分支点，并讨论 $z$ 绕一个分支点移动一周回到原处后多值函数值的变化. 如果同时绕两个、三个，乃至更多个分支点一周，多值函数的值又如何变化？

   (1) $\sqrt{(z-a)(z-b)},\quad a\ne b$；

   (2) $\sqrt[3]{(z-a)(z-b)},\quad a\ne b$；

   (3) $\sqrt{1-z^3}$；

   (4) $\sqrt[3]{1-z^3}$；

   (5) $\ln(z^2+1)$；

   (6) $\ln\cos z$.

10. 规定在图 2.10 所示割线上岸 $\arg(z-2)=0$，试求多值函数 $w=z\sqrt[3]{z-2}$ 在割线下岸 $z=3$ 处的数值.

又问：这个多值函数有几个单值分支？求出它在其他分支中割线下岸 $z=3$ 处的值.

![Figure 2.10](...)

图 2.10

11. 规定函数 $w=\sqrt{z^2-2z+2}$ 在 $z=0$ 时 $w=\sqrt2$，求当 $z$ 由原点出发沿圆 $|z-(1+i)|=\sqrt2$ 逆时针通过 $i$ 时的函数值.

当 $z$ 继续沿该圆移动而回到原点时，函数之值又如何？

12. 已知 $w=\ln(1-z^2)$，规定 $w(0)=0$，试讨论当 $z$ 限制在图 2.11 (a) 和 (b) 中的 $w(3)$ 值. 若作割线如图 2.11 (c)，则在割线上、下岸 $z=3$ 处 $w$ 又取何值？

![Figure 2.11](...)

图 2.11

13. 求下列多值函数在指定点的全部可能取值：

   (1) $\ln z,\ z=1,i,-1,1+i$；

   (2) $z^i,\ z=2,i,-1,(1+i)$.

14. 已知多值函数 $f(z)=z^{-p}(1-z)^p$，$p$ 为实数. 若在实轴上沿 0 到 1 作割线，规定在割线上岸 $\arg z=\arg(1-z)=0$，试求 $f(\pm i)$ 和 $f(\infty)$.

15. 反正切函数 $\arctan z$ 的定义为

$$
\arctan z\equiv\frac1{2i}\ln\frac{1+iz}{1-iz}.
$$

若作割线如图 2.12，并规定

$$
\left.\arctan z\right|_{z=0}=\pi,
$$

求函数在 $z=2$ 处的导数值.

![Figure 2.12](...)

图 2.12

16. 若函数 $f(z)$ 在区域 $G$ 内解析，且其模为常数，证明 $f(z)$ 本身也必为常数.
