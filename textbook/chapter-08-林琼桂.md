# 第八章 Fourier 变换法*

## 目录

§1 Fourier 变换 ................................................................ 168  
§1.1 复数形式的 Fourier 级数 .................................................... 168  
§1.2 Fourier 变换 ................................................................ 169  
§1.3 二维与三维 Fourier 变换 .................................................... 170  
§1.4 Fourier 变换的性质 .......................................................... 171  

§2 无界弦的自由振动 .............................................................. 172  
§2.1 Fourier 变换法 .............................................................. 172  
§2.2 行波法简介 .................................................................. 173  

§3 $\delta$ 函数 .................................................................. 174  
§3.1 定义一 ...................................................................... 174  
§3.2 定义二 ...................................................................... 175  
§3.3 $\delta$ 函数的基本性质 ...................................................... 176  
§3.4 $\delta$ 函数的几种表达式 .................................................... 176  
§3.5 *$\delta$ 函数的导数 ......................................................... 178  
§3.6 *复合 $\delta$ 函数 .......................................................... 178  

§4 无界杆的热传导问题 ............................................................ 179  
§4.1 无源问题 .................................................................... 180  
§4.2 *有源问题 ................................................................... 181  

§5 *无界弦的受迫振动 ............................................................. 181  

补充习题 .......................................................................... 183

上章研究了求解数理方程定解问题的分离变量法。由于所考虑的问题都只涉及两个自变量，所以情况比较简单。如果问题涉及多个自变量，情况就会比较复杂。在曲线坐标系中对数理方程分离变量，还常常导出一些特殊的常微分方程，它们的解就是所谓的特殊函数。所以分离变量法是一个很大的课题，今后还有较多的章节研究这一问题。

分离变量法适用于有界区间或有界区域，对于无界区间或无界区域上的数学物理方程定解问题，常常使用 Fourier 变换法。Fourier 变换法是积分变换法的一种。另一种常用的积分变换法是 Laplace 变换法。积分变换法的基本精神也是分离变量，只是形式上不甚明显。本章介绍 Fourier 变换法。对于其它的积分变换法，本书不作介绍。

## §1 Fourier 变换

### §1.1 复数形式的 Fourier 级数

函数族 $\{e^{in\pi x/l}\}_{n=-\infty}^{+\infty}$ 在区间 $[-l,l]$ 上是正交的：当 $m\ne n$，

$$
(e^{im\pi x/l},e^{in\pi x/l})\equiv\int_{-l}^{l}(e^{im\pi x/l})^*e^{in\pi x/l}\,dx
=\int_{-l}^{l}e^{i(n-m)\pi x/l}\,dx=0. \tag{1}
$$

注意内积定义中出现的复共轭，这是第二章关于实函数的内积定义的推广。

**注**　内积是普通三维空间中的矢量点积的推广。由于 $\boldsymbol a\cdot\boldsymbol a=|\boldsymbol a|^2$ 是矢量 $\boldsymbol a$ 的长度平方，作为它的推广，函数与自身的内积应该是实数。所以，上面的定义中出现复共轭是很容易理解的。

集合 $U$ 和 $V$ 的积集 $U\times V$ 定义为 $U\times V=\{(\alpha,\beta)|\alpha\in U,\beta\in V\}$，即其元素是由 $U$ 的一个元素和 $V$ 的一个元素构成的有序对。

在集合 $V$ 上定义了加法和数乘（设数域为 $\mathbb C$），$V$ 就成为线性空间。定义一种运算 $(\cdot,\cdot):V\times V\to\mathbb C$（这表示运算的定义域是 $V\times V$，值域是 $\mathbb C$），即 $(\alpha,\beta)\mapsto c\in\mathbb C$（这表示运算将有序对 $(\alpha,\beta)$ 映射为复数 $c$，其中 $\alpha,\beta\in V$），如果该运算满足下面三条规则，则称为内积。

1. 对称性：$(\beta,\alpha)=(\alpha,\beta)^*$；
2. 线性性：$(\gamma,k\alpha+l\beta)=k(\gamma,\alpha)+l(\gamma,\beta)$，其中 $k,l\in\mathbb C$；
3. 正定性：$(\alpha,\alpha)\ge0$，等号成立当且仅当 $\alpha=0$。

由第一条易知 $(\alpha,\alpha)$ 为实数；第二条和第一条结合可得 $(k\alpha+l\beta,\gamma)=k^*(\alpha,\gamma)+l^*(\beta,\gamma)$；第三条保证 $(\alpha,\alpha)$ 是矢量长度平方的适当推广。注意内积必须定义在线性空间中，因为如果集合 $V$ 中没有线性结构，则第二条中的 $k\alpha+l\beta$ 将没有意义，第三条中的 $\alpha=0$ 也没有明确意义。

**例 1**　考虑复数域 $\mathbb C$ 上的 $n$ 维矢量空间 $\mathbb C^n$（$n$ 个 $\mathbb C$ 的积集），记 $\alpha=(a_1,\cdots,a_n)\in\mathbb C^n$，$\beta=(b_1,\cdots,b_n)\in\mathbb C^n$，定义

$$
(\alpha,\beta)=\sum_{i=1}^{n}a_i^*b_i.
$$

容易验证，这一运算确实是内积（即满足上述三条规则）。如果将 $\mathbb C^n$ 换为 $\mathbb R^n$，容易看出，这就是普通三维矢量空间 $\mathbb R^3$ 中的点积的推广。

**例 2**　考虑三维空间的区域 $D$ 上所有平方可积的复值函数构成的集合 $\{\psi(\boldsymbol r)\}$，按常规的方式定义加法和数乘（仍设数域为 $\mathbb C$），即 $(\psi+\varphi)(\boldsymbol r)=\psi(\boldsymbol r)+\varphi(\boldsymbol r)$，$(k\psi)(\boldsymbol r)=k\psi(\boldsymbol r)$，则构成线性空间。容易验证，按下式定义的运算是内积

$$
(\varphi,\psi)=\int_D\varphi^*(\boldsymbol r)\psi(\boldsymbol r)\,d\boldsymbol r.
$$

上面式 (1) 中的内积是一维空间的一个特例。

可以证明，上述函数族在区间 $[-l,l]$ 上还是完备的。所以，区间 $[-l,l]$ 上任何解析性质良好的函数 $f(x)$ 都可以展开为

$$
f(x)=\sum_{n=-\infty}^{+\infty}f_ne^{in\pi x/l}, \tag{2}
$$

其中系数 $f_n$ 易由正交性求得为

$$
f_n=\frac{1}{2l}\int_{-l}^{l}f(\xi)e^{-in\pi\xi/l}\,d\xi,\qquad n\in\mathbb Z. \tag{3}
$$

如果 $f(x)$ 是以 $2l$ 为周期的周期函数，则式 (2) 对一切 $x$ 成立，它就是微积分中熟知的复形式的 Fourier 级数。

**注**　函数族 $\{e^{in\phi}\}_{n=-\infty}^{+\infty}$ 等价于 $\{\cos n\phi,\sin n\phi\}_{n=0}^{+\infty}$，它在区间 $[-\pi,\pi]$ 上是完备的，这是我们已经熟知的结论。由此容易知道函数族 $\{e^{in\pi x/l}\}_{n=-\infty}^{+\infty}$ 在区间 $[-l,l]$ 上是完备的。

### §1.2 Fourier 变换

现考虑展开式 (2) 当 $l\to\infty$ 时的极限。

将式 (3) 代入式 (2)，可得

$$
f(x)=\sum_{n=-\infty}^{+\infty}\left[\frac{1}{2l}\int_{-l}^{l}f(\xi)e^{-in\pi\xi/l}\,d\xi\right]e^{in\pi x/l}.
$$

令 $k_n=n\pi/l$，$\Delta k_n=k_{n+1}-k_n=\pi/l$，上式可以改写为

$$
f(x)=\sum_{k_n=-\infty}^{+\infty}\left[\frac{1}{2\pi}\int_{-l}^{l}f(\xi)e^{-ik_n\xi}\,d\xi\right]e^{ik_nx}\Delta k_n.
$$

当 $l\to\infty$ 时，$\Delta k_n\to0$，则 $k_n$ 的取值由分立变为连续，上式的求和也就变为积分，故得

$$
f(x)=\int_{-\infty}^{+\infty}\left[\frac{1}{2\pi}\int_{-\infty}^{+\infty}f(\xi)e^{-ik\xi}\,d\xi\right]e^{ikx}\,dk. \tag{4}
$$

由式 (4)，如果定义

$$
F(k)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f(x)e^{-ikx}\,dx, \tag{5a}
$$

则

$$
f(x)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}F(k)e^{ikx}\,dk. \tag{5b}
$$

$F(k)$ 称为 $f(x)$ 的 Fourier 变换，也称为 Fourier 变换的像函数；$f(x)$ 称为 $F(k)$ 的 Fourier 反变换，也称为 Fourier 变换的原函数。它们的关系常用以下一些符号来表示

$$
f(x)\leftrightarrow F(k); \tag{6a}
$$

或

$$
F(k)=\mathcal F[f(x)],\qquad f(x)=\mathcal F^{-1}[F(k)]; \tag{6b}
$$

或

$$
f(x)\coloneqq F(k),\qquad F(k)\eqqcolon f(x). \tag{6c}
$$

我们主要使用第一种写法，它最简单。有些书上将式 (5b) 称为 $f(x)$ 的 Fourier 积分。

以上通过对 Fourier 级数取 $l\to\infty$ 的极限形式上得到 Fourier 变换和反变换关系，可以证明，只要 $f(x)$ 具有良好的解析性质（参看下面的定理），则以上取极限的过程是合理的。

Fourier 变换和反变换关系可以这样理解：给定函数 $f(x)$，按式 (5a) 定义像函数 $F(k)$（如果 $f(x)$ 满足下面定理的条件，则像函数存在），然后将像函数 $F(k)$ 代入式 (5b) 的右边，积分后就可以重新得到原函数 $f(x)$。这里的关键问题是最后一句话是否真的成立。下面的定理回答了这个问题。

**定理**　设 $f(x)$ 定义在 $(-\infty,\infty)$ 上，且满足 (1) 在 $x$ 的任何有限区间上，$f(x)$ 只有有限个第一类间断点和有限个极值；(2) 积分 $\int_{-\infty}^{+\infty}|f(x)|\,dx$ 存在，即 $f(x)$ 在 $(-\infty,\infty)$ 上绝对可积，则像函数 $F(k)$ 存在且式 (5b) 右边的积分收敛于 $f(x)$（在 $f(x)$ 的连续点）或 $[f(x+0)+f(x-0)]/2$（在 $f(x)$ 的间断点）。

本定理证明从略。

### §1.3 二维与三维 Fourier 变换

对二元函数 $f(x,y)$，先把 $y$ 看作参数，对 $x$ 作 Fourier 变换，像函数记作 $\widetilde F(k,y)$，即变换关系为 $f(x,y)\leftrightarrow\widetilde F(k,y)$；然后将 $\widetilde F(k,y)$ 中的 $k$ 当作参数，对 $y$ 作 Fourier 变换，像函数记作 $F(k,l)$，即变换关系为 $\widetilde F(k,y)\leftrightarrow F(k,l)$。消去 $\widetilde F(k,y)$，即得下列 Fourier 变换与反变换关系

$$
F(k,l)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}f(x,y)e^{-i(kx+ly)}\,dx\,dy, \tag{7a}
$$

$$
f(x,y)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}\int_{-\infty}^{+\infty}F(k,l)e^{i(kx+ly)}\,dk\,dl. \tag{7b}
$$

类似地，对三元函数 $f(x,y,z)$，简记 $\boldsymbol r=(x,y,z)$，$\boldsymbol k=(k_1,k_2,k_3)$，有下列 Fourier 变换与反变换关系

$$
F(\boldsymbol k)=\frac{1}{(2\pi)^{3/2}}\int_{-\infty}^{+\infty} f(\boldsymbol r)e^{-i\boldsymbol k\cdot\boldsymbol r}\,d\boldsymbol r, \tag{8a}
$$

$$
f(\boldsymbol r)=\frac{1}{(2\pi)^{3/2}}\int_{-\infty}^{+\infty} F(\boldsymbol k)e^{i\boldsymbol k\cdot\boldsymbol r}\,d\boldsymbol k, \tag{8b}
$$

其中 $\boldsymbol k\cdot\boldsymbol r=k_1x+k_2y+k_3z$。对一般的 $n$ 元函数，也有类似的结果。

### §1.4 Fourier 变换的性质

关于 Fourier 变换的性质，读者在初学时只需关注最重要的几条，以免一时难以消化。在作 Fourier 变换和反变换的计算时，可以引用有关性质，也可以直接计算（这相当于把有关性质的证明包含在计算过程中）。如果经常使用 Fourier 变换，那么引用有关性质还是比直接计算来得方便，所以下面也列出较为次要的性质以供参考。

设 $f(x)\leftrightarrow F(k)$，$f_i(x)\leftrightarrow F_i(k)$（$i=1,2$），则有以下性质。

(1) 线性定理：$a_1f_1(x)+a_2f_2(x)\leftrightarrow a_1F_1(k)+a_2F_2(k)$。

这一性质是显而易见的，它直接来源于积分的线性性质。但是它的重要性也同样是显而易见的，不应受到忽视。

(2) 微分定理：如果 $f(\pm\infty)=f'(\pm\infty)=\cdots=f^{(n-1)}(\pm\infty)=0$，则 $f^{(n)}(x)\leftrightarrow(ik)^nF(k)$。

**证明**　由于 $f(\pm\infty)=0$，由定义

$$
f'(x)\leftrightarrow\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f'(x)e^{-ikx}\,dx
=\frac{1}{\sqrt{2\pi}}f(x)e^{-ikx}\Big|_{-\infty}^{\infty}
+\frac{ik}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f(x)e^{-ikx}\,dx
=ikF(k),
$$

其中利用了分部积分公式。同理可证 $n>1$ 的情况。证毕。

(3) 卷积定理：$f_1(x)*f_2(x)\leftrightarrow F_1(k)F_2(k)$，其中

$$
f_1(x)*f_2(x)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f_1(\xi)f_2(x-\xi)\,d\xi
=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f_1(x-\xi)f_2(\xi)\,d\xi \tag{9}
$$

称为两个函数 $f_1(x)$ 和 $f_2(x)$ 的卷积（convolution integral），其中两种表达式的等价性容易直接验证。

**注**　① 这一定理主要用于由像函数求原函数，基本前提是 $F_1(k)$ 和 $F_2(k)$ 的原函数已经知道或容易求出。② 定理可由定义直接证明，但证明略为繁琐，所以读者应该记得有这么一个定理（卷积的表达式可以不必记忆），以便需要时直接引用。③ 请读者尝试做出定理的证明。

(4) 延迟定理：$f(x-\xi)\leftrightarrow e^{-ik\xi}F(k)$，其中 $\xi\in\mathbb R$。

**证明**　由定义并作积分变量置换 $y=x-\xi$，即得

$$
f(x-\xi)\leftrightarrow\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f(x-\xi)e^{-ikx}\,dx
=e^{-ik\xi}\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}f(y)e^{-iky}\,dy
=e^{-ik\xi}F(k).
$$

证毕。

(5) 积分定理：如果 $\int_{-\infty}^{+\infty}f(x)\,dx=0$，则 $g(x)\equiv\int_{-\infty}^{x}f(\xi)\,d\xi\leftrightarrow F(k)/ik$。

**证明**　设 $g(x)\leftrightarrow G(k)$，由定义和已知条件显然有 $g(\pm\infty)=0$，故由微分定理有 $\mathcal F[g'(x)]=ikG(k)$。另一方面，$\mathcal F[g'(x)]=\mathcal F[f(x)]=F(k)$。两式结合即得 $ikG(k)=F(k)$，或 $G(k)=F(k)/ik$。证毕。

(6) 相似定理：$f(ax)\leftrightarrow (1/|a|)F(k/a)$，其中 $a\ne0$。

本定理证明简单，从略。

Fourier 变换还有其它一些性质，不一一介绍。

## §2 无界弦的自由振动

无界弦即无穷长的弦，实际上并不存在。但是，如果弦较长，而我们所关心的部分离边界较远，在较短的时间内，边界对所关心的部分尚未发生作用，则可以忽略边界的存在，从而将有界弦抽象为无界弦。

### §2.1 Fourier 变换法

考虑无界弦在初始激励下的振动。定解问题为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0,
\qquad (-\infty<x<+\infty,\quad t>0), \tag{10a}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad (-\infty<x<+\infty), \tag{10b}
$$

$$
\left.\frac{\partial u}{\partial t}\right|_{t=0}=\psi(x),
\qquad (-\infty<x<+\infty). \tag{10c}
$$

**注**　这里忽略了默认的边界条件 $u|_{x=\pm\infty}=0$ 和 $\partial u/\partial x|_{x=\pm\infty}=0$。

由于 $t$ 的变化范围是 $(0,+\infty)$，而 $x$ 的变化范围是 $(-\infty,+\infty)$，故若考虑用 Fourier 变换法求解，则应该对 $x$ 作 Fourier 变换，设

$$
u(x,t)\leftrightarrow U(k,t),\qquad \varphi(x)\leftrightarrow\Phi(k),\qquad \psi(x)\leftrightarrow\Psi(k). \tag{11}
$$

对定解问题 (10) 的各式作 Fourier 变换，即对每个方程两边同乘以 $(1/\sqrt{2\pi})e^{-ikx}$ 并对 $x$ 积分，利用微分定理，可得常微分方程初值问题

$$
\frac{d^2U}{dt^2}+k^2a^2U=0, \tag{12a}
$$

$$
U\big|_{t=0}=\Phi(k),\qquad \left.\frac{dU}{dt}\right|_{t=0}=\Psi(k). \tag{12b}
$$

容易求得上式的解为

$$
U(k,t)=\Phi(k)\cos kat+\Psi(k)\frac{\sin kat}{ka}. \tag{13}
$$

剩下的问题是对 $U(k,t)$ 作反变换求出原函数 $u(x,t)$。由反变换公式，

$$
u(x,t)=\mathcal F^{-1}[U(k,t)]
=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}U(k,t)e^{ikx}\,dk. \tag{14}
$$

将式 (13) 代入，两项分别计算，可得

$$
\mathcal F^{-1}[\Phi(k)\cos kat]
=\frac12\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}\Phi(k)
[e^{ik(x-at)}+e^{ik(x+at)}]dk
=\frac12[\varphi(x-at)+\varphi(x+at)], \tag{15a}
$$

以及

$$
\begin{aligned}
\mathcal F^{-1}\left[\Psi(k)\frac{\sin kat}{ka}\right]
&=\mathcal F^{-1}\left[\int_0^t\Psi(k)\cos ka\tau\,d\tau\right]
=\int_0^t\mathcal F^{-1}[\Psi(k)\cos ka\tau]d\tau\\
&=\frac12\int_0^t[\psi(x-a\tau)+\psi(x+a\tau)]d\tau\\
&=\frac{1}{2a}\int_{x-at}^{x}\psi(\xi)d\xi
+\frac{1}{2a}\int_x^{x+at}\psi(\xi)d\xi\\
&=\frac{1}{2a}\int_{x-at}^{x+at}\psi(\xi)d\xi. \tag{15b}
\end{aligned}
$$

其中第二步利用了 Fourier 反变换的线性性（如果具体写出反变换的形式，也可以说是交换了对 $\tau$ 积分和对 $k$ 积分的顺序），第三步利用了式 (15a) 的结果，第四步作了积分变量置换（第一项令 $\xi=x-a\tau$，第二项令 $\xi=x+a\tau$），最后一步利用了定积分的基本性质。综合两式的结果，得到最后结果为

$$
u(x,t)=\frac12[\varphi(x-at)+\varphi(x+at)]
+\frac{1}{2a}\int_{x-at}^{x+at}\psi(\xi)d\xi. \tag{16}
$$

这称为 d’Alembert 公式。

如果 $\psi(x)=0$，则

$$
u(x,t)=\frac12[\varphi(x-at)+\varphi(x+at)]. \tag{17}
$$

容易看出，初始位移的波形一半向右传播，一半向左传播，物理图像非常简单。$\psi(x)\ne0$ 时情况略为复杂，但整个波形仍然是由一列右行波和一列左行波叠加而成。

**注**　① 作 Fourier 变换，使得偏微分方程的定解问题 (10) 化为常微分方程的初值问题 (12)。对有界问题用分离变量法，也是把偏微分方程化作常微分方程来求解。两者的精神一致。② 由式 (14) 可以看出，$u(x,t)$ 是一系列“本征振动”的叠加（对 $k$ 积分），每一本征振动是两个因子的乘积，一个只是 $t$ 的函数，一个只是 $x$ 的函数，这与用分离变量法求解有界弦的振动所得结果形式上相同。可见 Fourier 变换法本质上也是分离变量法。③ 可以由延迟定理和线性定理直接得到式 (15a) 的结果，由延迟定理、线性定理和积分定理直接得到式 (15b) 的结果。

### §2.2 行波法简介

d’Alembert 公式也可以由以下介绍的行波法得到。其思路是先求出偏微分方程的通解，其中含有两个任意函数，然后用初始条件确定它们。对于无界弦的振动，这一方法非常简单。但是，这种方法不具有一般性，所以我们只在这里作一个简单的介绍。

令

$$
\xi=x-at,\qquad \eta=x+at. \tag{18}
$$

经过简单的计算，可以将式 (10a) 化为

$$
\frac{\partial^2u}{\partial\xi\partial\eta}=0. \tag{19}
$$

上式也就是 $(\partial/\partial\eta)(\partial u/\partial\xi)=0$，所以 $\partial u/\partial\xi$ 与 $\eta$ 无关，即

$$
\frac{\partial u}{\partial\xi}=f_0(\xi), \tag{20}
$$

其中 $f_0(\xi)$ 是 $\xi$ 的任意函数。再积分，即得

$$
u=f(\xi)+g(\eta), \tag{21}
$$

其中 $f(\xi)$ 和 $g(\eta)$ 都是任意函数。重新用变量 $x,t$ 写出来，上式就是

$$
u(x,t)=f(x-at)+g(x+at). \tag{22}
$$

这就是一维自由波动方程在未曾考虑任何定解条件时的最一般解，即真正的通解。可惜能够这样求出通解的方程甚少。另外，即使能够求出通解，如果定解条件比较复杂，要决定其中的任意函数也绝非易事。比如对于有限区间上的自由振动，要由定解条件确定上式中的任意函数就有一定的困难。不过，对于目前的定解问题，上述任意函数是不难确定的。

将初始条件 (10b) 和 (10c) 代入通解 (22)，得到

$$
f(x)+g(x)=\varphi(x),\qquad -af'(x)+ag'(x)=\psi(x). \tag{23}
$$

对第二式积分并与第一式结合可得

$$
f(x)=\frac12\varphi(x)-\frac{1}{2a}\int_{x_0}^{x}\psi(\xi)d\xi+c,
\qquad
g(x)=\frac12\varphi(x)+\frac{1}{2a}\int_{x_0}^{x}\psi(\xi)d\xi-c, \tag{24}
$$

其中 $x_0$ 可视方便取定，$c$ 是任意常数。代入 (22)，即得 d’Alembert 公式 (16)，结果与 $x_0$ 和 $c$ 无关。

**习题**　用 Fourier 变换法求解下列定解问题

$$
\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}=0,
\qquad (-\infty<x<+\infty,\quad y>0),
$$

$$
u\big|_{y=0}=\varphi(x),\qquad u\big|_{y=+\infty}=0,
$$

$$
u\big|_{x=\pm\infty}=0,
\qquad
\left.\frac{\partial u}{\partial x}\right|_{x=\pm\infty}=0.
$$

## §3 $\delta$ 函数

### §3.1 定义一

考虑位于 $x$ 轴上 $x=a$ 处具有单位质量的质点。为了描述它的质量密度（这里指的是线密度，下同），引入函数

$$
\delta_\varepsilon(x-a)=
\begin{cases}
\dfrac{1}{2\varepsilon},& |x-a|<\varepsilon,\\
0,& |x-a|>\varepsilon.
\end{cases} \tag{25}
$$

显然，它描述的是均匀分布在区间 $[a-\varepsilon,a+\varepsilon]$ 内的单位质量物体的质量密度。当 $\varepsilon\to0$，就得到上述质点的质量密度

$$
\delta(x-a)=\lim_{\varepsilon\to0}\delta_\varepsilon(x-a). \tag{26}
$$

它满足

$$
\delta(x-a)=
\begin{cases}
\infty,& x=a,\\
0,& x\ne a,
\end{cases} \tag{27a}
$$

和

$$
\int_{-\infty}^{+\infty}\delta(x-a)\,dx=1. \tag{27b}
$$

这就是 $\delta$ 函数的定义。如果上述质点的质量为 $m$，则其质量密度显然为

$$
\rho(x)=m\delta(x-a). \tag{28}
$$

$\delta$ 函数是数学物理中的重要概念，一般认为最早由物理学家 Dirac 引入。实际上，在 Dirac 之前，Heaviside 已经给出了明确的定义，至于这一函数的雏形，则具有更长的历史。不过，它在近代物理中才获得了广泛的应用。这一函数显然可以描述任何处于一点、总量为 1 的物理量的密度，比如单位质点的质量密度（如上所述），单位点电荷的电荷密度，集中在一点的单位磁通的磁感应强度，等等。如果自变量是时间 $t$ 而不是位置 $x$，那么它描述的就是存在于瞬时、总量为 1 的物理量的强度，比如具有单位冲量的瞬时力，具有单位电荷的瞬时电流，等等。

容易看出，$\delta$ 函数与我们过去所熟悉的函数颇为不同。后者通常是连续或分段连续的，起码在定义域内的每一点应该有确定的函数值，虽然也可能有奇点，但奇点通常不在定义域内。比如 $f(x)=1/x$，虽然 $x\to0$ 时，$f(x)\to\infty$，但 $x=0$ 并不在定义域内。而 $\delta$ 函数在 $x=a$ 处的函数值是 $\infty$，而且这一点对 $\delta$ 函数来说是最重要的。所以 $\delta$ 函数不是一个经典的函数，而是一个广义函数。严格的广义函数理论超出了本书的论题。我们总是从极限的意义上对 $\delta$ 函数作直观的理解，即把 $\delta$ 函数当作某种经典函数序列的极限（这不同于一个函数在某点的极限）。上面所述就是一个典型的例子，后面还会看到其它例子。

另一点值得注意的是，$\delta$ 函数的宗量和函数值都可以是有量纲的，由式 (27b) 容易看出，$\delta$ 函数的量纲是其宗量量纲的倒数。这也不同于我们所熟悉的函数，比如三角函数、指数函数等，它们的宗量必须是没有量纲的，其函数值当然也没有量纲。

### §3.2 定义二

由式 (27) 可以推出，对于任何连续函数 $f(x)$，有

$$
\int_{-\infty}^{+\infty}f(x)\delta(x-a)\,dx=f(a). \tag{29}
$$

证明如下：由于 $\delta(x-a)$ 只在 $a$ 点不为零，故有

$$
\int_{-\infty}^{+\infty}f(x)\delta(x-a)\,dx
=\int_{a-\varepsilon}^{a+\varepsilon}f(x)\delta(x-a)\,dx
=f(a+\theta\varepsilon)\int_{a-\varepsilon}^{a+\varepsilon}\delta(x-a)\,dx
=f(a+\theta\varepsilon). \tag{30}
$$

其中第二步用了积分中值定理，其中 $-1\le\theta\le1$。由于 $\varepsilon$ 的大小可以任意，令 $\varepsilon\to0$，就得到式 (29)。

上述证明实际上是会引起质疑的，因为其中引用了积分中值定理，而后者只适用于经典的函数（定理要求与 $f(x)$ 相乘的函数在积分区间上不变号且可积）。但是，如果我们用 $\delta_\varepsilon(x-a)$ 代替式 (30) 中的 $\delta(x-a)$，那么其中每一步都是严格成立的，于是就得到

$$
\int_{-\infty}^{+\infty}f(x)\delta_\varepsilon(x-a)\,dx=f(a+\theta\varepsilon). \tag{31}
$$

然后令 $\varepsilon\to0$，就得到式 (29)。这里读者可能会提出这样的问题：对上式左边取 $\varepsilon\to0$ 的极限，得到的是 $\lim_{\varepsilon\to0}\int_{-\infty}^{+\infty}f(x)\delta_\varepsilon(x-a)\,dx$，它是否就等于 $\int_{-\infty}^{+\infty}f(x)\delta(x-a)\,dx$？换句话说，取极限与积分是否可以交换次序？实际上，正确的理解应该是，$\int_{-\infty}^{+\infty}f(x)\delta(x-a)\,dx$ 正是用 $\lim_{\varepsilon\to0}\int_{-\infty}^{+\infty}f(x)\delta_\varepsilon(x-a)\,dx$ 来定义的。使用 $\delta$ 函数，使得我们可以省略用函数序列计算然后再取极限的过程。需要指出的是，这里的函数 $\delta_\varepsilon(x-a)$ 不一定是式 (25) 中定义的函数，只要它在 $|x-a|>\varepsilon$ 时为零，在 $|x-a|<\varepsilon$ 时为正，且积分为 1 即可。

今后我们就把式 (29) 当作 $\delta$ 函数的另一种定义，但仍然要从极限的意义上来理解。

很容易将上述定义推广到三维或其它维数的情况，比如三维 $\delta$ 函数 $\delta(\boldsymbol r-\boldsymbol a)$ 可以由下式定义

$$
\int f(\boldsymbol r)\delta(\boldsymbol r-\boldsymbol a)\,d\boldsymbol r=f(\boldsymbol a). \tag{32}
$$

记 $\boldsymbol r=(x,y,z)$，$\boldsymbol a=(a_1,a_2,a_3)$，我们也可以等价地定义

$$
\delta(\boldsymbol r-\boldsymbol a)
=\delta(x-a_1)\delta(y-a_2)\delta(z-a_3). \tag{33}
$$

### §3.3 $\delta$ 函数的基本性质

$\delta$ 函数的基本性质主要有以下几条。

1. 若 $f(x)$ 为连续函数，则 $f(x)\delta(x-a)=f(a)\delta(x-a)$。
2. $x\delta(x)=0$。
3. $\delta$ 函数是偶函数，即 $\delta(-x)=\delta(x)$。
4. $\delta(ax)=\delta(x)/|a|$，其中 $a\ne0$。

**注**　性质 2 可以写成 $\delta(x)/x^{-1}=0$，这说明 $\delta(x)$ 在 $x=0$ 处的奇性弱于 $1/x$。由此可见，$\delta$ 函数的奇性其实并不是很大。

**证明**　设 $\varphi(x)$ 为任意连续函数，根据 $\delta$ 函数的定义二，容易得到

$$
\begin{aligned}
&\int_{-\infty}^{+\infty}\varphi(x)[f(x)\delta(x-a)-f(a)\delta(x-a)]\,dx\\
&=\int_{-\infty}^{+\infty}[\varphi(x)f(x)]\delta(x-a)\,dx
-f(a)\int_{-\infty}^{+\infty}\varphi(x)\delta(x-a)\,dx\\
&=\varphi(a)f(a)-f(a)\varphi(a)=0.
\end{aligned}
$$

由于 $\varphi(x)$ 为任意连续函数，故左边方括号中的函数必须为零，于是就得到性质 1。其它性质可以用同样的方法加以证明，请读者自己完成。

### §3.4 $\delta$ 函数的几种表达式

1.
$$
\delta(x)=\frac{1}{2\pi}\int_{-\infty}^{+\infty}e^{ikx}\,dk.
$$

2.
$$
\delta(x)=\lim_{K\to+\infty}\frac{\sin Kx}{\pi x}.
$$

3.
$$
\delta(x)=\lim_{\sigma\to0}\rho_\sigma(x),\qquad
\rho_\sigma(x)=\frac{1}{\sqrt{2\pi}\sigma}e^{-x^2/2\sigma^2}
$$

是 Gauss 分布，$\sigma>0$。

4.
$$
\delta(x)=\lim_{b\to0}\rho_b(x),\qquad
\rho_b(x)=\frac{b}{\pi(x^2+b^2)}
$$

是 Lorentz 分布，$b>0$。

5.
$$
\delta(x)=\theta'(x),\qquad
\theta(x)=
\begin{cases}
1,& x>0,\\
0,& x<0.
\end{cases}
$$

将以上各式中的 $x$ 换为 $x-a$，即可得到 $\delta(x-a)$ 的相应表达式。

下面对各性质作论证、讨论或补充说明。

1. $\delta(x)$ 的 Fourier 变换为

$$
\Delta(k)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}\delta(x)e^{-ikx}\,dx
=\frac{1}{\sqrt{2\pi}}. \tag{34}
$$

由反变换公式，即得

$$
\delta(x)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}\Delta(k)e^{ikx}\,dk
=\frac{1}{2\pi}\int_{-\infty}^{+\infty}e^{ikx}\,dk. \tag{35}
$$

实际上，上式右边的积分并不存在，那么如何理解这一结果呢？如果我们用式 (25) 中的 $\delta_\varepsilon(x)$ 代替式 (34) 中的 $\delta(x)$，容易得到

$$
\Delta_\varepsilon(k)=\frac{1}{\sqrt{2\pi}}\frac{\sin k\varepsilon}{k\varepsilon}, \tag{36}
$$

由反变换公式，就得到

$$
\delta_\varepsilon(x)=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}\Delta_\varepsilon(k)e^{ikx}\,dk. \tag{37}
$$

由式 (36) 容易看出，$\lim_{\varepsilon\to0}\Delta_\varepsilon(k)\to1/\sqrt{2\pi}$，对式 (37) 取 $\varepsilon\to0$ 的极限，就得到式 (35)。这里可能仍会有读者问，对式 (37) 右边取 $\varepsilon\to0$ 的极限，是否真的可以得到式 (35) 的右边？也就是说，取极限与积分是否可以交换次序？答案仍然是，式 (35) 的右边正是用式 (37) 右边（或类似表达式）的极限来定义的。

2. 按式 (37)，有

$$
\begin{aligned}
\delta(x)
&=\frac{1}{\sqrt{2\pi}}\lim_{\varepsilon\to0}\int_{-\infty}^{+\infty}\Delta_\varepsilon(k)e^{ikx}\,dk\\
&=\frac{1}{\sqrt{2\pi}}\lim_{\varepsilon\to0}\lim_{K\to+\infty}\int_{-K}^{K}\Delta_\varepsilon(k)e^{ikx}\,dk\\
&=\frac{1}{\sqrt{2\pi}}\lim_{K\to+\infty}\lim_{\varepsilon\to0}\int_{-K}^{K}\Delta_\varepsilon(k)e^{ikx}\,dk,
\end{aligned}
$$

其中交换了两个极限的顺序，对于这一操作的合理性我们不作深究。最后的积分是普通的定积分，只要 $\delta_\varepsilon(x)$ 是 $x$ 和 $\varepsilon$ 的连续函数，则 $\Delta_\varepsilon(k)$ 是 $k$ 和 $\varepsilon$ 的连续函数，于是取 $\varepsilon\to0$ 的极限与积分可以交换次序，而得

$$
\delta(x)=\frac{1}{2\pi}\lim_{K\to+\infty}\int_{-K}^{K}e^{ikx}\,dk
=\frac{1}{\pi}\lim_{K\to+\infty}\int_0^K\cos kx\,dk
=\lim_{K\to+\infty}\frac{\sin Kx}{\pi x}.
$$

3. 当 $x\ne0$，$\lim_{\sigma\to0}\rho_\sigma(x)=0$；当 $x=0$，$\lim_{\sigma\to0}\rho_\sigma(x)=\lim_{\sigma\to0}(1/\sqrt{2\pi}\sigma)=\infty$；由于 $\rho_\sigma(x)$ 是概率分布，它已经归一化，即满足 $\int_{-\infty}^{+\infty}\rho_\sigma(x)\,dx=1$（这也可以直接验证）。因此它满足 $\delta$ 函数定义一中的所有条件。

4. 证明同上。事实上，任何概率分布取其宽度趋于零（从而高度趋于无穷大）的极限，都可以得到 $\delta$ 函数。

5. 由读者自己完成，证明方法可参照第三小节。

### §3.5 *$\delta$ 函数的导数

设函数 $f(x)$ 具有连续导数，则形式上可以计算以下积分

$$
\int_{-\infty}^{+\infty}f(x)\delta'(x-a)\,dx
=f(x)\delta(x-a)\Big|_{-\infty}^{+\infty}
-\int_{-\infty}^{+\infty}f'(x)\delta(x-a)\,dx
=-f'(a).
$$

正如前面用式 (29) 来定义 $\delta$ 函数，我们可以通过下列积分来定义 $\delta$ 函数的导数：

$$
\int_{-\infty}^{+\infty}f(x)\delta'(x-a)\,dx=-f'(a). \tag{38}
$$

类似地，还可以由下列积分定义 $\delta$ 函数的高阶导数

$$
\int_{-\infty}^{+\infty}f(x)\delta^{(n)}(x-a)\,dx=(-)^nf^{(n)}(a). \tag{39}
$$

其中函数 $f(x)$ 应该具有连续的 $n$ 阶导数。

$\delta$ 函数的导数具有下列性质：

1. $\delta'(x)$ 是奇函数，即 $\delta'(-x)=-\delta'(x)$。
2. $x\delta'(x)=-\delta(x)$，或 $\delta'(x)=-\delta(x)/x$。

这些性质的证明方法与第三小节类似，读者可以尝试自己做出证明。下面的证明供参考。

1. 设 $\varphi(x)$ 为具有连续导数的任意函数，根据 $\delta$ 函数的导数定义，容易得到

$$
\begin{aligned}
\int_{-\infty}^{+\infty}\varphi(x)[\delta'(-x)+\delta'(x)]\,dx
&=\int_{-\infty}^{+\infty}\varphi(-y)\delta'(y)\,dy-\varphi'(0)\\
&=-\left.\frac{d\varphi(-y)}{dy}\right|_{y=0}-\varphi'(0)
=\left.\frac{d\varphi(x)}{dx}\right|_{x=0}-\varphi'(0)\\
&=\varphi'(0)-\varphi'(0)=0.
\end{aligned}
$$

由于 $\varphi(x)$ 的任意性，左边方括号中的函数必须为零，于是就得到性质 1。

这里需要注意的是，$\delta'(-x)=[d\delta(y)/dy]_{y=-x}=-d\delta(-x)/dx\ne d\delta(-x)/dx$。

2. 同上可得

$$
\begin{aligned}
\int_{-\infty}^{+\infty}\varphi(x)[x\delta'(x)+\delta(x)]\,dx
&=\int_{-\infty}^{+\infty}[x\varphi(x)]\delta'(x)\,dx+\varphi(0)\\
&=-[x\varphi(x)]'_{x=0}+\varphi(0)=-\varphi(0)+\varphi(0)=0.
\end{aligned}
$$

于是得到性质 2。

### §3.6 *复合 $\delta$ 函数

复合 $\delta$ 函数指的是宗量为函数的 $\delta$ 函数，形如 $\delta[\varphi(x)]$。显然，它只在 $\varphi(x)$ 的零点处不为零，故可分解为一系列简单的 $\delta$ 函数。

设 $\varphi(x)$ 连续且只有离散的单根 $x_k$（$k=1,2,\cdots$），则有

$$
\delta[\varphi(x)]=\sum_k\frac{1}{|\varphi'(x_k)|}\delta(x-x_k). \tag{40}
$$

**证明**　设 $f(x)$ 是任意连续函数，由于 $\delta[\varphi(x)]$ 只在 $x_k$（$k=1,2,\cdots$）处不为零，故有

$$
\int_{-\infty}^{+\infty}f(x)\delta[\varphi(x)]\,dx
=\sum_k\int_{x_k-\varepsilon}^{x_k+\varepsilon}f(x)\delta[\varphi(x)]\,dx. \tag{41}
$$

只要 $\varepsilon$ 足够小，则在区间 $[x_k-\varepsilon,x_k+\varepsilon]$ 内，$\varphi(x)$ 是单调增或单调减的，从而在该区间内，其反函数存在，记 $u=\varphi(x)$，则 $x=\varphi^{-1}(u)$，且 $\varphi^{-1}(0)=x_k$。于是

$$
\int_{x_k-\varepsilon}^{x_k+\varepsilon}f(x)\delta[\varphi(x)]\,dx
=\int_{\varphi(x_k-\varepsilon)}^{\varphi(x_k+\varepsilon)}\frac{f(x)}{\varphi'(x)}\delta(u)\,du
=\epsilon[\varphi'(x_k)]\frac{f(x_k)}{\varphi'(x_k)}
=\frac{f(x_k)}{|\varphi'(x_k)|},
$$

其中 $\epsilon[\varphi'(x_k)]$ 是符号函数，即当 $\varphi'(x_k)>0$（$<0$）时，它等于 $1$（$-1$）。如果 $\varphi'(x_k)>0$，则 $u=\varphi(x)$ 在区间 $[x_k-\varepsilon,x_k+\varepsilon]$ 内单调增，上述积分结果就是 $f(x_k)/\varphi'(x_k)=f(x_k)/|\varphi'(x_k)|$；如果 $\varphi'(x_k)<0$，则 $u=\varphi(x)$ 在区间 $[x_k-\varepsilon,x_k+\varepsilon]$ 内单调减，此时应该将积分的上下限对调，然后才能用定义 (29)，从而积分结果就是 $-f(x_k)/\varphi'(x_k)=f(x_k)/|\varphi'(x_k)|$。将结果代入式 (41)，即得

$$
\int_{-\infty}^{+\infty}f(x)\delta[\varphi(x)]\,dx
=\sum_k\frac{f(x_k)}{|\varphi'(x_k)|}
=\int_{-\infty}^{+\infty}f(x)\sum_k\frac{\delta(x-x_k)}{|\varphi'(x_k)|}\,dx. \tag{42}
$$

由 $f(x)$ 的任意性，就得到式 (40)。证毕。

由式 (40) 容易得到以下几个特例：

1.
$$
\delta(ax)=\frac{\delta(x)}{|a|}.
$$

2.
$$
\delta(x^2-a^2)=\frac{\delta(x-a)+\delta(x+a)}{2|a|}
=\frac{\delta(x-a)+\delta(x+a)}{2|x|}.
$$

3.
$$
\delta(x^2)=\frac{\delta(x)}{|x|}.
$$

其中特例 3 可以在特例 2 中令 $a=0$ 得到。

**习题**　证明 $\theta'(x)=\delta(x)$。

## §4 无界杆的热传导问题

考虑无穷长导热细杆在热源和初始温度作用下的热传导问题。这是一个一维无界问题。定解问题如下：

$$
\frac{\partial u}{\partial t}-a^2\frac{\partial^2u}{\partial x^2}=f,
\qquad (-\infty<x<+\infty,\quad t>0), \tag{43a}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad (-\infty<x<+\infty). \tag{43b}
$$

**注**　这里忽略了默认的边界条件 $u|_{x=\pm\infty}=0$。关于有界问题抽象为无界问题的考虑与弦振动情况类似。

上述定解问题可以分解为两个定解问题，即由初始温度引起的热传导（无源问题）

$$
\frac{\partial u}{\partial t}-a^2\frac{\partial^2u}{\partial x^2}=0,
\qquad (-\infty<x<+\infty,\quad t>0), \tag{44a}
$$

$$
u\big|_{t=0}=\varphi(x),
\qquad (-\infty<x<+\infty), \tag{44b}
$$

以及由热源引起的热传导（有源问题）

$$
\frac{\partial u}{\partial t}-a^2\frac{\partial^2u}{\partial x^2}=f,
\qquad (-\infty<x<+\infty,\quad t>0), \tag{45a}
$$

$$
u\big|_{t=0}=0,
\qquad (-\infty<x<+\infty). \tag{45b}
$$

容易验证，把定解问题 (44) 和 (45) 的解加起来，就得到定解问题 (43) 的解。所以我们将两个问题分开考虑。

### §4.1 无源问题

先考虑定解问题 (44)。对 $x$ 作 Fourier 变换，设

$$
u(x,t)\leftrightarrow U(k,t),\qquad \varphi(x)\leftrightarrow\Phi(k). \tag{46}
$$

对定解问题 (44) 的各式作 Fourier 变换，利用微分定理，可得常微分方程初值问题

$$
\frac{dU}{dt}+k^2a^2U=0, \tag{47a}
$$

$$
U\big|_{t=0}=\Phi(k). \tag{47b}
$$

容易求得上式的解为

$$
U(k,t)=\Phi(k)e^{-k^2a^2t}. \tag{48}
$$

剩下的问题是对 $U(k,t)$ 作反变换求出原函数 $u(x,t)$。由卷积定理，只要求出 $e^{-k^2a^2t}$ 的原函数，并与 $\varphi(x)$ 作卷积即可。由反变换公式，

$$
\begin{aligned}
\mathcal F^{-1}(e^{-k^2a^2t})
&=\frac{1}{\sqrt{2\pi}}\int_{-\infty}^{+\infty}e^{-k^2a^2t}e^{ikx}\,dk
=\sqrt{\frac{2}{\pi}}\int_0^{+\infty}e^{-k^2a^2t}\cos kx\,dk\\
&=\sqrt{\frac{2}{\pi}}\cdot\frac12\sqrt{\frac{\pi}{a^2t}}e^{-x^2/4a^2t}
=\frac{1}{a\sqrt{2t}}e^{-x^2/4a^2t}. \tag{49}
\end{aligned}
$$

其中利用了积分公式

$$
\int_0^\infty e^{-ax^2}\cos bx\,dx
=\frac12\sqrt{\frac{\pi}{a}}e^{-b^2/4a},\qquad a>0. \tag{50}
$$

由卷积的表达式 (9) 即得

$$
u(x,t)=\frac{1}{2a\sqrt{\pi t}}\int_{-\infty}^{+\infty}\varphi(\xi)
\exp\left[-\frac{(x-\xi)^2}{4a^2t}\right]d\xi. \tag{51}
$$

今考虑特例 $\varphi(x)=\delta(x-x_0)$，即初始时刻有一定的热量集中在 $x_0$ 处，按上式，则 $t$ 时刻的温度分布为

$$
u(x,t)=\frac{1}{2a\sqrt{\pi t}}
\exp\left[-\frac{(x-x_0)^2}{4a^2t}\right]. \tag{52}
$$

上式为 Gauss 分布，其宽度随着 $t$ 的增大而增大，高度则随着 $t$ 的增大而减小，但对 $x$ 积分始终为 1。这表明热量不断向远处传播，但总能量守恒。但是，不管 $t$ 多么小，只要 $t>0$，杆上各处的温度均不为零，这表明温度的传播速度为无穷大。这显然不符合实际情况。引起这一结果的原因是推导热传导方程时没有考虑热传导过程的微观机制，即所用模型过于简化。不过，当 $t$ 较大以后，由热传导方程解得的结果还是符合实际情况的。

### §4.2 *有源问题

再考虑定解问题 (45)。对 $x$ 作 Fourier 变换，设

$$
u(x,t)\leftrightarrow U(k,t),\qquad f(x,t)\leftrightarrow F(k,t). \tag{53}
$$

对定解问题 (45) 的各式作 Fourier 变换，利用微分定理，可得常微分方程初值问题

$$
\frac{dU}{dt}+k^2a^2U=F, \tag{54a}
$$

$$
U\big|_{t=0}=0. \tag{54b}
$$

容易求得上式的解为

$$
U(k,t)=\int_0^tF(k,\tau)e^{-k^2a^2(t-\tau)}\,d\tau. \tag{55}
$$

于是

$$
u(x,t)=\mathcal F^{-1}[U(k,t)]
=\int_0^t\mathcal F^{-1}[F(k,\tau)e^{-k^2a^2(t-\tau)}]d\tau, \tag{56}
$$

其中利用了 Fourier 反变换的线性性，或者说，交换了对 $k$ 积分和对 $\tau$ 积分的次序。利用式 (49) 和卷积定理，即得

$$
u(x,t)=\int_0^t d\tau\int_{-\infty}^{+\infty}d\xi\,
\frac{f(\xi,\tau)}{2a\sqrt{\pi(t-\tau)}}
\exp\left[-\frac{(x-\xi)^2}{4a^2(t-\tau)}\right]. \tag{57}
$$

如果 $f(x,t)=\delta(x-x_0)\delta(t-t_0)$，即热源只在 $t_0$ 时刻作用在 $x_0$ 点，容易得到

$$
u(x,t)=\frac{\theta(t-t_0)}{2a\sqrt{\pi(t-t_0)}}
\exp\left[-\frac{(x-x_0)^2}{4a^2(t-t_0)}\right]. \tag{58}
$$

这一结果的物理意义是容易理解的：由于初始温度为零，而热源在 $t_0$ 时刻才发生作用，故 $t_0$ 时刻以前温度一直为零；$t_0$ 时刻由于热源的作用，使得杆上 $x_0$ 点获得一定的热量，所以 $t_0$ 时刻以后的温度分布与式 (52) 类似。实际上，只要将式 (52) 中的 $t$ 换为 $t-t_0$，即得式 (58)。

**习题**　半无界导热细杆，其坐标取为 $0<x<\infty$，$x=0$ 端温度保持零度，杆上初始温度分布是 $u|_{t=0}=\varphi(x)$，其中 $\varphi(0)=0$。求 $t>0$ 以后的温度分布。并具体计算特例 $\varphi(x)=\delta(x-x_0)$，其中 $x_0>0$。尝试画出几个时刻的温度分布草图。

## §5 *无界弦的受迫振动

类似于上节将定解问题 (43) 分解为定解问题 (44) 和 (45)，无界弦的受迫振动的一般定解问题也可以分解为两个定解问题来求解，即纯粹由初始激励引起的振动和纯粹由外力引起的振动。前者已经在 §2 解出，现在求解后者。其定解问题为

$$
\frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=f,
\qquad (-\infty<x<+\infty,\quad t>0), \tag{59a}
$$

$$
u\big|_{t=0}=0,
\qquad (-\infty<x<+\infty), \tag{59b}
$$

$$
\left.\frac{\partial u}{\partial t}\right|_{t=0}=0,
\qquad (-\infty<x<+\infty). \tag{59c}
$$

这里同样是忽略了默认的边界条件 $u|_{x=\pm\infty}=0$（因而 $\partial u/\partial x|_{x=\pm\infty}=0$）。

对 $x$ 作 Fourier 变换，设

$$
u(x,t)\leftrightarrow U(k,t),\qquad f(x,t)\leftrightarrow F(k,t). \tag{60}
$$

对以上定解问题的各式作 Fourier 变换，利用微分定理，可得常微分方程初值问题

$$
\frac{d^2U}{dt^2}+k^2a^2U=F, \tag{61a}
$$

$$
U\big|_{t=0}=0,
\qquad
\left.\frac{dU}{dt}\right|_{t=0}=0. \tag{61b}
$$

可以求得上式的解为

$$
U(k,t)=\int_0^tF(k,\tau)\frac{\sin ka(t-\tau)}{ka}\,d\tau. \tag{62}
$$

**注**　上式可以用常数变易法求出，但计算较为繁琐。如果用 Laplace 变换法，则计算可以大大简化。Laplace 变换法是求解常微分方程初值问题的最方便的方法。它还可以用于求解偏微分方程的定解问题。读者可能已经想到，Laplace 变换是对时间 $t$ 作变换。有兴趣的读者可以参考其它教材。

作 Fourier 反变换，利用其线性性（或者交换对 $k$ 积分和对 $\tau$ 积分的次序），可得

$$
u(x,t)=\mathcal F^{-1}[U(k,t)]
=\int_0^t\mathcal F^{-1}\left[F(k,\tau)\frac{\sin ka(t-\tau)}{ka}\right]d\tau. \tag{63}
$$

下面的计算与式 (15b) 的计算过程类似，

$$
\begin{aligned}
\mathcal F^{-1}\left[F(k,\tau)\frac{\sin ka(t-\tau)}{ka}\right]
&=\mathcal F^{-1}\left[\int_0^{t-\tau}F(k,\tau)\cos ka\eta\,d\eta\right]\\
&=\int_0^{t-\tau}\mathcal F^{-1}[F(k,\tau)\cos ka\eta]d\eta\\
&=\frac12\int_0^{t-\tau}\mathcal F^{-1}[F(k,\tau)e^{-ika\eta}+F(k,\tau)e^{ika\eta}]d\eta\\
&=\frac12\int_0^{t-\tau}[f(x-a\eta,\tau)+f(x+a\eta,\tau)]d\eta\\
&=\frac{1}{2a}\int_{x-a(t-\tau)}^{x+a(t-\tau)}f(\xi,\tau)d\xi. \tag{64}
\end{aligned}
$$

代入式 (63)，即得

$$
u(x,t)=\frac{1}{2a}\int_0^t d\tau\int_{x-a(t-\tau)}^{x+a(t-\tau)}d\xi\,f(\xi,\tau). \tag{65}
$$

由上式可以看出，只有 $t$ 时刻以前的外力才对 $t$ 时刻的振动有影响；$\tau$ 时刻作用于 $\xi$ 处的外力，在 $t$ 时刻所能影响到的位置 $x$ 满足 $|x-\xi|\le a(t-\tau)$，其中离 $\xi$ 最远的点距离 $\xi$ 为 $a(t-\tau)$，可见振动的传播速度是 $a$。

## 补充习题

1. 半无界导热细杆，其坐标取为 $0<x<\infty$，$x=0$ 端有已知热流 $q_0\sin\omega t$ 进入，其中 $q_0$ 和 $\omega$ 是已知常数。杆上初始温度分布是 $u|_{t=0}=\varphi(x)$。求 $t>0$ 以后的温度分布，并研究长时间以后的情况。

   提示：设温度分布为 $u(x,t)$，令 $u(x,t)=u_0(x,t)+v(x,t)$，其中

   $$
   u_0(x,t)=-\frac{q_0a}{k\sqrt{\omega}}
   \exp\left(-\sqrt{\frac{\omega}{2}}\frac{x}{a}\right)
   \cos\left(\omega t-\sqrt{\frac{\omega}{2}}\frac{x}{a}+\frac{\pi}{4}\right)
   $$

   则 $v(x,t)$ 的定解问题满足齐次方程和齐次边界条件。对 $v(x,t)$ 的初始条件作偶延拓，即可用现成的公式求解 $v(x,t)$。$u_0(x,t)$ 的解法如下。设

   $$
   U_0(x,t)=i\frac{q_0}{k}F(x)e^{i\omega t},\qquad F(x)=Ae^{\lambda x}.
   $$

   要求 $U_0(x,t)$ 满足 $u(x,t)$ 的方程，$F(x)$ 在 $x\to+\infty$ 时有界，并满足边界条件

   $$
   \left.\frac{\partial U_0}{\partial x}\right|_{x=0}=i\frac{q_0}{k}e^{i\omega t}.
   $$

   则可以定出

   $$
   F(x)=-e^{-i\pi/4}\frac{a}{\sqrt{\omega}}
   \exp\left(-e^{i\pi/4}\sqrt{\frac{\omega}{a}}\,x\right).
   $$

   从而得出 $U_0(x,t)$。然后取 $u_0(x,t)=\operatorname{Re}U_0(x,t)$ 即可。

2. 求解三维无界空间中的输运问题：

   $$
   \frac{\partial u}{\partial t}-a^2\nabla^2u=0,\qquad u\big|_{t=0}=\varphi(\boldsymbol r).
   $$

   答案：

   $$
   u(\boldsymbol r,t)=\frac{1}{(2a)^3(\pi t)^{3/2}}
   \int\varphi(\boldsymbol\xi)
   \exp\left[-\frac{(\boldsymbol r-\boldsymbol\xi)^2}{4a^2t}\right]d\boldsymbol\xi.
   $$

3. 求解半无界弦的振动问题：

   $$
   \frac{\partial^2u}{\partial t^2}-a^2\frac{\partial^2u}{\partial x^2}=0
   \qquad (0<x<\infty,\ t>0),
   $$

   $$
   u\big|_{t=0}=0,\qquad
   \left.\frac{\partial u}{\partial t}\right|_{t=0}=0,\qquad
   u\big|_{x=0}=f(t).
   $$

   （提示：利用零初始条件将问题延拓到 $-\infty<t<+\infty$，对 $t$ 作 Fourier 变换。求解时注意只有右行波。）
