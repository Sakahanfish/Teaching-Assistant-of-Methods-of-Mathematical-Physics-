# *第十二章 线性偏微分方程的通解

上一章我们示范性地推导了常用的几个典型的线性偏微分方程以及定解条件，下面我们讨论如何求解这些偏微分方程定解问题。

求解常微分方程的定解问题时，我们常常先求出常微分方程的通解，然后将通解代入定解条件，定出通解中的任意常数，从而确定常微分方程定解问题的解。按照常微分方程定解问题求解思路，我们也从偏微分方程的通解开始，讨论偏微分方程定解问题的求解问题。

## *§12.1 线性方程解的叠加性

在第十一章中已经导出的几种典型的二阶偏微分方程都是线性偏微分方程，也就是说，在方程中只出现对于未知函数的线性运算。在讨论偏微分方程的通解之前，我们先介绍线性方程共有的性质。下面为了叙述简洁起见，引进线性算符 $\hat L$ 的记号，把这些线性偏微分方程统一写成

$$
\hat L[u]=f
\tag{12.1}
$$

的形式，其中 $u$ 是未知函数，$f$ 是已知函数，称为方程的非齐次项。若 $f=0$，则称方程是齐次的，否则方程就是非齐次的。和上一章中得到的各个偏微分方程做比较，就可以看出线性算符 $\hat L$ 的具体形式。下面的表 12.1 中给出了几个典型的例子。而且，本书讨论的定解条件，也都是线性的，也可以把它们写成类似的算符形式。

表 12.1

| 方程类型 | 方程 | 线性算符 $\hat L$ |
|---|---|---|
| 波动方程 | $\displaystyle \frac{\partial^2u}{\partial t^2}-a^2\nabla^2u=f$ | $\displaystyle \hat L\equiv\frac{\partial^2}{\partial t^2}-a^2\nabla^2$ |
| 热传导方程 | $\displaystyle \frac{\partial u}{\partial t}-\kappa\nabla^2u=f$ | $\displaystyle \hat L\equiv\frac{\partial}{\partial t}-\kappa\nabla^2$ |
| Poisson 方程 | $\displaystyle \nabla^2u=f$ | $\displaystyle \hat L\equiv\nabla^2$ |
| Helmholtz 方程 | $\displaystyle \nabla^2u+k^2u=f$ | $\displaystyle \hat L\equiv\nabla^2+k^2$ |

根据线性算符的定义

$$
\hat L[c_1u_1+c_2u_2]
=
c_1\hat L[u_1]+c_2\hat L[u_2]
\qquad
(c_1,c_2\text{ 为常数}),
$$

立即可以导出线性方程的下列基本性质。

**性质 1**　若 $u_1$ 和 $u_2$ 都是齐次方程 $\hat L[u]=0$ 的解，

$$
\hat L[u_1]=0,
\qquad
\hat L[u_2]=0,
$$

则它们的线性组合 $c_1u_1+c_2u_2$ 也是该齐次方程的解，

$$
\hat L[c_1u_1+c_2u_2]=0.
\tag{12.2}
$$

**性质 2**　若 $u_1$ 和 $u_2$ 都是同一个非齐次方程 $\hat L[u]=f$ 的解，

$$
\hat L[u_1]=f,
\qquad
\hat L[u_2]=f,
$$

则它们的差 $u_1-u_2$ 一定是相应的齐次方程的解，

$$
\hat L[u_1-u_2]=0.
\tag{12.3}
$$

换言之，非齐次方程的一个特解加上相应齐次方程的解仍是原非齐次方程的解。

**性质 3**　若 $u_1$ 和 $u_2$ 分别满足非齐次方程

$$
\hat L[u_1]=f_1,
\qquad
\hat L[u_2]=f_2,
$$

则它们的线性组合 $c_1u_1+c_2u_2$ 满足非齐次方程

$$
\hat L[c_1u_1+c_2u_2]
=
c_1f_1+c_2f_2.
\tag{12.4}
$$

本节和以后几节均以两个自变量的线性偏微分方程为例，讨论方程的特解和通解。这类线性偏微分方程的普遍形式可以写为

$$
A_0\frac{\partial^nu}{\partial x^n}
+
A_1\frac{\partial^nu}{\partial x^{n-1}\partial y}
+\cdots+
A_n\frac{\partial^nu}{\partial y^n}
+
B_0\frac{\partial^{n-1}u}{\partial x^{n-1}}
+\cdots+
M\frac{\partial u}{\partial x}
+
N\frac{\partial u}{\partial y}
+
Pu
=
f(x,y),
\tag{12.5}
$$

或者引进简写符号 $\hat D_x\equiv\partial/\partial x,\ \hat D_y\equiv\partial/\partial y$，而将方程写成

$$
\hat L(\hat D_x,\hat D_y)u
\equiv
\left(
A_0\hat D_x^n
+
A_1\hat D_x^{n-1}\hat D_y
+\cdots+
A_n\hat D_y^n
+
B_0\hat D_x^{n-1}
+\cdots+
M\hat D_x
+
N\hat D_y
+
P
\right)u
=
f(x,y),
\tag{12.5'}
$$

其中 $A_0,A_1,\cdots,A_n,B_0,\cdots,M,N,P$ 都是 $x,y$ 的已知函数，称为方程的系数。我们只讨论最简单的情形，即常系数的线性偏微分方程（方程的系数均为常数），以及能化为常系数线性偏微分方程的方程。

## *§12.2 常系数线性齐次偏微分方程的通解

两个自变量的常系数线性齐次偏微分方程的普遍形式是

$$
A_0\frac{\partial^nu}{\partial x^n}
+
A_1\frac{\partial^nu}{\partial x^{n-1}\partial y}
+\cdots+
A_n\frac{\partial^nu}{\partial y^n}
+
B_0\frac{\partial^{n-1}u}{\partial x^{n-1}}
+\cdots+
M\frac{\partial u}{\partial x}
+
N\frac{\partial u}{\partial y}
+
Pu
=
0,
\tag{12.6}
$$

或者

$$
\hat L(\hat D_x,\hat D_y)u
\equiv
\left(
A_0\hat D_x^n
+
A_1\hat D_x^{n-1}\hat D_y
+\cdots+
A_n\hat D_y^n
+
B_0\hat D_x^{n-1}
+\cdots+
M\hat D_x
+
N\hat D_y
+
P
\right)u
=
0,
\tag{12.6'}
$$

其中方程的系数 $A_0,A_1,\cdots,A_n,B_0,\cdots,M,N,P$ 都是常数。求解这类方程，需要区分两种情形，即 $\hat L(\hat D_x,\hat D_y)$ 是否是 $\hat D_x,\hat D_y$ 的齐次式。

### 1. $\hat L(\hat D_x,\hat D_y)$ 是 $\hat D_x,\hat D_y$ 的齐次式

$$
\left(
A_0\hat D_x^n
+
A_1\hat D_x^{n-1}\hat D_y
+
A_2\hat D_x^{n-2}\hat D_y^2
+\cdots+
A_n\hat D_y^n
\right)u
=
0.
\tag{12.7}
$$

这时，线性算符 $\hat L(\hat D_x,\hat D_y)$ 可以分解成为 $n$ 个线性算符的乘积

$$
\hat L(\hat D_x,\hat D_y)
=
A_0
(\hat D_x-\alpha_1\hat D_y)
(\hat D_x-\alpha_2\hat D_y)
\cdots
(\hat D_x-\alpha_n\hat D_y),
\tag{12.8}
$$

其中 $\alpha_1,\alpha_2,\cdots,\alpha_n$ 也都是常数，因此这 $n$ 个因子的次序可以任意调换。

取试探解为 $u=\phi(y+\alpha x)$，因为

$$
\hat D_x^ku
=
\alpha^k\phi^{(k)}(y+\alpha x),
\qquad
\hat D_y^ku
=
\phi^{(k)}(y+\alpha x),
$$

$$
\hat D_x^r\hat D_y^su
=
\alpha^r\phi^{(r+s)}(y+\alpha x),
$$

代入方程即得

$$
\left(
A_0\alpha^n
+
A_1\alpha^{n-1}
+\cdots+
A_n
\right)
\phi^{(n)}(y+\alpha x)
=
0.
\tag{12.9}
$$

设代数方程（称为附加方程，auxiliary equation）

$$
A_0\alpha^n
+
A_1\alpha^{n-1}
+\cdots+
A_n
=
0
\tag{12.10}
$$

的解是 $\alpha_1,\alpha_2,\cdots,\alpha_n$，且互不相等，则方程 (12.7) 的通解为

$$
u
=
\phi_1(y+\alpha_1x)
+
\phi_2(y+\alpha_2x)
+\cdots+
\phi_n(y+\alpha_nx),
\tag{12.11}
$$

其中 $\phi_i,\ i=1,2,\cdots,n$ 是（互相独立的）任意（$n$ 次可微）函数。

**例 12.1**　求方程

$$
\frac{\partial^2u}{\partial x^2}
-
a^2\frac{\partial^2u}{\partial y^2}
=
0
$$

的通解，$a$ 为常数。

**解**　因附加方程 $\alpha^2-a^2=0$ 的解 $\alpha=\pm a$，故方程的通解为

$$
u
=
\phi_1(y+ax)
+
\phi_2(y-ax).
$$

若 $\alpha$ 是重根，例如是二重根，$(\hat D_x-\alpha\hat D_y)^2u=0$，则通解为

$$
u
=
x\phi_1(y+\alpha x)
+
\phi_2(y+\alpha x).
$$

若 $\alpha$ 为 $n$ 重根，即 $(\hat D_x-\alpha\hat D_y)^nu=0$，则方程的通解为

$$
u
=
x^{n-1}\phi_1(y+\alpha x)
+
x^{n-2}\phi_2(y+\alpha x)
+\cdots+
x\phi_{n-1}(y+\alpha x)
+
\phi_n(y+\alpha x).
$$

**例 12.2**　方程

$$
(\hat D_x^2-2\hat D_x\hat D_y+\hat D_y^2)u=0
$$

的通解为

$$
u
=
x\phi(x+y)+\psi(x+y).
$$

### 2. $\hat L(\hat D_x,\hat D_y)$ 不是 $\hat D_x,\hat D_y$ 的齐次式

首先考虑一阶偏微分方程

$$
(\hat D_x-\alpha\hat D_y-\beta)u=0.
\tag{12.12}
$$

前面已经求出此方程在 $\beta=0$ 时的通解 $u=\phi(y+\alpha x)$。当 $\beta\ne0$ 时可设解为

$$
u(x,y)
=
f(x)\phi(y+\alpha x).
$$

代入方程，有

$$
(\hat D_x-\alpha\hat D_y-\beta)
[f(x)\phi(y+\alpha x)]
=
f(x)(\hat D_x-\alpha\hat D_y)\phi(y+\alpha x)
+
\phi(y+\alpha x)(\hat D_x-\beta)f(x)
=
0.
$$

因为 $(\hat D_x-\alpha\hat D_y)\phi(y+\alpha x)=0$，就得到 $f(x)$ 满足的常微分方程

$$
f'(x)-\beta f(x)=0.
$$

解之得 $f(x)=e^{\beta x}$。因此，方程 (12.12) 的通解就是

$$
u(x,y)
=
e^{\beta x}\phi(y+\alpha x).
\tag{12.13}
$$

在上面的通解中，本来还会出现一个常系数，但它可以吸收进任意函数 $\phi(y+\alpha x)$ 中。

**练习 12.1**　将方程 (12.12) 的通解取为 $u(x,y)=g(y)\phi(y+\alpha x)$，试求解之。

显然，当 $\hat L(\hat D_x,\hat D_y)$ 不是 $\hat D_x$ 和 $\hat D_y$ 的齐次式，但能分解为 $\hat D_x$ 和 $\hat D_y$ 的线性式的乘积时，也可以容易地求出方程的通解。

**例 12.3**　求方程

$$
\frac{\partial^2u}{\partial x^2}
-
\frac{\partial^2u}{\partial x\partial y}
-
2\frac{\partial^2u}{\partial y^2}
+
2\frac{\partial u}{\partial x}
+
2\frac{\partial u}{\partial y}
=
0
$$

的通解。

**解**　容易看出，

$$
(\hat D_x^2-\hat D_x\hat D_y-2\hat D_y^2+2\hat D_x+2\hat D_y)u
=
(\hat D_x+\hat D_y)
(\hat D_x-2\hat D_y+2)u
=
0,
$$

故方程的通解为

$$
u
=
\phi(x-y)
+
e^{-2x}\psi(y+2x).
$$

注意：若有重复性因子，例如

$$
(\hat D_x-\alpha\hat D_y-\beta)^2z=0,
$$

则通解为

$$
z
=
xe^{\beta x}\phi(y+\alpha x)
+
e^{\beta x}\psi(y+\alpha x).
$$

## *§12.3 常系数线性非齐次偏微分方程的通解

对于线性非齐次偏微分方程，

$$
\hat L(\hat D_x,\hat D_y)u
=
f(x,y),
\tag{12.14}
$$

显然有

> 线性非齐次方程的通解 = 线性非齐次方程的任一特解 + 相应线性齐次方程的通解。

因此，问题便转化为只需求出线性非齐次方程的任意一个特解。

将方程 (12.14) 的特解形式地表示为

$$
u_0
=
\frac{1}{\hat L(\hat D_x,\hat D_y)}
f(x,y),
\tag{12.15}
$$

而后可按下列法则求出 $u_0(x,y)$。

(1) 若 $f(x,y)=e^{ax+by}$，则

$$
\boxed{
\frac{1}{\hat L(\hat D_x,\hat D_y)}
e^{ax+by}
=
\frac{1}{L(a,b)}
e^{ax+by},
\qquad
L(a,b)\ne0.
}
\tag{12.16}
$$

直接利用求偏导数的公式

$$
\hat D_xe^{ax+by}
=
ae^{ax+by},
\qquad
\hat D_ye^{ax+by}
=
be^{ax+by}
$$

就能够证明这个公式。

例外的情形是 (12.16) 式中 $L(a,b)=0$。这时，不妨假设

$$
\hat L(\hat D_x,\hat D_y)
=
b\hat D_x-a\hat D_y,
$$

于是便需要求解

$$
(b\hat D_x-a\hat D_y)u
=
e^{ax+by}.
\tag{12.17}
$$

可设特解为

$$
u_0(x,y)
=
f(x,y)e^{ax+by}.
$$

代入方程 (12.17)，可以得到 $f(x,y)$ 满足的微分方程

$$
(b\hat D_x-a\hat D_y)f(x,y)=1.
$$

不妨取此方程的解为

$$
f(x,y)
=
\alpha x+\beta y+\gamma,
$$

其中的系数满足

$$
b\alpha-a\beta=1.
$$

即可。特别是，可取 $\beta=\gamma=0,\ \alpha=1/b$，

$$
\boxed{
\frac{1}{b\hat D_x-a\hat D_y}
e^{ax+by}
=
\frac{x}{b}e^{ax+by},
}
\tag{12.18}
$$

或者取 $\alpha=\gamma=0,\ \beta=-1/a$，

$$
\boxed{
\frac{1}{b\hat D_x-a\hat D_y}
e^{ax+by}
=
-\frac{y}{a}e^{ax+by}.
}
\tag{12.18'}
$$

**例 12.4**　求解方程

$$
(2\hat D_x-3\hat D_y)(\hat D_x+\hat D_y)u
=
5e^{x-y}.
$$

**解**　显然，相应齐次方程的通解为

$$
\phi(y-x)+\psi(2y+3x).
$$

非齐次方程的特解可取为

$$
\begin{aligned}
u_0
&=
\frac{5}{(\hat D_x+\hat D_y)(2\hat D_x-3\hat D_y)}
e^{x-y}\\
&=
\frac1{\hat D_x+\hat D_y}
\left[
\frac5{2-3(-1)}
e^{x-y}
\right]\\
&=
\frac1{\hat D_x+\hat D_y}
e^{x-y}
=
xe^{x-y}.
\end{aligned}
$$

所以，非齐次方程的通解为

$$
u
=
xe^{x-y}
+
\phi(y-x)
+
\psi(2y+3x).
$$

(2) 若 $f(x,y)=e^{i(ax+by)}$，显然有

$$
\frac1{\hat L(\hat D_x,\hat D_y)}
e^{i(ax+by)}
=
\frac1{L(ia,ib)}
e^{i(ax+by)}.
$$

因此，当 $a$ 和 $b$ 为实数，且 $\hat L(\hat D_x,\hat D_y)$ 中的系数也为实数时，

$$
\boxed{
\frac1{\hat L(\hat D_x,\hat D_y)}
\sin(ax+by)
=
\operatorname{Im}
\left[
\frac1{L(ia,ib)}
e^{i(ax+by)}
\right],
}
\tag{12.19}
$$

$$
\boxed{
\frac1{\hat L(\hat D_x,\hat D_y)}
\cos(ax+by)
=
\operatorname{Re}
\left[
\frac1{L(ia,ib)}
e^{i(ax+by)}
\right].
}
\tag{12.20}
$$

(3) 若 $f(x,y)=e^{ax+by}g(x,y)$，则

$$
\boxed{
\frac1{\hat L(\hat D_x,\hat D_y)}
e^{ax+by}g(x,y)
=
e^{ax+by}
\frac1{\hat L(\hat D_x+a,\hat D_y+b)}
g(x,y).
}
\tag{12.21}
$$

要证明这个公式，只需注意到

$$
\hat D_x
[e^{ax+by}f(x,y)]
=
e^{ax+by}(\hat D_x+a)f(x,y),
$$

$$
\hat D_y
[e^{ax+by}f(x,y)]
=
e^{ax+by}(\hat D_y+b)f(x,y),
$$

因此就有

$$
\hat L(\hat D_x,\hat D_y)
[e^{ax+by}f(x,y)]
=
e^{ax+by}
\hat L(\hat D_x+a,\hat D_y+b)f(x,y).
$$

这样，用 $\hat L(\hat D_x,\hat D_y)$ 作用在公式 (12.21) 两端，显然，

$$
\begin{aligned}
&\hat L(\hat D_x,\hat D_y)
\left[
e^{ax+by}
\frac1{\hat L(\hat D_x+a,\hat D_y+b)}
g(x,y)
\right]\\
&\qquad=
e^{ax+by}
\hat L(\hat D_x+a,\hat D_y+b)
\left[
\frac1{\hat L(\hat D_x+a,\hat D_y+b)}
g(x,y)
\right]\\
&\qquad=
e^{ax+by}g(x,y).
\end{aligned}
$$

所以，公式 (12.21) 成立。

(4) 若 $f(x,y)=x^my^n$，则可将 $1/\hat L(\hat D_x,\hat D_y)$ 展开为 $\hat D_x,\hat D_y$ 的幂级数，而后求出特解。

**例 12.5**　求方程

$$
(\hat D_x^2-2\hat D_x\hat D_y+\hat D_y^2)u
=
12xy
$$

的通解。

**解**　方程的特解可取为

$$
\begin{aligned}
u_0
&=
\frac{12}
{\hat D_x^2-2\hat D_x\hat D_y+\hat D_y^2}
xy\\
&=
\frac{12}{(\hat D_x-\hat D_y)^2}xy\\
&=
\frac{12}{\hat D_x^2}
\left(
1-\frac{\hat D_y}{\hat D_x}
\right)^{-2}
xy\\
&=
\frac{12}{\hat D_x^2}
\left(
1+2\frac{\hat D_y}{\hat D_x}+\cdots
\right)xy\\
&=
\frac{12}{\hat D_x^2}
\left(
xy+\frac2{\hat D_x}x
\right)\\
&=
12
\left(
y\frac1{\hat D_x^2}x
+
2\frac1{\hat D_x^3}x
\right)\\
&=
12
\left(
\frac16x^3y+\frac1{12}x^4
\right)
=
x^4+2x^3y,
\end{aligned}
$$

其中利用了

$$
\frac1{\hat D_x^2}x
=
\frac16x^3
\qquad
\left(
\text{因 }
\frac{d^2}{dx^2}\frac{x^3}{6}=x
\right),
$$

$$
\frac1{\hat D_x^3}x
=
\frac1{24}x^4
\qquad
\left(
\text{因 }
\frac{d^3}{dx^3}\frac{x^4}{24}=x
\right).
$$

相应齐次方程的通解已在例 12.2 中求出，故非齐次方程的通解为

$$
u
=
x\phi(x+y)
+
\psi(x+y)
+
x^4+2x^3y.
$$

需要说明，在将 $1/\hat L(\hat D_x,\hat D_y)$ 展开时可以有不同的方法，因而得到不同的结果。例如，在上面的例题中，也可以得到

$$
\frac1{(\hat D_x-\hat D_y)^2}
=
\frac1{\hat D_y^2}
\left(
1-\frac{\hat D_x}{\hat D_y}
\right)^{-2}
=
\frac1{\hat D_y^2}
\left(
1+2\frac{\hat D_x}{\hat D_y}+\cdots
\right),
$$

因此，非齐次方程的特解也可以取为

$$
u_0
=
\frac{12}{(\hat D_x-\hat D_y)^2}xy
=
2xy^3+y^4.
$$

容易验证，这两种方法得到的特解之差

$$
x^4+2x^3y-2xy^3-y^4
=
(x-y)(x+y)^3
=
2x(x+y)^3-(x+y)^4
$$

正是相应齐次方程的解。

**推论**　若非齐次项为 $f(ax+by)$，且 $\hat L(\hat D_x,\hat D_y)$ 是 $\hat D_x,\hat D_y$ 的齐（$n$）次式，则

$$
\hat D_x^rg(ax+by)
=
a^rg^{(r)}(ax+by),
\qquad
\hat D_y^sg(ax+by)
=
b^sg^{(s)}(ax+by).
$$

所以

$$
\hat L(\hat D_x,\hat D_y)g(ax+by)
=
L(a,b)g^{(n)}(ax+by).
$$

因此，当 $L(a,b)\ne0$ 时，就有

$$
\boxed{
\frac1{\hat L(\hat D_x,\hat D_y)}
g^{(n)}(ax+by)
=
\frac1{L(a,b)}
g(ax+by).
}
\tag{12.22}
$$

**例 12.6**　求解方程

$$
\frac{\partial^2v}{\partial x^2}
+
\frac{\partial^2v}{\partial y^2}
=
12(x+y).
$$

**解**　先求特解。将方程写成

$$
(\hat D_x^2+\hat D_y^2)v
=
12(x+y),
$$

显然特解为

$$
v_0
=
\frac{12}{\hat D_x^2+\hat D_y^2}(x+y)
=
\frac{12}{1^2+1^2}
\frac1{3!}(x+y)^3
=
(x+y)^3.
$$

容易求出相应齐次方程的通解，从而得出非齐次方程的通解

$$
v
=
(x+y)^3
+
\phi(x+iy)
+
\psi(x-iy).
$$

此推论也有失效的情形，这就是 (12.22) 式中 $L(a,b)=0$。这时不妨先考虑一个相关的然而特殊的一阶非齐次偏微分方程

$$
(\hat D_x-\alpha\hat D_y)u
=
\phi(x)\psi(y+\alpha x).
\tag{12.23}
$$

同样可设

$$
u(x,y)
=
f(x)\psi(y+\alpha x).
$$

代入方程 (12.23)，可以得到 $f(x)$ 满足的常微分方程 $f'(x)=\phi(x)$，由此即可解得 $f(x)$。或者直接将方程 (12.23) 中的非齐次项写成 $f'(x)\psi(y+\alpha x)$，就有

$$
\boxed{
\frac1{\hat D_x-\alpha\hat D_y}
f'(x)\psi(y+\alpha x)
=
f(x)\psi(y+\alpha x).
}
\tag{12.24}
$$

反复利用 (12.24) 式的结果，还可以进一步得到

$$
\boxed{
\frac1{(\hat D_x-\alpha\hat D_y)^k}
f^{(k)}(x)\psi(y+\alpha x)
=
f(x)\psi(y+\alpha x).
}
\tag{12.25}
$$

**例 12.7**　求解

$$
(\hat D_x^2-6\hat D_x\hat D_y+9\hat D_y^2)u
=
6x+2y,
$$

即

$$
(\hat D_x-3\hat D_y)^2u
=
6x+2y.
$$

**解**　相应齐次方程的通解（见例 12.3）为

$$
x\phi(y+3x)+\psi(y+3x),
$$

而非齐次方程的特解为

$$
u_0
=
\frac1{(\hat D_x-3\hat D_y)^2}(6x+2y)
=
\frac2{(\hat D_x-3\hat D_y)^2}(3x+y)
=
x^2(y+3x).
$$

因此，非齐次方程的通解为

$$
u
=
x^2(y+3x)
+
x\phi(y+3x)
+
\psi(y+3x).
$$

对于方程 (12.23) 右端的非齐次项为 $\phi(x)\psi(y+\alpha x)$ 的更一般的情形，可类似地讨论。

## *§12.4 特殊的变系数线性齐次偏微分方程

先讨论

$$
x^my^n
\frac{\partial^{m+n}u}{\partial x^m\partial y^n}
\tag{12.26}
$$

形式的项。令

$$
x=e^t,
\qquad
y=e^s,
\tag{12.27}
$$

即 $t=\ln x,\ s=\ln y$，则有

$$
\hat D_t
\equiv
\frac{\partial}{\partial t}
=
x\frac{\partial}{\partial x},
\qquad
\hat D_s
\equiv
\frac{\partial}{\partial s}
=
y\frac{\partial}{\partial y}.
\tag{12.28}
$$

于是，

$$
x^2\frac{\partial^2}{\partial x^2}
=
\hat D_t(\hat D_t-1),
\qquad
y^2\frac{\partial^2}{\partial y^2}
=
\hat D_s(\hat D_s-1),
$$

$$
x^3\frac{\partial^3}{\partial x^3}
=
\hat D_t(\hat D_t-1)(\hat D_t-2),
\qquad
y^3\frac{\partial^3}{\partial y^3}
=
\hat D_s(\hat D_s-1)(\hat D_s-2),
$$

$$
\vdots
\qquad\qquad\qquad\qquad
\vdots
$$

更一般地，

$$
x^my^n
\frac{\partial^{m+n}}{\partial x^m\partial y^n}
=
\hat D_t(\hat D_t-1)\cdots(\hat D_t-m+1)
\hat D_s(\hat D_s-1)\cdots(\hat D_s-n+1).
\tag{12.29}
$$

所以，对于 $\hat L(\hat D_x,\hat D_y)$ 均由 (12.26) 式形式的项组成时，通过变换 (12.27) 可以化为常系数的微分方程。下面通过一个具体的例子来说明这类方程的解法。

**例 12.8**　求方程

$$
x^2\frac{\partial^2u}{\partial x^2}
-
y^2\frac{\partial^2u}{\partial y^2}
+
x\frac{\partial u}{\partial x}
-
y\frac{\partial u}{\partial y}
=
0
$$

的通解。

**解**　作变换 (12.27)，则方程化为

$$
[
\hat D_t(\hat D_t-1)
-
\hat D_s(\hat D_s-1)
+
\hat D_t
-
\hat D_s
]u
=
0,
$$

即

$$
(\hat D_t^2-\hat D_s^2)u=0.
$$

所以，方程的通解为

$$
\begin{aligned}
u
&=
\phi_1(t+s)+\psi_1(t-s)\\
&=
\phi_1(\ln x+\ln y)
+
\psi_1(\ln x-\ln y)\\
&=
\phi_1(\ln xy)
+
\psi_1\left(\ln\frac xy\right)\\
&=
\phi(xy)
+
\psi\left(\frac xy\right).
\end{aligned}
$$

## *§12.5 波动方程的行波解

在例 12.1 中，曾经讨论过波动方程

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0
\tag{12.30}
$$

的通解，这里，把它改写成

$$
u(x,t)
=
f(x-at)
+
g(x+at),
\tag{12.31}
$$

其中 $f$ 和 $g$ 是任意二阶可微函数。这个解式表明，波动方程 (12.30) 的通解由两个波组成：$f(x-at)$ 代表沿 $x$ 轴向右传播的波，当 $t=0$ 时，波形为 $f(x)$，而后以恒定速率 $a$ 向右传播，而保持波形不变；$g(x+at)$ 则代表沿 $x$ 轴向左传播的波，当 $t=0$ 时，波形为 $g(x)$，而后也以同样的恒定速率 $a$ 向左传播，保持波形不变。单独的 $f(x-at)$ 和 $g(x+at)$ 也是方程 (12.30) 的解，它们独立传播，互不干扰，这正是因为波动方程是线性齐次方程，其解具有叠加性。

原则上说，函数 $f$ 和 $g$ 应该由定解条件确定。如果把问题简化为一维无界弦上的波的传播问题，那么，$f$ 和 $g$ 当然便完全由初始条件决定。

**例 12.9**　求解定解问题①

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
-\infty<x<\infty,\ t>0,
\tag{12.32a}
$$

$$
\left.u(x,t)\right|_{t=0}
=
\phi(x),
\qquad
-\infty<x<\infty,
\tag{12.32b}
$$

$$
\left.
\frac{\partial u}{\partial t}
\right|_{t=0}
=
\psi(x),
\qquad
-\infty<x<\infty.
\tag{12.32c}
$$

① 这个定解问题中明显缺少了边界条件。严格说来，这里的确应该明确写出无穷远条件

$$
\left.u(x,t)\right|_{x\to\pm\infty}\to0
\qquad
\text{或}
\qquad
\left.u(x,t)\right|_{x\to\pm\infty}\text{ 有界}.
$$

但就具体问题而言，这个条件可以由 $\phi(x)$ 和 $\psi(x)$ 的具体形式来得到保证。$\phi(x)$ 和 $\psi(x)$ 总是会局限在一个有限的范围内。当 $|x|$ 增大时，$\phi(x)$ 和 $\psi(x)$ 都会足够快地趋于 0。因此，从解 (12.34) 就可以看出，在有限的时间内，$u(x,t)$ 总还是在一个有限的范围内才不为 0。从概念上说，所谓无穷长的弦当然只是一个理想化的抽象。它恰恰就是表示：在我们所考察的时间和空间范围内，端点的影响可以忽略不计。


**解**　方程 (12.32a) 的通解已由 (12.31) 式给出，现在的问题便是如何根据初始条件 (12.32b) 和 (12.32c) 确定函数 $f$ 和 $g$。为此，将解式 (12.31) 代入初始条件，得

$$
f(x)+g(x)=\phi(x),
\tag{12.33a}
$$

$$
a[f'(x)-g'(x)]
=
-\psi(x).
\tag{12.33b}
$$

将 (12.33b) 积分，可以得到

$$
f(x)-g(x)
=
-\frac1a
\int_0^x\psi(\xi)\,d\xi
+
C,
$$

其中 $C$ 是积分常数。将这个结果和 (12.33a) 联立，即可求得

$$
f(x)
=
\frac12\phi(x)
-
\frac1{2a}
\int_0^x\psi(\xi)\,d\xi
+
\frac C2,
$$

$$
g(x)
=
\frac12\phi(x)
+
\frac1{2a}
\int_0^x\psi(\xi)\,d\xi
-
\frac C2.
$$

再代回解式 (12.31)，即得

$$
\boxed{
u(x,t)
=
\frac12[\phi(x-at)+\phi(x+at)]
+
\frac1{2a}
\int_{x-at}^{x+at}
\psi(\xi)\,d\xi.
}
\tag{12.34}
$$

这样，就求得了一维无界区间上波动方程定解问题 (12.32) 的解。它称为一维波动方程定解问题的**行波解**，或 d’Alembert 解。这个解具有清楚的物理意义：第一项表示由初位移 $\left.u(x,t)\right|_{t=0}=\phi(x)$ 激发的行波，$t>0$ 以后分成相等的两部分，独立地向左、右传播，速率为 $a$；第二项表示由初速度 $\left.\partial u/\partial t\right|_{t=0}=\psi(x)$ 激发的行波，在 $t>0$ 时刻，它左右对称地扩展到 $[x-at,x+at]$ 的范围，所以传播速率也是 $a$。

无界弦自由振动的这种传播特性，可以利用几何图形更加形象地表现出来。显然，定解问题 (12.32) 的定义域是 $x-t$ 平面上的上半平面。如果弦上一点 $x_0$ 在 $t=0$ 时（即对应于 $x-t$ 平面 $x$ 轴上的 $(x_0,0)$ 点）受到激发，则此后一定只波及区域

$$
x_0-at
\le
x
\le
x_0+at.
$$

这个区域就是弦上 $x_0$ 点的影响区域，因为区域内任意一点一定会受到初始激发的影响，而区域外的任何一点一定不会受到初始激发的影响。这里的直线 $x=x_0-at$ 和 $x=x_0+at$ 称为波动方程的（过 $x_0$ 点的）特征线。类似地，由特征线 $x=x_1-at$ 和 $x=x_2+at$ 围成的平面区域

$$
x_1-at
\le
x
\le
x_2+at,
\qquad
x_2>x_1,\ t\ge0
$$

就是（$x$ 轴上）区间 $[x_1,x_2]$ 的影响区域。

![Figure 12.1](...)

(a) 弦上 $x_0$ 点的影响区域  
(b) 区间 $[x_1,x_2]$ 的影响区域  
(c) $(x,t)$ 点的依赖区间

图 12.1　一维无界空间波动方程的影响区域和依赖区间

还可以从相反的角度提出问题：$x-t$ 平面上一点 $(x,t)$（即弦上 $x$ 点在 $t$ 时刻）的位移到底与 $x$ 轴上哪些点的初始激发有关？解式 (12.34) 告诉我们，对于 $x-t$ 平面上的任意一点 $(x,t)$，其位移仅依赖于弦上 $[x-at,x+at]$ 中的初始激发，而与区间外任何点的初始激发无关。因此 $x$ 轴上的区间 $[x-at,x+at]$ 就是 $(x,t)$ 点的依赖区间。反之，由 $x$ 轴和 $x=x_1+at$ 及 $x=x_2-at$ 所围成的三角形区域，就是 $x$ 轴上区间 $[x_1,x_2]$ 的决定区域：决定区域内任意一点的位移，完全由区间 $[x_1,x_2]$ 上的初始激发决定。

简要回顾求解定解问题 (12.32) 的过程不难发现，我们是照搬求解二阶线性常微分方程定解问题的思路：首先求出方程的通解，然后将通解代入定解条件，确定通解中的任意常数（对常微分方程）或任意函数（对偏微分方程），从而得到定解问题的唯一解。这种求解偏微分方程定解问题的方法称为行波法。行波法思路清晰，得到的解物理图像清楚，但只适用于无界空间波动方程的初值问题，而且基本上只能限于一维无界空间。对于一般的二维或三维空间，我们就难以写出波动方程的通解，自然也就无法用行波法求解。为了克服这种局限性，我们在下一章将介绍另外一种求解偏微分方程定解问题的方法。它将不仅适用于波动问题，也适用于热传导问题和稳定问题，而且不受空间维数的限制。

## *§12.6 波的耗散和色散

本节在一维无界空间波动方程行波解的基础上，简要地讨论一下波动过程中的一些基本物理特征。

波动方程 (12.30) 描写的是以恒定速率 $a$ 传播的非衰减波。这是因为我们假定了在波动过程中不存在耗散，表现为无界弦的总能量①

$$
\begin{aligned}
E(t)
&=
\frac12
\int_{-\infty}^{\infty}
\left[
\rho
\left(
\frac{\partial u}{\partial t}
\right)^2
+
T
\left(
\frac{\partial u}{\partial x}
\right)^2
\right]dx\\
&=
\frac12
\int_{-\infty}^{\infty}
\left[
\left(
\frac{\partial u}{\partial t}
\right)^2
+
a^2
\left(
\frac{\partial u}{\partial x}
\right)^2
\right]\rho\,dx
\end{aligned}
\tag{12.35}
$$

守恒，

$$
\begin{aligned}
\frac{dE(t)}{dt}
&=
\int_{-\infty}^{\infty}
\left(
\frac{\partial u}{\partial t}
\frac{\partial^2u}{\partial t^2}
+
a^2
\frac{\partial u}{\partial x}
\frac{\partial^2u}{\partial x\partial t}
\right)\rho\,dx\\
&=
a^2\rho
\left.
\frac{\partial u}{\partial x}
\frac{\partial u}{\partial t}
\right|_{-\infty}^{\infty}
+
\int_{-\infty}^{\infty}
\frac{\partial u}{\partial t}
\left(
\frac{\partial^2u}{\partial t^2}
-
a^2
\frac{\partial^2u}{\partial x^2}
\right)
\rho\,dx
=
0.
\end{aligned}
\tag{12.36}
$$

相应地，反映在方程上，它在时间反演 $t\to-t$ 下是不变的。我们还做了小振动的假设，因此，涉及振动位移 $u(x,t)$ 的非线性效应可以忽略。与此相适应的是，传播速率为常数，不依赖于频率和波长。所有这些简化假设，就导致了方程 (12.30) 的形式与介质的具体性质无关。介质的具体性质仅仅体现在传播速率 $a$ 的大小上。

波动方程 (12.30) 的解可以分解为向左、右独立传播的两个波。我们可以跟踪其中的一个波，而完全忽略另一个波的存在。例如，不妨只考察向右传播的波 $u(x,t)=f(x-at)$，它是一阶波动方程

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}
=
0
\tag{12.37}
$$

的解。在此基础上，以适当方式进一步包括进方程建立过程中所略去的一些高级项，就可以明白地表现出波的耗散和色散。

首先，在方程 (12.37) 中加上 $\partial^2u/\partial x^2$ 项，例如

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}
-
\alpha\frac{\partial^2u}{\partial x^2}
=
0,
\tag{12.38}
$$

① 关于能量表达式的推导，可参阅钟锡华、周岳明所著《大学物理通用教程·力学》（北京：北京大学出版社，2000 年）的第 9.5 节。


$a$ 和 $\alpha$ 仍为常数。如果仍然要寻找谐波形式的解①

$$
u(x,t)
=
\int_{-\infty}^{\infty}
A(k)e^{i(kx-\omega t)}\,dk,
\tag{12.39}
$$

代入方程 (12.38)，假定可以和积分交换次序，那么就可以得到波数 $k$ 和“角频率” $\omega$ 所必须满足的关系式

$$
\omega
=
ka-i\alpha k^2,
$$

因此，方程 (12.38) 的解就是

$$
u(x,t)
=
\int_{-\infty}^{\infty}
A(k)e^{-\alpha k^2t}
e^{ik(x-at)}\,dk.
$$

这说明方程 (12.38) 所描写的波动过程仍是以恒定速率 $a$ 传播的波，然而振幅却随时间而指数地衰减（设常数 $\alpha>0$）。由于振幅的衰减因子与波数 $k$（或波长 $2\pi/k$）有关，不同波数的分量衰减速度不同，因此波的覆盖区间大小保持不变，但波形却将随着时间而改变。

如果在方程 (12.37) 中加上 $\partial^3u/\partial x^3$ 项，我们就会得到

$$
\frac{\partial u}{\partial t}
+
a\frac{\partial u}{\partial x}
+
\beta\frac{\partial^3u}{\partial x^3}
=
0,
\tag{12.40}
$$

$a$ 和 $\beta$ 仍为常数。倘若仍然要寻找 (12.39) 形式的解，则波数 $k$ 和角频率 $\omega$ 必须满足的关系式就是

$$
\omega
=
k(a-\beta k^2),
$$

所以

$$
kx-\omega t
=
k[x-(a-\beta k^2)t],
$$

方程 (12.40) 的解就是

$$
u(x,t)
=
\int_{-\infty}^{\infty}
A(k)e^{ik[x-(a-\beta k^2)t]}\,dk.
\tag{12.41}
$$

这样，波的传播速率（准确说，是相位的传播速率，即相速）

$$
v_p
=
\frac{\omega}{k}
=
a-\beta k^2
\tag{12.42}
$$

就是 $k$ 的函数。这说明，波数不同，传播速率不同。甚至对于 $k^2<a/\beta$ 和 $k^2>a/\beta$ 的分量，传播方向也不同，因此，随着时间的增大，波的覆盖区间将越来越大；在空间一点上，波的组分（不同波数的分量所占比例）也将随着时间而变化。这个现象，就称为波的色散。

还可以定义另一种传播速率

$$
v_g
=
\frac{d\omega}{dk}
=
a-3\beta k^2,
\tag{12.43}
$$

称为群速（见图 12.2）。它是波包的传播速率，因而也就是能量的传播速率。群速度和相速度不相等，是色散波的又一特点。

![Figure 12.2](...)

图 12.2　群速 $v_g$ 与相速 $v_p$

① 这是复数形式的解。真正的位移是它的实部或虚部。$k$ 必须是实数，因为在 $t=0$ 时解也必须是 $x$ 的振荡函数。

上面定性地讨论了波的耗散和色散。共同特点是方程仍保持为齐次，因而还具有解的可叠加性。如果在一阶波动方程中引进其他形式的线性修正项，也有可能产生耗散和色散，甚至会同时出现耗散和色散。回到原来的二阶波动方程，也可类似地讨论。

最后简要提一下非线性效应。若在一定条件下，在波动方程中出现非线性项，例如

$$
\frac{\partial u}{\partial t}
+
a(1+\gamma u)
\frac{\partial u}{\partial x}
=
0,
\tag{12.44}
$$

这可以理解为传播速率不只与介质的物理性质有关，还与位移的大小有关。可以直接验证，这个方程的解能够写成隐函数的形式：

$$
u(x,t)
=
f[x-a(1+\gamma u)t],
\tag{12.45}
$$

其中 $f$ 仍是任意可微函数，它就是 $u(x,t)$ 在 $t=0$ 时的形状 $f(x)$，

$$
\left.u(x,t)\right|_{t=0}
=
f(x).
$$

这里需要特别强调，由于在方程 (12.44) 中出现了非线性项

$$
u\frac{\partial u}{\partial x}
=
\frac12
\frac{\partial(u^2)}{\partial x},
$$

方程 (12.44) 的解不再具有可叠加性：如果 $u(x,t)$ 是方程 (12.44) 的解，$Au(x,t)$（$A$ 为常数，且 $\ne1$）也不再是该方程的解。同样，$u_1(x,t)$ 和 $u_2(x,t)$ 都是方程 (12.44) 的解，它们的线性组合一般也不会是 (12.44) 的解。这样必然会产生两个后果：一方面是方程更加难以求解，特别是难以求出全部解；另一方面，方程的解必然会表现出新的特点（相应地，物理上就表现为新的规律性）。研究各种非线性方程的求解以及解的特性，是非线性科学的课题。这里从略。

## *§12.7 热传导方程的定性讨论

对于热传导方程，最明显的特点便是存在耗散。热传导方程属于抛物型方程，在方程中含有未知函数对时间变量的一阶偏导数和对空间变量的二阶偏导数，因此，方程不具有时间反演不变性。换句话说，热传导过程是不可逆的。为了定性地了解热传导方程的特点，不妨讨论一下无穷长的一维介质上的热传导问题：

$$
\frac{\partial u}{\partial t}
-
\kappa
\frac{\partial^2u}{\partial x^2}
=
0,
\tag{12.46a}
$$

$$
\left.u\right|_{t=0}
=
\phi(x).
\tag{12.46b}
$$

与方程 (12.38) 相比较就可以看出，热传导方程 (12.46a) 的解可以表示成

$$
u(x,t)
=
\int_{-\infty}^{\infty}
A(k)e^{-\kappa k^2t}e^{ikx}\,dk,
\tag{12.47}
$$

其中的叠加系数 $A(k)$ 由初值函数 $\phi(x)$ 决定，

$$
\phi(x)
=
\int_{-\infty}^{\infty}
A(k)e^{ikx}\,dk,
\tag{12.48}
$$

$$
A(k)
=
\frac1{2\pi}
\int_{-\infty}^{\infty}
\phi(x)e^{-ikx}\,dx.
\tag{12.49}
$$

这里就明白无误地表现出了温度 $u(x,t)$ 随时间 $t$ 的衰减。

将 (12.49) 式代入 (12.47) 式，并交换积分次序，就得到

$$
\begin{aligned}
u(x,t)
&=
\frac1{2\pi}
\int_{-\infty}^{\infty}
e^{-\kappa k^2t}\,dk
\int_{-\infty}^{\infty}
e^{ikx}\phi(x')e^{-ikx'}\,dx'\\
&=
\sqrt{\frac1{4\pi\kappa t}}
\int_{-\infty}^{\infty}
\phi(x')
\exp
\left[
-\frac{(x-x')^2}{4\kappa t}
\right]dx'.
\end{aligned}
\tag{12.50}
$$

最后一步用到了 (4.29) 式的结果。(12.50) 式表明：$t=0$ 时在一维介质上任意一点的初始温度，都将立即传播到整个（无穷长）介质上的每一点。换言之，热传导方程隐含着传热速率为无穷。就介质上的任意一点而言，它的影响区域都是 $(x,t)$ 平面上的整个上半平面。反过来说，对于 $x-t$ 平面上的任意一点 $(x,t)$，它的依赖区间都是整个 $x$ 轴。

仿照 §10.3 的做法，还可定义热传导方程的 Green 函数 $g(x,t;x',t')$，它是定解问题

$$
\frac{\partial g}{\partial t}
-
\kappa
\frac{\partial^2g}{\partial x^2}
=
\delta(x-x')\delta(t-t'),
\qquad
t'>0,
\tag{12.51a}
$$

$$
\left.g(x,t;x',t')\right|_{t=0}
=
0
\tag{12.51b}
$$

的解。这个定解问题描述了瞬时（只存在于 $t'$ 时刻）点（集中在空间 $x'$ 一点）热源所产生的温度场。$t=0$ 时介质的温度为 0。这样，在热源出现之前（即 $t<t'$ 时），介质的温度显然仍一定维持为 0。当 $t>t'$ 时，可以求出（见 §19.6 的 (19.90) 式），

$$
g(x,t;x',t')
=
\frac1{2\sqrt{\kappa\pi(t-t')}}
\exp
\left[
-\frac{(x-x')^2}{4\kappa(t-t')}
\right].
\tag{12.52}
$$

这也明白无误地表现出热传导方程的无穷传热速率，因为只要 $t>t'$，则对于介质上的任何一点 $x$，$g(x,t;x',t')$ 总不为 0。这就是说，只要一旦热源在介质上某点出现，则不论距离多远处，总立即感受到它的影响。

无穷的传热速度，是热传导方程的基本特性之一。这个特性，对于（有界或无界的）二维或三维介质仍然成立。因此，如果 $t=0$ 时介质内某点温度达到极值，那么在 $t>0$ 时就会迅速弥散而消失。同样，如果介质的边界温度和初始温度都不超过某一常数 $M$，且介质内部无热源，则介质内部任意一点在任一时刻的温度也不会超过 $M$。以一维热传导方程

$$
\frac{\partial u}{\partial t}
-
\kappa
\frac{\partial^2u}{\partial x^2}
=
0
$$

为例，可以证明，若 $u(x,t)$ 在 $x-t$ 平面上的矩形区域 $a\le x\le b,\ 0\le t\le T$ 上连续，则它一定在矩形的两条侧边（$x=a$ 和 $x=b,\ 0\le t\le T$）及底边（$t=0,\ a\le x\le b$）达到最大值和最小值。这是热传导方程和波动方程的重要区别。波动方程具有有限的传播速度，在波的传播过程中位移（或其他相关的物理量）可以在区域内部达到最大和最小，甚至可以不断反复出现。

这种无穷的传热速率，当然只是一个近似。因此，只有在热源出现后的足够长的时间之后，热传导方程才可能足够好地描写实际的热传导过程。但由于这个方程形式简单，便于求解，只要物理问题的近似程度许可，我们总还是采用这种简化的热传导方程。

之所以出现无穷大的传热速度，是由于建立热传导方程时过分简化了热传导过程的微观机理。连续介质中出现的任何不均匀性，例如引起热传导过程的热运动的不均匀性，引起扩散过程的物质密度的不均匀性，总会由于微观粒子（分子、原子……）的运动（碰撞）而以有限的速度逐渐传播开来。这样，热传导方程便应该修改为

$$
\frac{\partial^2u}{\partial x^2}
=
\frac1\kappa
\frac{\partial u}{\partial t}
+
\frac1{a^2}
\frac{\partial^2u}{\partial t^2}.
\tag{12.53}
$$

方程的右端出现了两项。如果第一项占主导作用，第二项只是一个修正，方程的解仍然应当表现出热传导的基本特征①，但不容忽视，第二项的存在，的确可以导致传热速度为有限值。这从定解问题

① 如果第二项占主导作用，方程的解当然就表现出波动过程的基本特征，而第一项就可以看成耗散项。

$$
\frac{\partial^2\mathcal G}{\partial x^2}
-
\frac1\kappa
\frac{\partial\mathcal G}{\partial t}
-
\frac1{a^2}
\frac{\partial^2\mathcal G}{\partial t^2}
=
-\delta(x-x')\delta(t-t'),
\qquad
(t,t'>0,\ -\infty<x,x'<\infty),
\tag{12.54a}
$$

$$
\left.\mathcal G\right|_{t=0}
=
0,
\qquad
\left.
\frac{\partial\mathcal G}{\partial t}
\right|_{t=0}
=
0
\qquad
(-\infty<x,x'<\infty)
\tag{12.54b}
$$

的解

$$
\boxed{
\mathcal G(x,t)
=
\frac a2
J_0
\left(
\frac a{2\kappa}
\sqrt{
(x-x')^2-a^2(t-t')^2
}
\right)
\exp
\left[
-\frac{a^2}{2\kappa}(t-t')
\right]
\eta
\left(
t-t'
-
\frac{|x-x'|}{a}
\right)
}
\tag{12.55}
$$

可以看出，正如预料的那样，上面给出的 (12.52) 式，正是 $\mathcal G(x,t)$ 在 $a\to\infty$ 的极限下的结果。利用零阶 Bessel 函数 $J_0(z)$ 的渐近展开（见 (9.67) 或 (16.22) 式）就能够验证明这一点。

## *§12.8 Laplace 方程的定性讨论

在本书的复变函数部分已经指出，解析函数

$$
f(z)
=
f(x+iy)
=
u(x,y)+iv(x,y)
$$

的实部 $u(x,y)$（或虚部 $v(x,y)$），一定是二维 Laplace 方程

$$
\nabla^2u(x,y)
=
0
\tag{12.56}
$$

的解，称为二维调和函数。另外，还可以证明（见 §3.4 的 (3.18) 式）：在解析函数的解析区域内，任意一点 $a$ 的函数值一定等于以该点为圆心的圆周上各点函数值的平均值

$$
f(a)
=
\frac1{2\pi}
\oint_{|z-a|=R}
f(z)\,d\theta.
\tag{12.57}
$$

进一步比较等式两端的实部（或虚部）就可以看出，均值定理对解析函数的实部（或虚部）单独也成立：

$$
u(x_0,y_0)
=
\frac1{2\pi R}
\oint_C
u(x,y)\,dl,
\qquad
C:\ (x-x_0)^2+(y-y_0)^2=R^2.
\tag{12.58}
$$

由此还可推断出，只要 $u(x,y)\ne$ 常数，$u(x,y)$ 的最大值和最小值一定只能出现在圆周上。用反证法就可以证明这个结论（称为极值原理）。它可以看成是最大模原理在 $u(x,y)$ 上的具体表现。这样，均值定理和极值原理就是二维调和函数的最基本特性。

也可以建立三维调和函数的概念。如果在区域 $V$ 内函数的二阶偏导数存在，且满足三维 Laplace 方程，则称该函数为 $V$ 内的三维调和函数。可以证明，均值定理对于三维调和函数仍然成立，因此极值原理仍然成立。

作为三维调和函数的具体实例，不妨讨论一下（直角坐标系中）三维 Laplace 方程的多项式解（二维调和函数是它们的特殊情形）。因为 Laplace 方程中，对各自变量的偏导数阶数相同，所以，这样的多项式独立解可取为自变量的齐次函数。显然，常数（0 次齐次式）和一次函数（一次齐次式）$x,y,z$ 就是这样的独立解。对于二次函数，独立解有五个，可以取为

$$
2z^2-x^2-y^2,
\qquad
xz,
\qquad
yz,
\qquad
xy,
\qquad
x^2-y^2.
$$

三次函数的独立解有七个，可取为

$$
(4z^2-x^2-y^2)x,
\qquad
(4z^2-x^2-y^2)y,
\qquad
(2z^2-3x^2-3y^2)z,
$$

$$
xyz,
\qquad
(x^2-y^2)z,
\qquad
(y^2-3x^2)y,
\qquad
(x^2-3y^2)x.
$$

更一般地，$l$ 次函数的独立解有 $2l+1$ 个。它们的形式将在 §16.9 中给出，并且会看到，为什么 $l$ 次函数的独立解有 $2l+1$ 个。

## 习题

1. 求下列线性齐次偏微分方程的通解：

   (1)
   $$
   \frac{\partial^2u}{\partial x^2}
   -
   2\frac{\partial^2u}{\partial x\partial y}
   -
   3\frac{\partial^2u}{\partial y^2}
   =
   0;
   $$

   (2)
   $$
   \frac{\partial^2u}{\partial x^2}
   -
   2\frac{\partial^2u}{\partial x\partial y}
   +
   2\frac{\partial^2u}{\partial y^2}
   =
   0;
   $$

   (3)
   $$
   \frac{\partial^2u}{\partial x^2}
   -
   \frac{\partial^2u}{\partial x\partial y}
   =
   0;
   $$

   (4)
   $$
   \frac{\partial^2u}{\partial t^2}
   =
   \frac{c^2}{r^2}
   \frac{\partial}{\partial r}
   \left(
   r^2\frac{\partial u}{\partial r}
   \right),
   \qquad
   c\ne0;
   $$

   (5)
   $$
   (a^2-b^2)
   \frac{\partial^2u}{\partial x^2}
   +
   2a
   \frac{\partial^2u}{\partial x\partial t}
   +
   \frac{\partial^2u}{\partial t^2}
   =
   0,
   \qquad
   b\ne0;
   $$

   (6)
   $$
   \frac{\partial^4u}{\partial x^4}
   -
   \frac{\partial^4u}{\partial y^4}
   =
   0.
   $$

2. 求下列线性非齐次偏微分方程的通解：

   (1)
   $$
   \frac{\partial^2u}{\partial x^2}
   +
   \frac{\partial^2u}{\partial y^2}
   =
   x^2+xy;
   $$

   (2)
   $$
   \frac{\partial^2u}{\partial x^2}
   -
   \frac{\partial^2u}{\partial y^2}
   =
   xy-x;
   $$

   (3)
   $$
   \frac{\partial^2u}{\partial x^2}
   -
   2\frac{\partial^2u}{\partial x\partial y}
   +
   \frac{\partial^2u}{\partial y^2}
   =
   x^2+y.
   $$

3. 求解偏微分方程：

   (1)
   $$
   x^2\frac{\partial^2u}{\partial x^2}
   -
   2xy\frac{\partial^2u}{\partial x\partial y}
   +
   y^2\frac{\partial^2u}{\partial y^2}
   +
   x\frac{\partial u}{\partial x}
   +
   y\frac{\partial u}{\partial y}
   =
   0;
   $$

   (2)
   $$
   \frac{\partial^2u}{\partial x^2}
   -
   \frac{\partial^2u}{\partial y^2}
   =
   (x^2-y^2)\sin xy.
   $$

4. 用行波法求解一维半无界弦的波动问题

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
\qquad
x>0,\ t>0,
$$

$$
\left.u\right|_{x=0}
=
0,
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
x\ge0.
$$

5. 求解偏微分方程初值问题

$$
\frac{\partial}{\partial x}
\left[
\left(
1-\frac xl
\right)^2
\frac{\partial u}{\partial x}
\right]
-
\frac1{a^2}
\left(
1-\frac xl
\right)^2
\frac{\partial^2u}{\partial t^2}
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
