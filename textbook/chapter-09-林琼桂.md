# 第九章 正交曲线坐标系中的分离变量*

## 目录

§1 正交曲线坐标系中的微分算符 .................................................... 185  
§1.1 正交曲线坐标系 .............................................................. 185  
§1.2 正交曲线坐标系中的微分算符 .................................................. 187  

§2 球坐标系中的分离变量 .......................................................... 189  
§2.1 *Helmholtz* 方程 ............................................................. 189  
§2.2 *Helmholtz* 方程在球坐标系中的分离变量 ....................................... 189  

§3 柱坐标系中的分离变量 .......................................................... 191  

补充习题 .......................................................................... 194

本章开始把分离变量法推广到比较实际的三维问题。前已指出，分离变量时，应该根据边界的形状采用适当的坐标系。本章的目的就是研究如何在球坐标系和柱坐标系中对各类方程分离变量。我们所知道的几类方程都包含有 Laplace 算符，所以首先需要研究 Laplace 算符在曲线坐标系，尤其是球坐标系和柱坐标系中的形式。

## §1 正交曲线坐标系中的微分算符

### §1.1 正交曲线坐标系

直角坐标系 $(x,y,z)$，球坐标系 $(r,\theta,\phi)$，柱坐标系 $(\rho,\phi,z)$ 都是正交曲线坐标系。它们的三族坐标线处处相互正交。现考虑一般曲线坐标系，其坐标记作 $(q_1,q_2,q_3)$，它们与直角坐标之间的变换关系为

$$
x=x(q_1,q_2,q_3),\qquad y=y(q_1,q_2,q_3),\qquad z=z(q_1,q_2,q_3). \tag{1a}
$$

反之，$(q_1,q_2,q_3)$ 也可以表示为 $(x,y,z)$ 的函数。

**注**　一般来说，我们要求 Jacobi 行列式

$$
J=\frac{\partial(x,y,z)}{\partial(q_1,q_2,q_3)}=
\begin{vmatrix}
\partial x/\partial q_1 & \partial x/\partial q_2 & \partial x/\partial q_3\\
\partial y/\partial q_1 & \partial y/\partial q_2 & \partial y/\partial q_3\\
\partial z/\partial q_1 & \partial z/\partial q_2 & \partial z/\partial q_3
\end{vmatrix}\ne0.
$$

在 $J=0$ 处，给定 $(x,y,z)$，可能无法唯一确定 $(q_1,q_2,q_3)$，这些点就是坐标系的奇点。比如球坐标系，$J=r^2\sin\theta$，在 $r=0$ 处或 $\theta=0、\pi$ 处 $J=0$。此时给定 $(x,y,z)$，无法唯一确定 $(r,\theta,\phi)$。事实上，在 $\theta=0、\pi$ 处，$\phi$ 没有定义，而在 $r=0$ 处，$\theta、\phi$ 均没有定义。所以整个 $z$ 轴都是球坐标系的奇点。不过，这些点构成的集合测度（可以粗略地理解为体积）为零。这一般是允许的。

我们将式 (1a) 写成简洁的形式

$$
\boldsymbol r=\boldsymbol r(q_1,q_2,q_3), \tag{1b}
$$

并在下面的推导中使用类似的记号，这可以减少书写上的麻烦，也可以使有关结果显得更清晰。如果读者一时看不清楚式子的含义，可以在草稿纸上写出详细的分量形式加以对照。熟悉类似的记号，对于进一步学习后续的课程是非常有利的。

现在的问题是，什么样的曲线坐标系才算是正交的？

在直角坐标系中，相邻两点 $\boldsymbol r$ 和 $\boldsymbol r+d\boldsymbol r$ 之间的距离记作 $ds$，则

$$
(ds)^2=d\boldsymbol r\cdot d\boldsymbol r=(dx)^2+(dy)^2+(dz)^2. \tag{2}
$$

这一距离的定义是读者所熟悉的。但值得指出，这样的定义并不是天经地义或理所当然的。实际上我们还可以定义其它形式的距离。距离的定义不同，空间的性质可能也就不同。上述定义表明我们研究的是平直的 Euclid 空间。

现在看看上式在曲线坐标系中的形式。由于 $d\boldsymbol r=\sum_{i=1}^{3}(\partial\boldsymbol r/\partial q_i)dq_i$，故

$$
(ds)^2=d\boldsymbol r\cdot d\boldsymbol r
=\sum_{i,j=1}^{3}\frac{\partial\boldsymbol r}{\partial q_i}\cdot\frac{\partial\boldsymbol r}{\partial q_j}dq_i dq_j. \tag{3}
$$

如果写出详细的分量形式，则上式包含了 18 项。现在，我们可以给出正交曲线坐标系的定义：如果

$$
\frac{\partial\boldsymbol r}{\partial q_i}\cdot\frac{\partial\boldsymbol r}{\partial q_j}=0,\qquad i\ne j \tag{4}
$$

则曲线坐标系 $(q_1,q_2,q_3)$ 称为正交的。这时有

$$
(ds)^2=\sum_{i=1}^{3}\left(\frac{\partial\boldsymbol r}{\partial q_i}\right)^2(dq_i)^2
\equiv\sum_{i=1}^{3}(h_i dq_i)^2, \tag{5}
$$

其中 $h_i=|\partial\boldsymbol r/\partial q_i|$ 称为度规系数。这与直角坐标系中的形式 (2) 相似，只是 $dq_i$ 前面多了度规系数 $h_i$。所以，正交的关键就是 $(ds)^2$ 的表达式中不包含 $dq_1dq_2$ 这样的交叉项。

现在推导正交曲线坐标系 $(q_1,q_2,q_3)$ 中的单位矢量 $(\boldsymbol e_1,\boldsymbol e_2,\boldsymbol e_3)$。$\boldsymbol e_1$ 的方向就是 $q_1$ 坐标线的切线方向。沿着坐标线 $q_1$，有 $d\boldsymbol r=(\partial\boldsymbol r/\partial q_1)dq_1$，相应地 $ds=|d\boldsymbol r|=h_1dq_1$，故 $\boldsymbol e_1=d\boldsymbol r/ds=h_1^{-1}(\partial\boldsymbol r/\partial q_1)$。对 $\boldsymbol e_2$ 和 $\boldsymbol e_3$ 有类似结果。总结起来，就是

$$
\boldsymbol e_i=\frac{1}{h_i}\frac{\partial\boldsymbol r}{\partial q_i},\qquad i=1,2,3. \tag{6}
$$

由式 (4) 和度规系数的定义，易得

$$
\boldsymbol e_i\cdot\boldsymbol e_j=\delta_{ij}. \tag{7}
$$

可见三族坐标线确实是处处正交的。下面，我们用上述一般定义重新考察熟知的球坐标系和柱坐标系。

**例 1**　柱坐标系 $(\rho,\phi,z)$。由于 $\boldsymbol r=\rho\cos\phi\,\boldsymbol e_x+\rho\sin\phi\,\boldsymbol e_y+z\boldsymbol e_z$，故

$$
\frac{\partial\boldsymbol r}{\partial\rho}=\cos\phi\,\boldsymbol e_x+\sin\phi\,\boldsymbol e_y, \tag{8a}
$$

$$
\frac{\partial\boldsymbol r}{\partial\phi}=-\rho\sin\phi\,\boldsymbol e_x+\rho\cos\phi\,\boldsymbol e_y, \tag{8b}
$$

$$
\frac{\partial\boldsymbol r}{\partial z}=\boldsymbol e_z. \tag{8c}
$$

由此可以立得正交性，且有

$$
h_\rho=1,\qquad h_\phi=\rho,\qquad h_z=1. \tag{9}
$$

这些度规系数的几何意义是明显的，比如沿着 $\phi$ 坐标线由 $\phi$ 到 $\phi+d\phi$ 的距离（弧长）不是 $d\phi$，而是 $h_\phi d\phi=\rho d\phi$，这正是我们所熟知的。易得柱坐标系的正一单位矢量为

$$
\boldsymbol e_\rho=\cos\phi\,\boldsymbol e_x+\sin\phi\,\boldsymbol e_y, \tag{10a}
$$

$$
\boldsymbol e_\phi=-\sin\phi\,\boldsymbol e_x+\cos\phi\,\boldsymbol e_y, \tag{10b}
$$

$$
\boldsymbol e_z=\boldsymbol e_z. \tag{10c}
$$

**例 2**　球坐标系 $(r,\theta,\phi)$。由于 $\boldsymbol r=r\sin\theta\cos\phi\,\boldsymbol e_x+r\sin\theta\sin\phi\,\boldsymbol e_y+r\cos\theta\,\boldsymbol e_z$，故

$$
\frac{\partial\boldsymbol r}{\partial r}=\sin\theta\cos\phi\,\boldsymbol e_x+\sin\theta\sin\phi\,\boldsymbol e_y+\cos\theta\,\boldsymbol e_z, \tag{11a}
$$

$$
\frac{\partial\boldsymbol r}{\partial\theta}=r(\cos\theta\cos\phi\,\boldsymbol e_x+\cos\theta\sin\phi\,\boldsymbol e_y-\sin\theta\,\boldsymbol e_z), \tag{11b}
$$

$$
\frac{\partial\boldsymbol r}{\partial\phi}=r\sin\theta(-\sin\phi\,\boldsymbol e_x+\cos\phi\,\boldsymbol e_y). \tag{11c}
$$

由此可以立得正交性，且有

$$
h_r=1,\qquad h_\theta=r,\qquad h_\phi=r\sin\theta. \tag{12}
$$

这些度规系数的几何意义也是明显的。易得球坐标系的正一单位矢量为

$$
\boldsymbol e_r=\sin\theta\cos\phi\,\boldsymbol e_x+\sin\theta\sin\phi\,\boldsymbol e_y+\cos\theta\,\boldsymbol e_z, \tag{13a}
$$

$$
\boldsymbol e_\theta=\cos\theta\cos\phi\,\boldsymbol e_x+\cos\theta\sin\phi\,\boldsymbol e_y-\sin\theta\,\boldsymbol e_z, \tag{13b}
$$

$$
\boldsymbol e_\phi=-\sin\phi\,\boldsymbol e_x+\cos\phi\,\boldsymbol e_y. \tag{13c}
$$

### §1.2 正交曲线坐标系中的微分算符

微分算符主要有梯度、散度、旋度和 Laplace 算符，后者是梯度和散度的结合。首先考虑标量函数 $u(\boldsymbol r)$ 的梯度 $\nabla u(\boldsymbol r)$。

由于 $\nabla u$ 是矢量，它一定可以展开为 $\nabla u=\sum_{i=1}^{3}f_i\boldsymbol e_i$，注意其中 $f_i$ 是 $\boldsymbol r$ 的函数，$\boldsymbol e_i$ 的方向也是随着 $\boldsymbol r$ 变化的，这与直角坐标系的单位矢量不同。由正交关系 (7)，易得 $f_i=\boldsymbol e_i\cdot\nabla u$，另一方面，有

$$
\boldsymbol e_i\cdot\nabla u
=\frac{\partial u}{\partial s_i}
=\lim_{\Delta s_i\to0}\frac{\Delta u}{\Delta s_i}
=\lim_{\Delta q_i\to0}\frac{1}{h_i}\frac{\Delta u}{\Delta q_i}
=\frac{1}{h_i}\frac{\partial u}{\partial q_i},
$$

于是得到 $f_i=h_i^{-1}\partial u/\partial q_i$，从而

$$
\nabla u=\sum_{i=1}^{3}\frac{1}{h_i}\frac{\partial u}{\partial q_i}\boldsymbol e_i
=\frac{1}{h_1}\frac{\partial u}{\partial q_1}\boldsymbol e_1
+\frac{1}{h_2}\frac{\partial u}{\partial q_2}\boldsymbol e_2
+\frac{1}{h_3}\frac{\partial u}{\partial q_3}\boldsymbol e_3. \tag{14}
$$

与直角坐标系的相应表达式

$$
\nabla u=\frac{\partial u}{\partial x}\boldsymbol e_x+\frac{\partial u}{\partial y}\boldsymbol e_y+\frac{\partial u}{\partial z}\boldsymbol e_z
$$

比较，容易看出，不同之处在于曲线坐标系的结果中出现了度规系数。

由一般结果和度规系数 (9)，易得柱坐标系中的梯度算符为

$$
\nabla u=\frac{\partial u}{\partial\rho}\boldsymbol e_\rho+\frac{1}{\rho}\frac{\partial u}{\partial\phi}\boldsymbol e_\phi+\frac{\partial u}{\partial z}\boldsymbol e_z. \tag{15}
$$

又由度规系数 (12)，易得球坐标系中的梯度算符为

$$
\nabla u=\frac{\partial u}{\partial r}\boldsymbol e_r+\frac{1}{r}\frac{\partial u}{\partial\theta}\boldsymbol e_\theta+\frac{1}{r\sin\theta}\frac{\partial u}{\partial\phi}\boldsymbol e_\phi. \tag{16}
$$

本书主要用到 Laplace 算符，下面就给出两种常用的曲线坐标系、即柱坐标系和球坐标系中的 Laplace 算符形式。由于直角坐标系中的 Laplace 算符为

$$
\nabla^2u=\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}+\frac{\partial^2u}{\partial z^2},
$$

根据前面研究梯度算符所获得的经验，我们容易猜想，柱坐标系中的 Laplace 算符可能是

$$
\nabla^2u=\frac{\partial^2u}{\partial\rho^2}+\frac{1}{\rho^2}\frac{\partial^2u}{\partial\phi^2}+\frac{\partial^2u}{\partial z^2},
$$

而球坐标系中的 Laplace 算符可能是

$$
\nabla^2u=\frac{\partial^2u}{\partial r^2}+\frac{1}{r^2}\frac{\partial^2u}{\partial\theta^2}+\frac{1}{r^2\sin^2\theta}\frac{\partial^2u}{\partial\phi^2}.
$$

但这是不正确的。柱坐标系中的正确结果应该是

$$
\nabla^2u=\frac{1}{\rho}\frac{\partial}{\partial\rho}\left(\rho\frac{\partial u}{\partial\rho}\right)+\frac{1}{\rho^2}\frac{\partial^2u}{\partial\phi^2}+\frac{\partial^2u}{\partial z^2}, \tag{17}
$$

而球坐标系中的正确结果应该是

$$
\nabla^2u=\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial u}{\partial r}\right)
+\frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial u}{\partial\theta}\right)
+\frac{1}{r^2\sin^2\theta}\frac{\partial^2u}{\partial\phi^2}. \tag{18}
$$

如果想象括号中的内容与括号外的可以相消，就得到前面猜想的结果。由于微分算符的存在，相消的操作当然是不允许的。不过，这样的比较有助于我们把握这两个复杂的公式。

利用直角坐标系中的 Laplace 算符表达式和坐标变换，就可以得到以上两个公式。这是最简单的思路，但计算却比较繁琐，尤其是球坐标系的情况。下面给出的方法比较简单，而且适用于其它曲线坐标系，有兴趣的读者可以参考。

首先推导曲线坐标系中的散度表达式。矢量场 $\boldsymbol A(\boldsymbol r)$ 的散度定义为

$$
\nabla\cdot\boldsymbol A=\lim_{\Delta V\to0}\frac{\displaystyle\oint_{\partial(\Delta V)}\boldsymbol A\cdot d\boldsymbol S}{\Delta V}, \tag{19}
$$

其中 $\Delta V$ 是包含点 $\boldsymbol r$（欲求其散度值的一点，在以上定义中可以看作固定，但该点的选取是任意的）的区域，而 $\partial(\Delta V)$ 是 $\Delta V$ 的边界面，$d\boldsymbol S$ 是其面积元。今取 $\Delta V$ 为由坐标面 $q_1、q_1+\Delta q_1、q_2、q_2+\Delta q_2、q_3、q_3+\Delta q_3$ 所围成的六面体，计算上式分子中的积分，即矢量场在 $\Delta V$ 的边界面的通量。需要注意的是，六面体在 $q_1$ 处的面积是 $h_2h_3\Delta q_2\Delta q_3$，而不是 $\Delta q_2\Delta q_3$。于是

$$
\begin{aligned}
\oint_{\partial(\Delta V)}\boldsymbol A\cdot d\boldsymbol S
&=[(A_1h_2h_3)_{q_1+\Delta q_1}-(A_1h_2h_3)_{q_1}]\Delta q_2\Delta q_3+\cdots\\
&=\left.\frac{\partial(A_1h_2h_3)}{\partial q_1}\right|_{\boldsymbol r_0}\Delta q_1\Delta q_2\Delta q_3+\cdots.
\end{aligned}
$$

其中 $\boldsymbol r_0$ 是 $\Delta V$ 内一点，上式只写出了 $q_1$ 处和 $q_1+\Delta q_1$ 处两个面的贡献，“$\cdots$”表示其它四个面的贡献，容易写出其相应的表达式。由于 $\Delta V=h_1h_2h_3\Delta q_1\Delta q_2\Delta q_3$，而当 $\Delta q_1\to0、\Delta q_2\to0、\Delta q_3\to0$ 时，$\boldsymbol r_0\to\boldsymbol r$，故得

$$
\nabla\cdot\boldsymbol A=\frac{1}{h_1h_2h_3}\left[
\frac{\partial}{\partial q_1}(A_1h_2h_3)+
\frac{\partial}{\partial q_2}(A_2h_3h_1)+
\frac{\partial}{\partial q_3}(A_3h_1h_2)
\right]. \tag{20}
$$

这就是一般曲线坐标系中的散度表达式。

由于 $\nabla^2u=\nabla\cdot(\nabla u)$，结合式 (14) 和式 (20)，就得到一般曲线坐标系中的 Laplace 算符表达式：

$$
\nabla^2u=\frac{1}{h_1h_2h_3}\left[
\frac{\partial}{\partial q_1}\left(\frac{h_2h_3}{h_1}\frac{\partial u}{\partial q_1}\right)
+\frac{\partial}{\partial q_2}\left(\frac{h_3h_1}{h_2}\frac{\partial u}{\partial q_2}\right)
+\frac{\partial}{\partial q_3}\left(\frac{h_1h_2}{h_3}\frac{\partial u}{\partial q_3}\right)
\right]. \tag{21}
$$

将柱坐标系和球坐标系中的度规系数代入上式，就可以得到前面的结果 (17) 和 (18)。至于旋度，我们就不讨论了。

## §2 球坐标系中的分离变量

### §2.1 *Helmholtz* 方程

在三维空间，我们已经介绍过的数理方程有三种形式，即波动方程、输运方程和稳定场方程。首先看齐次的波动方程：

$$
\frac{\partial^2u}{\partial t^2}-a^2\nabla^2u=0.
$$

分离变量，令 $u(\boldsymbol r,t)=v(\boldsymbol r)T(t)$，代入上式，可得

$$
\nabla^2v+k^2v=0 \tag{22}
$$

和

$$
T''+k^2a^2T=0,
$$

其中 $k^2$ 是分离变量时引入的常数，记作 $k^2$ 是物理上的习惯。式 (22) 称为 Helmholtz 方程，加上适当的边界条件即可构成本征值问题，从而可以求解 $v$ 和常数 $k^2$（即本征值），合理的边界条件通常导致 $k^2\ge0$（所以将它记作 $k^2$ 是恰当的）。从物理上说，如果 $k^2<0$，则 $k$ 为虚数，于是有 $T(t)=e^{\pm|k|at}$，这显然不是波动，不符合问题的物理背景，由此也可推测 $k^2\ge0$。当然，如果推导方程或边界条件时作了不恰当的简化或近似，也可以导致不合理的结果。所以，可靠的结论还是应该来源于本征值问题的求解。目前我们承认 $k^2\ge0$。

其次看齐次的输运方程：

$$
\frac{\partial u}{\partial t}-a^2\nabla^2u=0.
$$

分离变量，令 $u(\boldsymbol r,t)=v(\boldsymbol r)T(t)$，代入上式，可得式 (22) 和

$$
T'+k^2a^2T=0.
$$

这里我们再次得到 Helmholtz 方程，适当的边界条件同样导致 $k^2\ge0$。从物理上说，如果 $k^2<0$，则 $k$ 为虚数，于是有 $T(t)=e^{|k|^2a^2t}$，这显然是不合理的，由此也可推测 $k^2\ge0$。

再次看齐次的稳定场方程：

$$
\nabla^2u=0.
$$

就其形式而言，这可以当作 Helmholtz 方程 (22) 当 $k=0$ 时的特殊情况。

综上所述，三类方程的齐次形式都可以归结为 Helmholtz 方程。这与坐标系无关，所以，研究各种曲线坐标系中的分离变量，都可以用 Helmholtz 方程作为代表。

所以下面我们就研究 Helmholtz 方程在球坐标系中的分离变量。

### §2.2 *Helmholtz* 方程在球坐标系中的分离变量

在球坐标系中，Helmholtz 方程的具体形式是

$$
\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial u}{\partial r}\right)
+\frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial u}{\partial\theta}\right)
+\frac{1}{r^2\sin^2\theta}\frac{\partial^2u}{\partial\phi^2}+k^2u=0. \tag{23}
$$

这里将函数写作 $u$。令 $u(\boldsymbol r)=R(r)Y(\theta,\phi)$，代入上式，两边同乘以 $r^2$，同除以 $R(r)Y(\theta,\phi)$，适当移项，可得

$$
\frac{1}{R}\left[\frac{d}{dr}\left(r^2\frac{dR}{dr}\right)+k^2r^2R\right]
=-\frac{1}{Y}\left[\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial Y}{\partial\theta}\right)+\frac{1}{\sin^2\theta}\frac{\partial^2Y}{\partial\phi^2}\right],
$$

上式左边与 $\theta、\phi$ 无关，右边与 $r$ 无关，两边相等，则与 $r、\theta、\phi$ 均无关，即为常数，记作 $\lambda$，于是得到角向方程

$$
\frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial Y}{\partial\theta}\right)
+\frac{1}{\sin^2\theta}\frac{\partial^2Y}{\partial\phi^2}+\lambda Y=0, \tag{24}
$$

和径向方程

$$
r^2R''+2rR'+(k^2r^2-\lambda)R=0. \tag{25}
$$

求解的顺序是：先由角向边界条件（见下）与角向方程 (24) 构成本征值问题求解 $Y(\theta,\phi)$ 并确定本征值 $\lambda$。然后将 $\lambda$ 代入径向方程 (25)，这时候需要分别讨论 Helmholtz 方程或 Laplace 方程。如果是 Helmholtz 方程，则由径向边界条件与径向方程构成本征值问题求解 $R(r)$ 并确定本征值 $k$，这时径向方程是球 Bessel 方程，需要用级数法求解；如果是 Laplace 方程，则 $k=0$ 是已知的，不需要由本征值问题确定，且径向方程简化为 Euler 方程，这很容易求解。Laplace 方程比较简单，后面会详细讨论。

角向方程是偏微分方程，需要进一步分离变量。令 $Y(\theta,\phi)=H(\theta)\Phi(\phi)$，代入式 (24)，两边同乘以 $\sin^2\theta$，同除以 $H(\theta)\Phi(\phi)$，适当移项，可得

$$
\frac{1}{H}\left[\sin\theta\frac{d}{d\theta}\left(\sin\theta\frac{dH}{d\theta}\right)+\lambda\sin^2\theta H\right]
=-\frac{\Phi''}{\Phi},
$$

上式左边与 $\phi$ 无关，右边与 $\theta$ 无关，两边相等，则与 $\theta、\phi$ 均无关，即为常数，记作 $\mu$，于是得到两个方程：

$$
\frac{1}{\sin\theta}\frac{d}{d\theta}\left(\sin\theta\frac{dH}{d\theta}\right)
+\left(\lambda-\frac{\mu}{\sin^2\theta}\right)H=0, \tag{26}
$$

和

$$
\Phi''+\mu\Phi=0. \tag{27}
$$

现考虑角向方程的边界条件和相应的本征值问题。由于 $(r,\theta,\phi+2\pi)$ 与 $(r,\theta,\phi)$ 表示同一几何点，而物理量在一个几何点的取值应该是确定的，不应该依赖于该点的数学描述，所以应有 $u(r,\theta,\phi+2\pi)=u(r,\theta,\phi)$，或 $R(r)H(\theta)\Phi(\phi+2\pi)=R(r)H(\theta)\Phi(\phi)$，但 $R(r)H(\theta)$ 不恒为零，故

$$
\Phi(\phi+2\pi)=\Phi(\phi). \tag{28}
$$

此即周期性边界条件，它与方程 (27) 构成本征值问题。我们已经在第七章求解了这一本征值问题，结果是

$$
\mu=m^2,\qquad \Phi(\phi)=\{e^{im\phi},e^{-im\phi}\},\qquad m=0,1,2,\cdots. \tag{29}
$$

在球坐标系，习惯上将本征函数写成指数形式，如上式所示，而在柱坐标系或平面极坐标系，习惯上将本征函数写成三角形式。它们是互相等价的。

再考虑 $\theta$ 方向的边界条件，现在 $Y(\theta,\phi)=H(\theta)e^{\pm im\phi}$。如果 $m\ne0$，应该要求 $H(0)=H(\pi)=0$，否则当 $\theta=0、\pi$ 时，$Y(\theta,\phi)$ 没有定义，因为此时 $\phi$ 没有定义；如果 $m=0$，则 $Y(\theta,\phi)=H(\theta)$，没有类似的问题，但 $H(0)、H(\pi)$ 应该有限。如前所述，原则上物理量处处都应该有限，但由于 $\theta=0、\pi$ 是球坐标系的奇点，微分方程的解在这些奇点处很容易出现奇性，所以对于物理量在这些地方的有限性需要特别加以强调。上述边界条件也是一种自然边界条件。将本征值 $\mu=m^2$ 代回方程 (26)，令 $x=\cos\theta$，$H(\theta)=P(x)$，将它化为 $P(x)$ 对 $x$ 的微分方程，考虑到以上自然边界条件，即得本征值问题

$$
\frac{d}{dx}\left[(1-x^2)\frac{dP}{dx}\right]+\left(\lambda-\frac{m^2}{1-x^2}\right)P=0, \tag{30a}
$$

$$
P(\pm1)=0\quad(m\ne0)\qquad\text{或}\qquad |P(\pm1)|<\infty\quad(m=0). \tag{30b}
$$

式 (30a) 称为连带 (associated) Legendre 方程。对于轴对称问题，取对称轴为球坐标的极轴，则问题的解 $u(r,\theta,\phi)$ 与 $\phi$ 无关，这时只需考虑 $m=0$，上述本征值问题就简化为

$$
\frac{d}{dx}\left[(1-x^2)\frac{dP}{dx}\right]+\lambda P=0, \tag{31a}
$$

$$
|P(\pm1)|<\infty. \tag{31b}
$$

即使问题并没有轴对称性，上式也是本征值问题 (30) 的一种特殊情况，也需要加以研究。总之，这是最重要的一种情况。式 (31a) 称为 Legendre 方程，需要用级数法求解，后面将会详细研究。至于连带 Legendre 方程，可以通过适当的变换将它与 Legendre 方程联系起来。

**习题**　微观粒子在中心力场 $V(r)$ 中的定态 Schrödinger 方程是

$$
-\frac{\hbar^2}{2\mu}\nabla^2u+V(r)u=Eu,
$$

其中 $u(r,\theta,\phi)$ 是波函数，$\mu$ 是粒子质量，$E$ 是能量，$\hbar=h/2\pi$，其中 $h$ 是 Planck 常数，$\hbar$ 亦称为 Planck 常数。试在球坐标系中对该方程分离变量。

## §3 柱坐标系中的分离变量

在柱坐标系中对 Helmholtz 方程分离变量，涉及的求解区域有许多不同情况，比如 $\rho$ 的取值范围可以是有限或半无限，$z$ 的取值范围可以是有限、无限或半无限，等等，对于确定的区域，边界条件的情况还可以不同（比如对有限的圆柱区域，可以在柱面上有齐次边界条件，也可以在上下底有齐次边界条件）。而在球坐标系，对于物理上比较常见的情况，角向坐标的取值范围都是固定的，边界条件就是自然边界条件（即上节所讨论的情况）。所以柱坐标系的情况比球坐标系要复杂。在柱坐标系中求解数理方程（即使是 Laplace 方程）的定解问题，最好是对每个具体问题都从分离变量开始做（“从猿到人”）。在计算经验不够丰富的时候，这样做比较不容易出错。有些书上列出了各种边界条件下的一般解形式，但是记忆这些结果有一定的困难，而且容易出错（比如遗漏掉一个本征值和相应的本征函数）。如果读者对这些一般结果有兴趣，也应该在熟练掌握求解过程的基础上自己加以总结，而不是简单地照搬书上的结果。

本节以圆柱内的 Laplace 方程定解问题为例进行分离变量，并假定在 $\rho=a$ 处有齐次边界条件，$z$ 的取值范围可以是有界的或半无界的，暂时不作明确限制。具体写出来是

$$
\frac{1}{\rho}\frac{\partial}{\partial\rho}\left(\rho\frac{\partial u}{\partial\rho}\right)
+\frac{1}{\rho^2}\frac{\partial^2u}{\partial\phi^2}
+\frac{\partial^2u}{\partial z^2}=0,\qquad (\rho<a), \tag{32a}
$$

$$
\left.\left(\alpha\frac{\partial u}{\partial\rho}+\beta u\right)\right|_{\rho=a}=0,
\qquad (\alpha\ge0,\ \beta\ge0,\ \alpha^2+\beta^2\ne0). \tag{32b}
$$

注意这里 $z$ 的取值范围可以是有界的或半无界的，因未具体指明，所以定解条件也不完整。不过这对于分离变量没有影响。令 $u(\rho,\phi,z)=R(\rho)\Phi(\phi)Z(z)$，代入式 (32a)，两边同乘以 $\rho^2$，同除以 $R(\rho)\Phi(\phi)Z(z)$，适当移项，可得

$$
\frac{\rho}{R}\frac{d}{d\rho}\left(\rho\frac{dR}{d\rho}\right)+\rho^2\frac{Z''}{Z}=-\frac{\Phi''}{\Phi},
$$

上式左边与 $\phi$ 无关，右边与 $\rho、z$ 无关，两边相等，则与 $\rho、\phi、z$ 均无关，即为常数，记作 $\mu$，于是得到两个方程：

$$
\frac{1}{R\rho}\frac{d}{d\rho}\left(\rho\frac{dR}{d\rho}\right)-\frac{\mu}{\rho^2}=-\frac{Z''}{Z} \tag{33}
$$

和

$$
\Phi''+\mu\Phi=0. \tag{34}
$$

先考虑角向方程 (34) 的边界条件和相应的本征值问题。由于 $(\rho,\phi+2\pi,z)$ 与 $(\rho,\phi,z)$ 表示同一几何点，与上节类似，有 $u(\rho,\phi+2\pi,z)=u(\rho,\phi,z)$，或 $R(\rho)\Phi(\phi+2\pi)Z(z)=R(\rho)\Phi(\phi)Z(z)$，但 $R(\rho)Z(z)$ 不恒为零，故得周期性边界条件

$$
\Phi(\phi+2\pi)=\Phi(\phi). \tag{35}
$$

它与方程 (34) 构成本征值问题。这是我们所熟悉的。我们将结果重新写在下面

$$
\mu=m^2,\qquad \Phi(\phi)=\{\cos m\phi,\sin m\phi\},\qquad m=0,1,2,\cdots. \tag{36}
$$

再看式 (33)，它的左边与 $z$ 无关，右边与 $\rho$ 无关，两边相等，则与 $\rho、z$ 均无关，即为常数，记作 $-\lambda$，于是得到两个方程：

$$
Z''-\lambda Z=0, \tag{37}
$$

和

$$
\frac{1}{\rho}\frac{d}{d\rho}\left(\rho\frac{dR}{d\rho}\right)+\left(\lambda-\frac{m^2}{\rho^2}\right)R=0. \tag{38a}
$$

边界条件 (32b) 在分离变量后成为

$$
\alpha R'(a)+\beta R(a)=0. \tag{38b}
$$

由于现在的解具有形式 $u(\rho,\phi,z)=R(\rho)\{\cos m\phi,\sin m\phi\}Z(z)$，而 $\rho=0$ 处 $\phi$ 没有定义，与上节类似，应该有自然边界条件

$$
R(0)=0\quad(m\ne0)\qquad\text{或}\qquad |R(0)|<\infty\quad(m=0). \tag{38c}
$$

式 (38a)–(38c) 构成本征值问题。以后可以证明，本征值 $\lambda\ge0$（Sturm–Liouville 本征值问题的一般结论）。

如果 $\lambda>0$，令 $x=\sqrt{\lambda}\rho$，$R(\rho)=y(x)$，则方程 (38a) 化为

$$
\frac{1}{x}\frac{d}{dx}\left(x\frac{dy}{dx}\right)+\left(1-\frac{m^2}{x^2}\right)y=0
\qquad\text{或}\qquad
y''+\frac{1}{x}y'+\left(1-\frac{m^2}{x^2}\right)y=0. \tag{39}
$$

这是 $m$ 阶 Bessel 方程的标准形式。Bessel 方程需要用级数法求解。本征值 $\lambda$ 的具体取值由边界条件 (38b)–(38c) 确定。确定本征值之后，不难求解式 (37)。

如果 $\lambda=0$，则方程 (38a) 化为

$$
\rho^2R''+\rho R'-m^2R=0. \tag{40}
$$

这是 Euler 方程，其解为 $R(\rho)=\{\rho^m,\rho^{-m}\}$（当 $m\ne0$），或 $R(\rho)=\{1,\ln\rho\}$（当 $m=0$）。当 $m\ne0$，其中只有 $R(\rho)=\rho^m$ 满足边界条件 (38c)，而边界条件 (38b) 化为 $m\alpha+\beta a=0$，由于 $m>0、a>0$，而 $\alpha、\beta$ 均非负且不同时为零，故不满足。当 $m=0$，只有 $R(\rho)=1$ 满足边界条件 (38c)，而边界条件 (38b) 化为 $\beta=0$，如果 $a$ 处为第二类齐次边界条件，则可以满足。即只当 $a$ 处为第二类齐次边界条件且 $m=0$ 时，$\lambda=0$ 才是一个本征值。对应于 $\lambda=0$，亦不难求解式 (37)。

各变量的常微分方程都解出以后，可以写出 $u(\rho,\phi,z)$ 的一般解，其中的任意常数由 $z$ 方向的边界条件决定。

从以上讨论可以看出，需要解决的关键问题是求解分离变量过程中出现的常微分方程及其本征值问题。有些常微分方程是我们所熟悉的，其解为初等函数。另一些则没有现成的解，需要加以研究。级数法是一种常用的解法，这将在下一章仔细讨论。但应该指出，能够用级数法求解的方程也是非常有限的。

将方程 (32a) 的求解区域改为 $0\le z\le h$，且 $z=0、h$ 处有齐次边界条件

$$
\left.\left(\alpha\frac{\partial u}{\partial z}-\beta u\right)\right|_{z=0}=0,
\qquad
\left.\left(\gamma\frac{\partial u}{\partial z}+\delta u\right)\right|_{z=h}=0, \tag{41}
$$

其中 $\alpha、\beta、\gamma、\delta\ge0$，但 $\alpha^2+\beta^2\ne0$，$\gamma^2+\delta^2\ne0$；而 $\rho$ 的取值范围暂时不作限制。则式 (36) 以下的讨论应该加以修改。此时式 (33) 两边仍为常数，记作 $\lambda$（由以下讨论可知 $\lambda\ge0$，所以这是方便的记法），于是得到两个方程：

$$
\frac{1}{\rho}\frac{d}{d\rho}\left(\rho\frac{dR}{d\rho}\right)-\left(\lambda+\frac{m^2}{\rho^2}\right)R=0, \tag{42}
$$

和

$$
Z''+\lambda Z=0. \tag{43}
$$

$z=0、h$ 处的齐次边界条件在分离变量后成为

$$
\alpha Z'(0)-\beta Z(0)=0,\qquad \gamma Z'(h)+\delta Z(h)=0. \tag{44}
$$

这与方程 (43) 构成本征值问题。以后可以证明，本征值 $\lambda\ge0$（如果 $z=0、h$ 处都不涉及第三类边界条件，则该结论是我们所熟悉的）。直接计算容易验证，只有 $z=0、h$ 处都是第二类齐次边界条件时，$\lambda=0$ 才是一个本征值。

如果 $\lambda=0$，则式 (42) 化为 Euler 方程 (40)，其解已在前面给出。

如果 $\lambda>0$，令 $x=\sqrt{\lambda}\rho$，$R(\rho)=y(x)$，则式 (42) 化为

$$
y''+\frac{1}{x}y'-\left(1+\frac{m^2}{x^2}\right)y=0. \tag{45}
$$

这称为 $m$ 阶虚宗量 Bessel 方程，因为只要令 $t=ix$，它就变成对 $t$ 的 $m$ 阶 Bessel 方程。它的解法与 Bessel 方程类似。

各变量的常微分方程都解出以后，可以写出 $u(\rho,\phi,z)$ 的一般解，其中的任意常数由 $\rho$ 方向的边界条件决定。

## 补充习题

1. 研究二维波动方程在平面极坐标系中的分离变量。

2. 研究二维输运方程在平面极坐标系中的分离变量。
