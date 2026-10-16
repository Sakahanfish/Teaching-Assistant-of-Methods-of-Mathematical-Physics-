# 第十一章 数学物理方程和定解条件

作为本书的“数学物理方程”部分的开始，在这一章里，我们先从一些物理问题中导出一些典型的二阶线性偏微分方程。以后再讨论这些方程的一般性质及解法。

从实际的物理问题出发，确定最方便用来描写该问题的物理量，经过抽象简化，抓住主要因素，作出合理的近似，建立起它所满足的微分方程。用数学语言描述物理问题，是我们认识物理世界的基础，是所有各门物理类课程都要讨论的内容，本课程只能就方程的推导举出几个经典例子，只能起到一个示范和引导的作用。

## §11.1 波动方程

### 弦的横振动方程

有一根完全柔软的均匀轻弦，绷紧，而后以某种方法激发，使弦在一个平面内作小振动。列出弦的横振动方程。

取弦的平衡位置为 $x$ 轴，两个端点的坐标分别为 $x=0$ 和 $x=l$。这样弦上每一点都可以用坐标 $x$ 标记。设 $u(x,t)$ 为弦上一点 $x$ 在 $t$ 时刻的横向位移。像在力学问题中常用的那样，采用微元分析法，在弦上隔离出长为 $\Delta x$ 的一小段（如图 11.1 所示），这一段弦是如此之小，以至于可以把它看成是质点。这一小段弦在两个端点 $x$ 及 $x+\Delta x$ 处受到弹性力的作用。因为弦完全柔软，故弹性力只有切向应力（张力），而无法向应力。$t$ 时刻弦上 $x$ 点的切线与 $x$ 正方向的夹角记为 $\theta(x,t)$，切线的斜率

$$
\tan\theta(x,t)
=
\left.\frac{\partial u}{\partial x}\right|_x.
$$

弦本身很轻，相对于张力而言，重力可以忽略。因此由 Newton 第二定律，可以写出 $t$ 时刻这一小段弦在 $x$ 方向及垂直方向上的动力学方程：

$$
\left.(T\cos\theta)\right|_{x+\Delta x}
-
\left.(T\cos\theta)\right|_x
=
0,
\tag{11.1}
$$

$$
\left.(T\sin\theta)\right|_{x+\Delta x}
-
\left.(T\sin\theta)\right|_x
=
\rho\Delta x\frac{\partial^2u}{\partial t^2}.
\tag{11.2}
$$

![Figure 11.1](...)

图 11.1　弦的横振动

$$
\tan\theta(x,t)
=
\left.\frac{\partial u}{\partial x}\right|_{(x,t)}.
$$

其中 $\left.T\right|_x$ 表示弦上 $x$ 点处切向应力的大小，$\rho$ 是弦的线密度，即单位长度的质量。弦均匀，不仅意味着弦的质量分布均匀，即线密度 $\rho$ 是常数，与 $x$ 无关，同时还意味着弦的弹性性质均匀。$\partial^2u/\partial t^2$ 表示微元的平均加速度。

在小振动近似下，$x+\Delta x$ 与 $x$ 两点间任一时刻横向位移之差 $u(x+\Delta x,t)-u(x,t)$ 与 $\Delta x$ 相比是一个小量，即 $|\partial u/\partial x|\ll1$。因此，在准确到小量 $\partial u/\partial x$ 的一级项的条件下，

$$
\sin\theta
\approx
\tan\theta
=
\frac{\partial u}{\partial x}
\qquad
\text{（略去了 $\partial u/\partial x$ 的三级项）},
\tag{11.3}
$$

$$
\cos\theta
\approx
1
\qquad
\text{（略去了 $\partial u/\partial x$ 的二级项）}.
\tag{11.4}
$$

将 (11.4) 式代入 (11.1) 式，就有

$$
\left.T\right|_{x+\Delta x}
-
\left.T\right|_x
=
0,
\qquad
\text{即}
\qquad
\left.T\right|_{x+\Delta x}
=
\left.T\right|_x,
\tag{11.5}
$$

表示 $T$ 与 $x$ 无关，即弦中各点的张力相等。将 (11.3) 和 (11.5) 式代入 (11.2) 式，即得

$$
\rho\Delta x\frac{\partial^2u}{\partial t^2}
=
T
\left(
\left.\frac{\partial u}{\partial x}\right|_{x+\Delta x}
-
\left.\frac{\partial u}{\partial x}\right|_x
\right).
$$

两端同时除以 $\Delta x$，并令 $\Delta x\to0$，得

$$
\rho\frac{\partial^2u}{\partial t^2}
-
T\frac{\partial^2u}{\partial x^2}
=
0.
$$

定义 $a=\sqrt{T/\rho}$①，则方程可以写成

$$
\boxed{
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0.
}
\tag{11.6}
$$

这就是自由振动的弦满足的方程。这里所谓“自由”表示弦上各点除端点外都不受外力作用。振动沿着弦（即沿 $x$ 方向）传播，而振动的方向与振动传播的方向垂直，即 $u$ 的方向与 $x$ 的方向互相垂直，这种振动称为横振动。

在小振动的条件下，还可以证明张力 $T$ 与时间 $t$ 无关。这是因为这一段弦的伸长

$$
\Delta s-\Delta x
=
\sqrt{\Delta u^2+\Delta x^2}
-
\Delta x
=
\left[
\sqrt{
1+
\left(
\frac{\partial u}{\partial x}
\right)^2
}
-
1
\right]\Delta x
=
O
\left[
\left(
\frac{\partial u}{\partial x}
\right)^2
\right],
$$

所以，弦的总长度不随时间变化。因此，由 Hooke 定律可知，张力 $T$ 不随时间变化。前面又已经证明，$T$ 也不随 $x$ 变化，所以 $T$ 是一个常量。

如果弦在位移 $u$ 的方向上还受到外力的作用，设单位长度所受的外力为 $f$，其正方向沿 $u$ 的正方向。仿照前面的推导，有

$$
\rho\Delta x\frac{\partial^2u}{\partial t^2}
=
T
\left(
\left.\frac{\partial u}{\partial x}\right|_{x+\Delta x}
-
\left.\frac{\partial u}{\partial x}\right|_x
\right)
+
f\Delta x.
$$

因此，受迫振动的弦满足非齐次方程

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
\frac f\rho,
\tag{11.7}
$$

其中的非齐次项 $f/\rho$ 是单位质量所受的外力。

① 通过量纲分析可知，$a$ 的量纲是速度量纲。以后将会看到，$a$ 就是弦的振动传播速率。

### 杆的纵振动方程

对于弹性杆的纵振动方程，可以类似地处理。

考虑一均匀细轻杆，沿杆长方向作小振动。列出杆的纵振动方程。

如图 11.2 所示，取杆长方向为 $x$ 轴方向。假设在垂直杆长方向的任一截面上各点的振动状况（即位移）完全相同，因此垂直于杆长方向的各截面均可用它的平衡位置 $x$ 标记。设在任一时刻 $t$，此截面相对于平衡位置的位移为 $u(x,t)$，约定沿 $x$ 正方向的位移为正。

仍然采用微元分析法。在杆中隔离出一小段 $(x,x+\Delta x)$，分析它所受的弹性力。通过截面 $x$，微元受到弹性力 $P(x,t)S$ 的作用（$P(x,t)$ 为 $x$ 处单位面积所受弹性力，即应力，规定沿 $x$ 方向为正），通过截面 $x+\Delta x$ 受到弹性力 $P(x+\Delta x,t)S$ 的作用。在 $t$ 时刻，两截面的位移分别为 $u(x,t)$ 与 $u(x+\Delta x,t)\approx u+\Delta u$（二者相差高阶小量）。对此小段应用 Newton 第二定律，即得

$$
\rho S\Delta x\frac{\partial^2u}{\partial t^2}
=
[P(x+\Delta x,t)-P(x,t)]S,
$$

其中 $\rho$ 为杆的体密度，$S$ 为垂直于杆长方向横截面的面积，$\partial^2u/\partial t^2$ 是这一小段杆的平均加速度。两端同时除以 $S\Delta x$，并令 $\Delta x\to0$，得

$$
\rho\frac{\partial^2u}{\partial t^2}
=
\frac{\partial P}{\partial x}.
\tag{11.8}
$$

如果略去垂直杆长方向的形变，根据 Hooke 定律，应力大小与应变 $\partial u/\partial x$ 成正比，

$$
\left.P\right|_x
=
E\frac{\partial u}{\partial x},
\tag{11.9}
$$

比例系数 $E$ 称为杨氏模量，是一个物质常数。将此式代入 (11.8) 式，就得到杆的纵振动方程

$$
\boxed{
\frac{\partial^2u}{\partial t^2}
-
a^2\frac{\partial^2u}{\partial x^2}
=
0,
}
\tag{11.10}
$$

其中 $a=\sqrt{E/\rho}$。这里，位移发生在杆的方向上，波动也是沿着杆的方向传播。这种位移与波动传播方向相同的振动，称为纵振动。我们看到，杆的纵振动与弦的横振动机理并不完全相同，但它们满足的偏微分方程 (11.10) 和 (11.6) 式的形式却完全一样。这一类方程统称为波动方程。

![Figure 11.2](...)

图 11.2　杆的纵振动用，应力与应变

更一般地，在三维空间中的波动方程是

$$
\boxed{
\frac{\partial^2u}{\partial t^2}
-
a^2\nabla^2u
=
0,
}
\tag{11.11}
$$

其中

$$
\nabla^2u
=
\nabla\cdot(\nabla u),
\qquad
\nabla^2
\equiv
\frac{\partial^2}{\partial x^2}
+
\frac{\partial^2}{\partial y^2}
+
\frac{\partial^2}{\partial z^2}
=
\nabla\cdot\nabla
\tag{11.12}
$$

称为 Laplace 算符。

## §11.2 热传导方程

推导热传导方程所用的数学方法和上面的完全相同。不同之处在于具体的物理规律不同。这里用到的是热学方面的两个基本规律，即能量守恒定律和热传导的 Fourier 定律。

扼要介绍一下 Fourier 定律。设有一块连续介质。取定坐标系，用 $u(x,y,z,t)$ 表示介质内空间坐标为 $(x,y,z)$ 的一点在 $t$ 时刻的温度。若沿 $x$ 方向存在温度差，则在 $x$ 方向就一定有热量的传递。当温度变化不大时，单位时间内通过垂直 $x$ 方向的单位面积的热量 $q_x$ 与温度沿 $x$ 方向的空间变化率成正比，

$$
\boxed{
q_x
=
-k\frac{\partial u}{\partial x},
}
\tag{11.13}
$$

$q_x$ 称为沿 $x$ 方向的热流密度，沿 $x$ 的正方向为正。上式中的负号表示热流方向和温度变化的方向相反，即热量由高温流向低温。$k$ 称为热导率，与介质的材料有关。严格说来，$k$ 与温度 $u$ 也有关系。但如果温度的变化幅度不大，则 $k$ 近似地与 $u$ 无关。

在三维各向同性的均匀介质中，$k$ 与热流方向以及空间位置均无关。因此，在 $x,y,z$ 三个方向上都存在温度差，则热流密度矢量 $\boldsymbol q$ 与温度梯度 $\nabla u$ 成正比，

$$
q_x
=
-k\frac{\partial u}{\partial x},
\qquad
q_y
=
-k\frac{\partial u}{\partial y},
\qquad
q_z
=
-k\frac{\partial u}{\partial z},
\tag{11.14}
$$

或者

$$
\boxed{
\boldsymbol q
=
-k\nabla u.
}
\tag{11.14'}
$$

(11.13)，(11.14) 或 (11.14') 诸式即为 Fourier 定律。如果讨论的介质各向异性，则 (11.14') 式应改写成

$$
\boldsymbol q
=
-\boldsymbol K\cdot\nabla u.
\tag{11.15}
$$

这里的 $\boldsymbol K$ 是一个 $3\times3$ 矩阵，它和 $\nabla u$ 按矩阵乘法的规则相乘。

下面根据 Fourier 定律和能量守恒定律来推导热传导方程。设想在介质内部隔离出一个长方体（见图 11.3），六个面都和坐标面重合。$\Delta t$ 时间内沿 $x$ 方向流入长方体的热量为

$$
\left(
\left.q_x\right|_x
-
\left.q_x\right|_{x+\Delta x}
\right)
\Delta y\Delta z\Delta t
=
-\frac{\partial q_x}{\partial x}
\Delta x\Delta y\Delta z\Delta t.
$$

同理，在 $\Delta t$ 时间内沿 $y$ 方向流入长方体的热量为

$$
\left(
\left.q_y\right|_y
-
\left.q_y\right|_{y+\Delta y}
\right)
\Delta x\Delta z\Delta t
=
-\frac{\partial q_y}{\partial y}
\Delta x\Delta y\Delta z\Delta t,
$$

在 $\Delta t$ 时间内沿 $z$ 方向流入长方体的热量为

$$
\left(
\left.q_z\right|_z
-
\left.q_z\right|_{z+\Delta z}
\right)
\Delta x\Delta y\Delta t
=
-\frac{\partial q_z}{\partial z}
\Delta x\Delta y\Delta z\Delta t.
$$

![Figure 11.3](...)

图 11.3　热传导方程的推导：位于 $(x,y,z)$ 处的小长方体

如果在介质内有热量产生或消耗（例如存在化学反应，或通有电流，等等），单位时间内单位体积中产生的热量为 $F(x,y,z,t)$，则应有

$$
-\nabla\cdot\boldsymbol q\,
\Delta x\Delta y\Delta z\Delta t
+
F(x,y,z,t)
\Delta x\Delta y\Delta z\Delta t
=
\rho\Delta x\Delta y\Delta z\cdot c\cdot\Delta u,
$$

其中 $\rho$ 是介质的体密度，$c$ 是比热容，即

$$
\boxed{
\frac{\partial(\rho cu)}{\partial t}
+
\nabla\cdot\boldsymbol q
=
F(x,y,z,t).
}
\tag{11.16}
$$

此方程常称为连续性方程。对各向异性介质，将 (11.15) 式代入方程 (11.16)，得

$$
\boxed{
\frac{\partial(\rho cu)}{\partial t}
-
\nabla\cdot
(\boldsymbol K\cdot\nabla u)
=
F(x,y,z,t).
}
\tag{11.17}
$$

对各向同性均匀介质，$\rho,c$ 和 $k$ 都是常数，将 (11.14') 式代入方程 (11.16)，得

$$
\boxed{
\frac{\partial u}{\partial t}
-
\kappa\nabla^2u
=
\frac{F}{\rho c}
=
f,
}
\tag{11.18}
$$

其中 $\kappa=k/(\rho c)$ 称为扩散率或温度传导率。如果介质内部没有其他热量来源或消耗，$F=0$，则

$$
\boxed{
\frac{\partial u}{\partial t}
-
\kappa\nabla^2u
=
0.
}
\tag{11.19}
$$

从分子运动的角度看，温度的高低是分子热运动激烈程度的反映。分子热运动的不平衡，通过碰撞交换能量，在宏观上就表现为热量的传递。可以设想，如果介质内存在别种不均匀状况，例如物质浓度的不均匀，通过分子的运动也会发生物质的交换，宏观上就表现为分子的扩散。微观机理上的这种相似性，就决定了扩散方程和热传导方程有相同的形式，

$$
\boxed{
\frac{\partial u}{\partial t}
-
D\nabla^2u
=
f(x,y,z,t),
}
\tag{11.20}
$$

其中的 $u(x,y,z,t)$ 代表分子浓度，$D$ 是扩散率，$f(x,y,z,t)$ 则是单位时间内在单位体积中该种分子的产率。

方程 (11.17)～(11.20) 式均称为热传导方程或者扩散方程，其中 (11.19) 式是齐次的，描写内部无热源的热传导过程（或者内部无源或汇的扩散过程），(11.17) 式、(11.18) 式和 (11.20) 式都是非齐次的，描写介质内部有热源的热传导过程（或者内部有源或汇的扩散过程）。

## §11.3 稳定问题

在一定条件下，物体温度达到稳定（温度不随时间变化，即 $\partial u/\partial t=0$）时，则温度分布满足 Poisson 方程

$$
\boxed{
\nabla^2u
=
-\frac f\kappa.
}
\tag{11.21}
$$

特别是，如果 $f=0$，则有 Laplace 方程

$$
\boxed{
\nabla^2u=0.
}
\tag{11.22}
$$

在电磁学中，将电场强度 $\boldsymbol E$ 与电势 $u$ 之间的关系式

$$
\boldsymbol E
=
-\nabla u,
$$

代入 Gauss 定理（的微分形式）

$$
\nabla\cdot\boldsymbol E
=
\frac{\rho}{\varepsilon_0},
$$

也能得到静电势 $u$ 满足的 Poisson 方程

$$
\nabla^2u
=
-\frac{\rho}{\varepsilon_0},
\tag{11.23}
$$

其中 $\rho$ 是电荷密度，$\varepsilon_0$ 为真空电容率（真空介电常数）。在电荷密度 $\rho=0$ 的区域内，静电势则满足 Laplace 方程。

如果齐次波动方程

$$
\frac{\partial^2u}{\partial t^2}
-
a^2\nabla^2u
=
0
$$

中，$u(x,y,z,t)$ 随时间周期地变化，频率为 $\omega$，

$$
u(x,y,z,t)
=
v(x,y,z)e^{-i\omega t},
\tag{11.24}
$$

则 $v(x,y,z)$ 满足 Helmholtz 方程

$$
\boxed{
\nabla^2v(x,y,z)
+
k^2v(x,y,z)
=
0,
}
\tag{11.25}
$$

其中 $k=\omega/a$ 称为波数。

以后我们还将看到，将波动方程或热传导方程分离变量，或者对它们作 Laplace 变换或是 Fourier 变换，也会得到 Helmholtz 方程。

以上介绍了几种基本的偏微分方程。从物理上看，有反映波动过程的波动方程，有反映热传导过程的热传导方程，也有反映稳恒状态的 Poisson 方程和 Laplace 方程。从数学上看，这也正好相应地分为三类：波动方程，在数学上属于双曲型方程；热传导方程，属于抛物型方程；而 Poisson 方程和 Laplace 方程，则属于椭圆型方程。求解这三类方程构成的定解问题，将是本课程的中心任务。

## §11.4 定解条件

上面几节建立的偏微分方程，并不能唯一地、确定地描写某一个具体的物理过程。正如只根据 Newton 第二定律列出的动力学方程（二阶常微分方程）并不能唯一地确定质点的运动（通解中含有两个任意常数）一样。要完全确定一个质点的运动，除了常微分方程之外，还必须有初始条件。

对于偏微分方程来说，情况还要更复杂一些。如果我们能求得二阶线性偏微分方程的通解，它应含有两个任意函数。这是容易理解的。例如，偏微分方程

$$
\frac{\partial^2u(x,y)}{\partial x^2}
=
0
$$

的通解就是

$$
u(x,y)
=
c_1(y)+xc_2(y),
$$

其中 $c_1(y)$ 和 $c_2(y)$ 是 $y$ 的任意函数。

仅有偏微分方程，解并不唯一，换句话说，还不足以唯一地决定介质中发生的具体物理过程，因为在推导偏微分方程时，只考虑了介质的内部，并没有考虑介质（通过边界面）与外界的相互作用。因此，微分方程只适用于介质内部。而且，如果问题与时间有关的话，在推导方程时也并没有考虑介质的历史状况。为了完全描写一个具有确定解的物理问题，在数学上就是要构成一个定解问题，这除了偏微分方程之外，还必须有定解条件，包括初始条件和边界条件。

先讨论初始条件。初始条件应该完全描写初始时刻（通常定为 $t=0$）介质内部以及边界上任意一点的状况。对于波动方程，就应该给出初始时刻的位移和速度（如果是力学问题的话），

$$
\boxed{
\left.u\right|_{t=0}
=
\phi(x,y,z),
\qquad
\left.\frac{\partial u}{\partial t}\right|_{t=0}
=
\psi(x,y,z),
\qquad
(x,y,z)\in\overline V.
}
\tag{11.26}
$$

对于热传导方程，由于方程中只出现未知函数 $u(x,y,z,t)$ 对 $t$ 的一阶偏微商，所以只需给出体系内部以及边界上任意一点在初始时刻的温度

$$
\boxed{
\left.u\right|_{t=0}
=
\phi(x,y,z),
\qquad
(x,y,z)\in\overline V.
}
\tag{11.27}
$$

边界条件的形式比较多样化，要由具体问题中描述的具体状况决定。总的原则是：边界条件应该完全描写边界上各点在任一时刻 $(t\ge0)$ 的状况。首先以弦的横振动为例。如果弦的“两端固定”，那么边界条件就是

$$
\boxed{
\left.u\right|_{x=0}=0,
\qquad
\left.u\right|_{x=l}=0,
\qquad
t\ge0.
}
\tag{11.28}
$$

对于杆的纵振动，如果 $x=0$ 端固定，则该处的边界条件仍是

$$
\left.u\right|_{x=0}
=
0.
\tag{11.29}
$$

如果杆的另一端 $(x=l)$ 受 $x$ 方向上的外力作用，单位面积上的拉力是 $F(t)$（见图 11.4，约定 $F(t)$ 沿 $x$ 正方向为正），那么这一端的边界条件并不能直接看出。我们应当模仿推导方程的办法，在端点 $x=l$ 处截取一小块介质，长度为 $\varepsilon$。根据 Newton 第二定律可知，这一小段介质所受的合力（外力以及介质的其余部分施加的内应力），应该等于这一小段介质的质量乘以介质的平均加速度，

$$
\rho\varepsilon S
\frac{\overline{\partial^2u}}{\partial t^2}
=
F(t)S
-
P(l-\varepsilon,t)S.
$$

![Figure 11.4](...)

图 11.4　端点所受外力与应力平衡

令 $\varepsilon\to0$，并利用 (11.8) 式，则有

$$
\boxed{
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
\frac1E F(t).
}
\tag{11.30}
$$

如果外力为 0，即 $x=l$ 端自由，则

$$
\boxed{
\left.
\frac{\partial u}{\partial x}
\right|_{x=l}
=
0.
}
\tag{11.31}
$$

如果外力 $F(t)$ 不是一个确定的已知函数，而是由弹簧提供的弹性力（约定沿 $x$ 正方向为正），则

$$
F(t)S
=
-k[u(l,t)-u_0],
$$

$k$ 是弹簧的劲度系数，于是，

$$
\boxed{
\left(
ku
+
ES\frac{\partial u}{\partial x}
\right)_{x=l}
=
ku_0.
}
\tag{11.32}
$$

对于热传导问题，常见的边界条件有下列几种类型。

第一种类型是边界上各点的温度已知，

$$
\boxed{
\left.u\right|_\Sigma
=
\phi(\Sigma,t).
}
\tag{11.33}
$$

这里我们用 $\Sigma$ 代表边界上的变点，同时也表示这些点的坐标。

第二种类型是单位时间内，通过单位面积的边界面散出的热量 $\psi(\Sigma,t)$ 已知。这时可在边界内侧截取一小薄层介质，它的一个底面在介质的表面，另一个底面在介质内部。两底面积相等，厚度趋于 0（见图 11.5）。根据能量守恒定律可知，介质从两个底面及侧面流入的热量之和，应该等于这一块介质温度升高所需要的热量。但是，当介质的厚度趋于 0 时，通过侧面流入的热量应该趋于 0（因为侧面积趋于 0），介质的热容量趋于 0（因为介质的质量趋于 0），因此，由介质内部流向薄层的热量应当全部通过介质表面散出。于是，可以写出边界条件

$$
\boxed{
-k
\left.
\frac{\partial u}{\partial n}
\right|_\Sigma
=
\psi(\Sigma,t).
}
\tag{11.34}
$$

![Figure 11.5](...)

图 11.5　边界面处的热流连续

其中 $\partial/\partial n$ 称为法向微商，它是温度梯度矢量在外法线方向上的投影，

$$
\frac{\partial}{\partial n}
=
\boldsymbol n\cdot\nabla,
\qquad
\text{即}
\qquad
\frac{\partial u}{\partial n}
=
\boldsymbol n\cdot\nabla u.
$$

如果边界绝热，则 $\psi=0$，

$$
\boxed{
\left.
\frac{\partial u}{\partial n}
\right|_\Sigma
=
0.
}
\tag{11.35}
$$

第三种类型是介质通过边界按 Newton 冷却定律散热：单位时间通过单位面积表面向外界散出的热量与介质表面温度 $\left.u\right|_\Sigma$ 和外界温度 $u_0$ 之差成正比。设比例系数为 $H$，则边界条件就是

$$
-k
\left.
\frac{\partial u}{\partial n}
\right|_\Sigma
=
H
\left(
\left.u\right|_\Sigma-u_0
\right),
\tag{11.36}
$$

或者写成

$$
\boxed{
\left(
Hu
+
k\frac{\partial u}{\partial n}
\right)_\Sigma
=
Hu_0.
}
$$

在上面的讨论中出现的边界条件有一个共同的特点：就未知函数而言，它们都是线性的。再进一步细分，还可以分为三类。

第一类边界条件：给出边界上各点的函数值；

第二类边界条件：给出边界上各点函数的法向微商值；

第三类边界条件：给出边界上各点的函数值与法向微商值之间的线性关系。

需要注意，就实际问题而言，在整个边界面上，各点的边界条件并不一定能有统一的表达式，也不见得属于同一种类型。其实上面讨论的弹性杆的边界条件就是如此。

### 无穷远条件

以上讨论的都是有界空间，都存在确定的几何边界。如果要讨论无界空间的问题，这时的边界条件当然就应当给出未知函数在无穷远处的极限行为，例如函数乃至它的微商在无穷远处有界。

### 有界条件

在有界空间的问题中，有时也要出现有界条件。例如，当我们采用平面极坐标系、柱坐标系或球坐标系时，偏微商 $\partial u/\partial r$ 在坐标原点失去意义，因而需要针对具体情况，在坐标原点补充上有界条件或其他条件。以后在有关章节中再做讨论。

### 内部界面上的连接条件

前面说到，微分方程只在空间区域的内部成立。这也是有条件的。如果在区域的内部出现结构（例如密度或其他相关的物理性质）上的跃变，处理这类问题，一种方法是在微分方程中引入间断函数乃至 $\delta$ 函数，另一种方法是认为微分方程在这些跃变点（线、面）上不成立。在这些点（线、面）上，需要补充上相应的条件，通常称为连接条件。

最简单的例子是由两种质料组成的介质，例如由两种不同材料的弹性细绳连接而成的弦。从波动方程的推导过程可以看出，质料不同，两段弦上的振动传播速度不同，方程当然不同：

对于第一段弦

$$
\frac{\partial^2u_1(x,t)}{\partial t^2}
-
a_1^2
\frac{\partial^2u_1(x,t)}{\partial x^2}
=
0,
\qquad
0<x<x_0,
\tag{11.37a}
$$

对于第二段弦

$$
\frac{\partial^2u_2(x,t)}{\partial t^2}
-
a_2^2
\frac{\partial^2u_2(x,t)}{\partial x^2}
=
0,
\qquad
x_0<x<l.
\tag{11.37b}
$$

而且，在连接点根本无法写出微分方程，而应该代之以连接条件。连接条件当然与连接的方式或状况有关。作为最理想的近似，如果这种连接（即跃变）只严格地发生于一点，而且连接得非常牢固光滑，那么，在连接点（设坐标为 $x_0$）的连接条件就是

$$
\left.u_1(x,t)\right|_{x=x_0-0}
=
\left.u_2(x,t)\right|_{x=x_0+0}
\qquad
\text{（即位移相等）},
\tag{11.38a}
$$

$$
\left.
\frac{\partial u_1(x,t)}{\partial x}
\right|_{x=x_0-0}
=
\left.
\frac{\partial u_2(x,t)}{\partial x}
\right|_{x=x_0+0}
\qquad
\text{（即张力相等）}.
\tag{11.38b}
$$

这里的 $\left.u(x,t)\right|_{x=x_0-0}$ 和 $\left.u(x,t)\right|_{x=x_0+0}$ 分别表示函数 $u(x,t)$ 在 $x=x_0$ 点的左、右极限，

$$
\left.u(x,t)\right|_{x=x_0\pm0}
=
\lim_{\varepsilon\to+0}
u(x_0\pm\varepsilon,t).
$$

即使是一种材料构成的介质，也可能出现连接条件。仍以弦的横振动为例。如果在均匀弦的某一点上受到有限大小的力（“集中外力”）$f(t)$ 的作用，方向与横向位移 $u$ 的正方向相反，那么，在这一点（仍记为 $x=x_0$）方程也不成立，也应该代之以连接条件

$$
\left.u(x,t)\right|_{x=x_0-0}
=
\left.u(x,t)\right|_{x=x_0+0}
\qquad
\text{（位移相等）},
\tag{11.39a}
$$

$$
T
\left[
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=x_0+0}
-
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=x_0-0}
\right]
=
f(t)
\qquad
\text{（张力与外力平衡）}.
\tag{11.39b}
$$

更进一步，如果这个集中外力是由一个重物 $M$ 提供的，而且，这个重物和弦同步地发生运动，重物和弦之间没有相对位移，这时，连接条件又应该变为

$$
\left.u(x,t)\right|_{x=x_0-0}
=
\left.u(x,t)\right|_{x=x_0+0},
\tag{11.40a}
$$

$$
T
\left[
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=x_0+0}
-
\left.
\frac{\partial u(x,t)}{\partial x}
\right|_{x=x_0-0}
\right]
=
Mg
+
M
\left.
\frac{\partial^2u}{\partial t^2}
\right|_{x=x_0}.
\tag{11.40b}
$$

在静电场的问题中，也常见到连接条件。例如，在两种电介质的界面 $\Sigma'$ 上，电势连续和电位移矢量的法向分量连续，

$$
\left.u_1\right|_{\Sigma'}
=
\left.u_2\right|_{\Sigma'},
\qquad
\varepsilon_1
\left.
\frac{\partial u_1}{\partial n}
\right|_{\Sigma'}
=
\varepsilon_2
\left.
\frac{\partial u_2}{\partial n}
\right|_{\Sigma'}.
\tag{11.41}
$$

其中 $u_1$ 和 $u_2$ 分别表示电介质 1 和电介质 2 中的电势，$\varepsilon_1$ 和 $\varepsilon_2$ 分别是这两种电介质的电容率。关于这两个条件的推导，见电磁学或电动力学教材，这里从略。

## §11.5 定解问题的适定性

通过以上几节的讨论，我们看到，在处理某些实际的物理问题时，可能会归结为在一定定解条件下求解偏微分方程。偏微分方程加上相应的定解条件就构成定解问题。

首先解释一下什么叫作定解问题的解。为了表述方便，不妨以有界空间内的热传导问题为例。设定解问题为

$$
\frac{\partial u}{\partial t}
-
\kappa\nabla^2u
=
f(x,y,z,t),
\qquad
(x,y,z)\in V,\ t>0,
\tag{11.42a}
$$

$$
\left.u\right|_\Sigma
=
\mu(\Sigma,t),
\qquad
t\ge0,
\tag{11.42b}
$$

$$
\left.u\right|_{t=0}
=
\phi(x,y,z),
\qquad
(x,y,z)\in\overline V,
\tag{11.42c}
$$

其中 $\overline V=V+\Sigma$。假设 $f(x,y,z,t),\mu(\Sigma,t)$ 和 $\phi(x,y,z)$ 均为连续函数，则此定解问题的解 $u(x,y,z,t)$ 应当满足：

(1) 是 $(x,y,z)\in\overline V,\ t>0$ 内的连续函数；

(2) 在 $(x,y,z)\in V,\ t>0$ 内有连续一阶偏导数 $\partial u/\partial t$ 和连续二阶偏导数 $\partial^2u/\partial x^2,\partial^2u/\partial y^2,\partial^2u/\partial z^2$；

(3) 是方程 (11.42a) 的解，即 $u(x,y,z,t)$ 使 (11.42a) 在 $(x,y,z)\in V,\ t>0$ 内为恒等式；

(4) 在边界面 $\Sigma$ 上，使 (11.42b) 式在 $t>0$ 的任一时刻恒成立，并且当 $t\to0+$ 时也成立；

(5) 在初始时刻 $t=0$，使 (11.42c) 式在体积 $V$ 内的任意一点恒成立，并且从体积内部逼近边界面时也成立。

现在的问题是：在什么条件下，定解问题的解是存在的，唯一的，并且是稳定的？

所谓解的存在性，是指在给定的定解条件下，偏微分方程是否有解。如果方程及定解条件中存在矛盾，这种定解问题的解就不可能存在。

所谓解的唯一性，是指在给定的定解条件下，偏微分方程的解是否只有一个。如果已知条件不足以保证实际问题解的唯一性，就说明还需要补充新的定解条件。

所谓解的稳定性，即如果定解问题中的已知条件（例如方程或定解条件中的已知函数）有微小改变时，相应地，解也只有微小的改变。由于在建立方程或者定解问题时，总会要做适当的简化与近似，解的稳定性就成为检验构建定解问题合理性的必然要求。

定解问题解的存在性、唯一性和稳定性，统称适定性。这当然需要从数学上加以研究。从物理上考虑，只要对实际物理问题的抽象是合理的，初始条件的确是完全、确定地描写了初始时刻（通常取为 $t=0$）体系内部以及边界面上任意一点的状况，边界条件的确是完全而且确定地描写了边界面上任意一点在 $t\ge0$ 的状况，那么，这样构成的定解问题就一定是适定的。

与此相关的问题是，按照上面对于初始条件和边界条件的要求，在这些条件中出现的已知函数必须满足一定的连续性要求。例如对于上面的定解问题 (11.42)，就应当有

$$
\left.
\mu(\Sigma,t)
\right|_{t=0}
=
\left.
\phi(x,y,z)
\right|_\Sigma.
\tag{11.43}
$$

有些定解问题不一定满足这个要求。可以设想，把初始温度分布为 $\phi(x,y,z)$ 的一块介质放到一个恒温环境（例如温度恒为 $u_0$）中，使介质表面的温度也迅速达到恒温 $u_0$，如果要求的精度许可，允许忽略介质表面冷却或升温细致过程的影响，那么，就可以简单地将边界条件写成

$$
\left.
u(x,y,z,t)
\right|_\Sigma
=
u_0.
\tag{11.44}
$$

这样做的结果，尽管和精确的边界条件还有差别，但只要这种差别足够小，那么，解的稳定性就告诉我们，由此引起的解的差异也是足够小的。当然，如果我们就是要研究这种冷却或升温细致过程的影响，这种近似就是不可取的。

## 习题

1. 在弦的横振动问题中，若弦受到一与速度成正比（比例系数为 $-\alpha$）的阻尼，试导出弦的有阻尼振动方程。又若除了阻尼力之外，弦还受到与弦的位移成正比（比例系数为 $-k$）的回复力，则此时弦的振动满足的方程是什么？

2. 一长为 $l$、横截面积为 $S$ 的均匀弹性杆，已知一端 $(x=0)$ 固定，另一端 $(x=l)$ 在杆轴方向上受拉力 $F$ 的作用而达到平衡（见图 11.6）。在 $t=0$ 时，撤去外力 $F$。试列出杆的纵振动所满足的方程、边界条件和初始条件。

![Figure 11.6](...)

图 11.6

3. 一长为 $l$ 的金属细杆（可近似地看成是一维的），通有稳定电流 $I$。如果杆的两端 $(x=0$ 和 $x=l)$ 均按 Newton 冷却定律与外界交换热量。外界温度为 $u_0$。初始时杆的温度分布为 $u_0(1-2x/l)^2$。试写出杆上温度场所满足的方程、边界条件和初始条件，设金属杆的电阻为 $R$。

4. 在铀块中，除了中子的扩散运动外，还存在中子的吸收和增殖过程。设在单位时间内、单位体积中吸收和增殖的中子数均正比于该时刻、该处的中子浓度 $u(\boldsymbol r,t)$，因而净增中子数可表为 $\alpha u(\boldsymbol r,t)$，$\alpha$ 为比例常数。试导出 $u(\boldsymbol r,t)$ 所满足的偏微分方程。

5. 有长为 $l$ 的均匀细杆，现通过其两端，在单位时间内经单位面积分别供给热量 $q_1$ 与 $q_2$。试写出相应的边界条件。

6. 有一半径为 $a$、表面涂黑的金属球，曝晒于日光下（见图 11.7），在垂直于光线的单位面积上，单位时间内吸收热量 $M$。同时，球面按 Newton 冷却定律散热（不妨取周围介质的温度为 0）。试在适当的坐标系中写出边界条件。

![Figure 11.7](...)

图 11.7

### 稳定问题中的三类定解问题

对于静电场或稳定的温度场，未知函数（静电势或温度）满足 Laplace 方程或 Poisson 方程。由于与时间无关，所以定解条件就应当只是边界条件，定解问题就称为相关偏微分方程的边值问题。依照边界条件纯粹是第一类（未知函数在边界面的数值已知），第二类（未知函数的法向微商在边界面上的数值已知），或者是第三类（未知函数及其法向微商的线性组合在边界面上的数值已知），定解问题更进一步称为偏微分方程的 Dirichlet 问题，Neumann 问题，或者 Robin 问题。
