# 第三章 复变积分

利用复变积分可以相对简明地证明解析函数的很多重要性质. 本章不刻意追求给出数学上完整的严格证明，一些结论会不加证明直接引用，但是由于这部分内容对于研究解析函数性质很重要，建议感兴趣的同学查阅相关的参考文献.

## §3.1 复变积分

复变积分是复平面 $\mathbb C$ 上的线积分. 设 $C$ 是 $\mathbb C$ 内的一条由 $A$ 点到 $B$ 点的曲线，函数 $f(z)$ 在 $C$ 上有定义. 如图 3.1，把曲线 $C$ 任意分割为 $n$ 段，分点为 $z_0(=A),z_1,z_2,\cdots,z_n(=B)$，$\zeta_k$ 是 $z_{k-1}\to z_k$ 段上的任意一点，作和数

$$
\sum_{k=1}^{n}f(\zeta_k)(z_k-z_{k-1})=\sum_{k=1}^{n}f(\zeta_k)\Delta z_k,
$$

其中 $\Delta z_k=z_k-z_{k-1}$. 若当 $n\to\infty,\max|\Delta z_k|\to0$ 时，此和数的极限存在，且极限值与 $\zeta_k$ 的选取无关，则称此极限值为函数 $f(z)$ 沿曲线 $C$ 的积分，记为

$$
\int_C f(z)\,dz=\lim_{\max|\Delta z_k|\to0}\sum_{k=1}^{n}f(\zeta_k)\Delta z_k. \tag{3.1}
$$

![Figure 3.1](...)

图 3.1

曲线 $C$ 也称为积分路径. 如果积分路径是闭合的，又常称为（积分）围道.

一个复变积分实际上是两个实变线积分的有序组合：

$$
\int_C f(z)\,dz=\int_C(u+iv)(dx+i\,dy)
=\int_C(u\,dx-v\,dy)+i\int_C(v\,dx+u\,dy). \tag{3.2}
$$

因此，根据实函数线积分的知识可以知道，如果 $C$ 是分段光滑曲线，$f(z)$ 是 $C$ 上的连续函数，则复变积分一定存在.

和实变积分相似，复变积分具有下列性质：

(1) 如果积分 $\int_C f_1(z)dz,\int_C f_2(z)dz,\cdots,\int_C f_n(z)dz$ 都存在，其中 $n$ 是正整数，则

$$
\int_C[f_1(z)+f_2(z)+\cdots+f_n(z)]dz
=\int_Cf_1(z)dz+\int_Cf_2(z)dz+\cdots+\int_Cf_n(z)dz;
$$

(2) 若 $C=C_1+C_2+\cdots+C_n$，其中 $n$ 是正整数，则

$$
\int_C f(z)dz=\int_{C_1}f(z)dz+\int_{C_2}f(z)dz+\cdots+\int_{C_n}f(z)dz;
$$

(3)

$$
\int_{C^-}f(z)dz=-\int_Cf(z)dz,
$$

其中 $C^-$ 表示 $C$ 的逆向；

(4)

$$
\int_C af(z)dz=a\int_C f(z)dz,
$$

其中 $a$ 为复常数；

(5)

$$
\left|\int_Cf(z)dz\right|\leq\int_C|f(z)|\,|dz|;
$$

(6)

$$
\left|\int_Cf(z)dz\right|\leq Ml,
$$

其中 $M$ 为 $|f(z)|$ 在 $C$ 上的上界，$l$ 为 $C$ 的长度.

显然，复变积分的值，不仅依赖于被积函数和积分路径的端点，还依赖于积分路径.

**例 3.1** 求 $\int_C\operatorname{Re}z\,dz$，其中 $C$ 为 (i) 沿实轴由 $0\to1$，再平行于虚轴由 $1\to1+i$；(ii) 沿虚轴由 $0\to i$，再平行于实轴由 $i\to1+i$；(iii) 沿直线由 $0\to1+i$.

**解** 对于 (i)，

$$
\int_C\operatorname{Re}z\,dz=\int_0^1x\,dx+\int_0^1i\,dy=\frac12+i;
$$

对于 (ii)，

$$
\int_C\operatorname{Re}z\,dz=\int_0^1x\,dx=\frac12;
$$

对于 (iii)，

$$
\int_C\operatorname{Re}z\,dz=\int_0^1(1+i)t\,dt=\frac12(1+i).
$$

## §3.2 Cauchy 定理

**定理 3.1（Cauchy 定理）** 如果函数 $f(z)$ 在有界闭区域 $\overline G$ 中解析，则沿 $\overline G$ 的边界 $C$（见图 3.2 与图 3.3(a)），有

$$
\boxed{\oint_C f(z)dz=0}. \tag{3.3}
$$

![Figure 3.2](...)

图 3.2 Cauchy 定理（单连通区域）

![Figure 3.3](...)

图 3.3 Cauchy 定理（有界多连通区域）

**证** 证明分两步走. 首先证明定理对于单连通区域成立，然后再推广到多连通区域.

对于单连通区域的情形，为简单起见，下面在更强的条件下证明这个定理. 附加的条件是 $f'(z)$ 在 $\overline G$ 中连续，因而可以对

$$
\oint_Cf(z)dz=\oint_C(u\,dx-v\,dy)+i\oint_C(v\,dx+u\,dy)
$$

应用 Green 公式

$$
\oint_C[P(x,y)dx+Q(x,y)dy]
=\iint_S\left(\frac{\partial Q}{\partial x}-\frac{\partial P}{\partial y}\right)dx\,dy,
$$

从而将闭合围道积分化为面积分：

$$
\oint_C(u\,dx-v\,dy)
=-\iint_S\left(\frac{\partial v}{\partial x}+\frac{\partial u}{\partial y}\right)dx\,dy,
$$

$$
\oint_C(v\,dx+u\,dy)
=\iint_S\left(\frac{\partial u}{\partial x}-\frac{\partial v}{\partial y}\right)dx\,dy.
$$

由 Cauchy-Riemann 方程，右端两个积分中的被积函数均为 0，故 (3.3) 式得证.

如果被积函数在区域内有奇点，为了保证被积函数的解析性就需要把奇点排除在外，这样 $f(z)$ 就是定义在多连通闭区域 $\overline G$ 中的解析函数. 设（有界）多连通闭区域 $\overline G$ 的边界是由分段光滑曲线 $C_0,C_1,C_2,\cdots,C_n$ 组成，并且 $C_1,C_2,\cdots,C_n$ 都包含在 $C_0$ 的内部. 于是，可以作适当的割线把 $C_0,C_1,C_2,\cdots,C_n$ 和 $C_0$ 连接起来（如图 3.3(b)），从而得到一个有界单连通闭区域 $\overline{G'}$，$f(z)$ 在 $\overline{G'}$ 中是解析的，因而可以应用有界单连通区域的 Cauchy 定理，得

$$
\oint_{C_0}f(z)dz+\int_{a_1}^{b_1}f(z)dz+\oint_{C_1}f(z)dz+\int_{b_1}^{a_1}f(z)dz
+\int_{a_2}^{b_2}f(z)dz+\oint_{C_2}f(z)dz
$$

$$
+\int_{b_2}^{a_2}f(z)dz+\cdots+\int_{a_n}^{b_n}f(z)dz+\oint_{C_n}f(z)dz+\int_{b_n}^{a_n}f(z)dz=0,
$$

这里积分路径均为沿边界的正向，即沿 $C_0$ 段为逆时针方向，而沿 $C_1,C_2,\cdots,C_n$ 等段则为顺时针方向. 由于 $f(z)$ 在 $\overline{G'}$ 中单值，故沿同一割线两岸的积分值互相抵消，

$$
\int_{a_i}^{b_i}f(z)dz+\int_{b_i}^{a_i}f(z)dz=0,
$$

所以

$$
\oint_{C_0}f(z)dz+\sum_{i=1}^n\oint_{C_i}f(z)dz=0.
$$

这样也就证明了 (3.3) 式. 注意这时的积分路径 $C$ 应当是沿多连通区域的边界（正向）.

就多连通区域而言，有时还将 Cauchy 定理改写成

$$
\oint_{C_0}f(z)dz=\sum_{i=1}^n\oint_{C_i^{(-)}}f(z)dz, \tag{3.3'}
$$

只是需要注意，这里 $C_0$ 与 $C_i^{(-)},\ i=1,2,3,\cdots,n$ 须保持为同向，即同为逆时针方向，或同为顺时针方向.

定理 3.1 的成立条件还可以减弱为：若 $f(z)$ 在有界单连通区域 $G$ 内解析，在闭区域 $\overline G$ 中连续. 证明从略.

**思考题** 对于任一解析函数的实部或虚部，Cauchy 定理仍成立吗？如果成立，试证明之. 如果不成立，试说明理由，并举一例.

还值得指出，以上讨论的，无论是单连通区域，或是多连通区域，均只限于有界区域，不能包含 $\infty$ 点在内. 从下面例 3.3 即可看到，即使 $f(z)$ 在 $\infty$ 点解析，例如 $f(z)=1/z$，它绕 $\infty$ 点一周的积分也可以不为 0.

Cauchy 定理还可以表述为复变积分的**变形定理**：若函数 $f(z)$ 在区域 $G$ 内解析，$C$ 为 $G$ 内的简单闭合曲线，如果能将曲线 $C$ 在 $G$ 内连续地变形为曲线 $C'$①，则必有

$$
\oint_Cf(z)dz=\oint_{C'}f(z)dz.
$$

Cauchy 定理，顾名思义，最早由 Cauchy 给出并证明. 但是当初 Cauchy 给出的解析函数定义是函数 $f(z)$ 在区域 $G$ 内每一点均可导且 $f'(z)$ 连续. 后来 Goursat 论证了，不需要“导函数连续”这个条件即可证明 Cauchy 定理. 从此解析函数的定义修改为区域 $G$ 内每一点都可导的函数. 所以，上面给出的 Cauchy 定理的证明，是解析函数理论发展历史中的产物. 现在，一些文献中也称 Cauchy 定理为 Cauchy-Goursat 定理. 下面我们不加证明地介绍几个引理，用 Goursat 的方法证明 Cauchy 定理②.

**引理 3.1** 设 $\overline{K_1},\overline{K_2},\overline{K_3},\cdots$ 是复平面 $\mathbb C$ 内的一系列有界闭区域，并且满足

$$
\overline{K_1}\supset\overline{K_2}\supset\overline{K_3}\supset\cdots,
$$

则至少存在一个点 $x_0$，满足

$$
x_0\in\bigcap_i^{\infty}\overline{K_i}.
$$

**引理 3.2** 设 $\alpha$ 和 $\beta$ 是复值常数，$\gamma$ 是复平面 $\mathbb C$ 内任意一个闭合的三角形围道，那么

$$
\oint_\gamma(\alpha+\beta z)dz=0.
$$

**引理 3.3** 设 $G$ 是复平面 $\mathbb C$ 内的一个有界单连通区域，$f$ 是定义在 $G$ 上的一个连续函数，$\gamma$ 是 $G$ 内的一条可求长曲线，则 $\forall\varepsilon>0$，存在内接于 $\gamma$ 且完全位于 $G$ 内的折线 $\Gamma$，使得

$$
\left|\int_\gamma fdz-\int_\Gamma fdz\right|<\varepsilon.
$$

根据引理 3.3，只需证明对于内接于 $C$ 的任何闭折线 $P$，均有

$$
\oint_P f(z)dz=0
$$

即可. 现在将 $P$ 所围成的多边形加以分割. 如果它不是凸多边形的话，首先总可以将它分割为若干个凸多边形，函数沿这些凸多边形的公共边上的积分互相抵消（因为积分路径方向相反），沿闭折线 $P$ 的积分就等于沿各凸多边形积分之和，因此我们只要证明函数沿凸多边形的积分为 0 即可. 而对于凸多边形，又可以从选定的一点出发，与多边形的各个顶点作连线，把它分解成若干个三角形区域. 于是，沿凸多边形的积分就是沿各个三角形积分之和. 这样，问题就又归结为只要证明函数沿任意一个三角形 $T$ 的积分 $\oint_Tf(z)dz=0$ 即可. 记

$$
\left|\oint_Tf(z)dz\right|=M.
$$

下面证明 $M=0$.

① 准确地说，$C$ 与 $C'$ 应当是同伦曲线.

② 严格的证明见参考书目 [22] 的第四章，第 2 节；或 L. A. Ahlfors, *Complex Analysis* (Third Edition), McGraw-Hill International Editions, 2004，第四章.

为此，把三角形 $T$ 三条边的中点连接起来，就得到 4 个与 $T$ 相似的、小的全等三角形 $T_1,T_2,T_3,T_4$（见图 3.4）. 这 4 个三角形相邻边的积分相互抵消，因此有

$$
M=\left|\oint_Tf(z)dz\right|
=\left|\sum_{i=1}^{4}\oint_{T_i}f(z)dz\right|
\leq\sum_{i=1}^{4}\left|\oint_{T_i}f(z)dz\right|,
$$

所以右端的 4 个积分中至少有一个积分的模不小于 $M/4$. 不妨设这个积分来自小三角形 $T^{(1)}$，即

$$
\left|\oint_{T^{(1)}}f(z)dz\right|\geq\frac M4.
$$

![Figure 3.4](...)

图 3.4 三角形闭合围道

再将这个小三角形分割为 4 个更小的三角形，重复上面的讨论，我们又能得到

$$
\left|\oint_{T^{(2)}}f(z)dz\right|\geq\frac M{4^2}.
$$

如此继续，将三角形 $T$ 分割 $n$ 次后，就得到

$$
\left|\oint_{T^{(n)}}f(z)dz\right|\geq\frac M{4^n}.
$$

若三角形 $T$ 的周长为 $l$，最大边长为 $d$，则 $T^{(n)}$ 的周长为 $l/2^n$，最大边长为 $d/2^n$. 将 $T^{(n)}$ 所包围的闭区域记为 $\Delta_n$，则应当有

$$
\Delta_{i+1}\subset\Delta_i,\qquad i=0,1,2,\cdots.
$$

按照上面的引理 3.1，就必然存在唯一的一点 $z_0$，$z_0\in\bigcap_i^\infty\Delta_i$. 因为 $z_0\in G$，$f(z)$ 在 $z_0$ 点可导，故在 $z_0$ 点的邻域 $V(z_0;\delta):|z-z_0|<\delta$ 内有

$$
f(z)-f(z_0)=f'(z_0)(z-z_0)+\rho(z,z_0)(z-z_0),
$$

其中

$$
\lim_{z\to z_0}\rho(z,z_0)=0.
$$

在此基础上，我们就能来考察 $\oint_{T^{(n)}}f(z)dz$. 因为当 $n$ 足够大时，$\Delta_n\subset V(z_0;\delta)$，所以

$$
\oint_{T^{(n)}}f(z)dz=f(z_0)\oint_{T^{(n)}}dz+f'(z_0)\oint_{T^{(n)}}(z-z_0)dz
+\oint_{T^{(n)}}\rho(z,z_0)(z-z_0)dz.
$$

由引理 3.2，

$$
\oint_{T^{(n)}}dz=0,\qquad \oint_{T^{(n)}}z\,dz=0,
$$

所以

$$
\oint_{T^{(n)}}f(z)dz=\oint_{T^{(n)}}\rho(z,z_0)(z-z_0)dz.
$$

这样就立即得到

$$
\left|\oint_{T^{(n)}}f(z)dz\right|
\leq\oint_{T^{(n)}}|\rho(z,z_0)(z-z_0)dz|
\leq\frac d{2^n}\cdot\frac l{2^n}\max_{z\in T^{(n)}}|\rho(z,z_0)|.
$$

上面已经证明 $M\leq4^n\left|\oint_{T^{(n)}}f(z)dz\right|$，因此

$$
M\leq dl\cdot\max_{z\in T^{(n)}}|\rho(z,z_0)|.
$$

因为 $\lim_{z\to z_0}\rho(z,z_0)=0$，亦即

$$
\lim_{n\to\infty}\max_{z\in T^{(n)}}|\rho(z,z_0)|=0,
$$

所以 $M=0$. 证毕.

**例 3.2** 计算 $\oint_Cz^n dz$ 值，其中 $n$ 为整数，$C$ 为复平面 $\mathbb C$ 内一条简单闭合围道.

**解** 当 $n$ 为自然数时，显然，按照（有界单连通区域的）Cauchy 定理

$$
\oint_Cz^n dz=0.
$$

当 $n$ 为负整数时，如果 $C$ 沿逆时针方向时所包含的区域内（简称 $C$ 内）不含 $z=0$，则也有

$$
\oint_Cz^n dz=0.
$$

如果 $C$ 内含有 $z=0$ 点，则按（有界多连通区域的）Cauchy 定理，有

$$
\oint_Cz^n dz=\oint_{|z|=\varepsilon}z^n dz
=\int_0^{2\pi}\varepsilon^{n+1}e^{i(n+1)\theta}i\,d\theta
=\begin{cases}
2\pi i,& n=-1,\\
0,& n=-2,-3,-4,\cdots.
\end{cases}
$$

上式中的 $\varepsilon$ 足够小，以保证圆 $|z|=\varepsilon$ 位于 $C$ 所包含的区域内.

总结上面的结果，就有

$$
\oint_Cz^n dz=
\begin{cases}
2\pi i,& n=-1,\text{且 }C\text{ 内含有 }z=0;\\
0,& \text{其他情形}.
\end{cases} \tag{3.4}
$$

或者，更一般地，

$$
\oint_C(z-a)^n dz=
\begin{cases}
2\pi i,& n=-1,\text{且 }C\text{ 内含有 }z=a;\\
0,& \text{其他情形}.
\end{cases} \tag{3.5}
$$

**例 3.3** 计算积分 $\int_0^{2\pi}\sin^{2n}\theta\,d\theta$.

**解** 令 $z=e^{i\theta}$，则

$$
\sin\theta=\frac{z^2-1}{2iz},\qquad d\theta=\frac{dz}{iz},
$$

而原积分即化为 $z$ 平面上沿单位圆的积分，

$$
\int_0^{2\pi}\sin^{2n}\theta\,d\theta
=\oint_{|z|=1}\left(\frac{z^2-1}{2iz}\right)^{2n}\frac{dz}{iz}
=\frac{(-1)^n}{2^{2n}i}\oint_{|z|=1}(z^2-1)^{2n}\frac{dz}{z^{2n+1}}.
$$

将 $(z^2-1)^{2n}$ 作二项式展开，

$$
(z^2-1)^{2n}=\sum_{l=0}^{2n}\frac{(2n)!}{l!(2n-l)!}(-1)^{2n-l}z^{2l},
$$

则其中只有 $z^{2n}$ 项（即 $l=n$ 项）对积分有贡献（为什么？），因此，

$$
\int_0^{2\pi}\sin^{2n}\theta\,d\theta
=2\pi i\times\frac1{2^{2n}i}\frac{(2n)!}{n!n!}
=\frac{\pi}{2^{2n-1}}\frac{(2n)!}{n!n!}.
$$

**例 3.4** 证明实系数多项式必有零点，除非该多项式恒为常数（零次多项式）.

**证** 用反证法. 设此多项式为

$$
P_n(z)=a_nz^n+a_{n-1}z^{n-1}+\cdots+a_0,\qquad a_n\neq0,
$$

$a_n,a_{n-1},\cdots,a_0$ 均为实数. 按所设，若 $P_n(z)$ 无零点，则当 $z$ 取实数值 $x$ 时，$P(x)$ 一定不为 0，且不变号，因此

$$
\int_0^{2\pi}\frac{d\theta}{P_n(2\cos\theta)}\neq0.
$$

令 $z=e^{i\theta}$，则上式化为

$$
\oint_{|z|=1}\frac1{P_n(z+z^{-1})}\frac{dz}{iz}
=\frac1i\oint_{|z|=1}\frac{z^{n-1}}{Q_{2n}(z)}dz\neq0, \tag{3.6}
$$

其中 $Q_{2n}(z)=z^nP_n(z+z^{-1})$ 为 $2n$ 次多项式. 因为 $P_n(z)\neq0$，故当 $z\neq0$ 时 $P_n(z+z^{-1})\neq0$，即 $Q_{2n}(z)\neq0$；又，$Q_{2n}(0)=a_n$ 亦不为 0，因此函数 $z^{n-1}/Q_{2n}(z)$（在全平面）解析. 根据 Cauchy 定理，一定有

$$
\oint_{|z|=1}\frac{z^{n-1}}{Q_{2n}(z)}dz=0,
$$

与 (3.6) 矛盾. 因此 $P_n(z)$ 必有零点. 命题得证.

**讨论** 此结论可推广到复系数多项式. 因为即使

$$
P_n(z)=a_nz^n+a_{n-1}z^{n-1}+\cdots+a_0
$$

为复系数多项式，也可以构造一个新的实系数多项式

$$
(a_nz^n+a_{n-1}z^{n-1}+\cdots+a_0)
(a_n^*z^n+a_{n-1}^*z^{n-1}+\cdots+a_0^*),
$$

因此上述结论仍然成立.

由 Cauchy 定理立即可得到下面的推论（其实它就是复变积分变形定理的另一种表述形式）：

**推论** 若 $f(z)$ 在有界单连通区域 $G$ 中解析，则复变积分 $\int_Cf(z)dz$ 与路径 $C$ 无关，其中 $C\subset G$.

既然在有界单连通区域中解析函数的积分值与路径无关，因此，如果固定起点 $z_0$，而令终点 $z$ 为变点，则作为积分上限的函数

$$
\int_{z_0}^{z}f(\zeta)d\zeta=F(z),\qquad z\in G \tag{3.7}
$$

是有界单连通区域 $G$ 内的单值函数，称为 $f(z)$ 的不定积分.

**定理 3.2** 如果函数 $f(z)$ 在有界单连通区域 $G$ 内解析，则 $f(z)$ 的不定积分

$$
F(z)=\int_{z_0}^{z}f(\zeta)d\zeta,\qquad z\in G \tag{3.8}
$$

也在 $G$ 内解析，并且

$$
F'(z)=\frac{d}{dz}\int_{z_0}^{z}f(\zeta)d\zeta=f(z),\qquad z\in G. \tag{3.9}
$$

**证** 只要直接求出 $F(z)$ 的导数即可. 为此，设 $z$ 是 $G$ 内一点，$z+\Delta z$ 是它的邻点，如图 3.5 所示，则

$$
F(z)=\int_{z_0}^{z}f(\zeta)d\zeta,\qquad
F(z+\Delta z)=\int_{z_0}^{z+\Delta z}f(\zeta)d\zeta.
$$

![Figure 3.5](...)

图 3.5 不定积分解析性的证明

因为积分与路径无关，所以

$$
\frac{\Delta F}{\Delta z}=\frac{F(z+\Delta z)-F(z)}{\Delta z}
=\frac1{\Delta z}\int_z^{z+\Delta z}f(\zeta)d\zeta.
$$

由此可得

$$
\left|\frac{\Delta F}{\Delta z}-f(z)\right|
=\left|\frac1{\Delta z}\int_z^{z+\Delta z}f(\zeta)d\zeta-f(z)\right|
=\left|\frac1{\Delta z}\int_z^{z+\Delta z}[f(\zeta)-f(z)]d\zeta\right|
$$

$$
\leq\frac1{|\Delta z|}\int_z^{z+\Delta z}|f(\zeta)-f(z)|\,|d\zeta|.
$$

由于 $f(z)$ 是连续的，$\forall\varepsilon>0,\exists\delta>0$，使当 $|\zeta-z|<\delta$ 时，$|f(\zeta)-f(z)|<\varepsilon$，所以

$$
\left|\frac{\Delta F}{\Delta z}-f(z)\right|
\leq\frac1{|\Delta z|}\cdot\varepsilon\cdot|\Delta z|=\varepsilon,
$$

即得

$$
F'(z)=\lim_{\Delta z\to0}\frac{\Delta F}{\Delta z}=f(z),\qquad z\in G.
$$

这样就证明了 $F(z)$ 在 $G$ 内可导，并且 $F'(z)=f(z)$.

在复变函数中，除了不定积分外，也还有原函数的概念.

如果函数 $\Phi(z)$ 的导数 $\Phi'(z)=f(z)$，则称 $\Phi(z)$ 为 $f(z)$ 的原函数. 上面 (3.7) 式定义的 $f(z)$ 的不定积分就是 $f(z)$ 的一个原函数. 对于给定的函数 $f(z)$，其原函数并不唯一，任意两个原函数之间相差一个常数. 这是因为，如果 $\Phi_1(z)$ 与 $\Phi_2(z)$ 都是 $f(z)$ 的原函数，则

$$
\Phi_1'(z)=f(z),\qquad \Phi_2'(z)=f(z).
$$

所以，$[\Phi_1(z)-\Phi_2(z)]'=0$，即

$$
\Phi_1(z)-\Phi_2(z)=C. \tag{3.10}
$$

知道了被积函数的原函数，可使复变积分的计算大为简化. 例如，设 $\Phi(z)$ 为 $f(z)$ 的一个原函数，则 $f(z)$ 的不定积分

$$
F(z)=\int_{z_0}^{z}f(\zeta)d\zeta=\Phi(z)+C.
$$

但是，显然有 $F(z_0)=\Phi(z_0)+C=0,C=-\Phi(z_0)$. 所以

$$
\int_{z_0}^{z}f(\zeta)d\zeta=\Phi(z)-\Phi(z_0). \tag{3.11}
$$

**例 3.5** 计算积分 $\int_a^b z^n dz$，$n$ 为整数.

**解** 当 $n$ 为自然数时，$z^n$ 在 $\mathbb C$ 内解析，$\frac1{n+1}z^{n+1}$ 是它的一个原函数. 因此，对于 $\mathbb C$ 内的任意一条积分路径，均有

$$
\int_a^b z^n dz=\frac1{n+1}(b^{n+1}-a^{n+1}). \tag{3.12}
$$

当 $n=-2,-3,-4,\cdots$ 时，$z^n$ 在 $\overline{\mathbb C}\setminus\{0\}$ 内解析，其原函数仍可取为 $\frac1{n+1}z^{n+1}$. 因此，在不包含 $z=0$ 的任一有界区域内，仍有

$$
\int_a^b z^n dz=\frac1{n+1}(b^{n+1}-a^{n+1}),\qquad n=-2,-3,-4,\cdots. \tag{3.12'}
$$

当 $n=-1$ 时，$z^{-1}$ 也是在 $\overline{\mathbb C}\setminus\{0\}$ 内解析，但原函数为 $\ln z$. 故在不包含 $z=0$ 的任一有界单连通区域内，

$$
\int_a^b\frac{dz}{z}=\ln b-\ln a. \tag{3.13}
$$

这里需要特别注意，一旦积分的上、下限给定，则在一个有界单连通区域内，上面的积分值一定与路径无关，但对于不同的有界单连通区域，可能会给出不同的积分值. 从计算过程看，(3.13) 式的原函数是多值函数，因此积分值与由 $a$ 变化到 $b$ 的方式有关. 当限制在不含 $z=0$ 的一个有界单连通区域内时，原函数的分支点 $z=0$ 和 $z=\infty$ 都不在区域内，就是把 $\ln z$ 限制在某一个单值分支内，故积分值 $\ln b-\ln a$ 唯一确定. 而对于不同的有界单连通区域，就可能对应于 $\ln z$ 的不同单值分支，因而积分值也就可能不同. 表面上看，积分的上、下限给定，即积分路径的终点与起点给定，但是实际上，它们在多叶 Riemann 面上的位置并不见得相同.

## §3.3 两个有用的引理

**引理 3.4（小圆弧引理）** 如果函数 $f(z)$ 在 $z=a$ 点的空心邻域内连续，并且在 $\theta_1\leq\arg(z-a)\leq\theta_2$ 中，当 $|z-a|\to0$ 时，$(z-a)f(z)$ 一致地趋近于 $k$，则

$$
\boxed{\lim_{\delta\to0}\int_{C_\delta}f(z)dz=ik(\theta_2-\theta_1)}. \tag{3.14}
$$

其中 $C_\delta$ 是以 $z=a$ 为圆心、$\delta$ 为半径、张角为 $\theta_2-\theta_1$ 的圆弧，$|z-a|=\delta,\theta_1\leq\arg(z-a)\leq\theta_2$，见图 3.6.

**证** 因为

$$
\int_{C_\delta}\frac{dz}{z-a}=i(\theta_2-\theta_1),
$$

所以

$$
\left|\int_{C_\delta}f(z)dz-ik(\theta_2-\theta_1)\right|
=\left|\int_{C_\delta}\left[f(z)-\frac{k}{z-a}\right]dz\right|
\leq\int_{C_\delta}|(z-a)f(z)-k|\frac{|dz|}{|z-a|}.
$$

由于在 $\theta_1\leq\arg(z-a)\leq\theta_2$ 中，当 $|z-a|\to0$ 时，$(z-a)f(z)$ 一致地趋近于 $k$，这就意味着 $\forall\varepsilon>0,\exists$（与 $\arg(z-a)$ 无关的）$r(\varepsilon)>0$，使当 $|z-a|=\delta<r$ 时，$|(z-a)f(z)-k|<\varepsilon$. 所以

$$
\left|\int_{C_\delta}f(z)dz-ik(\theta_2-\theta_1)\right|\leq\varepsilon(\theta_2-\theta_1),
$$

即

$$
\lim_{\delta\to0}\int_{C_\delta}f(z)dz=ik(\theta_2-\theta_1).
$$

![Figure 3.6](...)

图 3.6

**引理 3.5（大圆弧引理）** 设 $f(z)$ 在 $\infty$ 点的邻域内连续，在 $\theta_1\leq\arg z\leq\theta_2$ 中，当 $|z|\to\infty$ 时，$zf(z)$ 一致地趋近于 $K$，则

$$
\boxed{\lim_{R\to\infty}\int_{C_R}f(z)dz=iK(\theta_2-\theta_1)}. \tag{3.15}
$$

其中 $C_R$ 是以原点为心、$R$ 为半径、张角为 $\theta_2-\theta_1$ 的圆弧，$|z|=R,\theta_1\leq\arg z\leq\theta_2$，见图 3.7.

![Figure 3.7](...)

图 3.7

**证** 证明和引理 3.4 相仿. 因为

$$
\int_{C_R}\frac{dz}{z}=i(\theta_2-\theta_1),
$$

所以

$$
\left|\int_{C_R}f(z)dz-iK(\theta_2-\theta_1)\right|
=\left|\int_{C_R}\left[f(z)-\frac Kz\right]dz\right|
\leq\int_{C_R}|zf(z)-K|\frac{|dz|}{|z|}.
$$

由于在 $\theta_1\leq\arg z\leq\theta_2$ 中，当 $|z|\to\infty$ 时，$zf(z)$ 一致地趋近于 $K$，这意味着 $\forall\varepsilon>0,\exists$（与 $\arg z$ 无关的）$M(\varepsilon)>0$，使当 $|z|=R>M$ 时，$|zf(z)-K|<\varepsilon$. 所以

$$
\left|\int_{C_R}f(z)dz-iK(\theta_2-\theta_1)\right|\leq\varepsilon(\theta_2-\theta_1),
$$

即

$$
\lim_{R\to\infty}\int_{C_R}f(z)dz=iK(\theta_2-\theta_1).
$$

## §3.4 Cauchy 积分公式

Cauchy 定理从一个侧面反映了解析函数的基本特性：解析函数在其解析区域内各点的函数值是密切相关的. Cauchy-Riemann 方程可以说是这种关联的微分形式，而 Cauchy 定理则是它的积分形式. 这种关联性在 Cauchy 积分公式中也清楚地表现出来.

**有界区域的 Cauchy 积分公式** 设 $f(z)$ 是有界闭区域 $\overline G$ 中的单值解析函数，$\overline G$ 的边界 $C$ 是分段光滑曲线，$a$ 为 $G$ 内一点，则

$$
\boxed{f(a)=\frac1{2\pi i}\oint_C\frac{f(z)}{z-a}dz,\qquad a\in G}. \tag{3.16}
$$

其中积分路径沿 $C$ 的正向.

![Figure 3.8](...)

图 3.8 有界区域的 Cauchy 积分公式

**证** 在有界区域 $G$ 内作小圆 $|z-a|<r$（见图 3.8，$r$ 足够小以保证圆周 $|z-a|=r$ 位于 $G$ 内），则根据（有界多连通区域的）Cauchy 定理，有

$$
\oint_C\frac{f(z)}{z-a}dz
=\oint_{|z-a|=r}\frac{f(z)}{z-a}dz, \tag{3.17}
$$

其中 $C$ 沿区域边界正方向. 令 $r\to0$，因为

$$
\lim_{z\to a}(z-a)\frac{f(z)}{z-a}=f(a),
$$

(3.17) 式等号左边与 $r$ 的大小无关，故可令 $r\to0$，由上节的引理 3.4，就证得

$$
\frac1{2\pi i}\oint_C\frac{f(z)}{z-a}dz=f(a).
$$

作为 Cauchy 积分公式的特殊形式，取 $C$ 为以 $a$ 为圆心、$R$ 为半径的圆周，如果 $f(z)$ 在圆内解析，即可得到

$$
\boxed{f(a)=\frac1{2\pi}\int_0^{2\pi}f(a+Re^{i\theta})d\theta}. \tag{3.18}
$$

这个结果称为**均值定理**：解析函数 $f(z)$ 在解析区域 $G$ 内任意一点 $a$ 的函数值 $f(a)$，等于（完全位于 $G$ 内的）以该点为圆心的任一圆周上函数值的平均.

对于无界区域，需要假设 $f(z)$ 在简单闭合围道 $C$ 上及 $C$ 外（包括 $\infty$ 点）单值解析. 类似地计算 $\frac1{2\pi i}\oint_C\frac{f(z)}{z-a}dz$，其中 $a$ 为 $C$ 外一点，积分路径 $C$ 的走向是绕无穷远点的正向，即顺时针方向，如图 3.9. 在 $C$ 外再作一个以原点为圆心、$R$ 为半径的圆 $C_R$，对于 $C$ 和 $C_R$ 所包围的有界多连通区域，根据有界区域的 Cauchy 积分公式，就有

$$
\frac1{2\pi i}\left[\oint_{C_R}\frac{f(z)}{z-a}dz+\oint_C\frac{f(z)}{z-a}dz\right]=f(a),
$$

$C_R$ 的走向是逆时针方向. 等号两边令 $R\to\infty$，若

$$
\lim_{z\to\infty}z\frac{f(z)}{z-a}=\lim_{z\to\infty}f(z)=K, \tag{3.19}
$$

根据上一节的引理 3.5，即得

$$
\lim_{R\to\infty}\left[\frac1{2\pi i}\oint_{C_R}\frac{f(z)}{z-a}dz\right]=K,
$$

因此

$$
\frac1{2\pi i}\oint_C\frac{f(z)}{z-a}dz=f(a)-K. \tag{3.20}
$$

![Figure 3.9](...)

图 3.9 无界区域的 Cauchy 积分公式

当 $K=0$ 时，就得到了**无界区域的 Cauchy 积分公式**：如果 $f(z)$ 在简单闭合围道 $C$ 上及 $C$ 外解析，且当 $z\to\infty$ 时，$f(z)$ 趋于 0，则 Cauchy 积分公式

$$
f(a)=\frac1{2\pi i}\oint_C\frac{f(z)}{z-a}dz \tag{3.21}
$$

仍然成立，此处 $a$ 是 $C$ 为顺时针方向时所包含的区域内的一点.

## §3.5 解析函数的高阶导数

从 Cauchy 积分公式，可以推断出一个重要结论：如果 $f(z)$ 在有界闭区域 $\overline G$ 中解析，则在 $G$ 内 $f(z)$ 的任何阶导数 $f^{(n)}(z)$ 均存在，并且

$$
\boxed{f^{(n)}(z)=\frac{n!}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^{n+1}}d\zeta,\qquad z\in G}. \tag{3.22}
$$

其中 $C$ 是 $\overline G$ 的正向边界，$z$ 为 $G$ 内任意一点，如图 3.10 所示.

![Figure 3.10](...)

图 3.10 高阶导数公式

**证** 用数学归纳法证明. 首先求 $f'(z)$. 因为

$$
\frac{f(z+h)-f(z)}h
=\frac1{2\pi i}\frac1h\oint_C\left[\frac{f(\zeta)}{\zeta-z-h}-\frac{f(\zeta)}{\zeta-z}\right]d\zeta
$$

$$
=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z-h)(\zeta-z)}d\zeta,
$$

取极限 $h\to0$，左端即为 $f'(z)$，而右端被积函数的极限为 $f(\zeta)/(\zeta-z)^2$. 为了证明在积分号下求极限合法，不妨考察

$$
\oint_C\frac{f(\zeta)d\zeta}{(\zeta-z-h)(\zeta-z)}-\oint_C\frac{f(\zeta)d\zeta}{(\zeta-z)^2}
=h\oint_C\frac{f(\zeta)d\zeta}{(\zeta-z-h)(\zeta-z)^2}.
$$

由于 $f(\zeta)$ 在 $C$ 上连续，故在 $C$ 上有 $|f(\zeta)|\leq M$，

$$
\left|\oint_C\frac{f(\zeta)}{(\zeta-z-h)(\zeta-z)}d\zeta-\oint_C\frac{f(\zeta)}{(\zeta-z)^2}d\zeta\right|
\leq h\cdot\frac{Ml}{\delta^2(\delta-h)}\to0,
$$

其中 $\delta$ 为 $z$ 到 $C$ 的最短距离，$l$ 为 $C$ 的长度. 因此，

$$
f'(z)=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^2}d\zeta.
$$

现在假设

$$
f^{(k)}(z)=\frac{k!}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^{k+1}}d\zeta
$$

成立，于是

$$
\frac{f^{(k)}(z+h)-f^{(k)}(z)}h
=\frac{k!}{2\pi i}\frac1h\oint_C\left[\frac{f(\zeta)}{(\zeta-z-h)^{k+1}}-\frac{f(\zeta)}{(\zeta-z)^{k+1}}\right]d\zeta
$$

$$
=\frac{k!}{2\pi i}\oint_Cf(\zeta)\frac{(\zeta-z)^{k+1}-(\zeta-z-h)^{k+1}}{(\zeta-z)^{k+1}(\zeta-z-h)^{k+1}}d\zeta
$$

$$
=\frac{k!}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^{k+1}(\zeta-z-h)^{k+1}}
\left[(k+1)(\zeta-z-h)^k+O(h)\right]d\zeta,
$$

亦即

$$
\left|\frac{f^{(k)}(z+h)-f^{(k)}(z)}h-rac{(k+1)!}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^{k+2}}d\zeta\right|=O(h),
$$

由此即可证得

$$
f^{(k+1)}(z)=\frac{(k+1)!}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^{k+2}}d\zeta
$$

成立.

这个结果说明，一个复变函数，在其解析区域内任何阶导数都存在，并且都是这个区域内的解析函数. 而在实变函数中并非如此，并不能由 $f'(x)$ 的存在推断出 $f''(x)$ 的存在.

## §3.6 Cauchy 型积分和含参量积分的解析性

上一节关于解析函数高阶导数公式的证明过程中，$f(z)$ 的解析性只是用在：(1) $f(z)$ 可用 Cauchy 积分公式表示；(2) $f(z)$ 在 $C$ 上连续. 因此，重复上面的步骤就可以证明：在一条分段光滑的（闭合或不闭合）曲线 $C$ 上连续的函数 $\phi(\zeta)$ 所构成的积分

$$
f(z)=\frac1{2\pi i}\int_C\frac{\phi(\zeta)}{\zeta-z}d\zeta,\qquad z\notin C \tag{3.23}
$$

（称为 **Cauchy 型积分**）是曲线外点 $z$ 的解析函数，并且 $f(z)$ 的任意 $p$（$p$ 为自然数）阶导数 $f^{(p)}(z)$ 可通过积分号下求导而得到，

$$
f^{(p)}(z)=\frac{p!}{2\pi i}\int_C\frac{\phi(\zeta)}{(\zeta-z)^{p+1}}d\zeta,\qquad z\notin C. \tag{3.24}
$$

**例 3.6** 计算积分

$$
f(z)=\frac1{2\pi i}\oint_{|\zeta|=1}\frac{\zeta^*}{\zeta-z}d\zeta,
$$

其中 $|z|\neq1$.

**解** 这是一个 Cauchy 型积分. 因为在 $|\zeta|=1$ 上 $\zeta^*=1/\zeta$，故

$$
f(z)=\frac1{2\pi i}\oint_{|\zeta|=1}\frac1{\zeta(\zeta-z)}d\zeta.
$$

当 $|z|>1$ 时，$1/(\zeta-z)$ 在 $|\zeta|=1$ 内解析，因此可以用 Cauchy 积分公式计算此积分，

$$
f(z)=\frac1{2\pi i}\oint_{|\zeta|=1}\frac1\zeta\left(\frac1{\zeta-z}\right)d\zeta
=\left.\frac1{\zeta-z}\right|_{\zeta=0}=-\frac1z.
$$

当 $0<|z|<1$ 时，

$$
f(z)=\frac1{2\pi i}\oint_{|\zeta|=1}\frac1z\left(\frac1{\zeta-z}-\frac1\zeta\right)d\zeta=0.
$$

可以看出，此结果对于 $z=0$ 仍成立. 综合以上结果，就有

$$
f(z)=\frac1{2\pi i}\oint_{|\zeta|=1}\frac{\zeta^*}{\zeta-z}d\zeta
=\begin{cases}
-\dfrac1z,& |z|>1,\\[4pt]
0,& |z|<1.
\end{cases}
$$

由此可见，$f(z)$ 在 $|z|\neq1$ 处解析，尽管 $\zeta^*$ 在全平面不解析.

利用 Cauchy 型积分，就可以推出含参量积分的解析性.

**定理 3.3（含参量积分的解析性）** 设

1. $f(t,z)$ 是 $t$ 和 $z$ 的连续函数，$t\in[a,b],z\in\overline G$，其中 $\overline G$ 是有界闭区域；
2. 对于 $[a,b]$ 上的任何 $t$ 值，$f(t,z)$ 是 $\overline G$ 上的单值解析函数，

则 $F(z)=\int_a^bf(t,z)dt$ 在 $G$ 内是解析的，且

$$
F'(z)=\int_a^b\frac{\partial f(t,z)}{\partial z}dt,\qquad z\in G. \tag{3.25}
$$

**证** 因为 $f(t,z)$ 在 $\overline G$ 中解析，故对于 $G$ 内的任何一点 $z$，Cauchy 积分公式成立，

$$
f(t,z)=\frac1{2\pi i}\oint_C\frac{f(t,\zeta)}{\zeta-z}d\zeta,
$$

其中 $C$ 是 $\overline G$ 的边界，积分沿 $C$ 的正方向. 代入 $F(z)$ 的定义，并交换积分次序（因为 $f(t,z)$ 连续），得

$$
F(z)=\int_a^b\frac{dt}{2\pi i}\oint_C\frac{f(t,\zeta)}{\zeta-z}d\zeta
=\frac1{2\pi i}\oint_C\frac1{\zeta-z}\left[\int_a^bf(t,\zeta)dt\right]d\zeta.
$$

这是一个 Cauchy 型积分，$\int_a^bf(t,z)dt$ 在 $C$ 上连续，故 $F(z)$ 为 $G$ 内的解析函数，且

$$
F'(z)=\frac1{2\pi i}\oint_C\frac1{(\zeta-z)^2}\left[\int_a^bf(t,\zeta)dt\right]d\zeta
$$

$$
=\int_a^b\left[\frac1{2\pi i}\oint_C\frac{f(t,\zeta)}{(\zeta-z)^2}d\zeta\right]dt
=\int_a^b\frac{\partial f(t,z)}{\partial z}dt.
$$

显然，这个结论也适用于 $\int_{C'}f(t,z)dt$. 这时应当要求 $C'$ 是分段光滑曲线. 当 $t$ 在 $C'$ 上变动，$z\in\overline G$ 时，$f(t,z)$ 是 $t$ 和 $z$ 的连续函数. 证明的方法与上面相同.

## *§3.7 Poisson 公式

Cauchy 积分公式告诉我们，对于有界闭区域 $\overline G$ 中的解析函数 $f(z)$，边界上的数值就完全决定了函数在 $G$ 内任意一点的值. 作为它的特殊情形，当 $f(z)$ 在上半平面解析时，就有

$$
f(z)=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{\zeta-z}d\zeta,
$$

其中 $C$ 是上半平面内的简单闭合围道，$z$ 为 $C$ 内任意一点.

现在取一个特殊的围道，它由实轴上的线段 $[-R,R]$ 和以原点为圆心、$R$ 为半径的半圆弧 $C_R$ 组成（见图 3.11）. 只要 $R$ 足够大，则 $z$ 点必在围道内. 令 $\zeta=\xi+i\eta$，就有

$$
f(z)=\frac1{2\pi i}\int_{-R}^{R}\frac{f(\xi)}{\xi-z}d\xi
+\frac1{2\pi i}\int_{C_R}\frac{f(\zeta)}{\zeta-z}d\zeta.
$$

![Figure 3.11](...)

图 3.11 上半平面的 Poisson 公式

下面讨论 $R\to\infty$ 的极限情形. 因为 $f(z)$ 在上半平面解析，故只需 $z$ 在上半平面趋于 $\infty$ 时 $f(z)$ 趋于 0，则根据 §3.4 中的引理 3.5 可知，

$$
\lim_{R\to\infty}\int_{C_R}\frac{f(\zeta)}{\zeta-z}d\zeta=0,
$$

由此即得

$$
f(z)=\frac1{2\pi i}\int_{-\infty}^{\infty}\frac{f(\xi)}{\xi-z}d\xi. \tag{3.26}
$$

此结果说明，由 $f(z)$ 在实轴上的值，可以唯一地决定它在上半平面任意一点的数值. 更进一步，令 $z=x+iy,f(z)=u+iv$，还可以单独地写出 $f(z)$ 的实部与虚部：

$$
u(x,y)=\frac1{2\pi}\int_{-\infty}^{\infty}
\frac{(\xi-x)v(\xi,0)+yu(\xi,0)}{(\xi-x)^2+y^2}d\xi, \tag{3.27}
$$

$$
v(x,y)=-\frac1{2\pi}\int_{-\infty}^{\infty}
\frac{(\xi-x)u(\xi,0)-yv(\xi,0)}{(\xi-x)^2+y^2}d\xi. \tag{3.28}
$$

由于解析函数的实部和虚部并不互相独立，由解析函数的实部就可以求出虚部（反之亦然），从而确定解析函数本身. 因此，可以设想，如果只知道解析函数的实部在实轴上的数值，也有可能求出虚部在实轴上的数值，从而完全决定函数在上半平面内任意一点 $z$ 的值. 为此，考虑 $z$ 相对于实轴的反演点 $z^*$. 它一定位于下半平面，所以

$$
\lim_{R\to\infty}\left\{\frac1{2\pi i}\oint_C\frac{f(\zeta)}{\zeta-z^*}d\zeta\right\}
=\frac1{2\pi i}\int_{-\infty}^{\infty}\frac{f(\xi)}{\xi-x+iy}d\xi=0. \tag{3.29}
$$

或者也分别比较等式两端的实部与虚部，

$$
\frac1{2\pi}\int_{-\infty}^{\infty}
\frac{(\xi-x)u(\xi,0)+yv(\xi,0)}{(\xi-x)^2+y^2}d\xi=0, \tag{3.30}
$$

$$
\frac1{2\pi}\int_{-\infty}^{\infty}
\frac{(\xi-x)v(\xi,0)-yu(\xi,0)}{(\xi-x)^2+y^2}d\xi=0. \tag{3.31}
$$

将 (3.27) 和 (3.31) 式相加，就得到

$$
u(x,y)=\frac1\pi\int_{-\infty}^{\infty}
\frac{(\xi-x)v(\xi,0)}{(\xi-x)^2+y^2}d\xi \tag{3.32a}
$$

$$
=\frac1\pi\int_{-\infty}^{\infty}
\frac{yu(\xi,0)}{(\xi-x)^2+y^2}d\xi. \tag{3.32b}
$$

将 (3.28) 和 (3.30) 式相加，又能得到

$$
v(x,y)=-\frac1\pi\int_{-\infty}^{\infty}
\frac{(\xi-x)u(\xi,0)}{(\xi-x)^2+y^2}d\xi \tag{3.33a}
$$

$$
=\frac1\pi\int_{-\infty}^{\infty}
\frac{yv(\xi,0)}{(\xi-x)^2+y^2}d\xi. \tag{3.33b}
$$

(3.32) 和 (3.33) 诸式的特点是，只根据 $u(x,y)$ 或 $v(x,y)$ 在实轴上的数值，便可以完全确定 $u(x,y)$ 和 $v(x,y)$ 在上半平面内任意一点的值. 更进一步，也就给出 $f(z)$ 本身在上半平面内任意一点的值：

$$
f(z)=\frac1{\pi i}\int_{-\infty}^{\infty}\frac{u(\xi,0)}{\xi-(x+iy)}d\xi \tag{3.34a}
$$

$$
=\frac1\pi\int_{-\infty}^{\infty}\frac{v(\xi,0)}{\xi-(x+iy)}d\xi. \tag{3.34b}
$$

而且将 (3.32) 和 (3.33) 的各式适当组合，或者直接将 (3.26) 和 (3.29) 两式相减，还能得到

$$
f(z)=\frac y\pi\int_{-\infty}^{\infty}\frac{f(\xi)}{(\xi-x)^2+y^2}d\xi \tag{3.35a}
$$

$$
=\frac1{\pi i}\int_{-\infty}^{\infty}\frac{(\xi-x)f(\xi)}{(\xi-x)^2+y^2}d\xi. \tag{3.35b}
$$

(3.32)–(3.26) 诸式均称为上半平面的 Poisson 公式. 这些结果说明：如果 $f(z)$ 在上半平面解析，并且当 $z$ 在上半平面趋于 $\infty$ 时趋于 0，根据它（或者它的实部或虚部）在实轴上的数值，就可以唯一地决定它在上半平面内任意一点的数值.

**思考题** §2.3 中曾经指出，由解析函数的实部（虚部）就可以求出虚部（实部），但可差一个任意实常数. 可是在上面的 (3.32a)、(3.33a) 和 (3.34a)、(3.34b) 诸式中并不出现任意常数，为什么？

下面导出另一种 Poisson 公式：圆内区域内的 Poisson 公式. 设函数 $f(z)$ 在圆 $|z|\leq a$ 中解析，则对于圆内任意一点 $z=re^{i\phi},r<a$（见图 3.13），根据有界区域 Cauchy 积分公式，有

$$
f(z)=\frac1{2\pi i}\oint_{|\zeta|=a}\frac{f(\zeta)}{\zeta-z}d\zeta
=\frac a{2\pi}\int_0^{2\pi}\frac{f(ae^{i\theta})}{a-re^{i(\phi-\theta)}}d\theta
$$

$$
=\frac a{2\pi}\int_0^{2\pi}\frac{a-re^{-i(\phi-\theta)}}{a^2+r^2-2ar\cos(\phi-\theta)}f(ae^{i\theta})d\theta. \tag{3.36}
$$

![Figure 3.12](...)

图 3.12 圆内的 Poisson 公式

另一方面，对于圆外的点 $z_1=a^2/z^*=r_1e^{i\phi},r_1=a^2/r$，也应该有

$$
\frac1{2\pi i}\oint_{|\zeta|=a}\frac{f(\zeta)}{\zeta-z_1}d\zeta
=\frac a{2\pi}\int_0^{2\pi}\frac{a-r_1e^{-i(\phi-\theta)}}{a^2+r_1^2-2ar_1\cos(\phi-\theta)}f(ae^{i\theta})d\theta
$$

$$
=-\frac r{2\pi}\int_0^{2\pi}\frac{r-ae^{-i(\phi-\theta)}}{a^2+r^2-2ar\cos(\phi-\theta)}f(ae^{i\theta})d\theta=0. \tag{3.37}
$$

将 (3.36) 和 (3.37) 两式相减，就得到圆内区域的 Poisson 公式

$$
f(z)=\frac{a^2-r^2}{2\pi}\int_0^{2\pi}\frac{f(ae^{i\theta})}{a^2+r^2-2ar\cos(\phi-\theta)}d\theta. \tag{3.38}
$$

或者令 $f(z)$ 的实部和虚部分别是 $u(r,\theta)$ 和 $v(r,\theta)$，则对它们也有

$$
u(r,\phi)=\frac{a^2-r^2}{2\pi}\int_0^{2\pi}\frac{u(a,\theta)}{a^2+r^2-2ar\cos(\phi-\theta)}d\theta, \tag{3.39}
$$

$$
v(r,\phi)=\frac{a^2-r^2}{2\pi}\int_0^{2\pi}\frac{v(a,\theta)}{a^2+r^2-2ar\cos(\phi-\theta)}d\theta. \tag{3.40}
$$

## 习题

1. 试按给定的路径计算下列积分：

(1)

$$
\int_0^{2+i}\operatorname{Re}z\,dz,
$$

积分路径为：

(i) 线段 $[0,2]$ 和 $[2,2+i]$ 组成的折线，

(ii) 线段 $z=(2+i)t,0\leq t\leq1$.

(2)

$$
\int_C\frac{dz}{\sqrt{z}},
$$

规定 $\left.\sqrt z\right|_{z=1}=1$，积分路径为由 $z=1$ 出发的：

(i) 单位圆的上半周，

(ii) 单位圆的下半周.

2. 计算下列积分：

(1)

$$
\oint_{|z|=1}\frac{dz}{z};
$$

(2)

$$
\oint_{|z|=1}\frac{|dz|}{z};
$$

(3)

$$
\oint_{|z|=1}\frac{dz}{|z|};
$$

(4)

$$
\oint_{|z|=1}\left|\frac{dz}{z}\right|.
$$

3. 计算下列积分：

(1)

$$
\oint_C\frac1{z^2-1}\sin\frac{\pi z}{4}\,dz,
$$

$C$ 分别为：

(i) $|z|=1/2$，

(ii) $|z-1|=1$，

(iii) $|z|=3$，

(iv) $|z|=R,R\to\infty$.

(2)

$$
\oint_C\frac1{z^2+1}e^{iz}dz,
$$

$C$ 分别为：

(i) $|z-i|=1$，

(ii) $|z|=2$，

(iii) $|z+i|+|z-i|=2\sqrt2$，

(iv) 闭合曲线 $r=3-\sin^2\frac\theta4$.

4. 计算下列积分：

(1)

$$
\oint_{|z|=2}\frac{\cos z}{z}dz;
$$

(2)

$$
\oint_{|z|=2}\frac{z^2-1}{z^2+1}dz;
$$

(3)

$$
\oint_{|z|=2}\frac{\sin(e^z)}{z}dz;
$$

(4)

$$
\oint_{|z|=2}\frac{e^z}{\cosh z}dz;
$$

(5)

$$
\oint_{|z|=2}\frac{\sin z}{z^2}dz;
$$

(6)

$$
\oint_{|z|=2}\frac{|z|e^z}{z^2}dz;
$$

(7)

$$
\oint_{|z|=2}\frac{\sin z}{z^4}dz;
$$

(8)

$$
\oint_{|z|=2}\frac{dz}{z^2(z^2+16)}.
$$

5. (1) 计算积分

$$
\oint_{|z|=1}\frac{e^z}{z^3}dz;
$$

(2) $a$ 取何值时，函数

$$
F(z)=\int_{z_0}^{z}e^z\left(\frac1z+\frac a{z^3}\right)dz
$$

是单值的？

6. 计算下列积分：

(1)

$$
f(z)=\oint_{|\zeta|=2}\frac{(\zeta^*)^2e^\zeta}{\zeta-z}d\zeta
$$

（其中 $\zeta^*$ 是复数 $\zeta$ 的复共轭，$z$ 是复数，且 $|z|\neq2$）；

(2)

$$
\oint_{|z|=1}\ln\frac{2z+1}{2z-1}\frac{dz}{z+2},
$$

在实轴上沿 $-\frac12$ 到 $\frac12$ 作割线，规定割线上岸

$$
\arg\frac{2z+1}{2z-1}=\pi.
$$

### 无界区域的高阶导数公式

可以模仿有界区域高阶导数公式的证明方法. 根据 (3.20) 式，我们知道，如果 $f(\zeta)$ 在简单闭合围道 $C$ 上及 $C$ 外解析，且当 $z\to\infty$ 时 $f(z)$ 的极限存在，则对于围外一点 $z$，

$$
f(z)-f(\infty)=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{\zeta-z}d\zeta,
$$

其中 $C$ 为顺时针方向，即绕 $\infty$ 点的正向. 因此，

$$
\frac{f(z+\Delta z)-f(z)}{\Delta z}
=\frac1{2\pi i}\frac1{\Delta z}\oint_C\left[\frac{f(\zeta)}{\zeta-z-\Delta z}-\frac{f(\zeta)}{\zeta-z}\right]d\zeta
$$

$$
=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z-\Delta z)(\zeta-z)}d\zeta.
$$

仿照有界区域的做法，就能证得

$$
f'(z)\equiv\lim_{\Delta z\to0}\frac{f(z+\Delta z)-f(z)}{\Delta z}
=\frac1{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^2}d\zeta.
$$

重复以上计算，还可以推得更普遍的结果：

$$
f^{(n)}(z)\equiv\lim_{\Delta z\to0}\frac{f^{(n-1)}(z+\Delta z)-f^{(n-1)}(z)}{\Delta z}
=\frac{n!}{2\pi i}\oint_C\frac{f(\zeta)}{(\zeta-z)^{n+1}}d\zeta,
\qquad n=1,2,3,\cdots.
$$
