# 第三章　解析函数的幂级数展开\*

## 目录

§1　复常数项级数　47  
§1.1　基本概念与基本判敛法　47  
§1.2　\*Abel 判别法与 Dirichlet 判别法　49  
§1.3　\*更序级数　51  
§1.4　\*级数的乘积　53  
§1.5　\*Gauss 判别法　53

§2　复函数项级数　55  
§2.1　普通函数项级数　55  
§2.2　一致收敛的函数项级数　55  
§2.3　解析函数项级数　60  
§2.4　\*Riemann $\zeta$ 函数、Riemann 假设与千禧数学难题　62

§3　幂级数　63  
§3.1　幂级数的敛散性　63  
§3.2　收敛半径的求法　64  
§3.3　和函数的解析性　65

§4　解析函数的 Taylor 展开　66  
§4.1　Taylor 定理　66  
§4.2　展开式的收敛半径　67  
§4.3　初等单值函数的 Taylor 展开式　68  
§4.4　初等多值函数的 Taylor 展开式　69

§5　解析开拓的基本概念　70

补充习题　71

本章先介绍复数级数、尤其是幂级数的性质，然后利用 Cauchy 积分公式导出解析函数的 Taylor 展开式。

## §1　复常数项级数

### §1.1　基本概念与基本判敛法

复常数项级数具有形式

$$
\sum_{k=1}^{\infty}\alpha_k=\alpha_1+\alpha_2+\cdots+\alpha_k+\cdots, \tag{1}
$$

其中每一项都是复数，$s_n=\sum_{k=1}^n\alpha_k$ 称为其部分和，如果 $\lim_{n\to\infty}s_n=s$ 是有限复数，则称级数（1）收敛于 $s$，$s$ 称为其和，记作 $\sum_{k=1}^{\infty}\alpha_k=s$，否则称级数（1）发散。

**注**　发散包括两种情况：① $\lim_{n\to\infty}s_n=\infty$；② 数列 $\{s_n\}_{n=1}^{\infty}$ 没有极限。

设 $\alpha_k=a_k+ib_k$（其中 $a_k\in\mathbb R$，$b_k\in\mathbb R$，$k\in\mathbb N^+$，$\mathbb N^+$ 表示正整数的集合），则不难证明，级数（1）收敛的充要条件是级数 $\sum_{k=1}^{\infty}a_k$ 和 $\sum_{k=1}^{\infty}b_k$ 均收敛。这使得我们原则上可以用实级数的判敛方法来对复级数进行判敛，但在实用上述并不是很有效。比如级数 $\sum_{k=1}^{\infty}(p+iq)^k$（其中 $p\in\mathbb R$，$q\in\mathbb R$），其通项的实部和虚部都比较复杂，难以用上述方法来判敛。

根据上述充要条件，级数 $\sum_{k=1}^{\infty}a_k$ 和 $\sum_{k=1}^{\infty}b_k$ 中只要有一个发散，则级数（1）发散；另一方面，如果级数（1）发散，则级数 $\sum_{k=1}^{\infty}a_k$ 和 $\sum_{k=1}^{\infty}b_k$ 不可能都收敛，但其中有一个收敛显然是可能的。比如当所有的 $a_k=0$，则 $\sum_{k=1}^{\infty}a_k$ 显然收敛，不过这是平庸的情况。非平庸一点的例子如级数 $\sum_{k=1}^{\infty}(1/k^2+i/k)$，它显然发散，但其实部所构成的级数是收敛的。

对于数列 $\{s_n\}_{n=1}^{\infty}$ 是否收敛，即其极限是否存在的问题，类似于实数列，我们有 Cauchy 收敛原理：数列 $\{s_n\}_{n=1}^{\infty}$ 收敛的充要条件是，$\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，$\forall p\in\mathbb N^+$，均有 $|s_{n+p}-s_n|<\varepsilon$。（这一精确描述的大意是：从某一项开始，数列的各项均聚集在其附近。这样就容易理解存在极限了。）

**注**　Cauchy 收敛原理也称为 Cauchy 判敛准则，它的另一种更常见的表述是：数列 $\{s_n\}_{n=1}^{\infty}$ 收敛的充要条件是，$\forall\varepsilon>0$，$\exists n_0\in\mathbb N^+$，当 $n>n_0$ 时，$\forall p\in\mathbb N^+$，均有 $|s_{n+p}-s_n|<\varepsilon$。（这一精确描述的大意是：当数列的下标足够大以后，任何两项的距离都可以任意小。更形象地说，越往后，数列的各项就越是挤在一起。这样也容易理解存在极限。）两种表述的区别是：在第一种表述中，$n$ 足够大但固定；在第二种表述中，$n$ 足够大但也是任意的。容易证明，两种表述是等价的。第一种表述显得更简洁，对于多数具体应用情景也更便于操作。举例来说，如果数列 $\{s_n\}_{n=1}^{\infty}$ 的通项对于 $s_{3n}$、$s_{3n+1}$ 和 $s_{3n+2}$ 具有不同的形式，用 Cauchy 收敛原理研究其敛散性时，如果按照第一种表述，只需要考察三种情况，而按照第二种表述，则需要考察九种情况。今后我们较多地采用第一种表述，但视具体情况有时候也采用第二种。

根据数列的 Cauchy 收敛原理，立得关于级数（1）的 Cauchy 收敛原理：

**定理（Cauchy 收敛原理）**　级数（1）收敛的充要条件是，$\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，$\forall p\in\mathbb N^+$，均有

$$
\left|\sum_{k=1}^{p}\alpha_{n+k}\right|<\varepsilon. \tag{2}
$$

由上式和三角不等式，易得

$$
|\alpha_{n+p}|=\left|\sum_{k=1}^{p}\alpha_{n+k}-\sum_{k=1}^{p-1}\alpha_{n+k}\right|<2\varepsilon,
$$

所以级数（1）收敛的必要条件是

$$
\lim_{k\to\infty}\alpha_k=0. \tag{3}
$$

如果这一条件得不到满足，我们马上可以断定级数是发散的。

**注**　如果采用 Cauchy 收敛原理的第二种表述，则式（2）中的 $n$ 也是任意的（在足够大的前提下），此时只需取 $p=1$，则 $|\alpha_{n+1}|<\varepsilon$ 即是收敛的必要条件。可见此时用第二种表述较方便。

Cauchy 收敛原理无疑是一个最基本的判别法，但是使用这一判别法需要计算式（2）中的求和，或至少要估计其上限（更准确地说，是估计其模的上限），这只有在少数简单情况下可以做得到，比如几何级数 $\sum_{k=1}^{\infty}\alpha^k$。

根据 Cauchy 收敛原理，在级数（1）中添加或去掉有限项，并不改变级数的收敛或发散性质。

如果级数

$$
\sum_{k=1}^{\infty}|\alpha_k| \tag{4}
$$

收敛，则称级数（1）绝对收敛。收敛但不绝对收敛的级数称为条件收敛。

**定理**　绝对收敛的级数必定收敛。

这是因为

$$
\left|\sum_{k=1}^{p}\alpha_{n+k}\right|\le\sum_{k=1}^{p}|\alpha_{n+k}|.
$$

设 $\alpha_k=a_k+ib_k$（其中 $a_k\in\mathbb R$，$b_k\in\mathbb R$，$k\in\mathbb N^+$），显然 $|a_k|\le|\alpha_k|$，$|b_k|\le|\alpha_k|$，而 $|\alpha_k|\le|a_k|+|b_k|$，由正项级数收敛的比较判别法，易知级数（1）绝对收敛的充要条件是级数 $\sum_{k=1}^{\infty}a_k$ 和 $\sum_{k=1}^{\infty}b_k$ 均绝对收敛。所以，级数 $\sum_{k=1}^{\infty}a_k$ 和 $\sum_{k=1}^{\infty}b_k$ 当中如果有一个不绝对收敛，则级数（1）必不绝对收敛；另一方面，级数（1）条件收敛甚至发散，并不排除级数 $\sum_{k=1}^{\infty}a_k$ 和 $\sum_{k=1}^{\infty}b_k$ 当中可以有一个绝对收敛。比如当所有的 $a_k=0$，则 $\sum_{k=1}^{\infty}a_k$ 显然绝对收敛，不过这是平庸的情况。非平庸一点的例子如级数 $\sum_{k=1}^{\infty}(\cos k\theta/k^2+i/k)$ 和 $\sum_{k=1}^{\infty}[\cos k\theta/k^2+i(-)^{k-1}/k]$（其中 $0\le\theta<2\pi$），前者显然发散，而后者条件收敛，但它们的实部所构成的级数是绝对收敛的。

在实际操作上，给定一个级数，我们先看看式（3）是否满足，若不满足，我们马上可以断定级数发散。若式（3）得到满足，再看相应的级数（4）是否收敛，若收敛，则原级数绝对收敛。由于级数（4）是正项级数，我们可以用正项级数的种种判别法来判敛，比如比较判别法、Cauchy 判别法、d’Alembert 判别法、Cauchy 积分判别法等。如果式（3）得到满足，而相应的级数（4）又不收敛，则需要进一步判断原级数是发散还是条件收敛。不过，如果是由 Cauchy 判别法或 d’Alembert 判别法得到级数（4）发散，则原级数必定也发散，因为这时式（3）必不满足。其它情况下要进一步判断原级数是否条件收敛就比较困难。这时可以尝试计算式（2）中的求和，或估计其上限，如果成功，就可以用 Cauchy 收敛原理来判断，如果式（2）中的和式无法计算或估计上限，则可以尝试下面介绍的方法。就我们所知，并没有什么方法能够百分之百地解决任何级数的判敛问题。即使是正项级数，虽然有一些更敏锐的判别法可资利用，比如后面介绍的 Gauss 判别法，但也并不存在一劳永逸的解决方案。

**例 1**　几何级数 $\sum_{k=1}^{\infty}\alpha^k$。当 $|\alpha|<1$ 时，由于正项级数 $\sum_{k=1}^{\infty}|\alpha|^k$ 收敛，故原级数绝对收敛。当 $|\alpha|\ge1$ 时，式（3）不能成立，故原级数发散。

### §1.2　\*Abel 判别法与 Dirichlet 判别法

本小节介绍两个判别法，用于对形如 $\sum_{k=1}^{\infty}a_k\alpha_k$（其中 $a_k\in\mathbb R$，$\alpha_k\in\mathbb C$）的级数进行判敛。

对于和式 $\sum_{k=1}^n a_k\alpha_k$，其中 $a_k\in\mathbb R$，$\alpha_k\in\mathbb C$，令

$$
A_k=\sum_{i=1}^k\alpha_i, \tag{5}
$$

则 $\alpha_k=A_k-A_{k-1}$，其中 $A_0\equiv0$。于是

$$
\sum_{k=1}^n a_k\alpha_k
=\sum_{k=1}^n a_k(A_k-A_{k-1})
=\sum_{k=1}^n a_kA_k-\sum_{k=1}^n a_kA_{k-1}
=\sum_{k=1}^n a_kA_k-\sum_{k=1}^{n-1}a_{k+1}A_k
=a_nA_n+\sum_{k=1}^{n-1}(a_k-a_{k+1})A_k,
$$

即

$$
\sum_{k=1}^n a_k\alpha_k=a_nA_n+\sum_{k=1}^{n-1}(a_k-a_{k+1})A_k. \tag{6}
$$

上式称为 Abel 变换，它把一个求和变成另一个求和，这对于和式的计算未必有帮助，但在估计和式的上限时会很有用，下面马上就可以看到。如果 $a_k\in\mathbb C$，上式仍然成立，但下面用到的都是 $a_k\in\mathbb R$ 的情况。

**引理（Abel）**　若 $\{a_k\}_{k=1}^n$ 单调（增加或减少），而 $\{A_k\}_{k=1}^n$ 有界，即 $\exists M>0$，使得 $|A_k|\le M$，则

$$
\left|\sum_{k=1}^n a_k\alpha_k\right|\le M(|a_1|+2|a_n|). \tag{7}
$$

**注**　“$\{a_k\}_{k=1}^n$ 单调”就意味着 $a_k\in\mathbb R$，因为虚部不为 0 的复数是没有大小可言的。

**证明**　由 Abel 变换式（6），有

$$
\left|\sum_{k=1}^n a_k\alpha_k\right|
\le|a_n||A_n|+\sum_{k=1}^{n-1}|a_k-a_{k+1}||A_k|
\le M|a_n|+\sum_{k=1}^{n-1}M|a_k-a_{k+1}|.
$$

若 $\{a_k\}_{k=1}^n$ 单调减少，则 $a_k-a_{k+1}$ 为正数，这时 $|a_k-a_{k+1}|=a_k-a_{k+1}$，故

$$
\left|\sum_{k=1}^n a_k\alpha_k\right|
\le M|a_n|+\sum_{k=1}^{n-1}M(a_k-a_{k+1})
=M|a_n|+M(a_1-a_n)
\le M(|a_1|+2|a_n|). \tag{8a}
$$

$\{a_k\}_{k=1}^n$ 单调减少的一种特殊情况是 $a_1\ge a_2\ge\cdots\ge a_n\ge0$，这时容易看出

$$
\left|\sum_{k=1}^n a_k\alpha_k\right|
\le M|a_n|+\sum_{k=1}^{n-1}M(a_k-a_{k+1})
=Ma_n+M(a_1-a_n)\le Ma_1. \tag{8b}
$$

这是一个更强的结果。若 $\{a_k\}_{k=1}^n$ 单调增加，则 $a_k-a_{k+1}$ 为负数，这时 $|a_k-a_{k+1}|=a_{k+1}-a_k$，故

$$
\left|\sum_{k=1}^n a_k\alpha_k\right|
\le M|a_n|+\sum_{k=1}^{n-1}M(a_{k+1}-a_k)
=M|a_n|+M(a_n-a_1)
\le M(|a_1|+2|a_n|). \tag{8c}
$$

综合式（8a）和式（8c），即得式（7），证毕。

有了 Abel 引理，马上可以得到以下两个判别法。

**定理（Abel 判别法）**　若数列 $\{a_k\}_{k=1}^{\infty}$ 单调有界，而级数 $\sum_{k=1}^{\infty}\alpha_k$ 收敛，则级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 收敛。

**注**　① 如果 $\sum_{k=1}^{\infty}\alpha_k$ 是正项级数且收敛，则只要数列 $\{a_k\}_{k=1}^{\infty}$ 有界，级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 就收敛，而且是绝对收敛。在直观上，这是很容易理解的。证明也很容易。事实上，数列 $\{a_k\}_{k=1}^{\infty}$ 有界意味着 $\exists K>0$，使得 $|a_k|\le K$，$\forall k$。今

$$
\left|\sum_{k=n+1}^{n+p}a_k\alpha_k\right|
=\sum_{k=n+1}^{n+p}|a_k|\alpha_k
\le K\sum_{k=n+1}^{n+p}\alpha_k
=K\left|\sum_{k=n+1}^{n+p}\alpha_k\right|,
$$

故由正项级数 $\sum_{k=1}^{\infty}\alpha_k$ 收敛，立得级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 绝对收敛。② 但是，现在 $\alpha_k$ 可正可负，也可能是复数，这时，数列 $\{a_k\}_{k=1}^{\infty}$ 不仅要有界，而且必须单调，才能保证级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 收敛。这一单调的要求虽然不是很直观，但很容易举例说明这一要求是无法取消的。比如 $a_k=(-)^{k-1}$ 显然有界，但不单调，若 $\alpha_k=(-)^{k-1}/k$，则 $\sum_{k=1}^{\infty}\alpha_k=\sum_{k=1}^{\infty}(-)^{k-1}/k$ 收敛，但 $\sum_{k=1}^{\infty}a_k\alpha_k=\sum_{k=1}^{\infty}1/k$ 发散。有的读者可能会认为这是由于 $a_k$ 可正可负造成的，那么如果 $a_k$ 非负且有界，是否就可以不必单调呢？答案还是不可以。比如 $a_k=1+(-)^{k-1}/\sqrt{k}$ 显然有界而且是恒正的，但不单调，如果 $\alpha_k=(-)^{k-1}/\sqrt{k}$，则 $\sum_{k=1}^{\infty}\alpha_k=\sum_{k=1}^{\infty}(-)^{k-1}/\sqrt{k}$ 收敛，但这时 $\sum_{k=1}^{\infty}a_k\alpha_k=\sum_{k=1}^{\infty}[1/k+(-)^{k-1}/\sqrt{k}]$ 是发散的。由此可见，数列 $\{a_k\}_{k=1}^{\infty}$ 单调的要求是很重要的。

**证明**　已知数列 $\{a_k\}_{k=1}^{\infty}$ 有界，则 $\exists K>0$，使得 $|a_k|\le K$，$\forall k$。而级数 $\sum_{k=1}^{\infty}\alpha_k$ 收敛，故 $\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使得

$$
\left|\sum_{k=n+1}^{n+m}\alpha_k\right|\le\varepsilon,\qquad \forall m\in\mathbb N^+.
$$

又数列 $\{a_k\}_{k=1}^{\infty}$ 单调，由 Abel 引理，

$$
\left|\sum_{k=n+1}^{n+p}a_k\alpha_k\right|
\le\varepsilon(|a_{n+1}|+2|a_{n+p}|)\le3K\varepsilon,\qquad \forall p\in\mathbb N^+.
$$

由 Cauchy 收敛原理，级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 收敛。证毕。

**例 2**　如果级数 $\sum_{k=1}^{\infty}\alpha_k$ 收敛，则由 Abel 判别法易知以下几个级数均收敛：

1. $\displaystyle\sum_{k=1}^{\infty}\frac{\alpha_k}{k^s}$，其中 $s>0$。
2. $\displaystyle\sum_{k=1}^{\infty}\frac{ak+b}{ck+d}\alpha_k$，其中 $a,b,c,d\in\mathbb R$ 且 $cd>0$（思考 $ad-bc$ 为正、负和零三种情况）。
3. $\displaystyle\sum_{k=1}^{\infty}(\tanh k)\alpha_k$。
4. $\displaystyle\sum_{k=1}^{\infty}(\coth k)\alpha_k$。

**定理（Dirichlet 判别法）**　若数列 $\{a_k\}_{k=1}^{\infty}$ 单调趋于 0，而级数 $\sum_{k=1}^{\infty}\alpha_k$ 的部分和有界，则级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 收敛。

**注**　Abel 判别法可以看作 Dirichlet 判别法的一个推论。事实上，按 Abel 判别法的条件，数列 $\{a_k\}_{k=1}^{\infty}$ 单调有界，则必有极限，设极限为 $a$，则数列 $\{a_k-a\}_{k=1}^{\infty}$ 单调趋于 0，而级数 $\sum_{k=1}^{\infty}\alpha_k$ 因为收敛而部分和必定有界，于是，按照 Dirichlet 判别法，级数 $\sum_{k=1}^{\infty}(a_k-a)\alpha_k$ 收敛。另一方面，$a_k\alpha_k=(a_k-a)\alpha_k+a\alpha_k$，故原级数收敛。可见理论上只要有 Dirichlet 判别法就足够了。但在一些实际问题中直接引用 Abel 判别法是方便的。

**证明**　已知数列 $\{a_k\}_{k=1}^{\infty}$ 趋于 0，则 $\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使得当 $k>n$ 时，有 $|a_k|<\varepsilon$。而级数 $\sum_{k=1}^{\infty}\alpha_k$ 的部分和有界，故 $\exists M>0$，使得 $|s_k|\le M$，$\forall k\in\mathbb N^+$，于是，

$$
\left|\sum_{k=n+1}^{n+m}\alpha_k\right|=|s_{n+m}-s_n|\le|s_{n+m}|+|s_n|\le2M,\qquad \forall n,m\in\mathbb N^+.
$$

又数列 $\{a_k\}_{k=1}^{\infty}$ 单调，由 Abel 引理，

$$
\left|\sum_{k=n+1}^{n+p}a_k\alpha_k\right|
\le2M(|a_{n+1}|+2|a_{n+p}|)<6M\varepsilon,\qquad \forall p\in\mathbb N^+.
$$

由 Cauchy 收敛原理，级数 $\sum_{k=1}^{\infty}a_k\alpha_k$ 收敛。证毕。

**例 3**　研究级数 $\sum_{k=1}^{\infty}e^{ik\theta}/k^s$ 的敛散性，其中 $0\le\theta<2\pi$，$s>0$。

**解**　由于 $s>0$，收敛的必要条件总是满足的。

（1）如果 $s>1$，由于其模构成的级数 $\sum_{k=1}^{\infty}1/k^s$ 收敛，故原级数绝对收敛。由此可知，级数 $\sum_{k=1}^{\infty}\cos k\theta/k^s$ 和 $\sum_{k=1}^{\infty}\sin k\theta/k^s$ 当 $s>1$ 时绝对收敛。

（2）如果 $0<s\le1$，由于其模构成的级数 $\sum_{k=1}^{\infty}1/k^s$ 发散，故原级数不能绝对收敛。若 $\theta=0$，原级数亦成为 $\sum_{k=1}^{\infty}1/k^s$，显然发散。若 $\theta\ne0$，我们就有

$$
\left|\sum_{k=1}^{n}e^{ik\theta}\right|
=\left|\frac{e^{i\theta}(1-e^{in\theta})}{1-e^{i\theta}}\right|
=\left|\frac{\sin(n\theta/2)}{\sin(\theta/2)}\right|
\le\frac1{\sin(\theta/2)}, \tag{9}
$$

亦即级数 $\sum_{k=1}^{\infty}e^{ik\theta}$ 的部分和数列有界，而数列 $\{1/k^s\}_{k=1}^{\infty}$ 单调趋于 0，根据上面的 Dirichlet 判别法，原级数收敛，亦即条件收敛。由此可知，级数 $\sum_{k=1}^{\infty}\cos k\theta/k^s$ 和 $\sum_{k=1}^{\infty}\sin k\theta/k^s$ 当 $0<s\le1$ 且 $\theta\ne0$ 时收敛，不过对于后一复级数，$\theta\ne0$ 的限制是不必要的，因为 $\theta=0$ 时它的通项为 0，是一个平庸的级数。这两个级数之所以能够收敛，主要是因为三角函数的出现使得各项有正有负，互相抵消。

前面指出，当原级数条件收敛时，级数 $\sum_{k=1}^{\infty}\cos k\theta/k^s$ 和 $\sum_{k=1}^{\infty}\sin k\theta/k^s$ 当中有一个绝对收敛是可能的。不过，实际上这两个级数都不能绝对收敛。事实上，$|\cos k\theta|\ge\cos^2k\theta=(1+\cos2k\theta)/2$，由于 $\sum_{k=1}^{\infty}1/k^s$ 发散而 $\sum_{k=1}^{\infty}\cos2k\theta/k^s$ 收敛，故正项级数 $\sum_{k=1}^{\infty}(1+\cos2k\theta)/2k^s$ 发散，由正项级数的比较判别法，$\sum_{k=1}^{\infty}|\cos k\theta|/k^s$ 发散。类似可证 $\sum_{k=1}^{\infty}|\sin k\theta|/k^s$ 发散。

### §1.3　\*更序级数

在实级数理论中，我们已经知道，对收敛级数的相邻各项任意加上括号形成的新级数仍然收敛，因为新级数的部分和数列是原级数部分和数列的子列。这一结论对于复级数显然也成立。这可以理解为：加法的结合律对于收敛的级数（无穷多项加法）仍然成立。

本小节研究的问题可以直观理解为：加法交换律对级数是否成立？结论是：只对绝对收敛的级数成立。

设 $\varphi(k)$ 是 $\mathbb N^+\to\mathbb N^+$ 的一个双射（即一一对应），并且是非平庸的（即不是恒等映射），那么级数 $\sum_{k=1}^{\infty}\alpha_{\varphi(k)}$ 称为 $\sum_{k=1}^{\infty}\alpha_k$ 的更序级数（rearranged series）。

**例 4**　级数 $\alpha_2+\alpha_1+\alpha_4+\alpha_3+\cdots$ 就是 $\sum_{k=1}^{\infty}\alpha_k$ 的一个更序级数，此时 $\varphi(2k-1)=2k$，$\varphi(2k)=2k-1$，其中 $k\in\mathbb N^+$。

如果 $\sum_{k=1}^{\infty}\alpha_k$ 条件收敛，则 $\sum_{k=1}^{\infty}\alpha_{\varphi(k)}$ 可能收敛或发散，如果收敛，其和与前者也不一定相等。事实上，对于实级数，Riemann 证明了 $\forall l\in\mathbb R$ 或 $l=\pm\infty$，均可找到 $\varphi(k)$ 使得 $\sum_{k=1}^{\infty}\alpha_{\varphi(k)}$ 收敛于 $l$。

**例 5**　考虑级数

$$
\sum_{k=1}^{\infty}\frac{(-)^{k-1}}{k}=1-\frac12+\frac13-\frac14+\cdots,
$$

它是收敛的，其和记作 $s$（我们知道 $s=\ln2$），部分和记作 $s_n$。今作更序级数

$$
1+\frac13-\frac12+\frac15+\frac17-\frac14+\cdots+\frac1{4n-3}+\frac1{4n-1}-\frac1{2n}+\cdots.
$$

可以证明，这个级数仍然是收敛的（证明留给读者思考），其和记作 $\Sigma$，部分和记作 $\Sigma_n$。令 $\sigma_n=1+1/2+\cdots+1/n$，则不难证明

$$
s_{2n}=\sigma_{2n}-\sigma_n,\qquad
\Sigma_{3n}=s_{4n}+\frac12s_{2n}-\frac12\sigma_n.
$$

结合两式可得

$$
\Sigma_{3n}=s_{4n}+\frac12s_{2n}.
$$

令 $n\to\infty$，即得 $\Sigma=3s/2$，可见这个更序级数具有不同的和。

但是，对于绝对收敛的级数，我们有下述定理，它表明绝对收敛的级数可以像有限的求和一样随意交换各项的次序。

**定理**　如果 $\sum_{k=1}^{\infty}\alpha_k$ 绝对收敛于 $s$，则其更序级数 $\sum_{k=1}^{\infty}\alpha_{\varphi(k)}$ 亦绝对收敛于 $s$。

**注**　将 $\sum_{k=1}^{\infty}\alpha_k$ 的实部和虚部分开，利用实级数的类似定理，即可证明上述定理。下面的证明是独立的，不需要引用实级数的相关定理。

**证明**　首先证明 $\sum_{k=1}^{\infty}\alpha_{\varphi(k)}$ 绝对收敛。记 $\sum_{k=1}^{\infty}|\alpha_k|$ 的和为 $\sigma$，部分和为 $\sigma_n$。考虑其更序级数 $\sum_{k=1}^{\infty}|\alpha_{\varphi(k)}|$ 的部分和

$$
\sigma_n'=\sum_{k=1}^n|\alpha_{\varphi(k)}|.
$$

令 $N=\max\{\varphi(k)\}_{k=1}^n$，显然 $N\ge n$，则易知

$$
\sigma_n'\le\sum_{k=1}^{N}|\alpha_{\varphi(k)}|=\sigma_N\le\sigma,\qquad \forall n\in\mathbb N^+.
$$

可见部分和数列 $\{\sigma_n'\}_{n=1}^{\infty}$ 有界，故必有极限（因为正项级数的部分和数列单调增长），设极限为 $\sigma'$，对上式取极限可得 $\sigma'\le\sigma$。另一方面，$\sum_{k=1}^{\infty}|\alpha_k|$ 也可以看作 $\sum_{k=1}^{\infty}|\alpha_{\varphi(k)}|$ 的更序级数，故又有 $\sigma\le\sigma'$。综合两者即得 $\sigma'=\sigma$。这表明 $\sum_{k=1}^{\infty}\alpha_k$ 更序后仍然绝对收敛，而且 $\sum_{k=1}^{\infty}|\alpha_k|$ 不因更序而改变其和。后者实际上是关于正项级数更序的一个结论。

由于 $\sum_{k=1}^{\infty}\alpha_k$ 绝对收敛于 $s$，故 $\forall\varepsilon>0$，必能找到 $n\in\mathbb N^+$，使得

$$
|s_n-s|<\frac\varepsilon2,\qquad
\sum_{k=1}^{p}|\alpha_{n+k}|<\frac\varepsilon2,\quad \forall p\in\mathbb N^+.
$$

对于这样的 $n$，总可以找到 $m\ge n$，使得 $s_m'=\sum_{k=1}^m\alpha_{\varphi(k)}$ 包含 $s_n$。但 $s_m'$ 中可能还包含其他 $\alpha_k$ 项，其中 $k>n$。对于 $N>m$，$s_N'$ 中会有更多这样的项，可以表示为

$$
s_N'=s_n+\sum_{k>n}'\alpha_k,
$$

这里求和号上的撇号表示求和仅包含 $s_N'$ 中出现的 $k>n$ 的 $\alpha_k$ 项，而非所有 $k>n$ 的 $\alpha_k$ 项。利用上面的不等式和三角不等式，可得

$$
|s_N'-s|\le|s_n-s|+\sum_{k>n}'|\alpha_k|<\frac\varepsilon2+\frac\varepsilon2=\varepsilon,\qquad \forall N>m.
$$

因此

$$
\lim_{N\to\infty}s_N'=s,\qquad \text{或}\qquad \sum_{k=1}^{\infty}\alpha_{\varphi(k)}=s.
$$

证毕。

### §1.4　\*级数的乘积

本小节研究的问题可以大致理解为：乘法分配律对两个级数相乘是否成立？结论是：如果两者均绝对收敛则成立。

考虑两个收敛的级数 $\sum_{k=1}^{\infty}\alpha_k$ 和 $\sum_{k=1}^{\infty}\beta_k$，由所有形为 $\alpha_k\beta_l$ 的项按一定顺序构成级数，可以看作两个级数的某种乘积。我们希望后者收敛，并且其和与各项的顺序无关，否则谈论两个级数的乘积意义不大。这能否实现，显然与两个级数的性质（比如是绝对收敛还是条件收敛）有关。事实上，我们有下述定理。它给出了一个充分条件。

**定理**　设级数 $\sum_{k=1}^{\infty}\alpha_k$ 绝对收敛于 $s$，$\sum_{k=1}^{\infty}\beta_k$ 绝对收敛于 $t$，则由所有 $\alpha_k\beta_l$ 项按任意顺序构成的级数绝对收敛于 $st$。

**证明**　已知 $\sum_{k=1}^{\infty}\alpha_k$ 和 $\sum_{k=1}^{\infty}\beta_k$ 绝对收敛。记 $\sum_{k=1}^{\infty}|\alpha_k|=\sigma_1$，$\sum_{k=1}^{\infty}|\beta_k|=\sigma_2$。考虑乘积级数，设其通项为 $w_k=\alpha_{i_k}\beta_{j_k}$，考虑 $\sum_{k=1}^{\infty}|w_k|$，其部分和为

$$
\sigma_N=\sum_{k=1}^{N}|w_k|=\sum_{k=1}^{N}|\alpha_{i_k}||\beta_{j_k}|.
$$

令 $m=\max\{i_k\}_{k=1}^N$，$n=\max\{j_k\}_{k=1}^N$，则易知

$$
\sigma_N\le\sum_{k=1}^{m}|\alpha_k|\sum_{l=1}^{n}|\beta_l|\le\sigma_1\sigma_2,\qquad \forall N\in\mathbb N^+.
$$

故 $\sum_{k=1}^{\infty}w_k$ 绝对收敛。

根据上一小节的结论，绝对收敛的级数可以任意排序而不改变其和。另一方面，对于收敛级数，可以将相邻各项加括号而不改变其收敛性及其和。于是我们可以将乘积级数的各项适当排序并加括号形成级数

$$
\gamma_1=\alpha_1\beta_1,\qquad
\gamma_2=\alpha_2\beta_1+\alpha_2\beta_2+\alpha_1\beta_2,\qquad\cdots,\qquad
\gamma_n=s_nt_n-s_{n-1}t_{n-1},\qquad\cdots
$$

其中 $s_n=\sum_{k=1}^n\alpha_k$，$t_n=\sum_{k=1}^n\beta_k$。易知 $\sum_{k=1}^{\infty}\gamma_k$ 包含了所有 $\alpha_k\beta_l$ 项，又易得

$$
\sum_{k=1}^{n}\gamma_k=s_nt_n.
$$

取极限即得

$$
\sum_{k=1}^{\infty}\gamma_k=st,\qquad \text{或}\qquad
\sum_{k=1}^{\infty}\alpha_k\sum_{l=1}^{\infty}\beta_l=st.
$$

证毕。

### §1.5　\*Gauss 判别法

当 d’Alembert 判别法中的极限为 1 时，为了判断级数是否收敛，需要有更敏锐更细致（当然也就更复杂一些）的判别法。在这样的判别法中，Gauss 判别法是非常有用的一个。在数理方程部分，证明 Legendre 方程的级数解在求解区间端点发散，就需要引用这一判别法。Gauss 判别法实际上是关于正项级数是否收敛的一个判别法。对于复级数，可以用它来判断是否绝对收敛。但是需要注意，当判断结果是发散时，原来的级数可能是条件收敛的。

**定理（Gauss 判别法）**　对于正项级数 $\sum_{k=1}^{\infty}a_k$，设相邻两项之比可以表示为（至少对于足够大的 $k$）

$$
\frac{a_k}{a_{k+1}}=1+\frac{r}{k}+\frac{\theta_k}{k^2}, \tag{10}
$$

其中 $r$ 是常数，$|\theta_k|\le L$ 有界。则当 $r>1$ 时级数收敛，而 $r\le1$ 时级数发散。

**注**　按式（10），如果用 d’Alembert 判别法，因比值极限为 1，无法得出结论。按定理，级数收敛要求 $r>1$，这意味着 $a_k$ 随着 $k$ 下降要足够快。

**证明**　由所给条件易得

$$
\lim_{k\to\infty}k\left(\frac{a_k}{a_{k+1}}-1\right)=r.
$$

另一方面，考虑级数 $\sum_{k=1}^{\infty}c_k=\sum_{k=1}^{\infty}1/k^s$，其中 $s\in\mathbb R$，易得

$$
\lim_{k\to\infty}k\left(\frac{c_k}{c_{k+1}}-1\right)=s.
$$

先考虑 $r>1$。总可取 $s$ 使得 $1<s<r$，此时 $\sum_{k=1}^{\infty}c_k$ 收敛。由于 $r>s$，必存在 $N\in\mathbb N^+$，当 $k\ge N$，有

$$
k\left(\frac{a_k}{a_{k+1}}-1\right)>k\left(\frac{c_k}{c_{k+1}}-1\right),\qquad k\ge N,
$$

或

$$
\frac{a_{k+1}}{a_k}<\frac{c_{k+1}}{c_k},\qquad k\ge N.
$$

对 $k=N,N+1,\cdots,n-1$ 写下上式并将各式相乘，则可得

$$
a_n<\frac{a_N}{c_N}c_n,\qquad n>N.
$$

由此即知级数 $\sum_{k=1}^{\infty}a_k$ 收敛。

再考虑 $r<1$（允许 $r\le0$）。总可取 $s$ 使得 $r<s<1$，此时 $\sum_{k=1}^{\infty}c_k$ 发散。由于 $r<s$，类似可证存在 $N\in\mathbb N^+$，当 $k\ge N$，有

$$
\frac{a_{k+1}}{a_k}>\frac{c_{k+1}}{c_k},\qquad k\ge N.
$$

从而

$$
a_n>\frac{a_N}{c_N}c_n,\qquad n>N.
$$

于是可知级数 $\sum_{k=1}^{\infty}a_k$ 发散。

最后考虑 $r=1$。此时取参考级数为 $\sum_{k=1}^{\infty}c_k=\sum_{k=1}^{\infty}1/(k\ln k)$，它是发散的。根据已知条件可得

$$
\frac1{c_k}\frac{a_k}{a_{k+1}}-\frac1{c_{k+1}}
=-(k+1)\ln\left(1+\frac1k\right)+\frac{\theta_k\ln k}{k}\to-1,\qquad k\to\infty.
$$

由此可知存在 $N\in\mathbb N^+$，当 $k\ge N$，上式左边为负数，从而

$$
\frac{a_{k+1}}{a_k}>\frac{c_{k+1}}{c_k},\qquad k\ge N.
$$

类似上面第二种情况，可知级数 $\sum_{k=1}^{\infty}a_k$ 发散。证毕。

**例 6**　考虑以下两个级数的敛散性。（1）$\sum_{k=0}^{\infty}a_k$，其中

$$
\frac{a_{k+1}}{a_k}=\frac{2k(2k+1)-\lambda}{(2k+1)(2k+2)},\qquad \lambda\in\mathbb R.
$$

（2）$\sum_{k=0}^{\infty}b_k$，其中

$$
\frac{b_{k+1}}{b_k}=\frac{(2k+1)(2k+2)-\lambda}{(2k+2)(2k+3)},\qquad \lambda\in\mathbb R.
$$

先考虑（1）的级数。首先获知从 $k=0$ 开始，比前面的一般形式多了一项，但这不影响级数的敛散性，故可以忽略。其次，如果 $\lambda=2n(2n+1)$，$n\in\mathbb N$，则 $a_{n+1}=a_{n+2}=\cdots=0$，级数只有前 $n+1$ 项，即到 $a_n$ 为止，此时不存在收敛问题。最后，设 $\lambda$ 不取上述特殊值，则级数有无穷多项，当 $k$ 足够大时，各项同号，可按正项级数处理，不难证明，此时 $a_k/a_{k+1}$ 具有式（10）的形式，其中

$$
r=1,\qquad \theta_k=\frac{\lambda k(k+1)}{2k(2k+1)-\lambda}.
$$

显然 $\theta_k$ 有界，故根据 Gauss 判别法，级数发散。这个级数是 Legendre 方程的级数解 $y_0(x)$ 在 $x=1$ 处的形式。

再考虑（2）的级数。如果 $\lambda=(2n+1)(2n+2)$，$n\in\mathbb N$，则级数至于 $b_n$，不存在收敛问题。设 $\lambda$ 不取上述特殊值，则级数有无穷多项，此时 $b_k/b_{k+1}$ 具有式（10）的形式，其中

$$
r=1,\qquad \theta_k=\frac{(\lambda-2)k(k+1)}{(2k+1)(2k+2)-\lambda}.
$$

显然 $\theta_k$ 有界，级数发散。这个级数是 Legendre 方程的级数解 $y_1(x)$ 在 $x=1$ 处的形式。

## §2　复函数项级数

### §2.1　普通函数项级数

设级数

$$
\sum_{k=1}^{\infty}f_k(z) \tag{11}
$$

的每一项都是定义在点集 $E$ 上的复变函数，则该级数称为复函数项级数，$s_n(z)=\sum_{k=1}^n f_k(z)$ 称为其部分和。如果对 $E$ 中的一点 $z_0$，常数项级数 $\sum_{k=1}^{\infty}f_k(z_0)$ 收敛，则称级数（11）在 $z_0$ 收敛；如果 $\forall z\in E$，级数（11）都收敛，则称它在点集 $E$ 上收敛，这时级数在 $E$ 上每点有和，记作 $f(z)$，称为其和函数，即

$$
\sum_{k=1}^{\infty}f_k(z)=f(z).
$$

从上述定义可以看出，这是一个逐点收敛的概念。取定一点，则级数（11）成为一个常数项级数，其判敛方法与上节所述无异。

### §2.2　一致收敛的函数项级数

对于函数项级数的研究，不能只是满足于逐点收敛的概念。因为即使级数（11）在 $E$ 上处处收敛，我们也看不出它在 $E$ 上各点的行为有何共性。既然是函数项级数，我们自然更关心这种共性。更重要的是，即使级数（11）在 $E$ 上处处收敛，我们也无法由级数中各项均具有某种性质（比如连续、可导）而推断其和函数具有同种性质。为了研究和函数的性质，一个关键的概念是一致收敛，即级数在各点的收敛速度大致相同。这一概念的精确描述如下。

**定义（一致收敛）**　设在点集 $E$ 上存在函数 $f(z)$，如果 $\forall\varepsilon>0$，$\exists n_0\in\mathbb N^+$，当 $n>n_0$ 时，

$$
|s_n(z)-f(z)|<\varepsilon
$$

在 $E$ 上处处成立，则称级数（11）在 $E$ 上一致收敛于 $f(z)$。

按这一定义，一致收敛当然就收敛，但反之则不然。让我们看看收敛与一致收敛的区别在于何处。如果级数（11）在 $E$ 上收敛于 $f(z)$，则 $\forall\varepsilon>0$，对于每一点 $z\in E$，$\exists n_0\in\mathbb N^+$，当 $n>n_0$ 时，$|s_n(z)-f(z)|<\varepsilon$ 在该点 $z$ 成立，由此可以看出，对于同样的 $\varepsilon$，不同的点 $z$ 可能需要不同的 $n_0$，即 $n_0=n_0(\varepsilon,z)$。这样就存在一种可能性，当 $z$ 沿某一趋势变化时，比如 $z\to z_0$（其中 $z,z_0\in E$）时，所需的 $n_0(\varepsilon,z)$ 越来越大，以至于我们无法找到其有限的上界，亦即无法找到一个共同的 $n_0$，使得当 $n>n_0$ 时，$|s_n(z)-f(z)|<\varepsilon$ 在 $E$ 上一致成立。而一致收敛的定义正是要求存在这样一个共同的 $n_0=n_0(\varepsilon)$，它不依赖于 $z$，而只与 $\varepsilon$ 有关。

一致收敛的基本判别法是 Cauchy 一致收敛原理：

**定理（Cauchy 一致收敛原理）**　级数（11）在 $E$ 上一致收敛的充要条件是，$\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使得 $\forall p\in\mathbb N^+$ 和 $z\in E$，均有

$$
\left|\sum_{k=1}^{p}f_{n+k}(z)\right|<\varepsilon. \tag{12}
$$

**证明**　必要性。按定义，若级数（11）一致收敛，则存在函数 $f(z)$，$\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使得 $\forall p\in\mathbb N^+$ 和 $z\in E$，下面两式同时成立：

$$
|s_n(z)-f(z)|<\varepsilon/2,\qquad |s_{n+p}(z)-f(z)|<\varepsilon/2.
$$

于是，$\forall p\in\mathbb N^+$ 和 $z\in E$，有

$$
\left|\sum_{k=1}^{p}f_{n+k}(z)\right|
=|s_{n+p}(z)-s_n(z)|
=|[s_{n+p}(z)-f(z)]-[s_n(z)-f(z)]|
$$

$$
\le|s_{n+p}(z)-f(z)|+|s_n(z)-f(z)|<\varepsilon.
$$

充分性。按已知条件，级数（11）在 $E$ 上收敛，故在每点有和，记和函数为 $f(z)$，即 $\lim_{n\to\infty}s_n(z)=f(z)$。又按已知条件，$\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使得 $\forall p\in\mathbb N^+$ 和 $z\in E$，式（12）成立，将它改写为

$$
|s_n(z)-s_{n+p}(z)|<\varepsilon,
$$

令 $p\to\infty$，得到

$$
|s_n(z)-f(z)|\le\varepsilon,
$$

由于该式在 $E$ 上一致成立，故原级数一致收敛。证毕。

**例 1**　级数

$$
\sum_{k=1}^{\infty}\frac{(-)^{k-1}}{k+x^2}
$$

在 $x\in(-\infty,+\infty)$ 上一致收敛。事实上，当 $p=2q$，$q\in\mathbb N^+$，有

$$
\left|\sum_{k=n+1}^{n+2q}\frac{(-)^{k-1}}{k+x^2}\right|
=\frac1{(n+1)+x^2}-\left[\frac1{(n+2)+x^2}-\frac1{(n+3)+x^2}\right]-\cdots
$$

$$
-\left[\frac1{(n+2q-2)+x^2}-\frac1{(n+2q-1)+x^2}\right]-\frac1{(n+2q)+x^2}
<\frac1{(n+1)+x^2}<\frac1n;
$$

当 $p=2q-1$，$q\in\mathbb N^+$，有（第三行 $\le$ 号中的等号当 $q=1$，亦即 $p=1$ 时成立）

$$
\left|\sum_{k=n+1}^{n+2q-1}\frac{(-)^{k-1}}{k+x^2}\right|
=\frac1{(n+1)+x^2}-\left[\frac1{(n+2)+x^2}-\frac1{(n+3)+x^2}\right]-\cdots
$$

$$
-\left[\frac1{(n+2q-2)+x^2}-\frac1{(n+2q-1)+x^2}\right]
\le\frac1{(n+1)+x^2}<\frac1n.
$$

因此，只要取 $n_0>1/\varepsilon$，则当 $n>n_0$ 时，

$$
\left|\sum_{k=n+1}^{n+p}\frac{(-)^{k-1}}{k+x^2}\right|<\varepsilon,
$$

对一切 $x\in(-\infty,+\infty)$ 和 $p\in\mathbb N^+$ 成立。

但对于任何 $x\in(-\infty,+\infty)$，该级数不能绝对收敛。事实上，

$$
\sum_{k=1}^{\infty}\left|\frac{(-)^{k-1}}{k+x^2}\right|=\sum_{k=1}^{\infty}\frac1{k+x^2},
$$

而右边的正项级数发散，这是因为正项级数 $\sum_{k=1}^{\infty}1/k$ 发散，而

$$
\lim_{k\to\infty}\frac{1/(k+x^2)}{1/k}
=\lim_{k\to\infty}\frac{k}{k+x^2}
=\lim_{k\to\infty}\frac1{1+x^2/k}=1.
$$

**例 2**　级数 $\sum_{k=1}^{\infty}z^k$ 当 $|z|<1$ 时绝对收敛，因为正项级数 $\sum_{k=1}^{\infty}|z|^k$ 当 $|z|<1$ 时收敛。但该级数在单位圆 $|z|<1$ 上不能一致收敛，事实上，易得和函数 $f(z)$ 与部分和 $s_n(z)$ 为

$$
f(z)=\frac{z}{1-z},\qquad s_n(z)=\frac{z(1-z^n)}{1-z},
$$

故

$$
|s_n(z)-f(z)|=\left|\frac{z^{n+1}}{1-z}\right|=\frac{|z|^{n+1}}{|1-z|},
$$

只要取

$$
n_0=\left[\frac{\ln(\varepsilon|1-z|)}{\ln|z|}\right],
$$

则当 $n>n_0$ 时，有

$$
|s_n(z)-f(z)|<\varepsilon.
$$

然而，上述 $n_0$ 不仅依赖于 $\varepsilon$，而且依赖于 $z$，当 $|z|\to1$ 时，$n_0$ 会越来越大，以至于无法找到一个共同的 $n_0$，使上式在 $|z|<1$ 上一致成立，因此原级数不能一致收敛。

由以上二例可以看到，一致收敛与绝对收敛之间没有因果关系。

下面给出的一致收敛的 Weierstrass M-判别法是一个充分条件，它是一个非常有用的判别法，但象上面的例 1 是无法用它来判断的。

下面给出一个一致收敛的充分条件，它是一个常用的判别法。

**定理（Weierstrass M-判别法）**　如果存在正数列 $\{M_k\}_{k=1}^{\infty}$，使得

$$
|f_k(z)|\le M_k,\qquad \forall z\in E,\quad k\in\mathbb N^+,
$$

且正项级数 $\sum_{k=1}^{\infty}M_k$ 收敛，则级数（11）在 $E$ 上绝对收敛且一致收敛。

**证明**　由正项级数 $\sum_{k=1}^{\infty}M_k$ 收敛和比较判别法，级数 $\sum_{k=1}^{\infty}|f_k(z)|$ 在 $E$ 上收敛，亦即级数（11）在 $E$ 上绝对收敛。

又正项级数 $\sum_{k=1}^{\infty}M_k$ 收敛意味着 $\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使得 $\forall p\in\mathbb N^+$，有

$$
\sum_{k=1}^{p}M_{n+k}<\varepsilon.
$$

由已知条件，

$$
\left|\sum_{k=1}^{p}f_{n+k}(z)\right|
\le\sum_{k=1}^{p}|f_{n+k}(z)|
\le\sum_{k=1}^{p}M_{n+k}<\varepsilon
$$

在 $E$ 上一致成立，由 Cauchy 一致收敛原理，级数（11）在 $E$ 上一致收敛。证毕。

**例 3**　级数 $\sum_{k=1}^{\infty}z^k$ 在闭圆 $|z|\le r<1$ 上绝对收敛且一致收敛。事实上，存在正项级数 $\sum_{k=1}^{\infty}r^k$ 满足上述定理的要求。但该级数在单位圆 $|z|<1$ 上不能一致收敛，参看上面的例 2。

**例 4**　级数 $\sum_{k=1}^{\infty}e^{ik\theta}/k^s$ 当 $s>1$ 时在 $0\le\theta<2\pi$ 上绝对收敛且一致收敛。事实上，存在正项级数 $\sum_{k=1}^{\infty}1/k^s$ 满足上述定理的要求。等价地说，级数 $\sum_{k=1}^{\infty}z^k/k^s$ 当 $s>1$ 时在单位圆周 $|z|=1$ 上绝对且一致收敛。事实上此时它在闭圆 $|z|\le1$ 上绝对且一致收敛。

**例 5**　由上节例 2 已经知道，级数 $\sum_{k=1}^{\infty}e^{ik\theta}/k^s$ 当 $0<s\le1$ 且 $\theta\ne0$ 时条件收敛。下面进一步证明，它在 $\delta\le\theta\le2\pi-\delta$ 上一致收敛，其中 $0<\delta<\pi$。

事实上，当 $\theta\ne0$ 时，由式（9），

$$
\left|\sum_{k=n+1}^{n+m}e^{ik\theta}\right|
=\left|\sum_{k=1}^{m}e^{ik\theta}\right|\le\frac1{\sin(\theta/2)},
$$

当 $\delta\le\theta\le2\pi-\delta$ 时，$\delta/2\le\theta/2\le\pi-\delta/2$，从而 $\sin(\delta/2)\le\sin(\theta/2)\le1$，故

$$
\left|\sum_{k=n+1}^{n+m}e^{ik\theta}\right|\le\frac1{\sin(\delta/2)}.
$$

由于数列 $\{1/k^s\}_{k=1}^{\infty}$ 单调减少且恒正，故由式（8b）可得

$$
\left|\sum_{k=n+1}^{n+p}\frac{e^{ik\theta}}{k^s}\right|
\le\frac1{\sin(\delta/2)}\frac1{(n+1)^s}
<\frac1{n^s\sin(\delta/2)}.
$$

只要取 $n>1/[\varepsilon\sin(\delta/2)]^{1/s}$，则 $\forall p\in\mathbb N^+$，

$$
\left|\sum_{k=n+1}^{n+p}\frac{e^{ik\theta}}{k^s}\right|<\varepsilon
$$

在 $\delta\le\theta\le2\pi-\delta$ 上一致成立，由 Cauchy 一致收敛原理，原级数在 $\delta\le\theta\le2\pi-\delta$ 上一致收敛。注意本题无法用 Weierstrass M-判别法来判断。

一致收敛的重要性体现在下面的定理上。

**定理（和函数连续）**　若级数（11）的各项均在点集 $E$ 上连续，且级数在 $E$ 上一致收敛于函数 $f(z)$，则 $f(z)$ 亦在 $E$ 上连续。

**证明**　已知级数在 $E$ 上一致收敛于函数 $f(z)$，故 $\forall\varepsilon>0$，$\exists n\in\mathbb N^+$，使

$$
|s_n(z)-f(z)|<\varepsilon/3
$$

在 $E$ 上一致成立。另一方面，由于级数的各项均在点集 $E$ 上连续，故 $s_n(z)$ 亦在点集 $E$ 上连续，于是，$\forall\varepsilon>0$，$\exists\delta>0$，使得当 $|z-z_0|<\delta$ 时，有

$$
|s_n(z)-s_n(z_0)|<\varepsilon/3.
$$

于是，当 $|z-z_0|<\delta$ 时，

$$
|f(z)-f(z_0)|=|f(z)-s_n(z)+s_n(z)-s_n(z_0)+s_n(z_0)-f(z_0)|
$$

$$
\le|f(z)-s_n(z)|+|s_n(z)-s_n(z_0)|+|s_n(z_0)-f(z_0)|
<\varepsilon/3+\varepsilon/3+\varepsilon/3=\varepsilon,
$$

亦即和函数在 $E$ 上连续。证毕。

**例 6**　考虑级数（11），其中 $f_1(z)=z$，$f_k(z)=z^k-z^{k-1}$（$k=2,3,\cdots$），显然 $s_n(z)=z^n$，故级数在单位圆 $|z|<1$ 和点 $z=1$ 处收敛，特别地，它在 $z=x\in[0,1]$ 上收敛，在该区间上，其和函数为

$$
f(x)=
\begin{cases}
0,&0\le x<1,\\
1,&x=1.
\end{cases}
$$

显然和函数在 $x=1$ 处不连续，这是由于级数在 $z=x\in[0,1]$ 上不能一致收敛。事实上

$$
|s_n(x)-f(x)|=
\begin{cases}
x^n,&0\le x<1,\\
0,&x=1.
\end{cases}
$$

为使 $n>n_0$ 时，$|s_n(x)-f(x)|<\varepsilon$，需取 $n_0\ge\ln\varepsilon/\ln x$，当 $x\to1$，$n_0$ 越来越大，因此无法找到共同的 $n_0$，使 $|s_n(x)-f(x)|<\varepsilon$ 一致成立。

由上例可以看出，尽管级数的各项都具有良好的解析性质，并且级数在点集 $E$ 上处处收敛，但其和函数仍然可能在 $E$ 上不连续，由此可见一致收敛的重要性。

和函数在点 $z_0\in E$ 连续就是 $\lim_{z\to z_0}f(z)=f(z_0)$，也就是

$$
\lim_{z\to z_0}\sum_{k=1}^{\infty}f_k(z)
=\sum_{k=1}^{\infty}f_k(z_0)
=\sum_{k=1}^{\infty}\lim_{z\to z_0}f_k(z).
$$

所以，可以等价地说，取极限与求和可以交换次序，其前提是级数要一致收敛。

**定理（逐项积分）**　若级数（11）的各项均在曲线 $C$ 上连续，且级数在 $C$ 上一致收敛于函数 $f(z)$，则可以沿 $C$ 逐项积分，即

$$
\int_C f(z)\,dz=\sum_{k=1}^{\infty}\int_C f_k(z)\,dz. \tag{13}
$$

**证明**　记曲线 $C$ 的长度为 $L$。已知级数在 $C$ 上一致收敛于函数 $f(z)$，故 $\forall\varepsilon>0$，$\exists n_0\in\mathbb N^+$，当 $n>n_0$ 时，

$$
|s_n(z)-f(z)|<\varepsilon/L
$$

在 $C$ 上一致成立，于是，当 $n>n_0$ 时，

$$
\left|\sum_{k=1}^n\int_C f_k(z)\,dz-\int_C f(z)\,dz\right|
=\left|\int_C[s_n(z)-f(z)]\,dz\right|
\le\int_C|s_n(z)-f(z)||dz|
<\frac\varepsilon L\int_C|dz|=\varepsilon.
$$

证毕。

式（13）可以写作

$$
\sum_{k=1}^{\infty}\int_C f_k(z)\,dz
=\int_C\sum_{k=1}^{\infty}f_k(z)\,dz. \tag{13}
$$

亦即积分与求和可以交换次序，其前提是级数要一致收敛。

很容易想到的一个问题是，什么情况下求导与求和可以交换次序？当然，只有在级数各项均可导、或者说各项均解析的条件下，这一问题才有意义。所以这是一个关于解析函数项级数的问题，下一小节就会给出答案。

### §2.3　解析函数项级数

若级数（11）的各项均为解析函数，就称为解析函数项级数。关于解析函数项级数，有下面的重要定理。

**定理（Weierstrass）**　若级数（11）的各项均在区域 $D$ 内解析，且级数在 $D$ 内的任一有界闭集上一致收敛（称为在 $D$ 上内闭一致收敛）于函数 $f(z)$，则

1. 和函数 $f(z)$ 在 $D$ 内解析；
2. 级数 $\sum_{k=1}^{\infty}f_k^{(p)}(z)$ 在 $D$ 上内闭一致收敛到 $f^{(p)}(z)$，其中 $p\in\mathbb N^+$。

**注**　① 这是一个非常重要的定理，它是级数理论的基石，所以复变函数论中的级数理论也称为 Weierstrass 级数理论。② 这个定理是微积分中所没有的。在微积分中，要求级数 $\sum_{k=1}^{\infty}f_k(x)$ 在区间 $[a,b]$ 上收敛，各项都具有连续导数，且 $\sum_{k=1}^{\infty}f_k'(x)$ 在 $[a,b]$ 上一致收敛，才能逐项求导一次，并且不能保证可以进一步逐项求导。而这里只要求原级数一致收敛，且各项解析，那么就可以不断地逐项求导，且求导后的级数也一致收敛到和函数的同阶导数。③ 这里要求级数在 $D$ 上内闭一致收敛，这比要求级数在 $D$ 内一致收敛要弱。比如级数 $\sum_{k=1}^{\infty}z^k$ 虽然在单位圆 $|z|<1$ 上不能一致收敛，但内闭一致收敛（因为它在闭圆 $|z|\le r<1$ 上一致收敛，参看上面的例 3），所以其和函数 $f(z)$ 在单位圆 $|z|<1$ 内解析，且可以逐项求导。我们知道，其和函数是 $z/(1-z)$，它确实在单位圆 $|z|<1$ 内解析。

为了证明这一定理，需要引用 Morera 定理，后者是 Cauchy 积分定理的逆定理。这一定理本来应该放在第二章介绍，但因为在本书的其它章节不需要用到，所以放在此处，作为证明 Weierstrass 定理的预备。

**定理（Morera）**　设函数 $f(z)$ 在单通区域 $D$ 内连续，且对 $D$ 内任一围线 $C$，满足 $\int_Cf(z)\,dz=0$，则 $f(z)$ 在 $D$ 内解析。

**证明**　取固定点 $z_0\in D$，由于对 $D$ 内任一围线 $C$ 有 $\int_Cf(z)\,dz=0$，故积分 $\int_{z_0}^z f(\zeta)\,d\zeta$（其中 $z\in D$）只与上限 $z$ 有关，而与 $z_0$ 到 $z$ 的路径无关。定义 $F(z)=\int_{z_0}^z f(\zeta)\,d\zeta$，它是 $z$ 的单值函数。我们在第二章已经证明，只要 $f(z)$ 在 $D$ 内连续，就有 $F'(z)=f(z)$，亦即 $F(z)$ 解析。由于解析函数的导数也是解析函数，故 $f(z)$ 在 $D$ 内解析。证毕。

下面就给出 Weierstrass 定理的

**证明**　（1）$\forall z\in D$，总可以找到邻域 $N(z,r)\subset D$，按已知条件，级数（11）在 $N(z,r)$ 内的任一围线 $C$ 上一致收敛，故

$$
\int_Cf(z)\,dz=\sum_{k=1}^{\infty}\int_Cf_k(z)\,dz=0,
$$

其中最后一个等号用了 Cauchy 积分定理，因为级数各项在 $N(z,r)$ 内解析，而 $N(z,r)$ 为单通区域。由上式和 Morera 定理，易知 $f(z)$ 在 $N(z,r)$ 内解析，特别地，$f(z)$ 在点 $z$ 解析，由 $z$ 的任意性，可知 $f(z)$ 在 $D$ 内解析。

（2）$\forall z_0\in D$，总可以找到邻域 $N(z_0,2r)$，使得 $\overline{N(z_0,2r)}\subset D$，按已知条件，级数（11）在 $\overline{N(z_0,2r)}$ 上一致收敛，故 $\forall\varepsilon>0$，$\exists n_0\in\mathbb N^+$，当 $n>n_0$ 时，

$$
\left|\sum_{k=1}^{n}f_k(\zeta)-f(\zeta)\right|<\varepsilon
$$

在圆周 $C:|\zeta-z_0|=2r$ 上一致成立。今 $\forall z\in N(z_0,r)$，$\zeta\in C$，有 $|\zeta-z_0|=2r$，$|z-z_0|<r$，故

$$
|\zeta-z|=|(\zeta-z_0)-(z-z_0)|\ge|\zeta-z_0|-|z-z_0|>2r-r=r,
$$

故

$$
\left|\sum_{k=1}^{n}\frac{f_k(\zeta)}{(\zeta-z)^{p+1}}-\frac{f(\zeta)}{(\zeta-z)^{p+1}}\right|
=\frac{\left|\sum_{k=1}^{n}f_k(\zeta)-f(\zeta)\right|}{|\zeta-z|^{p+1}}
<\frac\varepsilon{r^{p+1}},\qquad p\in\mathbb N^+,
$$

对 $\zeta\in C$ 和 $z\in N(z_0,r)$ 一致成立。由于 $f(z)$ 和 $f_k(z)$ 均在 $\overline{N(z_0,2r)}$ 上解析，由 Cauchy 高阶导数公式，当 $z\in N(z_0,2r)$，有

$$
f_k^{(p)}(z)=\frac{p!}{2\pi i}\int_C\frac{f_k(\zeta)}{(\zeta-z)^{p+1}}\,d\zeta,
\qquad
f^{(p)}(z)=\frac{p!}{2\pi i}\int_C\frac{f(\zeta)}{(\zeta-z)^{p+1}}\,d\zeta.
$$

当 $z\in N(z_0,r)$，上式当然也成立，于是当 $n>n_0$ 时，有

$$
\left|\sum_{k=1}^{n}f_k^{(p)}(z)-f^{(p)}(z)\right|
=\frac{p!}{2\pi}\left|\sum_{k=1}^{n}\int_C\frac{f_k(\zeta)}{(\zeta-z)^{p+1}}\,d\zeta-\int_C\frac{f(\zeta)}{(\zeta-z)^{p+1}}\,d\zeta\right|
$$

$$
\le\frac{p!}{2\pi}\int_C\left|\sum_{k=1}^{n}\frac{f_k(\zeta)}{(\zeta-z)^{p+1}}-\frac{f(\zeta)}{(\zeta-z)^{p+1}}\right||d\zeta|
<\frac{p!}{2\pi}\frac\varepsilon{r^{p+1}}\int_C|d\zeta|
=\frac{p!}{2\pi}\frac\varepsilon{r^{p+1}}4\pi r
=\frac{2p!}{r^p}\varepsilon.
$$

上式在 $z\in N(z_0,r)$ 时一致成立，换句话说，级数 $\sum_{k=1}^{\infty}f_k^{(p)}(z)$ 在 $N(z_0,r)$ 上一致收敛到 $f^{(p)}(z)$，而 $z_0\in D$ 是任意的，所以对 $D$ 内任一闭集 $G$ 上的任一点，都存在该点的一个邻域，使 $\sum_{k=1}^{\infty}f_k^{(p)}(z)$ 在该邻域上一致收敛到 $f^{(p)}(z)$。这些邻域构成了闭集 $G$ 的一个开覆盖（就是说这些邻域（开集）覆盖了闭集 $G$），我们总可以从中选出有限个覆盖闭集 $G$（这里用到了一个定理，称为有限覆盖定理，又称 Heine–Borel 定理，不过，这一结论是比较直观的），由于 $\sum_{k=1}^{\infty}f_k^{(p)}(z)$ 在其中的每一个上都一致收敛到 $f^{(p)}(z)$，故在闭集 $G$ 上亦然。证毕。

### §2.4　\*Riemann $\zeta$ 函数、Riemann 假设与千禧数学难题

下面用 Weierstrass 定理来研究一个实例。

**例 7**　研究级数 $\sum_{k=1}^{\infty}1/k^z=\sum_{k=1}^{\infty}e^{-z\ln k}$。

显然，该级数的每一项都是复平面上的解析函数。只要找到级数的收敛区域，应用 Weierstrass 定理，其和函数的解析性也就清楚了。从级数的形式来看，我们猜想级数的敛散性主要取决于 $z$ 的实部。又已经知道，当 $z=x>1$ 时级数收敛，而 $z=x\le1$ 时级数发散，所以估计级数的收敛区域应该是 $\operatorname{Re}z>1$。

首先证明，如果级数在 $z_0$ 收敛，则当 $\operatorname{Re}z>\operatorname{Re}z_0$ 时亦收敛。事实上，由于级数在 $z_0$ 收敛，故 $\exists M>0$，使得 $\forall n\in\mathbb N^+$，均有 $|s_n(z_0)|<M$。显然

$$
e^{-z\ln k}=e^{-(z-z_0)\ln k}e^{-z_0\ln k}
=e^{-(z-z_0)\ln k}[s_k(z_0)-s_{k-1}(z_0)],\qquad k\in\mathbb N^+,
$$

其中 $s_0(z_0)\equiv0$，故有

$$
\sum_{k=n+1}^{n+p}e^{-z\ln k}
=e^{-(z-z_0)\ln(n+p)}s_{n+p}(z_0)
-e^{-(z-z_0)\ln(n+1)}s_n(z_0)
$$

$$
+\sum_{k=n+1}^{n+p-1}(z-z_0)s_k(z_0)\int_k^{k+1}\frac{e^{-(z-z_0)\ln t}}{t}\,dt,
$$

于是

$$
\left|\sum_{k=n+1}^{n+p}e^{-z\ln k}\right|
\le M[e^{-(x-x_0)\ln(n+p)}+e^{-(x-x_0)\ln(n+1)}]
+M|z-z_0|\sum_{k=n+1}^{n+p-1}\int_k^{k+1}\frac{e^{-(x-x_0)\ln t}}{t}\,dt
$$

$$
<2Me^{-(x-x_0)\ln n}
+M|z-z_0|\sum_{k=n+1}^{n+p-1}\int_k^{k+1}\frac{e^{-(x-x_0)\ln t}}{t}\,dt,
$$

其中 $x=\operatorname{Re}z$，$x_0=\operatorname{Re}z_0$。但

$$
\sum_{k=n+1}^{n+p-1}\int_k^{k+1}\frac{e^{-(x-x_0)\ln t}}{t}\,dt
=\int_{n+1}^{n+p}\frac{e^{-(x-x_0)\ln t}}{t}\,dt
$$

$$
=\frac{e^{-(x-x_0)\ln(n+1)}-e^{-(x-x_0)\ln(n+p)}}{x-x_0}
<\frac{e^{-(x-x_0)\ln n}}{x-x_0},
$$

于是

$$
\left|\sum_{k=n+1}^{n+p}e^{-z\ln k}\right|
<M\left(\frac{|z-z_0|}{x-x_0}+2\right)e^{-(x-x_0)\ln n}
=M\left(\frac{|z-z_0|}{x-x_0}+2\right)\frac1{n^{x-x_0}}.
$$

由于 $x>x_0$，容易看出，只要 $n$ 充分大，上式就可以任意小。事实上，$\forall\varepsilon>0$，只要取

$$
n>\left[\frac{M}{\varepsilon}\left(\frac{|z-z_0|}{x-x_0}+2\right)\right]^{1/(x-x_0)},
$$

则有

$$
\left|\sum_{k=n+1}^{n+p}e^{-z\ln k}\right|<\varepsilon,\qquad \forall p\in\mathbb N^+,
$$

由 Cauchy 收敛原理，级数当 $\operatorname{Re}z>\operatorname{Re}z_0$ 时收敛。但注意以上的 $n_0$ 不仅依赖于 $\varepsilon$，而且也依赖于 $z$。

我们已经知道，当 $z=1$ 时级数发散，由上述结论容易推知，级数当 $\operatorname{Re}z<1$ 时必然发散。

另一方面，我们也知道，当 $s>1$ 时，正项级数 $\sum_{k=1}^{\infty}1/k^s$ 收敛，而当 $\operatorname{Re}z\ge s$ 时，有 $|1/k^z|=1/k^x\le1/k^s$，由 Weierstrass M-判别法，级数当 $\operatorname{Re}z\ge s$ 时绝对且一致收敛。由此容易推知，级数在区域 $\operatorname{Re}z>1$ 内的任一有界闭集上一致收敛，亦即在区域 $\operatorname{Re}z>1$ 上内闭一致收敛，由 Weierstrass 定理，和函数在区域 $\operatorname{Re}z>1$ 内解析，并且逐项求导后的级数也在区域 $\operatorname{Re}z>1$ 上内闭一致收敛到和函数的同阶导数。

但是用上面的方法无法判断级数在 $z=1+iy$ 而 $y\ne0$ 时的敛散性。

由上面的级数定义的函数可以解析开拓到整个复平面上，所得的函数称为 Riemann $\zeta$ 函数，记作 $\zeta(z)$，$z=1$ 是它的唯一一个奇点，$z=-2n$（$n\in\mathbb N^+$）是它的简单零点（即一阶零点），其它的零点都位于带形区域 $0\le\operatorname{Re}z\le1$ 内。

Riemann 假设（Riemann’s hypothesis）：$\zeta(z)$ 在带形区域 $0\le\operatorname{Re}z\le1$ 内的所有零点都位于直线 $\operatorname{Re}z=1/2$ 上。

已经知道有一个不可数的零点集合位于该直线上。但 Riemann 假设的严格证明目前尚未完成。

1900 年，著名的数学家 D. Hilbert 在巴黎的世界数学家大会上提出了 23 个数学问题，认为对这些问题的研究会极大地影响二十世纪数学的发展，Riemann 假设就是这些问题之一。这些问题有的很快就得到解决，有的则直到目前尚未解决，Riemann 假设就是尚未解决的问题之一。关于 Hilbert 问题，参看 https://mathcs.clarku.edu/~djoyce/hilbert/toc.html。

2000 年，美国的 Clay 数学研究所在尚未解决的著名数学问题中选出七大难题，并设立了每题 100 万美元的奖金，Riemann 假设就在其中，参看 http://www.claymath.org/millennium-problems/。

## §3　幂级数

### §3.1　幂级数的敛散性

形如

$$
\sum_{n=0}^{\infty}c_n(z-a)^n\qquad (a\in\mathbb C,c_n\in\mathbb C) \tag{14}
$$

的级数称为幂级数。

**注**　从本节开始，我们改用 $n$ 来表示求和指标。并请注意现在求和多了 $n=0$ 一项，不过这对于级数的性质没有什么影响。

显然，级数（14）的每一项都是解析函数，事实上，幂级数是最简单的解析函数项级数。根据 Weierstrass 定理，只要把级数的收敛区域研究清楚，那么和函数的解析性质也就清楚了。关于这一幂级数的敛散性，Abel 给出了下述的重要定理。

**定理（Abel）**　如果级数（14）在某点 $z_0\ne a$ 收敛，则（1）它在圆 $K:|z-a|<|z_0-a|$ 内绝对收敛；（2）它在闭圆 $\overline K_\rho:|z-a|\le\rho<|z_0-a|$ 内一致收敛（即在圆 $K$ 上内闭一致收敛）。

**证明**　（1）由已知条件，级数 $\sum_{n=0}^{\infty}c_n(z_0-a)^n$ 收敛，故 $\exists M>0$，使得

$$
|c_n(z_0-a)^n|<M,\qquad n=0,1,2,\cdots,
$$

于是

$$
|c_n(z-a)^n|
=|c_n(z_0-a)^n|\left|\frac{z-a}{z_0-a}\right|^n
<M\left|\frac{z-a}{z_0-a}\right|^n,
$$

当 $|z-a|<|z_0-a|$，以右边为通项的级数收敛，由比较判别法，以左边为通项的级数亦收敛，即原级数在圆 $K$ 上绝对收敛。

（2）当 $|z-a|\le\rho<|z_0-a|$，

$$
|c_n(z-a)^n|<M\left|\frac{z-a}{z_0-a}\right|^n
\le M\left(\frac{\rho}{|z_0-a|}\right)^n,
$$

右边是与 $z$ 无关的正数，且以此为通项的级数显然收敛，由 Weierstrass M-判别法，原级数在闭圆 $\overline K_\rho$ 上一致收敛。证毕。

由 Abel 定理，立得以下

**推论 1**　如果级数（14）在某点 $z_0\ne a$ 发散，则它必在圆周 $|z-a|=|z_0-a|$ 外发散。

由 Abel 定理和推论 1 可知，级数（14）在某点 $z$ 收敛或发散只与该点到 $a$ 的距离有关，而与该点相对于 $a$ 的方向无关（即只与 $z-a$ 的模有关，而与其辐角无关），因此收敛点集与发散点集的分界必定是一个圆周，这就是下面的

**推论 2**　对于级数（14），$\exists R\ge0$，使它在圆 $|z-a|<R$ 上绝对收敛且内闭一致收敛，而在圆外处处发散。

上面的 $R$ 称为收敛半径，$|z-a|<R$ 称为收敛圆，而 $|z-a|=R$ 称为收敛圆周。

应当指出，对于级数（14）在收敛圆周 $|z-a|=R$ 上的敛散性，Abel 定理及其推论没有给出结论。

### §3.2　收敛半径的求法

上一小节的推论 2 可以换一种等价的方式来表达：级数 $\sum_{n=0}^{\infty}|c_n|r^n$ 当 $0\le r<R$ 时收敛，而当 $r>R$ 时发散。也就是说，$R$ 是该正项级数的收敛区间与发散区间的分界点。由正项级数收敛的 Cauchy 判别法和 d’Alembert 判别法，不难推出下面收敛半径的计算公式。

**定理（收敛半径的 Cauchy 计算公式）**　如果 $\lim_{n\to\infty}\sqrt[n]{|c_n|}=l$，则级数（14）的收敛半径为

$$
R=
\begin{cases}
1/l,&\text{当 }0<l<\infty;\\
0,&\text{当 }l=\infty;\\
\infty,&\text{当 }l=0.
\end{cases} \tag{15}
$$

实际上，即使 $\lim_{n\to\infty}\sqrt[n]{|c_n|}$ 不存在，只要 $\varlimsup_{n\to\infty}\sqrt[n]{|c_n|}=l$，式（15）仍然成立（相应的判别法也称为 Cauchy–Hadamard 判别法）。这里 $\varlimsup$ 表示上极限。一个实数列的上极限就是它的所有收敛子列的极限中最大的一个（虽然这不是上极限的原始定义）。

**例 1**　级数 $\sum_{n=0}^{\infty}[2+(-)^n]^n z^n$，$\sqrt[n]{|c_n|}=2+(-)^n$，显然 $\lim_{n\to\infty}\sqrt[n]{|c_n|}$ 不存在，但考虑子列 $\{c_{2k}\}_{k=0}^{\infty}$，易得 $\varlimsup_{n\to\infty}\sqrt[n]{|c_n|}=3$。因此，该级数的收敛半径为 $R=1/3$，即收敛圆为 $|z|<1/3$。

**例 2**　级数 $\sum_{k=0}^{\infty}z^{k^2}=1+z+z^4+z^9+z^{16}+\cdots$，当 $n=k^2$（$k=0,1,2,\cdots$），$c_n=1$，其它情况下 $c_n=0$，显然 $\lim_{n\to\infty}\sqrt[n]{|c_n|}$ 不存在，但考虑子列 $\{c_{k^2}\}_{k=0}^{\infty}$，易得 $\varlimsup_{n\to\infty}\sqrt[n]{|c_n|}=1$。因此，该级数的收敛半径为 $R=1$，即收敛圆为 $|z|<1$。

**定理（收敛半径的 d’Alembert 计算公式）**　如果 $\lim_{n\to\infty}|c_{n+1}/c_n|=l$，则级数（14）的收敛半径由式（15）给出。

如果极限 $\lim_{n\to\infty}|c_{n+1}/c_n|$ 不存在，则无法使用 d’Alembert 计算公式，所以它不如 Cauchy–Hadamard 计算公式来得普遍，甚至也不如 Cauchy 计算公式，但是，当极限 $\lim_{n\to\infty}|c_{n+1}/c_n|$ 存在时，d’Alembert 计算公式通常是比较方便的。

**例 3**　考虑级数 $\sum_{n=2}^{\infty}c_nz^n$，其中 $c_{2k}=2k$，$c_{2k+1}=k$（$k\in\mathbb N^+$）。易得 $\lim_{k\to\infty}|c_{2k+1}/c_{2k}|=1/2$，而 $\lim_{k\to\infty}|c_{2k+2}/c_{2k+1}|=2$，即 $\lim_{n\to\infty}|c_{n+1}/c_n|$ 不存在，根据 d’Alembert 判别法，易知当 $|z|<1/(\varlimsup_{n\to\infty}|c_{n+1}/c_n|)=1/2$ 时级数收敛，而当 $|z|>1/(\varliminf_{n\to\infty}|c_{n+1}/c_n|)=2$ 时级数发散，但是对于 $1/2<|z|<2$ 时级数的敛散性，我们无法由 d’Alembert 判别法得到进一步的信息。另一方面，不难求得 $\lim_{n\to\infty}\sqrt[n]{|c_n|}=1$，故由 Cauchy 计算公式可得该级数的收敛半径为 $R=1$，即收敛圆为 $|z|<1$。可见 Cauchy 计算公式更普遍一些。

### §3.3　和函数的解析性

知道了幂级数的收敛范围，由前面的有关定理，就可以得到以下结论：

**定理**　幂级数（14）在收敛圆 $K:|z-a|<R$ 上有如下性质：

1. 和函数 $f(z)$ 解析；
2. 可以逐项求导而收敛半径不变；
3. 可以逐项积分而收敛半径不变。

**证明**　考虑闭圆 $\overline K_\rho:|z-a|\le\rho<R$，总可以找到一点 $z_0$，使得 $\rho<|z_0-a|<R$。如此，$z_0$ 在收敛圆内，故级数（14）在 $z_0$ 收敛，由 Abel 定理，它必在 $\overline K_\rho$ 上一致收敛。也就是说，它在收敛圆 $K$ 上内闭一致收敛。显然，幂级数的每一项都是解析函数，由 Weierstrass 定理，和函数在 $K$ 内解析。

由 Weierstrass 定理，级数（14）逐项求导后的级数在 $K$ 上内闭一致收敛到 $f'(z)$，从而该级数在 $|z-a|<R$ 收敛，即

$$
\sum_{n=0}^{\infty}(n+1)c_{n+1}(z-a)^n=f'(z),\qquad |z-a|<R. \tag{16}
$$

这就是说，原级数可以逐项求导，而所得级数的收敛半径 $R_d\ge R$。

对于 $K$ 内的任何曲线 $C$，总可以找到 $\overline K_\rho:|z-a|\le\rho<R$，使得曲线 $C$ 在 $\overline K_\rho$ 内，由于级数（14）在 $\overline K_\rho$ 内一致收敛，故在 $C$ 上一致收敛，因此可以逐项积分。如果 $C$ 是由固定点 $b$ 到变点 $z$ 的曲线（$b$ 和 $z$ 均在收敛圆内），则逐项积分后所得级数仍为 $z-a$ 的幂级数，这就有收敛半径的问题。由于 $(z-a)^n$ 的积分与路径无关，且 $\int_b^z(\zeta-a)^n\,d\zeta=\int_a^z(\zeta-a)^n\,d\zeta-\int_a^b(\zeta-a)^n\,d\zeta$，后项所构成的级数是常数项级数，且收敛到 $\int_a^b f(z)\,dz$，即

$$
\sum_{n=0}^{\infty}\frac{c_n}{n+1}(b-a)^{n+1}=\int_a^b f(z)\,dz,
$$

故只需考虑前项所构成的级数，或者说，可以取 $b=a$ 而不失一般性，这时有

$$
\sum_{n=0}^{\infty}\frac{c_n}{n+1}(z-a)^{n+1}
=\int_a^z f(\zeta)\,d\zeta,\qquad |z-a|<R. \tag{17}
$$

这就是说，原级数可以逐项积分，而所得级数的收敛半径 $R_i\ge R$。

下面证明 $R_i=R_d=R$。如果对式（17）逐项求导，则所得级数的收敛半径为 $(R_i)_d\ge R_i$，但该级数即为级数（14），故 $(R_i)_d=R$，从而 $R\ge R_i$，但上面得到 $R_i\ge R$，两式综合即得 $R_i=R$。类似可证 $R_d=R$。证毕。

**习题**　求下列幂级数的收敛圆：（1）$\sum_{n=1}^{\infty}n^{\ln n}(z-2)^n$；（2）$\sum_{n=1}^{\infty}(z/n)^n$；（3）$\sum_{n=1}^{\infty}n!(z/n)^n$。

## §4　解析函数的 Taylor 展开

### §4.1　Taylor 定理

上节已经看到，一个幂级数在其收敛圆内具有解析的和函数，换句话说，它在收敛圆内代表一个解析函数。反过来，在圆内解析的函数是否可以展开为幂级数呢？下面的定理给出了肯定的答案。

**定理（Taylor）**　设函数 $f(z)$ 在圆 $K:|z-a|<R$ 内解析，则在 $K$ 内可以展开为幂级数：

$$
f(z)=\sum_{n=0}^{\infty}c_n(z-a)^n, \tag{18}
$$

其中

$$
c_n=\frac1{2\pi i}\int_{\Gamma_\rho}\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta
=\frac{f^{(n)}(a)}{n!},\qquad
\Gamma_\rho:|\zeta-a|=\rho<R,\quad n=0,1,2,\cdots, \tag{19}
$$

且展开式是唯一的。

**证明**　$\forall z\in K$，总可以找到正数 $\rho$ 使得 $|z-a|<\rho<R$，作圆周 $\Gamma_\rho:|\zeta-a|=\rho$，则 $z$ 在其内部。显然，$f(\zeta)$ 在闭圆 $\overline K_\rho:|\zeta-a|\le\rho$ 上解析，根据 Cauchy 积分公式，

$$
f(z)=\frac1{2\pi i}\int_{\Gamma_\rho}\frac{f(\zeta)}{\zeta-z}\,d\zeta.
$$

今利用恒等式（可用初等数学的方法求部分和并取极限得到）

$$
\frac1{1-u}=\sum_{n=0}^{\infty}u^n,\qquad |u|<1,
$$

则在 $\Gamma_\rho$ 上成立

$$
\frac{f(\zeta)}{\zeta-z}
=\frac{f(\zeta)}{(\zeta-a)-(z-a)}
=\frac{f(\zeta)}{\zeta-a}\frac1{1-(z-a)/(\zeta-a)}
=\frac{f(\zeta)}{\zeta-a}\sum_{n=0}^{\infty}\left(\frac{z-a}{\zeta-a}\right)^n. \tag{20}
$$

当 $\zeta\in\Gamma_\rho$，

$$
\left|\left(\frac{z-a}{\zeta-a}\right)^n\right|=\left(\frac{|z-a|}{\rho}\right)^n,
$$

右边与 $\zeta$ 无关，且以此为通项的级数显然收敛，故由 Weierstrass M-判别法，级数 $\sum_{n=0}^{\infty}[(z-a)/(\zeta-a)]^n$ 在 $\Gamma_\rho$ 上一致收敛。又函数 $f(\zeta)/(\zeta-a)$ 在 $\Gamma_\rho$ 上显然有界，因此式（20）右边的级数在 $\Gamma_\rho$ 上一致收敛，因而可以逐项积分，得

$$
f(z)=\frac1{2\pi i}\int_{\Gamma_\rho}\frac{f(\zeta)}{\zeta-z}\,d\zeta
=\frac1{2\pi i}\int_{\Gamma_\rho}\frac{f(\zeta)}{\zeta-a}\sum_{n=0}^{\infty}\left(\frac{z-a}{\zeta-a}\right)^n d\zeta
$$

$$
=\sum_{n=0}^{\infty}\left[\frac1{2\pi i}\int_{\Gamma_\rho}\frac{f(\zeta)}{(\zeta-a)^{n+1}}\,d\zeta\right](z-a)^n
=\sum_{n=0}^{\infty}\frac{f^{(n)}(a)}{n!}(z-a)^n,
$$

其中最后一步用到了高阶导数的 Cauchy 公式。上式即所求证。虽然开始时取 $\rho$ 满足 $|z-a|<\rho<R$，但最后出现在系数中的积分显然与 $\rho$ 的大小无关，只要 $\rho<R$ 即可。

为了证明展开式的唯一性，假设用某种方法找到在 $K$ 内成立的另一个展开式

$$
f(z)=\sum_{n=0}^{\infty}c_n'(z-a)^n,
$$

令 $z=a$，即得 $c_0'=f(a)=c_0$，又因为幂级数逐项求导后在 $K$ 内仍收敛，故可逐项求导 $n$ 次后令 $z=a$，而得 $c_n'=f^{(n)}(a)/n!=c_n$，$n=1,2,\cdots$，即展开式是唯一的。证毕。

上面的式（18）称为 $f(z)$ 在 $a$ 点的 Taylor 展开式，其右边称为 Taylor 级数，式（19）称为 Taylor 系数，而 $a$ 点称为展开中心。

因为只要 $|z-a|<R$，式（18）就成立，故该 Taylor 级数的收敛半径大于或等于 $R$。

### §4.2　展开式的收敛半径

如果求得了 Taylor 展开式的通项系数，就可以用上节的方法来计算该 Taylor 级数的收敛半径。然而，在 $f(z)$ 比较复杂的情况下，要求出展开式的通项系数是很困难的，因此，我们希望有简单的方法可以由被展开函数本身的性质来计算收敛半径。即使可以求出展开式的通项系数，这样的方法无疑也是有益的，因为两种方法求得的收敛半径可以互相印证。关于这一问题，我们有以下简单的结论。

以 $a$ 点为中心将解析函数 $f(z)$ 展开为 Taylor 级数，若 $f(z)$ 距离 $a$ 点最近的奇点为 $z_0$，则展开式的收敛半径为 $R=|z_0-a|$。

事实上，以 $a$ 点为中心作圆 $K$，只要 $f(z)$ 在其中解析，则 Taylor 展开式就在 $K$ 内成立，因此，展开圆的圆周至少可以推至 $z_0$ 处，从而收敛半径 $R\ge|z_0-a|$。

另一方面，以 $a$ 点为中心将 $f(z)$ 展开为 Taylor 级数，假设其收敛半径 $R>|z_0-a|$，由于幂级数的和函数在收敛圆内解析，故和函数在 $z_0$ 解析，这与 $z_0$ 是 $f(z)$ 的奇点这一已知条件矛盾，故 $R\le|z_0-a|$。

综上所述，即得 $R=|z_0-a|$。

应该注意，以上说的是展开式成立的半径，即 Taylor 级数收敛于 $f(z)$（收敛且等于 $f(z)$）的收敛半径，至于 Taylor 级数本身的收敛半径（收敛但不必等于 $f(z)$）则可以更大。例如函数

$$
f(z)=
\begin{cases}
1+z,&\text{当 }z\ne1;\\
0,&\text{当 }z=1
\end{cases}
$$

以 $z=1$ 为奇点。若以 $a=0$ 为中心作 Taylor 展开，则该 Taylor 级数收敛于 $f(z)$ 的收敛半径为 1。容易求得

$$
f(0)=1,\qquad f'(0)=1,\qquad f^{(n)}(0)=0\quad(n\ge2),
$$

故该 Taylor 展开式为

$$
f(z)=1+z,\qquad |z|<1.
$$

由于右边的 Taylor 级数为多项式，故其本身的收敛半径为 $\infty$，这与展开式成立的半径是不同的。

不过，如果 $\lim_{z\to z_0}f(z)=\infty$，则 Taylor 级数本身的收敛半径也只能是 $|z_0-a|$。事实上，设 Taylor 级数本身的收敛半径为 $R'>R=|z_0-a|$，该级数在圆 $K':|z-a|<R'$ 上有解析的和函数，记作 $F(z)$，而在圆 $K:|z-a|<R$ 上，$F(z)=f(z)$。令 $F(z)$ 在 $z_0$ 解析，故在 $z_0$ 连续，于是 $\lim_{z\to z_0}F(z)=F(z_0)\ne\infty$，另一方面，$\lim_{z\to z_0,z\in K}F(z)=\lim_{z\to z_0}f(z)=\infty$，矛盾。

### §4.3　初等单值函数的 Taylor 展开式

本小节给出若干常用的初等单值函数的 Taylor 展开式。

**例 1**　求 $f(z)=e^z$ 在 $a=0$ 处的 Taylor 展开式。

**解**　由于 $f^{(n)}(z)=e^z$（$n\in\mathbb N^+$），故 $f(0)=1$，$f^{(n)}(0)=1$（$n\in\mathbb N^+$），于是所求 Taylor 展开式为

$$
e^z=1+z+\frac{z^2}{2!}+\cdots+\frac{z^n}{n!}+\cdots
=\sum_{n=0}^{\infty}\frac{z^n}{n!}. \tag{21}
$$

由于 $e^z$ 在 $z$ 平面上解析，即在有限远处无奇点，故展开式的收敛半径应为 $R=\infty$，由 Taylor 系数计算收敛半径可得同样的结果。

**例 2**　求 $f(z)=1/(1-z)$ 在 $a=0$ 处的 Taylor 展开式。

**解**　$f(0)=1$，$c_0=1$，而

$$
f^{(n)}(z)=\frac{n!}{(1-z)^{n+1}},\qquad f^{(n)}(0)=n!,\qquad c_n=1,\quad n\in\mathbb N^+,
$$

故所求 Taylor 展开式为

$$
\frac1{1-z}=1+z+\cdots+z^n+\cdots=\sum_{n=0}^{\infty}z^n. \tag{22}
$$

由于 $f(z)$ 在 $z$ 平面上有奇点 $z_0=1$，而展开中心为 $a=0$，故收敛半径应为 $R=1$，由 Taylor 系数计算收敛半径可得同样的结果。这一展开式已经在前面当作恒等式用过，当 $|z|<1$ 时，可用初等数学方法求出右边级数的部分和，取极限即得左边的和函数。

上面两个结果经常被用来计算更复杂的函数的 Taylor 展开式，所以应该熟练掌握它们。

**例 3**　求 $\cos z$ 和 $\sin z$ 在 $a=0$ 处的 Taylor 展开式。

**解**　与上面两题类似，本题也可用系数公式计算展开系数而得到 Taylor 展开式。这里我们采用另一种方法。由例 1 的结果，可有

$$
e^{iz}=\sum_{n=0}^{\infty}\frac{(iz)^n}{n!},\qquad
e^{-iz}=\sum_{n=0}^{\infty}\frac{(-iz)^n}{n!},
$$

于是

$$
e^{iz}+e^{-iz}=2\sum_{k=0}^{\infty}\frac{(iz)^{2k}}{(2k)!}
=2\sum_{k=0}^{\infty}\frac{(-)^kz^{2k}}{(2k)!},
$$

$$
e^{iz}-e^{-iz}=2\sum_{k=0}^{\infty}\frac{(iz)^{2k+1}}{(2k+1)!}
=2i\sum_{k=0}^{\infty}\frac{(-)^kz^{2k+1}}{(2k+1)!},
$$

故所求 Taylor 展开式为

$$
\cos z=\sum_{k=0}^{\infty}\frac{(-)^kz^{2k}}{(2k)!},\qquad
\sin z=\sum_{k=0}^{\infty}\frac{(-)^kz^{2k+1}}{(2k+1)!}. \tag{23}
$$

由于 $\cos z$ 和 $\sin z$ 在 $z$ 平面上解析，即在有限远处无奇点，故两式的收敛半径均应为 $R=\infty$，这从上面的推导过程也可以看出。由 Taylor 系数计算收敛半径可得同样的结果。

### §4.4　初等多值函数的 Taylor 展开式

对于多值函数，首先要分出单值分支，然后才能作 Taylor 展开。

**例 4**　求 $f(z)=\ln(1+z)$ 在 $a=0$ 处的 Taylor 展开式，规定 $f(0)=0$。

**解**　$f(z)=\ln(1+z)$ 的支点为 $-1$ 和 $\infty$，将 $z$ 平面沿负实轴从 $-1$ 到 $\infty$ 割破，可以分出单值分支。规定 $f(0)=0$ 即规定 $\arg(1+z)|_{z=0}=\arg[z-(-1)]|_{z=0}=0$，就确定了单值分支。（为了使 Taylor 级数的收敛圆尽量大，割线不应该穿过单位圆 $|z|<1$。）按规定，$f(0)=0$，故 $c_0=0$。又

$$
f^{(n)}(z)=\frac{(-)^{n-1}(n-1)!}{(1+z)^n},\qquad
f^{(n)}(0)=(-)^{n-1}(n-1)!,\qquad
c_n=\frac{(-)^{n-1}}{n},\quad n\in\mathbb N^+,
$$

故所求 Taylor 展开式为

$$
\ln(1+z)=\sum_{n=1}^{\infty}\frac{(-)^{n-1}}{n}z^n. \tag{24}
$$

由于 $f(z)$ 距离展开中心 $a=0$ 最近的奇点为 $z_0=-1$，故展开式的收敛半径应为 $R=1$，由 Taylor 系数计算收敛半径可得同样的结果。

如果规定 $\arg(1+z)|_{z=0}=\arg[z-(-1)]|_{z=0}=2k\pi$（$k\in\mathbb Z$），则 $f(0)=i2k\pi$，$c_0=i2k\pi$，而其余系数不变，故所求 Taylor 展开式为

$$
\ln(1+z)=i2k\pi+\sum_{n=1}^{\infty}\frac{(-)^{n-1}}{n}z^n. \tag{25}
$$

**例 5**　求 $f(z)=(1+z)^\alpha$（$\alpha\in\mathbb C$）在 $a=0$ 处的 Taylor 展开式，规定 $f(0)=1$。

**解**　$f(z)=(1+z)^\alpha=e^{\alpha\ln(1+z)}$ 的支点为 $-1$ 和 $\infty$，类似于上题，将 $z$ 平面沿负实轴从 $-1$ 到 $\infty$ 割破，可以分出单值分支。规定 $f(0)=1$ 即规定 $\arg(1+z)|_{z=0}=\arg[z-(-1)]|_{z=0}=0$，就确定了单值分支。按规定，$f(0)=1$，故 $c_0=1$。又

$$
f'(z)=\frac{\alpha e^{\alpha\ln(1+z)}}{1+z}
=\alpha(1+z)^{\alpha-1}
=\frac\alpha{1+z}f(z),
$$

继续求导得

$$
f^{(n)}(z)=\frac{\alpha(\alpha-1)\cdots(\alpha-n+1)}{(1+z)^n}f(z),
$$

故

$$
c_n=\frac{f^{(n)}(0)}{n!}
=\frac{\alpha(\alpha-1)\cdots(\alpha-n+1)}{n!}
=C_\alpha^n,\qquad n\in\mathbb N^+,
$$

故所求 Taylor 展开式为

$$
(1+z)^\alpha=\sum_{n=0}^{\infty}C_\alpha^n z^n, \tag{26}
$$

其中规定 $C_\alpha^0=1$。由于 $f(z)$ 距离展开中心 $a=0$ 最近的奇点为 $z_0=-1$，故展开式的收敛半径应为 $R=1$。由 Taylor 系数，

$$
\left|\frac{c_{n+1}}{c_n}\right|
=\left|\frac{\alpha-n}{n+1}\right|\xrightarrow[n\to\infty]{}1,
$$

由此亦得收敛半径 $R=1$，如所期望。

如果规定 $\arg(1+z)|_{z=0}=\arg[z-(-1)]|_{z=0}=2k\pi$（$k\in\mathbb Z$），则 $f(0)=e^{i2k\pi\alpha}$，$c_n=e^{i2k\pi\alpha}C_\alpha^n$，故所求 Taylor 展开式为

$$
(1+z)^\alpha=e^{i2k\pi\alpha}\sum_{n=0}^{\infty}C_\alpha^n z^n. \tag{27}
$$

如果 $\alpha=K\in\mathbb N^+$，则由以上系数表达易知当 $n>K$ 时，$c_n=0$，上述结果即退化为二项式定理。

**习题**　在原点的邻域上将下列函数展开为 Taylor 级数并指出收敛半径。（1）$e^{1/(1-z)}$（计算五个非零项），（2）$\sin^2z$ 和 $\cos^2z$。

## §5　解析开拓的基本概念

**定义**　设函数 $f(z)$ 在区域 $D$ 内解析，考虑一个包含 $D$ 的更大的区域 $G$，如果存在函数 $F(z)$ 在 $G$ 内解析，且在 $D$ 内 $F(z)=f(z)$，则称函数 $f(z)$ 可以解析开拓到 $G$ 内，并称 $F(z)$ 为 $f(z)$ 在 $G$ 内的解析开拓。

如果这样定义的解析开拓存在的话，可以证明它一定是唯一的。

**例**　$f(z)=\sum_{n=0}^{\infty}z^n$ 在区域 $D:|z|<1$ 内解析，而 $F(z)=1/(1-z)$ 在区域 $G=\mathbb C\backslash\{1\}$（表示复平面去掉点 $z=1$）内解析。在 $D$ 内 $F(z)=f(z)$，故 $F(z)$ 是 $f(z)$ 在 $G$ 内的解析开拓。

## 补充习题

1. 已知级数 $\sum_{n=0}^{\infty}a_nz^n$ 和 $\sum_{n=0}^{\infty}b_nz^n$ 的收敛半径分别为 $R_1$ 和 $R_2$，研究下列幂级数的收敛半径。

   (a) $\displaystyle\sum_{n=0}^{\infty}(a_n+b_n)z^n$；

   (b) $\displaystyle\sum_{n=0}^{\infty}(a_n-b_n)z^n$；

   (c) $\displaystyle\sum_{n=0}^{\infty}a_nb_nz^n$；

   (d) $\displaystyle\sum_{n=0}^{\infty}(a_n/b_n)z^n$（其中 $b_n\ne0$）。

2. 将下列函数在指定的展开中心 $a$ 处展开为 Taylor 级数并指出其收敛半径。

   (a) $\cos z$，$a=m\pi$（$m\in\mathbb Z$）；

   答案：

   $$
   \sum_{k=0}^{\infty}\frac{(-)^{k+m}}{(2k)!}(z-m\pi)^{2k},\qquad R=\infty.
   $$

   (b) $\cos z$，$a=1$；

   答案：

   $$
   \sum_{k=0}^{\infty}\frac{(-)^k\cos1}{(2k)!}(z-1)^{2k}
   -\sum_{k=0}^{\infty}\frac{(-)^k\sin1}{(2k+1)!}(z-1)^{2k+1},\qquad R=\infty.
   $$

   (c) $\sin z/(1-z)$，$a=0$；

   答案：

   $$
   \sum_{n=1}^{\infty}c_nz^n,\qquad
   c_n=\sum_{k=0}^{[(n-1)/2]}\frac{(-)^k}{(2k+1)!},\qquad R=1.
   $$

   (d) $e^z\cos z$，$a=0$；

   答案：

   $$
   \sum_{k=0}^{\infty}\frac{2^{k/2}\cos(k\pi/4)}{k!}z^k,\qquad R=\infty.
   $$

   (e) $\displaystyle\int_0^z(\sin t/t)\,dt$，$a=0$；

   答案：

   $$
   \sum_{k=0}^{\infty}\frac{(-)^k}{(2k)!(2k+1)}z^{2k+1},\qquad R=\infty.
   $$

   (f) $\ln z$（规定 $0\le\arg z<2\pi$），$a=i$；

   答案：

   $$
   i\frac\pi2-\sum_{k=1}^{\infty}\frac{i^k}{k}(z-i)^k,\qquad R=1.
   $$

   (g) $z^\alpha$（其中 $\alpha\in\mathbb C$，规定 $-\pi<\arg z\le\pi$），$a=1$；

   答案：

   $$
   \sum_{k=0}^{\infty}C_\alpha^k(z-1)^k,\qquad R=1.
   $$

   (h) $\ln(1+e^z)$，$a=0$。计算前五个非零项。
