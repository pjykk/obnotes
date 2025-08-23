---
{"dg-publish":true,"permalink":"/oi/maths/poly/fft-and-ntt/","tags":["gardenEntry"]}
---


## 0.  符号约定

我们用大写字母表示多项式，小写字母表示多项式的系数，比如：

$$A(x)=a_0+a_1x+a_2x^2+\cdots+a_nx^n$$
$\omega_n^k$ 表示单位根，即
$$\omega_n^k=\cos\dfrac{2k\pi}{n}+i\sin\dfrac{2k\pi}{n}$$
## 1. FFT

### 0. 基本思路

考虑对两个多项式进行乘法：

$$C(x)=A(x)B(x)$$
显然系数满足
$$c_n=\sum\limits_{k=0}^na_kb_{n-k}$$
```ad-note
title: 卷积
形如 
$c_n=\sum\limits_{j\oplus k=n}a_jb_k$
的式子叫作卷积。

对于多项式乘法，所计算的就是加法卷积。
```

朴素计算显然是 $\Theta(n^2)$ 的，很不优秀。这种不优秀实质上来源于系数表示，我们需要用其他的方法表示多项式。

注意到对于 $(n-1)$ 次多项式，我们可以用 $n$ 个不同的点处的值来进行唯一表示。因此我们可以对多项式 $A$ 求出其点值 $A(x_0),A(x_1),\cdots,A(x_{n-1})$，对于多项式 $B$ 同理，求出 $B(x_0),B(x_1),\cdots,B(x_{n-1})$ 后，乘积 $C(x)$ 对应的点值显然有

$$C(x_k)=A(x_k)B(x_k)$$

然后再进行插值即可得到 $C$ 的系数表示。

```ad-attention 
title: 注意
上面的 $n$ 为 $C$ 的次数 $+1$。
```

如此我们将 $\Theta(n^2)$ 的朴素卷积转化为了 $\Theta(n)$ 的点值的乘法，关键在于求点值和插值两步的优化。

### 1. 求点值（DFT）

如果随意代入 $x_k$，复杂度仍为 $\Theta(n^2)$，没有优化。但是 $x_i$ 是可以任意选取的，我们可以选取更容易计算的点进行计算。这里的技巧是代入单位根 $x_k=\omega_n^k$。

设 $A(x)=a_0+a_1x+\cdots+a_{n-1}x^{n-1}$（次数不够时高次项系数视为 $0$ 即可），有

$$A(\omega_n^k)=\sum\limits_{j=0}^{n-1}a_j\omega_n^{jk}$$
利用单位根 $\omega_{nd}^{kd}=\omega_n^k$ 的性质，我们能够轻易处理数列的拉伸操作，进而通过分治解决问题。

```ad-note
title: 什么是拉伸操作？
collapse:
考虑数列 $a_0,a_1,a_2,\cdots$，拉伸变换将数列变为 $a_0,\underbrace{0,\cdots,0}_{(k-1)\text{个}},a_1,\underbrace{0,\cdots,0}_{(k-1)\text{个}},a_2,\cdots$
在生成函数（以数列的每一项作为系数的幂级数）的语境下，
设 $A(x)=a_0+a_1x+a_2x^2+\cdots$，拉伸后的序列对应的幂级数即为 $A(x^k)$。
```

这里以二分治为例。方便起见将 $n$ 补为 $2$ 的幂。

奇偶分组便于使用拉伸操作：

设 $A_0(x)=a_0+a_2x+\cdots+a_{n-2}x^{n/2-1}$，$A_1(x)=a_1+a_3x+\cdots+a_{n-1}x^{n/2-1}$，则有

$$A(x)=A_0(x^2)+xA_1(x^2)$$

代入单位根立得

$$A(\omega_n^k)=A_0(\omega_{n/2}^k)+\omega_{n}^kA_1(\omega_{n/2}^k)$$
容易发现 $A_0(\omega_{n/2}^k)$ 与 $A_1(\omega_{n/2}^k)$ 即为两个子问题，因此可以分治求解……吗？

这里有一个问题，子问题中 $k$ 的取值仅为 $0,1,\cdots,n/2-1$，还有另一半需要计算。

因此我们还需要代入 $\omega_n^{k+n/2}=-\omega_n^k$，进行计算得

$$A(\omega_n^{k+n/2})=A_0(\omega_{n/2}^k)-\omega_{n}^kA_1(\omega_{n/2}^k)$$
至此我们已经完成了快速计算点值的操作。由 Master 定理易知时间复杂度 $\Theta(n\log n)$。

实质上 DFT 是一个序列变换的过程，它将 $a_0,\cdots,a_{n-1}$ 变换为了 $A(x_0),\cdots,A(x_{n-1})$。
并且这是一个线性的过程：

$$\begin{pmatrix}1&1&1&\cdots&1\\1&\omega_n^1&\omega_n^2&\cdots&\omega_n^{n-1}\\ 1&\omega_n^2&\omega_n^4&\cdots&\omega_n^{2(n-1)}\\ \vdots&\vdots&\vdots&\ddots&\vdots\\1&\omega_n^{n-1}&\omega_n^{2(n-1)}&\cdots&\omega_n^{(n-1)(n-1)}\\ \end{pmatrix}\begin{pmatrix}a_0\\a_1\\a_2\\ \vdots\\a_{n-1}\end{pmatrix}=\begin{pmatrix}A(x_0)\\A(x_1)\\A(x_2)\\ \vdots\\A(x_{n-1})\end{pmatrix}$$

设 DFT 将序列 $X$ 变换为序列 $Y$，记作 $\mathsf{DFT}(X)=Y$，根据线性性得到：
- $\mathsf{DFT}(X_1+X_2)=\mathsf{DFT}(X_1)+\mathsf{DFT}(X_2)$
- $\mathsf{DFT}(kX)=k\mathsf{DFT}(X)$
### 2. 插值（IDFT）

利用上面的观点，容易发现 IDFT 只需把上面矩阵的逆矩阵求出即可。

但是我们有更直观的方法，就是直接使用拉格朗日插值：

$$A(x)=\sum\limits_{k=0}^{n-1}A(x_k)\prod\limits_{j\neq k}\dfrac{x-x_j}{x_k-x_j}$$

注意，利用单位根的性质，后面那个乘积项实际上是可以直接算出的，分子：

$$\prod\limits_{j\neq k}(x-x_j)=\dfrac{x^n-1}{x-\omega_n^k}$$

分母：

$$\prod\limits_{j\neq k}(x_k-x_j)=\lim_{x\rightarrow x_k}\dfrac{x^n-1}{x-\omega_n^k}=n\omega_n^{k(n-1)}=n\omega_n^{-k}$$

于是

$$\begin{aligned}\prod\limits_{j\neq k}\dfrac{x-x_j}{x_k-x_j}&=\dfrac{\omega_n^k}{n}\dfrac{x^n-1}{x-\omega_n^k}\\&=\dfrac{1}{n}\left(\dfrac{1}{1-\omega_n^{-k}x}-\dfrac{x^n}{1-\omega_n^{-k}x}\right)\\&=\dfrac{1}{n}\left(\sum\limits_{j=0}^{+\infty}\omega_n^{-jk}x^j-\sum\limits_{j=0}^{+\infty}\omega_n^{-jk}x^{j+n}\right)\\&=\dfrac{1}{n}\left(\sum\limits_{j=0}^{+\infty}\omega_n^{-jk}x^j-\sum\limits_{j=0}^{+\infty}\omega_n^{-(j+n)k}x^{j+n}\right)\\&=\dfrac{1}{n}\sum\limits_{j=0}^{n-1}\omega_n^{-jk}x^j\end{aligned}$$
代入整理得

$$\begin{aligned}A(x)&=\sum\limits_{k=0}^{n-1}A(x_k)\dfrac{1}{n}\sum\limits_{j=0}^{n-1}\omega_n^{-jk}x^j\\&=\sum\limits_{j=0}^{n-1}\left(\dfrac{1}{n}\sum\limits_{k=0}^{n-1}A(x_k)\omega_n^{-jk}\right)x^j\end{aligned}$$
于是（这里更换了下标以保证形式与 DFT 的统一）
$$a_k=\dfrac{1}{n}\sum\limits_{j=0}^{n-1}A(x_j)\omega_n^{-jk}$$
对比 DFT 的式子
$$A(\omega_n^k)=\sum\limits_{j=0}^{n-1}a_j\omega_n^{jk}$$
发现 IDFT 只不过是将单位根换成了 $\omega_{n}^{-j}$ 并多了一个 $\dfrac{1}{n}$ 的系数。
在性质上我们实际上无法区分 $\omega_n^j$ 和 $\omega_n^{-j}$，算法与先前 DFT 完全一致。

### 3. 代码实现

To be finished