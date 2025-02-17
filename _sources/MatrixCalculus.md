# Matrix Differentiation

( and some other stuff)

Randal J. Barnes,
Department of Civil Engineering, University of Minnesota.
Minneapolis, Minnesota, USA

## 1 Introduction

Throughout this presentation I have chosen to use a symbolic matrix notation. This choice was not made lightly. I am a strong advocate of index notation, when appropriate. For example, index notation greatly simplifies the presentation and manipulation of differential geometry. As a rule-of-thumb, if your work is going to primarily involve differentiation with respect to the spatial coordinates, then index notation is almost surely the appropriate choice.

In the present case, however, I will be manipulating large systems of equations in which the matrix calculus is relatively simply while the matrix algebra and matrix arithmetic is messy and more involved. Thus, I have chosen to use symbolic notation.

## 2 Notation and Nomenclature

Definition 1 Let $\mathfrak{a}_{\mathfrak{i j}} \in \mathfrak{R}, \mathfrak{i}=1,2, \ldots, \mathfrak{m}, \mathfrak{j}=1,2, \ldots, \mathfrak{n}$. Then the ordered rectangular array

$$
\mathbf{A}=\left[\begin{array}{cccc}
a_{11} & a_{12} & \cdots & a_{1 n}  \tag{1}\\
a_{21} & a_{22} & \cdots & a_{2 n} \\
\vdots & \vdots & & \vdots \\
a_{\mathfrak{m} 1} & a_{\mathfrak{m} 2} & \cdots & a_{\mathfrak{m} n}
\end{array}\right]
$$

is said to be a real matrix of dimension $\mathfrak{m} \times \mathfrak{n}$.
When writing a matrix I will occasionally write down its typical element as well as its dimension. Thus,

$$
\begin{equation*}
\mathbf{A}=\left[\mathfrak{a}_{\mathfrak{i j}}\right], \quad \mathfrak{i}=1,2, \ldots, \mathfrak{m} ; \mathfrak{j}=1,2, \ldots, \mathbf{n}, \tag{2}
\end{equation*}
$$

denotes a matrix with $m$ rows and $n$ columns, whose typical element is $a_{i j}$. Note, the first subscript locates the row in which the typical element lies while the second subscript locates the column. For example, $\mathfrak{a}_{\mathfrak{j k}}$ denotes the element lying in the $\mathfrak{j}$ th row and kth column of the matrix $\mathbf{A}$.

Definition 2 A vector is a matrix with only one column. Thus, all vectors are inherently column vectors.

## Convention 1

Multi-column matrices are denoted by boldface uppercase letters: for example, $\mathbf{A}, \mathbf{B}, \mathbf{X}$. Vectors (single-column matrices) are denoted by boldfaced lowercase letters: for example, $\mathbf{a}, \mathbf{b}, \mathbf{x}$. I will attempt to use letters from the beginning of the alphabet to designate known matrices, and letters from the end of the alphabet for unknown or variable matrices.

## Convention 2

When it is useful to explicitly attach the matrix dimensions to the symbolic notation, I will use an underscript. For example, $\underset{m \times n}{\mathbf{A}}$, indicates a known, multi-column matrix with $m$ rows and $n$ columns.

A superscript ${ }^{\top}$ denotes the matrix transpose operation; for example, $\mathbf{A}^{\top}$ denotes the transpose of $\mathbf{A}$. Similarly, if $\mathbf{A}$ has an inverse it will be denoted by $\mathbf{A}^{-1}$. The determinant of $\mathbf{A}$ will be denoted by either $|\mathbf{A}|$ or $\operatorname{det}(\mathbf{A})$. Similarly, the $\operatorname{rank}$ of a matrix $\mathbf{A}$ is denoted by $\operatorname{rank}(\mathbf{A})$. An identity matrix will be denoted by $\mathbf{I}$, and $\mathbf{0}$ will denote a null matrix.

## 3 Matrix Multiplication

Definition 3 Let $\mathbf{A}$ be $\mathfrak{m} \times \mathfrak{n}$, and $\mathbf{B}$ be $\mathfrak{n} \times p$, and let the product $\mathbf{A B}$ be

$$
\begin{equation*}
\mathbf{C}=\mathbf{A B} \tag{3}
\end{equation*}
$$

then $\mathbf{C}$ is a $\mathfrak{m} \times p$ matrix, with element $(\mathfrak{i}, \mathfrak{j})$ given by

$$
\begin{equation*}
c_{i j}=\sum_{k=1}^{n} a_{i k} b_{k j} \tag{4}
\end{equation*}
$$

for all $i=1,2, \ldots, m, \quad j=1,2, \ldots, p$.
Proposition 1 Let A be $\mathrm{m} \times \mathrm{n}$, and $\mathbf{x}$ be $\mathrm{n} \times 1$, then the typical element of the product

$$
\begin{equation*}
\mathbf{z}=\mathbf{A x} \tag{5}
\end{equation*}
$$

is given by

$$
\begin{equation*}
z_{i}=\sum_{k=1}^{n} a_{i k} x_{k} \tag{6}
\end{equation*}
$$

for all $\mathfrak{i}=1,2, \ldots, \boldsymbol{m}$. Similarly, let $\mathbf{y}$ be $\mathbf{m} \times 1$, then the typical element of the product

$$
\begin{equation*}
\mathbf{z}^{\top}=\mathbf{y}^{\top} \mathbf{A} \tag{7}
\end{equation*}
$$

is given by

$$
\begin{equation*}
z_{i}=\sum_{k=1}^{n} a_{k i} y_{k} \tag{8}
\end{equation*}
$$

for all $\mathbf{i}=1,2, \ldots, n$. Finally, the scalar resulting from the product

$$
\begin{equation*}
\alpha=\mathbf{y}^{\top} \mathbf{A} \mathbf{x} \tag{9}
\end{equation*}
$$

is given by

$$
\begin{equation*}
\alpha=\sum_{j=1}^{m} \sum_{k=1}^{n} a_{j k} y_{j} x_{k} \tag{10}
\end{equation*}
$$

Proof: These are merely direct applications of Definition 3. q.e.d.

Proposition 2 Let $\mathbf{A}$ be $\mathfrak{m} \times \mathfrak{n}$, and $\mathbf{B}$ be $\mathfrak{n} \times \mathrm{p}$, and let the product $\mathbf{A B}$ be

$$
\begin{equation*}
\mathbf{C}=\mathbf{A B} \tag{11}
\end{equation*}
$$

then

$$
\begin{equation*}
\mathbf{C}^{\top}=\mathbf{B}^{\top} \mathbf{A}^{\top} \tag{12}
\end{equation*}
$$

Proof: The typical element of $\mathbf{C}$ is given by

$$
\begin{equation*}
c_{i j}=\sum_{k=1}^{n} a_{i k} b_{k j} \tag{13}
\end{equation*}
$$

By definition, the typical element of $\mathbf{C}^{\top}$, say $\mathrm{d}_{\mathfrak{i j}}$, is given by

$$
\begin{equation*}
\mathrm{d}_{\mathfrak{i j}}=\mathrm{c}_{\mathfrak{j i}}=\sum_{\mathrm{k}=1}^{n} \mathrm{a}_{\mathfrak{j k}} \mathrm{b}_{\mathrm{ki}} \tag{14}
\end{equation*}
$$

Hence,

$$
\begin{equation*}
\mathbf{C}^{\top}=\mathbf{B}^{\top} \mathbf{A}^{\top} \tag{15}
\end{equation*}
$$

q.e.d.

Proposition 3 Let $\mathbf{A}$ and $\mathbf{B}$ be $\mathrm{n} \times \mathrm{n}$ and invertible matrices. Let the product $\mathbf{A B}$ be given by

$$
\begin{equation*}
\mathbf{C}=\mathbf{A B} \tag{16}
\end{equation*}
$$

then

$$
\begin{equation*}
\mathbf{C}^{-1}=\mathbf{B}^{-1} \mathbf{A}^{-1} \tag{17}
\end{equation*}
$$

Proof:

$$
\begin{equation*}
\mathbf{C B}^{-1} \mathbf{A}^{-1}=\mathbf{A B B}^{-1} \mathbf{A}^{-1}=\mathbf{I} \tag{18}
\end{equation*}
$$

q.e.d.

## 4 Partioned Matrices

Frequently, I will find it convenient to deal with partitioned matrices ${ }^{1}$. Such a representation, and the manipulation of this representation, are two of the relative advantages of the symbolic matrix notation.

Definition 4 Let $\mathbf{A}$ be $m \times n$ and write

$$
A=\left[\begin{array}{ll}
B & C  \tag{19}\\
D & E
\end{array}\right]
$$

where $\mathbf{B}$ is $\mathfrak{m}_{1} \times \mathfrak{n}_{1}, \mathbf{E}$ is $\mathfrak{m}_{2} \times \mathfrak{n}_{2}, \mathbf{C}$ is $\mathfrak{m}_{1} \times \mathfrak{n}_{2}, \mathbf{D}$ is $\mathfrak{m}_{2} \times \mathfrak{n}_{1}, \mathfrak{m}_{1}+\mathfrak{m}_{2}=\mathfrak{m}$, and $\mathfrak{n}_{1}+\mathfrak{n}_{2}=\mathfrak{n}$. The above is said to be a partition of the matrix $\mathbf{A}$.

[^0]Proposition 4 Let A be a square, nonsingular matrix of order m. Partition A as

$$
\mathbf{A}=\left[\begin{array}{ll}
\mathbf{A}_{11} & \mathbf{A}_{12}  \tag{20}\\
\mathbf{A}_{21} & \mathbf{A}_{22}
\end{array}\right]
$$

so that $\mathbf{A}_{11}$ is a nonsingular matrix of order $\mathrm{m}_{1}, \mathbf{A}_{22}$ is a nonsingular matrix of order $\mathrm{m}_{2}$, and $\mathrm{m}_{1}+\mathrm{m}_{2}=\mathrm{m}$. Then

$$
\mathbf{A}^{-1}=\left[\begin{array}{cc}
\left(\mathbf{A}_{11}-\mathbf{A}_{12} \mathbf{A}_{22}^{-1} \mathbf{A}_{21}\right)^{-1} & -\mathbf{A}_{11}^{-1} \mathbf{A}_{12}\left(\mathbf{A}_{22}-\mathbf{A}_{21} \mathbf{A}_{11}^{-1} \mathbf{A}_{12}\right)^{-1}  \tag{21}\\
-\mathbf{A}_{22}^{-1} \mathbf{A}_{21}\left(\mathbf{A}_{11}-\mathbf{A}_{12} \mathbf{A}_{22}^{-1} \mathbf{A}_{21}\right)^{-1} & \left(\mathbf{A}_{22}-\mathbf{A}_{21} \mathbf{A}_{11}^{-1} \mathbf{A}_{12}\right)^{-1}
\end{array}\right]
$$

Proof: Direct multiplication of the proposed $\mathbf{A}^{-1}$ and $\mathbf{A}$ yields

$$
\begin{equation*}
\mathbf{A}^{-1} \mathbf{A}=\mathbf{I} \tag{22}
\end{equation*}
$$

q.e.d.

## 5 Matrix Differentiation

In the following discussion I will differentiate matrix quantities with respect to the elements of the referenced matrices. Although no new concept is required to carry out such operations, the element-by-element calculations involve cumbersome manipulations and, thus, it is useful to derive the necessary results and have them readily available ${ }^{2}$.

## Convention 3

Let

$$
\begin{equation*}
\mathbf{y}=\psi(\mathbf{x}) \tag{23}
\end{equation*}
$$

where $\mathbf{y}$ is an $m$-element vector, and $\mathbf{x}$ is an $n$-element vector. The symbol

$$
\frac{\partial \mathbf{y}}{\partial \mathbf{x}}=\left[\begin{array}{cccc}
\frac{\partial y_{1}}{\partial x_{1}} & \frac{\partial y_{1}}{\partial x_{2}} & \cdots & \frac{\partial y_{1}}{\partial x_{n}}  \tag{24}\\
\frac{\partial y_{2}}{\partial x_{1}} & \frac{\partial y_{2}}{\partial x_{2}} & \cdots & \frac{\partial y_{2}}{\partial x_{n}} \\
\vdots & \vdots & & \vdots \\
\frac{\partial y_{m}}{\partial x_{1}} & \frac{\partial y_{m}}{\partial x_{2}} & \cdots & \frac{\partial y_{m}}{\partial x_{n}}
\end{array}\right]
$$

will denote the $m \times n$ matrix of first-order partial derivatives of the transformation from $\mathbf{x}$ to $\mathbf{y}$. Such a matrix is called the Jacobian matrix of the transformation $\psi()$.

Notice that if $\mathbf{x}$ is actually a scalar in Convention 3 then the resulting Jacobian matrix is a $m \times 1$ matrix; that is, a single column (a vector). On the other hand, if $\mathbf{y}$ is actually a scalar in Convention 3 then the resulting Jacobian matrix is a $1 \times \mathfrak{n}$ matrix; that is, a single row (the transpose of a vector).

Proposition 5 Let

$$
\begin{equation*}
y=A x \tag{25}
\end{equation*}
$$

[^1]where $\mathbf{y}$ is $\mathrm{m} \times 1$, $\mathbf{x}$ is $\mathfrak{n} \times 1, \mathbf{A}$ is $\mathfrak{m} \times \mathfrak{n}$, and $\mathbf{A}$ does not depend on $\mathbf{x}$, then
\$\$

$$
\begin{equation*}
\frac{\partial \mathbf{y}}{\partial \mathbf{x}}=\mathbf{A} \tag{26}
\end{equation*}
$$

\$\$

Proof: Since the $\mathbf{i}$ th element of $\mathbf{y}$ is given by

$$
\begin{equation*}
y_{i}=\sum_{k=1}^{n} a_{i k} x_{k} \tag{27}
\end{equation*}
$$

it follows that

$$
\begin{equation*}
\frac{\partial y_{i}}{\partial x_{j}}=a_{i j} \tag{28}
\end{equation*}
$$

for all $\mathfrak{i}=1,2, \ldots, \mathfrak{m}, \quad \mathfrak{j}=1,2, \ldots, \boldsymbol{n}$. Hence

$$
\begin{equation*}
\frac{\partial \mathbf{y}}{\partial \mathbf{x}}=\mathbf{A} \tag{29}
\end{equation*}
$$

q.e.d.

Proposition 6 Let

$$
\begin{equation*}
\mathbf{y}=\mathbf{A x} \tag{30}
\end{equation*}
$$

where $\mathbf{y}$ is $\mathrm{m} \times 1$, $\mathbf{x}$ is $\mathfrak{n} \times 1, \mathbf{A}$ is $\mathrm{m} \times \mathfrak{n}$, and $\mathbf{A}$ does not depend on $\mathbf{x}$, as in Proposition 5 . Suppose that $\mathbf{x}$ is a function of the vector $\mathbf{z}$, while $\mathbf{A}$ is independent of $\mathbf{z}$. Then

$$
\begin{equation*}
\frac{\partial \mathbf{y}}{\partial \mathbf{z}}=\mathbf{A} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{31}
\end{equation*}
$$

Proof: Since the $\mathbf{i}$ th element of $\mathbf{y}$ is given by

$$
\begin{equation*}
y_{i}=\sum_{k=1}^{n} a_{i k} x_{k} \tag{32}
\end{equation*}
$$

for all $\mathfrak{i}=1,2, \ldots, m$, it follows that

$$
\begin{equation*}
\frac{\partial y_{i}}{\partial z_{j}}=\sum_{k=1}^{n} a_{i k} \frac{\partial x_{k}}{\partial z_{j}} \tag{33}
\end{equation*}
$$

but the right hand side of the above is simply element $(\mathfrak{i}, \mathfrak{j})$ of $\mathbf{A} \frac{\partial \mathbf{x}}{\partial \mathrm{z}}$. Hence

$$
\begin{equation*}
\frac{\partial \mathbf{y}}{\partial z}=\frac{\partial \mathbf{y}}{\partial \mathbf{x}} \frac{\partial \mathbf{x}}{\partial \mathbf{z}}=\mathbf{A} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{34}
\end{equation*}
$$

q.e.d.

Proposition 7 Let the scalar $\alpha$ be defined by

$$
\begin{equation*}
\alpha=\mathbf{y}^{\top} \mathbf{A} \mathbf{x} \tag{35}
\end{equation*}
$$

where $\mathbf{y}$ is $\mathrm{m} \times 1, \mathbf{x}$ is $\mathrm{n} \times 1, \mathbf{A}$ is $\mathrm{m} \times \mathrm{n}$, and $\mathbf{A}$ is independent of $\mathbf{x}$ and $\mathbf{y}$, then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{x}}=\mathbf{y}^{\top} \mathbf{A} \tag{36}
\end{equation*}
$$

and

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{y}}=\mathbf{x}^{\top} \mathbf{A}^{\top} \tag{37}
\end{equation*}
$$

Proof: Define

$$
\begin{equation*}
\mathbf{w}^{\top}=\mathbf{y}^{\top} \mathbf{A} \tag{38}
\end{equation*}
$$

and note that

$$
\begin{equation*}
\alpha=\mathbf{w}^{\top} \mathbf{x} \tag{39}
\end{equation*}
$$

Hence, by Proposition 5 we have that

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{x}}=\mathbf{w}^{\top}=\mathbf{y}^{\top} \mathbf{A} \tag{40}
\end{equation*}
$$

which is the first result. Since $\alpha$ is a scalar, we can write

$$
\begin{equation*}
\alpha=\alpha^{\top}=\mathbf{x}^{\top} \mathbf{A}^{\top} \mathbf{y} \tag{41}
\end{equation*}
$$

and applying Proposition 5 as before we obtain

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{y}}=\mathbf{x}^{\top} \mathbf{A}^{\top} \tag{42}
\end{equation*}
$$

q.e.d.

Proposition 8 For the special case in which the scalar $\alpha$ is given by the quadratic form

$$
\begin{equation*}
\alpha=\mathbf{x}^{\top} \mathbf{A} \mathbf{x} \tag{43}
\end{equation*}
$$

where $\mathbf{x}$ is $\mathfrak{n} \times 1$, $\mathbf{A}$ is $\mathfrak{n} \times \mathfrak{n}$, and $\mathbf{A}$ does not depend on $\mathbf{x}$, then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{x}}=\mathbf{x}^{\top}\left(\mathbf{A}+\mathbf{A}^{\top}\right) \tag{44}
\end{equation*}
$$

Proof: By definition

$$
\begin{equation*}
\alpha=\sum_{j=1}^{n} \sum_{i=1}^{n} a_{i j} x_{i} x_{j} \tag{45}
\end{equation*}
$$

Differentiating with respect to the k th element of $\mathbf{x}$ we have

$$
\begin{equation*}
\frac{\partial \alpha}{\partial x_{k}}=\sum_{j=1}^{n} a_{k j} x_{j}+\sum_{i=1}^{n} a_{i k} x_{i} \tag{46}
\end{equation*}
$$

for all $\mathrm{k}=1,2, \ldots, \mathrm{n}$, and consequently,

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{x}}=\mathbf{x}^{\top} \mathbf{A}^{\top}+\mathbf{x}^{\top} \mathbf{A}=\mathbf{x}^{\top}\left(\mathbf{A}^{\top}+\mathbf{A}\right) \tag{47}
\end{equation*}
$$

q.e.d.

Proposition 9 For the special case where $\mathbf{A}$ is a symmetric matrix and

$$
\begin{equation*}
\alpha=\mathbf{x}^{\top} \mathbf{A} \mathbf{x} \tag{48}
\end{equation*}
$$

where $\mathbf{x}$ is $\mathrm{n} \times 1, \mathbf{A}$ is $\mathrm{n} \times \mathrm{n}$, and $\mathbf{A}$ does not depend on $\mathbf{x}$, then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{x}}=2 \mathbf{x}^{\top} \mathbf{A} \tag{49}
\end{equation*}
$$

Proof: This is an obvious application of Proposition 8. q.e.d.
Proposition 10 Let the scalar $\alpha$ be defined by

$$
\begin{equation*}
\alpha=\mathbf{y}^{\top} \mathbf{x} \tag{50}
\end{equation*}
$$

where $\mathbf{y}$ is $\mathfrak{n} \times 1$, $\mathbf{x}$ is $\mathfrak{n} \times 1$, and both $\mathbf{y}$ and $\mathbf{x}$ are functions of the vector $\mathbf{z}$. Then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=\mathbf{x}^{\top} \frac{\partial \mathbf{y}}{\partial \mathbf{z}}+\mathbf{y}^{\top} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{51}
\end{equation*}
$$

Proof: We have

$$
\begin{equation*}
\alpha=\sum_{j=1}^{n} x_{j} y_{j} \tag{52}
\end{equation*}
$$

Differentiating with respect to the k th element of $\mathbf{z}$ we have

$$
\begin{equation*}
\frac{\partial \alpha}{\partial z_{k}}=\sum_{j=1}^{n}\left(x_{j} \frac{\partial y_{j}}{\partial z_{k}}+y_{j} \frac{\partial x_{j}}{\partial z_{k}}\right) \tag{53}
\end{equation*}
$$

for all $\mathrm{k}=1,2, \ldots, \mathrm{n}$, and consequently,

$$
\begin{equation*}
\frac{\partial \alpha}{\partial z}=\frac{\partial \alpha}{\partial \mathbf{y}} \frac{\partial \mathbf{y}}{\partial \mathbf{z}}+\frac{\partial \alpha}{\partial \mathbf{x}} \frac{\partial \mathbf{x}}{\partial \mathbf{z}}=\mathbf{x}^{\top} \frac{\partial \mathbf{y}}{\partial \mathbf{z}}+\mathbf{y}^{\top} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{54}
\end{equation*}
$$

q.e.d.

Proposition 11 Let the scalar $\alpha$ be defined by

$$
\begin{equation*}
\alpha=\mathbf{x}^{\top} \mathbf{x} \tag{55}
\end{equation*}
$$

where $\mathbf{x}$ is $\mathrm{n} \times 1$, and $\mathbf{x}$ is a function of the vector $\mathbf{z}$. Then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=2 \mathbf{x}^{\top} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{56}
\end{equation*}
$$

Proof: This is an obvious application of Proposition 10. q.e.d.
Proposition 12 Let the scalar $\alpha$ be defined by

$$
\begin{equation*}
\alpha=\mathbf{y}^{\top} \mathbf{A} \mathbf{x} \tag{57}
\end{equation*}
$$

where $\mathbf{y}$ is $\mathrm{m} \times 1, \mathbf{x}$ is $\mathrm{n} \times 1, \mathbf{A}$ is $\mathrm{m} \times \mathfrak{n}$, and both $\mathbf{y}$ and $\mathbf{x}$ are functions of the vector $\mathbf{z}$, while $\mathbf{A}$ does not depend on $\mathbf{z}$. Then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=\mathbf{x}^{\top} \mathbf{A}^{\top} \frac{\partial \mathbf{y}}{\partial \mathbf{z}}+\mathbf{y}^{\top} \mathbf{A} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{58}
\end{equation*}
$$

Proof: Define

$$
\begin{equation*}
\mathbf{w}^{\top}=\mathbf{y}^{\top} \mathbf{A} \tag{59}
\end{equation*}
$$

and note that

$$
\begin{equation*}
\alpha=\mathbf{w}^{\top} \mathbf{x} \tag{60}
\end{equation*}
$$

Applying Propositon 10 we have

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=\mathbf{x}^{\top} \frac{\partial \mathbf{w}}{\partial \mathbf{z}}+\mathbf{w}^{\top} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{61}
\end{equation*}
$$

Substituting back in for $\mathbf{w}$ we arrive at

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=\frac{\partial \alpha}{\partial \mathbf{y}} \frac{\partial \mathbf{y}}{\partial \mathbf{z}}+\frac{\partial \alpha}{\partial \mathbf{x}} \frac{\partial \mathbf{x}}{\partial \mathbf{z}}=\mathbf{x}^{\top} \mathbf{A}^{\top} \frac{\partial \mathbf{y}}{\partial \mathbf{z}}+\mathbf{y}^{\top} \mathbf{A} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{62}
\end{equation*}
$$

q.e.d.

Proposition 13 Let the scalar $\alpha$ be defined by the quadratic form

$$
\begin{equation*}
\alpha=\mathbf{x}^{\top} \mathbf{A} \mathbf{x} \tag{63}
\end{equation*}
$$

where $\mathbf{x}$ is $\mathfrak{n} \times 1$, $\mathbf{A}$ is $\mathfrak{n} \times \mathfrak{n}$, and $\mathbf{x}$ is a function of the vector $\mathbf{z}$, while $\mathbf{A}$ does not depend on $\mathbf{z}$. Then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=\mathbf{x}^{\top}\left(\mathbf{A}+\mathbf{A}^{\top}\right) \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{64}
\end{equation*}
$$

Proof: This is an obvious application of Proposition 12. q.e.d.
Proposition 14 For the special case where $\mathbf{A}$ is a symmetric matrix and

$$
\begin{equation*}
\alpha=\mathbf{x}^{\top} \mathbf{A} \mathbf{x} \tag{65}
\end{equation*}
$$

where $\mathbf{x}$ is $\mathfrak{n} \times 1$, $\mathbf{A}$ is $\mathfrak{n} \times \mathfrak{n}$, and $\mathbf{x}$ is a function of the vector $\mathbf{z}$, while $\mathbf{A}$ does not depend on $\mathbf{z}$. Then

$$
\begin{equation*}
\frac{\partial \alpha}{\partial \mathbf{z}}=2 \mathbf{x}^{\top} \mathbf{A} \frac{\partial \mathbf{x}}{\partial \mathbf{z}} \tag{66}
\end{equation*}
$$

Proof: This is an obvious application of Proposition 13. q.e.d.
Definition 5 Let $\mathbf{A}$ be a $\mathfrak{m} \times \mathfrak{n}$ matrix whose elements are functions of the scalar parameter $\alpha$. Then the derivative of the matrix $\mathbf{A}$ with respect to the scalar parameter $\alpha$ is the $m \times n$ matrix of element-by-element derivatives:

$$
\frac{\partial \mathbf{A}}{\partial \alpha}=\left[\begin{array}{cccc}
\frac{\partial a_{11}}{\partial \alpha} & \frac{\partial a_{12}}{\partial \alpha} & \cdots & \frac{\partial a_{1 n}}{\partial \alpha}  \tag{67}\\
\frac{\partial a_{21}}{\partial \alpha} & \frac{\partial a_{22}}{\partial \alpha} & \cdots & \frac{\partial a_{2 n}}{\partial \alpha} \\
\vdots & \vdots & & \vdots \\
\frac{\partial a_{m 1}}{\partial \alpha} & \frac{\partial a_{m 2}}{\partial \alpha} & \cdots & \frac{\partial a_{m n}}{\partial \alpha}
\end{array}\right]
$$

Proposition 15 Let A be a nonsingular, $\mathrm{m} \times \mathrm{m}$ matrix whose elements are functions of the scalar parameter $\alpha$. Then

$$
\begin{equation*}
\frac{\partial \mathbf{A}^{-1}}{\partial \alpha}=-\mathbf{A}^{-1} \frac{\partial \mathbf{A}}{\partial \alpha} \mathbf{A}^{-1} \tag{68}
\end{equation*}
$$

Proof: Start with the definition of the inverse

$$
\begin{equation*}
\mathbf{A}^{-1} \mathbf{A}=\mathbf{I} \tag{69}
\end{equation*}
$$

and differentiate, yielding

$$
\begin{equation*}
\mathbf{A}^{-1} \frac{\partial \mathbf{A}}{\partial \alpha}+\frac{\partial \mathbf{A}^{-1}}{\partial \alpha} \mathbf{A}=\mathbf{0} \tag{70}
\end{equation*}
$$

rearranging the terms yields

$$
\begin{equation*}
\frac{\partial \mathbf{A}^{-1}}{\partial \alpha}=-\mathbf{A}^{-1} \frac{\partial \mathbf{A}}{\partial \alpha} \mathbf{A}^{-1} \tag{71}
\end{equation*}
$$

q.e.d.

## 6 References

- Dhrymes, Phoebus J., 1978, Mathematics for Econometrics, Springer-Verlag, New York, 136 pp.
- Golub, Gene H., and Charles F. Van Loan, 1983, Matrix Computations, Johns Hopkins University Press, Baltimore, Maryland, 476 pp.
- Graybill, Franklin A., 1983, Matrices with Applications in Statistics, 2nd Edition, Wadsworth International Group, Belmont, California, 461 pp .

[^0]: ${ }^{1}$ Much of the material in this section is extracted directly from Dhrymes (1978, Section 2.7).
[^1]: ${ }^{2}$ Much of the material in this section is extracted directly from Dhrymes (1978, Section 4.3). The interested reader is directed to this worthy reference to find additional results.
