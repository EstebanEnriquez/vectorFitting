# The latest revised Fast-Relaxed Vector Fitting for MATLAB

<p align="justify">
One of the most utilized rational approximation techniques is Vector Fitting (VF). In 1999, VF was introduced and proved to be a highly robust and efficient rational approximation method, applicable to both smooth and resonant responses. This caused VF to be rapidly adopted in several areas including power systems and macro-modeling systems. VF approximates a frequency response (generally an array) using a rational function approximation, as shown in <i>(1)</i>:
</p>

![Ecuación](https://quicklatex.com/cache3/a9/ql_cb5dcb83a99458ff7a8a1e0546ca49a9_l3.png)

<p align="justify">
where <i>f(s)</i> may represent the transfer function or the frequency response of a system, <i>c</i><sub><i>n</i></sub> and <i>a</i><sub><i>n</i></sub> represent the sets of residues and poles, respectively; and both <i>d</i> and <i>e</i> are real coefficients. Since the solution of <i>(1)</i> is a nonlinear problem, VF resolves this issue by performing the approximation in two stages, represented by two separated linear problems, as described next.
</p>

<p align="justify">
In the first stage, a set of predetermined poles <i>â</i><sub><i>n</i></sub> distributed over the frequency range of interest is used. If <i>f(s)</i>  is a smooth function, real poles are employed; otherwise, complex poles are assumed, as in <i>(2)</i> :
</p>

![Ecuación](https://quicklatex.com/cache3/10/ql_6fc08efe8be47e8c19e7d83439aeb910_l3.png)

![Ecuación](https://quicklatex.com/cache3/4a/ql_b44cef18954ea7e15ca7ad80dec8024a_l3.png)

<p align="justify">
Additionally, a frequency-dependent scaling function <i>σ(s)</i> is introduced with the same initial set of poles, as in <i>(4)</i>:
</p>

![Ecuación](https://quicklatex.com/cache3/0b/ql_d085311f8d68a9211a4dc0e638649c0b_l3.png)

<p align="justify">
Considering the multiplication between <i>f(s)</i> and <i>σ(s)</i>, the following applies:
</p>

![Ecuación](https://quicklatex.com/cache3/ef/ql_7ffc2176fc01c1e80ba614d23735f2ef_l3.png)

Substitution of <i>(4)</i> into <i>(5)</i> results in:

![Ecuación](https://quicklatex.com/cache3/b1/ql_e51e026ab43f63d2cc4a8647cc1879b1_l3.png)

<p align="justify">
Since the initial set of poles is known, <i>(6)</i> represents a linear problem which can be solved via linear least squares approximation and provides the set of residues <i>ĉ</i><sub><i>n</i></sub>, which are used to obtain the zeros of <i>σ(s)</i> through the eigenvalues computation of <i>(7)</i>:
</p>

![Ecuación](https://quicklatex.com/cache3/51/ql_5672bab639cbd4dd822f8d4603e5df51_l3.png)

<p align="justify">
where <i>A</i> is a diagonal matrix containing the set of starting poles <i>â</i><sub><i>n</i></sub>, <i>b</i> is a column vector of ones, and <i>ĉ</i><sup><i>T</i></sup> is a row vector containing the set of residues <i>ĉ</i><sub><i>n</i></sub>. Subsequently, <i>(6)</i> can be represented as <i>(8)</i> which leads to <i>(9)</i>. Hence, the zeros of <i>σ(s)</i> are considered as the new poles of <i>f(s)</i>.
</p>

<p align="justify">
In the second stage, VF uses the new set of calculated poles and <i>(1)</i> is solved in the least squares sense to calculate an updated set of residues <i>ĉ</i><sub><i>n</i></sub>. The aforementioned process can be performed iteratively until the error between the measured function and fitted function is minimum. The diagram summarizing this procedure is shown below:
</p>


