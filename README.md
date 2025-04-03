# The latest revised Fast-Relaxed Vector Fitting for MATLAB

## Vector Fitting Theory

### Classical Vector Fitting

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

![Ecuación](https://quicklatex.com/cache3/80/ql_6353db49280f60dc64bf357d2f367b80_l3.png)

![Ecuación](https://quicklatex.com/cache3/55/ql_d2cf3b48dc2179027c0ab5b2ad87cb55_l3.png)

<p align="justify">
where <i>ẑ</i><sub><i>n</i></sub> are the zeros of <i>σ</i><sub><i>fit</i></sub><i>(s)</i>, and <i>z</i><sub><i>n</i></sub> are the zeros of <i>(σf)</i><sub><i>fit</i></sub><i>(s)</i>.
</p>

<p align="justify">
In the second stage, VF uses the new set of calculated poles and <i>(1)</i> is solved in the least squares sense to obtain an updated set of residues <i>ĉ</i><sub><i>n</i></sub>. The aforementioned process can be performed iteratively until the error between the measured function and fitted function is minimum. The diagram summarizing this procedure is shown below:
</p>

### Improvement in pole relocation: Relaxed condition of scaling function

<p align="justify">
VF operates by iteratively adjusting an initial set of poles to more optimal locations. When fitting the frequency-domain response of a rational function with the correct order, the poles can often be positioned accurately in a single step. However, when a lower-order function is used for the fit, multiple iterations may be required. The process becomes more challenging when the frequency response includes nonrational elements, such as noise, which can hinder convergence and, in some cases, even cause it to stall. A significant improvement in VF’s convergence properties can be achieved through a minor modification. As shown above, the classical formulation of VF incorporates a scaling function <i>σ(s)</i> that approaches unity at high frequencies. However, this high-frequency asymptotic condition can negatively impact convergence. To address this issue, the asymptotic condition is replaced with a more flexible constraint that ensures a nontrivial solution to the least-squares problem while avoiding the adverse effects on convergence. Then equation <i>(4)</i> is replaced by <i>(10)</i>:
</p>

![Ecuación](https://quicklatex.com/cache3/5f/ql_b5c5ba191a254cb9f0bf7983de2dc95f_l3.png)

<p align="justify">
where <i>d̂</i> is real. In order to avoid the trivial solution. One equation is added to the resulting LS problem:
</p>

![Ecuación](https://quicklatex.com/cache3/f3/ql_900c54f8394d22d7cccc4cbd67ba42f3_l3.png)

<p align="justify">
Equation <i>(11)</i> ensures that the sum of the real part of <i>σ(s)</i> over the given frequency samples is a nonzero value, while leaving all free variables unrestricted. Since <i>σ(s)</i> during iterations does not approach unity at high frequencies, <i>(7)</i> must be replaced with:
</p>

![Ecuación](https://quicklatex.com/cache3/3d/ql_5a96abc92d502012111829a0c327a33d_l3.png)

<p align="justify">
The zero calculation in <i>(12)</i> is only valid when <i>d̂</i> is nonzero. If the absolute value of <i>d̂</i> is found to be smaller than <i>tol = 1x10<sup>-8</sup></i>, the solution is discarded, and the least-squares problem has to be solved again with a fixed value for <i>d̂</i> in equation <i>(10)</i>.
</p>

### QR algorithm for efficient implementation

## A comprehesive vecfitX.m tutorial

### New features

### Settings

<p align="justify">
Like the original function written by Bjørn Gustavsen, vecfitX.m can be configured according to the user's requirements. Below is a list of the available options for customizing the function. You can set all the parameters or just some of them. Numbers in parentheses indicate the possible values ​​that a parameter can take.
</p>

<ul>
    <li>opt.relax:
        <ul>
            <li>(0) Use classical nontriviality constraint.</li>
            <li>(1) Use relaxed nontriviality constraint.</li>
        </ul>
    </li>
    <li>opt.stable:
        <ul>
            <li>(0) Unstable poles are kept unchanged.</li>
            <li>(1) Unstable poles are made stable by flipping them into the left half-plane.</li>
        </ul>
    </li>
    <li>opt.asymp:
        <ul>
            <li>(1) D and E are omitted in fitting process.</li>
            <li>(2) Only E is omitted in fitting process.</li>
            <li>(3) D and E are taken into account.</li>
        </ul>
    </li>
    <li>opt.skip_pole:
        <ul>
            <li>(1) The pole identification part is skipped. C, D and E are identified using the initial poles as final poles.</li>
        </ul>
    </li>
    <li>opt.skip_res:
        <ul>
            <li>(1) The residue identification part is skipped. Only the poles are identified.</li>
        </ul>
    </li>
    <li>opt.repre:
        <ul>
            <li>(0) The returned state-space model has real elements only. Output variable A is square with <i>2x2</i> submatrices as diagonal elements.</li>
            <li>(1) The returned state-space model has real and complex conjugate parameters. Output variable A is diagonal and sparse.</li>
            <li>(2) The returned model has a residue-pole representation. Output variable A (poles) is a <i>Nx1</i> vector, variable C (residues) is a <i>Nrx(NxNc)</i> array. Variables D and E are <i>NrxNc</i> arrays.</li>
        </ul>
    </li>
    <li>opt.errplot:
        <ul>
            <li>(1) Include deviation in magnitude and phase angle plots.</li>
        </ul>
    </li>
    <li>opt.fitplot:
        <ul>
            <li>(1) Create plots of fitted function compared to the original function. Both magnitude and phase angle are shown.</li>
        </ul>
    </li>
    <li>opt.sigmaplot:
        <ul>
            <li>(1) Create plot of sigma function.</li>
        </ul>
    </li>
    <li>opt.savefig:
        <ul>
            <li>(0) Figures are not saved.</li>
            <li>(1) Save plots in PDF format.</li>
            <li>(2) Save plots in PNG format.</li>
            <li>(3) Save plots in JPEG format.</li>
            <li>(4) Save plots in SVG format.</li>
        </ul>
    </li>
</ul>

<p align="justify">
If you omit the vecfitX.m configuration, the function has the following default parameters:
</p>

```matlab
def.relax = 1;                            % Use vector fitting with relaxed non-triviality constraint
def.stable = 1;                           % Enforce stable poles
def.asymp = 3;                            % Include both D and E  
def.skip_pole = 0;                        % Do not skip pole identification
def.skip_res = 0;                         % Do not skip identification of residues (C,D,E) 
def.repre = 1;                            % Create complex state space representation
def.errplot = 1;                          % Include deviation in magnitude and phase angle plot
def.fitplot = 1;                          % Create plots of fitted and original functions
def.sigmaplot = 0;                        % Exclude plot of sigma function
def.savefig = 0;                          % Figures are not saved
```

### Input and output data

### Iterative implementation

### Test cases

## Contact info

## References
