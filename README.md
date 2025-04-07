# The latest revised Fast-Relaxed Vector Fitting for MATLAB

## Vector Fitting Theory

### Classical Vector Fitting

<p align="justify">
One of the most utilized rational approximation techniques is Vector Fitting (VF). In 1999, VF was introduced and proved to be a highly robust and efficient rational approximation method, applicable to both smooth and resonant responses [1]. This caused VF to be rapidly adopted in several areas including power systems and macro-modeling systems. VF approximates a frequency response (generally an array) using a rational function approximation, as shown in <i>(1)</i>:
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
VF operates by iteratively adjusting an initial set of poles to more optimal locations. When fitting the frequency-domain response of a rational function with the correct order, the poles can often be positioned accurately in a single step. However, when a lower-order function is used for the fit, multiple iterations may be required. The process becomes more challenging when the frequency response includes nonrational elements, such as noise, which can hinder convergence and, in some cases, even cause it to stall. A significant improvement in VF’s convergence properties can be achieved through a minor modification. As shown above, the classical formulation of VF incorporates a scaling function <i>σ(s)</i> that approaches unity at high frequencies. However, this high-frequency asymptotic condition can negatively impact convergence. To address this issue, the asymptotic condition is replaced with a more flexible constraint that ensures a nontrivial solution to the least-squares problem while avoiding the adverse effects on convergence [2]. Then equation <i>(4)</i> is replaced by <i>(10)</i>:
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

<p align="justify">
For multiport systems, the classical VF method requires solving large, sparse linear systems to estimate the rational function, as shown in equation (6). This becomes computationally expensive and memory-intensive as the number of ports grows. For example, a 60-port system with 101 frequency samples requires 54 GB of RAM in the standard implementation, making it impractical for typical workstations. Even if the structure under study has a moderate amount of ports, the size of the corresponding LS matrix may become prohibitively large. Furthermore, a lot of computational effort is wasted on the calculation of the residues <i>c</i><sub><i>n</i></sub>, which are discarded by the VF algorithm. In [3], the authors introduce a Fast VF method that leverages the QR decomposition to simplify the LS equations. This process leads to a simplified set of equations which depend only on <i>ĉ</i><sub><i>n</i></sub>. Instead of solving one large linear system for all ports, the method processes each matrix element (port) sequentially. For each element, the QR decomposition is applied to its corresponding linear system. This decomposes the problem into smaller, manageable parts. After decomposing all individual systems, the results are combined into a smaller, shared system.
</p>

<p align="justify">
As explained above, once the residues <i>ĉ</i><sub><i>n</i></sub> of <i>σ(s)</i> are computed, the common poles <i>â</i><sub><i>n</i></sub> of the transfer function are found by computing the eigenvalues of <i>(12)</i>. Note that the second stage of VF remains unchanged.
</p>

<p align="justify">
The modification of the VF with the two methods proposed in [2] and [3] have resulted in the third and last version of VF so far, better known as  <strong>Fast-Relaxed Vector Fitting (FRVF)</strong> [4].
</p>

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

### Function call

<p align="justify">
    The following are some valid calls to the vecfitX function. The inputs and outputs are explained in detail below.
</p>

```matlab
[SER,poles,rmserr,fit] = vectfitX(f,s,poles,weight,opt);
[SER,poles,rmserr,fit] = vectfitX(f,s,poles,weight);
[SER,poles,rmserr] = vectfitX(f,s,poles,weight);
```

### Input and output data

#### Description of the input data:

<ul>
    <li><strong><code>f</code>:</strong> Matrix <i>f(s)</i> function (3D array) to be fitted with dimensions <i>Nr X Nc X Ns</i>.
        <ul>
            <li><i>Nr</i>: Number of rows in array.</li>
            <li><i>Nc</i>: Number of columns in array.</li>
            <li><i>Ns</i>: Number of layers (frequency samples) in array.</li>
        </ul>
    </li>
    <li><strong><code>s</code>:</strong> Vector of frequency samples [rad/sec] with dimensions <i>1 X Ns</i>.</li>
    <li><strong><code>poles</code>:</strong> Vector of initial poles [rad/sec] with dimensions <i>1 X N</i>.  (Explain selection of inital poles)</li>
    <li><strong><code>weight</code>:</strong> The elements in the system matrix are weighted using this array. It can be used for achieving higher accuracy at desired frequency samples. If no weighting is desired, use unitary weights, i.e. weight array of ones, with dimensions <i>1 X Ns</i>. Otherwise, 1D and 2D arrays are allowed. (Explain computation of weight and stackM.m function).</li>
        <ul>
            <li><strong>1D array</strong>: Common weighting for all elements, weight array with dimensions <i>1 X Ns</i>.</li>
            <li><strong>2D array</strong>: Individual weighting, weight array with dimensions <i>(Nr*Nc) X Ns</i>.</li>
        </ul>
    <li><strong><code>opt</code>:</strong> Configuration options.</li>
</ul>

#### Description of the output data:

<ul>
    <li><strong><code>SER.A</code>:</strong> Array containing the poles of the fitted matrix function.
        <ul>
            <li> If <strong><code>opt.repre = 0</code>:</strong> A is sparse, square and real with <i>2 X 2</i> submatrices as diagonal elements, with dimensions <i>N X N</i>.</li>
            <li> If <strong><code>opt.repre = 1</code>:</strong> A is sparse, diagonal and complex, with dimensions <i>N X N</i>.</li>
            <li> If <strong><code>opt.repre = 2</code>:</strong> A is a complex vector with dimensions <i>N X 1</i>.</li>
        </ul>
    </li>
    <li><strong><code>SER.B</code>:</strong> Nx1 matrix only for real and complex state space representation.
        <ul>
            <li> If <strong><code>opt.repre = 0</code>:</strong> B is a column vector of ones.</li>
            <li> If <strong><code>opt.repre = 1</code>:</strong> B is a column vector of zeros, ones and twos.</li>
        </ul>
    </li>
    <li><strong><code>SER.C</code>:</strong> Array containing the residues of the fitted matrix function.
        <ul>
            <li> If <strong><code>opt.repre = 0</code>:</strong> C is real with dimensions <i>(Nr*Nc) X N</i>.</li>
            <li> If <strong><code>opt.repre = 1</code>:</strong> C is complex with dimensions <i>(Nr*Nc) X N</i>.</li>
            <li> If <strong><code>opt.repre = 2</code>:</strong> C is complex with dimensions <i>Nr X (N*Nc)</i>.</li>
        </ul>
    </li>
    <li><strong><code>SER.D</code>:</strong> <i>Nr X Nc</i> real constant term.</li>
    <li><strong><code>SER.E</code>:</strong> <i>Nr X Nc</i> real proportional term.</li>
    <li><strong><code>ord_zrs</code>:</strong> <i>1 X N</i> matrix containing the newly calculated poles.</li>
    <li><strong><code>rms</code>:</strong> root-mean-square error (scalar) of approximation for <i>f(s)</i>. 0 is returned if <code>skip_res = 1</code>.</li>
    <li><strong><code>fit</code>:</strong> <i>(Nr*Nc) X Ns</i> array containing the rational approximation of <i>f(s)</i>. 0 is returned if <code>skip_res = 1</code>. If <i>f(s)</i> is a symmetric matrix, then fit has <i>(Nr*(Nr+1)/2)</i> rows and <i>Ns columns</i>.</li>
</ul>


### Iterative implementation

<p align="justify">
    A common practice is to iteratively apply the FRVF to improve the relocation of the calculated poles and achieve higher accuracy. This is achieved as follows:
</p>

```matlab
opt.skip_res = 1;          % Residue identification part is skipped.
opt.repre = 2;             % Pole-residue representation.
opt.fitplot = 1;           % Include graphs of fitted function.

Niter = 5;                 % Number of the iterations.
for iter = 1:Niter
  if iter == Niter
      opt.skip_res = 0;    % In the last iteration residue identification part is NOT skipped.
  end
    [SER,poles,rmserr,fit] = vectfitX(bigY,s,poles,weight,opt);
end
```
<p align="justify">
    It is important to note that residue identification is skipped until the last iteration. This is to reduce computation time. Keep in mind that residue identification should be performed with the last set of poles calculated. When residue identification is disabled, graph creation is also disabled. Even if it has been manually enabled, as in the previous example. This way, only the plots corresponding to the last iteration are displayed.
</p>

### Test cases

<p align="justify">
    In order to demonstrate the use of the vecfitX function, test cases have been included in this repository. The examples and data have been taken and adapted from [4]. A brief description of these cases is shown below:
</p>

<ul>
  <li><strong><code>ex1.m</code>:</strong> In this program, an artificial scalar function <i>f(s)</i> is constructed. <i>f(s)</i> is a row vector whose elements correspond to two predefined partial fractions. In addition, 18 complex starting poles are used, vecfitX is called three times iteratively (new poles are used as starting poles in each iteration), and a real state-space representation is selected.</li>
    <li><strong><code>ex2.m</code>:</strong> In this example, the measured response of a transformer <i>f(s)</i> (scalar) is read from file <code>03pk10.txt</code>. In addition, 6 complex starting poles are used, vecfitX is called five times iteratively (new poles are used as starting poles in each iteration), a complex state-space representation is selected, and <i>f(s)</i> is weighted with the inverse of its magnitude.</li>
    <li><strong><code>ex3.m</code>:</strong> In this program, a calculated terminal admittance <i>6 X 6 X 300</i> matrix <i>Y(s)</i> of a power system is read from file <code>fdne.txt</code>. Only the first column of <i>Y(s)</i> is fitted, as follows: 50 complex starting poles are used, vecfitX is called five times iteratively (new poles are used as starting poles in each iteration), a complex state-space representation is selected, and <i>Y(s)</i> is weighted with the inverse of the square root of its magnitude.</li>
    <li><strong><code>ex4.m</code>:</strong> In this example, the same <i>6 X 6 X 300</i> matrix <i>Y(s)</i> is used. Now a complete complex state-space representation is calculated by fitting <i>Y(s)</i> column by column (each column has a different set of poles). In addition, 50 complex starting poles are used and improved on by fitting the element sum of the first column before fitting all columns, vecfitX is called 3 times iteratively (for each column),  and the columns of <i>Y(s)</i> are weighted with the inverse of the square root of its magnitudes.</li>
    <li><strong><code>ex5.m</code>:</strong> Again the same <i>6 X 6 X 300</i> matrix <i>Y(s)</i> is used. Now a residue-pole representation is calculated by fitting <i>Y(s)</i> (all elements has a common set of poles). The fitting is performed similarly to the previous example. After obtaining the residue-pole representation, the <code>res2fit</code> function is used to obtain the fitted matrix in its original shape. </li>
</ul>

## Contact info
<p align="justify">
If you would like to get in touch for professional collaboration, research proposals, or questions related to my work, feel free to reach out through any of the following channels:
</p>
<ul>
  <li><strong>Primary email:</strong> josejamangape@gmail.com</li>
  <li><strong>Secondary email:</strong> esteban.enriquez@cinvestav.com</li>
  <li><strong>LinkedIn:</strong> <a href="https://www.linkedin.com/in/jose-esteban-enriquez-jamangape-9a024131a/">Jose Esteban Enriquez</a></li>
  <li><strong>GitHub:</strong> <a href="https://github.com/EstebanEnriquez">Esteban Enriquez</a></li>
  <li><strong>ResearchGate:</strong> <a href="https://www.researchgate.net/profile/J-Enriquez-Jamangape?ev=hdr_xprf">J. E. Enríquez-Jamangape</a></li>
</ul>

## References
<p align="justify">
[1]  B. Gustavsen and A. Semlyen, "Rational approximation of frequency domain responses by Vector Fitting", IEEE Trans. Power Delivery, vol. 14, no. 3, pp. 1052-1061, July 1999.

[2] B. Gustavsen, "Improving the pole relocating properties of vector fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592, July 2006.

[3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter, “Macromodeling of Multiport Systems Using a Fast Implementation of the Vector Fitting Method”, IEEE Microwave and Wireless Components Letters, vol. 18, no. 6, pp. 383-385, June 2008.

[4] B. Gustavsen, "User's Guide for vectfit3.m (Fast, Relaxed Vector fitting)", SINTEF Energy Research, N-7465 Trondheim, Norway, 2008. Aviable online: <a href="https://www.sintef.no/en/software/vector-fitting/downloads/#menu">vecfit3.m for MATLAB</a>.
</p>

