function [SER, ord_zrs, rms, fit] = vectfitX(fm, s, poles, weight, opt)

%         ==============================================================
%         ||   Fast Relaxed Vector Fitting                            ||
%         ||   Version 2.0 in MATLAB                                  ||
%         ||   Last revised: 04-2025                                  ||
%         ||   Algorithm creaded by: Bjorn Gustavsen                  ||
%         ||   bjorn.gustavsen@sintef.no                              ||
%         ||   http://www.energy.sintef.no/Produkt/VECTFIT/index.asp  ||
%         ||   MATLAB code modified by: Esteban Enriquez-Jamangape    ||
%         ||   esteban.enriquez@cinvestav.mx                          ||
%         ||   https://github.com/EstebanEnriquez                     ||
%         ||   Note: RESTRICTED to NON-COMMERCIAL use                 ||
%         ==============================================================


%   =========================================================================
%   ||  + PURPOSE : Approximate f(s) with a state-space or pole-residue    ||
%   ||              model.                                                 ||
%   ||                                                                     ||
%   ||            f(s) = C*(s*I-A)^(-1)*B + D + s*E  (state space model)   ||
%   ||            f(s) = C_n/(s-a_n) + D + s*E       (pole-residue model)  ||
%   ||                                                                     ||
%   ||            where f(s) is a single element or a vector of elements.  ||  
%   ||            When f(s) is a vector, all elements become fitted with   ||
%   ||            a common set of poles.                                   ||
%   ||                                                                     ||
%   ||  + INPUT :                                                          ||
%   ||                                                                     ||
%   ||  f(s) :    Function (3D array) to be fitted; dimension Nr x Nc x Ns ||
%   ||            Nr : number of rows in array                             ||
%   ||            Nc : number of columns in array                          ||
%   ||            Ns : number of layers (freq. samples) in array           ||
%   ||                                                                     ||
%   ||  s :       Vector of frequency points [rad/sec]; dimension 1xNs     || 
%   ||                                                                     ||
%   ||  poles :   Vector of initial poles [rad/sec]; dimension 1xN         ||
%   ||                                                                     ||
%   ||  weight :  The elements in the system matrix are weighted using this||
%   ||            array. Can be used for achieving higher accuracy at      ||
%   ||            desired frequency samples. If no weighting is desired,   ||
%   ||            use unitary weights: weight = ones(1,Ns).                ||
%   ||                                                                     ||
%   ||            1D and 2D arrays are allowed:                            ||
%   ||            dimension: 1xNs -> Common weighting for all elements.    ||
%   ||            dimension: (NrxNc)xNs-> Individual weighting.            ||
%   ||                                                                     ||
%   || opt :      Configuration options                                    ||
%   ||                                                                     ||
%   || + OUTPUT :                                                          ||
%   ||                                                                     ||
%   || SER.A :    NxN sparse matrix A. If repre == 1, A is diagonal and    ||
%   ||            complex. If repre == 0, it is square and real with 2x2   ||
%   ||            submatrices as diagonal elements. If repre == 2, it is   ||
%   ||            a Nx1 vector.                                            ||
%   ||                                                                     ||
%   || SER.B :    Nx1 matrix B. If repre == 1: It is a column of 1's.      ||
%   ||            If repre == 0: It contains 0's, 1's and 2's.             ||
%   ||            If repre == 2: It is omitted.                            ||
%   ||                                                                     ||
%   || SER.C :    NrxN matrix C. If repre == 1: It is complex.             ||
%   ||            If repre == 0: Real-only.                                ||
%   ||            If repre == 2: It is a Nrx(NxNc) array.                  ||
%   ||                                                                     ||
%   || SERD.D :   NrxNc real constant term.                                ||
%   ||                                                                     ||
%   || SERE.E :   NrxNc real proportional term.                            ||
%   ||                                                                     ||
%   || ord_zrs :  1xN matrix, it contains the new poles.                   || 
%   ||                                                                     ||
%   || rms :      root-mean-square error (scalar) of approx. for f(s).     ||
%   ||            (0 is returned if skip_res == 1).                        ||
%   ||                                                                     ||
%   || fit :      (NrxNc)xNs matrix. Rational approximation of f(s).       ||
%   ||            (0 is returned if skip_res == 1). If f(s) is a symmetric ||
%   ||            matrix, then fit has (Nrx(Nr+1)/2) rows and Ns columns.  ||
%   =========================================================================

%   =========================================================================
%   || + CONFIGURATION OPTIONS                                             ||
%   ||                                                                     ||
%   || opt.relax == 1      -> Use relaxed nontriviality constraint.        ||
%   || opt.relax == 0      -> Use classical nontriviality constraint.      ||
%   ||                                                                     ||
%   || opt.stable == 0     -> Unstable poles are kept unchanged.           ||
%   || opt.stable == 1     -> Unstable poles are made stable by flipping   ||
%   ||                        them into the left half-plane.               ||
%   ||                                                                     ||
%   ||                                                                     ||
%   || opt.asymp == 1      -> D and E are omitted in fitting process.      ||
%   || opt.asymp == 2      -> E is omitted in fitting process.             ||
%   || opt.asymp == 3      -> D and E are taken into account.              ||
%   ||                                                                     ||
%   || opt.skip_pole == 1  -> The pole identification part is skipped.     ||
%   ||                        C, D and E are identified using the initial  ||
%   ||                        poles as final poles.                        ||
%   ||                                                                     ||
%   || opt.skip_res == 1   -> The residue identification part is skipped.  ||
%   ||                        Only the poles are identified.               ||
%   ||                                                                     ||
%   || opt.repre == 2      -> The returned model has a residue-pole        ||
%   ||                        representation. Output variable A (poles)    ||
%   ||                        is a Nx1 vector, variable C (residues) is a  ||
%   ||                        Nrx(NxNc) array. D, E are of dim: NrxNc      ||
%   || opt.repre == 1      -> The returned state-space model has real and  || 
%   ||                        complex conjugate parameters. Output         ||
%   ||                        variable A is diagonal and sparse.           ||
%   || opt.repre == 0      -> The returned state-space model has real      ||
%   ||                        elements only. Output variable A is square   ||
%   ||                        with 2x2 submatrices as diagonal elements.   ||
%   ||                                                                     ||
%   || opt.errplot == 1    -> Include deviation in magnitude and phase     ||
%   ||                        angle plots.                                 ||
%   ||                                                                     ||
%   || opt.fitplot == 1    -> Create plots of fitted function compared to  ||
%   ||                        the original function. Both magnitude and    ||
%   ||                        phase angle are shown.                       ||
%   ||                                                                     ||
%   || opt.sigmaplot == 1  -> Create plot of sigma function.               ||
%   ||                                                                     ||
%   || opt.savefig == 0    -> Figures are not saved.                       ||
%   || opt.savefig == 1    -> Save plots in PDF format.                    ||
%   || opt.savefig == 2    -> Save plots in PNG format.                    ||
%   || opt.savefig == 3    -> Save plots in JPEG format.                   ||
%   || opt.savefig == 4    -> Save plots in SVG format.                    ||
%   =========================================================================

%   ===========================================================================  %
%   APPROACH: The identification is done using the pole relocating method known  % 
%   as Vector Fitting [1], with relaxed non-triviality constraint for faster     %
%   convergence and smaller fitting errors [2], and utilization of matrix        %
%   structure for fast solution of the pole identification step [3].             %
%   ===========================================================================  %

%   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  % 
%   IMPORTANT NOTE: The use of this program is limited to NON-COMMERCIAL usage   % 
%   only. If the program code or a modified version is used in a scientific      % 
%   work, then reference should be made to the following:                        %
%                                                                                %
%   [1] B. Gustavsen and A. Semlyen, "Rational approximation of frequency        %      
%       domain responses by Vector Fitting", IEEE Trans. Power Delivery,         %      
%       vol. 14, no. 3, pp. 1052-1061, July 1999.                                %
%                                                                                %
%   [2] B. Gustavsen, "Improving the pole relocating properties of vector        %
%       fitting", IEEE Trans. Power Delivery, vol. 21, no. 3, pp. 1587-1592,     %
%       July 2006.                                                               %
%                                                                                %
%   [3] D. Deschrijver, M. Mrozowski, T. Dhaene, and D. De Zutter,               %
%       "Macromodeling of Multiport Systems Using a Fast Implementation of       %
%       the Vector Fitting Method", IEEE Microwave and Wireless Components       % 
%       Letters, vol. 18, no. 6, pp. 383-385, June 2008.                         %
%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  %


%---------------------------------------- Default parameters ---------------------------------------%

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


if nargin < 5                             % Use default values as opts
    opt = def;
else                                      % Merge default values into opts
    A = fieldnames(def);
    for m = 1:length(A)
        if ~isfield(opt, A{m})
            opt.(A{m}) = def.(A{m});
        end
    end
end

%--------------------------------------- Some sanity checks on data input --------------------------------%
if opt.relax ~= 0 && opt.relax ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.relax: ' num2str(opt.relax)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.stable ~= 0 && opt.stable ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opts.stable: ' num2str(opt.stable)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.asymp ~= 1 && opt.asymp ~= 2 && opt.asymp ~= 3
    disp(['ERROR in vectfitX.m: ==> Illegal value for opts.asymp: ' num2str(opt.asymp)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.skip_pole ~= 0 && opt.skip_pole ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.skip_pole: ' num2str(opt.skip_pole)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.skip_res ~= 0 && opt.skip_res ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.skip_res: ' num2str(opt.skip_res)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.repre ~= 0 && opt.repre ~= 1 && opt.repre ~= 2
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.repre: ' num2str(opt.repre)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.errplot ~= 0 && opt.errplot ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.errplot: ' num2str(opt.errplot)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.fitplot ~= 0 && opt.fitplot ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.fitplot: ' num2str(opt.fitplot)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.sigmaplot ~= 0 && opt.sigmaplot ~= 1
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.sigmaplot: ' num2str(opt.sigmaplot)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if opt.savefig < 0 || opt.savefig > 4
    disp(['ERROR in vectfitX.m: ==> Illegal value for opt.savefig: ' num2str(opt.savefig)]);
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end
if size(s,2) ~= size(fm,3)
    disp('Error in vectfitX.m: ==> Third dimension of f does not match length of s.'); 
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end  
if size(s,2) ~= size(weight,2)
    disp('Error in vectfitX.m: ==> Second dimension of weight does not match length of s.');
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end

%---------------------------------------- Some variables definition --------------------------------------%

Ns = length(s);
N = length(poles);
Nr = size(fm,1);
Nc = size(fm,2);
f = stackM(fm, Nr, Nc, Ns);          % 2D function from 3D array 
Nrr = size(f,1);

if size(weight,1) ~= 1 && size(weight,1) ~= size(f,1)
    disp('Error in vectfit3.m: ==> First dimension of weight is neither 1 nor matches the elements of f.');  
    SER = []; ord_zrs = []; rms = []; fit = [];
    return
end

if opt.asymp == 1
    offs = 0; 
elseif opt.asymp == 2
    offs = 1; 
elseif opt.asymp == 3
    offs = 2;
end

%% POLE IDENTIFICATION STAGE

if opt.skip_pole ~= 1 % Pole identification is not skipped
    
%--------------------------------------- Complex pole identification ---------------------------------%

    ind_poles1 = idComplexPoles(poles, N);
    
%------------------------------------------ Building matrix system -----------------------------------%
    
    % Poles definition
    Dk = polesDef(s, poles, N, Ns, ind_poles1);
    
    % Static blocks of array A
    A = zeros(Ns, (2*N) + offs);
    A(:,1:N) = Dk;
    if opt.asymp == 2 
        A(:,N+1) = [ones(Ns, 1)];
    elseif opt.asymp == 3
        A(:,N+1:N+offs) = [ones(Ns, 1), transpose(s)];
    end

%------------------------------------------ Method with relaxation -----------------------------------%

    if opt.relax == 1

        % Scaling for last row of LS-problem
        W = sum(sum(abs(weight.' .* f.').^2, 2));
        W = sqrt(W) / Ns;

        % Definition of some arrays
        A_red=zeros(Nrr*(N+1),N+1);
        b_red=zeros(Nrr*(N+1),1);
        Dk = [Dk, ones(Ns,1)];

        % Main loop of the pole identification stage
        for n = 1:Nrr
            % Weighting method
            if size(weight,1) == 1
                w = weight.';
            elseif size(weight,1) == Nrr
                w = weight(n,:).';
            end
            Aw = w.*A(:,1:N+offs);
            
            % Dynamic block of A (f dependence)
            for j = 1:N+1
                Aw(1:Ns,N+offs+j) = -w.*Dk(1:Ns,j).*f(n,1:Ns).';    
            end
            Ar = [real(Aw); imag(Aw)];
            
            % Integral criterion for sigma
            if n == Nrr
                Ar(2*Ns+1,N+offs+1:2*N+offs+1)=real(W*sum(Dk(:,1:N+1)));
            end
            
            % QR decomposition
            [Q, R] = qr(Ar,"econ");
            R22 = R(N+offs+1:end, N+offs+1:end);
            A_red((n-1)*(N+1)+1:n*(N+1),:) = R22;

            % Array b
            if n == Nrr
                b_red((n-1)*(N+1)+1:n*(N+1),1) = transpose(Q(end,N+offs+1:end))*Ns*W;
            end
        end

        % Scaling of matrix A
        [A_red, escale] = scaling(A_red);

        % Residue calculation of sigma function as an LS-problem
        x = A_red\b_red;
        x = x.*escale.';
    end
    
%--------------------------------------------- Method with no relaxation ---------------------------------------%
% If "relaxation" was selected and D of sigma was extremely small or large, pole identification is solved again % 
    
    % Tolerances used by relaxed version of vector fitting
    tmax = 1e18; tmin = 1e-18;
    
    if opt.relax == 0 || abs(x(end)) < tmin || abs(x(end)) > tmax
        % Definition of some arrays
        if opt.relax == 0
            dn = 1;
        else
            if x(end) == 0
                dn = 1;
            elseif abs(x(end)) < tmin
                dn = sign(x(end))*tmin;
            elseif abs(x(end)) > tmax
                dn = sign(x(end))*tmax;
            end
        end
        A_red = zeros(Nrr*N,N);
        b_red = zeros(Nrr*N,1);

        % Main loop of the pole identification stage
        for n = 1:Nrr
            % Weighting method
            if size(weight,1) == 1
                w = weight.';
            elseif size(weight,1) == Nrr
                w = weight(n,:).';
            end
            Aw = w.*A;
            
            % Dynamic block of A (f dependence)
            for j = 1:N
                Aw(:,N+offs+j) = -w.*Dk(:,j).*transpose(f(n,:));
            end
            Ar = [real(Aw); imag(Aw)];
            
            % Array b
            b = dn*w.*transpose(f(n, :));
            b = [real(b); imag(b)];
            
            % QR decomposition
            [Q, R] = qr(Ar,"econ");
            R22 = R(N+offs+1:end, N+offs+1:end);
            A_red((n-1)*N+1:n*N,:) = R22;
            b_red((n-1)*N+1:n*N,1) = transpose(Q(:,N+offs+1:end))*b;
        end
        
        % Scaling of matrix A
        [A_red, escale] = scaling(A_red);
        
        % Residue calculation of sigma function as an LS-problem
        x = A_red\b_red;
        x = x.*transpose(escale);
        x = [x;dn];
    end
    % Subtracting sigma d from solution
    dn = x(end);
    x = x(1:end-1);
    
%------------------------------------ Calculation of zeros of sigma function --------------------------------%
    b_sig = ones(N,1);
    d_poles = diag(poles); j = 0;
    for i = 1:N
        j = j+1;
        if j < N
            if imag(d_poles(j,j)) ~= 0
                d_poles(j,j+1) = imag(d_poles(j,j));
                d_poles(j+1,j) = -imag(d_poles(j,j));
                d_poles(j,j) = real(d_poles(j,j));
                d_poles(j+1,j+1) = d_poles(j,j);
                b_sig(j) = 2;
                b_sig(j+1) = 0;    
                j = j+1;
            end
        end
    end
    
    % Eigenvalue problem
    H = d_poles - b_sig*x.'/dn;
    zrs = eig(H);
    
%-------------------------------------- Forcing unstable poles to be stable ---------------------------------%
    if opt.stable == 1
        unstable = real(zrs) > 0;
        zrs(unstable) = zrs(unstable) - 2*real(zrs(unstable));
    end
    zrs = sort(zrs);
    
%------------------------------------------ Allocating real poles first -------------------------------------%
    reals = zrs(imag(zrs) == 0);
    imags = sort(zrs(imag(zrs) ~= 0));
    ord_zrs = [reals;imags];
    ord_zrs = conj(ord_zrs);

else % Pole identification is skipped
    
    ord_zrs = poles;
end

%% RESIDUE IDENTIFICATION STAGE

if opt.skip_res ~= 1 % Residue identification is not skipped

%------------------------------------------ Complex pole identification -------------------------------------%

    ind_poles = idComplexPoles(ord_zrs, N);
    
%-------------------------------------------- Building matrix system ----------------------------------------%
    
    % Poles definition
    Dk = polesDef(s, ord_zrs, N, Ns, ind_poles);
    
    % Static blocks of array A
    A = zeros(Ns, N + offs);
    A(:,1:N) = Dk;
    if opt.asymp == 2
        A(:,N+1:N+offs) = [ones(Ns, 1)];
    elseif opt.asymp == 3
        A(:,N+1:N+offs) = [ones(Ns, 1), transpose(s)];
    end
    
    % Common weighting
    if size(weight,1) == 1
        
        % Array A
        A = weight.'.*A;
        A = [real(A); imag(A)];
        
        % Array b
        b = weight.'.*f.';
        b = [real(b); imag(b)];

        % Scaling of matrix A
        [A, escale] = scaling(A);
        
        % Residue calculation of fitted function as an LS-problem
        X = A\b;
        X = X.*transpose(escale);
        
        % Subtracting C,D,E from solution
        if opt.asymp == 2
            D = X(N+1,:).';
            E = zeros(size(X,2), 1);
        elseif opt.asymp == 3
            D = X(N+1,:).';
            E = X(N+2,:).';
        elseif offs == 0
            D = zeros(size(X,2), 1);
            E = zeros(size(X,2), 1);
        end
        C = X(1:N,:);
    
    % Individual weighting
    elseif size(weight,1) == Nrr
        D = zeros(Nrr,1);
        E = zeros(Nrr,1);
        C = zeros(N,Nrr);
        % Main loop of the residue identification stage
        for i = 1:Nrr
            % Array A
            Aw = weight(i,:).'.*A;
            Ar = [real(Aw); imag(Aw)];
            
            % Array b
            b = weight(i,:).'.*f(i,:).';
            b = [real(b); imag(b)];
            
            % Scaling of matrix A
            [Ar, escale] = scaling(Ar);
            
            % Residue calculation of fitted function as an LS-problem
            X = Ar\b;
            X = X.*transpose(escale);
            
            % Subtracting C,D,E from solution
            if opt.asymp == 2
                D(i) = X(N+1);
            elseif opt.asymp == 3
                D(i) = X(N+1);
                E(i) = X(N+2);
            end
            C(:,i) = X(1:N);
        end
    end
    
    % Changing back to make C complex
    C = real2comp(C,ind_poles);
    
    % Function fitted from the obtained state space
    Dk = zeros(Ns,N);
    fit = zeros(Nrr,Ns);
    for i = 1:N
        Dk(:,i) = 1./(s-ord_zrs(i));
    end
    for i = 1:Nrr
        fit(i,:) = Dk*C(:,i);
        if opt.asymp == 2
            fit(i,:) = fit(i,:) + D(i);
        elseif opt.asymp == 3
            fit(i,:) = fit(i,:) + D(i) + E(i).*s;
        end
    end
    
    % RMS error
    dif = fit - f;
    rms = sqrt(sum(sum(abs(dif.^2))))/sqrt(Nrr*Ns);
    
    % Output as state space representation
    B = ones(N,1);
    if opt.repre == 1   % Complex state space
        A = diag(ord_zrs);
        C = C.';
    elseif opt.repre == 0  % Real state space
        A = diag(ord_zrs);
        [A, B, C] = real2complexSS(A, B, ord_zrs, N, C);
    elseif opt.repre == 2  % Residue-pole representation
        A = ord_zrs;
        D = vec2Mat(D, fm, Nr, Nc);
        E = vec2Mat(E, fm, Nr, Nc);
        C = vec2ResMat(C, fm, Nr, Nc, N);
    end
    SER.A = sparse(A); SER.B = B; SER.C = C; SER.D = D; SER.E = E;
    
%---------------------------------------------------- Plots --------------------------------------------------%
    if opt.fitplot == 1 || (opt.sigmaplot == 1 && opt.skip_pole ~= 1)
        close all;
        freq = s./(2*pi*1i);
    else
        return
    end
    
    if opt.fitplot == 1
        hfig1 = figure (1);
        if opt.errplot == 1
            subplot(2,1,1)
            h1 = loglog(freq, abs(f),"b-"); hold on;
            h2 = loglog(freq, abs(fit),"r-."); hold off; grid on;
            legend([h1(1) h2(1)], "Data", "VF"); xlim([freq(1) freq(end)])
            xlabel("Frequency [Hz]"); ylabel("Magnitude")
            subplot(2,1,2)
            loglog(freq, abs(f-fit), 'm'); grid on;
            xlabel("Frequency [Hz]"); ylabel("Deviation"); xlim([freq(1) freq(end)])
        else
            h1 = loglog(freq, abs(f),"b-"); hold on;
            h2 = loglog(freq, abs(fit),"r-."); hold off; grid on;
            legend([h1(1) h2(1)], "Data", "VF")
            xlabel("Frequency [Hz]"); ylabel("Magnitude"); xlim([freq(1) freq(end)])
        end
        
        hfig2 = figure (2);
        if opt.errplot == 1
            subplot(2,1,1)
            ph1 = 180*unwrap(angle(f), [], 2)/pi;
            ph2 = 180*unwrap(angle(fit),[], 2)/pi;
            h4 = plot(freq,ph1,'b'); hold on
            h5 = plot(freq,ph2,'r-.'); hold off; grid on;
            legend([h4(1) h5(1)], "Data", "VF"); xlim([freq(1) freq(end)])
            xlabel("Frequency [Hz]"); ylabel("Phase angle [deg]")
            subplot(2,1,2)
            semilogy(freq, abs(ph1-ph2), 'm'); grid on;
            xlabel("Frequency [Hz]"); ylabel("Deviation"); xlim([freq(1) freq(end)])
        else
            ph1 = 180*unwrap(angle(f), [], 2)/pi;
            ph2 = 180*unwrap(angle(fit),[], 2)/pi;
            h4 = plot(freq,ph1,'b'); hold on
            h5 = plot(freq,ph2,'r-.'); hold off; grid on;
            legend([h4(1) h5(1)], "Data", "VF")
            xlabel("Frequency [Hz]"); ylabel("Phase angle [deg]"); xlim([freq(1) freq(end)])
        end
    end
    
    if opt.sigmaplot == 1 && opt.skip_pole ~= 1
        Dk = zeros(Ns, N);
        for i = 1:N
            Dk(:,i) = 1./(s-poles(i));
        end
        x_c = real2comp(x,ind_poles1);
        sigma = dn + Dk*x_c;
        hfig3 = figure (3);
        h3 = loglog(freq, abs(sigma),"b-"); grid on;
        legend([h3(1)], "Sigma")
        xlabel("Frequency [Hz]"); ylabel("Magnitude"); xlim([freq(1) freq(end)])
    end
    
    if opt.fitplot == 1 && opt.sigmaplot == 1 && opt.skip_pole ~= 1
        hfig = [hfig1, hfig2, hfig3];
    elseif opt.fitplot == 1
        hfig = [hfig1, hfig2];
    elseif opt.sigmaplot == 1 && opt.skip_pole ~= 1
        hfig = hfig3;
    end
    
    nplots = cell(1, size(hfig,2));
    for k = 1:size(nplots,2)
        nplots{k} = sprintf("myfigure_%d", k);
        picturewidth = 30;   % Set this parameter and keep it forever (if you want consistent figures in a document)
        hw_ratio = 0.75;     % Feel free to play with this ratio
        set(findall(hfig(k),'-property','FontSize'),'FontSize', 17)
        set(findall(hfig(k),'-property','Interpreter'),'Interpreter','latex') 
        set(findall(hfig(k),'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
        set(hfig(k),'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
        pos = get(hfig(k),'Position');
        set(hfig(k),'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
        if opt.savefig == 1
            print(hfig(k),nplots{k},'-dpdf','-vector','-fillpage')
        elseif opt.savefig == 2
            print(hfig(k),nplots{k},'-dpng','-r300')
        elseif opt.savefig == 3
            print(hfig(k),nplots{k},'-djpeg','-r300')
        elseif opt.savefig == 4
            print(hfig(k),nplots{k},'-dsvg')
        end
    end

else % Residue identification is skipped

    % Output as state space representation
    B = ones(N,1);
    if opt.repre == 1 || opt.repre == 2      % Complex state space
        A = diag(ord_zrs);
        C = 0;
    elseif opt.repre == 0  % Real state space
        A = diag(ord_zrs);
        [A, B, C] = real2complexSS(A, B, ord_zrs, N);
    end
    SER.A = sparse(A); SER.B = B; SER.C = C; SER.D = 0; SER.E = 0;
    rms = 0; fit = zeros(Nrr, Ns);
end