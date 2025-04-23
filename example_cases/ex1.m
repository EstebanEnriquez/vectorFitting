%  ex1.m
% -Creating an 18th order frequency response f(s) of 2 elements (row vector). 
% -Initial poles: 9 linearly spaced complex pairs (N = 18).
% -3 iterations.
%
% This example script was taken from the vector fitting package "VFIT3.zip" 
% Created by: Bjorn Gustavsen.

clear; clc;

%%  Constructing the artificial scalar function
p = [-4500 -41000 (-100+1j*5e3) (-100-1j*5e3) (-120+1j*15e3) (-120-1j*15e3) (-3e3+1j*35e3) (-3e3-1j*35e3)];
r = [-3000 -83000 (-5+1j*7e3) (-5-1j*7e3) (-20+1j*18e3) (-20-1j*18e3) (6e3+1j*45e3) (6e3-1j*45e3)];
p = [p (-200+1j*45e3) (-200-1j*45e3) (-1500+1j*45e3) (-1500-1j*45e3)];
r = [r (40 +1j*60e3) (40-1j*60e3) (90 +1j*10e3) (90-1j*10e3)];
p = [p (-5e2+1j*70e3) (-5e2-1j*70e3) (-1e3+1j*73e3) (-1e3-1j*73e3) (-2e3+1j*90e3) (-2e3-1j*90e3)];
r = [r (5e4+1j*80e3) (5e4-1j*80e3) (1e3+1j*45e3) (1e3-1j*45e3) (-5e3+1j*92e3) (-5e3-1j*92e3)];

p = 2*pi*p; r = 2*pi*r;
p1 = p(1:10); r1 = r(1:10); N1 = length(p1);
p2 = p(9:18); r2 = r(9:18);
D = 0.2; E = 2e-5;

%% Linearly spaced frequency samples
w = 2*pi*linspace(1,1e5,100);               % Angular frequency array
Ns = length(w);                             % Number of frequency samples
s = 1j*w;                                   % Complex frequency array

%% Artificial function
f = zeros(1,2,Ns);
for k = 1:Ns
    for n = 1:N1
        f(1,1,k) = f(1,1,k) + r1(n)/(s(k)-p1(n));
        f(1,2,k) = f(1,2,k) + r2(n)/(s(k)-p2(n));
    end
    f(1,1,k) = f(1,1,k) + s(k)*E;
    f(1,2,k) = f(1,2,k) + s(k)*3*E;
end
f(1,1,:) = f(1,1,:) + D;
f(1,2,:) = f(1,2,:) + 2*D;

%% Starting poles computation
N = 18;                                  % Order of approximation
bet = linspace(w(1), w(end), N/2);
poles = zeros(1,N);
alf = -bet*1e-2;
poles(1,1:2:end) = complex(alf, -bet);   % Complex starting poles
poles(1,2:2:end) = complex(alf, bet);

%% vecfitX configuration (optional)
opt.skip_res = 1;
opt.savefig = 4;

%% Iterative implementation of vecfitX
weight = ones(1,Ns);              % No weighting method is used 
Niter = 3;                        % Number of iterations
for iter = 1:Niter
    if iter == Niter
      opt.skip_res = 0;
    end
    [SER, poles, rms, fit] = vectfitX(f, s, poles, weight, opt);
end