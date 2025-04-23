%  ex2.m             
% -Fitting a measured admittance function from a distribution transformer (single element).
% -Reading frequency response f(s) from 03pk10.txt file.
% -Initial poles: 15 linearly spaced complex pairs (N = 30).
% -5 iterations.
%
% This example script was taken from the vector fitting package "VFIT3.zip" 
% Created by: Bjorn Gustavsen.

clear; clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

%% Reading frequency response data from file
fid1 = fopen('03PK10.txt','r');
A1 = fscanf(fid1,'%f',1);
A2 = fscanf(fid1,'%f',1); 

f = zeros(1,1,160);
for k = 1:160
    A1 = fscanf(fid1,'%f',1);
    A2 = fscanf(fid1,'%f',1);
    f1 = real(A1*exp(1i*A2*pi/180));
    f2 = imag(A1*exp(1i*A2*pi/180));
    f(1,1,k) = f1 + 1i*f2;
end

%% Linearly spaced frequency samples
w = 2*pi*linspace(0,10e6,401);                 % Angular frequency array 
w(1) = []; w(161:400) = [];
Ns = length(w);                                % Number of frequency samples
s = 1i.*w;                                     % Complex frequency array

%% Starting poles computation
N = 30;                                        % Order of approximation 
bet = linspace(w(1), w(end), N/2);
poles = zeros(1,N);
alf = -bet*1e-2;
poles(1,1:2:end) = complex(alf, -bet);        % Complex starting poles
poles(1,2:2:end) = complex(alf, bet);

%% vecfitX configuration (optional)
opt.skip_res = 1;
opt.savefig = 4;

%% Iterative implementation of vecfitX
fw = stackM(f, 1, 1, Ns);
weight = 1./abs(fw);                          % Weighting with inverse of magnitude function
Niter = 5;                                    % Number of iterations
for iter = 1:Niter
  if iter == Niter
      opt.skip_res = 0;
  end
  [SER, poles, rms, fit] = vectfitX(f, s, poles, weight, opt);
end
