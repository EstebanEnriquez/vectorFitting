%  ex3.m
% -Fitting 1st column of the admittance matrix of 6-terminal system (power system distribution network)
% -Reading frequency admittance matrix Y(s) from fdne.txt file
% -Extracting 1st column: f(s) (contains 6 elements)
% -Initial poles: 25 linearly spaced complex pairs (N = 50)
% -5 iterations
%
% This example script was taken from the vector fitting package "VFIT3.zip" 
% Created by: Bjorn Gustavsen.

clear; clc;

%% Reading frequency response data from file
fid1 = fopen('fdne.txt','r');
Nr = fscanf(fid1,'%f',1);
Ns = fscanf(fid1,'%f',1);
bigY = zeros(Nr,Nr,Ns); w=zeros(1,Ns);
for k = 1:Ns
    w(k) = fscanf(fid1,'%e',1);
    for row = 1:Nr
        for col = 1:Nr
            dum1 = fscanf(fid1,'%e',1);
            dum2 = fscanf(fid1,'%e',1);   
            bigY(row,col,k) = dum1+1j*dum2;
        end
    end
end
fclose(fid1);

%% Extracting first column
f = zeros(Nr,1,Ns);
for n = 1:Nr
    f(n,1,:) = bigY(n,1,:);
end

%% Linearly spaced frequency samples
Ns = length(w);                                % Number of frequency samples
s = 1i.*w;                                     % Complex frequency array

%% Starting poles computation
N = 50;                                        % Order of approximation 
bet = linspace(w(1), w(end), N/2);
poles = zeros(1,N);
alf = -bet*1e-2;
poles(1,1:2:end) = complex(alf, -bet);         % Complex starting poles
poles(1,2:2:end) = complex(alf, bet);

%% vecfitX configuration (optional)
opt.skip_res = 1;
opt.savefig = 4;

%% Iterative implementation of vecfitX
fw = stackM(f, Nr, 1, Ns);
weight = 1./sqrt(abs(fw));                % Weighting with inverse of the square root of magnitude function
Niter = 5;                                % Number of iterations
for iter = 1:Niter
    if iter == Niter
        opt.skip_res = 0;
    end
    [SER,poles,rmserr,fit] = vectfitX(f,s,poles,weight,opt);
end