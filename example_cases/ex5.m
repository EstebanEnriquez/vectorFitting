%  ex5.m
% -Fitting all elements of the admittance matrix of a 6-terminal system (power system distribution network)
% -Reading frequency admittance matrix Y(s) from fdne.txt file
% -Initial poles: 25 linearly spaced complex pairs (N = 50)
% -Fitting g(s) = sum(f(s)) in order to get a good set of initial poles (5 iterations)
% -Fitting elements of vector f(s) using a common pole set (3 iterations)
% -Residue-pole representation
%
% This example script was taken from the vector fitting package "VFIT3.zip" 
% Created by: Bjorn Gustavsen.

clear,clc;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

%% Reading frequency response data from file
fid1 = fopen('fdne.txt','r');
Nr = fscanf(fid1,'%f',1); Nc = Nr;
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
s = 1i.*w;
fclose(fid1);

%% Starting poles computation
N = 50;                                        % Order of approximation 
bet = linspace(w(1), w(end), N/2);
poles = zeros(1,N);
alf = -bet*1e-2;
poles(1,1:2:end) = complex(alf, -bet);         % Complex starting poles
poles(1,2:2:end) = complex(alf, bet);

%% vecfitX configuration (optional)
opt.skip_res = 1;
opt.repre = 2;
opt.fitplot = 0;

%% Improvement of initial poles

% Forming weighted column sum
g = 0;
gn = zeros(1,1,Ns);
f = stackM(bigY, Nr, Nc, Ns);
for n = 1:Nc
    %g = g + f(n,:);               % Unweighted sum     
    g = g + f(n,:)/norm(f(n,:));
    %g = g + f(n,:)/sqrt(norm(f(n,:)));     
end
weight_g = 1./abs(g);
gn(1,1,:) = g(:);

Niter1 = 5; % Iterations for pole set improvement 
for iter = 1:Niter1
    [~,poles] = vectfitX(gn,s,poles,weight_g,opt);
end

%% Weighting techniques
% weight = ones(1,Ns);
weight = 1./abs(f);
%weight = 1./sqrt(abs(f(1,:)));

%% Iterative implementation of vecfitX
Niter2 = 3; % Iterations for main fitting
for iter = 1:Niter2
    if iter == Niter2
        opt.skip_res = 0; 
    end
    [SER,poles] = vectfitX(bigY,s,poles,weight,opt);
end

%% Fitted matrix from residue-pole model
pl = 1; % If pl = 1, res2fit creates plots of fitted function compared to the original f(s). Both magnitude and phase angle are shown.
[fs_fit, rms] = res2fit(s, SER, bigY, Nr, Nc, pl);
