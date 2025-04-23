%  ex4.m
% -Fitting all columns of the admittance matrix of a 6-terminal system (power system distribution network)
% -Reading frequency admittance matrix Y(s) from fdne.txt file
% -Extracting 1st column: f(s) (contains 6 elements)
% -Initial poles: 25 linearly spaced complex pairs (N = 50)
% -Fitting g(s) = sum(f(s)) in order to get a good set of initial poles (5 iterations)
% -Fitting the six-columns one-by-one (3 iterations)
% -Residue-pole representation
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
opt.fitplot = 0;
opt.repre = 2;

%% Arrays for storing the residue-pole representation
A = zeros(N, Nr);
C = zeros(Nr, N*Nr);
D = zeros(Nr, Nr);
E = zeros(Nr, Nr);

%% Main loop
Niter1 = 5; % Iterations for pole set improvement 
Niter2 = 3; % Iterations for main fitting
f = zeros(Nr,1,Ns);
for col = 1:Nr
    
    % Extracting elements in each column
    for n=1:Nr
        f(n,1,:) = bigY(n,col,:);
    end
    
    % Skip residue identification stage
    opt.skip_res = 1;
    
    %% Improvement of initial poles using the first column
    if col == 1
        g = zeros(1,1,Ns);
        for n = 1:Nr
            g = g + f(n,1,:)/norm(squeeze(f(n,1,:)));
        end
        weight_g = 1./abs(squeeze(g).');
        for iter = 1:Niter1
            [~,poles] = vectfitX(g,s,poles,weight_g,opt);
        end
    end

    %% Weighting techniques
    %weight = ones(1,Ns);
    %weight = 1./abs(squeeze(f));
    weight = 1./sqrt(abs(squeeze(f)));
    
    %% Iterative implementation of vecfitX
    for iter = 1:Niter2
        if iter == Niter2
            opt.skip_res = 0;
        end
        [SER,poles,rms,fit] = vectfitX(f,s,poles,weight,opt);
    end

    %% Residue-pole representation of the entire system
    A(:,col) = SER.A;
    D(:,col) = SER.D;
    E(:,col) = SER.E;
    C(:,((0:N-1)*Nr) + col) = SER.C;
end

%% Plots
nc = 1; % Select the column number you want to plot
SER.A = A(:,nc); SER.C = C(:,((0:N-1)*Nr) + nc); SER.D = D(:,nc); SER.E = E(:,nc);
pl = 1; % If pl = 1, res2fit creates plots of fitted function compared to the original f(s). Both magnitude and phase angle are shown.
[fs_fit, rms2] = res2fit(s, SER, bigY(:,nc,:), Nr, 1, pl);