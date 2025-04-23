function [ff, rms] = res2fit(s, SER, bigY, Nr, Nc, pl)

N = size(SER.A, 1);
Ns = size(s,2);

% Fitting function from residue-pole representation
faux = zeros(Nr, Nc, Ns); ff = faux;
for i = 1:N
    cn = SER.C(:,((i-1)*Nc)+1:i*Nc);
    sa = (s-SER.A(i));
    for j = 1:Ns
        faux(:,:,j) = cn./sa(j);
        if i == N
            faux(:,:,j) = faux(:,:,j) + SER.D + s(j).*SER.E; 
        end
    end
    ff = ff + faux;
end

% RMS error
dif = ff - bigY;
rms = sqrt(sum(sum(sum(abs(dif.^2)))))/sqrt(Nc*Nr*Ns);

% Plots
if pl == 1
    fitM = zeros(Nc*Nr, Ns);
    fM = zeros(Nc*Nr, Ns);
    for i = 1:Nc
        fitM(((i-1)*Nr)+1:i*Nr,:) = ff(:,i,:);
        fM(((i-1)*Nr)+1:i*Nr,:) = bigY(:,i,:);
    end
    
    figure (1)
    freq = s./(2*pi*1i);
    subplot(2,1,1)
    h1 = loglog(freq, abs(fM),"b-"); hold on;
    h2 = loglog(freq, abs(fitM),"r-."); hold off; grid on;
    legend([h1(1) h2(1)], "Data", "VF"); xlim([freq(1) freq(end)])
    xlabel("Frequency [Hz]"); ylabel("Magnitude")
    subplot(2,1,2)
    loglog(freq, abs(fM-fitM), 'm'); grid on;
    xlabel("Frequency [Hz]"); ylabel("Deviation"); xlim([freq(1) freq(end)])
    
    figure (2)
    subplot(2,1,1)
    ph1 = 180*unwrap(angle(fM), [], 2)/pi;
    ph2 = 180*unwrap(angle(fitM),[], 2)/pi;
    h4 = plot(freq,ph1,'b'); hold on
    h5 = plot(freq,ph2,'r-.'); hold off; grid on;
    legend([h4(1) h5(1)], "Data", "VF"); xlim([freq(1) freq(end)])
    xlabel("Frequency [Hz]"); ylabel("Phase angle [deg]")
    subplot(2,1,2)
    loglog(freq, abs(ph1-ph2), 'm'); grid on;
    xlabel("Frequency [Hz]"); ylabel("Deviation"); xlim([freq(1) freq(end)])
end