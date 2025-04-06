function faux = stackM(f, Nr, Nc, Ns)

if Ns ~= size(f,3)
    disp('Error in vectfit3.m: ==> Third dimension of f does not match length of s.');
    faux = [];
    return
end

if issymm(f(:,:,1)) && issymm(f(:,:,end))   % Symmetric matrix case
    faux = zeros(Nr*(Nr+1)/2, Ns);
    mask = tril(true(Nr));
    [r, ~] = find(mask);
    sm = 1; j = 0;
    for i = r'
        j = j + 1;
        faux(j,:) = squeeze(f(i,sm,:));
        if i == Nr, sm = sm + 1; end
    end
else                                        % Nonsymmetric matrix case

    if Nc ~= 1 && Nr ~= 1
        faux = zeros(Nr*Nc, Ns);
        k = 0;
        for i = 1:Nc
            for j = 1:Nr
                k = k + 1;
                faux(k,:) = squeeze(f(j,i,:));
            end
        end
    elseif Nc == 1 && Nr == 1
        faux = squeeze(f).';
    else
        faux = squeeze(f);
    end
end