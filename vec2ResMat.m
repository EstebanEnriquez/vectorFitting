function C = vec2ResMat(f, fm, Nr, Nc, N)

if Nc == 1, C = f.'; return; end
C = zeros(Nr,N*Nc);
if size(f,2) == numel(fm(:,:,1))   % Nonsymmetric matrix case  
    for k = 1:N
        F = zeros(Nr, Nc);
        l = 0;
        for j = 1:Nc
            for i = 1:Nr
                l = l + 1;
                F(i,j) = f(k,l);
            end
        end
        C(:,((k-1)*Nc)+1:k*Nc) = F;
    end
else                               % Symmetric matrix case
    for k = 1:N
        F = zeros(Nr, Nc);
        mask = tril(true(Nr));
        [r, ~] = find(mask);
        sm = 1; j = 0;
        for i = r'
            j = j + 1;
            F(i,sm) = f(k,j);
            if sm ~= i
                F(sm,i) = f(k,j);
            end
            if i == Nr, sm = sm + 1; end
        end
        C(:,((k-1)*Nc)+1:k*Nc) = F;
    end
end