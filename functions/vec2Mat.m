function F = vec2Mat(f, fm, Nr, Nc)

if Nc == 1, F = f; return; end

if size(f,1) == numel(fm(:,:,1))   % Nonsymmetric matrix case
    F = zeros(Nr, Nc);
    k = 0;
    for j = 1:Nc
        for i = 1:Nr
            k = k + 1;
            F(i,j) = f(k);
        end
    end
else                               % Symmetric matrix case
    F = zeros(Nr, Nc);
    mask = tril(true(Nr));
    [r, ~] = find(mask);
    sm = 1; j = 0;
    for i = r'
        j = j + 1;
        F(i,sm) = f(j);
        if sm ~= i
            F(sm,i) = f(j);
        end
        if i == Nr, sm = sm + 1; end
    end
end