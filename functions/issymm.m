function y = issymm(A)

if size(A,1) == size(A,2)
    if norm(A - A.', 'fro') < 1e-10
        y = true;
    else
        y = false;
    end
else
    y = false;
end