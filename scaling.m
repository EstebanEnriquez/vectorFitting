function [As, escale] = scaling(A)
As = A;
N = size(A,2);
escale = zeros(1,N);
for i = 1:N
    escale(i) = 1./norm(A(:,i));
    As(:,i) = escale(i).*A(:,i);
end