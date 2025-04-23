function C = real2comp(X,ind_poles)
N = length(ind_poles);
C = X(1:N,:);
for i = 1:N
    if ind_poles(i) == 1
        for j = 1:size(C,2)
            C(i,j) = complex(X(i,j),X(i+1,j));
            C(i+1,j) = complex(X(i,j),-X(i+1,j));
        end
    end
end
