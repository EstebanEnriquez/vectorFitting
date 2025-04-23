function Dk = polesDef(s, poles, N, Ns, ind_poles)

% Poles definition
Dk = zeros(Ns, N);
for i = 1:N
    if ind_poles(i) == 0
        Dk(:,i) = 1./(s-poles(i));
    elseif ind_poles(i) == 1
        Dk(:,i) = (1./(s-poles(i)))+(1./(s-poles(i)'));
        Dk(:,i+1) = (1i./(s-poles(i)))-(1i./(s-poles(i)'));
    end
end