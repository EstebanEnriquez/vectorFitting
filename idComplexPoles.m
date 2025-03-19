function ind_poles = idComplexPoles(poles, N)

ind_poles = zeros(1,N);
for i = 1:N
    if imag(poles(i)) ~= 0 && ind_poles(i) == 0  
        if abs(imag(poles(i))) == abs(imag(poles(i+1)))
            ind_poles(i) = 1;
            ind_poles(i+1) = 2;
        end
    end
end