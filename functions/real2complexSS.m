function [A, B, C] = real2complexSS(A, B, poles, N, C)

ind_poles = idComplexPoles(poles, N);
n = 0;

if nargin == 5
    for m = 1:N
        n = n + 1;
        if ind_poles(m) == 1
            a = poles(n); a1=real(a); a2=imag(a);
            c = C(n,:); c1=real(c); c2=imag(c);
            b = B(n,:); b1=2*real(b); b2=-2*imag(b);
            Ablock = [a1 a2;-a2 a1];
            A(n:n+1,n:n+1) = Ablock;
            C(n,:) = c1;
            C(n+1,:) = c2;
            B(n,:) = b1;
            B(n+1,:) = b2;
        end
    end 
    C = C.';
else
    for m = 1:N
        n = n + 1;
        if ind_poles(m) == 1
            a = poles(n); a1=real(a); a2=imag(a);
            b = B(n,:); b1=2*real(b); b2=-2*imag(b);
            Ablock = [a1 a2;-a2 a1];
            A(n:n+1,n:n+1) = Ablock;
            B(n,:) = b1;
            B(n+1,:) = b2;
        end
    end
    C = 0;
end