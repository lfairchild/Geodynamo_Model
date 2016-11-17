function dTdr = Derivative(r,T)
%
% returns the derivative of an arbitrarily spaced T field
%
%

 n = length(r);

 h1 = r(2) - r(1);
 h2 = r(3) - r(2);
 tmp = h1 * h2 * (h1 + h2);
 dTdr(1) = (-(h2^2 + 2*h1*h2)* T(1) + (h1+h2)^2*T(2) - h1^2*T(3))/tmp;

% evaluate temperature gradient
for k = 2 : n-1
    h1 = r(k) - r(k-1);
    h2 = r(k+1) - r(k);
    tmp = h1 * h2 * (h1 + h2);
    dTdr(k)= (-h2^2*T(k-1) + (h2^2-h1^2)*T(k) + h1^2*T(k+1))/tmp;

end

h1 = r(n-1) - r(n-2);
h2 = r(n) - r(n-1);
tmp = h1 * h2 * (h1 + h2);
dTdr(n) = ( (h1^2+2*h1*h2)*T(n) - (h1+h2)^2*T(n-1) + h2^2*T(n-2))/tmp;



end

    



 
 