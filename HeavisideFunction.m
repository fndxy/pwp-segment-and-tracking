% function H = Heaviside(phi,e)
% if x>e
%     H = 1.0;
% elseif x<-1.0*e
%     H = 0.0;
% else
%     H = 0.5*(1+2/pi*arctan(phi/e));
% end





function f = HeavisideFunction(phi,sigma)
[m n] = size(phi);
f = zeros(m,n);
for i = 1:m
    for j = 1:n
        temp = phi(i,j);
        if temp > sigma
            f(i,j) = 1;
        elseif temp < -1*sigma
            f(i,j) = 0;
        else
            f(i,j) = 0.5*(1+phi(i,j)/sigma+1/pi*sin(pi*phi(i,j)/sigma));
        end
    end
end


