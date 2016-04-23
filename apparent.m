function f = apparent(Img,phi,sigma,epislon)
% phi = NeumannBoundCond(phi);
H = Heavisidef(phi,sigma);
D=Dirac(phi,sigma); %Dirac函数

K = curvature(phi); % 曲率

[m,n] = size(phi);
L = Laplacian(phi);

nf = sum(H(:));
Neg_H = 1-H;
nb = sum(Neg_H(:));

[P_yi_Mf ,P_yi_Mb] = condition(Img,phi);
Pf = P_yi_Mf./(nf*P_yi_Mf + nb*P_yi_Mb+1e-10);
Pb = P_yi_Mb./(nf*P_yi_Mf + nb*P_yi_Mb+1e-10);
Px = H.*Pf + Neg_H.*Pb;


phi_1 = D.*(Pf-Pb)./(Px);
phi_2 = 1/epislon^2*(L-1* K);


f = phi + 0.05*(phi_1 - phi) ;







function [hist_f hist_b] = histForeBack(Img,phi)  %求 前景、背景 直方图
[m,n] = size(phi);
hist_f = zeros(m,n);
hist_b = zeros(m,n);
for i = 1:m
    for j = 1:n
        temp = Img(i,j);
        if phi(i,j)>0 %前景
            hist_f(temp+1) = hist_f(temp+1) + 1;
        elseif(phi(i,j)<0) %背景
            hist_b(temp+1) = hist_b(temp+1) + 1;
        end
    end
end

function [P_yi_Mf, P_yi_Mb] = condition(Img,phi)
[m,n] = size(phi);
P_yi_Mf = zeros(m,n);
P_yi_Mb = zeros(m,n);
[hist_f, hist_b] = histForeBack(Img,phi);
sum_pos = sum(hist_f(:));
sum_neg = sum(hist_b(:));
% for i = 1:m
%     for j = 1:n
%         temp = Img(i,j);
%         P_yi_Mf(i,j) = hist_f(temp+1)/sum(hist_f(:));
%         P_yi_Mb(i,j) = hist_b(temp+1)/sum(hist_b(:));
%     end
% end
P_yi_Mf = hist_f(Img+1)./sum_pos;
P_yi_Mb = hist_b(Img+1)./sum_neg;



function g = NeumannBoundCond(f)
% Make a function satisfy Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  

function H = Heavisidef(phi,sigma)
[m,n] = size(phi);
H = zeros(m,n);

for i =1:m
    for j = 1:n
        H(i,j) = 0.5*(1+2/pi*atan(phi(i,j)/sigma));
    end
end
function D = Dirac(phi,sigma)
[m,n] = size(phi);
D = zeros(m,n);

for i =1:m
    for j = 1:n
        D(i,j) = 1/pi*(sigma/(sigma^2+phi(i,j)^2));
    end
end    


function L = Laplacian(phi)

    
phi = double(phi);
[m,n] = size(phi);
P = padarray(phi,[1,1],1,'pre');
P = padarray(P,[1,1],1,'post');

% central difference
fy = P(3:end,2:n+1)-P(1:m,2:n+1);
fx = P(2:m+1,3:end)-P(2:m+1,1:n);
fyy = P(3:end,2:n+1)+P(1:m,2:n+1)-2*phi;
fxx = P(2:m+1,3:end)+P(2:m+1,1:n)-2*phi;
fxy = 0.25.*(P(3:end,3:end)-P(1:m,3:end)+P(3:end,1:n)-P(1:m,1:n));
G = (fx.^2+fy.^2).^(0.5);
L = fxx + fyy;



