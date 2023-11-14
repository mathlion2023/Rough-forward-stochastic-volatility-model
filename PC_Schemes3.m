function [ h ] = PC_Schemes3(x,a,alpha,T0 ,T )
% Solving rough ODE RICATTI classical case by Adam's (predictor-corrector) method
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za
kappa_0 = x(1:3);
kappa_1 = x(4:6);
eta = x(7:9);
kappa = x(10:12);
theta = x(13:15);
sigma = x(16:18);
rho = x(19:21);
V0=x(22:24);

t=linspace(0,T0,365*T0);
an=length(a);
N=length(t);

h = zeros(3*an,N+1);
dt=t(2);
C = (dt^alpha/gamma(alpha+1));


 k=1:N;
 bk = k.^alpha - (k-1).^alpha;
 ak = (k+1).^(alpha+1) - 2*k.^(alpha+1) + (k-1).^(alpha+1);
 
for i=0:an-1 
for j = 1:N
 for k=1:3
     phi=(kappa_0(k)+kappa_1(k).*(T-t(1:j))).*exp(-eta(k).*(T-t(1:j)));
    
    FCF = @(h,psi) 0.5*(-a(i+1).^2 - 1i*a(i+1)).*psi.^2 + (1i.*a(i+1).*rho(k).*sigma(k).*psi-kappa(k)).*h + 0.5*h.^2*(sigma(k))^2;
    p = C.*sum(bk(j:-1:1).*FCF(h(3*i+k,1:j),phi));
    h(3*i+k,j+1) = (dt^alpha/gamma(alpha+2))*(FCF(p,1) + ((j-1)^(alpha+1) - (j-1-alpha)*j^alpha)*FCF(h(3*i+k,1),phi(1)) + sum(ak((j-1):-1:1).*FCF(h(3*i+k,2:j),phi(2:j))));
end
end

end


end
