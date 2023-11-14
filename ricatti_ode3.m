function dy = ricatti_ode3(t,T,y,x,lambda)
% ODE RICATTI classical case
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za
kappa = x(10:12);
theta = x(13:15);

Nn = length(lambda);

dy = zeros(4*Nn,1);
for j=0:Nn-1
    for k=1:3
        dy(4*j+k) = a(t,k,x)*y(4*j+k)^2 + b(t,T,k,x,lambda(j+1))*y(4*j+k) + c(t,T,k,x,lambda(j+1));
        dy(4*j+4) = dy(4*j+4) + kappa(k)*theta(k)*y(4*j+k);
    end
end
end

function ar = a(t,i,x)
sigma = x(16:18);
ar = 0.5*sigma(i)^2;
end

function br = b(t,T,i,x,lambda)


kappa = x(10:12);
kappa_0 = x(1:3);
kappa_1 = x(4:6);
eta = x(7:9);
sigma = x(16:18);
rho = x(19:21);
phi = (kappa_0(i)+kappa_1(i)*(T-t))*exp(-eta(i)*(T-t));
br = 1i*sigma(i)*lambda*rho(i)*phi - kappa(i);

end

function cr = c(t,T,i,x,lambda)

kappa_0 = x(1:3);
kappa_1 = x(4:6);
eta = x(7:9);


phi = (kappa_0(i)+ kappa_1(i)*(T-t))*exp(-eta(i)*(T-t));
cr = -0.5*(lambda.^2+1i*lambda).*phi.^2;


end
