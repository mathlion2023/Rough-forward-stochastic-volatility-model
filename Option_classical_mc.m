function op_mc = Option_classical_mc(x,t,T0,T,F0,strike, NSim, NTime)
% Monte Carlo price for the classical case
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za

% Model parameters
kappa_0 = x(1:3);
kappa_1 = x(4:6);
eta = x(7:9);
kappa = x(10:12);
theta = x(13:15);
sigma = x(16:18);
rho = x(19:21);
V0=x(22:24);


dt = (T0-t)/NTime;
dw1 = cell(NSim,1);
dw2 = cell(NSim,1);
i0 = 0;

op_mc = 0.0;

parfor j=1:NSim
    Vt1 = V0;
    Ft1 = F0;
    dw1{j} = randn(NTime,3).*sqrt(dt);
    dw2{j} = randn(NTime,3).*sqrt(dt);
    %simulate variance process
    for i=1:NTime
        sig1 = sum((kappa_0 + kappa_1.*(T-t-i*dt)).^2.*exp(-2.*eta.*(T-t-i*dt)).*Vt1);
        sig2 = sum((kappa_0 + kappa_1.*(T-t-i*dt)).*exp(-eta.*(T-t-i*dt)).*sqrt(Vt1).*(rho.*dw1{j}(i,:)+sqrt(1-rho.^2).*dw2{j}(i,:)));
        Ft2 = Ft1.*exp(-0.5.*sig1*dt+sig2);
        Ft1 = Ft2;
        Vt2 = Vt1 + kappa.*(theta - Vt1).*dt+sigma.*sqrt(Vt1).*dw1{j}(i,:);
        Vt1 = Vt2;
        for iv=1:3
            if(Vt1(iv)<0)
                Vt1(iv) = 0.0;
                i0 = i0 + 1;
            end
        end
    end
    op_mc = op_mc + max(Ft1 - strike,0.0);
end
op_mc = op_mc/NSim;
end