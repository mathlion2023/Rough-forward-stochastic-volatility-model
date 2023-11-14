function price = option_rough_forward_mc(x,NTime,NSim,t, T0, T,K,F0)
% Monte Carlo price for the rough case
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za

% Parameters
kappa_0 = x(1:3);
kappa_1 = x(4:6);
eta = x(7:9);
kappa = x(10:12);
theta = x(13:15);
sigma = x(16:18);
rho = x(19:21);
V0=x(22:24);
H=x(25);

dt=T0/NTime;
alpha=0.5+H;
tv=linspace(0,T,NTime);
c = (1/gamma(alpha));

 op_mc = 0.0;

parfor j=1:NSim
    Vt1 = V0;
    Ft1 = F0;
    dw1= randn(NTime,3).*sqrt(dt);
    dw2 = randn(NTime,3).*sqrt(dt);
    V = zeros(3, NTime); V(:,1) = V0*ones(3,1);
    %simulate variance process
    for i=2:NTime-1
        sig1=(kappa_0(1)+kappa_1(1)*(T-t-i*dt)).^2*exp(-2.*eta(1)*(T-t-i*dt))*Vt1(1)+(kappa_0(2)+kappa_1(2)*(T-t-i*dt)).^2*exp(-2.*eta(2)*(T-t-i*dt))*Vt1(2)+(kappa_0(3)+kappa_1(3)*(T-t-i*dt)).^2*exp(-2.*eta(3)*(T-t-i*dt))*Vt1(3);
        sig2=(kappa_0(1)+kappa_1(1)*(T-t-i*dt))*exp(-eta(1)*(T-t-i*dt))*sqrt(max(Vt1(1),0))*(rho(1).*dw1(i+1,1)+sqrt(1-rho(1).^2).*dw2(i+1,1))+(kappa_0(2)+kappa_1(2)*(T-t-i*dt))*exp(-eta(2)*(T-t-i*dt))*sqrt(max(Vt1(2),0))*(rho(2).*dw1(i+1,2)+sqrt(1-rho(2).^2).*dw2(i+1,2))+(kappa_0(3)+kappa_1(3)*(T-t-i*dt))*exp(-eta(3)*(T-t-i*dt))*sqrt(max(Vt1(3),0))*(rho(3).*dw1(i+1,3)+sqrt(1-rho(3).^2).*dw2(i+1,3));
        Ft2 = Ft1*exp(-0.5.*sig1*dt+sig2); % Futures
        Ft1 = Ft2;
        V(:,i) = max(0,V0' + c *sum( ((tv(i)-(0:i-2).*dt).^(alpha-1) .* kappa' .*(theta'-V(:,1:i-1)))*dt,2) +sum( ((tv(i)-(0:i-2).*dt).^(alpha-1) .* sigma' .* sqrt(V(:,1:i-1))) .*dw1(1:i-1,:)',2 ));
        Vt1=V(:,i);
      
    end

    op_mc = op_mc + max(Ft1-K,0.0);
end
price = op_mc/NSim;
end