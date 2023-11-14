function [phi, delta] = phi_omega3(omega,x,T0,T,Ft,vt)
% Solving ODE RICATTI classical case by ODE45
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za


Nn = length(omega);
t_values=linspace(0,T0,365*T0); % the time vector you want to use, or use tspan type vector, [0 10]
initial_cond=zeros(4*Nn,1); %[0 ; 0 ; 0; 0];  
 [tv,Yv]=ode45(@(t,y) ricatti_ode3(t,T,y,x,omega),t_values,initial_cond);
 
 temp = reshape(Yv(end,:),4,Nn).';


 phi = exp(temp(:,4)+temp(:,1:3)*vt+1i.*omega.*log(Ft));

 delta = temp(:,1:3).*repmat(phi,1,3);

end
