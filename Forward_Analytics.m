function [call, put]=Forward_Analytics(x,T0,T,F,K,Nn,r)

% Following Duffie, Pan and Singleton (2000) pricing method
% Classical Case
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za
v=x(22:24)';
[x1,w]=fclencurt(Nn,10e-7,100);


[phi1,del1] = phi_omega3(x1,x,T0,T,F,v);

nk = length(K);
phi_all = repmat(phi1,nk,1);

K_all = reshape(repmat(K,Nn,1),Nn*nk,1);

x_all = repmat(x1,nk,1);
w_all = repmat(w,nk,1);

I11 = sum(reshape(w_all.*(imag(phi_all.*exp(-1i.*x_all.*log(K_all)))./x_all),Nn,nk));


[phi2,del2] = phi_omega3(-1i+x1,x,T0,T,F,v);
phi_all2 = repmat(phi2,nk,1);
I12 = sum(reshape(w_all.*(imag(phi_all2.*exp(-1i.*x_all.*log(K_all)))./x_all),Nn,nk));

disc = exp(-r*T); 
temp = K.*(0.5*phi_omega3(0,x,T0, T,F,v) - I11/pi) - (0.5*phi_omega3(-1i,x,T0,T,F,v)-I12/pi);

 put=abs(disc*temp);
 call = put + disc*(F - K);
 call=abs(call);
end