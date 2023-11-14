function [call, put]=Semi_Analytic3(x,alpha,T0,T,F,K,Nn,r)
% Following Duffie, Pan and Singleton (2000) pricing method
% Rough Case
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za
nk = length(K);
[u,w]=fclencurt(Nn,10e-7,70);

phi_1 = CF_RFM3(x,u, alpha,T0,T,F);
phi_2 = CF_RFM3(x,-1i+u, alpha,T0, T,F);
phi_0 = CF_RFM3(x,0, alpha,T0,T,F);
phi_i = CF_RFM3(x,-1i, alpha,T0,T,F);


phi_all = repmat(phi_1,nk,1);
phi_all2 = repmat(phi_2,nk,1);
K_all = reshape(repmat(K,Nn,1),Nn*nk,1);

x_all = repmat(u,nk,1);
w_all = repmat(w,nk,1);

I11 = sum(reshape(w_all.*(imag(phi_all.*exp(-1i.*x_all.*log(K_all)))./x_all),Nn,nk));
I12 = sum(reshape(w_all.*(imag(phi_all2.*exp(-1i.*x_all.*log(K_all)))./x_all),Nn,nk));

disc = exp(-r*T); % T0
temp = K.*(0.5*phi_0 - I11/pi) - (0.5*phi_i-I12/pi);
 put=disc*temp;
 call = put + disc*(F - K);
call=abs(call);

end