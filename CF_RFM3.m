function [cf] = CF_RFM3(x,u, alpha,T0,T,F)
% Characteristic function of Rough Forward Model
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za

kappa = x(10:12);

theta=x(13:15);

V0 = x(22:24)';



Nn=length(u);
 [ h ] = PC_Schemes3(x,u,alpha,T0 ,T );
 temp = reshape(h(:,end),3,Nn).';
 temp1=zeros(Nn,1);
 parfor i=0:Nn-1
     for k=1:3
         temp1(i+1)=temp1(i+1)+kappa(k)*theta(k)*h(3*i+k,end);
     end
 end

y = exp(temp1+temp*V0+1i.*u.*log(F));
cf=y;
end