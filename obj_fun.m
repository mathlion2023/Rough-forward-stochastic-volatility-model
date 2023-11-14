function Diff = obj_fun(x,id,TTM,futures, options,mod,Nn)
% Calibration Objective function
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za
diff = 0.0;

     F = futures{id};
     NI=0;
     TToM = TTM{id};% options time to maturity
     NT = min(109,length(TToM));
     if(id<4500)
        nop = 1;
    else
        nop = 7;
     end
    
for j=2:nop:NT 
         K = options{id,j}(:,2)';
         % filtering
        ind_c = find(K/F(j)>=1 & K/F(j)<1.2); % for call
        ind_p = find(K/F(j)<=1 & K/F(j)>0.8);
   if(isempty(ind_c) == 0 && isempty(ind_p) == 0)
            T = TToM(j)/365;
            [mdiff,indk] = max(options{id,j}(:,3) - options{id,j}(:,1) + F(j));
            r = abs(-log(mdiff/K(indk))/T);
            if mod==1
                alpha=x(25);
                  [call, put] = Semi_Analytic3(x,alpha,T,10,F(j),K,Nn,r);
            else
                 [call, put]=Forward_Analytics(x,T,10,F(j),K,Nn,r);
            end
   NI=NI+numel(options{id,j}(ind_c,1))+numel(options{id,j}(ind_p,3));

        diff = diff +mean((call(ind_c)' - options{id,j}(ind_c,1)).^2)+mean((put(ind_p)'-options{id,j}(ind_p,3)).^2); % Calibration to both calls and puts
   end
 
Diff = diff/NI;
end
