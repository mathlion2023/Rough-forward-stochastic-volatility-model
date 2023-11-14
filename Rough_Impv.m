function impv = Rough_Impv(x,id,TTM,futures, options,para_yield)

% Implied vol for the rough case
% Copyright
% Mesias Alfeus 2023
% Department of Statistics and Actuarial Science
% Stellenbosch University
% email: mesias@sun.ac.za

Tr = [30, 90, 180, 360];
ik = 1;Nn = 40;
Tt = zeros(4,2);
l = zeros(4,1);
va = zeros(4,1);
lo = zeros(4,1);

impv = zeros(1,4); % model
TF=2;
   
    for i=1:4
        l(i) = sum(TTM{id,1}<=Tr(i));
    end
    if(l(4)<length(TTM{id,1}))
        if(l(1)>0)
            for i=1:4
                Tt(i,1) = TTM{id,1}(1,l(i))/360;
                Tt(i,2) = TTM{id,1}(1,l(i)+1)/360;
                [va(i),lo(i)] = min(abs(options{id,l(i)}(:,2) - futures{id,1}(l(i))));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,1)))/Tt(i,1);
                
                
               call1 = Semi_Analytic3(x,x(25),Tt(i,1),TF,futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),Nn,abs(ra));
                ivm1 = blsimpv(futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),abs(ra),Tt(i,1),call1);
                
                [va(i),lo(i)] = min(abs(options{id,l(i)+1}(:,2) - futures{id,1}(l(i)+1)));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,2)))/Tt(i,2);
               
                
                call2 = Semi_Analytic3(x,x(25),Tt(i,2),TF,futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),Nn,abs(ra));
                ivm2 = blsimpv(futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),abs(ra),Tt(i,2),call2);
                
                
                if(isnan(ivm1)==0 && isnan(ivm2) ==0)
                    
                    impv(1,i) = ivm1+(ivm2-ivm1)*(Tr(i)/360 - Tt(i,1))/(Tt(i,2) - Tt(i,1));
                elseif(isnan(ivm1)==0 && isnan(ivm2) ==1)
                   
                    impv(1,i) = ivm1;
                elseif(isnan(ivm2)==0 && isnan(ivm1) ==1)
                 
                    impv(1,i) = ivm2;
                end
            end
        else
           for i=1:4
               if(i>1)
                Tt(i,1) = TTM{id,1}(1,l(i))/360;
                Tt(i,2) = TTM{id,1}(1,l(i)+1)/360;
                [va(i),lo(i)] = min(abs(options{id,l(i)}(:,2) - futures{id,1}(l(i))));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,1)))/Tt(i,1);
                
                call1 = Semi_Analytic3(x,x(25),Tt(i,1),TF,futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),Nn,abs(ra));
                ivm1 = blsimpv(futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),abs(ra),Tt(i,1),call1);
                
                
                [va(i),lo(i)] = min(abs(options{id,l(i)+1}(:,2) - futures{id,1}(l(i)+1)));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,2)))/Tt(i,2);
                
                
                call2 = Semi_Analytic3(x,x(25),Tt(i,2),TF,futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),Nn,abs(ra));
                ivm2 = blsimpv(futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),abs(ra),Tt(i,2),call2);
                
                if(isnan(ivm1)==0 && isnan(ivm2) ==0)
                    
                    impv(1,i) = ivm1+(ivm2-ivm1)*(Tr(i)/360 - Tt(i,1))/(Tt(i,2) - Tt(i,1));
                elseif(isnan(ivm1)==0 && isnan(ivm2) ==1)
                  
                    impv(1,i) = ivm1;
                elseif(isnan(ivm2)==0 && isnan(ivm1) ==1)
                    
                      impv(1,i) = ivm2;
                end
               else
                   Tt(i,2) = TTM{id,1}(1,l(i)+1)/360;
                [va(i),lo(i)] = min(abs(options{id,l(i)+1}(:,2) - futures{id,1}(l(i)+1)));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,2)))/Tt(i,2);
                
                
                 call2 = Semi_Analytic3(x,x(25),Tt(i,2),TF,futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),Nn,abs(ra));
                ivm2 = blsimpv(futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),abs(ra),Tt(i,2),call2);
            
                impv(1,i) = ivm2;
            end
           end 
        end
    else
        if(l(1)>0)
        for i=1:4
            if(i==4)
                Tt(i,1) = TTM{id,1}(1,l(i))/360;
                [va(i),lo(i)] = min(abs(options{id,l(i)}(:,2) - futures{id,1}(l(i))));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,1)))/Tt(i,1);
               
                call1 = Semi_Analytic3(x,x(25),Tt(i,1),1,futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),Nn,abs(ra));
                impv(1,i) = blsimpv(futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),abs(ra),Tt(i,1),call1);
                
            else
                Tt(i,1) = TTM{id,1}(1,l(i))/360;
                Tt(i,2) = TTM{id,1}(1,l(i)+1)/360;
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,1)))/Tt(i,1);
                [va(i),lo(i)] = min(abs(options{id,l(i)}(:,2) - futures{id,1}(l(i))));
                
         
                call1 = Semi_Analytic3(x,x(25),Tt(i,1),1,futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),Nn,abs(ra));
                ivm1 = blsimpv(futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),abs(ra),Tt(i,1),call1);
                
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,2)))/Tt(i,2);
                [va(i),lo(i)] = min(abs(options{id,l(i)+1}(:,2) - futures{id,1}(l(i)+1)));
                
                call2 = Semi_Analytic3(x,x(25),Tt(i,2),1,futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),Nn,abs(ra));
                ivm2 = blsimpv(futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),abs(ra),Tt(i,2),call2);
                
               if(isnan(ivm1)==0 && isnan(ivm2) ==0)
                   
                    impv(1,i) = ivm1+(ivm2-ivm1)*(Tr(i)/360 - Tt(i,1))/(Tt(i,2) - Tt(i,1));
               elseif(isnan(ivm1)==0 && isnan(ivm2) ==1)
                   
                   impv(1,i) = ivm1;
               elseif(isnan(ivm2)==0 && isnan(ivm1) ==1)
                  
                   impv(1,i) = ivm2;
               end
            end
        end
        else
            for i=1:4
            if(i==4)
                Tt(i,1) = TTM{id,1}(1,l(i))/360;
                [va(i),lo(i)] = min(abs(options{id,l(i)}(:,2) - futures{id,1}(l(i))));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,1)))/Tt(i,1);
                
                call1 = Semi_Analytic3(x,x(25),Tt(i,1),1,futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),Nn,abs(ra));
                impv(1,i) = blsimpv(futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),abs(ra),Tt(i,1),call1);
                
            elseif(i==1)
                Tt(i,2) = TTM{id,1}(1,l(i)+1)/360;
                [va(i),lo(i)] = min(abs(options{id,l(i)+1}(:,2) - futures{id,1}(l(i)+1)));
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,2)))/Tt(i,2);
                
                 call2 = Semi_Analytic3(x,x(25),Tt(i,2),1,futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),Nn,abs(ra));
                impv(1,i) = blsimpv(futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),abs(ra),Tt(i,2),call2);
                
            else
                Tt(i,1) = TTM{id,1}(1,l(i))/360;
                Tt(i,2) = TTM{id,1}(1,l(i)+1)/360;
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,1)))/Tt(i,1);
                [va(i),lo(i)] = min(abs(options{id,l(i)}(:,2) - futures{id,1}(l(i))));
              
                call1 = Semi_Analytic3(x,x(25),Tt(i,1),1,futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),Nn,abs(ra));
                ivm1 = blsimpv(futures{id,1}(l(i)),options{id,l(i)}(lo(i),2),abs(ra),Tt(i,1),call1);
                
                ra = -log(zero_coupon(para_yield(:,id-5000),0,Tt(i,2)))/Tt(i,2);
                [va(i),lo(i)] = min(abs(options{id,l(i)+1}(:,2) - futures{id,1}(l(i)+1)));
                
                call2 = Semi_Analytic3(x,x(25),Tt(i,2),1,futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),Nn,abs(ra));
                ivm2 = blsimpv(futures{id,1}(l(i)+1),options{id,l(i)+1}(lo(i),2),abs(ra),Tt(i,2),call2);
                
               if(isnan(ivm1)==0 && isnan(ivm2) ==0)
                  
                    impv(1,i) = ivm1+(ivm2-ivm1)*(Tr(i)/360 - Tt(i,1))/(Tt(i,2) - Tt(i,1));
               elseif(isnan(ivm1)==0 && isnan(ivm2) ==1)
                  
                   impv(1,i) = ivm1;
               elseif(isnan(ivm2)==0 && isnan(ivm1) ==1)
                   
                   impv(1,i) = ivm2;
               end
            end
        end
        end
    end
    
    

end
