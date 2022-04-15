function [OmegaAProb] = P_A_function(OmegaAProb,DerStr,PT,Po,T,Pmed,Xmed,dvector,a,b,n_k,OmegaA)
npatients=length(dvector);
c1 = 4*(Po-PT)/T^2;
c2 = -c1*T ;
c3 = Po;
for i=1:npatients
    for j=1:n_k
        
        if DerStr==1 %Linear
            
            if OmegaA(i,j)==1 && i>=2
                OmegaAProb(i,j)=((PT-Po)/T)*sum(dvector(1:i-1))+Po;
                
            elseif OmegaA(i,j)==0 && i>=2
                OmegaAProb(i,j)=1-((PT-Po)/T)*sum(dvector(1:i-1))-Po;
            else
                
                if OmegaA(i,j)==1 && i<2
                    OmegaAProb(i,j)=Po;
                else
                    OmegaAProb(i,j)=1-Po;
                    
                end
            end
            
            
        elseif DerStr==2 %quadratic
        
            if i==1
                p=Po;
            else

                x =   sum(dvector(1:i-1)); % arrival time of patient i 

                p= c1 * x^2 + c2*x +c3;

            end
         OmegaAProb(i,j)=p * OmegaA(i,j) + (1-p)*(1-OmegaA(i,j));
            
           
            
        elseif DerStr==3 %Specific Function
            
            %a+b defines the peak of the function, and b-a its lowest point
            a=(Po-PT)/2;
            b=(Po+PT)/2;
            
            if OmegaA(i,j)==1 && i>=2
                
                OmegaAProb(i,j)=a*cos(4*pi*sum(dvector(1:i-1))/T)+b;
                
            elseif OmegaA(i,j)==0 && i>=2
                OmegaAProb(i,j)=1-(a*cos(4*pi*sum(dvector(1:i-1))/T)+b);
            else
                
                if OmegaA(i,j)==1 && i<2
                    OmegaAProb(i,j)=Po;
                else
                    OmegaAProb(i,j)=1-Po;
                    
                end
            end
            
            
        end
    end
    
end


end

