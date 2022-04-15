function [ DerivateAProbReg ] = P_Der( DerStr,OmegaAReg,DerivateAProbReg,PT,Po,T,Pmed,Xmed,dvector,a,b,n_k)
npatients=length(dvector);

if DerStr==1 %Linear
    for k=1:(npatients-1) %Number of possible derivatives
        for i=1:npatients
            for j=1:n_k
                
                if OmegaAReg(i,j)==1 && i>=2 && k<=(i-1)
                    DerivateAProbReg(i,j,k)=((PT-Po)/T);
                    
                elseif OmegaAReg(i,j)==0 && i>=2 && k<=(i-1)
                    DerivateAProbReg(i,j,k)=-((PT-Po)/T);
                else
                    
                    DerivateAProbReg(i,j,k)=0;
                    
                end
            end
        end
    end

elseif DerStr==2 %Piecewise
    
    for k=1:(npatients-1) %Number of possible derivatives
        for i=1:npatients
            for j=1:n_k
                
                if OmegaAReg(i,j)==1 && i>=2 && k<=(i-1) && sum(dvector(1:i-1))<Xmed
                    DerivateAProbReg(i,j,k)=((Pmed-Po)/Xmed);
                    
                elseif OmegaAReg(i,j)==0 && i>=2 && k<=(i-1) && sum(dvector(1:i-1))<Xmed
                    DerivateAProbReg(i,j,k)=-((Pmed-Po)/Xmed);
                else
                    if OmegaAReg(i,j)==1 && i>=2 && k<=(i-1) && sum(dvector(1:i-1))>=Xmed
                        
                        DerivateAProbReg(i,j,k)=((PT-Pmed)/(T-Xmed));
                        
                    elseif OmegaAReg(i,j)==0 && i>=2 && k<=(i-1) && sum(dvector(1:i-1))>=Xmed
                        
                        DerivateAProbReg(i,j,k)=-((PT-Pmed)/(T-Xmed));
                        
                    else
                    
                    DerivateAProbReg(i,j,k)=0;
                    
                    end
                end
            end
        end
    end
    
elseif DerStr==3 %Specific Function
    %a+b defines the peak of the function, and b-a its lowest point
    a=(Po-PT)/2;
    b=(Po+PT)/2;
    
    for k=1:(npatients-1) %Number of possible derivatives
        for i=1:npatients
            for j=1:n_k
                
                if OmegaAReg(i,j)==1 && i>=2 && k<=(i-1)
                    DerivateAProbReg(i,j,k)=-a*(4*pi/10)*sin(4*pi*sum(dvector(1:i-1))/T);
                    
                elseif OmegaAReg(i,j)==0 && i>=2 && k<=(i-1)
                    DerivateAProbReg(i,j,k)=a*(4*pi/10)*sin(4*pi*sum(dvector(1:i-1))/T);
                else
                    
                    DerivateAProbReg(i,j,k)=0;
                    
                end
            end
        end
    end
    
end
        
        

end

