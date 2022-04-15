function [PvectorIS] = P(DerStr,PT,Po,T,Pmed,Xmed,dvectorIS,a,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
npatients=length(dvectorIS);
for i=1:npatients
    
    if DerStr==1 %Linear
        
        if i==1
            PvectorIS(i,1)=Po;
        else
            PvectorIS(i,1)=((PT-Po)/T)*sum(dvectorIS(1:i-1))+Po;
        end
        
    elseif DerStr==2 %quardratic
        c1 = 4*(Po-PT)/T^2;
        c2 = -c1*T;
        c3 = Po;
        
        if i==1
            PvectorIS(i,1)=Po;
        else
            x=   sum(dvectorIS(1:i-1)); % arrival time of patient i 
            PvectorIS(i,1)= c1*x^2 + c2*x +c3;
            
            
        end
        
    elseif DerStr==3 %Specific Function
        %a+b defines the peak of the function, and b-a its lowest point
        a=(Po-PT)/2;
        b=(Po+PT)/2;
        
        if i==1
            PvectorIS(i,1)=Po;
        else
            PvectorIS(i,1)=a*cos(4*pi*sum(dvectorIS(1:i-1))/T)+b;
        end
    end
    
end


end

