function [PvectorIS, dp] = DP(DerStr,PT,Po,T,Pmed,Xmed,dvector,omegas, a,b)
%UNTITLED3 Summary of this function goes here
% this code is updated on March 13, to calculate both the P_j and its
% derivatives with respect to x_i, and dummy patient is not considered 
%   Detailed explanation goes here
%this code is further updated on october 10 to use quadratic function.
npatients=length(dvector);

for j = 1: npatients-1
    r(j) = sum(dvector(1:j)); %x_1 is the arrival time of patient 1
end
dp = zeros(npatients-1, npatients-1);%p_i wrt x_0,x_1,...x_{n-1}
    
if DerStr==1 %Linear
    for j =1: npatients-1 %for patient 1 to 12
        PvectorIS(j,1)=((PT-Po)/T)*r(j)+Po;
        for k = 1:j % if j>k, then the derivative is zero
            dp(j,k) =( 2*omegas(j)-1)*(PT-Po)/T;
        end
    end
end
 if DerStr==2 %quadratic- updated on 2019-10-11
      c1 = 4*(Po-PT)/T^2;
      c2 = -c1*T;
      c3 = Po;
     
  % to be ocontiued...      
     for j = 1 : npatients-1
      
            %if r(j)<Xmed
                
                PvectorIS(j,1)=c1*r(j)^2 + c2*r(j) + c3; 
                for k =1: j
                    dp(j,k) = (2*omegas(j)-1)*(2*c1*r(j)+c2);
                end

                %else %sum(dvectorIS(1:i-1))>=Xmed

            
          
            
     end
 end
        
 if DerStr==3 %Specific Function
        %a+b defines the peak of the function, and b-a its lowest point
        a=(Po-PT)/2;
        b=(Po+PT)/2;
        
        
        for j = 1 : npatients-1
            PvectorIS(j,1)=a*cos(4*pi*r(j)/T)+b;
            for k =1: j
                    dp(j,k) =( 1- 2*omegas(j))*a*(4*pi/10)*sin(4*pi*r(j)/T);
            end
            
        end
 end
 PvectorIS = omegas.*PvectorIS + (1-omegas).*(1-PvectorIS);
 end
  
 
    





