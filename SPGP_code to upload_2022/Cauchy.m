function [Xk1,ProjGrad,Norm_ProjGrad, norm_sk] = Cauchy(Estim,gradient,Hessian,Deltak,State,Xk,T,eps)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters coming from STRONG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%State = 1;  %TEMPORARY
if State==1
    
    PolyApprox = @(x) Estim+transpose(gradient)*(x-Xk);
    
elseif State==2
    
    PolyApprox = @(x) Estim+transpose(gradient)*(x-Xk)+0.5*transpose(x-Xk)*Hessian*(x-Xk);
    
end

%Hessian

%SecOrderA = @(x) Estim+transpose(gradient)*(x-Xk)+0.5*transpose(x-Xk)*Hessian*(x-Xk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Precode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmin=0;
tmax_upper=999999;
tmax=tmax_upper;
k_epp=0.125;
k_lbs=0.95;
k_ubs=0.55;
k_frd=0.75;

stop=0;
j=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t0=Deltak/norm(gradient);

while stop==0
    
    if j==1
        tj=t0;
    else
        tj=tj_plus1;
    end
    
    Xk1_temporal=Xk-tj*gradient;
    
    [Xk1_temporal_P] = ProjFeasible(Xk1_temporal,T,eps);
    
    sk_tj=Xk1_temporal_P-Xk;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Checkig for stopping conditions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    norm_sk=norm(sk_tj);
    dotprod_gradsk=dot(gradient,sk_tj);
    
%    [ProjTC_TrasGrad] = ProjTangentCone(Xk1_temporal_P,(-gradient),T); %gives the point, not the vector
   ProjGrad=Xk1_temporal_P;  %%% This is just a trash value     ProjTC_TrasGrad-Xk1_temporal_P;
%    Norm_ProjGrad=norm(ProjGrad);

    %NOTE: The calculation below replaces the criterion based on the
    %projection of gradient onto the tangent cone.
    %The justification is in Thm 12.1.5(iii) of the Conn, Gould & Toint
    %book - the norm of projected gradient is replaced by a lower bound
    theta = norm_sk+1; %Any theta > norm_sk works
    [d,CM] = CriticalityMeasure(Xk,gradient,theta,T,eps);
    Norm_ProjGrad = (CM-abs(dotprod_gradsk))/(2*theta);
    
    C1=Deltak-norm_sk;
    C2=PolyApprox(Xk)+k_ubs*dotprod_gradsk-PolyApprox(Xk1_temporal_P);
    C3=norm_sk-k_frd*Deltak;
    C4=PolyApprox(Xk1_temporal_P)-PolyApprox(Xk)-k_lbs*dotprod_gradsk;
    C5=((k_epp*abs(dotprod_gradsk))*(1/Deltak))-Norm_ProjGrad;
    
    C1G=[C1;C2];
    C2G=[C3;C4;C5];
    
    
    if size(find(C1G<0),1)>=1
        
        tmax=tj;
        
        tj_plus1=0.5*(tmin+tmax);
        
    elseif size(find(C2G<0),1)==3

        tmin=tj;
        
        if tmax==tmax_upper
            tj_plus1=2*tj;
        else
            tj_plus1=0.5*(tmin+tmax);
            
        end
        
    else
        stop=1;
        Xk1=Xk1_temporal_P;
    end
    j=j+1;
end



 %%%%%%%%






end
