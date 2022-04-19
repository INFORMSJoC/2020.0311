function [Xk_Plus1,Deltak_Plus1,Innerk,Ro_k] = RCnSR_Tests(rejectHo,Estim,Estim_kOp,AproxEstim_Xk,AproxEstim_Xk_Op,Gamma1,Gamma2,Deltak,Xk,Xk_Op,Nu0,Nu1,State)
disp('Ratio Comparison Test Statistic')
Ro_k=(Estim-Estim_kOp)/(AproxEstim_Xk-AproxEstim_Xk_Op)

if State==1 % \delta is big enough
    
    if Ro_k<Nu0 % if fails the RC test 
        Deltak_Plus1=Gamma1*Deltak; % shrink the region
        Xk_Plus1=Xk;                % Don't update X_k
        Innerk=0;                   % stay in the outer loop
        
    elseif Nu0<=Ro_k && Ro_k<Nu1    % RC test is fair
        
        if rejectHo==1 %Passes the SR test
            
            Xk_Plus1=Xk_Op; %update x_k
            Deltak_Plus1=Deltak;% keep the trust region unchanged 
            Innerk=0;
        else          %  SR test is not passed
            Xk_Plus1=Xk;
            Deltak_Plus1=Gamma1*Deltak; %shrink the trust region
            Innerk=0;
        end
        
    else %Ro>Nu1 (Solution proposed is apparently better better...)
        
        if rejectHo==1 %Passes the test
            Xk_Plus1=Xk_Op;
            Deltak_Plus1=Gamma2*Deltak;%enlarge the  trust region
            Innerk=0;
            
        else
            Xk_Plus1=Xk;
            Deltak_Plus1=Gamma1*Deltak;
            Innerk=0;
        end
        
    end
    
else %State=2
    
    if Ro_k<Nu0 %Go to INNER LOOP
        Deltak_Plus1=Deltak;
        Xk_Plus1=Xk;
        Innerk=1;
        
    elseif Nu0<=Ro_k && Ro_k<Nu1
        
        if rejectHo==1 %Passes the test
            
            Xk_Plus1=Xk_Op;
            Deltak_Plus1=Deltak;
            Innerk=0;
        else
            Xk_Plus1=Xk;
            Deltak_Plus1=Deltak;
            Innerk=1;
        end
        
    else %Ro>Nu1 (Solution proposed is apparently better better...)
        
        if rejectHo==1 %Passes the test
            Xk_Plus1=Xk_Op;
            Deltak_Plus1=Gamma2*Deltak;
            Innerk=0;
            
        else
            Xk_Plus1=Xk;
            Deltak_Plus1=Deltak;
            Innerk=1;
        end
        
    end
    
    
    
end