function [Xk_Plus1,Deltak_Plus1,Estim_kplus1,Var_Xkplus1,NextGradient,Hessian_kplus1,VarGradientSim_kplus1,Ro_k,ProjGrad,Norm_ProjGrad,n_k, stopID,NewSeedServ,NewSeedShow]=STRONGITk(npatients,T,Epsilon,STlower,STupper,Po,PT,Cw,Cl,Cs,State,OptimMethod,Nu0,Nu1,Gamma1,Gamma2,n0,alpha_k,Xk,Estim,gradient,Hessian,Deltak,Var_Xk,kH,Var_gradient_Xk,OptimMethodInnerLoop,SimDistr,Pmed,Xmed,a,b,DerStr,k,sets,uVector,CurrSeedServ,CurrSeedShow)
Innerk=0;
n_k=n0;
%Epsilon=0.001;

stopID = 0;
%ProjGrad =
%Norm_ProjGrad = 


% Parameters for Strong

XkAux=Xk;
%max_nk =0.4*10^5 ; %maximum sample size in the inner loop

% Strong Code

x0=Xk; %Punto de partida para la b�squeda

FirstOrderA = @(x) Estim+transpose(gradient)*(x-Xk);
SecOrderA = @(x) Estim+transpose(gradient)*(x-Xk)+0.5*transpose(x-Xk)*Hessian*(x-Xk); %TEMPORARY

%Contrains Information
Apat=ones(1,(npatients-1));
b=T;
Aeq=[];
beq=[];
lb=zeros(npatients-1,1);
ub=[];

%What Method?

if OptimMethod==1
    nonlcon=@(x)circlecon(x,Deltak,Xk);
    if State==1
        [xOptimum] = fmincon(FirstOrderA,x0,Apat,b,Aeq,beq,lb,ub,nonlcon);
    elseif State==2
        [xOptimum] = fmincon(SecOrderA,x0,Apat,b,Aeq,beq,lb,ub,nonlcon);
    end
    Xk_Op=xOptimum;
elseif OptimMethod==2
    [xCauchy,ProjGrad,Norm_ProjGrad, norm_sk]=Cauchy(Estim,gradient,Hessian,Deltak,State,Xk,T,Epsilon);
    Xk_Op=xCauchy;
end
disp('Finished outer loop line search');

%Ratio Comparison Test

if State==1
    AproxEstim_Xk_Op=FirstOrderA(Xk_Op);
    AproxEstim_Xk=FirstOrderA(Xk);
elseif State==2
    AproxEstim_Xk_Op=SecOrderA(Xk_Op);
    AproxEstim_Xk=SecOrderA(Xk);
end

%Obtaining estimations for candidate
Xk=Xk_Op;
dLast = T - sum(Xk);
dvector=[Xk;dLast];
[OmegaAReg,~,~,~,~,OmegaAProbReg,CurrSeedShow]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr,CurrSeedShow);
%uVector used here is the same as in input parameter
[Estim_kOp,Var_kOp,gradientXk_Op,Var_gradientXk_Op]=GradEstim_DP(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg);
Xk=XkAux;


%Statistical Tests
Norm_ProjGrad
if State==1
    Zetak=Deltak*norm(gradient); %Original; Deltak*norm(gradient)
else
    Zetak=0.5*norm(gradient)*min([Deltak norm(gradient)/norm(Hessian)]); %Original: 0.5*norm(gradient)*min([Deltak norm(gradient)/norm(Hessian)])
end

[rejectHo] = SRTest_new(n0,n_k,Var_Xk,Var_kOp,Zetak,alpha_k,Nu0,k,Estim,Estim_kOp);

[Xk_Plus1,Deltak_Plus1,Innerk,Ro_k] = RCnSR_Tests(rejectHo,Estim,Estim_kOp,AproxEstim_Xk,AproxEstim_Xk_Op,Gamma1,Gamma2,Deltak,Xk,Xk_Op,Nu0,Nu1,State);

Innerk
if Innerk==1 %Begin the inner loop
    j=0;
    n_kj=n0;
    DeltakI=Deltak;
    State=2;
    disp('I am in Inner Loop');

    while Innerk==1
        j=j+1;
        
        Deltak=DeltakI;
        n_kaux=n_k;
        n_kj_1=n_kj*(1+ceil(Gamma1^(-2))); % m_ki, sample size of gradient estimator
        n_kj=n_kj_1;
%         n_k=n_kj_1; %Le saque el "aprovechar las simulaciones anteriores" n_k=n_kj_1-n_kaux;
        n_k = 2* n_k; % updated on 15-01, only double the sample size 
        n_k %Show Me

        

     
        %Increase simulations for center solution
        [uVector,CurrSeedServ]=Sampling(npatients,STlower,STupper,n_k,SimDistr,[]);
        Tdvector=transpose(Xk);
        dvectorAux=transpose(Tdvector);
        dLast=T-sum(dvectorAux);
        dvector=[dvectorAux;dLast];
        [OmegaAReg,~,~,~,~,OmegaAProbReg,CurrSeedShow]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr,[]);
        [EstimI,Var_I,gradientXkI,Var_gradient_I]=GradEstim_DP(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg);
        C_measure1 = norm(ProjFeasible(Xk-gradientXkI, T,Epsilon)-Xk)
        [dmin,C_measure2]=CriticalityMeasure(Xk,gradientXkI,1,T,Epsilon); %criticality measure 
        C_measure2
        %%%%THIS WASNT HERE, I ADDED AN UPDATE OF THE HESSIAN
        
%         Xepsi=Xk-Epsilon;
%         [~,~,GradientEpsi,~]=GradEstim_Uniform_ExplicitScenario(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets);
%         Hessian=BFGSHessian(Hessian,Xk,Xepsi,gradientXkI,GradientEpsi,kH);
        %%%%
        
        if OptimMethodInnerLoop==1
            Apat=ones(1,(npatients-1));
            b=T;
            Aeq=[];
            beq=[];
            lb=zeros(npatients-1,1);
            ub=[];
            x0=Xk;
            
            nonlconI=@(x)circlecon(x,Deltak,Xk);
            
            SecOrderAI = @(x) EstimI+transpose(gradientXkI)*(x-Xk)+0.5*transpose(x-Xk)*Hessian*(x-Xk); %TEMPORARY
            [xOptimumI] = fmincon(SecOrderAI,x0,Apat,b,Aeq,beq,lb,ub,nonlconI);
        elseif OptimMethodInnerLoop==2

            [xCauchy,ProjGrad,Norm_ProjGrad,norm_sk]=Cauchy(EstimI,gradientXkI,Hessian,Deltak,State,Xk,T,Epsilon);
           
            xOptimumI=xCauchy;
            SecOrderAI = @(x) EstimI+transpose(gradientXkI)*(x-Xk)+0.5*transpose(x-Xk)*Hessian*(x-Xk); %TEMPORARY
        end
        AproxEstim_XkI=SecOrderAI(Xk);
        AproxEstim_Xk_OpI=SecOrderAI(xOptimumI);
        
        Xk=xOptimumI;
        Tdvector=transpose(Xk);
        dvectorAux=transpose(Tdvector);
        dLast=T-sum(dvectorAux);
        dvector=[dvectorAux;dLast];
        [OmegaAReg,~,~,~,~,OmegaAProbReg]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr,CurrSeedShow);
        
        
        [EstimI_kOp,Var_IOp,gradientIXk_Op,Var_gradient_IOp]=GradEstim_DP(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg);

        Xk=XkAux;
        
        % start statistical test in the inner loop
        % ---------------------------------------------
        RoI_k=(EstimI-EstimI_kOp)/(AproxEstim_XkI-AproxEstim_Xk_OpI);
        
        if RoI_k<Nu0   % does not pass the RC test, shrink the trust region and continue 
            fprintf('Numerator of Ro_k=%f\n',EstimI-EstimI_kOp);
            fprintf('Denominator of Ro_k=%f\n',AproxEstim_XkI-AproxEstim_Xk_OpI);
            fprintf('Norm of s_k=%f\n',norm_sk);
            fprintf('Rho_k=%f\n',RoI_k);
            disp('Sorry, Still Here')
            DeltakI=Gamma1*DeltakI;
             fprintf('DeltakI=%f\n',DeltakI);
            Innerk=1;
            
        elseif Nu0<RoI_k && RoI_k<Nu1 % RC ratio is fair
            disp('Still Here but in the middle')
            Var_Xk=Var_I;
            Var_Xk_Op=Var_IOp;
            Zetak=0.5*norm(gradientXkI)*min([Deltak norm(gradientXkI)/norm(Hessian)]);
            
            [rejectHoI] = SRTest_new(n_k,n_k,Var_Xk,Var_Xk_Op,Zetak,alpha_k,Nu0,k,EstimI,EstimI_kOp);
            
            if rejectHoI==1  % passed the test, go to outer loop
                Xk_Plus1=xOptimumI;
                Deltak_Plus1=DeltakI;
                Estim_kplus1=EstimI_kOp;
                Var_Xkplus1=Var_Xk_Op;
                
                NextGradient=gradientIXk_Op;
                Hessian_kplus1=BFGSHessian(Hessian,Xk_Plus1,Xk,NextGradient,gradientXkI,kH);
                VarGradientSim_kplus1=Var_gradient_IOp;
                
                Innerk=0;
                
            elseif rejectHoI== 0 % continue with the innerloop
                Innerk=1;
                
            else % an local optimal is considered found 
                disp('Local optimal found, stop the program ')
                stopID = 1;
                Xk_Plus1=xOptimumI;
                Deltak_Plus1=DeltakI;
                Estim_kplus1=EstimI_kOp;
                Var_Xkplus1=Var_Xk_Op;
                disp('Criticality measure of current solution:')
                C_measure1 = norm(ProjFeasible(Xk_Plus1-gradientIXk_Op, T, Epsilon)-Xk_Plus1)
                [dmin,C_measure2]=CriticalityMeasure(Xk_Plus1,gradientIXk_Op,1,T,Epsilon);
                C_measure2
                
                NextGradient=gradientIXk_Op;
                Hessian_kplus1=BFGSHessian(Hessian,Xk_Plus1,Xk,NextGradient,gradientXkI,kH);
                VarGradientSim_kplus1=Var_gradient_IOp;
                
                break 
            end
            
        else %Ro bigger than Nu1
            disp('I advanced, but still need Hypothesis Test')
            Var_Xk=Var_I;
            Var_Xk_Op=Var_IOp;
            Zetak=0.5*norm(gradientXkI)*min([Deltak norm(gradientXkI)/norm(Hessian)]);

            [rejectHoI] = SRTest_new(n_k,n_k,Var_Xk,Var_Xk_Op,Zetak,alpha_k,Nu0,k,EstimI,EstimI_kOp);

            if rejectHoI==1
                disp('Got it! :)')
                Xk_Plus1=xOptimumI;
                Deltak_Plus1=DeltakI;
                Estim_kplus1=EstimI_kOp;
                Var_Xkplus1=Var_Xk_Op;
                NextGradient=gradientIXk_Op;
                Hessian_kplus1=BFGSHessian(Hessian,Xk_Plus1,Xk,NextGradient,gradientXkI,kH);
                VarGradientSim_kplus1=Var_gradient_IOp;
                
                Innerk=0;
                
            elseif rejectHoI==0
                Innerk=1;
                
            else 
                disp('Local optimal found, stop the program ')
                stopID = 1
                Xk_Plus1=xOptimumI;
                Deltak_Plus1=DeltakI;
                Estim_kplus1=EstimI_kOp;
                Var_Xkplus1=Var_Xk_Op;
                disp('Criticality measure of current solution:')
                C_measure1 = norm(ProjFeasible(Xk_Plus1-gradientIXk_Op, T,Epsilon)-Xk_Plus1)
                [dmin,C_measure2]=CriticalityMeasure(Xk_Plus1,gradientIXk_Op,1,T,Epsilon);
                C_measure2
                
                NextGradient=gradientIXk_Op;
                Hessian_kplus1=BFGSHessian(Hessian,Xk_Plus1,Xk,NextGradient,gradientXkI,kH);
                VarGradientSim_kplus1=Var_gradient_IOp;
                
                break
            end
            
            %nk=max([nk_1))])
        end
                %break while loop when n_k is too big (modifed on 19/01/08)

       if n_k> 20000000
            disp('Sample size is bigger than 20million, stop')
            stopID = 1

            Estim_kplus1 = Estim;
            Var_Xkplus1 =Var_Xk; 
            NextGradient = gradient;
            Hessian_kplus1 = Hessian; 
            VarGradientSim_kplus1 = Var_gradient_Xk;
            break
        end
    end
    
else %No Inner
    disp('no inner loop needed');
    
    if Xk==Xk_Plus1
        
        %Making the Hessian better
        Xepsi=Xk-Epsilon;
        
        [~,~,GradientEpsi,~]=GradEstim_DP(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg);
        %
        
        Estim_kplus1 = Estim;
        Var_Xkplus1=Var_Xk;
        NextGradient=gradient;
        VarGradientSim_kplus1=Var_gradient_Xk;
        
        Hessian_kplus1=BFGSHessian(Hessian,Xk,Xepsi,gradient,GradientEpsi,kH);
    else
        
        
        Estim_kplus1=Estim_kOp
        Var_Xkplus1=Var_kOp;
        NextGradient=gradientXk_Op;
        VarGradientSim_kplus1=Var_gradientXk_Op;
        
        Hessian_kplus1=BFGSHessian(Hessian,Xk_Plus1,Xk,NextGradient,gradient,kH);
    end
    
    
end

NewSeedServ=CurrSeedServ;
NewSeedShow=CurrSeedShow;




