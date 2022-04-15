%% note on 2020-11-20
%this code include a beta service distribution, and the code only
%applies when the mean is 1 and std is 1; the sampling code need to be
%changed if the mean/std change. 
% STlower used this as mean and   STupper=1 as std in this code 

%% Problem Parameters

T= SessionLength; %read from main file

realpatientn = NumberOfPatients; %read from main file
% 
% set a real patient number (8,10,12,14) because in the code npatient is real patient number +1

npatients =realpatientn+1;

SimDistr= TypeDistr; %1: beta distribution 2.5*beta(0.2,0.3), 2:Exponential, 3: Lognormal -updated on 2020-11-20
%Beta with mean=std.dev.=STlower
%Exponential: Mean is STLower
%Lognormal: Mean is STLower, Variance is STUpper

STlower=MeanDistr; % use this as mean in this code 
STupper=VarDistr; % use this a std in this code 

DerStr=ShapeShowUpProb; %1; %1:Linear; 2:quadratic (Po: starting point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
% note: static_063 Po=0.9, PT= 0.1, and DerStr =3
%quadratic covex Po= 0.9, PT=0.1; concave: Po=0.1, PT=0.9
Po=P1; %0.9;
PT=P2; %0.1;

%%%%
Pmed=0.5; %0.55; %In case piecewise probabilities is used
Xmed= T/2; %floor(npatients/2)+1; %In case piecewise probabilities is used
a=1; %In case quadratic probabilities is used
b=2; %In case quadratic probabilities is used
%%%%
% cost parameter
Cw=CostWait; % waiting time 
Cs=CostIdle; % idle time 
Cl=CostOver;  % overtime 

for round = 1:Trials
tic 

Epsilon=0.001;

%% Parameters for Strong
%This code has a Dummy Patient, so a 12 patient system will require 12
%variables
% X0 =Epsilon*ones(npatients-1,1);
%X0(1) = T-20*Epsilon;
%X0= T/npatients*ones(npatients-1, 1); %changes made on 19/01/08
% X0= [T/2; 0*ones(npatients-3,1); T/2];
% 
%X0=(T/2)/npatients*ones(npatients-1,1);
X0 =  rand (npatients-1,1)
% X0=[0.0010
% 0.2381
% 0.5582
% 0.4199
% 0.0778
% 0.0981
% 0.0082
% 0.1256
% 0.3462
% 0.4602
% 0.6050
% 0.7632
% ];
X0 = ProjFeasible(X0,T,Epsilon)  %Project to get feasible point

nSTRONGiterations=100; %changes made on 19/01/08
Xvector=[];

OptimMethod=2; %1:=Regular Optimization; 2:=Cauchy 
OptimMethodInnerLoop=2; %1:=Regular Optimization; 2:=Cauchy (FOR INNER LOOP)


n0=SampleSize;
Nu0=0.01; 
Nu1=0.30;
Gamma1=0.90;
Gamma2=1.10;
Delta0=1.5;
DeltaThreshold=0.95;

% n0=100;
sets=316; %Number of samples of A (matrix of no-show). OBS: For now, it doesnt matter. Its just 1.
%samples of service time will be n0/Sets

kH=npatients*1000; %Maximum norm for the hessian matrix
N_It_Max=300;
%Alpha for statistical tests parameters.
AlphaC1=0.6;
AlphaC2=0.95;

%% Strong Code
Innerk=0; %Just for the first iteration

%Hessian0=eye([npatients-1 npatients-1]); % initial of hessian to identical matrix
Hessian0=zeros(npatients-1,npatients-1); 
TCVector=[];%keep track of total cost value
RoVector=[];
ProjGradVector=[];
Nvector=[];
NormProjGradVector=[];
VarTC=[];
VarGradient=[];

k=1;
Norm_ProjGrad=9999;
TimeRecord=0;
Xk=X0; n_k=n0;
%while (k<N_It_Max) && ((TimeRecord)<700) && (Norm_ProjGrad>0.0005)
while (k<N_It_Max) && ((TimeRecord/(60*60)<7)) && (Norm_ProjGrad>0.0005)
  
   [uVector,CurrSeedServ]=Sampling(npatients,STlower,STupper,n0,SimDistr,[]);


%    Tdvector=transpose(Xk);
%    dvectorAux=transpose(Tdvector);
%    dLast=T-sum(dvectorAux);
   dLast = T - sum(Xk);
   dvector=[Xk;dLast];
   [OmegaAReg,~,~,~,~,OmegaAProbReg,CurrSeedShow]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr,[]);
   % OmegaAReg is the show-up status; OmegaAProbReg is the probability of
   % observing a certain show-up status

    if k==1
        
        alpha_k=AlphaC1*(AlphaC2^k); 
        
        Deltak=Delta0;
        Hessian=Hessian0;
        
        if Deltak>DeltaThreshold
            State=1; %1:=First Order; 2:=SecondOrder
        else
            State=2; %1:=First Order; 2:=SecondOrder
        end
        
 
 

        %[Estim1,Var_X1,gradient1,Var_gradient_X1]=GradEstim_Uniform_ExplicitScenario(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg);
        
        [Estim1,Var_X1,gradient1,Var_gradient_X1]=GradEstim_DP(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg);

        Estim=Estim1;
        Var_Xk=Var_X1;
        gradient=gradient1;
        Var_gradient_Xk=Var_gradient_X1;
        
        
        [Xk,Deltak,Estim,Var_Xk,gradient,Hessian,Var_gradient_Xk,Ro_k,ProjGrad,Norm_ProjGrad,Nsampled, stopID,CurrSeedServ,CurrSeedShow]=STRONGITk(npatients,T,Epsilon,STlower,STupper,Po,PT,Cw,Cl,Cs,State,OptimMethod,Nu0,Nu1,Gamma1,Gamma2,n0,alpha_k,Xk,Estim,gradient,Hessian,Deltak,Var_Xk,kH,Var_gradient_Xk,OptimMethodInnerLoop,SimDistr,Pmed,Xmed,a,b,DerStr,k,sets,uVector,CurrSeedServ,CurrSeedShow);
%         [Xk,Deltak,Estim,gradient,Hessian,Var_gradient_Xk,Ro_k,ProjGrad,Norm_ProjGrad,Nsampled, stopID]=STRONGITk(npatients,T,STlower,STupper,Po,PT,Cw,Cl,Cs,State,OptimMethod,Nu0,Nu1,Gamma1,Gamma2,n0,alpha_k,Xk,Estim,gradient,Hessian,Deltak,Var_Xk,kH,Var_gradient_Xk,OptimMethodInnerLoop,SimDistr,Pmed,Xmed,a,b,DerStr,k,sets,uVector,OmegaAReg,OmegaAProbReg);
       
        Xvector=Xk;
        DeltasVector=Deltak;
        Nvector=[Nvector Nsampled];
        RoVector=[RoVector Ro_k];
        %ProjGradVector=[ProjGradVector ProjGrad];
        NormProjGradVector=[NormProjGradVector Norm_ProjGrad];
%         VarTC=[VarTC Var_Xk];
        VarGradient=[VarGradient Var_gradient_Xk];
        k=k+1
    else
        
        alpha_k=AlphaC1*(AlphaC2^k); 
        
        if Deltak>DeltaThreshold %Stage I
            Xk=Xvector(:,k-1);

            State=1; %1:=First Order; 2:=SecondOrder
            
%             [Estim1,SimulationOutcome_X1,gradient1]=GradEstim_Uniform_ExplicitScenario(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,Epsilon,n0,n_k,Xk,State);
%             
%             Estim=Estim1;
%             SimulationOutcome_Xk=SimulationOutcome_X1;
%             gradient=gradient1;
%             
            
            [Xk,Deltak,Estim,Var_Xk,gradient,Hessian,Var_gradient_Xk,Ro_k,ProjGrad,Norm_ProjGrad,Nsampled, stopID,CurrSeedServ,CurrSeedShow]=STRONGITk(npatients,T,Epsilon,STlower,STupper,Po,PT,Cw,Cl,Cs,State,OptimMethod,Nu0,Nu1,Gamma1,Gamma2,n0,alpha_k,Xk,Estim,gradient,Hessian,Deltak,Var_Xk,kH,Var_gradient_Xk,OptimMethodInnerLoop,SimDistr,Pmed,Xmed,a,b,DerStr,k,sets,uVector,CurrSeedServ,CurrSeedShow);
            Nvector=[Nvector Nsampled];
            Xvector=[Xvector Xk];
            DeltasVector=[DeltasVector Deltak];
            RoVector=[RoVector Ro_k];
            %ProjGradVector=[ProjGradVector ProjGrad];
            NormProjGradVector=[NormProjGradVector Norm_ProjGrad];
            VarTC=[VarTC Var_Xk];
            VarGradient=[VarGradient Var_gradient_Xk];
            k=k+1;
        elseif Deltak<=DeltaThreshold %Stage II. We are going to put the innerloop inside STRONGITk
            
            
            Xk=Xvector(:,k-1);

            State=2; %1:=First Order; 2:=SecondOrder
            
%             [Estim1,SimulationOutcome_X1,gradient1]=GradEstim_Uniform_ExplicitScenario(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,Epsilon,n0,n_k,Xk,State);
%             
%             Estim=Estim1;
%             SimulationOutcome_Xk=SimulationOutcome_X1;
%             gradient=gradient1;
%             

            [Xk,Deltak,Estim,Var_Xk,gradient,Hessian,Var_gradient_Xk,Ro_k,ProjGrad,Norm_ProjGrad,Nsampled, stopID,CurrSeedServ,CurrSeedShow]=STRONGITk(npatients,T,Epsilon,STlower,STupper,Po,PT,Cw,Cl,Cs,State,OptimMethod,Nu0,Nu1,Gamma1,Gamma2,n0,alpha_k,Xk,Estim,gradient,Hessian,Deltak,Var_Xk,kH,Var_gradient_Xk,OptimMethodInnerLoop,SimDistr,Pmed,Xmed,a,b,DerStr,k,sets,uVector,CurrSeedServ,CurrSeedShow);
      
            Nvector=[Nvector Nsampled];
            Xvector=[Xvector Xk];
            DeltasVector=[DeltasVector Deltak];
            RoVector=[RoVector Ro_k];
            %ProjGradVector=[ProjGradVector ProjGrad];
            NormProjGradVector=[NormProjGradVector Norm_ProjGrad];
            VarTC=[VarTC Var_Xk];
            VarGradient=[VarGradient Var_gradient_Xk];
                 
            if stopID ==1
                TCVector=[TCVector Estim];
                disp('The Program Stops')
                 %k = N_It_Max;
                break 
            end  
            k=k+1
        end
        
    end
    TCVector=[TCVector Estim];
    k %Show me where you are
    
  %  Estim


%save('RESULTADOS.mat','X0','Xvector','TCVector','Nvector','DeltasVector','RoVector','gradient','ProjGrad','NormProjGradVector','VarTC','VarGradient')
 
end

cputimer = toc

 
X0= [X0; T-sum(X0)]
Xk= [Xk; T-sum(Xk)]
TCVector'
CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr)

namefile = strcat('schedule_',num2str(DerStr),'_',num2str(round),'.mat');
save (namefile,'X0','Xk','TCVector', 'CI', 'cputimer');

end



