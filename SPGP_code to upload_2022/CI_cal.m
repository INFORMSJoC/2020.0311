%load from excel 

% Xk= xlsread('Simulations_1309-2019.xlsx', 'Cosine', 'AD31:AD43' )
% CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr)
% X0= xlsread('Simulations_1309-2019.xlsx', 'Cosine', 'AC31:AC43' );
% TCVector = xlsread('Simulations_1309-2019.xlsx', 'Cosine', 'AE31:AE46' );
% save ('RESULTADOS_cos27.mat','X0','Xk','TCVector','CI')

%-------------------------------------------------------------------------------
%load from .mat 

%increasing case 
STlower=1;
STupper=1;

SimDistr=2; %1: Uniform, 2:Exponential, 3: Lognormal
%Uniform: STlower and STupper are the limits of the distr.
%Exponential: Mean is STLower
%Lognormal: Mean is STLower, Variance is STUpper

Po=0.9; %0.9;
PT=0.9; %0.1;
DerStr=2; %1; %1:Linear; 2:quadratic (Po: startin point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
% note: cosine Po=0.9, PT= 0.1, and DerStr =3
%%%%
Pmed=0.55; %0.55; %In case piecewise probabilities is used
Xmed= T/2; %floor(npatients/2)+1; %In case piecewise probabilities is used
for i = 1: 10
    if i ==1 
        a = load ('RESULTADOS_v21.mat');
        X0 = a.X0;
        Xk = a.Xvector;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v21.mat','X0','Xk','TCVector', 'CI')
    elseif i ==2
        a = load ('RESULTADOS_v22.mat');
        X0 = a.X0;
        Xk = a.Xvector;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v22.mat','X0','Xk','TCVector', 'CI')

    elseif i ==3
        a = load ('RESULTADOS_v23.mat');
        X0 = a.X0;
        Xk = a.Xvector;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v23.mat','X0','Xk','TCVector', 'CI')
    elseif i == 4
        a = load ('RESULTADOS_v24.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v24.mat','X0','Xk','TCVector', 'CI')
    elseif i ==5 
        a = load ('RESULTADOS_v25.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v25.mat','X0','Xk','TCVector', 'CI')
    elseif i ==6
        a = load ('RESULTADOS_v26.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v26.mat','X0','Xk','TCVector', 'CI')
    elseif i ==7
        a = load ('RESULTADOS_v27.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v27.mat','X0','Xk','TCVector', 'CI')
    elseif i == 8
        a = load ('RESULTADOS_v28.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v28.mat','X0','Xk','TCVector', 'CI')
    
    elseif i == 9
        a = load ('RESULTADOS_v29.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v29.mat','X0','Xk','TCVector', 'CI')
    elseif i ==10
        a = load ('RESULTADOS_v30.mat');
        X0 = a.X0;
        Xk = a.Xk;
        TCVector = a.TCVector;
        CI = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr);
        save ('RESULTADOS_v30.mat','X0','Xk','TCVector', 'CI')
end
end