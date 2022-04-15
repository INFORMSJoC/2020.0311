function [rejectHo,tstatistic,df,tcomparison] = SRTest(n0,n_k,Var_Xk,Var_kOp,Zetak,alpha_k,Nu0,k,Estim,Estim_kOp)

Sk2=Var_Xk;
Sk_op2=Var_kOp;

S2k=(Sk2/n_k)+(Sk_op2/n0);

df=(S2k^2)*((((Sk2/n_k)^2)/(n_k-1)+((Sk_op2/n0)^2)/(n0-1))^(-1));

tstatistic=((Estim-Estim_kOp)-Zetak*Nu0^2)/(S2k^(0.5));

tcomparison = qt(1-alpha_k,round(df));

%tcomparison = icdf('T',1-alpha_k,round(df));
%load('Matrix.mat')
%tcomparison=Matrix(round(df),k);

if tstatistic>tcomparison
    rejectHo=1;
    
else
    rejectHo=0;
    
end


    
end

