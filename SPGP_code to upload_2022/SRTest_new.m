function [rejectHo,tstatistic,df,tcomparison] = SRTest_new(n0,n_k,Var_Xk,Var_kOp,Zetak,alpha_k,Nu0,k,Estim,Estim_kOp)
% updated on Jan-14-2019 to include bound for type-II error

Sk2=Var_Xk;
Sk_op2=Var_kOp;

S2k=(Sk2/n_k)+(Sk_op2/n0);

df=(S2k^2)*((((Sk2/n_k)^2)/(n_k-1)+((Sk_op2/n0)^2)/(n0-1))^(-1));

tstatistic=((Estim-Estim_kOp)-Zetak*Nu0^2)/(S2k^(0.5));
tcomparison = qt(1-alpha_k,round(df));
v_alpha = tcomparison;
v_gama = tcomparison; % set the significance level for both types of errors to be the same 
zN = Estim;
zL = Estim_kOp + Zetak*Nu0^2;  
rN = zN/zL;
rhoN = rN -v_alpha *sqrt(S2k)/zL ;
delta = (v_alpha + v_gama) * sqrt (S2k) / zL 


if rhoN>1 % continue to the next iteration 
    rejectHo=1;
  
elseif delta <= 0.1
    rejectHo = 2; % stop the program because local optimal is reached
        
elseif delta >0.1
    rejectHo=0;  % continue with the inner loop
end

%tcomparison = icdf('T',1-alpha_k,round(df));
%load('Matrix.mat')
%tcomparison=Matrix(round(df),k);

% if tstatistic>tcomparison
%     rejectHo=1;
%     
% else
%     rejectHo=0;
%     
% end


    
end

