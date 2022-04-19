function [d,CM] = CriticalityMeasure(X,gradient,theta,T,eps)
% Compute criticality measure min {grad'd : X+d is feasible, ||d||<=theta} (see
% Conn et al. book, p. 452). 
%
%Input: X (feasible point), gradient (gradient of objective function at X),
%theta, T and eps (for teh feasible set  sum(X) <= T-eps, X >= eps

n=size(X,1);
fun=@(d) gradient'*d;

A = ones(1,n);
b = (T-eps)-sum(X); %Constraint sum(X+d) <= T-eps becomes sum(d)<=T-eps - sum(X)
Aeq = ones(1,n);
beq = 0; %feasibilty constriant becomes sum d_i=0
lb = (eps-X).*ones(n,1); %Constraint X+d>=eps becomes d>=eps-X
ub = theta*ones(n,1);  
opts = optimoptions('fmincon','Display','off','Algorithm','sqp');
[d,obj] = fmincon(fun,zeros(n,1),A,b,[],[],lb,ub,@(d) constraint_norm(d,theta),opts);
CM = abs(obj);

end