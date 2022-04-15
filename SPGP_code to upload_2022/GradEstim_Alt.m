function [meanTC,meanD] =GradEstim_Alt(npatients,n_k,Xk,uVector,OmegaAReg,T,Cw,Cs,Cl,DerStr,PT,Po,Pmed,Xmed,a,b)%Sample size input

TC = zeros(1,n_k);
dTC = zeros(npatients-1,n_k);

Tdvector=transpose(Xk);
dvectorAux=transpose(Tdvector);
dLast=T-sum(dvectorAux);
dvector=[dvectorAux;dLast];
for k = 1 : n_k
    
W = zeros(npatients,1);
S = zeros (npatients, 1);
%Compute beginning time of service for all patients, given the realization
%of service times and no-shows
 % first dummy patient is not included
uvs =uVector(2:npatients,k);
omegas = OmegaAReg(2:npatients,k);
 %Sum 1 to account for dummy patient at the end 
%----------comment made on March 07 --------------------------------------
%note that in the main code 'Intento" npatient is already added 1 and now
%npatient means patient 1, 2, ...n+1; but in the main code npatient means
%patient 0, 1, 2, ...n

vh=zeros(2*npatients-1,npatients); %vh records the optimal solution of the LPs,i.e., maximum flow route
r=zeros(npatients,1); % r is the arrival time of patient i

init=zeros(npatients,1);%



for j=1:npatients-1
    r(j)=sum(dvector(1:j)); % scheudled arrival time of patient j 
end
r(npatients) = T;
us= uvs .* omegas; %AU-actual service time 

init(1)=r(1);
vh(1,1) =1 ;
for ell=2:npatients
    f=[r(1:ell); us(1:ell-1)];
    clear Aeq
    Aeq(1,:)=[1 zeros(1,ell-1) -1 zeros(1,ell-2)];
    for j=2:ell-1
        Aeq=[Aeq; zeros(1,j-1) 1 zeros(1,ell-j) zeros(1,j-2) 1 -1 zeros(1,ell-j-1)];
    end
    Aeq=[Aeq; zeros(1,ell-1) 1 zeros(1,ell-2) 1];
    Aeq=[Aeq; ones(1,ell) zeros(1,ell-1)];
    lb=zeros(2*ell-1,1);
    ub=ones(2*ell-1,1);
    options = optimset('linprog');
    options.Display = 'off';
    [sol,opt]=linprog(-f,[],[],Aeq,[zeros(ell-1,1);1;1],lb,ub,options); %Later we can replace this LP by a recursion
    init(ell)=-opt;
    vh(1:ell,ell)=sol(1:ell);
    vh(npatients+1:npatients+ell-1,ell)=sol(ell+1:2*ell-1);
end

W=init-r; %waiting time 
S(1)=init(1); %idle time 
for j=2:npatients
   S(j)=init(j)-init(j-1)-us(j-1);
end

% calculate the derivatives of the objective function w.r.t. x_j
dw= zeros (npatients,npatients-1); %derivaties of w w.r.t. x_j
ds = zeros (npatients,npatients-1); %derivatives of s w.r.t. x_j
dw(1,1)=0;
ds(1,1)=1;
for i = 2: npatients-1 % derivatives of w_i wrt x_j
    for j = 1:i % if j>i, then dw(i,j)=0; ds(i,j)=0
        dw(i,j) = sum (vh(j:i ,i))-1;
        ds(i,j) = sum(vh(j:i,i))-sum(vh(j:i-1,i-1));
    end 
end

for j = 1:npatients-1
    dw(npatients,j) = sum (vh(j:npatients-1 ,npatients));
    ds(npatients,j) = sum(vh(j:npatients-1,npatients)) - sum(vh(j:npatients-1,npatients-1));
end

TC (k) = W(1:npatients-1)'*omegas*Cw + Cs*sum(S(1:npatients)) + Cl*W(npatients);
[PvectorIS, dp] = DP(DerStr,PT,Po,T,Pmed,Xmed,dvector,omegas, a,b);

%calculate the derivatives 
dTC1=zeros(npatients-1,1);
dTC2 = zeros(npatients-1,1);
dTC1(1) = Cl*dw(npatients,1) + Cs + Cw * omegas (2:npatients-1)'*dw(2:npatients-1,1) + Cs* sum(ds(2:npatients,1)); 
dTC2(1) = TC(k) *sum(dp(:,1)./PvectorIS(:));
for j = 2 : npatients-1 
    dTC1 (j) = Cl*dw(npatients,j)  + Cw * omegas (2:npatients-1)'*dw(2:npatients-1,j) + Cs* sum(ds(2:npatients,j)); 
    dTC2(j) = TC(k)*sum(dp(:,j)./PvectorIS(:));
end

dTC(:,k) = dTC1+dTC2;


end 
WT = mean (W');
ST =mean (S');
meanTC = mean(TC');
%now calculate the 2nd part of the derivative wrt show-up probability 

if n_k ==1
    meanD = dTC';
else
    meanD = mean(dTC');
end
%calculate gradients of the objective function value 



