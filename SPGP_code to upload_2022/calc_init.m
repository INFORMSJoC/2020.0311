function [vh,init,W,S] =calc_init(npatients,dvector,uvs,omegas,T)
%Compute beginning time of service for all patients, given the realization
%of service times and no-shows


npatients = npatients + 1; %Sum 1 to account for dummy patient at the end 
%----------comment made on March 07 --------------------------------------
%note that in the main code 'Intento" npatient is already added 1 and now
%npatient means patient 1, 2, ...n+1; but in the main code npatient means
%patient 0, 1, 2, ...n

vh=zeros(2*npatients-1,npatients); %vh records the optimal solution of the LPs,i.e., maximum flow rounte
r=zeros(npatients,1);

init=zeros(npatients,1);%

for j=1:npatients-1
    r(j)=sum(dvector(1:j)); % scheudled arrival time of patient j 
end
r(npatients) = T;
us= uvs .* omegas; %AU-actual service time 

init(1)=r(1);
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
    [sol,opt]=linprog(-f,[],[],Aeq,[zeros(ell-1,1);1;1],lb,ub); %Later we can replace this LP by a recursion
    init(ell)=-opt;
    vh(1:ell,ell)=sol(1:ell);
    vh(npatients+1:npatients+ell-1,ell)=sol(ell+1:2*ell-1);
end

W=init-r; %waiting time 
S(1)=init(1); %idle time 
for j=2:npatients
   S(j)=init(j)-init(j-1)-us(j-1);
end


end
