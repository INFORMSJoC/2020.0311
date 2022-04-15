function [meanTC,varTC,meanD,varD]=GradEstim_DP(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg)
        

%This code is used to calculate the expected totaal cost and derivatives
%for a set of given U and Omegas
%the derivative calculation uses dynamic programming 



TC = zeros(1,n_k);
dTC = zeros(npatients-1,n_k);
dTC1 = zeros(npatients-1,n_k);
dTC2 = zeros(npatients-1,n_k);

dLast = T - sum(Xk);
dvector=[Xk;dLast];

for k = 1:n_k
  
%clear W S  init d_init r dw ds 
W = zeros(npatients,1);
S = zeros (npatients, 1);
%Compute beginning time of service for all patients, given the realization
%of service times and no-shows
 % first dummy patient is not included
uvs =uVector(2:npatients,k);

omegas = OmegaAReg(2:npatients,k); %realized show-up status


 %Sum 1 to account for dummy patient at the end 
%----------comment made on March 07 --------------------------------------
%note that in the main code 'Intento" npatient is already added 1 and now
%npatient means patient 1, 2, ...n+1; but in the main code npatient means
%patient 0, 1, 2, ...n

%vh=zeros(2*npatients-1,npatients); %vh records the optimal solution of the LPs,i.e., maximum flow route for patient i
%vh(1: npatients)are the flow on the vectical arcs and vh(npatients+1: 2*npatients-1) is the flow on horizontal arcs 
r=zeros(npatients,1); % r is the arrival time of patient i

init=zeros(npatients,1);%
d_init = zeros (npatients, npatients-1); % start time of patient i w.r.t. x_j


for j=1:npatients-1
    r(j)=sum(Xk(1:j)); % scheudled arrival time of patient j 
end
r(npatients) = T;
us= uvs .* omegas; %AU-actual service time 

init(1)=r(1);
d_init(1,1) = 1;

for i = 2: npatients-1 % for each patient 
    if r(i) < init(i-1) + us(i-1) % when patient i arrives, service hasn't finished
        init(i) = init(i-1) + us(i-1) ;
        d_init(i,:) =d_init(i-1,:);
    else
        init(i) = r(i);
        d_init(i,1:i) =1;
    end
end
if T<=init(npatients-1) + us(npatients-1) % didn't finish at the end of the session 
    init(npatients) = init(npatients-1) + us(npatients-1);
    d_init(npatients,:) =d_init(npatients-1,:);
else
    init(npatients) = T ;
    %d_init(npatients, :) =0
end   


W=init-r; %waiting time 
dw= zeros (npatients,npatients-1); %derivaties of w w.r.t. x_j
for i = 2:npatients-1
    
      dw (i,1:i) =d_init(i,1:i) -1;
      dw(i,i+1:npatients-1) =d_init(i,i+1:npatients-1);
end

dw(npatients, :) = d_init(npatients,:); 

ds = zeros (npatients,npatients-1); %derivatives of s w.r.t. x_j
S(1)=init(1); %idle time 
ds(1,1)=1;
for i=2:npatients
   S(i)=init(i)-init(i-1)-us(i-1);
   ds(i,: ) = d_init(i,:) - d_init(i-1,:);
end



TC (k) = W(1:npatients-1)'*omegas*Cw + Cs*sum(S(1:npatients)) + Cl*W(npatients);
[PvectorIS, dp] = DP(DerStr,PT,Po,T,Pmed,Xmed,dvector,omegas, a,b);

%calculate the derivatives 
% dTC1=zeros(npatients-1,1);
% dTC2 = zeros(npatients-1,1);
%dTC1(1) = Cl*dw(npatients,1) + Cs + Cw * omegas (1:npatients-1)'*dw(1:npatients-1,1) + Cs* sum(ds(1:npatients,1)); 
%dTC2(1) = TC(k) *sum(dp(:,1)./PvectorIS(:));
for j = 1 : npatients-1 
    dTC1 (j,k) = Cl*dw(npatients,j)  + Cw * omegas (1:npatients-1)'*dw(1:npatients-1,j) + Cs* sum(ds(1:npatients,j)); 
    dTC2(j,k) = TC(k)*sum(dp(:,j)./PvectorIS(:));
end

dTC(:,k) = dTC1(:,k)+dTC2(:,k);


end 
% WT = mean (W');
% ST =mean (S');
meanTC = mean(TC');

%now calculate the 2nd part of the derivative wrt show-up probability 

if n_k ==1
    meanD = dTC';
    meanTC1 = dTC1'
else
    meanD = mean(dTC');
    meanTC1 = mean (dTC1');
   
end
meanD = meanD';
varTC = var (TC);
varD = var (dTC');
varD =varD';
%calculate gradients of the objective function value 



