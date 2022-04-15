function [CI]= out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk,a,b,Pmed,Xmed,DerStr)
%out-of-sample estimation 
%input the final scheduling vector and generated the expected total cost
%and confidence interval 
% first define the shape of show-up function 
%shape =1 - peak; 2, valley; 3,increase; 4 decrease, 5-consine; 6- static

n_sample = 1000000; %sample size 
dvector = Xk;
% 
% 
% Pmed=0.55; %0.55; %In case piecewise probabilities is used
% Xmed= T/2; %floor(npatients/2)+1; %In case piecewise probabilities is used
% a=1; %In case quadratic probabilities is used
% b=2; %In case quadratic probabilities is used
% 
% shape =1;
% if shape == 1 %peak
%     DerStr = 2 ;
%     Po = 0.1;
%     PT = 0.1;
% elseif shape ==2 
%     DerStr = 2 ;
%     Po = 0.9;
%     PT = 0.9;
% elseif shape ==3 %Increas 
%     DerStr = 1 ;
%     Po = 0.1;
%     PT = 0.9;
% elseif shape ==4 %decrease 
%     DerStr = 1 ;
%     Po = 0.9;
%     PT = 0.1;
% % elseif shape ==5 %TO BE Continued
% %     DerStr = 3 ;
% %     Po = 
% %     PT = 
% end
% if shape ==1
%     
% end
    

% 

%input dvector includes the allocated time to the last patient

p_th = 0.95 ; %CI threshold
t = norminv(1-(1-p_th)/2); 

%Step 1: generate service time and show-up status
[uVector,CurrSeedServ]=Sampling(npatients,STlower,STupper,n_sample,SimDistr,[]);
[OmegaAReg,~,~,~,~,OmegaAProbReg,CurrSeedShow]=IS_Epsilon_SAMPLING(npatients,n_sample,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr,[]);
     
TC = zeros(1,n_sample);



for k = 1:n_sample
  
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


for j=1:npatients-1
    r(j)=sum(Xk(1:j)); % scheudled arrival time of patient j 
end
r(npatients) = T;
us= uvs .* omegas; %AU-actual service time 

init(1)=r(1);


for i = 2: npatients-1 % for each patient 
    if r(i) < init(i-1) + us(i-1) % when patient i arrives, service hasn't finished
        init(i) = init(i-1) + us(i-1) ;
    else
        init(i) = r(i);
    end
end
if T<=init(npatients-1) + us(npatients-1) % didn't finish at the end of the session 
    init(npatients) = init(npatients-1) + us(npatients-1);
else
    init(npatients) = T ;
    %d_init(npatients, :) =0
end   


W=init-r; %waiting time 



ds = zeros (npatients,npatients-1); %derivatives of s w.r.t. x_j
S(1)=init(1); %idle time 
for i=2:npatients
   S(i)=init(i)-init(i-1)-us(i-1);
end


TC (k) = W(1:npatients-1)'*omegas*Cw + Cs*sum(S(1:npatients)) + Cl*W(npatients);



end 
% WT = mean (W');
% ST =mean (S');
meanTC = mean(TC');
varTC = var (TC);
stdTC= std(TC);
%confidence interval
CI_lower = meanTC- t*stdTC/sqrt(n_sample);
CI_upper = meanTC + t*stdTC/sqrt(n_sample);
CI =[CI_lower meanTC CI_upper];

