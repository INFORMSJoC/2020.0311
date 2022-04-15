function [uVector,NewSeedServ]=Sampling(npatients,STlower,STupper,n_k,SimDistr,CurrSeedServ)

if isempty(CurrSeedServ) %Generate new random numbers
    rng('shuffle');
else    %Keep sample
    rng(CurrSeedServ);
end
NewSeedServ= rng; 


if SimDistr==1 %Beta with mean=std.dev.=STlower
            uVectorDummy=(0.0001)*rand(1,n_k); %Including dummy Serv. Time
%            uVectorReal=(STupper-STlower)*rand(npatients-1,n_k)+STlower;
            uVectorReal=STlower*2.5*betarnd(0.2,0.3,npatients-1,n_k);
            uVector=[uVectorDummy;uVectorReal];
        elseif SimDistr==2 %Exponential with mean=variance=STlower
            uVectorDummy=(0.0001)*rand(1,n_k); %Including dummy Serv. Time
            uVectorReal=-log(rand(npatients-1,n_k))/(1/STlower);%The mean of exp. is STlower
            uVector=[uVectorDummy;uVectorReal];
elseif SimDistr==3 %Lognormal with mean=STlower, variance=STupper
            uVectorDummy=(0.0001)*rand(1,n_k); %Including dummy Serv. Time
            me=STlower; %Mean of Lognormal
            v=STupper; %Variance of Lognormal
            munorm= log((me^2)/sqrt(v+me^2));
            sdnorm= sqrt(log(1+v/(me^2)));

            uVectorReal=exp(randn(npatients-1,n_k)*sdnorm+munorm); %Lognormal
            uVector=[uVectorDummy;uVectorReal];
else % unform distribution 
            uVectorDummy=(0.0001)*rand(1,n_k); %Including dummy Serv. Time
            uVectorReal=5*rand(npatients-1,n_k)+2.5;
            uVector=[uVectorDummy;uVectorReal];
    
        end
        