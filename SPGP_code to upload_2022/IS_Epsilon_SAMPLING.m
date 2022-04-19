function [OmegaA,Pitatoria,PitatoriaIS,PitatioriaPrima,PitatoriaPrimaAux,OmegaAProb,NewSeed]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvectorIS,a,b,Pmed,Xmed,DerStr,CurrSeed)
%No hay arreglar el vector IS porque toma el vector d normal en GradEstim

if isempty(CurrSeed) %Generate new random numbers
    rng('shuffle');
else    %Keep sample
    rng(CurrSeed);
end
NewSeed= rng;

[PvectorIS] = P(DerStr,PT,Po,T,Pmed,Xmed,dvector,a,b); %Genera probabilidades de Show-Up


SimpleUnif=rand(npatients,n_k);

for k=1:n_k
    for i=1:npatients
        if SimpleUnif(i,k)>(1-PvectorIS(i,1))
            OmegaA(i,k)=1;
        elseif SimpleUnif(i,k)<=(1-PvectorIS(i,1))
            OmegaA(i,k)=0;
        end
    end
end

%Forma antigua del sampling
% for k=1:n_k
%     for i=1:npatients
%         OmegaA(i,k)=random('Binomial',1,PvectorIS(i,1),1);
%     end
% end


% Assigning Probabilities

OmegaAProb = zeros(npatients, n_k);
[OmegaAProb] = P_A_function(OmegaAProb,DerStr,PT,Po,T,Pmed,Xmed,dvector,a,b,n_k,OmegaA);


% Assigning Probabilities (Importance Sampling matrix of h(x))

OmegaAProbIS = zeros(npatients, n_k);
[OmegaAProbIS] = P_A_function(OmegaAProbIS,DerStr,PT,Po,T,Pmed,Xmed,dvectorIS,a,b,n_k,OmegaA);



Pitatoria=prod(OmegaAProb); %Resultado de pitatoria (Vector)
PitatoriaIS=prod(OmegaAProbIS);
PitatioriaPrima=Pitatoria./PitatoriaIS; %1,nAsim
PitatoriaPrimaAux=repmat(PitatioriaPrima,npatients-1,1);

end