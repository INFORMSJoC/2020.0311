function [MeanTC,VarTC,MeanGrad,VarGrad]=GradEstim_Uniform_ExplicitScenario(npatients,STlower,STupper,T,Po,PT,Cw,Cl,Cs,n_k,Xk,SimDistr,Pmed,Xmed,a,b,DerStr,sets,uVector,OmegaAReg,OmegaAProbReg)


%% Vector of decisions

Tdvector=transpose(Xk);

%% Pre-Code

dvectorAux=transpose(Tdvector);

dLast=T-sum(dvectorAux);
dvector=[dvectorAux;dLast];
dvectorReg=dvectorAux;

%% Code
%sets=round(n_k^0.15)+1; %Assigning root. If this is removed, it will take the value defined on the beginning
sets=n_k;  %% MUST USE sets=n_k TO KEEP SAMPLING OUTSIDE THE ROUTINE
m=floor(n_k/sets);
n_kAUX=n_k;

for i=1:m
    if i==1
        n_k=sets;
        
        OmegaA = zeros(npatients, n_k);
        OmegaAauxderivative=zeros(npatients, npatients-1, n_k);
        DerivateAProbReg = zeros(npatients, n_k, npatients-1);
        AuxTCrepReg=zeros(1, n_k, npatients-1);
        OmegaAProbAuxDerReg=zeros(npatients, n_k, npatients-1);
        OmegaAauxderivativeReg=zeros(npatients, npatients-1, n_k);
        
        %[OmegaAReg,~,~,~,~,OmegaAProbReg]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr);
        
        for k=1:n_k
            
            OmegaAauxderivative(:,:,k)=repmat(OmegaA(:,k),1,npatients-1);
            OmegaAauxderivativeReg(:,:,k)=repmat(OmegaAReg(:,k),1,npatients-1);
            
        end
        
        for h=1:(npatients-1)
            
            OmegaAProbAuxDerReg(:,:,h)=OmegaAProbReg;
        end
        %% Pi Probability Derivate Matrix
        
        [DerivateAProbReg]=P_Der(DerStr,OmegaAReg,DerivateAProbReg,PT,Po,T,Pmed,Xmed,dvector,a,b,n_k);
        
        
        %% ///STARTING WITH SIMULATION ITERATIONS///
       
        
        % CALCULATING TC Samples and TC derivative
        
        u=uVector(:,(i-1)*sets+1:(i-1)*sets+sets);
        [DerTCPerScenarioReg,TCPerScenarioReg,~]=IS_Epsilon_TCnDER(npatients,n_k,OmegaAReg,Tdvector,u,Cs,Cw,Cl,T);
        
        for h=1:(npatients-1)
            %     AuxTCrep(:,:,h)=TCPerScenario;
            AuxTCrepReg(:,:,h)=TCPerScenarioReg; %For non IS
        end
        
        AuxsumPReg=transpose(reshape((sum(DerivateAProbReg.*(OmegaAProbAuxDerReg.^(-1)),1)),[n_k,npatients-1]));
        AuxInnerReg=transpose(reshape(DerTCPerScenarioReg,[n_k,npatients-1]))+transpose(reshape(AuxTCrepReg,[n_k,npatients-1])).*AuxsumPReg;
        AuxTotalReg=AuxInnerReg;
        
       
        CurrMeanTC=sum(TCPerScenarioReg)/n_k;
        CurrVarTC=var(TCPerScenarioReg);
        
        CurrMeanGrad=sum(AuxTotalReg,2)/n_k;
        CurrVarGrad=var(AuxTotalReg,0,2);
        
    else
        
        OldMeanTC=CurrMeanTC;
        OldVarTC=CurrVarTC;
        OldMeanGrad=CurrMeanGrad;
        OldVarGrad=CurrVarGrad;
        
        n_k=sets;
        %[OmegaAReg,~,~,~,~,OmegaAProbReg]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr);
        
        for k=1:n_k
            
            OmegaAauxderivative(:,:,k)=repmat(OmegaA(:,k),1,npatients-1);
            OmegaAauxderivativeReg(:,:,k)=repmat(OmegaAReg(:,k),1,npatients-1);
            
        end
        
        for h=1:(npatients-1)
            
            OmegaAProbAuxDerReg(:,:,h)=OmegaAProbReg;
        end
        % Pi Probability Derivate Matrix
        
        [DerivateAProbReg]=P_Der(DerStr,OmegaAReg,DerivateAProbReg,PT,Po,T,Pmed,Xmed,dvector,a,b,n_k);
        
        
        % ///STARTING WITH SIMULATION ITERATIONS///
        
        
        % CALCULATING TC Samples and TC derivative
        u=uVector(:,(i-1)*sets+1:(i-1)*sets+sets);
        [DerTCPerScenarioReg,TCPerScenarioReg,~]=IS_Epsilon_TCnDER(npatients,n_k,OmegaAReg,Tdvector,u,Cs,Cw,Cl,T);
        
        for h=1:(npatients-1)
            %     AuxTCrep(:,:,h)=TCPerScenario;
            AuxTCrepReg(:,:,h)=TCPerScenarioReg; %For non IS
        end
        
        AuxsumPReg=transpose(reshape((sum(DerivateAProbReg.*(OmegaAProbAuxDerReg.^(-1)),1)),[n_k,npatients-1]));
        AuxInnerReg=transpose(reshape(DerTCPerScenarioReg,[n_k,npatients-1]))+transpose(reshape(AuxTCrepReg,[n_k,npatients-1])).*AuxsumPReg;
        AuxTotalReg=AuxInnerReg;
        
        
        % TotalSampleReg=sum(AuxTotalReg,2)/nAsim; %npatients-1,1
        %
        % GradientSamplesReg(:,y)=TotalSampleReg; %CORREGIR PROBLEMA CON Y, YA NO EXISTE
        % GradientSamplesRegEach(:,:,y)=AuxTotalReg; %CORREGIR PROBLEMA CON Y, YA NO EXISTE
        % TCSamplesReg(:,y)=sum(transpose(TCPerScenarioReg))/nAsim; %I.S. Expected Value Estimator, %CORREGIR PROBLEMA CON Y, YA NO EXISTE
        % RecopilandoGradientsE=[RecopilandoGradientsE AuxTotalReg];
        % RecopilandoTCE=[RecopilandoTCE TCPerScenarioReg];
        
        % for i=1:(npatients-1)
        %
        %     AuxVarDerCov(:,i)=reshape(GradientSamplesEach(i,:,:),nAsim*nUsim,1);
        %     AuxVarDerCovReg(:,i)=reshape(GradientSamplesRegEach(i,:,:),nAsim*nUsim,1);
        %
        % end
        
        % CovDerMatrix=cov(AuxVarDerCov)./(nAsim*nUsim);
        
        
        TempMeanTC=sum(TCPerScenarioReg)/n_k;
        TempVarTC=var(TCPerScenarioReg);
        TempMeanGrad=sum(AuxTotalReg,2)/n_k;
        TempVarGrad=var(AuxTotalReg,0,2);
        
      
        
        MeanTC=(sets*(i-1)*OldMeanTC+sets*TempMeanTC)/(sets*i);
        VarTC=(sets*(i-1)*OldVarTC+sets*TempVarTC+sets*(i-1)*((OldMeanTC-MeanTC)^2)+sets*((TempMeanTC-MeanTC)^2))/(sets*i);
        MeanGrad=(sets*(i-1)*OldMeanGrad+sets*TempMeanGrad)/(sets*i);
        VarGrad=(sets*(i-1)*OldVarGrad+sets*TempVarGrad+sets*(i-1)*((OldMeanTC-MeanTC)^2)+sets*((TempMeanTC-MeanTC)^2))/(sets*i);
        
        CurrMeanTC=MeanTC;
        CurrVarTC=VarTC;
        CurrMeanGrad=MeanGrad;
        CurrVarGrad=VarGrad;
    end
 
 end
    if (n_kAUX-sets*m)>0
        
        
        OldMeanTC=CurrMeanTC;
        OldVarTC=CurrVarTC;
        OldMeanGrad=CurrMeanGrad;
        OldVarGrad=CurrVarGrad;
        
        n_k=n_kAUX-m*sets;
        
        OmegaA = zeros(npatients, n_k);
        OmegaAauxderivative=zeros(npatients, npatients-1, n_k);
        DerivateAProbReg = zeros(npatients, n_k, npatients-1);
        AuxTCrepReg=zeros(1, n_k, npatients-1);
        OmegaAProbAuxDerReg=zeros(npatients, n_k, npatients-1);
        OmegaAauxderivativeReg=zeros(npatients, npatients-1, n_k);
        OmegaAReg = zeros(npatients, n_k);
        
        %[OmegaAReg,~,~,~,~,OmegaAProbReg]=IS_Epsilon_SAMPLING(npatients,n_k,Po,PT,T,dvector,dvector,a,b,Pmed,Xmed,DerStr);
        
        for k=1:n_k
            
            OmegaAauxderivative(:,:,k)=repmat(OmegaA(:,k),1,npatients-1);
            OmegaAauxderivativeReg(:,:,k)=repmat(OmegaAReg(:,k),1,npatients-1);
            
        end
        
        for h=1:(npatients-1)
            
            OmegaAProbAuxDerReg(:,:,h)=OmegaAProbReg;
        end
        % Pi Probability Derivate Matrix
        
        [DerivateAProbReg]=P_Der(DerStr,OmegaAReg,DerivateAProbReg,PT,Po,T,Pmed,Xmed,dvector,a,b,n_k);
       
        
          % CALCULATING TC Samples and TC derivative
        
        u=uVector(:,m*sets+1:m*sets+n_k);
        [DerTCPerScenarioReg,TCPerScenarioReg,~]=IS_Epsilon_TCnDER(npatients,n_k,OmegaAReg,Tdvector,u,Cs,Cw,Cl,T);
        
        for h=1:(npatients-1)
            %     AuxTCrep(:,:,h)=TCPerScenario;
            AuxTCrepReg(:,:,h)=TCPerScenarioReg; %For non IS
        end
        
        AuxsumPReg=transpose(reshape((sum(DerivateAProbReg.*(OmegaAProbAuxDerReg.^(-1)),1)),[n_k,npatients-1]));
        AuxInnerReg=transpose(reshape(DerTCPerScenarioReg,[n_k,npatients-1]))+transpose(reshape(AuxTCrepReg,[n_k,npatients-1])).*AuxsumPReg;
        AuxTotalReg=AuxInnerReg;
        
        TempMeanTC=sum(TCPerScenarioReg)/n_k;
        TempVarTC=var(TCPerScenarioReg);
        TempMeanGrad=sum(AuxTotalReg,2)/n_k;
        TempVarGrad=var(AuxTotalReg,0,2);
        
        MeanTC=(sets*m*OldMeanTC+(n_kAUX-sets*m)*TempMeanTC)/(n_kAUX);
        VarTC=(sets*m*OldVarTC+(n_kAUX-sets*m)*TempVarTC+sets*m*((OldMeanTC-MeanTC)^2)+(n_kAUX-sets*m)*((TempMeanTC-MeanTC)^2))/(n_kAUX);
        MeanGrad=(sets*m*(i-1)*OldMeanGrad+(n_kAUX-sets*m)*TempMeanGrad)/(n_kAUX);
        VarGrad=(sets*m*OldVarGrad+(n_kAUX-sets*m)*TempVarGrad+sets*m*((OldMeanTC-MeanTC)^2)+(n_kAUX-sets*m)*((TempMeanTC-MeanTC)^2))/(n_kAUX);
        
 
    else
        MeanTC=CurrMeanTC;
        VarTC=CurrVarTC;
        MeanGrad=CurrMeanGrad;
        VarGrad=CurrVarGrad;
       
      
    end
    
