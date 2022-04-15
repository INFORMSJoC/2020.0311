function [DerTCPerScenario,TCPerScenario,DerTCPerScenario_NoS]=IS_Epsilon_TCnDER(npatients,nAsim,OmegaA,Tdvector,uVector,Cs,Cw,Cl,T)

WMatrixNewFFF=zeros(npatients,nAsim, 1);
WMatrixNewFFFPositions=zeros(npatients,nAsim, npatients);

WDerivatives=zeros(npatients, npatients-1, nAsim); %W's, Derivatives, Scenarios
SDerivatives=zeros(npatients, npatients-1, nAsim); %S's, Derivatives, Scenarios
LDerivativesAUX=zeros(1,npatients-1,nAsim); % 1 (L), Derivatives, Scenarios
LMatrixNewFFFPositions=zeros(1,nAsim,npatients);
%% Pre-Code
OmegaAauxderivative=zeros(npatients, npatients-1, nAsim);
dvectorAux=transpose(Tdvector);

TNon=zeros(npatients, nAsim);
zerosnAsim=zeros(1,nAsim);

dLast=T-sum(dvectorAux);
dvector=[dvectorAux;dLast];

SMatrixNew=zeros(npatients,nAsim);
LMatrixNew=zeros(1,nAsim);

WDerivativesCompNew=zeros(npatients,nAsim,npatients-1);
SDerivativesCompNew=zeros(npatients,nAsim,npatients-1);
LDerivativesCompNew=zeros(1,nAsim,npatients-1);
OmegaAauxderivativeCompNew=zeros(npatients,nAsim,npatients-1);


for j=1:(npatients-1)
    
    OmegaAauxderivativeCompNew(:,:,j)=OmegaA;
    
end


for k=1:nAsim
    
    OmegaAauxderivative(:,:,k)=repmat(OmegaA(:,k),1,npatients-1);
    
end

for k=1:nAsim
    for i=1:npatients
        
        TNon(i,k)=OmegaA(i,k)*uVector(i,k)-dvector(i);
        
    end
end

TNonINV=flipud(TNon);


% Creating the array of sums of T's

for i=1:npatients   %Its correct, not neccesary to check
    if i==1
        TNonINViW=zeros(1,nAsim);
        TNonINViS=-1*TNon(1:i, :);
    elseif i<npatients
        TNonINViW=TNonINV((npatients-i+2):end, :); %Antes TNonINViW=TNonINV((npatients-i+2):end, :)
        TNonINViS=-1*TNon(1:i, :);
    else
        TNonINViW=TNonINV((npatients-i+2):end, :); %Antes TNonINViW=TNonINV((npatients-i+2):end, :)
        TNonINViS=-1*TNon(1:i, :);
        TNonINViL=TNonINV((npatients-i+1):end, :);
    end
    
    if i==1
        CumTNonINViW=cumsum(TNonINViW,1);
        CumTNonINViS=[zerosnAsim;cumsum(TNonINViS,1)]; %Does not work for S's. It is used on L
    elseif i<npatients
        
        CumTNonINViW=[zerosnAsim;cumsum(TNonINViW,1)]; %Temporal vector... (Ex: 0, T2, T2+T1)
        CumTNonINViS=[zerosnAsim;cumsum(TNonINViS,1)]; %Does not work for S's. It is used on L
        
    else
        CumTNonINViW=[zerosnAsim;cumsum(TNonINViW,1)]; %Temporal vector... (Ex: 0, T2, T2+T1)
        CumTNonINViS=[zerosnAsim;cumsum(TNonINViS,1)]; %Does not work for S's. It is used on L
        CumTNonINViL=[zerosnAsim;cumsum(TNonINViL,1)];
        LMatrixNewFFF(1,:,1)=max(CumTNonINViL);
    end
    
    WMatrixNewFFF(i,:,1)=max(CumTNonINViW); %Obtain Maximum
    
    
    for k=1:nAsim
        
        AuxPosition=find(CumTNonINViW(:,k)==WMatrixNewFFF(i,k,1));
        
        
        if length(AuxPosition)==npatients
            WMatrixNewFFFPositions(i,k,:)=AuxPosition;
        else
            WMatrixNewFFFPositions(i,k,:)=[AuxPosition;zeros((npatients-length(AuxPosition)),1)];
            
        end
        
        if i==npatients
            
            AuxPosition_L=find(CumTNonINViL(:,k)==LMatrixNewFFF(1,k,1));
            
            if length(AuxPosition)==npatients
                LMatrixNewFFFPositions(1,k,:)=AuxPosition_L;
            else
                LMatrixNewFFFPositions(1,k,:)=[AuxPosition_L;zeros((npatients-length(AuxPosition_L)),1)];
                
            end
        end
        
        if i==npatients
            
            if (WMatrixNewFFF(npatients,k,1)+TNon(npatients,k))>0
                LMatrixNew(1,k)=WMatrixNewFFF(npatients,k,1)+TNon(npatients,k);
            else
                LMatrixNew(1,k)=0;
            end
        end
        
        if (-WMatrixNewFFF(i,k,1)-TNon(i,k))>=0
            SMatrixNew(i,k)=-WMatrixNewFFF(i,k,1)-TNon(i,k);
        else
            SMatrixNew(i,k)=0;
        end
        
        for j=1:(npatients-1) %Number of derivatives
            
            %%%% W Derivatives
            if i==1
                WDerivatives(i,j,k)=0;
            else
                
                if j<=(i-1)
                    %WMatrixNew(i,k,2)<=(npatients-j)
                    
                    if any(WMatrixNewFFFPositions(i,k,:)>=(i-j+1))% && any(WMatrixNewFFFPositions(i,k,:)>0)
                        
                        WDerivatives(i,j,k)=-1;
                    else
                        WDerivatives(i,j,k)=0;
                    end
                    
                else
                    WDerivatives(i,j,k)=0;
                end
            end
            
            %%%% S Derivatives
            
            
            if (-WMatrixNewFFF(i,k,1)-TNon(i,k))<0
                SDerivatives(i,j,k)=0;
                
            else %-W-T>=0
                if i<=(npatients-1)
                    if j<i
                        SDerivatives(i,j,k)=-WDerivatives(i,j,k);
                    elseif j==i
                        SDerivatives(i,j,k)=1;
                    else
                        SDerivatives(i,j,k)=0;
                    end
                elseif i==npatients
                    SDerivatives(i,j,k)=-WDerivatives(i,j,k)-1;
                end
            end
            
            %%%% L Derivatives
            
            
            if i==npatients
                             
                % The second way to do it
                if any(LMatrixNewFFFPositions(1,k,:)==1) ||  any(LMatrixNewFFFPositions(1,k,:)>=(npatients-j+2))
                    
                    LDerivativesAUX(1,j,k)=0;
                else
                    LDerivativesAUX(1,j,k)=1;
                    
                end
                
             
            end
            
        end
    end
end


for i=1:npatients
    WDerivativesCompNew(i,:,:)= transpose(reshape(WDerivatives(i,:,:),[npatients-1,nAsim]));
    SDerivativesCompNew(i,:,:)= transpose(reshape(SDerivatives(i,:,:),[npatients-1,nAsim]));
    if i==npatients
        LDerivativesCompNew(1,:,:)=transpose(reshape(LDerivativesAUX(1,:,:),[npatients-1,nAsim]));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   reshape(LDerivativesAUX,[1,nAsim,npatients-1]);
    end
end

DerTCPerScenario=Cl*LDerivativesCompNew+Cw*sum(OmegaAauxderivativeCompNew.*WDerivativesCompNew,1)+Cs*sum(SDerivativesCompNew,1);
DerTCPerScenario_NoS=(Cs+Cl)*LDerivativesCompNew+Cw*sum(OmegaAauxderivativeCompNew.*WDerivativesCompNew,1);


TCPerScenario=Cl*LMatrixNewFFF+Cw*sum(WMatrixNewFFF.*OmegaA)+Cs*sum(SMatrixNew,1);

end
