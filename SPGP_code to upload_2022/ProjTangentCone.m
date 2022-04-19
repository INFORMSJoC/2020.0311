function [xGradP] = ProjTangentCone(x0,xGrad,T )

n=length(x0)+1;

Tras_Grad=xGrad+x0; %Traslated Gradient (this is the one that we need to check)

%Precode

AuxSum=sum(x0)-T;
AuxSumPos=sum(x0(find(x0>=0)))-T; %Function to see whether the assignment feasibility constraint is violated by the positive values
AuxSum_TrasGrad=sum(Tras_Grad)-T; % If its neg, the AF constraint is not violated

Aux_xP=x0;
Aux_xGrad=Tras_Grad;

% % Proj Matrixes
An=[-ones(1,n-2);eye(n-2)]; %For the AF
PAn=An*((transpose(An)*An)^(-1))*transpose(An); %Proj Matrix
PxGrad=PAn*Tras_Grad+(T/(n-1))*ones(n-1,1); %I have to traslate the projected vector

%Determine the active constraints. x0 is the point where the tangent cone
%is generated

A_plus=find(x0==0); %Positions of active non-neg constraints
NegCompGrad_Index=find(Tras_Grad<0); %Identify indexes of neg comp of traslated grad
NegCompProj_Index=find(PxGrad<0); %Identify the indexes of neg comp of the projected gradient

Coincidence=intersect(A_plus,NegCompGrad_Index); % Identify the non-neg active violiting indexes of gradient
Coincidence_Proj=intersect(A_plus,NegCompProj_Index);% Identify the non-neg active violiting indexes of projection

if (numel(A_plus)==0) && (AuxSum<0) %no active constraints
    
    xGradP=Tras_Grad;
    
else
    if (numel(A_plus)==0) && (round2(AuxSum,4)==0)  % Only Assignment F. active
        
        if AuxSum_TrasGrad<=0 % AF not violated by gradient
            xGradP=Tras_Grad;
            
        else %AF violated by the gradient
            
            xGradP=PxGrad; %The projection is has already been calculated
        end
        
    elseif (numel(A_plus)>0) && (AuxSum<0) %One or more Non-Neg constraints active (No AF)
        
        if numel(Coincidence)>0 % There is at least 1 violating component of the active constraints
            
            Aux_xGrad(Coincidence)=0; %Make the negative components zero.
            
            xGradP=Aux_xGrad;
            
            Aux_xGrad=Tras_Grad;
            
        else %There are non-neg active constraints, but the gradient does not violate them
            
            xGradP=Tras_Grad;
        end
        
    else % AF and one or more non neg active constraints
        
        if (numel(A_plus)>0) && (round2(AuxSum,4)==0) %AF and one or more Non-Neg are active
        
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% LONG CODE %%%%%%%%%%%%%%%%%%%%%%%
            
            Aux_xGrad(Coincidence)=0;
            AuxSumPosGrad_Coincidence=sum(Aux_xGrad)-T;
            Aux_xGrad=Tras_Grad;
            %AuxSumPosGrad_All=sum(xGrad)-T;

            
            if (any(Tras_Grad<0)) && (numel(Coincidence)==0)
            
            
                if AuxSum_TrasGrad<0 %AF is not violated by the traslated gradient
                    % (i)
                    xGradP=Tras_Grad;
                    
                else %AF is violated by the traslated gradient
                    xGradP=PxGrad;
                end
                
            elseif (any(Tras_Grad<0)) && (numel(Coincidence)>0)
                
                if AuxSumPosGrad_Coincidence<0 %The non-violating components don't violate the AF
                    %(ii)
                    
                    Aux_xGrad(Coincidence)=0; %Make the negative components zero.
                    xGradP=Aux_xGrad;
                    Aux_xGrad=Tras_Grad;
                    
                else %The non-violating components violate the AF
                    
                    %(v)
                    
                    Aux_xGrad(Coincidence)=0; %Make the negative components zero.
                    %%%##################################################################################################
                    
                    m=numel(Coincidence);
                    
                    if (n-m-2)==0
                        An_Sub=0;
                        AnAux=zeros(n-1,1);
                        
                        Flist=1:(n-1);
                        Flist(:,transpose(Coincidence))=[];
                        AnAux(Flist,:)=An_Sub;
                        
                        Ov=ones(n-1,1);
                        Ov(Coincidence,:)=0;
                        
                        if round(AuxSum,4)==0
                            xGradP=(T/(n-m-1))*Ov;
                        else
                            xGradP=AnAux*((transpose(AnAux)*AnAux)^(-1))*(transpose(AnAux))*Aux_xGrad+(T/(n-m-1))*Ov;
                        end
                        
                    else
                        An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
                        AnAux=zeros(n-1,n-m-2);
                        
                        
                        Flist=1:(n-1);
                        Flist(:,transpose(Coincidence))=[];
                        AnAux(Flist,:)=An_Sub;
                        
                        Ov=ones(n-1,1);
                        Ov(Coincidence,:)=0;
                        
                        xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Aux_xGrad+(T/(n-m-1))*Ov;
                        
                        Auch_Proj=find(xP_FirstProj<0); %%%&%& Last Change
                        ViolatingQuestion=intersect(A_plus,Auch_Proj); %%%&%& Last Change
                        
                        if size(ViolatingQuestion>0) %%%&%&
                            
                            while any(ViolatingQuestion>0) %%%&%&
                                
                                
                                Coincidence=[Coincidence;ViolatingQuestion]; %%%&%&
                                m=numel(Coincidence);
                                
                                if (n-m-2)==0
                                    An_Sub=0;
                                    AnAux=zeros(n-1,1);
                                    
                                    Flist=1:(n-1);
                                    Flist(:,transpose(Coincidence))=[];
                                    AnAux(Flist,:)=An_Sub;
                                    
                                    Ov=ones(n-1,1);
                                    Ov(Coincidence,:)=0;
                                    
                                    if round2(AuxSum,4)==0
                                        xP_FirstProj=(T/(n-m-1))*Ov;
                                        Auch_Proj=find(xP_FirstProj<0); %%%&%&
                                        ViolatingQuestion=intersect(A_plus,Auch_Proj); %%%&%&
                                    else
                                        xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*(transpose(AnAux))*Aux_xGrad+(T/(n-m-1))*Ov;
                                        Auch_Proj=find(xP_FirstProj<0); %%%&%&
                                        ViolatingQuestion=intersect(A_plus,Auch_Proj); %%%&%&
                                    end
                                    
                                else
                                    An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
                                    AnAux=zeros(n-1,n-m-2);
                                    
                                    
                                    Flist=1:(n-1);
                                    Flist(:,transpose(Coincidence))=[];
                                    AnAux(Flist,:)=An_Sub;
                                    
                                    Ov=ones(n-1,1);
                                    Ov(Coincidence,:)=0;
                                    
                                    xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Aux_xGrad+(T/(n-m-1))*Ov;
                                    Auch_Proj=find(xP_FirstProj<0); %%%&%&
                                    ViolatingQuestion=intersect(A_plus,Auch_Proj); %%%&%&
                                end
                                xGradP=xP_FirstProj;
                            end
                        else
                            xGradP=xP_FirstProj;
                            
                        end
                        
                        
                    end
                    
                    
                    %%%##################################################################################################
                end
                
                
                
            else
                
                if (numel(Coincidence)==0) && (AuxSum_TrasGrad>0) && (numel(Coincidence_Proj)>0) %The traslated gradient violates the AF and the projection violates AF
                    %(vi)
                    
                    %%%%%%%%%%%%%%%%####################################################################################
                    
                    
                    m=numel(Coincidence_Proj);
                    
                    if (n-m-2)==0
                        An_Sub=0;
                        AnAux=zeros(n-1,1);
                        
                        Flist=1:(n-1);
                        Flist(:,transpose(Coincidence_Proj))=[];
                        AnAux(Flist,:)=An_Sub;
                        
                        Ov=ones(n-1,1);
                        Ov(Coincidence_Proj,:)=0;
                        
                        if round(AuxSum,4)==0
                            xGradP=(T/(n-m-1))*Ov;
                        else
                            xGradP=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Tras_Grad+(T/(n-m-1))*Ov;
                        end
                        
                        
                    else
                        An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
                        AnAux=zeros(n-1,n-m-2);
                        
                        
                        Flist=1:(n-1);
                        Flist(:,transpose(Coincidence_Proj))=[];
                        AnAux(Flist,:)=An_Sub;
                        
                        Ov=ones(n-1,1);
                        Ov(Coincidence_Proj,:)=0;
                        %                         (T/(n-m-1))*Ov
                        
                        xGradP=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Tras_Grad+(T/(n-m-1))*Ov;
                        
                    end
                    
                    %%%%%%%%%%%%%%%%####################################################################################
                    
                elseif (AuxSum_TrasGrad>0) && (numel(Coincidence_Proj)==0) %The traslated gradient violates the AF and the projection doesn't violate AF
                    
                    
                    xGradP=PxGrad;
                    
                else
                    
                    if (AuxSum_TrasGrad<=0) && (numel(Coincidence_Proj)==0) && (numel(Coincidence)==0)
                        
                        xGradP=Tras_Grad;
                        
                    end
                end
            end
        end
    end
end

end

