function [xP] = ProjFeasible(x0,T,eps)
%SOlve min_y {||x0-y||^2 : sum y_i<=T-eps, y>=eps}

%n=size(x0,1)+1;

%Set up coefficents for QP
n=size(x0,1);
%eps=1e-3; 
f = -x0;
H = eye(n);
A = ones(1,n);
lb = eps*ones(n,1); %make sure we dont't hit the axes
H = eye(n);
opts = optimoptions('quadprog','Display','off');
xP = quadprog(H,f,A,T-eps,[],[],lb,[],[],opts);
% %Precode
% 
% AuxSum=sum(x0)-T;
% AuxSumPos=sum(x0(find(x0>=0)))-T;
% Aux_xP=x0;
% 
% %Code
% 
% if AuxSumPos<=0 %The the assignment feasibility constraint is not violated by the positive values
%     if any(x0<0)
%         Aux_xP(find(x0<0))=0;
%         xP=Aux_xP;
%         
%     else %Not violated and all the components are non-negative
%         
%         %Do nothing, none of the constraits are violated
%         xP=x0;
%     end
%     
% else %The the assignment feasibility constraint is violated by the positive values
%     
%     if any(x0<0)
%         NegLoc=find(x0<0);
%         Aux_xP(NegLoc)=0; %Make the negative components zero.
%         
%         m=numel(NegLoc);
%         
%         if (n-m-2)==0
%             An_Sub=0;
%             AnAux=zeros(n-1,1);
%             
%             Flist=1:(n-1);
%             Flist(:,transpose(NegLoc))=[];
%             AnAux(Flist,:)=An_Sub;
%             
%             Ov=ones(n-1,1);
%             Ov(NegLoc,:)=0;
%             
%             if AnAux==0
%                 xP=(T/(n-m-1))*Ov;
%             else
%                 xP=AnAux*((transpose(AnAux)*AnAux)^(-1))*(transpose(AnAux))*Aux_xP+(T/(n-m-1))*Ov;
%             end
%             
%         else
%             An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
%             AnAux=zeros(n-1,n-m-2);
%             
%             
%             Flist=1:(n-1);
%             Flist(:,transpose(NegLoc))=[];
%             AnAux(Flist,:)=An_Sub;
%             
%             Ov=ones(n-1,1);
%             Ov(NegLoc,:)=0;
%             
%             xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Aux_xP+(T/(n-m-1))*Ov;
%             
%             if any(xP_FirstProj<0)
%                 
%                 while any(xP_FirstProj<0)
%                     
%                     SecondNegLoc=find(xP_FirstProj<0);
%                     NegLoc=[NegLoc;SecondNegLoc];
%                     m=numel(NegLoc);
%                     
%                     if (n-m-2)==0
%                         An_Sub=0;
%                         AnAux=zeros(n-1,1);
%                         
%                         Flist=1:(n-1);
%                         Flist(:,transpose(NegLoc))=[];
%                         AnAux(Flist,:)=An_Sub;
%                         
%                         Ov=ones(n-1,1);
%                         Ov(NegLoc,:)=0;
%                         
%                         if AnAux==0
%                             xP_FirstProj=(T/(n-m-1))*Ov;
%                         else
%                             xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*(transpose(AnAux))*Aux_xP+(T/(n-m-1))*Ov;
%                         end
%                         
%                     else
%                         An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
%                         AnAux=zeros(n-1,n-m-2);
%                         
%                         
%                         Flist=1:(n-1);
%                         Flist(:,transpose(NegLoc))=[];
%                         AnAux(Flist,:)=An_Sub;
%                         
%                         Ov=ones(n-1,1);
%                         Ov(NegLoc,:)=0;
%                         
%                         xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Aux_xP+(T/(n-m-1))*Ov;
%                     end
%                     xP=xP_FirstProj;
%                 end
%             else
%                 xP=xP_FirstProj;
%                 
%             end
%             
%             
%         end
%         
%     else %All of them are positive
%         
%         An=[-ones(1,n-2);eye(n-2)];
%         PAn=An*((transpose(An)*An)^(-1))*transpose(An); %I have to traslate the projected vector
%         Px0=PAn*x0+(T/(n-1))*ones(n-1,1);
%         
%         if any(Px0<0) %Projection onto feasibility contraint plane Test has some negative components
%             NegLoc_Proj=find(Px0<0);
%             m=numel(NegLoc_Proj);
%             
%             if (n-m-2)==0
%                 An_Sub=0;
%                 AnAux=zeros(n-1,1);
%                 
%                 Flist=1:(n-1);
%                 Flist(:,transpose(NegLoc_Proj))=[];
%                 AnAux(Flist,:)=An_Sub;
%                 
%                 Ov=ones(n-1,1);
%                 Ov(NegLoc_Proj,:)=0;
%                 
%                 if AnAux==0
%                     xP=(T/(n-m-1))*Ov;
%                 else
%                     xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*x0+(T/(n-m-1))*Ov;
%                 end
%                 
%                 
%             else
%                 An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
%                 AnAux=zeros(n-1,n-m-2);
%                 
%                 
%                 Flist=1:(n-1);
%                 Flist(:,transpose(NegLoc_Proj))=[];
%                 AnAux(Flist,:)=An_Sub;
%                 
%                 Ov=ones(n-1,1);
%                 Ov(NegLoc_Proj,:)=0;
%                 
%                 xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*x0+(T/(n-m-1))*Ov;
%                 
%             end
%             
%                        
%             if any(xP_FirstProj<0)
%                 
%                 while any(xP_FirstProj<0)
%                     
%                     SecondNegLoc=find(xP_FirstProj<0);
%                     NegLoc=[NegLoc_Proj;SecondNegLoc];
%                     m=numel(NegLoc);
%                     
%                     if (n-m-2)==0
%                         An_Sub=0;
%                         AnAux=zeros(n-1,1);
%                         
%                         Flist=1:(n-1);
%                         Flist(:,transpose(NegLoc))=[];
%                         AnAux(Flist,:)=An_Sub;
%                         
%                         Ov=ones(n-1,1);
%                         Ov(NegLoc,:)=0;
%                         
%                         if AnAux==0
%                             xP_FirstProj=(T/(n-m-1))*Ov;
%                         else
%                             xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*(transpose(AnAux))*Aux_xP+(T/(n-m-1))*Ov;
%                         end
%                         
%                     else
%                         An_Sub=[-ones(1,n-m-2);eye(n-m-2)];
%                         AnAux=zeros(n-1,n-m-2);
%                         
%                         
%                         Flist=1:(n-1);
%                         Flist(:,transpose(NegLoc))=[];
%                         AnAux(Flist,:)=An_Sub;
%                         
%                         Ov=ones(n-1,1);
%                         Ov(NegLoc,:)=0;
%                         
%                         xP_FirstProj=AnAux*((transpose(AnAux)*AnAux)^(-1))*transpose(AnAux)*Aux_xP+(T/(n-m-1))*Ov;
%                     end
%                     xP=xP_FirstProj;
%                 end
%             else
%                 xP=xP_FirstProj;
%                 
%             end
%             
%         else % The projected vector has all its components positive
%             
%             xP=Px0; %Proj. onto the assignment feasibility plane
%             
%         end
%         
%     end
%     
% end

end

