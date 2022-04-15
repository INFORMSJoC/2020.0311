%This code runs the SPGP algorithm discussed in Homem-de-Mello, Kong and
%Godoy (Informs Journal on Computing, 2022)

clear

%%%%Problem parameters:

%Session length (in hours)
SessionLength= 6; 

%Number of patients
NumberOfPatients= 12; 

%Service time distribution: 1: beta distribution 2.5*beta(0.2,0.3), 2:Exponential, 3: Lognormal 
TypeDistr= 2; 

%Mean of the service time distribution, same units as session length (Used only when TypeDistr=2 or 3)
MeanDistr= 1;  

%Variance of the service time distribution (Used only when TypeDistr=3)
VarDistr= 1;   

%Shape of the show-up probability function
% 1:Linear (P1: starting point, P2: end point); 
% 2:Quadratic (P1: starting point, P2: middle point); 
% 3:Cosine Function (P1: Peak, P2: Lowest point)

ShapeShowUpProb= 1;   
P1= 0.1; 
P2= 0.9; 


%Cost parameters

%Cost per unit of waiting time 
CostWait= 0.1; 

%Cost per unit of idle time 
CostIdle= 1; 

%Cost per unit of overtime
CostOver= 1.5;   


%%%%%Code parameters

%Initial sample size
SampleSize= 5000;

%Number of trials for multi-start 
Trials= 5;


%%%%% Run the code
Intento_DP
