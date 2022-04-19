T=6;% session length
npatients = 13
STupper=1;
SimDistr= 3; %1: beta, 2:Exponential, 3: Lognormal
%Uniform: STlower and STupper are the limits of the distr.
%Exponential: Mean is STLower
%Lognormal: Mean is STLower, Variance is STUpper 
 

%%%%
Pmed=0.55; %0.55; %In case piecewise probabilities is used
Xmed= T/2; %floor(npatients/2)+1; %In case piecewise probabilities is used
a=1; %In case quadratic probabilities is used
b=2; %In case quadratic probabilities is used
%%%%
% cost parameter
Cw=0.1;
Cl=1.50;
Cs=1;

s = zeros (13, 8); % input schedule solutions for the six show-up prob.
at =zeros (12,8); % arrival time of patient 1 to 12;
%constant 0.5
s (:,1) = [0.0010
    0.0010
    0.0174
    0.2111
    0.3690
    0.4800
    0.3916
    0.7198
    0.4183
    0.6538
    0.3912
    0.5733
    1.7726];
%increasing
s(:,2) =[0.1842
0.1588
0.1789
0.1772
0.186
0.1987
0.2067
0.1945
0.4444
0.4375
0.5744
0.5965
2.4623
];
%decreasing
s(:,3) =[0.001
0.0718
0.6135
0.6793
0.6485
0.9703
0.728
0.7404
0.8858
0.6584
0.001
0.001
0.001
];
%concave
s(:,4) =[ 0.0010
    0.0010
    0.0010
    0.0010
    0.0010
    0.0010
    0.0010
    0.2433
    0.8823
    1.4420
    1.1836
    1.1969
    1.0449
];
%convex
s(:,5) = [0.0010
    0.0010
    0.2798
    0.4199
    0.4697
    0.4363
    0.5782
    0.4504
    0.1185
    0.2909
    0.3007
    0.0766
    2.5769];
%cosine
s(:,6) = [0.0010
    0.2013
    0.6928
    0.6406
    0.0787
    0.0352
    0.0010
    0.2018
    0.4430
    0.5761
    0.6883
    0.9601
    1.4803
];

% constant 0.6333 -concave
s(:,7) = [ 0.0010
    0.0010
    0.0191
    0.2451
    0.2838
    0.4395
    0.6907
    0.4108
    0.6628
    0.3146
    0.5394
    0.7635
    1.6288]; 
% constant 0.367 -convex
s(:,8) =[       0.0010
    0.0010
    0.0010
    0.0294
    0.1539
    0.1863
    0.4195
    0.4466
    0.2644
    0.5775
    0.3183
    0.5971
    3.0041
];

caseid = input ('Enter a show-up pattern:')
 switch caseid
     case 1
         disp('We are calculating an increasing case')
         Po=0.1;
         PT=0.9; 
         DerStr=1; %1; %1:Linear; 2:quadratic (Po: starting point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
         Xk_constant = s (1:13, 1); %constant solution in the increasing environment 
         Xk_optimal = s (1:13, 2); %increasing solution 
     case 2
         disp ('We are calculating a decreasing case')
         Po=0.9; 
         PT=0.1; 
         DerStr=1; %1; %1:Linear; 2:quadratic (Po: starting point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
         Xk_constant = s (1:13, 1); %constant solution in the decreasing environment 
         Xk_optimal = s (1:13, 3); %decreasing solution 
     case 3
         disp ('We are calculating a quadratic concave case')
         Po=0.1; 
         PT=0.9; 
         DerStr=2; %1; %1:Linear; 2:quadratic (Po: starting point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
         Xk_constant = s (1:13, 7); %constant solution in the concave environment 
         Xk_optimal = s (1:13, 4); %concave solution 
     case 4
         disp ('We are calculating a quadratic convex case')
         Po=0.9; 
         PT=0.1; 
         DerStr=2; %1; %1:Linear; 2:quadratic (Po: starting point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
         Xk_constant = s (1:13, 8); %constant solution in the convex environment 
         Xk_optimal = s (1:13, 5); %convex solution 
     case 5 
         disp ('We are calculating a cosine case')
         Po=0.9; 
         PT=0.1; 
         DerStr=3; %1; %1:Linear; 2:quadratic (Po: starting point, PT middle point)=2019-10-09; 3: Cosine Funct (Po: Peak, PT: Lowest point)
         Xk_constant = s (1:13, 1); %constant solution in the cosine environment 
         Xk_optimal = s (1:13, 6); %cosine solution 
 end

CI_constant = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk_constant,a,b,Pmed,Xmed,DerStr);
CI_optimal  = out_of_sample_etimation (Cw,Cl,Cs,npatients,Po,PT,T,STlower,STupper,SimDistr,Xk_optimal,a,b,Pmed,Xmed,DerStr);
% at (1,:) = s(1,:);
[round(CI_constant,3) round(CI_optimal,3)]
(CI_constant(2)- CI_optimal(2))/CI_constant(2)
% 
% for i = 1: 11 
%     at (i+1,:) =at(i,:) + s(i+1,:);  
% end