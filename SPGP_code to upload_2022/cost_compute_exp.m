T=6;% session length
npatients = 13
STupper=1;
SimDistr= 2; %1: beta, 2:Exponential, 3: Lognormal
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
s (:,1) = [0.001
0.001
0.001
0.1346
0.3045
0.3991
0.4214
0.4493
0.7833
0.3861
0.6246
0.5718
1.9224
];
%increasing
s(:,2) =[0.12
0.1024
0.1128
0.0885
0.2782
0.1198
0.1349
0.1592
0.4856
0.1878
0.8945
0.6895
2.6269
];
%decreasing
s(:,3) =[0.001
0.0027
0.5373
0.7031
0.8355
0.8989
0.7714
1.1962
1.0335
0.0175
0.001
0.001
0.001
];
%concave
s(:,4) =[0.0989
0.096
0.0756
0.065
0.0319
0.013
0.0263
0.4961
0.8772
1.1765
1.5651
1.4775
0.001
];
%convex
s(:,5) = [0.001
0.001
0.22
0.6117
0.546
0.5506
0.2929
0.4506
0.0947
0.1353
0.4076
0.0581
2.6305
];
%cosine
s(:,6) =[0.001
0.1126
0.5199
0.1998
0.513
0.046
0.0349
0.0509
0.1099
0.1324
0.1214
0.8921
3.2661
];

% constant 0.63 -concave
s(:,7) = [0.001
0.001
0.0411
0.2551
0.4606
0.4275
0.502
0.5579
0.8019
0.3938
0.5813
0.6698
1.3071
]; 
% constant 0.367 -convex
s(:,8) =[0.001
0.001
0.001
0.001
0.0732
0.3415
0.3627
0.3224
0.329
0.2814
0.6701
0.4666
3.1492
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