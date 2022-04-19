%Plot optimal schedules reported in the files schedule_r.mat (r indexes multiple rounds),
% together with the show-up probability function

clear

%Info anout session length, no. of patients
T=6; n=12;

%Info about show-up probability function used
% 1:Linear (P1: starting point, P2: end point); 
% 2:Quadratic (P1: starting point, P2: middle point); 
% 3:Cosine Function (P1: Peak, P2: Lowest point)
ShapeShowUpProb= 3;   
P1= 0.1; 
P2= 0.9; 

%Number of runs
Trials=5;

%Auxiliary variables
c1 = 4*(P1-P2)/T^2;
c2 = -c1*T;
c3 = P1;
a=(P1-P2)/2;
b=(P1+P2)/2;
lower=zeros(1,Trials);
upper=zeros(1,Trials);

k=ShapeShowUpProb;

    for r=1:Trials
        figure(r);
        namefile = strcat('schedule_',num2str(k),'_',num2str(r),'.mat');
        load(namefile);
        if (k==1)
            fplot(@(t) (P1 + (P2-P1)/T*t), [0,T]);  %linear
        elseif (k==2) 
            fplot(@(t) (c1*t^2 + c2*t +c3), [0,T]); %quadratic
        else 
            fplot(@(t) (a*cos(4*pi*t/T)+b), [0,T]); %cosine
        end
        hold on
        probs = P(k,P2,P1,T,0,0,Xk); probs=probs(2:n+1);
        arrtimes= cumsum(Xk); arrtimes=arrtimes(1:n);
        plot(arrtimes,probs,'o');
        hold off
        lower(r)=CI(1);
        upper(r)=CI(3);
    end

%Plot CI for the optimal costs
j=1:Trials;
figure(Trials+1);
plot(j,lower,'ro');
hold on
plot(j,upper,'bo');

