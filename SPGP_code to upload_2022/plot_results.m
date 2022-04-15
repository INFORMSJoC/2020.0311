%Plot optimal schedules reported in the files schedule_r.mat (r indexes multiple rounds),
% together with the show-up probability function

clear
nfun=2;
nrounds=5;

P1=0.1; P2=0.9;
T=6; n=12;

c1 = 4*(P1-P2)/T^2;
c2 = -c1*T;
c3 = P1;

lower=zeros(1,nfun*nrounds);
upper=zeros(1,nfun*nrounds);

for k=1:nfun
    for r=1:nrounds
        figure((k-1)*nrounds+r);
        namefile = strcat('schedule_',num2str(k),num2str(r),'.mat');
        load(namefile);
        if (k==1)
            fplot(@(t) (P1 + (P2-P1)/T*t), [0,T]);  %linear
        else fplot(@(t) (c1*t^2 + c2*t +c3), [0,T]); %quadratic
        end
        hold on
        probs = P(k,P2,P1,T,0,0,Xk); probs=probs(2:n+1);
        arrtimes= cumsum(Xk); arrtimes=arrtimes(1:n);
        plot(arrtimes,probs,'o');
        hold off
        lower((k-1)*5+r)=CI(1);
        upper((k-1)*5+r)=CI(3);
        cpu((k-1)*5+r)=cputimer;
    end
end

%Plot CI for the optimal costs
j=1:nfun*nrounds;
figure(nfun*nrounds+1);
plot(j,lower,'ro');
hold on
plot(j,upper,'bo');

%Plot CPU times
figure(nfun*nrounds+2);
plot(j,cpu/60,'o');