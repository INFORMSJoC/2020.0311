%load and plot data
% set up the paragrameter under different cases 

 
N= 20; % number of results
npatient = 14; % number of patients 
clear data_temp
for round = 1 : N
    if round ==1
        data_temp(round) =load ('RESULTADOS_inc_14_1.mat');
    elseif round ==2
        data_temp(round) =load ('RESULTADOS_inc_14_2.mat');
    elseif round ==3
        data_temp(round) =load ('RESULTADOS_inc_14_3.mat');
    elseif round ==4
         data_temp(round) =load ('RESULTADOS_inc_14_4.mat');
    elseif round ==5
         data_temp(round) =load ('RESULTADOS_inc_14_5.mat');
    elseif round ==6
         data_temp(round) =load ('RESULTADOS_inc_14_6.mat');
    elseif round ==7
         data_temp(round) =load ('RESULTADOS_inc_14_7.mat');
    elseif round ==8
         data_temp(round) =load ('RESULTADOS_inc_14_8.mat');
    elseif round ==9
         data_temp(round) =load ('RESULTADOS_inc_14_9.mat');
    elseif round ==10
         data_temp(round) =load ('RESULTADOS_inc_14_10.mat');
    elseif round ==11
        data_temp(round) =load ('RESULTADOS_inc_14_11.mat');
    elseif round ==12
        data_temp(round) =load ('RESULTADOS_inc_14_12.mat');
    elseif round ==13
        data_temp(round) =load ('RESULTADOS_inc_14_13.mat');
    elseif round ==14
         data_temp(round) =load ('RESULTADOS_inc_14_14.mat');
    elseif round ==15
         data_temp(round) =load ('RESULTADOS_inc_14_15.mat');
    elseif round ==16
         data_temp(round) =load ('RESULTADOS_inc_14_16.mat');
    elseif round ==17
         data_temp(round) =load ('RESULTADOS_inc_14_17.mat');
    elseif round ==18
         data_temp(round) =load ('RESULTADOS_inc_14_18.mat');
    elseif round ==19
         data_temp(round) =load ('RESULTADOS_inc_14_19.mat');
    elseif round ==20
         data_temp(round) =load ('RESULTADOS_inc_14_20.mat');
    elseif round ==21
        data_temp(round) =load ('RESULTADOS_inc_14_21.mat');% CI is not calculated 
    elseif round ==22
        data_temp(round) =load ('RESULTADOS_inc_14_22.mat');
    elseif round ==23
        data_temp(round) =load ('RESULTADOS_inc_14_23.mat');
    elseif round ==24
         data_temp(round) =load ('RESULTADOS_inc_14_24.mat');
    elseif round ==25
         data_temp(round) =load ('RESULTADOS_inc_14_25.mat');
    elseif round ==26
         data_temp(round) =load ('RESULTADOS_inc_14_26.mat');
    elseif round ==27
         data_temp(round) =load ('RESULTADOS_inc_14_27.mat');
    elseif round ==28
         data_temp(round) =load ('RESULTADOS_inc_14_28.mat');
    elseif round ==29
         data_temp(round) =load ('RESULTADOS_inc_14_29.mat');
    elseif round ==30
         data_temp(round) =load ('RESULTADOS_inc_14_30.mat');
    end
 
end
meanTC = zeros (N,3);
for i = 1 : N
    meanTC(i,:) = data_temp(i).CI;
end 

% find the minimum m meanTC
m =4; 
cost = zeros (m,1);
ind =  zeros (m,1);
Xk = zeros (npatient + 1, m);
CI = zeros (3,m);

for n = 1: m
    [cost(n) ind(n)] = min (meanTC(:,2));
    meanTC(ind(n),2) = 100;
    Xk(:, n) = data_temp(ind(n)).Xk;
    CI(:,n) = data_temp(ind(n)).CI;    
end
AT =zeros(npatient + 1,m);
AT(1,:) = Xk(1,:); % scheduled arrival time 
for n = 1:npatient
    AT(n+1, :) = AT (n, :) + Xk(n+1,:); 
end
% subplot (2, 1, 1)
% plot (Xk)
% xlabel ('Patients')
% ylabel ('Allocated Time')
% title('Allocated Service Time under Concave Case')
% subplot (2,1,2) 
% plot(AT) 
% xlabel ('Patients')
% ylabel ('Scheduled Arrival Time')
% title('Scheduled Arrival Time under Concave Case')


