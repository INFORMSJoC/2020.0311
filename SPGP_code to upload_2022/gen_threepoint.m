% generate uniform discrete distribution with mean 1 and std 1
npatients =6;
n_k= 100000; %sample size 
threepoint_v = [0.2 0.39 2.41];
            temp_rand = unidrnd(3,npatients-1,n_k);%generate discrete numbers 1, 2, 3 uniformly 
            uVector = zeros(1, (npatients-1)*n_k);
            for i = 1:3
                uVector(temp_rand==i) = threepoint_v(i);
            end
            uVectorReal = reshape (uVector, npatients-1, n_k);