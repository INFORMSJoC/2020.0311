%this code is used to generate lognormal service distribution 
npatients = 13;
n_k = 1000000;
STlower =1;
STupper = 1; 
me=STlower; %Mean of Lognormal
v=STupper; %Variance of Lognormal
munorm= log((me^2)/sqrt(v+me^2));
sdnorm= sqrt(log(1+v/(me^2)));

uVectorReal=exp(randn(npatients-1,n_k)*sdnorm+munorm); %Lognormal
%uVector=[uVectorDummy;uVectorReal];
uVector = reshape (uVectorReal, (npatients-1)*n_k, 1);
mean(uVector)
std(uVector)
xbins = [00000.1 0.2 0.4 0.6 0.8 1 2 3 4 5 6 7 8 9 10]
hist(uVector, xbins )