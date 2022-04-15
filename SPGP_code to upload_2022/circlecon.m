function [c,ceq] = circlecon(x,Deltak,Xk)

Aux=sum((x-Xk).^2);


c = Aux - (Deltak)^2;
ceq = [];

end