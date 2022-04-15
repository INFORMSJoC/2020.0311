function [c,ceq] =constraint_norm(d,theta)
    
c = norm(d)-theta;
ceq = [];

end