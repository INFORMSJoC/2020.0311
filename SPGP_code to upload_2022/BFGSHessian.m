function [Hessian_1] = BFGSHessian(Hessian,Xk,Xk_Op,gradient,gradientXk_Op,kH)


s=Xk_Op-Xk;
y=gradientXk_Op-gradient;

HessianAUX=Hessian-((Hessian*s*transpose(s)*Hessian)/(transpose(s)*Hessian*s))+((y*transpose(y))/(transpose(y)*s));


if norm(HessianAUX)<kH
    
    Hessian_1=HessianAUX;
    
else
    
    Hessian_1=kH*HessianAUX/norm(HessianAUX);
% ����Esta mal. Esta tomando la hessiana k-1. Hay que cambiarlo por:
% ����Hessian_1=kH*HessianAUX/norm(HessianAUX);

end

Hessian_1=Hessian-Hessian;
end
