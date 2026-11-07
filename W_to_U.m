function [ U ] = W_to_U(W)

rho = W(1,:);
u   = W(2,:);
p   = W(3,:);
xi   = W(4,:);
E0   = W(5,:);
e =E0./rho  + p./rho.*xi;

U  = [rho; rho.*u; rho.*e+rho.*(u.^2) /2;xi; E0;rho.*W(6,:);W(7,:)];
end

