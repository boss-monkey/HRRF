function [ xi,E0,gammaS, Pi ] = cal_xi_and_E0(c,e,T,p,rho)

h = e+p./rho;
gammaS = 1+c.^2./h;

xi = 1./(gammaS-1);
E0 = e.*rho - p.*xi;
Pi = E0./xi./gammaS;
end


