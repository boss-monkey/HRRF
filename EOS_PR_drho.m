function dpdrho = EOS_PR_drho(rho,T)

[a,b,M,~,~] = cal_PR(T);
V = M./rho;
Ru = 8.31443;
dpdV = -Ru.*T./(V-b).^2.*(1-2*a.*((Ru.*T.*(V+b).*((V.^2+2*V.*b-b.^2)./(V.^2-b.^2)).^2).^(-1)));
dpdrho = -M./rho.^2.*dpdV;
end
