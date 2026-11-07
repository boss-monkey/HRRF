function dpdT = EOS_PR_dT(rho,T)

[~,b,M,dadT,~] = cal_PR(T);
V = M./rho;
Ru =8.31443;
dpdT = Ru./(V-b) - dadT./(V.^2+2*V.*b-b.^2);

end

