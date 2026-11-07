function p = EOS_PR(rho,T)

Ru = 8.31443;
[a,b,M,~,~] = cal_PR(T);

V = M./rho;

p = Ru.*T./(V-b) - a./(V.^2+2*V.*b-b.^2);
end
