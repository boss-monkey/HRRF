function T = cal_T_EOS_PR(p,rho)

T0 = 300;
[ ~,b,M,~,~,~,c1,c2,c3 ] = cal_PR(T0);
V=M./rho;
Ru = 8.31443;
C1 = Ru./(V-b)-1./(V.^2+2.*b.*V-b.^2).*c1;
C2 = 1./(V.^2+2.*b.*V-b.^2).*c2;
C3 = -p-1./(V.^2+2.*b.*V-b.^2).*c3;

D = (C2.^2-4.*C1.*C3).^0.5;
T = ((-C2+D)./(2.*C1)).^2;
end

