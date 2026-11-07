function [ W ] = U_to_W(U)
rho = U(1,:);
u = U(2,:)./ U(1,:);
rho_et = U(3,:);
xi =  U(4,:);
E0 =  U(5,:);

p =1./xi.*(rho_et - E0 - 0.5*rho.*u.^2);

W  = [rho; u; p;xi;E0;U(6,:)./rho;U(7,:)];
end

