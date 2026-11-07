clear all
clc

p_0 =2e6; p_1 = 1e5;
rho_0 =400; rho_1 =2;
u_0 = 50; u_1 = 0;
% 
tau_M = 1e-5; % 
tau_T = 1e-4; %
%%
N = 1000;
x1 = 0.0; x2 = 1; x3 = 0.5; x4 = 1;
dx = (x2-x1)/N;
x = linspace(x1,x2,N);
epsilon =1.0*dx; %  Smooth initialization enhances initial convergence
if epsilon > 0
    phi_x3 = 0.5*(1-tanh((x-x3)/epsilon));
    smooth_cells = 2;
    smooth_width = smooth_cells*dx;
    phi_x3(x<x3-smooth_width) = 1.0;
    phi_x3(x>x3+smooth_width) = 0.0;
else
    phi_x3 = double(x <= x3);
end
rho0 = rho_0.*phi_x3+rho_1.*(1-phi_x3);
u0 = u_0.*phi_x3+u_1.*(1-phi_x3);
p0 = p_0.*phi_x3+p_1.*(1-phi_x3);

%%

% 
CFL = 0.5;
source_CFL = 0.5;
use_phase_merge = true; 

phase_tol = 1e-8; 
flag = 1;
time = 1e-3;
maxSteps = 5000000;
steps = 0;
current_time = 0;
energy_correction = 0;
phase_mass_correction = 0;
phase_mass_correction_abs = 0;
phase_volume_correction = 0;
run_error = [];
ids = 3:N-2;

result = [x',rho0',u0',p0',NaN(N,8)];
e_error = NaN(1,N);
tic
try
    [beta0,alpha0] = cal_initial_beta_alpha(p0,rho0);
    [c0,e0,T0,Tv0,rhol0,rhov0,Kl,Kv,Gamma_l,Gamma_v,cpl,cpv,hl,hv,c2l,c2v] = ...
        cal_eq_parameter_NN(p0,rho0,beta0,alpha0,tau_M);
    [xi,E0,gammaS,Pi] = cal_xi_and_E0(c0,e0,T0,p0,rho0);
    W_0 = [rho0;u0;p0;xi;E0;beta0;alpha0];
    U_0 = W_to_U(W_0);
    W_1 = W_0;
    c1 = c0; e1 = e0; u1 = u0; T1 = T0; Tv1 = Tv0;
    rhol1 = rhol0; rhov1 = rhov0;
    e2 = e0; ep = zeros(1,N);
    e_error = zeros(1,N);
    result = [x',rho0',u0',p0',T0',c0',T0',Tv0',beta0',alpha0',rhol0',rhov0'];
    disp('Initialization complete, starting calculation');

    while flag
        if current_time >= time, break; end
        if steps >= maxSteps
            error('HRM7:MaxSteps','Reached maxSteps, current physical time %.9g s.',current_time);
        end

        [R,Q,R_P,K,R_alpha,tau_M_e,tau_M_c,tau_T_local] = cal_hrm_source( ...
            W_1,T1,Tv1,rhol1,rhov1,Kl,Kv,Gamma_l,Gamma_v,cpl,cpv,hl,hv,c2l,c2v, ...
            tau_M,tau_T);
        F_U = cal_flux(W_1,N,K);
        S_U = zeros(7,N);
        S_U(5,:) = -W_1(4,:).*R_P;
        S_U(6,:) = W_1(1,:).*R;
        S_U(7,:) = R_alpha;
        rhs = -F_U/dx+S_U;

        dt = CFL*dx/max(max(abs(u1)+c1));
        % The source-term CFL uses the local, actually activated mass- and heat-relaxation times in the current step.
        tau_active=min(tau_M_e(ids),tau_M_c(ids));
        pick=R(ids)~=0 & isfinite(tau_active);
        if any(pick), dt=min(dt,source_CFL*min(tau_active(pick))); end
        tau_active=tau_T_local(ids);
        pick=Q(ids)~=0 & isfinite(tau_active);
        if any(pick), dt=min(dt,source_CFL*min(tau_active(pick))); end
        dt = min(dt,time-current_time); 
        if ~isfinite(dt) || dt <= 0 || current_time+dt == current_time
            error('HRM7:TimeStep','Invalid time step, current physical time %.9g s.',current_time);
        end

        U = U_0+dt*rhs;
        U(:,1:2) = [U(:,3),U(:,3)];
        U(:,N-1:N) = [U(:,N-2),U(:,N-2)];
        W = U_to_W(U);
        rho1 = W(1,:); u1 = W(2,:); p1 = W(3,:);
        beta1 = W(6,:); alpha1 = W(7,:);
        if any(~isfinite(W(:))) || any(rho1<=0 | p1<=0) ...
                || any(beta1<0 | beta1>1 | alpha1<0 | alpha1>1)
            error('HRM7:State','Predicted pressure, density, or phase fractions are unacceptable.');
        end
        if use_phase_merge
            tinyV = beta1<=phase_tol & alpha1<=phase_tol;
            tinyL = 1-beta1<=phase_tol & 1-alpha1<=phase_tol;
            beta1(tinyV)=0; alpha1(tinyV)=0;
            beta1(tinyL)=1; alpha1(tinyL)=1;
        end
        phase_mass_fix = rho1.*(beta1-W(6,:));
        phase_volume_fix = alpha1-W(7,:);

        [c1,e1,T1,Tv1,rhol1,rhov1,Kl,Kv,Gamma_l,Gamma_v,cpl,cpv,hl,hv,c2l,c2v] = ...
            cal_eq_parameter_NN(p1,rho1,beta1,alpha1,tau_M);
        e2 = (U(3,:)-0.5*rho1.*u1.^2)./rho1;
        ep = rho1.*(e1-e2);
        [xi,E0,gammaS,Pi] = cal_xi_and_E0(c1,e1,T1,p1,rho1);
        W_1 = [rho1;u1;p1;xi;E0;beta1;alpha1];

        U_0 = W_to_U(W_1);
        current_time = current_time+dt;
        steps = steps+1;
        result = [x',rho1',u1',p1',T1',c1',T1',Tv1',beta1',alpha1',rhol1',rhov1'];
        if steps==1 || mod(steps,1000)==0
            disp(['Updating flow field, current time ' num2str(current_time) ', current step ' num2str(steps)]);
        end
        if current_time >= time, flag=0; end
    end
catch ME
    run_error = ME;
    fprintf(2,'Calculation stopped, last completed time step t=%.9g s, steps=%d: %s\n', ...
        current_time,steps,ME.message);
end
toc
%%
rho1 = result(:,2)'; u1 = result(:,3)'; p1 = result(:,4)';
T1 = result(:,5)'; c1 = result(:,6)';
Tl1 = result(:,7)'; Tv1 = result(:,8)';
beta1 = result(:,9)'; alpha1 = result(:,10)';
rhol1 = result(:,11)'; rhov1 = result(:,12)';
T1(beta1==1)=Tv1(beta1==1);
Tl1(beta1==1)=NaN; rhol1(beta1==1)=NaN;
Tv1(beta1==0)=NaN; rhov1(beta1==0)=NaN;
sl1 = cal_entropy(rhol1,Tl1);  sl1(beta1==1)=0;
sv1 = cal_entropy(rhov1,Tv1); sv1(beta1==0)=0;
s1 = sv1.*beta1+sl1.*(1-beta1);

if isempty(run_error)
    disp(['Calculation complete, physical time ' num2str(current_time) 's']);
end


figure(1); clf;
subplot(2,3,1); plot(x,p1,'-r','linewidth',2);
ylabel('$p$','Interpreter','latex');
subplot(2,3,2); plot(x,u1,'-b','linewidth',2);
ylabel('$u$','Interpreter','latex');
subplot(2,3,3); plot(x,Tl1,'-k',x,Tv1,'--r','linewidth',2);
ylabel('$T$','Interpreter','latex');
legend({'$T_l$','$T_v$'},'Interpreter','latex');
subplot(2,3,4); plot(x,rho1,'-m','linewidth',2);
ylabel('$\rho$','Interpreter','latex');
subplot(2,3,5); plot(x,c1,'-b','linewidth',2);
ylabel('$c_W$','Interpreter','latex');
subplot(2,3,6); plot(x,ep,'-b','linewidth',2);
ylabel('$\epsilon$','Interpreter','latex');
% xlim([0.4,0.6])

for j=1:6
    subplot(2,3,j); xlabel('$x$','Interpreter','latex'); grid on;
    set(gca,'FontSize',14,'TickLabelInterpreter','latex');
end
if exist('sgtitle','file')
    sgtitle(sprintf(['RFQC-HRM, RK1, $t=%.6g\\,\\mathrm{s}$, ' ...
        '$N_{\\mathrm{step}}=%d$'],current_time,steps), ...
        'Interpreter','latex');
end
figure(2); clf;
subplot(2,2,1); plot(x,beta1,'-b','linewidth',2);
ylabel('$\beta$','Interpreter','latex');
subplot(2,2,2); plot(x,alpha1,'-r','linewidth',2);
ylabel('$\alpha$','Interpreter','latex');
subplot(2,2,3); plot(x,rhol1,'-b','linewidth',2);
ylabel('$\rho_l$','Interpreter','latex');
subplot(2,2,4); plot(x,rhov1,'-r','linewidth',2);
ylabel('$\rho_v$','Interpreter','latex');
for j=1:4
    subplot(2,2,j); xlabel('$x$','Interpreter','latex'); grid on;
    set(gca,'FontSize',14,'TickLabelInterpreter','latex');
end

  

%% 
