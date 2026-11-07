function [flux, g_flux, u_interface] = Numerical_flux_HLLC_Consistent(W_L, W_R)


    rho_L = W_L(1,:); rho_R = W_R(1,:);
    u_L   = W_L(2,:); u_R   = W_R(2,:);
    p_L   = W_L(3,:); p_R   = W_R(3,:);
    
    xi_L = W_L(4,:); xi_R = W_R(4,:);
    E0_L    = W_L(5,:); E0_R    = W_R(5,:); 

    er = 1e-12;
    rho_L = max(rho_L, er); rho_R = max(rho_R, er);

    h_L = (p_L .* (xi_L + 1) + E0_L) ./ rho_L;
    h_R = (p_R .* (xi_R + 1) + E0_R) ./ rho_R;

    rho_et_L = rho_L .* h_L - p_L + 0.5 * rho_L .* u_L.^2;
    rho_et_R = rho_R .* h_R - p_R + 0.5 * rho_R .* u_R.^2;
    
    c_L = sqrt(abs(h_L ./ xi_L));
    c_R = sqrt(abs(h_R ./ xi_R));

    % ----- Toro PVRS + qK with gamma from xi -----
    a_L = c_L; 
    a_R = c_R;

    rho_bar = 0.5 * (rho_L + rho_R);
    a_bar   = 0.5 * (a_L + a_R);
    
    p_pvrs = 0.5*(p_L + p_R) - 0.5*(u_R - u_L) .* (rho_bar .* a_bar);
    
    gamma_L = 1.0 + 1.0 ./ (xi_L + 1e-22);
    gamma_R = 1.0 + 1.0 ./ (xi_R + 1e-22);
    
    
    q_L = ones(size(p_L));
    q_R = ones(size(p_R));
    
    maskL = (p_pvrs > p_L);
    maskR = (p_pvrs > p_R);
    
    q_L(maskL) = sqrt( 1.0 + (gamma_L(maskL)+1)./(2*gamma_L(maskL)) .* (p_pvrs(maskL)./p_L(maskL) - 1.0) );
    q_R(maskR) = sqrt( 1.0 + (gamma_R(maskR)+1)./(2*gamma_R(maskR)) .* (p_pvrs(maskR)./p_R(maskR) - 1.0) );
    
    S_L = u_L - a_L .* q_L;
    S_R = u_R + a_R .* q_R;
   numer = p_R - p_L + rho_L.*u_L.*(S_L-u_L) - rho_R.*u_R.*(S_R-u_R);
    denom = rho_L.*(S_L-u_L) - rho_R.*(S_R-u_R) + er;
    S_star = numer ./ denom;

    % =========================================================================
    % 2. HLLC conservative flux calculation
    % =========================================================================
    len = length(rho_L);
    flux = zeros(3, len);
    
    factor_L = (S_L - u_L) ./ (S_L - S_star + er);
    factor_R = (S_R - u_R) ./ (S_R - S_star + er);
    
    rho_star_L = rho_L .* factor_L;
    rho_star_R = rho_R .* factor_R;
    
    m_dot_HLLC = zeros(1, len);
    rho_star_current = zeros(1, len); 

    for i = 1:len
        if S_L(i) >= 0

            flux(1, i) = rho_L(i) * u_L(i);
            flux(2, i) = rho_L(i) * u_L(i)^2 + p_L(i);
            flux(3, i) = (rho_et_L(i) + p_L(i)) * u_L(i);
            
            m_dot_HLLC(i) = flux(1, i);
            rho_star_current(i) = rho_L(i); 

        elseif S_R(i) <= 0
 
            flux(1, i) = rho_R(i) * u_R(i);
            flux(2, i) = rho_R(i) * u_R(i)^2 + p_R(i);
            flux(3, i) = (rho_et_R(i) + p_R(i)) * u_R(i);
            
            m_dot_HLLC(i) = flux(1, i);
            rho_star_current(i) = rho_R(i);

        elseif S_L(i) < 0 && S_star(i) >= 0

            U_L = [rho_L(i); rho_L(i)*u_L(i); rho_et_L(i)];
            F_L = [rho_L(i)*u_L(i); rho_L(i)*u_L(i)^2 + p_L(i); (rho_et_L(i)+p_L(i))*u_L(i)];

            e_t_L = rho_et_L(i) / rho_L(i);
            p_term = p_L(i) / (rho_L(i)*(S_L(i)-u_L(i)) + er);
            e_t_star = e_t_L + (S_star(i) - u_L(i)) * (S_star(i) + p_term);
            
            U_star_L = rho_star_L(i) * [1; S_star(i); e_t_star];
            flux(:, i) = F_L + S_L(i) * (U_star_L - U_L);
            m_dot_HLLC(i) = flux(1, i); 
            rho_star_current(i) = rho_star_L(i); 

        else

            
            U_R = [rho_R(i); rho_R(i)*u_R(i); rho_et_R(i)];
            F_R = [rho_R(i)*u_R(i); rho_R(i)*u_R(i)^2 + p_R(i); (rho_et_R(i)+p_R(i))*u_R(i)];
            
            e_t_R = rho_et_R(i) / rho_R(i);
            p_term = p_R(i) / (rho_R(i)*(S_R(i)-u_R(i)) + er);
            e_t_star = e_t_R + (S_star(i) - u_R(i)) * (S_star(i) + p_term);
            
            U_star_R = rho_star_R(i) * [1; S_star(i); e_t_star];
            
            flux(:, i) = F_R + S_R(i) * (U_star_R - U_R);
            
            m_dot_HLLC(i) = flux(1, i);
            rho_star_current(i) = rho_star_R(i);
        end
    end

    
    g_flux = zeros(2, len);
    u_interface = zeros(1, len);
    
    Phi_L = W_L(4:5, :);
    Phi_R = W_R(4:5, :);
    
    for i = 1:len
        m = m_dot_HLLC(i);
        
        if m >= 0
            u_consistent = m / rho_L(i);
            g_flux(:, i) = u_consistent * Phi_L(:, i);
        else
            u_consistent = m / rho_R(i);
            g_flux(:, i) = u_consistent * Phi_R(:, i);
        end  
        u_interface(i) = u_consistent;
    end
 end
 

 

