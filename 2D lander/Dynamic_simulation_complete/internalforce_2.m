function [f_I,varepsilon_s,sigma_s,varepsilon_b,sigma_b, ... 
    s_initiallength,b_initiallength, sigma_ss, sigma_si, ...  
    sigma_ss_diff, sigma_si_diff, sigma_b_c_diff, sigma_b_t_diff] = internalforce_2(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b, c_s,c_b,A_s,A_b,Yield_Nylon, Youngs_Titanium, Yield_Titanium)

 %  ----------------------- Strings force ------------------------------
    % Finding initial length matrix and current length matrix (both are n_s*1 matrix)
    
    s_nowlength = [];  
    s_initiallength = [];
    for i = 1:n_s
        s_initiallength = [s_initiallength; sqrt( s_0(1+D*(i-1))^2 + s_0(2+D*(i-1))^2 ) ];
        s_nowlength = [s_nowlength ; sqrt( s(1+D*(i-1))^2 + s(2+D*(i-1))^2 ) ]  ;
    end
       
    % d\varepsilon_s and \varepsilon_s (both are n_s*1 matrix)
    d_varepsilon_s = [];
    varepsilon_s = [];
    sT = s';
    for i = 1:n_s 
        d_varepsilon_s = [d_varepsilon_s ; 
        sT(1,1+D*(i-1):2+D*(i-1)) * ds(1+D*(i-1):2+D*(i-1),1) / ( s_initiallength(i) * s_nowlength(i) ) ];
        varepsilon_s = [varepsilon_s ; ( s_nowlength(i) - s_initiallength(i)) / s_initiallength(i) ];
    end
    % sigma_s matrix ( n_s*1 )
    sigma_s = [];
    for i = 1:n_s
        sigma_s = [sigma_s; E_s*varepsilon_s(i) + c_s*d_varepsilon_s(i) ];
    end

    % f_s matrix ( D*n_s * 1 )
    f_s = [];
    %f_s = [external ; internal]
    for i = 1:D*n_s/2
        j = ceil(i/D);
        f_s = [ f_s ; sigma_s(j) * A_s(1) * s(i) / s_nowlength(j) ];
    end
    
    for i = (D*n_s/2 + 1):D*n_s
        j = ceil(i/D);
        f_s = [ f_s ; sigma_s(j) * A_s(2) * s(i) / s_nowlength(j) ];
    end
    
    %  ----------------------- String Stress ------------------------------
   % Calculate the stresses in each of the strings and compare to max
   % For strings, only care if the values are positive. 
   % Negative value does not matter because string will not break in
   % compression
   
   sigma_ss = sigma_s(1:(n_s/2),:);
   sigma_si = sigma_s((n_s/2 + 1):n_s,:);
   
   sigma_s_diff = sigma_s - Yield_Nylon;
   sigma_s_check = sigma_s_diff > 0;
   %Get rid of all values of simga_s that are in range
   sigma_s_diff = sigma_s_diff.*sigma_s_check;
   sigma_ss_diff = sigma_s_diff(1:(n_s/2),:);
   sigma_si_diff = sigma_s_diff((n_s/2 + 1):n_s,:);   
    
    %  ------------------------ Bars force --------------------------------
   % Finding initial length matrix and current length matrix (both are n_b*1 matrix)
   
    b_nowlength = [];  
    b_initiallength = [];
    for i = 1:n_b
        b_initiallength = [b_initiallength; sqrt( b_0(1+D*(i-1))^2 + b_0(2+D*(i-1))^2 ) ];
        b_nowlength = [b_nowlength ; sqrt( b(1+D*(i-1))^2 + b(2+D*(i-1))^2 ) ];     
    end
  
    % d\varepsilon_b and \varepsilon_b (both are n_b*1 matrix)
    d_varepsilon_b = [];
    varepsilon_b = [];
    bT = b';
    for i = 1:n_b 
        d_varepsilon_b = [d_varepsilon_b ; 
        bT(1,1+D*(i-1):2+D*(i-1)) * db(1+D*(i-1):2+D*(i-1),1) / ( b_initiallength(i) * b_nowlength(i) ) ];
        varepsilon_b = [varepsilon_b ; ( b_nowlength(i) - b_initiallength(i) )/b_initiallength(i) ];
    end
    % sigma_b matrix ( n_b*1 )
    sigma_b = [];
    for i = 1:n_b
        sigma_b = [sigma_b; E_b*varepsilon_b(i) + c_b*d_varepsilon_b(i) ];
    end
    
    % f_b matrix ( D*n_b * 1 )
    
    f_b = [];
    for i = 1:D*n_b
        j = ceil(i/D);
        f_b = [ f_b ; sigma_b(j) * A_b * b(i) / b_nowlength(j) ];
    end
    
   %  ------------------------- Bar Stress --------------------------------
   % Compressive Stresses to prevent buckling
   r_bars = sqrt(A_b/pi);
   sigma_b_max_c = ((pi^2)*Youngs_Titanium*(r_bars^2))./((b_nowlength).^2);
   compressive_sigma_b = sigma_b > 0; %Get rid of all tensile values
   sigma_b_compressive = sigma_b;
   sigma_b_compressive(compressive_sigma_b) = 0;
   sigma_b_max_c(compressive_sigma_b) = 0;
   %Positive sigma_b_comp_diff means failure
   sigma_b_c_diff = abs(sigma_b_compressive) - sigma_b_max_c;
   sigma_b_c_check = sigma_b_c_diff > 0;
   sigma_b_c_diff = sigma_b_c_diff .* sigma_b_c_check;
   
   % Tensile Stresses to prevent yielding
   
   
   
   tensile_sigma_b = sigma_b > 0; % Get rid of compressive values
   tensile_sigma_b = sigma_b.*tensile_sigma_b;
   sigma_b_t_diff = tensile_sigma_b - Yield_Titanium;
   sigma_b_t_check = sigma_b_t_diff > 0;
   sigma_b_t_diff = sigma_b_t_diff.*sigma_b_t_check;
   
    %--------------------------  f_I --------------------------------------
    f_I = kron(C_sT,I_D)*f_s + kron(C_bT,I_D)*f_b ;

end

