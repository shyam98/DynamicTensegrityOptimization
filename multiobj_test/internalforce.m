function [f_I,varepsilon_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength, sigma_ss, sigma_si, sigma_ss_diff, sigma_si_diff, sigma_b_c_diff, sigma_b_t_diff] = internalforce(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b, c_s,c_b,A_s,A_b,Yield_Nylon, Youngs_Titanium, Yield_Titanium, n_ss, FoS)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

 %  ----------------------- Factor of Saftey ------------------------------
String_Yield = Yield_Nylon / FoS;
Bar_Yield = Yield_Titanium / FoS;
 %  ----------------------- Strings force ------------------------------
    % Finding initial length matrix and current length matrix (both are n_s*1 matrix)
    s_initiallength = [];
    s_nowlength = [];  
    for i = 1:n_s
        s_initiallength = [s_initiallength; sqrt( s_0(1+D*(i-1))^2 + s_0(2+D*(i-1))^2 + s_0(3+D*(i-1))^2  ) ] ;
        s_nowlength = [s_nowlength ; sqrt( s(1+D*(i-1))^2 + s(2+D*(i-1))^2 + s(3+D*(i-1))^2 ) ]  ;
    end
       
    % d\varepsilon_s and \varepsilon_s (both are n_s*1 matrix)
    d_varepsilon_s = [];
    varepsilon_s = [];
    sT = s';
    for i = 1:n_s 
        d_varepsilon_s = [d_varepsilon_s ; 
        sT(1,1+D*(i-1):3+D*(i-1)) * ds(1+D*(i-1):3+D*(i-1),1) / ( s_initiallength(i) * s_nowlength(i) ) ];
        varepsilon_s = [varepsilon_s ; ( s_nowlength(i) - s_initiallength(i)) / s_initiallength(i) ];
    end
    % sigma_s matrix ( n_s*1 )
    sigma_s = [];
    for i = 1:n_s
        string_buckle_check = E_s*varepsilon_s(i) + c_s*d_varepsilon_s(i);
        if string_buckle_check < 0 
            sigma_s = [sigma_s; 0]; 
        else
            sigma_s = [sigma_s; string_buckle_check];
        end
    end
    % Separate sigma into surface and internal stresses
    sigma_ss = sigma_s(1:(n_ss),:);
    sigma_si = sigma_s((n_ss + 1):n_s,:);    
    
   
    % f_s matrix ( D*n_s * 1 )
    %f_s = [external ; internal]
    f_s = [];
    for i = 1:D*n_ss
        j = ceil(i/D);
        f_s = [ f_s ; sigma_s(j) * A_s(1) * s(i) / s_nowlength(j) ];
    end
    for i = (D*n_ss + 1):D*n_s
        j = ceil(i/D);
        f_s = [ f_s ; sigma_s(j) * A_s(2) * s(i) / s_nowlength(j) ];
    end
    
    %  ------------------------ String Stress --------------------------------
   %Calculate the difference between string stress and yield stress
   sigma_s_diff = sigma_s - String_Yield;
   sigma_s_check = sigma_s_diff > 0;
   %Get rid of all values of simga_s that are in range
   sigma_s_diff = sigma_s_diff.*sigma_s_check;
   sigma_ss_diff = sigma_s_diff(1:(n_ss),:);
   sigma_si_diff = sigma_s_diff((n_ss + 1):n_s,:); 
    
    %  ------------------------ Bars force --------------------------------
   % Finding initial length matrix and current length matrix (both are n_b*1 matrix)
    b_initiallength = [];
    b_nowlength = [];  
    for i = 1:n_b
        b_initiallength = [b_initiallength; sqrt( b_0(1+D*(i-1))^2 +  b_0(2+D*(i-1))^2 + b_0(3+D*(i-1))^2 ) ];  % Only meets 3D system!
        b_nowlength = [b_nowlength ; sqrt( b(1+D*(i-1))^2 + b(2+D*(i-1))^2 + b(3+D*(i-1))^2 ) ];     
    end
  
    % d\varepsilon_b and \varepsilon_b (both are n_b*1 matrix)
    d_varepsilon_b = [];
    varepsilon_b = [];
    bT = b';
    for i = 1:n_b 
        d_varepsilon_b = [d_varepsilon_b ; 
        bT(1,1+D*(i-1):3+D*(i-1)) * db(1+D*(i-1):3+D*(i-1),1) / ( b_initiallength(i) * b_nowlength(i) ) ];
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
    
    %  ------------------------ Bars Stress --------------------------------
   %Compressive stresses to buckling
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
   
    %Tensile stresses for yield
   tensile_sigma_b = sigma_b > 0; % Get rid of compressive values
   tensile_sigma_b = sigma_b.*tensile_sigma_b;
   sigma_b_t_diff = tensile_sigma_b - Bar_Yield;
   sigma_b_t_check = sigma_b_t_diff > 0;
   sigma_b_t_diff = sigma_b_t_diff.*sigma_b_t_check;
 
   
    %--------------------------  f_I --------------------------------------
    f_I = kron(C_sT,I_D)*f_s + kron(C_bT,I_D)*f_b ;
end

