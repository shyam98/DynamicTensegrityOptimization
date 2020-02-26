function [f_I,varepsilon_s,sigma_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength, f_s] = internalforce(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b,c_s,c_b,A_s,A_b)

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
    for i = 1:D*n_s
        j = ceil(i/D);
        f_s = [ f_s ; sigma_s(j) * A_s * s(i) / s_nowlength(j) ];
    end
    
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
   
    %--------------------------  f_I --------------------------------------
    f_I = kron(C_sT,I_D)*f_s + kron(C_bT,I_D)*f_b ;



end

