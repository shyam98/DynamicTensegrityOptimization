function [f_e] = externalforce(D,nnodes,n,pc,dn,eta,f_g,Cc)
    %UNTITLED6 Summary of this function goes here
    %   Detailed explanation goes here

    % ---------------------------- Finding f_n ----------------------------
        % f_n is a (D*nnodes)*1 dimensional matrix
        % f_n = [f_1x ; f_1y ; ...... ];
        f_n = zeros(D*nnodes,1);
        % Find the normal force f_n (f_iy) on each nodes
        for i=1:nnodes
            if n(2+D*(i-1),1) < 0
                f_n(2+D*(i-1),1) = pc*n(2+D*(i-1),1)^2 + Cc*2*n(2+D*(i-1),1)*dn(2+D*(i-1),1); 
            end
        end
    %--------- Finding friction force f_f (f_ix, f_iz) on each nodes ------
        % calculate f_f
        f_f = zeros(D*nnodes,1);
      for i = 1:nnodes     
          if sqrt( dn(1+D*(i-1))^2 + dn(3+D*(i-1))^2 ) > 0.01  % If v_x^2 + v_z^2 < tolerance, we assume f_f = 0              
               f_f(1+D*(i-1),1) = - dn(1+D*(i-1)) / sqrt( dn(1+D*(i-1))^2 + dn(3+D*(i-1))^2 ) * eta*f_n(2+D*(i-1),1) ;
               f_f(3+D*(i-1),1) = - dn(3+D*(i-1)) / sqrt( dn(1+D*(i-1))^2 + dn(3+D*(i-1))^2 ) * eta*f_n(2+D*(i-1),1) ;
          end
      end
    f_c = f_f + f_n;
    % ---------------------- Finding f_e ---------------------------------
    f_e = f_g + f_c;
end

