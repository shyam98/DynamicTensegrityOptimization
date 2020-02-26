
function [f_e] = externalforce(D,nnodes,n,pc,dn,eta,f_g,Cc)

% ---------------------- Finding f_n ---------------------------------
    % f_n is a (D*nnodes)*1 dimensional matrix
    % f_n = [f_1x ; f_1y ; ...... ];
    f_n = zeros(D*nnodes,1);
    % Find the normal force f_n (f_iy) on each nodes
    for i=1:nnodes
        if n(2+D*(i-1),1) < 0
            f_n(2+D*(i-1),1) = pc*n(2+D*(i-1),1)^2 + Cc*2*n(2+D*(i-1),1)*dn(2+D*(i-1),1); 
        end
    end
    % Finding friction force f_f (f_ix) on each nodes
    % v_G = [v1_x;0;v2_x;0;... ];
    v_G = dn;
    for i = 1:D*nnodes
        if rem(i,2) == 0
            v_G(i,1) = 0;
        end
    end
    % calculate f_f
    f_f = zeros(D*nnodes,1);
  for i = 1:nnodes     
      if sqrt( v_G(1+D*(i-1))^2 ) < 0.01  % If v_x^2 + v_y^2 < tolerance, we assume f_f = 0 
                
      else     
           f_f(1+D*(i-1),1) = - v_G(1+D*(i-1)) / abs( v_G(1+D*(i-1)) ) * eta*f_n(2+D*(i-1),1) ;
      end
  end
    f_c = f_f + f_n;
    % ---------------------- Finding f_e ---------------------------------
    f_e = f_g + f_c;


end

