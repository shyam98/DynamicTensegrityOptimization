function [n,N] = nodematrix(N_norotation,height,nnodes)

% nodes position matrix n
   
    theta_0 = rand*2*pi;                 
    %theta_0 = 0/180*pi; 
    R = [  cos(theta_0)     -sin(theta_0)
            sin(theta_0)      cos(theta_0) ];
    
    
    heightmatrix = repmat([0;height],1,nnodes);
    
    N = heightmatrix + R*N_norotation;
    
    n = [];
    
    for i =1:length(N)
       n = [n; N(1,i); N(2,i)];
    end


end

