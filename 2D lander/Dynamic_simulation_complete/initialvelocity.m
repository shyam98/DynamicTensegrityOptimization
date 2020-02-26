function [dn] = initialvelocity(n,height,v_0,dtheta_0,nnodes)

    % Find out the center of the lander
    center = repmat([0;height],1,nnodes);
    % Calculate r matrix. It is a 3*nnodes matrix
    r = n - center;
    % Define initial translation velocity matrix v_0
    V_initialtranslation = repmat([0;v_0],1,nnodes);
    % Total initial velocity
    dN = V_initialtranslation - r*dtheta_0;
    
    dn = [];
    
    for i =1:length(dN)
       dn = [dn; dN(1,i); dN(2,i)];
    end


end

