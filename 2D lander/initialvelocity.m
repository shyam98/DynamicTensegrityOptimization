function [dn] = initialvelocity(n,height,v_0,dtheta_0,nnodes)
    
    %Randomize the velocities between x and y 
    % Max x dir velocity is between 0% and 5% of v_0
    v_x = v_0 * 0.05 * rand();
    v_y = sqrt((v_0^2)-(v_x^2));
    
    % Find out the center of the lander
    center = repmat([0;height],1,nnodes);
    % Calculate r matrix. It is a 3*nnodes matrix
    r = n - center;
    % Define initial translation velocity matrix v_0
    V_initialtranslation = repmat([v_x;-v_y],1,nnodes);
    % Total initial velocity
    dN = V_initialtranslation - r*dtheta_0;
    
    dn = [];
    
    for i =1:length(dN)
       dn = [dn; dN(1,i); dN(2,i)];
    end


end

