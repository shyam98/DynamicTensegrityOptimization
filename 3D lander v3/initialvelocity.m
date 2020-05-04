function [dn] = initialvelocity(N,height,v_0,nnodes)
%UNTITLED5 Summary of this function goes here
    v_y = v_0;
%   Detailed explanation goes here
    dtheta = rand(3,1)*2-1;
    dtheta_x0 = dtheta(1) / norm(dtheta);
    dtheta_y0 = dtheta(2) / norm(dtheta);
    dtheta_z0 = dtheta(3) / norm(dtheta);

 % Find out the center of the lander
    center = repmat([0;height;0],1,nnodes);
    % Calculate r matrix. It is a 3*nnodes matrix
    r = N - center;
    % Define initial translation velocity matrix v_0
    V_initialtranslation = repmat([0;v_y;0],1,nnodes);
    % Rotation matrix
    theta = repmat([0;0;0],1,nnodes);
    % Total initial velocity
    dN = V_initialtranslation;
    
    dn = [];
    
    for i =1:length(dN)
       dn = [dn; dN(1,i); dN(2,i);dN(3,i)];
    end
    
    
end

