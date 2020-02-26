function [dn] = initialvelocity(N,height,v_0,dtheta_max,nnodes)
%UNTITLED5 Summary of this function goes here
    % Separate the v_0 into random components of x y z
    v_x = rand() * 0.05 * v_0;
    v_z = rand() * 0.05 * v_0;
    v_y = sqrt((v_0^2)-(v_x^2)-(v_z^2));
%   Detailed explanation goes here
    
    dtheta = rand(3,1)*2-1;
    dtheta_x0 = dtheta_max * dtheta(1) / norm(dtheta);
    dtheta_y0 = dtheta_max * dtheta(2) / norm(dtheta);
    dtheta_z0 = dtheta_max * dtheta(3) / norm(dtheta);

 % Find out the center of the lander
    center = repmat([0;height;0],1,nnodes);
    % Calculate r matrix. It is a 3*nnodes matrix
    r = N - center;
    % Define initial translation velocity matrix v_0
    V_initialtranslation = repmat([v_x;-v_y;v_z],1,nnodes);
    % Rotation matrix
    theta = repmat([dtheta_x0;dtheta_y0;dtheta_z0],1,nnodes);
    % Total initial velocity
    dN = V_initialtranslation - cross(r,theta);
    
    dn = [];
    
    for i =1:length(dN)
       dn = [dn; dN(1,i); dN(2,i);dN(3,i)];
    end
    
    
end

