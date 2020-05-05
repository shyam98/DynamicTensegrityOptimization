function [n,N,theta_0x,theta_0y,theta_0z] = nodematrix(N_norotation,height,nnodes, theta_0x, theta_0y, theta_0z, L)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% nodes position matrix n
%    
%     theta_0x = rand*2*pi ;     
%     theta_0y = rand*2*pi ;
%     theta_0z = rand*2*pi ;

    
    
    Rx = [ 1       0                      0
           0     cos(theta_0x)     -sin(theta_0x)
           0     sin(theta_0x)      cos(theta_0x) ];
    
    Ry = [ cos(theta_0y)       0      sin(theta_0y)
           0                   1            0
          -sin(theta_0y)       0      cos(theta_0y) ];
        
    Rz = [ cos(theta_0z)    -sin(theta_0z)     0
           sin(theta_0z)     cos(theta_0z)     0
           0                      0            1 ];
    
    heightmatrix = repmat([0;height+(L/2);0],1,nnodes);
    N = heightmatrix + Rx*Ry*Rz*N_norotation;
    
    n = [];
    
     for i =1:length(N)
       n = [n; N(1,i); N(2,i);N(3,i)];
     end
end

