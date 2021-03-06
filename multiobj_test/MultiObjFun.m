function y = MultiObjFun(x)
%Structure 
p = x(1);
q = x(2);

%Radii
L = x(3);
r_ss = x(4);
r_si = x(5);
r_b = x(6);

%Internal Parameters
C_2 = x(7);
z_position = x(8);

%MISC
RL_ratio = x(9);
cyl = 'RCC';

%Call Function
[mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n,...
    sigma_b_c_diff_n, sigma_b_t_diff_n, volume_const,...
    node_flip] = Dynamic_simulation_lander_3D_fn(L, r_ss, r_si, r_b, p, q, RL_ratio, C_2, z_position, cyl);

if Max_g_of_different_orientation > 10
    lambda_g = ((Max_g_of_different_orientation - 10)^2)*1E12;
else
    lambda_g = 0;
end

if node_flip == 1
    lambda_node_flip = 1000;
else
    lambda_node_flip = 0;
end

lambda_ss = sigma_ss_diff_n^2;
lamba_si = sigma_si_diff_n^2;
lambda_bc = sigma_b_c_diff_n^2;
lambda_bt = sigma_b_t_diff_n^2;
lambda_vc = (volume_const^2)*1E14;
y = mass + lambda_g + lambda_ss + lamba_si + lambda_bc + lambda_bt + ...
    lambda_vc + lambda_node_flip;
fileID = fopen('GA_data.txt','a+');
fprintf(fileID, '%d %d %f %f %f %f %f %f %f %f\n',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9), y);


end