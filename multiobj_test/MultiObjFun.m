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
cyl = 'SP';

%Call Function
[mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n, volume_const] = Dynamic_simulation_lander_3D_fn(L, r_ss, r_si, r_b, p, q, RL_ratio, C_2, z_position, cyl);

if Max_g_of_different_orientation > 10
    lambda_g = (Max_g_of_different_orientation - 10)^2;
else
    lambda_g = 0;
end
if sigma_ss_diff_n > 0
    lambda_ss = sigma_ss_diff_n^2;
end
if sigma_si_diff_n > 0
    lamba_si = sigma_si_diff_n^2;
end
if sigma_bc_diff_n > 0
    lambda_bc = sigma_b_c_diff_n^2;
end
if sigma_bt_diff_n > 0
    lambda_bt = sigma_b_t_diff_n^2;
end
if volume_const == 1
    lambda_vc = 1000;
end
y = mass + lambda_g + lambda_ss + lamba_si + lambda_bc + lambda_bt + lambda_vc;

end