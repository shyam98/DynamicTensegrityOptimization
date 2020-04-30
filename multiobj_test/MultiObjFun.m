function y = MultiObjFun(x)
%Structure 
p = x(1);
q = x(2);

%Radii
r = x(3);
r_ss = x(4);
r_si = x(5);
r_b = x(6);

%Internal Parameters
C_2 = x(7);
z_position = x(8);

%MISC
RL_ratio = 1;
cyl = 'SP';

p
q
%Call Function
[mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n, volume_const] = Dynamic_simulation_lander_3D_fn(r, r_ss, r_si, r_b, p, q, RL_ratio, C_2, z_position, cyl);

f_Max_g = 0;
f_sigma_failure = 0;

%Check if node is out of bounds
if volume_const == 1
    f_Max_g = 2;
    f_sigma_failure = 2;
else
    if Max_g_of_different_orientation > 15 %Max-G constraint
        f_Max_g = 2;
    end
    if sigma_ss_diff_n ~= 0 || sigma_si_diff_n ~= 0 || ...
        sigma_b_c_diff_n~=0  || sigma_b_t_diff_n ~=0 %Stress Constraint
        f_sigma_failure = 2;
    end
end


y = 100*mass + Max_g_of_different_orientation^f_Max_g + ...
    sigma_ss_diff_n^f_sigma_failure + sigma_si_diff_n^f_sigma_failure ...
    + sigma_b_c_diff_n^f_sigma_failure + ...
    sigma_b_t_diff_n^f_sigma_failure;

end