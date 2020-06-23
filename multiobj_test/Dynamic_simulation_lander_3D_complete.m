clc; clear; close all;
r_ss = 0.0012606;
r_si = 0.00071109;
r_b = 0.0023568;
p = 5;
q = 4;

L = 0.93218;
RL_ratio = 0.36362;
C_2 = -0.0032284;
z_position = 0.3666;
cyl = 'RCC';

[mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n] = Dynamic_simulation_lander_3D_fn(L, r_ss, r_si, r_b, p, q, RL_ratio, C_2, z_position, cyl);
n_table = [1]';
VarNames = {'Run', 'Max G', 'B Fail Comp', 'B Fail Tens', 'Ss Fail Tens', 'Si Fail Tens'};
T = table( n_table , Max_g_of_different_orientation, sigma_b_c_diff_n, sigma_b_t_diff_n, sigma_ss_diff_n, sigma_si_diff_n, 'VariableNames', VarNames);