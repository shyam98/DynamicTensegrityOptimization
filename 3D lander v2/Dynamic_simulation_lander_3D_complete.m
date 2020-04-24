clc; clear;
r_ss = 0.01;
r_si = 0.005;
r_b = 0.02;
p = 3;
q = 4;

r = 1.5;
RL_ratio = 5

[mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n] = Dynamic_simulation_lander_3D_fn(r, r_ss, r_si, r_b, p, q, RL_ratio);
n_table = [1]';
VarNames = {'Run', 'Max G', 'Max B Stress', 'Min B Stress', 'Max Si Stress', 'Min Si Stress', 'Max Ss Stress', 'Min Ss Stress', 'B Fail Comp', 'B Fail Tens', 'Ss Fail Tens', 'Si Fail Tens'};
T = table( n_table , Max_g_of_different_orientation, sigma_bar_max_n, sigma_bar_min_n, sigma_si_max_n, sigma_si_min_n, sigma_ss_max_n, sigma_ss_min_n, sigma_b_c_diff_n, sigma_b_t_diff_n, sigma_ss_diff_n, sigma_si_diff_n, 'VariableNames', VarNames)