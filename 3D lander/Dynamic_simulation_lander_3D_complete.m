clc; clear;

%User input for testing before GA
r_ss = input('Input surface string radii, default 0.01 -- r_ss =  ');
if isempty(r_ss)
 r_ss = 0.01;
end

r_si = input('Input inner string radii, default 0.05 -- r_si = ');
if isempty(r_si)
 r_si = 0.005;
end

r_b = input('Input bar radii, default 0.02 -- r_b = ');
if isempty(r_b)
 r_b = 0.02;
end

p = input('Input bar configuration, default 3 -- p = ');
if isempty(p)
 p = 6;
end

q = input('Input string configuration, default 4 -- q = ');
if isempty(q)
 q = 4;
end

r = input('Input lattice radius, default 1.5 -- r = ');
if isempty(r)
 r = 1.5;
end

number_of_orientation = input('Input Number of times to run the code: ');
if isempty(number_of_orientation)
 number_of_orientation = 1;
end

[mass, Max_g_of_different_orientation, sigma_ss_max_n, sigma_ss_min_n, sigma_si_max_n, sigma_si_min_n, sigma_bar_max_n, sigma_bar_min_n, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n] = Dynamic_simulation_lander_3D_fn(r, r_ss, r_si, r_b, p, q, number_of_orientation);
n_table = [1:number_of_orientation]';
VarNames = {'Run', 'Max G', 'Max B Stress', 'Min B Stress', 'Max Si Stress', 'Min Si Stress', 'Max Ss Stress', 'Min Ss Stress', 'B Fail Comp', 'B Fail Tens', 'Ss Fail Tens', 'Si Fail Tens'};
T = table( n_table , Max_g_of_different_orientation, sigma_bar_max_n, sigma_bar_min_n, sigma_si_max_n, sigma_si_min_n, sigma_ss_max_n, sigma_ss_min_n, sigma_b_c_diff_n, sigma_b_t_diff_n, sigma_ss_diff_n, sigma_si_diff_n, 'VariableNames', VarNames)
