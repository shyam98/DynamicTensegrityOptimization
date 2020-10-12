clc; clear; close all;
p = 4;
q = 4;
L = 0.99028;
r_ss = 0.00079179;
r_si = 0.00071109;
r_b = 0.0023707;
C_2 = -0.068604;
z_position = 0.37304;
RL_ratio = 0.28606;
cyl = 'SP';

[mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n, ...
    sigma_b_c_diff_n, sigma_b_t_diff_n, volume_const, node_flip,...
    number_of_loop, n_bottom_min_track, n_bottom_max_track] ...
    = Dynamic_simulation_lander_3D_fn(L, r_ss, r_si, r_b, p, q, RL_ratio, C_2, z_position, cyl);
n_table = [1]';
VarNames = {'Run', 'Max G', 'B Fail Comp', 'B Fail Tens', 'Ss Fail Tens', 'Si Fail Tens', 'Flipped?'};
T = table( n_table , Max_g_of_different_orientation, sigma_b_c_diff_n, sigma_b_t_diff_n, sigma_ss_diff_n, sigma_si_diff_n, node_flip, 'VariableNames', VarNames);


%{
GA parameters
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

RCC: 5	4	0.93218	0.0012606	0.00071109	0.0023568	-0.0032284	0.3666	0.36362	4.7653
SP: 4	4	0.99028	0.00079179	0.00072517	0.0023707	-0.068604	0.37304	0.28606	4.6628

%}
