function [] = Tensegrity_simulation(r_ss, r_si, r_b, p, q, number_of_orientation, C_dist, R_range, RL_ratio, m_load)
tic
%% Input parameters
% ====================== Lander Initial Properties ========================
% D = Dimension; I_D = Identity Matrix; height = Initial height
% r = radius of the lander; m_load = mass of central node
% v_0 = Initial vertical speed

%Calculate Areas from Radius of input
A_ss = (r_ss)^2 * pi;
A_si = (r_si)^2 * pi;
A_b = (r_b)^2 * pi;
D = 2;
I_D = eye(D);
height = 10;
v_0 = -5;
A_s = [A_ss ; A_si];
dtheta_0 = 1*pi/180;

Yield_Nylon = 5e7;
Yield_Titanium = 1e9;
Youngs_Titanium = 115e9;


% ======================= Environment Properties ==========================
% pc = ground bounce constant; Cc = ground damping coefficient
% eta = ground frictional coefficient; % g_mars = Gravitational
% acceleration on Mars
pc = 5e05;  
Cc = 1e04;
eta = 0.5; 
g_earth = -9.81;

% ========================= Simulation Parameters =========================
% dt = timestep size ; total_time = totoal simulation time
% number_of_loop = the number of loops per simulation
% V_tol = Threshold value of average speed of all nodes
% Time_stop: If the average speed of all nodes remains below "V_tol" for
% "Time_stop" loops, we assume the lander stops. In this case, the lander
% is assume to stop when the average speed of all nodes remains below 
% 0.05 m/s for 2500 loops (2 sec).
% accerelation_tol = Maximum allowed value of the acceleration on the
% central node. If a > accerelation_tol, we believe that this configuration
% is not safe.
% Final_distance
dt = 4*10^(-4);
total_time = 10;
number_of_loop = total_time/dt;
V_tol = 0.05;
Time_stop = 2/dt;
acceleration_tol = 100*9.8;
toc
end

