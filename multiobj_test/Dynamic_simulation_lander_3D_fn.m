function [mass, Max_g_of_different_orientation, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n, volume_const] = Dynamic_simulation_lander_3D_fn(L, r_ss, r_si, r_b, p, q, RL_Ratio, C_2, z_position, cyl)

tic
%% Input parameters
% ====================== Lander Initial Properties ========================
% D = Dimension; I_D = Identity Matrix; height = Initial height
% r = radius of the lander; m_load = mass of central node
% m_s_total = Total mass of strings ; m_b_total = Total mass of bars
% v_0 = Initial vertical speed

%Calculate Areas from Radius Input
A_ss = (r_ss)^2 * pi;
A_si = (r_si)^2 * pi;
A_b = (r_b)^2 * pi;
A_s = [A_ss ; A_si];

node_mass = 0.010;

D = 3;
I_D = eye(D);
height = 1.01;
v_0 = -10;
m_load = 20;

dtheta_max = 1*pi/180;

Yield_Nylon = 9.4e7;
Yield_Titanium = 1e9;
Youngs_Titanium = 115e9;

% ======================= Environment Properties ==========================
% pc = ground bounce constant; Cc = ground damping coefficient
% eta = ground frictional coefficient; % g_earth = Gravitational
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
dt = (0.1*10^(-4))/0.8;
total_time = 0.8;
number_of_loop = total_time/dt;
number_of_loop = ceil(number_of_loop)
V_tol = 0.05;
Time_stop = 2/dt;
acceleration_tol = 5000*9.8;

%% ======================== Materials properties ===========================
% --------------------- bars -----------------------
rho_b = 5000;
E_b = 60e09;
c_b = 0;
% ------------------- strings ----------------------
rho_s = 2500;
E_s = 1e08;
c_s = 5e06;

%% ======================== Display parameters ============================
%number_of_configurations = 3;
%X=round(linspace(1,number_of_loop,number_of_configurations));
%----- Initializing arbitrary node position matrices to be displayed ------
%display_node = [];
%display_node_x_position = zeros(length(display_node),number_of_loop);
%display_node_y_position = zeros(length(display_node),number_of_loop);


%%  ============================================================================================
%   Start dynamic simulation
%   ============================================================================================
%% ==================== node matrix and C_s, C_b ======================
% 'RCC': right circular cylinder
% 'SP': sphere

%% Angle Adjustments
%Cylinders
%y rotation goes from 0 - pi
%x rotation goes from 0 - .222*pi (20d)
theta_0y = 0;
theta_0x = pi/6;
%Spheres
%y rotation goes from 0 - pi
%x rotation goes from 0 - pi
theta_0z = 0;

%Create 3D Lander
[N_norotation,C_b,C_s,nnodes,n_s,n_b, n_ss, V_c] = Lander_3D(q,p,L,cyl, C_2, z_position, RL_Ratio);
C_sT = C_s'; C_bT = C_b';

if V_c ~= 0
    mass = 1000;
    Max_g_of_different_orientation = 100;
    sigma_ss_diff_n = 1e6;
    sigma_si_diff_n = 1e6;
    sigma_b_c_diff_n = 1e6;
    sigma_b_t_diff_n = 1e6;
    volume_const = V_c;
    return
else
    volume_const = 0;
end
tenseg_plot(N_norotation, C_b, C_s);
axis on
%% Initializing Matrices or arrays
% ---------------------- Final distance arrays ----------------------------
Max_g_of_different_orientation = zeros(length(theta_0y), length(theta_0x)); % Store Maximum acceleration
deviation_distance = [];             % Store deviation distance
Falling_time = [];                   % Store stop time of all orientations
deepest_y = [];                      % Store the deepest penetration of all orientations
theta = [];                          % Record theta of one orientation

% Keep track of these parameters
sigma_ss_max_n = zeros(length(theta_0y), length(theta_0x));
sigma_ss_min_n = zeros(length(theta_0y), length(theta_0x));
sigma_si_max_n = zeros(length(theta_0y), length(theta_0x));
sigma_si_min_n = zeros(length(theta_0y), length(theta_0x));
sigma_bar_max_n = zeros(length(theta_0y), length(theta_0x));
sigma_bar_min_n = zeros(length(theta_0y), length(theta_0x));
sigma_si_diff_n = zeros(length(theta_0y), length(theta_0x));
sigma_ss_diff_n = zeros(length(theta_0y), length(theta_0x));
sigma_b_c_diff_n = zeros(length(theta_0y), length(theta_0x));
sigma_b_t_diff_n = zeros(length(theta_0y), length(theta_0x));

%% Iterate about the y-axis
for thetay_i = 1:length(theta_0y)
    deepest_y_of_each_orientation = 0;
    %% Iterate about the x-axis
    for thetax_i = 1:length(theta_0x)
        %% Node position matrix n
        [n,N] = nodematrix(N_norotation,height,nnodes,theta_0x(thetax_i),theta_0y(thetay_i),theta_0z, L);
        n_0 = n;
        N_0 = N;
        %% Initial velocity dn
        dn = initialvelocity(N,height,v_0,dtheta_max,nnodes);
        dn_0 = dn;
        %% Define strings matrix s and bars matrix b
        %  s and b indicates the length (2D - x and y direction) of each string and each bar
        s = kron(C_s,I_D)*n;
        b = kron(C_b,I_D)*n;
        ds = kron(C_s,I_D)*dn;
        db = kron(C_b,I_D)*dn;
        s_0 = s;
        b_0 = b;
        ds_0 = ds;
        db_0 = db;
        %% Mass (Only calculate once)
        if thetay_i == 1
            [M,m] = Mass(D,I_D,n_s,n_b,s,b,C_sT,C_bT,rho_s,rho_b,m_load, A_s, A_b, node_mass, n_ss);
            invM = inv(M);
        end
        %%  Finding f_g (Only calculate once)
        if thetay_i == 1
            g = repmat([0;g_earth;0],nnodes,1);
            f_g = M*g;
        end
        %% Initializing Energy matrices and velocity scaler
        U = [];
        E_k = [];
        U_po = [];
        V_scaler = [];
        center_node_g = zeros(number_of_loop,1);
        currenttime = [];
        number_of_configurations = 10;
        X=round(linspace(1,number_of_loop,number_of_configurations));
        loop_after_V_tol = 0;
        acceleration_excess = 0;
        %% Initializing arbitrary node position matrices to be displayed
        display_node = [1,6];
        display_node_x_position = zeros(length(display_node),number_of_loop);
        display_node_y_position = zeros(length(display_node),number_of_loop);
        display_node_z_position = zeros(length(display_node),number_of_loop);
        
        %% Initializing Stress matricies for tracking in each run
        sigma_ss_max = zeros(number_of_loop,1);
        sigma_ss_min = zeros(number_of_loop,1);
        
        sigma_si_max = zeros(number_of_loop,1);
        sigma_si_min = zeros(number_of_loop,1);
        
        sigma_bar_max = zeros(number_of_loop,1);
        sigma_bar_min = zeros(number_of_loop,1);
        
        sigma_ss_diff = zeros(number_of_loop,1);
        sigma_si_diff = zeros(number_of_loop,1);
        sigma_b_c_diff = zeros(number_of_loop,1);
        sigma_b_t_diff = zeros(number_of_loop,1);
        %% Initializing maximum and minimum value of x,y,z position
        % They are the maximum of minimum valve for every loops
        % They are used to confine the range of axis if we plot diagrams.
        Xmax = []; Xmin = [];
        Ymax = []; Ymin = [];
        Zmax = []; Zmin = [];
        %% Initialize the number of output images
        ii = 1;
        %% ====== Start simulation for each orientation =========
        for loop = 1:number_of_loop
            if loop == number_of_loop/10
                tenseg_plot(N,C_s,C_b)
            end
            
            
            
            %% Finding trajectory matrices of arbitrary nodes
            
            for j = 1:length(display_node)
                display_node_x_position(j,loop) = N(1,display_node(j));
                display_node_y_position(j,loop) = N(2,display_node(j));
                display_node_z_position(j,loop) = N(3,display_node(j));
            end
            %% Save the node position matrix every 0.1s to make gif
            
            
            % Finding the maximum and minimum value of x position
            Xmax = [Xmax max(N(1,:))]; Xmin = [Xmin min(N(1,:))];
            Ymax = [Ymax max(N(2,:))]; Ymin = [Ymin min(N(2,:))];
            Zmax = [Zmax max(N(3,:))]; Zmin = [Zmin min(N(3,:))];
            
            
            if rem(loop*dt,0.01) == 0 || loop == 1
                cell{ii} = N;
                ii = ii + 1;
            end
            
            
            %% Plotting an arbitrary number of lander configurations throughout the time frame
            
            figure(2);
            if ismember(loop,X) == 1
                tenseg_plot(N,C_b,C_s);
                hold on
                axis off
                ttext = string(find(X==loop));
                %text(N(1,1),N(2,1),ttext,'fontsize',15);
                hold on
            end
            
            
            %% External force
            f_e = externalforce(D,nnodes,n,pc,dn,eta,f_g,Cc);
            %% Internal force
            [f_I,varepsilon_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength, sigma_ss_dt, sigma_si_dt, sigma_ss_diff_dt, sigma_si_diff_dt, sigma_b_c_diff_dt, sigma_b_t_diff_dt] = internalforce(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b, c_s,c_b,A_s,A_b,Yield_Nylon, Youngs_Titanium, Yield_Titanium, n_ss);
            %% Keep track of stresses
            sigma_ss_max(loop) = max(sigma_ss_dt);
            sigma_ss_min(loop) = min(sigma_ss_dt);
            
            sigma_si_max(loop) = max(sigma_si_dt);
            sigma_si_min(loop) = min(sigma_si_dt);
            
            sigma_bar_max(loop) = max(sigma_b);
            sigma_bar_min(loop) = min(sigma_b);
            
            sigma_ss_diff(loop) = max(abs(sigma_ss_diff_dt));
            sigma_si_diff(loop) = max(abs(sigma_si_diff_dt));
            sigma_b_c_diff(loop) = max(abs(sigma_b_c_diff_dt));
            sigma_b_t_diff(loop) = max(abs(sigma_b_t_diff_dt));
            
            %% Finding acceleration ddn
            ddn = invM*(f_e - f_I);
            center_node_g(loop) = sqrt(ddn(end-3)^2 + ddn(end-1)^2 + ddn(end)^2);
            % Check the acceleration if it is excess the maximum value
            if center_node_g(end) >= acceleration_tol
                acceleration_excess = 1;
            end
            
            %% Finding n_-1 in the first loop
            if loop == 1
                n_minus1 = n - dn.*dt + 1/2.*ddn.*dt^2;
                n_iplus1 = 2.*n - n_minus1 + ddn*dt^2;
                dn_iplus1 = 1/(2*dt).*(2.*n - 2.*n_minus1 + 3.*ddn*dt^2);
            else
                n_iplus1 = 2*n - n_previous + ddn*dt^2;
                dn_iplus1 = 1/(2*dt)*(2*n - 2*n_previous + 3*ddn*dt^2);
            end
            
            %% Update position and velocity matrix
            n_previous = n;
            n = n_iplus1;
            dn = dn_iplus1;
            for i = 1:nnodes
                N(:,i) = n( 1+D*(i-1) : D*i ,1) ;
            end
            %% update strings matrix s and bars matrix b
            s = kron(C_s,I_D)*n;
            b = kron(C_b,I_D)*n;
            ds = kron(C_s,I_D)*dn;
            db = kron(C_b,I_D)*dn;
            %% Total energy (start from this first step not from the initial state)
            %{
            %  ------------------------- Potential energy ----------------------
            %  Gravity
            U_g = 0;
            for j = 1:nnodes
                U_g = U_g + (-g_earth) * m(j) * N(2,j);
            end
            % Elastic potential energy
            U_e = 0;
            for j = 1:n_s
                U_e = U_e + 1/2 * varepsilon_s(j) * sigma_s(j) * A_s * s_initiallength(j);
            end
            for j = 1:n_b
                U_e = U_e + 1/2 * varepsilon_b(j) * sigma_b(j) * A_b * b_initiallength(j);
            end
            U_po = [U_po, U_g + U_e];
            % ------------------------ Kinetic energy ---------------------------
            E_kj = 0;
            for j = 1:nnodes
               E_kj = E_kj + 1/2 * m(j) * ( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 + dn(3+D*(j-1))^2);
            end
            E_k = [E_k , E_kj ];
            % ------------------------ Total Energy ---------------------------
            U = [U , E_kj + U_g + U_e];
            %}
            %% currenttime
            currenttime = [currenttime , dt*loop ];
            %% deepest y position
            if min(N(2,:)) < deepest_y_of_each_orientation
                deepest_y_of_each_orientation = min(N(2,:));
            end
            
            %% Velocity Scaler
            
            V_each_step = 0;
            for j = 1:nnodes
                V_each_step = V_each_step + sqrt(( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 + dn(3+D*(j-1))^2));
            end
            V_each_step = V_each_step/nnodes;
            V_scaler = [V_scaler, V_each_step];
            % If V_scaler < V_tol for 'lopp_after_V_tol' loops, we assume it is stopped
            % loop_after_V_tol*dt = Time_stop
            if V_each_step <= V_tol
                loop_after_V_tol = loop_after_V_tol + 1;
            else
                loop_after_V_tol = 0;
            end
            if loop_after_V_tol >= Time_stop
                break
            end
            
            %%
            if acceleration_excess == 1
                %tenseg_plot(N,C_b,C_s);
                %axis off
                break
            end
            
        end
        % =========================================================================================================
        % The ending of one orientation
        % =========================================================================================================
        %% Record the max_g, depth, deviation distance and falling time of each orientation
        if acceleration_excess == 1
            Max_g_of_different_orientation(thetay_i,thetax_i) = max(center_node_g)/9.8;
            deviation_distance = [deviation_distance; -1];
            Falling_time = [Falling_time; -1];
            deepest_y  = [deepest_y ; deepest_y_of_each_orientation];
        else
            Max_g_of_different_orientation(thetay_i,thetax_i) = max(center_node_g)/9.8;
            deviation_distance = [deviation_distance; sqrt( N(1,nnodes)^2 + N(3,nnodes)^2) ];
            Falling_time = [Falling_time; dt*loop];
            deepest_y  = [deepest_y ; deepest_y_of_each_orientation];
        end
        %% Max and Min Stresses
        %Surface String Stresses
        sigma_ss_max_n(thetay_i,thetax_i) = max(sigma_ss_max);
        sigma_ss_min_n(thetay_i,thetax_i) = min(sigma_ss_min);
        %Internal String Stresses
        sigma_si_max_n(thetay_i,thetax_i) = max(sigma_si_max);
        sigma_si_min_n(thetay_i,thetax_i) = min(sigma_si_min);
        %Bar stresses
        sigma_bar_max_n(thetay_i,thetax_i) = max(sigma_bar_max);
        sigma_bar_min_n(thetay_i,thetax_i) = min(sigma_bar_min);
        %Difference between material limit and actual stress
        sigma_si_diff_n(thetay_i,thetax_i) = max(sigma_si_diff);
        sigma_ss_diff_n(thetay_i,thetax_i) = max(sigma_ss_diff);
        sigma_b_c_diff_n(thetay_i,thetax_i) = max(sigma_b_c_diff);
        sigma_b_t_diff_n(thetay_i,thetax_i) = max(sigma_b_t_diff);
        
        %Check if acceleration went over the given value
        if acceleration_excess == 1
            break
        end
        %% Output frames to make gif -- imwrite
        
        nImages = length(cell);
        
        for idx = 1:nImages
            fig = figure(idx+1);
            fill3([-30 30 30 -30],[0 0 0 0],[-30 -30 30 30],[0 0 0],'facealpha',0.3)
            hold on
            tenseg_plot_gif(cell{idx},C_b,C_s);
            set(fig, 'Position',  [600, 10, 1150, 1000])
            set(fig,'color',[1 1 1])
            xLabel = strcat('x',32,'[m]');
            yLabel = strcat('y',32,'[m]');
            xlabel(xLabel,'fontsize',28,'interpreter','latex')
            ylabel(yLabel,'fontsize',28,'interpreter','latex')
            set(gca,'FontSize',28,'TickLabelInterpreter','latex')
            axis([min(Xmin) max(Xmax) min(Ymin) max(Ymax) min(Zmin) max(Zmax)])
            drawnow
            set(fig,'visible','off')
            
            frame = getframe(fig);
            im{idx} = frame2im(frame);
            ax.Units = 'normalized';
        end
        close;
        
        filename = 'gif3D.gif'; % Specify the output file name
        for idx = 1:nImages
            [A,map] = rgb2ind(im{idx},256);
            if idx == 1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.06);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.06);
            end
        end
        
        
    end
end

%  ============================================================================================
%   End dynamic simulation
%   ============================================================================================
%axis([min(Xmin) max(Xmax) min(Ymin) max(Ymax) min(Zmin) max(Zmax)])
%% Total Mass
mass = sum(m,'all');

%mass_nodes = nnodes * node_mass;
%mass = mass_internal+mass_nodes;
Max_g_of_different_orientation = max(max(Max_g_of_different_orientation));


figure()
plot(currenttime,center_node_g/9.8)
xlabel('time [s]')
ylabel('Acceleration of the center node [g_earth]')

%re = reshape(Max_g_of_different_orientation,[20,42]);
sigma_ss_diff_n = max(max(abs(sigma_ss_diff_n)));
sigma_si_diff_n = max(max(abs(sigma_si_diff_n)));
sigma_b_c_diff_n = max(max(abs(sigma_b_c_diff_n)));
sigma_b_t_diff_n = max(max(abs(sigma_b_t_diff_n)));


toc
%
% sound(sin(2*pi*25*(1:4000)/100));
% pause(1);
% sound(sin(2*pi*25*(1:4000)/100));

end