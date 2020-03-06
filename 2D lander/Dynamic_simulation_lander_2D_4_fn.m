function [mass, Max_g_of_different_orientation, sigma_ss_max_n, sigma_ss_min_n, sigma_si_max_n, sigma_si_min_n, sigma_bar_max_n, sigma_bar_min_n, sigma_ss_diff_n, sigma_si_diff_n, sigma_b_c_diff_n, sigma_b_t_diff_n] = Dynamic_simulation_lander_2D_4_fn(r, r_ss, r_si, r_b, p_range, q_range, number_of_orientation)
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
m_load = 400;
v_0 = -5;
A_s = [A_ss ; A_si];
dtheta_0 = 1*pi/180;

stress_ss_max = [];
stress_ss_min = [];
stress_si_max = [];
stress_si_min = [];
stress_b_max = [];
stress_b_min = [];


sigma_ss_max_n = [];
sigma_ss_min_n = [];

sigma_si_max_n = [];
sigma_si_min_n = [];

sigma_bar_max_n = [];
sigma_bar_min_n = [];

sigma_si_diff_n = [];
sigma_ss_diff_n = [];
sigma_b_c_diff_n = [];
sigma_b_t_diff_n = [];


Yield_Nylon = 9.4e7;
Yield_Titanium = 1e9;
Youngs_Titanium = 115e9;
% ======================= Environment Properties ==========================
% pc = ground bounce constant; Cc = ground damping coefficient
% eta = ground frictional coefficient; % g_mars = Gravitational
% acceleration on Mars
pc = 5e05;  
Cc = 1e04;
eta = 0.5; 
g_mars = -3.711;
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
%number_of_orientation = 1; % The number of orientations for each configurations
% ======================== Materials properties ===========================
% --------------------- bars -----------------------
rho_b = 4705; % Density of Titanium Aloy
E_b = 11e10; % Elastic Modulus
c_b = 0;      % Strain damping coefficient
% ------------------- strings ----------------------
rho_s = 1123; % Density
E_s = 1e08;   % Elastic Modulus
c_s = 5e06;   % Strain damping coefficient 
%% ======================== Display parameters ============================
% number_of_configurations means how many configurations will be displayed.
number_of_configurations = 10;
X=round(linspace(1,number_of_loop,number_of_configurations));
%----- Initializing arbitrary node position matrices to be displayed ------   
% display_node matrix contains the nodes to be displayed, 
% e.g., if display_node = [1,6], the trajectory of the first and the fixth nodes will be displayed
display_node = [];    
display_node_x_position = zeros(length(display_node),number_of_loop);
display_node_y_position = zeros(length(display_node),number_of_loop);
%% Initializing Matrices or arrays
% ---------------------- Final distance arrays ----------------------------
Max_g_of_different_orientation = []; % Store Maximum acceleration of all orientations
deviation_distance = [];             % Store deviation distance of all orientations
Falling_time = [];                   % Store stop time of all orientations
deepest_y = [];                      % Store the deepest penetration of all orientations
%%  ============================================================================================
%   Start dynamic simulation 
%   ============================================================================================

for q = q_range
    %% ================ Initialize C_s and C_b ====================
    nnodes = q + 1;
    n_s = 2*q;  
    n_b = length(p_range)*q; 
    basic = linspace(0,2*pi,q+1);
    N_norotation = zeros(D, nnodes);
    for i_N = 1:nnodes-1
        N_norotation(1,i_N) = r*cos(basic(i_N));
        N_norotation(2,i_N) = r*sin(basic(i_N));
    end
    C_s = zeros(n_s,nnodes);
    C_b = zeros(n_b,nnodes);
    % -------------------------- C_s ------------------------------
    for i_cs = 1:q
        C_s(i_cs,i_cs) = -1;
        if i_cs ~= q
            C_s(i_cs,i_cs+1) = 1;
        else
            C_s(i_cs,1) = 1;
        end
    end
    for i_cs = (q+1):2*q
        C_s(i_cs,nnodes) = -1;
        C_s(i_cs,i_cs-q) = 1;
    end
    % -------------------------- C_b ------------------------------ 
    count = 0;
    for p = p_range
        for i_cb = 1:q
           C_b(i_cb+count,i_cb) = -1;
           if i_cb+p <= q
              C_b(i_cb+count,i_cb+p) = 1;
           else
              C_b(i_cb+count,i_cb+p-q) = 1;
           end
        end
        count = count + q;
    end
    C_sT = C_s'; C_bT = C_b';
    % loops: The number of orientations for each configuration
    
    % Keep Track of figures for each orientation
    max_g_n = [];
    for aa = 1:number_of_orientation
        %% Node position matrix n
        [n,N] = nodematrix(N_norotation,height,nnodes);
        n_0 = n ;
        %% Initial velocity dn
        dn = initialvelocity(N,height,v_0,dtheta_0,nnodes);
        %% Define strings matrix s and bars matrix b
        %  s and b indicates the length (2D : x and y direction) of each string and each bar
        s = kron(C_s,I_D)*n;
        b = kron(C_b,I_D)*n;
        ds = kron(C_s,I_D)*dn;
        db = kron(C_b,I_D)*dn;
        s_0 = s;
        b_0 = b; 
        %% Mass (Only calculate once for each configuration)
        % aa = 1 means This is the first orientation for one configuration, so
        % we calculate its Mass matrix. We do not have to calculate mass for
        % the remaining orientations.
        if aa == 1
            [M,m] = Mass_2(D,I_D,n_s,n_b,s,b,C_sT,C_bT,rho_s,rho_b,m_load, A_s, A_b);
            invM = inv(M);   
        end
        %%  Finding f_g (Only calculate once)
        if aa == 1
            g = repmat([0;g_mars],nnodes,1);
            f_g = M*g;
        end
        %% Initializing Energy matices and velocity scaler
        U = [];   % Total Energy
        E_k = [];  % Kinetic Energy
        U_po = [];  % Potential Energy 
        V_scaler = []; % Vector used to store the average speed of all nodes per simulation loops, and used to plot speed vs. time plot.
        center_node_g = []; % Store the acceleration of central nodes
        currenttime = []; % Vector used to store current simulation time.
        loop_after_V_tol = 0; % Used to record the number of loops when the speed is below 'V_tol'. If 'loop_after_V_tol' >= 'Time_stop', the lander stops.
        acceleration_excess = 0; %  acceleration_excess = 1 represents that the acceleration excesses the maximum allowed value, acceleration_excess = 0 represents that the acceleration is below the maximum allowed value
        deepest_y_of_each_orientation = 0;
        sigma_ss_max = [];
        sigma_ss_min = [];
        
        sigma_si_max = [];
        sigma_si_min = [];
        
        sigma_bar_max = [];
        sigma_bar_min = [];
        
        sigma_ss_diff = [];
        sigma_si_diff = [];
        sigma_b_c_diff = [];
        sigma_b_t_diff = [];
        
        %% ================================================
        %       Start simulation for each configuration
        % =================================================
        for loop = 1:number_of_loop
            %% Finding trajectory matrices of arbitrary nodes
            %{
            for j = 1:length(display_node)
                display_node_x_position(j,loop) = N(1,display_node(j));
                display_node_y_position(j,loop) = N(2,display_node(j));
            end
            %}
            %% Plotting an arbitrary number of lander configurations throughout the time frame
            %{
            figure(1);
            % Finding the maximum and minimum value of x position
            Xmax = [Xmax max(N(1,:))]; Xmin = [Xmin min(N(1,:))];
            Ymax = [Ymax max(N(2,:))]; Ymin = [Ymin min(N(2,:))];
            % Add z coordination
            N_z = zeros(1,nnodes);
            N_plot = [N;N_z];
            % Plot configurations
            if ismember(loop,X) == 1
                tenseg_plot(N_plot,C_b,C_s);
                hold on
                ttext = string(find(X==loop));
                %text(N(1,1),N(2,1),ttext,'fontsize',15);
                hold on
            end
            %}
            %% External force
            f_e = externalforce(D,nnodes,n,pc,dn,eta,f_g,Cc);
            %% Internal force                              
            [f_I,varepsilon_s,sigma_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength, sigma_ss_dt, sigma_si_dt, sigma_ss_diff_dt, sigma_si_diff_dt, sigma_b_c_diff_dt, sigma_b_t_diff_dt] = internalforce_2(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b,c_s,c_b,A_s,A_b,Yield_Nylon, Youngs_Titanium, Yield_Titanium);
            %% Keep track of max and min stresses
            sigma_ss_max = [sigma_ss_max; max(sigma_ss_dt)];
            sigma_ss_min = [sigma_ss_min; min(sigma_ss_dt)];
            
            sigma_si_max = [sigma_si_max; max(sigma_si_dt)];
            sigma_si_min = [sigma_si_min; min(sigma_si_dt)];

            
            sigma_bar_max = [sigma_bar_max; max(sigma_b)];
            sigma_bar_min = [sigma_bar_min; min(sigma_b)];
            
            
            sigma_ss_diff = [sigma_ss_diff; sigma_ss_diff_dt];
            sigma_si_diff = [sigma_si_diff; sigma_si_diff_dt];
            sigma_b_c_diff = [sigma_b_c_diff; sigma_b_c_diff_dt];
            sigma_b_t_diff = [sigma_b_t_diff; sigma_b_t_diff_dt];
            
            %% Finding acceleration ddn
            ddn = invM*(f_e - f_I);
            %% Check the acceleration of the central node
            center_node_g = [center_node_g, sqrt(ddn(1+D*(nnodes-1))^2 + ddn(2+D*(nnodes-1))^2) ];
            %  acceleration_excess = 1 represents that the acceleration excesses the maximum allowed
            %  value.
            if sqrt(ddn(1+D*(nnodes-1))^2 + ddn(2+D*(nnodes-1))^2) >= abs(acceleration_tol)
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
            %% Total energy 
            %  ------------------------- Potential energy ----------------------
            %  Gravity
            U_g = 0;
            for j = 1:nnodes
                U_g = U_g + (-g_mars) * m(j) * N(2,j);
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
               E_kj = E_kj + 1/2 * m(j) * ( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 );
            end
            E_k = [E_k , E_kj ];
            % ------------------------ Total Energy ---------------------------
            U = [U , E_kj + U_g + U_e];
            %% currenttime
            currenttime = [currenttime , dt*loop ];
            %% Finding deepest y position ( or finding the peneration of the lander)
            if min(N(2,:)) < deepest_y_of_each_orientation
                deepest_y_of_each_orientation = min(N(2,:));
            end

            %% Using velocity threshold to determine whether the lander stops or not
            V_each_step = 0;
            for j = 1:nnodes
                V_each_step = V_each_step + sqrt(( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 ));
            end
            V_each_step = V_each_step/nnodes;
            V_scaler = [V_scaler, V_each_step];
            % If V_scaler < V_tol for the number of 'lopp_after_V_tol' loops, we assume it stops
            % loop_after_V_tol*dt = Time_stop
            if V_each_step <= V_tol
                loop_after_V_tol = loop_after_V_tol + 1;
            else
                loop_after_V_tol = 0;
            end
            if loop_after_V_tol >= Time_stop
                break
            end
            %% Display the figure right after its acceleration excesses the allowed value
            %{
            if acceleration_excess == 1                    
                N_z = zeros(1,nnodes);
                N_plot = [N;N_z];
                tenseg_plot(N_plot,C_b,C_s);
                axis off
                break
            end
            %}
        end
       % =========================================================================================================
       % The ending of simulating one orientation
       % =========================================================================================================
        %% Record the max_g, depth, deviation distance and falling time of each orientation
        % If the acceleration excesses the maximum allowed value, 
        % we assume its deviation distance and stopping time are -1.
       max_g_n = mean(Max_g_of_different_orientation);
       if acceleration_excess == 1
            Max_g_of_different_orientation = [Max_g_of_different_orientation; max(center_node_g)/9.8];
            deviation_distance = [deviation_distance; -1];
            Falling_time = [Falling_time; -1];    
            deepest_y  = [deepest_y ; deepest_y_of_each_orientation];
       else
            Max_g_of_different_orientation = [Max_g_of_different_orientation; max(center_node_g)/9.8];
            deviation_distance = [deviation_distance; abs(N(1,nnodes))];
            Falling_time = [Falling_time; dt*loop];   
            deepest_y  = [deepest_y ; deepest_y_of_each_orientation];
       end
            %% Max and Min Stresses
        sigma_ss_max_n = [sigma_ss_max_n; max(sigma_ss_max)];
        sigma_ss_min_n = [sigma_ss_min_n; min(sigma_ss_min)];
        
        sigma_si_max_n = [sigma_si_max_n; max(sigma_si_max)];
        sigma_si_min_n = [sigma_si_min_n; min(sigma_si_min)];
        
        sigma_bar_max_n = [sigma_bar_max_n; max(sigma_bar_max)];
        sigma_bar_min_n = [sigma_bar_min_n; min(sigma_bar_min)];
        
        sigma_si_diff_n = [sigma_si_diff_n; max(sigma_si_diff)];
        sigma_ss_diff_n = [sigma_ss_diff_n; max(sigma_ss_diff)];
        sigma_b_c_diff_n = [sigma_b_c_diff_n; max(sigma_b_c_diff)];
        sigma_b_t_diff_n = [sigma_b_t_diff_n; max(sigma_b_t_diff)];
    %% Plot outputs to find averages
        mean_g(aa) = mean(Max_g_of_different_orientation)
        mean_sig_ss_max(aa) = mean(sigma_ss_max_n)/1e8;
        mean_sig_si_max(aa) = mean(sigma_si_max_n)/1e8;
        mean_sig_ss_min(aa) = mean(sigma_ss_min_n)/1e8;
        mean_sig_si_min(aa) = mean(sigma_si_min_n)/1e8;
        mean_sig_b_max(aa) = mean(sigma_bar_max_n)/1e8;
        mean_sig_b_min(aa) = mean(sigma_bar_min_n)/1e8;
        
        if aa>10
            std_mean_g = std(mean_g(:,end-5:end));
            std_mean_sig_ss_max = std(mean_sig_ss_max(:,end-5:end));
            std_mean_sig_si_max = std(mean_sig_si_max(:,end-5:end));
            std_mean_sig_ss_min = std(mean_sig_ss_min(:,end-5:end));
            std_mean_sig_si_min = std(mean_sig_si_min(:,end-5:end));
            std_mean_sig_b_max = std(mean_sig_b_max(:,end-5:end));
            std_mean_sig_b_min = std(mean_sig_b_min(:,end-5:end));
        end
        
        
        if aa>10 && std_mean_g<0.2 && std_mean_sig_ss_max<0.2 && std_mean_sig_si_max<0.2 && std_mean_sig_ss_min<0.2 && std_mean_sig_si_min<0.2 && std_mean_sig_b_max<0.2 && std_mean_sig_b_min<0.2
            plot(1:aa, mean_g)
            break;
        end
        
        
    end

    % =========================================================================================================
    % The ending of one configuration
    % =========================================================================================================

end

%  ============================================================================================
%   End dynamic simulation 
%   ============================================================================================
%% Total Mass
mass = sum(m,'all');      
%% Max Min Stress Vectors


%{
figure()
    plot(currenttime,center_node_g/9.8)
    xlabel('time [s]')
    ylabel('Acceleration of the center node [m/s^2]')

%}

toc
