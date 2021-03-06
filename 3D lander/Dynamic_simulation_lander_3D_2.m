clear all; close all; clc
% This code is to calculate distance, max acceleration, penetration and time

tic
%% Input parameters
% ====================== Lander Initial Properties ========================
% D = Dimension; I_D = Identity Matrix; height = Initial height
% r = radius of the lander; m_load = mass of central node
% m_s_total = Total mass of strings ; m_b_total = Total mass of bars
% v_0 = Initial vertical speed
D = 3;
I_D = eye(D);
height = 10;
r = 1.5;
v_0 = -5;
m_load = 500;
dtheta_max = 1*pi/180;
% ======================= Environment Properties ==========================
% pc = ground bounce constant; Cc = ground damping coefficient
% eta = ground frictional coefficient; % g_mars = Gravitational
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
dt = 1*10^(-4);
total_time = 20;
number_of_loop = total_time/dt;
V_tol = 0.05;
Time_stop = 2/dt;
acceleration_tol = 5000*9.8;
number_of_orientation = 1; % The number of orientations for each configurations
% ======================== Materials properties ===========================
% --------------------- bars -----------------------
rho_b = 5000;   
E_b = 60e09;
c_b = 0;
% ------------------- strings ----------------------
rho_s = 2500;
E_s = 1e08;
c_s = 5e06;
%% ======================== Display parameters ============================
number_of_configurations = 3;
X=round(linspace(1,number_of_loop,number_of_configurations));
%----- Initializing arbitrary node position matrices to be displayed ------   
display_node = [];    
%display_node_x_position = zeros(length(display_node),number_of_loop);
%display_node_y_position = zeros(length(display_node),number_of_loop);
%% Initializing Matrices or arrays
% ---------------------- Final distance arrays ----------------------------
Max_g_of_different_orientation = []; % Store Maximum acceleration of all orientations
deviation_distance = [];             % Store deviation distance of all orientations
Falling_time = [];                   % Store stop time of all orientations
deepest_y = [];                      % Store the deepest penetration of all orientations
theta = [];                          % Record theta of one orientation. We first simulate landers using a large timestep.
% if the maximum acceleration is greater than the allowed value, we need to re-simulate this rotation angle using a smaller timestep.
% We can save time in such way.
%% Initializing maximum and minimum value of x,y,z position
% They are the maximum of minimum valve for every loops
Xmax = []; Xmin = [];
Ymax = []; Ymin = [];
Zmax = []; Zmin = [];
%%  ============================================================================================
%   Start dynamic simulation 
%   ============================================================================================
% Loops for m_s (1/32*m_load,1/16*m_load,1/8*m_load,1/4*m_load,...)
for m_s_total = [1/16*m_load] 
    for p = 5:5 % p = 3, 4 ,5
        for q = 5:5 % q = 4 ,5
            %% ==================== node matrix and C_s, C_b ======================
            [N_norotation,C_b,C_s,nnodes,n_s,n_b] = Lander_3D(q,p);
            C_sT = C_s'; C_bT = C_b';
            %%  Loops for m_b (1/32*m_load,1/16*m_load,1/8*m_load,1/4*m_load,...)
                for m_b_total = [1/2*m_load]
                    % aa is the number of orientations for one configurations
                    for aa = 1:number_of_orientation
                        deepest_y_of_each_orientation = 0; 
                        %% Node position matrix n
                        [n,N,theta_0x,theta_0y,theta_0z] = nodematrix(N_norotation,height,nnodes);
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
                        if aa == 1
                            [M,m,A_s,A_b] = Mass(D,I_D,n_s,n_b,s,b,C_sT,C_bT,rho_s,rho_b,m_s_total,m_b_total,m_load);
                            invM = inv(M);
                        end
                        %%  Finding f_g (Only calculate once)
                        if aa == 1
                            g = repmat([0;g_mars;0],nnodes,1);
                            f_g = M*g;
                        end
                        %% Because 3D system runs very slow, so we prefer larger value of dt to decrease the running time.
                        % However, for some orientations, larger dt leads to high acceleration up to hundreds of g_earth.
                        % If we use lower dt to simulate those orientations, accelerations have a dramatical drop.
                        % Whereas, with regard to other orientations which have small accelerations under larger dt, 
                        % decreasing dt will not lead to further decrease of acceleration.
                        %
                        % Therefore, we first simulate an orientation using
                        % larger dt, e.g., 1*10^(-4), if it results in a very high acceleration, 
                        % we will re-simulate this orientation using lower dt.
                        % We only re-simulate this orientation once again.
                        resimulate = 0;
                        simulate = true;
                        while simulate                     
                            if resimulate == 1 % resimulate == 1 means we need to re-simulate this orientation again.
                                n = n_0 ;
                                N = N_0;
                                dn = dn_0;
                                s = s_0;
                                b = b_0; 
                                ds = ds_0;
                                db = db_0;
                                dt = 1*10^(-5);
                                number_of_loop = total_time/dt;
                            end
                            %% Initializing Energy matices and velocity scaler
                            U = [];
                            E_k = [];
                            U_po = [];
                            V_scaler = [];
                            center_node_g = [];
                            currenttime = [];
                            loop_after_V_tol = 0;
                            acceleration_excess = 0; 
                            %% ====== Start simulation for each orientation =========
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
                                figure(2);
                                % Finding the maximum and minimum value of x position
                                Xmax = [Xmax max(N(1,:))]; Xmin = [Xmin min(N(1,:))];
                                Ymax = [Ymax max(N(2,:))]; Ymin = [Ymin min(N(2,:))];
                                Zmax = [Zmax max(N(3,:))]; Zmin = [Zmin min(N(3,:))];
                                % Plot configurations
                                if ismember(loop,X) == 1
                                    tenseg_plot(N,C_b,C_s);
                                    hold on
                                    %ttext = string(find(X==loop));
                                    %text(N(1,1),N(2,1),N(3,1),ttext,'fontsize',15);
                                    %hold on
                                end
                                %}
                                %% External force
                                f_e = externalforce(D,nnodes,n,pc,dn,eta,f_g,Cc);
                                %% Internal force                              
                                [f_I,varepsilon_s,sigma_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength] = internalforce(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b,c_s,c_b,A_s,A_b);
                                %% Finding acceleration ddn
                                ddn = invM*(f_e - f_I);
                                center_node_g = [center_node_g, sqrt(ddn(end-3)^2 + ddn(end-1)^2 + ddn(end)^2) ];
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
                                Max_g_of_different_orientation = [Max_g_of_different_orientation; max(center_node_g)/9.8];
                                deviation_distance = [deviation_distance; -1]
                                Falling_time = [Falling_time; -1];    
                                deepest_y  = [deepest_y ; deepest_y_of_each_orientation]
                           else
                                Max_g_of_different_orientation = [Max_g_of_different_orientation; max(center_node_g)/9.8]
                                deviation_distance = [deviation_distance; sqrt( N(1,nnodes)^2 + N(3,nnodes)^2) ]
                                Falling_time = [Falling_time; dt*loop]   
                                deepest_y  = [deepest_y ; deepest_y_of_each_orientation]
                           end 
                            %% If the maximum acceleration is greater than 50*g_earth, we need to re-simulate this orientation again.
                            if max(center_node_g)/9.8 >= 50 && resimulate ~= 1
                                theta = [theta ; theta_0x  theta_0y  theta_0z ]
                                resimulate = 1;
                                Max_g_of_different_orientation(end) = [];
                                %deviation_distance(end) = [];
                                %Falling_time(end) = [];   
                                deepest_y(end)  = [];
                            else
                                simulate = false;
                                dt = 1*10^(-4);
                                number_of_loop = total_time/dt;
                            end

                        end
                        if acceleration_excess == 1
                            break
                        end
                    end 
                    % =========================================================================================================
                    % The ending of one configuration
                    % =========================================================================================================
                end 
            end
        end
end
%  ============================================================================================
%   End dynamic simulation 
%   ============================================================================================
%axis([min(Xmin) max(Xmax) min(Ymin) max(Ymax) min(Zmin) max(Zmax)])

figure()
plot(currenttime,center_node_g/9.8)
xlabel('time [s]')
ylabel('Acceleration of the center node [g_earth]')
    
re = reshape(Max_g_of_different_orientation,[20,42]);
    
toc

sound(sin(2*pi*25*(1:4000)/100));
pause(1);
sound(sin(2*pi*25*(1:4000)/100));



