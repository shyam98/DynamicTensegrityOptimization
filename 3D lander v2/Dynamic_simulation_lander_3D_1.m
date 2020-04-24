clear all; close all; clc
% This code is to export gif, plot energy figures, plot projectory of
% arbitrary nodes
% This is written by Liming Zhao, June 13, 2019
%% Input parameters
% ====================== Lander Initial Properties ========================
% D = Dimension; I_D = Identity Matrix; height = Initial height
% r = radius of the lander; m_load = mass of central node
% m_s_total = Total mass of strings ; m_b_total = Total mass of bars
% v_0 = Initial vertical speed
D = 3;
I_D = eye(D);
height = 10;
g_mars = -3.711;
r = 1.5;
dtheta_max = 1*pi/180;
accerelation_tol = 300*g_mars;
m_load = 400; m_s_total = 1/16*m_load; m_b_total = 1/2*m_load; 
v_0 = -5;

r_ss = 0.01;
r_si = 0.005;
r_b = 0.02;

% ======================= Environment Properties ==========================
% pc = ground bounce constant; Cc = ground damping coefficient
% eta = ground frictional coefficient; % g_mars = Gravitational
% acceleration on Mars
pc = 5e05; % f_n coefficient
Cc = 1e04;
eta = 0.5; % f_f coefficient
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
%% Materials properties
% --------------------- bars -----------------------
rho_b = 4705; %Density of Titanium Alloy
E_b = 11e10; % Elastic Modulus
c_b = 0; 
% ------------------- strings ----------------------
rho_s = 1123; % Density of Nylon wire
E_s = 1e08; % Elastic Modulus
c_s = 5e06;    % strain damping coefficient 
%% Choose one configuration
class = 1;
switch class
    case 1
        q = 4; p = 3;
    case 2
        q = 4; p = 4;
    case 3
        q = 4; p = 5;
    case 4
        q = 5; p = 3;
    case 5
        q = 5; p = 4;
    case 6
        q = 5; p = 5;
end
%% Initializing node matrix N and C_s, C_b
[N_norotation,C_b,C_s,nnodes,n_s,n_b] = Lander_3D(q,p);
C_sT = C_s'; C_bT = C_b';
%% Initializing maximum and minimum value of x,y,z position
% They are the maximum of minimum valve for every loops
% They are used to confine the range of axis if we plot diagrams.
Xmax = []; Xmin = [];
Ymax = []; Ymin = [];
Zmax = []; Zmin = [];
%% 
Max_g_of_different_orientation = []; % Store Maximum acceleration of all orientations
center_node_g = [];                  % Store the acceleration of the central node for every step
%% Initialize the number of output images
ii = 1;
%%  ============================================================================================
%  START DYNAMIC SIMULATION
%  ============================================================================================
number_of_orientation = 1; % Number of orientations per configuration
for aa = 1:number_of_orientation
    % Node position matrix n
	[n,N,theta_0x,theta_0y,theta_0z] = nodematrix(N_norotation,height,nnodes);
    n_0 = n ;
    %% Initial velocity dn
    dn = initialvelocity(N,height,v_0,dtheta_max,nnodes);
     %% Define strings matrix s and bars matrix b
     %  s and b indicates the length (2D - x and y direction) of each string and each bar
     s = kron(C_s,I_D)*n;
     b = kron(C_b,I_D)*n;
     ds = kron(C_s,I_D)*dn;
     db = kron(C_b,I_D)*dn;
     s_0 = s;
     b_0 = b; 
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
    %% Initializing Energy matices and velocity scaler
    U = [];
    E_k = [];
    U_po = [];
    V_scaler = [];
    %% Time controller and calculating the number of configurations to be displayed
    currenttime = [];
    number_of_configurations = 10;
    X=round(linspace(1,number_of_loop,number_of_configurations));
    loop_after_V_tol = 0;
    % velocity_tolerance = 0.001;
    %% Initializing arbitrary node position matrices to be displayed
    display_node = [1,6];
    display_node_x_position = zeros(length(display_node),number_of_loop);
    display_node_y_position = zeros(length(display_node),number_of_loop);
    display_node_z_position = zeros(length(display_node),number_of_loop);
    %% ============================================================================================

    %  ============================================================================================
    %  Start dynamic simulation for one lander
    %  ============================================================================================

    for loop = 1:number_of_loop
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
        [f_I,varepsilon_s,sigma_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength] = internalforce(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b,c_s,c_b,A_s,A_b);
        %% Finding acceleration ddn
        ddn = invM*(f_e - f_I);
        center_node_g = [center_node_g, sqrt(ddn(end-3)^2 + ddn(end-1)^2 + ddn(end)^2) ];
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
           E_kj = E_kj + 1/2 * m(j) * ( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 + dn(3+D*(j-1))^2 );
        end
        E_k = [E_k , E_kj ];
        % ------------------------ Total Energy ---------------------------
        U = [U , E_kj + U_g + U_e];
        %% currenttime
        currenttime = [currenttime , dt*loop ];
        %% Velocity Scaler
        V_each_step = 0;
        for j = 1:nnodes
            V_each_step = V_each_step + sqrt(( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 + dn(3+D*(j-1))^2 ));
        end
        V_each_step = V_each_step/nnodes;
        V_scaler = [V_scaler, V_each_step];
        % If the average velocity is smaller than V_tol for the number of 'lopp_after_V_tol' loops, we assume it stops
        % loop_after_V_tol*dt = Time_stop
        if V_each_step <= V_tol
            loop_after_V_tol = loop_after_V_tol + 1;
        else
            loop_after_V_tol = 0;
        end

        if loop_after_V_tol >= Time_stop
            break
        end
    end
    Max_g_of_different_orientation = [Max_g_of_different_orientation; max(center_node_g)/9.8];
    % Adjust the display parameters for configurations figure
    % view(90,0)
    axis([min(Xmin) max(Xmax) min(Ymin) max(Ymax) min(Zmin) max(Zmax) ])
    % ============================================================================================
    %  End dynamic simulation for one lander
    %  ============================================================================================
end

    %% Plotting Energy
    %{
    figure()
    plot(currenttime,U)
    xlabel('time [s]')
    ylabel('Total energy [J]')

    figure()
    plot(currenttime,U_po)
    xlabel('time [s]')
    ylabel('Potential energy [J]')

    figure()
    plot(currenttime,E_k)
    xlabel('time [s]')
    ylabel('Kinetic energy [J]')
    
    figure()
    plot(currenttime,V_scaler)
    xlabel('time [s]')
    ylabel('Average velocity [m/s]')
    %}
    %% Plot trajectory of arbitrary nodes
   %{
    figure();
    for i = 1:length(display_node)
        plot3(display_node_x_position(i,:),display_node_y_position(i,:),display_node_z_position(i,:))
        hold on
    end
    
    %axis equal
    xlabel('x')
    ylabel('y')
    
end
%}
%% Output frames to make gif -- imwrite
%{
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
%}
