clear all; close all; clc
% This code is to export gif, plot energy figures, plot projectory of
% arbitrary nodes for one configuration
% This is written by Liming Zhao, Jan 12, 2020
%% Input parameters
% ====================== Lander Initial Properties ========================
% D = Dimension; I_D = Identity Matrix; height = Initial height
% r = radius of the lander; m_load = mass of central node
% m_s_total = Total mass of strings ; m_b_total = Total mass of bars
% v_0 = Initial vertical speed
D = 2;
I_D = eye(D);
height = 10;
r = 1.5;
m_load = 400;
m_s_total = 1/8*m_load; m_b_total = 1/2*m_load; 
v_0 = -5;
% ======================= Environment Properties ==========================
% pc = ground bounce constant; Cc = ground damping coefficient
% eta = ground frictional coefficient; % g_mars = Gravitational
% acceleration on Mars
pc = 5e05; % f_n coefficient
Cc = 1e04;
eta = 0.5; % f_f coefficient
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
total_time = 40;
number_of_loop = total_time/dt;
V_tol = 0.05;
Time_stop = 5/dt;
accerelation_tol = 25*g_mars;
%% Materials properties
% --------------------- bars -----------------------
rho_b = 4705; % Density of Titanium Aloy
E_b = 11e10; % Elastic Modulus
c_b = 0;      % Strain damping coefficient
Y_b = 0;      %Max elastic yield strength
% ------------------- strings ----------------------
rho_s = 2500; % Density
E_s = 1e08;   % Elastic Modulus
c_s = 5e06;   % Strain damping coefficient 
Y_b = 0;      %Max elastic yield strength

%% Initializing node matrix and C_s, C_b
% class is from 1 to 22, which represents 22 different configurations in
% 2D situation
class = 2;
switch class
    case 1
        q = 5; p_range = 2;
    case 2 
        q = 6; p_range = 2;
    case 3
        q = 7; p_range = 2;
    case 4
        q = 8; p_range = 2;
    case 5
        q = 9; p_range = 2;
    case 6
        q = 10; p_range = 2;
    case 7 
        q = 7; p_range = 3;
    case 8
        q = 8; p_range = 3;
    case 9
        q = 9; p_range = 3;
    case 10
        q = 10; p_range = 3;
    case 11
        q = 9; p_range = 4;
    case 12
        q = 10; p_range = 4;
    case 13
        q = 7; p_range = [2,3];
    case 14
        q = 8; p_range = [2,3];
    case 15
        q = 9; p_range = [2,3];
    case 16
        q = 10; p_range = [2,3];
    case 17
        q = 9; p_range = [2,4];
    case 18
        q = 10; p_range = [2,4];
    case 19
        q = 9; p_range = [3,4];
    case 20
        q = 10; p_range = [3,4];
    case 21
        q = 9; p_range = [2,3,4];
    case 22
        q = 10; p_range = [2,3,4];
end
% -------------------------------------------------------------------------
nnodes = q + 1; % Number of nodes
n_s = 2*q; n_b = length(p_range)*q;  % n_s: Number of strings ; n_b: Number of bars

basic = linspace(0,2*pi,q+1);
N_norotation = zeros(D, nnodes);
for i_N = 1:nnodes-1
    N_norotation(1,i_N) = r*cos(basic(i_N));
    N_norotation(2,i_N) = r*sin(basic(i_N));
end
C_s = zeros(n_s,nnodes);
C_b = zeros(n_b,nnodes);
% ------------------------ C_s --------------------------------------------
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
% ------------------------ C_b --------------------------------------------
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

%% Initializing maximum and minimum value of x,y,z position
% They are the maximum of minimum valve for every loops.
% They are used to confine the range of axis if we plot diagrams.
Xmax = []; Xmin = [];
Ymax = []; Ymin = [];


%% Initialize the number of output images
%  Used for plotting gif image
ii = 1;
%%  ============================================================================================
%  START DYNAMIC SIMULATION
%  ============================================================================================
number_of_Orientation = 1; % Number of orientations per configuration
for aa = 1:number_of_Orientation
    % The function ‘nodematrix’ is used to produce node position matrix n and N
	[n,N] = nodematrix(N_norotation,height,nnodes);
    n_0 = n ;
    %% Initial velocity dn
    dtheta_0 = 30*pi/180;
    dn = initialvelocity(N,height,v_0,dtheta_0,nnodes);
     %% Define strings matrix s and bars matrix b
     %  s and b indicates the length (2D - x and y direction) of each string and each bar
     s = kron(C_s,I_D)*n;
     b = kron(C_b,I_D)*n;
     ds = kron(C_s,I_D)*dn;
     db = kron(C_b,I_D)*dn;
     s_0 = s;
     b_0 = b; 
    %% Mass (Only calculate once)
    % aa = 1 means This is the first orientation for one configuration, so
    % we calculate its Mass matrix. We do not have to calculate mass for
    % the remaining orientations.
    if aa == 1
        [M,m,A_s,A_b] = Mass(D,I_D,n_s,n_b,s_0,b_0,C_sT,C_bT,rho_s,rho_b,m_s_total,m_b_total,m_load,q,r);
        invM = inv(M);
    end
    
    %%  Finding f_g (Only calculate once)
    if aa == 1
        g = repmat([0;g_mars],nnodes,1);
        f_g = M*g;
    end
    %% Initializing Energy matices and time
    U = []; % Total Energy
    E_k = []; % Kinetic Energy
    U_po = []; % Potential Energy 
    V_scaler = []; % Vector used to store the average speed of all nodes per simulation loops, and used to plot speed vs. time plot.
    currenttime = []; % Vector used to store current simulation time.
    loop_after_V_tol = 0; % Used to record the number of loops when the speed is below 'V_tol'. If 'loop_after_V_tol' >= 'Time_stop', the lander stops.
    %% X represents whether plot the current lander or not during one simulation
    % number_of_configurations means how many configurations will be displayed.
    number_of_configurations = 10;
    X=round(linspace(1,number_of_loop,number_of_configurations));

    %% Initializing arbitrary node position matrices to be displayed
    % 'display_node matrix' contains the nodes to be displayed, 
    % e.g., if display_node = [1,6], the trajectory of the first and the fixth nodes will be
    % displayed
    display_node = [1,6];
    display_node_x_position = zeros(length(display_node),number_of_loop);
    display_node_y_position = zeros(length(display_node),number_of_loop);
    %% ============================================================================================
    %  START DYNAMIC SIMULATION FOR ONE ORIENTATION
    %  ============================================================================================
    for loop = 1:number_of_loop
        %% Finding trajectory matrices of arbitrary nodes
        for j = 1:length(display_node)
            display_node_x_position(j,loop) = N(1,display_node(j));
            display_node_y_position(j,loop) = N(2,display_node(j));
        end
         %% Save the node position matrix every 0.1s to make gif

        % Finding the maximum and minimum value of x position
        Xmax = [Xmax max(N(1,:))]; Xmin = [Xmin min(N(1,:))];
        Ymax = [Ymax max(N(2,:))]; Ymin = [Ymin min(N(2,:))];

        % Add z coordination
        N_z = zeros(1,nnodes);
        N_plot = [N;N_z];
        % Plot configurations
        if rem(loop*dt,0.1) == 0 || loop == 1
            cell{ii} = N_plot;
            ii = ii + 1;
        end
       

        %% Plotting an arbitrary number of lander configurations throughout the time frame

        figure(1);
        % Add z coordination
        N_z = zeros(1,nnodes);
        N_plot = [N;N_z];
        % Plot configurations is loop is in X
        if ismember(loop,X) == 1
            tenseg_plot(N_plot,C_b,C_s);
            axis off
            hold on
            N_plot;
            ttext = string(find(X==loop));
            %text(N(1,1),N(2,1),ttext,'fontsize',15);
            hold on       
        end

        
        %% External force
        f_e = externalforce(D,nnodes,n,pc,dn,eta,f_g,Cc);
        %% Internal force                              
        [f_I,varepsilon_s,sigma_s,varepsilon_b,sigma_b,s_initiallength,b_initiallength, f_s] = internalforce(D,I_D,C_sT,C_bT,s_0,b_0,s,b,n_s,n_b,ds,db,E_s,E_b,c_s,c_b,A_s,A_b);
        %% Finding acceleration ddn
        ddn = invM*(f_e - f_I);
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
           E_kj = E_kj + 1/2 * m(j) * ( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 );
        end
        E_k = [E_k , E_kj ];
        % ------------------------ Total Energy ---------------------------
        U = [U , E_kj + U_g + U_e];
        %% currenttime
        currenttime = [currenttime , dt*loop ];

        %% Using velocity threshold to determine whether the lander stops or not
        V_each_step = 0;
        for j = 1:nnodes
           V_each_step = V_each_step + sqrt(( dn(1+D*(j-1))^2 + dn(2+D*(j-1))^2 ));
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
end
    % ============================================================================================
    %  End dynamic simulation for one lander
    %  ============================================================================================   
    %% Plot trajectory of arbitrary nodes
   %{
    figure();
    for i = 1:length(display_node)
        plot(display_node_x_position(i,:),display_node_y_position(i,:))
        hold on
    end
    
    %axis equal
    xlabel('x')
    ylabel('y')
    
end
%}

%% Output frames to make gif
%{
nImages = length(cell);

for idx = 1:nImages
    fig = figure(idx);
    fill3([-30 30 30 -30],[0 0 0 0],[-30 -30 30 30],[0 0 0])
    hold on
    tenseg_plot(cell{idx},C_b,C_s);
    set(fig, 'Position',  [600, 10, 1150, 1000])
    set(fig,'color',[1 1 1]) 
    xLabel = strcat('x',32,'[m]');
    yLabel = strcat('y',32,'[m]');
    xlabel(xLabel,'fontsize',28,'interpreter','latex')
    ylabel(yLabel,'fontsize',28,'interpreter','latex')
    set(gca,'FontSize',28,'TickLabelInterpreter','latex')
    axis([min(Xmin) max(Xmax) min(Ymin) max(Ymax)])
    drawnow
    set(fig,'visible','off')
    
    frame = getframe(fig);
    im{idx} = frame2im(frame);
    ax.Units = 'normalized';
end
close;

filename = 'gif2D.gif'; % Specify the output file name
for idx = 1:nImages
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
%}

  %% Plotting Energy figures
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