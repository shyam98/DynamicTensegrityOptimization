clear all; close all; clc

%[N,Cb,Cs,nnodes,n_s,n_b] = Lander_3D(4,3, 1.5);

q = 6;
p = 6;
r_ss = 0.001;
%r_ss = 0;
%r_b = 0;
r_si = 0.0005;
r_b = 0.003;


L = 1;
cyl = 'RCC';
z_position = 0;
C_2 = 0;
RL_Ratio = 0.5;
[N,Cb,Cs,nnodes,n_s,n_b, n_ss, V_c] = Lander_3D(q,p,L,cyl, C_2, z_position, RL_Ratio);

central_node_r = 0.02;
lander_node_r = 0.008;
R3Ddata.Bradius = r_b*ones(size(Cb,1),1);  
R3Ddata.Sradius = r_ss*ones(size(Cs,1),1);
R3Ddata.Nradius = [ones((nnodes-1),1)*lander_node_r ; central_node_r];
%tenseg_plot( N,Cb,Cs, [], [], [0,90,0])

%highlightNode = [1];

tenseg_plot4( N,Cb,Cs,[],[],[0,90,0], [], []);
%view(62.2344, 25.3369)
 az = 10;
 el = 20;
 view(az, el);
cameratoolbar('SetCoordSys','z')
axis off