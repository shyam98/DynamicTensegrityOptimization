clear all; close all; clc

%[N,Cb,Cs,nnodes,n_s,n_b] = Lander_3D(4,3, 1.5);

q = 6;
p = 4;
%r = 1.5;
L = 4;
cyl = 'SP';
z_position = 0;
C_2 = 0.1;
RL_Ratio = 0.5;
[N,Cb,Cs,nnodes,n_s,n_b, n_ss, V_c] = Lander_3D(q,p,L,cyl, C_2, z_position, RL_Ratio);

R3Ddata.Bradius = 0.02*ones(size(Cb,1),1);  
R3Ddata.Sradius = 0.00*ones(size(Cs,1),1);

tenseg_plot( N,Cb,Cs, [], [], [0,90,0])
%tenseg_plot( N,Cb,Cs,[],[],[], [], R3Ddata)
axis on