clear all; close all; clc

%[N,Cb,Cs,nnodes,n_s,n_b] = Lander_3D(4,3, 1.5);

q = 5;
p = 7;
r = 1.5;
L = 4;
cyl = 'SP';
z_position = 0;
C_2 = 0;
[N,Cb,Cs,nnodes,n_s,n_b, zl_i] = Lander_3D(q,p,r,L,cyl, C_2, z_position);

R3Ddata.Bradius = 0.02*ones(size(Cb,1),1);
R3Ddata.Sradius = 0.00*ones(size(Cs,1),1);

tenseg_plot( N,Cb,Cs)
%tenseg_plot( N,Cb,Cs,[],[],[], [], R3Ddata)
axis on