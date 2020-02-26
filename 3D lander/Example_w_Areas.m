clear all; close all; clc

[N,Cb,Cs,nnodes,n_s,n_b] = Lander_3D(4,3);


R3Ddata.Bradius = 0.02*ones(size(Cb,1),1);
R3Ddata.Sradius = 0.0075*ones(size(Cs,1),1);

tenseg_plot( N,Cb,Cs)
tenseg_plot( N,Cb,Cs,[],[],[], [], R3Ddata)
axis off