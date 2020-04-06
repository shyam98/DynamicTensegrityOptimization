clear all; close all; clc

[N,Cb,Cs,nnodes,n_s,n_b] = Lander_3D(10, 4); %Input (p,q)

R3Ddata.Bradius = 0.01*ones(size(Cb,1),1);
R3Ddata.Sradius = 0.0025*ones(size(Cs,1),1);

%tenseg_plot( N,C_b,C_s,fig_handle,highlight_nodes,view_vec, PlotTitle, R3Ddata)
tenseg_plot( N,Cb,Cs)
tenseg_plot( N,Cb,Cs,[],[],[], [], R3Ddata)
axis off