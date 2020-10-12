clc; clear all; close all; 
IntCon = [1,2];
fun = @MultiObjFun;
A = [];
b = [];
Aeq = [];
beq = [];
% [p;q;L;r_ss;r_si;r_b;C_2;z_position;RL_ratio]
options = optimoptions('ga','PlotFcn',"gaplotbestf",'PopulationSize', 270, 'UseParallel', true, 'MaxGenerations', 90, 'MaxStallGenerations', 5);


lb = [4;4;0.5;0.0005;0.0005;0.001;-0.01;-0.4;0.25];
ub = [8;8;1;0.002;0.002;0.004;0.01;0.4;1];
nonlcon = [];
rng default % For reproducibility
[x,fval,exitflag,output,population,scores] = ga(fun,9,A,b,Aeq,beq,lb,ub,nonlcon,IntCon, options);

