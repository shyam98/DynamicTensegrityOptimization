IntCon = [1,2];
rng default % For reproducibility
fun = @MultiObjFun;
A = [];
b = [];
Aeq = [];
beq = [];
% [p;q;L;r_ss;r_si;r_b;C_2;z_position;RL_ratio]
options = optimoptions('ga','Display', 'iter');


lb = [4;4;0.5;0.0005;0.0005;0.0005;-0.01;-0.4;0.25];
ub = [8;8;1;0.0025;0.0025;0.005;0.01;0.4;1.5];
nonlcon = [];
[x,fval,exitflag,output,population,scores] = ga(fun,9,A,b,Aeq,beq,lb,ub,nonlcon,IntCon, options);


