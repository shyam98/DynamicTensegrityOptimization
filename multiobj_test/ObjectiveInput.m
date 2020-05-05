IntCon = [1,2];
rng default % For reproducibility
fun = @MultiObjFun;
A = [];
b = [];
Aeq = [];
beq = [];
% [p;q;L;r_ss;r_si;r_b;C_2;z_position;RL_ratio]

cyl = 'SP';
lb = [4;4;0.5;0.005;0.005;0.01;-0.01;-0.8;0.25,cyl];
ub = [10;10;1.5;0.02;0.02;0.03;0.01;0.8;1.5,cyl];
nonlcon = [];
x = ga(fun,10,A,b,Aeq,beq,lb,ub,nonlcon,IntCon);
