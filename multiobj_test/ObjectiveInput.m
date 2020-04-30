IntCon = [1,2];
rng default % For reproducibility
fun = @MultiObjFun;
A = [];
b = [];
Aeq = [];
beq = [];
% [p;q;r;r_ss;r_si;r_b;C_2;z_position;RL_ratio(if RCC)]
lb = [4;4;0.25;0.005;0.005;0.01;-0.01;-0.75];
ub = [5;5;1.5;0.02;0.02;0.03;0.01;0.75];
nonlcon = [];
x = ga(fun,8,A,b,Aeq,beq,lb,ub,nonlcon,IntCon)
