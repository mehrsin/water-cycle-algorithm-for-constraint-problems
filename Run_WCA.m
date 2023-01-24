clear all;
clc;
close all;

% % % %  Problem 1
% nvars=7;
% LB=ones(1,nvars)*-10;
% UB=ones(1,nvars)*10;
% constraints=@Const_Problem_1;
% objective_function=@Problem_1;

% % % %  Problem 3
% nvars=10;
% LB=ones(1,nvars)*0;
% UB=ones(1,nvars)*1;
% constraints=@Const_Problem_3;
% objective_function=@Problem_3;

% % % %  Problem 7
% nvars=3;
% LB=ones(1,nvars)*0;
% UB=ones(1,nvars)*100;
% constraints=@Const_Problem_7;
% objective_function=@Problem_7;

% % % %  Problem 15
nvars=2;
LB=ones(1,nvars)*0;
UB=ones(1,nvars)*4;
constraints=@Const_Problem_15;
objective_function=@Problem_15;

% % % %  Problem 24
% nvars=6;
% LB=ones(1,nvars)*0;
% UB=ones(1,nvars)*1;
% constraints=@Const_Problem_24;
% objective_function=@Problem_24;


[Xmin,Fmin,Sum_Const,NFEs,Elapsed_Time]=Const_WCA(objective_function,constraints,LB,UB,nvars)

