clc;
clear;

%x = ga(fun,nvars,A,b,Aeq,beq,lb,ub)

%x = coordinates = times
%fun to find delta V given the t (dep, ga, arr)
% lb and ub to put in place the domain

% domains: nb of Julian days! because inequality 
% lower et upper bounds of the domains

[td, tga, ta] = ga (@f, 3,[],[],[],[],lb,ub);