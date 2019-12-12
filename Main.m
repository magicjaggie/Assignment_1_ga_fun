clc;
clear;

%x = ga(fun,nvars,A,b,Aeq,beq,lb,ub)

%x = coordinates = times in MJD
%fun to find delta V given the t (dep, ga, arr)
% lb and ub to put in place the domain

% domains: nb of Julian days! because inequality 
% lower et upper bounds of the domains

% Departure domain - from Mercury
min_dep_date = [2020 1 1 0 0 0];
max_dep_date = [2020 8 1 0 0 0];

% GA domain - Venus
min_ga_date = [2020 6 1 0 0 0];
max_ga_date = [2021 2 1 0 0 0];

% Arrival domain - to Mars
min_arr_date = [2021 1 1 0 0 0];
max_arr_date = [2021 8 1 0 0 0];

lb = [date2mjd2000(min_dep_date), date2mjd2000(min_ga_date), date2mjd2000(min_arr_date)];
ub = [date2mjd2000(max_dep_date), date2mjd2000(max_ga_date), date2mjd2000(max_arr_date)];


% Genetic algorithms
[td, tga, ta] = ga (@f, 3,[],[],[],[],lb,ub);

