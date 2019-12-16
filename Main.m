clc;
clear;

%x = ga(fun,nvars,A,b,Aeq,beq,lb,ub)
%x = coordinates = times in MJD
%fun to find delta V given the t (dep, ga, arr)


% Departure domain - from Mercury
min_dep_date = [2020 1 1 0 0 0];
max_dep_date = [2040 8 1 0 0 0];

% GA domain - Venus
min_ga_date = [2020 6 1 0 0 0];
max_ga_date = [2041 2 1 0 0 0];

% Arrival domain - to Mars
min_arr_date = [2021 1 1 0 0 0];
max_arr_date = [2041 8 1 0 0 0];


% Inequality constraints
A = [1 -1 0;
    0 1 -1;
    0 0 0];

duree_min_transfer1 = date2mjd2000([2000 6 1 0 0 0]) - date2mjd2000([2000 1 1 0 0 0]); % 5 mois
duree_min_transfer2 = date2mjd2000([2000 6 1 0 0 0]) - date2mjd2000([2000 1 1 0 0 0]); % 5 mois

b = [duree_min_transfer1;
    duree_min_transfer2;
    0];


% Boundaries for dates (in Julian days because inequalities)
lb = [date2mjd2000(min_dep_date), date2mjd2000(min_ga_date), date2mjd2000(min_arr_date)];
ub = [date2mjd2000(max_dep_date), date2mjd2000(max_ga_date), date2mjd2000(max_arr_date)];


% Set options
%options = optimoptions('ga','FitnessScalingFcn', {@fitscalingtop,1/10});
%options=gaoptimset('PlotFcn'{@gaplotbestf,@gaplotscores},'TolFun',0,'PopulationSize',4000,'Generations',100000,'Display','iter');
%options = optimoptions('ga','FitnessLimit', {-1e-100});

options = optimoptions('ga', 'FunctionTolerance', 1e-6, 'Display', 'off');
%options = optimoptions('ga', 'PlotFcn', @gaplotbestf);


% Genetic algorithms
[t, fval, exitflag, output, population, scores] = ga (@f, 3,A,b,[],[],lb,ub);


td = mjd20002date(t(1));
tga = mjd20002date(t(2));
ta = mjd20002date(t(3));

fprintf(['\ntd = \t', num2str(td), '\ntga = \t', num2str(tga), '\nta = ', num2str(ta), '\n']);

