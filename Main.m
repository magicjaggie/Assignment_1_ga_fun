clc;
clear;

%x = ga(fun,nvars,A,b,Aeq,beq,lb,ub)
%x = coordinates = times in MJD
%fun to find delta V given the t (dep, ga, arr)

%% INPUTS

% Planets
dep_Id = 1; % Mercury
ga_Id = 2; % Venus
arr_Id = 4; % Mars

% Departure domain - from Mercury
min_dep_date = [2020 1 1 0 0 0];
max_dep_date = [2020 8 1 0 0 0];

% GA domain - Venus
min_ga_date = [2020 6 1 0 0 0];
max_ga_date = [2021 2 1 0 0 0];

% Arrival domain - to Mars
min_arr_date = [2021 1 1 0 0 0];
max_arr_date = [2021 8 1 0 0 0];

% Other inputs
mu = astroConstants(4);
R_S = 6955100; %km  0 to remove
R_Me = 2439000.7; %km  000 to remove
R_V = 6051000.8; %km  000 to remove
R_Ma = 3389000.5; %km  000 to remove

%% GA algorithm

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


%% Transfer caracteristics

% Mercury
kep_dep = uplanet(t(1), dep_Id);
[rr_dep, ~] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu);
% Venus
kep_ga = uplanet(t(2), ga_Id);
[rr_ga,~] = kep2car(kep_ga(1), kep_ga(2), kep_ga(3), kep_ga(4), kep_ga(5), kep_ga(6), mu);
% Mars
kep_arr = uplanet(t(3), arr_Id);
[rr_arr,~] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu);

ToF_1 = (t(2) - t(1))*(24*3600);
ToF_2 = (t(3) - t(2))*(24*3600);

% First transfer
[a1,~,~,~,vv_i1, vv_f1,~,~]=lambertMR(rr_dep,rr_ga,ToF_1,mu,0,0,0,2);
% Second transfer
[a2,~,~,~,vv_i2, vv_f2,~,~]=lambertMR(rr_ga,rr_arr,ToF_2,mu,0,0,0,2);


%% Plot

% Time
%T = 2*pi()*sqrt(a1^3/mu); % period first transfer
tfin = ToF_1; % time of parabolic flight first transfer
t0 = 0;
%tspan1 = linspace(0,ToF_1, 800)'; % faster
tspan1 = [t0,tfin]; % slower
%tspan1 = sort(tspan1(:),'descend');

tfin = ToF_2; % time of parabolic flight first transfer
t0 = 0;
tspan2 = [t0,tfin];

% Initial conditions
Y0_1 = [rr_dep(1) rr_dep(2) rr_dep(3) vv_i1(1) vv_i1(2) vv_i1(3)]; % First transfer
Y0_2 = [rr_ga(1) rr_ga(2) rr_ga(3) vv_i2(1) vv_i2(2) vv_i2(3)]; % First transfer

% Set options
options = odeset( 'RelTol', 1e-12, 'AbsTol', 1e-12 );

% Perform the integration
[ T1, Y1 ] = ode113( @(t,y)ode_keplerian_orbit(t,y,mu), tspan1, Y0_1, options);
r1 = [Y1(:,1),Y1(:,2),Y1(:,3)];
v1 = [Y1(:,4),Y1(:,5),Y1(:,6)];

[ T2, Y2 ] = ode113( @(t,y)ode_keplerian_orbit(t,y,mu), tspan2, Y0_2, options);
r2 = [Y2(:,1),Y2(:,2),Y2(:,3)];
v2 = [Y2(:,4),Y2(:,5),Y2(:,6)];

% Plot 
figure(1)
plot(T1,r1(:,1))
title('Trajectory of 1st transfer')


% % Draw planets
figure(2)
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha = 1;

% Sun
image_file_0 = '/Users/morgane/Desktop/Laboratories Orb Mech/Assignment/Assignment_1_ga_fun/textures/Sun.jpg';
[x0, y0, z0] = ellipsoid(0, 0, 0, R_S, R_S, R_S, npanels);
globe1 = surf(x0, y0, z0, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata0 = imread(image_file_0);
set(globe1, 'FaceColor', 'texturemap', 'CData', cdata0, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% Mercury
hold on
image_file_1 = '/Users/morgane/Desktop/Laboratories Orb Mech/Assignment/Assignment_1_ga_fun/textures/Mercury.jpg';
[x1, y1, z1] = ellipsoid(rr_dep(1), rr_dep(2), rr_dep(3), R_Me, R_Me, R_Me, npanels);
globe1 = surf(x1, y1, z1, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata1 = imread(image_file_1);
set(globe1, 'FaceColor', 'texturemap', 'CData', cdata1, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% Venus
hold on
image_file_2 = '/Users/morgane/Desktop/Laboratories Orb Mech/Assignment/Assignment_1_ga_fun/textures/Venus.jpg';
[x2, y2, z2] = ellipsoid(rr_ga(1), rr_ga(2), rr_ga(3), R_V, R_V, R_V, npanels);
globe2 = surf(x2, y2, z2, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata2 = imread(image_file_2);
set(globe2, 'FaceColor', 'texturemap', 'CData', cdata2, 'FaceAlpha', alpha, 'EdgeColor', 'none');

% Mars
hold on
image_file_3 = '/Users/morgane/Desktop/Laboratories Orb Mech/Assignment/Assignment_1_ga_fun/textures/Mars.jpg';
[x3, y3, z3] = ellipsoid(rr_arr(1), rr_arr(2), rr_arr(3), R_Ma, R_Ma, R_Ma, npanels);
globe3 = surf(x3, y3, z3, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
cdata3 = imread(image_file_3);
set(globe3, 'FaceColor', 'texturemap', 'CData', cdata3, 'FaceAlpha', alpha, 'EdgeColor', 'none');


% Aspect
hold on
axis equal
grid on
axis([-2*a1,2*a1,-2*a1,2*a1,-2*a1,2*a1])
xlabel('rx [L]');
ylabel('ry [L]');
zlabel('rz [L]');
title('Orbit');

% % Plot the planets trajectory around the Sun

% Mercury
hold on
T_Me = 2*pi()*sqrt(kep_dep(1)^3/mu);
tMe = linspace(0,T_Me, 800)';
rr_Me = zeros(length(tMe),3);
for i = 1:length(tMe)
    kepMe = uplanet(i, dep_Id);
    [rrMe,~] = kep2car(kepMe(1), kepMe(2), kepMe(3), kepMe(4), kepMe(5), kepMe(6), mu);
    rr_Me(i,:) = transpose(rrMe);
end
plot3(rr_Me(:,1), rr_Me(:,2), rr_Me(:,3));

% Venus
hold on
T_V = 2*pi()*sqrt(kep_ga(1)^3/mu);
tV = linspace(0,T_V, 800)';
rr_V = zeros(length(tV),3);
for i = 1:length(tV)
    kepV = uplanet(i, ga_Id);
    [rrV,~] = kep2car(kepV(1), kepV(2), kepV(3), kepV(4), kepV(5), kepV(6), mu);
    rr_V(i,:) = transpose(rrV);
end
plot3(rr_V(:,1), rr_V(:,2), rr_V(:,3));

% Mars
hold on
T_Ma = 2*pi()*sqrt(kep_arr(1)^3/mu);
tMa = linspace(0,T_Ma, 800)';
rr_Ma = zeros(length(tMa),3);
for i = 1:length(tMa)
    kepMa = uplanet(i, arr_Id);
    [rrMa,~] = kep2car(kepMa(1), kepMa(2), kepMa(3), kepMa(4), kepMa(5), kepMa(6), mu);
    rr_Ma(i,:) = transpose(rrMa);
end
plot3(rr_Ma(:,1), rr_Ma(:,2), rr_Ma(:,3));

% Plot the results
hold on
comet3(r1(:,1),r1(:,2),r1(:,3))
hold on
comet3(r2(:,1),r2(:,2),r2(:,3))
P = [rr_dep,rr_ga];
% p1 = plot3(P(1,1),P(2,1),P(3,1),'.', 'Color', '#0072BD', 'MarkerSize', 20); % blue
% p2 = plot3(P(1,2),P(2,2),P(3,2),'.', 'Color', '#77AC30', 'MarkerSize', 20); % green


