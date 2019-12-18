function [deltaV_tot] = f (t)

t_dep = t(1);
t_ga = t(2);
t_arr = t(3);

mu = astroConstants(4);
dep_Id = 1; % Mercury
ga_Id = 2; % Venus
arr_Id = 4; % Mars

kep_dep = uplanet(t_dep, dep_Id);
[rr_dep, vv_dep] = kep2car(kep_dep(1), kep_dep(2), kep_dep(3), kep_dep(4), kep_dep(5), kep_dep(6), mu);


kep_ga = uplanet(t_ga, ga_Id);
[rr_ga,vv_ga] = kep2car(kep_ga(1), kep_ga(2), kep_ga(3), kep_ga(4), kep_ga(5), kep_ga(6), mu);


kep_arr = uplanet(t_arr, arr_Id);
[rr_arr,vv_arr] = kep2car(kep_arr(1), kep_arr(2), kep_arr(3), kep_arr(4), kep_arr(5), kep_arr(6), mu);


ToF_1 = (t_ga - t_dep)*86400;
ToF_2 = (t_arr - t_ga)*86400;


% DeltaV first transfer - Mercury to Venus
[~,~,~,~,vv_i1, vv_f1,~,~]=lambertMR(rr_dep,rr_ga,ToF_1,mu,0,0,0,2);

% DeltaV second transfer - Venus to Mars
[~,~,~,~,vv_i2, vv_f2,~,~]=lambertMR(rr_ga,rr_arr,ToF_2,mu,0,0,0,2);


% Tp1 = 2*pi()*sqrt(a1^3/mu);
% Tp2 = 2*pi()*sqrt(a2^3/mu);


% Fly By - Venus

vv_min = vv_ga + vv_f1;   % v_inf_min = vv_f1;
vv_plus = vv_ga + vv_i2;  % v_inf_plus = vv_i2;

deltaV_ga = norm(vv_plus - vv_min); 


deltaV_1 = norm(vv_i1 - vv_dep) + norm(vv_ga - vv_f1);    % the velocity minus the velocity of the planet
deltaV_2 = norm(vv_i2 - vv_ga) + norm(vv_arr - vv_f2);  

% if ToF_1 > (Tp1 - ToF_1)
%     deltaV_tot = Nan;
% elseif ToF_2 > (Tp2 - ToF_2)
%     deltaV_tot = Nan;
% else
%     deltaV_tot = deltaV_1 + deltaV_ga + deltaV_2;
% end

deltaV_tot = deltaV_1 + deltaV_ga + deltaV_2;

return
