function [deltaV_tot] = f (t_dep, t_ga, t_arr)

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


ToF_1 = (t_ga - t_dep)*(24*3600);
ToF_2 = (t_arr - t_ga)*(24*3600);

[~,~,~,~,vv_i1, vv_f1,~,~]=lambertMR(rr_dep,rr_ga,ToF_1,mu,0,0,0,2);

deltaV_1 = (norm(vv_f1 - vv_ga) + norm(vv_i1 - vv_dep));

[~,~,~,~,vv_i2, vv_f2,~,~]=lambertMR(rr_ga,rr_arr,ToF_2,mu,0,0,0,2);

deltaV_2 = (norm(vv_f2 - vv_arr) + norm(vv_i2 - vv_ga));

deltaV_tot = deltaV_1 + deltaV_2;

