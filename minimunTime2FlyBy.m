%% Minimum and Maximum Time to FlyBy
clear all
%Given a maximum DeltaV this function calculates the minimum and mazimum 
%time required to arrive to Venus as the two intersections of the maximum
%transfer ellipse from Mercury.

%All orbits are assumed coplanar, Mercury eccentricity is taken into
%account, impulse is assumed at Mercury's periapsis, Venus is considered at
%its minimum distance from the Sun at all times

date=[2020 1 1 0 0 0];
dateMJD=date2mjd2000(date);

[kep_hg,ksun] = uplanet(dateMJD, 1);

a=kep_hg(1); %Semi-Major Axis of Mercury
e=kep_hg(2); %Eccentricity of Mercury

Rp_hg = a*(1-e); %Periapsis radius of Mercury
Vp_hg = sqrt(ksun*(1+e)/Rp_hg); %Periapsis speed of Mercury

[DV, DVtot, DV_circ] = homannHg2Mrs; %DV is maximum DeltaV possible
%First guess is the Homann transfer between Mercury and Mars

a_tr = 1/( 2/Rp_hg - (DV+Vp_hg)^2/ksun ); %Semi-Major Axis of transfer orbit
e_tr = 1 - Rp_hg/a_tr; %Eccentricity of transfer orbit

[kep_ve,~] = uplanet(dateMJD, 2);
r_ve = kep_ve(1)*(1-kep_ve(2)); %Venus minimum distance from Sun

%% Minimum time to FlyBy

th = acos( ( a_tr*(1-e_tr^2)/r_ve - 1) / e_tr ); %True anomaly of maximum ellipse when it intersects minimum Venus distance
E = 2*atan( tan(th/2) / sqrt((1+e_tr)/(1-e_tr)) ); %Eccentric anomaly
t = (E - e_tr*sin(E)) / sqrt( ksun/(a_tr^3)); %Time to that anomaly
    
secondsAday=23*3600+56*60+4.1;
tdays= t/secondsAday

%% Maximum time to FlyBy

tmax = sqrt((a_tr^3)/ksun) - 2*t;
tmax_days = tmax / secondsAday





