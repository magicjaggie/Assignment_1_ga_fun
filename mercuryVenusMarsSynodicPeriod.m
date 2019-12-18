%% Mercury-Venus-Mars Synodic Period
%Calculating total average error on the difference of the three repective 
%true anomalies for the three planets in function of time 

date0=[2020 1 1 0 0 0];
dateMJD0=date2mjd2000(date0);

datef=[2060 1 1 0 0 0];
dateMJDf=date2mjd2000(datef);

[kep_hg0,~] = uplanet(dateMJD0, 1);
[kep_ve0,~] = uplanet(dateMJD0, 2);
[kep_ma0,~] = uplanet(dateMJD0, 4);

dth10 = (wrapToPi(kep_hg0(6)- kep_ve0(6)));
dth20 = (wrapToPi(kep_ve0(6)- kep_ma0(6)));
%dth30 = (wrapToPi(kep_hg0(6)- kep_ma0(6)));

period = dateMJDf-dateMJD0;
err=zeros(1,(period+1));
err1=zeros(1,(period+1));
err2=zeros(1,(period+1));

% abs(wtp(th1))-abs(wtp(th2))

for i=0:period
    
    [kep_hg,~] = uplanet((dateMJD0+i), 1);
    [kep_ve,~] = uplanet((dateMJD0+i), 2);
    [kep_ma,~] = uplanet((dateMJD0+i), 4);
    
    dth1 = (wrapToPi(kep_hg(6)- kep_ve(6)));
    dth2 = (wrapToPi(kep_ve(6)- kep_ma(6)));
    
    err1(i+1)= abs(wrapToPi(dth1-dth10));
    err2(i+1)= abs(wrapToPi(dth2-dth20));
    
    err(i+1) = sqrt( (err1(i+1))^2 + (err2(i+1))^2 );
    
end

hold on
grid MINOR
plot((0:period),err,'k')
%plot((0:period),err,'k',(0:period),err1,'y',(0:period),err2,'r')