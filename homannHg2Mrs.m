%% Homann transfer Mercury to Mars
%Both orbits considered circular and coplanar
function [DV, DVtot, DV_circ] = homannHg2Mrs
    date=[2020 1 1 0 0 0];
    dateMJD=date2mjd2000(date);
    
    %% Total Homann DV

    [kep_hg,~] = uplanet(dateMJD, 1);
    [kep_mrs,ksun] = uplanet(dateMJD, 4);
    %Must change date to 1/1/2020

    R_hg=kep_hg(1); %Mercury orbital radius
    R_mrs=kep_mrs(1); %Mars orbital radius
    a_tran = (R_hg+R_mrs)/2; %Semi major axis of transfer orbit

    Vhg = sqrt(ksun/R_hg); %Mercury orbital speed
    Vmrs= sqrt(ksun/R_mrs); %Mars orbital radius

    VtransP = sqrt(ksun * (2/R_hg-1/a_tran)); %Final speed at periapsis
    VtransA = sqrt(ksun * (2/R_mrs-1/a_tran)); %Initial speed at apoapsis

    DV1= VtransP - Vhg; %First burn DeltaV
    DV2= Vmrs - VtransA; %Second burn DeltaV
    DVtot=DV1+DV2; %Total DeltaV
    
    %% Minimum DV to circularize after Venus FlyBy
    %After the Venus FlyBy the orbit must be again modified with a final
    %insertion burn in the Martian orbit, the minimum DV possible is given by
    %a burn at the apocentre from an orbit that has the pericentre at the Venus
    %radius and apocentre at Mars radius.
    %This DV must be subtracted from the maximum initial DV

    [kep_v,~] = uplanet(dateMJD, 2);
    [kep_m,ksun] = uplanet(dateMJD, 4);
    %Planetary orbits considered planar and circular

    Rv=kep_v(1);
    Rm=kep_m(1);

    Vm = sqrt(ksun/Rm);

    a1 = (Rv + Rm) / 2;

    DV_circ = Vm - sqrt(ksun * (2/Rm - 1/a1));

    DV=DVtot-DV_circ;

end
