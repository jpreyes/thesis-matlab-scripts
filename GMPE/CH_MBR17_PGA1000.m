%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: PGA1000 = CH_MBR17_PGA1000(R, M, I, Zh, F_FABA)
%
% SCRIPT NAME: CH_MBR17_PGA1000.m
%
% DESCRIPTION:
%   Computes the spectral acceleration at T = 0.01 s (PGA equivalent)
%   assuming a reference shear-wave velocity of Vs30 = 1000 m/s, based on
%   the GMPE model of Montalva et al. (2017) for Chilean subduction earthquakes.
%   This function is used to estimate site amplification in cases where 
%   Vs30 is below the linear threshold (vlin), following a two-stage approach.
%
%   The function loads the appropriate coefficient set depending on whether 
%   the event is classified as interface or intraslab.
%
%   This subroutine supports the full CH-MBR17 ground motion prediction
%   framework, and is called internally by CH_MBR17_CoefModel.
%
%   Developed as part of the doctoral thesis:
%   "Lessons Learned from the Survey of Damage to School Buildings by the 
%    Mw = 8.4 Illapel Earthquake (Chile, September 2015)"
%    Author: Juan Patricio Reyes Cancino
%    Doctoral Program in Earthquake Engineering and Structural Dynamics
%    Universitat Politècnica de Catalunya (UPC) – Universidad Austral de Chile (UACh)
%
% DESCRIPTION:
%
% INPUT:
%   R        - Source-to-site distance [km]
%   M        - Moment magnitude
%   I        - Event type (0 = interface, 1 = intraslab)
%   Zh       - Hypocentral depth [km]
%   F_FABA   - Tectonic region indicator (0 = forearc, 1 = backarc)
%
% OUTPUT:
%   PGA1000  - Pseudo-acceleration [g] at T = 0.01 s for Vs30 = 1000 m/s
%
% DEPENDENCIES:
%   Requires coefficient tables:
%   - CH_MBR17_CoefTable.mat (intraslab)
%   - CH_MBR17_CoefTableInter.mat (interface)
%
% LICENSE: MIT License – See LICENSE_MIT.txt
%
% AUTHOR: Juan Patricio Reyes
% DATE CREATED: 2025-05-17
% LAST MODIFIED: 2025-05-19
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CH_MBR17_PGA1000.m - Calcula PGA (Sa en T=0.01s) con Vs30=1000 m/s
function PGA1000 = CH_MBR17_PGA1000(R, M, I, Zh, F_FABA)

% Cargar coeficientes
if I==0
load('CH_MBR17_CoefTableInter.mat','coef_table');
else
load('CH_MBR17_CoefTable.mat','coef_table');
end

T = coef_table.T;
c4   = coef_table.c4(1);
C1   = coef_table.C1(1);
n    = coef_table.n(1);
c    = coef_table.c(1);

Vs = 1000;
F_event = I;
FFABA = F_FABA;

Sa = zeros(size(T));

for i = 1:length(T)
    theta1  = coef_table.theta1(i); theta2  = coef_table.theta2(i); theta3  = coef_table.theta3(i); theta4  = coef_table.theta4(i);
    theta5  = coef_table.theta5(i); theta6  = coef_table.theta6(i); theta7  = coef_table.theta7(i); theta8  = coef_table.theta8(i);
    theta9  = coef_table.theta9(i); theta10 = coef_table.theta10(i); theta11 = coef_table.theta11(i); theta12 = coef_table.theta12(i);
    theta13 = coef_table.theta13(i); theta14 = coef_table.theta14(i); theta15 = coef_table.theta15(i); theta16 = coef_table.theta16(i);
    DC1  = coef_table.DC1(i);
    vlin = coef_table.vlin(i);
    b    = coef_table.b(i);
    if M <= C1 + DC1
        fmag = theta4 * (M - (C1 + DC1));
    else
        fmag = theta5 * (M - (C1 + DC1));
    end

    fsource=theta4*DC1+fmag;

    Rmod = max(R, 1);
    fpath = (theta2 + theta14 * F_event + theta3 * (M - 7.2)) * log(Rmod + c4 * exp(theta9 * (M - 6))) + theta6 * Rmod;

    fdepth = (theta10 + theta11 * (min(Zh, 120) - 60)) * F_event;
    fsite = theta12 * log(Vs / vlin)+b*n*log(Vs/vlin);

    if F_event == 1
        fFABA = (theta7 + theta8 * log(max(Rmod,85)/40)) * FFABA;
    else
        fFABA = (theta15 + theta16 * log(max(Rmod,100)/40)) * FFABA;
    end

    lnSa = theta1 + fsource + fpath + fdepth + fsite + fFABA;
    Sa(i) = exp(lnSa);
end

% Devolver PGA para T = 0.01 s
PGA1000 = interp1(T, Sa, 0.01, 'linear', 'extrap');

end
