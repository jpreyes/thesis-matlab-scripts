% CH_MBR17_CoefModel.m - Modelo GMPE depurado basado en Montalva et al. (2017)
function [T,Sa,sig,phi,phiSS,phiS2S,tau] = CH_MBR17_CoefModel(R,M,I,Zh,Vs30,F_FABA,epsilon)

% Cargar tabla de coeficientes (matriz coef_table debe estar en workspace)
if I==0
%Se debe cargar previamente con:
%coef_table=readtable('CH_MBR17_CoefTable_full_inter.csv');
%save('CH_MBR17_CoefTableInter.mat','coef_table');
load('CH_MBR17_CoefTableInter.mat','coef_table');
else
%coef_table=readtable('CH_MBR17_CoefTable_full.csv');
%save('CH_MBR17_CoefTable.mat','coef_table');
load('CH_MBR17_CoefTable.mat','coef_table');
end

T       = coef_table.T;
c4      = coef_table.c4(1);
C1      = coef_table.C1(1);
n       = coef_table.n(1);
c       = coef_table.c(1);

F_event = I;
FFABA = F_FABA;

% Truncamiento Vs30
Vs = min(Vs30, 1000);

% Llamar función auxiliar para estimar PGA1000 si Vs30 < vlin
vlin = coef_table.vlin(1);
if Vs30 < vlin
    PGA1000 = CH_MBR17_PGA1000(R, M, I, Zh, F_FABA);
else
    PGA1000 = NaN; % no se usa
end

Sa = zeros(size(T));
sig = Sa; phi = Sa; phiSS = Sa; phiS2S = Sa; tau = Sa;

for i = 1:length(T)
    theta1  = coef_table.theta1(i); theta2  = coef_table.theta2(i); theta3  = coef_table.theta3(i); theta4  = coef_table.theta4(i);
    theta5  = coef_table.theta5(i); theta6  = coef_table.theta6(i); theta7  = coef_table.theta7(i); theta8  = coef_table.theta8(i);
    theta9  = coef_table.theta9(i); theta10 = coef_table.theta10(i); theta11 = coef_table.theta11(i); theta12 = coef_table.theta12(i);
    theta13 = coef_table.theta13(i); theta14 = coef_table.theta14(i); theta15 = coef_table.theta15(i); theta16 = coef_table.theta16(i);
    DC1     = coef_table.DC1(i);
    vlin = coef_table.vlin(i);
    b    = coef_table.b(i);
    % fmag
    if M <= C1 + DC1
        fmag = theta4 * (M - (C1 + DC1));
    else
        fmag = theta5 * (M - (C1 + DC1));
    end

    % fsource
    fsource=theta4*DC1+fmag;

    % fpath
    Rmod = max(R, 1);
    fpath = (theta2 + theta14 * F_event + theta3 * (M - 7.2)) * log(Rmod + c4 * exp(theta9 * (M - 6))) + theta6 * Rmod;

    % fdepth
    fdepth = (theta10 + theta11 * (min(Zh, 120) - 60)) * F_event;

    % fsite usando PGA1000
    if Vs30 < vlin
        fsite = theta12 * log(Vs / vlin) - b * log(PGA1000 + c) + b * log(PGA1000 + c * (Vs / vlin)^n);
    else
        fsite = theta12 * log(Vs / vlin)+b*n*log(Vs/vlin);
    end

    % fFABA
    if F_event == 1
        fFABA = (theta7 + theta8 * log(max(Rmod,85)/40)) * FFABA;
    else
        fFABA = (theta15 + theta16 * log(max(Rmod,100)/40)) * FFABA;
    end

    % Mediana ln(Sa)
    lnSa = theta1 + fsource + fpath + fdepth + fsite + fFABA;

    % Desviaciones
    phi(i)     = coef_table.phi(i);
    tau(i)     = coef_table.tau(i);
    sig(i)     = coef_table.sigma(i);
    phiS2S(i)  = coef_table.phiS2S(i);
    sig1=sqrt(phi(i)^2 + tau(i)^2);
    % Sa mas desviaciones
    Sa(i) = exp(lnSa + epsilon*sig1);
end

phiSS = sqrt(phi.^2 - phiS2S.^2);

end
