%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: [Sa, T, errFlag, sig, phi, phiSS, phiS2S, tau] = GMPE_MBR17(...)
%
% SCRIPT NAME: GMPE_MBR17.m
%
% DESCRIPTION:
%   Computes spectral acceleration values (Sa) using the Ground Motion 
%   Prediction Equation (GMPE) proposed by Montalva et al. (2017), calibrated 
%   for Chilean subduction interface and intraslab events. 
%   The implementation supports scalar and vector inputs for distance, 
%   magnitude, and period, and performs log-log interpolation of spectral values.
%
%   This script was developed as part of the doctoral thesis:
%   "Lessons Learned from the Survey of Damage to School Buildings by the 
%    Mw = 8.4 Illapel Earthquake (Chile, September 2015)"
%    Author: Juan Patricio Reyes Cancino
%    Doctoral Program in Earthquake Engineering and Structural Dynamics
%    Universitat Politècnica de Catalunya (UPC) – Universidad Austral de Chile (UACh)
%
% DESCRIPTION:
%
% INPUT:
%   R        - Source-to-site distance [km] (scalar or vector)
%   M        - Moment magnitude (scalar or vector)
%   I        - Event type: 0 = interface, 1 = intraslab
%   Zh       - Hypocentral depth [km]
%   Vs30     - Average shear-wave velocity of top 30 m [m/s]
%   F_FABA   - Forearc/backarc flag: 0 = forearc, 1 = backarc
%   epsilon  - Number of standard deviations (typically 0)
%   per      - Period vector [s]
%
% OUTPUT:
%   Sa       - Spectral acceleration [g]
%   T        - Periods [s]
%   errFlag  - Error flag (0 = OK, -1 = invalid input combination)
%   sig      - Total standard deviation
%   phi      - Intra-event variability
%   phiSS    - Site-to-site variability
%   phiS2S   - Inter-site variability
%   tau      - Inter-event variability
%
%
% AUTHOR: Juan Patricio Reyes
% DATE CREATED: 2024-12-10
% LAST MODIFIED: 2025-04-28
%
% LICENSE: MIT License – See LICENSE_MIT.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sa, T, errFlag, sig, phi, phiSS, phiS2S, tau] = GMPE_MBR17(R, M, I, Zh, Vs30, F_FABA, epsilon, per)

%load('CH_MBR17_CoefTable.mat', 'coef_table'); % asegurarse que coef_table esté en el entorno

errFlag = 0;

% Validaciones básicas
if nargin < 8 || isempty(per)
    per = logspace(log10(0.01), log10(10), 100);
end

% Inicialización
tVec = ~isscalar(per);
rVec = ~isscalar(R);
mVec = ~isscalar(M);

% Error si se dan múltiples vectores simultáneamente
if tVec + rVec + mVec > 2
    Sa = 0;
    T = per;
    errFlag = -1;
    return;
end

% Vector de periodos base
T = per;
rows = length(T);

if tVec && ~rVec && ~mVec
    [T_native, Sa_tmp, sig, phi, phiSS, phiS2S, tau] = CH_MBR17_CoefModel(R, M, I, Zh, Vs30, F_FABA, epsilon);
    Sa = exp(interp1(log(T_native), log(Sa_tmp), log(per), 'linear', 'extrap'));

elseif tVec && rVec && ~mVec
    cols = length(R);
    Sa = zeros(rows, cols);
    for j = 1:cols
        [T_native, Sa_tmp, sig_tmp, phi_tmp, ~, phiS2S_tmp, tau_tmp] = CH_MBR17_CoefModel(R(j), M, I, Zh, Vs30, F_FABA, epsilon);
        Sa(:, j) = exp(interp1(log(T_native), log(Sa_tmp), log(per), 'linear', 'extrap'));
    end

elseif tVec && ~rVec && mVec
    cols = length(M);
    Sa = zeros(rows, cols);
    for j = 1:cols
        [T_native, Sa_tmp, sig_tmp, phi_tmp, ~, phiS2S_tmp, tau_tmp] = CH_MBR17_CoefModel(R, M(j), I, Zh, Vs30, F_FABA, epsilon);
        Sa(:, j) = exp(interp1(log(T_native), log(Sa_tmp), log(per), 'linear', 'extrap'));
    end

else
    rows = length(M);
    cols = length(R);
    Sa = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            [T_native, Sa_tmp, sig_tmp, phi_tmp, ~, phiS2S_tmp, tau_tmp] = CH_MBR17_CoefModel(R(j), M(i), I, Zh, Vs30, F_FABA, epsilon);
            Sa_tmp_interp = exp(interp1(log(T_native), log(Sa_tmp), log(per), 'linear', 'extrap'));
            Sa(i, j) = Sa_tmp_interp;
        end
    end
end

T = per;

% Visualización básica (solo si se llama sin salida)
if nargout == 0
    figure;
    subplot(1,2,1)
    [Sa1,~] = GMPE_MBR17(50, 7.5, 1, 80, 1200, 0, 0, per); % Vs30 alto, R bajo
    semilogx(per, Sa1, 'r', 'LineWidth', 1.5); hold on;
    [Sa2,~] = GMPE_MBR17(50, 7.5, 1, 80, 600, 0, 0, per);
    semilogx(per, Sa2, 'b', 'LineWidth', 1.5);
    xlabel('Periodo (s)'); ylabel('Sa (g)');
    title('Comparación Vs30 alto vs bajo');
    legend('Vs30=1200, R=5km', 'Vs30=600, R=50km');
    grid on

    subplot(1,2,2)
    [Sa3,~] = GMPE_MBR17(5, 7.5, 0, 25, 1100, 1, 0, per);
    semilogx(per, Sa3, 'm', 'LineWidth', 1.5); hold on;
    [Sa4,~] = GMPE_MBR17(5, 7.5, 1, 80, 1100, 1, 0, per);
    semilogx(per, Sa4, 'g', 'LineWidth', 1.5);
    xlabel('Periodo (s)'); ylabel('Sa (g)');
    title('Comparación tipo de evento');
    legend('Interplaca', 'Intraplaca');
    grid on

end
end
