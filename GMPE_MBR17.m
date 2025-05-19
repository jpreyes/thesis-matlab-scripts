% GMPE_MBR17.m - Script principal
function [Sa, T, errFlag, sig, phi, phiSS, phiS2S, tau] = GMPE_MBR17(R, M, I, Zh, Vs30, F_FABA, epsilon, per)
% R        = distancia (puede ser escalar o vector)
% M        = magnitud (puede ser escalar o vector)
% I        = tipo de evento (1=intraslab, 0=interface)
% Zh       = profundidad hipocentral
% Vs30     = velocidad de la onda de corte promedio de los 30m
% F_FABA   = forearc/backarc (0=forearc, 1=backarc)
% epsilon  = desviaciones estándar
% per      = vector de periodos (s)

% Carga de coeficientes embebidos 
% Extraído de Montalva et al., (2017)
% Cada fila representa un periodo, columnas: [T, c4, C1, DC1, n, c, vlin, b, theta(1-16), phi, tau, sigma, phiS2S]
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
