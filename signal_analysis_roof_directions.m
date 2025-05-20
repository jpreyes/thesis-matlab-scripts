
% signal_analysis_roof_directions_commented.m
%
% This script was developed as part of the doctoral thesis:
% "Lessons Learned from the Survey of Damage to School Buildings by the Mw = 8.4 Illapel Earthquake (Chile, September 2015)"
% Author: Juan Patricio Reyes Cancino
% PhD Program in Earthquake Engineering and Structural Dynamics
% UPC – Universidad Austral de Chile
%
% DESCRIPTION:
% This script reads acceleration response signals measured at the roof level of a school building
% in three directions: X (longitudinal), Y (transverse), and RZ (rotational).
% It performs:
%   - Time-history plotting
%   - Short-Time Fourier Transform (spectrograms)
%   - Parametric spectral estimation using Yule-Walker method
%   - Classical power spectrum via periodograms
%
% REQUIRED INPUT FILES (ASCII or CSV format):
%   - 'roof_x.txt' or 'roof_x.csv'   → two columns: [time, acceleration_x]
%   - 'roof_y.txt' or 'roof_y.csv'   → two columns: [time, acceleration_y]
%   - 'roof_rz.txt' or 'roof_rz.csv' → two columns: [time, acceleration_rz]
%
% The analysis supports the dynamic assessment of modal energy distribution and damage-related frequency shifts.

% Sampling frequency [Hz]
fm = 10;

% Load seismic response in X direction (roof)
fid = fopen('roof_x.txt');
p1 = fscanf(fid, '%g %g', [2 inf]);
fclose(fid);

% Load seismic response in Y direction (roof)
fid = fopen('roof_y.txt');
p2 = fscanf(fid, '%g %g', [2 inf]);
fclose(fid);

% Load seismic response in RZ direction (roof)
fid = fopen('roof_rz.txt');
p3 = fscanf(fid, '%g %g', [2 inf]);
fclose(fid);

% Format: transpose and extract acceleration components
p1 = p1'; p2 = p2'; p3 = p3';
p1s = p1(:,2); p2s = p2(:,2); p3s = p3(:,2);

%% Time-domain plot of the X direction
figure(14)
plot(p1(:,1), p1(:,2))
title('Response Signal in X Direction')
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');

%% Spectrogram (Short-Time Fourier Transform) for X direction
[y,f,t,p] = spectrogram(p1s,500,490,[],fm,'yaxis');
figure(15)
contour(t,f,p,200)
title('Spectrogram of the Response Signal in X Direction')
xlabel('Time [s]');
ylabel('Frequency [Hz]');

figure(30)
surf(t,f,10*log10(abs(p)),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);
title('Spectrogram of the Response Signal in X Direction')
xlabel('Time [s]');
ylabel('Frequency [Hz]');

%% Parametric spectral estimation (Yule-Walker) for X direction
figure(16)
pyulear(p1s,50,[],fm)

%% Periodogram (power spectral density) for X direction
[ps2, pf2] = periodogram(p1s, [], [], fm);
figure(17)
plot(pf2, ps2)
title('Periodogram Response Signal in X Direction')
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');

%% Time-domain plot of the Y direction
figure(18)
plot(p2(:,1), p2(:,2))
title('Response Signal in Y Direction')
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');

%% Spectrogram for Y direction
[y,f,t,p] = spectrogram(p2s,500,490,[],fm,'yaxis');
figure(19)
contour(t,f,p,200)
title('Spectrogram of the Response Signal in Y Direction')
xlabel('Time [s]');
ylabel('Frequency [Hz]');

figure(31)
surf(t,f,10*log10(abs(p)),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);
title('Spectrogram of the Response Signal in Y Direction')
xlabel('Time [s]');
ylabel('Frequency [Hz]');

%% Yule-Walker for Y direction
figure(20)
pyulear(p2s,90,[],fm)

%% Periodogram for Y direction
figure(21)
[ps3, pf3] = periodogram(p2s,[],[],fm);
plot(pf3,ps3)
title('Periodogram Response Signal in Y Direction')
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');

%% Time-domain plot of the RZ (rotational) direction
figure(22)
plot(p3(:,1), p3(:,2))
title('Response Signal in RZ Direction')
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');

%% Spectrogram for RZ direction
[y,f,t,p] = spectrogram(p3s,500,490,[],fm,'yaxis');
figure(23)
contour(t,f,p,200)
title('Spectrogram of the Response Signal in RZ Direction')
xlabel('Time [s]');
ylabel('Frequency [Hz]');

figure(32)
surf(t,f,10*log10(abs(p)),'EdgeColor','none');
axis xy; axis tight; colormap(jet); view(0,90);
title('Spectrogram of the Response Signal in RZ Direction')
xlabel('Time [s]');
ylabel('Frequency [Hz]');

%% Yule-Walker for RZ direction
figure(24)
pyulear(p3s,90,[],fm)

%% Periodogram for RZ direction
figure(25)
[ps4, pf4] = periodogram(p3s,[],[],fm);
plot(pf4,ps4)
title('Periodogram Response Signal in RZ Direction')
xlabel('Frequency [Hz]');
ylabel('Amplitude [dB]');

%% Additional comparisons (early vs late sections of signal)
% Short-term segment Yule-Walker comparison
figure(40); pyulear(p1s(1500:2001),5,[],fm)
figure(41); pyulear(p2s(1500:2001),5,[],fm)
figure(42); pyulear(p3s(1500:2001),5,[],fm)

% Overlay Yule-Walker spectra: initial vs final segments
figure(51); pyulear(p1s(1:500),50,[],fm); hold on; pyulear(p1s(1500:2001),50,[],fm);
figure(52); pyulear(p2s(1:500),50,[],fm); hold on; pyulear(p2s(1500:2001),50,[],fm);
figure(53); pyulear(p3s(1:500),50,[],fm); hold on; pyulear(p3s(1500:2001),50,[],fm);
