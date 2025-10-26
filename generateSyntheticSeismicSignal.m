% generateSyntheticSeismicSignal.m
%
% MATLAB function developed as part of the doctoral thesis:
% "Lessons Learned from the Survey of Damage to School Buildings by the Mw = 8.4 Illapel Earthquake (Chile, September 2015)"
% Author: Juan Patricio Reyes Cancino
% Doctoral Program in Earthquake Engineering and Structural Dynamics
% Universitat Politècnica de Catalunya (UPC) – Universidad Austral de Chile (UACh)
% Date: june, 2025 – Version 1.0
%
% DESCRIPTION:
% This function generates a synthetic seismic acceleration signal with a Gaussian-modulated 
% envelope and superimposed white noise, intended for sensitivity analyses and 
% seismic input generation in nonlinear dynamic simulations. The function allows for 
% parametric control over signal duration, peak ground acceleration (PGA), dominant frequency, 
% damping ratio, noise amplitude, and the time location of the envelope peak.
%
% The signal represents a damped harmonic motion modulated by a smooth pulse centered 
% at a specified fraction of the total time, allowing for physically realistic temporal concentration 
% of seismic energy. White Gaussian noise is added to approximate high-frequency content and 
% signal variability.
%
% INPUT PARAMETERS:
% - filename        : String, name of the output CSV file
% - duration        : Total signal duration [s]
% - a_max           : Peak acceleration [m/s²]
% - f_dom           : Dominant frequency of the signal [Hz]
% - zeta            : Critical damping ratio (dimensionless, usually 0.03–0.05)
% - fs              : Sampling frequency [Hz]
% - noise_level     : Relative amplitude of added white noise (0 to 1, scaled to a_max)
% - peak_fraction   : Relative time position of the signal peak (e.g., 0.2 for 20% of total duration)
%
% OUTPUT:
% - A CSV file with two columns: time [s] and acceleration [m/s²]
% - A time-history plot showing the generated signal
%
% APPLICATION:
% The generated signals can be used in sensitivity studies, parametric response analyses,
% and for validation of seismic energy-based methodologies. In this thesis, they are
% used to isolate the influence of input characteristics on the nonlinear seismic response
% of reinforced concrete frames.
%
% See Appendix D of the thesis for further methodological details.

function generateSyntheticSeismicSignal(filename, duration, a_max, f_dom, zeta, fs, noise_level, peak_fraction)
% Time vector
dt = 1/fs;
t = 0:dt:duration;
omega = 2 * pi * f_dom;

% Define Gaussian-like pulse envelope centered at desired peak
center = peak_fraction * duration;
width = 0.4 * duration;
sigma = width / 6;  % standard deviation for smooth pulse shape

% Envelope function: Gaussian bell centered at 'center'
envelope = exp(-0.5 * ((t - center) / sigma).^2);
envelope = envelope / max(envelope);  % normalize to peak = 1

% Damped harmonic wave modulated by envelope
signal = a_max * envelope .* sin(omega * t) .* exp(-zeta * omega * t);

% Add white Gaussian noise
noise = noise_level * a_max * randn(size(t));
acc = signal + noise;

% Export to CSV
data = table(t', acc', 'VariableNames', {'Time_s', 'Acceleration_mps2'});
writetable(data, filename);

% Plot
figure;
plot(t, acc, 'LineWidth', 1.2); grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');
title(['Synthetic Seismic Signal – Pulse Envelope (Peak at ', num2str(peak_fraction * 100), '%)']);
legend('Signal + Noise');
end

