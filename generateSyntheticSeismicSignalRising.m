function generateSyntheticSeismicSignal(filename, duration, a_max, f_dom, zeta, fs, noise_level)
% generateSyntheticSeismicSignal generates a synthetic seismic signal with a rising envelope
%
% INPUT PARAMETERS:
% - filename: output CSV file name (e.g., 'signal.csv')
% - duration: signal duration [s]
% - a_max: peak acceleration [m/s²]
% - f_dom: dominant frequency [Hz]
% - zeta: damping ratio (e.g. 0.05 to 0.3)
% - fs: sampling frequency [Hz]
% - noise_level: noise amplitude as a fraction of a_max (0 to 1)

% Time vector
dt = 1/fs;
t = 0:dt:duration;
omega = 2 * pi * f_dom;

% Rising envelope using a half Hanning window or a logistic function
ramp = sin(pi * t / duration).^2;  % smooth rising from 0 to 1

% Damped oscillation
damped_sine = sin(omega * t) .* exp(-zeta * omega * t);

% Final modulated signal
signal = a_max * ramp .* damped_sine;

% Additive Gaussian noise
noise = noise_level * a_max * randn(size(t));
acc = signal + noise;

% Export to CSV
data = table(t', acc', 'VariableNames', {'Time_s', 'Acceleration_mps2'});
writetable(data, filename);

% Plot
figure;
plot(t, acc, 'b-', 'LineWidth', 1.2); grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s²]');
title('Synthetic Seismic Signal with Rising Envelope and Noise');
legend('Signal + Noise');

end