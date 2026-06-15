% VIBRATION DATA ANALYSIS SCRIPT
% Reads 3-axis accelerometer data from a text file and performs time and frequency domain analysis.
clear all;
close all;

%% PARAMETERS
%data_dir = ['C:\Users\User\OneDrive - ITCC\02 IIA Projects\14 HCT Upgrade\11 MCS Modelling and SImulation\11 Vibration identification\HCT August2025 Test Data\data\T1'];
data_dir = '';

filename = fullfile(data_dir, 'S5_long2.txt');

sampleTime = 0.00040816;        % Sampling interval in seconds C1= C2=2.45Khz
fs = 1/sampleTime;              % Sampling frequency in Hz
g_per_unit = 1/256000.0;        % Scaling factor for ADXL355 raw data

%% DATA LOADING
% Assumes data format: [Ax Ay Az]
data = load(filename);
Ax = data(:, 1) * g_per_unit;
Ay = data(:, 2) * g_per_unit;
Az = data(:, 3) * g_per_unit;

N = length(Ax);                 % Number of samples
t = (0:N-1) * sampleTime;       % Time vector

%% RESULTANT ACCELERATION
A_res = sqrt(Ax.^2 + Ay.^2 + Az.^2);

%% TIME DOMAIN PLOTS
figure;
subplot(4,1,1);
plot(t, Ax, 'r');
ylabel('Ax [g]');
title('Time Domain Vibration Data');
grid on;

subplot(4,1,2);
plot(t, Ay, 'g');
ylabel('Ay [g]');
grid on;

subplot(4,1,3);
plot(t, Az, 'b');
ylabel('Az [g]');
grid on;

subplot(4,1,4);
plot(t, A_res, 'k');
ylabel('A_{res} [g]');
xlabel('Time [s]');
grid on;

%% EFFECTIVE (RMS) ACCELERATION
A_eff = sqrt(mean(Ax.^2 + Ay.^2 + Az.^2));
A_res_eff = sqrt(mean(A_res.^2));
fprintf('Effective (RMS) Vector Sum Acceleration = %.4f g\n', A_eff);
fprintf('Effective (RMS) Resultant Acceleration = %.4f g\n', A_res_eff);

%% FREQUENCY DOMAIN ANALYSIS
f = (0:N-1)*(fs/N);   % Frequency vector

% Apply windowing
window = hann(N);
Ax_win = Ax .* window;
Ay_win = Ay .* window;
Az_win = Az .* window;
Ares_win = A_res .* window;

% Compute FFTs
Ax_fft = abs(fft(Ax_win)) / (N/2);
Ay_fft = abs(fft(Ay_win)) / (N/2);
Az_fft = abs(fft(Az_win)) / (N/2);
Ares_fft = abs(fft(Ares_win)) / (N/2);

% Keep only positive frequencies
f_plot = f(1:N/2);
Ax_fft = Ax_fft(1:N/2);
Ay_fft = Ay_fft(1:N/2);
Az_fft = Az_fft(1:N/2);
Ares_fft = Ares_fft(1:N/2);

%% FREQUENCY DOMAIN PLOTS
fs=1428;
figure;
subplot(4,1,1);
semilogx(f_plot, Ax_fft, 'r');
xlim([1 fs/2]);
ylabel('|Ax(f)| [g]');
title('Frequency Domain (Log-Log) Spectrum');
grid on;

subplot(4,1,2);
semilogx(f_plot, Ay_fft, 'g');
xlim([1 fs/2]);
ylabel('|Ay(f)| [g]');
grid on;

subplot(4,1,3);
semilogx(f_plot, Az_fft, 'b');
xlim([1 fs/2]);
ylabel('|Az(f)| [g]');
grid on;

subplot(4,1,4);
semilogx(f_plot, Ares_fft, 'k');
xlim([1 fs/2]);
ylabel('|A_{res}(f)| [g]');
xlabel('Frequency [Hz]');
grid on;


%% SAVE PLOTS (300 dpi, PNG)
[~, name, ~] = fileparts(filename);   % extract file name without extension
out_dir = data_dir;                   % save in same folder

% Time domain plot
exportgraphics(figure(1), fullfile(out_dir, [name '_time.png']), 'Resolution', 300);

% Frequency domain plot
exportgraphics(figure(2), fullfile(out_dir, [name '_freq.png']), 'Resolution', 300);