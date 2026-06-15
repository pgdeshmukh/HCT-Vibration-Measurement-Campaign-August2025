% VIBRATION DATA ANALYSIS SCRIPT
clear all;
close all;

%% PARAMETERS
data_dir = 'data\T13\Set1';
filename = fullfile(data_dir, 'logfile_11.08.2025, 18.55.27.txt');

g_per_unit = 1 / 256000.0;  % ADXL355 raw to g conversion
g_to_ms2 = 9.80665;

%% LOAD DATA
lines = readlines(filename);
lines = strtrim(lines);
lines = lines(~strcmp(lines, ""));
lines = lines(2:end-1);

data = cellfun(@(L) sscanf(L, '%f,%f,%f,%f').', lines, 'UniformOutput', false);
data = vertcat(data{:});

timestamp_us = data(:, 1);
Ax = data(:, 2) * g_per_unit;
Ay = data(:, 3) * g_per_unit;
Az = data(:, 4) * g_per_unit;

t = (timestamp_us - timestamp_us(1)) / 1e6;
N = length(t);
dt = mean(diff(t));
fs = 1 / dt;

%% RESULTANT
A_res = sqrt(Ax.^2 + Ay.^2 + Az.^2);

%% TIME DOMAIN
figure(1);
subplot(4,1,1); plot(t, Ax, 'r'); ylabel('Ax [g]'); title('Time Domain'); grid on;
subplot(4,1,2); plot(t, Ay, 'g'); ylabel('Ay [g]'); grid on;
subplot(4,1,3); plot(t, Az, 'b'); ylabel('Az [g]'); grid on;
subplot(4,1,4); plot(t, A_res, 'k'); ylabel('A_{res} [g]'); xlabel('Time [s]'); grid on;

%% RMS
A_eff = sqrt(mean(Ax.^2 + Ay.^2 + Az.^2));
A_res_eff = sqrt(mean(A_res.^2));
fprintf('RMS Vector Sum = %.4f g\n', A_eff);
fprintf('RMS Resultant  = %.4f g\n', A_res_eff);

%% FFT
f = (0:N-1)*(fs/N);
window = hann(N);

Ax_fft = abs(fft(Ax .* window)) / (N/2);
Ay_fft = abs(fft(Ay .* window)) / (N/2);
Az_fft = abs(fft(Az .* window)) / (N/2);
Ares_fft = abs(fft(A_res .* window)) / (N/2);

f_plot = f(1:N/2);
Ax_fft = Ax_fft(1:N/2);
Ay_fft = Ay_fft(1:N/2);
Az_fft = Az_fft(1:N/2);
Ares_fft = Ares_fft(1:N/2);

%% FREQUENCY PLOTS
figure(2);
subplot(4,1,1); semilogx(f_plot, Ax_fft, 'r'); ylabel('|Ax|'); xlim([1 fs/2]); grid on;
subplot(4,1,2); semilogx(f_plot, Ay_fft, 'g'); ylabel('|Ay|'); xlim([1 fs/2]); grid on;
subplot(4,1,3); semilogx(f_plot, Az_fft, 'b'); ylabel('|Az|'); xlim([1 fs/2]); grid on;
subplot(4,1,4); semilogx(f_plot, Ares_fft, 'k'); ylabel('|Ares|'); xlabel('Hz'); xlim([1 fs/2]); grid on;

%% DOMINANT FREQUENCIES
valid_idx = f_plot > 1;
f_valid = f_plot(valid_idx);

topN = 3;

% RESULTANT
[Ares_sorted, idx] = sort(Ares_fft(valid_idx), 'descend');
res_freq = f_valid(idx(1:topN));
res_amp  = Ares_sorted(1:topN);

% X
[Ax_sorted, idx] = sort(Ax_fft(valid_idx), 'descend');
Ax_freq = f_valid(idx(1:topN));
Ax_amp  = Ax_sorted(1:topN);

% Y
[Ay_sorted, idx] = sort(Ay_fft(valid_idx), 'descend');
Ay_freq = f_valid(idx(1:topN));
Ay_amp  = Ay_sorted(1:topN);

% Z
[Az_sorted, idx] = sort(Az_fft(valid_idx), 'descend');
Az_freq = f_valid(idx(1:topN));
Az_amp  = Az_sorted(1:topN);

%% PRINT RESULTS (FORMATTED)

fprintf('\nDominant Frequencies and Amplitudes\n');
fprintf('===================================\n\n');

fprintf('RESULTANT:\n');
for i=1:topN
    fprintf('%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, res_freq(i), res_amp(i), res_amp(i)*g_to_ms2);
end

fprintf('\nX AXIS:\n');
for i=1:topN
    fprintf('%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, Ax_freq(i), Ax_amp(i), Ax_amp(i)*g_to_ms2);
end

fprintf('\nY AXIS:\n');
for i=1:topN
    fprintf('%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, Ay_freq(i), Ay_amp(i), Ay_amp(i)*g_to_ms2);
end

fprintf('\nZ AXIS:\n');
for i=1:topN
    fprintf('%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, Az_freq(i), Az_amp(i), Az_amp(i)*g_to_ms2);
end

%% SAVE PLOTS
[~, name, ~] = fileparts(filename);
out_dir = data_dir;

exportgraphics(figure(1), fullfile(out_dir, [name '_time.png']), 'Resolution', 300);
exportgraphics(figure(2), fullfile(out_dir, [name '_freq.png']), 'Resolution', 300);

%% SAVE TEXT FILE
txt_filename = fullfile(out_dir, ['Dominant_Modes.txt']);
fid = fopen(txt_filename, 'w');

fprintf(fid, 'Dominant Frequencies and Amplitudes\n');
fprintf(fid, '===================================\n\n');

fprintf(fid, 'RESULTANT:\n');
for i=1:topN
    fprintf(fid, '%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, res_freq(i), res_amp(i), res_amp(i)*g_to_ms2);
end

fprintf(fid, '\nX AXIS:\n');
for i=1:topN
    fprintf(fid, '%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, Ax_freq(i), Ax_amp(i), Ax_amp(i)*g_to_ms2);
end

fprintf(fid, '\nY AXIS:\n');
for i=1:topN
    fprintf(fid, '%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, Ay_freq(i), Ay_amp(i), Ay_amp(i)*g_to_ms2);
end

fprintf(fid, '\nZ AXIS:\n');
for i=1:topN
    fprintf(fid, '%d) %.2f Hz | %.6f g | %.6f m/s^2\n', ...
        i, Az_freq(i), Az_amp(i), Az_amp(i)*g_to_ms2);
end

fclose(fid);

fprintf('\nSaved TEXT + plots in: %s\n', out_dir);