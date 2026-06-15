% VIBRATION DATA ANALYSIS - BATCH MODE WITH DATA QUALITY CHECK
clear all;
close all;
clc;

%% PARAMETERS
data_dir = 'data\T13\set2';
g_per_unit = 1 / 256000.0;
g_to_ms2 = 9.80665;
topN = 5;

files = dir(fullfile(data_dir, 'logfile_*.txt'));

txt_filename = fullfile(data_dir, 'Dominant_Modes.txt');
fid = fopen(txt_filename, 'w');

fprintf(fid, 'Dominant Frequencies and Amplitudes\n');
fprintf(fid, '===================================\n\n');

%% LOOP THROUGH FILES
for k = 1:length(files)

    filename = fullfile(data_dir, files(k).name);
    fprintf('\nProcessing file: %s\n', files(k).name);

    %% LOAD DATA
    lines = readlines(filename);
    lines = strtrim(lines);
    lines = lines(~strcmp(lines, ""));
    lines = lines(2:end-1);

    data = cellfun(@(L) sscanf(L, '%f,%f,%f,%f').', lines, 'UniformOutput', false);
    data = vertcat(data{:});

    timestamp_us = data(:,1);
    Ax = data(:,2) * g_per_unit;
    Ay = data(:,3) * g_per_unit;
    Az = data(:,4) * g_per_unit;

    %% =========================
    %% DATA QUALITY CHECK
    %% =========================
    dt_all = diff(timestamp_us);
    dt_median = median(dt_all);
    threshold = 5 * dt_median;

    jump_idx = find(dt_all > threshold);

    if isempty(jump_idx)
        quality_flag = 'GOOD';
    else
        quality_flag = sprintf('DISCONTINUOUS (%d gaps)', length(jump_idx));
    end

    %% SEGMENTATION (use longest continuous segment)
    segment_edges = [1; jump_idx+1; length(timestamp_us)];
    num_segments = length(segment_edges)-1;

    max_len = 0;
    best_seg = 1;

    for s = 1:num_segments
        seg_len = segment_edges(s+1) - segment_edges(s);
        if seg_len > max_len
            max_len = seg_len;
            best_seg = s;
        end
    end

    idx_start = segment_edges(best_seg);
    idx_end   = segment_edges(best_seg+1)-1;

    timestamp_us = timestamp_us(idx_start:idx_end);
    Ax = Ax(idx_start:idx_end);
    Ay = Ay(idx_start:idx_end);
    Az = Az(idx_start:idx_end);

    %% TIME VECTOR
    t = (timestamp_us - timestamp_us(1)) / 1e6;
    N = length(t);
    dt = mean(diff(t));
    fs = 1 / dt;

    %% RESULTANT
    A_res = sqrt(Ax.^2 + Ay.^2 + Az.^2);

    %% TIME DOMAIN PLOT
    figure(1); clf;
    subplot(4,1,1); plot(t, Ax, 'r'); ylabel('Ax [g]'); grid on;
    subplot(4,1,2); plot(t, Ay, 'g'); ylabel('Ay [g]'); grid on;
    subplot(4,1,3); plot(t, Az, 'b'); ylabel('Az [g]'); grid on;
    subplot(4,1,4); plot(t, A_res, 'k'); ylabel('A_{res} [g]'); xlabel('Time [s]'); grid on;
    sgtitle([files(k).name ' | ' quality_flag], 'Interpreter', 'none');

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

    %% FREQUENCY PLOT
    figure(2); clf;
    subplot(4,1,1); semilogx(f_plot, Ax_fft, 'r'); ylabel('|Ax|'); xlim([1 fs/2]); grid on;
    subplot(4,1,2); semilogx(f_plot, Ay_fft, 'g'); ylabel('|Ay|'); xlim([1 fs/2]); grid on;
    subplot(4,1,3); semilogx(f_plot, Az_fft, 'b'); ylabel('|Az|'); xlim([1 fs/2]); grid on;
    subplot(4,1,4); semilogx(f_plot, Ares_fft, 'k'); ylabel('|Ares|'); xlabel('Hz'); xlim([1 fs/2]); grid on;
    sgtitle([files(k).name ' | ' quality_flag], 'Interpreter', 'none');

    %% DOMINANT FREQUENCIES
    valid_idx = f_plot > 1;
    f_valid = f_plot(valid_idx);

    [Ares_sorted, idx] = sort(Ares_fft(valid_idx), 'descend');
    res_freq = f_valid(idx(1:topN));
    res_amp  = Ares_sorted(1:topN);

    [Ax_sorted, idx] = sort(Ax_fft(valid_idx), 'descend');
    Ax_freq = f_valid(idx(1:topN));
    Ax_amp  = Ax_sorted(1:topN);

    [Ay_sorted, idx] = sort(Ay_fft(valid_idx), 'descend');
    Ay_freq = f_valid(idx(1:topN));
    Ay_amp  = Ay_sorted(1:topN);

    [Az_sorted, idx] = sort(Az_fft(valid_idx), 'descend');
    Az_freq = f_valid(idx(1:topN));
    Az_amp  = Az_sorted(1:topN);

    %% PRINT RESULTS
    fprintf('\nFile: %s\n', files(k).name);
    fprintf('Data Quality: %s\n', quality_flag);

    print_axis('RESULTANT', res_freq, res_amp, g_to_ms2, topN);
    print_axis('X AXIS', Ax_freq, Ax_amp, g_to_ms2, topN);
    print_axis('Y AXIS', Ay_freq, Ay_amp, g_to_ms2, topN);
    print_axis('Z AXIS', Az_freq, Az_amp, g_to_ms2, topN);

    %% SAVE PLOTS
    [~, name, ~] = fileparts(files(k).name);
    exportgraphics(figure(1), fullfile(data_dir, [name '_time.png']), 'Resolution', 300);
    exportgraphics(figure(2), fullfile(data_dir, [name '_freq.png']), 'Resolution', 300);

    %% WRITE TO TEXT FILE
    fprintf(fid, 'File: %s\n', files(k).name);
    fprintf(fid, 'Data Quality: %s\n', quality_flag);
    fprintf(fid, '-----------------------------------\n');

    write_axis(fid, 'RESULTANT', res_freq, res_amp, g_to_ms2, topN);
    write_axis(fid, 'X AXIS', Ax_freq, Ax_amp, g_to_ms2, topN);
    write_axis(fid, 'Y AXIS', Ay_freq, Ay_amp, g_to_ms2, topN);
    write_axis(fid, 'Z AXIS', Az_freq, Az_amp, g_to_ms2, topN);

    fprintf(fid, '\n===================================\n\n');

end

fclose(fid);

fprintf('\n✅ Batch processing complete. Results saved in: %s\n', txt_filename);

%% =========================
%% HELPER FUNCTIONS (MUST BE AT END)
%% =========================

function label = classify_amp(amp)
    if amp > 1e-3
        label = 'HIGH';
    elseif amp > 1e-4
        label = 'MODERATE';
    elseif amp > 2e-5
        label = 'LOW';
    else
        label = 'Noise Floor';
    end
end

function print_axis(name, freq, amp, g_to_ms2, topN)
    fprintf('\n%s:\n', name);
    for i=1:topN
        label = classify_amp(amp(i));
        fprintf('%d) %.2f Hz | %.6f g | %.6f m/s^2 | %s\n', ...
            i, freq(i), amp(i), amp(i)*g_to_ms2, label);
    end
end

function write_axis(fid, name, freq, amp, g_to_ms2, topN)
    fprintf(fid, '\n%s:\n', name);
    for i=1:topN
        label = classify_amp(amp(i));
        fprintf(fid, '%d) %.2f Hz | %.6f g | %.6f m/s^2 | %s\n', ...
            i, freq(i), amp(i), amp(i)*g_to_ms2, label);
    end
end