% clear; clc; close all;

%% Read file with timestamps
fname = '70cm_0.00Hz_120sec_2_260518_mouse.txt';   % Change to your filename
L = readlines(fname);
L = L(strlength(L)>0);  % Remove empty lines

% Header
hdr_RX1_I = uint8(hex2dec(split('55 AA 00 91 01 04 00 00 00 00')));
hdr_RX1_Q = uint8(hex2dec(split('55 AA 00 91 01 04 00 01 00 00')));
FRAME_LEN  = 10 + 128;  
% Parameters
N   = 32;                 % Fast-time samples per frame
T_meas = 120;            
B   = 200e6;              % 200 MHz
c   = 3e8;
range_resolution = c/(2*B);   % 0.75 m
REMOVE_MEAN = true;

% Buffers
I_RX1_all = []; Q_RX1_all = [];
t_I = datetime.empty(0,1); t_Q = datetime.empty(0,1);
last_ts = NaT;

% Little-endian 4 bytes -> single
convert_bytes_to_float = @(b4) typecast( ...
  uint32( bitor( bitor(bitshift(uint32(b4(4)),24), bitshift(uint32(b4(3)),16)), ...
                 bitor(bitshift(uint32(b4(2)), 8), uint32(b4(1))) ) ), 'single');

i = 1;
while i <= numel(L)
    line = strtrim(L(i));

    % 1) Timestamp line
    if startsWith(line,"==== Time:")
        ts_text = extractBetween(line, "==== Time:", "====");
        if ~isempty(ts_text)
            ts_text = strtrim(ts_text{1});
            try
                last_ts = datetime(ts_text,'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            catch
                last_ts = datetime(ts_text,'InputFormat','yyyy-MM-dd HH:mm:ss');
            end
        else
            last_ts = NaT;
        end
        i = i+1;
        continue
    end

    % 2) Frame data line
    toks = split(line);
    toks = toks(strlength(toks) == 2);
    if numel(toks) ~= FRAME_LEN
        i = i+1;  
        continue
    end
    bytes = uint8(hex2dec(toks));

    header  = bytes(1:10);
    payload = bytes(11:end);   

    % 3) 128 bytes -> 32 singles
    vals = zeros(N,1,'single');
    for k = 1:N
        b4 = payload((4*k-3):(4*k));
        vals(k) = convert_bytes_to_float(b4);
    end

    % 4) Route to I or Q channel
    if isequal(header, hdr_RX1_I)
        I_RX1_all(:, end+1) = vals; 
        t_I(end+1,1) = last_ts;     
    elseif isequal(header, hdr_RX1_Q)
        Q_RX1_all(:, end+1) = vals; 
        t_Q(end+1,1) = last_ts;    
    end

    i = i+1;
end

% 5) Align I/Q to common length
T_I = size(I_RX1_all,2);
T_Q = size(Q_RX1_all,2);
T   = min(T_I, T_Q);
if T == 0
    error('No complete 138-byte frames found in file (or I/Q count is zero for one channel).');
end
if T_I ~= T_Q
    warning('RX1 I/Q mismatch: I=%d, Q=%d -> trim to %d', T_I, T_Q, T);
end
I_RX1_all = I_RX1_all(:,1:T);
Q_RX1_all = Q_RX1_all(:,1:T);
t_I = t_I(1:T); t_Q = t_Q(1:T);

%% DC removal
if REMOVE_MEAN
    I_RX1_dc = I_RX1_all - mean(I_RX1_all, 1);
    Q_RX1_dc = Q_RX1_all - mean(Q_RX1_all, 1);
else
    I_RX1_dc = I_RX1_all;
    Q_RX1_dc = Q_RX1_all;
end

%% Build complex signal, apply Range FFT
IQ_RX1 = I_RX1_dc - 1j*Q_RX1_dc;    % N x T
RFFT   = fft(IQ_RX1, N, 1);         % N x T

%% Find most stable range-bin
range_axis = (0:(N-1)) * range_resolution;
r_min = 0.5; r_max = 1;                     
mag_stat = median(abs(RFFT), 2);                

rid = find(range_axis >= r_min & range_axis <= r_max);
if isempty(rid) || all(~isfinite(mag_stat))
    error('No valid range-bin found, or mag_stat is invalid. Check data and r_min/r_max settings.');
end
[pk_val, rid_local] = max(mag_stat(rid));
peak_idx = rid(rid_local);
DUT_range = range_axis(peak_idx);
fprintf('Peak bin = %d (Range approx. %.2f m)\n', peak_idx, DUT_range);

% Range profile plot
figure; stem(range_axis, mag_stat, 'filled','^'); grid on; hold on;
plot(DUT_range, pk_val, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
text(DUT_range, pk_val*1.05, sprintf('R=%.2f m', DUT_range), ...
    'HorizontalAlignment','center','Color','r');
xlabel('Range (m)'); ylabel('Magnitude (median over frames)');
title('Range Profile (Median Across Frames)'); xlim([0,3]);
ymax = max(mag_stat); if ~isfinite(ymax) || ymax<=0, ymax = 1; end
ylim([0, ymax*1.2]); hold off;

%% Extract slow-time complex sequence at peak range-bin 
y = RFFT(peak_idx, :);                 % 1 x T
num_frames = length(y);
Fs_slow = num_frames / T_meas;         % Slow-time sampling rate based on total duration
fprintf('Frames=%d, T=%.1fs, Fs_slow=%.3f Hz, Nyquist=%.3f Hz\n',...
        num_frames, T_meas, Fs_slow, Fs_slow/2);

%% Phase unwrapping -> detrend -> displacement 
phase_rad = unwrap(angle(y));
fc = 24e9;  lambda = c/fc;
displacement_mm = (lambda/(4*pi)) * phase_rad * 1e3;

t = (0:num_frames-1) / Fs_slow;
figure; plot(t, displacement_mm, 'LineWidth', 1.5); grid on;
xlabel('Time (s)'); ylabel('Displacement (mm)'); xlim([0, t(end)]);
title('Displacement vs Time');

%% Respiration / Heartbeat filtering + segmented peak detection
resp_band  = [0.1, 0.5];
heart_band = [0.6, 1.2];
heart_band(2) = min(heart_band(2), 0.45*Fs_slow);

bp_resp  = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',resp_band(1),'HalfPowerFrequency2',resp_band(2), ...
    'SampleRate',Fs_slow);
bp_heart = designfilt('bandpassiir','FilterOrder',4, ...
    'HalfPowerFrequency1',heart_band(1),'HalfPowerFrequency2',heart_band(2), ...
    'SampleRate',Fs_slow);

win_resp   = round(Fs_slow * 30);
win_heart  = round(Fs_slow * 30);
overlap    = 0.5;

% Respiratory rate trend
step_resp = max(1, round(win_resp*(1-overlap)));
seg_start_resp = 1:step_resp:(num_frames-win_resp);
resp_freqs = nan(1, numel(seg_start_resp));
t_resp     = nan(1, numel(seg_start_resp));
for k = 1:numel(seg_start_resp)
    s = seg_start_resp(k);
    seg = displacement_mm(s:s+win_resp-1);
    seg_f = filtfilt(bp_resp, seg);
    w = hann(length(seg_f)); xw = seg_f(:).*w;
    Nfft_seg = 2^nextpow2(length(xw));
    F = (0:Nfft_seg-1)*(Fs_slow/Nfft_seg);
    S = abs(fft(xw, Nfft_seg));
    rid = find(F>=resp_band(1) & F<=resp_band(2));
    if numel(rid) > 2
        rid = rid(2:end-1);
        [~,ix] = max(S(rid));
        f_peak = F(rid(ix));
        resp_freqs(k) = f_peak;
        t_resp(k) = mean(t(s:s+win_resp-1));
    end
end
resp_bpm_plot = resp_freqs*60;
figure; plot(t_resp, resp_bpm_plot, 'b-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Respiratory Rate (BPM)');
title('Respiratory Rate Trend'); grid on; xlim([0, t(end)]);

% Heart rate trend
step_heart = max(1, round(win_heart*(1-overlap)));
seg_start_heart = 1:step_heart:(num_frames-win_heart);
heart_freqs = nan(1, numel(seg_start_heart));
t_heart     = nan(1, numel(seg_start_heart));
for k = 1:numel(seg_start_heart)
    s = seg_start_heart(k);
    seg = displacement_mm(s:s+win_heart-1);   % [FIX] was s+win_heart (off-by-one)
    seg_f = filtfilt(bp_heart, seg);
    w = hann(length(seg_f)); xw = seg_f(:).*w;
    Nfft_seg = 2^nextpow2(length(xw));
    F = (0:Nfft_seg-1)*(Fs_slow/Nfft_seg);
    S = abs(fft(xw, Nfft_seg));
    hid = find(F>=heart_band(1) & F<=heart_band(2));
    if numel(hid) > 2
        hid = hid(2:end-1);
        [~,jx] = max(S(hid));
        f_peak = F(hid(jx));
        heart_freqs(k) = f_peak;
        t_heart(k) = mean(t(s:s+win_heart-1));
    end
end
heart_bpm_plot = heart_freqs*60;

figure; plot(t_heart, heart_bpm_plot, 'r-o','LineWidth',1.5);
xlabel('Time (s)'); ylabel('Heart Rate (BPM)');
title('Heart Rate Trend'); grid on; xlim([0, t(end)]);

%% First 8 Resp Rate Spectra
figure;
K = min(8, numel(seg_start_resp));
for k = 1:K
    s = seg_start_resp(k);
    seg   = displacement_mm(s:s+win_resp-1);
    seg_f = filtfilt(bp_resp, seg);
    w  = hann(length(seg_f));
    xw = seg_f(:).*w;
    Nfft_seg = 2^nextpow2(length(xw));
    F  = (0:Nfft_seg-1) * (Fs_slow/Nfft_seg);
    S  = abs(fft(xw, Nfft_seg));
    rid = find(F>=resp_band(1) & F<=resp_band(2));
    [~,ix_local] = max(S(rid));
    idx_full = rid(ix_local);
    f_peak = F(idx_full);
    pk     = S(idx_full);
    subplot(4,2,k);
    plot(F, S); grid on; xlim([0, 4]);
    hold on; plot(f_peak, pk, 'ro','MarkerFaceColor','r');
    title(sprintf('Resp %d (Peak=%.2f Hz, %.1f BPM)', k, f_peak, f_peak*60));
    xlabel('Hz'); ylabel('|FFT|');
end
sgtitle('Resp Spectra');

%% First 8 Heart Rate Spectra 
figure;
K = min(8, numel(seg_start_heart));
for k = 1:K
    s = seg_start_heart(k);
    seg   = displacement_mm(s:s+win_heart-1);
    seg_f = filtfilt(bp_heart, seg);
    w  = hann(length(seg_f));
    xw = seg_f(:).*w;
    Nfft_seg = 2^nextpow2(length(xw));
    F  = (0:Nfft_seg-1) * (Fs_slow/Nfft_seg);
    S  = abs(fft(xw, Nfft_seg));
    hid = find(F>=heart_band(1) & F<=heart_band(2));
    [~,ix_local] = max(S(hid));
    idx_full = hid(ix_local);
    f_peak = F(idx_full);
    pk     = S(idx_full);
    subplot(4,2,k);
    plot(F, S); grid on; xlim([0, 2]);
    hold on; plot(f_peak, pk, 'ro','MarkerFaceColor','r');
    title(sprintf('Heart %d (Peak=%.2f Hz, %.1f BPM)', k, f_peak, f_peak*60));
    xlabel('Hz'); ylabel('|FFT|');
end
sgtitle('Heart Rate Spectra');

%% Raw Displacement Segments
seconds_per_plot = 10;
total_duration = t(end);
num_plots = ceil(total_duration / seconds_per_plot);

figure('Name', 'Raw Displacement - Local Inspection (10s per segment)', 'NumberTitle', 'off');
for p = 1:num_plots
    t_start = (p-1) * seconds_per_plot;
    t_end = min(p * seconds_per_plot, total_duration);
    idx = find(t >= t_start & t < t_end);
    if ~isempty(idx)
        subplot(ceil(num_plots/2), 2, p);
        plot(t(idx), displacement_mm(idx), 'LineWidth', 1);
        grid on;
        title(sprintf('Time: %.0f - %.0f s', t_start, t_end));
        xlabel('Time (s)'); ylabel('Displacement (mm)');
        axis tight;
        curr_y = displacement_mm(idx);
        ylim([min(curr_y)-0.1, max(curr_y)+0.1]);
    end
end
sgtitle('Raw Displacement - Local Detail (10s Segments)');

%% Resample to 250 Hz Multirate DSP
% Upsample displacement + filtered signals to 250 Hz so they can be synchronized with ECG/BP reference signal
% I do radar signals resampling only.

Fs_target = 250; 
[p_rs, q_rs] = rat(Fs_target / Fs_slow, 1e-4);
fprintf('Resampling: Fs_slow=%.3f Hz -> %d Hz  (ratio %d/%d)\n', ...
        Fs_slow, Fs_target, p_rs, q_rs);
disp_resp_native   = filtfilt(bp_resp,  displacement_mm(:));
disp_heart_native  = filtfilt(bp_heart, displacement_mm(:));

% Resample to 250 Hz
disp_rs       = resample(displacement_mm(:), p_rs, q_rs);
resp_rs       = resample(disp_resp_native,   p_rs, q_rs);
cardiac_rs    = resample(disp_heart_native,  p_rs, q_rs);

N_rs = length(disp_rs);
t_rs = (0:N_rs-1) / Fs_target;

fprintf('After resampling: %d samples at %d Hz (%.1f s)\n', N_rs, Fs_target, t_rs(end));

% Plot resampled signals
figure('Name','Resampled Signals @ 250 Hz');
subplot(3,1,1);
plot(t_rs, disp_rs, 'k', 'LineWidth', 1);
title('Displacement (resampled @ 250 Hz)');
xlabel('Time (s)'); ylabel('mm'); grid on; xlim([0, t_rs(end)]);

subplot(3,1,2);
plot(t_rs, resp_rs, 'b', 'LineWidth', 1.2);
title(sprintf('Respiratory signal (%.1f–%.1f Hz) @ 250 Hz', resp_band(1), resp_band(2)));
xlabel('Time (s)'); ylabel('mm'); grid on; xlim([0, t_rs(end)]);

subplot(3,1,3);
plot(t_rs, cardiac_rs, 'r', 'LineWidth', 1.2);
title(sprintf('Cardiac signal (%.1f–%.1f Hz) @ 250 Hz', heart_band(1), heart_band(2)));
xlabel('Time (s)'); ylabel('mm'); grid on; xlim([0, t_rs(end)]);
sgtitle('Resampled Signals');


%%  R-Peak Detection 
% In the paper, R-peak detection is performed on the ECG signal. But I apply an algorithm to the radar-derived cardiac signal as an approximation.

fprintf('\n--- R-Peak Detection on radar cardiac signal ---\n');
 
x_card = cardiac_rs(:);   % BPF cardiac signal at 250 Hz
pos_vals = x_card(x_card > 0);
if isempty(pos_vals)
    threshold_pk = 0;
else
    threshold_pk = prctile(pos_vals, 40);
end
 

min_dist_smp = round(0.5 * Fs_target);
prom_thresh = 0.20 * (max(x_card) - min(x_card));
 
[~, r_locs] = findpeaks(x_card, ...
    'MinPeakHeight',      threshold_pk, ...
    'MinPeakDistance',    min_dist_smp, ...
    'MinPeakProminence',  prom_thresh);
 
fprintf('Detected %d R-peaks over %.1f s\n', numel(r_locs), t_rs(end));
 
if numel(r_locs) >= 2
    rr_intervals_s = diff(r_locs) / Fs_target;   % R-R in seconds
    hr_bpm_rr      = 60 ./ rr_intervals_s;        % BPM per beat
    t_rr           = t_rs(r_locs(2:end));
    fprintf('Mean HR from R-R: %.1f BPM  (std: %.1f BPM)\n', ...
            mean(hr_bpm_rr), std(hr_bpm_rr));
else
    warning('Too few R-peaks detected. Try adjusting threshold_pk or prom_thresh.');
    hr_bpm_rr = []; t_rr = [];
end
 
% Plot result
figure('Name','R-Peak Detection');
subplot(2,1,1);
plot(t_rs, x_card, 'r', 'LineWidth', 1); hold on;
if ~isempty(r_locs)
    plot(t_rs(r_locs), x_card(r_locs), 'k^', ...
         'MarkerFaceColor','k', 'MarkerSize', 7);
end
yline(threshold_pk, 'b--', 'LineWidth', 1, 'Label', 'threshold');
title('Cardiac Signal with Detected Peaks');
xlabel('Time (s)'); ylabel('mm'); grid on; xlim([0, t_rs(end)]);
legend('Cardiac signal', 'Detected peaks', 'Location', 'northeast');
 
subplot(2,1,2);
if ~isempty(hr_bpm_rr)
    plot(t_rr, hr_bpm_rr, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 5);
    yline(mean(hr_bpm_rr), 'k--', sprintf('Mean=%.1f BPM', mean(hr_bpm_rr)));
    title('Beat-to-Beat Heart Rate (from peak intervals)');
    xlabel('Time (s)'); ylabel('BPM'); grid on; xlim([0, t_rs(end)]);
    ylim([0, 200]);
end
sgtitle('Peak Detection on Radar Cardiac Signal');


%% Beat Segmentation 

SEG_LEN = 256;  
half    = floor(SEG_LEN / 2);
valid_locs = r_locs(r_locs > half & r_locs <= (N_rs - half));
N_beats    = numel(valid_locs);
fprintf('\nBeat segmentation: %d valid windows (SEG_LEN=%d, Fs=%d Hz)\n', ...
        N_beats, SEG_LEN, Fs_target);

if N_beats > 0
    segments = zeros(N_beats, SEG_LEN, 2);
    for b = 1:N_beats
        idx_win = (valid_locs(b) - half) : (valid_locs(b) + half - 1);
        segments(b, :, 1) = resp_rs(idx_win);
        segments(b, :, 2) = cardiac_rs(idx_win);
    end
    fprintf('Segment tensor shape: [%d x %d x 2]\n', N_beats, SEG_LEN);

    % Plot first 4 beat windows (cardiac channel)
    figure('Name','Beat Segmentation — first 4 windows');
    K_plot = min(4, N_beats);
    t_win  = (0:SEG_LEN-1) / Fs_target * 1000;   % ms
    for b = 1:K_plot
        subplot(K_plot, 2, 2*b-1);
        plot(t_win, segments(b,:,1), 'b', 'LineWidth', 1.2);
        title(sprintf('Beat %d — Respiratory', b));
        xlabel('ms'); ylabel('mm'); grid on;

        subplot(K_plot, 2, 2*b);
        plot(t_win, segments(b,:,2), 'r', 'LineWidth', 1.2);
        title(sprintf('Beat %d — Cardiac', b));
        xlabel('ms'); ylabel('mm'); grid on;
    end
    sgtitle(sprintf('Beat Segments (anchor=R-peak, %d samples @ %d Hz)', SEG_LEN, Fs_target));
else
    warning('No valid beat windows — check R-peak detection results.');
end