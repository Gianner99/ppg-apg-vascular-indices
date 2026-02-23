% =========================================================================
% PPG & APG COMPLETE ANALYSIS
% =========================================================================
%
% Description:
%   Automated pipeline for photoplethysmography (PPG) signal processing
%   and acceleration plethysmogram (APG) fiducial point extraction.
%   Computes APG-derived vascular indices (b/a, c/a, d/a, e/a, AGI) from
%   the PPG-BP database and analyses their correlation with chronological
%   age across clinical groups (Healthy, Hypertensive, Diabetic).
%
% Requirements:
%   - MATLAB R2022b or later
%   - Signal Processing Toolbox (cheby2, filtfilt, sgolay, findpeaks)
%
% Data Availability Statement: 
%   The PPG-BP dataset analysed in this study is publicly available 
%   on Figshare at https://doi.org/10.6084/m9.figshare.5459299 
%   (Liang, Y., Liu, G., Chen, Z. & Elgendi, M. PPG-BP Database). 
%
% Note:
%   This is a preliminary version of the code. For questions, bug reports,
%   or suggestions, please contact: gianluca.diana@unipa.it
%
% License:
%   MIT License
%
%   Copyright (c) 2026 Gianluca Diana
%
%   Permission is hereby granted, free of charge, to any person obtaining
%   a copy of this software and associated documentation files (the
%   "Software"), to deal in the Software without restriction, including
%   without limitation the rights to use, copy, modify, merge, publish,
%   distribute, sublicense, and/or sell copies of the Software, and to
%   permit persons to whom the Software is furnished to do so, subject to
%   the following conditions:
%
%   The above copyright notice and this permission notice shall be
%   included in all copies or substantial portions of the Software.
%
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
%   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
%   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
%   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
%   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
%   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
%   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
% =========================================================================

clearvars; close all; clc;

%% === CONFIGURATION ===
config = struct();
config.Fs = 1000;  % Sampling frequency (Hz)

% Bandpass filter: 4th-order Chebyshev type II, 0.5-8 Hz, 20 dB stopband
config.filter_type = 'cheby2';
config.filter_order = 4;
config.filter_ripple = 20;       % Stopband attenuation (dB)
config.bandpass = [0.5 8];       % Passband (Hz)

% Savitzky-Golay differentiation filter parameters
config.sgolay_order = 3;               % Polynomial order
config.sgolay_framelen_sec = 0.15;     % Window length (seconds)

% Analysis parameters
config.subjects_per_fig = 3;      % Subjects per figure
config.ssqi_threshold = 0.41;     % Minimum SSQI for segment inclusion

% Peak detection parameters
config.peak_min_distance = 0.4;   % Minimum inter-peak distance (s)
config.peak_min_prominence = 0.2; % Minimum peak prominence (normalised)

% Valid peaks temporal range (analysis window within each segment)
config.valid_peak_start = 0.2;    % seconds
config.valid_peak_end = 1.75;     % seconds

% APG normalisation
config.normalize_apg = true;

% --- Fiducial point detection thresholds ---
% These are orientative values tuned for reliable identification of APG
% fiducial points on the PPG-BP dataset. Minor adjustments may be needed
% for signals acquired with different hardware or sampling rates.
config.notch_min_time = 0.1;      % Min time after systolic peak for dicrotic notch (s)
config.notch_max_time = 0.3;      % Max time after systolic peak for dicrotic notch (s)
config.apg_c_max_time = 0.15;     % Max time after 'b' for point 'c' (s)
config.apg_d_max_time = 0.22;     % Max time after 'b' for point 'd' (s)
config.apg_e_max_time = 0.30;     % Max time after 'b' for point 'e' (s)
config.min_time_after_b = 0.05;   % Min time between b and c/d (s)
config.min_time_before_e = 0.05;  % Min time between c/d and e (s)
config.fallback_a_window = 0.3;   % Window before b to search for 'a' in fallback (s)
config.fallback_window = 0.4;     % Window around systolic peak for third fallback (s)
config.apg_peak_min_dist = 0.03;  % Min distance between APG peaks (s)
config.apg_peak_min_prom = 0.03;  % Min prominence for APG peaks (normalised)
config.val_a_threshold = 0.25;    % Minimum amplitude for valid 'a' wave (normalised)
config.val_b_threshold = -0.1;    % Maximum amplitude for valid 'b' wave (normalised)
config.val_b_fallback_thr = -0.05;% Maximum amplitude for 'b' in third fallback (normalised)

% Outlier removal for correlation plots (mean +/- threshold * std)
config.outlier_std_threshold = 2.0;

% File paths
config.excel_file = fullfile(pwd, 'PPG-BP dataset.xlsx');
config.ppg_folder = fullfile(pwd, '0_subject');

% Create output folder with timestamp
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
config.output_folder = fullfile(pwd, sprintf('PPG_Analysis_%s', timestamp));
config.screen_folder = fullfile(config.output_folder, 'screen');
config.matlab_folder = fullfile(config.output_folder, 'matlab_figures');
config.correlation_folder = fullfile(config.output_folder, 'correlation_plots');
mkdir(config.output_folder);
mkdir(config.screen_folder);
mkdir(config.matlab_folder);
mkdir(config.correlation_folder);

fprintf('=== PPG & APG ANALYSIS ===\n');
fprintf('Output folder: %s\n\n', config.output_folder);

%% === LOAD EXCEL DATA ===
fprintf('Loading Excel file...\n');
T = readtable(config.excel_file, "VariableNamingRule", "modify");

% Find required columns
id_var = findColumn(T, ["subject_id", "id"], ["subject", "id"]);
[ssqi_vars, ok] = findSSQIColumns(T);
if ~ok
    error('Cannot find SSQI columns in Excel file');
end

% Find optional metadata columns
meta = struct();
meta.sex_var = findColumnOptional(T, "sex", "sex");
meta.age_var = findColumnOptional(T, "age", "age");
meta.hyper_var = findColumnOptional(T, "hypertension", "hypertension");
meta.diab_var = findColumnOptional(T, "diabetes", "diabetes");
meta.cinf_var = findColumnOptional(T, "cerebral_infarction", ["cerebral", "infar"]);

fprintf('Configuration loaded successfully\n\n');

%% === SELECT ALL VALID SUBJECTS ===
fprintf('Selecting subjects with SSQI >= %.2f...\n', config.ssqi_threshold);

n_total = height(T);
valid_subjects_temp = zeros(1, n_total);
valid_ages = zeros(1, n_total);
valid_ids = zeros(1, n_total);
count = 0;

for i = 1:n_total
    % Check if at least one segment has valid SSQI
    has_valid_segment = false;
    for seg = 1:3
        ssqi_val = convertToNum(T.(ssqi_vars(seg))(i));
        if ssqi_val >= config.ssqi_threshold
            has_valid_segment = true;
            break;
        end
    end

    if has_valid_segment && meta.age_var ~= ""
        age = convertToNum(T.(meta.age_var)(i));
        if ~isnan(age)
            count = count + 1;
            valid_subjects_temp(count) = i;
            valid_ages(count) = age;
            valid_ids(count) = convertToNum(T.(id_var)(i));
        end
    end
end

% Trim preallocated arrays
valid_subjects_temp = valid_subjects_temp(1:count);
valid_ages = valid_ages(1:count);
valid_ids = valid_ids(1:count);

fprintf('Found %d valid subjects with age data\n', count);

% Sort by ID for figures
[sorted_ids, sort_idx] = sort(valid_ids);
sorted_subjects = valid_subjects_temp(sort_idx);
sorted_ages = valid_ages(sort_idx);

fprintf('ID range: %d - %d\n', min(sorted_ids), max(sorted_ids));
fprintf('Age range: %d - %d years\n\n', min(sorted_ages), max(sorted_ages));

%% === INITIALIZE RESULTS STRUCTURE ===
results = struct();
for i = 1:n_total
    results(i).subject_id = convertToNum(T.(id_var)(i));
    results(i).age = NaN;
    if meta.age_var ~= ""
        results(i).age = convertToNum(T.(meta.age_var)(i));
    end

    for seg = 1:3
        seg_name = sprintf('seg%d', seg);
        results(i).(seg_name).processed = false;
        results(i).(seg_name).onset_times = [];
        results(i).(seg_name).peak_times = [];
        results(i).(seg_name).notch_times = [];
        results(i).(seg_name).t_a = [];
        results(i).(seg_name).a = [];
        results(i).(seg_name).t_b = [];
        results(i).(seg_name).b = [];
        results(i).(seg_name).t_c = [];
        results(i).(seg_name).c = [];
        results(i).(seg_name).t_d = [];
        results(i).(seg_name).d = [];
        results(i).(seg_name).t_e = [];
        results(i).(seg_name).e = [];
    end
end

%% === PROCESS ALL SUBJECTS AND GENERATE FIGURES ===
n_figures = ceil(length(sorted_subjects) / config.subjects_per_fig);

fprintf('Processing %d subjects across %d figures (ordered by ID)...\n\n', ...
        length(sorted_subjects), n_figures);

for fig_idx = 1:n_figures
    % Get subjects for this figure
    start_idx = (fig_idx - 1) * config.subjects_per_fig + 1;
    end_idx = min(fig_idx * config.subjects_per_fig, length(sorted_subjects));
    subjects_in_fig = sorted_subjects(start_idx:end_idx);

    fig = figure('Color', 'w', 'Position', [50 50 1800 1400], ...
                 'Name', sprintf('PPG & APG Analysis - Fig %d', fig_idx));

    subject_ids_str = '';

    for s_idx = 1:length(subjects_in_fig)
        row_idx = subjects_in_fig(s_idx);
        subj_id = convertToNum(T.(id_var)(row_idx));

        if s_idx == 1
            subject_ids_str = sprintf('%d', subj_id);
        else
            subject_ids_str = sprintf('%s_%d', subject_ids_str, subj_id);
        end

        fprintf('Processing Subject ID %d...\n', subj_id);
        metadata = extractMetadata(T, row_idx, meta);
        base_row = (s_idx - 1) * 2;

        for seg = 1:3
            fname = fullfile(config.ppg_folder, sprintf('%d_%d.txt', subj_id, seg));
            raw_ppg = loadPPGFile(fname);

            if isempty(raw_ppg)
                warning('File missing: %s', fname);
                continue;
            end

            ssqi_val = convertToNum(T.(ssqi_vars(seg))(row_idx));

            % --- Signal Processing Pipeline ---
            % Step 1: Bandpass filter (Chebyshev type II)
            filtered_ppg = applyBandpassFilter(raw_ppg, config);

            % Step 2: Normalise PPG to [0, 1]
            filtered_ppg = (filtered_ppg - min(filtered_ppg)) / ...
                          (max(filtered_ppg) - min(filtered_ppg));

            % Step 3: Find systolic peaks
            [peaks_val, peaks_idx, valid_peaks_idx, valid_peaks_time, valid_peaks_val] = ...
                findSystolicPeaks(filtered_ppg, config);

            t = (0:length(filtered_ppg)-1) / config.Fs;
            peaks_time = t(peaks_idx);

            % Step 4: Find PPG onset (foot points)
            [onset_idx, onset_val] = findPPGOnsets(filtered_ppg, peaks_idx, config.Fs);

            % Step 5: Find dicrotic notch
            [notch_idx, notch_val] = findDicroticNotch(filtered_ppg, peaks_idx, config);

            % Step 6: Calculate APG via Savitzky-Golay differentiation
            [apg, t_apg] = calculateAPG(filtered_ppg, config);

            % Step 7: Extract APG fiducial points (a, b, c, d, e)
            apg_points = extractAPGPoints(t, apg, valid_peaks_time, config);

            % --- Store results ---
            seg_name = sprintf('seg%d', seg);
            results(row_idx).(seg_name).processed = true;

            n_valid_peaks = length(valid_peaks_idx);
            selected_peak_indices = selectBestPeaks(n_valid_peaks, apg_points);

            valid_mask = ismember(peaks_idx, valid_peaks_idx);

            if ~isempty(onset_idx)
                onset_times_all = t(onset_idx(valid_mask));
                results(row_idx).(seg_name).onset_times = onset_times_all(selected_peak_indices);
            end

            results(row_idx).(seg_name).peak_times = valid_peaks_time(selected_peak_indices);

            if ~isempty(notch_idx)
                notch_times_valid = matchNotchesToPeaks(notch_idx, valid_peaks_idx, t, config.Fs);
                results(row_idx).(seg_name).notch_times = notch_times_valid(selected_peak_indices);
            end

            if ~isempty(apg_points)
                for lbl = {'a','b','c','d','e'}
                    t_field = ['t_' lbl{1}];
                    results(row_idx).(seg_name).(t_field) = selectPeakData(apg_points.(t_field), selected_peak_indices);
                    results(row_idx).(seg_name).(lbl{1}) = selectPeakData(apg_points.(lbl{1}), selected_peak_indices);
                end
            end

            % --- Plot ---
            plotPPG(base_row, seg, t, filtered_ppg, peaks_idx, peaks_val, ...
                   valid_peaks_idx, valid_peaks_val, onset_idx, onset_val, ...
                   notch_idx, notch_val, subj_id, ssqi_val, metadata);

            plotAPG(base_row, seg, t_apg, apg, apg_points, subj_id, ...
                   ssqi_val, metadata, config.normalize_apg);
        end
    end

    % Figure title and save
    ids_in_fig = sorted_ids(start_idx:end_idx);
    ages_in_fig = sorted_ages(start_idx:end_idx);
    sgtitle(sprintf('PPG & APG Analysis - Subjects: %s  |  Age: %d-%d years', ...
                    strrep(subject_ids_str, '_', ', '), ...
                    min(ages_in_fig), max(ages_in_fig)), ...
            'FontSize', 18, 'FontWeight', 'bold');

    fig_name = sprintf('Fig_%03d_IDs_%s', fig_idx, subject_ids_str);
    saveas(fig, fullfile(config.screen_folder, [fig_name '.png']));
    savefig(fig, fullfile(config.matlab_folder, [fig_name '.fig']));
    fprintf('Saved figure: %s (PNG + FIG)\n\n', fig_name);
    close(fig);
end

%% === EXPORT RESULTS TO EXCEL ===
fprintf('Exporting results to Excel...\n');
exportResultsToExcel(T, results, config, meta, id_var, ssqi_vars);

%% === CREATE CORRELATION PLOTS ===
fprintf('\nGenerating correlation plots (Age vs Vascular Parameters)...\n');
createCorrelationPlots(T, results, config, meta, id_var);

%% === SAVE MATLAB CODE COPY ===
fprintf('\nSaving MATLAB code copy...\n');
try
    current_file = mfilename('fullpath');
    if ~isempty(current_file)
        if ~endsWith(current_file, '.m')
            current_file = [current_file '.m'];
        end
        if exist(current_file, 'file')
            [~, name, ext] = fileparts(current_file);
            copyfile(current_file, fullfile(config.output_folder, [name ext]));
            fprintf('MATLAB code saved successfully\n');
        end
    end
catch ME
    fprintf('Warning: Could not copy MATLAB code (%s)\n', ME.message);
end

fprintf('\n=== ANALYSIS COMPLETED ===\n');
fprintf('Total subjects processed: %d\n', length(sorted_subjects));
fprintf('Total figures saved: %d (ordered by ID)\n', n_figures);
fprintf('Correlation plots: 10 parameters\n');
fprintf('Output folder: %s\n', config.output_folder);

%% =====================================================================
%%  HELPER FUNCTIONS
%% =====================================================================

% --- Signal Processing ---

function filtered = applyBandpassFilter(signal, config)
    % Apply 4th-order Chebyshev type II bandpass filter (0.5-8 Hz).
    Fs = config.Fs;
    Wp = config.bandpass / (Fs/2);
    [z, p, k] = cheby2(config.filter_order, config.filter_ripple, Wp, 'bandpass');
    [sos, g] = zp2sos(z, p, k);
    filtered = filtfilt(sos, g, signal);
end

function [apg, t_apg] = calculateAPG(ppg, config)
    % Calculate APG using Savitzky-Golay differentiation.
    %
    % The VPG (first derivative) and APG (second derivative) are computed
    % by applying the Savitzky-Golay filter with first-order derivation
    % twice in succession, as described in the paper:
    %   - VPG = SG 1st derivative of PPG
    %   - APG = SG 1st derivative of VPG
    %
    % The SG filter performs simultaneous smoothing and differentiation,
    % intrinsically attenuating high-frequency noise amplification without
    % requiring additional filtering.

    frameLen = max(11, round(config.sgolay_framelen_sec * config.Fs));
    if mod(frameLen, 2) == 0
        frameLen = frameLen + 1;
    end
    if frameLen > length(ppg)
        frameLen = length(ppg);
        if mod(frameLen, 2) == 0
            frameLen = frameLen - 1;
        end
    end

    dt = 1 / config.Fs;

    % Compute SG differentiation filter coefficients
    [~, g] = sgolay(config.sgolay_order, frameLen);
    sg_diff_coeff = flipud(g(:, 2));  % 1st derivative coefficients

    % VPG: 1st derivative of PPG via SG differentiation
    vpg = conv(ppg(:), sg_diff_coeff, 'same') / dt;

    % APG: 1st derivative of VPG via SG differentiation (= 2nd derivative of PPG)
    apg = conv(vpg, sg_diff_coeff, 'same') / dt;

    % Normalise APG to [-1, 1]
    if config.normalize_apg
        max_abs_val = max(abs(apg));
        if max_abs_val > 0
            apg = apg / max_abs_val;
        end
    end

    apg = apg(:)';
    t_apg = (0:length(apg)-1) / config.Fs;
end

% --- Peak Detection ---

function [peaks_val, peaks_idx, valid_peaks_idx, valid_peaks_time, valid_peaks_val] = findSystolicPeaks(ppg, config)
    % Find systolic peaks and filter those within the valid time range.
    min_peak_distance = round(config.peak_min_distance * config.Fs);
    [peaks_val, peaks_idx] = findpeaks(ppg, ...
        'MinPeakDistance', min_peak_distance, ...
        'MinPeakProminence', config.peak_min_prominence);

    t = (0:length(ppg)-1) / config.Fs;
    peaks_time = t(peaks_idx);

    valid_mask = (peaks_time >= config.valid_peak_start) & ...
                 (peaks_time <= config.valid_peak_end);
    valid_peaks_idx = peaks_idx(valid_mask);
    valid_peaks_time = peaks_time(valid_mask);
    valid_peaks_val = peaks_val(valid_mask);
end

function [onset_idx, onset_val] = findPPGOnsets(ppg, peaks_idx, Fs)
    % Find PPG onset (foot point) before each systolic peak.
    n_peaks = length(peaks_idx);
    onset_idx = NaN(1, n_peaks);
    onset_val = NaN(1, n_peaks);

    for p = 1:n_peaks
        peak_idx = peaks_idx(p);

        if p == 1
            search_start = max(1, peak_idx - round(0.8*Fs));
        else
            search_start = peaks_idx(p-1) + round(0.1*Fs);
        end
        search_end = peak_idx - round(0.05*Fs);

        if search_end > search_start
            [min_val, min_rel_idx] = min(ppg(search_start:search_end));
            onset_idx(p) = search_start + min_rel_idx - 1;
            onset_val(p) = min_val;
        end
    end
end

function [notch_idx, notch_val] = findDicroticNotch(ppg, peaks_idx, config)
    % Find dicrotic notch after each systolic peak using VPG-based method.
    n_peaks = length(peaks_idx);
    notch_idx = NaN(1, n_peaks);
    notch_val = NaN(1, n_peaks);
    vpg = gradient(ppg) * config.Fs;

    for p = 1:n_peaks
        peak_idx = peaks_idx(p);

        if p < n_peaks
            next_peak_idx = peaks_idx(p+1);
        else
            next_peak_idx = min(peak_idx + round(0.6*config.Fs), length(ppg));
        end

        search_valley_start = peak_idx + round(0.1*config.Fs);
        search_valley_end = next_peak_idx - round(0.05*config.Fs);

        if search_valley_end <= search_valley_start
            continue;
        end

        [~, valley_rel_idx] = min(ppg(search_valley_start:search_valley_end));
        valley_idx = search_valley_start + valley_rel_idx - 1;

        % Notch search window: proportional to peak-valley distance,
        % constrained by absolute temporal limits
        dist_to_valley = valley_idx - peak_idx;
        notch_start = peak_idx + round(0.35 * dist_to_valley);
        notch_end = peak_idx + round(0.50 * dist_to_valley);

        abs_min_idx = peak_idx + round(config.notch_min_time * config.Fs);
        abs_max_idx = peak_idx + round(config.notch_max_time * config.Fs);

        notch_start = max(notch_start, abs_min_idx);
        notch_end = min(notch_end, abs_max_idx);
        notch_start = max(notch_start, peak_idx + 20);
        notch_end = min(notch_end, valley_idx - 10);

        % Require minimum 30 samples in search window
        if notch_end <= notch_start + 30
            continue;
        end

        vpg_seg = vpg(notch_start:notch_end);
        prom_thresh = 0.025 * abs(min(vpg_seg));

        [~, locs_vpg] = findpeaks(-vpg_seg, ...
            'MinPeakProminence', prom_thresh, ...
            'SortStr', 'descend', 'NPeaks', 3);

        if isempty(locs_vpg)
            continue;
        end

        if length(locs_vpg) == 1
            notch_rel_idx = locs_vpg(1);
        else
            notch_rel_idx = min(locs_vpg(1:min(2, length(locs_vpg))));
        end

        notch_idx_cand = notch_start + notch_rel_idx - 1;
        time_after_peak = (notch_idx_cand - peak_idx) / config.Fs;

        if time_after_peak >= config.notch_min_time && ...
           time_after_peak <= config.notch_max_time
            % Refine notch position by finding local minimum
            refine_start = max(1, notch_idx_cand - 10);
            refine_end = min(length(ppg), notch_idx_cand + 10);
            [refined_val, refined_rel] = min(ppg(refine_start:refine_end));
            notch_idx(p) = refine_start + refined_rel - 1;
            notch_val(p) = refined_val;
        end
    end
end

% --- APG Fiducial Point Extraction ---

function apg_points = extractAPGPoints(t, apg, valid_peaks_time, config)
    % Extract APG fiducial points (a, b, c, d, e) with temporal constraints.
    %
    % Three methods are applied in order:
    %   1. Primary: detect 'a' as first prominent maximum near the systolic
    %      peak, then 'b' as first prominent minimum after 'a'.
    %   2. Fallback 1: if 'b' is found but not 'a', search for 'a' as the
    %      maximum in a window before 'b'.
    %   3. Fallback 2: if neither 'a' nor 'b' are found, use the absolute
    %      minimum in a window around the systolic peak as 'b' and derive
    %      the remaining points.
    %
    % Points c, d, e are extracted via a shared subroutine (findCDE).

    n_peaks = length(valid_peaks_time);

    % Preallocate
    apg_points = struct();
    for lbl = {'a','b','c','d','e'}
        apg_points.(['t_' lbl{1}]) = NaN(1, n_peaks);
        apg_points.(lbl{1}) = NaN(1, n_peaks);
    end

    if n_peaks == 0
        return;
    end

    for i = 1:n_peaks
        tPeak = valid_peaks_time(i);

        % Define analysis window around the systolic peak
        tStart = max(tPeak - 0.2, t(1));
        tEnd = min(tPeak + 0.5, t(end));
        if i < n_peaks
            tEnd = min(tEnd, valid_peaks_time(i+1) - 0.1);
        end

        idx = (t >= tStart) & (t <= tEnd);
        if sum(idx) < 20
            continue;
        end

        tSeg = t(idx);
        apgSeg = apg(idx);

        [maxVals, maxLocs] = findpeaks(apgSeg, tSeg, ...
            'MinPeakDistance', config.apg_peak_min_dist, ...
            'MinPeakProminence', config.apg_peak_min_prom);
        [minVals_inv, minLocs] = findpeaks(-apgSeg, tSeg, ...
            'MinPeakDistance', config.apg_peak_min_dist, ...
            'MinPeakProminence', config.apg_peak_min_prom);
        minVals = -minVals_inv;

        if isempty(maxLocs) || isempty(minLocs)
            continue;
        end

        % --- Primary method: find 'a' then 'b' ---
        idx_a = find(maxLocs >= tStart & maxLocs <= tPeak+0.1, 1, 'first');

        if ~isempty(idx_a) && maxVals(idx_a) > config.val_a_threshold
            t_a = maxLocs(idx_a);
            val_a = maxVals(idx_a);

            idx_b = find(minLocs > t_a, 1, 'first');
            if ~isempty(idx_b) && minVals(idx_b) < config.val_b_threshold
                t_b = minLocs(idx_b);
                val_b = minVals(idx_b);

                % Find c, d, e
                [t_c, val_c, t_d, val_d, t_e, val_e] = ...
                    findCDE(maxLocs, maxVals, minLocs, minVals, t_b, config);

                % Second fallback: if a, b, e found but not c, d
                if ~isnan(t_e) && isnan(t_c) && isnan(t_d)
                    [t_c, val_c, t_d, val_d] = ...
                        findCD_between_b_and_e(maxLocs, maxVals, minLocs, minVals, t_b, t_e, config);
                end

                storePoint(i);
                continue;
            end
        end

        % --- Fallback 1: 'b' found but not 'a' ---
        idx_b = find(minLocs > tStart, 1, 'first');
        if ~isempty(idx_b) && minVals(idx_b) < config.val_b_threshold
            t_b = minLocs(idx_b);
            val_b = minVals(idx_b);

            % Find 'a' as maximum in window before 'b'
            t_search_start = max(t_b - config.fallback_a_window, t(1));
            idx_search = (t >= t_search_start) & (t <= t_b);
            if sum(idx_search) > 0
                [val_a, idx_max] = max(apg(idx_search));
                t_search = t(idx_search);
                t_a = t_search(idx_max);

                [t_c, val_c, t_d, val_d, t_e, val_e] = ...
                    findCDE(maxLocs, maxVals, minLocs, minVals, t_b, config);

                storePoint(i);
                continue;
            end
        end

        % --- Fallback 2: use absolute minimum as 'b' ---
        tryFallbackAbsMin(i);
    end

    % === Nested functions (share workspace with extractAPGPoints) ===

    function storePoint(pk_idx)
        apg_points.t_a(pk_idx) = t_a;
        apg_points.a(pk_idx) = val_a;
        apg_points.t_b(pk_idx) = t_b;
        apg_points.b(pk_idx) = val_b;
        apg_points.t_c(pk_idx) = t_c;
        apg_points.c(pk_idx) = val_c;
        apg_points.t_d(pk_idx) = t_d;
        apg_points.d(pk_idx) = val_d;
        apg_points.t_e(pk_idx) = t_e;
        apg_points.e(pk_idx) = val_e;
    end

    function tryFallbackAbsMin(pk_idx)
        tPk = valid_peaks_time(pk_idx);
        tStart_fb = max(tPk - config.fallback_window/2, t(1));
        tEnd_fb = min(tPk + config.fallback_window/2, t(end));

        idx_fb = (t >= tStart_fb) & (t <= tEnd_fb);
        if sum(idx_fb) < 20
            return;
        end

        tSeg_fb = t(idx_fb);
        apgSeg_fb = apg(idx_fb);

        [val_b_local, idx_b_local] = min(apgSeg_fb);
        if val_b_local >= config.val_b_fallback_thr
            return;
        end

        t_b = tSeg_fb(idx_b_local);
        val_b = val_b_local;

        % Find 'a' as maximum before 'b'
        t_a = NaN; val_a = NaN;
        t_search_start_fb = max(t_b - config.fallback_a_window, t(1));
        idx_a_search = (t >= t_search_start_fb) & (t <= t_b);
        if sum(idx_a_search) > 0
            apg_a_search = apg(idx_a_search);
            t_a_search = t(idx_a_search);
            [val_a, idx_max_a] = max(apg_a_search);
            t_a = t_a_search(idx_max_a);
        end

        % Find peaks/valleys in fallback window
        [mxV, mxL] = findpeaks(apgSeg_fb, tSeg_fb, 'MinPeakDistance', 0.02);
        [mnV_inv, mnL] = findpeaks(-apgSeg_fb, tSeg_fb, 'MinPeakDistance', 0.02);
        mnV = -mnV_inv;

        % Find c
        t_c = NaN; val_c = NaN;
        idx_c_fb = find(mxL > t_b, 1, 'first');
        if ~isempty(idx_c_fb) && (mxL(idx_c_fb) - t_b) >= config.min_time_after_b
            t_c = mxL(idx_c_fb);
            val_c = mxV(idx_c_fb);
        end

        % Find d
        t_d = NaN; val_d = NaN;
        if ~isnan(t_c)
            idx_d_fb = find(mnL > t_c, 1, 'first');
        else
            idx_d_fb = find(mnL > t_b, 1, 'first');
        end
        if ~isempty(idx_d_fb) && (mnL(idx_d_fb) - t_b) >= config.min_time_after_b
            t_d = mnL(idx_d_fb);
            val_d = mnV(idx_d_fb);
        end

        % Find e (maximum among candidates after d or c or b)
        t_e = NaN; val_e = NaN;
        if ~isnan(t_d)
            idx_e_all = find(mxL > t_d);
        elseif ~isnan(t_c)
            idx_e_all = find(mxL > t_c);
        else
            idx_e_all = find(mxL > t_b);
        end

        if ~isempty(idx_e_all)
            [val_e_cand, rel_idx] = max(mxV(idx_e_all));
            t_e_cand = mxL(idx_e_all(rel_idx));

            valid_e = true;
            if ~isnan(t_c) && (t_e_cand - t_c) < config.min_time_before_e
                valid_e = false;
            end
            if ~isnan(t_d) && (t_e_cand - t_d) < config.min_time_before_e
                valid_e = false;
            end
            if (t_e_cand - t_b) < config.min_time_after_b
                valid_e = false;
            end
            if valid_e
                t_e = t_e_cand;
                val_e = val_e_cand;
            end
        end

        % Discard c or d if they coincide with e
        if ~isnan(t_c) && ~isnan(t_e) && abs(t_c - t_e) < 0.01
            t_c = NaN; val_c = NaN;
        end
        if ~isnan(t_d) && ~isnan(t_e) && abs(t_d - t_e) < 0.01
            t_d = NaN; val_d = NaN;
        end

        storePoint(pk_idx);
    end
end

function [t_c, val_c, t_d, val_d, t_e, val_e] = findCDE(maxLocs, maxVals, minLocs, minVals, t_b, config)
    % Find APG points c, d, e after point b using temporal constraints.
    t_c = NaN; val_c = NaN;
    t_d = NaN; val_d = NaN;
    t_e = NaN; val_e = NaN;

    % Point 'c': first local maximum after b, within c_max_time
    t_max_c = t_b + config.apg_c_max_time;
    idx_c = find(maxLocs > t_b & maxLocs <= t_max_c, 1, 'first');
    if ~isempty(idx_c)
        t_c = maxLocs(idx_c);
        val_c = maxVals(idx_c);
    end

    % Point 'd': first local minimum after c (or b), within d_max_time
    t_max_d = t_b + config.apg_d_max_time;
    if ~isnan(t_c)
        idx_d = find(minLocs > t_c & minLocs <= t_max_d, 1, 'first');
    else
        idx_d = find(minLocs > t_b & minLocs <= t_max_d, 1, 'first');
    end
    if ~isempty(idx_d)
        t_d = minLocs(idx_d);
        val_d = minVals(idx_d);
    end

    % Point 'e': first local maximum after d (or c, or b), within e_max_time
    t_max_e = t_b + config.apg_e_max_time;
    if ~isnan(t_d)
        idx_e = find(maxLocs > t_d & maxLocs <= t_max_e, 1, 'first');
    elseif ~isnan(t_c)
        idx_e = find(maxLocs > t_c & maxLocs <= t_max_e, 1, 'first');
    else
        idx_e = find(maxLocs > t_b & maxLocs <= t_max_e, 1, 'first');
    end
    if ~isempty(idx_e)
        t_e = maxLocs(idx_e);
        val_e = maxVals(idx_e);
    end
end

function [t_c, val_c, t_d, val_d] = findCD_between_b_and_e(maxLocs, maxVals, minLocs, minVals, t_b, t_e, config)
    % Second fallback: find c and d between b and e when primary method
    % found a, b, e but missed c, d.
    t_c = NaN; val_c = NaN;
    t_d = NaN; val_d = NaN;

    idx_c = find(maxLocs > t_b & maxLocs < t_e, 1, 'first');
    if ~isempty(idx_c)
        t_c_cand = maxLocs(idx_c);
        if (t_c_cand - t_b) >= config.min_time_after_b && ...
           (t_e - t_c_cand) >= config.min_time_before_e
            t_c = t_c_cand;
            val_c = maxVals(idx_c);

            idx_d = find(minLocs > t_c & minLocs < t_e, 1, 'first');
            if ~isempty(idx_d)
                t_d_cand = minLocs(idx_d);
                if (t_d_cand - t_b) >= config.min_time_after_b && ...
                   (t_e - t_d_cand) >= config.min_time_before_e
                    t_d = t_d_cand;
                    val_d = minVals(idx_d);
                end
            end
        end
    end
end

% --- Peak Selection ---

function selected_indices = selectBestPeaks(n_valid_peaks, apg_points)
    % Select up to 2 peaks with the most complete APG fiducial points.
    if n_valid_peaks <= 2
        selected_indices = 1:n_valid_peaks;
        return;
    end

    if ~isempty(apg_points)
        completeness = zeros(1, n_valid_peaks);
        for pk = 1:n_valid_peaks
            for lbl = {'a','b','c','d','e'}
                t_field = ['t_' lbl{1}];
                if pk <= length(apg_points.(t_field)) && ~isnan(apg_points.(t_field)(pk))
                    completeness(pk) = completeness(pk) + 1;
                end
            end
        end

        peaks_with_5 = find(completeness == 5);
        if length(peaks_with_5) >= 2
            selected_indices = peaks_with_5(1:2);
        elseif length(peaks_with_5) == 1
            other_peaks = setdiff(1:n_valid_peaks, peaks_with_5);
            selected_indices = [peaks_with_5(1), other_peaks(1)];
        else
            selected_indices = [1, 2];
        end
    else
        selected_indices = [1, 2];
    end
end

function notch_times = matchNotchesToPeaks(notch_idx, valid_peaks_idx, t, Fs)
    % Match each dicrotic notch to its corresponding valid peak.
    n_vp = length(valid_peaks_idx);
    notch_times = NaN(1, n_vp);

    for vp = 1:n_vp
        vp_idx = valid_peaks_idx(vp);
        notch_match = notch_idx(notch_idx > vp_idx & ~isnan(notch_idx));
        if ~isempty(notch_match) && notch_match(1) < vp_idx + 0.5*Fs
            notch_times(vp) = t(notch_match(1));
        end
    end
end

function data = selectPeakData(source_data, selected_indices)
    % Extract data for selected peak indices.
    n = length(selected_indices);
    data = NaN(1, n);
    for k = 1:n
        idx = selected_indices(k);
        if idx <= length(source_data)
            data(k) = source_data(idx);
        end
    end
end

% --- Plotting ---

function plotPPG(base_row, seg, t, ppg, peaks_idx, peaks_val, valid_idx, valid_val, ...
                 onset_idx, onset_val, notch_idx, notch_val, id, ssqi, metadata)
    subplot(6, 3, base_row*3 + seg);
    plot(t, ppg, 'b', 'LineWidth', 1.2);
    hold on;

    if ~isempty(peaks_idx)
        plot(t(peaks_idx), peaks_val, 'ro', 'MarkerSize', 7, ...
             'MarkerFaceColor', 'r', 'DisplayName', 'All peaks');
    end
    if ~isempty(valid_idx)
        plot(t(valid_idx), valid_val, 'go', 'MarkerSize', 9, ...
             'LineWidth', 2, 'DisplayName', 'Valid peaks');
    end
    if ~isempty(onset_idx)
        valid_onsets = ~isnan(onset_idx);
        if any(valid_onsets)
            plot(t(onset_idx(valid_onsets)), onset_val(valid_onsets), 'ks', ...
                 'MarkerSize', 6, 'MarkerFaceColor', 'k', 'DisplayName', 'Onset');
        end
    end
    if ~isempty(notch_idx)
        valid_notches = ~isnan(notch_idx);
        if any(valid_notches)
            plot(t(notch_idx(valid_notches)), notch_val(valid_notches), 'md', ...
                 'MarkerSize', 6, 'MarkerFaceColor', 'm', 'DisplayName', 'Notch');
        end
    end
    hold off;

    grid on; box on;
    xlim([0 max(t)]);
    ylim([-0.1 1.1]);
    xlabel('Time (s)', 'FontSize', 9);
    ylabel('PPG (norm)', 'FontSize', 9);
    title(sprintf('ID %d | Seg %d | PPG | SSQI=%.2f | %s %dy', ...
                  id, seg, ssqi, metadata.sex, metadata.age), ...
          'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
    set(gca, 'FontSize', 9);

    if base_row == 0 && seg == 1
        leg = legend('PPG', 'All peaks', 'Valid', 'Onset', 'Notch', ...
              'FontSize', 7, 'Location', 'northoutside', 'Orientation', 'horizontal');
        leg.Position(2) = leg.Position(2) + 0.02;
    end

    % Background colour based on SSQI quality
    if ssqi >= 0.6
        set(gca, 'Color', [0.93 1.00 0.93]);  % Green (excellent)
    elseif ssqi >= 0.41
        set(gca, 'Color', [1.00 1.00 0.92]);  % Yellow (acceptable)
    else
        set(gca, 'Color', [1.00 0.95 0.95]);  % Red (poor)
    end
end

function plotAPG(base_row, seg, t_apg, apg, apg_points, id, ssqi, metadata, normalize)
    subplot(6, 3, (base_row+1)*3 + seg);
    plot(t_apg, apg, 'Color', [0.8 0 0], 'LineWidth', 1.2);
    hold on;
    plot([0 max(t_apg)], [0 0], 'k--', 'LineWidth', 0.5);

    if ~isempty(apg_points)
        colors_map = struct('a', [1 0 0], 'b', [0 0 1], 'c', [0 0.7 0], ...
                       'd', [1 0.5 0], 'e', [0.5 0 0.5]);
        markers_map = struct('a', '^', 'b', 'v', 'c', 'o', 'd', 's', 'e', 'd');
        labels = {'a', 'b', 'c', 'd', 'e'};

        for lbl_idx = 1:length(labels)
            lbl = labels{lbl_idx};
            times = apg_points.(['t_' lbl]);
            vals = apg_points.(lbl);
            valid_mask = ~isnan(times) & ~isnan(vals);
            if any(valid_mask)
                scatter(times(valid_mask), vals(valid_mask), 70, ...
                       colors_map.(lbl), markers_map.(lbl), 'filled', ...
                       'LineWidth', 1.5, 'MarkerEdgeColor', 'k', ...
                       'DisplayName', lbl);
            end
        end
    end
    hold off;

    grid on; box on;
    xlim([0 max(t_apg)]);
    xlabel('Time (s)', 'FontSize', 9);
    if normalize
        ylabel('APG (norm)', 'FontSize', 9);
    else
        ylabel('APG (a.u.)', 'FontSize', 9);
    end
    title(sprintf('ID %d | Seg %d | APG | SSQI=%.2f | %s %dy', ...
                  id, seg, ssqi, metadata.sex, metadata.age), ...
          'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'none');
    set(gca, 'FontSize', 9);

    if base_row == 0 && seg == 1 && ~isempty(apg_points)
        leg = legend('APG', 'Zero', 'a', 'b', 'c', 'd', 'e', ...
              'FontSize', 7, 'Location', 'northoutside', 'Orientation', 'horizontal', ...
              'NumColumns', 7);
        leg.Position(2) = leg.Position(2) + 0.02;
    end
end

% --- Excel Export ---

function exportResultsToExcel(T, results, config, meta, id_var, ssqi_vars)
    output_file = fullfile(config.output_folder, 'PPG_Analysis_Results.xlsx');

    % Sheet 1: Raw Data
    T_out = T;
    max_peaks = 2;

    for seg = 1:3
        seg_name = sprintf('seg%d', seg);
        for peak = 1:max_peaks
            prefix = sprintf('Seg%d_Peak%d', seg, peak);
            col_names = {[prefix '_Onset_s'], [prefix '_Systolic_s'], [prefix '_Notch_s'], ...
                         [prefix '_t_a_s'], [prefix '_a_norm'], ...
                         [prefix '_t_b_s'], [prefix '_b_norm'], ...
                         [prefix '_t_c_s'], [prefix '_c_norm'], ...
                         [prefix '_t_d_s'], [prefix '_d_norm'], ...
                         [prefix '_t_e_s'], [prefix '_e_norm']};
            for c = 1:length(col_names)
                T_out.(col_names{c}) = NaN(height(T), 1);
            end

            for i = 1:height(T)
                if ~results(i).(seg_name).processed
                    continue;
                end
                if peak <= length(results(i).(seg_name).onset_times)
                    T_out.([prefix '_Onset_s'])(i) = results(i).(seg_name).onset_times(peak);
                end
                if peak <= length(results(i).(seg_name).peak_times)
                    T_out.([prefix '_Systolic_s'])(i) = results(i).(seg_name).peak_times(peak);
                end
                if peak <= length(results(i).(seg_name).notch_times)
                    T_out.([prefix '_Notch_s'])(i) = results(i).(seg_name).notch_times(peak);
                end
                for lbl = {'a','b','c','d','e'}
                    t_field = ['t_' lbl{1}];
                    if peak <= length(results(i).(seg_name).(t_field))
                        T_out.([prefix '_t_' lbl{1} '_s'])(i) = results(i).(seg_name).(t_field)(peak);
                        T_out.([prefix '_' lbl{1} '_norm'])(i) = results(i).(seg_name).(lbl{1})(peak);
                    end
                end
            end
        end

        if seg < 3
            T_out.(sprintf('Empty_After_Seg%d', seg)) = repmat({''}, height(T), 1);
        end
    end

    writetable(T_out, output_file, 'Sheet', 'Raw_Data');
    fprintf('  Sheet 1 (Raw_Data): Saved\n');

    % Sheet 2: Vascular Indices (individual measurements)
    T_indices = createVascularIndicesTable(T, results, id_var, meta);
    writetable(T_indices, output_file, 'Sheet', 'Vascular_Indices');
    fprintf('  Sheet 2 (Vascular_Indices): %d rows\n', height(T_indices));

    % Sheet 3: Subject Averages
    T_averaged = createAveragedDataTable(T, results, id_var, meta);
    writetable(T_averaged, output_file, 'Sheet', 'Subject_Averages');
    fprintf('  Sheet 3 (Subject_Averages): %d rows\n', height(T_averaged));
end

function T_indices = createVascularIndicesTable(T, results, id_var, meta)
    % Create vascular indices table (one row per measurement).

    % Estimate upper bound for preallocation (3 segments * 2 peaks * n_subjects)
    max_rows = height(T) * 3 * 2;

    subject_ids = NaN(max_rows, 1);
    ages = NaN(max_rows, 1);
    sexes = cell(max_rows, 1);
    segments = NaN(max_rows, 1);
    peaks = NaN(max_rows, 1);
    ratio_b_a = NaN(max_rows, 1);
    ratio_c_a = NaN(max_rows, 1);
    ratio_d_a = NaN(max_rows, 1);
    ratio_e_a = NaN(max_rows, 1);
    agi_index = NaN(max_rows, 1);
    time_onset_to_peak = NaN(max_rows, 1);
    time_peak_to_notch = NaN(max_rows, 1);
    time_b_to_c = NaN(max_rows, 1);
    time_b_to_d = NaN(max_rows, 1);
    time_b_to_e = NaN(max_rows, 1);

    row_count = 0;

    for i = 1:height(T)
        subj_id = convertToNum(T.(id_var)(i));
        age_val = getMetaNum(T, i, meta.age_var);
        sex_val = getMetaSex(T, i, meta.sex_var);

        for seg = 1:3
            seg_name = sprintf('seg%d', seg);
            if ~results(i).(seg_name).processed
                continue;
            end

            n_pk = min(length(results(i).(seg_name).peak_times), 2);
            for pk = 1:n_pk
                row_count = row_count + 1;
                subject_ids(row_count) = subj_id;
                ages(row_count) = age_val;
                sexes{row_count} = char(sex_val);
                segments(row_count) = seg;
                peaks(row_count) = pk;

                [a_val, b_val, c_val, d_val, e_val] = getAPGValues(results(i).(seg_name), pk);

                if ~isnan(a_val) && a_val ~= 0
                    ratio_b_a(row_count) = b_val / a_val;
                    ratio_c_a(row_count) = c_val / a_val;
                    ratio_d_a(row_count) = d_val / a_val;
                    ratio_e_a(row_count) = e_val / a_val;
                    agi_index(row_count) = (b_val - c_val - d_val - e_val) / a_val;
                end

                [onset_t, peak_t, notch_t, t_b, t_c, t_d, t_e] = getTimingValues(results(i).(seg_name), pk);

                if ~isnan(onset_t) && ~isnan(peak_t)
                    time_onset_to_peak(row_count) = peak_t - onset_t;
                end
                if ~isnan(peak_t) && ~isnan(notch_t)
                    time_peak_to_notch(row_count) = notch_t - peak_t;
                end
                if ~isnan(t_b) && ~isnan(t_c)
                    time_b_to_c(row_count) = t_c - t_b;
                end
                if ~isnan(t_b) && ~isnan(t_d)
                    time_b_to_d(row_count) = t_d - t_b;
                end
                if ~isnan(t_b) && ~isnan(t_e)
                    time_b_to_e(row_count) = t_e - t_b;
                end
            end
        end
    end

    % Trim to actual size
    T_indices = table(subject_ids(1:row_count), ages(1:row_count), ...
                     sexes(1:row_count), segments(1:row_count), peaks(1:row_count), ...
                     ratio_b_a(1:row_count), ratio_c_a(1:row_count), ...
                     ratio_d_a(1:row_count), ratio_e_a(1:row_count), ...
                     agi_index(1:row_count), ...
                     time_onset_to_peak(1:row_count), time_peak_to_notch(1:row_count), ...
                     time_b_to_c(1:row_count), time_b_to_d(1:row_count), ...
                     time_b_to_e(1:row_count), ...
                     'VariableNames', {'Subject_ID', 'Age', 'Sex', 'Segment', 'Peak', ...
                                      'b_a_ratio', 'c_a_ratio', 'd_a_ratio', 'e_a_ratio', ...
                                      'AGI', ...
                                      'Time_Onset_to_Peak_s', 'Time_Peak_to_Notch_s', ...
                                      'Time_b_to_c_s', 'Time_b_to_d_s', 'Time_b_to_e_s'});
end

function T_averaged = createAveragedDataTable(T, results, id_var, meta)
    % Create averaged data table (one row per subject).
    n_subjects = height(T);

    subject_ids = NaN(n_subjects, 1);
    ages = NaN(n_subjects, 1);
    sexes = cell(n_subjects, 1);
    hypertension = cell(n_subjects, 1);
    diabetes = cell(n_subjects, 1);
    cerebral_infarction = cell(n_subjects, 1);
    n_measurements = zeros(n_subjects, 1);

    avg_b_a = NaN(n_subjects, 1);
    avg_c_a = NaN(n_subjects, 1);
    avg_d_a = NaN(n_subjects, 1);
    avg_e_a = NaN(n_subjects, 1);
    avg_agi = NaN(n_subjects, 1);
    avg_time_onset_peak = NaN(n_subjects, 1);
    avg_time_peak_notch = NaN(n_subjects, 1);
    avg_time_b_c = NaN(n_subjects, 1);
    avg_time_b_d = NaN(n_subjects, 1);
    avg_time_b_e = NaN(n_subjects, 1);

    for i = 1:n_subjects
        subject_ids(i) = convertToNum(T.(id_var)(i));
        ages(i) = getMetaNum(T, i, meta.age_var);
        sexes{i} = char(getMetaSex(T, i, meta.sex_var));
        hypertension{i} = getMetaStr(T, i, meta.hyper_var);
        diabetes{i} = getMetaStr(T, i, meta.diab_var);
        cerebral_infarction{i} = getMetaStr(T, i, meta.cinf_var);

        % Collect all measurements for this subject
        all_b_a = []; all_c_a = []; all_d_a = []; all_e_a = []; all_agi = [];
        all_time_op = []; all_time_pn = [];
        all_time_bc = []; all_time_bd = []; all_time_be = [];

        for seg = 1:3
            seg_name = sprintf('seg%d', seg);
            if ~results(i).(seg_name).processed
                continue;
            end
            n_pk = min(length(results(i).(seg_name).peak_times), 2);
            for pk = 1:n_pk
                [a_val, b_val, c_val, d_val, e_val] = getAPGValues(results(i).(seg_name), pk);
                if ~isnan(a_val) && a_val ~= 0
                    all_b_a(end+1) = b_val / a_val;
                    all_c_a(end+1) = c_val / a_val;
                    all_d_a(end+1) = d_val / a_val;
                    all_e_a(end+1) = e_val / a_val;
                    all_agi(end+1) = (b_val - c_val - d_val - e_val) / a_val;
                end

                [onset_t, peak_t, notch_t, t_b, t_c, t_d, t_e] = getTimingValues(results(i).(seg_name), pk);
                if ~isnan(onset_t) && ~isnan(peak_t)
                    all_time_op(end+1) = peak_t - onset_t;
                end
                if ~isnan(peak_t) && ~isnan(notch_t)
                    all_time_pn(end+1) = notch_t - peak_t;
                end
                if ~isnan(t_b) && ~isnan(t_c)
                    all_time_bc(end+1) = t_c - t_b;
                end
                if ~isnan(t_b) && ~isnan(t_d)
                    all_time_bd(end+1) = t_d - t_b;
                end
                if ~isnan(t_b) && ~isnan(t_e)
                    all_time_be(end+1) = t_e - t_b;
                end
            end
        end

        n_measurements(i) = length(all_b_a);
        avg_b_a(i) = mean(all_b_a, 'omitnan');
        avg_c_a(i) = mean(all_c_a, 'omitnan');
        avg_d_a(i) = mean(all_d_a, 'omitnan');
        avg_e_a(i) = mean(all_e_a, 'omitnan');
        avg_agi(i) = mean(all_agi, 'omitnan');
        avg_time_onset_peak(i) = mean(all_time_op, 'omitnan');
        avg_time_peak_notch(i) = mean(all_time_pn, 'omitnan');
        avg_time_b_c(i) = mean(all_time_bc, 'omitnan');
        avg_time_b_d(i) = mean(all_time_bd, 'omitnan');
        avg_time_b_e(i) = mean(all_time_be, 'omitnan');
    end

    T_averaged = table(subject_ids, ages, sexes, hypertension, diabetes, ...
                      cerebral_infarction, n_measurements, ...
                      avg_b_a, avg_c_a, avg_d_a, avg_e_a, avg_agi, ...
                      avg_time_onset_peak, avg_time_peak_notch, ...
                      avg_time_b_c, avg_time_b_d, avg_time_b_e, ...
                      'VariableNames', {'Subject_ID', 'Age', 'Sex', 'Hypertension', ...
                                       'Diabetes', 'Cerebral_Infarction', 'N_Measurements', ...
                                       'Mean_b_a', 'Mean_c_a', 'Mean_d_a', 'Mean_e_a', ...
                                       'Mean_AGI', ...
                                       'Mean_Time_Onset_Peak_s', 'Mean_Time_Peak_Notch_s', ...
                                       'Mean_Time_b_c_s', 'Mean_Time_b_d_s', 'Mean_Time_b_e_s'});
end

% --- Correlation Plots ---

function createCorrelationPlots(T, results, config, meta, id_var)
    fprintf('  Creating correlation plots...\n');

    T_indices = createVascularIndicesTable(T, results, id_var, meta);
    [clinical_categories, category_names, colors] = getClinicalCategories(T, meta);

    params_to_track = {'d_a_ratio', 'AGI', 'b_a_ratio', 'c_a_ratio', 'e_a_ratio'};

    params = {
        'b_a_ratio',             'b_a_Ratio',        'b/a Ratio',                  'b/a'
        'c_a_ratio',             'c_a_Ratio',        'c/a Ratio',                  'c/a'
        'd_a_ratio',             'd_a_Ratio',        'd/a Ratio',                  'd/a'
        'e_a_ratio',             'e_a_Ratio',        'e/a Ratio',                  'e/a'
        'AGI',                   'AGI',              'Aging Index',                'AGI = (b-c-d-e)/a'
        'Time_Onset_to_Peak_s',  'Time_Onset_Peak',  'Time Onset to Peak (s)',     'Upstroke Time'
        'Time_Peak_to_Notch_s',  'Time_Peak_Notch',  'Time Peak to Notch (s)',     'Systolic Duration'
        'Time_b_to_c_s',         'Time_b_c',         'Time b to c (s)',            'Early Diastolic Timing'
        'Time_b_to_d_s',         'Time_b_d',         'Time b to d (s)',            'Dicrotic Wave Timing'
        'Time_b_to_e_s',         'Time_b_e',         'Time b to e (s)',            'Late Diastolic Timing'
    };

    for p = 1:size(params, 1)
        param_col = params{p, 1};
        param_short = params{p, 2};
        param_full = params{p, 3};
        param_formula = params{p, 4};

        age_data = T_indices.Age(:);
        param_data = T_indices.(param_col)(:);

        valid_idx = ~isnan(age_data) & ~isnan(param_data);
        age_clean = age_data(valid_idx);
        param_clean = param_data(valid_idx);

        if length(age_clean) < 3
            fprintf('    Skipped %s: insufficient data\n', param_short);
            continue;
        end

        % Outlier removal (mean +/- threshold * std)
        param_mean = mean(param_clean);
        param_std = std(param_clean);
        subject_ids_valid = T_indices.Subject_ID(valid_idx);

        if param_std > 0
            outlier_mask = (param_clean >= param_mean - config.outlier_std_threshold * param_std) & ...
                          (param_clean <= param_mean + config.outlier_std_threshold * param_std);

            n_outliers = sum(~outlier_mask);
            if n_outliers > 0 && any(strcmp(param_col, params_to_track))
                outlier_ids = unique(subject_ids_valid(~outlier_mask));
                fprintf('    %s: Outliers removed (IDs): ', param_short);
                fprintf('%d ', outlier_ids);
                fprintf('\n');
            end

            age_clean = age_clean(outlier_mask);
            param_clean = param_clean(outlier_mask);
            subject_ids_valid = subject_ids_valid(outlier_mask);
        end

        if length(age_clean) < 3 || std(age_clean) == 0 || std(param_clean) == 0
            fprintf('    Skipped %s: insufficient data after outlier removal\n', param_short);
            continue;
        end

        % Map subjects to clinical categories
        categories_valid = zeros(size(subject_ids_valid));
        for i = 1:length(subject_ids_valid)
            subj_row = find(convertToNum(T.(id_var)) == subject_ids_valid(i), 1);
            if ~isempty(subj_row)
                categories_valid(i) = clinical_categories(subj_row);
            end
        end

        % Global statistics
        p_fit_global = polyfit(age_clean, param_clean, 1);
        [r_pearson, p_value] = corr(age_clean, param_clean);

        fig = figure('Position', [100, 100, 800, 600], 'Color', 'w');
        hold on;

        for cat = 1:length(category_names)
            cat_idx = categories_valid == cat;
            if any(cat_idx)
                scatter(age_clean(cat_idx), param_clean(cat_idx), 80, ...
                       colors(cat, :), 'filled', 'MarkerEdgeColor', 'k', ...
                       'LineWidth', 0.5, 'DisplayName', category_names{cat});

                age_cat = age_clean(cat_idx);
                param_cat = param_clean(cat_idx);
                if length(age_cat) >= 2 && std(age_cat) > 0
                    p_fit_cat = polyfit(age_cat, param_cat, 1);
                    age_range = linspace(min(age_cat), max(age_cat), 100);
                    plot(age_range, polyval(p_fit_cat, age_range), '-', ...
                         'Color', colors(cat, :) * 0.7, 'LineWidth', 2.5, ...
                         'DisplayName', [category_names{cat} ' - fit']);
                end
            end
        end
        hold off;

        xlabel('Age (years)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(param_full, 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('%s vs Age', param_full), 'FontSize', 14, 'FontWeight', 'bold');

        if p_value < 0.001
            p_str = 'p < 0.001 ***';
        elseif p_value < 0.01
            p_str = sprintf('p = %.3f **', p_value);
        elseif p_value < 0.05
            p_str = sprintf('p = %.3f *', p_value);
        else
            p_str = sprintf('p = %.4f (n.s.)', p_value);
        end

        stats_text = sprintf(['Formula: %s\n' ...
                             'Global linear fit: y = %.4fx + %.4f\n' ...
                             'Pearson r = %.3f\n' ...
                             '%s\n' ...
                             'n = %d'], ...
                           param_formula, p_fit_global(1), p_fit_global(2), ...
                           r_pearson, p_str, length(age_clean));
        text(0.05, 0.95, stats_text, 'Units', 'normalized', ...
             'VerticalAlignment', 'top', 'FontSize', 9, ...
             'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 5, ...
             'FontWeight', 'bold');

        legend('Location', 'best', 'FontSize', 8);
        grid on;
        set(gca, 'FontSize', 11, 'LineWidth', 1);
        box on;

        fig_name = sprintf('Correlation_%s_vs_Age', param_short);
        print(fig, fullfile(config.correlation_folder, [fig_name '.png']), '-dpng', '-r300');
        savefig(fig, fullfile(config.correlation_folder, [fig_name '.fig']));

        fprintf('    Plotted: %s (r=%.3f, p=%.4f)\n', param_short, r_pearson, p_value);
    end
end

% --- Clinical Categories ---

function [categories, category_names, colors] = getClinicalCategories(T, meta)
    n_subjects = height(T);
    categories = ones(n_subjects, 1);

    for i = 1:n_subjects
        has_hyper = false;
        has_diab = false;
        has_cinf = false;

        if meta.hyper_var ~= ""
            hyper_str = lower(string(T.(meta.hyper_var)(i)));
            if ~ismissing(hyper_str) && hyper_str ~= "" && ...
               ~contains(hyper_str, "normal") && ~contains(hyper_str, "normo")
                has_hyper = true;
            end
        end

        if meta.diab_var ~= ""
            diab_str = lower(string(T.(meta.diab_var)(i)));
            if ~ismissing(diab_str) && diab_str ~= "" && ...
               (contains(diab_str, "diabetes") || contains(diab_str, "type 2"))
                has_diab = true;
            end
        end

        if meta.cinf_var ~= ""
            cinf_str = string(T.(meta.cinf_var)(i));
            if ~ismissing(cinf_str) && cinf_str ~= ""
                has_cinf = true;
            end
        end

        if ~has_hyper && ~has_diab && ~has_cinf
            categories(i) = 1;  % Healthy
        elseif has_diab
            categories(i) = 3;  % Diabetes
        elseif has_hyper && has_cinf
            categories(i) = 4;  % Multiple conditions
        elseif has_hyper
            categories(i) = 2;  % Hypertension only
        else
            categories(i) = 4;  % Other
        end
    end

    category_names = {'Healthy', 'Hypertension', 'Diabetes', 'Multiple Conditions'};
    colors = [
        0.00, 0.45, 0.74;  % Blue
        0.85, 0.33, 0.10;  % Orange
        0.93, 0.11, 0.14;  % Red
        0.49, 0.18, 0.56   % Purple
    ];
end

% --- Metadata & Utility Functions ---

function metadata = extractMetadata(T, row, meta)
    metadata = struct();
    if meta.sex_var ~= ""
        metadata.sex = normalizeSex(string(T.(meta.sex_var)(row)));
    else
        metadata.sex = "";
    end
    if meta.age_var ~= ""
        metadata.age = convertToNum(T.(meta.age_var)(row));
    else
        metadata.age = NaN;
    end
end

function val = getMetaNum(T, row, var_name)
    if var_name ~= ""
        val = convertToNum(T.(var_name)(row));
    else
        val = NaN;
    end
end

function val = getMetaSex(T, row, var_name)
    if var_name ~= ""
        val = normalizeSex(string(T.(var_name)(row)));
    else
        val = "";
    end
end

function val = getMetaStr(T, row, var_name)
    if var_name ~= ""
        val = char(string(T.(var_name)(row)));
    else
        val = '';
    end
end

function [a_val, b_val, c_val, d_val, e_val] = getAPGValues(seg_data, pk)
    % Extract APG amplitude values for a given peak.
    labels = {'a', 'b', 'c', 'd', 'e'};
    vals = NaN(1, 5);
    for k = 1:5
        if pk <= length(seg_data.(labels{k}))
            vals(k) = seg_data.(labels{k})(pk);
        end
    end
    a_val = vals(1); b_val = vals(2); c_val = vals(3);
    d_val = vals(4); e_val = vals(5);
end

function [onset_t, peak_t, notch_t, t_b, t_c, t_d, t_e] = getTimingValues(seg_data, pk)
    % Extract timing values for a given peak.
    onset_t = NaN; peak_t = NaN; notch_t = NaN;
    t_b = NaN; t_c = NaN; t_d = NaN; t_e = NaN;

    if pk <= length(seg_data.onset_times), onset_t = seg_data.onset_times(pk); end
    if pk <= length(seg_data.peak_times),  peak_t = seg_data.peak_times(pk);   end
    if pk <= length(seg_data.notch_times), notch_t = seg_data.notch_times(pk); end
    if pk <= length(seg_data.t_b), t_b = seg_data.t_b(pk); end
    if pk <= length(seg_data.t_c), t_c = seg_data.t_c(pk); end
    if pk <= length(seg_data.t_d), t_d = seg_data.t_d(pk); end
    if pk <= length(seg_data.t_e), t_e = seg_data.t_e(pk); end
end

function s = normalizeSex(sex)
    x = lower(strtrim(sex));
    if contains(x, "female")
        s = "F";
    elseif contains(x, "male")
        s = "M";
    else
        s = sex;
    end
end

function ppg = loadPPGFile(fname)
    if ~exist(fname, 'file')
        ppg = [];
        return;
    end
    m = readmatrix(fname);
    m = m(:);
    m = m(~isnan(m));
    ppg = m;
end

function vname = findColumn(T, exact_preferred, contains_all)
    vars = string(T.Properties.VariableNames);
    vlow = lower(vars);
    for k = 1:numel(exact_preferred)
        idx = find(vlow == exact_preferred(k), 1);
        if ~isempty(idx)
            vname = vars(idx);
            return;
        end
    end
    idx = 1:numel(vars);
    for k = 1:numel(contains_all)
        idx = idx(contains(vlow(idx), contains_all(k)));
    end
    if ~isempty(idx)
        vname = vars(idx(1));
        return;
    end
    error('Required column not found.');
end

function vname = findColumnOptional(T, exact_preferred, contains_all)
    try
        vname = findColumn(T, exact_preferred, contains_all);
    catch
        vname = "";
    end
end

function [ssqi_vars, ok] = findSSQIColumns(T)
    vars = string(T.Properties.VariableNames);
    vlow = lower(vars);
    ssqi_vars = strings(1, 3);
    ok = false;

    cand = vars(contains(vlow, "ssqi") | contains(vlow, "segment"));
    if numel(cand) >= 3
        ssqi_vars = cand(1:3);
        ok = true;
        return;
    end

    s1 = find(contains(vlow, "segment1") | contains(vlow, "segment_1") | ...
              contains(vlow, "segment 1"), 1);
    s2 = find(contains(vlow, "segment2") | contains(vlow, "segment_2") | ...
              contains(vlow, "segment 2"), 1);
    s3 = find(contains(vlow, "segment3") | contains(vlow, "segment_3") | ...
              contains(vlow, "segment 3"), 1);

    if ~isempty(s1) && ~isempty(s2) && ~isempty(s3)
        ssqi_vars = [vars(s1), vars(s2), vars(s3)];
        ok = true;
    end
end

function y = convertToNum(x)
    if iscell(x)
        x = x{1};
    end
    if isstring(x) || ischar(x)
        s = string(x);
        s = replace(s, ",", ".");
        y = str2double(s);
    elseif isnumeric(x)
        y = double(x);
    else
        y = NaN;
    end
end
