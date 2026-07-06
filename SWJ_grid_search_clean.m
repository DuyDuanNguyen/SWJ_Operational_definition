%% SWJ_GRID_SEARCH_CLEAN.M
% =========================================================================
%  Square-wave jerk (SWJ) grid-search pipeline
%
%  Purpose
%  -------
%  This script detects right-eye saccades during attempted fixation and then
%  performs a grid sweep of candidate SWJ operational definitions. It is
%  designed for transparent, reproducible sharing with a manuscript/GitHub
%  repository.
%
%  Core workflow
%  -------------
%    1. Read each exported eye-tracking CSV file.
%    2. Extract gaze, velocity, and acceleration traces.
%    3. Mark blink/signal-loss/artifact periods.
%    4. Detect saccades and export Saccades_right_main.
%    5. Sweep candidate SWJ criteria:
%          - saccade amplitude range
%          - intra-SWJ interval range
%          - directional opposition tolerance
%          - amplitude similarity threshold
%    6. Export:
%          - per-file threshold sweep CSV
%          - master threshold sweep CSV across all files
%          - event-level CSV for the default SWJ definition
%          - participant/file-level summary for the default SWJ definition
%          - optional Saccades_right_main CSV for quality control/reuse
%
%  Saccades_right_main columns
%  ---------------------------
%    1  StartSample
%    2  EndSample
%    3  StartX_px
%    4  EndX_px
%    5  StartY_px
%    6  EndY_px
%    7  PeakVelSample
%    8  PeakVel_deg_s
%    9  MaxAcc_deg_s2
%    10 MaxDec_deg_s2
%    11 Duration_ms
%    12 AmpX_deg
%    13 AmpY_deg
%    14 AmpTotal_deg
%    15 Direction_deg
%    16 SacPart
%    17 ClusterID
%    18 VerticalComponent
%    19 Trial
%    20 Reserved
%    21 Reserved
%
%  Vertical component
%  ------------------
%    VerticalComponent = abs(vertical amplitude) / total amplitude
%    0 = purely horizontal; 1 = purely vertical.
%
%
%  Required input columns
%  ----------------------
%    GazeX, GazeY,
%    RIGHT_VELOCITY_X, RIGHT_VELOCITY_Y,
%    RIGHT_ACCELERATION_X, RIGHT_ACCELERATION_Y
%
%  Optional input columns
%  ----------------------
%    TRIAL or Trial: used to prevent SWJ pairs crossing trial boundaries.
% =========================================================================

clear; close all; clc;

%% ========================================================================
% USER SETTINGS
% ========================================================================

cfg = struct();

% ---- Paths ----
cfg.inputFolder  = '...';  % folder containing input *.csv files
cfg.outputFolder = '...';  % folder for output CSV files

if ~exist(cfg.outputFolder, 'dir')
    mkdir(cfg.outputFolder);
end

% ---- Eye-tracker / screen setup ----
% Change these values if applying the same code to another acquisition setup.
cfg.fs        = 500;   % sampling frequency, Hz
cfg.widthPix  = 1024;
cfg.heightPix = 768;
cfg.widthCm   = 41.6;
cfg.heightCm  = 23.5;
cfg.viewDistCm = 60;

% ---- Blink/artifact handling ----
cfg.preBlinkBuffer_ms  = 100;
cfg.postBlinkBuffer_ms = 250;
cfg.artifactBuffer_ms  = 10;

cfg.velocityArtifact_deg_s = 900;
cfg.accArtifact_deg_s2     = 100000;

% ---- Saccade detection ----
cfg.accelerationSDMultiplier = 2.576;
cfg.mergeGap_ms              = 20;
cfg.refineWindow_ms          = 500;
cfg.minSaccadeDuration_ms    = 8;
cfg.minSaccadeAmplitude_deg  = 0.5;
cfg.velocityOffset_deg_s     = 5;
cfg.maxDirectionDeviation_deg = 60;
cfg.maxSampleDirectionChange_deg = 20;

% ---- SWJ pairing ----
% 0 = only consecutive detected saccades can form an SWJ pair.
% 1 = allow one intervening detected saccade between the first and return
%     component. Use the setting that matches the final analysis.
cfg.maxInterveningSaccades = 0;

% ---- Export options ----
cfg.exportPerFileSweep     = true;
cfg.exportDefaultEvents    = true;
cfg.exportSaccadesRightMain = true;

%% ========================================================================
% GRID SEARCH DEFINITION
% ========================================================================

GRID = struct();

GRID.ampMinList = [0.5 1];
GRID.ampMaxList = 8;

GRID.isiWindows_ms = [ ...
    20 300; ...
    20 400; ...
    20 500; ...
    20 600];

GRID.dirTolList_deg = [15 45];

% Current manuscript grid uses the dissimilarity index:
% |A1 - A2| / |A1 + A2|
GRID.simTypeList = "dissim";
GRID.simThrList  = [0.25 0.5];

% Default definition used for event-level and default-summary exports.
% Set this to the final operational definition selected from the grid search.
DEFAULT = struct();
DEFAULT.amp_min  = 1;
DEFAULT.amp_max  = 8;
DEFAULT.isi_min  = 20;
DEFAULT.isi_max  = 600;
DEFAULT.dir_tol  = 45;
DEFAULT.sim_type = "dissim";
DEFAULT.sim_thr  = 0.25;

%% ========================================================================
% MAIN LOOP
% ========================================================================

fileList = dir(fullfile(cfg.inputFolder, '*.csv'));

AllThresholdResults = table();
Summary_SWJ_default = table();

for fileIdx = 1:numel(fileList)

    inputFile = fullfile(cfg.inputFolder, fileList(fileIdx).name);
    [~, filename, ~] = fileparts(inputFile);

    fprintf('\n[%d/%d] Processing %s\n', fileIdx, numel(fileList), filename);

    % ---------------------------------------------------------------------
    % 1. Read file
    % ---------------------------------------------------------------------
    T = readEyeTrackingCSV(inputFile);

    % ---------------------------------------------------------------------
    % 2. Detect right-eye saccades
    % ---------------------------------------------------------------------
    [Saccades_right_main, validGaze] = detectRightMainSaccades(T, cfg);

    durSec = sum(validGaze) / cfg.fs;
    durMin = durSec / 60;

    % Optional export of the detected main saccades.
    if cfg.exportSaccadesRightMain
        outSaccades = fullfile(cfg.outputFolder, sprintf('%s_Saccades_right_main.csv', filename));
        writetable(saccadesMainToTable(Saccades_right_main), outSaccades);
    end

    % ---------------------------------------------------------------------
    % 3. Grid sweep
    % ---------------------------------------------------------------------
    FileThresholdRows = sweepSWJGrid( ...
        filename, fileIdx, durSec, durMin, ...
        Saccades_right_main, validGaze, GRID, cfg);

    AllThresholdResults = [AllThresholdResults; FileThresholdRows]; %#ok<AGROW>

    if cfg.exportPerFileSweep
        outFileSweep = fullfile(cfg.outputFolder, sprintf('%s_SWJ_threshold_sweep.csv', filename));
        writetable(FileThresholdRows, outFileSweep);
    end

    % ---------------------------------------------------------------------
    % 4. Default-threshold event-level and summary export
    % ---------------------------------------------------------------------
    [SWJ_pairs_def, SWJ_details_def] = detectSWJPairs( ...
        Saccades_right_main, DEFAULT, validGaze, cfg);

    if cfg.exportDefaultEvents
        outEvents = fullfile(cfg.outputFolder, sprintf('%s_SWJ_default_events.csv', filename));
        writetable(swjDetailsToTable(SWJ_details_def), outEvents);
    end

    defaultMetrics = summariseSWJDetails(SWJ_details_def, durMin);

    thisSummary = table( ...
        string(filename), fileIdx, durSec, ...
        defaultMetrics.SWJ_Count, ...
        defaultMetrics.SWJ_Rate_per_min, ...
        defaultMetrics.MeanPairAmp_deg, ...
        defaultMetrics.Total_Amp, ...
        defaultMetrics.MeanISI_ms, ...
        defaultMetrics.MeanPairVerticalComp, ...
        'VariableNames', { ...
            'File', 'FileIndex', 'Dur_sec', ...
            'SWJ_Count', 'SWJ_Rate_per_min', ...
            'MeanPairAmp_deg', 'Total_Amp', ...
            'MeanISI_ms', 'MeanPairVerticalComp'} ...
        );

    Summary_SWJ_default = [Summary_SWJ_default; thisSummary]; %#ok<AGROW>
end

%% ========================================================================
% MASTER EXPORTS
% ========================================================================

outMaster = fullfile(cfg.outputFolder, 'SWJ_threshold_sweep_all_files.csv');
writetable(AllThresholdResults, outMaster);
fprintf('\nSaved master threshold sweep:\n%s\n', outMaster);

outDefaultSummary = fullfile(cfg.outputFolder, 'SWJ_default_threshold_summary_all_files.csv');
writetable(Summary_SWJ_default, outDefaultSummary);
fprintf('Saved default-threshold summary:\n%s\n', outDefaultSummary);


%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function T = readEyeTrackingCSV(inputFile)
% Read an eye-tracking CSV robustly. Some exports contain row names and some
% do not, so this function first tries the row-name format and then falls
% back to the standard readtable call.

    try
        T = readtable(inputFile, 'TreatAsMissing', '.', 'ReadRowNames', true);
    catch
        T = readtable(inputFile, 'TreatAsMissing', '.');
    end
end


function [Saccades_right_main, validGaze] = detectRightMainSaccades(T, cfg)
% Detect right-eye main saccades using the DEMoNS-style approach:
% acceleration-based approximate intervals, gap merging, peak-velocity
% candidate identification, and onset/offset refinement using velocity and
% direction criteria.

    % ---------------------------------------------------------------------
    % Required traces
    % ---------------------------------------------------------------------
    gazeX = getRequiredColumn(T, "GazeX");
    gazeY = getRequiredColumn(T, "GazeY");

    velX = getRequiredColumn(T, "RIGHT_VELOCITY_X");
    velY = getRequiredColumn(T, "RIGHT_VELOCITY_Y");

    accX = getRequiredColumn(T, "RIGHT_ACCELERATION_X");
    accY = getRequiredColumn(T, "RIGHT_ACCELERATION_Y");

    n = height(T);

    velTotal = hypot(velX, velY);
    accTotal = [NaN; diff(velTotal)] * cfg.fs;

    accXpos = accX;
    accYpos = accY;
    accXpos(velX < 0) = -accXpos(velX < 0);
    accYpos(velY < 0) = -accYpos(velY < 0);

    % ---------------------------------------------------------------------
    % Blink and artifact exclusion
    % ---------------------------------------------------------------------
    blinkMask = isnan(gazeX) | isnan(gazeY);
    blinkMask = expandMask( ...
        blinkMask, ...
        ms2samples(cfg.preBlinkBuffer_ms, cfg.fs), ...
        ms2samples(cfg.postBlinkBuffer_ms, cfg.fs));

    artifactMask = abs(velX) >= cfg.velocityArtifact_deg_s | ...
                   abs(velY) >= cfg.velocityArtifact_deg_s | ...
                   abs(accXpos) >= cfg.accArtifact_deg_s2 | ...
                   abs(accYpos) >= cfg.accArtifact_deg_s2;

    artifactMask = expandMask( ...
        artifactMask, ...
        ms2samples(cfg.artifactBuffer_ms, cfg.fs), ...
        ms2samples(cfg.artifactBuffer_ms, cfg.fs));

    invalidMask = blinkMask | artifactMask | ...
                  isnan(gazeX) | isnan(gazeY) | ...
                  isnan(velX)  | isnan(velY)  | ...
                  (velX == 0 & accX == 0);

    gazeX_clean = gazeX;
    gazeY_clean = gazeY;

    gazeX_clean(invalidMask) = NaN;
    gazeY_clean(invalidMask) = NaN;

    velX(invalidMask) = NaN;
    velY(invalidMask) = NaN;
    velTotal(invalidMask) = NaN;

    accX(invalidMask) = NaN;
    accY(invalidMask) = NaN;
    accTotal(invalidMask) = NaN;

    validGaze = ~isnan(gazeX_clean) & ~isnan(gazeY_clean);

    % ---------------------------------------------------------------------
    % Approximate saccade intervals
    % ---------------------------------------------------------------------
    sdAccX = std(accX, 'omitnan');
    sdAccY = std(accY, 'omitnan');

    appSac = abs(accX) > cfg.accelerationSDMultiplier * sdAccX | ...
             abs(accY) > cfg.accelerationSDMultiplier * sdAccY;

    appSac(~validGaze) = false;
    appSac = mergeShortGaps(appSac, ms2samples(cfg.mergeGap_ms, cfg.fs));

    AS = buildApproxSaccadeParameterMatrix(appSac, velTotal, velX, velY);

    if isempty(AS)
        Saccades_right_main = zeros(0, 21);
        return;
    end

    % ---------------------------------------------------------------------
    % Refine onset/offset and construct candidate saccades
    % ---------------------------------------------------------------------
    Saccades_right = refineSaccadesFromApproxIntervals( ...
        AS, gazeX_clean, gazeY_clean, velX, velY, velTotal, accTotal, cfg);

    if isempty(Saccades_right)
        Saccades_right_main = zeros(0, 21);
        return;
    end

    % ---------------------------------------------------------------------
    % Add kinematic features and select main saccades
    % ---------------------------------------------------------------------
    Saccades_right = addSaccadeMetrics( ...
        Saccades_right, T, velX, velY, cfg);

    keep = Saccades_right(:, 11) >= cfg.minSaccadeDuration_ms & ...
           Saccades_right(:, 14) >= cfg.minSaccadeAmplitude_deg;

    Saccades_right_long = Saccades_right(keep, :);

    if isempty(Saccades_right_long)
        Saccades_right_main = zeros(0, 21);
        return;
    end

    Saccades_right_main = selectMainSaccadePerCluster(Saccades_right_long, cfg);
end


function AS = buildApproxSaccadeParameterMatrix(appSac, velTotal, velX, velY)
% Construct approximate saccade intervals and store up to 12 local
% peak-velocity candidates per interval. This preserves the structure of
% the original analysis while making the code easier to read.

    n = numel(appSac);
    candidateCols = [3, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58];

    AS = zeros(5000, 62);
    nIntervals = 0;
    maxVelSac = zeros(5, 12);
    maxVelSac(5, :) = candidateCols;

    for i = 2:(n - 1)

        startsNewInterval = appSac(i) && ~appSac(i - 1);

        if startsNewInterval
            nIntervals = nIntervals + 1;
            AS(nIntervals, 1) = i;

            maxVelSac = zeros(5, 12);
            maxVelSac(5, :) = candidateCols;
            maxVelSac(1:4, 1) = [velTotal(i); i; velX(i); velY(i)];
        end

        if appSac(i) && nIntervals > 0 && AS(nIntervals, 1) ~= 0

            if AS(nIntervals, 3) == 0 && velTotal(i) > maxVelSac(1, 1)
                maxVelSac(1:4, 1) = [velTotal(i); i; velX(i); velY(i)];
            end

            isIntervalEnd = appSac(i) && ~appSac(i + 1);

            if AS(nIntervals, 3) == 0 && ...
               (velTotal(i) < maxVelSac(1, 1) || isIntervalEnd)

                AS(nIntervals, 3:6) = maxVelSac(1:4, 1);
                AS(nIntervals, 7) = directionDeg( ...
                    velX(maxVelSac(2, 1)), velY(maxVelSac(2, 1)));
            end

            for k = 2:12
                previousCol = maxVelSac(5, k - 1);
                currentCol  = maxVelSac(5, k);

                if AS(nIntervals, previousCol) == 0
                    break;
                end

                foundRisingPeak = ...
                    AS(nIntervals, previousCol) ~= 0 && ...
                    AS(nIntervals, currentCol) == 0 && ...
                    i > AS(nIntervals, previousCol + 1) && ...
                    velTotal(i) > velTotal(i - 1);

                if foundRisingPeak
                    maxVelSac(1:4, k) = [velTotal(i); i; velX(i); velY(i)];
                end

                if maxVelSac(1, k) ~= 0 && ...
                   AS(nIntervals, currentCol) == 0 && ...
                   velTotal(i) > maxVelSac(1, k)

                    maxVelSac(1:4, k) = [velTotal(i); i; velX(i); velY(i)];
                end

                isIntervalEnd = appSac(i) && ~appSac(i + 1);

                if maxVelSac(1, k) ~= 0 && ...
                   AS(nIntervals, currentCol) == 0 && ...
                   (velTotal(i) < maxVelSac(1, k) || isIntervalEnd)

                    AS(nIntervals, currentCol:(currentCol + 3)) = maxVelSac(1:4, k);
                    AS(nIntervals, currentCol + 4) = directionDeg( ...
                        velX(maxVelSac(2, k)), velY(maxVelSac(2, k)));
                end
            end
        end

        if appSac(i) && ~appSac(i + 1) && nIntervals > 0 && AS(nIntervals, 1) ~= 0
            AS(nIntervals, 2) = i;
        end
    end

    AS((AS(:, 1) == 0 | AS(:, 2) == 0), :) = [];
end


function Saccades_right = refineSaccadesFromApproxIntervals( ...
    AS, gazeX, gazeY, velX, velY, velTotal, accTotal, cfg)
% Refine the onset/offset of each candidate saccade from its peak-velocity
% sample.

    candidateCols = [3, 8, 13, 18, 23, 28, 33, 38, 43, 48, 53, 58];

    nMax = max(5000, size(AS, 1) * numel(candidateCols));
    Saccades_right = zeros(nMax, 21);

    sac = 0;
    nSamples = numel(gazeX);
    refineWindow = ms2samples(cfg.refineWindow_ms, cfg.fs);

    for i = 1:size(AS, 1)

        sacPart = 9;

        for col = candidateCols

            if AS(i, col) == 0
                break;
            end

            peakSample = AS(i, col + 1);
            mainDir    = AS(i, col + 4);

            if peakSample < 2 || peakSample > (nSamples - 1)
                continue;
            end

            maxVel = 0;
            maxAcc = 0;
            maxDec = 0;
            tMaxVel = 0;

            sac = sac + 1;
            sacPart = sacPart + 1;

            % ---- Onset search ----
            startSearch = max(2, peakSample - refineWindow);

            for k = peakSample:-1:startSearch

                if isnan(velX(k))
                    break;
                end

                sampleDir       = directionDeg(velX(k), velY(k));
                previousSampleDir = directionDeg(velX(k - 1), velY(k - 1));

                if velTotal(k) > maxVel
                    maxVel = velTotal(k);
                    tMaxVel = k;
                end

                if accTotal(k) > maxAcc
                    maxAcc = accTotal(k);
                end

                if accTotal(k) < maxDec
                    maxDec = accTotal(k);
                end

                onsetCriterion = ...
                    velTotal(k - 1) < cfg.velocityOffset_deg_s || ...
                    angleDiffDeg(previousSampleDir, mainDir) > cfg.maxDirectionDeviation_deg || ...
                    angleDiffDeg(sampleDir, previousSampleDir) > cfg.maxSampleDirectionChange_deg;

                if onsetCriterion
                    Saccades_right(sac, 1) = k;
                    Saccades_right(sac, 3) = gazeX(k);
                    Saccades_right(sac, 5) = gazeY(k);
                    break;
                end
            end

            if Saccades_right(sac, 1) == 0
                continue;
            end

            % ---- Offset search ----
            endSearch = min(nSamples - 1, peakSample + refineWindow);

            for k = peakSample:endSearch

                if isnan(velX(k))
                    break;
                end

                sampleDir     = directionDeg(velX(k), velY(k));
                nextSampleDir = directionDeg(velX(k + 1), velY(k + 1));

                if velTotal(k) > maxVel
                    maxVel = velTotal(k);
                    tMaxVel = k;
                end

                if accTotal(k) > maxAcc
                    maxAcc = accTotal(k);
                end

                if accTotal(k) < maxDec
                    maxDec = accTotal(k);
                end

                offsetCriterion = ...
                    velTotal(k + 1) < cfg.velocityOffset_deg_s || ...
                    angleDiffDeg(nextSampleDir, mainDir) > cfg.maxDirectionDeviation_deg || ...
                    angleDiffDeg(sampleDir, nextSampleDir) > cfg.maxSampleDirectionChange_deg;

                if offsetCriterion
                    Saccades_right(sac, 2) = k;
                    Saccades_right(sac, 4) = gazeX(k);
                    Saccades_right(sac, 6) = gazeY(k);
                    break;
                end
            end

            Saccades_right(sac, 7)  = tMaxVel;
            Saccades_right(sac, 8)  = maxVel;
            Saccades_right(sac, 9)  = maxAcc;
            Saccades_right(sac, 10) = maxDec;
            Saccades_right(sac, 16) = sacPart;
            Saccades_right(sac, 17) = i;

            % Remove duplicate candidates with the same peak-velocity sample.
            if sac > 1 && Saccades_right(sac, 7) == Saccades_right(sac - 1, 7)
                durationCurrent = Saccades_right(sac, 2)     - Saccades_right(sac, 1);
                durationPrevious = Saccades_right(sac - 1, 2) - Saccades_right(sac - 1, 1);

                if durationCurrent <= durationPrevious
                    Saccades_right(sac, :) = [];
                    sac = sac - 1;
                    sacPart = sacPart - 1;
                else
                    Saccades_right(sac - 1, :) = [];
                    sac = sac - 1;
                    sacPart = sacPart - 1;
                end
            end
        end
    end

    Saccades_right = Saccades_right(Saccades_right(:, 1) ~= 0 & Saccades_right(:, 2) ~= 0, :);
end


function Saccades_right = addSaccadeMetrics(Saccades_right, T, velX, velY, cfg)
% Calculate duration, amplitude, direction, vertical component, and trial ID.

    if isempty(Saccades_right)
        return;
    end

    trialVec = getOptionalColumn(T, ["TRIAL", "Trial"], NaN(height(T), 1));

    for i = 1:size(Saccades_right, 1)

        startSample = round(Saccades_right(i, 1));
        endSample   = round(Saccades_right(i, 2));
        peakSample  = round(Saccades_right(i, 7));

        Saccades_right(i, 11) = (endSample - startSample + 1) / cfg.fs * 1000;

        xDelta_px = abs(Saccades_right(i, 4) - Saccades_right(i, 3));
        yDelta_px = abs(Saccades_right(i, 6) - Saccades_right(i, 5));

        Saccades_right(i, 12) = pixToDegX(xDelta_px, cfg);
        Saccades_right(i, 13) = pixToDegY(yDelta_px, cfg);
        Saccades_right(i, 14) = hypot(Saccades_right(i, 12), Saccades_right(i, 13));

        Saccades_right(i, 15) = directionDeg(velX(peakSample), velY(peakSample));

        if Saccades_right(i, 14) > 0
            Saccades_right(i, 18) = Saccades_right(i, 13) / Saccades_right(i, 14);
        else
            Saccades_right(i, 18) = NaN;
        end

        if all(~isnan(trialVec(startSample:endSample)))
            Saccades_right(i, 19) = mode(trialVec(startSample:endSample));
        else
            Saccades_right(i, 19) = NaN;
        end
    end
end


function Saccades_right_main = selectMainSaccadePerCluster(Saccades_right_long, cfg)
% Select one main saccade per cluster. Saccades separated by less than the
% merge gap are considered part of the same cluster; the saccade with the
% highest peak velocity is kept.

    if isempty(Saccades_right_long)
        Saccades_right_main = zeros(0, 21);
        return;
    end

    gapSamples = ms2samples(cfg.mergeGap_ms, cfg.fs);
    n = size(Saccades_right_long, 1);

    isClusterStart = true(n, 1);
    isClusterStart(2:end) = ...
        Saccades_right_long(2:end, 1) - Saccades_right_long(1:end-1, 2) >= gapSamples;

    clusterID = cumsum(isClusterStart);

    nClusters = max(clusterID);
    Saccades_right_main = zeros(nClusters, 21);

    for c = 1:nClusters
        rows = find(clusterID == c);
        [~, localMaxRow] = max(Saccades_right_long(rows, 8));
        row = rows(localMaxRow);

        Saccades_right_main(c, :) = Saccades_right_long(row, :);
    end

    Saccades_right_main = Saccades_right_main(Saccades_right_main(:, 1) ~= 0, :);
end


function FileThresholdRows = sweepSWJGrid( ...
    filename, fileIdx, durSec, durMin, Saccades_right_main, validGaze, GRID, cfg)
% Run all SWJ criteria combinations for one file.

    FileThresholdRows = table();

    for aMin = GRID.ampMinList
        for aMax = GRID.ampMaxList
            for w = 1:size(GRID.isiWindows_ms, 1)

                isiMin = GRID.isiWindows_ms(w, 1);
                isiMax = GRID.isiWindows_ms(w, 2);

                for dirTol = GRID.dirTolList_deg
                    for simType = GRID.simTypeList
                        for simThr = GRID.simThrList

                            params = struct();
                            params.amp_min  = aMin;
                            params.amp_max  = aMax;
                            params.isi_min  = isiMin;
                            params.isi_max  = isiMax;
                            params.dir_tol  = dirTol;
                            params.sim_type = simType;
                            params.sim_thr  = simThr;

                            [~, details] = detectSWJPairs( ...
                                Saccades_right_main, params, validGaze, cfg);

                            metrics = summariseSWJDetails(details, durMin);

                            newRow = table( ...
                                string(filename), fileIdx, durSec, ...
                                aMin, aMax, isiMin, isiMax, ...
                                dirTol, simType, simThr, ...
                                metrics.SWJ_Count, ...
                                metrics.SWJ_Rate_per_min, ...
                                metrics.MeanPairAmp_deg, ...
                                metrics.Total_Amp, ...
                                metrics.MeanISI_ms, ...
                                metrics.MeanPairVerticalComp, ...
                                'VariableNames', { ...
                                    'File', 'FileIndex', 'Dur_sec', ...
                                    'AmpMin', 'AmpMax', 'ISImin', 'ISImax', ...
                                    'DirTol_deg', 'SimType', 'SimThr', ...
                                    'SWJ_Count', 'SWJ_Rate_per_min', ...
                                    'MeanPairAmp_deg', 'Total_Amp', ...
                                    'MeanISI_ms', 'MeanPairVerticalComp'} ...
                                );

                            FileThresholdRows = [FileThresholdRows; newRow]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end


function [pairs, details] = detectSWJPairs(Saccades, params, validGaze, cfg)
% Detect SWJ pairs from Saccades_right_main.
%
% Pair definition:
%   - both components satisfy the amplitude range
%   - intra-SWJ interval is inside [isi_min, isi_max]
%   - direction difference is approximately 180 degrees
%   - amplitudes are sufficiently similar
%   - no invalid gaze samples between the first component onset and the
%     return component offset
%   - if trial IDs are available, both components must occur in the same trial
%
% Output details columns:
%   1  Saccade1_Index
%   2  Saccade2_Index
%   3  Amplitude_Sac1
%   4  Amplitude_Sac2
%   5  Direction_Sac1
%   6  Direction_Sac2
%   7  PeakVel_Sac1
%   8  PeakVel_Sac2
%   9  Angle_Difference
%   10 InterSaccadic_Interval_ms
%   11 Amplitude_Dissimilarity
%   12 VerticalComp_Sac1
%   13 VerticalComp_Sac2
%   14 PairMeanVerticalComp

    pairs = zeros(0, 2);
    details = zeros(0, 14);

    n = size(Saccades, 1);

    if n < 2
        return;
    end

    used = false(n, 1);
    i = 1;

    while i < n

        if used(i)
            i = i + 1;
            continue;
        end

        pairFound = false;
        maxJ = min(n, i + 1 + cfg.maxInterveningSaccades);

        for j = (i + 1):maxJ

            if used(j)
                continue;
            end

            sac1 = Saccades(i, :);
            sac2 = Saccades(j, :);

            start1 = round(sac1(1));
            end1   = round(sac1(2));
            start2 = round(sac2(1));
            end2   = round(sac2(2));

            if start1 < 1 || end2 > numel(validGaze) || start2 <= end1
                continue;
            end

            if any(~validGaze(start1:end2))
                continue;
            end

            trial1 = sac1(19);
            trial2 = sac2(19);

            if ~isnan(trial1) && ~isnan(trial2) && trial1 ~= trial2
                continue;
            end

            amp1 = sac1(14);
            amp2 = sac2(14);

            ampOK = amp1 >= params.amp_min && amp1 <= params.amp_max && ...
                    amp2 >= params.amp_min && amp2 <= params.amp_max;

            if ~ampOK
                continue;
            end

            isi_ms = (start2 - end1) / cfg.fs * 1000;

            if isi_ms < params.isi_min || isi_ms > params.isi_max
                continue;
            end

            dir1 = sac1(15);
            dir2 = sac2(15);

            angleDiff = angleDiffDeg(dir1, dir2);
            directionOK = abs(angleDiff - 180) <= params.dir_tol;

            if ~directionOK
                continue;
            end

            ampDissim = amplitudeDissimilarity(amp1, amp2, params.sim_type);

            if ampDissim > params.sim_thr
                continue;
            end

            vc1 = sac1(18);
            vc2 = sac2(18);
            pairMeanVC = mean([vc1 vc2], 'omitnan');

            pairs(end + 1, :) = [i j]; %#ok<AGROW>
            details(end + 1, :) = [ ...
                i, j, ...
                amp1, amp2, ...
                dir1, dir2, ...
                sac1(8), sac2(8), ...
                angleDiff, isi_ms, ampDissim, ...
                vc1, vc2, pairMeanVC]; %#ok<AGROW>

            used(i) = true;
            used(j) = true;
            pairFound = true;
            break;
        end

        if pairFound
            i = i + 2;
        else
            i = i + 1;
        end
    end
end


function metrics = summariseSWJDetails(details, durMin)
% Calculate final summary metrics for a set of SWJ pairs.

    metrics = struct();

    metrics.SWJ_Count = size(details, 1);
    metrics.SWJ_Rate_per_min = metrics.SWJ_Count / (durMin + eps);

    if isempty(details)
        metrics.MeanPairAmp_deg = NaN;
        metrics.Total_Amp = NaN;
        metrics.MeanISI_ms = NaN;
        metrics.MeanPairVerticalComp = NaN;
        return;
    end

    pairAmp = mean(details(:, 3:4), 2, 'omitnan');

    metrics.MeanPairAmp_deg = mean(pairAmp, 'omitnan');
    metrics.Total_Amp = sum(pairAmp, 'omitnan');
    metrics.MeanISI_ms = mean(details(:, 10), 'omitnan');
    metrics.MeanPairVerticalComp = mean(details(:, 14), 'omitnan');
end


function outTable = swjDetailsToTable(details)
% Convert event-level SWJ details to a labelled table.

    varNames = { ...
        'Saccade1_Index', 'Saccade2_Index', ...
        'Amplitude_Sac1', 'Amplitude_Sac2', ...
        'Direction_Sac1', 'Direction_Sac2', ...
        'PeakVel_Sac1', 'PeakVel_Sac2', ...
        'Angle_Difference', 'InterSaccadic_Interval_ms', ...
        'Amplitude_Dissimilarity', ...
        'VerticalComp_Sac1', 'VerticalComp_Sac2', ...
        'PairMeanVerticalComp'};

    if isempty(details)
        outTable = array2table(zeros(0, numel(varNames)), 'VariableNames', varNames);
    else
        outTable = array2table(details, 'VariableNames', varNames);
    end
end


function outTable = saccadesMainToTable(S)
% Convert Saccades_right_main matrix to a labelled table for QC/reuse.

    varNames = { ...
        'StartSample', 'EndSample', ...
        'StartX_px', 'EndX_px', ...
        'StartY_px', 'EndY_px', ...
        'PeakVelSample', 'PeakVel_deg_s', ...
        'MaxAcc_deg_s2', 'MaxDec_deg_s2', ...
        'Duration_ms', ...
        'AmpX_deg', 'AmpY_deg', 'AmpTotal_deg', ...
        'Direction_deg', ...
        'SacPart', 'ClusterID', ...
        'VerticalComponent', ...
        'Trial', 'Reserved20', 'Reserved21'};

    if isempty(S)
        outTable = array2table(zeros(0, numel(varNames)), 'VariableNames', varNames);
    else
        outTable = array2table(S, 'VariableNames', varNames);
    end
end


function value = amplitudeDissimilarity(amp1, amp2, simType)
% Return amplitude dissimilarity. Lower values indicate more similar
% amplitudes.

    switch lower(string(simType))
        case "dissim"
            value = abs(amp1 - amp2) / (abs(amp1 + amp2) + eps);

        case "absdiff"
            value = abs(amp1 - amp2);

        case "relmax"
            value = abs(amp1 - amp2) / (max([amp1 amp2]) + eps);

        case "none"
            value = 0;

        otherwise
            error("Unknown amplitude similarity type: %s", simType);
    end
end


function x = getRequiredColumn(T, name)
% Return a required column from a table.

    name = string(name);

    if ~ismember(name, string(T.Properties.VariableNames))
        error("Missing required input column: %s", name);
    end

    x = T.(name);
    x = double(x(:));
end


function x = getOptionalColumn(T, names, defaultValue)
% Return the first available optional column, or a default vector.

    names = string(names);

    for i = 1:numel(names)
        if ismember(names(i), string(T.Properties.VariableNames))
            x = T.(names(i));
            x = double(x(:));
            return;
        end
    end

    x = defaultValue;
end


function maskOut = expandMask(maskIn, nBefore, nAfter)
% Expand true samples in a logical mask by nBefore and nAfter samples.

    maskIn = logical(maskIn(:));
    n = numel(maskIn);
    maskOut = false(n, 1);

    idx = find(maskIn);

    for k = 1:numel(idx)
        lo = max(1, idx(k) - nBefore);
        hi = min(n, idx(k) + nAfter);
        maskOut(lo:hi) = true;
    end
end


function maskOut = mergeShortGaps(maskIn, maxGapSamples)
% Merge two true segments if the gap between them is <= maxGapSamples.

    maskOut = logical(maskIn(:));

    if ~any(maskOut)
        return;
    end

    d = diff([false; maskOut; false]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;

    for k = 1:(numel(ends) - 1)
        gapStart = ends(k) + 1;
        gapEnd   = starts(k + 1) - 1;
        gapLength = gapEnd - gapStart + 1;

        if gapLength <= maxGapSamples
            maskOut(gapStart:gapEnd) = true;
        end
    end
end


function nSamples = ms2samples(ms, fs)
% Convert milliseconds to samples.

    nSamples = round(ms / 1000 * fs);
end


function deg = pixToDegX(deltaPix, cfg)
% Convert horizontal pixel displacement to degrees of visual angle.

    totalAngleDeg = 2 * atan2d(cfg.widthCm / 2, cfg.viewDistCm);
    degPerPix = totalAngleDeg / cfg.widthPix;
    deg = deltaPix * degPerPix;
end


function deg = pixToDegY(deltaPix, cfg)
% Convert vertical pixel displacement to degrees of visual angle.

    totalAngleDeg = 2 * atan2d(cfg.heightCm / 2, cfg.viewDistCm);
    degPerPix = totalAngleDeg / cfg.heightPix;
    deg = deltaPix * degPerPix;
end


function dir = directionDeg(vx, vy)
% Direction of the velocity vector in degrees, wrapped to [0, 360).

    dir = mod(atan2d(vy, vx), 360);

    if isnan(vx) || isnan(vy)
        dir = NaN;
    end
end


function d = angleDiffDeg(a, b)
% Smallest absolute angle difference between two directions in degrees.
% Output range: [0, 180].

    d = abs(mod((a - b) + 180, 360) - 180);
end
