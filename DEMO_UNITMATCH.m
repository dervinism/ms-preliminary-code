%% DEMO UNIT MATCH 


%% READ ME
% UnitMatch contains a script (i.e. ExtractAndSaveAverageWaveform.m) to extract average waveforms from raw SpikeGLX data. 
% If you do not use SpikeGLX, you have multiple options:
% 1. Convert your data to SpikeGLX format (e.g. Use OpenEphys?)
% 2. Use Bombcell frist, the extracted waveforms with the defaults for UnitMatch, raw waveforms will be stored in the KS folder
% 3. Extract the waveforms yourself. 
% Either one of the options above should result in 'KiloSortPaths' containing a subfolder called 'RawWaveforms'. There
% should be a NPY file for every cluster with the dimensions of
% UMparam.spikeWidth X nRecordingChannels X 2 (1 for each half of a
% recording). This should contain the average waveform (recommended of at
% least 500 spikes) for every recording channel for every half of a
% recording for that cluster.


%% User input:
UMparam.DecompressLocal = 0; %if 1, uncompress data first if it's currently compressed
UMparam.RunQualityMetrics = 0;
UMparam.SaveDir = 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\unitmatch'; % Recommended to use end this path with \Probe0\IMRO_1\ if more probes/IMRO tables were used or \AllProbes\AllIMRO\ otherwise
UMparam.KSDir = {'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-25_15-54-22', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-29_13-15-57', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-30_12-27-33', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-01_13-07-16', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-02_12-41-03', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-08_17-24-19', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-09_13-22-13', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-14_15-58-44', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-15_15-17-43', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-21_13-00-21', ...
                 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-23_16-22-54'};  % This is a cell array with a path, in the path there should be a subfolder called 'RawWaveforms'. 
% N.B. if you want to use the functional score evaluation of UnitMatch, 'KSDir' should also contain typical 'Kilosort output', (e.g. spike times etc.)

%% N.B. the following user input can also be automatically extracted and prepared/cleaned up using UMparam = ExtractKilosortData(KiloSortPaths, UMparam) for Kilosorted data of SpikeGLX recorded data (see next section);
UMparam.RawDataPaths = {'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-25_15-54-22\PP01_2020-06-25_15-54-22.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-29_13-15-57\PP01_2020-06-29_13-15-57.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-30_12-27-33\PP01_2020-06-30_12-27-33.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-01_13-07-16\PP01_2020-07-01_13-07-16.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-02_12-41-03\PP01_2020-07-02_12-41-03.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-08_17-24-19\PP01_2020-07-08_17-24-19.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-09_13-22-13\PP01_2020-07-09_13-22-13.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-14_15-58-44\PP01_2020-07-14_15-58-44.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-15_15-17-43\PP01_2020-07-15_15-17-43.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-21_13-00-21\PP01_2020-07-21_13-00-21.dat', ...
                        'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-23_16-22-54\PP01_2020-07-23_16-22-54.dat'};  % This is a cell array with info on where to find the compressed recording (.cbin files OR .bin files)
UMparam.AllDecompPaths = {'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-25_15-54-22\PP01_2020-06-25_15-54-22.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-29_13-15-57\PP01_2020-06-29_13-15-57.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-06-30_12-27-33\PP01_2020-06-30_12-27-33.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-01_13-07-16\PP01_2020-07-01_13-07-16.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-02_12-41-03\PP01_2020-07-02_12-41-03.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-08_17-24-19\PP01_2020-07-08_17-24-19.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-09_13-22-13\PP01_2020-07-09_13-22-13.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-14_15-58-44\PP01_2020-07-14_15-58-44.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-15_15-17-43\PP01_2020-07-15_15-17-43.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-21_13-00-21\PP01_2020-07-21_13-00-21.dat', ...
                          'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-23_16-22-54\PP01_2020-07-23_16-22-54.dat'};  % This is a cell array with info on where to find the decompressed recording (.bin files) --> Necessary when you want UnitMatch to do waveform extraction
UMparam.AllChannelPos = [(1:384)' (1:384)']; % These are coordinates of every recording channel on the probe (e.g. nRecordingChannels x 2)
spikesRepoPath = 'C:\Users\Martynas\Matlab_libraries\EnnyvanBeest_spikes';
loadKSdirPath = 'C:\Users\Martynas\Matlab_libraries\UnitMatch\MATLAB\DataPreparation';
clusinfo = struct; % Note, this can be kilosort input, 
% - clusinfo (this is a struct that contains per unit the following information):
% * cluster_id (e.g. kilosort output clus_id)
% * Good_ID: ones for units that should be included in the analysis
% * RecSesID: Recording Session ID
% * Probe: Which probe (if just 1, ones of numel cluster_id)
% * Depth: depth on probe (optional)
% * Shank: Which shank (optional)
% * Coordinates: Typically 3D Allen Common Coordinate framework coordinates per unit (optional)


% N.B. clusinfo can also be automatically extracted using clusinfo =
% getClusinfo


%% Add paths and subpaths
mfilePath = mfilename('fullpath');
if contains(mfilePath,'LiveEditorEvaluationHelper')
    mfilePath = matlab.desktop.editor.getActiveFilename;
end
Components = strsplit(mfilePath,filesep);
addpath(genpath(fullfile(Components{1:end-1})));
addpath(genpath(spikesRepoPath));
addpath(loadKSdirPath);

%% Optional (for Kilosort + SpikeGLX users) --- see ExampleAnalysisPipelines for more detail!!
UMparam = ExtractKilosortData(UMparam.KSDir, UMparam); % Extract KS data and do some noise removal, optionally decompresses cbin to bin data and uses BOMBCELL quality metric to define good single units
clusinfo = getClusinfo(UMparam.KSDir); % prepare clusinfo struct

%% Load default parameters
UMparam = DefaultParametersUnitMatch(UMparam);

%% UnitMatch algorithm:
[UniqueIDConversion, MatchTable, WaveformInfo, UMparam] = UnitMatch(clusinfo, UMparam);
if UMparam.AssignUniqueID
    AssignUniqueID(UMparam.SaveDir);
end

%%% N.B. From here it is all evaluation, you don't need this to use UnitMatch
%%% results in your further analysis
%% Automatic evaluation:
EvaluatingUnitMatch(UMparam.SaveDir); % Within session cross-validation
QualityMetricsROCs(UMparam.SaveDir); % Only works in combination with BOMBCELL
ComputeFunctionalScores(UMparam.SaveDir) % Only works when having access to Kilosort output (e.g. spike times etc.) 

%% Curation:
if UMparam.MakePlotsOfPairs
    DrawPairsUnitMatch(UMparam.SaveDir);
    if UMparam.GUI
        FigureFlick(UMparam.SaveDir)
        pause
    end
end