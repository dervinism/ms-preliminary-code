% This script performs principal component analyses for ms-preliminary dataset

%% Initialise parameters
% Data setup parameters
params

% Analysis parameters
samplingIntervalCoh = 0.02;
samplingIntervalPCA = 0.002;
convolutionPoints = 25;
frequencyRange = [1 20];
minIntervalLength = 1;
parallelise = true;
nPhaseHistogramBins = 9;
nPCs = 3;
filterType = 'wideband';
fitCircle = true;
fitCircleColour = true;
fitCircleDensity = true;

% Unit quality criteria
firingRateCutoff = 200; % At least that many spikes per hour or at least a single 20-minute window (i.e., firingRateCutoff/3)
refractoryContaminationTh = 1; % Not more than this percentage of spiking ACG refractory period contamination.
                               % Will take the ratio mean(+/- 1 ms relative to 0 lag of ACG) /
                               % mean(1000-1500ms ACG away from zero lag (ACG shoulder)).
coherenceTh = 0.0;

% Other unit selection criteria
oscTh = 0; %oscSignificanceCutoff;

% Script run parameters
loadThetaPhaseData = true;
loadSpikeTimesData = true;
runThetaMazeAnalyses_maze = true;
createThetaMazeFigs_maze = true;
runThetaMazeAnalyses_trials = true;
createThetaMazeFigs_trials = true;
runThetaMazeAnalyses_speed = true;
createThetaMazeFigs_speed = true;
runSleepAnalyses_rem = true;
createThetaMazeFigs_rem = true;
runSleepAnalyses_rem1 = true;
createThetaMazeFigs_rem1 = true;
runSleepAnalyses_rem2 = true;
createThetaMazeFigs_rem2 = true;
runSleepAnalyses_remTheta = true;
createThetaMazeFigs_remTheta = true;
runSleepAnalyses_nrem = true;
createThetaMazeFigs_nrem = true;
runSleepAnalyses_nrem1 = true;
createThetaMazeFigs_nrem1 = true;
runSleepAnalyses_nrem2 = true;
createThetaMazeFigs_nrem2 = true;
runSleepAnalyses_nremTheta = true;
createThetaMazeFigs_nremTheta = true;
runModeratePowerAnalyses = false;
createModeratePowerFigs = false;
runHighPowerAnalyses = false;
createHighPowerFigs = false;
runModerateAmplitudeAnalyses = false;
createModerateAmplitudeFigs = false;
runHighAmplitudeAnalyses = false;
createHighAmplitudeFigs = false;


%% Get population rate theta phase values
if ~exist('thetaPhaseData','var') || loadThetaPhaseData
  [thetaPhaseData, thetaUnwrappedPhaseData, thetaFrequencyData, ...
    thetaAmplitudeData, thetaPowerData] = thetaWrapper(dataFiles, filterType);
end
if ~exist('spikeTimes','var') || loadSpikeTimesData
  [spikeTimes, ~, populationRates, populationRatesTimesBins, ...
    filtPopulationRates, unitIDs] = spikeTimesAggregator(dataFiles, ...
    samplingIntervalPCA, convolutionPoints, frequencyRange);
end


%% Correlate PCA coefficients with various measures
if runThetaMazeAnalyses_maze

  % Define periods of interest: ThetaMaze_AlternativeRunning
  epochsOfInterest = 'ThetaMaze_AlternativeRunning'; %#ok<*UNRCH>
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_maze, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_maze, pcaOutFilt_maze, scoreCorr_maze] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: ThetaMaze_AlternativeRunning
if createThetaMazeFigs_maze
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'ThetaMaze_AlternativeRunning');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_maze, pcaOutFilt_maze, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_maze, pcaOutFilt_maze, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_maze, pcaOutFilt_maze, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_maze, pcaOutFilt_maze, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_maze, pcaOutFilt_maze, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_maze, pcaOutFilt_maze, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_maze, pcaOutFilt_maze, ...
    fullInterpCoherence_maze, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_maze, pcaOutFilt_maze, scoreCorr_maze, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: Behaviour
  [radiiSpeed_rho_maze, radiiSpeed_p_maze, radiiSpeed_rhoFilt_maze, radiiSpeed_pFilt_maze, ...
    radiiAcceleration_rho_maze, radiiAcceleration_p_maze, radiiAcceleration_rhoFilt_maze, radiiAcceleration_pFilt_maze, ...
    radiiPGDI_rho_maze, radiiPGDI_p_maze, radiiPGDI_rhoFilt_maze, radiiPGDI_pFilt_maze, ...
    radiiPGDISmoothed_rho_maze, radiiPGDISmoothed_p_maze, radiiPGDISmoothed_rhoFilt_maze, radiiPGDISmoothed_pFilt_maze, ...
    radiiAmplitude_rho_maze, radiiAmplitude_p_maze, radiiAmplitude_rhoFilt_maze, radiiAmplitude_pFilt_maze, ...
    radiiPower_rho_maze, radiiPower_p_maze, radiiPower_rhoFilt_maze, radiiPower_pFilt_maze, ...
    radiiFrequency_rho_maze, radiiFrequency_p_maze, radiiFrequency_rhoFilt_maze, radiiFrequency_pFilt_maze, ...
    speedAmplitude_rho_maze, speedAmplitude_p_maze, speedFrequency_rho_maze, speedFrequency_p_maze] = ...
    cycleBehavCorrWrapper(dataFiles, pcaOut_maze, pcaOutFilt_maze, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, drawFigs=true, saveFigs=true);
end


if runThetaMazeAnalyses_trials

  % Define periods of interest: ThetaMaze_AlternativeRunning_trials
  epochsOfInterest = 'ThetaMaze_AlternativeRunning';
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = true;
  onlyHighSpeed = false;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_trials, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_trials, pcaOutFilt_trials, scoreCorr_trials] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: ThetaMaze_AlternativeRunning_trials
if createThetaMazeFigs_trials
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'ThetaMaze_AlternativeRunning_trials');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_trials, pcaOutFilt_trials, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_trials, pcaOutFilt_trials, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_trials, pcaOutFilt_trials, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_trials, pcaOutFilt_trials, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_trials, pcaOutFilt_trials, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_trials, pcaOutFilt_trials, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_trials, pcaOutFilt_trials, ...
    fullInterpCoherence_trials, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_trials, pcaOutFilt_trials, scoreCorr_trials, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: Behaviour
  [radiiSpeed_rho_trials, radiiSpeed_p_trials, radiiSpeed_rhoFilt_trials, radiiSpeed_pFilt_trials, ...
    radiiAcceleration_rho_trials, radiiAcceleration_p_trials, radiiAcceleration_rhoFilt_trials, radiiAcceleration_pFilt_trials, ...
    radiiPGDI_rho_trials, radiiPGDI_p_trials, radiiPGDI_rhoFilt_trials, radiiPGDI_pFilt_trials, ...
    radiiPGDISmoothed_rho_trials, radiiPGDISmoothed_p_trials, radiiPGDISmoothed_rhoFilt_trials, radiiPGDISmoothed_pFilt_trials, ...
    radiiAmplitude_rho_trials, radiiAmplitude_p_trials, radiiAmplitude_rhoFilt_trials, radiiAmplitude_pFilt_trials, ...
    radiiPower_rho_trials, radiiPower_p_trials, radiiPower_rhoFilt_trials, radiiPower_pFilt_trials, ...
    radiiFrequency_rho_trials, radiiFrequency_p_trials, radiiFrequency_rhoFilt_trials, radiiFrequency_pFilt_trials, ...
    speedAmplitude_rho_trials, speedAmplitude_p_trials, speedFrequency_rho_trials, speedFrequency_p_trials] = ...
    cycleBehavCorrWrapper(dataFiles, pcaOut_trials, pcaOutFilt_trials, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, drawFigs=true, saveFigs=true);
end


if runThetaMazeAnalyses_speed
  % Define periods of interest: ThetaMaze_AlternativeRunning_trials_highSpeed
  epochsOfInterest = 'ThetaMaze_AlternativeRunning';
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = true;
  onlyHighSpeed = true;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_speed, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_speed, pcaOutFilt_speed, scoreCorr_speed] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: ThetaMaze_AlternativeRunning_speed
if createThetaMazeFigs_speed
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'ThetaMaze_AlternativeRunning_speed');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_speed, pcaOutFilt_speed, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_speed, pcaOutFilt_speed, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_speed, pcaOutFilt_speed, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_speed, pcaOutFilt_speed, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_speed, pcaOutFilt_speed, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_speed, pcaOutFilt_speed, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_speed, pcaOutFilt_speed, ...
    fullInterpCoherence_speed, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_speed, pcaOutFilt_speed, scoreCorr_speed, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: Behaviour
  [radiiSpeed_rho_speed, radiiSpeed_p_speed, radiiSpeed_rhoFilt_speed, radiiSpeed_pFilt_speed, ...
    radiiAcceleration_rho_speed, radiiAcceleration_p_speed, radiiAcceleration_rhoFilt_speed, radiiAcceleration_pFilt_speed, ...
    radiiPGDI_rho_speed, radiiPGDI_p_speed, radiiPGDI_rhoFilt_speed, radiiPGDI_pFilt_speed, ...
    radiiPGDISmoothed_rho_speed, radiiPGDISmoothed_p_speed, radiiPGDISmoothed_rhoFilt_speed, radiiPGDISmoothed_pFilt_speed, ...
    radiiAmplitude_rho_speed, radiiAmplitude_p_speed, radiiAmplitude_rhoFilt_speed, radiiAmplitude_pFilt_speed, ...
    radiiPower_rho_speed, radiiPower_p_speed, radiiPower_rhoFilt_speed, radiiPower_pFilt_speed, ...
    radiiFrequency_rho_speed, radiiFrequency_p_speed, radiiFrequency_rhoFilt_speed, radiiFrequency_pFilt_speed, ...
    speedAmplitude_rho_speed, speedAmplitude_p_speed, speedFrequency_rho_speed, speedFrequency_p_speed] = ...
    cycleBehavCorrWrapper(dataFiles, pcaOut_speed, pcaOutFilt_speed, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_rem

  % Define periods of interest: REM
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_rem, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_rem, pcaOutFilt_rem, scoreCorr_rem] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM
if createThetaMazeFigs_rem
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'rem');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_rem, pcaOutFilt_rem, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_rem, pcaOutFilt_rem, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_rem, pcaOutFilt_rem, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_rem, pcaOutFilt_rem, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_rem, pcaOutFilt_rem, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_rem, pcaOutFilt_rem, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_rem, pcaOutFilt_rem, ...
    fullInterpCoherence_rem, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_rem, pcaOutFilt_rem, scoreCorr_rem, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_rem, radiiPGDI_p_rem, radiiPGDI_rhoFilt_rem, radiiPGDI_pFilt_rem, ...
    radiiPGDISmoothed_rho_rem, radiiPGDISmoothed_p_rem, radiiPGDISmoothed_rhoFilt_rem, radiiPGDISmoothed_pFilt_rem, ...
    radiiAmplitude_rho_rem, radiiAmplitude_p_rem, radiiAmplitude_rhoFilt_rem, radiiAmplitude_pFilt_rem, ...
    radiiPower_rho_rem, radiiPower_p_rem, radiiPower_rhoFilt_rem, radiiPower_pFilt_rem, ...
    radiiFrequency_rho_rem, radiiFrequency_p_rem, radiiFrequency_rhoFilt_rem, radiiFrequency_pFilt_rem] = ...
    cycleCorrWrapper(dataFiles, pcaOut_rem, pcaOutFilt_rem, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_rem1

  % Define periods of interest: REM prior to behaviour
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem1';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_rem1, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_rem1, pcaOutFilt_rem1, scoreCorr_rem1] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM prior to behaviour
if createThetaMazeFigs_rem1
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'rem1');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, ...
    fullInterpCoherence_rem1, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, scoreCorr_rem1, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_rem1, radiiPGDI_p_rem1, radiiPGDI_rhoFilt_rem1, radiiPGDI_pFilt_rem1, ...
    radiiPGDISmoothed_rho_rem1, radiiPGDISmoothed_p_rem1, radiiPGDISmoothed_rhoFilt_rem1, radiiPGDISmoothed_pFilt_rem1, ...
    radiiAmplitude_rho_rem1, radiiAmplitude_p_rem1, radiiAmplitude_rhoFilt_rem1, radiiAmplitude_pFilt_rem1, ...
    radiiPower_rho_rem1, radiiPower_p_rem1, radiiPower_rhoFilt_rem1, radiiPower_pFilt_rem1, ...
    radiiFrequency_rho_rem1, radiiFrequency_p_rem1, radiiFrequency_rhoFilt_rem1, radiiFrequency_pFilt_rem1] = ...
    cycleCorrWrapper(dataFiles, pcaOut_rem1, pcaOutFilt_rem1, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_rem2

  % Define periods of interest: REM after to behaviour
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem2';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_rem2, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_rem2, pcaOutFilt_rem2, scoreCorr_rem2] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM after to behaviour
if createThetaMazeFigs_rem2
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'rem2');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, ...
    fullInterpCoherence_rem2, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, scoreCorr_rem2, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_rem2, radiiPGDI_p_rem2, radiiPGDI_rhoFilt_rem2, radiiPGDI_pFilt_rem2, ...
    radiiPGDISmoothed_rho_rem2, radiiPGDISmoothed_p_rem2, radiiPGDISmoothed_rhoFilt_rem2, radiiPGDISmoothed_pFilt_rem2, ...
    radiiAmplitude_rho_rem2, radiiAmplitude_p_rem2, radiiAmplitude_rhoFilt_rem2, radiiAmplitude_pFilt_rem2, ...
    radiiPower_rho_rem2, radiiPower_p_rem2, radiiPower_rhoFilt_rem2, radiiPower_pFilt_rem2, ...
    radiiFrequency_rho_rem2, radiiFrequency_p_rem2, radiiFrequency_rhoFilt_rem2, radiiFrequency_pFilt_rem2] = ...
    cycleCorrWrapper(dataFiles, pcaOut_rem2, pcaOutFilt_rem2, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_remTheta

  % Define periods of interest: REM with high theta/delta ratio
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = 'high';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_remTheta, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_remTheta, pcaOutFilt_remTheta, scoreCorr_remTheta] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM
if createThetaMazeFigs_remTheta
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'remTheta');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, ...
    fullInterpCoherence_remTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, scoreCorr_remTheta, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_remTheta, radiiPGDI_p_remTheta, radiiPGDI_rhoFilt_remTheta, radiiPGDI_pFilt_remTheta, ...
    radiiPGDISmoothed_rho_remTheta, radiiPGDISmoothed_p_remTheta, radiiPGDISmoothed_rhoFilt_remTheta, radiiPGDISmoothed_pFilt_remTheta, ...
    radiiAmplitude_rho_remTheta, radiiAmplitude_p_remTheta, radiiAmplitude_rhoFilt_remTheta, radiiAmplitude_pFilt_remTheta, ...
    radiiPower_rho_remTheta, radiiPower_p_remTheta, radiiPower_rhoFilt_remTheta, radiiPower_pFilt_remTheta, ...
    radiiFrequency_rho_remTheta, radiiFrequency_p_remTheta, radiiFrequency_rhoFilt_remTheta, radiiFrequency_pFilt_remTheta] = ...
    cycleCorrWrapper(dataFiles, pcaOut_remTheta, pcaOutFilt_remTheta, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_nrem

  % Define periods of interest: NREM
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'nrem';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_nrem, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_nrem, pcaOutFilt_nrem, scoreCorr_nrem] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: NREM
if createThetaMazeFigs_nrem
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'nrem');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, ...
    fullInterpCoherence_nrem, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, scoreCorr_nrem, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_nrem, radiiPGDI_p_nrem, radiiPGDI_rhoFilt_nrem, radiiPGDI_pFilt_nrem, ...
    radiiPGDISmoothed_rho_nrem, radiiPGDISmoothed_p_nrem, radiiPGDISmoothed_rhoFilt_nrem, radiiPGDISmoothed_pFilt_nrem, ...
    radiiAmplitude_rho_nrem, radiiAmplitude_p_nrem, radiiAmplitude_rhoFilt_nrem, radiiAmplitude_pFilt_nrem, ...
    radiiPower_rho_nrem, radiiPower_p_nrem, radiiPower_rhoFilt_nrem, radiiPower_pFilt_nrem, ...
    radiiFrequency_rho_nrem, radiiFrequency_p_nrem, radiiFrequency_rhoFilt_nrem, radiiFrequency_pFilt_nrem] = ...
    cycleCorrWrapper(dataFiles, pcaOut_nrem, pcaOutFilt_nrem, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_nrem1

  % Define periods of interest: REM prior to behaviour
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'nrem1';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_nrem1, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_nrem1, pcaOutFilt_nrem1, scoreCorr_nrem1] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM prior to behaviour
if createThetaMazeFigs_nrem1
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'nrem1');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, ...
    fullInterpCoherence_nrem1, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, scoreCorr_nrem1, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_nrem1, radiiPGDI_p_nrem1, radiiPGDI_rhoFilt_nrem1, radiiPGDI_pFilt_nrem1, ...
    radiiPGDISmoothed_rho_nrem1, radiiPGDISmoothed_p_nrem1, radiiPGDISmoothed_rhoFilt_nrem1, radiiPGDISmoothed_pFilt_nrem1, ...
    radiiAmplitude_rho_nrem1, radiiAmplitude_p_nrem1, radiiAmplitude_rhoFilt_nrem1, radiiAmplitude_pFilt_nrem1, ...
    radiiPower_rho_nrem1, radiiPower_p_nrem1, radiiPower_rhoFilt_nrem1, radiiPower_pFilt_nrem1, ...
    radiiFrequency_rho_nrem1, radiiFrequency_p_nrem1, radiiFrequency_rhoFilt_nrem1, radiiFrequency_pFilt_nrem1] = ...
    cycleCorrWrapper(dataFiles, pcaOut_nrem1, pcaOutFilt_nrem1, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_nrem2

  % Define periods of interest: REM after to behaviour
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'nrem2';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_nrem2, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_nrem2, pcaOutFilt_nrem2, scoreCorr_nrem2] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM after to behaviour
if createThetaMazeFigs_nrem2
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'nrem2');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, ...
    fullInterpCoherence_nrem2, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, scoreCorr_nrem2, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_nrem2, radiiPGDI_p_nrem2, radiiPGDI_rhoFilt_nrem2, radiiPGDI_pFilt_nrem2, ...
    radiiPGDISmoothed_rho_nrem2, radiiPGDISmoothed_p_nrem2, radiiPGDISmoothed_rhoFilt_nrem2, radiiPGDISmoothed_pFilt_nrem2, ...
    radiiAmplitude_rho_nrem2, radiiAmplitude_p_nrem2, radiiAmplitude_rhoFilt_nrem2, radiiAmplitude_pFilt_nrem2, ...
    radiiPower_rho_nrem2, radiiPower_p_nrem2, radiiPower_rhoFilt_nrem2, radiiPower_pFilt_nrem2, ...
    radiiFrequency_rho_nrem2, radiiFrequency_p_nrem2, radiiFrequency_rhoFilt_nrem2, radiiFrequency_pFilt_nrem2] = ...
    cycleCorrWrapper(dataFiles, pcaOut_nrem2, pcaOutFilt_nrem2, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runSleepAnalyses_nremTheta

  % Define periods of interest: NREM with high theta/delta ratio
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  thetaAmplitude = '';
  theta2deltaRatio = 'high';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'nrem';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_nremTheta, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_nremTheta, pcaOutFilt_nremTheta, scoreCorr_nremTheta] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: REM
if createThetaMazeFigs_nremTheta
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'nremTheta');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, ...
    fullInterpCoherence_nremTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);

  % Correlate PCs
  pcCorrFigWrapper(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, scoreCorr_nremTheta, ...
    nPCs, pcaFigFolderLocal)

  % Correlate theta cycle-averaged measures: No behaviour
  [radiiPGDI_rho_nremTheta, radiiPGDI_p_nremTheta, radiiPGDI_rhoFilt_nremTheta, radiiPGDI_pFilt_nremTheta, ...
    radiiPGDISmoothed_rho_nremTheta, radiiPGDISmoothed_p_nremTheta, radiiPGDISmoothed_rhoFilt_nremTheta, radiiPGDISmoothed_pFilt_nremTheta, ...
    radiiAmplitude_rho_nremTheta, radiiAmplitude_p_nremTheta, radiiAmplitude_rhoFilt_nremTheta, radiiAmplitude_pFilt_nremTheta, ...
    radiiPower_rho_nremTheta, radiiPower_p_nremTheta, radiiPower_rhoFilt_nremTheta, radiiPower_pFilt_nremTheta, ...
    radiiFrequency_rho_nremTheta, radiiFrequency_p_nremTheta, radiiFrequency_rhoFilt_nremTheta, radiiFrequency_pFilt_nremTheta] = ...
    cycleCorrWrapper(dataFiles, pcaOut_nremTheta, pcaOutFilt_nremTheta, ...
    thetaUnwrappedPhaseData, thetaFrequencyData, thetaAmplitudeData, ...
    thetaPowerData, intervals, nPCs, figFolder=pcaFigFolderLocal, ...
    drawFigs=true, saveFigs=true);
end


if runModeratePowerAnalyses

  % Define periods of interest: Moderate theta power
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = 'moderate';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_modTheta, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_modTheta, pcaOutFilt_modTheta, scoreCorr_modTheta] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: Moderate theta power
if createModeratePowerFigs
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'moderateThetaPower');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, ...
    fullInterpCoherence_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);
end


if runHighPowerAnalyses

  % Define periods of interest: High theta power
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = 'high';
  thetaAmplitude = '';
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_hiTheta, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_hiTheta, pcaOutFilt_hiTheta, scoreCorr_hiTheta] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: High theta power
if createHighPowerFigs
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'highThetaPower');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, ...
    fullInterpCoherence_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);
end


if runModerateAmplitudeAnalyses

  % Define periods of interest: Moderate theta amplitude
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  if strcmpi(filterType, 'narrowband')
    thetaAmplitude = 'moderate';
  elseif strcmpi(filterType, 'wideband')
    thetaAmplitude = 'widebandModerate';
  end
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_modTheta, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_modTheta, pcaOutFilt_modTheta, scoreCorr_modTheta] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: Moderate theta amplitude
if createModerateAmplitudeFigs
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'moderateThetaAmplitude');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_modTheta, pcaOutFilt_modTheta, ...
    fullInterpCoherence_modTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);
end


if runHighAmplitudeAnalyses

  % Define periods of interest: High theta amplitude
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunning','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  if strcmpi(filterType, 'narrowband')
    thetaAmplitude = 'high';
  elseif strcmpi(filterType, 'wideband')
    thetaAmplitude = 'widebandHigh';
  end
  theta2deltaRatio = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = '';

  % Find the time intervals of interest
  intervals = intervalWrapper(dataFiles, epochsOfInterest, thetaPower, ...
    thetaAmplitude, theta2deltaRatio, minIntervalLength, onlyTrials, ...
    onlyHighSpeed, sleepState, true);

  % Carry out coherence analyses
  [~, fullInterpCoherence_hiTheta, includeUnits] = ...
    coherenceWrapper(dataFiles, intervals, samplingIntervalCoh, ...
    frequencyRange, parallelise, firingRateCutoff, ...
    refractoryContaminationTh, oscTh, coherenceTh);

  % Carry out PCAs for all sessions
  [pcaOut_hiTheta, pcaOutFilt_hiTheta, scoreCorr_hiTheta] = pcaWrapper( ...
    spikeTimes, intervals, includeUnits, samplingInterval=samplingIntervalPCA, ...
    freqRange=frequencyRange, convPoints=convolutionPoints, normalise=true);
end

% Produce various graphs: High theta amplitude
if createHighAmplitudeFigs
  pcaFigFolderLocal = fullfile(pcaFigFolder, 'highThetaAmplitude');

  % Cumulative explained variance bar graphs
  barExplained(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, pcaFigFolderLocal);

  % Visualisation of combined 2 PCs
  combine2PCs(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, intervals, pcaFigFolderLocal, fitCircle);

  % Visualisation of combined 2 PCs in colour
  combine2PCsColour(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, pcaFigFolderLocal, ...
    fitCircleColour);

  % 2PC density map
  density2PCs(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, ...
    pcaFigFolderLocal, fitCircleDensity);

  % Visualisation of combined 3 PCs
  combine3PCs(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, pcaFigFolderLocal);

  % Visualisation of combined 3 PCs in colour
  combine3PCsColour(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, thetaPhaseData, ...
    intervals, pcaFigFolderLocal);

  % Unit phase distributions separated based on loadings for individual PCs
  splitPhaseDistros(dataFiles, pcaOut_hiTheta, pcaOutFilt_hiTheta, ...
    fullInterpCoherence_hiTheta, filtPopulationRates, ...
    populationRatesTimesBins, thetaFrequencyData, intervals, ...
    nPhaseHistogramBins, pcaFigFolderLocal);
end




%% Local functions
function unitPhases = getUnitThetaPhases(dataFiles, filterType) %#ok<*DEFNU> 
% unitPhases = getUnitThetaPhases(dataFiles, filterType)
%
% Function estimates unit theta phases with respect to the population
% firing rate for all animals and sessions in the ms-preliminary dataset.
%
% Args:
%   dataFiles
%   filterType
%
% Returns:
%   unitPhases
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  filterType
end

unitPhases = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})

    % Load unit spiking data
    spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
    if ~exist(spikesFile, 'file')
      continue
    end
    load(spikesFile); %#ok<*LOAD>
    
    % Load theta phase data
    if strcmpi(filterType, 'narrowband')
      thetaPhaseFile = strrep(dataFiles{animal}{session}, '*', 'instantThetaPhase.timeseriesCollection');
    elseif strcmpi(filterType, 'wideband')
      thetaPhaseFile = strrep(dataFiles{animal}{session}, '*', 'widebandInstantPhase.timeseries');
    end
    try
      load(thetaPhaseFile);
    catch
      if strcmpi(filterType, 'narrowband')
        thetaPhaseFile = strrep(dataFiles{animal}{session}, '*', 'instantThetaPhase.timeseriesCollection');
        load(thetaPhaseFile);
        instantThetaPhase.data = instantThetaPhase.data(:,1);
      elseif strcmpi(filterType, 'wideband')
        thetaPhaseFile = strrep(dataFiles{animal}{session}, '*', 'widebandInstantPhase.timeseriesCollection');
        load(thetaPhaseFile);
        widebandInstantPhase.data = widebandInstantPhase.data(:,1);
      end
    end
    if strcmpi(filterType, 'narrowband')
      samplingInterval = instantThetaPhase.timestamps(2) - instantThetaPhase.timestamps(1);
    elseif strcmpi(filterType, 'wideband')
      samplingInterval = widebandInstantPhase.timestamps(2) - widebandInstantPhase.timestamps(1);
    end

    % Get spike theta phases
    unitPhases{animal}{session} = cell(spikes.numcells,1);
    for unit = 1:spikes.numcells
      spikeIndices = round(spikes.times{unit}./samplingInterval);
      if strcmpi(filterType, 'narrowband')
        unitPhases{animal}{session}{unit} = instantThetaPhase.data(spikeIndices);
      elseif strcmpi(filterType, 'wideband')
        unitPhases{animal}{session}{unit} = widebandInstantPhase.data(spikeIndices);
      end
    end
  end
end
end


function [thetaPhaseData, thetaUnwrappedPhaseData, thetaFrequencyData, ...
  thetaAmplitudeData, thetaPowerData] = thetaWrapper(dataFiles, filterType)
% [thetaPhaseData, thetaUnwrappedPhaseData, thetaFrequencyData, ...
%   thetaAmplitudeData, thetaPowerData] = thetaPhaseWrapper(dataFiles, type)
%
% Local wrapper function.
%
% Args:
%   dataFiles
%   type
%
% Returns:
%   thetaPhaseData
%   thetaUnwrappedPhaseData
%   thetaFrequencyData
%   thetaAmplitudeData
%   thetaPowerData
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

thetaPhaseData = cell(numel(dataFiles),1);
thetaUnwrappedPhaseData = cell(numel(dataFiles),1);
thetaFrequencyData = cell(numel(dataFiles),1);
thetaAmplitudeData = cell(numel(dataFiles),1);
thetaPowerData = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if strcmpi(filterType, 'narrowband')
      thetaPhaseFile = strrep(dataFiles{animal}{session}, '*', 'instantThetaPhase.timeseriesCollection');
      thetaFrequencyFile = strrep(dataFiles{animal}{session}, '*', 'instantThetaFrequency.timeseries');
    elseif strcmpi(filterType, 'wideband')
      thetaPhaseFile = strrep(dataFiles{animal}{session}, '*', 'widebandInstantPhase.timeseries');
      thetaFrequencyFile = strrep(dataFiles{animal}{session}, '*', 'widebandInstantFrequency.timeseries');
    end
    thetaAmplitudeFile = strrep(dataFiles{animal}{session}, '*', 'thetaAmplitude.timeseries');
    thetaPowerFile = strrep(dataFiles{animal}{session}, '*', 'thetaPower.timeseriesCollection');
    if exist(thetaPhaseFile,'file')
      load(thetaPhaseFile); %#ok<*LOAD>
      if strcmpi(filterType, 'narrowband')
        thetaPhaseData{animal}{session}.data = instantThetaPhase.data(:,1)';
        thetaPhaseData{animal}{session}.timestamps = instantThetaPhase.timestamps';
        thetaUnwrappedPhaseData{animal}{session}.data = instantThetaPhase.data(:,2)';
        thetaUnwrappedPhaseData{animal}{session}.timestamps = instantThetaPhase.timestamps';
      elseif strcmpi(filterType, 'wideband')
        thetaPhaseData{animal}{session}.data = widebandInstantPhase.data';
        thetaPhaseData{animal}{session}.timestamps = widebandInstantPhase.timestamps';
        thetaUnwrappedPhaseData{animal}{session}.data = unwrap(widebandInstantPhase.data');
        thetaUnwrappedPhaseData{animal}{session}.timestamps = widebandInstantPhase.timestamps';
      end
      load(thetaFrequencyFile);
      if strcmpi(filterType, 'narrowband')
        thetaFrequencyData{animal}{session}.data = instantThetaFrequency.data(:,1)';
        thetaFrequencyData{animal}{session}.timestamps = instantThetaFrequency.timestamps';
      elseif strcmpi(filterType, 'wideband')
        thetaFrequencyData{animal}{session}.data = widebandInstantFrequency.data';
        thetaFrequencyData{animal}{session}.timestamps = widebandInstantFrequency.timestamps';
      end
      load(thetaAmplitudeFile);
      thetaAmplitudeData{animal}{session}.data = thetaAmplitude.data(:,1)';
      thetaAmplitudeData{animal}{session}.timestamps = thetaAmplitude.timestamps';
      load(thetaPowerFile);
      thetaPowerData{animal}{session}.data = thetaPower.data(:,2)';
      thetaPowerData{animal}{session}.timestamps = thetaPower.timestamps';
    else
      thetaPhaseData{animal}{session} = [];
      warning(['Theta phase file ' thetaPhaseFile ' does not exist. Continuing...']);
    end
  end
end
end


function intervals = intervalWrapper(dataFiles, epochsOfInterest, ...
  thetaPower, thetaAmplitude, theta2deltaRatio, minIntervalLength, ...
  onlyTrials, onlyHighSpeed, sleepState, excludeNoise)
% intervals = intervals = intervalWrapper(dataFiles, epochsOfInterest, ...
%   thetaPower, thetaAmplitude, theta2deltaRatio, minIntervalLength, ...
%   onlyTrials, onlyHighSpeed, sleepState, excludeNoise)
%
% Local wrapper function.
%
% Args:
%   dataFiles
%   epochsOfInterest
%   thetaPower
%   thetaAmplitude
%   theta2deltaRatio
%   minIntervalLength
%   onlyTrials
%   onlyHighSpeed
%   sleepState
%   excludeNoise
%
% Returns:
%   intervals
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  epochsOfInterest
  thetaPower
  thetaAmplitude
  theta2deltaRatio
  minIntervalLength
  onlyTrials
  onlyHighSpeed
  sleepState
  excludeNoise
end

% Cycle through animal IDs and sessions
intervals = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    intervals{animal}{session} = getTimeIntervals( ...
      dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
      thetaAmplitude=thetaAmplitude, theta2deltaRatio=theta2deltaRatio, ...
      minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
      onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, ...
      excludeNoise=excludeNoise);
  end
end
end


function [pcaOut, pcaOutFilt, scoreCorr] = pcaWrapper(spikeTimes, ...
  intervals, includeUnits, options)
% [pcaOut, pcaOutFiltscoreCorr] = pcaWrapper(spikeTimes, intervals, ...
%   includeUnits, <options>)
%
% Args:
%   dataFiles
%   intervals
%   includeUnits
%   <samplingInterval>
%   <freqRange>
%   <convPoints>
%   <normalise>
%
% Returns:
%   pcaOut
%   pcaOutFilt
%   scoreCorr
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  spikeTimes
  intervals
  includeUnits
  options.samplingInterval = 0.005
  options.freqRange = [4 12]
  options.convPoints = 5
  options.normalise = false;
end

% Create the output containers
pcaOut = cell(numel(spikeTimes),1);
pcaOutFilt = cell(numel(spikeTimes),1);
scoreCorr = cell(numel(spikeTimes),1);
for animal = 1:numel(spikeTimes)
  for session = 1:numel(spikeTimes{animal})
    if ~isempty(intervals{animal}{session})
      % Perform the PCA
      [pcaOut{animal}{session}, pcaOutFilt{animal}{session}, ...
        scoreCorr{animal}{session}] = pcaConvSpikes( ...
        spikeTimes{animal}{session}, intervals=intervals{animal}{session}, ...
        includeUnits=includeUnits{animal}{session}, ...
        samplingInterval=options.samplingInterval, ...
        convPoints=options.convPoints, freqRange=options.freqRange);
    else
      pcaOut{animal}{session} = [];
      pcaOutFilt{animal}{session} = [];
      scoreCorr{animal}{session} = [];
    end
  end
end
end


function [rho, p, rhoFilt, pFilt, rho_interp, p_interp, rhoFilt_interp, ...
  pFilt_interp] = corrWrapper(fullCoherence, fullInterpCoherence, ...
  pcaOut, pcaOutFilt) %#ok<*DEFNU>
% [rho, p, rhoFilt, pFilt, rho_interp, p_interp, rhoFilt_interp, ...
%   pFilt_interp] = corrWrapper(fullCoherence, fullInterpCoherence, ...
%   pcaOut, pcaOutFilt)
%
% Local wrapper function.
%
% Args:
%   fullCoherence
%   fullInterpCoherence
%   pcaOut
%   pcaOutFilt
%
% Returns:
%   rho
%   p
%   rhoFilt
%   pFilt
%   rho_interp
%   p_interp
%   rhoFilt_interp
%   pFilt_interp
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  fullCoherence
  fullInterpCoherence
  pcaOut
  pcaOutFilt
end

% Create the output containers
rho = cell(numel(fullCoherence),1);
p = cell(numel(fullCoherence),1);
rhoFilt = cell(numel(fullCoherence),1);
pFilt = cell(numel(fullCoherence),1);
rho_interp = cell(numel(fullCoherence),1);
p_interp = cell(numel(fullCoherence),1);
rhoFilt_interp = cell(numel(fullCoherence),1);
pFilt_interp = cell(numel(fullCoherence),1);
for animal = 1:numel(fullCoherence)
  for session = 1:numel(fullCoherence{animal})
    if ~isempty(fullCoherence{animal}{session})

      % Regular coherence
      nFreq = size(fullCoherence{animal}{session}.rateAdjustedCoherence,2);
      nComponents = size(pcaOut{animal}{session}.score,2);
      rho{animal}{session} = zeros(nComponents, nFreq);
      p{animal}{session} = zeros(nComponents, nFreq);
      for f = 1:nFreq
        for comp = 1:nComponents
          [rho{animal}{session}(comp,f), p{animal}{session}(comp,f)] = corr( ...
            fullCoherence{animal}{session}.rateAdjustedCoherence(fullCoherence{animal}{session}.validUnits,f), ...
            pcaOut{animal}{session}.coeff(:,comp), 'type','Spearman');
          [rhoFilt{animal}{session}(comp,f), pFilt{animal}{session}(comp,f)] = corr( ...
            fullCoherence{animal}{session}.rateAdjustedCoherence(fullCoherence{animal}{session}.validUnits,f), ...
            pcaOutFilt{animal}{session}.coeff(:,comp), 'type','Spearman');
        end
      end

      % Interpolated coherence
      nFreq = size(fullInterpCoherence{animal}{session}.rateAdjustedCoherence,2);
      nComponents = size(pcaOut{animal}{session}.score,2);
      rho_interp{animal}{session} = zeros(nComponents, nFreq);
      p_interp{animal}{session} = zeros(nComponents, nFreq);
      for f = 1:nFreq
        for comp = 1:nComponents
          [rho_interp{animal}{session}(comp,f), p_interp{animal}{session}(comp,f)] = corr( ...
            fullInterpCoherence{animal}{session}.rateAdjustedCoherence(fullInterpCoherence{animal}{session}.validUnits,f), ...
            pcaOut{animal}{session}.coeff(:,comp), 'type','Spearman');
          [rhoFilt_interp{animal}{session}(comp,f), pFilt_interp{animal}{session}(comp,f)] = corr( ...
            fullInterpCoherence{animal}{session}.rateAdjustedCoherence(fullInterpCoherence{animal}{session}.validUnits,f), ...
            pcaOutFilt{animal}{session}.coeff(:,comp), 'type','Spearman');
        end
      end

    else
      rho{animal}{session} = [];
      p{animal}{session} = [];
      rhoFilt{animal}{session} = [];
      pFilt{animal}{session} = [];
      rho_interp{animal}{session} = [];
      p_interp{animal}{session} = [];
      rhoFilt_interp{animal}{session} = [];
      pFilt_interp{animal}{session} = [];
    end
  end
end
end


function [radiiSpeed_rho, radiiSpeed_p, radiiSpeed_rhoFilt, radiiSpeed_pFilt, ...
  radiiAcceleration_rho, radiiAcceleration_p, radiiAcceleration_rhoFilt, radiiAcceleration_pFilt, ...
  radiiPGDI_rho, radiiPGDI_p, radiiPGDI_rhoFilt, radiiPGDI_pFilt, ...
  radiiPGDISmoothed_rho, radiiPGDISmoothed_p, radiiPGDISmoothed_rhoFilt, radiiPGDISmoothed_pFilt, ...
  radiiAmplitude_rho, radiiAmplitude_p, radiiAmplitude_rhoFilt, radiiAmplitude_pFilt, ...
  radiiPower_rho, radiiPower_p, radiiPower_rhoFilt, radiiPower_pFilt, ...
  radiiFrequency_rho, radiiFrequency_p, radiiFrequency_rhoFilt, radiiFrequency_pFilt, ...
  speedAmplitude_rho, speedAmplitude_p, speedFrequency_rho, speedFrequency_p] = ...
  cycleBehavCorrWrapper(dataFiles, pcaOut, pcaOutFilt, thetaUnwrappedPhaseData, ...
  thetaFrequencyData, thetaAmplitudeData, thetaPowerData, intervals, nPCs, options)
% [radiiSpeed_rho, radiiSpeed_p, radiiSpeed_rhoFilt, radiiSpeed_pFilt, ...
%   radiiAcceleration_rho, radiiAcceleration_p, radiiAcceleration_rhoFilt, radiiAcceleration_pFilt, ...
%   radiiPGDI_rho, radiiPGDI_p, radiiPGDI_rhoFilt, radiiPGDI_pFilt, ...
%   radiiPGDISmoothed_rho, radiiPGDISmoothed_p, radiiPGDISmoothed_rhoFilt, radiiPGDISmoothed_pFilt, ...
%   radiiAmplitude_rho, radiiAmplitude_p, radiiAmplitude_rhoFilt, radiiAmplitude_pFilt, ...
%   radiiPower_rho, radiiPower_p, radiiPower_rhoFilt, radiiPower_pFilt, ...
%   radiiFrequency_rho, radiiFrequency_p, radiiFrequency_rhoFilt, radiiFrequency_pFilt, ...
%   speedAmplitude_rho, speedAmplitude_p, speedFrequency_rho, speedFrequency_p] = ...
%   cycleBehavCorrWrapper(dataFiles, pcaOut, pcaOutFilt, thetaUnwrappedPhaseData, ...
%   thetaFrequencyData, thetaAmplitudeData, thetaPowerData, intervals, nPCs, <options>)
%
% Local wrapper function.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   thetaUnwrappedPhaseData
%   thetaFrequencyData
%   thetaAmplitudeData
%   thetaPowerData
%   intervals
%   nPCs
%   <figFolder>
%   <drawFigs>
%   <saveFigs>
%
% Returns:
%   radiiSpeed_rho
%   radiiSpeed_p
%   radiiSpeed_rhoFilt
%   radiiSpeed_pFilt
%   radiiPGDI_rho
%   radiiPGDI_p
%   radiiPGDI_rhoFilt
%   radiiPGDI_pFilt
%   radiiPGDISmoothed_rho
%   radiiPGDISmoothed_p
%   radiiPGDISmoothed_rhoFilt
%   radiiPGDISmoothed_pFilt
%   radiiAmplitude_rho
%   radiiAmplitude_p
%   radiiAmplitude_rhoFilt
%   radiiAmplitude_pFilt
%   radiiFrequency_rho
%   radiiFrequency_p
%   radiiFrequency_rhoFilt
%   radiiFrequency_pFilt
%   speedAmplitude_rho
%   speedAmplitude_p
%   speedFrequency_rho
%   speedFrequency_p
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  thetaUnwrappedPhaseData
  thetaFrequencyData
  thetaAmplitudeData
  thetaPowerData
  intervals
  nPCs
  options.figFolder = [];
  options.drawFigs = true
  options.saveFigs = true
end

% Create data containers
radiiSpeed_rho = cell(numel(dataFiles),1);
radiiSpeed_p = cell(numel(dataFiles),1);
radiiSpeed_rhoFilt = cell(numel(dataFiles),1);
radiiSpeed_pFilt = cell(numel(dataFiles),1);
radiiAcceleration_rho = cell(numel(dataFiles),1);
radiiAcceleration_p = cell(numel(dataFiles),1);
radiiAcceleration_rhoFilt = cell(numel(dataFiles),1);
radiiAcceleration_pFilt = cell(numel(dataFiles),1);
radiiPGDI_rho = cell(numel(dataFiles),1);
radiiPGDI_p = cell(numel(dataFiles),1);
radiiPGDI_rhoFilt = cell(numel(dataFiles),1);
radiiPGDI_pFilt = cell(numel(dataFiles),1);
radiiPGDISmoothed_rho = cell(numel(dataFiles),1);
radiiPGDISmoothed_p = cell(numel(dataFiles),1);
radiiPGDISmoothed_rhoFilt = cell(numel(dataFiles),1);
radiiPGDISmoothed_pFilt = cell(numel(dataFiles),1);
radiiAmplitude_rho = cell(numel(dataFiles),1);
radiiAmplitude_p = cell(numel(dataFiles),1);
radiiAmplitude_rhoFilt = cell(numel(dataFiles),1);
radiiAmplitude_pFilt = cell(numel(dataFiles),1);
radiiFrequency_rho = cell(numel(dataFiles),1);
radiiFrequency_p = cell(numel(dataFiles),1);
radiiFrequency_rhoFilt = cell(numel(dataFiles),1);
radiiFrequency_pFilt = cell(numel(dataFiles),1);
radiiPower_rho = cell(numel(dataFiles),1);
radiiPower_p = cell(numel(dataFiles),1);
radiiPower_rhoFilt = cell(numel(dataFiles),1);
radiiPower_pFilt = cell(numel(dataFiles),1);
speedAmplitude_rho = cell(numel(dataFiles),1);
speedAmplitude_p = cell(numel(dataFiles),1);
speedFrequency_rho = cell(numel(dataFiles),1);
speedFrequency_p = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})

    % Load data
    % Behaviour
    circularTrackFile = strrep(dataFiles{animal}{session}, '*', 'circular_track.behavior');
    if exist(circularTrackFile,'file') && ~isempty(intervals{animal}{session})
      load(circularTrackFile);
      % Theta phase
      [phaseTimestamps, inds] = selectArrayValues( ...
        thetaUnwrappedPhaseData{animal}{session}.timestamps, intervals{animal}{session});

      % Interpolate PCA scores
      pcaOut{animal}{session}.score = interp1( ...
        pcaOut{animal}{session}.timestamps, ...
        pcaOut{animal}{session}.score(:,1:nPCs), phaseTimestamps);
      pcaOutFilt{animal}{session}.score = interp1( ...
        pcaOutFilt{animal}{session}.timestamps, ...
        pcaOutFilt{animal}{session}.score(:,1:nPCs), phaseTimestamps);

      cycles = ceil((thetaUnwrappedPhaseData{animal}{session}.data(inds)+pi)./(2*pi));
      allCycles = unique(cycles);
      % PGDI
      pgdFile = strrep(dataFiles{animal}{session}, '*', 'pgdIndex_travellingSpikingThetaWave.timeseries');
      load(pgdFile);
      pgd = pgdIndex_travellingSpikingThetaWave.data(inds);
      pgdFile = strrep(dataFiles{animal}{session}, '*', 'pgdIndexSmoothed_travellingSpikingThetaWave.timeseries');
      load(pgdFile);
      pgdSmoothed = pgdIndexSmoothed.data(inds);
      % Interpolate theta frequency
      thetaFrequency = interp1(thetaFrequencyData{animal}{session}.timestamps, ...
        thetaFrequencyData{animal}{session}.data, phaseTimestamps);
      % Theta amplitude
      thetaAmplitude = thetaAmplitudeData{animal}{session}.data(inds);
      % Theta power
      thetaPower = interp1(thetaPowerData{animal}{session}.timestamps, thetaPowerData{animal}{session}.data, phaseTimestamps);
      % Interpolate speed
      speed = interp1(circular_track.timestamps, circular_track.speed, phaseTimestamps);


      % Calculate radii
      radii{animal}{session} = zeros(1,numel(allCycles)); %#ok<*AGROW>
      centres{animal}{session} = zeros(nPCs,numel(allCycles));
      radiiFilt{animal}{session} = zeros(1,numel(allCycles));
      centresFilt{animal}{session} = zeros(nPCs,numel(allCycles));
      speeds{animal}{session} = zeros(1,numel(allCycles));
      accelerations{animal}{session} = zeros(1,numel(allCycles));
      pgds{animal}{session} = zeros(1,numel(allCycles));
      pgdsSmoothed{animal}{session} = zeros(1,numel(allCycles));
      thetaAmplitudes{animal}{session} = zeros(1,numel(allCycles));
      thetaPowers{animal}{session} = zeros(1,numel(allCycles));
      thetaFrequencies{animal}{session} = zeros(1,numel(allCycles));
      for iCycle = 1:numel(allCycles)
        dataInds = cycles == allCycles(iCycle);
        [radii{animal}{session}(iCycle), centres{animal}{session}(1:nPCs,iCycle)] = ...
          effectiveRotationRadius(pcaOut{animal}{session}.score(dataInds,1:nPCs), nPCs=nPCs);
        speeds{animal}{session}(iCycle) = mean(speed(dataInds),'omitnan');
        acceleration = diff(speed(dataInds));
        if isempty(acceleration)
          accelerations{animal}{session}(iCycle) = 0;
        else
          accelerations{animal}{session}(iCycle) = mean([acceleration(1) acceleration]);
        end
        pgds{animal}{session}(iCycle) = mean(pgd(dataInds),'omitnan');
        pgdsSmoothed{animal}{session}(iCycle) = mean(pgdSmoothed(dataInds),'omitnan');
        thetaAmplitudes{animal}{session}(iCycle) = mean(thetaAmplitude(dataInds),'omitnan');
        thetaPowers{animal}{session}(iCycle) = mean(thetaPower(dataInds),'omitnan');
        thetaFrequencies{animal}{session}(iCycle) = mean(thetaFrequency(dataInds),'omitnan');
        %figure; plot3(pcaOut{animal}{session}.score(dataInds,1), pcaOut{animal}{session}.score(:,2), pcaOut{animal}{session}.score(:,3)); hold on
        %plot3([rotationCentre(1,1) rotationCentre(1,1)], [rotationCentre(1,2) rotationCentre(1,2)+radii{animal}{session}(iCycle)], [rotationCentre(1,3) rotationCentre(1,3)]); hold off

        [radiiFilt{animal}{session}(iCycle), centresFilt{animal}{session}(1:nPCs,iCycle)] = ...
          effectiveRotationRadius(pcaOutFilt{animal}{session}.score(dataInds,1:nPCs), nPCs=nPCs);
        %figure; plot3(pcaOutFilt{animal}{session}.score(dataInds,1), pcaOutFilt{animal}{session}.score(:,2), pcaOutFilt{animal}{session}.score(:,3)); hold on
        %plot3([rotationCentre(1,1) rotationCentre(1,1)], [rotationCentre(1,2) rotationCentre(1,2)+radii{animal}{session}(iCycle)], [rotationCentre(1,3) rotationCentre(1,3)]); hold off
      end

      % Correlate rotation radii with speed
      inds = ~isnan(radii{animal}{session}) & ~isnan(speeds{animal}{session});
      if sum(inds)
        [radiiSpeed_rho{animal}{session}, radiiSpeed_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', speeds{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), speeds{animal}{session}(inds), ...
            radiiSpeed_rho{animal}{session}, radiiSpeed_p{animal}{session}, ...
            'Radius (a.u.)', 'Speed (cm/s)', options.figFolder, ...
            'Convolved 3PC radius v speed', saveFig=options.saveFigs);
          close(fH1);
        end
      else
        radiiSpeed_rho{animal}{session} = [];
        radiiSpeed_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(speeds{animal}{session});
      if sum(inds)
        [radiiSpeed_rhoFilt{animal}{session}, radiiSpeed_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', speeds{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), speeds{animal}{session}(inds), ...
            radiiSpeed_rhoFilt{animal}{session}, radiiSpeed_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'Speed (cm/s)', options.figFolder, ...
            'Filtered 3PC radius v speed', saveFig=options.saveFigs);
          close(fH2);
        end
      else
        radiiSpeed_rhoFilt{animal}{session} = [];
        radiiSpeed_pFilt{animal}{session} = [];
      end

      % Correlate rotation radii with acceleration
      inds = ~isnan(radii{animal}{session}) & ~isnan(accelerations{animal}{session});
      if sum(inds)
        [radiiAcceleration_rho{animal}{session}, radiiAcceleration_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', accelerations{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), accelerations{animal}{session}(inds), ...
            radiiAcceleration_rho{animal}{session}, radiiAcceleration_p{animal}{session}, ...
            'Radius (a.u.)', 'Acceleration (cm/s^2)', options.figFolder, ...
            'Convolved 3PC radius v acceleration', saveFig=options.saveFigs);
          close(fH1);
        end
      else
        radiiAcceleration_rho{animal}{session} = [];
        radiiAcceleration_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(accelerations{animal}{session});
      if sum(inds)
        [radiiAcceleration_rhoFilt{animal}{session}, radiiAcceleration_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', accelerations{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), accelerations{animal}{session}(inds), ...
            radiiAcceleration_rhoFilt{animal}{session}, radiiAcceleration_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'Acceleration (cm/s^2)', options.figFolder, ...
            'Filtered 3PC radius v acceleration', saveFig=options.saveFigs);
          close(fH2);
        end
      else
        radiiAcceleration_rhoFilt{animal}{session} = [];
        radiiAcceleration_pFilt{animal}{session} = [];
      end

      % Correlate rotation radii with PGDI
      inds = ~isnan(radii{animal}{session}) & ~isnan(pgds{animal}{session});
      if sum(inds)
        [radiiPGDI_rho{animal}{session}, radiiPGDI_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', pgds{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), pgds{animal}{session}(inds), ...
            radiiPGDI_rho{animal}{session}, radiiPGDI_p{animal}{session}, ...
            'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
            'Convolved 3PC radius v PGDI', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        radiiPGDI_rho{animal}{session} = [];
        radiiPGDI_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(pgds{animal}{session});
      if sum(inds)
        [radiiPGDI_rhoFilt{animal}{session}, radiiPGDI_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', pgds{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), pgds{animal}{session}(inds), ...
            radiiPGDI_rhoFilt{animal}{session}, radiiPGDI_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
            'Filtered 3PC radius v PGDI', saveFig=options.saveFigs);
          close(fH2)
        end
      else
        radiiPGDI_rhoFilt{animal}{session} = [];
        radiiPGDI_pFilt{animal}{session} = [];
      end

      % Correlate rotation radii with smoothed PGDI
      inds = ~isnan(radii{animal}{session}) & ~isnan(pgdsSmoothed{animal}{session});
      if sum(inds)
        [radiiPGDISmoothed_rho{animal}{session}, radiiPGDISmoothed_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', pgdsSmoothed{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), pgdsSmoothed{animal}{session}(inds), ...
            radiiPGDISmoothed_rho{animal}{session}, radiiPGDISmoothed_p{animal}{session}, ...
            'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
            'Convolved 3PC radius v smoothed PGDI', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        radiiPGDISmoothed_rho{animal}{session} = [];
        radiiPGDISmoothed_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(pgdsSmoothed{animal}{session});
      if sum(inds)
        [radiiPGDISmoothed_rhoFilt{animal}{session}, radiiPGDISmoothed_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', pgdsSmoothed{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), pgdsSmoothed{animal}{session}(inds), ...
            radiiPGDISmoothed_rhoFilt{animal}{session}, radiiPGDISmoothed_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
            'Filtered 3PC radius v smoothed PGDI', saveFig=options.saveFigs);
          close(fH2)
        end
      else
        radiiPGDISmoothed_rhoFilt{animal}{session} = [];
        radiiPGDISmoothed_pFilt{animal}{session} = [];
      end

      % Correlate rotation radii with theta oscillation amplitude
      inds = ~isnan(radii{animal}{session}) & ~isnan(thetaAmplitudes{animal}{session});
      if sum(inds)
        [radiiAmplitude_rho{animal}{session}, radiiAmplitude_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', thetaAmplitudes{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), thetaAmplitudes{animal}{session}(inds), ...
            radiiAmplitude_rho{animal}{session}, radiiAmplitude_p{animal}{session}, ...
            'Radius (a.u.)', 'Theta oscillation amplitude (a.u.)', options.figFolder, ...
            'Convolved 3PC radius v osc amplitude', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        radiiAmplitude_rho{animal}{session} = [];
        radiiAmplitude_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(thetaAmplitudes{animal}{session});
      if sum(inds)
        [radiiAmplitude_rhoFilt{animal}{session}, radiiAmplitude_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', thetaAmplitudes{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), thetaAmplitudes{animal}{session}(inds), ...
            radiiAmplitude_rhoFilt{animal}{session}, radiiAmplitude_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'Theta oscillation amplitude (a.u.)', options.figFolder, ...
            'Filtered 3PC radius v osc amplitude', saveFig=options.saveFigs);
          close(fH2)
        end
      else
        radiiAmplitude_rhoFilt{animal}{session} = [];
        radiiAmplitude_pFilt{animal}{session} = [];
      end

      % Correlate rotation radii with theta power
      inds = ~isnan(radii{animal}{session}) & ~isnan(thetaPowers{animal}{session});
      if sum(inds)
        [radiiPower_rho{animal}{session}, radiiPower_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', thetaPowers{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), thetaPowers{animal}{session}(inds), ...
            radiiPower_rho{animal}{session}, radiiPower_p{animal}{session}, ...
            'Radius (a.u.)', 'Theta power (a.u.)', options.figFolder, ...
            'Convolved 3PC radius v theta power', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        radiiPower_rho{animal}{session} = [];
        radiiPower_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(thetaPowers{animal}{session});
      if sum(inds)
        [radiiPower_rhoFilt{animal}{session}, radiiPower_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', thetaPowers{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), thetaPowers{animal}{session}(inds), ...
            radiiPower_rhoFilt{animal}{session}, radiiPower_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'Theta power (a.u.)', options.figFolder, ...
            'Filtered 3PC radius v theta power', saveFig=options.saveFigs);
          close(fH2)
        end
      else
        radiiPower_rhoFilt{animal}{session} = [];
        radiiPower_pFilt{animal}{session} = [];
      end

      % Correlate rotation radii with theta oscillation frequency
      inds = ~isnan(radii{animal}{session}) & ~isnan(thetaFrequencies{animal}{session});
      if sum(inds)
        [radiiFrequency_rho{animal}{session}, radiiFrequency_p{animal}{session}] = corr( ...
          radii{animal}{session}(inds)', thetaFrequencies{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radii{animal}{session}(inds), thetaFrequencies{animal}{session}(inds), ...
            radiiFrequency_rho{animal}{session}, radiiFrequency_p{animal}{session}, ...
            'Radius (a.u.)', 'Theta oscillation frequency (Hz)', options.figFolder, ...
            'Convolved 3PC radius v osc frequency', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        radiiFrequency_rho{animal}{session} = [];
        radiiFrequency_p{animal}{session} = [];
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(thetaFrequencies{animal}{session});
      if sum(inds)
        [radiiFrequency_rhoFilt{animal}{session}, radiiFrequency_pFilt{animal}{session}] = corr( ...
          radiiFilt{animal}{session}(inds)', thetaFrequencies{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
            radiiFilt{animal}{session}(inds), thetaFrequencies{animal}{session}(inds), ...
            radiiFrequency_rhoFilt{animal}{session}, radiiFrequency_pFilt{animal}{session}, ...
            'Radius (a.u.)', 'Theta oscillation frequency (Hz)', options.figFolder, ...
            'Filtered 3PC radius v osc frequency', saveFig=options.saveFigs);
          close(fH2)
        end
      else
        radiiFrequency_rhoFilt{animal}{session} = [];
        radiiFrequency_pFilt{animal}{session} = [];
      end

      % Correlate speed with theta oscillation amplitude
      inds = ~isnan(speeds{animal}{session}) & ~isnan(thetaAmplitudes{animal}{session});
      if sum(inds)
        [speedsAmplitude_rho{animal}{session}, speedsAmplitude_p{animal}{session}] = corr( ...
          speeds{animal}{session}(inds)', thetaAmplitudes{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            speeds{animal}{session}(inds), thetaAmplitudes{animal}{session}(inds), ...
            speedsAmplitude_rho{animal}{session}, speedsAmplitude_p{animal}{session}, ...
            'Speed (cm/s)', 'Theta oscillation amplitude (a.u)', options.figFolder, ...
            'Speed v osc amplitude', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        speedsAmplitude_rho{animal}{session} = [];
        speedsAmplitude_p{animal}{session} = [];
      end

      % Correlate speed with theta oscillation frequency
      inds = ~isnan(speeds{animal}{session}) & ~isnan(thetaFrequencies{animal}{session});
      if sum(inds)
        [speedsFrequency_rho{animal}{session}, speedsFrequency_p{animal}{session}] = corr( ...
          speeds{animal}{session}(inds)', thetaFrequencies{animal}{session}(inds)', 'type','Spearman');
        if options.drawFigs
          fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
            speeds{animal}{session}(inds), thetaFrequencies{animal}{session}(inds), ...
            speedsFrequency_rho{animal}{session}, speedsFrequency_p{animal}{session}, ...
            'Speed (cm/s)', '{Theta} oscillation frequency (Hz)', options.figFolder, ...
            'Speed v osc frequency', saveFig=options.saveFigs);
          close(fH1)
        end
      else
        speedsFrequency_rho{animal}{session} = [];
        speedsFrequency_p{animal}{session} = [];
      end
    else
      radiiSpeed_rho{animal}{session} = [];
      radiiSpeed_p{animal}{session} = [];
      radiiSpeed_rhoFilt{animal}{session} = [];
      radiiSpeed_pFilt{animal}{session} = [];
      radiiAcceleration_rho{animal}{session} = [];
      radiiAcceleration_p{animal}{session} = [];
      radiiAcceleration_rhoFilt{animal}{session} = [];
      radiiAcceleration_pFilt{animal}{session} = [];
      radiiPGDI_rho{animal}{session} = [];
      radiiPGDI_p{animal}{session} = [];
      radiiPGDI_rhoFilt{animal}{session} = [];
      radiiPGDI_pFilt{animal}{session} = [];
      radiiPGDISmoothed_rho{animal}{session} = [];
      radiiPGDISmoothed_p{animal}{session} = [];
      radiiPGDISmoothed_rhoFilt{animal}{session} = [];
      radiiPGDISmoothed_pFilt{animal}{session} = [];
      radiiAmplitude_rho{animal}{session} = [];
      radiiAmplitude_p{animal}{session} = [];
      radiiAmplitude_rhoFilt{animal}{session} = [];
      radiiAmplitude_pFilt{animal}{session} = [];
      radiiPower_rho{animal}{session} = [];
      radiiPower_p{animal}{session} = [];
      radiiPower_rhoFilt{animal}{session} = [];
      radiiPower_pFilt{animal}{session} = [];
      radiiFrequency_rho{animal}{session} = [];
      radiiFrequency_p{animal}{session} = [];
      radiiFrequency_rhoFilt{animal}{session} = [];
      radiiFrequency_pFilt{animal}{session} = [];
      speedAmplitude_rho{animal}{session} = [];
      speedAmplitude_p{animal}{session} = [];
      speedFrequency_rho{animal}{session} = [];
      speedFrequency_p{animal}{session} = [];
    end
  end
end
end


function [radiiPGDI_rho, radiiPGDI_p, radiiPGDI_rhoFilt, radiiPGDI_pFilt, ...
  radiiPGDISmoothed_rho, radiiPGDISmoothed_p, radiiPGDISmoothed_rhoFilt, radiiPGDISmoothed_pFilt, ...
  radiiAmplitude_rho, radiiAmplitude_p, radiiAmplitude_rhoFilt, radiiAmplitude_pFilt, ...
  radiiPower_rho, radiiPower_p, radiiPower_rhoFilt, radiiPower_pFilt, ...
  radiiFrequency_rho, radiiFrequency_p, radiiFrequency_rhoFilt, radiiFrequency_pFilt] = ...
  cycleCorrWrapper(dataFiles, pcaOut, pcaOutFilt, thetaUnwrappedPhaseData, ...
  thetaFrequencyData, thetaAmplitudeData, thetaPowerData, intervals, nPCs, options)
% [radiiPGDI_rho, radiiPGDI_p, radiiPGDI_rhoFilt, radiiPGDI_pFilt, ...
%   radiiPGDISmoothed_rho, radiiPGDISmoothed_p, radiiPGDISmoothed_rhoFilt, radiiPGDISmoothed_pFilt, ...
%   radiiAmplitude_rho, radiiAmplitude_p, radiiAmplitude_rhoFilt, radiiAmplitude_pFilt, ...
%   radiiPower_rho, radiiPower_p, radiiPower_rhoFilt, radiiPower_pFilt, ...
%   radiiFrequency_rho, radiiFrequency_p, radiiFrequency_rhoFilt, radiiFrequency_pFilt] = ...
%   cycleCorrWrapper(dataFiles, pcaOut, pcaOutFilt, thetaUnwrappedPhaseData, ...
%   thetaFrequencyData, thetaAmplitudeData, thetaPowerData, intervals, nPCs, <options>)
%
% Local wrapper function.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   thetaUnwrappedPhaseData
%   thetaFrequencyData
%   thetaAmplitudeData
%   thetaPowerData
%   intervals
%   nPCs
%   <figFolder>
%   <drawFigs>
%   <saveFigs>
%
% Returns:
%   radiiPGDI_rho
%   radiiPGDI_p
%   radiiPGDI_rhoFilt
%   radiiPGDI_pFilt
%   radiiPGDISmoothed_rho
%   radiiPGDISmoothed_p
%   radiiPGDISmoothed_rhoFilt
%   radiiPGDISmoothed_pFilt
%   radiiAmplitude_rho
%   radiiAmplitude_p
%   radiiAmplitude_rhoFilt
%   radiiAmplitude_pFilt
%   radiiPower_rho
%   radiiPower_p
%   radiiPower_rhoFilt
%   radiiPower_pFilt
%   radiiFrequency_rho
%   radiiFrequency_p
%   radiiFrequency_rhoFilt
%   radiiFrequency_pFilt
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  thetaUnwrappedPhaseData
  thetaFrequencyData
  thetaAmplitudeData
  thetaPowerData
  intervals
  nPCs
  options.figFolder = [];
  options.drawFigs = true
  options.saveFigs = true
end

% Create data containers
radiiPGDI_rho = cell(numel(dataFiles),1);
radiiPGDI_p = cell(numel(dataFiles),1);
radiiPGDI_rhoFilt = cell(numel(dataFiles),1);
radiiPGDI_pFilt = cell(numel(dataFiles),1);
radiiPGDISmoothed_rho = cell(numel(dataFiles),1);
radiiPGDISmoothed_p = cell(numel(dataFiles),1);
radiiPGDISmoothed_rhoFilt = cell(numel(dataFiles),1);
radiiPGDISmoothed_pFilt = cell(numel(dataFiles),1);
radiiAmplitude_rho = cell(numel(dataFiles),1);
radiiAmplitude_p = cell(numel(dataFiles),1);
radiiAmplitude_rhoFilt = cell(numel(dataFiles),1);
radiiAmplitude_pFilt = cell(numel(dataFiles),1);
radiiPower_rho = cell(numel(dataFiles),1);
radiiPower_p = cell(numel(dataFiles),1);
radiiPower_rhoFilt = cell(numel(dataFiles),1);
radiiPower_pFilt = cell(numel(dataFiles),1);
radiiFrequency_rho = cell(numel(dataFiles),1);
radiiFrequency_p = cell(numel(dataFiles),1);
radiiFrequency_rhoFilt = cell(numel(dataFiles),1);
radiiFrequency_pFilt = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})

    % Load data
    if ~isempty(thetaUnwrappedPhaseData{animal}{session}) && ~isempty(intervals{animal}{session})
      % Theta phase
      [phaseTimestamps, inds] = selectArrayValues( ...
        thetaUnwrappedPhaseData{animal}{session}.timestamps, intervals{animal}{session});

      % Interpolate PCA scores
      pcaOut{animal}{session}.score = interp1( ...
        pcaOut{animal}{session}.timestamps, ...
        pcaOut{animal}{session}.score(:,1:nPCs), phaseTimestamps);
      pcaOutFilt{animal}{session}.score = interp1( ...
        pcaOutFilt{animal}{session}.timestamps, ...
        pcaOutFilt{animal}{session}.score(:,1:nPCs), phaseTimestamps);

      cycles = ceil((thetaUnwrappedPhaseData{animal}{session}.data(inds)+pi)./(2*pi));
      allCycles = unique(cycles);
      % PGDI
      pgdFile = strrep(dataFiles{animal}{session}, '*', 'pgdIndex_travellingSpikingThetaWave.timeseries');
      load(pgdFile);
      pgd = pgdIndex_travellingSpikingThetaWave.data(inds);
      pgdFile = strrep(dataFiles{animal}{session}, '*', 'pgdIndexSmoothed_travellingSpikingThetaWave.timeseries');
      load(pgdFile);
      pgdSmoothed = pgdIndexSmoothed.data(inds);
      % Interpolate theta frequency
      thetaFrequency = interp1(thetaFrequencyData{animal}{session}.timestamps, ...
        thetaFrequencyData{animal}{session}.data, phaseTimestamps);
      % Theta amplitude
      thetaAmplitude = thetaAmplitudeData{animal}{session}.data(inds);
      % Theta power
      thetaPower = interp1(thetaPowerData{animal}{session}.timestamps, thetaPowerData{animal}{session}.data, phaseTimestamps);

      % Calculate radii
      radii{animal}{session} = zeros(1,numel(allCycles)); %#ok<*AGROW>
      centres{animal}{session} = zeros(nPCs,numel(allCycles));
      radiiFilt{animal}{session} = zeros(1,numel(allCycles));
      centresFilt{animal}{session} = zeros(nPCs,numel(allCycles));
      pgds{animal}{session} = zeros(1,numel(allCycles));
      pgdsSmoothed{animal}{session} = zeros(1,numel(allCycles));
      thetaAmplitudes{animal}{session} = zeros(1,numel(allCycles));
      thetaFrequencies{animal}{session} = zeros(1,numel(allCycles));
      for iCycle = 1:numel(allCycles)
        dataInds = cycles == allCycles(iCycle);
        [radii{animal}{session}(iCycle), centres{animal}{session}(1:nPCs,iCycle)] = ...
          effectiveRotationRadius(pcaOut{animal}{session}.score(dataInds,1:nPCs), nPCs=nPCs);
        pgds{animal}{session}(iCycle) = mean(pgd(dataInds),'omitnan');
        pgdsSmoothed{animal}{session}(iCycle) = mean(pgdSmoothed(dataInds),'omitnan');
        thetaAmplitudes{animal}{session}(iCycle) = mean(thetaAmplitude(dataInds),'omitnan');
        thetaPowers{animal}{session}(iCycle) = mean(thetaPower(dataInds),'omitnan');
        thetaFrequencies{animal}{session}(iCycle) = mean(thetaFrequency(dataInds),'omitnan');
        %figure; plot3(pcaOut{animal}{session}.score(dataInds,1), pcaOut{animal}{session}.score(:,2), pcaOut{animal}{session}.score(:,3)); hold on
        %plot3([rotationCentre(1,1) rotationCentre(1,1)], [rotationCentre(1,2) rotationCentre(1,2)+radii{animal}{session}(iCycle)], [rotationCentre(1,3) rotationCentre(1,3)]); hold off

        [radiiFilt{animal}{session}(iCycle), centresFilt{animal}{session}(1:nPCs,iCycle)] = ...
          effectiveRotationRadius(pcaOutFilt{animal}{session}.score(dataInds,1:nPCs), nPCs=nPCs);
        %figure; plot3(pcaOutFilt{animal}{session}.score(dataInds,1), pcaOutFilt{animal}{session}.score(:,2), pcaOutFilt{animal}{session}.score(:,3)); hold on
        %plot3([rotationCentre(1,1) rotationCentre(1,1)], [rotationCentre(1,2) rotationCentre(1,2)+radii{animal}{session}(iCycle)], [rotationCentre(1,3) rotationCentre(1,3)]); hold off
      end

      % Correlate rotation radii with PGDI
      inds = ~isnan(radii{animal}{session}) & ~isnan(pgds{animal}{session});
      [radiiPGDI_rho{animal}{session}, radiiPGDI_p{animal}{session}] = corr( ...
        radii{animal}{session}(inds)', pgds{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radii{animal}{session}(inds), pgds{animal}{session}(inds), ...
          radiiPGDI_rho{animal}{session}, radiiPGDI_p{animal}{session}, ...
          'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
          'Convolved 3PC radius v PGDI', saveFig=options.saveFigs);
        close(fH1)
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(pgds{animal}{session});
      [radiiPGDI_rhoFilt{animal}{session}, radiiPGDI_pFilt{animal}{session}] = corr( ...
        radiiFilt{animal}{session}(inds)', pgds{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radiiFilt{animal}{session}(inds), pgds{animal}{session}(inds), ...
          radiiPGDI_rhoFilt{animal}{session}, radiiPGDI_pFilt{animal}{session}, ...
          'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
          'Filtered 3PC radius v PGDI', saveFig=options.saveFigs);
        close(fH2)
      end

      % Correlate rotation radii with smoothed PGDI
      inds = ~isnan(radii{animal}{session}) & ~isnan(pgdsSmoothed{animal}{session});
      [radiiPGDISmoothed_rho{animal}{session}, radiiPGDISmoothed_p{animal}{session}] = corr( ...
        radii{animal}{session}(inds)', pgdsSmoothed{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radii{animal}{session}(inds), pgdsSmoothed{animal}{session}(inds), ...
          radiiPGDISmoothed_rho{animal}{session}, radiiPGDISmoothed_p{animal}{session}, ...
          'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
          'Convolved 3PC radius v smoothed PGDI', saveFig=options.saveFigs);
        close(fH1)
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(pgdsSmoothed{animal}{session});
      [radiiPGDISmoothed_rhoFilt{animal}{session}, radiiPGDISmoothed_pFilt{animal}{session}] = corr( ...
        radiiFilt{animal}{session}(inds)', pgdsSmoothed{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radiiFilt{animal}{session}(inds), pgdsSmoothed{animal}{session}(inds), ...
          radiiPGDISmoothed_rhoFilt{animal}{session}, radiiPGDISmoothed_pFilt{animal}{session}, ...
          'Radius (a.u.)', 'PGDI (a.u.)', options.figFolder, ...
          'Filtered 3PC radius v smoothed PGDI', saveFig=options.saveFigs);
        close(fH2)
      end

      % Correlate rotation radii with theta oscillation amplitude
      inds = ~isnan(radii{animal}{session}) & ~isnan(thetaAmplitudes{animal}{session});
      [radiiAmplitude_rho{animal}{session}, radiiAmplitude_p{animal}{session}] = corr( ...
        radii{animal}{session}(inds)', thetaAmplitudes{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radii{animal}{session}(inds), thetaAmplitudes{animal}{session}(inds), ...
          radiiAmplitude_rho{animal}{session}, radiiAmplitude_p{animal}{session}, ...
          'Radius (a.u.)', 'Theta oscillation amplitude (a.u.)', options.figFolder, ...
          'Convolved 3PC radius v osc amplitude', saveFig=options.saveFigs);
        close(fH1)
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(thetaAmplitudes{animal}{session});
      [radiiAmplitude_rhoFilt{animal}{session}, radiiAmplitude_pFilt{animal}{session}] = corr( ...
        radiiFilt{animal}{session}(inds)', thetaAmplitudes{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radiiFilt{animal}{session}(inds), thetaAmplitudes{animal}{session}(inds), ...
          radiiAmplitude_rhoFilt{animal}{session}, radiiAmplitude_pFilt{animal}{session}, ...
          'Radius (a.u.)', 'Theta oscillation amplitude (a.u.)', options.figFolder, ...
          'Filtered 3PC radius v osc amplitude', saveFig=options.saveFigs);
        close(fH2)
      end

      % Correlate rotation radii with theta power
      inds = ~isnan(radii{animal}{session}) & ~isnan(thetaPowers{animal}{session});
      [radiiAmplitude_rho{animal}{session}, radiiAmplitude_p{animal}{session}] = corr( ...
        radii{animal}{session}(inds)', thetaPowers{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radii{animal}{session}(inds), thetaPowers{animal}{session}(inds), ...
          radiiAmplitude_rho{animal}{session}, radiiAmplitude_p{animal}{session}, ...
          'Radius (a.u.)', 'Theta power (a.u.)', options.figFolder, ...
          'Convolved 3PC radius v theta power', saveFig=options.saveFigs);
        close(fH1)
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(thetaPowers{animal}{session});
      [radiiAmplitude_rhoFilt{animal}{session}, radiiAmplitude_pFilt{animal}{session}] = corr( ...
        radiiFilt{animal}{session}(inds)', thetaPowers{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radiiFilt{animal}{session}(inds), thetaPowers{animal}{session}(inds), ...
          radiiAmplitude_rhoFilt{animal}{session}, radiiAmplitude_pFilt{animal}{session}, ...
          'Radius (a.u.)', 'Theta power (a.u.)', options.figFolder, ...
          'Filtered 3PC radius v theta power', saveFig=options.saveFigs);
        close(fH2)
      end

      % Correlate rotation radii with theta oscillation frequency
      inds = ~isnan(radii{animal}{session}) & ~isnan(thetaFrequencies{animal}{session});
      [radiiFrequency_rho{animal}{session}, radiiFrequency_p{animal}{session}] = corr( ...
        radii{animal}{session}(inds)', thetaFrequencies{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH1 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radii{animal}{session}(inds), thetaFrequencies{animal}{session}(inds), ...
          radiiFrequency_rho{animal}{session}, radiiFrequency_p{animal}{session}, ...
          'Radius (a.u.)', 'Theta oscillation frequency (Hz)', options.figFolder, ...
          'Convolved 3PC radius v osc frequency', saveFig=options.saveFigs);
        close(fH1)
      end
      inds = ~isnan(radiiFilt{animal}{session}) & ~isnan(thetaFrequencies{animal}{session});
      [radiiFrequency_rhoFilt{animal}{session}, radiiFrequency_pFilt{animal}{session}] = corr( ...
        radiiFilt{animal}{session}(inds)', thetaFrequencies{animal}{session}(inds)', 'type','Spearman');
      if options.drawFigs
        fH2 = cycleFigWrapper(dataFiles{animal}{session}, ...
          radiiFilt{animal}{session}(inds), thetaFrequencies{animal}{session}(inds), ...
          radiiFrequency_rhoFilt{animal}{session}, radiiFrequency_pFilt{animal}{session}, ...
          'Radius (a.u.)', 'Theta oscillation frequency (Hz)', options.figFolder, ...
          'Filtered 3PC radius v osc frequency', saveFig=options.saveFigs);
        close(fH2)
      end
    else
      radiiPGDI_rho{animal}{session} = [];
      radiiPGDI_p{animal}{session} = [];
      radiiPGDI_rhoFilt{animal}{session} = [];
      radiiPGDI_pFilt{animal}{session} = [];
      radiiPGDISmoothed_rho{animal}{session} = [];
      radiiPGDISmoothed_p{animal}{session} = [];
      radiiPGDISmoothed_rhoFilt{animal}{session} = [];
      radiiPGDISmoothed_pFilt{animal}{session} = [];
      radiiAmplitude_rho{animal}{session} = [];
      radiiAmplitude_p{animal}{session} = [];
      radiiAmplitude_rhoFilt{animal}{session} = [];
      radiiAmplitude_pFilt{animal}{session} = [];
      radiiPower_rho{animal}{session} = [];
      radiiPower_p{animal}{session} = [];
      radiiPower_rhoFilt{animal}{session} = [];
      radiiPower_pFilt{animal}{session} = [];
      radiiFrequency_rho{animal}{session} = [];
      radiiFrequency_p{animal}{session} = [];
      radiiFrequency_rhoFilt{animal}{session} = [];
      radiiFrequency_pFilt{animal}{session} = [];
    end
  end
end
end


function fH = cycleFigWrapper(dataFile, vec1, vec2, rho, p, xLabel, yLabel, figFolder, figName, options)
% fH = cycleFigWrapper(dataFile, vec1, vec2, rho, p, xLabel, yLabel, figFolder, figName, <saveFig>)
%
% A wrapper function for pcaScript.
%
% Args:
%   dataFile
%   vec1
%   vec2
%   rho
%   p
%   xLabel
%   yLabel
%   figFolder
%   figName
%   <saveFig>
%
% Returns:
%   fH

arguments
  dataFile
  vec1
  vec2
  rho
  p
  xLabel
  yLabel
  figFolder
  figName
  options.saveFig = true
end

% Create the figure folder
[~, animalFolder] = fileparts(fileparts(fileparts(dataFile)));
[~, sessionFolder] = fileparts(fileparts(dataFile));
fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
if ~exist(fullSessionFolder, 'dir')
  mkdir(fullSessionFolder);
end

% Draw the figure
fH = figure;
plot(vec1, vec2, '.', 'MarkerSize',5);

% Text
yLim = ylim;
yDist = yLim(2) - yLim(1);
xLim = xlim;
xDist = xLim(2) - xLim(1);
textStr = ['r=' num2str(rho) ' p=' num2str(p)];
text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);

% Line fit
xAxisStep = xDist/10000;
x = xLim(1):xAxisStep:xLim(2);
[~, slope, coefficients] = fitLine(vec1, vec2, type='linear-linear');
yFit = x.*slope + coefficients(2);
hold on; plot(x, yFit, 'k--'); hold off;
ylim(yLim);
xlim(xLim);

% Labels
xlabel(xLabel);
ylabel(yLabel);
title(figName, 'Interpreter','none');
set(fH, 'Name',figName);

% Save the figure
if options.saveFig
  if ~exist(fullSessionFolder, 'dir')
    mdir(fullSessionFolder);
  end
  figName = strrep(figName, ' ', '_');
  figName = strrep(figName, '(', '_');
  figName = strrep(figName, ')', '_');
  figName = fullfile(fullSessionFolder, figName);
  savefig(fH, figName,'compact');
  title('');
  saveas(fH, figName, 'png');
  saveas(fH, figName, 'pdf');
end
end


function barExplained(dataFiles, pcaOut, pcaOutFilt, figFolder)
% barExplained(dataFiles, pcaOut, pcaOutFilt, figFolder)
%
% Function produces and saves graphs showing cumulative explained variance
% by principal components.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   pcaFigFolder
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  figFolder
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session})

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Calculate a change in the angle of the line fitted to the explained cummulative variance curve
      cumSumDataPoints = numel(pcaOut{animal}{session}.explained);
      angleChangePoints = 2:(cumSumDataPoints - 1);
      csPCA = cumsum(pcaOut{animal}{session}.explained);
      csPCAFilt = cumsum(pcaOutFilt{animal}{session}.explained);
      slopePCA = diff(csPCA);
      slopePCAFilt = diff(csPCAFilt);
      anglePCA = atand(slopePCA);
      anglePCAFilt = atand(slopePCAFilt);
      angleChangePCA = deg2rad(abs(diff(anglePCA)));
      angleChangePCAFilt = deg2rad(abs(diff(anglePCAFilt)));

      % Convolved data only figure
      fH1 = figure;
      bar(pcaOut{animal}{session}.explained); hold on
      yyaxis left
      plot(cumsum(pcaOut{animal}{session}.explained), '.-', 'MarkerSize',10);
      xLim = xlim;
      xlim([0 xLim(2)]);
      ylim([0 100])

      % Text
      yLim = ylim;
      yDist = yLim(2) - yLim(1);
      xLim = xlim;
      xDist = xLim(2) - xLim(1);
      if numel(csPCA) >= 2
        textStr = ['R^2_2=' num2str(round(csPCA(2),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.965*yDist, textStr);
      end
      if numel(csPCA) >= 3
        textStr = ['R^2_3=' num2str(round(csPCA(3),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.89*yDist, textStr);
      end
      if numel(csPCA) >= 6
        textStr = ['R^2_6=' num2str(round(csPCA(6),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.815*yDist, textStr);
      end
      if numel(csPCA) >= 10
        textStr = ['R^2_{10}=' num2str(round(csPCA(10),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.74*yDist, textStr);
      end
      if numel(csPCA) >= 20
        textStr = ['R^2_{20}=' num2str(round(csPCA(20),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.665*yDist, textStr);
      end
      
      % Labels and saving
      xlabel('PC #', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('(Cumulative) explained variance %', ...
        'FontSize',fontSize, 'FontWeight','bold');
      yyaxis right
      plot(angleChangePoints, angleChangePCA); hold off
      ylabel('Tangent line inclination change (rad)', ...
        'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Explained variance by PCs (convolved-only data)';
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');

      % Convolved and band-pass filtered data figure
      fH2 = figure;
      bar(pcaOutFilt{animal}{session}.explained); hold on
      yyaxis left
      plot(cumsum(pcaOutFilt{animal}{session}.explained), '.-', 'MarkerSize',10); hold off
      xLim = xlim;
      xlim([0 xLim(2)]);
      ylim([0 100])

      % Text
      yLim = ylim;
      yDist = yLim(2) - yLim(1);
      xLim = xlim;
      xDist = xLim(2) - xLim(1);
      if numel(csPCAFilt) >= 2
        textStr = ['R^2_2=' num2str(round(csPCAFilt(2),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.965*yDist, textStr);
      end
      if numel(csPCAFilt) >= 3
        textStr = ['R^2_3=' num2str(round(csPCAFilt(3),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.89*yDist, textStr);
      end
      if numel(csPCAFilt) >= 6
        textStr = ['R^2_6=' num2str(round(csPCAFilt(6),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.815*yDist, textStr);
      end
      if numel(csPCAFilt) >= 10
        textStr = ['R^2_{10}=' num2str(round(csPCAFilt(10),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.74*yDist, textStr);
      end
      if numel(csPCAFilt) >= 20
        textStr = ['R^2_{20}=' num2str(round(csPCAFilt(20),2)) '%'];
        text(xLim(1)+0.025*xDist, yLim(1)+0.665*yDist, textStr);
      end
      
      % Labels and saving
      xlabel('PC #', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('(Cumulative) explained variance %', ...
        'FontSize',fontSize, 'FontWeight','bold');
      yyaxis right
      plot(angleChangePoints, angleChangePCAFilt);
      ylabel('Tangent line inclination change (rad)', ...
        'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Explained variance by PCs (filtered data)';
      title(titleStr);
      set(fH2, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH2, figName,'compact');
      title('');
      saveas(fH2, figName, 'png');
      saveas(fH2, figName, 'pdf');

      close([fH1 fH2]);
    end
  end
end
end


function combine2PCs(dataFiles, pcaOut, pcaOutFilt, populationRates, ...
  populationRatesTimesBins, intervals, figFolder, fitCircle)
% combine2PCs(dataFiles, pcaOut, pcaOutFilt, populationRates, ...
%   populationRatesTimesBins, intervals, figFolder, fitCircle)
%
% Function produces and saves graphs showing 2 PCs drawn against each other.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   populationRates
%   populationRatesTimesBins
%   intervals
%   figFolder
%   fitCircle
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  populationRates
  populationRatesTimesBins
  intervals
  figFolder
  fitCircle
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session}) && ~isempty(pcaOut{animal}{session}.score) && ...
        size(pcaOut{animal}{session}.score, 2) > 1
      
      % Order PCs
      [pcInds, pcIndsFilt] = orderPCs(pcaOut, pcaOutFilt, ...
        populationRates, populationRatesTimesBins, intervals, animal, session);

      % Fit a circle
      try
        if fitCircle
          [~, ~, pcaR, ellipse] = fitPlotEllipse([pcaOut{animal}{session}.score(:,pcInds(1)) ...
            pcaOut{animal}{session}.score(:,pcInds(2))], type='ellipseLinearTrace', draw=false);
          [~, ~, pcaFiltR, ellipseFilt] = fitPlotEllipse([pcaOutFilt{animal}{session}.score(:,pcInds(1)) ...
            pcaOutFilt{animal}{session}.score(:,pcInds(2))], type='ellipseLinearTrace', draw=false);
          %[radius, xCentreCoord, yCentreCoord] = circfit( ...
          %  pcaOut{animal}{session}.score(:,pcInds(1)), ...
          %  pcaOut{animal}{session}.score(:,pcInds(2)));
          %pcaR = explainedVarCircleFit(pcaOut{animal}{session}.score(:,pcInds(1)), ...
          %  pcaOut{animal}{session}.score(:,pcInds(1)), xCentreCoord, yCentreCoord, radius);
          %[radiusFilt, xCentreCoordFilt, yCentreCoordFilt] = circfit( ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)));
          %pcaFiltR = explainedVarCircleFit(pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)), xCentreCoordFilt, yCentreCoordFilt, radiusFilt);
        end
      catch
        continue
      end

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Convolved data only figure
      fH1 = figure;
      plot(pcaOut{animal}{session}.score(:,pcInds(1)), ...
        pcaOut{animal}{session}.score(:,pcInds(2)), '.', 'MarkerSize',1);
      if fitCircle % Plot circle fit
        hold on; plot(ellipse(1,:),ellipse(2,:),'k'); hold off;
        %rectangle('Curvature',[1 1],'Position',[xCentreCoord-radius ...
        %  yCentreCoord-radius 2*radius 2*radius]);
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['R^2=' num2str(pcaR)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
      end
      xlabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      ylabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 2 PCs (convolved-only data)';
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');

      % Convolved and band-pass filtered data figure
      fH2 = figure;
      plot(pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
        pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)), ...
        '.', 'MarkerSize',1);
      if fitCircle % Plot circle fit
        hold on; plot(ellipseFilt(1,:),ellipseFilt(2,:),'k'); hold off;
        %rectangle('Curvature',[1 1],'Position',[xCentreCoordFilt-radiusFilt ...
        %  yCentreCoordFilt-radiusFilt 2*radiusFilt 2*radiusFilt]);
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['R^2=' num2str(pcaFiltR)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
      end
      xlabel(['PC' num2str(pcIndsFilt(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      ylabel(['PC' num2str(pcIndsFilt(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 2 PCs (filtered data)';
      title(titleStr);
      set(fH2, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH2, figName,'compact');
      title('');
      saveas(fH2, figName, 'png');
      saveas(fH2, figName, 'pdf');

      close([fH1 fH2]);
    end
  end
end
end


function combine3PCs(dataFiles, pcaOut, pcaOutFilt, figFolder)
% combine3PCs(dataFiles, pcaOut, pcaOutFilt, figFolder)
%
% Function produces and saves graphs showing 3 PCs drawn against each other.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   pcaFigFolder
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  figFolder
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session}) && ...
        size(pcaOut{animal}{session}.score, 2) > 2

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Convolved data only figure
      fH1 = figure;
      plot3(pcaOut{animal}{session}.score(:,1), pcaOut{animal}{session}.score(:,2), ...
        pcaOut{animal}{session}.score(:,3), '.', 'MarkerSize',1);
      xlabel('PC1 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('PC2 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      zlabel('PC3 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 3 PCs (convolved-only data)';
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');

      % Convolved and band-pass filtered data figure
      fH2 = figure;
      plot3(pcaOutFilt{animal}{session}.score(:,1), pcaOutFilt{animal}{session}.score(:,2), ...
        pcaOutFilt{animal}{session}.score(:,3), '.', 'MarkerSize',1);
      xlabel('PC1 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('PC2 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      zlabel('PC3 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 3 PCs (filtered data)';
      title(titleStr);
      set(fH2, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH2, figName,'compact');
      title('');
      saveas(fH2, figName, 'png');
      saveas(fH2, figName, 'pdf');

      close([fH1 fH2]);
    end
  end
end
end


function combine2PCsColour(dataFiles, pcaOut, pcaOutFilt, populationRates, ...
    populationRatesTimesBins, thetaPhaseData, intervals, figFolder, fitCircle)
% combine2PCsColour(dataFiles, pcaOut, pcaOutFilt, populationRates, ...
%   populationRatesTimesBins, thetaPhaseData, intervals, figFolder, fitCircle)
%
% Function produces and saves graphs showing 2 PCs drawn against each other
% with data points colour-coding the theta phase.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   populationRates
%   populationRatesTimesBins
%   thetaPhaseData
%   intervals
%   figFolder
%   fitCircle
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  populationRates
  populationRatesTimesBins
  thetaPhaseData
  intervals
  figFolder
  fitCircle
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session}) && ~isempty(pcaOut{animal}{session}.score) && ...
        size(pcaOut{animal}{session}.score, 2) > 1

      % Order PCs
      [pcInds, pcIndsFilt] = orderPCs(pcaOut, pcaOutFilt, ...
        populationRates, populationRatesTimesBins, intervals, animal, session);

      try
        if fitCircle
          [~, ellipseParams, pcaR, ellipse] = fitPlotEllipse([pcaOut{animal}{session}.score(:,pcInds(1)) ...
            pcaOut{animal}{session}.score(:,pcInds(2))], type='ellipseLinearTrace', draw=false);
          xCentreCoord = ellipseParams.x0(1);
          yCentreCoord = ellipseParams.x0(2);
          [~, ellipseParamsFilt, pcaFiltR, ellipseFilt] = fitPlotEllipse([pcaOutFilt{animal}{session}.score(:,pcInds(1)) ...
            pcaOutFilt{animal}{session}.score(:,pcInds(2))], type='ellipseLinearTrace', draw=false);
          xCentreCoordFilt = ellipseParamsFilt.x0(1);
          yCentreCoordFilt = ellipseParamsFilt.x0(2);
          %[radius, xCentreCoord, yCentreCoord] = circfit( ...
          %  pcaOut{animal}{session}.score(:,pcInds(1)), ...
          %  pcaOut{animal}{session}.score(:,pcInds(2)));
          %pcaR = explainedVarCircleFit(pcaOut{animal}{session}.score(:,pcInds(1)), ...
          %  pcaOut{animal}{session}.score(:,pcInds(1)), xCentreCoord, yCentreCoord, radius);
          %[radiusFilt, xCentreCoordFilt, yCentreCoordFilt] = circfit( ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)));
          %pcaFiltR = explainedVarCircleFit(pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)), xCentreCoordFilt, yCentreCoordFilt, radiusFilt);
        end
      catch
        continue
      end

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Get phase labels
      [~, thetaPhaseInds] = selectArrayValues( ...
        thetaPhaseData{animal}{session}.timestamps, intervals{animal}{session});
      thetaPhaseSession = thetaPhaseData{animal}{session}.data(thetaPhaseInds);
      %colours = hsv;
      %colourInds = ceil(((thetaPhaseSession + pi)./(2*pi)).*size(colours,1));
      %dataColours = colours(colourInds,:);
      %colormap(hsv);
      if numel(thetaPhaseSession) ~= numel(pcaOut{animal}{session}.score(:,1))
        phaseTimestamps = thetaPhaseData{animal}{session}.timestamps(thetaPhaseInds);
        pcaTimestamps = pcaOut{animal}{session}.timestamps;
        [~, phaseInds, pcaInds] = intersect(round(phaseTimestamps.*1000), ...
          round(pcaTimestamps.*1000));
        thetaPhaseSession = thetaPhaseSession(phaseInds);
      else
        pcaInds = 1:numel(pcaOut{animal}{session}.score(:,1));
      end
      if isempty(pcaInds)
        continue
      end

      % Convolved data only figure: PCs
      fH1 = figure;
      scatter(pcaOut{animal}{session}.score(pcaInds,pcInds(1)), ...
        pcaOut{animal}{session}.score(pcaInds,pcInds(2)), ...
        ones(numel(thetaPhaseSession),1), thetaPhaseSession, '.');
      colormap(hsv);
      c = colorbar;
      c.Label.String = 'Theta phase (rad)';
      c.Label.FontSize = fontSize;
      c.Label.FontWeight = 'bold';
      if fitCircle % Plot circle fit
        hold on; plot(ellipse(1,:),ellipse(2,:),'k'); hold off;
        %rectangle('Curvature',[1 1],'Position',[xCentreCoord-radius ...
        %  yCentreCoord-radius 2*radius 2*radius]);
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['R^2=' num2str(pcaR)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
      end
      xlabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      ylabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 2 PCs with phase (convolved-only data)';
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');

      % Convolved data only figure: PC space phase vs theta phase
      if fitCircle
        centreRelativePC1 = pcaOut{animal}{session}.score(pcaInds,pcInds(1)) - xCentreCoord;
        centreRelativePC2 = pcaOut{animal}{session}.score(pcaInds,pcInds(2)) - yCentreCoord;
        pcSpacePhase = atan2(centreRelativePC2, centreRelativePC1);
        %cycleOrigin = datamean(pcSpacePhase(find(diff(ceil(unwrap(thetaPhaseSession)+pi/(2*pi))))+1), 'circularNP');
        %pcSpacePhaseShift = -pi - cycleOrigin;
        %pcSpacePhase = recentrePhase(pcSpacePhase - pcSpacePhaseShift, 0);
        fH2 = figure; plot(thetaPhaseSession, pcSpacePhase, '.', 'MarkerSize',1);

        [r, pval] = corrLinearCircular(thetaPhaseSession, pcSpacePhase, ...
          type='circularnp');
        yLim = [-pi pi];
        yDist = yLim(2) - yLim(1);
        xLim = [-pi pi];
        xDist = xLim(2) - xLim(1);
        textStr = ['r=' num2str(r) ' p=' num2str(pval)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);

        [~, slope, coefficients] = fitLine(thetaPhaseSession, pcSpacePhase, ...
          type='linear-circular-fma');
        xAxisStep = xDist/10000;
        x = xLim(1):xAxisStep:xLim(2);
        yFit = x.*slope + coefficients(2);
        hold on; plot(x, yFit, 'k--'); hold off;
        xlim([-pi pi]); ylim([-pi pi]);
        
        xlabel('Theta phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        ylabel('PC space phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        titleStr = 'PC space phase vs theta phase correlation (convolved-only data)';
        title(titleStr);
        set(fH2, 'Name',titleStr);
        titleStr = strrep(titleStr, ' ', '_');
        titleStr = strrep(titleStr, '(', '_');
        titleStr = strrep(titleStr, ')', '_');
        figName = fullfile(fullSessionFolder, titleStr);
        savefig(fH2, figName,'compact');
        title('');
        saveas(fH2, figName, 'png');
        saveas(fH2, figName, 'pdf');

        nbins=[50 50];
        [phaseMatrix, phaseGridCoords] = hist3([pcSpacePhase(:) ...
          thetaPhaseSession(:)], nbins);

        fH3 = figure;
        s = pcolor(phaseGridCoords{1}, phaseGridCoords{2}, phaseMatrix);
        s.FaceColor = 'interp';
        s.EdgeColor = 'none';
        c = colorbar;
        c.Label.String = 'Count (samples)';
        c.Label.FontSize = fontSize;
        c.Label.FontWeight = 'bold';
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
        xLim = xlim; yLim = ylim;
        hold on; plot(x, yFit, 'k--'); hold off;
        xlim(xLim); ylim(yLim);
        xlabel('Theta phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        ylabel('PC space phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        titleStr = 'PC space phase vs theta phase density map (convolved-only data)';
        title(titleStr);
        set(fH3, 'Name',titleStr);
        titleStr = strrep(titleStr, ' ', '_');
        titleStr = strrep(titleStr, '(', '_');
        titleStr = strrep(titleStr, ')', '_');
        figName = fullfile(fullSessionFolder, titleStr);
        savefig(fH3, figName,'compact');
        title('');
        saveas(fH3, figName, 'png');
        saveas(fH3, figName, 'pdf');
      end

      % Convolved and band-pass filtered data figure
      fH4 = figure;
      scatter(pcaOutFilt{animal}{session}.score(pcaInds,pcIndsFilt(1)), ...
        pcaOutFilt{animal}{session}.score(pcaInds,pcIndsFilt(2)), ...
        ones(numel(thetaPhaseSession),1), thetaPhaseSession, '.');
      colormap(hsv);
      c = colorbar;
      c.Label.String = 'Theta phase (rad)';
      c.Label.FontSize = fontSize;
      c.Label.FontWeight = 'bold';
      if fitCircle % Plot circle fit
        hold on; plot(ellipseFilt(1,:),ellipseFilt(2,:),'k'); hold off;
        %rectangle('Curvature',[1 1],'Position',[xCentreCoordFilt-radiusFilt ...
        %  yCentreCoordFilt-radiusFilt 2*radiusFilt 2*radiusFilt]);
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['R^2=' num2str(pcaFiltR)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
      end
      xlabel(['PC' num2str(pcIndsFilt(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      ylabel(['PC' num2str(pcIndsFilt(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 2 PCs with phase (filtered data)';
      title(titleStr);
      set(fH4, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH4, figName,'compact');
      title('');
      saveas(fH4, figName, 'png');
      saveas(fH4, figName, 'pdf');

      % Convolved data only figure: PC space phase vs theta phase
      if fitCircle
        centreRelativePC1 = pcaOutFilt{animal}{session}.score(pcaInds,pcIndsFilt(1)) - xCentreCoordFilt;
        centreRelativePC2 = pcaOutFilt{animal}{session}.score(pcaInds,pcIndsFilt(2)) - yCentreCoordFilt;
        pcSpacePhase = atan2(centreRelativePC2, centreRelativePC1);
        %cycleOrigin = datamean(pcSpacePhase(find(diff(ceil(unwrap(thetaPhaseSession)+pi/(2*pi))))+1), 'circularNP');
        %pcSpacePhaseShift = -pi - cycleOrigin;
        %pcSpacePhase = recentrePhase(pcSpacePhase - pcSpacePhaseShift, 0);
        fH5 = figure; plot(thetaPhaseSession, pcSpacePhase, '.', 'MarkerSize',1);

        [r, pval] = corrLinearCircular(thetaPhaseSession, pcSpacePhase, ...
          type='circularnp');
        yLim = [-pi pi];
        yDist = yLim(2) - yLim(1);
        xLim = [-pi pi];
        xDist = xLim(2) - xLim(1);
        textStr = ['r=' num2str(r) ' p=' num2str(pval)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);

        [~, slope, coefficients] = fitLine(thetaPhaseSession, pcSpacePhase, ...
          type='linear-circular-fma');
        xAxisStep = xDist/10000;
        x = xLim(1):xAxisStep:xLim(2);
        yFit = x.*slope + coefficients(2);
        hold on; plot(x, yFit, 'k--'); hold off;
        xlim([-pi pi]); ylim([-pi pi]);
        
        xlabel('Theta phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        ylabel('PC space phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        titleStr = 'PC space phase vs theta phase correlation (filtered data)';
        title(titleStr);
        set(fH5, 'Name',titleStr);
        titleStr = strrep(titleStr, ' ', '_');
        titleStr = strrep(titleStr, '(', '_');
        titleStr = strrep(titleStr, ')', '_');
        figName = fullfile(fullSessionFolder, titleStr);
        savefig(fH5, figName,'compact');
        title('');
        saveas(fH5, figName, 'png');
        saveas(fH5, figName, 'pdf');

        nbins=[50 50];
        [phaseMatrix, phaseGridCoords] = hist3([pcSpacePhase(:) ...
          thetaPhaseSession(:)], nbins);

        fH6 = figure;
        s = pcolor(phaseGridCoords{1}, phaseGridCoords{2}, phaseMatrix);
        s.FaceColor = 'interp';
        s.EdgeColor = 'none';
        c = colorbar;
        c.Label.String = 'Count (samples)';
        c.Label.FontSize = fontSize;
        c.Label.FontWeight = 'bold';
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
        xLim = xlim; yLim = ylim;
        hold on; plot(x, yFit, 'k--'); hold off;
        xlim(xLim); ylim(yLim);
        xlabel('Theta phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        ylabel('PC space phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        titleStr = 'PC space phase vs theta phase density map (filtered data)';
        title(titleStr);
        set(fH6, 'Name',titleStr);
        titleStr = strrep(titleStr, ' ', '_');
        titleStr = strrep(titleStr, '(', '_');
        titleStr = strrep(titleStr, ')', '_');
        figName = fullfile(fullSessionFolder, titleStr);
        savefig(fH6, figName,'compact');
        title('');
        saveas(fH6, figName, 'png');
        saveas(fH6, figName, 'pdf');
        close([fH1 fH2 fH3 fH4 fH5 fH6]);
      else
        close([fH1 fH4]);
      end
    end
  end
end
end


function combine3PCsColour(dataFiles, pcaOut, pcaOutFilt, thetaPhaseData, ...
  intervals, figFolder)
% combine3PCsColour(dataFiles, pcaOut, pcaOutFilt, thetaPhaseData, ...
%   intervals, figFolder)
%
% Function produces and saves graphs showing 3 PCs drawn against each other
% with data points colour-coding the theta phase.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   thetaPhaseData, intervals, pcaFigFolder
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  thetaPhaseData
  intervals
  figFolder
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session}) && ...
        size(pcaOut{animal}{session}.score, 2) > 2

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Get phase labels
      [~, thetaPhaseInds] = selectArrayValues( ...
        thetaPhaseData{animal}{session}.timestamps, intervals{animal}{session});
      thetaPhaseSession = thetaPhaseData{animal}{session}.data(thetaPhaseInds);
      %colours = hsv;
      %colourInds = ceil(((thetaPhaseSession + pi)./(2*pi)).*size(colours,1));
      %dataColours = colours(colourInds,:);
      %colormap(hsv);
      if numel(thetaPhaseSession) ~= numel(pcaOut{animal}{session}.score(:,1))
        phaseTimestamps = thetaPhaseData{animal}{session}.timestamps(thetaPhaseInds);
        pcaTimestamps = pcaOut{animal}{session}.timestamps;
        [~, phaseInds, pcaInds] = intersect(round(phaseTimestamps.*1000), ...
          round(pcaTimestamps.*1000));
        thetaPhaseSession = thetaPhaseSession(phaseInds);
      else
        pcaInds = 1:numel(pcaOut{animal}{session}.score(:,1));
      end

      % Convolved data only figure
      fH1 = figure;
      scatter3(pcaOut{animal}{session}.score(pcaInds,1), pcaOut{animal}{session}.score(pcaInds,2), ...
        pcaOut{animal}{session}.score(pcaInds,3), ones(numel(thetaPhaseSession),1), thetaPhaseSession, '.');
      colormap(hsv);
      c = colorbar;
      c.Label.String = 'Theta phase (rad)';
      c.Label.FontSize = fontSize;
      c.Label.FontWeight = 'bold';
      xlabel('PC1 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('PC2 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      zlabel('PC3 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 3 PCs with phase (convolved-only data)';
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');

      % Convolved and band-pass filtered data figure
      fH2 = figure;
      scatter3(pcaOutFilt{animal}{session}.score(pcaInds,1), pcaOutFilt{animal}{session}.score(pcaInds,2), ...
        pcaOutFilt{animal}{session}.score(pcaInds,3), ones(numel(thetaPhaseSession),1), thetaPhaseSession, '.');
      colormap(hsv);
      c = colorbar;
      c.Label.String = 'Theta phase (rad)';
      c.Label.FontSize = fontSize;
      c.Label.FontWeight = 'bold';
      xlabel('PC1 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('PC2 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      zlabel('PC3 (a.u.)', 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Combined 3 PCs with phase (filtered data)';
      title(titleStr);
      set(fH2, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH2, figName,'compact');
      title('');
      saveas(fH2, figName, 'png');
      saveas(fH2, figName, 'pdf');

      close([fH1 fH2]);
    end
  end
end
end


function density2PCs(dataFiles, pcaOut, pcaOutFilt, populationRates, ...
    populationRatesTimesBins, spikeTimes, unitIDs, intervals, figFolder, ...
    fitCircle)
% density2PCs(dataFiles, pcaOut, pcaOutFilt, populationRates, ...
%   populationRatesTimesBins, spikeTimes, unitIDs, intervals, figFolder, ...
%   fitCircle)
%
% Function produces and saves graphs showing 2 PC density maps.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   populationRates
%   populationRatesTimesBins
%   spikeTimes
%   unitIDs
%   intervals
%   figFolder
%   fitCircle
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  populationRates
  populationRatesTimesBins
  spikeTimes
  unitIDs
  intervals
  figFolder
  fitCircle
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session}) && ~isempty(pcaOut{animal}{session}.score) && ...
        size(pcaOut{animal}{session}.score, 2) > 1

      % Order PCs
      [pcInds, pcIndsFilt] = orderPCs(pcaOut, pcaOutFilt, ...
        populationRates, populationRatesTimesBins, intervals, animal, session);

      try
        if fitCircle
          [~, ~, pcaR, ellipse] = fitPlotEllipse([pcaOut{animal}{session}.score(:,pcInds(1)) ...
            pcaOut{animal}{session}.score(:,pcInds(2))], type='ellipseLinearTrace', draw=false);
          [~, ~, pcaFiltR, ellipseFilt] = fitPlotEllipse([pcaOutFilt{animal}{session}.score(:,pcInds(1)) ...
            pcaOutFilt{animal}{session}.score(:,pcInds(2))], type='ellipseLinearTrace', draw=false);
          %[radius, xCentreCoord, yCentreCoord] = circfit( ...
          %  pcaOut{animal}{session}.score(:,pcInds(1)), ...
          %  pcaOut{animal}{session}.score(:,pcInds(2)));
          %pcaR = explainedVarCircleFit(pcaOut{animal}{session}.score(:,pcInds(1)), ...
          %  pcaOut{animal}{session}.score(:,pcInds(1)), xCentreCoord, yCentreCoord, radius);
          %[radiusFilt, xCentreCoordFilt, yCentreCoordFilt] = circfit( ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)));
          %pcaFiltR = explainedVarCircleFit(pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1)), ...
          %  pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)), xCentreCoordFilt, yCentreCoordFilt, radiusFilt);
        end
      catch
        continue
      end

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Bin the data
      nbins=[50 50];
      [pcaMatrix, pcaGridCoords] = hist3([ ...
        pcaOut{animal}{session}.score(:,pcInds(2)), ...
        pcaOut{animal}{session}.score(:,pcInds(1))], nbins);
      [pcaFiltMatrix, pcaFiltGridCoords] = hist3([ ...
        pcaOutFilt{animal}{session}.score(:,pcIndsFilt(2)), ...
        pcaOutFilt{animal}{session}.score(:,pcIndsFilt(1))], nbins);

      % Convolved data only figure
      fH1 = figure;
      s = pcolor(pcaGridCoords{1}, pcaGridCoords{2}, pcaMatrix');
      s.FaceColor = 'interp';
      s.EdgeColor = 'none';
      c = colorbar;
      c.Label.String = 'Count (samples)';
      c.Label.FontSize = fontSize;
      c.Label.FontWeight = 'bold';
      if fitCircle % Plot circle fit
        hold on; plot(ellipse(2,:),ellipse(1,:),'k'); hold off;
        %rectangle('Curvature',[1 1],'Position',[xCentreCoord-radius ...
        %  yCentreCoord-radius 2*radius 2*radius]);
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['R^2=' num2str(pcaR)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
      end
      view([90 -90])
      ylabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      xlabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Density map 2 PCs (convolved-only data)';
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');
      
      % Convolved and band-pass filtered data figure
      fH2 = figure;
      s = pcolor(pcaFiltGridCoords{1}, pcaFiltGridCoords{2}, pcaFiltMatrix');
      s.FaceColor = 'interp';
      s.EdgeColor = 'none';
      c = colorbar;
      c.Label.String = 'Count (samples)';
      c.Label.FontSize = fontSize;
      c.Label.FontWeight = 'bold';
      if fitCircle % Plot circle fit
        hold on; plot(ellipseFilt(2,:),ellipseFilt(1,:),'k'); hold off;
        %rectangle('Curvature',[1 1],'Position',[xCentreCoordFilt-radiusFilt ...
        %  yCentreCoordFilt-radiusFilt 2*radiusFilt 2*radiusFilt]);
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['R^2=' num2str(pcaFiltR)];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);
      end
      view([90 -90])
      ylabel(['PC' num2str(pcIndsFilt(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      xlabel(['PC' num2str(pcIndsFilt(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = 'Density map 2 PCs (filtered data)';
      title(titleStr);
      set(fH2, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH2, figName,'compact');
      title('');
      saveas(fH2, figName, 'png');
      saveas(fH2, figName, 'pdf');

      close([fH1 fH2]);

      % Individual unit figures
      samplingInterval = populationRatesTimesBins{animal}{session}(2) - ...
        populationRatesTimesBins{animal}{session}(1);
      nbins=[50 50];
      smoothingWindow = 9;
      for u = 1:numel(spikeTimes{animal}{session})
        % Bin the unit data
        unitSpikeTimes = selectArrayValues(spikeTimes{animal}{session}{u}, ...
          intervals{animal}{session});
        spikeInds = round(unitSpikeTimes./samplingInterval);
        pcaInds = round(pcaOut{animal}{session}.timestamps./samplingInterval);
        selectSpikeInds = ismember(pcaInds, spikeInds);
        [pcaMatrixUnit, pcaGridCoords] = hist3([ ...
          pcaOut{animal}{session}.score(selectSpikeInds,pcInds(2)), ...
          pcaOut{animal}{session}.score(selectSpikeInds,pcInds(1))], nbins);
        pcaMatrixUnit = imgaussfilt(pcaMatrixUnit,'FilterSize',smoothingWindow);
        [pcaFiltMatrixUnit, pcaFiltGridCoords] = hist3([ ...
          pcaOutFilt{animal}{session}.score(selectSpikeInds,pcIndsFilt(2)), ...
          pcaOutFilt{animal}{session}.score(selectSpikeInds,pcIndsFilt(1))], nbins);
        pcaFiltMatrixUnit = imgaussfilt(pcaFiltMatrixUnit,'FilterSize',smoothingWindow);

        % Regular figure for convolved data
        fH3 = figure;
        s = pcolor(pcaGridCoords{1}, pcaGridCoords{2}, pcaMatrixUnit');
        s.FaceColor = 'interp';
        s.EdgeColor = 'none';
        c = colorbar;
        c.Label.String = 'Count (samples)';
        c.Label.FontSize = fontSize;
        c.Label.FontWeight = 'bold';
        view([90 -90])
        ylabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
        xlabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
        titleStr = ['Unit ' num2str(unitIDs{animal}{session}(u)) ...
          ' density map 2 PCs (convolved-only data)'];
        title(titleStr);
        set(fH3, 'Name',titleStr);
        titleStr = strrep(titleStr, ' ', '_');
        titleStr = strrep(titleStr, '(', '_');
        titleStr = strrep(titleStr, ')', '_');
        figName = fullfile(fullSessionFolder, titleStr);
        savefig(fH3, figName,'compact');
        title('');
        saveas(fH3, figName, 'png');
        saveas(fH3, figName, 'pdf');

        % Normalised figure for convolved data
%         pcaMatrixNormalised = (pcaMatrixUnit./sum(sum(pcaMatrixUnit)))./ ...
%           (imgaussfilt(pcaMatrix,'FilterSize',smoothingWindow)./sum(sum(pcaMatrix)));
%         pcaMatrixNormalised(isnan(pcaMatrixNormalised)) = 0;
%         pcaMatrixNormalised(isinf(pcaMatrixNormalised)) = 0;
%         pcaMatrixNormalised(pcaMatrixNormalised >= sum(sum(pcaMatrixUnit))/sum(sum(pcaMatrix))) = 0;
%         pcaMatrixNormalised = imgaussfilt(pcaMatrixNormalised,'FilterSize',smoothingWindow);
%         fH4 = figure;
%         s = pcolor(pcaGridCoords{1}, pcaGridCoords{2}, pcaMatrixNormalised');
%         s.FaceColor = 'interp';
%         s.EdgeColor = 'none';
%         c = colorbar;
%         c.Label.String = 'Count (samples)';
%         c.Label.FontSize = fontSize;
%         c.Label.FontWeight = 'bold';
%         view([90 -90])
%         ylabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
%         xlabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
%         titleStr = ['Unit ' num2str(unitIDs{animal}{session}(u)) ...
%          ' normalised density map 2 PCs (convolved-only data)'];
%         title(titleStr);
%         set(fH4, 'Name',titleStr);
%         titleStr = strrep(titleStr, ' ', '_');
%         titleStr = strrep(titleStr, '(', '_');
%         titleStr = strrep(titleStr, ')', '_');
%         figName = fullfile(fullSessionFolder, titleStr);
%         savefig(fH4, figName,'compact');
%         title('');
%         saveas(fH4, figName, 'png');
%         saveas(fH4, figName, 'pdf');

        % Regular figure for filtered data
        fH5 = figure;
        s = pcolor(pcaFiltGridCoords{1}, pcaFiltGridCoords{2}, pcaFiltMatrixUnit');
        s.FaceColor = 'interp';
        s.EdgeColor = 'none';
        c = colorbar;
        c.Label.String = 'Count (samples)';
        c.Label.FontSize = fontSize;
        c.Label.FontWeight = 'bold';
        view([90 -90])
        ylabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
        xlabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
        titleStr = ['Unit ' num2str(unitIDs{animal}{session}(u)) ...
          ' density map 2 PCs (filtered data)'];
        title(titleStr);
        set(fH5, 'Name',titleStr);
        titleStr = strrep(titleStr, ' ', '_');
        titleStr = strrep(titleStr, '(', '_');
        titleStr = strrep(titleStr, ')', '_');
        figName = fullfile(fullSessionFolder, titleStr);
        savefig(fH5, figName,'compact');
        title('');
        saveas(fH5, figName, 'png');
        saveas(fH5, figName, 'pdf');

        % Normalised figure for filtered data
%         pcaMatrixNormalised = (pcaFiltMatrixUnit./sum(sum(pcaFiltMatrixUnit)))./ ...
%           (imgaussfilt(pcaFiltMatrix,'FilterSize',smoothingWindow)./sum(sum(pcaFiltMatrix)));
%         pcaMatrixNormalised(isnan(pcaMatrixNormalised)) = 0;
%         pcaMatrixNormalised(isinf(pcaMatrixNormalised)) = 0;
%         %pcaMatrixNormalised(pcaMatrixNormalised > 0.15) = 0.15;
%         pcaMatrixNormalised = imgaussfilt(pcaMatrixNormalised,'FilterSize',smoothingWindow);
%         fH6 = figure;
%         s = pcolor(pcaFiltGridCoords{1}, pcaFiltGridCoords{2}, pcaMatrixNormalised');
%         s.FaceColor = 'interp';
%         s.EdgeColor = 'none';
%         c = colorbar;
%         c.Label.String = 'Count (samples)';
%         c.Label.FontSize = fontSize;
%         c.Label.FontWeight = 'bold';
%         view([90 -90])
%         ylabel(['PC' num2str(pcInds(1)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
%         xlabel(['PC' num2str(pcInds(2)) ' (a.u.)'], 'FontSize',fontSize, 'FontWeight','bold');
%         titleStr = ['Unit ' num2str(unitIDs{animal}{session}(u)) ...
%          ' normalised density map 2 PCs (filtered data)'];
%         title(titleStr);
%         set(fH6, 'Name',titleStr);
%         titleStr = strrep(titleStr, ' ', '_');
%         titleStr = strrep(titleStr, '(', '_');
%         titleStr = strrep(titleStr, ')', '_');
%         figName = fullfile(fullSessionFolder, titleStr);
%         savefig(fH6, figName,'compact');
%         title('');
%         saveas(fH6, figName, 'png');
%         saveas(fH6, figName, 'pdf');
%         close([fH3 fH4 fH5 fH6]);
        close([fH3 fH5]);
      end
    end
  end
end
end


function splitPhaseDistros(dataFiles, pcaOut, pcaOutFilt, ...
  fullInterpCoherence, populationRates, populationRatesTimesBins, ...
  thetaFrequencyData, intervals, nPhaseHistogramBins, ...
  figFolder)
% splitPhaseDistros(dataFiles, pcaOut, pcaOutFilt, ...
%   fullInterpCoherence, populationRates, populationRatesTimesBins, ...
%   thetaFrequencyData, intervals, nPhaseHistogramBins, ...
%   figFolder)
%
% Function produces and saves unit phase distributions split on the basis
% of their PC loadings.
%
% Args:
%   dataFiles
%   pcaOut
%   pcaOutFilt
%   fullInterpCoherence
%   populationRates
%   populationRatesTimesBins
%   thetaFrequencyData
%   intervals
%   nPhaseHistogramBins
%   figFolder
%
% Returns:
%   None.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  fullInterpCoherence
  populationRates
  populationRatesTimesBins
  thetaFrequencyData
  intervals
  nPhaseHistogramBins
  figFolder
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session}) && ...
        ~isempty(fullInterpCoherence{animal}{session}) && ...
        ~isempty(fullInterpCoherence{animal}{session}.frequency)

      % Order PCs
      [pcInds, pcIndsFilt] = orderPCs(pcaOut, pcaOutFilt, ...
        populationRates, populationRatesTimesBins, intervals, animal, session);

      % Create the figure folder
      [~, animalFolder] = fileparts(fileparts(fileparts(dataFiles{animal}{session})));
      [~, sessionFolder] = fileparts(fileparts(dataFiles{animal}{session}));
      fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
      if ~exist(fullSessionFolder, 'dir')
        mkdir(fullSessionFolder);
      end

      % Find the mean instantaneous theta frequency
      [~, instantThetaFrequencyInds] = selectArrayValues( ...
        thetaFrequencyData{animal}{session}.timestamps, intervals{animal}{session});
      instantThetaFrequency = thetaFrequencyData{animal}{session}.data(instantThetaFrequencyInds);
      meanInstantThetaFrequency = mean(instantThetaFrequency, 'omitnan');
      [~,fIndMeanInstant] = min(abs(fullInterpCoherence{animal}{session}.frequency(1,:) - ...
        meanInstantThetaFrequency));

      % Split units into PC1 and PC2 groups
      unitsPC1 = pcaOut{animal}{session}.coeff(:,pcInds(1)) > pcaOut{animal}{session}.coeff(:,pcInds(2));
      unitsPC2 = ~unitsPC1;
      unitsPC1Filt = pcaOutFilt{animal}{session}.coeff(:,pcIndsFilt(1)) > pcaOutFilt{animal}{session}.coeff(:,pcIndsFilt(2));
      unitsPC2Filt = ~unitsPC1Filt;

      % Draw the histogram: Convolved-only data
      phase = fullInterpCoherence{animal}{session}.phase(fullInterpCoherence{animal}{session}.validUnits,fIndMeanInstant);
      phasePC1 = phase(unitsPC1);
      [phaseHistogramsPC1, phaseHistogramBinsPC1, ~, meanPhasePC1] = ...
        phaseHistrogram(phasePC1, centre=0, nBins=nPhaseHistogramBins);
      phasePC2 = phase(unitsPC2);
      [phaseHistogramsPC2, phaseHistogramBinsPC2, ~, meanPhasePC2] = ...
        phaseHistrogram(phasePC2, centre=0, nBins=nPhaseHistogramBins);
      fH1 = figure;
      p1 = plot(-phaseHistogramBinsPC1, ...
        phaseHistogramsPC1, 'g', 'LineWidth',2); hold on
      p2 = plot(-phaseHistogramBinsPC2, ...
        phaseHistogramsPC2, 'm', 'LineWidth',2);
      yLim = ylim;
      plot([-meanPhasePC1 -meanPhasePC1], yLim, 'g:');
      plot([-meanPhasePC2 -meanPhasePC2], yLim, 'm:'); hold off
      legend([p1, p2], {['PC' num2str(pcInds(1)) '>PC' num2str(pcInds(2))], ...
        ['PC' num2str(pcInds(1)) '<PC' num2str(pcInds(2))]});
      legend('boxoff');
      xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('Unit count', 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = ['Unit phase distros PC' num2str(pcInds(1)) '>PC' num2str(pcInds(2)) ...
        ' v PC' num2str(pcInds(1)) '<PC' num2str(pcInds(2)) ' (convolved-only data)'];
      title(titleStr);
      set(fH1, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      titleStr = strrep(titleStr, '>', '_');
      titleStr = strrep(titleStr, '<', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH1, figName,'compact');
      title('');
      saveas(fH1, figName, 'png');
      saveas(fH1, figName, 'pdf');

      % Draw the histogram: Filtered data
      phasePC1 = phase(unitsPC1Filt);
      [phaseHistogramsPC1, phaseHistogramBinsPC1, ~, meanPhasePC1] = ...
        phaseHistrogram(phasePC1, centre=0, nBins=nPhaseHistogramBins);
      phasePC2 = phase(unitsPC2Filt);
      [phaseHistogramsPC2, phaseHistogramBinsPC2, ~, meanPhasePC2] = ...
        phaseHistrogram(phasePC2, centre=0, nBins=nPhaseHistogramBins);
      fH2 = figure;
      p1 = plot(-phaseHistogramBinsPC1, ...
        phaseHistogramsPC1, 'g', 'LineWidth',2); hold on
      p2 = plot(-phaseHistogramBinsPC2, ...
        phaseHistogramsPC2, 'm', 'LineWidth',2);
      yLim = ylim;
      plot([-meanPhasePC1 -meanPhasePC1], yLim, 'g:');
      plot([-meanPhasePC2 -meanPhasePC2], yLim, 'm:'); hold off
      legend([p1, p2], {['PC' num2str(pcIndsFilt(1)) '>PC' num2str(pcIndsFilt(2))], ...
        ['PC' num2str(pcIndsFilt(1)) '<PC' num2str(pcIndsFilt(2))]});
      legend('boxoff');
      xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('Unit count', 'FontSize',fontSize, 'FontWeight','bold');
      titleStr = ['Unit phase distros PC' num2str(pcIndsFilt(1)) '>PC' num2str(pcIndsFilt(2)) ...
        ' v PC' num2str(pcIndsFilt(1)) '<PC' num2str(pcIndsFilt(2)) ' (filtered data)'];
      title(titleStr);
      set(fH2, 'Name',titleStr);
      titleStr = strrep(titleStr, ' ', '_');
      titleStr = strrep(titleStr, '(', '_');
      titleStr = strrep(titleStr, ')', '_');
      titleStr = strrep(titleStr, '>', '_');
      titleStr = strrep(titleStr, '<', '_');
      figName = fullfile(fullSessionFolder, titleStr);
      savefig(fH2, figName,'compact');
      title('');
      saveas(fH2, figName, 'png');
      saveas(fH2, figName, 'pdf');

      close([fH1 fH2]);
    end
  end
end
end


function pcCorrFigWrapper(dataFiles, pcaOut, pcaOutFilt, scoreCorr, nPCs, figFolder)
% pcCorrFigWrapper(dataFiles, pcaOut, pcaOutFilt, scoreCorr, nPC, figFolder)

arguments
  dataFiles
  pcaOut
  pcaOutFilt
  scoreCorr
  nPCs
  figFolder
end

% Drawing parameters
fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(pcaOut{animal}{session})
      for iPC = 1:nPCs

        % Create the figure folder
        dataFile = dataFiles{animal}{session};
        [~, animalFolder] = fileparts(fileparts(fileparts(dataFile)));
        [~, sessionFolder] = fileparts(fileparts(dataFile));
        fullSessionFolder = fullfile(figFolder, animalFolder, sessionFolder);
        if ~exist(fullSessionFolder, 'dir')
          mkdir(fullSessionFolder);
        end

        % Draw the figure
        fH = figure;
        plot(pcaOut{animal}{session}.score(:,iPC), ...
          pcaOutFilt{animal}{session}.score(:,iPC), '.', 'MarkerSize',1);

        % Text
        yLim = ylim;
        yDist = yLim(2) - yLim(1);
        xLim = xlim;
        xDist = xLim(2) - xLim(1);
        textStr = ['r=' num2str(scoreCorr{animal}{session}.rho(iPC)) ...
          ' p=' num2str(scoreCorr{animal}{session}.p(iPC))];
        text(xLim(1)+0.025*xDist, yLim(1)+0.95*yDist, textStr);

        % Line fit
        xAxisStep = xDist/10000;
        x = xLim(1):xAxisStep:xLim(2);
        [~, slope, coefficients] = fitLine(pcaOut{animal}{session}.score(:,iPC), ...
          pcaOutFilt{animal}{session}.score(:,iPC), type='linear-linear');
        yFit = x.*slope + coefficients(2);
        hold on; plot(x, yFit, 'k--'); hold off;
        ylim(yLim);
        xlim(xLim);

        % Labels
        xlabel(['Convolved PC' num2str(iPC) ' score (a.u.)'], ...
          'FontSize',fontSize, 'FontWeight','bold');
        ylabel(['Band-pass filtered PC' num2str(iPC) ' score (a.u.)'], ...
          'FontSize',fontSize, 'FontWeight','bold');
        figName = ['Correlation between filtered and nonfiltered PC' num2str(iPC)];
        title(figName, 'Interpreter','none');
        set(fH, 'Name',figName);

        % Save the figure
        if ~exist(fullSessionFolder, 'dir')
          mdir(fullSessionFolder);
        end
        figName = strrep(figName, ' ', '_');
        figName = strrep(figName, '(', '_');
        figName = strrep(figName, ')', '_');
        figName = fullfile(fullSessionFolder, figName);
        savefig(fH, figName,'compact');
        title('');
        saveas(fH, figName, 'png');
        saveas(fH, figName, 'pdf');
        close(fH);
      end
    end
  end
end
end


function [regularPCOrder, filteredPCOrder] = orderPCs(pcaOut, pcaOutFilt, ...
  populationRates, populationRatesTimesBins, intervals, animal, session)
% [regularPCOrder, filteredPCOrder] = orderPCs(pcaOut, pcaOutFilt, ...
%   populationRates, populationRatesTimesBins, intervals, animal, session)

arguments
  pcaOut
  pcaOutFilt
  populationRates
  populationRatesTimesBins
  intervals
  animal
  session
end

% Correlate PCs
nPCs = 3;
[~, populationRateInds] = selectArrayValues( ...
  populationRatesTimesBins{animal}{session}, intervals{animal}{session});
r = zeros(1,nPCs);
rFilt = zeros(1,nPCs);
for pc = 1:nPCs
  r(pc) = corrLinearCircular(pcaOut{animal}{session}.score(:,pc), ...
    populationRates{animal}{session}(populationRateInds), type='Spearman');
  rFilt(pc) = corrLinearCircular(pcaOutFilt{animal}{session}.score(:,pc), ...
    populationRates{animal}{session}(populationRateInds), type='Spearman');
end

% Order regular PCs
if r(3) < 0
  [~, maxSecondaryPC] = max(abs(r(2:3)));
else
  maxSecondaryPC = 1;
end
regularPCOrder = [1 maxSecondaryPC+1];
[~, pcOrder] = sort(r(regularPCOrder), 'descend');
regularPCOrder = regularPCOrder(pcOrder);

% Order filtered PCs
if rFilt(3) < 0
  [~, maxSecondaryPC] = max(abs(rFilt(2:3)));
else
  maxSecondaryPC = 1;
end
filteredPCOrder = [1 maxSecondaryPC+1];
[~, pcOrder] = sort(rFilt(filteredPCOrder), 'descend');
filteredPCOrder = filteredPCOrder(pcOrder);
end