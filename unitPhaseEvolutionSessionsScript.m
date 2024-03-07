% This script analyses firing patterns of individual units with respect
% to the MS theta activity within individual recording sessions

%% Initialise parameters
% Data setup parameters
params

% Analysis parameters
samplingInterval = 0.002;
minIntervalLength = 1;
nPhaseHistogramBins = 9;
filterType = 'wideband';

% Unit quality criteria
firingRateCutoff = 200; % At least that many spikes per hour or at least a single 20-minute window (i.e., firingRateCutoff/3)
refractoryContaminationTh = 5; % Not more than this percentage of spiking ACG refractory period contamination.
                               % Will take the ratio mean(+/- 1 ms relative to 0 lag of ACG) /
                               % mean(1000-1500ms ACG away from zero lag (ACG shoulder)).

% Script run parameters
loadPhases = false;
runThetaMazeAnalyses_speed = false;
runThetaMazeAnalyses = false;
runThetaMazeAnalyses_trials = false;
runSleepAnalyses = false;
runCombinedAnalyses = true;


%% Estimate theta phase of unit spikes
if loadPhases
  unitPhases = getUnitThetaPhases(dataFiles, filterType);
end


%% Plot unit phase histograms and phase evolution
if runThetaMazeAnalyses_speed
  % Define periods of interest: ThetaMaze_AlternativeRunning_trials_highSpeed
  epochsOfInterest = 'ThetaMaze_AlternativeRunning'; %#ok<*UNRCH> 
  thetaPower = '';
  onlyTrials = true;
  onlyHighSpeed = true;

  % Find the time intervals of interest
  intervals = cell(numel(dataFiles),1);
  for animal = 1:numel(dataFiles)
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate unit firing rates
  firingRates = getFiringRates(dataFiles, intervals);

  % plot unit phase
  figSubfolder = [epochsOfInterest '_trials_highSpeed'];
  [phaseSortingIndsSpeed, phaseHistosSpeed] = plotUnitPhase(dataFiles, ...
    unitPhases, firingRates, firingRateCutoff, refractoryContaminationTh, ...
    intervals, channelOrder, nPhaseHistogramBins, figSubfolder, ...
    unitPhaseFigFolder);


  % Define periods of interest: ThetaMaze_AlternativeRunning_highSpeed
  onlyTrials = false;
  onlyHighSpeed = true;

  % Find the time intervals of interest
  intervals = cell(numel(dataFiles),1);
  for animal = 1:numel(dataFiles)
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate unit firing rates
  firingRates = getFiringRates(dataFiles, intervals);

  % plot unit phase
  figSubfolder = [epochsOfInterest '_highSpeed'];
  plotUnitPhase(dataFiles, unitPhases, firingRates, firingRateCutoff, ...
    refractoryContaminationTh, intervals, channelOrder, ...
    nPhaseHistogramBins, figSubfolder, unitPhaseFigFolder, ...
    secondarySorting=phaseSortingIndsSpeed);
end


if runThetaMazeAnalyses
  % Define periods of interest: ThetaMaze_AlternativeRunning
  onlyTrials = false;
  onlyHighSpeed = false;

  % Find the time intervals of interest
  intervals = cell(numel(dataFiles),1);
  for animal = 1:numel(dataFiles)
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate unit firing rates
  firingRates = getFiringRates(dataFiles, intervals);

  % plot unit phase
  [~, phaseHistosBehav] = plotUnitPhase(dataFiles, unitPhases, firingRates, ...
    firingRateCutoff, refractoryContaminationTh, intervals, channelOrder, ...
    nPhaseHistogramBins, epochsOfInterest, unitPhaseFigFolder, ...
    secondarySorting=phaseSortingIndsSpeed);
end


if runThetaMazeAnalyses_trials
  % Define periods of interest: ThetaMaze_AlternativeRunning_trials
  onlyTrials = true;
  onlyHighSpeed = false;

  % Find the time intervals of interest
  intervals = cell(numel(dataFiles),1);
  for animal = 1:numel(dataFiles)
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate unit firing rates
  firingRates = getFiringRates(dataFiles, intervals);

  % plot unit phase
  figSubfolder = [epochsOfInterest '_trials'];
  plotUnitPhase(dataFiles, unitPhases, firingRates, firingRateCutoff, ...
    refractoryContaminationTh, intervals, channelOrder, ...
    nPhaseHistogramBins, figSubfolder, unitPhaseFigFolder, ...
    secondarySorting=phaseSortingIndsSpeed);
end


if runSleepAnalyses
  % Define periods of interest: REM
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunnig','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem';

  % Find the time intervals of interest
  intervals = cell(numel(dataFiles),1);
  for animal = 1:numel(dataFiles)
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate unit firing rates
  firingRates = getFiringRates(dataFiles, intervals);

  % plot unit phase
  [~, phaseHistosREM] = plotUnitPhase(dataFiles, unitPhases, firingRates, ...
    firingRateCutoff, refractoryContaminationTh, intervals, channelOrder, ...
    nPhaseHistogramBins, sleepState, unitPhaseFigFolder, ...
    secondarySorting=phaseSortingIndsSpeed);


  % Define periods of interest: NREM
  sleepState = 'nrem';

  % Find the time intervals of interest
  intervals = cell(numel(dataFiles),1);
  for animal = 1:numel(dataFiles)
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate unit firing rates
  firingRates = getFiringRates(dataFiles, intervals);

  % plot unit phase
  [~, phaseHistosNREM] = plotUnitPhase(dataFiles, unitPhases, firingRates, ...
    firingRateCutoff, refractoryContaminationTh, intervals, channelOrder, ...
    nPhaseHistogramBins, sleepState, unitPhaseFigFolder, ...
    secondarySorting=phaseSortingIndsSpeed);


%   % Define periods of interest: All sessions when animal is awake in its homecage
%   epochsOfInterest = {'HomeCage_Sleep','HomeCage_Wake'};
%   sleepState = 'wake';
% 
%   % Find the time intervals of interest
%   intervals = cell(numel(dataFiles),1);
%   for animal = 1:numel(dataFiles)
%     for session = 1:numel(dataFiles{animal})
%       intervals{animal}{session} = getTimeIntervals( ...
%         dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
%         minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
%         onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
%     end
%   end
% 
%   % Calculate unit firing rates
%   firingRates = getFiringRates(dataFiles, intervals);
% 
%   % plot unit phase
%   plotUnitPhase(dataFiles, unitPhases, firingRates, firingRateCutoff, ...
%     refractoryContaminationTh, intervals, channelOrder, ...
%     nPhaseHistogramBins, sleepState, unitPhaseFigFolder, ...
%     secondarySorting=phaseSortingIndsBehav);
end


if runCombinedAnalyses
  plotUnitPhaseAcrossSessions(dataFiles, phaseHistosSpeed, phaseHistosBehav, ...
    phaseHistosREM, phaseHistosNREM, fullfile(unitPhaseFigFolder, 'combined'));
end



%% Local functions
function sessionIntervals = getTimeIntervals(dataFile, epochsOfInterest, options)
% sessionIntervals = getTimeIntervals(dataFile, epochsOfInterest, <options>)
%
% Function extracts time intervals of interest for a single recording
% session located in the ms-preliminary repository.
%
% Args:
%   dataFile
%   epochsOfInterest
%   <thetaPower>
%   <minIntervalLength>
%   <onlyTrials>
%   <onlyHighSpeed>
%   <sleepState>
%   <excludeNoise>
%
% Returns:
%   sessionIntervals
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFile
  epochsOfInterest {mustBeCharOrListedType(epochsOfInterest,'cell')}
  options.thetaPower = '';
  options.minIntervalLength = 0;
  options.onlyTrials = false
  options.onlyHighSpeed = false
  options.sleepState {mustBeMember(options.sleepState, {'wake','ma','nrem','rem',''})} = ''
  options.excludeNoise (1,1) {mustBeA(options.excludeNoise,'logical')} = false;
end

% Parse input
if ischar(epochsOfInterest)
  epochsOfInterest = {epochsOfInterest};
end

% Load data
sessionFile = strrep(dataFile, '*', 'session');
if exist(sessionFile,'file')
  load(sessionFile); %#ok<*LOAD>
else
  sessionIntervals = [];
  warning(['Session file ' sessionFile ' does not exist. Terminating function call.'])
  return
end
if options.excludeNoise
  lfpNoiseFile = strrep(dataFile, '*', 'lfpNoise.events');
  if exist(lfpNoiseFile,'file')
    load(lfpNoiseFile);
  else
    options.excludeNoise = false;
    warning(['LFP noise file ' lfpNoiseFile ' does not exist.'])
  end
end
if strcmpi(options.thetaPower,'moderate') || strcmpi(options.thetaPower,'high')
  thetaPowerFile = strrep(dataFile, '*', 'thetaPower.timeseriesCollection');
  if exist(thetaPowerFile,'file')
    load(thetaPowerFile);
  else
    options.thetaPower = '';
    warning(['Theta power file ' thetaPowerFile ' does not exist.'])
  end
end
if options.onlyTrials || options.onlyHighSpeed
  circularTrackFile = strrep(dataFile, '*', 'circular_track.behavior');
  if exist(circularTrackFile,'file')
    load(circularTrackFile);
  else
    options.onlyTrials = false;
    options.onlyHighSpeed = false;
    warning(['Circular track behaviour file ' circularTrackFile ' does not exist. Behaviour data is missing.'])
  end
end
if ~isempty(options.sleepState)
  eegStatesFile = strrep(dataFile, '*', 'eegstates');
  if exist(eegStatesFile,'file')
    load(eegStatesFile);
  else
    warning(['EEG states file ' eegStatesFile ' does not exist. EEG states data is missing.'])
  end
  sleepStateFile = strrep(dataFile, '*', 'SleepState.states');
  if exist(sleepStateFile,'file')
    load(sleepStateFile);
  else
    warning(['Sleep state file ' sleepStateFile ' does not exist. Sleep state data is missing.'])
  end
end

% Proceed with interval selection
sessionIntervals = [];
for epoch = 1:numel(session.epochs)
  if isfield(session.epochs{epoch},'behavioralParadigm') && ...
      ismember(session.epochs{epoch}.behavioralParadigm, epochsOfInterest)

    % Select time intervals corresponding to the sessions of interest
    interval = [session.epochs{epoch}.startTime session.epochs{epoch}.stopTime];

    % Exclude intervals with noisy (saturated) LFPs
    if options.excludeNoise && exist('lfpNoise','var')
      cleanIntervals = invertIntervals(lfpNoise.timestamps, ...
        session.epochs{epoch}.stopTime, startTime=session.epochs{epoch}.startTime);
      if ~isempty(cleanIntervals)
        interval = intervalOverlap(interval, cleanIntervals);
      end
    end

    % Select time intervals corresponding to increased theta power
    if (strcmpi(options.thetaPower,'moderate') || strcmpi(options.thetaPower,'high')) && ~isempty(interval)
      if strcmpi(options.thetaPower,'moderate')
        increasedThetaPowerIntervals = thetaPower.data(:,3)' | thetaPower.data(:,4)';
      elseif strcmpi(options.thetaPower,'high')
        increasedThetaPowerIntervals = logical(thetaPower.data(:,4)');
      end
      increasedThetaPowerIntervals = logical2intervals(increasedThetaPowerIntervals);
      increasedThetaPowerIntervals(:,2) = increasedThetaPowerIntervals(:,2) + 1;
      increasedThetaPowerIntervals(end,2) = min([increasedThetaPowerIntervals(end,2) numel(thetaPower.timestamps)]);
      increasedThetaPowerIntervals = thetaPower.timestamps(increasedThetaPowerIntervals);
      interval = intervalOverlap(interval, increasedThetaPowerIntervals);
    end

    % Select time intervals corresponding to trials only
    if options.onlyTrials && ~isempty(interval)
      if exist('circular_track','var') && isfield(circular_track,'trials')
        trialIntervals = [circular_track.trials.alternation.start ...
          circular_track.trials.alternation.end];
      else
        trialIntervals = [];
      end
      interval = intervalOverlap(interval, trialIntervals);
    end

    % Select time intervals corresponding to high speeds
    if options.onlyHighSpeed && ~isempty(interval)
      if exist('circular_track','var') && isfield(circular_track,'speed')
        highSpeedIntervals = circular_track.speed >= circular_track.speed_th;
        highSpeedIntervals = logical2intervals(highSpeedIntervals);
        highSpeedIntervals(:,2) = highSpeedIntervals(:,2) + 1;
        highSpeedIntervals(end,2) = min([highSpeedIntervals(end,2) numel(circular_track.timestamps)]);
        highSpeedIntervals = circular_track.timestamps(highSpeedIntervals);
      else
        highSpeedIntervals = [];
      end
      interval = intervalOverlap(interval, highSpeedIntervals);
    end

    % Select intervals corresponding to a brain state of interest
    if ~isempty(options.sleepState)
      if exist('SleepState','var') && isfield(SleepState,'ints')
        switch options.sleepState
          case 'wake'
            interval = intervalOverlap(interval, SleepState.ints.WAKEstate);
          case 'ma'
            interval = intervalOverlap(interval, SleepState.ints.MAstate);
          case 'nrem'
            interval = intervalOverlap(interval, SleepState.ints.NREMstate);
          case 'rem'
            interval = intervalOverlap(interval, SleepState.ints.REMstate);
          otherwise
            interval = [];
        end
      else
        interval = [];
      end
    end

    % Concatenate all relevant intervals for the session
    if ~isempty(interval) && sum(interval(:,2) - interval(:,1) >= options.minIntervalLength)
      sessionIntervals = [sessionIntervals; ...
        interval(interval(:,2) - interval(:,1) >= options.minIntervalLength,:)]; %#ok<*AGROW>
    end
  end
end
end


function maxFiringRates = getFiringRates(dataFiles, intervals)
% maxFiringRates = getFiringRates(dataFiles, intervals)
%
% Function calculates unit firing rates over sliding windows of 20 minutes
% and outputs window values with the largest estimates for all units in the
% ms-preliminary dataset.
%
% Args:
%   dataFiles
%   intervals
%
% Returns:
%   maxFiringRates
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  intervals
end

maxFiringRates = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})

    % Load unit spiking data
    spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
    if ~exist(spikesFile, 'file')
      continue
    end
    load(spikesFile); %#ok<*LOAD>

    % Loop through units
    maxFiringRates{animal}{session} = zeros(spikes.numcells,1);
    if isempty(intervals{animal}{session})
      continue
    end
    for unit = 1:spikes.numcells

      % Exclude spike times outside time intervals of interest
      spikeTimes = selectArrayValues( ...
        spikes.times{unit}', intervals{animal}{session}([1 end]));

      % Calculate maximal firing rates
      if isempty(spikeTimes)
        maxFiringRates{animal}{session}(unit) = 0;
      else
        maxFiringRates{animal}{session}(unit) = firingRateWindows( ...
          spikeTimes, stepSize=30, startTime=intervals{animal}{session}(1,1));
      end
    end
  end
end
end


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


function [phaseSortingInds, phaseHistos] = plotUnitPhase(dataFiles, ...
  unitPhases, firingRates, firingRateCutoff, refractoryContaminationTh, ...
  intervals, channelOrder, nPhaseHistogramBins, subFolderName, ...
  outputFolder, options)
% [phaseSortingInds, phaseHistos] = plotUnitPhase(dataFiles, ...
%   unitPhases, firingRates, firingRateCutoff, refractoryContaminationTh, ...
%   intervals, channelOrder, nPhaseHistogramBins, subFolderName, ...
%   outputFolder, <secondarySorting>)
%
% Function displays unit spiking phase histograms and spiking phase
% evolution graphs.
%
% Args:
%   dataFiles
%   unitPhases
%   firingRates
%   firingRateCutoff
%   refractoryContaminationTh
%   intervals
%   channelOrder
%   nPhaseHistogramBins
%   subFolderName
%   outputFolder
%   <secondarySorting>
%
% Returns:
%   phaseSortingInds
%   phaseHistos
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  unitPhases
  firingRates
  firingRateCutoff
  refractoryContaminationTh
  intervals
  channelOrder
  nPhaseHistogramBins
  subFolderName
  outputFolder
  options.secondarySorting = {};
end

phaseSortingInds = cell(numel(dataFiles),1);
phaseHistos.histograms = cell(numel(dataFiles),1);
phaseHistos.bins = cell(numel(dataFiles),1);
phaseHistos.unitIDs = cell(numel(dataFiles),1);
phaseHistos.channelIDs = cell(numel(dataFiles),1);
phaseHistos.validUnits = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    % Load unit spiking data
    spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
    if ~exist(spikesFile, 'file')
      phaseSortingInds{animal}{session} = [];
      phaseHistos.histograms{animal}{session} = [];
      phaseHistos.bins{animal}{session} = [];
      phaseHistos.unitIDs{animal}{session} = [];
      phaseHistos.channelIDs{animal}{session} = [];
      phaseHistos.validUnits{animal}{session} = [];
      continue
    end
    load(spikesFile); %#ok<*LOAD>

    % Load unit refractory period contamination percentages
    refractoryContaminationFile = strrep(dataFiles{animal}{session}, '*', 'contaminationPercent.cellinfo');
    load(refractoryContaminationFile);

    % Eliminate subquality units
    includeUnits = firingRates{animal}{session} >= firingRateCutoff/3600;
    includeUnits = includeUnits & contaminationPercent.data <= refractoryContaminationTh;
    phaseHistos.validUnits{animal}{session} = includeUnits;

    % Loop through units
    phaseHistograms = zeros(spikes.numcells,nPhaseHistogramBins);
    phaseHistogramBins = zeros(spikes.numcells,nPhaseHistogramBins);
    phaseMeans = NaN(spikes.numcells,1);
    phaseSDs = NaN(spikes.numcells,1);
    for unit = 1:spikes.numcells
      if includeUnits(unit)

        % Plot unit phase histogram
        [selectSpikeTimes, selectSpikeInds] = selectArrayValues( ...
          spikes.times{unit}, intervals{animal}{session}); %#ok<*ASGLU> 
        unitPhaseVec = unitPhases{animal}{session}{unit}(selectSpikeInds);
        if ~isempty(unitPhaseVec)
          [phaseHistograms(unit,:), phaseHistogramBins(unit,:), ~, phaseMeans(unit), ...
            phaseSDs(unit)] = phaseHistrogram(unitPhaseVec(:), centre=0, ...
            nBins=nPhaseHistogramBins);

          histTitle = ['uid ' num2str(spikes.cluID(unit)) ' ch ' ...
            num2str(channelOrder{animal}(session,spikes.maxWaveformCh1(unit))) ...
            ' spiking phase histogram'];
          figPath = fullfile(outputFolder,subFolderName,spikes.basename);
          fH1 = phaseHistogramPlot(phaseHistograms(unit,:)', binLocs=phaseHistogramBins(unit,:), ...
            dataMean=phaseMeans(unit), figTitle=histTitle, figPath=figPath);

          % Plot unit phase evolution
          fH2 = figure; plot(selectSpikeTimes, unitPhaseVec, '.', 'MarkerSize',2)
          ylim([-pi pi]);
          xlabel('Time (s)');
          ylabel('Phase wrt population rate (rad)');
          titleStr = ['uid ' num2str(spikes.cluID(unit)) ' ch ' ...
            num2str(channelOrder{animal}(session,spikes.maxWaveformCh1(unit))) ...
            ' spiking phase evolution'];
          title(titleStr, 'Interpreter','none');
          set(fH2, 'Name',titleStr);
          figFilename = strrep(titleStr, ' ', '_');
          figFilename = fullfile(figPath, [figFilename '.fig']);
          savefig(fH2, figFilename, 'compact');

          close([fH1 fH2])
        end
      end
    end
    phaseHistos.histograms{animal}{session} = phaseHistograms;
    phaseHistos.bins{animal}{session} = phaseHistogramBins;
    phaseHistos.unitIDs{animal}{session} = spikes.cluID;
    phaseHistos.channelIDs{animal}{session} = channelOrder{animal}(session,spikes.maxWaveformCh1);
    [sortedMeans, sortedInds] = sort(phaseMeans);
    sortedSDs = phaseSDs(sortedInds);
    phaseSortingInds{animal}{session} = sortedInds;

    % Plot all unit phase histograms in a single figure
    if any(includeUnits)
      fontSize = 14;
      spikeCounts = sum(phaseHistograms(includeUnits,:)'); %#ok<UDIM> 
      
      % Non-normalised figures
      fH1 = figure;
      plot(phaseHistogramBins(includeUnits,:)', phaseHistograms(includeUnits,:)');

      % Label axes
      xTickLocs = [-pi, -pi/2, 0, pi/2, pi];
      xTickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
      xticks(xTickLocs);
      xticklabels(xTickLabels);
      xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('Spike count', 'FontSize',fontSize, 'FontWeight','bold');

      % Add the figure title
      tH = title('Phase Histograms', 'Interpreter','none');

      % Save the figure
      if ~exist(figPath,'dir')
        mkdir(figPath);
      end
      filename = strrep(tH.String,' ','_');
      filename = strrep(filename, ',','_');
      filename = strrep(filename, '-','_');
      filename = strrep(filename, ':','' );
      filename = strrep(filename, '^','' );
      filename = strrep(filename, '{','' );
      filename = strrep(filename, '}','' );
      filename = strrep(filename, '.','p');
      filename = [figPath filesep filename '.fig']; %#ok<*AGROW>
      savefig(fH1,filename,'compact');
      title('');
      saveas(fH1,filename(1:end-4),'png');
      saveas(fH1,filename(1:end-4),'pdf');
      close(fH1);

      % Normalised figures
      fH2 = figure;
      plot(phaseHistogramBins(includeUnits,:)', phaseHistograms(includeUnits,:)'./spikeCounts);

      % Label axes
      xTickLocs = [-pi, -pi/2, 0, pi/2, pi];
      xTickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
      xticks(xTickLocs);
      xticklabels(xTickLabels);
      xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
      ylabel('Normalised spike count', 'FontSize',fontSize, 'FontWeight','bold');

      % Add the figure title
      tH = title('Normalised Phase Histograms', 'Interpreter','none');

      % Save the figure
      if ~exist(figPath,'dir')
        mkdir(figPath);
      end
      filename = strrep(tH.String,' ','_');
      filename = strrep(filename, ',','_');
      filename = strrep(filename, '-','_');
      filename = strrep(filename, ':','' );
      filename = strrep(filename, '^','' );
      filename = strrep(filename, '{','' );
      filename = strrep(filename, '}','' );
      filename = strrep(filename, '.','p');
      filename = [figPath filesep filename '.fig']; %#ok<*AGROW>
      savefig(fH2,filename,'compact');
      title('');
      saveas(fH2,filename(1:end-4),'png');
      saveas(fH2,filename(1:end-4),'pdf');
      close(fH2);

      % Phase distribution SDs
      fH3 = figure; hold on
      for iMean = 1:numel(sortedMeans)
        if ~isnan(sortedMeans(iMean))
          plot(sortedMeans(iMean),-iMean, 'b.', 'MarkerSize',10);
          plot([sortedMeans(iMean)-sortedSDs(iMean)/2 ...
            sortedMeans(iMean)+sortedSDs(iMean)/2],[-iMean -iMean], 'b', ...
            'LineWidth',1);
        else
          break
        end
      end

      % Label axes
      xTickLocs = [-pi, -pi/2, 0, pi/2, pi];
      xTickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
      xticks(xTickLocs);
      xticklabels(xTickLabels);
      xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
      ax1 = gca;
      ax1.YAxis.Visible = 'off';
      ylim([-iMean-0.5 0])

      % Add the figure title
      tH = title('Sorted Unit Phases', 'Interpreter','none');

      % Save the figure
      if ~exist(figPath,'dir')
        mkdir(figPath);
      end
      filename = strrep(tH.String,' ','_');
      filename = strrep(filename, ',','_');
      filename = strrep(filename, '-','_');
      filename = strrep(filename, ':','' );
      filename = strrep(filename, '^','' );
      filename = strrep(filename, '{','' );
      filename = strrep(filename, '}','' );
      filename = strrep(filename, '.','p');
      filename = [figPath filesep filename '.fig']; %#ok<*AGROW>
      savefig(fH3,filename,'compact');
      title('');
      saveas(fH3,filename(1:end-4),'png');
      saveas(fH3,filename(1:end-4),'pdf');
      close(fH3);

      % Phase distribution SDs resorted based on behaviour
      if ~isempty(options.secondarySorting) && ...
          ~isempty(options.secondarySorting{animal}{session})
        sortedMeans = phaseMeans(options.secondarySorting{animal}{session});
        sortedSDs = phaseSDs(options.secondarySorting{animal}{session});

        fH4 = figure; hold on
        for iMean = 1:numel(sortedMeans)
          if ~isnan(sortedMeans(iMean))
            plot(sortedMeans(iMean),-iMean, 'b.', 'MarkerSize',10);
            plot([sortedMeans(iMean)-sortedSDs(iMean)/2 ...
              sortedMeans(iMean)+sortedSDs(iMean)/2],[-iMean -iMean], 'b', ...
              'LineWidth',1);
          else
            break
          end
        end

        % Label axes
        xTickLocs = [-pi, -pi/2, 0, pi/2, pi];
        xTickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
        xticks(xTickLocs);
        xticklabels(xTickLabels);
        xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
        ax1 = gca;
        ax1.YAxis.Visible = 'off';
        ylim([-iMean-0.5 0])

        % Add the figure title
        tH = title('Sorted Unit Phases based on behaviour', 'Interpreter','none');

        % Save the figure
        if ~exist(figPath,'dir')
          mkdir(figPath);
        end
        filename = strrep(tH.String,' ','_');
        filename = strrep(filename, ',','_');
        filename = strrep(filename, '-','_');
        filename = strrep(filename, ':','' );
        filename = strrep(filename, '^','' );
        filename = strrep(filename, '{','' );
        filename = strrep(filename, '}','' );
        filename = strrep(filename, '.','p');
        filename = [figPath filesep filename '.fig']; %#ok<*AGROW>
        savefig(fH4,filename,'compact');
        title('');
        saveas(fH4,filename(1:end-4),'png');
        saveas(fH4,filename(1:end-4),'pdf');
        close(fH4);
      end
    end
  end
end
end


function plotUnitPhaseAcrossSessions(dataFiles, phaseHistosSpeed, ...
  phaseHistosBehav, phaseHistosREM, phaseHistosNREM, figFolder)
% plotUnitPhaseAcrossSessions(dataFiles, phaseHistosSpeed, ...
%   phaseHistosBehav, phaseHistosREM, phaseHistosNREM, figFolder);
%
% Function displays unit spiking phase histograms across same day recording
% sessions.
%
% Args:
%   dataFiles
%   phaseHistosSpeed
%   phaseHistosBehav
%   phaseHistosREM
%   phaseHistosNREM
%   figFolder
%
% Returns:
%   None.
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  phaseHistosSpeed
  phaseHistosBehav
  phaseHistosREM
  phaseHistosNREM
  figFolder
end

fontSize = 14;

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    nUnits = max([numel(phaseHistosSpeed.unitIDs{animal}{session}) ...
      numel(phaseHistosBehav.unitIDs{animal}{session}) ...
      numel(phaseHistosREM.unitIDs{animal}{session}) ...
      numel(phaseHistosNREM.unitIDs{animal}{session})]);
    for unit = 1:nUnits

      % Loop through units
      for condition = 1:4
        if condition == 1
          phaseHistos = phaseHistosSpeed;
        elseif condition == 2
          phaseHistos = phaseHistosBehav;
        elseif condition == 3
          phaseHistos = phaseHistosREM;
        elseif condition == 4
          phaseHistos = phaseHistosNREM;
        end
        if ~isempty(phaseHistos.histograms{animal}{session}) && ...
            phaseHistos.validUnits{animal}{session}(unit) && ...
            sum(phaseHistos.histograms{animal}{session}(unit,:))

          if ~exist('fH1', 'var')
            fH1 = figure; hold on
            legendStr = {};
          end
          plot(phaseHistos.bins{animal}{session}(unit,:), ...
            phaseHistos.histograms{animal}{session}(unit,:)./ ...
            sum(phaseHistos.histograms{animal}{session}(unit,:)));

          if condition == 1
            legendStr{numel(legendStr)+1} = 'Run';
          elseif condition == 2
            legendStr{numel(legendStr)+1} = 'Behave';
          elseif condition == 3
            legendStr{numel(legendStr)+1} = 'REM';
          elseif condition == 4
            legendStr{numel(legendStr)+1} = 'NREM';
          end
        end

        if condition == 4 && exist('fH1','var')
          % Label axes
          xTickLocs = [-pi, -pi/2, 0, pi/2, pi];
          xTickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
          xticks(xTickLocs);
          xticklabels(xTickLabels);
          xlabel('Phase (rad)', 'FontSize',fontSize, 'FontWeight','bold');
          ylabel('Normalised spike count', 'FontSize',fontSize, 'FontWeight','bold');

          legend(legendStr);
          legend('boxoff');

          % Add the figure title
          histTitle = ['uid ' num2str(phaseHistos.unitIDs{animal}{session}(unit)) ...
            ' ch ' num2str(phaseHistos.channelIDs{animal}{session}(unit)) ...
            ' spiking phase histograms'];
          tH = title(histTitle, 'Interpreter','none');

          % Save the figure
          [~, basename] = fileparts(dataFiles{animal}{session});
          basename = strrep(basename,'.*','');
          figPath = fullfile(figFolder, basename);
          if ~exist(figPath,'dir')
            mkdir(figPath);
          end
          filename = strrep(tH.String,' ','_');
          filename = strrep(filename, ',','_');
          filename = strrep(filename, '-','_');
          filename = strrep(filename, ':','' );
          filename = strrep(filename, '^','' );
          filename = strrep(filename, '{','' );
          filename = strrep(filename, '}','' );
          filename = strrep(filename, '.','p');
          filename = [figPath filesep filename '.fig']; %#ok<*AGROW>
          savefig(fH1,filename,'compact');
          title('');
          saveas(fH1,filename(1:end-4),'png');
          saveas(fH1,filename(1:end-4),'pdf');
          close(fH1); clear fH1
        end
      end
    end
  end
end
end

% phaseHistos.histograms = cell(numel(dataFiles),1);
% phaseHistos.bins = cell(numel(dataFiles),1);
% phaseHistos.unitIDs = cell(numel(dataFiles),1);
% phaseHistos.channelIDs = cell(numel(dataFiles),1);
% phaseHistos.validUnits = cell(numel(dataFiles),1);