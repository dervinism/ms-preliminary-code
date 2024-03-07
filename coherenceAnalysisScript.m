% This script performs coherence analysis for ms-preliminary dataset


%% Initialise parameters
% Data setup parameters
params

% Analysis parameters
resamplingInterval = 0.02;
frequencyRange = [4 12];
minIntervalLength = 1;
parallelise = true;
nPhaseHistogramBins = 9;

% Unit quality criteria
firingRateCutoff = 200; % At least that many spikes per hour or at least a single 20-minute window (i.e., firingRateCutoff/3)
refractoryContaminationTh = 5; % Not more than this percentage of spiking ACG refractory period contamination.
                               % Will take the ratio mean(+/- 1 ms relative to 0 lag of ACG) /
                               % mean(1000-1500ms ACG away from zero lag (ACG shoulder)).
coherenceTh = 0.0; % Exclude units with coherence in the theta frequency range below this value

% Other unit selection criteria
oscTh = 0; %oscSignificanceCutoff;

% Script run parameters
loadData = false;
runThetaMazeAnalyses = true;
runHighThetaAnalyses = false;
runSleepAnalyses = true;


%% Load all preprocessed data
if loadData || (~exist('data', 'var') || isempty(data))
  data = getData(animals, dataFiles);
end


%% Measure theta phase topography in the medial septum
if runThetaMazeAnalyses

  % Define periods of interest: ThetaMaze_AlternativeRunning
  epochsOfInterest = 'ThetaMaze_AlternativeRunning';
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  figTitleInit = epochsOfInterest;

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed, excludeNoise=true);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitleInit, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins, saveResults=false);


  % Define periods of interest: ThetaMaze_AlternativeRunning onlyTrials
  onlyTrials = true;
  onlyHighSpeed = false;
  figTitle = [figTitleInit '_trials'];

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);


  % Define periods of interest: ThetaMaze_AlternativeRunning onlyHighSpeed
  onlyTrials = false;
  onlyHighSpeed = true;
  figTitle = [figTitleInit '_highSpeed'];

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);


  % Define periods of interest: ThetaMaze_AlternativeRunning onlyTrials onlyHighSpeed
  onlyTrials = true;
  onlyHighSpeed = true;
  figTitle = [figTitleInit '_trials_highSpeed'];

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);
end


if runHighThetaAnalyses

  % Define periods of interest: HomeCage_Sleep
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake'}; %#ok<*UNRCH> 
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  figTitleInit = 'Homecage';

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitleInit, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);


  % Define periods of interest: HomeCage_Sleep moderateTheta
  thetaPower = 'moderate';
  figTitle = [figTitleInit '_moderateTheta'];

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);


  % Define periods of interest: HomeCage_Sleep highTheta
  thetaPower = 'high';
  figTitle = [figTitleInit '_strongTheta'];

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);
end


if runSleepAnalyses
  % Define periods of interest: All sessions with REM periods
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunnig','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem';
  figTitle = 'REM';

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed, sleepState=sleepState);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);


  % Define periods of interest: All sessions with NREM periods
  sleepState = 'nrem';
  figTitle = 'NREM';

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed, sleepState=sleepState);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);


  % Define periods of interest: All sessions when animal is awake in its homecage
  epochsOfInterest = {'HomeCage_Sleep','HomeCage_Wake'};
  sleepState = 'wake';
  figTitle = 'wake';

  % Find the time intervals of interest
  intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, ...
    thetaPower=thetaPower, minIntervalLength=minIntervalLength, ...
    onlyTrials=onlyTrials, onlyHighSpeed=onlyHighSpeed, sleepState=sleepState);

  % Probe the theta phase topography
  thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
    cohFigFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
    firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
    nPhaseHistogramBins);
end



%% Local functions
function data = getData(animals, dataFiles)
% data = loadData(animals, dataFiles)
%
% Function loads preprocessed data for multiple recording sessions located
% in the ms-preliminary repository and returns it stored in a single data
% structure variable. It's a helper function of coherenceAnalysisScript.
%
% Args:
%   animals
%   dataFiles
%
% Returns:
%   data
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

data = cell(numel(dataFiles),1);
for animal = 1:numel(animals)
  for session = 1:numel(dataFiles{animal})
    dataFilename = dataFiles{animal}{session};
    if exist(dataFilename, 'file') % If there is an exact file match, load it
      data{animal}{session} = load(dataFilename);
    else
      folderContents = dir(fileparts(dataFilename));
      data{animal}{session} = [];
      if ~isempty(folderContents)
        for file = 1:numel(folderContents) % If there are multiple matching files, load them all
          checkedFolder = fileparts(dataFilename);
          checkedFilename = fullfile(checkedFolder, folderContents(file).name);
          if regexp(checkedFilename, regexptranslate('wildcard', dataFilename))
            fileContents = load(checkedFilename);
            loadedFieldNames = fieldnames(fileContents);
            for field = 1:numel(loadedFieldNames)
              data{animal}{session}.(loadedFieldNames{field}) = fileContents.(loadedFieldNames{field});
            end
          end
        end
        if ~isempty(data{animal}{session}) % Check if the loaded content is meaningful
          loadedFieldNames = fieldnames(data{animal}{session});
          if numel(loadedFieldNames) == 1
            data{animal}{session} = [];
          end
        end
      end
    end
  end
end
end


function data = getRefractoryContamination(data, dataFiles) %#ok<*DEFNU> 
% data = getRefractoryContamination(data, dataFiles)
%
% Function estimates refractory period contamination of individual unit
% autocorellograms for multiple recording sessions located in the
% ms-preliminary repository and returns it stored in a single data
% structure variable. It's a helper function of coherenceAnalysisScript.
%
% Args:
%   data
%   dataFiles
%
% Returns:
%   data
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(data{animal}{session})
      sessionData = data{animal}{session};
      nUnits = numel(sessionData.spikes.ids);
      sessionData.contaminationPercentage = zeros(nUnits,1);
      for unit = 1:nUnits
        sessionData.contaminationPercentage(unit) = ...
          refractoryContamination(sessionData.spikes.times{unit}');
      end
      data{animal}{session} = sessionData;
    end
  end
end
end


function intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, options)
% intervals = getTimeIntervals(data, dataFiles, epochsOfInterest, <options>)
%
% Function extracts time intervals of interest for multiple recording
% sessions located in the ms-preliminary repository.
%
% Args:
%   data
%   dataFiles
%   epochsOfInterest
%   <thetaPower>
%   <minIntervalLength>
%   <onlyTrials>
%   <onlyHighSpeed>
%   <sleepState>
%   <excludeNoise>
%
% Returns:
%   intervals
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  data
  dataFiles
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

% Proceed with interval selection
intervals = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(data{animal}{session})
      sessionData = data{animal}{session};
      sessionInfo = sessionData.session;
      sessionIntervals = [];
      for epoch = 1:numel(sessionInfo.epochs)
        if isfield(sessionInfo.epochs{epoch},'behavioralParadigm') && ...
            ismember(sessionInfo.epochs{epoch}.behavioralParadigm, epochsOfInterest)

          % Select time intervals corresponding to the sessions of interest
          interval = [sessionInfo.epochs{epoch}.startTime sessionInfo.epochs{epoch}.stopTime];

          % Exclude intervals with noisy (saturated) LFPs
          if options.excludeNoise && isfield(sessionData,'lfpNoise')
            cleanIntervals = invertIntervals(sessionData.lfpNoise.timestamps, ...
              sessionInfo.epochs{epoch}.stopTime, startTime=sessionInfo.epochs{epoch}.startTime);
            if ~isempty(cleanIntervals)
              interval = intervalOverlap(interval, cleanIntervals);
            else
              pause(1)
            end
          end

          % Select time intervals corresponding to increased theta power
          if (strcmpi(options.thetaPower,'moderate') || strcmpi(options.thetaPower,'high')) && ~isempty(interval)
            if strcmpi(options.thetaPower,'moderate')
              increasedThetaPowerIntervals = sessionData.thetaPower.data(:,3)' | sessionData.thetaPower.data(:,4)';
            elseif strcmpi(options.thetaPower,'high')
              increasedThetaPowerIntervals = logical(sessionData.thetaPower.data(:,4)');
            end
            increasedThetaPowerIntervals = logical2intervals(increasedThetaPowerIntervals);
            increasedThetaPowerIntervals(:,2) = increasedThetaPowerIntervals(:,2) + 1;
            increasedThetaPowerIntervals(end,2) = min([increasedThetaPowerIntervals(end,2) numel(sessionData.thetaPower.timestamps)]);
            increasedThetaPowerIntervals = sessionData.thetaPower.timestamps(increasedThetaPowerIntervals);
            interval = intervalOverlap(interval, increasedThetaPowerIntervals);
          end

          % Select time intervals corresponding to trials only
          if options.onlyTrials && ~isempty(interval)
            if isfield(sessionData,'circular_track') && isfield(sessionData.circular_track,'trials')
              trialIntervals = [sessionData.circular_track.trials.alternation.start ...
                sessionData.circular_track.trials.alternation.end];
            else
              trialIntervals = [];
            end
            interval = intervalOverlap(interval, trialIntervals);
          end

          % Select time intervals corresponding to high speeds
          if options.onlyHighSpeed && ~isempty(interval)
            if isfield(sessionData,'circular_track') && isfield(sessionData.circular_track,'speed')
              highSpeedIntervals = sessionData.circular_track.speed >= sessionData.circular_track.speed_th;
              highSpeedIntervals = logical2intervals(highSpeedIntervals);
              highSpeedIntervals(:,2) = highSpeedIntervals(:,2) + 1;
              highSpeedIntervals(end,2) = min([highSpeedIntervals(end,2) numel(sessionData.circular_track.timestamps)]);
              highSpeedIntervals = sessionData.circular_track.timestamps(highSpeedIntervals);
            else
              highSpeedIntervals = [];
            end
            interval = intervalOverlap(interval, highSpeedIntervals);
          end

          % Select intervals corresponding to a brain state of interest
          if ~isempty(options.sleepState)
            if isfield(sessionData,'SleepState') && isfield(sessionData.SleepState,'ints')
              switch options.sleepState
                case 'wake'
                  interval = intervalOverlap(interval, sessionData.SleepState.ints.WAKEstate);
                case 'ma'
                  interval = intervalOverlap(interval, sessionData.SleepState.ints.MAstate);
                case 'nrem'
                  interval = intervalOverlap(interval, sessionData.SleepState.ints.NREMstate);
                case 'rem'
                  interval = intervalOverlap(interval, sessionData.SleepState.ints.REMstate);
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
      intervals{animal}{session} = sessionIntervals;
    end
  end
end
end