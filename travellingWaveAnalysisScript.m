% This script detects travelling theta waves and their direction for ms-preliminary dataset


%% Initialise parameters
% Data setup parameters
params
figureFolder = travellingWavesFigFolder;

% Analysis parameters
samplingInterval = 0.002;
frequencyRange = [4 12];
minIntervalLength = 1;
electrodeAxis = 'all';
envTh = [];
minTravelDuration = 0.05;

% Spiking only
detectSpikingWaves = false;
firingRateTh = 0.1;
oscillationThSpikes = 0;
pgdThSpikes = [];
aggregateSpikingWaves = true;
runThetaMazeAnalysesSpiking = true;
runSleepAnalysesSpiking = true;

% LFP only
detectLFPWaves = false;
lfpSamplingInterval = 0.0004;
oscillationThLFP = false;
pgdThLFP = [];
aggregateLFPWaves = true;
runThetaMazeAnalysesLFP = true;
runSleepAnalysesLFP = true;


%% Preprocessing steps for spiking data
% Detect and save spiking travelling waves
if detectSpikingWaves
  analyseTravellingChannelWaves(dataFiles, samplingInterval, probeLayouts, ...
    channelsOI, freqRange=frequencyRange, axis=electrodeAxis, ...
    omitnans=true, firingRateTh=firingRateTh, ...
    oscillationTh=oscillationThSpikes, pgdTh=pgdThSpikes, envTh=envTh); %#ok<*UNRCH>
end

% Aggregate travelling theta spiking wave direction data
if aggregateSpikingWaves
  [travellingWaveDirSpk, travellingWaveDirSpk2, travellingWaveDirSpk3, ...
    travellingWaveTimesSpk, travellingWaveDirSpk4, travellingWaveTimesSpk4, ...
    timestampsSpk, pgdSignificanceCutoffSpk] = aggregateWaveData( ...
    dataFiles, minTravelDuration=minTravelDuration, dataType='spikes', ...
    saveData=true);
end


%% Summary analysis steps for spiking data
% Calculate and plot travelling theta spiking wave direction incidence rates
if runThetaMazeAnalysesSpiking

  % Define periods of interest: ThetaMaze_AlternativeRunning
  epochsOfInterest = 'ThetaMaze_AlternativeRunning';
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' epochsOfInterest];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(epochsOfInterest);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: ThetaMaze_AlternativeRunning_trials
  onlyTrials = true;
  onlyHighSpeed = false;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' epochsOfInterest '_trials'];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp([epochsOfInterest '_trials']);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: ThetaMaze_AlternativeRunning_highSpeed
  onlyTrials = false;
  onlyHighSpeed = true;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' epochsOfInterest '_highSpeed'];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp([epochsOfInterest '_highSpeed']);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: ThetaMaze_AlternativeRunning_trials_highSpeed
  onlyTrials = true;
  onlyHighSpeed = true;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' epochsOfInterest '_trials_highSpeed'];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp([epochsOfInterest '_trials_highSpeed']);
  disp(recSessions{1});
  disp(recSessions{2});
end


if runSleepAnalysesSpiking
  % Define periods of interest: REM
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunnig','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem';

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' sleepState];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(sleepState);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: NREM
  sleepState = 'nrem';

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' sleepState];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(sleepState);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: All sessions when animal is awake in its homecage
  epochsOfInterest = {'HomeCage_Sleep','HomeCage_Wake'};
  sleepState = 'wake';

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate travelling theta spiking wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirSpk, ...
    travellingWaveDirSpk2, travellingWaveDirSpk3, travellingWaveTimesSpk, ...
    travellingWaveDirSpk4, travellingWaveTimesSpk4, timestampsSpk, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_spiking_' sleepState];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(sleepState);
  disp(recSessions{1});
  disp(recSessions{2});
end


%% Preprocessing steps for LFP data
% Detect and save LFP travelling waves
pgdThLFP = 0.379773934741062; %pgdSignificanceCutoffSpk(1);
if detectLFPWaves
  analyseTravellingLFPWaves(lfpFiles, lfpSamplingInterval, samplingInterval, ...
    probeLayouts, channelsOI, freqRange=frequencyRange, axis=electrodeAxis, ...
    omitnans=true, oscillationTh=oscillationThLFP, pgdTh=pgdThLFP, envTh=envTh);
end

% Aggregate travelling theta LFP wave direction data
if aggregateLFPWaves
  [travellingWaveDirLFP, travellingWaveDirLFP2, travellingWaveDirLFP3, ...
    travellingWaveTimesLFP, travellingWaveDirLFP4, travellingWaveTimesLFP4, ...
    timestampsLFP, pgdSignificanceCutoffLFP] = aggregateWaveData( ...
    dataFiles, minTravelDuration=minTravelDuration, dataType='lfp', ...
    saveData=true);
end


%% Summary analysis steps for spiking data
% Calculate and plot travelling theta LFP wave direction incidence rates
if runThetaMazeAnalysesLFP

  % Define periods of interest: ThetaMaze_AlternativeRunning
  epochsOfInterest = 'ThetaMaze_AlternativeRunning';
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' epochsOfInterest];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(epochsOfInterest);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: ThetaMaze_AlternativeRunning_trials
  onlyTrials = true;
  onlyHighSpeed = false;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' epochsOfInterest '_trials'];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp([epochsOfInterest '_trials']);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: ThetaMaze_AlternativeRunning_highSpeed
  onlyTrials = false;
  onlyHighSpeed = true;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' epochsOfInterest '_highSpeed'];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp([epochsOfInterest '_highSpeed']);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: ThetaMaze_AlternativeRunning_trials_highSpeed
  onlyTrials = true;
  onlyHighSpeed = true;

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' epochsOfInterest '_trials_highSpeed'];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp([epochsOfInterest '_trials_highSpeed']);
  disp(recSessions{1});
  disp(recSessions{2});
end


if runSleepAnalysesLFP
  % Define periods of interest: REM
  epochsOfInterest = {'Homecage_Sleep','Homecage_Wake', ...
    'ThetaMaze_AlternativeRunnig','ThetaMaze_FreeRunning','LinearTrack_EndToEnd'};
  thetaPower = '';
  onlyTrials = false;
  onlyHighSpeed = false;
  sleepState = 'rem';

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' sleepState];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(sleepState);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: NREM
  sleepState = 'nrem';

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' sleepState];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(sleepState);
  disp(recSessions{1});
  disp(recSessions{2});


  % Define periods of interest: All sessions when animal is awake in its homecage
  epochsOfInterest = {'HomeCage_Sleep','HomeCage_Wake'};
  sleepState = 'wake';

  % Find the time intervals of interest
  nAnimals = numel(dataFiles);
  intervals = cell(nAnimals,1);
  for animal = 1:nAnimals
    for session = 1:numel(dataFiles{animal})
      intervals{animal}{session} = getTimeIntervals( ...
        dataFiles{animal}{session}, epochsOfInterest, thetaPower=thetaPower, ...
        minIntervalLength=minIntervalLength, onlyTrials=onlyTrials, ...
        onlyHighSpeed=onlyHighSpeed, sleepState=sleepState, excludeNoise=true);
    end
  end

  % Calculate travelling theta LFP wave incidence rates
  [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
    recSessions] = waveIncidenceRate(dataFiles, intervals, travellingWaveDirLFP, ...
    travellingWaveDirLFP2, travellingWaveDirLFP3, travellingWaveTimesLFP, ...
    travellingWaveDirLFP4, travellingWaveTimesLFP4, timestampsLFP, ...
    minIntervalLength=minIntervalLength);

  % Plot and save incidence rates
  figFilenameSuffix = ['_lfp_' sleepState];
  plotIncidenceRates(figureFolder, figFilenameSuffix, ...
    travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
    dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
    dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
    travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4);
  disp(sleepState);
  disp(recSessions{1});
  disp(recSessions{2});
end



%% Local functions
function analyseTravellingChannelWaves(dataFiles, samplingInterval, ...
  channelOrder, channelsOI, options) %#ok<*DEFNU>
% analyseTravellingChannelWaves(dataFiles, samplingInterval, ...
%   channelOrder, channelsOI, <options>)
%
% Function detects travelling recording channel spiking waves, their
% direction, saves this data in accordance with the CellExplorer format,
% and performs further analyses of travelling waves.
%
% Args:
%   dataFiles
%   samplingInterval
%   channelOrder
%   channelsOI
%   options.freqRange
%   options.axis ('vertical','horizontal','all')
%   options.omitnans (true,false)
%   options.firingRateTh (default=0)
%   options.oscillationTh (default=0)
%   options.pgdTh (default=0)
%   options.envTh (default=0)
%
% Returns:
%   None.
%
% Dependencies
%   CellExplorer (https://cellexplorer.org/)
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   The Oscillation Score (https://www.raulmuresan.ro/sources/oscore/).
%   erfanzabeh/WaveMonk (https://github.com/erfanzabeh/WaveMonk).
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  samplingInterval
  channelOrder
  channelsOI
  options.freqRange
  options.axis
  options.omitnans
  options.firingRateTh = 0
  options.oscillationTh = 0
  options.pgdTh = 0
  options.envTh = 0
end

% Detect and save theta travelling waves
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session})
      disp(dataFiles{animal}{session})
      saveTravellingThetaChannelWaves(dataFiles{animal}{session}, ...
        samplingInterval=samplingInterval, ...
        channelOrder=channelOrder{animal}(session,:), ...
        channelsOI=channelsOI{animal}(session,:), ...
        freqRange=options.freqRange, axis=options.axis, omitnans=true, ...
        firingRateTh=options.firingRateTh, ...
        oscillationTh=options.oscillationTh, pgdTh=options.pgdTh, ...
        envTh=options.envTh, verbose=true);
    end
  end
end
end


function analyseTravellingLFPWaves(lfpFiles, lfpSamplingInterval, ...
  samplingInterval, channelOrder, channelsOI, options)
% analyseTravellingLFPWaves(lfpFiles, lfpSamplingInterval, ...
%   samplingInterval, channelOrder, channelsOI, <options>)
%
% Function detects travelling recording channel LFP waves, their direction,
% saves this data in accordance with the CellExplorer format, and performs
% further analyses of travelling waves.
%
% Args:
%   lfpFiles
%   lfpSamplingInterval
%   samplingInterval
%   channelOrder
%   channelsOI
%   options.freqRange
%   options.axis ('vertical','horizontal','all')
%   options.omitnans (true,false)
%   options.oscillationTh (default=0)
%   options.pgdTh (default=0)
%   options.envTh (default=0)
%
% Returns:
%   None.
%
% Dependencies
%   CellExplorer (https://cellexplorer.org/)
%   petersen-lab/petersen-lab-matlab
%     (https://github.com/petersen-lab/petersen-lab-matlab/).
%   erfanzabeh/WaveMonk (https://github.com/erfanzabeh/WaveMonk).
%   Circular Statistics Toolbox
%     (https://uk.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics).
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  lfpFiles
  lfpSamplingInterval
  samplingInterval
  channelOrder
  channelsOI
  options.freqRange
  options.axis
  options.omitnans
  options.oscillationTh = 0
  options.pgdTh = 0
  options.envTh = 0
end

% Detect and save theta travelling waves
for animal = 1:numel(lfpFiles)
  for session = 1:numel(lfpFiles{animal})
    if ~isempty(lfpFiles{animal}{session})
      disp(lfpFiles{animal}{session})
      saveTravellingThetaLFPWaves(lfpFiles{animal}{session}, ...
        lfpSamplingInterval=lfpSamplingInterval, ...
        samplingInterval=samplingInterval, ...
        channelOrder=channelOrder{animal}(session,:), ...
        channelsOI=channelsOI{animal}(session,:), ...
        freqRange=options.freqRange, axis=options.axis, omitnans=true, ...
        oscillationTh=options.oscillationTh, pgdTh=options.pgdTh, ...
        envTh=options.envTh, verbose=true);
    end
  end
end
end


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


function [travellingWaveDir, travellingWaveDir2, travellingWaveDir3, ...
  travellingWaveTimes, travellingWaveDir4, travellingWaveTimes4, ...
  timestamps, pgdSignificanceCutoff] = aggregateWaveData(dataFiles, options)
% [travellingWaveDir, travellingWaveDir2, travellingWaveDir3, ...
%   travellingWaveTimes, travellingWaveDir4, travellingWaveTimes4, ...
%   timestamps, pgdSignificanceCutoff] = aggregateWaveData(dataFiles, <options>)
%
% Function aggregates travelling theta wave detection data for multiple
% recording sessions located in the ms-preliminary repository.
%
% Args:
%   dataFile
%   <smoothingWindowSize>
%   <minTravelDuration>
%   <dataType>
%   <saveData>
%
% Returns:
%   travellingWaveDir
%   travellingWaveDir2
%   travellingWaveDir3
%   travellingWaveTimes
%   travellingWaveDir4
%   travellingWaveTimes4
%   timestamps
%   pgdSignificanceCutoff
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  dataFiles
  options.smoothingWindowSize = 0.1
  options.minTravelDuration = 0.05
  options.dataType (1,:) {mustBeMember(options.dataType,{'spikes','lfp'})} = 'spikes'
  options.saveData (1,1) {mustBeA(options.saveData,'logical')} = false
end

nAnimals = numel(dataFiles);
travellingWaveDir = cell(nAnimals,1);
travellingWaveDir2 = cell(nAnimals,1);
travellingWaveDir3 = cell(nAnimals,1);
travellingWaveTimes = cell(nAnimals,1);
travellingWaveDir4 = cell(nAnimals,1);
travellingWaveTimes4 = cell(nAnimals,1);
travellingWaveLocs4_postprocess = cell(nAnimals,1);
timestamps = cell(nAnimals,1);
pgdIndexShuffled = cell(nAnimals,1);
for animal = 1:nAnimals
  for session = 1:numel(dataFiles{animal})
    disp(dataFiles{animal}{session});

    % Load data
    if strcmpi(options.dataType,'spikes')
      sessionDataFile = strrep(dataFiles{animal}{session}, '*', ...
        'travellingSpikingThetaWave.timeseriesCollection');
    elseif strcmpi(options.dataType,'lfp')
      sessionDataFile = strrep(dataFiles{animal}{session}, '*', ...
        'travellingLFPThetaWave.timeseriesCollection');
    end
    if exist(sessionDataFile,'file')
      load(sessionDataFile); %#ok<*LOAD>
    else
      travellingWaveDir{animal}{numel(travellingWaveDir{animal})+1} = {};
      travellingWaveDir2{animal}{numel(travellingWaveDir2{animal})+1} = {};
      travellingWaveDir3{animal}{numel(travellingWaveDir3{animal})+1} = {};
      travellingWaveTimes{animal}{numel(travellingWaveTimes{animal})+1} = {};
      travellingWaveDir4{animal}{numel(travellingWaveDir4{animal})+1} = {};
      travellingWaveTimes4{animal}{numel(travellingWaveTimes4{animal})+1} = {};
      travellingWaveLocs4_postprocess{animal}{numel(travellingWaveLocs4_postprocess{animal})+1} = {};
      timestamps{animal}{numel(timestamps{animal})+1} = {};
      continue
    end
    if strcmpi(options.dataType,'spikes')
      travellingWave = travellingSpikingThetaWave;
    elseif strcmpi(options.dataType,'lfp')
      travellingWave = travellingLFPThetaWave;
    end
    ind = startsWith(travellingWave.channelNames,'netWaveDir') ...
      & endsWith(travellingWave.channelNames,'netWaveDir');
    netWaveDir = travellingWave.data(:,ind);
    ind = startsWith(travellingWave.channelNames,'travellingWaveLocs') ...
      & endsWith(travellingWave.channelNames,'travellingWaveLocs');
    travellingWaveLocs = logical(travellingWave.data(:,ind));
    ind = startsWith(travellingWave.channelNames,'pgdIndex') ...
      & endsWith(travellingWave.channelNames,'pgdIndex');
    pgdIndex = travellingWave.data(:,ind);
    if strcmpi(options.dataType,'spikes')
      ind = startsWith(travellingWave.channelNames,'pgdIndexShuffled') ...
        & endsWith(travellingWave.channelNames,'pgdIndexShuffled');
      pgdIndexShuffled_session = travellingWave.data(:,ind)';
    elseif strcmpi(options.dataType,'lfp')
      pgdIndexShuffled_session = travellingWave.data(:,ind)';
    end
    pgdIndexShuffled{animal} = [pgdIndexShuffled{animal} pgdIndexShuffled_session];
    ind = startsWith(travellingWave.channelNames,'travellingWaveCycleLocs') ...
      & endsWith(travellingWave.channelNames,'travellingWaveCycleLocs');
    travellingWaveCycleLocs = logical(travellingWave.data(:,ind))';
    ind = startsWith(travellingWave.channelNames,'cycleNumbers') ...
      & endsWith(travellingWave.channelNames,'cycleNumbers');
    cycleNumbers = travellingWave.data(:,ind);

    % Smooth data and estimate smoothed PGDI significance cutoff
    samplingInterval = travellingWave.timestamps(2) - travellingWave.timestamps(1);
    pgdIndexSmooth = movmean(pgdIndex, round(options.smoothingWindowSize/samplingInterval));
    pgdIndexShuffledSmooth = movmean(pgdIndexShuffled_session, ...
      round(options.smoothingWindowSize/samplingInterval));
    pgdIndexCutoff = prctile(pgdIndexShuffledSmooth,95);
    travellingWaveLocs4 = pgdIndexSmooth' > pgdIndexCutoff;

    % Eliminate short travel periods
    waveStartInds = find([0 diff(travellingWaveLocs4)] > 0);
    waveEndInds = find([0 diff(travellingWaveLocs4)] < 0);
    if waveStartInds(1) > waveEndInds(1)
      waveEndInds(1) = [];
      waveStartInds = waveStartInds(1:numel(waveEndInds));
    elseif waveStartInds(end) > waveEndInds(end)
      waveStartInds(end) = [];
      waveEndInds = waveEndInds(1:numel(waveStartInds));
    end
    assert(numel(waveStartInds) == numel(waveEndInds));
    acceptableWaves = (waveEndInds - waveStartInds)*samplingInterval >= options.minTravelDuration;
    waveStartInds = waveStartInds(acceptableWaves);
    waveEndInds = waveEndInds(acceptableWaves);
    travellingWaveLocs4 = false(size(travellingWaveLocs4));
    for iWave = 1:sum(acceptableWaves)
      travellingWaveLocs4(waveStartInds(iWave):waveEndInds(iWave)) = true;
    end

    %figure; plot(travellingWave.timestamps,pgdIndexSmooth); hold on
    %plot(travellingWave.timestamps,pgdIndexCutoff.*ones(size(pgdIndexSmooth)));
    %plot(travellingWave.timestamps,travellingWaveLocs4); hold off
    %ylim([-0.1 1.1]);

    % Aggregate data: Legacy method
    p = gcp('nocreate');
    if isempty(p)
      parpool('local', feature('numcores'));
    elseif p.NumWorkers < feature('numcores')
      delete(gcp('nocreate'));
      parpool('local', feature('numcores'));
    end
    nTravellingWaves = numel(unique(cycleNumbers(travellingWaveCycleLocs)));
    if nTravellingWaves
      travellingWaveCycles = sort(unique(cycleNumbers(travellingWaveCycleLocs)));
      travellingWaveDir_session = zeros(1,nTravellingWaves);
      travellingWaveDir2_session = zeros(1,nTravellingWaves);
      travellingWaveDir3_session = zeros(1,nTravellingWaves);
      travellingWaveTimes_session = zeros(1,nTravellingWaves);
      parfor iWave = 1:nTravellingWaves
        %disp(iCycle/nTravellingWaves);
        cycleLocs = cycleNumbers == travellingWaveCycles(iWave);
        selectPGDI = zeros(size(pgdIndex));
        selectPGDI(cycleLocs) = pgdIndex(cycleLocs);
        [~, maxCyclePGDI] = max(selectPGDI);
        travellingWaveDir_session(iWave) = netWaveDir(maxCyclePGDI); %#ok<*PFBNS>
        travellingWaveDir2_session(iWave) = ...
          (pi/2)*(sum(netWaveDir(cycleLocs))/abs(sum(netWaveDir(cycleLocs))));
        travellingWaveDir3_session(iWave) = ...
          (pi/2)*(sum(netWaveDir(cycleLocs&travellingWaveLocs)) / ...
          abs(sum(netWaveDir(cycleLocs&travellingWaveLocs))));
        travellingWaveTimes_session(iWave) = ...
          travellingWave.timestamps(maxCyclePGDI);
      end
      iCell = numel(travellingWaveDir{animal})+1;
      travellingWaveDir{animal}{iCell} = travellingWaveDir_session;
      travellingWaveDir2{animal}{iCell} = travellingWaveDir2_session;
      travellingWaveDir3{animal}{iCell} = travellingWaveDir3_session;
      travellingWaveTimes{animal}{iCell} = travellingWaveTimes_session;
    else
      iCell = numel(travellingWaveDir{animal})+1;
      travellingWaveDir{animal}{iCell} = {};
      travellingWaveDir2{animal}{iCell} = {};
      travellingWaveDir3{animal}{iCell} = {};
      travellingWaveTimes{animal}{iCell} = {};
    end

    % Aggregate data: New method
    travellingWaveInds = logical2intervals(travellingWaveLocs4);
    travellingWaveInds = travellingWaveInds( ...
      travellingWaveInds(:,2) - travellingWaveInds(:,1) >= options.minTravelDuration, :);
    nTravellingWaves4 = size(travellingWaveInds, 1);
    if nTravellingWaves4
      travellingWaveDir4_session = nan(1,nTravellingWaves4);
      travellingWaveTimes4_session = nan(1,nTravellingWaves4);
      parfor iWave = 1:nTravellingWaves4
        indsOI = travellingWaveInds(iWave,1):travellingWaveInds(iWave,2);
        if sum(netWaveDir(indsOI))
          travellingWaveDir4_session(iWave) = ...
            (pi/2)*(sum(netWaveDir(indsOI)) / ...
            abs(sum(netWaveDir(indsOI))));
          travellingWaveTimes4_session(iWave) = ...
            mean(travellingWave.timestamps(indsOI));
        end
      end
      travellingWaveDir4{animal}{numel(travellingWaveDir4{animal})+1} = travellingWaveDir4_session;
      travellingWaveTimes4{animal}{numel(travellingWaveTimes4{animal})+1} = travellingWaveTimes4_session;
      travellingWaveDir4{animal}{end}(isnan(travellingWaveDir4{animal}{end})) = [];
      travellingWaveTimes4{animal}{end}(isnan(travellingWaveTimes4{animal}{end})) = [];
      assert(numel(travellingWaveDir4{animal}{end}) == numel(travellingWaveTimes4{animal}{end}));

      travellingWaveLocs4_postprocess{animal}{numel(travellingWaveLocs4_postprocess{animal})+1} = false(size(travellingWaveLocs4'));
      for iWave = 1:nTravellingWaves4
        indsOI = travellingWaveInds(iWave,1):travellingWaveInds(iWave,2);
        if sum(netWaveDir(indsOI))
          travellingWaveLocs4_postprocess{animal}{end}(indsOI) = true;
        end
      end

      %figure; plot(travellingWave.timestamps,pgdIndexSmooth); hold on
      %plot(travellingWave.timestamps,pgdIndexCutoff.*ones(size(pgdIndexSmooth)));
      %plot(travellingWave.timestamps,travellingWaveLocs4_postprocess{animal}{end}); hold off
      %ylim([-0.1 1.1]);
    else
      travellingWaveDir4{animal}{numel(travellingWaveDir4{animal})+1} = {};
      travellingWaveTimes4{animal}{numel(travellingWaveTimes4{animal})+1} = {};
      travellingWaveTimes4{animal}{numel(travellingWaveTimes4{animal})+1} = {};
      travellingWaveLocs4_postprocess{animal}{numel(travellingWaveLocs4_postprocess{animal})+1} = {};
    end

    if nTravellingWaves || nTravellingWaves4
      timestamps{animal}{numel(timestamps{animal})+1} = travellingWave.timestamps;
    else
      timestamps{animal}{numel(timestamps{animal})+1} = {};
    end

    % Save smoothed PGD Index and travelling wave locations based on this measure
    if options.saveData

      % PGD Index
      if strcmpi(options.dataType,'spikes')
        timeseriesFile = strrep(dataFiles{animal}{session}, '*', ...
          'pgdIndexSmoothed_travellingSpikingThetaWave.timeseries');
      elseif strcmpi(options.dataType,'lfp')
        timeseriesFile = strrep(dataFiles{animal}{session}, '*', ...
          'pgdIndexSmoothed_travellingLFPThetaWave.timeseries');
      end
      pgdIndexSmoothed.data = pgdIndexSmooth;
      pgdIndexSmoothed.timestamps = timestamps{animal}{end};
      pgdIndexSmoothed.precision = class(pgdIndexSmoothed.data);
      pgdIndexSmoothed.units = 'a.u.';
      pgdIndexSmoothed.nChannels = 1;
      pgdIndexSmoothed.channelNames = 'pgdIndexSmoothed';
      pgdIndexSmoothed.sr = round(1/samplingInterval);
      pgdIndexSmoothed.nSamples = numel(pgdIndexSmoothed.data);
      pgdIndexSmoothed.description = ...
        'Smoothed PGD Index produced using Matlab''s movmean function.';
      pgdIndexSmoothed.processingInfo.params = round(options.smoothingWindowSize/samplingInterval);
      pgdIndexSmoothed.processingInfo.function = 'movmean';
      pgdIndexSmoothed.processingInfo.date = datetime;
      pgdIndexSmoothed.processingInfo.username = getenv('username');
      pgdIndexSmoothed.processingInfo.hostname = getenv('computername');
      save(timeseriesFile, 'pgdIndexSmoothed', '-v7.3');

      % Travelling wave locations
      if strcmpi(options.dataType,'spikes')
        timeseriesFile = strrep(dataFiles{animal}{session}, '*', ...
          'pgdIndexSmoothedWaveLocs_travellingSpikingThetaWave.timeseries');
      elseif strcmpi(options.dataType,'lfp')
        timeseriesFile = strrep(dataFiles{animal}{session}, '*', ...
          'pgdIndexSmoothedWaveLocs_travellingLFPThetaWave.timeseries');
      end
      pgdIndexSmoothedWaveLocs.data = travellingWaveLocs4_postprocess{animal}{end};
      pgdIndexSmoothedWaveLocs.timestamps = timestamps{animal}{end};
      pgdIndexSmoothedWaveLocs.precision = class(pgdIndexSmoothedWaveLocs.data);
      pgdIndexSmoothedWaveLocs.units = 'boolean';
      pgdIndexSmoothedWaveLocs.nChannels = 1;
      pgdIndexSmoothedWaveLocs.channelNames = 'pgdIndexSmoothedWaveLocs';
      pgdIndexSmoothedWaveLocs.sr = round(1/samplingInterval);
      pgdIndexSmoothedWaveLocs.nSamples = numel(pgdIndexSmoothedWaveLocs.data);
      pgdIndexSmoothedWaveLocs.description = ...
        'Detected travelling wave locations based on smoothed PGD Index.';
      pgdIndexSmoothedWaveLocs.processingInfo.params = options;
      pgdIndexSmoothedWaveLocs.processingInfo.function = 'aggregateWaveData';
      pgdIndexSmoothedWaveLocs.processingInfo.date = datetime;
      pgdIndexSmoothedWaveLocs.processingInfo.username = getenv('username');
      pgdIndexSmoothedWaveLocs.processingInfo.hostname = getenv('computername');
      save(timeseriesFile, 'pgdIndexSmoothedWaveLocs', '-v7.3');
    end
  end
end

% Calculate PGD Index significance cutoff
pgdSignificanceCutoff = [0; 0];
for animal = 1:nAnimals
  pgdSignificanceCutoff(animal) = prctile(pgdIndexShuffled{animal},95);
end
end


function [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
  dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
  dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
  travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
  recSessions] = waveIncidenceRate(dataFiles, intervals, ...
  travellingWaveDir, travellingWaveDir2, travellingWaveDir3, ...
  travellingWaveTimes, travellingWaveDir4, travellingWaveTimes4, ...
  timestamps, options)
% [travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
%   dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
%   dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
%   travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4, ...
%   recSessions] = waveIncidenceRate(dataFiles, intervals, ...
%   travellingWaveDir, travellingWaveDir2, travellingWaveDir3, ...
%   travellingWaveTimes, travellingWaveDir4, travellingWaveTimes4, ...
%   timestamps, <options>)
%
% Function counts travelling waves with the same spread direction for
%   multiple recording sessions located in the ms-preliminary repository
%   and aggregates them over individual animals.
%
% Args:
%   dataFiles
%   intervals
%   travellingWaveDir
%   travellingWaveDir2
%   travellingWaveDir3
%   travellingWaveTimes
%   travellingWaveDir4
%   travellingWaveTimes4,
%   timestamps
%   <minIntervalLength>
%
% Returns:
%   travellingWaveIncidence
%   dorsoVentralWaveIncidence
%   ventroDorsalWaveIncidence
%   dorsoVentralWaveIncidence2
%   ventroDorsalWaveIncidence2
%   dorsoVentralWaveIncidence3
%   ventroDorsalWaveIncidence3
%   travellingWaveIncidence4
%   dorsoVentralWaveIncidence4
%   ventroDorsalWaveIncidence4
%   recSessions
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
  travellingWaveDir
  travellingWaveDir2
  travellingWaveDir3
  travellingWaveTimes
  travellingWaveDir4
  travellingWaveTimes4
  timestamps
  options.minIntervalLength = 0;
end

nAnimals = numel(dataFiles);
travellingWaveIncidence = cell(nAnimals,1);
dorsoVentralWaveIncidence = cell(nAnimals,1);
ventroDorsalWaveIncidence = cell(nAnimals,1);
dorsoVentralWaveIncidence2 = cell(nAnimals,1);
ventroDorsalWaveIncidence2 = cell(nAnimals,1);
dorsoVentralWaveIncidence3 = cell(nAnimals,1);
ventroDorsalWaveIncidence3 = cell(nAnimals,1);
travellingWaveIncidence4 = cell(nAnimals,1);
dorsoVentralWaveIncidence4 = cell(nAnimals,1);
ventroDorsalWaveIncidence4 = cell(nAnimals,1);
recSessions = cell(nAnimals,1);
for animal = 1:nAnimals
  for session = 1:numel(dataFiles{animal})
    if ~isempty(travellingWaveTimes{animal}{session}) && ~isempty(intervals{animal}{session})

      % Exclude data outside the periods of interest
      [~, indsOI] = selectArrayValues( ...
        travellingWaveTimes{animal}{session}, intervals{animal}{session});
      samplingInterval = timestamps{animal}{session}(2) - timestamps{animal}{session}(1);
      if numel(indsOI)*samplingInterval < options.minIntervalLength
        continue
      end
      selectData = travellingWaveDir{animal}{session}(indsOI);
      selectData2 = travellingWaveDir2{animal}{session}(indsOI);
      selectData3 = travellingWaveDir3{animal}{session}(indsOI);

      [~, indsOI4] = selectArrayValues( ...
        travellingWaveTimes4{animal}{session}, intervals{animal}{session});
      selectData4 = travellingWaveDir4{animal}{session}(indsOI4);

      % Count direction data
      netWaveDirHist = hist(selectData, sort(unique(selectData(~isnan(selectData))))); %#ok<*HIST>
      if isempty(selectData2(~isnan(selectData2)))
        netWaveDirHist2 = [0 0];
      elseif numel(unique(selectData2(~isnan(selectData2)))) == 1
        netWaveDirHist2 = hist(selectData2, ...
          [-abs(unique(selectData2(~isnan(selectData2)))) abs(unique(selectData2(~isnan(selectData2))))]);
      else
        netWaveDirHist2 = hist(selectData2, sort(unique(selectData2(~isnan(selectData2)))));
      end
      if isempty(selectData3(~isnan(selectData3)))
        netWaveDirHist3 = [0 0];
      elseif numel(unique(selectData3(~isnan(selectData3)))) == 1
        netWaveDirHist3 = hist(selectData3, ...
          [-abs(unique(selectData3(~isnan(selectData3)))) abs(unique(selectData3(~isnan(selectData3))))]);
      else
        netWaveDirHist3 = hist(selectData3, sort(unique(selectData3(~isnan(selectData3)))));
      end
      if isempty(selectData4(~isnan(selectData4)))
        netWaveDirHist4 = [0 0];
      elseif numel(unique(selectData4(~isnan(selectData4)))) == 1
        netWaveDirHist4 = hist(selectData4, ...
          [-abs(unique(selectData4(~isnan(selectData4)))) abs(unique(selectData4(~isnan(selectData4))))]);
      else
        netWaveDirHist4 = hist(selectData4, sort(unique(selectData4(~isnan(selectData4)))));
      end

      % Aggregate incidence data
      nTravellingWaves = numel(selectData);
      nTravellingWaves4 = numel(selectData4);
      timestampsOI = ...
        selectArrayValues(timestamps{animal}{session}, intervals{animal}{session});
      dataDuration = numel(timestampsOI)*samplingInterval;
      travellingWaveIncidence{animal} = [travellingWaveIncidence{animal}; ...
        nTravellingWaves/dataDuration];
      travellingWaveIncidence4{animal} = [travellingWaveIncidence4{animal}; ...
        nTravellingWaves4/dataDuration];
      % type 1
      if ~isempty(netWaveDirHist(1))
        dorsoVentralWaveIncidence{animal} = [dorsoVentralWaveIncidence{animal}; ...
          netWaveDirHist(1)/dataDuration];
      else
        dorsoVentralWaveIncidence{animal} = [dorsoVentralWaveIncidence{animal}; 0];
      end
      if ~isempty(netWaveDirHist(2))
        ventroDorsalWaveIncidence{animal} = [ventroDorsalWaveIncidence{animal}; ...
          netWaveDirHist(2)/dataDuration];
      else
        dorsoVentralWaveIncidence{animal} = [dorsoVentralWaveIncidence{animal}; 0];
      end
      % type 2
      if ~isempty(netWaveDirHist2(1))
        dorsoVentralWaveIncidence2{animal} = [dorsoVentralWaveIncidence2{animal}; ...
          netWaveDirHist2(1)/dataDuration];
      else
        dorsoVentralWaveIncidence2{animal} = [dorsoVentralWaveIncidence2{animal}; 0];
      end
      if ~isempty(netWaveDirHist2(2))
        ventroDorsalWaveIncidence2{animal} = [ventroDorsalWaveIncidence2{animal}; ...
          netWaveDirHist2(2)/dataDuration];
      else
        dorsoVentralWaveIncidence2{animal} = [dorsoVentralWaveIncidence2{animal}; 0];
      end
      % type 3
      if ~isempty(netWaveDirHist3(1))
        dorsoVentralWaveIncidence3{animal} = [dorsoVentralWaveIncidence3{animal}; ...
          netWaveDirHist3(1)/dataDuration];
      else
        dorsoVentralWaveIncidence3{animal} = [dorsoVentralWaveIncidence3{animal}; 0];
      end
      if ~isempty(netWaveDirHist3(2))
        ventroDorsalWaveIncidence3{animal} = [ventroDorsalWaveIncidence3{animal}; ...
          netWaveDirHist3(2)/dataDuration];
      else
        dorsoVentralWaveIncidence3{animal} = [dorsoVentralWaveIncidence3{animal}; 0];
      end
      % type 4
      if ~isempty(netWaveDirHist4(1))
        dorsoVentralWaveIncidence4{animal} = [dorsoVentralWaveIncidence4{animal}; ...
          netWaveDirHist4(1)/dataDuration];
      else
        dorsoVentralWaveIncidence4{animal} = [dorsoVentralWaveIncidence4{animal}; 0];
      end
      if ~isempty(netWaveDirHist4(2))
        ventroDorsalWaveIncidence4{animal} = [ventroDorsalWaveIncidence4{animal}; ...
          netWaveDirHist4(2)/dataDuration];
      else
        dorsoVentralWaveIncidence4{animal} = [dorsoVentralWaveIncidence4{animal}; 0];
      end
      recSessions{animal} = [recSessions{animal}; fileparts(dataFiles{animal}{session})];
    end
  end
end
end


function plotIncidenceRates(figFolder, figFilenameSuffix, ...
  travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
  dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
  dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
  travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4)
% plotIncidenceRates(figFolder, figFilenameSuffix, ...
%   travellingWaveIncidence, dorsoVentralWaveIncidence, ventroDorsalWaveIncidence, ...
%   dorsoVentralWaveIncidence2, ventroDorsalWaveIncidence2, ...
%   dorsoVentralWaveIncidence3, ventroDorsalWaveIncidence3, ...
%   travellingWaveIncidence4, dorsoVentralWaveIncidence4, ventroDorsalWaveIncidence4)
%
% Function plots and saves travelling wave incidence rates across sessions
%   for the ms-preliminary repository.
%
% Args:
%   figFolder
%   figFilenameSuffix
%   travellingWaveIncidence
%   dorsoVentralWaveIncidence
%   ventroDorsalWaveIncidence
%   dorsoVentralWaveIncidence2
%   ventroDorsalWaveIncidence2
%   dorsoVentralWaveIncidence3
%   ventroDorsalWaveIncidence3
%   travellingWaveIncidence4
%   dorsoVentralWaveIncidence4
%   ventroDorsalWaveIncidence4
%
% Returns:
%   None
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

arguments
  figFolder
  figFilenameSuffix
  travellingWaveIncidence
  dorsoVentralWaveIncidence
  ventroDorsalWaveIncidence
  dorsoVentralWaveIncidence2
  ventroDorsalWaveIncidence2
  dorsoVentralWaveIncidence3
  ventroDorsalWaveIncidence3
  travellingWaveIncidence4
  dorsoVentralWaveIncidence4
  ventroDorsalWaveIncidence4
end

for animal = 1:numel(travellingWaveIncidence)
  % Output folder
  animalID = ['PP0' num2str(animal)];
  figFolder_animal = fullfile(figFolder, animalID);
  if ~exist(figFolder_animal,'dir')
    mkdir(figFolder_animal);
  end

  % Type 1: Direction at max PGD Index values only (legacy)
  markerSize = 20;
  lineWidth = 2;
  fH = figure;
  if numel(dorsoVentralWaveIncidence{animal}) == 1
    plot(dorsoVentralWaveIncidence{animal}, '.', 'MarkerSize',markerSize); hold on
    plot(ventroDorsalWaveIncidence{animal}, '.', 'MarkerSize',markerSize);
    plot(travellingWaveIncidence{animal}, '.', 'MarkerSize',markerSize); hold off
  else
    plot(dorsoVentralWaveIncidence{animal}, 'LineWidth',lineWidth); hold on
    plot(ventroDorsalWaveIncidence{animal}, 'LineWidth',lineWidth);
    plot(travellingWaveIncidence{animal}, 'LineWidth',lineWidth); hold off
  end
  yLim = ylim;
  ylim([0 yLim(2)]);
  xlabel('Recording session #');
  ylabel('Incidence rate (per second)');
  legend('Ventro-dorsal','Dorso-ventral','Total');
  legend('boxoff');
  title(['Trvl Wave Inc (max PGDI,' figFilenameSuffix '): ' animalID], ...
    'Interpreter','none');
  set(fH, 'Name',['travelling_wave_incidence_rates_' animalID]);
  figFilename = fullfile(figFolder_animal, ...
    ['travellingWaveIncidenceRatesMax' figFilenameSuffix '.fig']);
  savefig(fH,figFilename,'compact');
  title('');
  saveas(fH,figFilename(1:end-4),'png');
  close(fH);

  % Type 2: Net direction over full cycle (legacy)
  fH = figure;
  if numel(dorsoVentralWaveIncidence2{animal}) == 1
    plot(dorsoVentralWaveIncidence2{animal}, '.', 'MarkerSize',markerSize); hold on
    plot(ventroDorsalWaveIncidence2{animal}, '.', 'MarkerSize',markerSize);
    plot(travellingWaveIncidence{animal}, '.', 'MarkerSize',markerSize); hold off
  else
    plot(dorsoVentralWaveIncidence2{animal}, 'LineWidth',lineWidth); hold on
    plot(ventroDorsalWaveIncidence2{animal}, 'LineWidth',lineWidth);
    plot(travellingWaveIncidence{animal}, 'LineWidth',lineWidth); hold off
  end
  yLim = ylim;
  ylim([0 yLim(2)]);
  xlabel('Recording session #');
  ylabel('Incidence rate (per second)');
  legend('Ventro-dorsal','Dorso-ventral','Total');
  legend('boxoff');
  title(['Trvl Wave Inc (hiPGDI cycle,' figFilenameSuffix '): ' animalID], ...
    'Interpreter','none');
  set(fH, 'Name',['travelling_wave_incidence_rates_' animalID]);
  figFilename = fullfile(figFolder_animal, ...
    ['travellingWaveIncidenceRatesCycle' figFilenameSuffix '.fig']);
  savefig(fH,figFilename,'compact');
  title('');
  saveas(fH,figFilename(1:end-4),'png');
  close(fH);

  % Type 3: Net direction over the travelling part of the cycle only (legacy)
  fH = figure;
  if numel(dorsoVentralWaveIncidence3{animal}) == 1
    plot(dorsoVentralWaveIncidence3{animal}, '.', 'MarkerSize',markerSize); hold on
    plot(ventroDorsalWaveIncidence3{animal}, '.', 'MarkerSize',markerSize);
    plot(travellingWaveIncidence{animal}, '.', 'MarkerSize',markerSize); hold off
  else
    plot(dorsoVentralWaveIncidence3{animal}, 'LineWidth',lineWidth); hold on
    plot(ventroDorsalWaveIncidence3{animal}, 'LineWidth',lineWidth);
    plot(travellingWaveIncidence{animal}, 'LineWidth',lineWidth); hold off
  end
  yLim = ylim;
  ylim([0 yLim(2)]);
  xlabel('Recording session #');
  ylabel('Incidence rate (per second)');
  legend('Ventro-dorsal','Dorso-ventral','Total');
  legend('boxoff');
  title(['Trvl Wave Inc (hiPGDI,' figFilenameSuffix '): ' animalID], ...
    'Interpreter','none');
  set(fH, 'Name',['travelling_wave_incidence_rates_' animalID]);
  figFilename = fullfile(figFolder_animal, ...
    ['travellingWaveIncidenceRatesWave' figFilenameSuffix '.fig']);
  savefig(fH,figFilename,'compact');
  title('');
  saveas(fH,figFilename(1:end-4),'png');
  close(fH);

  % Type 4: Net direction over the travelling part of the cycle only
  fH = figure;
  if numel(dorsoVentralWaveIncidence4{animal}) == 1
    plot(dorsoVentralWaveIncidence4{animal}, '.', 'MarkerSize',markerSize); hold on
    plot(ventroDorsalWaveIncidence4{animal}, '.', 'MarkerSize',markerSize);
    plot(travellingWaveIncidence4{animal}, '.', 'MarkerSize',markerSize); hold off
  else
    plot(dorsoVentralWaveIncidence4{animal}, 'LineWidth',lineWidth); hold on
    plot(ventroDorsalWaveIncidence4{animal}, 'LineWidth',lineWidth);
    plot(travellingWaveIncidence4{animal}, 'LineWidth',lineWidth); hold off
  end
  yLim = ylim;
  ylim([0 yLim(2)]);
  xlabel('Recording session #');
  ylabel('Incidence rate (per second)');
  legend('Ventro-dorsal','Dorso-ventral','Total');
  legend('boxoff');
  title(['Smth Trvl Wave Inc (hiPGDI,' figFilenameSuffix '): ' animalID], ...
    'Interpreter','none');
  set(fH, 'Name',['smoothed_travelling_wave_incidence_rates_' animalID]);
  figFilename = fullfile(figFolder_animal, ...
    ['smoothedTravellingWaveIncidenceRatesWave' figFilenameSuffix '.fig']);
  savefig(fH,figFilename,'compact');
  title('');
  saveas(fH,figFilename(1:end-4),'png');
  close(fH);
end
end