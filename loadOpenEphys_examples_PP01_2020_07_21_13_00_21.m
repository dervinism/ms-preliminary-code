basepath = 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data\PP01\PP01_2020-07-21_13-00-21';

basename = basenameFromBasepath(basepath);
cd(basepath)
session = loadSession(basepath,basename,'showGUI',false); % Loading session info


%% Detect open ephys data

session = preprocessOpenEphysData('session',session,'showGUI',false);


%% Load OpenEphys Settings File

session = loadOpenEphysSettingsFile(session);
% saveStruct(session);


% %% Load digital pulses
% digSeriesFilename = fullfile(basepath,'openephysDig.digitalseries.mat');
% if exist(digSeriesFilename, 'file')
%   load(digSeriesFilename);
% else
%   openephysDig = loadOpenEphysDigital(session); %#ok<NASGU> 
% end
% 
% 
% %% Behavior processing
% scaling_factor = 0.5;
% offset_origin = [5,-5,0];
% 
% offset_rigid_body = [5,-5,0]; % Not implemented yet
% circular_track = loadOptitrack('session',session,'dataName','circular_track','offset_origin',offset_origin,'scaling_factor',scaling_factor);
% 
% % linear_track = loadOptitrack('session',session,'dataName','linear_track','offset',offset,'scaling_factor',scaling_factor);
% 
% 
% %% maze parameters
% maze = {};
% maze.type = 'theta';
% maze.radius_in = 96.5/2; % From the centre of the central arm to the inner circle wall of the maze
% maze.radius_out =  116.5/2; % From the centre of the central arm to the outer circle wall of the maze
% maze.arm_half_width = 4; % Central arm half-width
% maze.cross_radii = 47.9; % ???
% maze.rim_buffer = 10; % rim width
% maze.polar_rho_limits = [44,75]; % rim limits for assigning trials
% maze.polar_theta_limits = [15,2.8*maze.radius_in]; % In units of cm % ???
% maze.pos_x_limits = [-12,12]; % arm delimiters in x-coords
% maze.pos_y_limits = [-36,36]; % arm delimiters in y-coords
% 
% subplot(1,2,1)
% plot_ThetaMaze(maze)
% 
% %% Align optitrack behavior data with TTL pulses
% 
% % Loading pulses
% openephysDig = loadStruct('openephysDig','digitalseries','session',session);
% 
% if session.behavioralTracking{1}.epoch==1
%     TTL_offset = 0;
% else
%     TTL_offset = sum(openephysDig.nOnPrFile(1:session.behavioralTracking{1}.epoch-1));    
% end
% circular_track.timestamps = openephysDig.on{1}((1:circular_track.nSamples)+TTL_offset);
% circular_track.timestamps_reference = 'ephys';
% saveStruct(circular_track,'behavior','session',session);
% 
% 
% %% Get trials from behavior
% 
% % Define the trials struct:
% circular_track = getTrials_thetamaze(circular_track,maze, 1);
% 
% 
% %% Linearizing and defining boundaries
% 
% circular_track = linearize_theta_maze(circular_track,maze);
% 
% % Setting a minimum speed threshold
% circular_track.speed_th = 10;
% 
% % Generating left_right states data
% circular_track.states.left_right = nan(size(circular_track.trials));
% for i = 1:circular_track.trials.alternation.nTrials
%     circular_track.states.left_right(circular_track.trials.alternation.trials==i) = circular_track.states.left_right(i);
% end
% circular_track.stateNames.left_right = {'Left','Right'};
% 
% saveStruct(circular_track,'behavior','session',session);
% trials = circular_track.trials;
% saveStruct(trials,'behavior','session',session);
% 
% % After this you can load the generated files:
% % circular_track = loadStruct('circular_track','behavior','session',session);
% % trials = loadStruct('trials','behavior','session',session);


%% Plot spikes
%spikesUnits = loadSpikes('session',session, 'format','CellExplorer', 'labelsToRead',{'good'}, 'saveMat',true);
spikesUnits = load(fullfile(session.general.basePath, [session.general.name '.spikes.cellinfo.mat']));
spikesUnits = spikesUnits.spikes;

% % A figure for each unit showing the X,Y position of the spikes
% for unit1 = 1:spikesUnits.numcells
%     spikesUnits.position_x{unit1} = interp1(circular_track.timestamps,circular_track.position.x,spikesUnits.times{unit1},'linear');
%     spikesUnits.position_y{unit1} = interp1(circular_track.timestamps,circular_track.position.y,spikesUnits.times{unit1},'linear');
%     
%     figure
%     plot(circular_track.position.x,circular_track.position.y,'-k'), hold on
%     plot(spikesUnits.position_x{unit1},spikesUnits.position_y{unit1},'.r'), hold on
%     plot_ThetaMaze(maze), title(['Unit ' num2str(unit1)]), xlabel('X'), ylabel('Y')
% end
% 
% for unit1 = 1:spikesUnits.numcells
%     spikesUnits.position_linearized{unit1} = interp1(circular_track.timestamps,circular_track.position.linearized,spikesUnits.times{unit1},'linear');
%     spikesUnits.position_trials{unit1} = interp1(circular_track.timestamps,circular_track.trials.alternation.trials,spikesUnits.times{unit1},'linear');
%     
%     figure
%     plot(circular_track.position.linearized,circular_track.trials.alternation.trials,'.k'), hold on
%     plot(spikesUnits.position_linearized{unit1},spikesUnits.position_trials{unit1},'.r'), hold on
%     title(['Unit ' num2str(unit1)]), xlabel('Linearized position'), ylabel('Trials')
% end

spikesUnits.spindices = generateSpinDices(spikesUnits.times);
% figure, plot(spikesUnits.spindices(:,1),spikesUnits.spindices(:,2),'|')
% 
% for i = 1:circular_track.trials.alternation.nTrials
%     trial1 = i;
%     figure
%     idx = circular_track.trials.alternation.trials == trial1;
%     plot(circular_track.position.x(idx),circular_track.position.y(idx),'-k'), hold on
%     for unit1 = 1:spikesUnits.numcells
%         idx2 = spikesUnits.position_trials{unit1} == trial1;
%         plot(spikesUnits.position_x{unit1}(idx2),spikesUnits.position_y{unit1}(idx2),'.'), hold on
%     end
%     title(['Trial ' num2str(trial1)]), xlabel('X'), ylabel('Y')
% end


%% Get convolved population rate
%muaSpikes = loadSpikes('session',session, 'format','Phy', 'labelsToRead',{'mua'}, 'saveMat',false);
load(fullfile(session.general.basePath, [session.general.name '.muaSpikes.cellinfo.mat']));
%save([basename '.muaSpikes.cellinfo.mat'], 'muaSpikes', '-v7.3');

populationRate.numcells = 1;
populationRate.times = [spikesUnits.times, muaSpikes.times];
populationRate.total = numel(populationRate.times);
populationRate.times = {sort(concatenateCells(populationRate.times, 1))};
populationRate.basename = basename;
populationRate.description = 'Spike times across the entire recording session including all units and MUAs.';
% save([basename '.populationRate.cellinfo.mat'], 'populationRate', '-v7.3');
% 
% [spikes_presentation,time_bins,parameters] = spikes_convolution(populationRate, 0.002, 25);
% figure; plot(time_bins, spikes_presentation);
% convolvedPopulationRate.data = spikes_presentation';
% convolvedPopulationRate.timestamps = time_bins';
% convolvedPopulationRate.precision = class(spikes_presentation);
% convolvedPopulationRate.units = 'spikes';
% convolvedPopulationRate.nChannels = 1;
% convolvedPopulationRate.sr = parameters.stepsize;
% convolvedPopulationRate.nSamples = numel(spikes_presentation);
% convolvedPopulationRate.description = 'Population firing rate convolved with a Gaussian function.';
% convolvedPopulationRate.processingInfo.params = parameters;
% convolvedPopulationRate.processingInfo.function = 'petersen-lab-matlab/spikes/spikes_convolution';
% convolvedPopulationRate.processingInfo.date = datetime;
% convolvedPopulationRate.processingInfo.username = getenv('username');
% convolvedPopulationRate.processingInfo.hostname = getenv('computername');
% save([basename '.convolvedPopulationRate.timeseries.mat'], 'convolvedPopulationRate', '-v7.3');


% Get population rate delta phase and power
[instantDelta, parameters] = freqBandPropertiesForPointProcess(populationRate.times{1}, [1 4], pad=3);

% delta frequency range spectrogram
deltaSpectrogram.data = instantDelta.spectrogram';
deltaSpectrogram.timestamps = instantDelta.spectrogramTimestamps';
deltaSpectrogram.precision = class(instantDelta.spectrogram);
deltaSpectrogram.units = 'a.u.';
deltaSpectrogram.nChannels = numel(instantDelta.spectrogramFrequencies);
deltaSpectrogram.channelNames = instantDelta.spectrogramFrequencies';
deltaSpectrogram.sr = round(1/(instantDelta.spectrogramTimestamps(2)-instantDelta.spectrogramTimestamps(1)));
deltaSpectrogram.nSamples = numel(instantDelta.spectrogramTimestamps);
deltaSpectrogram.description = ['A population rate delta frequency range (1-4 Hz) spectrogram. ' ...
                                'Data columns correspond to individual frequencies given by channelNames.'];
deltaSpectrogram.processingInfo.params = parameters;
deltaSpectrogram.processingInfo.function = 'petersen-lab-matlab/waves/freqBandPropertiesForPointProcess';
deltaSpectrogram.processingInfo.date = datetime;
deltaSpectrogram.processingInfo.username = getenv('username');
deltaSpectrogram.processingInfo.hostname = getenv('computername');
save([basename '.deltaSpectrogram.timeseries.mat'], 'deltaSpectrogram', '-v7.3');

% Delta frequency range power
deltaPower.data = [instantDelta.spectrogramMaxFrequency; instantDelta.spectrogramPower; ...
                   instantDelta.moderatePowerPeriods; instantDelta.highPowerPeriods]';
deltaPower.timestamps = instantDelta.spectrogramTimestamps';
deltaPower.precision = class(instantDelta.spectrogramPower);
deltaPower.units = {'Hz', 'a.u.', 'None', 'None'};
deltaPower.nChannels = size(deltaPower.data,2);
deltaPower.channelNames = {'spectrogramMaxFrequency' 'spectrogramPower' ...
                           'moderatePowerPeriods' 'highPowerPeriods'};
deltaPower.sr = round(1/(instantDelta.spectrogramTimestamps(2)-instantDelta.spectrogramTimestamps(1)));
deltaPower.nSamples = numel(instantDelta.spectrogramTimestamps);
deltaPower.description = ['A population rate delta frequency range (1-4 Hz) maximum power (the 2nd column of the data matrix). ' ...
                          'The data of this timeseries collection also contains other information, like maximum delta frequency ' ...
                          '(the 1st data column), a logical vector with true values representing points when the delta power is moderate ' ...
                          '(above the mean power but below mean + norminv(0.95)*SD; 3rd data column) and a logical vector with true values ' ...
                          'representing points when the delta power is high (above mean + norminv(0.95)*SD; 4th data column).'];
deltaPower.processingInfo.params = parameters;
deltaPower.processingInfo.function = 'petersen-lab-matlab/waves/freqBandPropertiesForPointProcess';
deltaPower.processingInfo.date = datetime;
deltaPower.processingInfo.username = getenv('username');
deltaPower.processingInfo.hostname = getenv('computername');
save([basename '.deltaPower.timeseriesCollection.mat'], 'deltaPower', '-v7.3');

% Instantaneous delta frequency
instantDeltaFrequency.data = instantDelta.instFrequency';
instantDeltaFrequency.timestamps = instantDelta.instTimestamps';
instantDeltaFrequency.precision = class(instantDelta.instFrequency);
instantDeltaFrequency.units = 'Hz';
instantDeltaFrequency.nChannels = 1;
instantDeltaFrequency.sr = round(1/(instantDelta.instTimestamps(2)-instantDelta.instTimestamps(1)));
instantDeltaFrequency.nSamples = numel(instantDelta.instTimestamps);
instantDeltaFrequency.description = ['A vector of instantaneous delta frequencies. The number of samples ' ...
                                     'corresponds to the number of delta oscillation cycles derived using Hilbert transform.'];
instantDeltaFrequency.processingInfo.params = parameters;
instantDeltaFrequency.processingInfo.function = 'petersen-lab-matlab/waves/freqBandPropertiesForPointProcess';
instantDeltaFrequency.processingInfo.date = datetime;
instantDeltaFrequency.processingInfo.username = getenv('username');
instantDeltaFrequency.processingInfo.hostname = getenv('computername');
save([basename '.instantDeltaFrequency.timeseries.mat'], 'instantDeltaFrequency', '-v7.3');

% Instantaneous delta phase
instantDeltaPhase.data = [instantDelta.instPhase; instantDelta.instPhaseUnwrapped; instantDelta.spikingRateFiltered]';
instantDeltaPhase.timestamps = instantDelta.instPhaseTimestamps';
instantDeltaPhase.precision = class(instantDelta.instPhase);
instantDeltaPhase.units = {'rad', 'rad', 'spikes/s'};
instantDeltaPhase.nChannels = size(instantDeltaPhase.data,2);
instantDeltaPhase.channelNames = {'instantPhase' 'instantUnwrappedPhase' 'convolvedFilteredSpikes'};
instantDeltaPhase.sr = round(1/(instantDelta.spectrogramTimestamps(2)-instantDelta.spectrogramTimestamps(1)));
instantDeltaPhase.nSamples = numel(instantDelta.spectrogramTimestamps);
instantDeltaPhase.description = ['The 1st data column contains instantaneous delta oscillation phase estimated using Hilbert transform. ' ...
                                 'The number of samples corresponds to the length of the convolvedFilteredSpikes vector (3rd data column).' ...
                                 'The 2nd data column contains the output of unwrap(1st-data-column).' ...
                                 'The 3rd data column contains the convolved population firing rate band-passed filterd at 1-4 Hz.'];
instantDeltaPhase.processingInfo.params = parameters;
instantDeltaPhase.processingInfo.function = 'petersen-lab-matlab/waves/freqBandPropertiesForPointProcess';
instantDeltaPhase.processingInfo.date = datetime;
instantDeltaPhase.processingInfo.username = getenv('username');
instantDeltaPhase.processingInfo.hostname = getenv('computername');
save([basename '.instantDeltaPhase.timeseriesCollection.mat'], 'instantDeltaPhase', '-v7.3');

% Delta amplitude
deltaAmplitude.data = instantDelta.amplitude';
deltaAmplitude.timestamps = instantDelta.instPhaseTimestamps';
deltaAmplitude.precision = class(instantDelta.amplitude);
deltaAmplitude.units = 'Hz';
deltaAmplitude.nChannels = 1;
deltaAmplitude.sr = round(1/(instantDelta.instPhaseTimestamps(2)-instantDelta.instPhaseTimestamps(1)));
deltaAmplitude.nSamples = numel(instantDelta.instPhaseTimestamps);
deltaAmplitude.description = 'Delta oscillation amplitude vector: Derived using Hilbert transform.';
deltaAmplitude.processingInfo.params = parameters;
deltaAmplitude.processingInfo.function = 'petersen-lab-matlab/waves/freqBandPropertiesForPointProcess';
deltaAmplitude.processingInfo.date = datetime;
deltaAmplitude.processingInfo.username = getenv('username');
deltaAmplitude.processingInfo.hostname = getenv('computername');
save([basename '.deltaAmplitude.timeseries.mat'], 'deltaAmplitude', '-v7.3');

% % Get population rate theta phase and power
% [instantTheta, parameters] = instThetaForPointProcess(populationRate.times{1}, pad=3, showPower=true);
% 
% % % Theta frequency range spectrogram
% % thetaSpectrogram.data = instantTheta.spectrogram';
% % thetaSpectrogram.timestamps = instantTheta.spectrogramTimestamps';
% % thetaSpectrogram.precision = class(instantTheta.spectrogram);
% % thetaSpectrogram.units = 'a.u.';
% % thetaSpectrogram.nChannels = numel(instantTheta.spectrogramFrequencies);
% % thetaSpectrogram.channelNames = instantTheta.spectrogramFrequencies';
% % thetaSpectrogram.sr = round(1/(instantTheta.spectrogramTimestamps(2)-instantTheta.spectrogramTimestamps(1)));
% % thetaSpectrogram.nSamples = numel(instantTheta.spectrogramTimestamps);
% % thetaSpectrogram.description = ['A population rate theta frequency range (4-12 Hz) spectrogram. ' ...
% %                                 'Data columns correspond to individual frequencies given by channelNames.'];
% % thetaSpectrogram.processingInfo.params = parameters;
% % thetaSpectrogram.processingInfo.function = 'petersen-lab-matlab/waves/instThetaForPointProcess';
% % thetaSpectrogram.processingInfo.date = datetime;
% % thetaSpectrogram.processingInfo.username = getenv('username');
% % thetaSpectrogram.processingInfo.hostname = getenv('computername');
% % save([basename '.thetaSpectrogram.timeseries.mat'], 'thetaSpectrogram', '-v7.3');
% % 
% % % Theta frequency range power
% % thetaPower.data = [instantTheta.spectrogramMaxFrequency; instantTheta.spectrogramPower; ...
% %                    instantTheta.moderatePowerPeriods; instantTheta.highPowerPeriods]';
% % thetaPower.timestamps = instantTheta.spectrogramTimestamps';
% % thetaPower.precision = class(instantTheta.spectrogramPower);
% % thetaPower.units = {'Hz', 'a.u.', 'None', 'None'};
% % thetaPower.nChannels = size(thetaPower.data,2);
% % thetaPower.channelNames = {'spectrogramMaxFrequency' 'spectrogramPower' ...
% %                            'moderatePowerPeriods' 'highPowerPeriods'};
% % thetaPower.sr = round(1/(instantTheta.spectrogramTimestamps(2)-instantTheta.spectrogramTimestamps(1)));
% % thetaPower.nSamples = numel(instantTheta.spectrogramTimestamps);
% % thetaPower.description = ['A population rate theta frequency range (4-12 Hz) maximum power (the 2nd column of the data matrix). ' ...
% %                           'The data of this timeseries collection also contains other information, like maximum theta frequency ' ...
% %                           '(the 1st data column), a logical vector with true values representing points when the theta power is moderate ' ...
% %                           '(above the mean power but below mean + norminv(0.95)*SD; 3rd data column) and a logical vector with true values ' ...
% %                           'representing points when the theta power is high (above mean + norminv(0.95)*SD; 4th data column).'];
% % thetaPower.processingInfo.params = parameters;
% % thetaPower.processingInfo.function = 'petersen-lab-matlab/waves/instThetaForPointProcess';
% % thetaPower.processingInfo.date = datetime;
% % thetaPower.processingInfo.username = getenv('username');
% % thetaPower.processingInfo.hostname = getenv('computername');
% % save([basename '.thetaPower.timeseriesCollection.mat'], 'thetaPower', '-v7.3');
% 
% % Instantaneous theta frequency
% instantThetaFrequency.data = instantTheta.instFrequency';
% instantThetaFrequency.timestamps = instantTheta.instTimestamps';
% instantThetaFrequency.precision = class(instantTheta.instFrequency);
% instantThetaFrequency.units = 'Hz';
% instantThetaFrequency.nChannels = 1;
% instantThetaFrequency.sr = round(1/(instantTheta.instTimestamps(2)-instantTheta.instTimestamps(1)));
% instantThetaFrequency.nSamples = numel(instantTheta.instTimestamps);
% instantThetaFrequency.description = ['A vector of instantaneous theta frequencies. The number of samples ' ...
%                                      'corresponds to the number of theta oscillation cycles derived using Hilbert transform.'];
% instantThetaFrequency.processingInfo.params = parameters;
% instantThetaFrequency.processingInfo.function = 'petersen-lab-matlab/waves/instThetaForPointProcess';
% instantThetaFrequency.processingInfo.date = datetime;
% instantThetaFrequency.processingInfo.username = getenv('username');
% instantThetaFrequency.processingInfo.hostname = getenv('computername');
% save([basename '.instantThetaFrequency.timeseries.mat'], 'instantThetaFrequency', '-v7.3');
% 
% % Instantaneous theta phase
% instantThetaPhase.data = [instantTheta.instPhase; instantTheta.instPhaseUnwrapped; instantTheta.spikingRateFiltered]';
% instantThetaPhase.timestamps = instantTheta.instPhaseTimestamps';
% instantThetaPhase.precision = class(instantTheta.instPhase);
% instantThetaPhase.units = {'rad', 'rad', 'spikes/s'};
% instantThetaPhase.nChannels = size(instantThetaPhase.data,2);
% instantThetaPhase.channelNames = {'instantPhase' 'instantUnwrappedPhase' 'convolvedFilteredSpikes'};
% instantThetaPhase.sr = round(1/(instantTheta.spectrogramTimestamps(2)-instantTheta.spectrogramTimestamps(1)));
% instantThetaPhase.nSamples = numel(instantTheta.spectrogramTimestamps);
% instantThetaPhase.description = ['The 1st data column contains instantaneous theta oscillation phase estimated using Hilbert transform. ' ...
%                                  'The number of samples corresponds to the length of the convolvedFilteredSpikes vector (3rd data column).' ...
%                                  'The 2nd data column contains the output of unwrap(1st-data-column).' ...
%                                  'The 3rd data column contains the convolved population firing rate band-passed filterd at 4-12 Hz.'];
% instantThetaPhase.processingInfo.params = parameters;
% instantThetaPhase.processingInfo.function = 'petersen-lab-matlab/waves/instThetaForPointProcess';
% instantThetaPhase.processingInfo.date = datetime;
% instantThetaPhase.processingInfo.username = getenv('username');
% instantThetaPhase.processingInfo.hostname = getenv('computername');
% save([basename '.instantThetaPhase.timeseriesCollection.mat'], 'instantThetaPhase', '-v7.3');
% 
% % Theta amplitude
% thetaAmplitude.data = instantTheta.amplitude';
% thetaAmplitude.timestamps = instantTheta.instPhaseTimestamps';
% thetaAmplitude.precision = class(instantTheta.amplitude);
% thetaAmplitude.units = 'Hz';
% thetaAmplitude.nChannels = 1;
% thetaAmplitude.sr = round(1/(instantTheta.instPhaseTimestamps(2)-instantTheta.instPhaseTimestamps(1)));
% thetaAmplitude.nSamples = numel(instantTheta.instPhaseTimestamps);
% thetaAmplitude.description = 'Theta oscillation amplitude vector: Derived using Hilbert transform.';
% thetaAmplitude.processingInfo.params = parameters;
% thetaAmplitude.processingInfo.function = 'petersen-lab-matlab/waves/instThetaForPointProcess';
% thetaAmplitude.processingInfo.date = datetime;
% thetaAmplitude.processingInfo.username = getenv('username');
% thetaAmplitude.processingInfo.hostname = getenv('computername');
% save([basename '.thetaAmplitude.timeseries.mat'], 'thetaAmplitude', '-v7.3');