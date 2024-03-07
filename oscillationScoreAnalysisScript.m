% This script creates oscillation score figures for ms-preliminary dataset


%% Initialise parameters
% Data setup parameters
params


%% Plot oscillation score and firing rate recording channel profiles
nAnimals = numel(dataFiles);
oscScoreAgregate = cell(nAnimals,1);
oscScoreShuffledAgregate = cell(nAnimals,1);
maxWaveformChAgregate = cell(nAnimals,1);
firingRatesAgregate = cell(nAnimals,1);
fH_OSvCh = zeros(nAnimals,1);
fH_uFRvCh = zeros(nAnimals,1);
fH_chFRvCh = zeros(nAnimals,1);
fH_OSvFR = zeros(nAnimals,1);
p_OSvCh = cell(nAnimals,1);
p_uFRvCh = cell(nAnimals,1);
p_chFRvCh = cell(nAnimals,1);
p_OSvFR = cell(nAnimals,1);
l_OSvCh = cell(nAnimals,1);
l_uFRvCh = cell(nAnimals,1);
l_chFRvCh = cell(nAnimals,1);
l_OSvFR = cell(nAnimals,1);
for animal = 1:nAnimals
  for session = 1:numel(dataFiles{animal})
    if session == 1
      fH_OSvCh(animal) = figure;
      fH_uFRvCh(animal) = figure;
      fH_chFRvCh(animal) = figure;
      fH_OSvFR(animal) = figure;
    end

    % Estimate and draw the oscillation score significance threshold
    if session == numel(dataFiles{animal}) && ~isempty(oscScoreAgregate{animal})
      figure(fH_OSvCh(animal)); hold on
      significanceCutoff = prctile(oscScoreShuffledAgregate{animal},95);
      xLim = xlim;
      p = plot(xLim, [significanceCutoff significanceCutoff], 'r:'); hold off
      p_OSvCh{animal} = [p_OSvCh{animal}; p];
      l_OSvCh{animal}{numel(l_OSvCh{animal})+1} = 'cutoff';
    end

    % Load spiking data
    spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
    if ~exist(spikesFile, 'file')
      continue
    end
    load(spikesFile); %#ok<*LOAD>
    spikesFile = strrep(dataFiles{animal}{session}, '*', 'muaSpikes.cellinfo');
    load(spikesFile);

    % Load oscillation score data
    oscScoreFile = strrep(dataFiles{animal}{session}, '*', 'oscillationScore.cellinfo');
    load(oscScoreFile);
    oscScoreFile = strrep(dataFiles{animal}{session}, '*', 'oscillationScoreShuffled.cellinfo');
    load(oscScoreFile);

    % Calculate unit firing rates
    lastSpikeTime = max(cellfun(@(x) getMaxSpikeTime(x), spikes.times));
    unitFiringRates = spikes.total/lastSpikeTime;

    % Calculate channel firing rates
    chSpikeTimes_units = getChannelSpikeTimes(spikes.times, spikes.maxWaveformCh1, 384);
    chSpikeTimes_muas = getChannelSpikeTimes(muaSpikes.times, muaSpikes.maxWaveformCh1, 384);
    chSpikeCounts_units = cellfun(@(x) numel(x), chSpikeTimes_units);
    chSpikeCounts_muas = cellfun(@(x) numel(x), chSpikeTimes_muas);
    chFiringRates_units = chSpikeCounts_units./lastSpikeTime;
    chSpikeCounts_muas = chSpikeCounts_muas./lastSpikeTime;
    chFiringRates = chFiringRates_units + chSpikeCounts_muas;

    % Agregate the data
    oscScoreAgregate{animal} = [oscScoreAgregate{animal}; oscillationScore.data(:,1)];
    oscScoreShuffledAgregate{animal} = [oscScoreShuffledAgregate{animal}; ...
      oscillationScoreShuffled.data(:,1)];
    chOrder = channelOrder{animal}(session, spikes.maxWaveformCh1);
    maxWaveformChAgregate{animal} = [maxWaveformChAgregate{animal}; chOrder'];
    firingRatesAgregate{animal} = [firingRatesAgregate{animal}; unitFiringRates'];

    % Plot the oscillation score channel profile
    figure(fH_OSvCh(animal)); hold on
    p = plot(chOrder, oscillationScore.data(:,1), '.', 'MarkerSize',10); hold off
    p_OSvCh{animal} = [p_OSvCh{animal}; p];
    if isempty(l_OSvCh{animal})
      l_OSvCh{animal} = {['s' num2str(session)]};
    else
      l_OSvCh{animal}{numel(l_OSvCh{animal})+1} = ['s' num2str(session)];
    end

    % Plot the unit firing rate profile
    figure(fH_uFRvCh(animal));
    if ~isempty(get(fH_uFRvCh(animal),'Children'))
      hold on
    end
    p = semilogy(chOrder, unitFiringRates, '.', 'MarkerSize',10); hold off
    p_uFRvCh{animal} = [p_uFRvCh{animal}; p];
    if isempty(l_uFRvCh{animal})
      l_uFRvCh{animal} = {['s' num2str(session)]};
    else
      l_uFRvCh{animal}{numel(l_uFRvCh{animal})+1} = ['s' num2str(session)];
    end

    % Plot the channel firing rate profile
    figure(fH_chFRvCh(animal));
    if ~isempty(get(fH_chFRvCh(animal),'Children'))
      hold on
    end
    p = semilogy(channelOrder{animal}(session,logical(chFiringRates)), ...
      chFiringRates(logical(chFiringRates)), '.', 'MarkerSize',10); hold off
    p_chFRvCh{animal} = [p_chFRvCh{animal}; p];
    if isempty(l_chFRvCh{animal})
      l_chFRvCh{animal} = {['s' num2str(session)]};
    else
      l_chFRvCh{animal}{numel(l_chFRvCh{animal})+1} = ['s' num2str(session)];
    end
    
    % Plot oscillation score vs firing rate
    figure(fH_OSvFR(animal)); hold on
    p = plot(log10(unitFiringRates), oscillationScore.data(:,1), '.', 'MarkerSize',10); hold off
    p_OSvFR{animal} = [p_OSvFR{animal}; p];
    if isempty(l_OSvFR{animal})
      l_OSvFR{animal} = {['s' num2str(session)]};
    else
      l_OSvFR{animal}{numel(l_OSvFR{animal})+1} = ['s' num2str(session)];
    end
  end
end

% Perform oscillation score vs firing rate correlation analysis and fit the line
r = zeros(nAnimals,1);
pval = zeros(nAnimals,1);
coefficients = zeros(nAnimals,2);
for animal = 1:nAnimals

  % Correlation analysis
  [r(animal), pval(animal)] = corrLinearCircular(log10(firingRatesAgregate{animal}), ...
    oscScoreAgregate{animal}, type='Spearman');

  % Line fitting
  figure(fH_OSvFR(animal));
  xLim = xlim;
  xAxisLength = xLim(2) - xLim(1);
  xAxisStep = xAxisLength/10000;
  x = xLim(1):xAxisStep:xLim(2);
  yLim = ylim;
  [~, slope, coefficients(animal,:)] = fitLine(log10(firingRatesAgregate{animal}), ...
    oscScoreAgregate{animal}, type='linear-linear');
  yFit = x.*slope + coefficients(animal,2);
  hold on; plot(x, yFit, 'k--'); hold off;
  xlim(xLim); ylim(yLim);
end

% Label and save the figures
for animal = 1:nAnimals
  animalID = ['PP0' num2str(animal)];
  figFolder = fullfile(generalFigFolder, animalID);
  if ~exist(figFolder,'dir')
    mkdir(figFolder);
  end

  % Unit oscillation score profiles
  figure(fH_OSvCh(animal));
  xlabel('Channel #');
  ylabel('Oscillation score');
  if animal == 1
    legend([p_OSvCh{animal}(1:end-2);p_OSvCh{animal}(end);p_OSvCh{animal}(end-1)], ...
      [l_OSvCh{animal}(1:end-2) l_OSvCh{animal}(end) l_OSvCh{animal}(end-1)], ...
      'Location','northwest');
  else
    legend([p_OSvCh{animal}(1:end-2);p_OSvCh{animal}(end);p_OSvCh{animal}(end-1)], ...
      [l_OSvCh{animal}(1:end-2) l_OSvCh{animal}(end) l_OSvCh{animal}(end-1)], ...
      'Location','northeast');
  end
  legend('boxoff');
  title(['Unit Oscillation Score Profile: ' animalID]);
  set(fH_OSvCh(animal), 'Name',['unit_oscillation_score_profile_' animalID]);
  figFilename = fullfile(figFolder, 'unitOscScoreProfile.fig');
  savefig(fH_OSvCh(animal),figFilename,'compact');

  % Unit firing rate profiles
  figure(fH_uFRvCh(animal));
  xlabel('Channel #');
  ylabel('Firing rate (APs/s)');
  if animal == 1
    legend(p_uFRvCh{animal},l_uFRvCh{animal}, 'Location','northwest');
  else
    legend(p_uFRvCh{animal},l_uFRvCh{animal}, 'Location','northeast');
  end
  legend('boxoff');
  title(['Unit Firing Rate Profile: ' animalID]);
  set(fH_uFRvCh(animal), 'Name',['unit_firing_rate_profile_' animalID]);
  figFilename = fullfile(figFolder, 'unitFiringRateProfile.fig');
  savefig(fH_uFRvCh(animal),figFilename,'compact');

  % Channel firing rate profiles
  figure(fH_chFRvCh(animal));
  xlabel('Channel #');
  ylabel('Firing rate (APs/s)');
  if animal == 1
    legend(p_chFRvCh{animal},l_chFRvCh{animal}, 'Location','northwest');
  else
    legend(p_chFRvCh{animal},l_chFRvCh{animal}, 'Location','northeast');
  end
  legend('boxoff');
  title(['Channel Firing Rate Profile: ' animalID]);
  set(fH_chFRvCh(animal), 'Name',['channel_firing_rate_profile_' animalID]);
  figFilename = fullfile(figFolder, 'channelFiringRateProfile.fig');
  savefig(fH_chFRvCh(animal),figFilename,'compact');

  % Oscillation score vs firing rate
  figure(fH_OSvFR(animal));
  str = ['r=' num2str(r(animal)), ', p=' num2str(pval(animal))];
  yLim = ylim;
  yAxisLength = yLim(2) - yLim(1);
  text(xLim(1)+0.2*xAxisLength, yLim(1)+0.95*yAxisLength, str)
  xlabel('Log firing rate');
  ylabel('Oscillation score');
  if animal == 1
    xlim([-3 2.1])
    legend(p_OSvFR{animal},l_OSvFR{animal}, 'Location','northwest');
  else
    legend(p_OSvFR{animal},l_OSvFR{animal}, 'Location','northeast');
  end
  legend('boxoff');
  title(['Unit Firing Rate vs Oscillation Score: ' animalID]);
  set(fH_OSvFR(animal), 'Name',['unit_firing_rate_vs_oscillation_score_' animalID]);
  figFilename = fullfile(figFolder, 'unitFiringRateVOscillationScore.fig');
  savefig(fH_OSvFR(animal),figFilename,'compact');
end


%% Plot unit count histograms
fH_OSH = zeros(nAnimals,1);
for animal = 1:nAnimals
  animalID = ['PP0' num2str(animal)];
  figFolder = fullfile(generalFigFolder, animalID);
  if ~exist(figFolder,'dir')
    mkdir(figFolder);
  end

  fH_OSH(animal) = figure;
  histogram(maxWaveformChAgregate{animal},max(maxWaveformChAgregate{animal}), ...
    'FaceColor',[.7 .7 .7], 'EdgeColor',[.7 .7 .7]);
  xlabel('Channel #');
  ylabel('Unit count');
  title(['Unit Count across Recording Channels: ' animalID]);
  set(fH_OSH(animal), 'Name',['unit_count_across_recording_channels_' animalID]);
  figFilename = fullfile(figFolder, 'unitCountAcrossRecordingChannels.fig');
  savefig(fH_OSH(animal),figFilename,'compact');
end