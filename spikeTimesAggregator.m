function [spikeTimesAggregate, populationSpikeTimesAggregate, ...
  populationRateAggregate, populationRateTimeBins, ...
  filtPopulationRateAggregate, unitIDAggregate] = ...
  spikeTimesAggregator(dataFiles, samplingInterval, convPoints, frequencyRange)
% [spikeTimesAggregate, populationSpikeTimesAggregate, ...
%   populationRate, populationRateTimeBins, filtPopulationRateAggregate, ...
%   unitIDAggregate] = spikeTimesAggregator(dataFiles, samplingInterval, ...
%   convPoints, frequencyRange)
%
% A helper function of pcaScript: A wrapper for loading spiking data.
%
% Args:
%   dataFiles
%   samplingInterval
%   convPoints
%
% Returns:
%   spikeTimesAggregate
%   populationSpikeTimesAggregate
%   populationRate
%   populationRateTimeBins
%   filtPopulationRateAggregate
%   unitIDAggregate
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
  convPoints
  frequencyRange
end

% Estimate unit firing rates over 20-minute windows
spikeTimesAggregate = cell(numel(dataFiles),1);
populationSpikeTimesAggregate = cell(numel(dataFiles),1);
populationRateAggregate = cell(numel(dataFiles),1);
populationRateTimeBins = cell(numel(dataFiles),1);
filtPopulationRateAggregate = cell(numel(dataFiles),1);
unitIDAggregate = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session})

      % Load spiking data
      spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
      if ~exist(spikesFile, 'file')
        continue
      end
      load(spikesFile); %#ok<*LOAD>
      spikeTimesAggregate{animal}{session} = spikes.times; %#ok<*NODEF>
      unitIDAggregate{animal}{session} = spikes.cluID;

      % Load the population spike times
      populationRateFile = strrep(dataFiles{animal}{session}, '*', 'populationRate.cellinfo');
      load(populationRateFile);
      populationSpikeTimesAggregate{animal}{session} = populationRate.times{1};

      % Convolve the population spike times
      [populationRateAggregate{animal}{session}, populationRateTimeBins{animal}{session}] = ...
        convolveSpikes(populationSpikeTimesAggregate{animal}{session}, ...
        stepSize=samplingInterval, convolutionPoints=convPoints, ...
        endTime=max(populationSpikeTimesAggregate{animal}{session}));

      % Filter the population spiking rate
      filtPopulationRateAggregate{animal}{session} = bandpassFilterTimeSeries( ...
        populationRateAggregate{animal}{session}, ...
        sampleRate=round(1/samplingInterval), frequencyRange=frequencyRange);
    end
  end
end