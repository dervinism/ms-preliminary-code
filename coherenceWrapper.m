function [fullCoherence, fullInterpCoherence, includeUnits, spikeTimesAggregate, ...
  populationRateAggregate] = coherenceWrapper(dataFiles, intervals, ...
  resamplingInterval, frequencyRange, parallelise, firingRateCutoff, ...
  refractoryContaminationTh, oscTh, coherenceTh)
% [fullCoherence, fullInterpCoherence, includeUnits, spikeTimesAggregate, ...
%   populationRateAggregate] = coherenceWrapper(dataFiles, intervals, ...
%   resamplingInterval, frequencyRange, parallelise, firingRateCutoff, ...
%   refractoryContaminationTh, oscTh, coherenceTh)
%
% A helper function of pcaScript: A wrapper for coherence analysis.
%
% Args:
%   dataFiles
%   intervals
%   resamplingInterval
%   frequencyRange
%   parallelise
%   firingRateCutoff
%   refractoryContaminationTh
%   oscTh
%   coherenceTh
%
% Returns:
%   fullCoherence
%   fullInterpCoherence
%   includeUnits
%   spikeTimesAggregate
%   populationRateAggregate
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
  resamplingInterval
  frequencyRange
  parallelise
  firingRateCutoff
  refractoryContaminationTh
  oscTh
  coherenceTh
end

% Estimate unit firing rates over 20-minute windows
spikeTimesAggregate = cell(numel(dataFiles),1);
maxFiringRate = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session}) && ~isempty(intervals{animal}{session})

      % Load spiking data
      spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
      if ~exist(spikesFile, 'file')
        continue
      end
      load(spikesFile); %#ok<*LOAD>
      spikeTimesAggregate{animal}{session} = spikes.times; %#ok<*NODEF>

      % count spikes
      nUnits = numel(spikes.ids);
      maxFiringRate{animal}{session} = zeros(nUnits,1);
      for unit = 1:nUnits
        spikeTimes = selectArrayValues( ...
          spikes.times{unit}', intervals{animal}{session}([1 end]));
        if isempty(spikeTimes)
          maxFiringRate{animal}{session}(unit) = 0;
        else
          maxFiringRate{animal}{session}(unit) = firingRateWindows( ...
            spikeTimes, stepSize=30, startTime=intervals{animal}{session}(1,1));
        end
      end
    else
      spikeTimesAggregate{animal}{session} = [];
      maxFiringRate{animal}{session} = [];
    end
  end
end

% Estimate coherence for individual sessions
populationRateAggregate = cell(numel(dataFiles),1);
includeUnits = cell(numel(dataFiles),1);
fullCoherence = cell(numel(dataFiles),1);
fullInterpCoherence = cell(numel(dataFiles),1);
contaminationPercentData = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session}) && ~isempty(intervals{animal}{session})

      % Load unit quality data
      qualityFile = strrep(dataFiles{animal}{session}, '*', 'contaminationPercent.cellinfo');
      load(qualityFile);
      contaminationPercentData{animal}{session} = contaminationPercent.data;

      includeUnits{animal}{session} = logical( ...
        maxFiringRate{animal}{session} >= firingRateCutoff/3600 & ...
        contaminationPercent.data <= refractoryContaminationTh/100);

      % Oscillation score calculations
      if oscTh > 0
        if ~isempty(dataFiles{animal}{session})
          if ~isempty(intervals{animal}{session})
            nUnits = numel(spikeTimesAggregate{animal}{session});
            oscScore = zeros(nUnits,1);
            oscFreq = zeros(nUnits,1);
            for unit = 1:nUnits
              spikeTimes = selectArrayValues( ...
                spikeTimesAggregate{animal}{session}{unit}', intervals{animal}{session}([1 end]));
              if isempty(spikeTimes)
                oscScore(unit) = 0;
              else
                sr = 500;
                [oscScore(unit), ~, oscFreq(unit)] = OScoreSpikes( ...
                  {round(spikeTimes*sr)}, round(spikeTimes(end)*sr), 5, 11, sr);
              end
            end

            if oscTh > 0
              includeUnits{animal}{session} = logical(includeUnits{animal}{session} ...
                & oscScore > oscTh);
            elseif oscTh < 0
              includeUnits{animal}{session} = logical(includeUnits{animal}{session} ...
                & oscScore <= abs(oscTh));
            end
          end
        end
      end

      % Load the population firing rate
      populationRateFile = strrep(dataFiles{animal}{session}, '*', 'populationRate.cellinfo');
      load(populationRateFile);
      populationRateAggregate{animal}{session} = selectArrayValues( ...
        populationRate.times{1}, intervals{animal}{session});

      % Calculate coherence
      unitSpikeTimes = selectCellValues(spikeTimesAggregate{animal}{session}, intervals{animal}{session});
      freqGrid = frequencyRange(1):0.01:frequencyRange(2);
      [fullCoherence{animal}{session}, ~, ~, fullInterpCoherence{animal}{session}] = ...
        coherence(unitSpikeTimes, populationRateAggregate{animal}{session}, ...
        intervals=intervals{animal}{session}, stepsize=resamplingInterval, ...
        freqRange=frequencyRange, freqGrid=freqGrid, parallelise=parallelise);
    else
      includeUnits{animal}{session} = [];
      fullCoherence{animal}{session} = [];
      fullInterpCoherence{animal}{session} = [];
      contaminationPercentData{animal}{session} = [];
    end
  end
end

% Mark valid unit coherence values for further analyses
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session}) && ~isempty(fullCoherence{animal}{session}) ...
        && ~isempty(fullCoherence{animal}{session}.rateAdjustedCoherence)
      % Regular coherence
      maxCoherence = max(fullCoherence{animal}{session}.rateAdjustedCoherence,[],2);
      fullCoherence{animal}{session}.validCoherenceValues = ...
        false(size(fullCoherence{animal}{session}.rateAdjustedCoherence));
      for f = 1:size(fullCoherence{animal}{session}.frequency,2)
        fullCoherence{animal}{session}.validCoherenceValues(:,f) = logical( ...
          includeUnits{animal}{session} & ...
          fullCoherence{animal}{session}.rateAdjustedCoherence(:,f) >= coherenceTh & ...
          fullCoherence{animal}{session}.rateAdjustedCoherence(:,f) == maxCoherence);
      end
      fullCoherence{animal}{session}.validUnits = logical( ...
        sum(fullCoherence{animal}{session}.validCoherenceValues,2));
      includeUnits{animal}{session} = logical( ...
        sum(fullCoherence{animal}{session}.validCoherenceValues,2));
      % Interpolated coherence
      maxCoherence = max(fullInterpCoherence{animal}{session}.rateAdjustedCoherence,[],2);
      fullInterpCoherence{animal}{session}.validCoherenceValues = ...
        false(size(fullInterpCoherence{animal}{session}.rateAdjustedCoherence));
      for f = 1:size(fullInterpCoherence{animal}{session}.frequency,2)
        fullInterpCoherence{animal}{session}.validCoherenceValues(:,f) = logical( ...
          includeUnits{animal}{session} & ...
          fullInterpCoherence{animal}{session}.rateAdjustedCoherence(:,f) >= coherenceTh & ...
          fullInterpCoherence{animal}{session}.rateAdjustedCoherence(:,f) == maxCoherence);
      end
      fullInterpCoherence{animal}{session}.validUnits = logical( ...
        sum(fullInterpCoherence{animal}{session}.validCoherenceValues,2));
    end
  end
end