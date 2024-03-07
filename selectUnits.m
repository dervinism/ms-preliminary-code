function includeUnits = selectUnits(dataFiles, intervals, firingRateCutoff, ...
  refractoryContaminationTh, oscTh, options)
% includeUnits = selectUnits(dataFiles, intervals, firingRateCutoff, ...
%   refractoryContaminationTh, oscTh, <spikeTimes>)
%
% A helper function of pcaScript: Include/exclude units from analysis.
%
% Args:
%   dataFiles
%   intervals
%   firingRateCutoff
%   refractoryContaminationTh
%   oscTh
%   <spikeTimes>
%
% Returns:
%   includeUnits
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
  firingRateCutoff
  refractoryContaminationTh
  oscTh
  options.spikeTimes = [];
end

% Estimate unit firing rates over 20-minute windows
if isempty(options.spikeTimes)
  spikeTimes = cell(numel(dataFiles),1);
else
  spikeTimes = options.spikeTimes;
end
maxFiringRate = cell(numel(dataFiles),1);
for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session}) && ~isempty(intervals{animal}{session})

      % Load spiking data
      if isempty(options.spikeTimes)
        spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
        if ~exist(spikesFile, 'file')
          continue
        end
        load(spikesFile); %#ok<*LOAD>
        spikeTimes{animal}{session} = spikes.times;
      end

      % count spikes
      nUnits = numel(spikeTimes{animal}{session});
      maxFiringRate{animal}{session} = zeros(nUnits,1);
      for unit = 1:nUnits
        unitSpikeTimes = selectArrayValues(spikeTimes{animal}{session}{unit}', ...
          intervals{animal}{session}([1 end]));
        if isempty(unitSpikeTimes)
          maxFiringRate{animal}{session}(unit) = 0;
        else
          maxFiringRate{animal}{session}(unit) = firingRateWindows( ...
            unitSpikeTimes, stepSize=30, startTime=intervals{animal}{session}(1,1));
        end
      end
    else
      if isempty(options.spikeTimes)
        spikeTimes{animal}{session} = [];
      end
      maxFiringRate{animal}{session} = [];
    end
  end
end

% Select units for individual sessions
includeUnits = cell(numel(dataFiles),1);
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
        nUnits = numel(spikeTimes{animal}{session});
        oscScore = zeros(nUnits,1);
        oscFreq = zeros(nUnits,1);
        for unit = 1:nUnits
          spikeTimes = selectArrayValues( ...
            spikeTimes{animal}{session}{unit}', intervals{animal}{session}([1 end]));
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
    else
      includeUnits{animal}{session} = [];
    end
  end
end