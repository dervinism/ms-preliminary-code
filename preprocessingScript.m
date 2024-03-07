% This script performs data preprocessing for the ms-preliminary dataset


%% Initialise parameters
% Data setup parameters
params

% Analysis parameters
runRefractContamAnalysis = false;
saveLFPNoisePeriods = false;
saveOscillationScore = true;


%% Preprocessing steps
% Refractory period contamination analysis
if runRefractContamAnalysis
  saveRefractoryContamination(dataFiles); %#ok<*UNRCH>
end

% Mark noise periods
if saveLFPNoisePeriods
  saveLFPNoise(lfpNoise, dataFiles);
end

% Save oscillation score
if saveOscillationScore
  saveOScore(dataFiles)
end




%% Local functions
function refractoryContaminationData = saveRefractoryContamination(dataFiles)
% refractoryContaminationData = saveRefractoryContamination(dataFiles)
%
% Function estimates refractory period contamination of individual unit
% autocorellograms for multiple recording sessions located in the
% ms-preliminary repository, saves it, and returns it stored in a single
% data structure variable. It's a helper function of preprocessingScript.
%
% Args:
%   dataFiles
%
% Returns:
%   refractoryContaminationData
%
% Comments:
%   The function is not intended for wide use. If you intend to use it, get
%   in touch with the author.
%
% Authors:
%   Martynas Dervinis (martynas.dervinis@gmail.com).

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    if ~isempty(dataFiles{animal}{session})

      % Load spiking data
      spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
      if ~exist(spikesFile, 'file')
        continue
      end
      load(spikesFile); %#ok<*LOAD>
      nUnits = numel(spikes.ids);

      % Calculate refractory period contamination percentage
      contaminationPercentage = zeros(nUnits,1);
      for unit = 1:nUnits
        contaminationPercentage(unit) = ...
          refractoryContamination(spikes.times{unit}');
      end

      % Set output
      refractoryContaminationData{animal}{session} = contaminationPercentage; %#ok<*AGROW>

      % Save calculations
      contaminationPercent.data = contaminationPercentage;
      contaminationPercent.definition = ['Stores contamination percentage of ' ...
        'the unit''s spiking autocorrelagram (ACG) 1 ms refractory period around zero lag ' ...
        'relative to the ACG shoulder (1000-1500 ms away from the zero lag). ' ...
        'The measure is calculated for each unit.'];
      contaminationPercent.processingInfo.params = [];
      contaminationPercent.processingInfo.function = 'petersen-lab-matlab/spikes/refractoryContamination';
      contaminationPercent.processingInfo.date = datetime;
      contaminationPercent.processingInfo.username = getenv('username');
      contaminationPercent.processingInfo.hostname = getenv('computername');
      qualityFile = strrep(spikesFile, 'spikes', 'contaminationPercent');
      save(qualityFile, 'contaminationPercent', '-v7.3');
    end
  end
end
end


function saveLFPNoise(lfpNoisePeriods, dataFiles)
% saveLFPNoise(lfpNoisePeriods, dataFiles)
%
% Function saves LFP noise periods in CellExplorer format for the
% ms-preliminary repoistory data. It's a helper function of
% preprocessingScript.
%
% Args:
%   lfpNoisePeriods
%   dataFiles
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

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    for fn = 1:numel(lfpNoisePeriods)
      if contains(dataFiles{animal}{session}, lfpNoisePeriods(fn).basename)
        saveLFPNoisePeriods(lfpNoisePeriods(fn).periods, ...
          outputFilename=fileparts(dataFiles{animal}{session}));
      end
    end
  end
end
end


function saveOScore(dataFiles)
% saveOScore(dataFiles)
%
% Function saves unit oscillation scores in CellExplorer format for the
% ms-preliminary dataset. It's a helper function of preprocessingScript.
%
% Args:
%   dataFiles
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

for animal = 1:numel(dataFiles)
  for session = 1:numel(dataFiles{animal})
    % Load spiking data
    spikesFile = strrep(dataFiles{animal}{session}, '*', 'spikes.cellinfo');
    if ~exist(spikesFile, 'file')
      continue
    end
    saveOscScore(spikesFile, shuffle=true);
  end
end
end