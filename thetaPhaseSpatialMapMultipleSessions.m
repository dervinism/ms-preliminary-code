function thetaPhaseSpatialMapMultipleSessions(data, channelOrder, ...
  intervals, figureFolder, figTitle, resamplingInterval, frequencyRange, ...
  parallelise, firingRateCutoff, refractoryContaminationTh, oscTh, ...
  coherenceTh, nPhaseHistogramBins, options)
% thetaPhaseSpatialMapMultipleSessions(data, channelOrder, intervals, ...
%   figureFolder, figTitle, resamplingInterval, frequencyRange, parallelise, ...
%   firingRateCutoff, refractoryContaminationTh, oscTh, coherenceTh, ...
%   nPhaseHistogramBins, <saveResults>)
%
% A helper function of coherenceAnalysisScript.
%
% Args:
%   data
%   channelOrder
%   intervals
%   figureFolder
%   figTitle
%   resamplingInterval
%   frequencyRange
%   parallelise
%   firingRateCutoff
%   refractoryContaminationTh
%   coherenceTh
%   oscTh
%   nPhaseHistogramBins
%   <saveResults>
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
  data
  channelOrder
  intervals
  figureFolder
  figTitle
  resamplingInterval
  frequencyRange
  parallelise
  firingRateCutoff
  refractoryContaminationTh
  oscTh
  coherenceTh
  nPhaseHistogramBins
  options.saveResults = false
end

% Estimate unit firing rates over 20-minute windows
for animal = 1:numel(data)
  for session = 1:numel(data{animal})
    if ~isempty(data{animal}{session})
      sessionData = data{animal}{session};
      if ~isempty(intervals{animal}{session})
        nUnits = numel(sessionData.spikes.ids);
        sessionData.maxFiringRate = zeros(nUnits,1);
        for unit = 1:nUnits
          spikeTimes = selectArrayValues( ...
            sessionData.spikes.times{unit}', intervals{animal}{session}([1 end]));
          if isempty(spikeTimes)
            sessionData.maxFiringRate(unit) = 0;
          else
            sessionData.maxFiringRate(unit) = firingRateWindows( ...
              spikeTimes, stepSize=30, startTime=intervals{animal}{session}(1,1));
          end
        end
        data{animal}{session} = sessionData;
      end
    end
  end
end

% Estimate coherence and perform correlation analyses for individual sessions
for animal = 1:numel(data)
  for session = 1:numel(data{animal})
    if ~isempty(data{animal}{session})
      sessionData = data{animal}{session};
      if ~isempty(intervals{animal}{session})

        includeUnits = logical( ...
          sessionData.maxFiringRate >= firingRateCutoff/3600 & ...
          sessionData.contaminationPercent.data <= refractoryContaminationTh/100);
        maxThetaFrequency = [sessionData.thetaPower.timestamps'; ...
          sessionData.thetaPower.data(:,1)'];
        sessionFigFolder = fullfile(figureFolder, sessionData.sessionInfo.FileName);

        % Oscillation score calculations and histograms
        if oscTh > 0
          if ~isempty(data{animal}{session})
            if ~isempty(intervals{animal}{session})
              nUnits = numel(sessionData.spikes.ids);
              sessionData.oscScore = zeros(nUnits,1);
              sessionData.oscFreq = zeros(nUnits,1);
              for unit = 1:nUnits
                spikeTimes = selectArrayValues( ...
                  sessionData.spikes.times{unit}', intervals{animal}{session}([1 end]));
                if isempty(spikeTimes)
                  sessionData.oscScore(unit) = 0;
                else
                  sr = 500;
                  [sessionData.oscScore(unit), ~, sessionData.oscFreq(unit)] = OScoreSpikes( ...
                    {round(spikeTimes*sr)}, round(spikeTimes(end)*sr), 5, 11, sr);
                end
              end
              data{animal}{session} = sessionData;

              % Plot histograms
              fH = figure;
              histogram(sessionData.oscScore,round(max(sessionData.oscScore))+1, ...
                'FaceColor',[.7 .7 .7], 'EdgeColor',[.7 .7 .7]);
              hold on
              yLim = ylim;
              p = plot([oscTh oscTh],yLim, ':r');
              ylim(yLim);
              hold off
              xlabel('Oscillation score');
              ylabel('Unit count');
              legend(p,{'cutoff'}); legend('boxoff');
              title('Oscillation Score Distribution');

              % Save histograms
              if ~exist(sessionFigFolder,'dir')
                mkdir(sessionFigFolder);
              end
              filename = [sessionFigFolder filesep 'osc_score_hist.fig']; %#ok<*AGROW>
              savefig(fH,filename,'compact');
              title('');
              saveas(fH,filename(1:end-4),'png');
              close(fH);
              
              if oscTh > 0
                includeUnits = logical(includeUnits & sessionData.oscScore > oscTh);
              elseif oscTh < 0
                includeUnits = logical(includeUnits & sessionData.oscScore <= abs(oscTh));
              end
            end
          end
        end

        % Phase scattergrams
        [sessionData.fullCoherence, sessionData.thetaPhaseTopography, ...
          sessionData.fullInterpCoherence, sessionData.interpThetaPhaseTopography] = ...
          thetaPhaseSpatialMapInterval(intervals{animal}{session}, ...
          sessionData.spikes.times', ...
          sessionData.populationRate.times{1}', ...
          channelOrder{animal}(session,sessionData.spikes.maxWaveformCh+1)', ...
          figTitle, resamplingInterval, frequencyRange, ...
          sessionData.instantThetaFrequency.data, ...
          sessionData.instantThetaFrequency.timestamps, ...
          coherenceRange='full', parallelise=parallelise, ...
          include=includeUnits, maxThetaFrequency=maxThetaFrequency, ...
          figPath=sessionFigFolder);

        % Phase histograms
        % Regular coherence
        [sessionData.fullCoherence.phaseHistograms, ...
          sessionData.fullCoherence.phaseHistogramBins, ...
          sessionData.fullCoherence.significantPhaseCounts, ...
          sessionData.fullCoherence.meanPhase] = ...
          phaseHistrogram(sessionData.fullCoherence.phase(includeUnits,:), ...
          centre=0, nBins=nPhaseHistogramBins);
        % Interpolated coherence
        [sessionData.fullInterpCoherence.phaseHistograms, ...
          sessionData.fullInterpCoherence.phaseHistogramBins, ...
          sessionData.fullInterpCoherence.significantPhaseCounts, ...
          sessionData.fullInterpCoherence.meanPhase] = ...
          phaseHistrogram(sessionData.fullInterpCoherence.phase(includeUnits,:), ...
          centre=0, nBins=nPhaseHistogramBins);

        % Phase histogram plot: Most coherent theta frequency
        fInd = sessionData.thetaPhaseTopography.fIndMostCoh;
        figText = {[num2str(sessionData.fullCoherence.significantPhaseCounts(fInd)) ...
          '/' num2str(sum(includeUnits))]};
        histTitle = ['hist_' figTitle ' (most coherent freq)'];
        phaseHistogramPlot(sessionData.fullCoherence.phaseHistograms(:,fInd)', ...
          binLocs=sessionData.fullCoherence.phaseHistogramBins, ...
          dataMean=sessionData.fullCoherence.meanPhase(fInd), ...
          figText=figText, figTitle=histTitle, figPath=sessionFigFolder);

        % Phase histogram plot: Most significant theta frequency
        fInd = sessionData.thetaPhaseTopography.fIndMostSignificant;
        figText = {[num2str(sessionData.fullCoherence.significantPhaseCounts(fInd)) ...
          '/' num2str(sum(includeUnits))]};
        histTitle = ['hist_' figTitle ' (most significant freq)'];
        phaseHistogramPlot(sessionData.fullCoherence.phaseHistograms(:,fInd)', ...
          binLocs=sessionData.fullCoherence.phaseHistogramBins, ...
          dataMean=sessionData.fullCoherence.meanPhase(fInd), ...
          figText=figText, figTitle=histTitle, figPath=sessionFigFolder);

        % Phase histogram plot: Most powerful theta frequency
        fInd = sessionData.interpThetaPhaseTopography.fIndMostPower;
        figText = {[num2str(sessionData.fullInterpCoherence.significantPhaseCounts(fInd)) ...
          '/' num2str(sum(includeUnits))]};
        histTitle = ['hist_' figTitle ' (most powerful freq)'];
        phaseHistogramPlot(sessionData.fullInterpCoherence.phaseHistograms(:,fInd)', ...
          binLocs=sessionData.fullInterpCoherence.phaseHistogramBins, ...
          dataMean=sessionData.fullInterpCoherence.meanPhase(fInd), ...
          figText=figText, figTitle=histTitle, figPath=sessionFigFolder);

        % Phase histogram plot: Mean instantaneous theta frequency
        fInd = sessionData.interpThetaPhaseTopography.fIndMeanInstant;
        figText = {[num2str(sessionData.fullInterpCoherence.significantPhaseCounts(fInd)) ...
          '/' num2str(sum(includeUnits))]};
        histTitle = ['hist_' figTitle ' (mean instant freq)'];
        phaseHistogramPlot(sessionData.fullInterpCoherence.phaseHistograms(:,fInd)', ...
          binLocs=sessionData.fullInterpCoherence.phaseHistogramBins, ...
          dataMean=sessionData.fullInterpCoherence.meanPhase(fInd), ...
          figText=figText, figTitle=histTitle, figPath=sessionFigFolder);

        % Save unit phase
        if options.saveResults
          spikes = sessionData.spikes;
          spikes.phase_coherence = sessionData.interpThetaPhaseTopography.phase(:,sessionData.interpThetaPhaseTopography.fIndMeanInstant)';
          spikesFilename = [fullfile(sessionData.sessionInfo.session.path, ...
            sessionData.sessionInfo.session.name) '.spikes.cellinfo.mat'];
          save(spikesFilename, 'spikes', '-v7.3');
        end

        data{animal}{session} = sessionData;
      end
    end
  end
end

% Mark valid unit coherence values for further analyses
for animal = 1:numel(data)
  for session = 1:numel(data{animal})
    if ~isempty(data{animal}{session})
      sessionData = data{animal}{session};
      if isfield(sessionData, 'fullCoherence')
        % Regular coherence
        maxCoherence = max(sessionData.fullCoherence.rateAdjustedCoherence,[],2);
        sessionData.fullCoherence.validCoherenceValues = ...
          false(size(sessionData.fullCoherence.rateAdjustedCoherence));
        for f = 1:size(sessionData.fullCoherence.frequency,2)
          sessionData.fullCoherence.validCoherenceValues(:,f) = logical( ...
            sessionData.maxFiringRate >= firingRateCutoff/3600 & ...
            sessionData.contaminationPercent.data <= refractoryContaminationTh/100 & ...
            sessionData.fullCoherence.rateAdjustedCoherence(:,f) >= coherenceTh & ...
            sessionData.fullCoherence.rateAdjustedCoherence(:,f) == maxCoherence);
        end
        sessionData.fullCoherence.validUnits = logical( ...
          sum(sessionData.fullCoherence.validCoherenceValues,2));
        % Interpolated coherence
        maxCoherence = max(sessionData.fullInterpCoherence.rateAdjustedCoherence,[],2);
        sessionData.fullInterpCoherence.validCoherenceValues = ...
          false(size(sessionData.fullInterpCoherence.rateAdjustedCoherence));
        for f = 1:size(sessionData.fullInterpCoherence.frequency,2)
          sessionData.fullInterpCoherence.validCoherenceValues(:,f) = logical( ...
            sessionData.maxFiringRate >= firingRateCutoff/3600 & ...
            sessionData.contaminationPercent.data <= refractoryContaminationTh/100 & ...
            sessionData.fullInterpCoherence.rateAdjustedCoherence(:,f) >= coherenceTh & ...
            sessionData.fullInterpCoherence.rateAdjustedCoherence(:,f) == maxCoherence);
        end
        sessionData.fullInterpCoherence.validUnits = logical( ...
          sum(sessionData.fullInterpCoherence.validCoherenceValues,2));
        data{animal}{session} = sessionData;
      end
    end
  end
end

% Concatenate phase and regression slope values
%phase = [];
%frequency = [];
%maxChan = [];
%spikeTimes = [];
mostCoherentSlopes = cell(numel(data),1);
mostCoherentPvalues = cell(numel(data),1);
mostCoherentExplainedVar = cell(numel(data),1);
mostPowerfulSlopes = cell(numel(data),1);
mostPowerfulPvalues = cell(numel(data),1);
mostPowerfulExplainedVar = cell(numel(data),1);
mostSignificantSlopes = cell(numel(data),1);
mostSignificantPvalues = cell(numel(data),1);
mostSignificantExplainedVar = cell(numel(data),1);
meanInstantSlopes = cell(numel(data),1);
meanInstantPvalues = cell(numel(data),1);
meanInstantExplainedVar = cell(numel(data),1);
for animal = 1:numel(data)
  for session = 1:numel(data{animal})
    if ~isempty(data{animal}{session})
      sessionData = data{animal}{session};
      if isfield(sessionData, 'fullCoherence')
        % Phase
        %validValuesIdx = sessionData.fullCoherence.validCoherenceValues;
        %validUnitsIdx = sessionData.fullCoherence.validUnits;
        %phase = [phase; sessionData.fullCoherence.phase(validValuesIdx)]; %#ok<*AGROW> 
        %frequency = [frequency; sessionData.fullCoherence.frequency(validValuesIdx)];
        %maxChan = [maxChan; channelOrder{animal}(session,sessionData.spikes.maxWaveformCh(validUnitsIdx)+1)'];
        %spikeTimes = [spikeTimes; sessionData.spikes.times(validUnitsIdx)'];
        % Slopes
        % Most coherent theta frequency
        fInd = sessionData.thetaPhaseTopography.fIndMostCoh;
        mostCoherentSlopes{animal} = [mostCoherentSlopes{animal}; ...
          sessionData.thetaPhaseTopography.coefficients(fInd,1)];
        mostCoherentPvalues{animal} = [mostCoherentPvalues{animal}; ...
          sessionData.thetaPhaseTopography.pval(fInd)];
        mostCoherentExplainedVar{animal} = [mostCoherentExplainedVar{animal}; ...
          sessionData.thetaPhaseTopography.r(fInd)^2];
        % Most significant theta frequency
        fInd = sessionData.thetaPhaseTopography.fIndMostSignificant;
        mostSignificantSlopes{animal} = [mostSignificantSlopes{animal}; ...
          sessionData.thetaPhaseTopography.coefficients(fInd,1)];
        mostSignificantPvalues{animal} = [mostSignificantPvalues{animal}; ...
          sessionData.thetaPhaseTopography.pval(fInd)];
        mostSignificantExplainedVar{animal} = [mostSignificantExplainedVar{animal}; ...
          sessionData.thetaPhaseTopography.r(fInd)^2];
        % Most powerful theta frequency
        mostPowerfulSlopes{animal} = [mostPowerfulSlopes{animal}; ...
          sessionData.interpThetaPhaseTopography.coefficients(1,1)];
        mostPowerfulPvalues{animal} = [mostPowerfulPvalues{animal}; ...
          sessionData.interpThetaPhaseTopography.pval(1)];
        mostPowerfulExplainedVar{animal} = [mostPowerfulExplainedVar{animal}; ...
          sessionData.interpThetaPhaseTopography.r(1)^2];
        % Mean Instantaneous theta frequency
        meanInstantSlopes{animal} = [meanInstantSlopes{animal}; ...
          sessionData.interpThetaPhaseTopography.coefficients(2,1)];
        meanInstantPvalues{animal} = [meanInstantPvalues{animal}; ...
          sessionData.interpThetaPhaseTopography.pval(2)];
        meanInstantExplainedVar{animal} = [meanInstantExplainedVar{animal}; ...
          sessionData.interpThetaPhaseTopography.r(2)^2];
      end
    end
  end
end

% Perform phase correlation analysis and display the results
% thetaPhaseSpatialMap(phase, maxChan, mean(frequency), spikeTimes, ...
%   figTitle, include=true(size(phase)), figPath=figureFolder);

% Display correlation slopes for individual animals
colours = {'b','r'};
for animal = 1:numel(data)
  % Most coherent theta frequency
  slopeFigTitle = [figTitle '_PP0' num2str(animal) '_most_coherent_slopes'];
  slopePlot(mostCoherentSlopes{animal}, mostCoherentPvalues{animal}, ...
    colours{animal}, slopeFigTitle, figureFolder);
  significanceFigTitle = [figTitle '_PP0' num2str(animal) '_most_coherent_p-values'];
  significancePlot(mostCoherentPvalues{animal}, ...
    colours{animal}, significanceFigTitle, figureFolder);
  explainedVarFigTitle = [figTitle '_PP0' num2str(animal) '_most_coherent_R^2'];
  explainedVarPlot(mostCoherentExplainedVar{animal}, mostCoherentPvalues{animal}, ...
    colours{animal}, explainedVarFigTitle, figureFolder);
  % Most significant theta frequency
  slopeFigTitle = [figTitle '_PP0' num2str(animal) '_most_significant_slopes'];
  slopePlot(mostSignificantSlopes{animal}, mostSignificantPvalues{animal}, ...
    colours{animal}, slopeFigTitle, figureFolder);
  significanceFigTitle = [figTitle '_PP0' num2str(animal) '_most_significant_p-values'];
  significancePlot(mostSignificantPvalues{animal}, ...
    colours{animal}, significanceFigTitle, figureFolder);
  explainedVarFigTitle = [figTitle '_PP0' num2str(animal) '_most_significant_R^2'];
  explainedVarPlot(mostSignificantExplainedVar{animal}, mostSignificantPvalues{animal}, ...
    colours{animal}, explainedVarFigTitle, figureFolder);
  % Most powerful theta frequency
  slopeFigTitle = [figTitle '_PP0' num2str(animal) '_most_powerful_slopes'];
  slopePlot(mostPowerfulSlopes{animal}, mostPowerfulPvalues{animal}, ...
    colours{animal}, slopeFigTitle, figureFolder);
  significanceFigTitle = [figTitle '_PP0' num2str(animal) '_most_powerful_p-values'];
  significancePlot(mostPowerfulPvalues{animal}, ...
    colours{animal}, significanceFigTitle, figureFolder);
  explainedVarFigTitle = [figTitle '_PP0' num2str(animal) '_most_powerful_R^2'];
  explainedVarPlot(mostPowerfulExplainedVar{animal}, mostPowerfulPvalues{animal}, ...
    colours{animal}, explainedVarFigTitle, figureFolder);
  % Mean Instantaneous theta frequency
  slopeFigTitle = [figTitle '_PP0' num2str(animal) '_mean_instant_slopes'];
  slopePlot(meanInstantSlopes{animal}, meanInstantPvalues{animal}, ...
    colours{animal}, slopeFigTitle, figureFolder);
  significanceFigTitle = [figTitle '_PP0' num2str(animal) '_mean_instant_p-values'];
  significancePlot(meanInstantPvalues{animal}, ...
    colours{animal}, significanceFigTitle, figureFolder);
  explainedVarFigTitle = [figTitle '_PP0' num2str(animal) '_mean_instant_R^2'];
  explainedVarPlot(meanInstantExplainedVar{animal}, meanInstantPvalues{animal}, ...
    colours{animal}, explainedVarFigTitle, figureFolder);
end

close all
end



%% Local functions
function fH = slopePlot(slopes, pval, colour, figTitle, figPath)

% Plot slopes
nonSignificantIdx = pval >= 0.05;
session = 1:numel(slopes);
fH = figure; plot(session(nonSignificantIdx), slopes(nonSignificantIdx), ...
  'o', 'Color',colour); hold on
plot(session(~nonSignificantIdx), slopes(~nonSignificantIdx), ...
  '.', 'Color',colour, 'MarkerSize',20); hold off

% Create the legend
if sum(nonSignificantIdx)
  legend('p\geq.05','p<.05', 'Location','SouthEast')
else
  legend('p<.05', 'Location','SouthEast')
end

% Get axes' dimensions
xLim = xlim;
xAxisLength = xLim(2) - xLim(1);
yLim = ylim;
yAxisLength = yLim(2) - yLim(1);

% Adjust axes' limits
xlim([xLim(1)-0.05*xAxisLength xLim(2)+0.05*xAxisLength])
ylim([yLim(1)-0.05*yAxisLength yLim(2)+0.05*yAxisLength])

% Axes labels
xlabel('session')
ylabel('Regression slope (channels/rad)')

% Figure title
tH = title(figTitle, 'Interpreter','none');

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
filename = [figPath filesep filename '.fig'];
savefig(fH,filename,'compact');
title('');
saveas(fH,filename(1:end-4),'png');
close(fH);
end


function fH = significancePlot(pval, colour, figTitle, figPath)

% Plot p-values
nonSignificantIdx = pval >= 0.05;
session = 1:numel(pval);
fH = figure; plot(session(nonSignificantIdx), pval(nonSignificantIdx), ...
  'o', 'Color',colour); hold on
plot(session(~nonSignificantIdx), pval(~nonSignificantIdx), ...
  '.', 'Color',colour, 'MarkerSize',20); hold off

% Create the legend
if sum(nonSignificantIdx)
  legend('p\geq.05','p<.05', 'Location','NorthEast')
else
  legend('p<.05', 'Location','SouthEast')
end

% Get axes' dimensions
xLim = xlim;
xAxisLength = xLim(2) - xLim(1);
yLim = ylim;
yAxisLength = yLim(2) - yLim(1);

% Adjust axes' limits
xlim([xLim(1)-0.05*xAxisLength xLim(2)+0.05*xAxisLength])
ylim([yLim(1)-0.05*yAxisLength yLim(2)+0.05*yAxisLength])
set(gca, 'YScale', 'log')

% Axes labels
xlabel('session')
ylabel('p-values')

% Figure title
tH = title(figTitle, 'Interpreter','none');

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
filename = [figPath filesep filename '.fig'];
savefig(fH,filename,'compact');
title('');
saveas(fH,filename(1:end-4),'png');
close(fH);
end


function fH = explainedVarPlot(RSquared, pval, colour, figTitle, figPath)

% Plot R2-values
nonSignificantIdx = pval >= 0.05;
session = 1:numel(RSquared);
fH = figure; plot(session(nonSignificantIdx), RSquared(nonSignificantIdx), ...
  'o', 'Color',colour); hold on
plot(session(~nonSignificantIdx), RSquared(~nonSignificantIdx), ...
  '.', 'Color',colour, 'MarkerSize',20); hold off

% Create the legend
if sum(nonSignificantIdx)
  legend('p\geq.05','p<.05', 'Location','SouthEast')
else
  legend('p<.05', 'Location','SouthEast')
end

% Get axes' dimensions
xLim = xlim;
xAxisLength = xLim(2) - xLim(1);
yLim = ylim;
yAxisLength = yLim(2) - yLim(1);

% Adjust axes' limits
xlim([xLim(1)-0.05*xAxisLength xLim(2)+0.05*xAxisLength])
ylim([yLim(1)-0.05*yAxisLength yLim(2)+0.05*yAxisLength])
set(gca, 'YScale', 'log')

% Axes labels
xlabel('session')
ylabel('R^2')

% Figure title
tH = title(figTitle, 'Interpreter','none');

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
filename = [figPath filesep filename '.fig'];
savefig(fH,filename,'compact');
title('');
saveas(fH,filename(1:end-4),'png');
close(fH);
end