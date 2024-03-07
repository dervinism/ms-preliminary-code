%% I/O parameters
% Data path info
dataFolder = 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\03_data';
animals = {'PP01'}; %{'PP01';'PP02'};
figureFolder = 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\05_figures\phase_topography\individual_figures';
generalFigFolder = fullfile(figureFolder, 'general');
cohFigFolder = fullfile(figureFolder, 'coherence_analyses\all_units');
travellingWavesFigFolder = fullfile(figureFolder, 'travelling_waves');
unitPhaseFigFolder = fullfile(figureFolder,'unit_phases_within_sessions\wideband');

figureFolder = 'Z:\SUN-IN-Petersen-lab\EphysData\MartynasDervinis\ms-preliminary\05_figures\dimensionality_reduction\individual_figures\pca\latestOutput';
pcaFigFolder = figureFolder;

% Lists of recording sessions and raw/preprocessed data files
recSessions = cell(numel(animals),1);
dataFiles = cell(numel(animals),1);
lfpFiles = cell(numel(animals),1);
for animal = 1:numel(animals)
  animalDataFolder = fullfile(dataFolder, animals{animal});
  animalDataFolderContents = {dir(animalDataFolder).name}';
  animalDataFolderContents = animalDataFolderContents(3:end);
  recSessions{animal} = animalDataFolderContents;
  dataFiles{animal} = cellfun(@(x) fullfile(dataFolder, animals{animal}, ...
    x, [x '.*.mat']), animalDataFolderContents, 'UniformOutput', false);
  lfpFiles{animal} = cellfun(@(x) fullfile(dataFolder, animals{animal}, ...
    x, [x '.lfp']), animalDataFolderContents, 'UniformOutput', false);
end


%% Probe recording channel order
regularLayout = 1:384;
disjointLayout2_1 = [93:384 1:92];
disjointLayout2_2 = [117:384 1:116];
msChannels = {1:384; 1:384}; %{201:330; 1:384};
probeLayouts = cell(numel(animals),1);
channelOrder = cell(numel(animals),1);
channelsOI = cell(numel(animals),1);
for animal = 1:numel(animals)
  for session = 1:numel(dataFiles{animal})
    if strcmpi(animals{animal},'PP01')
      dataFilename = dataFiles{animal}{session};
      if contains(dataFilename, {'PP01_2020-06-23_14-00-45'})
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      elseif contains(dataFilename, {'PP01_2020-06-25_15-54-22'})
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      elseif contains(dataFilename, {'PP01_2020-06-29_13-15-57'})
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      elseif contains(dataFilename, {'PP01_2020-06-30_12-27-33'})
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      elseif contains(dataFilename, {'PP01_2020-07-01_13-07-16'})
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      elseif contains(dataFilename, {'PP01_2020-07-02_12-41-03'})
        probeLayouts{animal} = [probeLayouts{animal}; [83:384 1:82]];
      elseif contains(dataFilename, {'PP01_2020-07-08_17-24-19'})
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      elseif contains(dataFilename, {'PP01_2020-07-09_13-22-13'})
        probeLayouts{animal} = [probeLayouts{animal}; [83:384 1:82]];
      elseif contains(dataFilename, {'PP01_2020-07-14_15-58-44'})
        probeLayouts{animal} = [probeLayouts{animal}; [101:384 1:100]];
      elseif contains(dataFilename, {'PP01_2020-07-15_15-17-43'})
        probeLayouts{animal} = [probeLayouts{animal}; [101:384 1:100]];
      elseif contains(dataFilename, {'PP01_2020-07-21_13-00-21'})
        probeLayouts{animal} = [probeLayouts{animal}; [101:384 1:100]];
      elseif contains(dataFilename, {'PP01_2020-07-23_16-22-54'})
        probeLayouts{animal} = [probeLayouts{animal}; [101:384 1:100]];
      else
        probeLayouts{animal} = [probeLayouts{animal}; regularLayout];
      end
    elseif strcmpi(animals{animal},'PP02')
      if contains(dataFilename, '2020-07-07')
        probeLayouts{animal} = [probeLayouts{animal}; disjointLayout2_1];
      else
        probeLayouts{animal} = [probeLayouts{animal}; disjointLayout2_2];
      end
    end
    [~, layoutSorting] = sort(probeLayouts{animal}(session,:));
    channelOrder{animal} = [channelOrder{animal}; layoutSorting];
    channelsOI{animal} = [channelsOI{animal}; msChannels{animal}];
  end
end


%% LFP saturations
lfpNoise(1).basename = 'PP01_2020-06-23_14-00-45';
lfpNoise(1).periods =  [   98.858    98.867
                          567.500   587.436
                         1157.220  1157.230
                         1257.190  1262.320
                         1654.490  1655.520
                         1742.370  1763.090
                         2255.080  2255.190];

lfpNoise(2).basename = 'PP01_2020-06-25_15-54-22';
lfpNoise(2).periods =  [  162.824   208.937
                          862.011   862.014
                         1599.800  1607.100
                         4186.130  4188.730];

lfpNoise(3).basename = 'PP01_2020-06-29_13-15-57';
lfpNoise(3).periods =  [    0.000   405.121
                          574.593   580.738
                          776.577   783.793
                         1249.190  1259.480
                         1820.000  1834.550
                         2408.920  2431.980
                         2749.720  2750.170
                         2928.870  2933.310
                         3332.170  3363.330
                         4265.670  4300.380
                         4628.700  4631.810
                         4822.360  4833.520
                         5308.080  5322.810
                         7578.120  7604.120
                         8566.930  8612.710];

lfpNoise(4).basename = 'PP01_2020-06-30_12-27-33';
lfpNoise(4).periods =  [  233.150   269.487
                          605.654   753.600
                          950.984  1010.640
                         1558.930  1564.280
                         2245.670  2376.790
                         3090.710  3090.730
                         3170.470  3197.560
                         3276.060  3280.570
                         3572.120  3595.260
                         3666.350  3671.810
                         3699.420  3702.300
                         3716.650  3719.700
                         5788.690  5797.260
                         6763.780  6770.430
                         7812.150  7853.020
                         7976.020  7995.580];

lfpNoise(5).basename = 'PP01_2020-07-01_13-07-16';
lfpNoise(5).periods =  [   44.920    50.396
                          629.048   634.993
                          5431.060  5437.130
                          6450.950  6495.960];

lfpNoise(6).basename = 'PP01_2020-07-02_12-41-03';
lfpNoise(6).periods =  [  197.068   269.604
                         2133.090  2149.520
                         2377.440  2388.320
                         2644.480  2669.250
                         3167.560  3171.010
                         3212.460  3215.080
                         7582.780  7584.020
                         7747.840  7751.400
                         7862.110  7865.890
                         8472.610  8475.660
                         9370.400  9390.040
                        10072.900 10084.900
                        10513.600 10534.600
                        10616.400 10620.200];

lfpNoise(7).basename = 'PP01_2020-07-08_17-24-19';
lfpNoise(7).periods =  [   88.112    91.980
                          297.288   310.314
                          672.329   682.526
                          794.038   841.886
                         1417.730  1431.540];

lfpNoise(8).basename = 'PP01_2020-07-09_13-22-13';
lfpNoise(8).periods =  [   34.138    42.106
                          171.340   243.480
                          385.493   389.479
                          404.148   407.314
                          439.042   531.939
                         2359.200  2381.280
                         3450.460  3460.010
                         4031.290  4054.530
                         4521.000  4580.800
                         7815.020  7828.510
                         8253.180  8304.620
                         9703.090  9706.330
                        13953.600 13958.500
                        14211.600 14241.700];

lfpNoise(9).basename = 'PP01_2020-07-14_15-58-44';
lfpNoise(9).periods =  [ 1333.390  1333.820
                         2124.200  2128.180
                         3561.340  3564.380];

lfpNoise(10).basename = 'PP01_2020-07-15_15-17-43';
lfpNoise(10).periods = [  100.771   146.393
                         1215.050  1235.610
                         1878.920  1907.940
                         2983.810  3172.060
                         3338.100  3368.810
                         6617.620  6641.070
                         7132.440  7475.300];

lfpNoise(11).basename = 'PP01_2020-07-21_13-00-21';
lfpNoise(11).periods = [  102.648   107.396
                          342.067   346.680
                          986.742  1003.800
                         2700.660  2727.690
                         6717.440  6742.110
                         8578.120  8758.240
                        10445.000 10458.600
                        11099.600 11119.900];

lfpNoise(12).basename = 'PP01_2020-07-23_16-22-54';
lfpNoise(12).periods = [ 1022.030  1028.530
                         1138.220  1140.270
                         1491.500  1508.650
                         2604.310  2636.930
                         3758.810  3798.190];

lfpNoise(13).basename = 'PP02_2020-07-07_14-02-12';
lfpNoise(13).periods = [  812.442  1236.050
                         1323.950  1873.950
                         2021.300  2185.810
                         2252.080  2257.050
                         2389.630  2407.790
                         3069.930  3412.620
                         4565.840  4800.950
                         5037.070  5238.000
                         5436.770  5440.890
                         8127.830  8280.030
                         8422.220  8462.260
                         8537.460  8577.350
                         9467.820  9536.480
                        10853.300 10856.500];

lfpNoise(14).basename = 'PP02_2020-07-14_12-50-12';
lfpNoise(14).periods = [  896.613   906.902
                         1975.640  1985.820
                         3438.770  3441.740
                         3811.960  3817.450
                         4149.650  4230.030
                         5111.560  5114.670
                         5287.540  5311.420
                         5701.690  5706.310
                         5806.880  5871.060
                         6274.710  6335.300
                         6743.580  6766.870
                         6891.600  6929.470];

lfpNoise(15).basename = 'PP02_2020-07-17_14-18-20';
lfpNoise(15).periods = [    0.000   809.262
                         1652.950  1655.490
                         1926.480  1960.630
                         6727.610  6900.390
                         7565.610  7574.140
                         8306.040  8369.740
                         8708.220  8820.850
                         9226.230  9394.780
                        10929.500 10939.500];

lfpNoise(16).basename = 'PP02_2020-07-20_13-15-54';
lfpNoise(16).periods = [   45.785    77.283
                          698.077   754.592
                         2656.840  2741.530
                         2998.400  3004.690
                         3155.830  3180.930
                         8250.650  8297.510
                         8766.810  8849.110];

lfpNoise(17).basename = 'PP02_2020-07-22_14-38-24';
lfpNoise(17).periods = [   39.830   144.170
                          450.884   456.352
                          748.636   801.614
                          846.070   848.605
                         1002.530  1013.040
                         1097.050  1107.120
                         6103.320  6403.920
                         6526.800  6539.610
                         8870.740  8879.650
                         9922.390  9963.400
                        11345.500 11519.400
                        12570.000 12604.300
                        12809.900 12814.000
                        14080.200 14089.300];

lfpNoise(18).basename = 'PP02_2020-07-23_12-12-45';
lfpNoise(18).periods = [  346.056   474.647
                         1188.200  1421.460
                         1910.300  1955.190
                         2352.830  2414.530
                        10407.800 10418.700
                        10888.600 10951.300
                        11680.600 11714.800
                        12561.200 12576.500
                        12698.700 12739.200];

%% Scalar parameters
% PGD significance cutoff
%oscSignificanceCutoff = [5.74266821989813; 5.94095602315232];
oscSignificanceCutoff = [3; 3];