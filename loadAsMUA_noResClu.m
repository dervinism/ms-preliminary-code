% A helper function to AnPSD_load.
% Inputs: dirname - data folder
%         basefilename - binary data file
%         probeFile - probe configuration file (forPRB*)
%         shanks - which shanks/tetrodes to load, e.g. 2:2:6
%         sr - sampling rate (in Hz), 30000 by default
%         binSize - binSize of the raster (in seconds, 1 by default)
%         opt -
%             .removeLastBin - true by default (because such bin typically spans beyond the end of recording
%             .SUAonly - false by default, if true only SUAs (i.e. cluster no. >=2 will be used)
% Output: raster - a shanks by time (2D) matrix
function raster = loadAsMUA_noResClu(dirname, basefilename, probeFile, shanks, shCh, chOI, sr, binSize, opt)
if nargin < 7
  sr = 3e4;
end
if nargin < 8
  binSize = 1;
end
if nargin < 9 || ~isfield(opt, 'removeLastBin') 
  opt.removeLastBin = true;
end


maxSample = 0;
for s = shanks
  if exist([basefilename, '.res.', num2str(s)], 'file')
    res = load([basefilename, '.res.', num2str(s)]);
  else
    [~, res] = resCluFromKilosort(dirname, s, shCh, chOI, probeFile);
  end
  if ~isempty(res) % Any units on this shank?
    maxSample = max(maxSample, res(end));
  end
end
binSize = round(sr*binSize); % translate to samples
assert(binSize >= 1, 'binSize less than sampling rate allows?!')
raster = zeros(ceil(maxSample/binSize)*binSize, numel(shanks));

for s = 1:numel(shanks)
  if ~opt.templates
    if exist([basefilename, '.clu.', num2str(s)], 'file') && exist([basefilename, '.res.', num2str(s)], 'file')
      clu = load([basefilename, '.clu.', num2str(s)]);
      clu = clu(2:end);
      res = load([basefilename, '.res.', num2str(s)]);
      positions = unitPos(dirname, shCh, probeFile);
      unitsOI = positions(logical(sum(positions(:,2) == chOI,2)),1)';
      res = res(logical(sum(clu == unitsOI,2)));
      clu = clu(logical(sum(clu == unitsOI,2)));
    else
      [clu, res] = resCluFromKilosort(dirname, s, shCh, chOI, probeFile);
      clu = clu(2:end);
    end
  else
    if exist([basefilename, '.tmpl.', num2str(s)], 'file') && exist([basefilename, '.res.', num2str(s)], 'file')
      clu = load([basefilename, '.tmpl.', num2str(s)]);
      uClu = double(unique(clu));
      res = load([basefilename, '.res.', num2str(s)]);
      positions = unitPos(dirname, shCh, probeFile);
      unitsOI = positions(logical(sum(positions(:,1) == uClu',2)) & logical(sum(positions(:,2) == chOI,2)),1)';
      res = res(logical(sum(clu == unitsOI,2)));
      clu = clu(logical(sum(clu == unitsOI,2)));
    else
      [~, res, clu] = resCluFromKilosort(dirname, s, shCh, chOI, probeFile);
    end
  end
  assert(numel(res) == numel(clu))  
  res = res(clu > 0); % throw out 'noise' spikes
  n = numel(res);
  res = sort(res, 'ascend'); % just a precaution, this is typically the case 
  raster(res, s) = 1:numel(res); % Even if res contains repeats, raster(res, s) will show the right cumsum of spikes till that point.
  res = unique(res);
  raster(res, s) = diff([0; raster(res, s)]);
  assert(n == sum(raster(res, s)))
end
 
if binSize > 1
  if numel(shanks) == 1
    raster = squeeze(sum(reshape(raster, binSize, [])));
  else
    raster = squeeze(sum(reshape(raster, binSize, [], numel(shanks))))';
  end
else
  raster = raster';
end

if opt.removeLastBin
  raster = raster(:, 1:end-1);
end