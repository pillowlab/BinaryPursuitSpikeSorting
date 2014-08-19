function ydat = loadElecDatWin(twin,eefilename,NE);
% ydat = loadElecDatWin(twin,eefilename,NE);
%
% Loads electrode data for all electrodes in given time window
%
% Input:
%    twin = [t0,t1].  Loads all samples from number t0+1 to t1;
%    eefilename = filename to load.  (Should have a %d in filename for
%        block number).
%
% Note: relies on value of global variable nSampsPerFile, which specifies
% the number of time samples stored in each ydat file.

global nSampsPerFile;

slen = diff(twin);
blk1 = ceil((twin(1)+1)/nSampsPerFile);  % 1st block to load
blkn = ceil((twin(2))/nSampsPerFile);    % last block to load

ydat = zeros(slen,NE);
icum = 0;  % cumulative length of ydat
istrt = mod(twin(1),nSampsPerFile)+1; %starting index in electrode data
for j = blk1:blkn
    istp = min(j*nSampsPerFile,twin(2)) - (j-1)*nSampsPerFile;
    ii = istrt:istp;
    ilen = length(ii);
    
    % Load file and extract relevant dat
    fload = sprintf(eefilename, j);
    xx = load(fload);
    ydat(icum+1:icum+ilen,:) = xx.ydat(ii,:);

    % Set for next block
    icum = icum+ilen;
    istrt = 1;
end

