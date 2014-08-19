function Ydat = loadRawElecDatWin(twin,filenamestring)
% Ydat = loadRawElecDatWin(twin,filenamestring)
%
% Loads electrode data for all electrodes in given time window
%
% Input:
%    twin = [t0,t1].  Loads all samples from index (t0+1) to (t1);
%    filenamestring = filename to load
%
% Note: this simple version loads a single file, but for longer experiments, users should
% write their own function so that the user passes in the desired time index range, and the function is
% clever about loading the relevant data files and stitching them together (if needed).

Ydat = struct2array(load(filenamestring));
Ydat = Ydat(twin(1)+1:twin(2),:);
