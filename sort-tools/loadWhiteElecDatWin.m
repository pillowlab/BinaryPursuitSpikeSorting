function Ydat = loadWhiteElecDatWin(twin,filenamestring,nsampsperchunk)
% Ydat = loadRawElecDatWin(twin,filenamestring,nsampsperchunk)
%
% Loads whitened electrode data for all electrodes in given time window
%
% Input:
%    twin = [t0,t1].  Loads all samples from index (t0+1) to (t1);
%    filenamestring = filename to load (contains '%d' for chunk number)
%    nsampsperchunk = number of samples per file
%
% Note: this simple version loads a single file, but for longer experiments, users should
% write their own function so that the user passes in the desired time index range, and the function is
% clever about loading the relevant data files and stitching them together (if needed).
%
% jw pillow 8/18/2014

slen = diff(twin); % number of time bins (rows in matrix)
ichunk1 = floor(twin(1)/nsampsperchunk)+1; % index for first chunk
ichunkN = ceil(twin(2)/nsampsperchunk); % index for last chunk
i1 = twin(1)-((ichunk1-1)*nsampsperchunk)+1; % index of first sample in 1st chunk
iN = twin(2)-((ichunkN-1)*nsampsperchunk); % index of last sample in last chunk

% Load first chunk
Ydat = struct2array(load(sprintf(filenamestring,ichunk1)));

% --- Determine if we need multiple chunks -------
if ichunk1 == ichunkN
        
    % Remove rows if necessary
    if (i1>1) || (iN < nsampsperchunk)
        Ydat = Ydat(i1:iN,:);
    end
    
else   %  ---- Concatenate multiple chunks together ----
  
    % Remove initial rows, if necessary
    if i1>1
        Ydat = Ydat(i1:end,:);
    end
    
    % Allocate space for remaining samples
    [n1,nelec] = size(Ydat);  % number of samples and number of electrodes
    Ydat = [Ydat; zeros(slen-n1,nelec)];
    
    % Load additional chunks
    for ichunk = ichunk1+1:ichunkN-1
        Ychnk = struct2array(load(sprintf(filenamestring,ichunk)));
        Ydat(n1+(ichunk-ichunk1-1)*nsampsperchunk+(1:nsampsperchunk),:) = Ychnk;
    end

    % Load last chunk
    Ychnk = struct2array(load(sprintf(filenamestring,ichunkN)));
    Ydat(n1+(ichunkN-ichunk1-1)*nsampsperchunk+1:end,:) = Ychnk(1:iN,:);
    
end
