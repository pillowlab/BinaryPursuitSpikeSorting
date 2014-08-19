function [xx,nrefviols] = runBinaryPursuit(xx,yy,W,pspike,wproj,wnorm,minISI)
% [xx,nrefviols] = runBinaryPursuit(xx,yy,W,pspike,wproj,wnorm)
%
% Binary programming solution for 
%    argmin_xx  ||W*xx- Y||^2
% via greedy binary updates to xx.
%
% Inputs:  xx = initial guess at spike trains
%          yy = raw voltages
%          WW = waveform tensor
%          pspike = prior probability of a spike in a bin for each cell
%          wproj = projection of w onto itself
%          wnorm = vector norm of each waveform
%          minISI = minimum ISI (in number of samples)
%
% Outputs:  
%   xx - estimated binary spike train
%   nrefviols [1 x ncells] - number of spikes removed for violating refractory period
%
% This is the workhorse function that does the actual BP (generally called from
% estimSps_BinaryPursuit.m, which handles passing in and out chunks of data).
%
% Also includes (kludgey) final step of removing spikes that violate minimum ISI. 
%
% jw pillow 8/18/2014

maxSteps = 2e5;  % Max # of passes (i.e., # spikes added/removed)

% initialize params
[nw,~,nc] = size(W);  % size of waveform (time bins) and number of cells
slen = size(xx,1);  % number of time bins
nw2 = nw/2;  % half the time-length of the waveform

spkthresh = -log(pspike)+log(1-pspike);  % diff in pior log-li for adding sp

% Compute residual errors between electrode data and prediction from initial spike train
rr = yy-compVpredictionSprse(xx,W);  % residual errors
    
% -----------------------------------------------------------------
% 1. Compute MSE reduction for each spike bin position

%fprintf('estimSps_binary_chunk: initializing dlogli matrix...\n');
rwcnv = validfilt_sprse(rr,W);
dlogli = rwcnv - repmat((.5*wnorm'+spkthresh),slen-nw+1,1);
for jcell = 1:nc
    ispk = find(xx(nw2+1:slen-nw2+1,jcell));
    dlogli(ispk,jcell) = -dlogli(ispk,jcell)-wnorm(jcell);
end
dlogli = [-100*ones(nw2,nc); dlogli];


% -----------------------------------------------------------------
% 2. Now begin maximizing posterior by greedily inserting/removing spikes

% fprintf('estimSps_binary_chunk: beginning to add/subtract spikes\n');
% modreport = 100;
nstp = 1; % Counts # of passes through data
iirel = (-nw+1:nw-1);  % indices whose dlogli affected by inserting a spike

[mxvals,iimx] = max(dlogli);  % Find maximum of dlogli in each column
[mx,cellsp] = max(mxvals);  % Max of the maxima (get cell #)
ispk = iimx(cellsp);  % The time bin to flip (add/remove spike)

while (mx > 0) && (nstp <= maxSteps);

%     % Uncomment to print reports on progress
%     % --------------------------------------
%     if mod(nstp,modreport) == 0     % Report output only every 100 bins
%         if (xx(ispk,cellsp) == 1)
%             fprintf('step %d: REMOVED sp, cell %d, bin %d, Dlogli=%.3f\n',...
%                 nstp,cellsp,ispk,dlogli(ispk,cellsp));
%         else
%             fprintf(1, 'step %d: inserted cell %d, bin %d, Dlogli=%.3f\n',...
%                 nstp,cellsp,ispk,dlogli(ispk,cellsp));
%         end
%     end % ----------------------------------

    % ----------------------------------------
    % 2A.  Insert or remove spike in location where logli is most improved
    if xx(ispk,cellsp) == 0   % Insert spike --------------
        
        xx(ispk,cellsp) = 1;
        dloglictrbin = -dlogli(ispk,cellsp); % dLogli for this bin
        inds = ispk-nw2:ispk+nw2-1;
        rr(inds,:) = rr(inds,:)- W(:,:,cellsp); % remove waveform

        % update dlogli for all bins within +/- nw
        if ((ispk+iirel(1)) >= nw2+1) && ((ispk+iirel(end)) <= slen-nw2+1)
            iichg = ispk+iirel;
            dlogli(iichg,:)=  dlogli(iichg,:)-wproj(:,:,cellsp);

        else % update only for relevant range of indices
            iichg = ispk+iirel;
            ii = find((iichg>=(nw2+1)) & (iichg<=(slen-nw2+1)));
            iichg = iichg(ii);
            dlogli(iichg,:)= dlogli(iichg,:)-wproj(ii,:,cellsp);
        end
        dlogli(ispk,cellsp) = dloglictrbin; % set for center bin 

    else   % Remove spike ----------------------------------

        xx(ispk,cellsp) = 0;
        dloglictrbin = -dlogli(ispk,cellsp); % dLogli for this bin
        inds = ispk-nw2:ispk+nw2-1;
        rr(inds,:) = rr(inds,:) + W(:,:,cellsp); % add waveform back to residuals

        % update dlogli for all bins within +/- nw
        if ((ispk+iirel(1)) >= nw2+1) && ((ispk+iirel(end))<=slen-nw2+1)
            iichg = ispk+iirel;
            dlogli(iichg,:)= dlogli(iichg,:)+wproj(:,:,cellsp);
        else % update only for relevant range of indices
            iichg = ispk+iirel;
            ii = find((iichg>=(nw2+1)) & (iichg<=(slen-nw2+1)));
            iichg = iichg(ii);
            dlogli(iichg,:)=  dlogli(iichg,:)+wproj(ii,:,cellsp);
        end
        dlogli(ispk,cellsp) = dloglictrbin; % set for center bin
     
    end

    % ----------------------------------------
    % 2B. Do some index arithmetic to max maximum dlogli for each cell
    %  (Big speedup from searching for the max over all bins for each cell).
    
    % Find any cells whose prev max was in the region just changed
    iimaxchg = find((iimx>=iichg(1)) & (iimx<=iichg(end)));
    [mx0,iinw] = max(dlogli(:,iimaxchg));
    mxvals(iimaxchg) = mx0;
    iimx(iimaxchg) = iinw;
    
    % Now see if new maxima arose in region of dlogli that was just altered
    [mxvals0,iimx0] = max(dlogli(iichg,:));

    % Combine to get new maximum & position
    [mxvals,OneOrTwo] = max([mxvals0; mxvals]); % max of [changed-bin ; other-bin];
    iflip = find(OneOrTwo == 1);
    iimx(iflip) = iichg(iimx0(iflip));
    
    % Find next bin to adjust
    [mx,cellsp] = max(mxvals);
    ispk = iimx(cellsp);

    nstp = nstp+1;
end

% Notify if MaxSteps exceeded
if nstp > maxSteps
    fprintf('estimSps_binary_chun: max # passes exceeded (dlogli=%.3f)\n', dlogli(ispk,cellsp));
end

% ----------------------------------------------------------
% 3. Finally, remove spikes that violate refractory period 
nrefviols = zeros(1,nc);  % initialize counter
for jcell = 1:nc
    isis = diff(find(xx(:,jcell)));
    while any(isis<minISI)
        tsp = find(xx(:,jcell)); % spike indices
        nsp = length(tsp); % number of spikes
        badii = find(isis<minISI);
        badii = union(badii,badii(badii<nsp)+1); % indices of problem spikes
        badtsp = tsp(badii);
        [~,ii] = max(dlogli(badtsp,jcell)); % Find which spike is least likely
        ispk = badtsp(ii);
        
        % Remove this spike and update
        nrefviols(jcell) = nrefviols(jcell)+1; % count the ISI removal
        xx(ispk,jcell) = 0; % remove from spike train
        dloglictrbin = -dlogli(ispk,jcell); % dLogli for this bin
        inds = ispk-nw2:ispk+nw2-1;
        rr(inds,:) = rr(inds,:)+ W(:,:,jcell); % add waveform back to residuals

        % update dlogli for all bins within +/- nw
        if ((ispk+iirel(1)) >= nw2+1) && ((ispk+iirel(end))<=slen-nw2+1)
            iichg = ispk+iirel;
            dlogli(iichg,:)= dlogli(iichg,:)+wproj(:,:,jcell);
        else % update only for relevant range of indices
            iichg = ispk+iirel;
            ii = find((iichg>=(nw2+1)) & (iichg<=(slen-nw2+1)));
            iichg = iichg(ii);
            dlogli(iichg,:)=  dlogli(iichg,:)+wproj(ii,:,jcell);
        end
        dlogli(ispk,jcell) = dloglictrbin; % set for center bin

        % Recompute ISIs
        isis = diff(find(xx(:,jcell)));
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function G = validfilt_sprse(A,B)
% f = validfilt_sprse(A,B);
%
% Convolve sparse A with each frame of B and return only "valid" part

am = size(A,1);
[bm,~,nflt] = size(B);
nn = am+bm-1;
npre = bm-1;
npost = bm-1;

% Do convolution
G = zeros(nn,nflt);
for i = 1:nflt
    yy = A*B(:,:,i)';
    for j = 1:bm
        G(j:j+am-1,i) = G(j:j+am-1,i) + yy(:,bm-j+1);
    end
end
G = G(npre+1:nn-npost,:);
