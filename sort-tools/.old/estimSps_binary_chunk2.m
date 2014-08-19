function [xx,rr,dlogli] = estimSps_binary_chunk(xx,yy,WW,pspike,eenums,wproj,wnorm);
%  [xx,rr,dlogli] = estimSps_binary_chunk(xx,yy,WW,pspike,eenums,eeIntersct,wproj,wnorm);
%
%  Binary programming solution for xx that solves
%  W*xx = Y;     
%  via greedy binary updates to xx.
%
%  Inputs:  xx = initial guess at spike trains
%           yy = raw voltages
%           WW = waveform tensor
%           pspike = prior probability of a spike in a bin for each cell
%           eenums = cell array of electrode numbers for each waveform
%           wproj = projection of w onto itself
%           wnorm = vector norm of each waveform.
%
% Outputs:  
%   xx = estimated binary spike train
%   rr = electrode residuals for estimated spike train
%   dlogli = increase in log-likelihood for spike added/removed at each bin
%
%  Generally called from estimSps_binary.m


maxSteps = 2e5;  % Max # of passes (i.e., # spikes added/removed)

% initialize params
[nw,ne,nc] = size(WW);
slen = size(xx,1);
nw2 = nw/2;

spkthresh = -log(pspike)+log(1-pspike);  % diff in pior log-li for adding sp

% Compute voltage predicted by initial guess of spike train
rr = compVpredictionSprseCols(xx,WW,eenums);  % voltage predicted by known spikes
rr = yy-rr;  % initial residuals

% Compute projection of WW onto all shifted copies of WW (if necessary) ---------
if nargin <= 5  % Means wproj and wnorm *not* passed in as params
    fprintf(['Computing Waveform auto-convolutions: better to compute ' ...
             'in advance!\n']);
    wproj = zeros(2*nw-1,nc,nc);
    for jcell = 1:nc
        wproj(:,:,jcell) = fullfilt(WW(:,:,jcell), WW);
    end
    wnorm = diag(squeeze(wproj(nw,:,:)));  % dot prod of each waveform with itself
end
    

% -----------------------------------------------------------------
% 1. Compute MSE reduction for each spike bin position

%fprintf('estimSps_binary_chunk: initializing dlogli matrix...\n');
rwcnv = validfilt_sprse(rr,WW,eenums);
dlogli = rwcnv - repmat((.5*wnorm'+spkthresh),slen-nw+1,1);
for jcell = 1:nc
    ispk = find(xx(nw2+1:slen-nw2+1,jcell));
    dlogli(ispk,jcell) = -dlogli(ispk,jcell)-wnorm(jcell);
end
dlogli = [-100*ones(nw2,nc); dlogli];


% -----------------------------------------------------------------
% 2. Now begin maximizing posterior by greedily inserting/removing spikes

%fprintf('estimSps_binary_chunk: beginning to add/subtract spikes\n');
modreport = 100;
nstp = 1; % Counts # of passes through data
iirel = [-nw+1:nw-1];  % indices whose dlogli affected by inserting a spike

[mxvals,iimx] = max(dlogli);  % Find maximum of dlogli
[mx,cellsp] = max(mxvals);
ispk = iimx(cellsp);
while (mx > 0) & (nstp <= maxSteps);

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

    if xx(ispk,cellsp) == 0   % Insert spike --------------
        
        xx(ispk,cellsp) = 1;  
        dloglictrbin = -dlogli(ispk,cellsp); % dLogli for this bin
        inds = ispk-nw2:ispk+nw2-1;
        rr(inds,eenums{cellsp}) = rr(inds,eenums{cellsp}) ...
            - WW(:,eenums{cellsp},cellsp); % remove waveform

        % update dlogli for nearby bins
        if ((ispk+iirel(1)) >= nw2+1) & ((ispk+iirel(end)) <= slen-nw2+1)
            dlogli(ispk+iirel,:)=  dlogli(ispk+iirel,:)-wproj(:,:,cellsp);
        else
            % do index math to set only relevant indices
            iichg = ispk+iirel;
            ii = find((iichg>=(nw2+1)) & (iichg<=(slen-nw2+1)));
            dlogli(iichg(ii),:)= dlogli(iichg(ii),:)-wproj(ii,:,cellsp);
        end
            
        dlogli(ispk,cellsp) = dloglictrbin; % set for this bin 

    else   % Remove spike ----------------------------------

        xx(ispk,cellsp) = 0;
        dloglictrbin = -dlogli(ispk,cellsp); % dLogli for this bin
        inds = ispk-nw2:ispk+nw2-1;
        rr(inds,eenums{cellsp}) = rr(inds,eenums{cellsp}) ...
            + WW(:,eenums{cellsp},cellsp); % add waveform back to residuals

        % update dlogli for nearby bins
        if ((ispk+iirel(1)) >= nw2+1) & ((ispk+iirel(end))<=slen-nw2+1)
            dlogli(ispk+iirel,:)= dlogli(ispk+iirel,:)+wproj(:,:,cellsp);
        else
            % do index math to set only relevant indices
            iichg = ispk+iirel;
            ii = find((iichg>=(nw2+1)) & (iichg<=(slen-nw2+1)));
            dlogli(iichg(ii),:)=  dlogli(iichg(ii),:)+wproj(ii,:,cellsp);
        end

        dlogli(ispk,cellsp) = dloglictrbin; % set for this bin
     
    end

%% NOTE (jp): this step here (finding the maximum) is a bottleneck
[mxvals,iimx] = max(dlogli);  % Find maximum of dlogli
[mx,cellsp] = max(mxvals);
ispk = iimx(cellsp);
nstp = nstp+1;

end

% Notify if MaxSteps exceeded
if nstp > maxSteps
    fprintf('estimSps_binary_chun: max # passes exceeded (dlogli=%.3f)\n', dlogli(ispk,cellsp));
end

nstp

% ==========================================
function G = validfilt_sprse(A,B,cnx);
% f = validfilt_cnx(A,B,cnx);
%
% Analagous to samefilt, but returns only "valid" part of
% convolution, and uses column-sparsity information cnx

[am, an] = size(A);
[bm, bn,nflt] = size(B);
nn = am+bm-1;
npre = bm-1;
npost = bm-1;

% Do convolution
G = zeros(nn,nflt);
for i = 1:nflt
    yy = A(:,cnx{i})*B(:,cnx{i},i)';
    for j = 1:bm
        G(j:j+am-1,i) = G(j:j+am-1,i) + yy(:,bm-j+1);
    end
end

G = G(npre+1:nn-npost,:);
