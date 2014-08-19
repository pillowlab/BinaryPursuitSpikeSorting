function [xxEst, YW]  = estimSps_binarizeScalars(xx,YW,Wproj,WW);
%  [xxEst, YW]  = estimSps_binarizeScalars(xx,YW,Wproj,WW);
%
%  Binarizes the scalar solution obtained by estimation with a Laplace
%  prior.  Works by greedily flipping bins to 0 or 1 based on which moves
%  incur the smallest likelihood penalty.
%
%  Inputs:  
%    xx = initial guess at spike trains
%    YW = central piece needed to compute posterior (includes prior)
%    Wproj = inner product of waveform for each cell with all
%             shifted copies of other waveforms  
%    WW = waveform tensor
%
% Outputs:  
%   xxEst = estimated (scalar-valued) spike train
%   YW = updated values

% initialize params
[nw,ne,nc] = size(WW);
slen = size(xx,1);
nw2 = nw/2;

% Shrink xx to relevant window
xx = xx(nw2+1:end-nw2+1,:);
xlen = size(xx,1);

wnorm = diag(squeeze(Wproj(nw,:,:)));  % dot product of each waveform with itself
iirel = [-nw+1:nw-1];  % indices whose xhat affected by inserting a spikeiirel = [-nw+1:nw-1];  % indices whose xhat affected by inserting a spike
wnormmat = repmat(wnorm',length(iirel),1);

% Compute reduction in logli from setting relevant bins to 0 or 1
dxsp = YW.*(1-xx) - .5*repmat(wnorm',xlen,1).*(1-xx.^2);
dxnosp = YW.*(-xx) - .5*repmat(wnorm',xlen,1).*(-xx.^2);
maxdlogli = max(dxsp,dxnosp);

xxnotset = sparse(xlen,nc);
xxnotset((xx>0)&(xx<1)) = 1;
iinotset = find(xxnotset);

fprintf('Binarizing scalar spike train estimate...\n');

nremaining = length(iinotset);
while (nremaining>0) 
    [mx,iimx] = max(maxdlogli(iinotset));
    cellsp = floor((iinotset(iimx)-1)/xlen)+1;
    ispk = iinotset(iimx)-((cellsp-1)*xlen);

    if dxsp(ispk,cellsp) > dxnosp(ispk,cellsp)
        xnew = 1;
    else
        xnew = 0;
    end
    dx = xnew -xx(ispk,cellsp);
    xx(ispk,cellsp) = xnew;

    iichg = ispk+iirel;
    ywctr = YW(ispk,cellsp);
    if ((iichg(1)) >= nw) & ((iichg(end)) <= xlen)
        % update xhat for all nearby bins
        YW(iichg,:) =  YW(iichg,:)-Wproj(:,:,cellsp)*dx;
        YW(ispk,cellsp) = ywctr;

        dxsp(iichg,:) = YW(iichg,:).*(1-xx(iichg,:)) - ...
            .5*wnormmat.*(1-xx(iichg,:).^2);
        dxnosp(iichg,:) = YW(iichg,:).*(-xx(iichg,:)) - ...
            .5*wnormmat.*(-xx(iichg,:).^2);        
        maxdlogli(iichg,:) = max(dxsp(iichg,:),dxnosp(iichg,:));

    else   % do index math to set only relevant indices
        ii = find((iichg>=1) & (iichg<=xlen));
        YW(iichg(ii),:)= YW(iichg(ii),:)-Wproj(ii,:,cellsp)*dx;
        YW(ispk,cellsp) = ywctr;

        dxsp(iichg(ii),:) = YW(iichg(ii),:).*(1-xx(iichg(ii),:)) - ...
            .5*wnormmat(ii,:).*(1-xx(iichg(ii),:).^2);
        dxnosp(iichg(ii),:) = YW(iichg(ii),:).*(-xx(iichg(ii),:)) - ...
            .5*wnormmat(ii,:).*(-xx(iichg(ii),:).^2);
        maxdlogli(iichg(ii),:) = max(dxsp(iichg(ii),:),dxnosp(iichg(ii),:));

    end
    
    iinotset(iimx) = [];
    nremaining = nremaining-1;

end

% -----------------------------------------------------------------
% Finally, check that no bins need to be flipped after binarization 

nflipped = 0;
while  max(max(maxdlogli)) > 0
    [mxvals,iimx] = max(maxdlogli);
    [mx,cellsp] = max(mxvals);
    ispk = iimx(cellsp);
    if dxsp(ispk,cellsp) > dxnosp(ispk,cellsp)
        xnew = 1;
    else
        xnew = 0;
    end
    dx = xnew -xx(ispk,cellsp);
    xx(ispk,cellsp) = xnew;

    iichg = ispk+iirel;
    ywctr = YW(ispk,cellsp);
    ii = find((iichg>=1) & (iichg<=xlen));
    YW(iichg(ii),:)= YW(iichg(ii),:)-Wproj(ii,:,cellsp)*dx;
    YW(ispk,cellsp) = ywctr;

    dxsp(iichg(ii),:) = YW(iichg(ii),:).*(1-xx(iichg(ii),:)) - ...
        .5*wnormmat(ii,:).*(1-xx(iichg(ii),:).^2);
    dxnosp(iichg(ii),:) = YW(iichg(ii),:).*(-xx(iichg(ii),:)) - ...
        .5*wnormmat(ii,:).*(-xx(iichg(ii),:).^2);
    maxdlogli(iichg(ii),:) = max(dxsp(iichg(ii),:),dxnosp(iichg(ii),:));

    nflipped = nflipped+1;
end
 
xxEst = [sparse(nw2,nc); xx; sparse(nw2-1,nc)];

