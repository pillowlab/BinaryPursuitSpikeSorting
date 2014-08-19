function xc = circxcorr(x1, x2, maxlag, str)
% xc = circxcorr(x1, x2, maxlag, str);  %% cross-correlation of x1,x2
%   or
% xc = circxcorr(x1,maxlag,str);  %% auto-correlation
%
% Performs convolution of x1 and x2 with circular boundary conditions (in Fourier domain) 
%
% jw pillow 8/19/2014


% parse inputs
slen = length(x1);
switch nargin
    case 1,  % Full auto-correlation
        nargs=1; LAGflag=0; str=[];
    case 2,
        if isnumeric(x2) & (length(x2)==1), % 1 arg, with maxlag
            nargs=1; LAGflag=1; maxlag=x2; str=[];
        elseif ischar(x2) % 1 arg, no lag, scaling
            nargs=1; LAGflag=0; str=x2;
        elseif length(x2) == slen % 2 args
            nargs=2; LAGflag=0; str=[];
        else
            error('unrecognized inputs');
        end
    case 3,
        if isnumeric(x2) & (length(x2)==1), % 1 arg, with maxlag
            nargs = 1; LAGflag=1; str=maxlag; maxlag=x2;
        elseif isnumeric(maxlag) % 1 arg, no lag, scaling
            nargs=2; LAGflag=1; str=[];
        elseif ischar(maxlag)
            nargs=2; LAGflag=0; str=maxlag;
        else
            error('unrecognized inputs');
        end
    case 4,
        nargs = 2; LAGflag=1;
end

if (nargs == 1)
    xc = ifft(abs(fft(x1)).^2);
else
    xc = ifft(conj(fft(x1)).*fft(x2));
end

if (LAGflag)
    ii = [slen-maxlag+1:slen, 1:maxlag+1];
    xc = xc(ii);
end

if ~isempty(str)
    if strcmp(str,'unbiased');
        xc = xc./(slen-1);
    elseif strcmp(str,'biased')
        xc = xc./slen;
    elseif strcmp(str,'none');
    else
        str
        error('unrecognized SCALE OPTION (str) for circxcorr.m');
    end
end

