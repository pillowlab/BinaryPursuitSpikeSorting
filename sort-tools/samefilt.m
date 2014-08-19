function [G,G1,G2] = samefilt(A, B, str)
%  [G,G1,G2] = samefilt(A, B, flag)
%   
%  Filters matrix A with each (same-width) filters contained in B.
%  
%  Input: A = tall matrix
%         B = 3-tensor with filters B(:,:,j) the same width as A
%         flag = set to 'conv' for B to be vertically flipped;
%                default is 'filt'  (optional)
%
%  Output: 
%    G = same height as A and width is size(B,3).  (Vertical
%        dimension clipped as with 'same' flag under conv2.m)
%    G1,G2 = pieces before and after if 'full' temporal conv2 desired
%
%  Notes: 
%   - Convolution performed using matrix multiplication:
%       best when B is short but fat
%
% jw pillow 8/19/2014

if (nargin <= 2)
    str = 'filt';
end

if strcmp(str, 'conv')
    B = flipdim(B,1);
elseif ~strcmp(str, 'filt')
    error('Unrecognized string: samefilt.m');
end

[am, an] = size(A);
[bm, bn,nflt] = size(B);
nn = am+bm-1;
npre = ceil((bm-1)/2);
npost = floor((bm-1)/2);

% Do convolution
G = zeros(nn,nflt);
for i = 1:nflt
    yy = A*B(:,:,i)';
    for j = 1:bm
        G(j:j+am-1,i) = G(j:j+am-1,i) + yy(:,bm-j+1);
    end
end

if nargout > 1
    G1 = G(1:npre,:);
end
if nargout > 2
    G2 = G(nn-npost+1:nn,:);
end

G = G(npre+1:nn-npost,:);
