function G = green(m)
%GREEN   Green Colormap
%   GREEN(M) is an M-by-3 matrix colormap for increasing red intensity.
%   GREEN, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   See also RED, JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

G = zeros(m,3);
G(:,2) = (0:(1/(m-1)):1);
