function R = red(m)
%RED   Red Colormap
%   RED(M) is an M-by-3 matrix colormap for increasing red intensity.
%   RED, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   See also JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

R = zeros(m,3);
R(:,1) = (0:(1/(m-1)):1);
