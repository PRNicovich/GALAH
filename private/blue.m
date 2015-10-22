function B = blue(m)
%BLUE   Blue Colormap
%   BLUE(M) is an M-by-3 matrix colormap for increasing red intensity.
%   BLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   See also GREEN, RED, JET, HSV, HOT, PINK, FLAG, COLORMAP, RGBPLOT.


if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

B = zeros(m,3);
B(:,3) = (0:(1/(m-1)):1);