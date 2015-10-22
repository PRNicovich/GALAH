%Vector2Colormap returns an M-by-N-by-3 matrix of colormap values corresponding
%to values in M-by-N-by-1 input_vector.  
%Call as map = Vector2Colormap(input_vector, input_map), where
%input_vector is vector data and input_map is string corresponding to
%desired colormap.  This allows for multiple colormaps to be used in a
%single figure by expressing each as a M-by-N-by-3 matrix in the desired
%colormap.  The resulting group of images can each be expressed as color
%images in their respective colormaps.
% Altered from Vector2Colormap to allow for NON-autoscaled output.
% Vector2Colormap can be used if output range is always 0-1.  Here a third argument
% has to be supplied for the CLim property.  Values in input_vector in
% excess of max(scaling) will be set to max(scaling).  Converse true with
% min(scaling).  

% Set scaling to [] to leave default scaling of [min(input) max(input)] -> [0 1] in place.



function map = Vector2Colormap_setscale(input_vector, input_map, scaling, varargin)


if size(varargin) == 0;

    
    N_steps = 256;
    
else

    
    N_steps = varargin{1};
    

end

cm = feval(input_map, N_steps);

%%%%% Find which bins each data bit goes into.

% Normalize input_vector from 1 to N_steps, in integer values only

if ~isempty(scaling)
    
    input_vector(input_vector > max(scaling)) = max(scaling);
    input_vector(input_vector < min(scaling)) = min(scaling);
    


    scaled = round((N_steps-1)*((input_vector - min(scaling))/(max(scaling) - min(scaling)))) + 1;
    
else
    
    scaled = round((N_steps-1)*((input_vector - min(input_vector(:)))/(max(input_vector(:)) - min(input_vector(:))))) + 1;

end

if isvector(input_vector) == 1;
    
    map = [cm(scaled,1) cm(scaled,2) cm(scaled,3)]; % Map scaled across [0 1]x3
    
       
else

	m = cm(scaled, :);
	map = reshape(m, size(input_vector, 1), size(input_vector, 2), 3);

  
end
    
    