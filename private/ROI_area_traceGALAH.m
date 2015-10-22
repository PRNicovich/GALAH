%%%%%%%%%% ROI_area_traceGALAH %%%%%%%%%%%%%%%
% Call as ROI_area_traceGALAH(Image_stack, WhichChannels, PixelTraceCoordinates,
% ROIWidth, NumSteps, StepSmoothingFactor, EndCapsYorN, ShowPlotsYorN);

% Originally guts from roi.m from
%
% Shanrong Zhang
% Department of Radiology
% University of Washington
% 02/09/2004
%
% email: zhangs@u.washington.edu
%
% Changed significantly to be 
% streamlined to work better inside GALAH by PRN and Toan Nguyen

function ROI_struct = ROI_area_traceGALAH(varargin)

if nargin == 7
    
    Image_stack = varargin{1};
    channels = varargin{2};
    points = varargin{3};
    wide = varargin{4};
    num_steps = varargin{5};
    SmoothingFactor.span = 0;
    SmoothingFactor.degree = 0;
    end_caps = varargin{6};
    display_plots = varargin{7};
    
elseif nargin == 8
    
    Image_stack = varargin{1};
    channels = varargin{2};
    points = varargin{3};
    wide = varargin{4};
    num_steps = varargin{5};
    SmoothingFactor = varargin{6};
    end_caps = varargin{7};
    display_plots = varargin{8};
    
else
    disp(' ')
    disp('  Number of arguments is incorrect !!!')
    disp(' ')
    help ROI_area_traceGALAH
    
    return
end


clear Segs;

img = Image_stack(:,:,1,channels(1));

imh = findobj(0, 'Type', 'Image');

if ischar(wide)
    wide = str2double(wide);
end

if ischar(num_steps)
    num_steps = str2double(num_steps);
end

if iscell(wide)
    wide = double(wide{1});
end

if iscell(num_steps)
    num_steps = double(num_steps{1});
    
end

if display_plots == 1
    figure(2);
    hold off
    imagesc(img)
    set(gca,'XDir','normal');   % set as 'reverse' or 'normal' to change
    set(gca, 'YDir', 'reverse');% orientation of image to match other plots
    hold on
    plot(points(:,1), points(:,2), 'w');
end


nroi = 1;
outfn = 0;
% generate a jet colormap according to nroi
cmap = jet(1);
rndp = randperm(1);


croi = 1;

ROI_stack = zeros(size(img, 1), size(img, 2), 1);

Linedist = 0;
% slope = [];
% slope_x = [];
% slope_y = [];

if SmoothingFactor.span > 0
    
    % Smooth out using Savitzky-Golay method
    points(:,1) = smooth(points(:,1), SmoothingFactor.span, 'sgolay', SmoothingFactor.degree);
    points(:,2) = smooth(points(:,2), SmoothingFactor.span, 'sgolay', SmoothingFactor.degree);
    
end


for k = 1:(size(points, 1)-1)
    
    Linedist = [Linedist; sqrt( (points(k+1,1) - points(k,1)).^2 + (points(k+1,2) - points(k,2)).^2)];
%     slope = [slope ;( (points(k+1,2) - points(k,2)) / (points(k+1,1) - points(k,1)) )];
%     
%     slope_x = [slope_x; (points(k+1,2) - points(k,2))];
%     slope_y = [slope_y; (points(k+1,1) - points(k,1))];
    
end


step = sum(Linedist)/num_steps;


% If end_caps == 1, add end and extra point off either end of the
% points list.  This point is a straight shot off the end of distance = step_dist with slope
% equal to the slope between the last two (or first two) points.

if end_caps == 1;
    
    caboose_slope = (points(end,2) - points(end-1,2)) / (points(end,1) - points(end-1,1));
    
    engine_slope = (points(2,2) - points(1,2)) / (points(2,1) - points(1,1));
    
    engine_position = [(points(1,1) + cos(atan2((points(1,2) - points(2,2)),(points(1,1) - points(2,1))))*step), ...
        (points(1,2) + sin(atan2((points(1,2) - points(2,2)),(points(1,1) - points(2,1))))*step)];
    
    caboose_position = [(points(end,1) + cos(atan2((points(end,2) - points(end-1,2)), (points(end,1) - points(end-1,1))))*step), ...
        (points(end,2) + sin(atan2((points(end,2) - points(end-1,2)), (points(end,1) - points(end-1,1))))*step)];
    
    points = [engine_position; points; caboose_position];
    
    %         if size(points,1) > num_steps
    %             num_steps = size(points,1);
    %         end
    
    %plot(points(:,1), points(:,2), 'w-o');
    
end



%dist = zeros((length(x) - 1), numel(img));



x = points(:,1);
y = points(:,2);

Linedist = 0;
slope = [];
slope_x = [];
slope_y = [];


for k = 1:(size(points, 1)-1)
    
    Linedist = [Linedist; sqrt( (points(k+1,1) - points(k,1)).^2 + (points(k+1,2) - points(k,2)).^2)];
    slope = [slope ;( (points(k+1,2) - points(k,2)) / (points(k+1,1) - points(k,1)) )];
    
    slope_x = [slope_x; (points(k+1,2) - points(k,2))];
    slope_y = [slope_y; (points(k+1,1) - points(k,1))];
    
end



if display_plots == 1
    figure(2)
    hold on
    plot(engine_position(1), engine_position(2), 'ws');
    plot(caboose_position(1), caboose_position(2), 'wx');
end





theta = atan2(slope_x, slope_y);
theta = [theta; theta(end)]; % Pad with a trailing value for later
width = wide;


step = sum(Linedist)/num_steps;
k_steps = 0:step:sum(Linedist);
cum_dist = cumsum(Linedist);

inv_slope = -1./(slope);
inv_theta = -atan2(slope_y, slope_x);
inv_theta = [inv_theta; inv_theta(end)];

Segs_tot = zeros(numel(k_steps), 7);
%Segs_tot(:,3) = k_steps;



for k = 1:numel(k_steps)
    
    if k == 1
        
        x_center = x(1);
        y_center = y(1);
        
        past_point = find(le(cum_dist, k_steps(k)), 1, 'last');
        dist_past = k_steps(k) - cum_dist(past_point);
        
    else
        
        past_point = find(le(cum_dist, k_steps(k)), 1, 'last');
        dist_past = k_steps(k) - cum_dist(past_point);
        
        x_center = x(past_point) + dist_past*cos(theta(past_point));
        y_center = y(past_point) + dist_past*sin(theta(past_point));
        
    end
    
    width_y = width*sin(inv_theta(past_point));
    width_x = width*cos(inv_theta(past_point));
    
    step_down = [(x_center + width_x), (y_center + width_y)];
    
    step_up = [(x_center - width_x), (y_center - width_y)];
    
    Segs_tot(k, 1) = x_center;
    Segs_tot(k, 2) = y_center;
    Segs_tot(k, 3) = dist_past;
    Segs_tot(k, [4 5]) = step_up;
    Segs_tot(k, [6 7]) = step_down;
    Segs_tot(k, 8) = k_steps(k);
    
end

for chan = min(channels):max(channels);
    
    img = Image_stack(:,:,1,chan);
    
    mean_mask = zeros((size(Segs_tot, 1) - 1), 1);
    std_mask = zeros((size(Segs_tot, 1) - 1), 1);
    x_center = zeros((size(Segs_tot, 1) - 1), 1);
    y_center = zeros((size(Segs_tot, 1) - 1), 1);
    area_mask = zeros((size(Segs_tot, 1) - 1), 1);
    
    [X_mesh, Y_mesh] = meshgrid(1:size(img, 1), 1:size(img, 2));
    %ROI_stack = zeros(size(img, 1), size(img, 2), (size(Segs_tot, 1)-1));
    
    if display_plots == 1
        figure(2)
        hold on
    end
    
    mean_mask = zeros(size(Segs_tot,1)-1, size(Image_stack, 3));
    std_mask = zeros(size(Segs_tot,1)-1, size(Image_stack, 3));
    
    for pol = 1:(size(Segs_tot, 1) - 1);
        
        x_mask = [Segs_tot(pol, 4) Segs_tot(pol+1, 4) Segs_tot(pol+1, 6) Segs_tot(pol, 6) Segs_tot(pol, 4)];
        y_mask = [Segs_tot(pol, 5) Segs_tot(pol+1, 5) Segs_tot(pol+1, 7) Segs_tot(pol, 7) Segs_tot(pol, 5)];
        
        mask = double(poly2mask(x_mask, y_mask, size(img, 1), size(img, 2)));
        mask(mask == 0) = NaN;
        
        filter_I = mask.*double(img);
        mean_mask(pol, :) = mean(filter_I(~isnan(filter_I))) * ones(1, size(Image_stack, 3));
        std_mask(pol, :) = std(filter_I(~isnan(filter_I))) * ones(1, size(Image_stack, 3));
        
        %{
            for m = 1:size(Image_stack, 3)
                
                filter_I = mask.*double(img);
                mean_mask(pol, m) = mean(filter_I(~isnan(filter_I)));
                std_mask(pol, m) = std(filter_I(~isnan(filter_I)));
                
            end
        %}
        
        [geom, ~, ~] = polygeom(x_mask, y_mask);
        
        x_center(pol) = geom(2);
        y_center(pol) = geom(3);
        area_mask(pol) = geom(1);
        
        %ROI_stack(:,:,pol) = inpolygon(X_mesh, Y_mesh, x_mask, y_mask);
        
        if display_plots == 1
            plot(x_mask', y_mask', 'w');
            plot(x_center, y_center, 'w.');
        end
        
    end
    
    center_dist = sqrt((repmat(x_center,1, length(x_center)) - repmat(x_center',length(x_center), 1)).^2 + ...
        (repmat(y_center,1, length(y_center)) - repmat(y_center',length(y_center), 1)).^2);
    
    ROI_out(chan).x_center = x_center;
    ROI_out(chan).y_center = y_center;
    ROI_out(chan).center_dist = cumsum([0; diag(center_dist,1)]);
    ROI_out(chan).center_spacing = [0; diag(center_dist,1)];
    ROI_out(chan).Segs_tot = Segs_tot;
    ROI_out(chan).areas = area_mask;
    %ROI_out.ROI_stack = ROI_stack;
    ROI_out(chan).Mean_mask = mean_mask;
    ROI_out(chan).std_mask = std_mask;
    ROI_out(chan).points = points;
    
end


ROI_struct = ROI_out;