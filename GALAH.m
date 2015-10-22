% General Automated Localization Area Handler (GALAH)
% Display TILL Photonics-exported TIFF images and extract ROIs from single-
% and dual-channel data.
% 
% PRN UNSW BMIF March 2013
% Updated July 2013-2014

% Fixes
% Begun 26082013
% 26/08/2013 - Now GALAHv05
% Resolved issue in second color of kymographs
% Now supports odd-shaped (cropped) images in backbone trace and main
% screen.
% Small issue with second axis in right axis panel - now only one.
% 07/03/2014 - Now GALAHv07
% Previous version dealt with some issues requested by PO.
% Here the exponential decay fitter feature added
% Cleaned up to allow for 2-channel data as well.
% 04042014
% Cosmetic fixes for decay fitter - added end fit frame parameter and
% changed boxes to be in engineering notation.
% 11092014 
% GALAHv09
% Fixing error in ROI generation, pseudocoloring
% GALAHv10
% Revise workflow for kinetics of filament binding + unbinding
% Kinetics button removed (associated code now vestigial)
% GALAHv11
% Added two-channel alignment module
% GALAHv12
% Added drift correction
% Changed ROI selection in Exponential Curve Fit to ROI poly
% Added BIG TIFF (*.tif8) import option
% Fixed issue with autoscale in pop-up windows
% Changed name to 'GALAH' and started adding to Bitbucket repo
% 070715 Using updated ROI_area_trace function
% 180915 Modify segmentation code so that;
%     - Trace boxes persist until Threshold change event, letting you see
%       boxes in through all frames + channels
%     - Recalc on substantial change of backbone trace parameters (width,
%     length, smoothing filter, end caps)
%     - Option and parameter for smoothing filter on positions added
%       Don't make the degree 0 or funny things happen.  
%     - Option + checkbox for segmenting based on average of whole image stack added
% 011015 Error fixes and addition of slope sensing in kymographs
%   - Update, clarify, and remove excess comments throughout.  This
%       includes 'Kinetics ROI' vestigial code where I can find it. 
%   - Fixed issue with display of ROI area in Backbone Trace plots display.
%       Error had field mis-assigned in ROI_area_traceGALAH.m
%   - Removed automatic recalc of regions on change in backbone trace
%       parameters (width, length, smoothing filter, end caps).  User must
%       now set 
%   - Error when a file has a colormap applied in ImageJ, saved as TIFF,
%       and then tried to load here.  Issue is with base MATLAB  and LibTIFF and not with
%       GALAH.  Workaround is re-open file in ImageJ, change colormap to
%       'grays', re-save, and then it should load here. 
%   - Loading speed of Exponential Curve Data Fitter window improved by
%       being smarter about how to segment data from image stack
%   - Multiple exponential fit windows can be launched.  This is done by
%       transferring all of guidata for each into the window itself instead of
%       drawing from the main guidata.  ROIs for active windows are drawn in
%       parent Exponential Decay Fit window until closed.
%   - Max intensity projection option for frame stack.  Acts analagous to
%       Use mean stack option
%   - Inverted orientation of kymograph in image such that distance down->
%       up in kymograph corresponds to left->right in ROI plot
%   - Kymographs use actual mean values in image and do not interpolate
%       data to yield smoother distance points in kymograph
%   - Structure field added to make it easier to access kymograph data


function GALAH

% Close out previous windows so no two are open at same time
close(findobj('Tag', 'TIFF viewer'));


scrsz = get(0,'ScreenSize');

Window_size_big = [60 190 1151 678]; % Desired window size in pixels
Window_size_small = [20 50 1036 614];
if scrsz(3) < 1300 || scrsz(4) < 800
    Window_size = Window_size_small;
else
    Window_size = Window_size_big;
end


fig1 = figure('Name','GALAH Image Viewer', 'Tag', 'TIFF viewer', 'Units', ...
    'normalized','Position',[Window_size(1)/scrsz(3) Window_size(2)/scrsz(4) Window_size(3)/scrsz(3) Window_size(4)/scrsz(4)], ...
    'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure');
set(fig1, 'Color',[0.9 0.9 0.9]);

%%%%%%%%%%%%
% Set up toolbar
hToolbar = findall(fig1,'tag','FigureToolBar');
AllToolHandles = findall(hToolbar);
ToolBarTags = get(AllToolHandles,'Tag');
ToolsToKeep = {'FigureToolBar'; 'Exploration.DataCursor'; 'Exploration.Pan'; 'Exploration.ZoomOut'; 'Exploration.ZoomIn'};
WhichTools = ~ismember(ToolBarTags, ToolsToKeep);
delete(AllToolHandles(WhichTools));

% Yields figure position in form [left bottom width height].
fig1_size = get(fig1, 'Position');
set(fig1, 'DeleteFcn', @GUI_close_fcn);

handles.handles.fig1 = fig1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize GUI data
handles.Load_file = [];
handles.N_frames = 2;
handles.N_channels = 2; % 1 or 2 for single, double channel data
handles.Primary_channel = 1;
handles.Img_stack = [];
handles.Left_color = 'green';
handles.Right_color = 'red';
handles.Load_file = [];
handles.Left_invert = 0;
handles.Right_invert = 0;
handles.scrsz_pixels = get(0, 'ScreenSize');
handles.Autoscale_left = 0;
handles.Autoscale_right = 0;
handles.Min_max_left = [0 1];
handles.Min_max_right = [0 1];
handles.Display_range_left = [0 1];
handles.Display_range_right = [0 1];
handles.Display_range_ROI = [0 1];
handles.handles.border_colors = jet(64);
handles.handles.Border_randomize = randperm(64);
handles.handles.Unselected_ROI_color = [.4 1 .2];
handles.handles.Selected_ROI_color = [.4 .2 1];
handles.pkfnd_radius = 2;
handles.Dilate_radius = 4;
handles.Dilate_shape = 'disk';
handles.End_caps = 1;
handles.ROI_length = 10;
handles.ROI_width = 3;
handles.ROI_smoothing.span = 0;
handles.ROI_smoothing.degree = 0;
handles.ROI_zoom_border = 50;
handles.ROIsAreSet = false;
handles.FilterType = 'disk';
handles.FilterBkgdRadius = [15 5];
handles.UseFlatKinetics = 0;
handles.Min_Traj_length = 3;
handles.StartFrame = 1;
handles.AnticipatedFrameSize = [512 512]; % Native size of image file

handles.TransformType = 'rigid';
handles.TransformFixedChannel = 1;
handles.AffineMatrix = eye(3,3);
handles.TransformMobileChannel = 2;
handles.TransformRefFrame = 1;
handles.TransformRefFrameCustom = 0;

handles.polyROI = [];

handles.ROI_boundaries.B = [];
handles.ROI_boundaries.L = [];
handles.ROI_boundaries.bound_square = [];
handles.ROI_output = [];
handles.ROI_parameters = [];
handles.ROI_Mask = [];

handles.ExpFits.A(1,1) = 1;
handles.ExpFits.A(2,1) = 0.3;
handles.ExpFits.A(3,1) = 0.1;
handles.ExpFits.tau(1,1) = 1;
handles.ExpFits.tau(2,1) = 0.1;
handles.ExpFits.tau(3,1) = 0.01;
handles.ExpFits.Nterms(1) = 2;
handles.ExpFits.FloatBkgd{1} = 1;
handles.ExpFits.ExpFitBkgdValue(1) = 0;
handles.ExpFits.FrameTime = 0.05;
handles.ExpFits.FixFit = zeros(3, 2);


handles.ExpFits.A(1,2) = 1;
handles.ExpFits.A(2,2) = 0.3;
handles.ExpFits.A(3,2) = 0.1;
handles.ExpFits.tau(1,2) = 1;
handles.ExpFits.tau(2,2) = 0.1;
handles.ExpFits.tau(3,2) = 0.01;
handles.ExpFits.Nterms(2) = 2;
handles.ExpFits.FloatBkgd{2} = 1;
handles.ExpFits.ExpFitBkgdValue(2) = 0;
handles.ExpFits.WindowCount = 0;

handles.DriftRefFrame = 1;
handles.DriftSmoothing = 1;
handles.DriftFrameRange = [1 2];
handles.DriftMatrix = [0, 0];
handles.DriftMatrixSmoothed = [0, 0];

guidata(fig1, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define panels.  

fig1_size_pixels = fig1_size.*scrsz;

panel_border = fig1_size_pixels(4)/max(fig1_size_pixels);

butt_panel = uipanel(fig1, 'Units', 'normalized', 'Position', [0 .95, 1, .05], ...
    'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'button_panel');

ax_panel1 = uipanel(fig1, 'Units', 'normalized', 'Position', [0 .1 .5 .85], ...
    'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'axes_panel1');

ax_panel2 = uipanel(fig1, 'Units', 'normalized', 'Position', [.5 .1 .5 .85], ...
    'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'axes_panel2');

slider_panel = uipanel(fig1, 'Units', 'normalized', 'Position', [0 0 1 .1], ...
    'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'slider_panel');

handles.handles.butt_panel = butt_panel;
handles.handles.ax_panel1 = ax_panel1;
handles.handles.ax_panel2 = ax_panel2;
handles.handles.slider_panel = slider_panel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define axes positions

ax1 = axes('Parent', ax_panel1, 'Position', [0.002 .005 .994 .994]);
set(ax1, 'Tag', 'Left axis');
path_here = mfilename('fullpath');

% Find logo file

if isdeployed
        logo_1 = BMIFLogoGenerate;
        fill_image = imagesc(Vector2Colormap(-logo_1,handles.Left_color), 'Parent', ax1);
        set(fill_image, 'Tag', 'fill_image_left', 'HitTest', 'on');
    
else
    logo_file = fullfile(fileparts(path_here), 'private', 'BMIF_logo.jpg');


    if exist(logo_file, 'file') == 2;

        logo_hold = single(imread(logo_file));
        logo_1 = logo_hold(:,:,1);
        clear logo_hold  
        fill_image = imagesc(Vector2Colormap(-logo_1,handles.Left_color), 'Parent', ax1);
        set(fill_image, 'Tag', 'fill_image_left', 'HitTest', 'on');

    else

        % Dummy data to put into the axes on startup
        z=peaks(1000);
        z = z./max(abs(z(:)));
        fill_image = imshow(z, 'Parent', ax1, 'ColorMap', jet, 'DisplayRange', [min(z(:)) max(z(:))]);
        set(fill_image, 'Tag', 'fill_image_left', 'HitTest', 'on');
        freezeColors(ax1);

    end
end

% Get rid of tick labels
set(ax1, 'xtick', [], 'ytick', [])
axis image % Freezes axis aspect ratio to that of the initial image - disallows skewing due to figure reshaping.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = axes('Parent', ax_panel2, 'Position', [0 0 1 1]);
set(ax2, 'Tag', 'Axis2');

if isdeployed
    
        logo_1 = BMIFLogoGenerate;
        fill_image = imagesc(Vector2Colormap(-logo_1,handles.Right_color), 'Parent', ax2);
        set(fill_image, 'Tag', 'fill_image_right', 'HitTest', 'on');
        
else

    if exist(logo_file, 'file') == 2;

        logo_hold = single(imread(logo_file));
        logo_1 = logo_hold(:,:,1);
        clear logo_hold  
        fill_image = imagesc(Vector2Colormap(-logo_1, handles.Right_color), 'Parent', ax2);
        set(fill_image, 'Tag', 'fill_image_right', 'HitTest', 'on');

    else

        % Dummy data to put into the axes on startup
        z=peaks(1000);
        z = z./max(abs(z(:)));
        fill_image = imshow(z, 'Parent', ax2, 'ColorMap', jet, 'DisplayRange', [min(z(:)) max(z(:))]);
        set(fill_image, 'Tag', 'fill_image_right', 'HitTest', 'on');
        freezeColors(ax2);

    end
end

% Get rid of tick labels
set(ax2, 'xtick', [], 'ytick', []);
axis image % Freezes axis aspect ratio to that of the initial image - disallows skewing due to figure reshaping.

handles.handles.ax1 = ax1;
handles.handles.ax2 = ax2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define button positions

%%%%%%%%%%%%%%%%%%%%%%
% Top Button panel buttons

% Button
Load_out =     uicontrol(butt_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Load Path',...
        'Position', [0 .05 .1 .9],...
        'Callback', @Load_pts, 'Tag', 'Load Path');

% Button %%%%% 
width = .1;
Image_preferences_out =     uicontrol(butt_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Image Preferences',...
        'Position', [(1 - width) .05 width .9],...
        'Callback', @Image_prefs, 'Tag', 'Image_prefs');  
    
handles.handles.Load_out = Load_out;
handles.handles.Image_preferences_out = Image_preferences_out;

%%%%%%%%%%%%%%%%%%%%%%
% Slider panel buttons

% Button
ROI_finder_out =     uicontrol(slider_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Backbone Trace',...
        'Position', [.01 .05 .1 .45], 'Enable', 'off', ...
        'Callback', @ROI_launch, 'Tag', 'ROI_launch');

    % Button
ExpFit_out =     uicontrol(slider_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Fit Exponential Decay',...
        'Position', [.23 .05 .1 .45], 'Enable', 'off',...
        'Callback', @ExpFit_launch, 'Tag', 'ExpFit_launch');
    
    % Button
Flatten_out =     uicontrol(slider_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Flatten Background',...
        'Position', [.98 - 2*width - 0.01 .05 width .45], 'Enable', 'off',...
        'Callback', @Flatten_launch, 'Tag', 'Flatten_launch');
    
        % Button
Align_out =     uicontrol(slider_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Align Channels',...
        'Position', [.98 - 3*width - 0.02 .05 width .45], 'Enable', 'off',...
        'Callback', @Align_launch, 'Tag', 'Align_launch');
    
        % Button
Drift_out =     uicontrol(slider_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Correct Drift',...
        'Position', [.98 - 4*width - 0.03 .05 width .45], 'Enable', 'off',...
        'Callback', @Drift_launch, 'Tag', 'Drift_launch');


% Button %%%%%
width = .1;
Export_out =     uicontrol(slider_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Export Results',...
        'Position', [(.98 - width) .05 width .45], 'Enable', 'off',...
        'Callback', @Export_data, 'Tag', 'Export_data'); 
    
    
handles.handles.ROI_finder_out = ROI_finder_out;
handles.handles.Flatten_out = Flatten_out;
handles.handles.Export_out = Export_out;
handles.handles.ExpFit_out = ExpFit_out;
handles.handles.Align_out = Align_out;
handles.handles.Drift_out = Drift_out;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define text box positions


Load_text = uicontrol(butt_panel, 'Style', 'edit', 'Units', 'normalized', ...
	'Position',[.11 .15 .39 .7], 'BackgroundColor', [1 1 1], ...
	'String', 'File', 'Callback', @Load_edit, 'Tag', 'Load_textbox');

handles.handles.Load_text = Load_text;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define sider positions


slider_value = 0;
slider_step = 1/(handles.N_frames-1);

slide_hand = uicontrol(slider_panel, 'Style', 'slider', 'Units', 'normalized',...  
	'SliderStep', [slider_step slider_step], 'Min', 0, 'Max', 1, 'Value', slider_value, 'Position', [.01 .67 .85 .25],...
	'Callback', @slider_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'Slider handle');
        
addlistener(slide_hand, 'Value', 'PostSet', @slider_listener);

slide_box = uicontrol(slider_panel, 'Style', 'edit', 'Units', 'normalized', ...
	'Position', [.88 .63 .1 .35], 'BackgroundColor', [1 1 1], ...
	'String', 'Frame', 'Callback', @edit_call);

set(slide_hand, 'Enable', 'off');
set(slide_box, 'Enable', 'off');



handles.handles.slide_hand = slide_hand;
handles.handles.slide_box = slide_box;

guidata(fig1, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Callback functions

%%%%%%%%%%%%%%%%%%%%%%
% Frame slider update functions

    function slider_call(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));

        set(slide_box, 'String', (1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')))));
        
        Display_images_in_axes;
        
        
    end

    function slider_listener(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        set(slide_box, 'String', (1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')))));
        
        Display_images_in_axes;
        
    end

    function edit_call(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        slide_string = str2double(get(slide_box, 'String'));
        
        % Make sure the string fed is actually a string
        
        if length(slide_string) ~= 1
            slide_set = get(slide_hand, 'Value');
            slide_str2 = round(1+slide_set*(handles.N_frames-1));
            set(slide_box, 'String', slide_str2);
            
        else
        
            slide_set = ((slide_string - 1)/(handles.N_frames - 1));
            slide_range = [get(slide_hand, 'Min') get(slide_hand, 'Max')];

            if slide_set > slide_range(2)

                slide_set = slide_range(2);
                slide_str2 = 1+slide_range(2)*(handles.N_frames-1);
                set(slide_box, 'String', num2str(slide_str2));

            elseif slide_set < slide_range(1)

                slide_set = slide_range(1);
                slide_str2 = 1+slide_range(1)*(handles.N_frames-1);
                set(slide_box, 'String', num2str(slide_range(1)));

            end
        
        end
            
        
        set(slide_hand, 'Value', slide_set);
        
        Display_images_in_axes;
        
    end




%%%%%%%%%%%%%%%%%%%%%%
% Use uigetfile to load up a file

    function Load_pts(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));

        [fname, pathname, filterindex] = uigetfile({'*.tif; *.tiff', 'TIFF file (*.tif, *.tiff)';...
                                                  '*.tf8', 'Big TIFF file (*.tf8)'});
        
        if ismember(filterindex, [1 2]);
           
        if ~strcmp(fullfile(pathname, fname), handles.Load_file)
            % Reset stuff now that there is a new file being loaded (as long as
            % it's actually new).
        
            handles.ROI_boundaries = [];
            handles.ROI_output = [];
            handles.ROI_parameters = [];
            
            
            close(findobj('Tag', 'GALAH_ROI_calc'));
            close(findobj('Tag', 'GALAH_plot_display'));

        end
        
        if ~isequal(fname, 0) && ~isequal(pathname, 0)
            
            set(handles.handles.ROI_finder_out, 'Enable', 'on');
            
            set(handles.handles.Flatten_out, 'Enable', 'on');
            set(handles.handles.Export_out, 'Enable', 'on');
 
            set(Load_text, 'String', fullfile(pathname, fname));
            handles.Load_file = fullfile(pathname, fname);
            
            load_wait = waitbar(0, 'Loading File');

            fname = handles.Load_file;
            info = imfinfo(fname);
            num_images = numel(info);
            
            Img_stack = zeros(info(1).Height, info(1).Width, num_images);
            
            for k = 1:num_images
                Img_stack(:,:,k) = imread(fname, k, 'Info', info);
                
                waitbar(k/num_images, load_wait);
            end

            waitbar(1, load_wait);
            close(load_wait)
            
            if num_images == 1;
                handles.N_channels = 1;
            end
            
            if mod(num_images, 2) == 1
                handles.N_channels = 1;
            end
                    
            handles.ROI_Mask = ones(info(1).Height, info(1).Width);

            if handles.N_channels == 1
                


                [fa, fb, fc, fd] = size(Img_stack);

                handles.Img_stack = zeros(fa, fb, fc, 1);

                handles.Img_stack = double(Img_stack);
                clear Img_stack;
                
                if fc > 1
                    
                    set(slide_hand, 'Enable', 'on');
                    set(slide_box, 'Enable', 'on');
                    
                end
                
                handles.N_frames = fc;
                
                % Pull slider value
                slide_frame = 1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')));
                
                if handles.N_frames == 1;
                    set(slide_hand, 'SliderStep', [1 1]);
                    
                else
                
                    set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                    
                end
                
                slide_frame = min([handles.N_frames slide_frame]);
                set(slide_box, 'String', num2str(slide_frame));
                
             OldXLimits = [0.5 size(handles.Img_stack, 2)+0.5];
             OldYLimits = [0.5 size(handles.Img_stack, 1)+0.5];
             set(ax1, 'XLim', OldXLimits, 'YLim', OldYLimits);
             %set(ax2, 'XLim', OldXLimits, 'YLim', OldYLimits);
             
             set(handles.handles.Align_out, 'Enable', 'off');
             
 

            elseif handles.N_channels == 2
                


                [fa, fb, fc, ~] = size(Img_stack);

                % Reshape 'cause input image is always single-channel by
                % default.
                
                handles.Img_stack = zeros(fa, fb, fc/2, 2);

                red_vect = 1:2:(fc);
                green_vect = 2:2:(fc);

                handles.Img_stack(:,:,:,1) = double(Img_stack(:,:,red_vect, 1));
                handles.Img_stack(:,:,:,2) =  double(Img_stack(:,:,green_vect, 1));
                
               if fc > 1
                    
                    set(slide_hand, 'Enable', 'on');
                    set(slide_box, 'Enable', 'on');
                    
               end
                
               handles.N_frames = fc/2;
                
                % Pull slider value
                slide_frame = 1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')));
                
                if handles.N_frames == 1;
                    set(slide_hand, 'SliderStep', [1 1]);
                    
                else
                
                    set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                    
                end
                
                slide_frame = min([handles.N_frames slide_frame]);
                set(slide_box, 'String', num2str(slide_frame));

             OldXLimits = [0.5 size(handles.Img_stack, 2)+0.5];
             OldYLimits = [0.5 size(handles.Img_stack, 1)+0.5];
             set(ax1, 'XLim', OldXLimits, 'YLim', OldYLimits);
             set(ax2, 'XLim', OldXLimits, 'YLim', OldYLimits);
             
             set(handles.handles.Align_out, 'Enable', 'on');

            end
        
        end
        
        temp_left = reshape(handles.Img_stack(:,:,:,1), [], handles.N_frames);
        handles.Min_max_left = [min(temp_left)' max(temp_left)'];
        handles.Display_range_left = [min(temp_left(:)) max(temp_left(:))];

        if handles.N_channels == 2;
            temp_right = reshape(handles.Img_stack(:,:,:,2), [], handles.N_frames);
            handles.Min_max_right = [min(temp_right)' max(temp_right)'];
            handles.Display_range_right = [min(temp_right(:)) max(temp_right(:))];
        end
        
        clear temp_left temp_right
        
        handles.ImgHold = handles.Img_stack;
        handles.MeanImg = squeeze(mean(handles.Img_stack, 3));
        handles.MaxImg = squeeze(max(handles.Img_stack, [], 3));
        handles.UseFlatKinetics = 0;
        
            if handles.N_frames > 1
%                 set(handles.handles.Unbind_out, 'Enable', 'on');
                set(handles.handles.ExpFit_out, 'Enable', 'on');
                set(handles.handles.Drift_out, 'Enable', 'on');
                handles.DriftFrameRange = [1 handles.N_frames];
                handles.DriftMatrix = zeros(handles.N_frames, 2);
                handles.DriftMatrixSmoothed = zeros(handles.N_frames, 2);
            else
%                 set(handles.handles.Unbind_out, 'Enable', 'off');
                set(handles.handles.ExpFit_out, 'Enable', 'off');
                set(handles.handles.Drift_out, 'Enable', 'off');
                set(slide_hand, 'Enable', 'off')
                set(slide_box, 'Enable', 'off')
            end
        handles.TransformRefFrame = max([handles.TransformRefFrame, size(handles.Img_stack, 3)]);

        guidata(findobj('Tag', 'TIFF viewer'), handles);
                
        Display_images_in_axes;
        
        end

	end
    
%%%%%%%%%%%%%%%%%%%%%%
% Edit text box to load file

    function Load_edit(varargin)
        
        text_input = get(Load_text, 'String');
        %disp(text_input)
        
        if exist(text_input, 'file') == 2;
            
            % If valid filename is entered...
            
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        [pathname, fname, ext] = fileparts(text_input);

        set(Load_text, 'String', fullfile(pathname, strcat(fname, ext)));
        
            set(handles.handles.ROI_finder_out, 'Enable', 'on');
            
            set(handles.handles.Flatten_out, 'Enable', 'on');
            set(handles.handles.Export_out, 'Enable', 'on');
        
           
        if ~strcmp(fullfile(pathname, fname), handles.Load_file)
            % Reset stuff now that there is a new file being loaded (as long as
            % it's actually new).
        
            handles.ROI_boundaries = [];
            handles.ROI_output = [];
            handles.ROI_parameters = [];
            
            
            close(findobj('Tag', 'GALAH_ROI_calc'));
            close(findobj('Tag', 'GALAH_plot_display'));

        end
        
        handles.Load_file = fullfile(pathname, strcat(fname, ext));
            
        load_wait = waitbar(0, 'Loading File');

        fname = handles.Load_file;
        info = imfinfo(fname);
        num_images = numel(info);
            
        Img_stack = zeros(info(1).Height, info(1).Width, num_images);
            
        for k = 1:num_images
            Img_stack(:,:,k) = imread(fname, k, 'Info', info);
                
            waitbar(num_images/k, load_wait);
        end

        waitbar(1, load_wait);
        close(load_wait)
        
         handles.ROI_Mask = ones(info(1).Height, info(1).Width);
                    
            if num_images == 1;
                handles.N_channels = 1;
            end
        
           if handles.N_channels == 1

                [fa, fb, fc, ~] = size(Img_stack);

                handles.Img_stack = zeros(fa, fb, fc, 1);

                handles.Img_stack = double(Img_stack);
                clear Img_stack;
                
                if fc > 1
                    
                    set(slide_hand, 'Enable', 'on');
                    set(slide_box, 'Enable', 'on');
                    
                end
                
                handles.N_frames = fc;
                
                % Pull slider value
                slide_frame = 1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')));
                
                %disp(slide_frame);

                
                if handles.N_frames == 1;
                    set(slide_hand, 'SliderStep', [1 1]);
                    
                else
                
                    set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                    
                end
                
                slide_frame = min([handles.N_frames slide_frame]);
                set(slide_box, 'String', num2str(slide_frame));
                
             OldXLimits = [0.5 size(handles.Img_stack, 2)+0.5];
             OldYLimits = [0.5 size(handles.Img_stack, 1)+0.5];
             set(ax1, 'XLim', OldXLimits, 'YLim', OldYLimits);
             set(handles.handles.Align_out, 'Enable', 'off');


            elseif handles.N_channels == 2
                


                [fa, fb, fc, ~] = size(Img_stack);

                % Reshape 'cause input image is always single-channel by
                % default.
                
                handles.Img_stack = zeros(fa, fb, fc/2, 2);

                red_vect = 1:2:(fc);
                green_vect = 2:2:(fc);

                handles.Img_stack(:,:,:,1) = double(Img_stack(:,:,red_vect, 1));
                handles.Img_stack(:,:,:,2) =  double(Img_stack(:,:,green_vect, 1));
                
               if fc > 1
                    
                    set(slide_hand, 'Enable', 'on');
                    set(slide_box, 'Enable', 'on');
                    
                end
                
                % Pull slider value
                slide_frame = round(handles.N_frames*get(slide_hand, 'Value'));
                
                
                %disp(slide_frame);

                handles.N_frames = fc/2;
                %set(slide_hand, 'Max', handles.N_frames);
                %set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                
                if handles.N_frames == 1;
                    set(slide_hand, 'SliderStep', [1 1]);
                    
                else
                
                    set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                    
                end
                
                slide_frame = min([handles.N_frames slide_frame]);
                set(slide_box, 'String', num2str(slide_frame));
                
              OldXLimits = [0.5 size(handles.Img_stack, 2)+0.5];
             OldYLimits = [0.5 size(handles.Img_stack, 1)+0.5];
             set(ax1, 'XLim', OldXLimits, 'YLim', OldYLimits);
             set(ax2, 'XLim', OldXLimits, 'YLim', OldYLimits);
             set(handles.handles.Align_out, 'Enable', 'on');
                
           end
           
        temp_left = reshape(handles.Img_stack(:,:,:,1), [], handles.N_frames);
        handles.Min_max_left = [min(temp_left)' max(temp_left)'];
        handles.Display_range_left = [min(temp_left(:)) max(temp_left(:))];

        if handles.N_channels == 2;
            temp_right = reshape(handles.Img_stack(:,:,:,2), [], handles.N_frames);
            handles.Min_max_right = [min(temp_right)' max(temp_right)'];
            handles.Display_range_right = [min(temp_right(:)) max(temp_right(:))];
        end
        
        clear temp_left temp_right
        handles.ImgHold = handles.Img_stack;   
        handles.UseFlatKinetics = 0;
        
            if handles.N_frames > 1
                set(handles.handles.Unbind_out, 'Enable', 'on');
                set(handles.handles.ExpFit_out, 'Enable', 'on');
                set(handles.handles.Drift_out, 'Enable', 'on');
                handles.DriftFrameRange = [1 handles.N_frames];
                handles.DriftMatrix = zeros(handles.N_frames, 2);
                handles.DriftMatrixSmoothed = zeros(handles.N_frames, 2);
            else
                set(handles.handles.Unbind_out, 'Enable', 'off')
                set(handles.handles.ExpFit_out, 'Enable', 'off');
                set(handles.handles.Drift_out, 'Enable', 'off');
                set(slide_hand, 'Enable', 'off')
                set(slide_box, 'Enable', 'off')
            end
        handles.TransformRefFrame = max([handles.TransformRefFrame, size(handles.Img_stack, 3)]);
        guidata(findobj('Tag', 'TIFF viewer'), handles);
        Display_images_in_axes;
        
        end
            

    end

%%%%%%%%%%%%%%%%%%%%%%
% Display images in axes.  Used by multiple calls in GUI.

    function Display_images_in_axes(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        %disp(handles);
        
        if isempty(handles.Load_file) % No data loaded, just dummy images
            
            ax1 = handles.handles.ax1;
            ax2 = handles.handles.ax2;

            path_here = mfilename('fullpath');

            if isdeployed
                    logo_1 = BMIFLogoGenerate;
                    fill_image = imagesc(Vector2Colormap(-logo_1,handles.Left_color), 'Parent', ax1);
                    fill_image2 = imagesc(Vector2Colormap(-logo_1,handles.Right_color), 'Parent', ax2);
                    set(fill_image, 'Tag', 'fill_image_left', 'HitTest', 'on');
                    set(fill_image2, 'Tag', 'fill_image_right', 'HitTest', 'on');
            else
                logo_file = fullfile(fileparts(path_here), 'BMIF_logo.jpg');

                if exist(logo_file, 'file') == 2;

                    logo_hold = single(imread(logo_file));
                    logo_1 = logo_hold(:,:,1);
                    clear logo_hold  
                    fill_image = imagesc(Vector2Colormap(-logo_1,handles.Left_color), 'Parent', ax1);
                    fill_image2 = imagesc(Vector2Colormap(-logo_1,handles.Right_color), 'Parent', ax2);
                    set(fill_image, 'Tag', 'fill_image_left', 'HitTest', 'on');
                    set(fill_image2, 'Tag', 'fill_image_right', 'HitTest', 'on');

                else

                    % Dummy data to put into the axes on startup
                    z=peaks(1000);
                    z = z./max(abs(z(:)));
                    fill_image = imshow(z, 'Parent', ax1, 'ColorMap', jet, 'DisplayRange', [min(z(:)) max(z(:))]);
                    set(fill_image, 'Tag', 'fill_image_left', 'HitTest', 'on');
                    freezeColors(ax1);

                end
            end
                
        else        
        
            if handles.N_channels == 1;

                % Pull slider value
                slide_frame = 1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')));


                    if handles.Autoscale_left == 0;
                        min_max_left = handles.Display_range_left;
                    else
                        min_max_left = handles.Min_max_left(slide_frame, :);
                    end
                    
                    OldXLimits = get(ax1, 'XLim');
                    OldYLimits = get(ax1, 'YLim');

                    % Set left axis to that frame
                    image(Vector2Colormap_setscale(handles.Img_stack(:,:,slide_frame,1), handles.Left_color, min_max_left), ...
                        'Parent', ax1, 'Tag', 'Left Image');
                    set(ax1, 'xtick', [], 'ytick', []);
                    axis(ax1, 'image');
                    set(ax1, 'XLim', OldXLimits, 'YLim', OldYLimits);
                    


            elseif handles.N_channels == 2;

                    % Pull slider value
                    slide_frame = 1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')));


                    % Set both axes to that frame


                    if handles.Autoscale_left == 0;
                        min_max_left = handles.Display_range_left;
                    else
                        min_max_left = handles.Min_max_right(slide_frame, :);
                    end

                    if handles.Autoscale_right == 0;
                        min_max_right = handles.Display_range_right;
                    else
                        min_max_right = handles.Min_max_right(slide_frame, :);
                    end
                    OldXLimits = get(ax1, 'XLim');
                    OldYLimits = get(ax1, 'YLim');
                    
                    image(Vector2Colormap_setscale(handles.Img_stack(:,:,slide_frame,1), handles.Left_color, min_max_left), ...
                        'Parent', ax1, 'Tag', 'Left Image');
                    set(ax1, 'xtick', [], 'ytick', []);
                    axis(ax1, 'image');
                    set(ax1, 'XLim', OldXLimits, 'YLim', OldYLimits);
                    
                    OldXLimits = get(ax2, 'XLim');
                    OldYLimits = get(ax2, 'YLim');
                    %                     disp(handles.Right_color)
                    image(Vector2Colormap_setscale(handles.Img_stack(:,:,slide_frame,2), handles.Right_color, min_max_right), ...
                        'Parent', ax2, 'Tag', 'Right Image');
                    set(ax2, 'xtick', [], 'ytick', []);
                    axis(ax2, 'image');
                    set(ax2, 'XLim', OldXLimits, 'YLim', OldYLimits);
            end
        end
        
        if handles.Left_invert == 1;
            
            axis_handle = get(findobj('Tag', 'axes_panel1'), 'Children');
            set(axis_handle, 'XDir', 'reverse');
            
        end
        
        if handles.Right_invert == 1;
            
            axis_handle = get(findobj('Tag', 'axes_panel2'), 'Children');
            set(axis_handle, 'XDir', 'reverse');
            
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%
% Set parameters for image display

    function Image_prefs(varargin)
        
        % Make sure there isn't another one of these already open.  If so,
        % bring it to the front.  
        
        if ~isempty(findobj('Tag', 'GALAH_Image_prefs'))
        
            uistack(findobj('Tag', 'GALAH_Image_prefs'), 'top');
            
        else
            
            %fig1 = findobj('Tag', 'TIFF viewer');
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);      
            fig2_size = [400 300];
            fig2_position = [(mf_post(1) + (mf_post(3) - fig2_size(1))/2) (mf_post(2) + (mf_post(4) - fig2_size(2))/2)];
            fig2 = figure('Name','Image Preferences', 'Tag', 'GALAH_Image_prefs', 'Units', 'pixels',...
                'Position',[fig2_position fig2_size], 'NumberTitle', 'off', 'Toolbar', 'none', 'Menu', 'none');
            set(fig2, 'Color',[0.9 0.9 0.9]);

            fig2_green = uipanel(fig2, 'Units', 'normalized', 'Position', [0 .45, 1, .44], ...
                'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'green_panel', 'Title', 'Channel 1');

            fig2_red = uipanel(fig2, 'Units', 'normalized', 'Position', [0 0, 1, .44], ...
                'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'red_panel', 'Title', 'Channel 2');

            fig2_top = uipanel(fig2, 'Units', 'normalized', 'Position', [0 .89, 1, .11], ...
                'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'none', 'Tag', 'top_panel');

            handles.handles.fig2_green = fig2_green;
            handles.handles.fig2_red = fig2_red;
            handles.handles.fig2_top = fig2_top;

            %%%%%%%%%%%%%%%%%%
            % Single/dual channel toggle

            dual_single_radio = uibuttongroup('visible', 'off', 'Parent', fig2_top, 'Units', 'normalized', ...
                'Position', [0 0 1 1], 'BorderType', 'none', 'BackgroundColor', [.9 .9 .9]);
            ds1 = uicontrol('Style', 'togglebutton', 'String', 'Single Channel', 'Parent', dual_single_radio, ...
                'Units', 'normalized', 'Position', [.05 .05 .4 .9]);
            ds2 = uicontrol('Style', 'togglebutton', 'String', 'Dual Channel', 'Parent', dual_single_radio, ...
                'Units', 'normalized', 'Position', [.55 .05 .4 .9]);
            set(dual_single_radio, 'SelectionChangeFcn', @dual_single_push);
            radio_handles = [ds1 ds2];
            set(dual_single_radio, 'SelectedObject', radio_handles(handles.N_channels));
            set(dual_single_radio, 'Visible', 'on');
            
            handles.handles.dual_single_radio.Single = ds1;
            handles.handles.dual_single_radio.Dual = ds2;
            handles.handles.dual_slingle_radio = dual_single_radio;

            %%%%%%%%%%%%%%%%%%
            % Channel 1 (green channel) sliders and such
            
             if isempty(handles.Load_file)
                 
                slider_step = 1;
                green_range = [0 1];
%                 green_max_slider_display = 1;
%                 green_min_slider_display = 0;
                slider_value_green_max = 1;
                slider_value_green_min = 0;
                
             elseif ~isempty(handles.Load_file)

                green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];

                slider_value_green_max = (handles.Display_range_left(2) - green_range(1))/(green_range(2) - green_range(1));
                slider_step = 1/((green_range(2)-green_range(1))-1);
                
                slider_value_green_min = (handles.Display_range_left(1) - green_range(1))/(green_range(2) - green_range(1));
                slider_step = 1/((green_range(2)-green_range(1))-1);
                
             end

            green_max_slide_hand = uicontrol(fig2_green, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [slider_step slider_step], 'Min', 0, 'Max', 1, 'Value', slider_value_green_max, 'Position', [.30 .77 .68 .1],...
                'Callback', @slider_green_max_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'Green max');

            addlistener(green_max_slide_hand, 'Value', 'PostSet', @slider_green_max_listener);

            green_max_slide_box = uicontrol(fig2_green, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.18 .71 .1 .25], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.Display_range_left(2)), 'Callback', @edit_green_max_call);

            green_max_slide_text = uicontrol(fig2_green, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.01 .75 .16 .14], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Display Max:');

            set(green_max_slide_hand, 'Enable', 'off');
            set(green_max_slide_box, 'Enable', 'off');

            green_min_slide_hand = uicontrol(fig2_green, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [slider_step slider_step], 'Min', 0, 'Max', 1, 'Value', slider_value_green_min, 'Position', [.3 .46 .68 .1],...
                'Callback', @slider_green_min_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'Green max');

            addlistener(green_min_slide_hand, 'Value', 'PostSet', @slider_green_min_listener);

            green_min_slide_box = uicontrol(fig2_green, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.18 .39 .1 .25], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.Display_range_left(1)), 'Callback', @edit_green_min_call);

            green_min_slide_text = uicontrol(fig2_green, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.01 .43 .16 .14], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Display Min:');

            Colormap_strings = {'Gray'; 'Jet'; 'Green'; 'Red'; 'Hot'; 'Cool'; 'Spring'; 'Summer'; 'Autumn'; 'Winter'};
            handles.Colormap_strings = Colormap_strings;
            left_value = find(strcmpi(handles.Left_color, Colormap_strings));

            green_colormap_listbox = uicontrol(fig2_green, 'Style', 'popupmenu', 'Units', 'normalized', ...
                'Position', [.18 .095 .22 .2], 'String', Colormap_strings, 'Value', left_value, 'Callback', @popup_green_colormap);

            green_colormap_text = uicontrol(fig2_green, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.01 .05 .16 .2], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Colormap:');

            green_autoscale = uicontrol('Style', 'checkbox', 'String', 'Autoscale', 'Parent', fig2_green, ...
                'Units', 'normalized', 'Position', [.50 .06 .2 .25], 'BackgroundColor', [.9 .9 .9], ...
                'Value', handles.Autoscale_left, 'Callback', @autoscale_green);

            green_invert = uicontrol('Style', 'checkbox', 'String', 'Invert Image', 'Parent', fig2_green, ...
                'Units', 'normalized', 'Position', [.76 .06 .2 .25], 'BackgroundColor', [.9 .9 .9], ...
                'Value', handles.Left_invert, 'Callback', @invert_green);

            if handles.Autoscale_left == 1
                set(green_max_slide_hand, 'Enable', 'off');
                set(green_max_slide_box, 'Enable', 'off');
                set(green_min_slide_hand, 'Enable', 'off');
                set(green_min_slide_box, 'Enable', 'off');
            else
                set(green_max_slide_hand, 'Enable', 'on');
                set(green_max_slide_box, 'Enable', 'on');
                set(green_min_slide_hand, 'Enable', 'on');
                set(green_min_slide_box, 'Enable', 'on');
            end
            set(green_colormap_listbox, 'Enable', 'on');

            %%%%%%%%%%%%%%%%%%
            % Channel 2 (red channel) sliders and such
            
            if isempty(handles.Load_file)
                red_range = [0 1];
                slider_step = 1;
                red_max_slider_display = 1;
                red_min_slider_display = 0;
                slider_value_red_max = 1;
                slider_value_red_min = 0;
                
            else

                if handles.N_channels == 2
                    red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
                    slider_step = 1/((red_range(2)-red_range(1))-1);
                    red_max_slider_display = handles.Display_range_right(2);
                    red_min_slider_display = handles.Display_range_right(1);
                    slider_value_red_max = (handles.Display_range_right(2) - red_range(1))/(red_range(2) - red_range(1));
                    slider_value_red_min = (handles.Display_range_right(1) - red_range(1))/(red_range(2) - red_range(1));
                else
                    slider_step = 1;
                    red_range = [0 1];
                    red_max_slider_display = 1;
                    red_min_slider_display = 0;
                    slider_value_red_max = 1;
                    slider_value_red_min = 0;
                end
                
            end

            red_max_slide_hand = uicontrol(fig2_red, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [slider_step slider_step], 'Min', 0, 'Max', 1, 'Value', slider_value_red_max, 'Position', [.30 .77 .68 .1],...
                'Callback', @slider_red_max_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'red max');

            addlistener(red_max_slide_hand, 'Value', 'PostSet', @slider_red_max_listener);

            red_max_slide_box = uicontrol(fig2_red, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.18 .71 .1 .25], 'BackgroundColor', [1 1 1], ...
                'String', num2str(red_max_slider_display), 'Callback', @edit_red_max_call);

            red_max_slide_text = uicontrol(fig2_red, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.01 .75 .16 .14], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Display Max:');

            set(red_max_slide_hand, 'Enable', 'off');
            set(red_max_slide_box, 'Enable', 'off');

            

            red_min_slide_hand = uicontrol(fig2_red, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [slider_step slider_step], 'Min', 0, 'Max', 1, 'Value', slider_value_red_min, 'Position', [.3 .46 .68 .1],...
                'Callback', @slider_red_min_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'red max');

            addlistener(red_min_slide_hand, 'Value', 'PostSet', @slider_red_min_listener);

            red_min_slide_box = uicontrol(fig2_red, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.18 .39 .1 .25], 'BackgroundColor', [1 1 1], ...
                'String', num2str(red_min_slider_display), 'Callback', @edit_red_min_call);

            red_min_slide_text = uicontrol(fig2_red, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.01 .43 .16 .14], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Display Min:');

            Colormap_strings = {'Gray'; 'Jet'; 'Green'; 'Red'; 'Hot'; 'Cool'; 'Spring'; 'Summer'; 'Autumn'; 'Winter'};
            right_value = find(strcmpi(handles.Right_color, Colormap_strings));

            red_colormap_listbox = uicontrol(fig2_red, 'Style', 'popupmenu', 'Units', 'normalized', ...
                'Position', [.18 .095 .22 .2], 'String', Colormap_strings, 'Value', right_value, 'Callback', @popup_red_colormap);

            red_colormap_text = uicontrol(fig2_red, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.01 .05 .16 .2], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Colormap:');

            red_autoscale = uicontrol('Style', 'checkbox', 'String', 'Autoscale', 'Parent', fig2_red, ...
                'Units', 'normalized', 'Position', [.50 .06 .2 .25], 'BackgroundColor', [.9 .9 .9], ...
                'Value', handles.Autoscale_right, 'Callback', @autoscale_red);

            red_invert = uicontrol('Style', 'checkbox', 'String', 'Invert Image', 'Parent', fig2_red, ...
                'Units', 'normalized', 'Position', [.76 .06 .2 .25], 'BackgroundColor', [.9 .9 .9], ...
                'Value', handles.Right_invert, 'Callback', @invert_red);
            


            if handles.Autoscale_left == 1;
                set(green_max_slide_hand, 'Enable', 'off');
                set(green_max_slide_box, 'Enable', 'off');
                set(green_min_slide_hand, 'Enable', 'off');
                set(green_min_slide_box, 'Enable', 'off');
                set(green_min_slide_text, 'Enable', 'off');
                set(green_max_slide_text, 'Enable', 'off');
            else
                set(green_max_slide_hand, 'Enable', 'on');
                set(green_max_slide_box, 'Enable', 'on');
                set(green_min_slide_hand, 'Enable', 'on');
                set(green_min_slide_box, 'Enable', 'on');
                set(green_min_slide_text, 'Enable', 'on');
                set(green_max_slide_text, 'Enable', 'on');
            end

            set(red_colormap_listbox, 'Enable', 'on');

            if handles.Autoscale_right == 1;
                set(red_max_slide_hand, 'Enable', 'off');
                set(red_max_slide_box, 'Enable', 'off');
                set(red_min_slide_hand, 'Enable', 'off');
                set(red_min_slide_box, 'Enable', 'off');
                set(red_min_slide_text, 'Enable', 'off');
                set(red_max_slide_text, 'Enable', 'off');
            else
                set(red_max_slide_hand, 'Enable', 'on');
                set(red_max_slide_box, 'Enable', 'on');
                set(red_min_slide_hand, 'Enable', 'on');
                set(red_min_slide_box, 'Enable', 'on');
                set(red_min_slide_text, 'Enable', 'on');
                set(red_max_slide_text, 'Enable', 'on');
            end
            
            if handles.N_channels == 1;
                
                set(red_max_slide_hand, 'Enable', 'off');
                set(red_max_slide_box, 'Enable', 'off', 'String', []);
                set(red_min_slide_hand, 'Enable', 'off');
                set(red_min_slide_box, 'Enable', 'off', 'String', []);
                set(red_min_slide_text, 'Enable', 'off');
                set(red_max_slide_text, 'Enable', 'off');
                set(red_autoscale, 'Enable', 'off');
                set(red_invert, 'Enable', 'off');
                set(red_colormap_listbox, 'Enable', 'off');
                set(red_colormap_text, 'Enable', 'off');
                
            end
            
            if isempty(handles.Load_file)
                
                set(green_max_slide_hand, 'Enable', 'off');
                set(green_max_slide_box, 'Enable', 'off', 'String', []);
                set(green_min_slide_hand, 'Enable', 'off');
                set(green_min_slide_box, 'Enable', 'off', 'String', []);
                set(green_min_slide_text, 'Enable', 'off');
                set(green_max_slide_text, 'Enable', 'off');
                set(green_autoscale, 'Enable', 'off');

                set(red_max_slide_hand, 'Enable', 'off');
                set(red_max_slide_box, 'Enable', 'off', 'String', []);
                set(red_min_slide_hand, 'Enable', 'off');
                set(red_min_slide_box, 'Enable', 'off', 'String', []);
                set(red_min_slide_text, 'Enable', 'off');
                set(red_max_slide_text, 'Enable', 'off');
                set(red_autoscale, 'Enable', 'off');
                
            end
            
            if mod(handles.N_frames*handles.N_channels,2) == 1
                set(handles.handles.dual_single_radio.Single, 'Enable', 'off')
                set(handles.handles.dual_single_radio.Dual, 'Enable', 'off')
            end
            
            handles.handles.green_max_slide_hand = green_max_slide_hand;
            handles.handles.green_max_slide_box = green_max_slide_box;
            handles.handles.green_min_slide_hand = green_min_slide_hand;
            handles.handles.green_min_slide_box = green_min_slide_box;
            handles.handles.green_colormap_text = green_colormap_text;
            handles.handles.green_colormap_listbox = green_colormap_listbox;
            handles.handles.green_autoscale = green_autoscale;
            handles.handles.green_invert = green_invert;
            
            handles.handles.red_max_slide_hand = red_max_slide_hand;
            handles.handles.red_max_slide_box = red_max_slide_box;
            handles.handles.red_min_slide_hand = red_min_slide_box;
            handles.handles.red_min_slide_box = red_max_slide_box;
            handles.handles.red_colormap_listbox = red_colormap_listbox;
            handles.handles.red_autoscale = red_autoscale;
            handles.handles.red_invert = red_invert;
            

            guidata(findobj('Tag', 'TIFF viewer'), handles);
        end
            %%%% Big pile of callback functions for display

        function dual_single_push(~, eventdata, ~)
            % Set as single or dual-channel data
            % Dual channel available only if data has even number of frames
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            channels_now = find(eventdata.NewValue == [handles.handles.dual_single_radio.Single handles.handles.dual_single_radio.Dual]);
            
            if handles.N_frames*handles.N_channels == 1;
                
                % If there is only one frame, it can only be a
                % single-channel data set.  This forces that fact. 
                % There shouldn't ever be anything to change as the
                % single-frame/single-channel issue is addressed upon
                % loading.
                
                channels_now = 3;
                set(handles.handles.dual_slingle_radio, 'SelectedObject', handles.handles.dual_single_radio.Single);
                
                
            end
            
            if channels_now == 1;
            	handles.N_channels = 1;
            
                % Disable all of right channel
                
                set(red_max_slide_hand, 'Enable', 'off');
                set(red_max_slide_box, 'Enable', 'off', 'String', []);
                set(red_min_slide_hand, 'Enable', 'off');
                set(red_min_slide_box, 'Enable', 'off', 'String', []);
                set(red_min_slide_text, 'Enable', 'off');
                set(red_max_slide_text, 'Enable', 'off');
                set(red_autoscale, 'Enable', 'off');
                set(red_invert, 'Enable', 'off');
                set(red_colormap_listbox, 'Enable', 'off');
                set(red_colormap_text, 'Enable', 'off');

                if ~isempty(handles.Load_file)
                    
                    % Collapse Img_stack, Min_max_XXX down to a single dimension

                    green_frames = 1:2:(handles.N_frames*2);
                    red_frames = 2:2:(handles.N_frames*2);

                    Img_hold = zeros(size(handles.Img_stack,1), size(handles.Img_stack,2), 2*handles.N_channels, 1);
                    Img_hold(:,:,green_frames) = handles.Img_stack(:,:,:,1);
                    Img_hold(:,:,red_frames) = handles.Img_stack(:,:,:,2);

                    Min_max_hold = zeros(handles.N_frames, 2);
                    Min_max_hold(green_frames, :) = handles.Min_max_left;
                    Min_max_hold(red_frames, :) = handles.Min_max_right;

                    handles.Img_stack = Img_hold;
                    handles.Min_max_left = Min_max_hold;
                    handles.Min_max_right = [];
                    clear Img_hold Min_max_hold;

                    handles.N_frames = size(handles.Img_stack, 3);

                    if handles.Primary_channel > handles.N_channels
                        handles.Primary_channel = handles.N_channels;
                    end


                    % Figure out where slider should be with new N_channels

                    slider_step = 1/(handles.N_frames-1);

                    if handles.N_frames == 1;
                        set(slide_hand, 'SliderStep', [1 1]);
                    else
                        set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                    end 
                    
                    set(slide_box, 'String', (1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')))));

                    green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
                    slider_step_green = 1/((green_range(2)-green_range(1))-1);


                    set(green_max_slide_hand, 'SliderStep', [slider_step_green slider_step_green]);
                    set(green_min_slide_hand, 'SliderStep', [slider_step_green slider_step_green]); 

                    slide_string_max = str2double(get(green_max_slide_box, 'String'));
                    slide_set_max = ((slide_string_max - green_range(1))/(green_range(2) - green_range(1)));
                    slide_set_max = min([slide_set_max 1]); 
                    slider_value_max = (green_range(1) + slide_set_max*(green_range(2) - green_range(1)));
                    set(green_max_slide_box, 'String', num2str(slider_value_max));
                    set(green_max_slide_hand, 'Value', slide_set_max);

                    slide_string_min = str2double(get(green_min_slide_box, 'String'));
                    slide_set_min = ((slide_string_min - green_range(1))/(green_range(2) - green_range(1)));
                    slide_set_min = max([slide_set_min 0]);
                    slider_value_min = (green_range(1) + slide_set_min*(green_range(2) - green_range(1)));
                    set(green_min_slide_box, 'String', num2str(slider_value_min));
                    set(green_min_slide_hand, 'Value', slide_set_min);


                    % Fill in red channel with dummy image
                    path_here = mfilename('fullpath');
                    logo_file = fullfile(fileparts(path_here), 'BMIF_logo.jpg');

                    %disp(logo_file);

                    ax2 = handles.handles.ax2;

                    if exist(logo_file, 'file') == 2;

                        logo_hold = single(imread(logo_file));
                        logo_2 = logo_hold(:,:,1);
                        clear logo_hold
                        %disp(size(logo_2));
                        fill_image = imagesc(Vector2Colormap(-logo_2,handles.Right_color), 'Parent', ax2);
                        set(fill_image, 'Tag', 'fill_image_right', 'HitTest', 'on');

                    else

                        % Dummy data to put into the axes on startup
                        z=peaks(1000);
                        z = z./max(abs(z(:)));
                        fill_image = imshow(z, 'Parent', ax2, 'ColorMap', jet, 'DisplayRange', [min(z(:)) max(z(:))]);
                        set(fill_image, 'Tag', 'fill_image_right', 'HitTest', 'on');
                        freezeColors(ax2);

                    end

                    % Get rid of tick labels
                    set(ax2, 'xtick', [], 'ytick', []);

                    guidata(findobj('Tag', 'TIFF viewer'), handles);

                    Display_images_in_axes;
                    
                else
                    
                    guidata(findobj('Tag', 'TIFF viewer'), handles);
                
                end
                
            elseif channels_now == 2;
                handles.N_channels = 2;
                
                if ~isempty(handles.Load_file)
                
                % Enable right channel
                
                set(red_max_slide_hand, 'Enable', 'on');
                set(red_max_slide_box, 'Enable', 'on');
                set(red_min_slide_hand, 'Enable', 'on');
                set(red_min_slide_box, 'Enable', 'on');
                set(red_min_slide_text, 'Enable', 'on');
                set(red_max_slide_text, 'Enable', 'on');
                set(red_autoscale, 'Enable', 'on');
                set(red_invert, 'Enable', 'on');
                set(red_colormap_listbox, 'Enable', 'on');
                set(red_colormap_text, 'Enable', 'on');
                
                set(handles.handles.Align_out, 'Enable', 'on');
                
                % Expand Img_stack to two channels
        
                    green_frames = 1:2:(handles.N_frames);
                    red_frames = 2:2:(handles.N_frames);

                    Img_hold = zeros(size(handles.Img_stack,1), size(handles.Img_stack,2), handles.N_frames/2, 2);
                    Img_hold(:,:,:,1) = handles.Img_stack(:,:,green_frames);
                    Img_hold(:,:,:,2) = handles.Img_stack(:,:,red_frames);

%                     Min_max_hold_left = zeros(handles.N_frames, 2);
%                     Min_max_hold_right = zeros(handles.N_frames, 2);
                    Min_max_hold_left = handles.Min_max_left(green_frames, :);
                    Min_max_hold_right = handles.Min_max_left(red_frames, :);

                    handles.Min_max_left = Min_max_hold_left;
                    handles.Min_max_right = Min_max_hold_right;
                    handles.Img_stack = Img_hold;
                    clear Img_hold Min_max_hold_left Min_max_hold_right

                    handles.N_frames = size(handles.Img_stack, 3);

                    if handles.Primary_channel > handles.N_channels
                        handles.Primary_channel = handles.N_channels;
                    end

                    % Figure out where sliders should be with new N_channels

                    slider_step = 1/(handles.N_frames-1);
                    slide_hand = handles.handles.slide_hand;
                    slide_box = handles.handles.slide_box;
                    
                    if handles.N_frames == 1;
                        set(slide_hand, 'SliderStep', [1 1]);
                    
                    else
                        set(slide_hand, 'SliderStep', [1/(handles.N_frames-1) 1/(handles.N_frames-1)]);
                    end
                    
                    set(slide_box, 'String', (1 + round((handles.N_frames - 1)*(get(slide_hand, 'Value')))));

                    green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
                    slider_step_green = 1/((green_range(2)-green_range(1))-1);

                    set(green_max_slide_hand, 'SliderStep', [slider_step_green slider_step_green]);
                    set(green_min_slide_hand, 'SliderStep', [slider_step_green slider_step_green]);

                    slide_string_max = str2double(get(green_max_slide_box, 'String'));
                    slide_set_max = ((slide_string_max - green_range(1))/(green_range(2) - green_range(1)));
                    slide_set_max = min([slide_set_max 1]); 
                    slider_value_max = (green_range(1) + slide_set_max*(green_range(2) - green_range(1)));
                    set(green_max_slide_box, 'String', num2str(slider_value_max));
                    set(green_max_slide_hand, 'Value', slide_set_max);

                    slide_string_min = str2double(get(green_min_slide_box, 'String'));
                    slide_set_min = ((slide_string_min - green_range(1))/(green_range(2) - green_range(1)));
                    slide_set_min = max([slide_set_min 0]); 
                    slider_value_min = (green_range(1) + slide_set_min*(green_range(2) - green_range(1)));
                    set(green_min_slide_box, 'String', num2str(slider_value_min));
                    set(green_min_slide_hand, 'Value', slide_set_min);


                    red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
                    slider_step_red = 1/((red_range(2)-red_range(1))-1);
                    set(red_max_slide_hand, 'SliderStep', [slider_step_red slider_step_red]);
                    set(red_min_slide_hand, 'SliderStep', [slider_step_red slider_step_red]);
                    set(red_max_slide_box, 'String', num2str(handles.Display_range_right(2)));
                    set(red_min_slide_box, 'String', num2str(handles.Display_range_right(1)));

                    % Replot channels
                    
                    NewXLim = [0.5 size(handles.Img_stack, 2)+0.5];
                    NewYLim = [0.5 size(handles.Img_stack, 1)+0.5];
                    set(handles.handles.ax2, 'XLim', NewXLim, 'YLim', NewYLim);


                    guidata(findobj('Tag', 'TIFF viewer'), handles);
                    Display_images_in_axes;
                    
                else
                    
                    set(red_invert, 'Enable', 'on');
                    set(red_colormap_listbox, 'Enable', 'on');
                    set(red_colormap_text, 'Enable', 'on');
                    
                   guidata(findobj('Tag', 'TIFF viewer'), handles); 
                
                end
            
            end
            
                    if handles.N_frames == 1
                        set(handles.handles.Unbind_out, 'Enable', 'off');
                        set(handles.handles.ExpFit_out, 'Enable', 'off');
                        set(slide_hand, 'Enable', 'off')
                        set(slide_box, 'Enable', 'off')
                    else
                        set(handles.handles.Unbind_out, 'Enable', 'on');
                        set(handles.handles.ExpFit_out, 'Enable', 'on');
                        set(slide_hand, 'Enable', 'on')
                        set(slide_box, 'Enable', 'on')
                        
                    end

        end

        function slider_green_max_call(varargin)
            % Max intenstiy green channel
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_green_here = get(green_max_slide_hand, 'Value');
            slider_check_here = get(green_min_slide_hand, 'Value');
            slider_step = get(green_max_slide_hand, 'SliderStep');
            
            if le(slider_green_here, slider_check_here)
                %disp('slider_check');
                slider_green_here = slider_check_here + slider_step(1);
                set(green_max_slide_hand, 'Value', slider_green_here);
            end    
            
            green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
            slider_value = round(slider_green_here*(green_range(2) - green_range(1)) + green_range(1));
            
            set(green_max_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_left(2) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
        end

        function slider_green_max_listener(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_green_here = get(green_max_slide_hand, 'Value');
            slider_check_here = get(green_min_slide_hand, 'Value');
            slider_step = get(green_max_slide_hand, 'SliderStep');
            
            if le(slider_green_here, slider_check_here)
                slider_green_here = slider_check_here + slider_step(1);
                set(green_max_slide_hand, 'Value', slider_green_here);
            end 
            
            green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
            slider_value = round(slider_green_here*(green_range(2) - green_range(1)) + green_range(1));
            
            set(green_max_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_left(2) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function edit_green_max_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slide_string = str2double(get(green_max_slide_box, 'String'));
            green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
            
            if length(slide_string) ~= 1
                slide_set = get(green_max_slide_hand, 'Value');
                slide_str2 = round(green_range(2)+slide_set*(green_range(2) - green_range(1)));
            
            else
        
            slide_set = ((slide_string - green_range(1))/(green_range(2) - green_range(1)));
            slide_range = [get(green_max_slide_hand, 'Min') get(green_max_slide_hand, 'Max')];

                if slide_set > slide_range(2)

                    slide_set = slide_range(2);
                    slide_str2 = (green_range(1) + slide_set*(green_range(2) - green_range(1)));

                elseif slide_set < slide_range(1)

                    slide_set = slide_range(1);
                    slide_str2 = (green_range(1) + slide_set*(green_range(2) - green_range(1)));
                    
                else 
                    
                    slide_str2 = (green_range(1) + slide_set*(green_range(2) - green_range(1)));

                end
        
            end
            
            
            slider_value = slide_str2;
            
            set(green_max_slide_box, 'String', num2str(slider_value));
            set(green_max_slide_hand, 'Value', slide_set);
            
            handles.Display_range_left(2) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust')

        end

        function slider_green_min_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_green_here = get(green_min_slide_hand, 'Value');
            slider_check_here = get(green_max_slide_hand, 'Value');
            slider_step = get(green_min_slide_hand, 'SliderStep');
            
            if ge(slider_green_here, slider_check_here)
                %disp('slider_check');
                slider_green_here = slider_check_here - slider_step(1);
                set(green_min_slide_hand, 'Value', slider_green_here);
            end 
            

            green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
            slider_value = round(slider_green_here*(green_range(2) - green_range(1)) + green_range(1));
            
            set(green_min_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_left(1) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            

        end

        function slider_green_min_listener(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_green_here = get(green_min_slide_hand, 'Value');
            slider_check_here = get(green_max_slide_hand, 'Value');
            slider_step = get(green_min_slide_hand, 'SliderStep');
            
            if ge(slider_green_here, slider_check_here)
                %disp('slider_check');
                slider_green_here = slider_check_here - slider_step(1);
                set(green_min_slide_hand, 'Value', slider_green_here);
            end 
            
            green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
            slider_value = round(slider_green_here*(green_range(2) - green_range(1)) + green_range(1));
            
            set(green_min_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_left(1) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust')

        end

        function edit_green_min_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slide_string = str2double(get(green_min_slide_box, 'String'));
            green_range = [min(handles.Min_max_left(:,1)), max(handles.Min_max_left(:,2))];
            
            if length(slide_string) ~= 1
                slide_set = get(green_max_slide_hand, 'Value');
                slide_str2 = round(green_range(2)+slide_set*(green_range(2) - green_range(1)));
            
            else
        
            slide_set = ((slide_string - green_range(1))/(green_range(2) - green_range(1)));
            slide_range = [get(green_max_slide_hand, 'Min') get(green_max_slide_hand, 'Max')];

                if slide_set > slide_range(2)

                    slide_set = slide_range(2);
                    slide_str2 = (green_range(1) + slide_set*(green_range(2) - green_range(1)));

                elseif slide_set < slide_range(1)

                    slide_set = slide_range(1);
                    slide_str2 = (green_range(1) + slide_set*(green_range(2) - green_range(1)));
                    
                else 
                    
                    slide_str2 = (green_range(1) + slide_set*(green_range(2) - green_range(1)));

                end
        
            end
            
            
            slider_value = slide_str2;
            
            set(green_min_slide_box, 'String', num2str(slider_value));
            set(green_min_slide_hand, 'Value', slide_set);
            
            handles.Display_range_left(1) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function popup_green_colormap(varargin) %%%% Change colormap in green channel
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            string_here = get(green_colormap_listbox, 'Value');
            handles.Left_color = lower(handles.Colormap_strings{string_here});

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;

        end

        function autoscale_green(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Autoscale_left = get(green_autoscale, 'Value');
            
            if handles.Autoscale_left == 1;
                set(green_max_slide_hand, 'Enable', 'off');
                set(green_max_slide_box, 'Enable', 'off');
                set(green_min_slide_hand, 'Enable', 'off');
                set(green_min_slide_box, 'Enable', 'off');
                set(green_min_slide_text, 'Enable', 'off');
                set(green_max_slide_text, 'Enable', 'off');
            else
                set(green_max_slide_hand, 'Enable', 'on');
                set(green_max_slide_box, 'Enable', 'on');
                set(green_min_slide_hand, 'Enable', 'on');
                set(green_min_slide_box, 'Enable', 'on');
                set(green_min_slide_text, 'Enable', 'on');
                set(green_max_slide_text, 'Enable', 'on');
            end

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function invert_green(varargin)
            % Set value (1 or 0) if green channel should be displayed as LR flip
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Left_invert = get(green_invert, 'Value');

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function slider_red_max_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_right_here = get(red_max_slide_hand, 'Value');
            slider_check_here = get(red_min_slide_hand, 'Value');
            slider_step = get(red_max_slide_hand, 'SliderStep');
            
            if le(slider_right_here, slider_check_here)
                %disp('slider_check');
                slider_right_here = slider_check_here + slider_step(1);
                set(red_max_slide_hand, 'Value', slider_right_here);
            end 

            red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
            slider_value = round(slider_right_here*(red_range(2) - red_range(1)) + red_range(1));
            
            set(red_max_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_right(2) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;

        end

        function slider_red_max_listener(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_right_here = get(red_max_slide_hand, 'Value');
            slider_check_here = get(red_min_slide_hand, 'Value');
            slider_step = get(red_max_slide_hand, 'SliderStep');
            
            if le(slider_right_here, slider_check_here)
                %disp('slider_check');
                slider_right_here = slider_check_here + slider_step(1);
                set(red_max_slide_hand, 'Value', slider_right_here);
            end 
            
            
            red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
            slider_value = round(slider_right_here*(red_range(2) - red_range(1)) + red_range(1));
            
            set(red_max_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_right(2) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function edit_red_max_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slide_string = str2double(get(red_max_slide_box, 'String'));
            red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
            
            if length(slide_string) ~= 1
                slide_set = get(red_max_slide_hand, 'Value');
                slide_str2 = round(red_range(2)+slide_set*(red_range(2) - red_range(1)));
            
            else
        
            slide_set = ((slide_string - red_range(1))/(red_range(2) - red_range(1)));
            slide_range = [get(red_max_slide_hand, 'Min') get(red_max_slide_hand, 'Max')];

                if slide_set > slide_range(2)

                    slide_set = slide_range(2);
                    slide_str2 = (red_range(1) + slide_set*(red_range(2) - red_range(1)));

                elseif slide_set < slide_range(1)

                    slide_set = slide_range(1);
                    slide_str2 = (red_range(1) + slide_set*(red_range(2) - red_range(1)));
                    
                else 
                    
                    slide_str2 = (red_range(1) + slide_set*(red_range(2) - red_range(1)));

                end
        
            end
            
            
            slider_value = slide_str2;
            
            set(red_max_slide_box, 'String', num2str(slider_value));
            set(red_max_slide_hand, 'Value', slide_set);
            
            handles.Display_range_right(2) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');
            


        end

        function slider_red_min_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_right_here = get(red_min_slide_hand, 'Value');
            slider_check_here = get(red_max_slide_hand, 'Value');
            slider_step = get(red_min_slide_hand, 'SliderStep');
            
            if ge(slider_right_here, slider_check_here)
                %disp('slider_check');
                slider_right_here = slider_check_here - slider_step(1);
                set(red_min_slide_hand, 'Value', slider_right_here);
            end 
            
            red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
            slider_value = round(slider_right_here*(red_range(2) - red_range(1)) + red_range(1));
            
            set(red_min_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_right(1) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;

        end

        function slider_red_min_listener(varargin)
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_right_here = get(red_min_slide_hand, 'Value');
            slider_check_here = get(red_max_slide_hand, 'Value');
            slider_step = get(red_min_slide_hand, 'SliderStep');
            
            if ge(slider_right_here, slider_check_here)
                %disp('slider_check');
                slider_right_here = slider_check_here - slider_step(1);
                set(red_min_slide_hand, 'Value', slider_right_here);
            end 

            red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
            slider_value = round(slider_right_here*(red_range(2) - red_range(1)) + red_range(1));
            
            set(red_min_slide_box, 'String', num2str(slider_value));
            
            handles.Display_range_right(1) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function edit_red_min_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slide_string = str2double(get(red_min_slide_box, 'String'));
            red_range = [min(handles.Min_max_right(:,1)), max(handles.Min_max_right(:,2))];
            
            if length(slide_string) ~= 1
                slide_set = get(red_min_slide_hand, 'Value');
                slide_str2 = round(red_range(2)+slide_set*(red_range(2) - red_range(1)));
            
            else
        
            slide_set = ((slide_string - red_range(1))/(red_range(2) - red_range(1)));
            slide_range = [get(red_min_slide_hand, 'Min') get(red_min_slide_hand, 'Max')];

                if slide_set > slide_range(2)

                    slide_set = slide_range(2);
                    slide_str2 = (red_range(1) + slide_set*(red_range(2) - red_range(1)));

                elseif slide_set < slide_range(1)

                    slide_set = slide_range(1);
                    slide_str2 = (red_range(1) + slide_set*(red_range(2) - red_range(1)));
                    
                else 
                    
                    slide_str2 = (red_range(1) + slide_set*(red_range(2) - red_range(1)));

                end
        
            end
            
            
            slider_value = slide_str2;
            
            set(red_min_slide_box, 'String', num2str(slider_value));
            set(red_min_slide_hand, 'Value', slide_set);
            
            handles.Display_range_right(1) = slider_value;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end

        function popup_red_colormap(varargin) % Change colormap in red channel
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            string_here = get(red_colormap_listbox, 'Value');
            handles.Right_color = lower(handles.Colormap_strings{string_here});

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;

        end

        function autoscale_red(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Autoscale_right = get(red_autoscale, 'Value');

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            if handles.Autoscale_right == 1;
                set(red_max_slide_hand, 'Enable', 'off');
                set(red_max_slide_box, 'Enable', 'off');
                set(red_min_slide_hand, 'Enable', 'off');
                set(red_min_slide_box, 'Enable', 'off');
                set(red_min_slide_text, 'Enable', 'off');
                set(red_max_slide_text, 'Enable', 'off');
            else
                set(red_max_slide_hand, 'Enable', 'on');
                set(red_max_slide_box, 'Enable', 'on');
                set(red_min_slide_hand, 'Enable', 'on');
                set(red_min_slide_box, 'Enable', 'on');
                set(red_min_slide_text, 'Enable', 'on');
                set(red_max_slide_text, 'Enable', 'on');
            end
            
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end
        
        function invert_red(varargin)
            
            % Set value (1 or 0) if red channel should be shown with flip
            % LR
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Right_invert = get(red_invert, 'Value');

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            Display_ROI_img('adjust');

        end
        
        
        
    end

%%%%%%%%%%%%%%%%%%%%%%
% Launch Background Trace window and utility

    function ROI_launch(varargin)
       
        if ~isempty(findobj('Tag', 'GALAH_ROI_calc'))
        
            uistack(findobj('Tag', 'GALAH_ROI_calc'), 'top');
            
        else
            
            fig1 = findobj('Tag', 'TIFF viewer');
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            fig1_post = get(fig1, 'Position');
            fig1_outer = get(fig1, 'OuterPos');
            up_pad = (fig1_outer(4) - fig1_post(4))*handles.scrsz_pixels(4);
            right_pad = (fig1_outer(3) - fig1_post(3))*handles.scrsz_pixels(3);
            
            mf_post = fig1_post.*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);  
            mf_post(1) = mf_post(1) - fig1_post(3)*handles.scrsz_pixels(3)*.2;
            %fig3_size = [560 560+190];
            
            scrsz = get(0, 'ScreenSize');
            Window_size_big = [560 560+190]; % Desired window size in pixels
            Window_size_small = [477 637];
            if scrsz(3) < 1300 || scrsz(4) < 800
                fig3_size = Window_size_small;
            else
                fig3_size = Window_size_big;
            end
            
            
            fig3_position = floor([(mf_post(1) + (mf_post(3) - fig3_size(1))/2) (mf_post(2) + (mf_post(4) - fig3_size(2))/2)]);
            
            fig3_loc = [fig3_size(1)+fig3_position(1)+right_pad, fig3_size(2)+fig3_position(2)+up_pad];
            
            if fig3_loc(1) > handles.scrsz_pixels(3)
                fig3_position(1) = fig3_position(1) - (fig3_loc(1) - handles.scrsz_pixels(3) - 1);
            end
            
            if fig3_loc(2) > handles.scrsz_pixels(4)
                fig3_position(2) = fig3_position(2) - (fig3_loc(2) - handles.scrsz_pixels(4) - 1);
            end
            
            fig3 = figure('Name','ROI Generator', 'Tag', 'GALAH_ROI_calc', 'Units', 'pixels',...
                'Position',[fig3_position fig3_size], 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure');
            set(fig3, 'Color',[0.9 0.9 0.9]);
            axis image;
            set(fig3, 'DeleteFcn', @CloseROIFcn);
            
            %%%%%%%%%%%%
            %Set up toolbar
%             hToolbar = findall(fig3,'tag','FigureToolBar');
%             AllToolHandles = findall(hToolbar);
%             ToolBarTags = get(AllToolHandles,'Tag');
%             ToolsToKeep = {'FigureToolBar'; 'Exploration.DataCursor'; 'Exploration.Pan'; 'Exploration.ZoomOut'; 'Exploration.ZoomIn'};
%             WhichTools = ~ismember(ToolBarTags, ToolsToKeep);
%             delete(AllToolHandles(WhichTools));
%             

            % Temporary menubar for dev purposes
            set(fig3, 'MenuBar', 'figure');
            
            %%%%%%%%%%%%%%%%%%%
            % Initialize panels
            

            fig3_img_panel = uipanel(fig3, 'Units', 'normalized', 'Position', [0, 1-(fig3_size(1)/fig3_size(2)), 1, (fig3_size(1)/fig3_size(2))], ...
                'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'ROI_img_panel');

            fig3_button_panel = uipanel(fig3, 'Units', 'normalized', 'Position', [0, 0, 1, 1-(fig3_size(1)/fig3_size(2))], ...
                'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'ROI_button_panel');
            
            %%%%%%%%%%%%%%%%%%%
            % Initialize axes
             OldXLimits = [0.5 size(handles.Img_stack, 2)+0.5];
             OldYLimits = [0.5 size(handles.Img_stack, 1)+0.5];
            ax_ROI = axes('Parent', fig3_img_panel, 'Position', [0 0 1 1]);
            set(ax_ROI, 'Tag', 'Axis_ROI');
            axis(ax_ROI, 'image');
            set(ax_ROI, 'xtick', [], 'ytick', []);
            set(ax_ROI, 'XLim', OldXLimits, 'YLim', OldYLimits);
            %axis(ax_ROI, 'image');
             %opengl software
            
            %%%%%%%%%%%%%%%%%%%
            % UI Components
            
            ROI_frame_slider_value = get(handles.handles.slide_hand, 'Value');
            ROI_frame_display = get(handles.handles.slide_box, 'String');
            
            if handles.N_frames == 1;
                ROI_frame_step = 1;
            else
                ROI_frame_step = 1/(handles.N_frames - 1);
            end
            
            ROI_frame_slide_hand = uicontrol(fig3_button_panel, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [ROI_frame_step ROI_frame_step], 'Min', 0, 'Max', 1, 'Value', ROI_frame_slider_value, ...
                'Position', [.36 .8146 .452 .07],...
                'Callback', @ROI_frame_slide_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'ROI_frame');

            addlistener(ROI_frame_slide_hand, 'Value', 'PostSet', @ROI_frame_listener_call);

            ROI_frame_slide_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.2646 .795 .08 .12], 'BackgroundColor', [1 1 1], ...
                'String', ROI_frame_display, 'Callback', @ROI_slide_edit_call);

            ROI_frame_slide_text = uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.182 .825 .08 .07], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Frame:');
            
            %%%%
            
            
            ROI_frame_UseAverage = uicontrol(fig3_button_panel, 'Style', 'checkbox', 'Units', 'normalized', ...
                'Position', [.82 .884 .176 .08], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Use Stack Mean', 'Value', 0, 'Callback', @ROI_frame_UseStackMean);
            
            ROI_frame_UseMaxProjection = uicontrol(fig3_button_panel, 'Style', 'checkbox', 'Units', 'normalized', ...
                'Position', [.82 .755 .176 .08], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Use Max Proj''n', 'Value', 0, 'Callback', @ROI_frame_UseMaxProjCallback);
            
            
            %%%%
            
            if handles.Primary_channel == 1;
                min_max_here = handles.Min_max_left;
%                 display_here = handles.Display_range_left;
            elseif handles.Primary_channel == 2;
                min_max_here = handles.Min_max_right;
%                 display_here = handles.Display_range_right;
            end
            
            ROI_thresh_step = 1/(max(min_max_here(:,2) - min(min_max_here(:,1))));
            ROI_thresh_big_step = min([1e-1 100*ROI_thresh_step]); 
            
            frame_here = str2double(ROI_frame_display);
            scale_frame = (handles.Img_stack(:, :, frame_here, handles.Primary_channel) - min_max_here(frame_here,1))./ ...
                (min_max_here(frame_here, 2) - min_max_here(frame_here, 1));
            
            ROI_thresh_slider_value = graythresh(scale_frame);
            ROI_thresh_display = num2str(round(ROI_thresh_slider_value*(max(min_max_here(:,2) - min(min_max_here(:,1)))) + min(min_max_here(:,1))));
            
            ROI_thresh_slide_hand = uicontrol(fig3_button_panel, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [ROI_thresh_step ROI_thresh_big_step], 'Min', 0, 'Max', 1, 'Value', ROI_thresh_slider_value, ...
                'Position', [.42 .59 .57 .07],...
                'Callback', @ROI_thresh_slide_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'ROI_thresh');

            addlistener(ROI_thresh_slide_hand, 'Value', 'PostSet', @ROI_thresh_listener_call);

            ROI_thresh_slide_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.31 .57 .09 .12], 'BackgroundColor', [1 1 1], ...
                'String', ROI_thresh_display, 'Callback', @ROI_thresh_edit_call);

            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.195 .59 .11 .08], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Threshold:');
            
           %%%%%
            
            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.0074 .775 .1 .16], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Primary Channel:');

            Primary_strings = cellstr(num2str((1:size(handles.Img_stack, 4))'));
                        
            Active_pulldown = uicontrol(fig3_button_panel, 'Style', 'popupmenu', 'Units', 'normalized', ...
                'Position', [.1156 .78 .06 .14], 'String', Primary_strings, 'Value', handles.Primary_channel, 'Callback', @Active_ROI_call);
            
            %%%%%
            
            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.0074 .55 .095 .16], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Dilate Radius:');
            
            Dilate_rad_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.1156 .57 .06 .12], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.Dilate_radius), 'Callback', @Dilate_edit_call);

            %%%%%
            
            
            End_caps_hand = uicontrol('Style', 'checkbox', 'String', 'Extend Ends', 'Parent', fig3_button_panel, ...
                'Units', 'normalized', 'Position', [.798 .33 .155 .09], 'BackgroundColor', [.9 .9 .9], ...
                'Value', handles.End_caps, 'Callback', @End_caps_callback);
            
            %%%%%
            
            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.005 .42 .2 .09], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'ROI Length (pixels):');
            
            ROI_length_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.21 .405 .09 .12], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.ROI_length), 'Callback', @ROI_length_edit_call);
            
            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.009 .27 .2 .09], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'ROI Width (pixels):');
            
            ROI_width_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.21 .255 .09 .12], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.ROI_width), 'Callback', @ROI_width_edit_call);
            
            
            %%%%
            
            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.305 .42 .225 .09], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Smoothing span (pixels):');
            
            ROI_smoothingSpan_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.53 .405 .09 .12], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.ROI_smoothing.span), 'Callback', @ROI_smoothingSpan_edit_call);
            
            uicontrol(fig3_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.3444 .27 .2 .09], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'Smoothing degree:');
            
            ROI_smoothingDegree_box = uicontrol(fig3_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.53 .255 .09 .12], 'BackgroundColor', [1 1 1], ...
                'String', num2str(handles.ROI_smoothing.degree), 'Callback', @ROI_smoothingDegree_edit_call);
            
            
            %%%%%
            
            ROI_choose_ROI = uicontrol(fig3_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Add ROI',...
                'Position', [.01 .04 .22 .175],...
                'Callback', @ROI_ROI, 'Tag', 'ROI_choose_ROI_button');

            %%%%%
            
            ROI_apply_hand = uicontrol(fig3_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Apply Settings',...
                'Position', [.27 .04 .22 .175],...
                'Callback', @ROI_apply, 'Tag', 'ROI_apply_botton');
            
           %%%%%
            
%             ROI_review_hand = uicontrol(fig3_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Review ROIs',...
%                 'Position', [.37 .04 .26 .175],...
%                 'Callback', @ROI_review, 'Tag', 'ROI_review_button');
            
            ROI_Plots_launch_hand = uicontrol(fig3_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Display Plots',...
                'Position', [.52 .04 .22 .175],...
                'Callback', @Plots_launch, 'Tag', 'Plots_launch_from_');
            
            %%%%%
            
            ROI_close_hand = uicontrol(fig3_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Set and Close',...
                'Position', [.77 .04 .22 .175],...
                'Callback', @ROI_close, 'Tag', 'ROI_close_button');
            
            %%%%%%%%%%%%%%%%%%%
            % Collect handles
            
            handles.handles.fig3 = fig3;
            handles.handles.fig3_img_panel = fig3_img_panel;
            handles.handles.fig3_button_panel = fig3_button_panel;
            handles.handles.ax_ROI = ax_ROI;
            handles.handles.ROI_frame_slide_hand = ROI_frame_slide_hand;
            handles.handles.ROI_frame_slide_box = ROI_frame_slide_box;
            handles.handles.ROI_thresh_slide_hand = ROI_thresh_slide_hand;
            handles.handles.ROI_thresh_slide_box = ROI_thresh_slide_box;
            handles.handles.Active_pulldown = Active_pulldown;
            handles.handles.Dilate_rad_box = Dilate_rad_box;
            handles.handles.End_caps_hand = End_caps_hand;
            handles.handles.ROI_length_box = ROI_length_box;
            handles.handles.ROI_width_box = ROI_width_box;
            handles.handles.ROI_apply_hand = ROI_apply_hand;
            %handles.handles.ROI_review_hand = ROI_review_hand;
            handles.handles.ROI_Plots_launch_hand = ROI_Plots_launch_hand;
            handles.handles.ROI_close_hand = ROI_close_hand;
            handles.handles.ROI_choose_ROI = ROI_choose_ROI;
            handles.handles.ROI_frame_UseAverage = ROI_frame_UseAverage;
            handles.handles.ROI_smoothingSpan_box = ROI_smoothingSpan_box;
            handles.handles.ROI_smoothingDegree_box = ROI_smoothingDegree_box;
            handles.handles.ROI_frame_UseMaxProjection = ROI_frame_UseMaxProjection;
            
            if handles.N_frames == 1;
                
                set(ROI_frame_slide_hand, 'Enable', 'off');
                set(ROI_frame_slide_box, 'Enable', 'off');
                set(ROI_frame_slide_text, 'Enable', 'off');
                
            end
            
            if handles.ROI_smoothing.span == 0
                set(handles.handles.ROI_smoothingDegree_box, 'Enable', 'off');
            end
            
            
            % Make ROI reviewer disabled until it is fixed.
            %set(ROI_review_hand, 'Enable', 'off');
            
            handles.ROIsAreSet = false;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            Display_ROI_img('initialize');
            Threshold_ROI_display;
            

        end
        
            %%%%%%%%%%%%%%%%%%%
            % UI Callbacks
            
        function ROI_frame_slide_call(varargin)
            
            % Listener handles all functions that would go here.

        end
            
        function ROI_frame_listener_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            slider_value = get(handles.handles.ROI_frame_slide_hand, 'Value');
            frame_num = 1 + slider_value*(handles.N_frames - 1);
            
            set(handles.handles.ROI_frame_slide_box, 'String', num2str(round(frame_num)));
            
            Display_ROI_img('framechange');
            if ~handles.ROIsAreSet
                Threshold_ROI_display;
            end
            
        end
            
        function ROI_slide_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            box_value = round(str2double(get(handles.handles.ROI_frame_slide_box, 'String')));
            
            %disp(box_value);
            
            if box_value < 1
                
                box_value = 1;
                
            elseif box_value > handles.N_frames
                
                box_value = handles.N_frames;
                
            end
            
            frame_num = (box_value - 1)./(handles.N_frames - 1);
            
            
            set(handles.handles.ROI_frame_slide_box, 'String', num2str(box_value));
            set(handles.handles.ROI_frame_slide_hand, 'Value', frame_num);
            
            Display_ROI_img('framechange');
   
            if ~handles.ROIsAreSet
                Threshold_ROI_display;
            end
            
        end
        
        function ROI_frame_UseStackMean(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            if get(varargin{1}, 'Value') == 1
                % Use stack mean active
                set(handles.handles.ROI_frame_slide_box, 'Enable', 'off');
                set(handles.handles.ROI_frame_slide_hand, 'Enable', 'off');
                
                % Unselect the MaxImg option
                set(handles.handles.ROI_frame_UseMaxProjection, 'Value', 0);
                
                Display_ROI_img('MeanImg');
            else
                % Use individual frames
               
                if handles.N_frames > 1
                    set(handles.handles.ROI_frame_slide_box, 'Enable', 'on');
                    set(handles.handles.ROI_frame_slide_hand, 'Enable', 'on');
                end
                
                Display_ROI_img('framechange');
            end
            
            
            if ~handles.ROIsAreSet
                Threshold_ROI_display;
            end
            
        end
        
        function ROI_frame_UseMaxProjCallback(varargin)
            % Use maximum intensity projection instead of a single frame
            % for segmentation
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            if get(varargin{1}, 'Value') == 1
                % Use stack mean active
                set(handles.handles.ROI_frame_slide_box, 'Enable', 'off');
                set(handles.handles.ROI_frame_slide_hand, 'Enable', 'off');
                
                % Unselect the MeanImg option
                set(handles.handles.ROI_frame_UseAverage, 'Value', 0);
                
                Display_ROI_img('MaxImg');
            else
                % Use individual frames
               
                if handles.N_frames > 1
                    set(handles.handles.ROI_frame_slide_box, 'Enable', 'on');
                    set(handles.handles.ROI_frame_slide_hand, 'Enable', 'on');
                end
                
                Display_ROI_img('framechange');
            end
            
            
            if ~handles.ROIsAreSet
                Threshold_ROI_display;
            end
            
        end
        
        
        
        function ROI_thresh_slide_call(varargin)
            
            % Covered by listener function below
            
        end
            
        function ROI_thresh_listener_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            ROI_thresh_slider_value = get(handles.handles.ROI_thresh_slide_hand, 'Value');
            
            if handles.Primary_channel == 1;
                min_max_here = handles.Min_max_left;
            elseif handles.Primary_channel == 2;
                min_max_here = handles.Min_max_right;
            end
            ROI_thresh_display = num2str(round(ROI_thresh_slider_value*(max(min_max_here(:,2) - min(min_max_here(:,1)))) + min(min_max_here(:,1))));
            
            set(handles.handles.ROI_thresh_slide_box, 'String', ROI_thresh_display);
            

          
            Threshold_ROI_display;
            
        end
            
        function ROI_thresh_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            box_value = round(str2double(get(handles.handles.ROI_thresh_slide_box, 'String')));
            
            if handles.Primary_channel == 1;
                min_max_here = handles.Min_max_left;
            elseif handles.Primary_channel == 2;
                min_max_here = handles.Min_max_right;
            end

            if box_value < min(min_max_here(:,1))
                
                box_value = min(min_max_here(:,1));
                
            elseif box_value > max(min_max_here(:,2))
                
                box_value = max(min_max_here(:,2));
                
            end
            
            frame_val = (box_value - (min(min_max_here(:,1))))/(max(min_max_here(:,2)) - min(min_max_here(:,1)));
            
            
            set(handles.handles.ROI_thresh_slide_box, 'String', num2str(box_value));
            set(handles.handles.ROI_thresh_slide_hand, 'Value', frame_val);
            

            
            Display_ROI_img('threshold');
            Threshold_ROI_display;
            
        end
        
        function Active_ROI_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            active_possible = get(handles.handles.Active_pulldown, 'String');
            active_now = get(handles.handles.Active_pulldown, 'Value');
            
            handles.Primary_channel = str2double(active_possible{active_now});
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off'); 
            
            if (get(handles.handles.ROI_frame_UseAverage, 'Value') == 1) 
                Display_ROI_img('MeanImg');
            elseif (get(handles.handles.ROI_frame_UseAverage, 'Value') == 1)
                Display_ROI_img('MaxImg');
            else
                Display_ROI_img('framechange');
            end
            
            if ~handles.ROIsAreSet
                Threshold_ROI_display;
            end
            
        end

        function Dilate_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Dilate_radius = str2double(get(handles.handles.Dilate_rad_box, 'String'));
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off'); 

            Threshold_ROI_display;
            
        end
        
        function End_caps_callback(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.End_caps = get(handles.handles.End_caps_hand, 'Value');
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off'); 
%             Calc_some_ROIs;
            
        end
        
        function ROI_length_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            tempString = (get(handles.handles.ROI_length_box, 'String'));
            
            if all(isstrprop(tempString, 'digit')) && str2double(tempString) > 0

                handles.ROI_length = str2double(tempString);
                
            else
                
                set(handles.handles.ROI_length_box, 'String', num2str(handles.ROI_length));
                
            end
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off');
%             Calc_some_ROIs;
            
        end
        
        function ROI_width_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            tempString = (get(handles.handles.ROI_width_box, 'String'));
            
            if all(isstrprop(tempString, 'digit')) && str2double(tempString) > 0

                handles.ROI_width = str2double(tempString);
                
            else
                
                set(handles.handles.ROI_width_box, 'String', num2str(handles.ROI_width));
                
            end
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off'); 
%             Calc_some_ROIs;
            
        end
        
        function ROI_smoothingSpan_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));            
            tempString = (get(handles.handles.ROI_smoothingSpan_box, 'String'));
            
            
            if all(isstrprop(tempString, 'digit') | tempString == '.') 
                % Entry is at least a number
                
                if (str2double(tempString) > 0)
                    % Number is greater than 0

                    handles.ROI_smoothing.span = round(str2double(tempString));
                    set(handles.handles.ROI_smoothingDegree_box, 'Enable', 'on');
                    guidata(findobj('Tag', 'TIFF viewer'), handles);
                    CheckDegreeValue(handles.ROI_smoothing.span);

                elseif (str2double(tempString) == 0)
                    % Number is equal to 0
                    % Disable smoothing
                    
                    handles.ROI_smoothing.span = 0;
                    set(handles.handles.ROI_smoothingDegree_box, 'Enable', 'off');
                    guidata(findobj('Tag', 'TIFF viewer'), handles);
%                     CheckDegreeValue(handles.ROI_smoothing.span);
                    
                else
                
                    set(handles.handles.ROI_smoothingSpan_box, 'String', num2str(handles.ROI_smoothing.span));
                    
                end

            else
               
                set(handles.handles.ROI_smoothingSpan_box, 'String', num2str(handles.ROI_smoothing.span));
                
            end

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off'); 
%             Calc_some_ROIs;
            
        end
        
        function CheckDegreeValue(valIn)
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            degValue = str2double(get(handles.handles.ROI_smoothingDegree_box, 'String'));
                       
            if (valIn >= degValue) && (valIn > 1)
                
                % Nothing to do - values are OK
                
            else
                
                if (valIn == 0) || (valIn == 1)
                    
                    
                    handles.ROI_smoothing.degree = 0;
                    
                else
                    
                    handles.ROI_smoothing.degree = handles.ROI_smoothing.span - 1;
                    
                end
                
            end
            
            set(handles.handles.ROI_smoothingDegree_box, 'String', num2str(handles.ROI_smoothing.degree));
            guidata(findobj('Tag', 'TIFF viewer'), handles);
        end
                
        
        function ROI_smoothingDegree_edit_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));            
            tempString = (get(handles.handles.ROI_smoothingDegree_box, 'String'));
            
            if all(isstrprop(tempString, 'digit') | tempString == '.') && (str2double(tempString) >= 0) && (str2double(tempString) < handles.ROI_smoothing.span)

                handles.ROI_smoothing.degree = round(str2double(tempString));
                               
            else
                
               %%% 
                
            end
            
            set(handles.handles.ROI_smoothingDegree_box, 'String', num2str(handles.ROI_smoothing.degree));
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            EnableDisablePlot('off'); 
%             Calc_some_ROIs;
            
        end
        
        function ROI_ROI(varargin)
            % Add a new ROI using impoly
           
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            polyHold = impoly(handles.handles.ax_ROI);

            % New position callbacks
            polyHold.addNewPositionCallback(@UpdateROIPosition);
            
            handles.polyROI = [handles.polyROI; polyHold];
            
            % Add area in polygon to Mask
            MaskHold = false(size(handles.ROI_Mask, 1), size(handles.ROI_Mask, 2));
            for k = 1:numel(handles.polyROI)
                MaskHold = (MaskHold | handles.polyROI(k).createMask(handles.handles.ROI_image));
            end
            handles.ROI_Mask = MaskHold;

            % Incorporate mask into threshold
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Threshold_ROI_display;
            
        end
        
        function UpdateROIPosition(varargin)
            % Callback for resetting ROI position
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            handles.polyROI(~handles.polyROI.isvalid()) = [];

            MaskHold = false(size(handles.ROI_Mask, 1), size(handles.ROI_Mask, 2));
            for k = 1:numel(handles.polyROI)
                MaskHold = (MaskHold | handles.polyROI(k).createMask(handles.handles.ROI_image));
            end
            handles.ROI_Mask = MaskHold;

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Threshold_ROI_display;
            
        end
        
        function ROI_apply(varargin)
            
            % Apply settings and calculate backbone traces
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            set(handles.handles.fig3, 'pointer', 'watch');
            %disp('Start delay')
            % Find everything in relevant panel that is turned on. Turn it
            % off.
            On_handles = findobj('Parent', handles.handles.fig3_button_panel, 'Enable', 'on');
            set(On_handles, 'Enable', 'off');
            drawnow;
 
            Calc_some_ROIs
            
            waitfor(handles.handles.fig3, 'UserData', 'ROI_calced');
            set(handles.handles.fig3, 'pointer', 'arrow');
            set(On_handles, 'Enable', 'on');
            set(handles.handles.fig3, 'UserData', []);
            
            EnableDisablePlot('on'); 
            
            drawnow;
            %disp('End Delay');
                    
        end
                
        function ROI_close(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            set(handles.handles.fig3, 'pointer', 'watch');
            %disp('Start delay')
            On_handles = findobj('Parent', handles.handles.fig3_button_panel, 'Enable', 'on');
            set(On_handles, 'Enable', 'off');
            drawnow;
 
            Calc_some_ROIs
            
            %Plot_some_ROI_boxes;
            
            waitfor(handles.handles.fig3, 'UserData', 'ROI_calced');
            set(handles.handles.fig3, 'pointer', 'arrow');
            set(On_handles, 'Enable', 'on');
            set(handles.handles.fig3, 'UserData', []);
            
            EnableDisablePlot('on'); 
            
            drawnow;
            %disp('End Delay');

            close(handles.handles.fig3);
            
            
        end
        
        function CloseROIFcn(varargin)  %Function called when Background trace window is closed
            
            handles = guidata(handles.handles.fig1);
            
            handles.polyROI = [];
            
            handles.ROI_Mask = ones(size(handles.Img_stack, 1), size(handles.Img_stack, 2));
           
            guidata(handles.handles.fig1, handles);
        end
            

    end

%%%%%%%%%%%%%%%%%%%%%%
% ROI plot display

    function Display_ROI_img(Reason) %Display image in Background Trace window.  
        % Function varies depending on reason for callback and what needs
        % to be drawn in window due to changes elsewhere
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        if ~isempty(findobj('Tag', 'GALAH_ROI_calc'));
            
            set(handles.handles.ax_ROI, 'NextPlot', 'add');

            if strcmp(Reason, 'threshold')
                % Clear out old lines, patches if this was a change in frame
                % number to call Display_ROI_img

                delete(findobj('Parent', handles.handles.ax_ROI, 'Type', 'line'));
                delete(findobj('Parent', handles.handles.ax_ROI, 'Type', 'patch'));
                
                handles.ROIsAreSet = false;
            end
            
            % Clear out old image
            delete(findobj('Tag', 'ROI_image'));

            ax_ROI = handles.handles.ax_ROI;
            frame_slider = str2double(get(handles.handles.ROI_frame_slide_box, 'String')); % Set this to ROI_slider

            if handles.Primary_channel == 1;
                min_max_here = handles.Min_max_left;
                display_here = handles.Display_range_left;
                invert_here = handles.Left_invert;
                Autoscale_here = handles.Autoscale_left;
                yaxis_direction = get(handles.handles.ax1, 'YDir');
            elseif handles.Primary_channel == 2;
                min_max_here = handles.Min_max_right;
                display_here = handles.Display_range_right;
                invert_here = handles.Right_invert;
                Autoscale_here = handles.Autoscale_right;
                yaxis_direction = get(handles.handles.ax2, 'YDir');
            end

            % Check if Autoscale
            
            if Autoscale_here == 1;
                display_here = min_max_here(frame_slider,:);
            end
            
            OldXLimits = get(ax_ROI, 'XLim');
            OldYLimits = get(ax_ROI, 'YLim');
            
            if strcmp(Reason, 'MeanImg')
                
                ROI_image = image(Vector2Colormap_setscale(handles.MeanImg(:,:,handles.Primary_channel), 'gray', display_here), ...
                    'Parent', ax_ROI, 'Tag', 'ROI_image');
            elseif strcmp(Reason, 'MaxImg')
                ROI_image = image(Vector2Colormap_setscale(handles.MaxImg(:,:,handles.Primary_channel), 'gray', display_here), ...
                    'Parent', ax_ROI, 'Tag', 'ROI_image');
            else 
                ROI_image = image(Vector2Colormap_setscale(handles.Img_stack(:,:,frame_slider,handles.Primary_channel), 'gray', display_here), ...
                    'Parent', ax_ROI, 'Tag', 'ROI_image');
                set(ax_ROI, 'xtick', [], 'ytick', []);
                set(ax_ROI, 'XLim', OldXLimits, 'YLim', OldYLimits);
            end

            %axis image
            
            % Invert if needed
            if invert_here == 1;
                set(ax_ROI, 'XDir', 'reverse');
            else
                set(ax_ROI, 'XDir', 'normal');
            end
            
            set(ax_ROI, 'YDir', yaxis_direction);

            % Move image to back
            uistack(ROI_image, 'bottom');

            handles.handles.ROI_image = ROI_image;
            
            set(handles.handles.ax_ROI, 'NextPlot', 'replace');
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
    end

%%%%%%%%%%%%%%%%%%%%%%
% ROI Threshold Display

    function Threshold_ROI_display(varargin)
        % Calculate threshold mask and display threshold region outlines in
        % Background Trace window
     
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        handles.ROIsAreSet = false;
        set(handles.handles.ROI_Plots_launch_hand, 'Enable', 'off');

            OldXLimits = get(handles.handles.ax_ROI, 'XLim');
            OldYLimits = get(handles.handles.ax_ROI, 'YLim');
        
        % Clear out old lines, ROI values
        delete(findobj('Parent', handles.handles.ax_ROI, 'Type', 'line'));
        delete(findobj('Parent', handles.handles.ax_ROI, 'Type', 'patch'));
        handles.handles.ROI_line = [];
        handles.handles.ROI_box ={};
        handles.handles.ROI_color = [];
        handles.ROI_output = {};
        handles.ROI_parameters = [];
        
        guidata(handles.handles.fig1, handles);
        
        EnableDisablePlot('off'); 
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
       
        frame_here = str2double(get(handles.handles.ROI_frame_slide_box, 'String'));
        
        if handles.Primary_channel == 1;
            min_max_here = handles.Min_max_left;
%             display_here = handles.Display_range_left;
        elseif handles.Primary_channel == 2;
            min_max_here = handles.Min_max_right;
%             display_here = handles.Display_range_right;
        end
        
        min_max_here = [min(min_max_here(:,1)), max(min_max_here(:,2))];
        
        if get(handles.handles.ROI_frame_UseAverage, 'Value') == 1
            
            scale_frame = (handles.MeanImg(:,:,handles.Primary_channel) - min_max_here(1))./ ...
                (min_max_here(2) - min_max_here(1));
        elseif get(handles.handles.ROI_frame_UseMaxProjection, 'Value') == 1
            scale_frame = (handles.MaxImg(:,:,handles.Primary_channel) - min_max_here(1))./ ...
                (min_max_here(2) - min_max_here(1));
        else
        
            scale_frame = (handles.Img_stack(:, :, frame_here, handles.Primary_channel) - min_max_here(1))./ ...
                (min_max_here(2) - min_max_here(1));
        end
        
        thresh_here = get(handles.handles.ROI_thresh_slide_hand, 'Value');
        
        dilate_radius = str2double(get(handles.handles.Dilate_rad_box, 'String'));
            
        Dilate_kernel = strel('disk', dilate_radius);

        BW = im2bw(scale_frame, thresh_here);
        
        % Apply ROI Mask 
        
        BW = (BW & handles.ROI_Mask);
        
        BW2 = bwmorph(BW, 'clean');
        BW2 = imdilate(BW2, Dilate_kernel);

        
        [B,L] = bwboundaries(BW2,'noholes');

        border_colors = handles.handles.border_colors;
        Border_randomize = handles.handles.Border_randomize;
        
        
        set(handles.handles.ax_ROI, 'NextPlot', 'add');
        
        Border_plot = zeros(length(B), 1);
        if ~isempty(Border_plot)

            for k = 1:length(B);
                
                color_mod = mod(k, size(border_colors,1)) + 1;
                
                boundary = B{(k)};
                Border_plot(k) = plot(handles.handles.ax_ROI, boundary(:,2), boundary(:,1), 'Color', ...
                    border_colors(Border_randomize(color_mod),:), 'LineWidth', 1, 'HitTest', 'off');
                
                % Determine range for zoom image later.  Should be a square
                % that encloses ROI with borders of ROI_zoom_border;

                range_here = [min(B{k},[],1) max(B{k},[],1)];
                range_here = range_here + (repmat(handles.ROI_zoom_border, 1, 4).*[-1 -1 1 1]);
                ROI_center = [mean(range_here([1,3])) mean(range_here([2 4]))];
                square_size = max([range_here(3) - range_here(1), range_here(4) - range_here(2)]);
                ROI_range{k} = [ROI_center(1) - square_size/2, ROI_center(1) + square_size/2, ROI_center(2) - square_size/2, ROI_center(2) + square_size/2];
                
                
            end

            handles.handles.Border_plot = Border_plot;
            handles.ROI_boundaries.B = B;
            handles.ROI_boundaries.L = L;
            handles.ROI_boundaries.bound_square = ROI_range;
        else 
            handles.handles.Border_plot = [];
            handles.ROI_boundaries.B = [];
            handles.ROI_boundaries.L = [];
        end
        
        set(handles.handles.ax_ROI, 'NextPlot', 'replace');
        axis(handles.handles.ax_ROI, 'image');
        set(handles.handles.ax_ROI, 'XLim', OldXLimits, 'YLim', OldYLimits);
        
%         disp(handles.ROIsAreSet)
        
        guidata(findobj('Tag', 'TIFF viewer'), handles);
        
        
    end


%%%%%%%%%%%%%%%%%%%%%%
% Launch ROI plot viewer

    function Plots_launch(varargin)
        % Show plots for currently-selected backbone trace
        % Lets you scroll through generated traces and provides launcher
        % for Kymographs

        % This gets reset if traces have been re-calc'd, too.
       
       if ~isempty(findobj('Tag', 'GALAH_plot_display'))
        
            uistack(findobj('Tag', 'GALAH_plot_display'), 'top');
            
       else
                        
            fig1 = findobj('Tag', 'TIFF viewer');
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            if ~isempty(handles.ROI_output) && handles.ROIsAreSet;
            
                mf_post = get(fig1, 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);      
                fig8_size = [700 400];
                fig8_position = [(mf_post(1) + (mf_post(3) - fig8_size(1))/2) (mf_post(2) + (mf_post(4) - fig8_size(2))/2)];
                fig8 = figure('Name','ROI Property Plotter', 'Tag', 'GALAH_plot_display', 'Units', 'pixels', 'MenuBar', 'none',...
                    'Position',[fig8_position fig8_size]);
                set(fig8, 'Color',[0.9 0.9 0.9]);
                set(fig8, 'DeleteFcn', @ROI_plot_close_fcn);
                set(fig8, 'NumberTitle', 'off');

                fig8_img_panel = uipanel(fig8, 'Units', 'normalized', 'Position', [0, .1, 1, .9], ...
                    'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'none', 'Tag', 'Plot_plot_panel');

                fig8_button_panel = uipanel(fig8, 'Units', 'normalized', 'Position', [0, 0, 1, .1], ...
                    'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'none', 'Tag', 'Plot_button_panel');



                % Set up Axis

                ax_ROI_plot = axes('Parent', fig8_img_panel, 'Position', [.10 .14 .83 .80]);
                set(ax_ROI_plot, 'Tag', 'Axis_Review');

                handles.handles.fig8_img_panel = fig8_img_panel;
                handles.handles.fig8_button_panel = fig8_button_panel;
                handles.handles.ax_ROI_plot = ax_ROI_plot;

                guidata(findobj('Tag', 'TIFF viewer'), handles);

                % Call function to update or initialize figure contents

                FillInPlotWindow;

            else
                errordlg(sprintf('No ROI settings have been applied. \nClick ''Apply Settings'' to Create ROIs.'), 'No ROIs Applied');
            end
            
       end
       
       function ROI_plot_close_fcn(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));

            % Reset purple boxes to green 
        
            purple_objects = findobj('Parent', handles.handles.ax_ROI, 'Color', handles.handles.Selected_ROI_color);
            set(purple_objects, 'Color', handles.handles.Unselected_ROI_color);
            
            close(findobj('Tag', 'Kymograph Figure'));
            close(findobj('Tag', 'KymographCurve'));
            
        end
   
        
    end


            
	function FillInPlotWindow(varargin)
        % Fill in values in the plot window based on number of channels
        % available
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        % Check to make sure plots window is open
        if ~isempty(findobj('Tag', 'GALAH_plot_display'))
        
            % Set up ROI slider
            
            ROI_plot_slider_value = 1;
            ROI_plot_display = num2str(1);
                        
            if length(handles.ROI_output) > 1;
                plot_slide_max = length(handles.ROI_output);
                popup_plot_options = [{' '};{'X Position'}; {'Y Position'}; {'ROI Area'}; {'ROI Mean'}; {'ROI StdDev'}];
                popup_start = 5;
                %popup_plot_options = fieldnames(handles.ROI_output);
                ROI_plot_step = 1/(length(handles.ROI_output) - 1);
                disable_slider = 'on';
            elseif length(handles.ROI_output) == 1
                plot_slide_max = 2;
                popup_plot_options = [{' '};{'X Position'}; {'Y Position'}; {'ROI Area'}; {'ROI Mean'}; {'ROI StdDev'}];
                popup_start = 5;
                ROI_plot_step = .8;
                disable_slider = 'off';
            elseif isempty(handles.ROI_output)
                plot_slide_max = 2;
                popup_plot_options = {' '};
                popup_start = 1;
                ROI_plot_step = 1;
                disable_slider = 'off';
            end

            ROI_plot_slide_hand = uicontrol(handles.handles.fig8_button_panel, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [ROI_plot_step ROI_plot_step], 'Min', 1, 'Max', plot_slide_max, ...
                'Value', ROI_plot_slider_value, 'Position', [.19 .90 .44 .325],...
                'Callback', @ROI_plot_slider_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'ROI_num');

            addlistener(ROI_plot_slide_hand, 'Value', 'PostSet', @ROI_plot_listener_call);

            ROI_plot_slide_box = uicontrol(handles.handles.fig8_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.11 .825 .06 .5], 'BackgroundColor', [1 1 1], ...
                'String', ROI_plot_display, 'Callback', @ROI_plot_edit_call);

            uicontrol(handles.handles.fig8_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.016 .75 .09 .5], 'BackgroundColor', [.9 .9 .9], ...
                'String', 'ROI Number:');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if handles.N_frames == 1;
            
                ROI_plot_Framestep = .5;
                ROI_plot_Frameslider_value = 1;
                ROI_FrameMinMax = [0 1];
                ROI_Frameplot_display = num2str(ROI_plot_Frameslider_value);
                EnableDisableFrame = 'off';

            else
                
                ROI_plot_Framestep = 1/(handles.N_frames - 1);
                ROI_plot_Frameslider_value = str2double(get(handles.handles.ROI_frame_slide_box, 'String'));
                ROI_Frameplot_display = get(handles.handles.ROI_frame_slide_box, 'String');
                ROI_FrameMinMax = [1 handles.N_frames];
                EnableDisableFrame = 'on';
                
            end
                
            
            ROI_plot_Frameslide_hand = uicontrol(handles.handles.fig8_button_panel, 'Style', 'slider', 'Units', 'normalized',...  
                'SliderStep', [ROI_plot_Framestep ROI_plot_Framestep], 'Min', ROI_FrameMinMax(1), 'Max', ROI_FrameMinMax(2), ...
                'Value', ROI_plot_Frameslider_value, 'Position', [.19 .25 .44 .325], 'Enable', EnableDisableFrame,...
                'Callback', @ROI_plot_Frameslider_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'ROI_frame');

            addlistener(ROI_plot_Frameslide_hand, 'Value', 'PostSet', @ROI_plot_Framelistener_call);

            ROI_plot_Frameslide_box = uicontrol(handles.handles.fig8_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
                'Position', [.11 .175 .06 .5], 'BackgroundColor', [1 1 1], 'Enable', EnableDisableFrame,...
                'String', ROI_Frameplot_display, 'Callback', @ROI_plot_Frameedit_call);

            uicontrol(handles.handles.fig8_button_panel, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [.016 .1 .09 .5], 'BackgroundColor', [.9 .9 .9],'Enable', EnableDisableFrame, ...
                'String', 'Frame:');
            
            %%%%
            
%             ROI_plot_popup_text = uicontrol(handles.handles.fig8_button_panel, 'Style', 'text', 'Units', 'normalized', ...
%                 'Position', [.65 .155 .07 .4], 'BackgroundColor', [.9 .9 .9], ...
%                 'String', 'Plot Field:');
            
            
             ROI_plot_popup = uicontrol(handles.handles.fig8_button_panel, 'Style', 'popupmenu', 'Units', 'normalized', ...
                'Position', [.66 .45 .12 .5], 'String', popup_plot_options, 'Value', popup_start, 'Callback', @popup_plot_option);

            
            %%%%
            
            ROI_plot_2D = uicontrol(handles.handles.fig8_button_panel, 'Style', 'pushbutton', 'Units', 'normalized',...
                'Position', [.851 .8 .14 .625], 'String', 'Kymograph', 'Callback', @ROI_plot_2D_launch);
            %%%%
            
            ROI_plot_close = uicontrol(handles.handles.fig8_button_panel, 'Style', 'pushbutton', 'Units', 'normalized',...
                'Position', [.851 .05 .14 .625], 'String', 'OK', 'Callback', @ROI_plot_close_push);
            
            
            % Set up close button


            handles.handles.ROI_plot_slide_hand = ROI_plot_slide_hand;
            handles.handles.ROI_plot_slide_box = ROI_plot_slide_box;
            handles.handles.ROI_plot_close = ROI_plot_close;
            handles.handles.ROI_plot_popup = ROI_plot_popup;
            handles.handles.ROI_plot_Frameslide_hand = ROI_plot_Frameslide_hand;
            handles.handles.ROI_plot_Frameslide_box = ROI_plot_Frameslide_box;
            handles.handles.ROI_plot_2D = ROI_plot_2D;
            
            
            set(handles.handles.ROI_plot_slide_hand, 'Enable', disable_slider);
            set(handles.handles.ROI_plot_slide_box, 'Enable', disable_slider);
            
            
            
            if size(handles.Img_stack, 3)==1
                set(handles.handles.ROI_plot_2D, 'Enable', 'off')
            end
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            Display_ROI_plot;
            ROI_highlight;
        
        end
            
            function ROI_plot_slider_call(varargin)
                
                % Necessary function but does nothing. 
                % Operations covered by listener.

            end
       
            function ROI_plot_listener_call(varargin)

                handles = guidata(findobj('Tag', 'TIFF viewer'));
                ROI_value = round(get(handles.handles.ROI_plot_slide_hand, 'Value'));
                set(handles.handles.ROI_plot_slide_box, 'String', num2str(ROI_value));

                Display_ROI_plot;

                ROI_highlight;


            end
        
            function ROI_plot_edit_call(varargin)

                handles = guidata(findobj('Tag', 'TIFF viewer'));
                box_val = round(str2double(get(handles.handles.ROI_plot_slide_box, 'String')));

                if box_val > get(handles.handles.ROI_plot_slide_hand, 'Max')
                    box_val = get(handles.handles.ROI_plot_slide_hand, 'Max');
                elseif box_val < get(handles.handles.ROI_plot_slide_hand, 'Min');
                    box_val = get(handles.handles.ROI_plot_slide_hand, 'Min');
                end

                set(handles.handles.ROI_plot_slide_box, 'String', num2str(box_val));
                set(handles.handles.ROI_plot_slide_hand, 'Value', box_val);


            end
            
            %%%%%%%%%%%%%%
            
            function ROI_plot_Frameslider_call(varargin)
                
                % Necessary function that can remain empty
                % Listener does everything

            end

            function ROI_plot_Framelistener_call(varargin)
                
                handles = guidata(findobj('Tag', 'TIFF viewer'));
                ROI_value = round(get(handles.handles.ROI_plot_Frameslide_hand, 'Value'));
                set(handles.handles.ROI_plot_Frameslide_box, 'String', num2str(ROI_value));

                Display_ROI_plot;

            end

            function ROI_plot_Frameedit_call(varargin)
                
                handles = guidata(findobj('Tag', 'TIFF viewer'));
                box_val = round(str2double(get(handles.handles.ROI_plot_Frameslide_box, 'String')));

                if box_val > get(handles.handles.ROI_plot_Frameslide_hand, 'Max')
                    box_val = get(handles.handles.ROI_plot_Frameslide_hand, 'Max');
                elseif box_val < get(handles.handles.ROI_plot_Frameslide_hand, 'Min');
                    box_val = get(handles.handles.ROI_plot_Frameslide_hand, 'Min');
                end

                set(handles.handles.ROI_plot_Frameslide_box, 'String', num2str(box_val));
                set(handles.handles.ROI_plot_Frameslide_hand, 'Value', box_val);

            end
            
            %%%%%%%%%%%%%%
        
        
            function popup_plot_option(varargin)

                Display_ROI_plot

            end
        
            function ROI_plot_close_push(varargin)

                close(findobj('Tag', 'GALAH_plot_display'));

            end
        

            
    end


%%%%%%%%%%%%%%%%%%%%%%
% Display ROI plot

    function Display_ROI_plot(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
            frame_here = round(get(handles.handles.ROI_plot_Frameslide_hand, 'Value'));
       
            ROI_num = str2double(get(handles.handles.ROI_plot_slide_box, 'String'));
            ROI_field_num = get(handles.handles.ROI_plot_popup, 'Value');
            ROI_strings = get(handles.handles.ROI_plot_popup, 'String');

            if ROI_field_num == 1;
                % Clear plot
                delete(findobj('Parent', handles.handles.ax_ROI_plot));
                ylabel('');
                
                set(handles.handles.ROI_plot_2D, 'Enable', 'off');

            else
                
                if size(handles.Img_stack, 3) == 1
                    set(handles.handles.ROI_plot_2D, 'Enable', 'off');
                else
                    set(handles.handles.ROI_plot_2D, 'Enable', 'on');
                end
                ROI_field_names = [{'x_center'};{'y_center'};{'areas'};{'Mean_mask'};{'std_mask'}];
                data_here = handles.ROI_output{ROI_num}(1, frame_here).(ROI_field_names{ROI_field_num-1});
                x_post = handles.ROI_output{ROI_num}(1).x_center;
                y_post = handles.ROI_output{ROI_num}(1).y_center;
                domain = sqrt((x_post-circshift(x_post, 1)).^2 + (y_post-circshift(y_post, 1)).^2);
                domain(1) = 0;
                domain = cumsum(domain);

            ROI_plot_line = plot(handles.handles.ax_ROI_plot, domain, data_here, '--o');
            set(ROI_plot_line, 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', [.4 .2 1], ...
                'MarkerFaceColor', [.6 .6 1], 'Color', [.4 .2 1]); 
            
            if handles.N_channels == 2;
                set(handles.handles.ax_ROI_plot, 'NextPlot', 'add')
                
                data_here = handles.ROI_output{ROI_num}(2, frame_here).(ROI_field_names{ROI_field_num-1});
                
                ROI_plot_line(2) = plot(handles.handles.ax_ROI_plot, domain, data_here, '--o');
                set(ROI_plot_line(2), 'LineWidth', 2, 'MarkerSize', 8, 'MarkerEdgeColor', [1 .4 .2], ...
                'MarkerFaceColor', [1 .8 .6], 'Color', [1 .4 .2]); 
                set(handles.handles.ax_ROI_plot, 'NextPlot', 'replace');
            end
            
            xlabel('Distance (pixels)'); 
            ylabel(ROI_strings{ROI_field_num});
            
            

            end

            %disp(sprintf('Plot ROI %d, Field %d', ROI_num, ROI_field_num));
            guidata(findobj('Tag', 'TIFF viewer'), handles);

        
    end


%%%%%%%%%%%%%%%%%%%%%%%
% Pop-out figure for dist vs time (kymograph) with current data displayed

 function ROI_plot_2D_launch(varargin)
     
    handles = guidata(findobj('Tag', 'TIFF viewer'));
    
%     frame_here = round(get(handles.handles.ROI_plot_Frameslide_hand, 'Value'));
    ROI_num = str2double(get(handles.handles.ROI_plot_slide_box, 'String'));
    ROI_field_num = get(handles.handles.ROI_plot_popup, 'Value');
    ROI_strings = get(handles.handles.ROI_plot_popup, 'String');
     
	ROI_field_names = [{'x_center'};{'y_center'};{'areas'};{'Mean_mask'};{'std_mask'}];
	data_green = horzcat(handles.ROI_output{ROI_num}(1, :).(ROI_field_names{ROI_field_num-1}));
	x_post = handles.ROI_output{ROI_num}(1).x_center;
	y_post = handles.ROI_output{ROI_num}(1).y_center;
	domain = sqrt((x_post-circshift(x_post, 1)).^2 + (y_post-circshift(y_post, 1)).^2);
	domain(1) = 0;
	interpDomain = cumsum(domain);
    
%     interpDomainStep = mean(diff(domain));
%     interpDomain = 0:interpDomainStep:max(domain);
%     [X, Y] = meshgrid(domain, 1:handles.N_frames);
%     [Xq, Yq] = meshgrid(interpDomain, 1:handles.N_frames);
    
%     interpGreen = interp2(X, Y, data_green', Xq, Yq)';
    % Changed here to show true data, not interpolated data
%     interpGreenMap = Vector2Colormap(interpGreen, 'green');
    
    interpGreenMap = repmat(flipud(data_green), [1, 1, 3]).*reshape(repmat([.6 .6 1], [numel(data_green), 1]), [size(data_green, 1), size(data_green, 2), 3]);
    interpGreenMap = (interpGreenMap - min(interpGreenMap(:)))/(max(interpGreenMap(:)) - min(interpGreenMap(:)));
    
    handles.Kymograph(ROI_num).PixelDistance = (interpDomain);
    handles.Kymograph(ROI_num).Channel1 = flipud(data_green);
    
    KymographSize = [600 250];

    if handles.N_channels == 2;
        
        data_red = horzcat(handles.ROI_output{ROI_num}(2, :).(ROI_field_names{ROI_field_num-1}));
%         interpRed = interp2(X, Y, data_red', Xq, Yq)';
        % Changed here to show true data, not interpolated data
%         interpRedMap = Vector2Colormap(interpRed, 'red');
        interpRedMap = repmat(flipud(data_red), [1, 1, 3]).*reshape(repmat([1 .8 .6], [numel(data_red), 1]), [size(data_red, 1), size(data_red, 2), 3]);
        interpRedMap = (interpRedMap - min(interpRedMap(:)))/(max(interpRedMap(:)) - min(interpRedMap(:)));
        
        KymographSize = [600 450];
        
        handles.Kymograph(ROI_num).Channel2 = flipud(data_red);

    end
    
    % Make figure
    handles.handles.lastKymograph = figure();
    
    FigPost = get(handles.handles.lastKymograph, 'Position');
    set(handles.handles.lastKymograph, 'Position', [FigPost(1) FigPost(2)-100 KymographSize(1) KymographSize(2)]);
    KymographTitle = sprintf('Kymograph - ROI #%d - %s', ROI_num, ROI_strings{ROI_field_num});
    set(handles.handles.lastKymograph, 'Tag', 'Kymograph Figure', 'Name', KymographTitle, 'Units', 'pixels', ...
        'Toolbar', 'none', 'MenuBar', 'none');
    set(handles.handles.lastKymograph, 'Color',[0.9 0.9 0.9]);
    set(handles.handles.lastKymograph, 'NumberTitle', 'off');
    
    if handles.N_channels == 1
        handles.handles.lastKymoAxis = axes('Parent', handles.handles.lastKymograph, 'Position', [.05 .234 .94 .85]);
        
        handles.handles.lastKymoImg(1) = image(interpDomain, 1:handles.N_frames, interpGreenMap, 'Parent', handles.handles.lastKymoAxis(1));
        set(handles.handles.lastKymoAxis, 'xTickLabel', [], 'yTickLabel', []);
        set(handles.handles.lastKymoAxis, 'XTick', [], 'YTick', []);
        ylabel('Distance (pixels)'); xlabel('Time (frames)')
        
    elseif handles.N_channels == 2
        
        handles.handles.lastKymoAxis(1) = axes('Parent', handles.handles.lastKymograph, 'Position', [.05 .13 .94 .43]);
        handles.handles.lastKymoAxis(2) = axes('Parent', handles.handles.lastKymograph, 'Position', [.05 .56 .94 .43]);
        
        handles.handles.lastKymoImg(1) = image(interpDomain, 1:handles.N_frames, interpGreenMap, 'Parent', handles.handles.lastKymoAxis(1));
        handles.handles.lastKymoImg(2) = image(interpDomain, 1:handles.N_frames, interpRedMap, 'Parent', handles.handles.lastKymoAxis(2));
        
        set(handles.handles.lastKymoAxis, 'xTickLabel', [], 'yTickLabel', []);
        set(handles.handles.lastKymoAxis, 'XTick', [], 'YTick', []);
        ylabel(handles.handles.lastKymoAxis(1), 'Distance (pixels)'); xlabel(handles.handles.lastKymoAxis(1), 'Time (frames)')
        ylabel(handles.handles.lastKymoAxis(2), 'Distance (pixels)');
    end
    buttonHeight(1) = .1;
    buttonHeight(2) = .06;
    
    
    % Add buttons
    handles.handles.ExportKymographData = uicontrol(handles.handles.lastKymograph, 'Style', 'pushbutton', 'Units', 'normalized',...
                'Position', [.85 .015 .14 buttonHeight(handles.N_channels)], 'String', 'Export Data', 'Callback', @KymoDataOut);
            
    handles.handles.DisplayKymographCurves = uicontrol(handles.handles.lastKymograph, 'Style', 'pushbutton', 'Units', 'normalized',...
                'Position', [.70 .015 .14 buttonHeight(handles.N_channels)], 'String', 'Curve Stack', 'Callback', @KymoCurveStack);
    
            
            
	guidata(findobj('Tag', 'TIFF viewer'), handles);

    
     function KymoDataOut(varargin)
         
         handles = guidata(handles.handles.fig1);
         
         [fname, pname, filtIndex] = uiputfile({'*.csv','CSV File (*.csv)'; ...
            '*.xls','Excel File (*.xls)';       
            '*.tif',  'TIFF File (*.tif)';
            '*.csv',  'Summary Data (*.csv)'}, ...
            'Pick a file', 'DataExport');
        
        switch filtIndex
            case 1
                
                % Export kymograph as CSV file
                % Data bookended by blank lines, rows of ------'s
                
                expObj = fopen(fullfile(pname, fname), 'wt');

                fwrite(expObj, sprintf('Data File : %s\r\n', handles.Load_file));
                fwrite(expObj, sprintf('ROI Number : %.f\r\n', ROI_num));
                fwrite(expObj, sprintf('Data Type : %s\r\n', ROI_strings{ROI_field_num}));
                
                
                for k = 1:handles.N_channels
                    fwrite(expObj, sprintf('\r\n'));
                    fwrite(expObj, sprintf('--------------------------------------------------------\r\n'));
                    fwrite(expObj, sprintf('Channel %.f\r\n', k));
                    fwrite(expObj, sprintf('First row is interpolated distance along trace from first ROI center in units of pixels \r\n'));
                    fwrite(expObj, sprintf('Subsequent rows correspond to one ROI through time \r\n'));
                    fwrite(expObj, sprintf('--------------------------------------------------------\r\n'));
                    fwrite(expObj, sprintf('\r\n'));
                    
                    if k == 1
                        BlockToWrite = ([interpDomain', data_green])';
                    else
                        BlockToWrite = ([interpDomain', data_red])';
                    end
                    
                    formatString = repmat('%.2f,', 1, size(BlockToWrite, 2));
                    formatString = strcat(formatString(1:end-1), '\r\n');
                    
                    for m = 1:size(BlockToWrite, 1)
                        fwrite(expObj, sprintf(formatString, BlockToWrite(m, :)));
                    end

                end
                
                fclose(expObj);


            case 2
                
                % Export kymograph as Excel file
                % Data begins on A5
                
                for k = 1:handles.N_channels

                    if k == 1
                        BlockToWrite = ([interpDomain', data_green])';
                    else
                        BlockToWrite = ([interpDomain', data_red])';
                    end
                    xlswrite(fullfile(pname, fname), {'Data File', handles.Load_file}, k, 'A1');
                    xlswrite(fullfile(pname, fname), {'ROI Number', ROI_num}, k, 'A2');
                    xlswrite(fullfile(pname, fname), {'Data Type', ROI_strings{ROI_field_num}}, k, 'A3');
                    xlswrite(fullfile(pname, fname), {'Channel', k}, k, 'A4');
                    xlswrite(fullfile(pname, fname), BlockToWrite, k, 'A5');
                    
                end


            case 3

                % Export figure as TIFF image
                % Two-channel data results in two-channel TIFF stack
                
                for k = 1:handles.N_channels
                
                    if k == 1
                        imwrite(uint16(data_green), fullfile(pname, fname), 'TIFF');
                    elseif k == 2
                        imwrite(uint16(data_red), fullfile(pname, fname), 'TIFF', 'writemode', 'append');
                    end
                    
                end
                
            case 4
                
                % Export summary data from kymograph as CSV file
                % Format is one column time (frames), then next two
                % the mean values in kymograph for channel 1 and, if
                % needed, channel 2
                % Data bookended by blank lines, rows of ------'s
                
                expObj = fopen(fullfile(pname, fname), 'wt');

                fwrite(expObj, sprintf('Data File : %s\r\n', handles.Load_file));
                fwrite(expObj, sprintf('ROI Number : %.f\r\n', ROI_num));
                fwrite(expObj, sprintf('Data Type : %s\r\n', ROI_strings{ROI_field_num}));
                fwrite(expObj, sprintf('\r\n'));
                fwrite(expObj, sprintf('--------------------------------------------------------\r\n'));
                
                switch handles.N_channels
                
                    case 1
                        fwrite(expObj, sprintf('Time (frames),Mean Chan 1\r\n'));
                        for k = 1:handles.N_frames
                            fprintf(expObj, '%.3f,%.3f\r\n', k, mean(data_green(:,k)));
                        end
                    case 2
                        fwrite(expObj, sprintf('Time (frames),Mean Chan 1,Mean Chan 2\r\n'));
                        for k = 1:handles.N_frames
                            fprintf(expObj, '%.3f,%.3f,%.3f\r\n', k, mean(data_green(:,k)), mean(data_red(:,k)));
                        end

                end
                                
                fclose(expObj);

        end
  
         
     end
     
     function KymoCurveStack(varargin)
         
         handles.handles.KymoCurveFig = figure('Name', sprintf('Kymograph Curves, ROI %.f', ROI_num),...
             'Tag', 'KymographCurve', 'NumberTitle', 'off');
         postNow = get(handles.handles.KymoCurveFig, 'Position');
         
         
         
         if handles.N_channels == 2
             set(handles.handles.KymoCurveFig, 'Position', [postNow(1:2), 900, 400]);
             handles.handles.KymoCurveAxes(1) = axes('Parent', handles.handles.KymoCurveFig, ...
                 'Units', 'normalized', 'Position', [.055 .1 .44 .83]);
             handles.handles.KymoCurveAxes(2) = axes('Parent', handles.handles.KymoCurveFig, ...
                 'Units', 'normalized', 'Position', [.55 .1 .44 .83]);

             handles.handles.KymoCurveSingles{1} = plot(handles.handles.KymoCurveAxes(1), data_green', 'Color', [.82 .74 1]);
             handles.handles.KymoCurveSingles{2} = plot(handles.handles.KymoCurveAxes(2), data_red', 'Color', [1 .8 .68]);
             set(handles.handles.KymoCurveAxes, 'NextPlot', 'add');

             handles.handles.KymoCurveMean{1} = plot(handles.handles.KymoCurveAxes(1), mean(data_green, 1), 'Color', [.4 .2 1], 'LineWidth', 2);
             handles.handles.KymoCurveMean{2} = plot(handles.handles.KymoCurveAxes(2), mean(data_red, 1), 'Color', [1 .4 .2], 'LineWidth', 2);
             set(handles.handles.KymoCurveAxes, 'NextPlot', 'replace');

             set(handles.handles.KymoCurveAxes, 'XLim', [1 handles.N_frames]);
             xlabel(handles.handles.KymoCurveAxes(1), 'Time (frame)');
             ylabel(handles.handles.KymoCurveAxes(1), ROI_strings{ROI_field_num});
             xlabel(handles.handles.KymoCurveAxes(2), 'Time (frame)');
             ylabel(handles.handles.KymoCurveAxes(2), ROI_strings{ROI_field_num});
             title(handles.handles.KymoCurveAxes(1), 'Channel 1');
             title(handles.handles.KymoCurveAxes(2), 'Channel 2');
             
             leghandles1 = [handles.handles.KymoCurveSingles{1}(1); handles.handles.KymoCurveMean{1}(1)];
             leghandles2 = [handles.handles.KymoCurveSingles{2}(1); handles.handles.KymoCurveMean{2}(1)];
             
             legend(handles.handles.KymoCurveAxes(1), leghandles1, {'ROI Trace', 'Mean Trace'});
             legend(handles.handles.KymoCurveAxes(2), leghandles2, {'ROI Trace', 'Mean Trace'});
             
         else
             set(handles.handles.KymoCurveFig, 'Position', [postNow(1:2), 500, 400]);
             handles.handles.KymoCurveAxes(1) = axes('Parent', handles.handles.KymoCurveFig, ...
                 'Units', 'normalized', 'Position', [.055 .1 .86 .83]);


             handles.handles.KymoCurveSingles{1} = plot(handles.handles.KymoCurveAxes(1), data_green', 'Color', [.82 .74 1]);

             set(handles.handles.KymoCurveAxes, 'NextPlot', 'add');

             handles.handles.KymoCurveMean{1} = plot(handles.handles.KymoCurveAxes(1), mean(data_green, 1), 'Color', [.4 .2 1], 'LineWidth', 2);

             set(handles.handles.KymoCurveAxes, 'NextPlot', 'replace');

             set(handles.handles.KymoCurveAxes, 'XLim', [1 handles.N_frames]);
             xlabel(handles.handles.KymoCurveAxes(1), 'Time (frame)');
             ylabel(handles.handles.KymoCurveAxes(1), ROI_strings{ROI_field_num});
                title(handles.handles.KymoCurveAxes(1), 'Channel 1');
                
                leghandles1 = [handles.handles.KymoCurveSingles{1}(1); handles.handles.KymoCurveMean{1}(1)];
                legend(handles.handles.KymoCurveAxes(1), leghandles1, {'ROI Trace', 'Mean Trace'});

         end
             
         guidata(handles.handles.fig1, handles);
         
     end
     
 end

%%%%%%%%%%%%%%%%%%%%%%
% Highlight selected ROI

    function ROI_highlight(varargin)

        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        ROI_num = str2double(get(handles.handles.ROI_plot_slide_box, 'String'));
        
        if sum(ishandle(handles.handles.ROI_line))
            ROI_box_hands = handles.handles.ROI_box{ROI_num};
            ROI_line_hands = handles.handles.ROI_line(ROI_num);


            % Reset purple patches, lines to green 

            purple_lines = findobj('Parent', handles.handles.ax_ROI, 'Color', handles.handles.Selected_ROI_color);
            purple_patches = findobj('Parent', handles.handles.ax_ROI, 'FaceColor', handles.handles.Selected_ROI_color);
            set(purple_lines, 'Color', handles.handles.Unselected_ROI_color);
            %set(purple_patches, 'FaceColor', handles.handles.Unselected_ROI_color);%FaceAlpha
            set(purple_patches, 'Color', handles.handles.Unselected_ROI_color);

            % Change proper green box to purple

%             set(ROI_box_hands, 'FaceColor', handles.handles.Selected_ROI_color); %FaceAlpha
            set(ROI_box_hands, 'Color', handles.handles.Selected_ROI_color);
            set(ROI_line_hands, 'Color', handles.handles.Selected_ROI_color);
        end

        
    end

%%%%%%%%%%%%%%%%%%%%%%
% Click on ROI in Editor Window

    function PatchClickFcn(hObj, ~, ~)

        handles = guidata(findobj('Tag', 'TIFF viewer'));
        click_handle = hObj;
        ROI_num = find(cell2mat(handles.handles.ROI_box) == click_handle);
        
        clickType = get(handles.handles.fig3, 'SelectionType');
        
        if strcmp(clickType, 'normal')
        
            if ~isempty(findobj('Tag', 'GALAH_plot_display'))

                set(handles.handles.ROI_plot_slide_box, 'String', num2str(ROI_num));
                set(handles.handles.ROI_plot_slide_hand, 'Value', ROI_num);

                label_string = get(handles.handles.ROI_plot_popup, 'String');
                xlabel(handles.handles.ax_ROI_plot, 'Distance (pixels)');
                ylabel(handles.handles.ax_ROI_plot, label_string{get(handles.handles.ROI_plot_popup,'Value')});


                ROI_highlight;


            end
            
        else
                        
            WhichROI = find(cell2mat(handles.handles.ROI_box) == hObj);

            %disp(WhichROI);
            delete(handles.handles.ROI_line(WhichROI));
            delete(handles.handles.ROI_box{WhichROI});
            
            handles.handles.ROI_line(WhichROI) = [];
            handles.handles.ROI_box(WhichROI) = [];
            handles.ROI_output(WhichROI) = [];
            handles.ROI_parameters(WhichROI) = [];
            
            guidata(handles.handles.fig1, handles);
            
            if (~isempty(handles.ROI_parameters))
                if (ROI_num > size(handles.ROI_parameters, 2))  
                    ROI_num = size(handles.ROI_parameters, 2);
                end
                
                if ~isempty(findobj('Tag', 'GALAH_plot_display'))
                    
                    ROI_plot_step = 1/(length(handles.ROI_output) - 1);

                    set(handles.handles.ROI_plot_slide_box, 'String', num2str(ROI_num));
                    set(handles.handles.ROI_plot_slide_hand, 'Value', ROI_num, ...
                        'SliderStep', [ROI_plot_step ROI_plot_step], 'Max', length(handles.ROI_output));

                    label_string = get(handles.handles.ROI_plot_popup, 'String');
                    xlabel(handles.handles.ax_ROI_plot, 'Distance (pixels)');
                    ylabel(handles.handles.ax_ROI_plot, label_string{get(handles.handles.ROI_plot_popup,'Value')});

                    ROI_highlight;
                end
 
            else
                EnableDisablePlot('off')
            end

        end
            
            

    end




%%%%%%%%%%%%%%%%%%%%%%
% Enable or Disable Plots window

    function EnableDisablePlot(on_or_off)

        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        if ishandle((findobj('Tag', 'GALAH_plot_display')))
            
                set(findobj('Parent', handles.handles.fig8_img_panel), 'Visible', on_or_off);
                set(findobj('Parent', handles.handles.ax_ROI_plot), 'Visible', on_or_off);
                set(findobj('Parent', handles.handles.fig8_button_panel), 'Enable', on_or_off);
                
                set(handles.handles.ROI_Plots_launch_hand, 'Enable', 'off');
                
        end
       
        
        if strcmp(on_or_off, 'on');
        
            FillInPlotWindow;
            set(handles.handles.ROI_Plots_launch_hand, 'Enable', 'on');
            
        end

        
    end

        

%%%%%%%%%%%%%%%%%%%%%%
% Format and export data to .mat file

    function Export_data(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        GALAH_Export = [];
        
        if ~isempty(handles.Load_file);
            
            GALAH_Export.File_properties.Filename = handles.Load_file;
            %GALAH_Export.File_properties.Img_stack = handles.Img_stack;
            GALAH_Export.File_properties.Primary_channel = handles.Primary_channel;
            
        end
        

        
        if ~isempty(handles.ROI_boundaries)
            
            for k = 1:length(handles.ROI_output)
                
                Measure = [];
                
                GALAH_Export.ROI(k).Settings = handles.ROI_parameters(k);
                GALAH_Export.ROI(k).Boundary = handles.ROI_boundaries.B{k};
                
                L_vals = unique(handles.ROI_boundaries.L);
                L_vals(L_vals == 0) = [];
                GALAH_Export.ROI(k).Mask = (handles.ROI_boundaries.L == L_vals(k));
                
                x_post = handles.ROI_output{k}.x_center;
                y_post = handles.ROI_output{k}.y_center;
                domain = sqrt((x_post-circshift(x_post, 1)).^2 + (y_post-circshift(y_post, 1)).^2);
                domain(1) = 0;
                domain = cumsum(domain);
                
                Measure.Distance = domain;
                Measure.X_center = handles.ROI_output{k}(1).x_center;
                Measure.Y_center = handles.ROI_output{k}(1).y_center;
                Measure.Area = handles.ROI_output{k}(1).areas;
                
                if handles.N_channels == 1
                    
                    Measure.Mean = horzcat(handles.ROI_output{k}.Mean_mask);
                    Measure.StdDev = horzcat(handles.ROI_output{k}.std_mask);
                    
                elseif handles.N_channels == 2
                    
                    Measure.Mean(:,:,1) = horzcat(handles.ROI_output{k}(1,:).Mean_mask);
                    Measure.Mean(:,:,2) = horzcat(handles.ROI_output{k}(2,:).Mean_mask);
                    Measure.StdDev(:,:,1) = horzcat(handles.ROI_output{k}(1,:).std_mask);
                    Measure.StdDev(:,:,2) = horzcat(handles.ROI_output{k}(2,:).std_mask);
                    
                end

                
                GALAH_Export.ROI(k).Measure = Measure;
                
            end
            
        end

        if ~isempty(GALAH_Export)
            
            set(get(handles.handles.slider_panel, 'Children'), 'Enable', 'off');
            set(handles.handles.fig1, 'Pointer', 'watch');
            drawnow;
            uisave('GALAH_Export', 'Results.mat');
            set(get(handles.handles.slider_panel, 'Children'), 'Enable', 'on');
            set(handles.handles.fig1, 'Pointer', 'arrow');
            drawnow;
            
            
        end
        
        
    end


%%%%%%%%%%%%%%%%%%%%%%
% ROI Calc function

    function Calc_some_ROIs(varargin)
        tic;
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
                               
            % Generate ROI structure
            
            %disp('doing calc')
            
        chanHere = cell(1,handles.N_channels);
        
        for chan = 1:handles.N_channels;
                    
            chanHere{chan} = reshape(handles.Img_stack(:,:,:,chan), [size(handles.Img_stack, 1)*size(handles.Img_stack, 2), handles.N_frames]);
        end
            
        threshold_val = get(handles.handles.ROI_thresh_slide_hand, 'Value');
        frame_here = str2double(get(handles.handles.ROI_frame_slide_box, 'String'));
        
        if get(handles.handles.ROI_frame_UseAverage, 'Value') == 1  
            image_here = handles.MeanImg(:,:,handles.Primary_channel);
        elseif get(handles.handles.ROI_frame_UseMaxProjection, 'Value') == 1
            image_here = handles.MaxImg(:,:,handles.Primary_channel);
        else
            image_here = handles.Img_stack(:,:, frame_here, handles.Primary_channel);
        end
        
        
        if handles.Primary_channel == 1;
            min_max_here = handles.Min_max_left;
%             display_here = handles.Display_range_left;
        elseif handles.Primary_channel == 2;
            min_max_here = handles.Min_max_right;
%             display_here = handles.Display_range_right;
        end
        
        min_max_here = [min(min_max_here(:,1)), max(min_max_here(:,2))];
        
        scale_frame = (image_here - min_max_here(1))./(min_max_here(2) - min_max_here(1));

        [sizx, sizy] = size(image_here);
        
        bw = im2bw(scale_frame, threshold_val);  
        
        % Incorporate mask values
        bw = bw & handles.ROI_Mask;
        
        bw = bwmorph(bw, 'clean');
        
        dilate_radius = str2double(get(handles.handles.Dilate_rad_box, 'String'));

        Dilate_kernel = strel('disk', dilate_radius);
        dil = imdilate(bw, Dilate_kernel);


        sk = bwmorph(dil, 'thin', Inf);
        sk = bwmorph(sk, 'spur');
        sk = bwmorph(sk, 'clean');
        %bd = bwmorph(dil, 'remove');

        % Check for branches and for loops.  If either exists, run them through
        % removal function
        
        handles.binary_mask_initial = sk;
        

        
        handles.binary_mask_final = CleanUpSk(sk);
        
        % Shoot line from end of sk to center of caps on dil
        
        outLines = bwmorph(dil, 'remove');
        lineEnds = bwmorph(sk, 'endpoints');
        
        % Remove lineEnds that might be next to one another.  This cause
        % some serious issues later.  
        lineEndProps = regionprops(lineEnds, 'Area', 'PixelIdxList');
        lineRegSize = vertcat(lineEndProps.Area);
        for k = 1:numel(lineRegSize)
            if lineRegSize(k) > 1
                lineEnds(lineEndProps(k).PixelIdxList) = 0;
            end
        end

        % If there is a dilated region without a lineEnds inside of it,
        % clear the thing out and don't consider it further.
        
        labelRegs = bwlabel(dil, 4);
        lineEndsCovered = unique(labelRegs(lineEnds));
        notCovered = ismember(labelRegs, lineEndsCovered);
        
        lineEnds = lineEnds & notCovered;
        outLines = outLines & notCovered;
        sk = sk & notCovered;
        
        % Find caps to ROIs
        
        lineMinusEnds = sk - lineEnds;
        epDist = bwdist(lineEnds);
        minusDist = bwdist(lineMinusEnds);
        caps = outLines & (epDist < minusDist);
        capPoints = bwmorph(caps, 'shrink', Inf);

        [capPostX, capPostY] = find(capPoints);
        [lineEndX, lineEndY] = find(lineEnds);
        
        % Connect caps to ends to extend out ROI
        lineList = zeros(numel(lineEndX), 1);
        for m = 1:numel(lineEndX)
            lineList(m) = labelRegs(lineEndX(m), lineEndY(m));
        end
        
        for k = 1:sum(capPoints(:))
            
            % Find nearest end point
            RegHere = labelRegs(capPostX(k), capPostY(k));
            
            possibleEndY = lineEndY(lineList == RegHere);
            possibleEndX = lineEndX(lineList == RegHere);
            capEndDist = sqrt((possibleEndX - capPostX(k)).^2 + (possibleEndY - capPostY(k)).^2);
            whichEnd = find(capEndDist(:) == min(capEndDist(:)));
            
            bw1 = false(size(handles.Img_stack, 1), size(handles.Img_stack, 2));
            bw2 = bw1;
            
            bw1(capPostX(k), capPostY(k)) = 1;
            bw2(possibleEndX(whichEnd), possibleEndY(whichEnd)) = 1;
            
            D1 = bwdist(bw1, 'quasi-euclidean');
            D2 = bwdist(bw2, 'quasi-euclidean');
            D = D1 + D2;
            D = round(D * 32) / 32;

            paths = imregionalmin(D);

            paths_thinned_many = bwmorph(paths, 'thin', inf);
            
            sk = sk | paths_thinned_many;
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gotta repeat the sk opening/closing from before.
        sk = CleanUpSk(sk);
        
        %
        CCfin = bwconncomp(sk);
        handles.ROI_parameters = [];

        for k = 1:CCfin.NumObjects

            % Figure pair-pair distances along route
            sorted = zeros(numel(CCfin.PixelIdxList{k}), 2);
            distSoFar = zeros(numel(CCfin.PixelIdxList{k}), 1);
            end_find = false(CCfin.ImageSize(1), CCfin.ImageSize(2));
            end_find(CCfin.PixelIdxList{k}) = 1;
            end_points = bwmorph(end_find, 'endpoints');
            [endy, endx] = find(end_points);

            kx = [-1 0 1; -1 0 1; -1 0 1];
            ky = [-1 -1 -1; 0 0 0; 1 1 1];
            kdist = [sqrt(2) 1 sqrt(2); 1 0 1; sqrt(2) 1 sqrt(2)];

            % Pick and end point and start marching along
            end_start = [endx(1) endy(1)];
            %unsorted(ismember(unsorted, end_start, 'rows'), :) = [];
            sorted(1,:) = end_start;
%             point_now = end_start;
            sorted_num = 2;

            while sum(sum(end_find)) > 0

                 end_find(sorted((sorted_num-1),2), sorted((sorted_num-1),1)) = 0;
                 
                 if sum(end_find(:) > 0)

                    % Pull out local matrix to current point

                    kernx1 = max([1 (sorted((sorted_num-1),2)-1)]);
                    kernx2 = min([sizx (sorted((sorted_num-1),2)+1)]);
                    kerny1 = max([1 (sorted((sorted_num-1),1)-1)]);
                    kerny2 = min([sizy (sorted((sorted_num-1),1)+1)]);

                    kern = end_find((kernx1:kernx2), (kerny1:kerny2));

                    % Special case padding when the kernel is against a border
                    if (sorted((sorted_num-1),2)-1) < 1

                       kern = [0 0 0; kern]; % Pad top

                    end

                    if (sorted((sorted_num-1),2)+1) > sizx;

                        kern = [kern; 0 0 0]; % Pad bottom

                    end

                    if (sorted((sorted_num-1),1)-1) < 1

                       kern = [0 0 0; kern']'; % Pad left

                    end

                    if (sorted((sorted_num-1),1)-1) > sizy

                       kern = [kern'; 0 0 0]'; % Pad right

                    end



                    ka = sum(sum(kern.*kx));
                    kb = sum(sum(kern.*ky));
                    kd = sum(sum(kern.*kdist));

                    sorted(sorted_num, :) = sorted((sorted_num-1), :) + [ka kb];
                    distSoFar(sorted_num) = kd; 
                    
                 end

                sorted_num = sorted_num + 1;

            end

            CCfin.dist{k} = distSoFar;
            CCfin.cumdist{k} = cumsum(distSoFar);
            CCfin.SortedPixelList{k} = sorted; 
            
            handles.ROI_parameters(k).End_caps = handles.End_caps;
            handles.ROI_parameters(k).ROI_length = handles.ROI_length;
            handles.ROI_parameters(k).ROI_width = handles.ROI_width;
            handles.ROI_parameters(k).Dilate_radius = handles.Dilate_radius;
            handles.ROI_parameters(k).Threshold_absolute = str2double(get(handles.handles.ROI_thresh_slide_box, 'String'));
            handles.ROI_parameters(k).Threshold_relative = (get(handles.handles.ROI_thresh_slide_hand, 'Value'));

        end
        
        ROI_line = zeros(size(CCfin.PixelIdxList, 2), 1);
        ROI_color = zeros(size(CCfin.PixelIdxList,2), 3);
        ROI_box = cell(size(CCfin.PixelIdxList, 2), 1);
        forget_list = [];
        
        set(handles.handles.ax_ROI, 'NextPlot', 'add');
        
        % Clear out the borders, old boxes
        delete(findobj('Parent', handles.handles.ax_ROI, 'Type', 'line'));

        handles.ROI_output = {};

        for k = 1:size(CCfin.PixelIdxList, 2)
            
                       ptraw = CCfin.SortedPixelList{k};            
            pt = interparc(round(max(CCfin.cumdist{k})/handles.ROI_length), ptraw(:,1), ptraw(:,2) , 'linear');

            CCfin.InterpPoints{k} = reshape(pt(~isnan(pt)), [], 2);


            points_here = CCfin.InterpPoints{k};
 

            if size(points_here, 1) > 1
                
                ROIstructHold = ROI_area_traceGALAH(handles.Img_stack(:,:,1,:), 1:handles.N_channels,...
                    points_here, handles.ROI_width, length(points_here), handles.ROI_smoothing, handles.End_caps, 0);
                
                ROI_line(k) = plot(handles.handles.ax_ROI, ROIstructHold(1).points(:,1), ROIstructHold(1).points(:,2), ':', 'Color', ...
                    handles.handles.Unselected_ROI_color);
                
                Segs_tot = ROIstructHold.Segs_tot;

%%%%%%%%%%%%%%%%%%
% This section rewritten 070715 for increased speed
                
%                 mask = false(size(handles.Img_stack, 1), size(handles.Img_stack, 2), (size(Segs_tot, 1) - 1));
                PixelIdxHold = cell(1, (size(Segs_tot, 1) - 1));
                for pol = 1:(size(Segs_tot, 1) - 1);

                   % disp(pol)
                        x_mask = [Segs_tot(pol, 4) Segs_tot(pol+1, 4) Segs_tot(pol+1, 6) Segs_tot(pol, 6) Segs_tot(pol, 4)];
                        y_mask = [Segs_tot(pol, 5) Segs_tot(pol+1, 5) Segs_tot(pol+1, 7) Segs_tot(pol, 7) Segs_tot(pol, 5)];

                        maskHere = logical(poly2mask(x_mask, y_mask, size(handles.Img_stack, 1), size(handles.Img_stack, 2)));
                        
                        PixelIdxHold{pol} = find(maskHere);
                         
                end
                
                FakeCC.Connectivity = 8;
                FakeCC.ImageSize = [size(handles.Img_stack, 1), size(handles.Img_stack, 2)];
                FakeCC.NumObjects = (size(Segs_tot, 1) - 1);
                FakeCC.PixelIdxList = PixelIdxHold;
                
                onesVect = ones(handles.N_frames, 1);
                
                structStart = struct('x_center', ROIstructHold(1).x_center, ...
                                     'y_center', ROIstructHold(1).y_center, ...
                                     'center_dist', ROIstructHold(1).center_dist, ...
                                     'center_spacing', ROIstructHold(1).center_spacing, ... 
                                     'areas', ROIstructHold(1).areas, ... 
                                     'Mean_mask', onesVect, ... 
                                     'std_mask', onesVect, ...
                                     'Segs_tot', Segs_tot);
                                 
                ROI_struct = repmat(structStart, [handles.N_channels, handles.N_frames]);
                
                for chan = 1:handles.N_channels;
                    
%                     chanHere = reshape(handles.Img_stack(:,:,:,chan), [size(handles.Img_stack, 1)*size(handles.Img_stack, 2), handles.N_frames]);
                    
                    mean_mask = zeros((size(Segs_tot, 1) - 1), handles.N_frames);
                    std_mask = zeros((size(Segs_tot, 1) - 1), handles.N_frames);
                    
                    %                     for fk = 1:handles.N_frames
                    
                    %                         regStruct = regionprops(FakeCC, handles.Img_stack(:,:,fk,chan), 'MeanIntensity', 'PixelValues');
                    
                    for pol = 1:(size(Segs_tot, 1) - 1)
                        
                        %                             frameHere = handles.Img_stack(:,:,fk, chan);
                        %                             frameHere = frameHere(:);
                        %                             mskHere = mask(:,:,pol);
                        
                        mean_mask(pol, :) = mean(chanHere{chan}(PixelIdxHold{pol}, :), 1);
                        std_mask(pol, :) = std(chanHere{chan}(PixelIdxHold{pol}, :), [], 1);
                        %                         polHold(polHold == 0) = NaN;
                        %
                        %                         mean_mask = nanmean(reshape(polHold, [size(handles.Img_stack, 1)*size(handles.Img_stack, 2), handles.N_frames]));
                        
                    end
                    
                    
                    
                    for fk = 1:handles.N_frames
                        ROI_struct(chan, fk).Mean_mask = mean_mask(:, fk);
                        ROI_struct(chan, fk).std_mask = std_mask(:, fk);                   
                    end
                    
                    
                end
                
% End 070715 edit section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %box_handles = zeros((size(Segs_tot, 1) - 1), 1);

%                 for pol = 1:(size(Segs_tot, 1) - 1);
% 
%                     x_mask = [Segs_tot(pol, 4) Segs_tot(pol+1, 4) Segs_tot(pol+1, 6) Segs_tot(pol, 6) Segs_tot(pol, 4)];
%                     y_mask = [Segs_tot(pol, 5) Segs_tot(pol+1, 5) Segs_tot(pol+1, 7) Segs_tot(pol, 7) Segs_tot(pol, 5)];
% 
% 
%                     box_handles(pol) = plot(handles.handles.ax_ROI, x_mask', y_mask', 'Color', ...
%                         handles.handles.Unselected_ROI_color);
% 
%                 end

                xx = [Segs_tot(:,4); flipud(Segs_tot(:,6)); Segs_tot(1,4)];
                yy = [Segs_tot(:,5); flipud(Segs_tot(:,7)); Segs_tot(1,5)];
                
                %box_handles = patch(xx, yy, handles.handles.Unselected_ROI_color, 'EdgeColor', 'none');%FaceAlpha
                box_handles = plot(xx, yy, 'Color', handles.handles.Unselected_ROI_color);
                %set(box_handles, 'Color', handles.handles.Unselected_ROI_color);
                
                ROI_box{k} = box_handles;
                if k > numel(handles.handles.Border_randomize)
                	ROI_color(k,:) = handles.handles.border_colors(handles.handles.Border_randomize(mod(k, nume(handles.handles.Border_randomize))), :);
                else
                    ROI_color(k,:) = handles.handles.border_colors(handles.handles.Border_randomize(k), :);
                end
                handles.ROI_output{k} = ROI_struct;
                
            else
                
                forget_list = [forget_list; k];
                handles.ROI_output{k} = [];

            end

        end
        
%         toc
        
        handles.ROI_parameters(cellfun(@isempty, handles.ROI_output)) = [];
        handles.handles.ROI_line = ROI_line(~cellfun(@isempty, handles.ROI_output));
        handles.handles.ROI_color = ROI_color(~cellfun(@isempty, handles.ROI_output), :);
        handles.handles.ROI_box = ROI_box(~cellfun(@isempty, handles.ROI_output)); % Not sure if needed, best way to do this.
        handles.ROI_output(cellfun(@isempty, handles.ROI_output)) = [];

        handles.ROIsAreSet = true;
        set(handles.handles.ROI_Plots_launch_hand, 'Enable', 'on');
        
        set(handles.handles.ax_ROI, 'NextPlot', 'replace');
        set(cell2mat(handles.handles.ROI_box), 'ButtonDownFcn', @PatchClickFcn);


            set(handles.handles.fig3, 'UserData', 'ROI_calced');
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            clear chanHere
            
    end


%%%%%%%%%%%%%%%%%%%%%%
% Actions on GUI close

    function GUI_close_fcn(varargin)
        
        close(findobj('Tag', 'GALAH_Image_prefs'));
        close(findobj('Tag', 'GALAH_ROI_review'));
        close(findobj('Tag', 'GALAH_ROI_calc'));
        close(findobj('Tag', 'GALAH_FilterBkgd'));
        close(findobj('Tag', 'GALAH_ROI_calc'));
        close(findobj('Tag', 'GALAH_Regions'));
        close(findobj('Tag', 'GALAH_ExpRegions'));
        close(findobj('Tag', 'GALAH_ExpFit'));
        close(findobj('Tag', 'GALAH_ExpRegions'));
        close(findobj('Tag', 'GALAH_AlignChannels'));
        close(findobj('Tag', 'GALAH_DriftCorrect'));
        
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Align_launch(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
   
        % Launch window to enter type, size of flatten object
        TransformOptions = {'translation', 'rigid', 'similarity', 'affine', 'custom'};
        WhichOption = find(strcmp(TransformOptions, handles.TransformType));
        
        if handles.TransformRefFrameCustom == 1
            startFrame = handles.TransformRefFrame;
        else
            startFrame = calcStartFrame(handles.TransformFixedChannel);
            handles.TransformRefFrame = startFrame;
        end
        
        TransformFixed = {'1', '2'};

        if ~isempty(findobj('Tag', 'GALAH_AlignChannels'))

            uistack(findobj('Tag', 'GALAH_AlignChannels'), 'top');

        else

            mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
            fig14_size = [300 100];
            fig14_position = [(mf_post(1) + (mf_post(3) - fig14_size(1))/2) (mf_post(2) + (mf_post(4) - fig14_size(2))/2)];
            handles.handles.fig14 = figure('Name','Align Channels', 'Tag', 'GALAH_AlignChannels', 'Units', 'pixels',...
                'Position',[fig14_position fig14_size], 'NumberTitle', 'off', 'Toolbar', 'none', 'Menu', 'none');
            set(handles.handles.fig14, 'Color',[0.9 0.9 0.9]);

            handles.handles.Align_pulldown = uicontrol(handles.handles.fig14, 'Style', 'popupmenu', 'Units', 'normalized', ...
                    'Position', [0.3100 0.8100 0.2567 0.1400], 'String', TransformOptions, 'Value', WhichOption, 'Callback', @Align_pulldown_call);

            handles.handles.Align_text = uicontrol(handles.handles.fig14, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.0217 0.600 0.2550 0.4000], 'String', 'Transformation Type:', 'BackgroundColor', [0.9 0.9 0.9]);
            
            handles.handles.Align_Fixed_pulldown = uicontrol(handles.handles.fig14, 'Style', 'popupmenu', 'Units', 'normalized', ...
                    'Position', [0.3100 0.4800 0.11 0.1400], 'String', TransformFixed, 'Value', handles.TransformFixedChannel, ...
                    'Callback', @Align_FixedSelection_call);
            
            handles.handles.Align_FixedChan_text = uicontrol(handles.handles.fig14, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.0117 .190 0.2550 0.4000], 'String', 'Fixed Channel:', 'BackgroundColor', [0.9 0.9 0.9]);
            
            handles.handles.AlignReferenceFrame = uicontrol(handles.handles.fig14, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', [0.3100 0.0400 0.11 0.2300], 'String', num2str(startFrame), ...
                    'Callback', @Align_ReferenceFrame_call);
            
            handles.handles.Align_ReferenceFrame_text = uicontrol(handles.handles.fig14, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.0117 -.080 0.2550 0.4000], 'String', 'Reference Frame:', 'BackgroundColor', [0.9 0.9 0.9]);
% 
            handles.handles.ClearTransButton = uicontrol(handles.handles.fig14, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [.60 .68 .38 .28], 'String', 'Reset to Original', 'Callback', @Reset_trans_push, ...
                'TooltipString', 'Clear all applied image transformations');
% 
            handles.handles.AlignChannelsButton = uicontrol(handles.handles.fig14, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [.6 .35 .38 .28], 'String', 'Align Channels', 'Callback', @ApplyTrans_push, ...
                'TooltipString', 'Align Channels via set parameters');
% 
            handles.handles.SetandCloseTransButton = uicontrol(handles.handles.fig14, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [.6 .02 .38 .28], 'String', 'Save Aligned Stack', 'Callback', @SetTransPush, ...
                'TooltipString', 'Save aligned image stack');


            guidata(findobj('Tag', 'TIFF viewer'), handles);

        end
        
        function Align_pulldown_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.TransformType = TransformOptions{get(handles.handles.Align_pulldown, 'Value')};
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            if strcmp(handles.TransformType, 'custom')
                % Launch custom transform input window
                
                if ~isempty(findobj('Tag', 'GALAH_CustomTransformation'))

                    uistack(findobj('Tag', 'GALAH_CustomTransformation'), 'top');

                else
                
                
                    mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
                    fig15_size = [300 250];
                    fig15_position = [(mf_post(1) + (mf_post(3) - fig15_size(1))/2) (mf_post(2) + (mf_post(4) - fig15_size(2))/2)];
                    handles.handles.fig15 = figure('Name','Custom Transformation', 'Tag', 'GALAH_CustomTransformation', 'Units', 'pixels',...
                        'Position',[fig15_position fig15_size], 'NumberTitle', 'off', 'Toolbar', 'none', 'Menu', 'none');
                    set(handles.handles.fig15, 'Color',[0.9 0.9 0.9]);

                    startCorner = [.15 .7];
                    tAspacing = [.05 .05];
                    tABoxSize = [.2 .1];

                    transBoxArrayDims = [3 3];
                    
                    handles.handles.CustomMatrixText = uicontrol(handles.handles.fig15, 'Style', 'text', 'Units', 'normalized', ...
                        'Position', [0.0267    0.8633    0.4100    0.0667], 'String', 'Transformation Matrix:', 'BackgroundColor', [0.9 0.9 0.9]);

                    for cx = 1:transBoxArrayDims(1)
                        for cy = 1:transBoxArrayDims(2)

                            handles.handles.CustomTransBox(cy, cx) = uicontrol(handles.handles.fig15, 'Style', 'edit', 'Units', 'normalized', ...
                                'Position', [startCorner(1)+(cx-1)*tAspacing(1)+(cx-1)*tABoxSize(1) ...
                                             startCorner(2)-(cy-1)*tAspacing(2)-(cy-1)*tABoxSize(2) ...
                                             tABoxSize(1) tABoxSize(2)], ...
                                'String', sprintf('%.4f', handles.AffineMatrix(cy, cx)), ...
                                'Callback', @CustomBoxEditCall);
                        end
                    end
                    
                    
                    handles.handles.ImportTransMatrix = uicontrol(handles.handles.fig15, 'Style', 'pushbutton', 'Units', 'normalized', ...
                        'Position', [0.0767 0.0467 0.3867 0.1500], 'String', 'Import Transform', 'Callback', @ImportTransMatrixPush, ...
                        'TooltipString', 'Import Transformation Matrix');
% 
                    handles.handles.CustomTransOK = uicontrol(handles.handles.fig15, 'Style', 'pushbutton', 'Units', 'normalized', ...
                        'Position', [0.56 0.0467 0.3867 0.1500], 'String', 'OK', 'Callback', @CustomTransOKPush, ...
                        'TooltipString', 'Set and Close');
                    
                    
                 

                end
                

                
            end
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
        function CustomBoxEditCall(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            isEdited = (handles.handles.CustomTransBox == varargin{1});
            newVal = get(handles.handles.CustomTransBox(isEdited), 'String');
            if all(isstrprop(newVal, 'digit') | newVal == '.')
                handles.AffineMatrix(isEdited) = str2double(newVal);
                set(handles.handles.CustomTransBox(isEdited), 'String', ...
                    sprintf('%.4f', str2double(newVal)));
            else
                set(handles.handles.CustomTransBox(isEdited), 'String', ...
                    sprintf('%.4f', handles.AffineMatrix(isEdited)));
            end
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
                
        end
        
        function ImportTransMatrixPush(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            [~, LoadString] = fileparts(handles.Load_file);
            LoadString = strcat(LoadString, '_aligned');

            % Load transform matrix file
            [filename, pathname, filterspec] = uigetfile(...
                {'*.dat', 'Transform Matrix (*.dat)'; ...
                '*.txt', 'ASCII Text File (*.txt)';}, 'Select Transform Matrix File', LoadString);
            
            if filterspec == 1 || filterspec == 2
                % Previously-saved transform matrix file.  Assume any
                % non-digit character that starts a line indicates a
                % comment
                Infile = fopen(fullfile(pathname, filename), 'r');
                dataIn = zeros(3,3);
                linesIn = 1;
                while true
                    dataMaybe = fgetl(Infile);
                    if dataMaybe == -1
                        break
                    else
                    
                        if ~isempty(regexp(dataMaybe, '^[-.\d]', 'once'))
                            dataIn(linesIn,:) = cell2mat(textscan(dataMaybe, '%f %f %f', 'delimiter', '\t'));
                            linesIn = linesIn+1;
                        end
                    end
                end
                    
                fclose(Infile);
                
%                 disp(dataIn)
                
                % Add loaded-in data to the right spots
                handles.AffineMatrix = dataIn;
                handles.tformObject = affine2d(dataIn);
                
                for bxNum = 1:numel(handles.handles.CustomTransBox)
                    set(handles.handles.CustomTransBox(bxNum), 'String', sprintf('%.4f', dataIn(bxNum)));
                end
                
            end
                
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
        function CustomTransOKPush(varargin)
            
            close(findobj('Tag', 'GALAH_CustomTransformation'));
            
        end
        
        function Align_FixedSelection_call(varargin)
            
%             disp('Fixed channel selection')
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.TransformFixedChannel = get(varargin{1}, 'Value');
            
            switch handles.TransformFixedChannel
                case 1
                    handles.TransformMobileChannel = 2;
                case 2 
                    handles.TransformMobileChannel = 1;
            end
            
%             disp(handles.TransformFixedChannel)
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
        function Align_ReferenceFrame_call(varargin)
           % Set a new alignment reference frame if desired
           oldRefFrame = round(str2double(get(varargin{1}, 'String')));
           if (oldRefFrame > 0) && le(oldRefFrame, size(handles.Img_stack, 3))
               handles.TransformRefFrame = oldRefFrame;
               handles.TransformRefFrameCustom = 1;
           else
               % No change
           end
           set(handles.handles.AlignReferenceFrame, 'String', num2str(handles.TransformRefFrame));
           guidata(findobj('Tag', 'TIFF viewer'), handles);
            
            
        end
        
        function Reset_trans_push(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Img_stack = handles.ImgHold;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            
        end
        
        function ApplyTrans_push(varargin)
            
            set(findobj('Parent', handles.handles.fig14), 'Enable', 'off')
            set(handles.handles.fig14, 'Pointer', 'watch')
            set(handles.handles.fig14, 'UserData', 'Calculating')
            drawnow;

            CalcTransformation;
            waitfor(handles.handles.fig14, 'UserData', 'Finished')
            set(findobj('Parent', handles.handles.fig14), 'Enable', 'on')
            set(handles.handles.fig14, 'Pointer', 'arrow')
            set(handles.handles.fig14, 'UserData', [])
            drawnow;

            Display_images_in_axes;
            
        end
        
        function SetTransPush(varargin)
            
            handles = guidata(findobj('tag', 'TIFF viewer'));
            
%             set(findobj('Parent', handles.handles.fig14), 'Enable', 'off')
%             set(handles.handles.fig14, 'Pointer', 'watch')
%             set(handles.handles.fig14, 'UserData', 'Calculating')
%             drawnow;
%             
%             CalcTransformation;
%             waitfor(handles.handles.fig14, 'UserData', 'Finished')
%             set(findobj('Parent', handles.handles.fig14), 'Enable', 'on')
%             set(handles.handles.fig14, 'Pointer', 'arrow')
%             set(handles.handles.fig14, 'UserData', [])
%             drawnow;

%             handles = guidata(findobj('tag', 'TIFF viewer'));
            
            [~, SaveString] = fileparts(handles.Load_file);
            SaveString = strcat(SaveString, '_aligned');
            
            % Save as a 16-bit TIFF file
            [filename, pathname, filterspec] = uiputfile(...
                {'*.tif', '8-bit TIFF File (*.tif)'; ...
                '*.tif', '16-bit TIFF File (*.tif)'; ...
                '*.mat', 'MATLAB Data File (*.mat)'; ...
                '*.dat', 'Transform Matrix (*.dat)'}, 'Select Save Filename', SaveString);

            if (filterspec > 0) && (filterspec < 3)

                fnamehere = filename;

                if filterspec == 1;
                    ImgOut = uint8(255*handles.Img_stack/max(handles.Img_stack(:)));
                elseif filterspec == 2;
                    ImgOut = uint16(65535*handles.Img_stack/max(handles.Img_stack(:)));
                end

                %oldObj = Tiff(handles.Load_file);
                tiffObj = Tiff(fullfile(pathname, fnamehere), 'w');

                clear tagstruct;
                tagstruct.ImageLength = size(ImgOut,1);
                tagstruct.ImageWidth = size(ImgOut,2);
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample = 8*filterspec;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.RowsPerStrip = 16;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Software = 'MATLAB';
%                 tagstruct.SubIFD = handles.N_channels*size(handles.Img_stack, 3);
    %             
    %             tagstruct.Orientation =  Tiff.Orientation.TopLeft;
    %             tagstruct.ExtraSamples =  Tiff.ExtraSamples.AssociatedAlpha;
                tagstruct.Compression = Tiff.Compression.None;
                tiffObj.setTag(tagstruct);

%                 tagstruct = rmfield(tagstruct, 'SubIFD');
                w = waitbar(0, 'Saving Image Stack');
                


                for j = 1:size(ImgOut, 3)

                    if j > 1
                        tiffObj.writeDirectory();
                        tiffObj.setTag(tagstruct);
                    end
                    tiffObj.write(ImgOut(:,:,j,1));
                    if handles.N_channels == 2
                        tiffObj.writeDirectory();
                        tiffObj.setTag(tagstruct);
                        tiffObj.write(ImgOut(:,:,j,2));
                    end

                    waitbar(j/size(handles.Img_stack, 3), w);
                end

                tiffObj.close;
                close(w);
                
            elseif filterspec == 3
                %.mat file
                save(fullfile(pathname, filename), '-struct', 'handles', 'Img_stack');
            elseif filterspec == 4
                %.dat transformation matrix
                % Outputs affine matrix used for registration as a
                % tab-delimited text file.  Can be loaded when using
                % 'custom' option
                tformFile = fopen(fullfile(pathname, filename), 'w');
                fprintf(tformFile, '# Transform object saved %s\r\n', datestr(now));
                fprintf(tformFile, '# Transform type %s\r\n', handles.TransformType);
                fprintf(tformFile, '# Fixed channel %.f\r\n', handles.TransformFixedChannel);
                fprintf(tformFile, '# Reference frame %.f\r\n', handles.TransformRefFrame);
                fprintf(tformFile, '%.5f\t%.5f\t%.5f\r\n', handles.AffineMatrix(1,:));
                fprintf(tformFile, '%.5f\t%.5f\t%.5f\r\n', handles.AffineMatrix(2,:));
                fprintf(tformFile, '%.5f\t%.5f\t%.5f\r\n', handles.AffineMatrix(3,:));
                fclose(tformFile);
                
            end
            
        end
        
        function CalcTransformation(varargin)
           
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
%             disp(handles.TransformFixedChannel)
%             disp(handles.TransformMobileChannel)
             
%             disp('CalcTransformation')
            [optimizer,metric] = imregconfig('monomodal');
%             optimizer.InitialRadius = optimizer.InitialRadius/3.5;
            % Loop through each frame to do the alignment
            
            if ~strcmp(handles.TransformType, 'custom')
            
                handles.tformObject = imregtform(handles.Img_stack(:,:, handles.TransformRefFrame, handles.TransformMobileChannel), ...
                    handles.Img_stack(:,:,handles.TransformRefFrame, handles.TransformFixedChannel), handles.TransformType, optimizer, metric);
            else
               % Apply custom transformation
               handles.tformObject = affine2d(handles.AffineMatrix);
               
            end
            
            w = waitbar(0, 'Aligning Channels');
            for frm = 1:size(handles.Img_stack, 3)
                handles.Img_stack(:,:,frm, handles.TransformMobileChannel) = imwarp_same(handles.Img_stack(:,:,frm, handles.TransformMobileChannel), ...
                    handles.tformObject);
                waitbar(frm/size(handles.Img_stack, 3), w);
            end
            close(w);
            
            handles.AffineMatrix = handles.tformObject.T;
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
%             Display_images_in_axes;
            
            set(handles.handles.fig14, 'UserData', 'Finished')
            
        end
        
        
    end

    function Drift_launch(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        if ~isempty(findobj('Tag', 'GALAH_DriftCorrect'))

            uistack(findobj('Tag', 'GALAH_DriftCorrect'), 'top');

        else

            DriftChannel = {'1', '2'};

            mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
            fig15_size = [300 100];
            fig15_position = [(mf_post(1) + (mf_post(3) - fig15_size(1))/2) (mf_post(2) + (mf_post(4) - fig15_size(2))/2)];
            handles.handles.fig15 = figure('Name','Drift Correction', 'Tag', 'GALAH_DriftCorrect', 'Units', 'pixels',...
                'Position',[fig15_position fig15_size], 'NumberTitle', 'off', 'Toolbar', 'none', 'Menu', 'none');
            set(handles.handles.fig15, 'Color',[0.9 0.9 0.9]);
            

            handles.handles.Drift_Fixed_pulldown = uicontrol(handles.handles.fig15, 'Style', 'popupmenu', 'Units', 'normalized', ...
                    'Position', [0.3100 0.8100 0.11 0.1400], 'String', DriftChannel, 'Value', handles.DriftRefFrame, ...
                    'Callback', @Drift_FixedSelection_call);
            
            handles.handles.Drift_FixedChan_text = uicontrol(handles.handles.fig15, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.0217 0.600 0.2550 0.4000], 'String', 'Reference Channel:', 'BackgroundColor', [0.9 0.9 0.9]);
            
            handles.handles.DriftRangeBox(1) = uicontrol(handles.handles.fig15, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', [0.2700    0.3700    0.1100    0.2300], 'String', num2str(handles.DriftFrameRange(1)), ...
                    'Callback', @DriftEditBoxCall);
  
            handles.handles.DriftStartFrameText = uicontrol(handles.handles.fig15, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.0075    0.3900    0.2425    0.1729], 'String', 'Frame Range:', 'BackgroundColor', [0.9 0.9 0.9]);
            
            handles.handles.DriftRangeBox(2) = uicontrol(handles.handles.fig15, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', [0.4667    0.3700    0.1100    0.2300], 'String', num2str(handles.DriftFrameRange(2)), ...
                    'Callback', @DriftEditBoxCall);
            
            handles.handles.DriftEndFrameText = uicontrol(handles.handles.fig15, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.3855    0.4000    0.0745    0.1565], 'String', '-', 'BackgroundColor', [0.9 0.9 0.9]);
            
            handles.handles.DriftSmoothBox = uicontrol(handles.handles.fig15, 'Style', 'edit', 'Units', 'normalized', ...
                    'Position', [0.3067    0.0600    0.1100    0.2300], 'String', num2str(handles.DriftSmoothing), ...
                    'Callback', @DriftSmoothBoxCall);
            
            handles.handles.DriftSmoothFrameText = uicontrol(handles.handles.fig15, 'Style', 'text', 'Units', 'normalized', ...
                'Position', [0.0121    0.0700    0.2345    0.1800], 'String', 'Smoothing:', 'BackgroundColor', [0.9 0.9 0.9]);
            
            handles.handles.ClearTransButton = uicontrol(handles.handles.fig15, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [.60 .68 .38 .28], 'String', 'Reset to Original', 'Callback', @Reset_Drift_push, ...
                'TooltipString', 'Clear all applied image transformations');
% 
            handles.handles.CorrectDriftButton = uicontrol(handles.handles.fig15, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [.6 .35 .38 .28], 'String', 'Apply Correction', 'Callback', @DriftCorrect_push, ...
                'TooltipString', 'Correct Drift via set parameters');
% 
            handles.handles.SetandCloseTransButton = uicontrol(handles.handles.fig15, 'Style', 'pushbutton', 'Units', 'normalized', ...
                'Position', [.6 .02 .38 .28], 'String', 'Save Aligned Stack', 'Callback', @SaveDriftPush, ...
                'TooltipString', 'Save aligned image stack');


            guidata(findobj('Tag', 'TIFF viewer'), handles);

        end
        
        function Drift_FixedSelection_call(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.DriftRefFrame = get(varargin{1}, 'Value');
            
%             disp(handles.TransformFixedChannel)
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
        function DriftEditBoxCall(varargin)
            
%             disp(varargin{1})
           % Set a new alignment reference frame if desired
           oldRefFrame = round(str2double(get(varargin{1}, 'String')));
           whichBox = find((handles.handles.DriftRangeBox  == varargin{1}));
           switch whichBox
               case 1
                   
                   if (oldRefFrame > 0) && le(oldRefFrame, size(handles.Img_stack, 3)) && ...
                           oldRefFrame < round(str2double(get(handles.handles.DriftRangeBox(2), 'String')));
                       
                       handles.DriftFrameRange(whichBox) = (oldRefFrame);
                       set(handles.handles.DriftRangeBox(whichBox), 'String', num2str(oldRefFrame));
                   else
                       % No change
                   end
                   
               case 2
                   
                    if (oldRefFrame > 0) && le(oldRefFrame, size(handles.Img_stack, 3)) && ...
                        oldRefFrame > round(str2double(get(handles.handles.DriftRangeBox(1), 'String')));
                    
                       handles.DriftFrameRange(whichBox) = (oldRefFrame);
                       set(handles.handles.DriftRangeBox(whichBox), 'String', num2str(oldRefFrame));
                   else
                       % No change
                    end
                   
           end
                   
                   guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
        function DriftSmoothBoxCall(varargin)
            
            oldRefFrame = (str2double(get(varargin{1}, 'String')));
%             disp(oldRefFrame)
            isset = 0;
                if (oldRefFrame) > 0
                    handles.DriftSmoothing = (oldRefFrame);
                    isset = 1;
                end
            
            if isset == 0
                set(varargin{1}, 'String', num2str(handles.DriftSmoothing));
            end
            
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            
        end
        
        function Reset_Drift_push(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            handles.Img_stack = handles.ImgHold;
            guidata(findobj('Tag', 'TIFF viewer'), handles);
            Display_images_in_axes;
            
        end
        
        function DriftCorrect_push(varargin)
            
            % pop out plot of drift correct with bits labeled, then make
            % push in pop-out apply correction
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            
            
            handles.DriftMatrix = zeros(handles.N_frames, 2);
            handles.DriftMatrixSmoothed = handles.DriftMatrix;
            
            [optimizer,metric] = imregconfig('monomodal');
            w = waitbar(0, 'Calculating Drift');

            for k = (handles.DriftFrameRange(1)+1):(handles.DriftFrameRange(2))
            
                
                handles.tformObject = imregtform(handles.Img_stack(:,:, k, handles.DriftRefFrame), ...
                    handles.Img_stack(:,:,k-1, handles.DriftRefFrame), 'translation', optimizer, metric);
                
                waitbar((k - handles.DriftFrameRange(1))/diff(handles.DriftFrameRange), w);
                
                handles.DriftMatrix(k, :) = handles.tformObject.T(3, 1:2);
                
            end
            close(w);
            
            handles.DriftMatrix = cumsum(handles.DriftMatrix);
            
            % Apply smoothing as Lowess filter
            for m = 1:2
                handles.DriftMatrixSmoothed((handles.DriftFrameRange(1):handles.DriftFrameRange(2)),m) = ...
                    smooth(handles.DriftMatrix((handles.DriftFrameRange(1):handles.DriftFrameRange(2)),m), ...
                    handles.DriftSmoothing, 'lowess');
            end
            
            mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
            fig16_size = [700 400];
            fig16_position = [(mf_post(1) + (mf_post(3) - fig16_size(1))/2) (mf_post(2) + (mf_post(4) - fig16_size(2))/2)];
            handles.handles.fig16 = figure('Name','Drift Display', 'Tag', 'GALAH_DriftDisplay', 'Units', 'pixels',...
                'Position',[fig16_position fig16_size], 'NumberTitle', 'off', 'Toolbar', 'figure', 'Menu', 'none');
            set(handles.handles.fig16, 'Color',[0.9 0.9 0.9]);
            
            % Axes
            
            handles.handles.DriftAxes = axes('Parent', handles.handles.fig16, 'Position', [0.0871    0.2000    0.8700    0.7500]);
            plot(handles.handles.DriftAxes, handles.DriftFrameRange(1):handles.DriftFrameRange(2),...
                handles.DriftMatrix(handles.DriftFrameRange(1):handles.DriftFrameRange(2), 1), ...
                'Color', sqrt([1 .4 .2]), 'LineWidth', 1);
            set(handles.handles.DriftAxes, 'NextPlot', 'add');
            plot(handles.handles.DriftAxes, handles.DriftFrameRange(1):handles.DriftFrameRange(2),...
                handles.DriftMatrixSmoothed(handles.DriftFrameRange(1):handles.DriftFrameRange(2), 1), ...
                'Color', [1 .4 .2], 'LineWidth', 2);
            
            plot(handles.handles.DriftAxes, handles.DriftFrameRange(1):handles.DriftFrameRange(2),...
                handles.DriftMatrix(handles.DriftFrameRange(1):handles.DriftFrameRange(2), 2), ...
                'Color', sqrt([.4 .2 1]), 'LineWidth', 1);
            plot(handles.handles.DriftAxes, handles.DriftFrameRange(1):handles.DriftFrameRange(2),...
                handles.DriftMatrixSmoothed(handles.DriftFrameRange(1):handles.DriftFrameRange(2), 2), ...
                'Color', [.4 .2 1], 'LineWidth', 2);
            set(handles.handles.DriftAxes, 'NextPlot', 'replace', 'XLim', handles.DriftFrameRange);
            xlabel(handles.handles.DriftAxes, 'Frame Number'); ylabel(handles.handles.DriftAxes, 'Drift (pixels)');
            
            % Button push
            
            handles.handles.SetAndCloseDriftCorrection = uicontrol(handles.handles.fig16, 'Style', 'pushbutton', ...
                'Units', 'normalized', 'Position', [0.8243    0.0275    0.1286    0.0950], ...
                'String', 'Apply', 'Callback', @ApplyDriftCorrection, ...
                'TooltipString', 'Apply calculated drift correction');
            
            handles.handles.CancelDriftCorrection = uicontrol(handles.handles.fig16, 'Style', 'pushbutton', ...
                'Units', 'normalized', 'Position', [0.6543    0.0275    0.1286    0.0950], ...
                'String', 'Cancel', 'Callback', @CancelDriftCorrection, ...
                'TooltipString', 'Close without applying correction');
            
            
            function ApplyDriftCorrection(varargin)
                
                w = waitbar(0, 'Apply Drift Correction');
            
                for frm = (handles.DriftFrameRange(1)+1):handles.DriftFrameRange(2)
                    
                        tformMatrix = eye(3,3);
                        tformMatrix(3, 1:2) = handles.DriftMatrixSmoothed(frm, [1 2]);
                        handles.tformObject = affine2d(tformMatrix);
                        
                        for chan = 1:2
                            
                            handles.Img_stack(:, :, frm, chan) = imwarp_same(handles.Img_stack(:,:,frm, chan), ...
                                handles.tformObject);
                        end
                        
                        waitbar((frm - handles.DriftFrameRange(1))/diff(handles.DriftFrameRange), w);

                end
                close(w);
                
                guidata(findobj('tag', 'TIFF viewer'), handles);
                Display_images_in_axes;
                
                close(handles.handles.fig16);
                
            end
            
            function CancelDriftCorrection(varargin)
            
                close(handles.handles.fig16);
                
            end
            
            
            
            
            guidata(findobj('tag', 'TIFF viewer'), handles);
            
        end
        
        function SaveDriftPush(varargin)
            
            handles = guidata(findobj('tag', 'TIFF viewer'));
            
            
            [~, SaveString] = fileparts(handles.Load_file);
            SaveString = strcat(SaveString, '_aligned');
            
            % Save as a 16-bit TIFF file, .mat file, or .dat file
            [filename, pathname, filterspec] = uiputfile(...
                {'*.tif', '8-bit TIFF File (*.tif)'; ...
                '*.tif', '16-bit TIFF File (*.tif)'; ...
                '*.mat', 'MATLAB Data File (*.mat)'; ...
                '*.dat', 'Transform Matrix (*.dat)'}, 'Select Save Filename', SaveString);

            if (filterspec > 0) && (filterspec < 3)

                fnamehere = filename;

                if filterspec == 1;
                    ImgOut = uint8(255*handles.Img_stack/max(handles.Img_stack(:)));
                elseif filterspec == 2;
                    ImgOut = uint16(65535*handles.Img_stack/max(handles.Img_stack(:)));
                end

                %oldObj = Tiff(handles.Load_file);
                tiffObj = Tiff(fullfile(pathname, fnamehere), 'w');

                clear tagstruct;
                tagstruct.ImageLength = size(ImgOut,1);
                tagstruct.ImageWidth = size(ImgOut,2);
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample = 8*filterspec;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.RowsPerStrip = 16;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Software = 'MATLAB';
%                 tagstruct.SubIFD = handles.N_channels*size(handles.Img_stack, 3);
    %             
    %             tagstruct.Orientation =  Tiff.Orientation.TopLeft;
    %             tagstruct.ExtraSamples =  Tiff.ExtraSamples.AssociatedAlpha;
                tagstruct.Compression = Tiff.Compression.None;
                tiffObj.setTag(tagstruct);

%                 tagstruct = rmfield(tagstruct, 'SubIFD');
                w = waitbar(0, 'Saving Image Stack');
                


                for j = 1:size(ImgOut, 3)

                    if j > 1
                        tiffObj.writeDirectory();
                        tiffObj.setTag(tagstruct);
                    end
                    tiffObj.write(ImgOut(:,:,j,1));
                    if handles.N_channels == 2
                        tiffObj.writeDirectory();
                        tiffObj.setTag(tagstruct);
                        tiffObj.write(ImgOut(:,:,j,2));
                    end

                    waitbar(j/size(handles.Img_stack, 3), w);
                end

                tiffObj.close;
                close(w);
                
            elseif filterspec == 3
                %.mat file
                save(fullfile(pathname, filename), '-struct', 'handles', 'Img_stack');
            elseif filterspec == 4
                %.dat transformation matrix
                % Outputs affine matrix used for drift correction as a
                % tab-delimited text file.  Each row corresponds to a new
                % frame between $Firstframe and $Lastframe

                tformFile = fopen(fullfile(pathname, filename), 'w');
                fprintf(tformFile, '# Drift correction saved %s\r\n', datestr(now));
                fprintf(tformFile, '# Reference channel %.f\r\n', handles.DriftRefFrame);
                fprintf(tformFile, '# Frame range %.f\t%.f\r\n', handles.DriftFrameRange(1), handles.DriftFrameRange(2));
                fprintf(tformFile, '# Drift smoothing %.2f\r\n', handles.DriftSmoothing);
                for m = handles.DriftMatrix(handles.DriftFrameRange(1)):handles.DriftMatrix(handles.DriftFrameRange(2))
                    fprintf(tformFile, '%.5f\t%.5f\r\n', handles.DriftMatrix(m,1), handles.DriftMatrix(m,2));
                end
                fclose(tformFile);
                
            end
            
        end
        
        
    end
        
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Flatten_launch(varargin)

    handles = guidata(findobj('Tag', 'TIFF viewer'));
   

    % Launch window to enter type, size of flatten object
    FilterOptions = {'disk', 'line', 'square', 'ball', 'diamond', 'rectangle', 'Gaussian'};
    WhichOption = find(strcmp(FilterOptions, handles.FilterType));
    
    singleRadius = {'diamond', 'disk', 'square', 'octagon', 'Gaussian'};
%     dualRadius = {'line', 'ball', 'rectangle'};


    if ~isempty(findobj('Tag', 'GALAH_FilterBkgd'))

        uistack(findobj('Tag', 'GALAH_FilterBkgd'), 'top');

    else
        
        mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
        fig4_size = [300 100];
        fig4_position = [(mf_post(1) + (mf_post(3) - fig4_size(1))/2) (mf_post(2) + (mf_post(4) - fig4_size(2))/2)];
        handles.handles.fig4 = figure('Name','Flatten Background', 'Tag', 'GALAH_FilterBkgd', 'Units', 'pixels',...
            'Position',[fig4_position fig4_size], 'NumberTitle', 'off', 'Toolbar', 'none', 'Menu', 'none');
        set(handles.handles.fig4, 'Color',[0.9 0.9 0.9]);
        
        handles.handles.Filter_pulldown = uicontrol(handles.handles.fig4, 'Style', 'popupmenu', 'Units', 'normalized', ...
                'Position', [.25 .72 .23 .14], 'String', FilterOptions, 'Value', WhichOption, 'Callback', @Filter_pulldown_call);
            
        handles.handles.Filter_text = uicontrol(handles.handles.fig4, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [.01 .62 .20 .2], 'String', 'Filter Type:', 'BackgroundColor', [0.9 0.9 0.9]);
            
        handles.handles.Radius_text = uicontrol(handles.handles.fig4, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [.03 .19 .13 .2], 'String', 'Radius:', 'BackgroundColor', [0.9 0.9 0.9]);
        
        handles.handles.Radius_box(1) = uicontrol(handles.handles.fig4, 'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.20 .16 .13 .30], 'String', num2str(handles.FilterBkgdRadius(1)), 'Callback', @Filter_edit_call);
        
        handles.handles.Radius_box(2) = uicontrol(handles.handles.fig4, 'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.4 .16 .13 .30], 'String', num2str(handles.FilterBkgdRadius(2)), 'Callback', @Filter_edit_call);
        
        handles.handles.ClearFlatButton = uicontrol(handles.handles.fig4, 'Style', 'pushbutton', 'Units', 'normalized', ...
            'Position', [.60 .68 .38 .28], 'String', 'Reset to Original', 'Callback', @Reset_button_push, ...
            'TooltipString', 'Clear all applied background corrections');
        
        handles.handles.RemoveBkgdButton = uicontrol(handles.handles.fig4, 'Style', 'pushbutton', 'Units', 'normalized', ...
            'Position', [.6 .35 .38 .28], 'String', 'Remove Background', 'Callback', @Remove_button_push, ...
            'TooltipString', 'Remove background via set parameters');
        
        handles.handles.SetandCloseBkgdButton = uicontrol(handles.handles.fig4, 'Style', 'pushbutton', 'Units', 'normalized', ...
            'Position', [.6 .02 .38 .28], 'String', 'Save Flat Stack', 'Callback', @SetBkgdPush, ...
            'TooltipString', 'Save flattened image stack');
        
        
        if ismember(FilterOptions{WhichOption}, singleRadius)
            
            set(handles.handles.Radius_box(2), 'Visible', 'off')
            
        end
       
        
        guidata(findobj('Tag', 'TIFF viewer'), handles);

    end
    
    function Filter_pulldown_call(varargin)

        handles.FilterType = FilterOptions{get(handles.handles.Filter_pulldown, 'Value')};
        
        % Change display of edit boxes, text if needed
        
        if ismember(handles.FilterType, singleRadius)
           set(handles.handles.Radius_box(2), 'Visible', 'off')
        else
           set(handles.handles.Radius_box(2), 'Visible', 'on') 
        end

        
        guidata(findobj('Tag', 'TIFF viewer'), handles);
    end

    function Filter_edit_call(varargin)
        
        if str2double(get(varargin{1}, 'String')) > 0
            handles.FilterBkgdRadius(handles.handles.Radius_box == varargin{1}) = str2double(get(varargin{1}, 'String'));
        end
        
        guidata(findobj('Tag', 'TIFF viewer'), handles);
    end
    
    function Reset_button_push(varargin)
    
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        handles.Img_stack = handles.ImgHold;
        handles.UseFlatKinetics = 0;
        guidata(findobj('Tag', 'TIFF viewer'), handles);
        Display_images_in_axes;

    end
    

    function Remove_button_push(varargin)
        
        set(findobj('Parent', handles.handles.fig4), 'Enable', 'off')
        set(handles.handles.fig4, 'Pointer', 'watch')
        set(handles.handles.fig4, 'UserData', 'Calculating')
        drawnow;
        
        CalcRemoveBkgd
        waitfor(handles.handles.fig4, 'UserData', 'Finished')
        set(findobj('Parent', handles.handles.fig4), 'Enable', 'on')
        set(handles.handles.fig4, 'Pointer', 'arrow')
        set(handles.handles.fig4, 'UserData', [])
        drawnow;
        
        Display_images_in_axes;
        
    end

    function SetBkgdPush(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        
            set(findobj('Parent', handles.handles.fig4), 'Enable', 'off')
            set(handles.handles.fig4, 'Pointer', 'watch')
            set(handles.handles.fig4, 'UserData', 'Calculating')
            drawnow;
        
        if handles.UseFlatKinetics == 0;
        

            CalcRemoveBkgd
            Display_images_in_axes;
            
        end
        
        [~, SaveString] = fileparts(handles.Load_file);
        SaveString = strcat(SaveString, '_flat');
        
            % Save as a 16-bit TIFF file
            [filename, pathname, filterspec] = uiputfile(...
                {'*.tif', '8-bit TIFF File (*.tif)'; ...
                '*.tif', '16-bit TIFF File (*.tif)'; ...
                '*.mat', 'MATLAB Data File (*.mat)'}, ...
                'Select Save Filename', SaveString);

            if (filterspec > 0) && (filterspec < 3)

                fnamehere = filename;

                if filterspec == 1;
                    ImgOut = uint8(255*handles.Img_stack/max(handles.Img_stack(:)));
                elseif filterspec == 2;
                    ImgOut = uint16(65535*handles.Img_stack/max(handles.Img_stack(:)));
                end

                %oldObj = Tiff(handles.Load_file);
                tiffObj = Tiff(fullfile(pathname, fnamehere), 'w');

                clear tagstruct;
                tagstruct.ImageLength = size(ImgOut,1);
                tagstruct.ImageWidth = size(ImgOut,2);
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                tagstruct.BitsPerSample = 8*filterspec;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.RowsPerStrip = 16;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Software = 'MATLAB';
%                 tagstruct.SubIFD = handles.N_channels*size(handles.Img_stack, 3);
    %             
    %             tagstruct.Orientation =  Tiff.Orientation.TopLeft;
    %             tagstruct.ExtraSamples =  Tiff.ExtraSamples.AssociatedAlpha;
                tagstruct.Compression = Tiff.Compression.None;
                tiffObj.setTag(tagstruct);

%                 tagstruct = rmfield(tagstruct, 'SubIFD');
                w = waitbar(0, 'Saving Image Stack');

                for j = 1:size(ImgOut, 3)

                    if j > 1
                        tiffObj.writeDirectory();
                        tiffObj.setTag(tagstruct);
                    end
                    tiffObj.write(ImgOut(:,:,j,1));
                    if handles.N_channels == 2
                        tiffObj.writeDirectory();
                        tiffObj.setTag(tagstruct);
                        tiffObj.write(ImgOut(:,:,j,2));
                    end

                    waitbar(j/size(handles.Img_stack, 3), w);
                end

                tiffObj.close;
                close(w);
                set(handles.handles.fig4, 'UserData', 'Finished');
                
            elseif filterspec == 3
                %.mat file
                save(fullfile(pathname, filename), '-struct', 'handles', 'Img_stack');
                set(handles.handles.fig4, 'UserData', 'Finished');
            else
                % Do nothing
                set(findobj('Parent', handles.handles.fig4), 'Enable', 'on')
                set(handles.handles.fig4, 'Pointer', 'arrow')
                set(handles.handles.fig4, 'UserData', [])
                drawnow;
            end

            waitfor(handles.handles.fig4, 'UserData', 'Finished')
            set(findobj('Parent', handles.handles.fig4), 'Enable', 'on')
            if ishandle(handles.handles.fig4)
                set(handles.handles.fig4, 'Pointer', 'arrow')
                set(handles.handles.fig4, 'UserData', [])
                drawnow;
                close(handles.handles.fig4);
            end
        
    end
    
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function CalcRemoveBkgd(varargin)
        
        handles = guidata(findobj('Tag', 'TIFF viewer'));

        % Correct for background inhomogeneities
        
%         disp(FilterOptions{WhichOption})
        
        if strcmp(handles.FilterType, 'Gaussian')
            
            BkgdFilt = fspecial('Gaussian', handles.FilterBkgdRadius(1)*3, handles.FilterBkgdRadius(1));
            
            w = waitbar(0, 'Flattening Background...');

            for k = 1:handles.N_frames
                for m = 1:handles.N_channels


                    SubtractImg = imfilter(handles.Img_stack(:,:,k, m), BkgdFilt, 'symmetric');
                    handles.Img_stack(:,:,k, m) = (handles.Img_stack(:,:,k,m) - SubtractImg);

                    if m == 1
                        handles.Min_max_left(k, 1) = min(min(handles.Img_stack(:,:,k,m)));
                        handles.Min_max_left(k, 2) = max(max(handles.Img_stack(:,:,k,m)));
                    elseif m == 2
                        handles.Min_max_right(k, 1) = min(min(handles.Img_stack(:,:,k,m)));
                        handles.Min_max_right(k, 2) = max(max(handles.Img_stack(:,:,k,m)));
                    end

                end

                waitbar(k/handles.N_frames, w);

            end
            
        else
        
            if ismember(FilterOptions{WhichOption}, singleRadius)
               BkgdFilt = strel(handles.FilterType, handles.FilterBkgdRadius(1));
            else
               BkgdFilt = strel(handles.FilterType, handles.FilterBkgdRadius(1), handles.FilterBkgdRadius(2));
            end

    %         disp(BkgdFilt)

            w = waitbar(0, 'Flattening Background...');

            for k = 1:handles.N_frames
                for m = 1:handles.N_channels


                    SubtractImg = imopen(handles.Img_stack(:,:,k, m), BkgdFilt);
                    handles.Img_stack(:,:,k, m) = (handles.Img_stack(:,:,k,m) - SubtractImg);

                    if m == 1
                        handles.Min_max_left(k, 1) = min(min(handles.Img_stack(:,:,k,m)));
                        handles.Min_max_left(k, 2) = max(max(handles.Img_stack(:,:,k,m)));
                    elseif m == 2
                        handles.Min_max_right(k, 1) = min(min(handles.Img_stack(:,:,k,m)));
                        handles.Min_max_right(k, 2) = max(max(handles.Img_stack(:,:,k,m)));
                    end

                end

                waitbar(k/handles.N_frames, w);

            end

        end
        
        if min(handles.Min_max_left(:,1)) < 1
            handles.Img_stack(:,:,:,1) = handles.Img_stack(:,:,:,1) - (min(handles.Min_max_left(:,1)) + 1);
            handles.Min_max_left(handles.Min_max_left(:,1) == min(handles.Min_max_left(:,1)),1) = 1;
        end
        
        if (min(handles.Min_max_right(:,1)) < 1) && (handles.N_channels == 2)
            handles.Img_stack(:,:,:,2) = handles.Img_stack(:,:,:,2) - (min(handles.Min_max_right(:,1)) + 1);
            handles.Min_max_right(handles.Min_max_right(:,1) == min(handles.Min_max_right(:,1)),1) = 1;
        end
        
        handles.Display_range_left = [min(handles.Min_max_left(:,1)) max(handles.Min_max_left(:,2))];
        handles.Display_range_right = [min(handles.Min_max_right(:,1)) max(handles.Min_max_right(:,2))];
        
        
        close(w);
        
        
        %handles.UseFlatKinetics = 1;
        
        
        % Display in axes, if open
        
        
        handles.UseFlatKinetics = 1;
        
        % Close it out
        guidata(findobj('Tag', 'TIFF viewer'), handles);
        set(handles.handles.fig4, 'UserData', 'Finished');

    end

 end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unbind_launch functionality deleted 01102015
% See earlier versions for this code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExpFit_launch(varargin)

    if ~isempty(findobj('Tag', 'GALAH_ExpFit'))
        
        uistack(findobj('Tag', 'GALAH_ExpFit'), 'top');
        
    else
        
        opengl('software'); % Fixes issues with transparent overlay
        
        
        fig1 = findobj('Tag', 'TIFF viewer');
        handles = guidata(findobj('Tag', 'TIFF viewer'));

        fig1_post = get(fig1, 'Position');
        fig1_outer = get(fig1, 'OuterPos');
        up_pad = (fig1_outer(4) - fig1_post(4))*handles.scrsz_pixels(4);
        right_pad = (fig1_outer(3) - fig1_post(3))*handles.scrsz_pixels(3);
        
        mf_post = fig1_post.*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
        mf_post(1) = mf_post(1) - fig1_post(3)*handles.scrsz_pixels(3)*.2;
        %fig3_size = [560 560+190];
        
        scrsz = get(0, 'ScreenSize');
        
        Window_size_big = [560 560+190]; % Desired window size in pixels
        Window_size_small = [477 637];
        if scrsz(3) < 1300 || scrsz(4) < 800
            fig5_size = Window_size_small;
        else
            fig5_size = Window_size_big;
        end
        
        
        fig5_position = floor([(mf_post(1) + (mf_post(3) - fig5_size(1))/2) (mf_post(2) + (mf_post(4) - fig5_size(2))/2)]);
        
        fig5_loc = [fig5_size(1)+fig5_position(1)+right_pad, fig5_size(2)+fig5_position(2)+up_pad];
        
        if fig5_loc(1) > handles.scrsz_pixels(3)
            fig5_position(1) = fig5_position(1) - (fig5_loc(1) - handles.scrsz_pixels(3) - 1);
        end
        
        if fig5_loc(2) > handles.scrsz_pixels(4)
            fig5_position(2) = fig5_position(2) - (fig5_loc(2) - handles.scrsz_pixels(4) - 1);
        end
        
        fig6 = figure('Name','Select Region for Exponential Fits', 'Tag', 'GALAH_ExpFit', 'Units', 'pixels',...
            'Position',[fig5_position fig5_size], 'NumberTitle', 'off', 'MenuBar', 'none', 'Toolbar', 'figure');
        
        %%%%%%%%%%%%
        % Set up toolbar
        hToolbar = findall(fig6,'tag','FigureToolBar');
        AllToolHandles = findall(hToolbar);
        ToolBarTags = get(AllToolHandles,'Tag');
        ToolsToKeep = {'FigureToolBar'; 'Exploration.DataCursor'; 'Exploration.Pan'; 'Exploration.ZoomOut'; 'Exploration.ZoomIn'};
        WhichTools = ~ismember(ToolBarTags, ToolsToKeep);
        delete(AllToolHandles(WhichTools));
        
        set(fig6, 'Color',[0.9 0.9 0.9]);
        %axis image;
        
        %%%%%%%%%%%%%%%%%%%
        % Initialize panels
        
        
        fig6_img_panel = uipanel(fig6, 'Units', 'normalized', 'Position', [0, 1-(fig5_size(1)/fig5_size(2)), 1, (fig5_size(1)/fig5_size(2))], ...
            'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'ExpFit_img_panel');
        
        fig6_button_panel = uipanel(fig6, 'Units', 'normalized', 'Position', [0, 0, 1, 1-(fig5_size(1)/fig5_size(2))], ...
            'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'ExpFit_button_panel');
        
        %%%%%%%%%%%%%%%%%%%
        % Initialize axes
        
        ax_ExpFit = axes('Parent', fig6_img_panel, 'Position', [0 0 1 1]);
        set(ax_ExpFit, 'Tag', 'Axis_Kinetics');
        set(ax_ExpFit, 'xtick', [], 'ytick', []);
        axis(ax_ExpFit, 'image');
        
        %opengl software

        %%%%%%%%%%%%%%%%%%%
        % UI Components
        
        if handles.StartFrame == 1
            FrameToStart = calcStartFrame(handles.Primary_channel);
        else
            FrameToStart = handles.StartFrame;
        end

        ExpFit_frame_slider_value = (FrameToStart - 1)/(handles.N_frames-1);
        ExpFit_frame_display = num2str(FrameToStart);
        
        if handles.N_frames == 1;
            ExpFit_frame_step = 1;
        else
            ExpFit_frame_step = 1/(handles.N_frames - 1);
        end
        
        ExpFit_frame_slide_hand = uicontrol(fig6_button_panel, 'Style', 'slider', 'Units', 'normalized',...
            'SliderStep', [ExpFit_frame_step ExpFit_frame_step], 'Min', 0, 'Max', 1, 'Value', ExpFit_frame_slider_value, ...
            'Position', [.42 .82 .57 .07],...
            'Callback', @ExpFit_frame_slide_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'ROI_frame');
        
        addlistener(ExpFit_frame_slide_hand, 'Value', 'PostSet', @ExpFit_frame_listener_call);
        
        ExpFit_frame_slide_box = uicontrol(fig6_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.31 .80 .09 .12], 'BackgroundColor', [1 1 1], ...
            'String', ExpFit_frame_display, 'Callback', @ExpFit_slide_edit_call);
        
        ExpFit_frame_slide_text = uicontrol(fig6_button_panel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [.22 .825 .08 .07], 'BackgroundColor', [.9 .9 .9], ...
            'String', 'Frame:');
        
        %%%%
        
        if handles.Primary_channel == 1;
            min_max_here = handles.Min_max_left;
            display_here = handles.Display_range_left;
        elseif handles.Primary_channel == 2;
            min_max_here = handles.Min_max_right;
            display_here = handles.Display_range_right;
        end
%         
        ExpFit_thresh_step = 1/(max(min_max_here(:,2) - min(min_max_here(:,1))));
        ExpFit_thresh_big_step = min([1e-1 100*ExpFit_thresh_step]);
        
        frame_here = FrameToStart;
        scale_frame = (handles.Img_stack(:, :, frame_here, handles.Primary_channel));
        
        ExpFit_thresh_slider_value = graythresh(scale_frame);
        ExpFit_thresh_display = num2str(round(ExpFit_thresh_slider_value*(max(min_max_here(:,2) - min(min_max_here(:,1)))) + min(min_max_here(:,1))));
        
        ExpFit_thresh_slide_hand = uicontrol(fig6_button_panel, 'Style', 'slider', 'Units', 'normalized',...
            'SliderStep', [ExpFit_thresh_step ExpFit_thresh_big_step], 'Min', 0, 'Max', 1, 'Value', ExpFit_thresh_slider_value, ...
            'Position', [.42 .59 .57 .07],...
            'Callback', @ExpFit_thresh_slide_call, 'BackgroundColor', [.6 .6 .6], 'Tag', 'ExpFit_thresh');
        
        addlistener(ExpFit_thresh_slide_hand, 'Value', 'PostSet', @ExpFit_thresh_listener_call);
        
        ExpFit_thresh_slide_box = uicontrol(fig6_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.31 .57 .09 .12], 'BackgroundColor', [1 1 1], ...
            'String', ExpFit_thresh_display, 'Callback', @ExpFit_thresh_edit_call);
        
        uicontrol(fig6_button_panel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [.195 .59 .11 .08], 'BackgroundColor', [.9 .9 .9], ...
            'String', 'Threshold:');
        
        %%%%%
        
        uicontrol(fig6_button_panel, 'Style', 'text', 'Units', 'normalized', ...
            'Position', [.02 .775 .1 .16], 'BackgroundColor', [.9 .9 .9], ...
            'String', 'Primary Channel:');
        
        Primary_strings = cellstr(num2str((1:size(handles.Img_stack, 4))'));
        
        Channel_pulldown = uicontrol(fig6_button_panel, 'Style', 'popupmenu', 'Units', 'normalized', ...
            'Position', [.13 .78 .06 .14], 'String', Primary_strings, 'Value', handles.Primary_channel, 'Callback', @Active_ExpFit_call);
        
        
        Choose_ROI_hand = uicontrol(fig6_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Set ROI',...
            'Position', [.01 .04 .22 .175],...
            'Callback', @ROI_choose, 'Tag', 'Choose_ROI_botton');
        
        
        ExpFit_Regions_launch_hand = uicontrol(fig6_button_panel, 'Units', 'normalized', ...
            'Style', 'pushbutton', 'String', 'Fit Region Decay',...
            'Position', [.27 .04 .22 .175],...
            'Callback', @ExpRegions_launch, 'Tag', 'ExpFit_regions_launch');

        
        ExpFit_close_hand = uicontrol(fig6_button_panel, 'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Close',...
            'Position', [.77 .04 .22 .175],...
            'Callback', @ExpFit_close, 'Tag', 'ExpFit_close_button');
        
        %%%%%%%%%%%%%%%%%%%
        % Collect handles
        
        handles.handles.fig6 = fig6;
        handles.handles.fig6_img_panel = fig6_img_panel;
        handles.handles.fig6_button_panel = fig6_button_panel;
        handles.handles.ax_ExpFit = ax_ExpFit;
        handles.handles.ExpFit_frame_slide_hand = ExpFit_frame_slide_hand;
        handles.handles.ExpFit_frame_slide_box = ExpFit_frame_slide_box;
        handles.handles.ExpFit_thresh_slide_hand = ExpFit_thresh_slide_hand;
        handles.handles.ExpFit_thresh_slide_box = ExpFit_thresh_slide_box;
        handles.handles.ExpFit_channel_pulldown = Channel_pulldown;
        handles.handles.Choose_ROI_hand = Choose_ROI_hand;
        handles.handles.ExpFit_Regions_launch_hand = ExpFit_Regions_launch_hand;
        handles.handles.ExpFit_close_hand = ExpFit_close_hand;
        
        
        if handles.N_frames == 1;
            
            set(ExpFit_frame_slide_hand, 'Enable', 'off');
            set(ExpFit_frame_slide_box, 'Enable', 'off');
            set(ExpFit_frame_slide_text, 'Enable', 'off');
            
        end
        
        % These start as disabled until an ROI is defined.
        set(handles.handles.ExpFit_thresh_slide_hand, 'Enable', 'off');
        set(handles.handles.ExpFit_thresh_slide_box, 'Enable', 'off');

        % These start as disabled until trajectories are calculated
        set(handles.handles.ExpFit_Regions_launch_hand, 'Enable', 'off');

        guidata(findobj('Tag', 'TIFF viewer'), handles);
        Display_ExpFit_img('initialize');

    end
    %%%%%%%%%%%%%%%%%%%
    % UI Callbacks
    
    function ExpFit_close(varargin)
        close(findobj('Tag', 'GALAH_ExpFit'));
        close(findobj('Tag', 'GALAH_ExpRegions'));
    end
    
    function ExpFit_frame_slide_call(varargin)
    
    % Listener supercedes, replaces all of this

    end
    
    function ExpFit_frame_listener_call(varargin)
    
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        slider_value = get(handles.handles.ExpFit_frame_slide_hand, 'Value');
        frame_num = 1 + slider_value*(handles.N_frames - 1);

        set(handles.handles.ExpFit_frame_slide_box, 'String', num2str(round(frame_num)));

        Display_ExpFit_img('framechange');
        Threshold_ExpFit_display
    
    end
    
    function ExpFit_slide_edit_call(varargin)
    
    handles = guidata(findobj('Tag', 'TIFF viewer'));
    box_value = round(str2double(get(handles.handles.ExpFit_frame_slide_box, 'String')));
    
    %disp(box_value);
    
    if box_value < 1
        
        box_value = 1;
        
    elseif box_value > handles.N_frames
        
        box_value = handles.N_frames;
        
    end
    
    frame_num = (box_value - 1)./(handles.N_frames - 1);
    
    
    set(handles.handles.ExpFit_frame_slide_box, 'String', num2str(box_value));
    set(handles.handles.ExpFit_frame_slide_hand, 'Value', frame_num);
    
    Display_ExpFit_img('framechange');
    Threshold_ExpFit_display
    
    end
    
    function ExpFit_thresh_slide_call(varargin)
    
    % Covered by listener function below
    
    
    end
    
    function ExpFit_thresh_listener_call(varargin)
    
    handles = guidata(findobj('Tag', 'TIFF viewer'));
    ExpFit_thresh_slider_value = get(handles.handles.ExpFit_thresh_slide_hand, 'Value');
    
    if handles.Primary_channel == 1;
        min_max_here = handles.Min_max_left;
    elseif handles.Primary_channel == 2;
        min_max_here = handles.Min_max_right;
    end
    ExpFit_thresh_display = num2str(round(ExpFit_thresh_slider_value*(max(min_max_here(:,2) - min(min_max_here(:,1)))) + min(min_max_here(:,1))));
    
    set(handles.handles.ExpFit_thresh_slide_box, 'String', ExpFit_thresh_display);
    
    
    
    Threshold_ExpFit_display
    
    end
    
    function ExpFit_thresh_edit_call(varargin)
    
    handles = guidata(findobj('Tag', 'TIFF viewer'));
    box_value = round(str2double(get(handles.handles.ExpFit_thresh_slide_box, 'String')));
    
    if handles.Primary_channel == 1;
        min_max_here = handles.Min_max_left;
    elseif handles.Primary_channel == 2;
        min_max_here = handles.Min_max_right;
    end
    
    if box_value < min(min_max_here(:,1))
        
        box_value = min(min_max_here(:,1));
        
    elseif box_value > max(min_max_here(:,2))
        
        box_value = max(min_max_here(:,2));
        
    end
    
    frame_val = (box_value - (min(min_max_here(:,1))))/(max(min_max_here(:,2)) - min(min_max_here(:,1)));
    
    
    set(handles.handles.ExpFit_thresh_slide_box, 'String', num2str(box_value));
    set(handles.handles.ExpFit_thresh_slide_hand, 'Value', frame_val);
    
    
    
                 Display_ExpFit_img('threshold');
                 Threshold_ExpFit_display;
    %
    end
    
    
    function Active_ExpFit_call(varargin)
    
    handles = guidata(findobj('Tag', 'TIFF viewer'));
    active_possible = get(handles.handles.ExpFit_channel_pulldown, 'String');
    active_now = get(handles.handles.ExpFit_channel_pulldown, 'Value');
    
    handles.Primary_channel = str2double(active_possible{active_now});
    guidata(findobj('Tag', 'TIFF viewer'), handles);
    
    %             EnableDisablePlot('off');
    
                  Display_ExpFit_img('framechange');
                  Threshold_ExpFit_display;
    
    end

    function ROI_choose(varargin)
        % Add/update ROI.  Here you're only allowed to have one drawn
        % region in the entire image.
    
        handles = guidata(findobj('Tag', 'TIFF viewer'));

        zoom(handles.handles.fig6, 'off');
        pan(handles.handles.fig6, 'off');

        ExpRectROIHandle = findobj(handles.handles.ax_ExpFit, 'Type', 'hggroup', 'Tag', 'impoly');
        

        if ~isempty(ExpRectROIHandle)
            RectROI = iptgetapi(ExpRectROIHandle);
            RectROI.setVerticesDraggable(1);

        else
            handles.handles.RectROI = impoly(handles.handles.ax_ExpFit);
            handles.handles.RectROI.addNewPositionCallback(@UpdateExpROIPosition);

            % Enable threshold slider, buttons now that there's an ROI
            set(handles.handles.ExpFit_thresh_slide_hand, 'Enable', 'on');
            set(handles.handles.ExpFit_thresh_slide_box, 'Enable', 'on');
            set(handles.handles.ExpFit_Regions_launch_hand, 'Enable', 'on');
    %        set(handles.handles.ExpFit_extract_hand, 'Enable', 'on');

        end

        handles.ExpFitROI = handles.handles.RectROI.getPosition();

        guidata(findobj('Tag', 'TIFF viewer'), handles);
        Threshold_ExpFit_display
    
    end

    function UpdateExpROIPosition(varargin)
        % Function called when ExpFit ROI position is changed
        handles = guidata(findobj('Tag', 'TIFF viewer'));
        
        handles.ExpFitROI = handles.handles.RectROI.getPosition();
        guidata(handles.handles.fig1, handles);
        Threshold_ExpFit_display;
%         get(handles.handles.fig6)
    end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    function ExpRegions_launch(varargin)
     % Launches Exponential Curve Data Fitter window
     
     % Proposed functionality:
     % New window will be launched each time the button in 'Select Region
     % for Exponential Fits' is clicked.  
     % Each time new window is made, the ROI will be drawn in 'SRfEF'
     % window as dashed orange line with number in the middle to show which
     % window that ROI belongs to.  The polygon, frame, and threshold value
     % are held in the guidata for that window only.  It is then saved when
     % the 'Export Results' button is pushed in that window.  
     % This lets you draw one ROI, launch a window, move the ROI, launch
     % another window....
     % 
     

    handles = guidata(findobj('Tag', 'TIFF viewer'));
    
    WindowNumber = handles.ExpFits.WindowCount + 1;
    handles.ExpFits.WindowCount = handles.ExpFits.WindowCount + 1;

    if ~isfield(handles, 'EndFrame')
        handles.EndFrame = handles.N_frames;
    end
    
        ExpFitROIParameters.ExpFits = handles.ExpFits; % Clone default ExpFits data for this window
    
        ExpFitROIParameters.Frame = str2double(get(handles.handles.ExpFit_frame_slide_box, 'String'));
        ExpFitROIParameters.Threshold = str2double(get(handles.handles.ExpFit_thresh_slide_box, 'String'));
        ExpRectROIHandle = findobj(handles.handles.ax_ExpFit, 'Type', 'hggroup', 'Tag', 'impoly');
        ExpFitROIParameters.ROIHandle = iptgetapi(ExpRectROIHandle);
        ExpFitROIParameters.ROIMask = ExpFitROIParameters.ROIHandle.getPosition();
        
        ExpFitROIParameters.StartFrame = handles.StartFrame;
        ExpFitROIParameters.EndFrame = handles.EndFrame;
        
    % Draw ROI in ROI selection window with number it the middle.
    % The number corresponds to the WindowNumber for this panel so window
    % and ROI can be kept straight.  
    
    set(handles.handles.ax_ExpFit, 'NextPlot', 'add');
    
    ExpFitROIParameters.ROIPlotInCallWindow.Line = plot(handles.handles.ax_ExpFit, [ExpFitROIParameters.ROIMask(:,1); ExpFitROIParameters.ROIMask(1,1)], ...
        [ExpFitROIParameters.ROIMask(:,2); ExpFitROIParameters.ROIMask(1,2)], ':', 'LineWidth', 1, 'Color', [1 0.4 .4], 'Hittest', 'off');
    ExpFitROIParameters.ROIPlotInCallWindow.Number = text(mean(ExpFitROIParameters.ROIMask(:,1)), mean(ExpFitROIParameters.ROIMask(:,2)),...
        num2str(WindowNumber), 'Color', [1 .4 .4], 'Parent', handles.handles.ax_ExpFit, 'Hittest', 'off');
    
    set(handles.handles.ax_ExpFit, 'NextPlot', 'replace');

            
    %%%%%%%%%%%%
    % Make new window showing Exponential curve + fit UI
    
    mf_post = get(findobj('Tag', 'TIFF viewer'), 'Position').*([handles.scrsz_pixels(3) handles.scrsz_pixels(4) handles.scrsz_pixels(3) handles.scrsz_pixels(4)]);
    fig7_size = [1050 600];
    fig7_position = [(mf_post(1) + (mf_post(3) - fig7_size(1))/2) (mf_post(2) + (mf_post(4) - fig7_size(2))/2)];
    
    ExpFitROIParameters.handles.ExpRegionsFig = figure('Name',sprintf('Exponential Curve Data Fitter Region #%.f', WindowNumber), ...
        'Tag', 'GALAH_ExpRegions', 'Units', 'pixels',...
        'Position',[fig7_position fig7_size], 'NumberTitle', 'off', 'Toolbar', 'figure', 'Menu', 'none', 'DeleteFcn', @ExpRegion_close);
    set(ExpFitROIParameters.handles.ExpRegionsFig, 'Color',[0.9 0.9 0.9]);
    
    %%%%%%%%%%%%%%%%%%%
    % Initialize panels
    
    
    ExpFitROIParameters.handles.fig7_img_panel = uipanel(ExpFitROIParameters.handles.ExpRegionsFig, 'Units', 'normalized', ...
        'Position', [0.3, 0, 0.7, 1], ...
        'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'ExpRegions_img_panel');
    
    ExpFitROIParameters.handles.fig7_button_panel = uipanel(ExpFitROIParameters.handles.ExpRegionsFig, 'Units', 'normalized', ...
        'Position', [0, 0, 0.3, 1], ...
        'BackgroundColor', [0.9 0.9 0.9], 'BorderType', 'etchedin', 'Tag', 'ExpRegions_button_panel');
    
    %%%%%%%%%%%%%%%%%%%
    % Initialize axes
    
    YRange = [handles.ExpFitROI(1) sum(handles.ExpFitROI([1 3]))];
    XRange = [handles.ExpFitROI(2) sum(handles.ExpFitROI([2 4]))];
    
    ExpFitROIParameters.handles.ax_ExpRegions = axes('Parent', ExpFitROIParameters.handles.fig7_img_panel, 'OuterPosition', [-.05 0.3 1.15 0.74]);
    set(ExpFitROIParameters.handles.ax_ExpRegions, 'Tag', 'Axis_ExpFit');
    set(ExpFitROIParameters.handles.ax_ExpRegions, 'XLim', [0 handles.ExpFits.FrameTime*size(handles.Img_stack, 3)]);
    set(ExpFitROIParameters.handles.ax_ExpRegions, 'YLim', [0 1]);
    xlabel(ExpFitROIParameters.handles.ax_ExpRegions, 'Time (s)'); ylabel(ExpFitROIParameters.handles.ax_ExpRegions, 'Intensity'); 
    
    ExpFitROIParameters.handles.ax_ExpResiduals = axes('Parent', ExpFitROIParameters.handles.fig7_img_panel, 'OuterPosition', [-.05 0.01 1.15 0.31]);
    set(ExpFitROIParameters.handles.ax_ExpResiduals, 'Tag', 'Axis_ExpResiduals');
    set(ExpFitROIParameters.handles.ax_ExpResiduals, 'XLim', [0 handles.ExpFits.FrameTime*size(handles.Img_stack, 3)]);
    set(ExpFitROIParameters.handles.ax_ExpResiduals, 'YLim', [0 1]);
    xlabel(ExpFitROIParameters.handles.ax_ExpResiduals, 'Time (s)'); ylabel(ExpFitROIParameters.handles.ax_ExpResiduals, 'Residuals'); 
    
    ExpFitROIParameters.handles.ax_ExpTextLabels = axes('Parent', ExpFitROIParameters.handles.fig7_button_panel,...
        'Position', [0 0 1 1]);
    set(ExpFitROIParameters.handles.ax_ExpTextLabels, 'Tag', 'Axis_ExpTextLabels');
    
    %%%%%%%%%%%%%%%%%
    % Define areas for two channel parameters as pseudo-panels
    
    ExpFitROIParameters.handles.ChannelBox(1) = plot(ExpFitROIParameters.handles.ax_ExpTextLabels, ...
        [0.01 0.01 0.995 .995 0.01], [0.7 0.995 0.995 0.7 0.7], 'Color', [0.4 0.2 1], 'LineWidth', 2);
    
    
    if handles.N_channels == 2
        set(ExpFitROIParameters.handles.ax_ExpTextLabels, 'NextPlot', 'add');
        ExpFitROIParameters.handles.ChannelBox(2) = plot(ExpFitROIParameters.handles.ax_ExpTextLabels, ...
            [0.01 0.01 0.995 .995 0.01], [0.4 0.69 0.69 0.4 0.4], 'Color', [1 0.4 .4], 'LineWidth', 2);
        
    end
    
    set(ExpFitROIParameters.handles.ax_ExpTextLabels, 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
    
	ExpFitROIParameters.handles.ExpFitChannelText(1) = text(0.06, 0.97, 'Channel 1', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels, 'Color', [0.4 0.2 1]);
    
    if handles.N_channels == 2
	ExpFitROIParameters.handles.ExpFitChannelText(2) = text(0.06, 0.67, 'Channel 2', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels, 'Color', [1 0.4 .4]);
    end
    
    set(ExpFitROIParameters.handles.ax_ExpTextLabels, 'NextPlot', 'replace');
    %%%%%%%%%%%%%
    % Intialize buttons
    
    ExpFitROIParameters.handles.ExpRegionFit = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, 'Units', ...
        'normalized', 'Style', 'pushbutton', 'String', 'Fit',...
        'Position', [.55 .13 .38 .046],...
        'Callback', @ExpRegion_fit, 'Tag', 'ExpRegion_fit_button');
    
    ExpFitROIParameters.handles.ExpRegionExport = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, 'Units', ...
        'normalized', 'Style', 'pushbutton', 'String', 'Export Results',...
        'Position', [.025 .01 .45 .07],...
        'Callback', @ExpRegion_export, 'Tag', 'ExpRegion_close_button', 'Enable', 'off');
    
    
    ExpFitROIParameters.handles.ExpRegionClose = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Units', 'normalized', 'Style', 'pushbutton', 'String', 'Close',...
        'Position', [.525 .01 .45 .07],...
        'Callback', @ExpRegion_close, 'Tag', 'ExpRegion_close_button');
    

    
    
    %%%%%%%%%%%%%
    % Intialize UI
    
    ExpFitROIParameters.handles.TermsText(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'text', 'Units', 'normalized', ...
        'Position', [.39 .95 .29 .03], 'BackgroundColor', [.9 .9 .9], ...
        'String', 'Fit Terms:');
    
    ExpFitROIParameters.handles.TermsPulldown(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'popupmenu', 'Units', 'normalized', ...
        'Position', [.706 .935 .245 .05], 'String', {'1' '2' '3'}, ...
        'Value', ExpFitROIParameters.ExpFits.Nterms(1), 'Callback', @NumFitTermsCallback);
    
	ExpFitROIParameters.handles.ExpFitCheckboxText(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'text', 'Units', 'normalized', ...
        'Position', [.06 .76 .56 .03], 'BackgroundColor', [.9 .9 .9], ...
        'String', 'Float Background:');
  
	ExpFitROIParameters.handles.ExpFitCheckboxBox(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [.75 .7534 .15 .05], 'Value', (ExpFitROIParameters.ExpFits.FloatBkgd{1}), ...
        'Callback', @floatBkgdBoxAction, 'BackgroundColor', [.9 .9 .9]);
    
    %%% Tau and Amp text boxes
    
    ExpFitROIParameters.handles.ExpFitFitText(1,1,1) = text(0.04, 0.92, 'Amp_1', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
    ExpFitROIParameters.handles.ExpFitParamsEdit(1,1,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.165 .90 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.A(1)), ...
        'Callback', @ExpFitParamsEdit);
    
    ExpFitROIParameters.handles.ExpFitFitText(1,2,1) = text(0.48, 0.92, '\tau_1', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
    ExpFitROIParameters.handles.ExpFitParamsEdit(1,2,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.54 .90 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.tau(1)), ...
        'Callback', @ExpFitParamsEdit);
    
	ExpFitROIParameters.handles.ExpFitFitText(1,3,1) = text(0.84, 0.92, 'Fix', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
	ExpFitROIParameters.handles.ExpFixCheckboxBox(1,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'checkbox', 'Units', 'normalized', ...
        'Position', [.913 .8965 .05 .05], 'Value', (ExpFitROIParameters.ExpFits.FixFit(1,1)), ...
        'Callback', @FixValsCheckbox, 'BackgroundColor', [.9 .9 .9]);
    
    %%%
	ExpFitROIParameters.handles.ExpFitFitText(2,1,1) = text(0.04, 0.870, 'Amp_2', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
    ExpFitROIParameters.handles.ExpFitParamsEdit(2,1,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.165 .8525 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.A(2)), ...
        'Callback', @ExpFitParamsEdit);
    
    ExpFitROIParameters.handles.ExpFitFitText(2,2,1) = text(0.48, 0.870, '\tau_2', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
    ExpFitROIParameters.handles.ExpFitParamsEdit(2,2,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.54 .8525 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.tau(2)), ...
        'Callback', @ExpFitParamsEdit);
    
	ExpFitROIParameters.handles.ExpFitFitText(2,3,1) = text(0.84, 0.8734, 'Fix', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
        
	ExpFitROIParameters.handles.ExpFixCheckboxBox(2,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'checkbox', 'Units', 'normalized', ...
            'Position', [.913 .8466 .05 .05], 'Value', (ExpFitROIParameters.ExpFits.FixFit(2,1)), ...
            'Callback', @FixValsCheckbox, 'BackgroundColor', [.9 .9 .9]);
    
    %%%
	ExpFitROIParameters.handles.ExpFitFitText(3,1,1) = text(0.04, 0.825, 'Amp_3', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
    ExpFitROIParameters.handles.ExpFitParamsEdit(3,1,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.165 .805 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.A(3)), ...
        'Callback', @ExpFitParamsEdit);
    
    ExpFitROIParameters.handles.ExpFitFitText(3,2,1) = text(0.48, 0.825, '\tau_3', ...
        'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
    
    ExpFitROIParameters.handles.ExpFitParamsEdit(3,2,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.54 .805 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.tau(3)), ...
        'Callback', @ExpFitParamsEdit);
    
	ExpFitROIParameters.handles.ExpFitFitText(3,3,1) = text(0.84, 0.825, 'Fix', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);
        
	ExpFitROIParameters.handles.ExpFixCheckboxBox(3,1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'checkbox', 'Units', 'normalized', ...
            'Position', [.913 .80 .05 .05], 'Value', (ExpFitROIParameters.ExpFits.FixFit(3,1)), ...
            'Callback', @FixValsCheckbox, 'BackgroundColor', [.9 .9 .9]);
    
    %%%%
    
	ExpFitROIParameters.handles.ExpFitCheckboxText(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'text', 'Units', 'normalized', ...
        'Position', [.07 .715 .49 .03], 'BackgroundColor', [.9 .9 .9], ...
        'String', 'Background Value:');
    
    
    ExpFitROIParameters.handles.ExpFitBackgroundValEdit(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.66 .711 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.ExpFitBkgdValue), ...
        'Callback', @ExpFitBkgdEdit, 'BackgroundColor', [1 1 1]);
    

    %%%%%%%%%%
    % Second channel of UI
    
    
    if handles.N_channels == 2
            ExpFitROIParameters.handles.TermsText(2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'text', 'Units', 'normalized', ...
            'Position', [.39 .65 .29 .03], 'BackgroundColor', [.9 .9 .9], ...
            'String', 'Fit Terms:');

        ExpFitROIParameters.handles.TermsPulldown(2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'popupmenu', 'Units', 'normalized', ...
            'Position', [.706 .635 .245 .05], 'String', {'1' '2' '3'}, ...
            'Value', ExpFitROIParameters.ExpFits.Nterms(2), 'Callback', @NumFitTermsCallback);

        ExpFitROIParameters.handles.ExpFitCheckboxText(2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'text', 'Units', 'normalized', ...
            'Position', [.06 .46 .56 .03], 'BackgroundColor', [.9 .9 .9], ...
            'String', 'Float Background:');

        ExpFitROIParameters.handles.ExpFitCheckboxBox(2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'checkbox', 'Units', 'normalized', ...
            'Position', [.75 .4534 .15 .05], 'Value', (ExpFitROIParameters.ExpFits.FloatBkgd{2}), ...
            'Callback', @floatBkgdBoxAction, 'BackgroundColor', [.9 .9 .9]);

        %%% Tau and Amp text boxes

        ExpFitROIParameters.handles.ExpFitFitText(1,1,2) = text(0.04, 0.62, 'Amp_1', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFitParamsEdit(1,1,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.165 .60 .26 .04], 'String', sprintf('%.2d', ExpFitROIParameters.ExpFits.A(1)), ...
            'Callback', @ExpFitParamsEdit);

        handles.handles.ExpFitFitText(1,2,2) = text(0.48, 0.62, '\tau_1', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFitParamsEdit(1,2,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.54 .60 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.tau(1)), ...
            'Callback', @ExpFitParamsEdit);
            
        ExpFitROIParameters.handles.ExpFitFitText(1,3,3) = text(0.84, 0.619, 'Fix', ...
                'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFixCheckboxBox(1,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
                'Style', 'checkbox', 'Units', 'normalized', ...
                'Position', [.913 .5937 .05 .05], 'Value', (ExpFitROIParameters.ExpFits.FixFit(1,2)), ...
                'Callback', @FixValsCheckbox, 'BackgroundColor', [.9 .9 .9]);

        %%%
        ExpFitROIParameters.handles.ExpFitFitText(2,1,2) = text(0.04, 0.570, 'Amp_2', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFitParamsEdit(2,1,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.165 .5525 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.A(2)), ...
            'Callback', @ExpFitParamsEdit);

        ExpFitROIParameters.handles.ExpFitFitText(2,2,2) = text(0.48, 0.570, '\tau_2', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFitParamsEdit(2,2,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.54 .5525 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.tau(2)), ...
            'Callback', @ExpFitParamsEdit);
                    
        ExpFitROIParameters.handles.ExpFitFitText(2,3,3) = text(0.84, 0.57, 'Fix', ...
                'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFixCheckboxBox(2,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
                'Style', 'checkbox', 'Units', 'normalized', ...
                'Position', [.913 .545 .05 .05], 'Value', (ExpFitROIParameters.ExpFits.FixFit(2,2)), ...
                'Callback', @FixValsCheckbox, 'BackgroundColor', [.9 .9 .9]);

        %%%
        ExpFitROIParameters.handles.ExpFitFitText(3,1,2) = text(0.04, 0.525, 'Amp_3', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFitParamsEdit(3,1,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.165 .505 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.A(3)), ...
            'Callback', @ExpFitParamsEdit);

        ExpFitROIParameters.handles.ExpFitFitText(3,2,2) = text(0.48, 0.525, '\tau_3', ...
            'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFitParamsEdit(3,2,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.54 .505 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.tau(3)), ...
            'Callback', @ExpFitParamsEdit);
                    
        ExpFitROIParameters.handles.ExpFitFitText(3,3,3) = text(0.84, 0.52, 'Fix', ...
                'Parent', ExpFitROIParameters.handles.ax_ExpTextLabels);

        ExpFitROIParameters.handles.ExpFixCheckboxBox(3,2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
                'Style', 'checkbox', 'Units', 'normalized', ...
                'Position', [.913 .495 .05 .05], 'Value', (ExpFitROIParameters.ExpFits.FixFit(3,2)), ...
                'Callback', @FixValsCheckbox, 'BackgroundColor', [.9 .9 .9]);
        %%%%

        ExpFitROIParameters.handles.ExpFitCheckboxText(2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'text', 'Units', 'normalized', ...
            'Position', [.07 .415 .49 .03], 'BackgroundColor', [.9 .9 .9], ...
            'String', 'Background Value:');


        ExpFitROIParameters.handles.ExpFitBackgroundValEdit(2) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
            'Style', 'edit', 'Units', 'normalized', ...
            'Position', [.66 .411 .26 .04], 'String', sprintf('%.2d',ExpFitROIParameters.ExpFits.ExpFitBkgdValue), ...
            'Callback', @ExpFitBkgdEdit, 'BackgroundColor', [1 1 1]);

    end
    
    
            
    ExpFitROIParameters.handles.Start_frame_text = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.2780    0.3134    0.3683    0.0300], 'BackgroundColor', [.9 .9 .9], ...
        'String', 'Analysis Start Frame:');
    
    ExpFitROIParameters.handles.Start_frame_box = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.66 .31 .26 .04], 'BackgroundColor', [1 1 1], ...
        'String', ExpFitROIParameters.StartFrame, 'Callback', @Start_frame_edit_call);
    
    ExpFitROIParameters.handles.End_frame_text = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.3055    0.2634    0.3183    0.0300], 'BackgroundColor', [.9 .9 .9], ...
        'String', 'Analysis End Frame:');
    
    ExpFitROIParameters.handles.End_frame_box = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, 'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.66 .26 .26 .04], 'BackgroundColor', [1 1 1], ...
        'String', num2str(ExpFitROIParameters.EndFrame), 'Callback', @End_frame_edit_call);
    
	ExpFitROIParameters.handles.ExpFitFrameTimeText(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'text', 'Units', 'normalized', ...
        'Position', [0.3666    0.2135    0.2667    0.0300], 'BackgroundColor', [.9 .9 .9], ...
        'String', 'Frame Time (s):');
    
    
    ExpFitROIParameters.handles.ExpFitFrameTimeEdit(1) = uicontrol(ExpFitROIParameters.handles.fig7_button_panel, ...
        'Style', 'edit', 'Units', 'normalized', ...
        'Position', [.66 .21 .26 .04], 'String', num2str(ExpFitROIParameters.ExpFits.FrameTime), ...
        'Callback', @ExpFitFrameEdit, 'BackgroundColor', [1 1 1]);
    
    % Finish setting up Exponential Fitter window
    
    guidata(ExpFitROIParameters.handles.ExpRegionsFig, ExpFitROIParameters)
    guidata(findobj('Tag', 'TIFF viewer'), handles);
    
    ExpParametersCheckEnabled(ExpFitROIParameters);
    FloatBackgroundCallbox(ExpFitROIParameters);
    DisplayExpData(ExpFitROIParameters);
    DisplayStartEndFrame(ExpFitROIParameters);
    
    
    
    
    end
    
    %%%%%%%%%%%%%%
    % Callback functions
    
        function ExpParametersCheckEnabled(ExpFitROIParameters)
            
            %disp('check enabled')
            for k = 1:handles.N_channels
                NTerms = ExpFitROIParameters.ExpFits.Nterms(k);

                set(ExpFitROIParameters.handles.ExpFitParamsEdit(1:3 > NTerms, :, k), ...
                    'Enable', 'off', 'String', '');

                for m = 1:NTerms;

                    set(ExpFitROIParameters.handles.ExpFitParamsEdit(m, 1, k), ...
                        'Enable', 'on', 'String', sprintf('%.2d', (ExpFitROIParameters.ExpFits.A(m, k))),...
                        'BackgroundColor', [1 1 1]);

                    set(ExpFitROIParameters.handles.ExpFitParamsEdit(m, 2, k), ...
                        'Enable', 'on', 'String', sprintf('%.2d', (ExpFitROIParameters.ExpFits.tau(m, k))), ...
                        'BackgroundColor', [1 1 1]);

                end

            end

        end
        
        %%%%%%%%%%%%%%
        
        function floatBkgdBoxAction(varargin)

            FloatBackgroundCallbox(guidata(ancestor(varargin{1}, 'figure', 'toplevel')));

        end

        
        function FloatBackgroundCallbox(ExpFitROIParameters)

            ExpFitROIParameters = guidata(ExpFitROIParameters.handles.ExpRegionsFig);
            CheckboxValues = get(ExpFitROIParameters.handles.ExpFitCheckboxBox, 'Value');
            
            if iscell(CheckboxValues)
                CheckboxValues = cell2mat(CheckboxValues);
            end

            for k = 1:handles.N_channels
                
                ExpFitROIParameters.ExpFits.FloatBkgd{k} = CheckboxValues(k);
                
                if ~get(ExpFitROIParameters.handles.ExpFitCheckboxBox(k), 'Value')
                    set(ExpFitROIParameters.handles.ExpFitBackgroundValEdit(k), 'Enable', 'off', ...
                        'String', '');
                else
                    set(ExpFitROIParameters.handles.ExpFitBackgroundValEdit(k), 'Enable', 'on', ...
                        'String', sprintf('%.2d', (ExpFitROIParameters.ExpFits.ExpFitBkgdValue(k))), ...
                        'BackgroundColor', [1 1 1]);
                end
            end

            guidata(ExpFitROIParameters.handles.ExpRegionsFig, ExpFitROIParameters);
            
        end
        
        function DisplayStartEndFrame(ExpFitROIParameters)
            % Indicate start and end frame for Exponential function fit

            ExpFitROIParameters = guidata(ExpFitROIParameters.handles.ExpRegionsFig);
            
            set(ExpFitROIParameters.handles.ax_ExpRegions, 'NextPlot', 'add');
            
            if isfield(ExpFitROIParameters.handles, 'StartEndLines')
                if ishandle(ExpFitROIParameters.handles.StartEndLines)
                    delete(ExpFitROIParameters.handles.StartEndLines)
                end
            end

            startFrame = str2double(get(ExpFitROIParameters.handles.Start_frame_box, 'String'))-1;
            endFrame = str2double(get(ExpFitROIParameters.handles.End_frame_box, 'String')) - 1;
            
            frameTime = str2double(get(ExpFitROIParameters.handles.ExpFitFrameTimeEdit, 'String'));
            
            yRange = get(ExpFitROIParameters.handles.ax_ExpRegions, 'YLim');
            
            ExpFitROIParameters.handles.StartEndLines(1) = plot(ExpFitROIParameters.handles.ax_ExpRegions, ...
                frameTime*([startFrame startFrame]), yRange, ':', 'Color', [.4 1 .2], 'LineWidth', 2); 
            ExpFitROIParameters.handles.StartEndLines(2) = plot(ExpFitROIParameters.handles.ax_ExpRegions, ...
                frameTime*([endFrame endFrame]), yRange, ':', 'Color', [.4 1 .2], 'LineWidth', 2); 
            
             uistack(ExpFitROIParameters.handles.StartEndLines(1), 'bottom');
             uistack(ExpFitROIParameters.handles.StartEndLines(2), 'bottom');
            
            set(ExpFitROIParameters.handles.ax_ExpRegions, 'NextPlot', 'replace');

             guidata(ExpFitROIParameters.handles.ExpRegionsFig, ExpFitROIParameters);
            
        end
        

    
        function DisplayExpData(ExpFitROIParameters)
            % Display data in Exponential fit function

            handles = guidata(findobj('Tag', 'TIFF viewer'));
            ExpFitROIParameters = guidata(ExpFitROIParameters.handles.ExpRegionsFig);
            
            % Method to extract data values improved ~4x in speed
            
            InMaskValues = reshape(handles.Img_stack((repmat(handles.ExpFitBkgdMask|handles.ExpFitMask,[1 1 size(handles.Img_stack, 3) size(handles.Img_stack, 4)]))), ...
                [sum(handles.ExpFitBkgdMask(:))+sum(handles.ExpFitMask(:)), handles.N_frames, handles.N_channels]);
            
            bkgdInMask = handles.ExpFitBkgdMask(handles.ExpFitBkgdMask|handles.ExpFitMask);
            maskInMask = handles.ExpFitMask(handles.ExpFitBkgdMask|handles.ExpFitMask);
            
            ExpFitROIParameters.ExpBkgdFitData = squeeze(mean(InMaskValues(bkgdInMask, :,:), 1));
            ExpFitROIParameters.ExpFitData = squeeze(mean(InMaskValues(maskInMask, :,:), 1));

            % Older, slow way.  Left to ensure new way matches in resulting
            % data
%             handles.ExpFitData = squeeze(sum(sum(handles.Img_stack.*repmat(handles.ExpFitMask, ...
%                 [1 1 size(handles.Img_stack, 3) size(handles.Img_stack, 4)]), 1), 2))./sum(handles.ExpFitMask(:));
%             
%             
            ExpFitROIParameters.ExpFitData = ExpFitROIParameters.ExpFitData - ExpFitROIParameters.ExpBkgdFitData;
            % Display data in axes
            
            ExpFitROIParameters.timeAxis = 0:ExpFitROIParameters.ExpFits.FrameTime:(size(handles.Img_stack, 3)-1)*ExpFitROIParameters.ExpFits.FrameTime;
                
                ExpFitROIParameters.handles.ExpFitDataPlot(1) = plot(ExpFitROIParameters.handles.ax_ExpRegions,...
                    ExpFitROIParameters.timeAxis, ...
                    ExpFitROIParameters.ExpFitData(:,1), 'o', 'MarkerEdgeColor', [0.6 0.6 1], 'LineWidth', 1, ...
                    'MarkerSize', 4);
                
             if size(handles.Img_stack, 4)  == 2
                 set(ExpFitROIParameters.handles.ax_ExpRegions, 'NextPlot', 'add');
                 ExpFitROIParameters.handles.ExpFitDataPlot(2) = plot(ExpFitROIParameters.handles.ax_ExpRegions,...
                    ExpFitROIParameters.timeAxis, ...
                    ExpFitROIParameters.ExpFitData(:,2), '.', 'MarkerEdgeColor', [1 0.4 0.2], 'LineWidth', 2);
                
                set(ExpFitROIParameters.handles.ax_ExpRegions, 'NextPlot', 'replace');
             end
             
            set(ExpFitROIParameters.handles.ax_ExpRegions, 'XLim', [0 max(ExpFitROIParameters.timeAxis)]);
            xlabel(ExpFitROIParameters.handles.ax_ExpRegions, 'Time (s)'); ylabel(ExpFitROIParameters.handles.ax_ExpRegions, 'Intensity'); 
            
            guidata(ExpFitROIParameters.handles.ExpRegionsFig, ExpFitROIParameters);

        end
        
        %%%%%%%%%%%%%%%%%%
        % Execute fit routine and display results
        function ExpRegion_fit(varargin)
            
            handles = guidata(findobj('Tag', 'TIFF viewer'));
            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            
            %disp('Fit Results')
            set(varargin{1}, 'Enable', 'off')
            set(ExpFitROIParameters.handles.ExpRegionsFig, 'Pointer', 'watch');
            drawnow;
            
            curveFitOptions = optimset('Display', 'off');
            
            % Clear out old curves
            if isfield(ExpFitROIParameters.handles, 'ExpFitPlot')
                if ishandle(ExpFitROIParameters.handles.ExpFitPlot)
                    delete(ExpFitROIParameters.handles.ExpFitPlot);
                end
            end
            
            if isfield(ExpFitROIParameters.handles, 'ExpResidPlot')
                if ishandle(ExpFitROIParameters.handles.ExpResidPlot)
                    delete(ExpFitROIParameters.handles.ExpResidPlot);
                end
            end
            
            if isfield(ExpFitROIParameters.handles, 'eqnText')
                if ishandle(ExpFitROIParameters.handles.eqnText)
                    delete(ExpFitROIParameters.handles.eqnText);
                end
            end
            
            if isfield(ExpFitROIParameters.handles, 'RsquaredText')
                if ishandle(ExpFitROIParameters.handles.RsquaredText)
                    delete(ExpFitROIParameters.handles.RsquaredText);
                end
            end

            ExpFitROIParameters.resids = zeros(ExpFitROIParameters.EndFrame - ExpFitROIParameters.StartFrame + 1, ...
                size(handles.Img_stack, 4));
                
            for fN = 1:handles.N_channels

                % Enter into fit routine
                xData = ExpFitROIParameters.timeAxis(ExpFitROIParameters.StartFrame:ExpFitROIParameters.EndFrame)';
                yData = ExpFitROIParameters.ExpFitData(ExpFitROIParameters.StartFrame:ExpFitROIParameters.EndFrame, fN);


                initial_guesses = [ExpFitROIParameters.ExpFits.A(1:ExpFitROIParameters.ExpFits.Nterms(fN))'; ...
                    ExpFitROIParameters.ExpFits.tau(ExpFitROIParameters.ExpFits.FixFit(:,fN) == 0 & ismember(1:3, 1:ExpFitROIParameters.ExpFits.Nterms(fN))')];
                                
                if (ExpFitROIParameters.ExpFits.FloatBkgd{fN})
                    initial_guesses = [initial_guesses; ExpFitROIParameters.ExpFits.ExpFitBkgdValue(fN)];
                end

                fitIndVars.xdata = xData - ExpFitROIParameters.timeAxis(ExpFitROIParameters.StartFrame);
                fitIndVars.fixedParams = zeros(3,2);
                fitIndVars.fixedParams = [zeros(3,1),...
                    ExpFitROIParameters.ExpFits.tau(:,fN).*ExpFitROIParameters.ExpFits.FixFit(:,fN)];
                
                lowerBound = zeros(numel(initial_guesses), 1);
                
                if numel(initial_guesses) > 0
                
                    fitsOut = lsqcurvefit(@ExpFitFunction, initial_guesses, ...
                        fitIndVars, yData, lowerBound, [], curveFitOptions);
                    
                else
                    fitsOut = [];
                end
                
                ExpFitROIParameters.resids(:,fN) = yData - ExpFitFunction(fitsOut, fitIndVars);
                
            	% Calc and display R^2 value
                SStot = sum((yData - mean(yData)).^2);
                SSres = sum(ExpFitROIParameters.resids(:,fN).^2);
                ExpFitROIParameters.ExpFits.Rsquared(fN) = 1 - SSres/SStot;
                
                if fN == 1
                    ExpFitROIParameters.handles.RsquaredText(fN) = text(0.2*(max(xData(:)) - min(xData(:))), ...
                        0.750*max(yData(:)), ...
                        sprintf('R^2 = %.4f', ExpFitROIParameters.ExpFits.Rsquared(fN)), ...
                        'Color', [0.4 0.2 1], 'Parent', ExpFitROIParameters.handles.ax_ExpRegions);
                    
                elseif fN == 2
                    
                    ExpFitROIParameters.handles.RsquaredText(fN) = text(0.25*(max(xData(:)) - min(xData(:))), ...
                        0.62*max(yData(:)), ...
                        sprintf('R^2 = %.4f', ExpFitROIParameters.ExpFits.Rsquared(fN)), ...
                        'Color', [1 0.4 0.2], 'Parent', ExpFitROIParameters.handles.ax_ExpRegions);
                    
                end
                
                % Collate all parameters that were free in fit and change
                % these in the necessary boxes
                changeParams = zeros(sum(fitIndVars.fixedParams(:,fN) == 1) + numel(fitsOut), 1);
                changeParams(1:floor(numel(changeParams)/2)) = fitsOut(1:floor(numel(changeParams)/2));

                
                numPulled = 1;
                for k = 1:floor(numel(changeParams)/2)
                    if fitIndVars.fixedParams(k, 2) ~= 0
                        changeParams(floor(numel(changeParams)/2) + k) = fitIndVars.fixedParams(k, 2);
                    else
                        changeParams(floor(numel(changeParams)/2) + k) = fitsOut(floor(numel(changeParams)/2) + numPulled);
                        numPulled = numPulled + 1;
                    end
                end

                if mod(numel(changeParams), 2) == 1
                    changeParams(end) = fitsOut(end);
                end
 
                if ExpFitROIParameters.ExpFits.FloatBkgd{fN} == 1
                    ExpFitROIParameters.ExpFits.expFits{fN} = reshape(changeParams(1:end-1), [], 2);
                else
                    ExpFitROIParameters.ExpFits.expFits{fN} = reshape(changeParams, [], 2);
                end
                
                guidata(handles.handles.fig1, handles);
                
                guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
                
                ExpFitROIParameters = ChangeEditBoxesAfterFit(changeParams, fN, ExpFitROIParameters);
                guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
                handles = guidata(handles.handles.fig1);
                
            end
            
            DisplayFitExpCurves(ExpFitROIParameters);
            
            
            set(ExpFitROIParameters.handles.ExpRegionsFig, 'Pointer', 'arrow');
            set(varargin{1}, 'Enable', 'on')
            set(ExpFitROIParameters.handles.ExpRegionExport, 'Enable', 'on');
            drawnow;

        end
        
        function DisplayFitExpCurves(ExpFitROIParameters)
            
            handles = guidata(handles.handles.fig1);
            
            xData = ExpFitROIParameters.timeAxis(ExpFitROIParameters.StartFrame:ExpFitROIParameters.EndFrame)';
            
            % Display fit curves + residuals
            
            fN = 1;
            
            plotData.xdata = xData - ExpFitROIParameters.timeAxis(ExpFitROIParameters.StartFrame);
            plotData.fixedParams = zeros(3,2);
            
            fitIndVars.fixedParams = zeros(3,2);
            fitIndVars.fixedParams = [ExpFitROIParameters.ExpFits.A(:,fN).*ExpFitROIParameters.ExpFits.FixFit(:,fN),...
                ExpFitROIParameters.ExpFits.tau(:,fN).*ExpFitROIParameters.ExpFits.FixFit(:,fN)];

            plotParams = ExpFitROIParameters.ExpFits.expFits{fN}(:);
            
            if ExpFitROIParameters.ExpFits.FloatBkgd{fN} == 1
                plotParams = [plotParams(:); ExpFitROIParameters.ExpFits.ExpFitBkgdValue(fN)];
            end
            
            yData = ExpFitROIParameters.ExpFitData(ExpFitROIParameters.StartFrame:ExpFitROIParameters.EndFrame, fN);

            set(ExpFitROIParameters.handles.ax_ExpRegions, 'NextPlot', 'add')
            ExpFitROIParameters.handles.ExpFitPlot(1) = plot(ExpFitROIParameters.handles.ax_ExpRegions, xData, ...
                ExpFitFunction(plotParams, plotData), ...
                'Color', [0.4 0.2 1], 'LineWidth', 2);

            set(ExpFitROIParameters.handles.ax_ExpResiduals, 'NextPlot', 'add')
            ExpFitROIParameters.handles.ExpResidPlot(1) = plot(ExpFitROIParameters.handles.ax_ExpResiduals, xData, ...
                ExpFitROIParameters.resids(:,1), ...
                'Color', [0.4 0.2 1], 'LineWidth', 1);
            
            if handles.N_channels == 2
                
                fN = 2;
            
                fitIndVars.xdata = xData - ExpFitROIParameters.timeAxis(ExpFitROIParameters.StartFrame);
                fitIndVars.fixedParams = zeros(3,2);
                fitIndVars.fixedParams = [ExpFitROIParameters.ExpFits.A(:,fN).*ExpFitROIParameters.ExpFits.FixFit(:,fN),...
                    ExpFitROIParameters.ExpFits.tau(:,fN).*ExpFitROIParameters.ExpFits.FixFit(:,fN)];
                
                plotParams = ExpFitROIParameters.ExpFits.expFits{fN}((ExpFitROIParameters.ExpFits.FixFit(1:size(ExpFitROIParameters.ExpFits.expFits{fN}, 1), fN) == 0), :);
                
                if ExpFitROIParameters.ExpFits.FloatBkgd{fN} == 1
                    plotParams = [plotParams(:); ExpFitROIParameters.ExpFits.ExpFitBkgdValue(fN)];
                end
                
                yData = ExpFitROIParameters.ExpFitData(ExpFitROIParameters.StartFrame:ExpFitROIParameters.EndFrame, fN);
                
                ExpFitROIParameters.handles.ExpFitPlot(2) = plot(ExpFitROIParameters.handles.ax_ExpRegions, xData, ...
                    ExpFitFunction(plotParams(:), plotData), ...
                    'Color', [1 0.8 0.6], 'LineWidth', 2);
                
                
                ExpFitROIParameters.handles.ExpResidPlot(2) = plot(ExpFitROIParameters.handles.ax_ExpResiduals, xData, ...
                    ExpFitROIParameters.resids(:,2), ...
                    'Color', [1 0.4 0.2], 'LineWidth', 1);
                
            end
            set(ExpFitROIParameters.handles.ax_ExpResiduals, 'XLim', [0 max(ExpFitROIParameters.timeAxis(:))]);
            set(ExpFitROIParameters.handles.ax_ExpResiduals, 'YLim', [1.05*min(ExpFitROIParameters.resids(:)) 1.05*max(ExpFitROIParameters.resids(:))]);
            set(ExpFitROIParameters.handles.ax_ExpRegions, 'NextPlot', 'replace')
            set(ExpFitROIParameters.handles.ax_ExpResiduals, 'NextPlot', 'replace')
            
            
            % Format and display equations            

            for fN = 1:handles.N_channels
                stringOut = ExpFitROIParameters.ExpFits.expFits{fN}(:);
                if ExpFitROIParameters.ExpFits.FloatBkgd{fN} == 1
                    stringOut = [stringOut(:); ExpFitROIParameters.ExpFits.ExpFitBkgdValue(fN)];
                end
            
                if fN == 1
                    colorHere = [.4 .2 1];
                elseif fN == 2
                    colorHere = [1 .4 .2];
                end
                
                stringHere = FormatStringForExpFitLegend(stringOut);
                %textPosition = 0.9*max(yData(:)), 0.1*(max(xData(:)) - min(xData(:)));
                ExpFitROIParameters.handles.eqnText(fN) = text(0.2*(max(xData(:)) - min(xData(:))), max(yData(:))*(0.90-(0.2*(fN-1))), ...
                    stringHere, 'Color', colorHere, 'Parent', ExpFitROIParameters.handles.ax_ExpRegions);

            end

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            guidata(ExpFitROIParameters.handles.ExpRegionsFig, ExpFitROIParameters);
            
        end
        
        
        
        %%%%%%%%%%%%%%%%%%
        % Reset number of fit terms
        function NumFitTermsCallback(varargin)
            
            termsNow = get(varargin{1}, 'Value');
            
            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            
            ExpFitROIParameters.ExpFits.Nterms(ExpFitROIParameters.handles.TermsPulldown == varargin{1}) = termsNow;
            
            guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
            
            ExpParametersCheckEnabled(ExpFitROIParameters)
        end
        
        %%%%%%%%%%%%%%%
        % Edit parameters boxes for Background float value
        
        function ExpFitBkgdEdit(varargin)
            
            
            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            valInput = str2double(get(varargin{1}, 'String'));
            
             valOld = ExpFitROIParameters.ExpFits.ExpFitBkgdValue(ExpFitROIParameters.handles.ExpFitBackgroundValEdit == varargin{1});

            if valInput < 0
                valInput = valOld;
                
            end
            ExpFitROIParameters.ExpFits.ExpFitBkgdValue(ExpFitROIParameters.handles.ExpFitBackgroundValEdit == varargin{1}) = valInput;
            set(varargin{1}, 'String', sprintf('%.2d',valInput));
            guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);

            
        end
        
        %%%%%%%%%%%%%%%
        % Edit parameters boxes for Amp, tau
        
        function ExpFitParamsEdit(varargin)

            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            
            IndexVal = find(ExpFitROIParameters.handles.ExpFitParamsEdit == varargin{1});
            [whichValX, whichValY, whichValZ] = ind2sub(size(ExpFitROIParameters.handles.ExpFitParamsEdit), IndexVal);
            valInput = str2double(get(varargin{1}, 'String'));
            
            if whichValY == 1
                % Amp value
                valOld = ExpFitROIParameters.ExpFits.A(whichValX, whichValZ);

            elseif whichValY == 2
                % tau value
                valOld = ExpFitROIParameters.ExpFits.tau(whichValX, whichValZ);
            end
            
            if valInput < 0
                valInput = valOld;
                
            else
                if whichValY == 1
                    % Amp value
                    ExpFitROIParameters.ExpFits.A(whichValX, whichValZ) = valInput;

                elseif whichValY == 2
                    % tau value 
                    ExpFitROIParameters.ExpFits.tau(whichValX, whichValZ) = valInput;
                end 
            end
            
            set(varargin{1}, 'String', num2str(valInput));
            
            guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
            
            
            
        end
        
            
        
        %%%%%%%%%%%%%%%
        % Edit parameters boxes for Frame time
        
        function ExpFitFrameEdit(varargin)
            
            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            valInput = str2double(get(varargin{1}, 'String'));

             valOld = ExpFitROIParameters.ExpFits.FrameTime;

            if valInput < 0
                valInput = valOld;
                
            end
            
            set(varargin{1}, 'String', num2str(valInput));
            
            % Reset plot axes
            ExpFitROIParameters.timeAxis = 0:valInput:((size(handles.Img_stack, 3)-1)*valInput);
            ExpFitROIParameters.ExpFits.FrameTime = valInput;
            
            % Clear out old curves
            if isfield(ExpFitROIParameters.handles, 'ExpFitPlot')
                if ishandle(ExpFitROIParameters.handles.ExpFitPlot)
                    delete(ExpFitROIParameters.handles.ExpFitPlot);
                end
            end
            
            if isfield(ExpFitROIParameters.handles, 'ExpResidPlot')
                if ishandle(ExpFitROIParameters.handles.ExpResidPlot)
                    delete(ExpFitROIParameters.handles.ExpResidPlot);
                end
            end
            
            if isfield(ExpFitROIParameters.handles, 'eqnText')
                if ishandle(ExpFitROIParameters.handles.eqnText)
                    delete(ExpFitROIParameters.handles.eqnText);
                end
            end
            
            
            guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
            DisplayExpData(ExpFitROIParameters);
            
        end
        
        %%%%%%%%%%%%%%%
        % Edit start frame
            
        function Start_frame_edit_call(varargin)

            handles = guidata(findobj('Tag', 'TIFF viewer'));
            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            ValHere = str2double(get(ExpFitROIParameters.handles.Start_frame_box, 'String'));

            if (ValHere > 0) && (ValHere < handles.N_frames)
                ExpFitROIParameters.StartFrame = ValHere;
            else
                set(ExpFitROIParameters.handles.Start_frame_box, 'String', num2str(ValHere));
            end

            guidata(findobj('Tag', 'TIFF viewer'), handles);
            guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
            DisplayStartEndFrame(ExpFitROIParameters);

        end
        
        %%%%%%%%%%%%%%%
        % Edit end frame
            
        function End_frame_edit_call(varargin)

        handles = guidata(findobj('Tag', 'TIFF viewer'));
        ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
        ValHere = str2double(get(ExpFitROIParameters.handles.End_frame_box, 'String'));

        if (ValHere > ExpFitROIParameters.StartFrame) && (ValHere <= handles.N_frames)
            ExpFitROIParameters.EndFrame = ValHere;
        else
            set(ExpFitROIParameters.handles.End_frame_box, 'String', num2str(ValHere));
        end

        guidata(findobj('Tag', 'TIFF viewer'), handles);
        guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);
        DisplayStartEndFrame(ExpFitROIParameters);

        end
        
        %%%%%%%%%%%%%%%
        % Callback for fix fit values checkboxes
        function FixValsCheckbox(varargin)

            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            
            % Format of ExpFitROIParameters.ExpFits.FixFit matrix is
            % [Param1_channel1_isFit, Param1_channel2_isFit;...
            ExpFitROIParameters.ExpFits.FixFit(ExpFitROIParameters.handles.ExpFixCheckboxBox == varargin{1}) = get(varargin{1}, 'Value');

            guidata(ancestor(varargin{1}, 'figure', 'toplevel'), ExpFitROIParameters);

        end
        
        
        %%%%%%%%%%%%%%%
        % Export Data
        
        function ExpRegion_export(varargin)

            handles = guidata(findobj('Tag', 'TIFF viewer'));
            expData = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            
            StartFrame = str2double(get(expData.handles.Start_frame_box, 'String'));
            EndFrame = str2double(get(expData.handles.End_frame_box, 'String'));
            
            [putfilename, pathname, ~] = uiputfile([fileparts(handles.Load_file), '\*.txt'], ...
                'TXT file');
            
            if isequal(putfilename,0) || isequal(pathname,0)
                % Cancel button was pressed, do nothing.
            else
                expObj = fopen(fullfile(pathname, putfilename), 'wt');
                fwrite(expObj, sprintf('Data File : %s\n', handles.Load_file));
                fwrite(expObj, sprintf('Frame Time (s): %.4f\n', expData.ExpFits.FrameTime));
                fwrite(expObj, sprintf('Primary Channel: %.4f\n', handles.Primary_channel));
                fwrite(expObj, sprintf('Start Frame : %.f\n', str2double(get(expData.handles.Start_frame_box, 'String'))));
                fwrite(expObj, sprintf('End Frame : %.f\n', str2double(get(expData.handles.End_frame_box, 'String'))));
                fwrite(expObj, sprintf('Threshold Level : %.f\n', expData.Threshold));
                fwrite(expObj, sprintf('ROI Boundaries : [%.2f %.2f;\n', expData.ROIMask(1, 1:2)));
                for k = 2:(size(expData.ROIMask, 1)-1)
                    fwrite(expObj, sprintf('\t\t%.2f %.2f;\n', expData.ROIMask(k, 1:2)));
                end
                fwrite(expObj, sprintf('\t\t%.2f %.2f];\n', expData.ROIMask(end, 1:2)));
                
                for k = 1:handles.N_channels
                    fwrite(expObj, sprintf('\n'));
                    fwrite(expObj, sprintf('--------------------------------------------------------\n'));
                    fwrite(expObj, sprintf('\n'));
                    
                    fwrite(expObj, sprintf('Channel %.f\n', k));
                    fwrite(expObj, sprintf('N Fit Terms : %.f\n', expData.ExpFits.Nterms(k)));
                    
                    fixedString = '';
                    for m = 1:expData.ExpFits.Nterms(k)
                        if expData.ExpFits.FixFit(m,k) == 1
                            fixedString = sprintf('%s\t%s', fixedString, 'True');
                        else
                            fixedString = sprintf('%s\t%s', fixedString, 'False');
                        end
                    end
                    
                    fwrite(expObj, sprintf('Tau Fixed : %s\n', fixedString));

                    fwrite(expObj, sprintf('Amp1 : %.6f\ttau1 : %.6f\n', expData.ExpFits.A(1, k), expData.ExpFits.tau(1, k)));
                    if expData.ExpFits.Nterms(k) > 1
                        fwrite(expObj, sprintf('Amp2 : %.6f\ttau2 : %.6f\n', expData.ExpFits.A(2, k), expData.ExpFits.tau(2, k)));
                    end
                    if expData.ExpFits.Nterms(k) > 2
                        fwrite(expObj, sprintf('Amp3 : %.6f\ttau3 : %.6f\n', expData.ExpFits.A(3, k), expData.ExpFits.tau(3, k)));
                    end

                    if expData.ExpFits.FloatBkgd{k} 
                        bkgdString = 'True';
                    else
                        bkgdString = 'False';
                    end

                    fwrite(expObj, sprintf('Float Background : %s\n', bkgdString));
                    if expData.ExpFits.FloatBkgd{k}
                        fwrite(expObj, sprintf('Background Fit : %.4f\n', expData.ExpFits.ExpFitBkgdValue(k)));
                    end

                    fwrite(expObj, sprintf('R squared : %.4f\n', expData.ExpFits.Rsquared(k)));

                    fwrite(expObj, sprintf('\n'));
                    fwrite(expObj, sprintf('Time(s)\tData\tFit\tResiduals\n'));

                    % Block write data in format [time(s) Data Fit resids]

                    inData.xdata = expData.timeAxis(StartFrame:EndFrame) - expData.timeAxis(StartFrame);
                    inData.fixedParams = zeros(3,2);
                    
                    plotData = expData.ExpFits.expFits{k}(:);
                    if expData.ExpFits.FloatBkgd{k}
                        plotData = [plotData; expData.ExpFits.ExpFitBkgdValue(k)];
                    end
                    
                    PadFit = [nan(1, StartFrame-1), ...
                        ExpFitFunction(plotData, inData), ...
                        nan(1, (handles.N_frames - EndFrame))];
                    PadResids = [nan(StartFrame-1, 1); expData.resids(:,k); nan((handles.N_frames - EndFrame), 1)];

                    BlockToWrite = [expData.timeAxis' expData.ExpFitData(:,k) PadFit' PadResids];

                    for m = 1:size(BlockToWrite, 1)
                        fwrite(expObj, sprintf('%.4f\t%.4f\t%.4f\t%.4f\n', BlockToWrite(m, :)));
                    end
                    fwrite(expObj, sprintf('\n'));
                    
                end

                fclose(expObj);
                
                guidata(findobj('Tag', 'TIFF viewer'), handles);
                guidata(ancestor(varargin{1}, 'figure', 'toplevel'), expData);
            end
 
            
        end

        
        %%%%%%%%%%%%%%%
        % Close Fitter 
        
        function ExpRegion_close(varargin)
            
            ExpFitROIParameters = guidata(ancestor(varargin{1}, 'figure', 'toplevel'));
            
            % Delete ROI drawn in the ROI selection window
            if ishandle(ExpFitROIParameters.ROIPlotInCallWindow.Number)
                delete(ExpFitROIParameters.ROIPlotInCallWindow.Number);
            end
            
            if ishandle(ExpFitROIParameters.ROIPlotInCallWindow.Line)
                delete(ExpFitROIParameters.ROIPlotInCallWindow.Line);
            end
            
            % Close the figure window itself
            close(ancestor(varargin{1}, 'figure', 'toplevel'));

        end
        
    end
    
    %%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function Display_ExpFit_img(Reason) %
    
        handles = guidata(findobj('Tag', 'TIFF viewer'));

        if ~isempty(findobj('Tag', 'GALAH_ExpFit'));

            set(handles.handles.ax_ExpFit, 'NextPlot', 'add');

            if ~strcmp(Reason, 'adjust')
                % Clear out old lines, patches if this was a change in frame
                % number to call Display_ExpFit_img

                delete(findobj('Parent', handles.handles.ax_ExpFit, 'Type', 'line'));
                delete(findobj('Parent', handles.handles.ax_ExpFit, 'Type', 'patch'));
            end

            % Clear out old image
            delete(findobj('Tag', 'ExpFit_image'));

            ax_ExpFit = handles.handles.ax_ExpFit;
            frame_slider = str2double(get(handles.handles.ExpFit_frame_slide_box, 'String')); % Set this to ExpFit_slider

            if handles.Primary_channel == 1;
                min_max_here = handles.Min_max_left;
                display_here = handles.Display_range_left;
                invert_here = handles.Left_invert;
                Autoscale_here = handles.Autoscale_left;
                yaxis_direction = get(handles.handles.ax1, 'YDir');
            elseif handles.Primary_channel == 2;
                min_max_here = handles.Min_max_right;
                display_here = handles.Display_range_right;
                invert_here = handles.Right_invert;
                Autoscale_here = handles.Autoscale_right;
                yaxis_direction = get(handles.handles.ax2, 'YDir');
            end

            % Check if Autoscale

            if Autoscale_here == 1;
                display_here = min_max_here(frame_slider, :);
            end

            ExpFit_image = image(Vector2Colormap_setscale(handles.Img_stack(:,:,frame_slider,handles.Primary_channel), 'gray', display_here), ...
                'Parent', ax_ExpFit, 'Tag', 'ExpFit_image');
            set(ax_ExpFit, 'xtick', [], 'ytick', []);

            if strcmp(Reason, 'initialize')
                set(handles.handles.ax_ExpFit, 'XLim', [0.5 size(handles.Img_stack, 1)+0.5]);
                set(handles.handles.ax_ExpFit, 'YLim', [0.5 size(handles.Img_stack, 1)+0.5]);
            end

            % Invert if needed
            if invert_here == 1;
                set(ax_ExpFit, 'XDir', 'reverse');
            else
                set(ax_ExpFit, 'XDir', 'normal');
            end

            set(ax_ExpFit, 'YDir', yaxis_direction);

            % Move image to back
            uistack(ExpFit_image, 'bottom');

            handles.handles.ExpFit_image = ExpFit_image;
            
            %%%%%%%%%%%%
            % Check for already-open Exponential Fit windows
            alreadyOpenExpRegions = findobj('Tag', 'GALAH_ExpRegions');
            if ~isempty(alreadyOpenExpRegions)
                
                for m = 1:numel(alreadyOpenExpRegions)
                    
                    ExpFitOpen = guidata(alreadyOpenExpRegions(m));
                    plot(handles.handles.ax_ExpFit, [ExpFitOpen.ROIMask(:,1); ExpFitOpen.ROIMask(1,1)], ...
                        [ExpFitOpen.ROIMask(:,2); ExpFitOpen.ROIMask(1,2)], ':', 'LineWidth', 1, 'Color', [1 0.4 .4], 'Hittest', 'off');
                    text(mean(ExpFitOpen.ROIMask(:,1)), mean(ExpFitOpen.ROIMask(:,2)),...
                        num2str(ExpFitOpen.ExpFits.WindowCount), 'Color', [1 .4 .4], 'Parent', handles.handles.ax_ExpFit, 'Hittest', 'off');
                    
                end
            end

            set(handles.handles.ax_ExpFit, 'NextPlot', 'replace');

            guidata(findobj('Tag', 'TIFF viewer'), handles);

        end
    
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function Threshold_ExpFit_display(varargin)
    
    
    handles = guidata(findobj('Tag', 'TIFF viewer'));
    
    if strcmp(get(handles.handles.ExpFit_thresh_slide_hand, 'Enable'), 'on')

        XaxisLim = get(handles.handles.ax_ExpFit, 'XLim');
        YaxisLim = get(handles.handles.ax_ExpFit, 'YLim');
        
        % Clear out old lines, ROI values
        delete(findobj('Parent', handles.handles.ax_ExpFit, 'Type', 'line', 'LineStyle', '-'));
        delete(findobj('Parent', handles.handles.ax_ExpFit, 'Type', 'patch'));
        delete(findobj('Parent', handles.handles.ax_ExpFit, 'Tag', 'OverlayImg'));
        handles.handles.ExpFit_line = [];
        handles.handles.ExpFit_box ={};
        handles.handles.ExpFit_color = [];
        handles.ExpFit_output = {};
        handles.ExpFit_parameters = [];

        frame_here = str2double(get(handles.handles.ExpFit_frame_slide_box, 'String'));
        
        if handles.Primary_channel == 1;
            min_max_here = handles.Min_max_left;
            display_here = handles.Display_range_left;
        elseif handles.Primary_channel == 2;
            min_max_here = handles.Min_max_right;
            display_here = handles.Display_range_right;
        end
        
        min_max_here = min_max_here(handles.StartFrame, :);
        
        scale_frame = (handles.Img_stack(:, :, frame_here, handles.Primary_channel));
        
        
        thresh_here = get(handles.handles.ExpFit_thresh_slide_hand, 'Value');

        NowBW = scale_frame;
        
        NowBW(NowBW < thresh_here*max(min_max_here(:))) = 0;
        NowBW(NowBW > 0) = 1;

        % impoly to select ROI
        BoxMask = poly2mask(round(handles.ExpFitROI(:,1)), round(handles.ExpFitROI(:,2)), size(NowBW, 1), size(NowBW, 2));
        NowBW = NowBW.*BoxMask;
        
        NowBW = bwmorph(NowBW, 'clean');
        NowBW = bwmorph(NowBW, 'thicken');

        [B,L] = bwboundaries(NowBW,'noholes');
        
        
        handles.ExpFitMask = NowBW;
        handles.ExpFitBkgdMask = handles.handles.RectROI.createMask() & ~ NowBW;
        set(handles.handles.ax_ExpFit, 'NextPlot', 'add');
        
        Border_plot = zeros(length(B), 1);
        if ~isempty(Border_plot)
            
            for k = 1:length(B);

                boundary = B{(k)};
                Border_plot(k) = plot(handles.handles.ax_ExpFit, boundary(:,2), boundary(:,1), ...
                    'Color', handles.handles.Unselected_ROI_color, 'LineWidth', 1);
                
            end
            
            handles.handles.ExpFit_Border_plot = Border_plot;
            handles.ExpFit_boundaries.B = B;
            handles.ExpFit_boundaries.L = L;
        else
            handles.handles.ExpFit_Border_plot = [];
            handles.ExpFit_boundaries.B = [];
            handles.ExpFit_boundaries.L = [];
        end
        
        
        
        set(handles.handles.ax_ExpFit, 'NextPlot', 'replace');
        set(handles.handles.ExpFit_Border_plot, 'hittest', 'off');
        guidata(findobj('Tag', 'TIFF viewer'), handles);

    else
        
        % Threshold is disabled - do nothing
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frameStart = calcStartFrame(whichChannel)

    handles = guidata(findobj('Tag', 'TIFF viewer'));
    % Figure out the frame to start watching for

    VarInTime = var(reshape(handles.Img_stack(:,:,:,whichChannel), ...
        [size(handles.Img_stack, 1)*size(handles.Img_stack, 2), size(handles.Img_stack, 3)]), [], 1);

    % Find first frame that is at least 50% of the peak.
    frameStart = find(VarInTime > 0.5*max(VarInTime(:)), 1, 'first');
    handles.StartFrame = frameStart;
    guidata(findobj('Tag', 'TIFF viewer'), handles);

end

function ExpFcnOutFinal = ExpFitFunction(inParams, indIn)

% Based on number of parameters specified, determine functional form of the
% Exponential fit function

xdata = indIn.xdata;
fixedParams = indIn.fixedParams;

parameters = zeros(sum(indIn.fixedParams(:,2) > 1) + numel(inParams), 1);
parameters(1:floor(numel(parameters)/2)) = inParams(1:floor(numel(parameters)/2));

numPulled = 1;
for k = 1:floor(numel(parameters)/2)
    if fixedParams(k, 2) ~= 0
        parameters(floor(numel(parameters)/2) + k) = fixedParams(k, 2);
    else
        parameters(floor(numel(parameters)/2) + k) = inParams(floor(numel(parameters)/2) + numPulled);
        numPulled = numPulled + 1;
    end
end

if mod(numel(parameters), 2) == 1
    parameters(end) = inParams(end);
end

switch numel(parameters)
    
    case 0
        
        ExpFcnOut = 0;
        
    case 1
        
        bkgd = parameters(1);
        
        ExpFcnOut = bkgd;
           
    
    case 2
        
        A1 = parameters(1);
        tau1 = parameters(2);
        
        ExpFcnOut = A1.*exp(-tau1*xdata);
        
    case 3
        
        A1 = parameters(1);
        tau1 = parameters(2);
        bkgd = parameters(3);
        
        ExpFcnOut = A1.*exp(-tau1*xdata) + bkgd;
        
    case 4
        
        A1 = parameters(1);
        tau1 = parameters(3);
        A2 = parameters(2);
        tau2 = parameters(4);
        
        ExpFcnOut = A1.*exp(-tau1*xdata) + A2.*exp(-tau2*xdata);
        
    case 5
        
        A1 = parameters(1);
        tau1 = parameters(3);
        A2 = parameters(2);
        tau2 = parameters(4);
        bkgd = parameters(5);
        
        ExpFcnOut = A1.*exp(-tau1*xdata) + A2.*exp(-tau2*xdata) + bkgd;
        
    case 6
        
        A1 = parameters(1);
        tau1 = parameters(4);
        A2 = parameters(2);
        tau2 = parameters(5);
        A3 = parameters(3);
        tau3 = parameters(6);
        
        ExpFcnOut = A1.*exp(-tau1*xdata) + A2.*exp(-tau2*xdata) + A3.*exp(-tau3*xdata);
        
    case 7
        A1 = parameters(1);
        tau1 = parameters(4);
        A2 = parameters(2);
        tau2 = parameters(5);
        A3 = parameters(3);
        tau3 = parameters(6);
        bkgd = parameters(7);
        
        ExpFcnOut = A1.*exp(-tau1*xdata) + A2.*exp(-tau2*xdata) + A3.*exp(-tau3*xdata) + bkgd;
end

for k = 1:size(fixedParams, 1)
    
    ExpFcnOut = ExpFcnOut + fixedParams(k,1).*exp(-fixedParams(k,2)*xdata);
end

ExpFcnOutFinal = ExpFcnOut;

end

function stringOut = FormatStringForExpFitLegend(parameters)
% Format output exponential fit function based on specified number of
% parameters
% Number of parameters parsed into number of components included in fit

    switch numel(parameters)
        
        case 1 
            
            bkgd = parameters(1);

            stringOut = sprintf('Int = %.3f', bkgd);
        
        case 2

            A1 = parameters(1);
            tau1 = parameters(2);

            stringOut = sprintf('Int = %.3f * e^{-%.3ft}', A1, tau1);

        case 3

            A1 = parameters(1);
            tau1 = parameters(2);
            bkgd = parameters(3);

            stringOut = sprintf('Int = %.3f * e^{-%.3ft} + %.3f', A1, tau1, bkgd);

        case 4

            A1 = parameters(1);
            tau1 = parameters(3);
            A2 = parameters(2);
            tau2 = parameters(4);

            stringOut = sprintf('Int = %.3f * e^{-%.3ft} + %.3f * e^{-%.3ft}', ...
                A1, tau1, A2, tau2);

        case 5

            A1 = parameters(1);
            tau1 = parameters(3);
            A2 = parameters(2);
            tau2 = parameters(4);
            bkgd = parameters(5);
            
            stringOut = sprintf('Int = %.3f * e^{-%.3ft} + %.3f * e^{-%.3ft} + %.3f', ...
                A1, tau1, A2, tau2, bkgd);

        case 6

            A1 = parameters(1);
            tau1 = parameters(4);
            A2 = parameters(2);
            tau2 = parameters(5);
            A3 = parameters(3);
            tau3 = parameters(6);

            stringOut = sprintf('Int = %.3f * e^{-%.3ft} + %.3f * e^{-%.3ft} + %.3f * e^{-%.3ft}', ...
                A1, tau1, A2, tau2, A3, tau3);

        case 7
            A1 = parameters(1);
            tau1 = parameters(4);
            A2 = parameters(2);
            tau2 = parameters(5);
            A3 = parameters(3);
            tau3 = parameters(6);
            bkgd = parameters(7);

            stringOut = sprintf('Int = %.3f * e^{-%.3ft} + %.3f * e^{-%.3ft} + %.3f * e^{-%.3ft} + %.3f', ...
                A1, tau1, A2, tau2, A3, tau3, bkgd);
    end

end

function NewHandles = ChangeEditBoxesAfterFit(parameters, k, handles)
% Given number of parameters for exponential fit, change the number of
% visible edit boxes in the GUI window
% handles = guidata(findobj('Tag', 'TIFF viewer'));

    switch numel(parameters)
        case 2
            
            handles.ExpFits.A(1,k) = parameters(1);
            set(handles.handles.ExpFitParamsEdit(1,1,k), 'String', sprintf('%.2d', parameters(1)));
            handles.ExpFits.tau(1,k) = parameters(2);
            set(handles.handles.ExpFitParamsEdit(1,2,k), 'String', sprintf('%.2d', parameters(2)));

        case 3
            
            handles.ExpFits.A(1,k) = parameters(1);
            set(handles.handles.ExpFitParamsEdit(1,1,k), 'String', sprintf('%.2d', parameters(1)));
            handles.ExpFits.tau(1,k) = parameters(2);
            set(handles.handles.ExpFitParamsEdit(1,2,k), 'String', sprintf('%.2d', parameters(2)));
            handles.ExpFits.ExpFitBkgdValue(k) = parameters(3);
            set(handles.handles.ExpFitBackgroundValEdit(k), 'String', sprintf('%.2d', parameters(3)));

        case 4

            handles.ExpFits.A(1,k) = parameters(1);
            set(handles.handles.ExpFitParamsEdit(1,1,k), 'String', sprintf('%.2d', parameters(1)));
            handles.ExpFits.tau(1,k) = parameters(3);
            set(handles.handles.ExpFitParamsEdit(1,2,k), 'String', sprintf('%.2d', parameters(3)));
            handles.ExpFits.A(2,k) = parameters(2);
            set(handles.handles.ExpFitParamsEdit(2,1,k), 'String', sprintf('%.2d', parameters(2)));
            handles.ExpFits.tau(2,k) = parameters(4);
            set(handles.handles.ExpFitParamsEdit(2,2,k), 'String', sprintf('%.2d', parameters(4)));            

        case 5

            handles.ExpFits.A(1,k) = parameters(1);
            set(handles.handles.ExpFitParamsEdit(1,1,k), 'String', sprintf('%.2d', parameters(1)));
            handles.ExpFits.tau(1,k) = parameters(3);
            set(handles.handles.ExpFitParamsEdit(1,2,k), 'String', sprintf('%.2d', parameters(3)));
            handles.ExpFits.A(2,k) = parameters(2);
            set(handles.handles.ExpFitParamsEdit(2,1,k), 'String', sprintf('%.2d', parameters(2)));
            handles.ExpFits.tau(2,k) = parameters(4);
            set(handles.handles.ExpFitParamsEdit(2,2,k), 'String', sprintf('%.2d', parameters(4))); 
            handles.ExpFits.ExpFitBkgdValue(k) = parameters(5);
            set(handles.handles.ExpFitBackgroundValEdit(k), 'String', sprintf('%.2d', parameters(5)));

        case 6

            handles.ExpFits.A(1,k) = parameters(1);
            set(handles.handles.ExpFitParamsEdit(1,1,k), 'String', sprintf('%.2d', parameters(1)));
            handles.ExpFits.tau(1,k) = parameters(4);
            set(handles.handles.ExpFitParamsEdit(1,2,k), 'String', sprintf('%.2d', parameters(4)));
            handles.ExpFits.A(2,k) = parameters(2);
            set(handles.handles.ExpFitParamsEdit(2,1,k), 'String', sprintf('%.2d', parameters(2)));
            handles.ExpFits.tau(2,k) = parameters(5);
            set(handles.handles.ExpFitParamsEdit(2,2,k), 'String', sprintf('%.2d', parameters(5)));  
            handles.ExpFits.A(3,k) = parameters(3);
            set(handles.handles.ExpFitParamsEdit(3,1,k), 'String', sprintf('%.2d', parameters(3)));
            handles.ExpFits.tau(3,k) = parameters(6);
            set(handles.handles.ExpFitParamsEdit(3,2,k), 'String', sprintf('%.2d', parameters(6)));  

        case 7
            
            handles.ExpFits.A(1,k) = parameters(1);
            set(handles.handles.ExpFitParamsEdit(1,1,k), 'String', sprintf('%.2d', parameters(1)));
            handles.ExpFits.tau(1,k) = parameters(4);
            set(handles.handles.ExpFitParamsEdit(1,2,k), 'String', sprintf('%.2d', parameters(4)));
            handles.ExpFits.A(2,k) = parameters(2);
            set(handles.handles.ExpFitParamsEdit(2,1,k), 'String', sprintf('%.2d', parameters(2)));
            handles.ExpFits.tau(2,k) = parameters(5);
            set(handles.handles.ExpFitParamsEdit(2,2,k), 'String', sprintf('%.2d', parameters(5)));  
            handles.ExpFits.A(3,k) = parameters(3);
            set(handles.handles.ExpFitParamsEdit(3,1,k), 'String', sprintf('%.2d', parameters(3)));
            handles.ExpFits.tau(3,k) = parameters(6);
            set(handles.handles.ExpFitParamsEdit(3,2,k), 'String', sprintf('%.2d', parameters(6))); 
            handles.ExpFits.ExpFitBkgdValue(k) = parameters(7);
            set(handles.handles.ExpFitBackgroundValEdit(k), 'String', sprintf('%.2d', parameters(7)));
    end

    NewHandles = handles;
end


