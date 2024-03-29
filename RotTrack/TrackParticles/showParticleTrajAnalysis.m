function processedTraj = showParticleTrajAnalysis(trajXlsPath,image_data,analysedAllTraj,n_traj,start_frame,tsamp,pixelsize_nm) 
%
% ========================================
% RotTrack.
% Copyright (c) 2017. Isabel Llorente-Garcia, Dept. of Physics and Astronomy, University College London, United Kingdom.
% Released and licensed under a BSD 2-Clause License:
% https://github.com/illg-ucl/rotTrack/blob/master/LICENSE
% This program is free software: you can redistribute it and/or modify it under the terms of the BSD 2-Clause License.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the BSD 2-Clause License for more details. You should have received a copy of the BSD 2-Clause License along with this program.
% Citation: If you use this software for your data analysis please acknowledge it in your publications and cite as follows.
% -	Citation example 1: 
% RotTrack software. (Version). 2017. Isabel Llorente-Garcia, 
% Dept. of Physics and Astronomy, University College London, United Kingdom.
% https://github.com/illg-ucl/rotTrack. (Download date).
% 
% -	Citation example 2:
% @Manual{... ,
% title  = {RotTrack software. (Version).},
% author       = {{Isabel Llorente-Garcia}},
% organization = { Dept. of Physics and Astronomy, University College London, United Kingdom.},
% address      = {Gower Place, London, UK.},
% year         = 2017,
% url          = {https://github.com/illg-ucl/rotTrack}}
% ========================================
%
%
% Analyse and process data for a single trajectory/track (for trajectory number n_traj), for a given image sequence:
% - Plot orientation angle (degrees) versus time; fit all points to a line to get
% angular velocity.
% - plot x and y (referenced to first point in track) versus time, 
% - plot trajectory on the x-y plane.
% - plot mean square displacement (msd) vs Delta t with error bars first
% all points, then only points with relative error < 150%. Fit the latter
% (tries fitting to a line and to a confined curve and saves results).
% - show the first frame with the trajectory overlaid.
% - print the trajectory number and other info, including angular velocity from fit.
% Saves all plots in a .PNG image file and saves also an excel file with
% all results.
% IMPORTANT: note that the postprocessing of the raw angle assummes anti-clockwise rotation!!
%  
% NOTE: this function works for short tracks, it only fails if no. of
% points in track is < 3.
%
% 
% INPUTS (these relate to function showManyParticleTrajAnalysis.m): 
% - 'trajXlsPath' is the path to an excel file with the full trajectories as returned by
% function "linkTrajSegmentsParticles.m".
% - 'image_data': structure array with each element being a 2D matrix with 
% the image data for each frame in the image sequence. 
% - 'analysedAllTraj': trajectory information. This is generated within function showManyParticleTrajAnalysis.m. 
% It is a structure array with as many elements as
% analysed trajectories (the ones with enough points in them), and
% with fields:
% 'trajNumber','AngleDegrees','majorAxisLength','minorAxisLength','xvalues','yvalues','mean_xvalue','mean_yvalue','xvalues_offset',
% 'yvalues_offset','msd_unavg','frame','timeabs','timerel','numel','minNumPointsInTraj',
% 'deltaTime','msd','errorMsd','errorMsdRelPercent','disp'.
% - 'n_traj' refers to the number of trajectory out of all the analysed
% trajectories for which the results are shown.
% - 'start_frame' is the number of frame considered as the origin of time (as t=0), in frames. 
% - 'tsamp' is the sampling time, used to calibrate the absolute time, to go
% from frames to time in seconds. 
% It is the time between frames in seconds. Use tsamp = 1 for the time to be in units of frames. A proper calibration
% have tsamp = 40*10^(-3), i.e., 40ms per frame, for example.
% start_frame*tsamp is therefore the absolute time origin in seconds.
% - 'pixelsize_nm': pixel size in nm (35.333nm for OXPHOS data).
%
% OUTPUT: 
% Saves all plots in a .PNG image file and saves also an excel file with
% all results.
% 'processedTraj' is a cell array 6 elements:
% The first element, {1}, and "Track info" tab in output Excel file,
% contains summary results and is a structure with fields:
%             track_with_jumps_flag
%                good_tracking_flag
%                  short_track_flag
%                   long_track_flag
%              very_long_track_flag
%        max_NumFramesForShortTrack
%         min_NumFramesForLongTrack
%     min_NumFramesForVeryLongTrack
%                           XLSpath
%                minNumPointsInTraj
%                           TrajNum
%                   OriginalTrajNum
%                    FirstTrajFrame
%                     LastTrajFrame
%                     NumDataPoints
%                     TrajStartTime
%                       TrajEndTime
%                      TrajDuration
%                 TimeBetweenFrames
%                FrameForTimeOrigin
%                     AbsTimeOrigin
%                      pixelsize_nm
%                       Track_meanX
%                       Track_meanY
%                Track_meanX_offset
%                Track_meanY_offset
%
% The next element, {2}, and "Track data" tab in output Excel file, is a
% structure with the following fields, each containing a vector (columns in
% Excel file):
%     'frame'
%     'timeabs'
%     'angleDegrees' % original angle found, between -90 and 90 degrees.
%     'angleDegreesPos' % positive, cyclic, processed angle, can be >360 degrees.
%     'majorAxisLength' % length of mayor axis of ellipse fitted to particle.
%     'minorAxisLength' % length of minor axis of ellipse fitted to particle.
%     'xvalues' % x position of found centre of mass
%     'yvalues' % y position of found centre of mass
%     'xvalues_offset' % same as above but relative to first point in track.
%     'yvalues_offset'
%
% The next element, {3}, and "MSD results error<150%" tab in output Excel
% file, is a structure with the following fields, each containing a vector,
% for MSD data with relative error bars < 150%:
%     'deltaTime'
%     'msd'
%     'errorMsd'
%     'errorMsdRelPercent'
%
% The next element, {4}, and "MSD results all" tab in output Excel
% file, is a structure with the following fields, each containing a vector, for all MSD data: 
%     'deltaTime'
%     'msd'
%     'errorMsd'
%     'errorMsdRelPercent'
%
% The next element, {5}, and "Mobility results" tab in output Excel
% file, is a structure with the field 'results_mobility'
% (results of MSD fitting to line and to confined curve). 
% Contains numbers:
%     lengthMsdVector
%     fit_msd_line_offset
%     fit_msd_line_offset_stDev
%     fit_msd_line_slope
%     fit_msd_line_slope_stDev
%     fit_msd_line_rsq
%     fit_msd_conf_limit
%     fit_msd_conf_limit_stDev
%     fit_msd_conf_timeconst
%     fit_msd_conf_timeconst_stDev
%     fit_msd_conf_offset
%     fit_msd_conf_offset_stDev
%     fit_msd_conf_rsq
%     diffusion1D_coeff_micronsSqrdPerSec
%     diffusion1D_coeff_error
%     size_1Dconfin_nm
%     size_1Dconfin_nm_error
%     Brownian_flag
%     Confined_flag
%     OtherMobility_flag
% 
% % The next element, {6}, and "Angular veloc results" tab in output Excel
% file, is a structure with the following fields, each containing a number,
% resulting from a linear fit of the angle versus time:
%     fit_angle_offset, in degrees
%     fit_angle_offset_stDev, in degrees
%     fit_angle_slope, angular velocity in degrees/s
%     fit_angle_slope_stDev, error of angular velocity in degrees/s
%     fit_angle_rsq, r-squared of linear fit.
%     
%
% Note: function "showManyParticleTrajAnalysis.m" calls this function "showParticleTrajAnalysis.m".
%
%

%% PARAMETERS

% PARAMETER mobility_rsq_limit: value of rsquare of fit above which we
% accept that the msd versus delta-time fit is either Brownian diffusion
% (linear) or confined diffusion:
mobility_rsq_limit = 0.4;
% Maximum timeconstant in seconds from confined-trajectory fit (if time
% constant is too large, the fit is actually linear...) for trajectory to
% be labelled as confined diffusion:
max_conf_timeconst = 100; 
% Note: the guesses for the fits of the msd to a line or to a saturating
% curve are given later and their values can make the fits fair or not.

% PARAMETER nbins: number of bins in histogram of intensity pair-wise differences.
% See fullPwD.m.
nbins = 500; % (50, 200)

% PARAMETERS for flags:
% Max number of frames in track for it to be flagged as "short" track:
max_NumFramesForShortTrack = 20; 
% Number of frames in track for it to be flagged as "long" track:
min_NumFramesForLongTrack = 50; 
% Number of frames in track for it to be flagged as "very long" track:
min_NumFramesForVeryLongTrack = 120; 

% % ------------------
% % Alternatively, one can save all parameter values in a file, e.g., "paramsForShowTrajAnalysis2.m"
% % in the current directory. Just calling the name of that file loads the
% % parameter values into the workspace:
% paramsForShowTrajAnalysis2
% % In this way, we can save different parameter sets for different data sets
% % in an easy manner and run two Matlabs at the same time working with different parameter sets.
 

%% Display some info first:

% Some needed track parameters-info:
firstFrameInTraj = analysedAllTraj(n_traj).frame(1); % first frame in trajectory.
lastFrameInTraj = analysedAllTraj(n_traj).frame(end); % last frame in trajectory.
traj_duration = tsamp*(lastFrameInTraj-firstFrameInTraj); % duration of track in seconds.
nPointsInTrack = analysedAllTraj(n_traj).numel;

% t_origin is the absolute time origin in seconds:
t_origin = start_frame*tsamp;

disp(' ') % empty line
disp(['The analysed file is: ',trajXlsPath]) 
disp(['The number of long enough trajectories in the file is: ',num2str(length(analysedAllTraj))]) 
disp(['The time origin for this image sequence is frame no.: ',num2str(start_frame)]) 
disp(' ') % empty line
disp(['TRAJECTORY number: ',num2str(n_traj)]) 
disp(['The number of data points in this trajectory is: ',num2str(nPointsInTrack)]) 


% if n_traj==50
%     n_traj
% end

%% Create new directory for saving trajectory-analysis results for this image sequence:

% Make new folder (new directory) to save trajectory analysis results:
pos1 = strfind(trajXlsPath,'fullTrajs.xls'); % position of the start of the string 'fullTraj.xls' in the xls input file name.
new_folder_name = trajXlsPath(1:(pos1-1)); % Take the name of the input excel file (with the end bit 'fullTraj.xls' removed) as the new folder name.
warning('off','MATLAB:MKDIR:DirectoryExists'); % Turn off warning: "Warning: Directory already exists." .
mkdir(new_folder_name); % make new directory.


%% Read in the image-sequence data:

% Read image-sequence file: 
%[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.
% --------------------------------------------------------------

%% Start multifigure here to plot/display all results:

% Create figure. 'position' vector is [left, bottom, width, height]. 
h1 = figure('Tag','Trajectory results','color','white','units','inches','position',[4 2.7 16 9]); 
% left, bottom control the position at which the window appears when it
% pops. 

time_abs = analysedAllTraj(n_traj).timeabs;
x_values_offset = analysedAllTraj(n_traj).xvalues_offset;
y_values_offset = analysedAllTraj(n_traj).yvalues_offset;

% Plot xvalues and yvalues (origin is first point in trajectory) versus time:
subplot(2,3,1);
plot(time_abs,x_values_offset,'.-b'); 
hold on;
plot(time_abs,y_values_offset,'.-g');
xlabel('t_{abs} (s)'); 
ylabel('x-x0 (b),  y-y0 (g) (pix)');
xlim([0 max(analysedAllTraj(n_traj).timeabs)]);
% legend('x-x0 (pix)','y-y0 (pix)');
legend('hide');
% display image file path:
axes_limits = axis; % get axes limits: axes_limits = [xmin xmax ymin ymax].
text(axes_limits(1),1.2*axes_limits(4),trajXlsPath,'Interpreter','none'); % display path of image sequence at top left of figure at position (xmin,ymax) of axes.
hold off;

% Plot trajectory on x-y plane (origin is first point in trajectory), axes in nm now:
subplot(2,3,2);
plot(pixelsize_nm.*analysedAllTraj(n_traj).xvalues_offset,pixelsize_nm.*analysedAllTraj(n_traj).yvalues_offset,'.-k');
xlabel('x-x0 (nm)'); 
ylabel('y-y0 (nm)'); 


%% Orientation angle versus time
angleDeg = analysedAllTraj(n_traj).AngleDegrees;
angleDegPos = angleDegToPos(angleDeg); % orientation angle by fitting ellipse.
% angleDeg is the angle originally found, between -90 and 90 degrees.
% angleDegPos is a postprocessed positive cyclic angle that can be larger than
% 360 degrees if there are many turns of the particle. The postprocessing
% within function angleDegToPos assummes anti-clockwise rotation.

% Fit orientation angle versus time to a LINE to get angular velocity of
% particle rotation:
linearFunction = fittype('A + B*x','independent','x'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% CHECK: values of guesses:
guess_A = 0; % offset guess in degrees.
guess_B = 30; % slope guess in degrees/s when correct frame rate has been entered.
options.StartPoint = [guess_A guess_B];
try % error control
    [fit_angle gof] = fit(time_abs,angleDegPos,linearFunction,options);
    fit_param_values = coeffvalues(fit_angle); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_msd_line) to find out order of parameters.
    fit_angle_offset = fit_param_values(1); % offset from fit in deg.
    fit_angle_slope = fit_param_values(2); % slope from fit in deg/s.
    fit_angle_rsq = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_angle,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    fit_angle_offset_stDev = errorSTDEV(1);
    fit_angle_slope_stDev = errorSTDEV(2);
catch ME1 % catch error. If there was an error:
    fit_angle = [0]; % needed for plot later.
    fit_angle_offset = []; 
    fit_angle_offset_stDev = []; 
    fit_angle_slope = []; 
    fit_angle_slope_stDev = []; 
    fit_angle_rsq = []; 
end

disp(' ') % empty line
disp('Fit for orientation angle versus time to a line:') 
disp([' fit_angle_offset = ',num2str(fit_angle_offset),' +- ',num2str(fit_angle_offset_stDev),';   fit_angle_slope = ',num2str(fit_angle_slope),' +- ',num2str(fit_angle_slope_stDev)]) 
disp([' fit_angle_rsq = ',num2str(fit_angle_rsq)]) % r squared of fit.

% Orientation results for angle vs time (numbers):
results_angle.fit_angle_offset = fit_angle_offset; 
results_angle.fit_angle_offset_stDev = fit_angle_offset_stDev;
results_angle.fit_angle_slope = fit_angle_slope;
results_angle.fit_angle_slope_stDev = fit_angle_slope_stDev;
results_angle.fit_angle_rsq = fit_angle_rsq;

% Angular velocity RESULTS from fits (numbers):
n_angle = 6;
processedTraj{n_angle} = results_angle;

% Plot angle versus time:
subplot(2,3,3);
plot(time_abs,angleDeg,'.-b')
hold on;
plot(time_abs,angleDegPos,'.-r')
% Plot angle versus time linear fit:
try % error control, if fit failed, do not plot
    plot(fit_angle,'k'); % plot linear fit to msd data as a black line.
catch ME3
end
xlabel('t_{abs} (s)'); 
ylabel('orientation angle (deg)');
xlim([0 max(analysedAllTraj(n_traj).timeabs)]);
% legend('angle(deg)','cyclicAngle(deg)');
legend('hide');
hold off;


%% Msd (mean square displacement) versus delta-time:

% Original msd data:
% insert break point here for figures for paper
xdata_msd_0 = double(analysedAllTraj(n_traj).deltaTime); % delta-time is in seconds.
ydata_msd_0 = double(analysedAllTraj(n_traj).msd); % fits fail later if data is not of class double. msd is in pixels^2 here.
% Error bars for msd versus Delta t plot:
halfErrorBars_0 = double(analysedAllTraj(n_traj).errorMsd); % These are the lower and upper bounds of the error bars around the msd result points for the graph, the standard deviation.
msd_relative_errors_0 = double(analysedAllTraj(n_traj).errorMsdRelPercent);

% % -----------
% Aid plot:
% Plot all Msd data:
% subplot(2,3,3);
% errorbar(xdata_msd_0,ydata_msd_0,halfErrorBars_0,halfErrorBars_0,'.r'); 
% xlabel('\Deltat (s)'); 
% ylabel('msd (pix^2)'); 
% xlim([0 max(xdata_msd_0)]); 
% ylim([0 1.5*max(ydata_msd_0)]);  
% title('All msd (2D, x-y) data with error bars')

% Only fit data points with relative error < 150%:
% % Use warning('query','last') to find out error identifier:
% warning('off','MATLAB:ignoreImagPart'); % turn off warning message for Warning: "Imaginary part ignored in comparison operation". 
% Select only real numbers (for large errors, msd_relative_errors is Inf and values are imaginary...)
pos_acceptedError_0 = []; % initialise empty vector.
for i=1:length(msd_relative_errors_0)
    if isreal(msd_relative_errors_0(i))==1        
        pos_acceptedError_0 = [pos_acceptedError_0  i]; % positions of real numbers in vector.
    end
end
msd_relative_errors_bis = msd_relative_errors_0(pos_acceptedError_0); % vector of relative errors excluding Inf and imaginary numbers. 
xdata_msd_bis = xdata_msd_0(pos_acceptedError_0);
ydata_msd_bis = ydata_msd_0(pos_acceptedError_0);
halfErrorBars_bis = halfErrorBars_0(pos_acceptedError_0);
% Now take only values with relative error < 150% :
pos_acceptedError = find(msd_relative_errors_bis < 150);
xdata_msd = xdata_msd_bis(pos_acceptedError);
ydata_msd = ydata_msd_bis(pos_acceptedError);
halfErrorBars = halfErrorBars_bis(pos_acceptedError);
msd_relative_errors = msd_relative_errors_bis(pos_acceptedError);

% -----------------------
% Fit the MSD versus delta-time to a LINE (Brownian-linear mobility):
linearFunction = fittype('A + B*x','independent','x'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% CHECK: values of guesses:
guess_A = 0; % offset guess in pix^2.
guess_B = 5; % slope guess in pix^2/s.
options.StartPoint = [guess_A guess_B];
try % error control
    [fit_msd_line gof] = fit(xdata_msd,ydata_msd,linearFunction,options);
    fit_param_values = coeffvalues(fit_msd_line); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_msd_line) to find out order of parameters.
    fit_msd_line_offset = fit_param_values(1); % offset from fit in pix^2.
    fit_msd_line_slope = fit_param_values(2); % slope from fit in pix^2/s.
    fit_msd_line_rsq = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_msd_line,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    fit_msd_line_offset_stDev = errorSTDEV(1);
    fit_msd_line_slope_stDev = errorSTDEV(2);
catch ME1 % catch error. If there was an error:
    fit_msd_line = [0]; % needed for plot later.
    fit_msd_line_offset = []; 
    fit_msd_line_offset_stDev = []; 
    fit_msd_line_slope = []; 
    fit_msd_line_slope_stDev = []; 
    fit_msd_line_rsq = []; 
end

disp(' ') % empty line
disp('Fit for msd versus delta-time to a line:') 
disp([' fit_msd_line_offset = ',num2str(fit_msd_line_offset),' +- ',num2str(fit_msd_line_offset_stDev),';   fit_msd_line_slope = ',num2str(fit_msd_line_slope),' +- ',num2str(fit_msd_line_slope_stDev)]) 
disp([' fit_msd_line_rsq = ',num2str(fit_msd_line_rsq)]) % r squared of fit.

% Mobility results for msd vs delta-time (numbers):
results_mobility.lengthMsdVector = length(ydata_msd); % no. of points in msd vector with error < 150%.
results_mobility.fit_msd_line_offset = fit_msd_line_offset; 
results_mobility.fit_msd_line_offset_stDev = fit_msd_line_offset_stDev;
results_mobility.fit_msd_line_slope = fit_msd_line_slope;
results_mobility.fit_msd_line_slope_stDev = fit_msd_line_slope_stDev;
results_mobility.fit_msd_line_rsq = fit_msd_line_rsq;

% -----------------------
% Fit the MSD versus delta-time to a SATURATING CURVE (confined mobility):
% As before, only fit data points with relative error < 150%, i.e., xdata_msd and ydata_msd:
% Note: min. no. of points in track must be at least 6 to get at least 4
% points in msd so that a fit with 3 params works ok, otherwise, it fails.
saturating_function = fittype('C+A*(1-exp(-t/B))','independent','t'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% CHECK: values of guesses:
guess_A = 3; % saturation limit guess (pix^2).
guess_B = 0.2; % time-constant guess in seconds.
guess_C = 0; % offset guess (pix^2).
options.StartPoint = [guess_A guess_B guess_C];
% Error control: catch error message to avoid exiting the whole function:
try
    [fit_msd_conf gof] = fit(xdata_msd,ydata_msd,saturating_function,options);
    fit_param_values = coeffvalues(fit_msd_conf); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_msd_conf) to find out order of parameters.
    fit_msd_conf_limit = fit_param_values(1); 
    fit_msd_conf_timeconst = fit_param_values(2); 
    fit_msd_conf_offset = fit_param_values(3);
    fit_msd_conf_rsq = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_msd_conf,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    fit_msd_conf_limit_stDev = errorSTDEV(1);
    fit_msd_conf_timeconst_stDev = errorSTDEV(2);
    fit_msd_conf_offset_stDev = errorSTDEV(3);
catch ME1 % MessageError1. % If fit fails:
    try
        % CHECK: second guesses:
        options.StartPoint = [max(ydata_msd) guess_B guess_C]; % use a different guess for the saturating limit.
        [fit_msd_conf gof] = fit(xdata_msd,ydata_msd,saturating_function,options);
        fit_param_values = coeffvalues(fit_msd_conf); % parameter values resulting from fit. First one is 'A', second one is 'B'.
        fit_msd_conf_limit = fit_param_values(1); % offset from fit. This is the initial intensity, I at t_origin or tstart.
        fit_msd_conf_timeconst = fit_param_values(2); % slope from fit.
        fit_msd_conf_offset = fit_param_values(3);
        fit_msd_conf_rsq = gof.rsquare; % rsquare coefficient of fit.
        errors = confint(fit_msd_conf,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
        errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
        fit_msd_conf_limit_stDev = errorSTDEV(1);
        fit_msd_conf_timeconst_stDev = errorSTDEV(2);
        fit_msd_conf_offset_stDev = errorSTDEV(3);
    catch ME2 % if second try of fit fails too, set all results to zero:
        disp('Fit of msd versus Delta-time to a saturating curve failed.');
        
        fit_msd_conf_limit = 0;
        fit_msd_conf_timeconst = 0;
        fit_msd_conf_offset = 0;
        fit_msd_conf_rsq = 0; % this is safe for later mobility flags too.
        fit_msd_conf_limit_stDev = 0;
        fit_msd_conf_timeconst_stDev = 0;
        fit_msd_conf_offset_stDev = 0;
        
        disp(ME2.message);
        ME2 = addCause(ME2, ME1);
        % rethrow(ME2)
    end
end

disp(' ') % empty line
disp('Fit for msd versus delta-time to a saturating curve:') 
disp([' fit_msd_conf_limit = ',num2str(fit_msd_conf_limit),' +- ',num2str(fit_msd_conf_limit_stDev),';   fit_msd_conf_timeconst = ',num2str(fit_msd_conf_timeconst),' +- ',num2str(fit_msd_conf_timeconst_stDev),';   fit_msd_conf_offset = ',num2str(fit_msd_conf_offset),' +- ',num2str(fit_msd_conf_offset_stDev)]) 
disp([' fit_msd_conf_rsq = ',num2str(fit_msd_conf_rsq)]) % r squared of fit.

% Mobility results for msd vs delta-time (numbers):
results_mobility.fit_msd_conf_limit = fit_msd_conf_limit; 
results_mobility.fit_msd_conf_limit_stDev = fit_msd_conf_limit_stDev;
results_mobility.fit_msd_conf_timeconst = fit_msd_conf_timeconst;
results_mobility.fit_msd_conf_timeconst_stDev = fit_msd_conf_timeconst_stDev;
results_mobility.fit_msd_conf_offset = fit_msd_conf_offset;
results_mobility.fit_msd_conf_offset_stDev = fit_msd_conf_offset_stDev;
results_mobility.fit_msd_conf_rsq = fit_msd_conf_rsq;

% ----------------
% Plot only msd data which has relative error < 150% :
subplot(2,3,4);
errorbar(xdata_msd,ydata_msd,halfErrorBars,halfErrorBars,'.b'); 
try % error control, if setting limits fails, do not set:
    xlim([0 max(xdata_msd)]);
    ylim([-0.5*max(ydata_msd) 2.5*max(ydata_msd)]);
catch ME4
    disp(ME4.message);
end
hold on;
% Plot msd fits:
plot(fit_msd_line,'k'); % plot linear fit to msd data as a black line.
try % error control, if fit to saturating curve failed, do not plot:
    plot(fit_msd_conf,'r'); % plot saturating-curve fit to msd data as a red line.
catch ME5 % catch error to avoid getting out of the entire function which is runnings
    disp(ME5.message);
end
legend('hide');
xlabel('\Deltat (s)'); 
ylabel('msd (pix^2)'); 
title('msd data with error bars < 150%')
hold off;

% multifigure continues to be filled in later...

% -------------
% Get values of size of 1D confinement region and Diffusion coefficient in
% meaningful units.
% For the "Brownian" trajectories: from the slope of msd vs Delta-t I get the
% 1D Diffusion coefficient as slope/2 and use the pixel size to convert to
% microns^2/s. 
% The slope from the fit of msd vs deltat-time is in pix^2/s.
diffusion1D_coeff_micronsSqrdPerSec = (fit_msd_line_slope/2)*(pixelsize_nm/1000)^2; % 1D diffusion coeff in microns^2/s.
diffusion1D_coeff_error = (fit_msd_line_slope_stDev/2)*(pixelsize_nm/1000)^2; % error or standard deviation of the 1D diffusion coeff in microns^2/s.
disp(' ') % empty line
disp(['Diffusion1D_coeff_micronsSqrdPerSec = ',num2str(diffusion1D_coeff_micronsSqrdPerSec),' +- ',num2str(diffusion1D_coeff_error)])

% For the "confined" trajectories, I assume confined diffusion in 1D to a
% segment of length L. The msd saturates to a value of (L^2)/6. I use the
% pixel size to convert to nm. From the fit, the value the msd saturates to
% is (A) fit_msd_conf_limit (pixels^2). 
size_1Dconfin_nm = pixelsize_nm*sqrt(6*fit_msd_conf_limit); % size of 1D confinement region in nm (L).
% L = pixtonm*sqrt(6)*sqrt(A), where A =fit_msd_conf_limit and L = size_1Dconfin_nm. 
% Propagating errors we have that error_L = (L/2)*(error_A/A):
size_1Dconfin_nm_error = (size_1Dconfin_nm/2)*(fit_msd_conf_limit_stDev/fit_msd_conf_limit); % size of 1D confinement region in nm.
disp(' ') % empty line
disp(['size_1Dconfin_nm = ',num2str(size_1Dconfin_nm),' +- ',num2str(size_1Dconfin_nm_error)])

% ----------------------
% Mobility flags:
if fit_msd_line_rsq > mobility_rsq_limit
    Brownian_flag = 1;
else 
    Brownian_flag = 0;
end

if (fit_msd_conf_rsq > mobility_rsq_limit) && (fit_msd_conf_timeconst < max_conf_timeconst) && (fit_msd_conf_timeconst < 1.5*traj_duration)
    % Reject fits to confined trajectory if timeconstant of fit is >= max_conf_timeconst
    % seconds (that is really a linear fit):
    Confined_flag = 1;
else
    Confined_flag = 0;
end

% (Brownian_flag == 1 && Confined_flag == 1) ||
if  Brownian_flag == 0 && Confined_flag == 0
    OtherMobility_flag = 1; % this results in mobility flags [1 1 1] or [0 0 1]
else
    OtherMobility_flag = 0;
end

mobility_flags = [Brownian_flag Confined_flag OtherMobility_flag];

disp(' '); % empty line.
disp('Mobility flags [Brownian_flag  Confined_flag  OtherMobility_flag] are now: ');
disp(mobility_flags);
% CHECK: give user option to override Mobility Flags or not:
% Give user a chance to override mobility flags:
% agree_mobility_flags = input('Do you agree with previous mobility_flags? (1 for yes, 0 for no): '); % request user input
agree_mobility_flags = 1;
if agree_mobility_flags == 0 % if user does not agree with mobility flags:
    mobility_flags = input('Enter new mobility_flags as [traj_Brownian  traj_confined  traj_other] (ones and zeros):'); % request user input.
    Brownian_flag = mobility_flags(1);
    Confined_flag = mobility_flags(2);
    OtherMobility_flag = mobility_flags(3);
end

% ------------------------
% Add to mobility results (numbers):
results_mobility.diffusion1D_coeff_micronsSqrdPerSec = diffusion1D_coeff_micronsSqrdPerSec;
results_mobility.diffusion1D_coeff_error = diffusion1D_coeff_error;
results_mobility.size_1Dconfin_nm = size_1Dconfin_nm;
results_mobility.size_1Dconfin_nm_error = size_1Dconfin_nm_error;
results_mobility.Brownian_flag = Brownian_flag;
results_mobility.Confined_flag = Confined_flag;
results_mobility.OtherMobility_flag = OtherMobility_flag;



%% Plot frame 1 with all trajectory points:

% Display frame 1 and all trajectory positions overlaid on it:
subplot(2,3,5);
frame1 = image_data(1).frame_data;
imshow(frame1,[],'Border','tight'); % display frame average scaled between its min and max values ([]).
hold on;

x_values = analysedAllTraj(n_traj).xvalues;
y_values = analysedAllTraj(n_traj).yvalues;

for k = 1:length(x_values)
    plot(x_values(k),y_values(k),'.-','Color','g','MarkerSize',5) % plot track positions in green.
end
hold off;


%% Continue adding results to figure:

% % Reminder:
% firstFrameInTraj = analysedAllTraj(n_traj).frame(1); % first frame in trajectory.
% lastFrameInTraj = analysedAllTraj(n_traj).frame(end); % last frame in trajectory.
% traj_duration = tsamp*(lastFrameInTraj-firstFrameInTraj); % duration of track in seconds.
% nPointsInTrack = analysedAllTraj(n_traj).numel;

% Display track info:
subplot(2,3,6); axis off; 
% text(-0.07,1,trajXlsPath,'Interpreter','none'); % display path of image sequence at position (0,1) of the axes (when there is no specific plot, axes go from 0 to 1).
str1(1) = {['Plots shown are for TRAJECTORY number ',num2str(n_traj)]};
str1(2) = {['No. of data points in this trajectory: ',num2str(nPointsInTrack)]};
str1(3) = {['First frame in this trajectory is:  ',num2str(firstFrameInTraj)]};
str1(4) = {['Time origin for sequence is frame no. ',num2str(start_frame)]};
str1(5) = {['Time between frames in seconds is: ',num2str(tsamp)]};
str1(6) = {['No. trajectories analysed: ',num2str(length(analysedAllTraj))]};
str1(7) = {['original TRAJECTORY no.: ',num2str(analysedAllTraj(n_traj).trajNumber)]};
str1(8) = {['Angular velocity of particle from fit: ',num2str(fit_angle_slope),' +- ',num2str(fit_angle_slope_stDev),' degrees/s.']};


text(0,0.5,str1)



%% ADD FLAGS:
% Add flags to first element of cell array of results:

n_trackInfo = 1; % cell element in which to save track info (single numbers).
figure(h1); % bring figure to front.

% Flag indicating if all points in track are equally spaced or if there are
% jumps in track (due to linking segments up to a certain no. of frames away):
if traj_duration > tsamp*(nPointsInTrack-1)
    track_with_jumps_flag = 1;
else
    track_with_jumps_flag = 0;
end

% Good-tracking flag:
% only files which show "good tracking" (on video) in showManyTrajAnalysis.m are analysed, and their results saved.
good_tracking_flag = 1;

% Flag trajectories as "short" (see PARAMETERS section):
if nPointsInTrack < max_NumFramesForShortTrack
    short_track_flag = 1;
else
    short_track_flag = 0;
end

% Flag trajectories as "long" or "very long" (see PARAMETERS section):
if nPointsInTrack < min_NumFramesForLongTrack
    long_track_flag = 0;
    very_long_track_flag = 0;
else
    if nPointsInTrack < min_NumFramesForVeryLongTrack
        long_track_flag = 1;
        very_long_track_flag = 0;
    else
        long_track_flag = 0;
        very_long_track_flag = 1;
    end    
end


% Save flags into TRACK INFO results:
processedTraj{n_trackInfo}.track_with_jumps_flag = track_with_jumps_flag;
processedTraj{n_trackInfo}.good_tracking_flag = good_tracking_flag; 
processedTraj{n_trackInfo}.short_track_flag = short_track_flag; 
processedTraj{n_trackInfo}.long_track_flag = long_track_flag;
processedTraj{n_trackInfo}.very_long_track_flag = very_long_track_flag;

% Save deciding criteria for flags into TRACK INFO results:
processedTraj{n_trackInfo}.max_NumFramesForShortTrack = max_NumFramesForShortTrack; 
processedTraj{n_trackInfo}.min_NumFramesForLongTrack = min_NumFramesForLongTrack; 
processedTraj{n_trackInfo}.min_NumFramesForVeryLongTrack = min_NumFramesForVeryLongTrack;


% Display flags on command window:
disp('Trajectory flags are now: ');
disp(processedTraj{n_trackInfo}); 
figure(h1); % bring figure to front.
% numbers, first column is fieldnames of element {n_trackInfo}, second one is values.
% CHECK: give user chance to override TRAJECTORY FLAGS or not:
% Give user a chance to override trajectory flags:
% agree_track_flags = input('Do you agree with previous Trajectory Flags? (1 for yes, 0 for no): '); % request user input
agree_track_flags = 1;
% if agree_track_flags ~= 1 % if user does not agree with track flags:
%     track_flags = input('Enter new Track Flags as [good_tracking  ] (1 or 0 each):'); % request user input.
%     good_tracking_flag = track_flags(1);
%     track_closeToTimeOrigin_flag = track_flags(4);
%     goodExpFit_flag = track_flags(5);
%     goodExpFit_wo_flag = track_flags(6);
%     
%     % Save new flags into TRACK INFO results:
%     processedTraj{n_trackInfo}.good_tracking_flag = good_tracking_flag;
%     processedTraj{n_trackInfo}.track_closeToTimeOrigin_flag = track_closeToTimeOrigin_flag;
%     processedTraj{n_trackInfo}.goodExpFit_flag = goodExpFit_flag;
%     processedTraj{n_trackInfo}.goodExpFit_wo_flag = goodExpFit_wo_flag;
% end



%% Save graphical analysis results (last figure):

% figName = strcat('figResultsTraj',num2str(n_traj)); % name of figure file to save to.
% set(gcf, 'PaperPositionMode', 'auto')  % Use screen size.
% print('-dpng',figName)  % add -r200 after -dpng to control resolution.

% Move into folder previously created to save traj analysis results and
% SAVE current FIGURE as a .png (at screen size) and then close the figure window:
% Result figure files are saved in a new folder within the same directory where the input excel file was.
figName = strcat(new_folder_name,'_figResults_traj',num2str(n_traj),'.png'); % name of figure file to save to. 
saveFigurePNG(new_folder_name,figName); % See saveFigurePNG.m.


%% Output, analysis results:
% Prepare result data to save them later.

% TRACK INFO continued (numbers):
% First element in output cell array contains general info about the image and the trajectory.
% Output result: 
% n_trackInfo = 1; % cell element in which to save track info. See above.
processedTraj{n_trackInfo}.XLSpath = trajXlsPath;
processedTraj{n_trackInfo}.minNumPointsInTraj = analysedAllTraj(n_traj).minNumPointsInTraj; % min no. of points in trajectory for it to be analysed.
processedTraj{n_trackInfo}.TrajNum = n_traj; % numbers of trajectories analysed (long enough) as chosen in analyseTraj.m.
processedTraj{n_trackInfo}.OriginalTrajNum = analysedAllTraj(n_traj).trajNumber; % Original traj number from all those tracked initially.
processedTraj{n_trackInfo}.FirstTrajFrame = firstFrameInTraj; % First frame in this trajectory.
processedTraj{n_trackInfo}.LastTrajFrame = lastFrameInTraj; % Last frame in this trajectory.
processedTraj{n_trackInfo}.NumDataPoints = analysedAllTraj(n_traj).numel; % no. of points in track.
processedTraj{n_trackInfo}.TrajStartTime = tsamp*analysedAllTraj(n_traj).frame(1); % absolute start time of trajectory (s).
processedTraj{n_trackInfo}.TrajEndTime = tsamp*analysedAllTraj(n_traj).frame(end); % absolute end time of trajectory (s).
processedTraj{n_trackInfo}.TrajDuration = traj_duration; % duration of track in seconds.
processedTraj{n_trackInfo}.TimeBetweenFrames = tsamp;
processedTraj{n_trackInfo}.FrameForTimeOrigin = start_frame; % time origin (as frame number) for this whole image sequence.
processedTraj{n_trackInfo}.AbsTimeOrigin = t_origin;
processedTraj{n_trackInfo}.pixelsize_nm = pixelsize_nm; % pixel size in nm.
processedTraj{n_trackInfo}.Track_meanX = mean(x_values); % mean x position (original) of track on image.
processedTraj{n_trackInfo}.Track_meanY = mean(y_values); % mean y position (original) of track on image.
processedTraj{n_trackInfo}.Track_meanX_offset = mean(x_values_offset); % mean x-x_com position of track with respect to first point in track.
processedTraj{n_trackInfo}.Track_meanY_offset = mean(y_values_offset); % mean y-y_com position of track with respect to first point in track.

% -----------------
% TRACK DATA (vectors):
% Next element of results cell array contains track data:
n_trackData = 2; % cell element in which to save track info.
processedTraj{n_trackData}.frame = analysedAllTraj(n_traj).frame; % vectors.
processedTraj{n_trackData}.timeabs = analysedAllTraj(n_traj).timeabs;
processedTraj{n_trackData}.angleDegrees = analysedAllTraj(n_traj).AngleDegrees; % orientation angle by fitting ellipse.
processedTraj{n_trackData}.angleDegreesPos = angleDegPos; % positive orientation angle in degrees, between 0 and 180 deg.   
processedTraj{n_trackData}.majorAxisLength = analysedAllTraj(n_traj).majorAxisLength; % orientation angle by fitting ellipse.
processedTraj{n_trackData}.minorAxisLength = analysedAllTraj(n_traj).minorAxisLength; % orientation angle by fitting ellipse.
processedTraj{n_trackData}.xvalues = analysedAllTraj(n_traj).xvalues; % original position on image.
processedTraj{n_trackData}.yvalues = analysedAllTraj(n_traj).yvalues; % original position on image.
processedTraj{n_trackData}.xvalues_offset = analysedAllTraj(n_traj).xvalues_offset; % position with respect to 1st point in track.
processedTraj{n_trackData}.yvalues_offset = analysedAllTraj(n_traj).yvalues_offset; % position with respect to 1st point in track.
processedTraj{n_trackData}.area = analysedAllTraj(n_traj).area; % particle area (number of pixels).

% Others:

% -----------------
% MSD RESULTS (vectors):
% Next element of results cell array contains MSD results (points with error < 150%):
n_MSD = 3;
% Msd data points for which error < 150% :
processedTraj{n_MSD}.deltaTime = xdata_msd; % vectors.
processedTraj{n_MSD}.msd = ydata_msd;
processedTraj{n_MSD}.errorMsd = halfErrorBars;
processedTraj{n_MSD}.errorMsdRelPercent = msd_relative_errors;

% Next element of results cell array contains original MSD data (all points):
n_MSD_all = 4;
processedTraj{n_MSD_all}.deltaTime_all = analysedAllTraj(n_traj).deltaTime; % vectors.
processedTraj{n_MSD_all}.msd_all = analysedAllTraj(n_traj).msd;
processedTraj{n_MSD_all}.errorMsd_all = analysedAllTraj(n_traj).errorMsd;
processedTraj{n_MSD_all}.errorMsdRelPercent_all = analysedAllTraj(n_traj).errorMsdRelPercent;

% -----------------
% MOBILITY RESULTS from fits (numbers):
n_mobility = 5;
processedTraj{n_mobility} = results_mobility;
% from fits of data with error < 150%.


%% SAVE excel file with trajectory results:

% Move into folder previously created to save traj analysis results:
cd(new_folder_name); % move into that directory.

% Save results in another excel file in the same folder where the
% analysis plot .png file has been saved, new_folder_name:
output_filename = strcat(new_folder_name,'_traj',num2str(n_traj),'.xls'); % name of excel file to save to.
dataForSheet1 = [fieldnames(processedTraj{n_trackInfo}) struct2cell(processedTraj{n_trackInfo})]; % numbers, first column is fieldnames of element {n_trackInfo}, second one is values.
dataForSheet2 = [fieldnames(processedTraj{n_trackData})'; num2cell(cell2mat(struct2cell(processedTraj{n_trackData})'))]; % vectors
dataForSheet3 = [fieldnames(processedTraj{n_angle}) struct2cell(processedTraj{n_angle})]; % numbers
dataForSheet4 = [fieldnames(processedTraj{n_MSD})'; num2cell(cell2mat(struct2cell(processedTraj{n_MSD})'))]; % vectors
dataForSheet5 = [fieldnames(processedTraj{n_MSD_all})'; num2cell(cell2mat(struct2cell(processedTraj{n_MSD_all})'))]; % vectors
dataForSheet6 = [fieldnames(processedTraj{n_mobility}) struct2cell(processedTraj{n_mobility})]; % number

warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.

xlswrite(output_filename,dataForSheet1,'Track info'); % write data to sheet 'Track info' in excel file.
xlswrite(output_filename,dataForSheet2,'Track data'); % write data to sheet 'Track data' in excel file.
xlswrite(output_filename,dataForSheet3,'Angular veloc results'); % write data to sheet 'Angular veloc results' in excel file.
xlswrite(output_filename,dataForSheet4,'MSD results error<150%'); % write data to sheet 'MSD results error<150%' in excel file.
xlswrite(output_filename,dataForSheet5,'MSD results all'); % write data to sheet 'MSD results all' in excel file.
xlswrite(output_filename,dataForSheet6,'Mobility results'); % write data to sheet 'Mobility results' in excel file.

cd('..'); % go back to previous directory.

% CHECK: skip this pause or not.
% Give time to user to prepare before looking at next track:
% input('Press Enter when ready to look at next track: '); 








