function processedManyTrajs = showManyParticleTrajAnalysis(image_label,n_traj_start,n_traj_end,start_frame,tsamp,pixelsize_nm,showVideo,saveAvi,minPointsTraj) 
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
% This function uses showParticleTrajAnalysis.m, getDisplacement.m, etc.
%
% Analyse all trajectory data for a given image sequence (), from trajectory
% 'n_traj_start' to 'n_traj_end', only for the trajectory numbers in "good_track_nums_image_label.mat", i.e., selected as
% "good ones" either going through videos of all tracks by eye using function
% goThroughTracksVideo(image_label,n_traj_start,n_traj_end,minPointsTraj)
% or generating the mat file by hand.
% The input list of "good track" numbers must be contained in a .mat file in the current directory (see below).
%
% Shows results of trajectory analysis and saves them, can show the
% trajectory overlayed on the image sequence on a video to check or not,
% saves results of trajectory analysis to an excel file (done within
% showParticleTrajAnalysis.m), saves output to a .mat file with name
% ('procManyTraj'+image_label) within a new folder named data_set_label +
% image_label in current directory.
% 
% NOTE: before running this function you should move into a directory which
% contains both the video sequence (labelled by 'image_label') which has previously been analysed with
% 'FindTrajectsParticles.m' and 'linkTrajSegmentsParticles.m' to produce the .xls file which
% contains the trajectory results, which should be in the same directory as the image sequence file.
%
% INPUTS: 
% - 'image_label' string that labels the image sequence under analysis, e.g. '101'.
% - 'n_traj_start': first trajectory we want to analyse and check.
% - 'n_traj_end': last trajectory we want to analyse and check. If the
% string 'end' is entered, we go through to the last analysed trajectory.
% - 'start_frame' is the number of frame considered as the origin of time
% (as t=0), in frames. It is the first frame for which the shutter is fully open and detecting images.
% - 'tsamp' is the sampling time (time between frames in seconds), used to calibrate the absolute time, to go from frames to time in seconds. 
% Use tsamp = 1 for the time to be in units of frames. A proper calibration
% would have tsamp = 40*10^(-3), i.e., 40ms per frame, for example.
% start_frame*tsamp is therefore the absolute time origin in seconds.
% - 'pixelsize_nm': pixel size in nm (e.g. 35.333nm for fluorescence data).
% - 'showVideo' is an input parameter to show a video of the trajectory
% overlaid on the image sequence or not (1 or 0).
% - 'saveAvi' to save or not an .avi video file of the trajectory with
% overlaid particle centre and orientation, shows subarray, not full image (1 or 0).
% - minPointsTraj: Minimum number of data points that a trajectory must have in order to be
% analised (default minPointsTraj = 10, at least 3). 
%
%
% OUTPUT: 
% 'processedManyTrajs' is a cell array with as many elements (processedManyTrajs{n}) as
% analysed trajectories (those with good tracking only).
% The first element within trajectory {n} in the cell array, {n}{1}, is a structure with fields:
% fieldnames(processedManyTrajs{1}{1}):
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
% fieldnames(processedManyTrajs{1}{2}):
%     'frame'
%     'timeabs'
%     'xvalues'
%     'yvalues'
%     'xvalues_offset'
%     'yvalues_offset'
%
%
% fieldnames(processedManyTrajs{1}{3}): MSD data with error bar <150%:
%     'deltaTime'
%     'msd'
%     'errorMsd'
%     'errorMsdRelPercent'
%
% fieldnames(processedManyTrajs{1}{4}): all MSD data:
%     'deltaTime'
%     'msd'
%     'errorMsd'
%     'errorMsdRelPercent'
%
% fieldnames(processedManyTrajs{1}{5}): results of MSD fitting
% structure with the field 'results_mobility'
%
% Eg.  T0 = showManyParticleTrajAnalysis('101',5,7,5,0.04,35.333,1);  
% image '101' in folder, look at analysed trajectories 5 to 7, with frame 5
% being the time origin and 40ms between frames.
% Eg. to show only one trajectory (no. 8 of ATPase-GFP_101fullTrajs.xls, eg) do:
% T0 = showManyParticleTrajAnalysis('101',8,8,5,0.04,35.333,1);
% eg. analyse only traj number 2 for image 500 in data set 'cybD-mCherry-ATPase-GFp':  T500 = showManyParticleTrajAnalysis('500',2,2,5,0.04,35.333,1);
% To get results, do T500{1}{i},   for i from 1 to 8.  
% ------------------------------------


%% PARAMETERS: 

% Set the "quickLook" parameter as 1 if you are just looking quickly at
% "good" trajectories and you don't want to be asked for user input as to
% whether the trajectory is "good" or not (to flag them), and you don't
% want to save the result structure as a .mat in a user specified folder
% either...
quickLook = 0;

% Frame rate to show video and for the saved .avi video file:
framesPerSecond = 10;
% e.g., 25 frames per second corresponds to 40ms between frames.

% Extra factor multiplier to allow for particle movement when displaying
% and saving video file of particle with overlaid centre-of-mass and
% orientation. A subarray is saved around average centre-of-mass
% position of particle, with half size equal to 1.2 times subarray-halfsize
% (parameter of function findTrajectsParticles.m). Corrected later if particle too close
% to edge of image.
extraSubarraySizeFactor = 1.2;

initial_folder_path = cd; 

%% Get path for trajectory data (excel file):

% You need to be in the correct directory before running the function!!!!
% Find paths in current folder which contain 'image_label' string:
trajXlsPath0 = dir(strcat('*',image_label,'*.xls')); % Trajectory data path (excel file with the full trajectories as returned by function "linkTrajSegments.m").
% Error control:
if isempty(trajXlsPath0) % If there is no .xls trajectory data file for such image number, show error and exit function:
    error('Check you are in the correct directory and run again. No .xls file found for that image number. Make sure image number is in between quotes ''.'); 
end
trajXlsPath = trajXlsPath0.name;
% trajXlsPath0 is a structure and the file names is stored in the field 'name'.


%% Create new directory for saving trajectory-analysis result structure for this image sequence:

% Make new folder (new directory) to save trajectory analysis results:
pos1 = strfind(trajXlsPath,'fullTrajs.xls'); % position of the start of the string 'fullTraj.xls' in the xls input file name.
new_folder_name = trajXlsPath(1:(pos1-1)); % Output folder path. Take the name of the input excel file (with the end bit 'fullTraj.xls' removed) as the new folder name.
% Note that the directory "new_folder_name" is created by function
% showTrajAnalysis2.m when called from this function.

%% Get path for .mat file with the numbers of the good tracks:

% A file with a name "good_track_nums_1757.mat" (eg) is produced by
% function "goThroughTracksVideo.m":
good_tracks_path = dir(strcat('*good_track_nums','*',image_label,'*.mat'));
load(good_tracks_path.name); % This loads the structure "good_tracks" onto the workspace.
% good_tracks is a structure with fields:
% good_tracks.image_label = image_label; input value.
% good_tracks.n_traj_start = n_traj_start; input value.
% good_tracks.n_traj_end = n_traj_end; input value.
% good_tracks.minPointsTraj = minPointsTraj, input value.
% good_tracks.track_numbers: is a row vector containing the numbers of tracks considered as
% "good" by the user after seing the videos of the track overlaid on the
% image sequence.


%% Read in all trajectory data (excel file with all tracks): 
% (.xls file previously generated with functions 'FindTrajectsParticles' and 'linkTrajSegmentsParticles'):

% error control:
if 2~=exist('xlsread') 
    error('Check file dependencies - you need to install xlsread'); 
end

% Open Excel file and read the data:

[NUMERIC0,TXT0,RAW0]=xlsread(trajXlsPath,'params FindTrajects'); % import the data in the sheet named 'params FindTrajectsParticles'.
pos_subarraySize = find(strcmp(TXT0,'subarray_halfwidth')); % row number corresponding to 'subarray_halfwidth' value.
subarray_halfwidth = RAW0{pos_subarraySize,2}; % extract subarray_halfwidth to later use to plot/save video of only a subarray.

[NUMERIC,TXT,RAW]=xlsread(trajXlsPath,'Track results'); % import the data in the sheet named 'Track results'.
% Import the column heads and assign ID
colheads = TXT;
% The column titles are: CentreX, CentreY, ClipFlag, rsqFitX, rsqFitY,
% FrameNumber, ParticleNumber,	TrajNumber.

% Generate ID: ID is a structure with fiels with the same names as the
% column titles, and each has an ID value of 1, 2, 3, etc (see below).
for i=1:numel(colheads) 
    ID.(colheads{i}) = find(strcmp(TXT,colheads{i})); 
end
% eg. ID = 
%     estimateXcentre: 98
%     estimateYcentre: 391
%                Xcom: 98.4190
%                Ycom: 391.4190
%        AngleDegrees: 63.1974
%     majorAxisLength: 13.0164
%     minorAxisLength: 10.4070
%            ClipFlag: 0
%      TooCloseToEdge: 0
%         FrameNumber: 1
%      ParticleNumber: 1
%          TrajNumber: 1

% The trajectory number column:
traj = NUMERIC(:,ID.TrajNumber); % NUMERIC is the numeric data read from the excel file (without the row of column titles).
disp('Excel file read successfully.');

% Get individual tracks:

% List the points at which the trajectory number first appears:
[A,I,J] = unique(traj,'first'); % [A,I,J] = UNIQUE(traj,'first') returns the vector I to index the first occurrence of each unique value in traj.  
% A has the same values as in traj but with no repetitions. A will also be sorted.
% List the points at which the trajectory number last appears:
[A,Y,Z] = unique(traj,'last'); % UNIQUE(traj,'last'), returns the vector Y to index the last occurrence of each unique value in traj.

% Get the number of tracks (no. of different trajectories):
numtracks = numel(A);

% Create tracks structure:
tracks(1:numtracks) = struct('trajNumber',[],'angleDegrees',[],'AngleDegreesPos',[],...
    'majorAxisLength',[],'minorAxisLength',[],...
    'xvalues',[],'yvalues',[],'mean_xvalue',[],'mean_yvalue',[],...
    'xvalues_offset',[],'yvalues_offset',[],...
    'msd_unavg',[],'frame',[],'timeabs',[],'timerel',[],'numel',[],...
    'minNumPointsInTraj',[],'deltaTime',[],'msd',[],'errorMsd',[],...
    'errorMsdRelPercent',[],'disp',[]);

del = []; % initialise for later.

for i=1:numtracks 
    % i
    a = I(i); % index for starting point in trajectory.
    b = Y(i); % index for ending point in trajectory.
    
    % Delete tracks that are less than minPointsTraj data points long:
    if b-a+1 >= minPointsTraj  % Only analyse tracks which have at least "minPointsTraj" points (frames) in them (5, or 15, e.g.).
    
    data{i} = NUMERIC(a:b,:);
    % tracks(i).XLS.track_index = A(i);
    tracks(i).trajNumber = A(i);
    tracks(i).AngleDegrees = data{i}(1:end,ID.AngleDegrees); % orientation angle in degrees. Original angle found, between -90 and 90 degrees.   
    tracks(i).AngleDegreesPos = angleDegToPos(tracks(i).AngleDegrees); % positive orientation angle in degrees, between 0 and 180 deg.   
    tracks(i).majorAxisLength = data{i}(1:end,ID.majorAxisLength); % length of major axis of ellipsoid fitted to particle shape.
    tracks(i).minorAxisLength = data{i}(1:end,ID.minorAxisLength); % length of major axis of ellipsoid fitted to particle shape.
    tracks(i).numel = b-a+1; % Number of points in the track.
    tracks(i).frame = data{i}(1:end,ID.FrameNumber); % frame numbers.
    tracks(i).timeabs = data{i}(1:end,ID.FrameNumber).*tsamp; % Absolute time in seconds. tsamp is the time between frames in s.
    tracks(i).minNumPointsInTraj = minPointsTraj;
    % All absolute position values in pixels:
    tracks(i).xvalues = data{i}(1:end,ID.Xcom); % xvalues for centre of mass of particle on image (used later for plotting traj on image).
    tracks(i).yvalues = data{i}(1:end,ID.Ycom); % yvalues for centre of mass of particle on image.           
    tracks(i).mean_xvalue = mean(data{i}(1:end,ID.Xcom)); % mean x value. 
    tracks(i).mean_yvalue = mean(data{i}(1:end,ID.Ycom)); % mean y value.
    % Set spatial origin to position of first point in track:
    tracks(i).xvalues_offset = tracks(i).xvalues - (tracks(i).xvalues(1)); % xvalues relative to the first point in the trajectory.
    tracks(i).yvalues_offset = tracks(i).yvalues - (tracks(i).yvalues(1)); % % yvalues relative to the first point in the trajectory.
    tracks(i).msd_unavg = tracks(i).xvalues_offset.^2+tracks(i).yvalues_offset.^2; % absolute squared displacement from the origin (0,0): x^2 + y^2.
    % Set time origin to first point in track:
    tracks(i).timerel = tracks(i).timeabs-tracks(i).timeabs(1); % Time relative to first point in track. Set the first frame analysed as time zero reference (not used for now). 
    % Calculate and add to structure the 2D mean square displacement (msd):
    tracks(i) = getDisplacement(tracks(i),tsamp); % calculate msd and its error and add it to result structure.
    % getDisplacement.m adds the following fields to the tracks structure:
    % 'deltaTime','msd','errorMsd','errorMsdRelPercent' and 'disp'.
    else
        % save indices to delete later:
        del(i) = i;     
    end
    
end

% Delete tracks which were too short: 
tracks(find(del))=[];

% ========================

%% Analyse all trajectory data (get msd): 
analysedAllTraj = tracks; 

if isempty(analysedAllTraj)
    disp('The total number of long enough trajectories (to analyse) for this file is zero.');
    disp('Exiting program');
    return % exits function.
end

%'The total number of trajectories analysed (long enough) in the file is: 
n_trajs_analysed = length(analysedAllTraj);


if strcmp(n_traj_end,'end') % string compare, true if equal
    % go through all analysed trajectories til the last one: 
   n_traj_end = n_trajs_analysed;
end


%% Error control: IMPORTANT!!
% Check that the "minPointsTraj" (min no. of points in track for it to be
% analysed) and all other inputs for function "goThroughTracksVideo.m" were
% the same as for this function:
if (good_tracks.minPointsTraj ~= minPointsTraj) || ...
        (strcmp(good_tracks.image_label,image_label)~=1) || ...
        (good_tracks.n_traj_start ~= n_traj_start) || ...
        (good_tracks.n_traj_end ~= n_traj_end)
    disp('ERROR: parameters minPointsTraj, image_label, n_traj_start, n_traj_end must be the same as those used to create list of "good track" numbers with function goThroughTracksVideo.m or goThroughTracksVideo2.m. Exiting function...')
    return % exit function.
end
% row vector with  numbers of "good" tracks:
good_track_numbers = good_tracks.good_track_numbers;


%% Read in the image-sequence data:

% Read image-sequence file:
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.



%% Loop selected trajectories:

n_good_tracking = 1; % initialise index for trajs with good tracking which are saved within loop.

for n = n_traj_start:n_traj_end
    
    % Check if track number n is one of the "good" ones:
    B = ismember(good_track_numbers,n); % result is a vector with zeros at all positions except at the position of n in vector good_track_numbers, if it is a "good" one.
    % sum(B) is equal to 0 if "n" is not a "good" track, and equal to 1 if
    % "n" is on the list of "good" track numbers.
    
    if sum(B) == 1 % If track number "n" is a "good" one:
        
        % Close any pre-existing figures:
        close(findobj('Tag','Trajectory results'));
        % close(all); % close all figures.
        
        frames_list = analysedAllTraj(n).frame; % list of frame numbers in trajectory n.
        x_values = analysedAllTraj(n).xvalues; % list of original x centres of particles in trajectory n.
        y_values = analysedAllTraj(n).yvalues; % list of original y centres of particles in trajectory n.
        
        mean_xvalue = analysedAllTraj(n).mean_xvalue; % average x position of particle in track.
        mean_yvalue = analysedAllTraj(n).mean_yvalue; % average y position of particle in track.
        
        angles_deg = analysedAllTraj(n).AngleDegrees; % list of orientation angles in degrees.
        majAxisLength = analysedAllTraj(n).majorAxisLength; % list of lengths of major axis in pixels.
        
        % Show video of trajectory overlaid on actual image:
        if showVideo == 1 % 'showVideo' is input parameter.
                                              
            % Loop through frames in each trajectory analysed:
            figure('Tag','Data video','units','inches','position',[12 4 6 6]); % Figure number 2.
            % 'position' vector is [left, bottom, width, height].
            % left, bottom control the position at which the window appears when it pops.
            
            if saveAvi == 1
                % clear mex % close all open .avi files (this avoids future errors).
                % Move into output folder to save videos:
                cd(new_folder_name);
                video_filename = strcat(image_label,'_Video','Traj_',num2str(n),'.avi'); % filename of video.
                % Create and open avi file for saving frames onto it later:
                mov = VideoWriter(video_filename);
                mov.FrameRate = framesPerSecond; % frame rate for the saved video.
                % Open video file for writing:
                open(mov);
            end
           
            % SUBARRAY: Initial size of subarray image around average centre-of-mass
            % position of particle. Corrected later if particle too close
            % to edge of image. An extra factor extraSubarraySizeFactor is
            % multiplied to allow for particle movement (see PARAMETERS):
            d_top = subarray_halfwidth*extraSubarraySizeFactor;
            d_bottom = subarray_halfwidth*extraSubarraySizeFactor;
            d_left = subarray_halfwidth*extraSubarraySizeFactor;
            d_right = subarray_halfwidth*extraSubarraySizeFactor;
            % If the particle is at edge of image, take a smaller subarray around
            % it but as large as possible until the edge is reached,
            % so reasign the d values:
            if (round(mean_yvalue)-d_top) < 1
                d_top = round(mean_yvalue) - 1;
            end
            if (round(mean_yvalue)+d_bottom) > frame_Ysize
                d_bottom = frame_Ysize - round(mean_yvalue);
            end
            if (round(mean_xvalue)-d_left) < 1
                d_left = round(mean_xvalue) - 1;
            end
            if (round(mean_xvalue)+d_right) > frame_Xsize
                d_right = frame_Xsize - round(mean_xvalue);
            end
            % Chose the minimum distance:
            d = min([d_top d_bottom d_left d_right]);
            
            for k = 1:length(frames_list) % loop through frames in track
                
                frame0 = image_data(frames_list(k)).frame_data; % extract frame data.
                frame0 = double(frame0);
                
                %% Create image subarray (frame) around particle average centre of mass:
                % Create squared image subarray of size
                % (2*d+1)x(2*d+1) centered on (mean_xvalue,mean_yvalue):
                frame = frame0(round(mean_yvalue)-d:round(mean_yvalue)+d,round(mean_xvalue)-d:round(mean_xvalue)+d);
                
                imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
                hold on;
                
                % Transform centre-of-mass to coordinate system referred to
                % subarray. x_values(k) is position on full original frame.
                % x_valuesB is position on subarray:
                x_valuesB = x_values(k)-(mean_xvalue-d)+1;
                y_valuesB = y_values(k)-(mean_yvalue-d)+1;
                
                plot(x_valuesB,y_valuesB,'x','Color','g','MarkerSize',10) % plot found centre-of-mass positions in green.
                
                % Plot major axis of ellipse to indicate particle orientation:
                xpointsMajorAxis = [
                    x_valuesB - 0.5*majAxisLength(k)*cos(angles_deg(k)/180*pi)
                    x_valuesB
                    x_valuesB + 0.5*majAxisLength(k)*cos(angles_deg(k)/180*pi)
                    ];
                ypointsMajorAxis = [
                    y_valuesB + 0.5*majAxisLength(k)*sin(angles_deg(k)/180*pi)
                    y_valuesB
                    y_valuesB - 0.5*majAxisLength(k)*sin(angles_deg(k)/180*pi)
                    ];
                plot(xpointsMajorAxis,ypointsMajorAxis,'y','LineWidth',1);
                
                pause(0.1); % this pause is needed to give time for the plot to appear (0.1 to 0.3 default)
                hold off;
                
                % ---
                if saveAvi == 1
                    % Save each frame plotted to the video mov:
                    F = getframe(gca);
                    writeVideo(mov,F)
                end
                % ---
            end
            
            if saveAvi == 1
                cd('..') % go back to previous folder after saving video file in new_folder_name
            end
        
        end
        
        % For a quick analysis of "good" trajectories, quickLook = 1, no user
        % input requested and trajectories flagged as "good tracking":
        if quickLook ==1
            good_tracking_flag = 1;
        else
            % CHECK: skip the following visual check of every track or not.
            good_tracking_flag = 1;
            % good_tracking_flag = input('Is the tracking "good" for this trajectory? (1 for "yes", anything else for "no"): '); % request user input.
            % flag saying if trajectory is a good one or not (bgnd point, not good tracking, etc).
        end
        
        close(findobj('Tag','Data video')); % close video figure;
        
        if good_tracking_flag == 1
            % only analyse n-th trajectory if tracking is good (and folder created for good-tracking trajectories only).
            % "good_tracking_flag" is added to result structure in
            % "showParticleTrajAnalysis.m", always 1, because if tracking is no good, the
            % trajectory is not analysed.
            
            % Analyse n-th trajectory data, produce result plots and save them:
            processedTraj = showParticleTrajAnalysis(trajXlsPath,image_data,analysedAllTraj,n,start_frame,tsamp,pixelsize_nm);
            
            % output: structure array, starting at element 1, only save to results trajs with good tracking.
            processedManyTrajs{n_good_tracking} = processedTraj; % Function output: cell array, starting at element 1.
            
            % Note that saving is done within function "showParticleTrajAnalysis.m".
            
            n_good_tracking = n_good_tracking + 1;
        end
        
    end
end



%% Save result (as .mat) in output folder which contains track results for this image sequence:

cd(new_folder_name) % move into folder corresponding to track results for this image sequence.
output_filename = strcat('procManyTraj',image_label);
save(output_filename,'processedManyTrajs'); % save variable processedManyTrajs to a .mat file
cd('..') % go back to previous folder.

% Save result (as .mat) in a folder specified by user input:
% % Do not save if we are just having a quick look at good trajectories
% % (quickLook = 1):
% if quickLook ~=1
%     uisave('processedManyTrajs','procManyTraj')
% end