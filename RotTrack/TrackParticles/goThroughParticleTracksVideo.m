function good_tracks = goThroughParticleTracksVideo(image_label,data_set_label,n_traj_start,n_traj_end,minPointsTraj,maxMajorAxisLength,showVideos) 
%
% good_tracks = goThroughParticleTracksVideo(image_label,data_set_label,n_traj_start,n_traj_end,minPointsTraj,maxMajorAxisLength)  
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
% Plot video of image with tracks overlaid on top for a given image sequence 'image_label', from trajectory
% 'n_traj_start' to 'n_traj_end'.
% The user is requested for input after each track is shown to say if it is
% a "good" one (1) or not (0).
% 
% NOTE: before running this function you should move into a directory which
% contains both the video file (labelled by 'image_label') which has previously been analysed with
% 'FindTrajectsParticles.m' and 'linkTrajSegmentsParticles.m' to produce the .xls file which
% contains the trajectory results, which should be in the same directory as the video file.
%
% INPUTS: 
% - image_label: string that labels a given image sequence found in current
% folder. The code finds the path of the image file automatically based on a string label 
% that is equal to the file name without the file extension. For example,
% for image video file "210217r25.tif", an appropriate label would be the
% string '210217r25'. This will be used throughout the entire RotTrack code.
% - data_set_label: string that labels set of data or parameters. Use
% same as in input to linkTrajSegmentsParticles.m.
% - n_traj_start: first trajectory we want to analyse and check.
% - n_traj_end: last trajectory we want to analyse and check. If the
% string 'end' is entered, we go through to the last analysed trajectory.
% - minPointsTraj: minimum number of data points that a trajectory must have in order to be
% analised.
% - maxMajorAxisLength: maximum length in pixels that a particle can have
% (major axis of fitted ellipsoid) to be accepted as valid.
% (A value of at least 6 needs to be used for all methods in
% "showTrajAnalysis.m" (and therefore "showManyTrajAnalysis.m") to work
% well).
% - showVideos: 1 if user wants to see one video for each track and provide
% manual input as to whether each trajectory is good or not. This can take
% a long time if there are many tracks. Set to 0 to not show any videos and 
% accept all tracks with more than "minPointsTraj" points and a max length 
% below "maxMajorAxisLength" as good tracks without requiring user input for each track.

% -----------------
% IMPORTANT NOTE!!!: The values of all inputs: "image_label",
% "n_traj_start", "n_traj_end" and "minPointsTraj" here need to be the same
% as those used later on for function "showManyParticleTrajAnalysis.m", 
% otherwise the trajectory numbers will be different!.
% -----------------
%
% OUTPUT: 
% - good_tracks is a structure with fields:
% good_tracks.image_label = image_label; input value.
% good_tracks.n_traj_start = n_traj_start; input value.
% good_tracks.n_traj_end = n_traj_end; input value.
% good_tracks.minPointsTraj = minPointsTraj, input value.
% good_tracks.track_numbers: is a row vector containing the numbers of tracks considered as
% "good" by the user after seing the videos of the track overlaid on the
% image sequence.
% The output is saved as a .mat file.
%
% Example of how to call this function:
% gt = goThroughParticleTracksVideo('1757_1',1,'end',3);
% gt = goThroughParticleTracksVideo('498',1,'end',6);
% ------------------------------------


%% Get path for trajectory data (excel file):

% You need to be in the correct directory before running the function!!!!
% Find paths in current folder which contain 'image_label' string:
trajXlsPath0 = dir(strcat('*',data_set_label,'_',image_label,'_fullTrajs.xls')); % Trajectory data path (excel file with the full trajectories as returned by function "linkTrajSegmentsParticles.m").

% Error control:
if isempty(trajXlsPath0) % If there is no .xls trajectory data file for such image number, show error and exit function:
    error('Check you are in the correct directory and run again. No .xls file found for that image number. Make sure image number is in between quotes ''.'); 
end
trajXlsPath = trajXlsPath0.name;
% trajXlsPath0 is a structure and the file names is stored in the field 'name'.


%% Open and analyse all trajectory data (excel file): 
% (.xls file previously generated with functions 'FindTrajects' and 'linkTrajSegmentsParticles'):
% ======================
% The following derives from old function analyseTraj(file,tsamp,minPointsTraj):

% Error control:
if 2~=exist('xlsread') 
    error('Check file dependencies - you need to install xlsread'); 
end

% Open excel file and read the data:
[NUMERIC,TXT,RAW]=xlsread(trajXlsPath,'Track results'); % import the data in the sheet named 'Track results'.
% Import the column heads and assign ID
colheads = TXT;
% The column titles are: Xcom, Ycom, AngleDegrees, ClipFlag, FrameNumber, ParticleNumber,
% TrajNumber, etc.

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
tracks(1:numtracks) = struct('TrajNumber',[],'Xcom',[],'Ycom',[],'FrameNumber',[], ...
    'AngleDegrees',[],'majorAxisLength',[],'numel',[],'minNumPointsInTraj',[]);
del = []; % initialise for later.

for i=1:numtracks 
    % i
    a = I(i); % index for starting point in trajectory.
    b = Y(i); % index for ending point in trajectory.
    
    % Delete tracks that are less than minPointsTraj data points long:
    if b-a+1 >= minPointsTraj  % Only analyse tracks which have at least "minPointsTraj" points (frames) in them (5, or 15, e.g.).
        
        data{i} = NUMERIC(a:b,:);
        % tracks(i).XLS.track_index = A(i);
        tracks(i).TrajNumber = A(i);
        tracks(i).Xcom = data{i}(1:end,ID.Xcom); % centre of mass x values on image (used later for plotting).
        tracks(i).Ycom = data{i}(1:end,ID.Ycom); % centre of mass y values on image (used later for plotting).
        tracks(i).FrameNumber = data{i}(1:end,ID.FrameNumber); % frame number.
        tracks(i).AngleDegrees = data{i}(1:end,ID.AngleDegrees); % orientation angle in degrees.
        tracks(i).majorAxisLength = data{i}(1:end,ID.majorAxisLength); % length of major axis of ellipsoid fitted to particle shape.
        tracks(i).numel = b-a+1; % Number of points in the track.
        tracks(i).minNumPointsInTraj = minPointsTraj;
        
    else
        % save indices to delete later:
        del(i) = i;
    end
    
end

% Delete tracks which were too short: 
tracks(find(del))=[];

analysedAllTraj = tracks;

% ======================

if isempty(analysedAllTraj)
    disp('The total number of long enough trajectories (to analyse) for this file is zero.');
    disp('Exiting program');
    return % exits function.
end

%'The total number of trajectories analysed (long enough) in the file is: 
n_trajs_analysed = length(analysedAllTraj);
disp(['The number of long enough tracks to go through for this image sequence is: ' num2str(n_trajs_analysed)])

if strcmp(n_traj_end,'end') % string compare, true if equal
    % go through all analysed trajectories til the last one: 
   n_traj_end = n_trajs_analysed;
end


%% Read in the image-sequence data:

% Read image-sequence file:
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.



%% Loop through trajectories:

good_track_numbers = []; % initialise empty vector where I will store numbers of tracks labelled as "good" by the user.

for n = n_traj_start:n_traj_end
     
    % Close any pre-existing figures:
    close(findobj('Tag','Trajectory results'));
    % close(all); % close all figures.
   
    frames_list = analysedAllTraj(n).FrameNumber; % list of frame numbers in trajectory n.
    x_values = analysedAllTraj(n).Xcom; % list of x centres of particle in trajectory n.
    y_values = analysedAllTraj(n).Ycom; % list of y centres of particle in trajectory n. 
    
    angles_deg = analysedAllTraj(n).AngleDegrees; % list of orientation angles in degrees.
    majAxisLength = analysedAllTraj(n).majorAxisLength; % list of lengths of major axis in pixels.
    
    if mean(majAxisLength) < maxMajorAxisLength
        
        if showVideos == 1
            
            % Show video of trajectory overlaid on actual image:
            % Loop through frames in each trajectory analysed:
            figure('Tag','Data video','Units','normalized','Position',[0 1 0.2 0.2]); % Figure number 2.
            % 'position' vector is [left, bottom, width, height].
            % left, bottom control the position at which the window appears when it
            % pops.
            
            for k = 1:length(frames_list) % loop through frames in track
                
                frame = image_data(frames_list(k)).frame_data; % extract frame data which is stored in field 'frame_data'.
                frame = double(frame);
                
                imshow(frame,[],'Border','tight'); % show image scaled between its min and max values ([]).
                hold on;
                
                plot(x_values(k),y_values(k),'x','Color','g','MarkerSize',5) % plot accepted particle centres in green.
                
                % Plot major axis of ellipse to indicate particle orientation:
                xpointsMajorAxis = [
                    x_values(k) - 0.5*majAxisLength(k)*cos(angles_deg(k)/180*pi)
                    x_values(k)
                    x_values(k) + 0.5*majAxisLength(k)*cos(angles_deg(k)/180*pi)
                    ];
                ypointsMajorAxis = [
                    y_values(k) + 0.5*majAxisLength(k)*sin(angles_deg(k)/180*pi)
                    y_values(k)
                    y_values(k) - 0.5*majAxisLength(k)*sin(angles_deg(k)/180*pi)
                    ];
                plot(xpointsMajorAxis,ypointsMajorAxis,'y','LineWidth',1);
                
                pause(0.1); % this pause is needed to give time for the plot to appear (0.1 to 0.3 default)
                hold off;
                
            end
            
            close(findobj('Tag','Data video')); % close video figure;
            
            disp(['Track number ',num2str(n),' out of ' num2str(n_trajs_analysed) ' tracks:'])
            
            % REQUEST USER INPUT: for GOOD tracking or not:
            good_tracking_flag = input('Is the tracking "good" for this trajectory? (1 for "yes", anything else for "no"): ');
            % flag saying if trajectory is a good one or not (bgnd point, not good tracking, etc.).
            
            if good_tracking_flag == 1
                good_track_numbers = [good_track_numbers n]; % append to good_tracks (track numbers) vector.
            end           
            
        else % when showVideos is different from 1, all tracks with more than
            % "minPointsTraj" points and a max length below "maxMajorAxisLength"
            % are stored as good ones           
            good_track_numbers = [good_track_numbers n]; % append to good_tracks (track numbers) vector.            
        end
    end    
    
end

disp(['The number of tracks that fullfill all input conditions are: ' num2str(length(good_track_numbers))])


% OUTPUT: structure with two fields;
good_tracks.image_label = image_label;
good_tracks.n_traj_start = n_traj_start;
good_tracks.n_traj_end = n_traj_end;
good_tracks.minPointsTraj = minPointsTraj;
good_tracks.maxMajorAxisLength = maxMajorAxisLength;
good_tracks.good_track_numbers = good_track_numbers;

% Save result (as .mat):
output_filename = strcat('good_tracks_',image_label,'.mat');
save(output_filename,'good_tracks') % save variable good_tracks.

% --------------------
