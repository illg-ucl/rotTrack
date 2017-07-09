function plotParticleTrajNumbers(image_label,minPointsTraj) 
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
% This function plots Trajectory Numbers (as they appear on the Excel file
% generated by function linkTrajSegmentsParticles) next to each corresponding
% particle overlaid on the first frame of the chosen image sequence. This
% allows easy further analysis, excluding beads that are too close, etc.
% Only trajectories with more than a certain number, minPointsTraj, of
% points in them are looked at.
%
% INPUTS: 
% - 'image_label' string that labels the image sequence under analysis, e.g. '101'.
% - minPointsTraj: minimum number of data points that a trajectory must have in order to be
% analised.
%
% OUTPUTS:
% This function plots and saves a .png image with the overlaid trajectory
% (track) numbers in the current directory.
% 
% NOTE: to run this function, you need to be in the directory that contains
% the image sequence and the excel file (both with same label, e.g., '25').

%% Read in the image-sequence data:

% Read image-sequence file:
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);


%% Get path for trajectory data (excel file):

% You need to be in the correct directory before running the function!!!!
% Find paths in current folder which contain the 'image_label' string and are .xls files:
trajXlsPath0 = dir(strcat('*',image_label,'*.xls')); % Trajectory data path (excel file with the full trajectories as returned by function "linkTrajSegmentsParticles.m").
% Error control:
if isempty(trajXlsPath0) % If there is no .xls trajectory data file for such image number, show error and exit function:
    error('Check you are in the correct directory and run again. No .xls file found for that image label. Make sure image number is in between quotes ''.'); 
end
trajXlsPath = trajXlsPath0.name;
% trajXlsPath0 is a structure and the file names is stored in the field 'name'.


%% Open all trajectory data (excel file): 
% (.xls file previously generated with functions 'FindTrajectsParticles' and 'linkTrajSegmentsParticles'):

% error control:
if 2~=exist('xlsread') 
    error('Check file dependencies - you need to install xlsread'); 
end

% Open excel file and read the data:
[NUMERIC,TXT,RAW]=xlsread(trajXlsPath,'Track results'); % import the data in the sheet named 'Track results'.
% Import the column heads and assign ID
colheads = TXT;
% The column titles are: Xcom, Ycom, ClipFlag,
% FrameNumber, ParticleNumber,	TrajNumber, etc.

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
tracks(1:numtracks) = struct('trajNumber',[],'Xcom',[],'Ycom',[],'mean_Xcom',[],'mean_Ycom',[]);
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
    % all values in pixels.
    tracks(i).Xcom = data{i}(1:end,ID.Xcom); % original Xcom (centre of mass) in image.
    tracks(i).Ycom = data{i}(1:end,ID.Ycom); % original Ycom in image.        
    tracks(i).mean_Xcom = mean(data{i}(1:end,ID.Xcom)); % mean Xcom value at which the Traj Number will be displayed. 
    tracks(i).mean_Ycom = mean(data{i}(1:end,ID.Ycom)); % mean Ycom value at which the Traj Number will be displayed.
    else
        % save indices to delete later:
        del(i) = i;     
    end
    
end

% Delete tracks which were too short: 
tracks(find(del))=[];


%% Plot first frame and overlay trajectory numbers:

traj_number = [tracks(:).trajNumber]'; % column vector with track numbers as in excel file.
xpos_vector = [tracks(:).mean_Xcom]'; % column vector with mean x positions to display track number.
ypos_vector = [tracks(:).mean_Ycom]'; % column vector with mean y positions to display track number.

% Plot numbers of accepted particle tracks (index j) overlaid 
% on top of first frame to be able to visually identify them:
frame1 = image_data(1).frame_data; % extract frame data which is stored in field 'frame_data'.
figure;
imshow(frame1,[],'Border','tight'); % show image scaled between its min and max values ([]).
hold on;

for j=1:length(traj_number)
    text(xpos_vector(j),ypos_vector(j),num2str(traj_number(j)),'Color',[1 1 0],'FontSize',8); % number in yellow.
end
pause(0.1); % this pause is needed to give time for the plot to appear
hold off;
% SAVE current FIGURE as a .png and then close the figure window:
figName = strcat('Traj_Numbers_',image_label); % choose name of figure file for saving .png.
saveFigurePNG(cd,figName); % Save file in current directory. See saveFigurePNG.m.