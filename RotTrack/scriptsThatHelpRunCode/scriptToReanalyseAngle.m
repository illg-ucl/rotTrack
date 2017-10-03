function scriptToReanalyseAngle
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
% Reanalyse a set of excel files containing track data for frame number, time, angle,
% etc, to obtain angular velocity and frequency of rotation. The raw angle
% is postprocessed assuming anti-clockwise rotation. Then the linear
% sections with a large enough slope (above input thresh_slope) in the 
% the angle vs time plot are automatically detected and fitted to a line. 
% The input frame rate is used to convert frame number into time (s).
%
% IMPORTANT: calculations assume anti-clockwise rotation.
%
% Before running this function, make sure Matlab's current directory is the
% directory containing the analysis folder that in turn contain the analysis excel files.

%% PARAMETERS:
frameRate = 30; % frame rate in frames per second.
minSectionPoints = 5; % minimum number of points in a linear section (in angle vs time plot)
% for it to be fitted to a line to obtain the slope (angular velocity).
thresh_slope_1Hz = 130; % minimum slope in a linear section for it to be fitted to a line. 
thresh_slope_5Hz = 500;
% Value in degrees/s. 360deg/s corresponds to a frequency of 1Hz.
% E.g., 250-300 deg/s is a good threshold for 10Hz rotating field. 
% For 5Hz field, ~200 deg/s is good.
% For 1Hz field, ~130 deg/s is good.

% Name used for analysis folders (in scriptToAnalyseManyVideos.m):
data_set_label = 'analysis';

%% Make list of folder names:
% There is one folder per video file.
% Choose video file extension:
videoFile_extension = '.m4v';
listVideoNames0 = dir(strcat('*',videoFile_extension));
% Error control:
if isempty(listVideoNames0) % If there are no video files of chosen extension, show error and exit function:
    error('Check you are in the correct directory. No video files of chosen extension found in folder.');
end
% List of video file names:
listVideoNames = {listVideoNames0.name}; % cell array of strings with video file names.
% Generate list of folder names:
for k = 1:length(listVideoNames)
    fullName = listVideoNames{k};
    pos = strfind(fullName,videoFile_extension); % position of the start of the string videoFile_extension in the file name.
    folderName{k} = strcat(data_set_label,'_',fullName(1:(pos-1)),'_'); 
end


%% Loop through analysis folders:

for k = 1:length(folderName)
    
    folderName{k} % print to command window
    cd(folderName{k}) % move into folder
    
    % Choose appropriate slope threshold for reanalysis:
    if ~isempty(strfind(folderName{k},'1Hz'))
        thresh_slope = thresh_slope_1Hz;
    elseif ~isempty(strfind(folderName{k},'5Hz'))
        thresh_slope = thresh_slope_5Hz;
    end
    
    % Find excel files in current directory:
    list_excelFiles = dir('*.xls');
    
    for i = 1:length(list_excelFiles)
        reanalyseAngle(list_excelFiles(i).name,frameRate,thresh_slope,minSectionPoints);
    end
    
end
    


