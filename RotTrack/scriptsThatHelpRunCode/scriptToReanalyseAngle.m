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
% etc, to obtain angular velocity and frequency of rotation, postprocessing
% the angle and fitting the appropriate linear regions. All valid for clockwise rotation.
% Before running this function, make sure Matlab's current directory is the
% directory containing the excel files.

% PARAMETERS:
frameRate = 15; % frame rate in frames per second.
thresh_slope = 130; % minimum slope in a linear section for it to be fitted to a line. 
% Value in degrees/s. 360deg/s corresponds to a frequency of 1Hz.
% E.g., 250-300 deg/s is a good threshold for 10Hz rotating field. 
% For 5Hz field, ~200 deg/s is good.
% For 1Hz field, ~130 deg/s is good.
minSectionPoints = 5; % minimum number of points in a linear section (in angle vs time plot)
% for it to be fitted to a line to obtain the slope (angular velocity).

% Find excel files in current directory:
list_excelFiles = dir('*.xls');

for i = 1:length(list_excelFiles)
   reanalyseAngle(list_excelFiles(i).name,frameRate,thresh_slope,minSectionPoints); 
end

