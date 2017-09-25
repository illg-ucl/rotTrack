function [time, angleRaw] = reanalyseAngle(excelFileName)
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
% INPUTS:
% - excelFileName: string with the name of the Excel file containing the original extracted
% angle and the postprocessed angle (result of showManyParticleTrajAnalysis.m).

% PARAMETERS:
% Actual frame rate for the data in frames per second:
frameRateReal = 15; 

% Open Excel file and read the data:

[NUMERIC0,TXT0,RAW0]=xlsread(excelFileName,'Track info'); % import the data in the sheet named 'Track info'.
pos_subarraySize = find(strcmp(TXT0,'TimeBetweenFrames')); % row number corresponding to 'TimeBetweenFrames' value.
timeBetweenFrames = RAW0{pos_subarraySize,2}; % extract TimeBetweenFrames (in seconds).
frameRate = 1/timeBetweenFrames;

[NUMERIC,TXT,RAW]=xlsread(excelFileName,'Track data'); % import the data in the sheet named 'Track data'.
% Import the column heads and assign ID
colheads = TXT;

% Generate ID: ID is a structure with fiels with the same names as the
% column titles.
for i=1:numel(colheads) 
    ID.(colheads{i}) = find(strcmp(TXT,colheads{i})); 
end

% Obtain the relevant vectors:
frameNumber = NUMERIC(:,ID.frame); % NUMERIC is the numeric data read from the excel file (without the row of column titles).
% timeabs = NUMERIC(:,ID.timeabs); % absolute time in seconds using previous frameRate value.
timeabsReal = frameNumber/frameRateReal; % corrected time in seconds for the frameRateReal in PARAMETERS.
angleDegrees = NUMERIC(:,ID.angleDegrees); % raw angle data between -90 and 90 degrees.
angleDegreesPos = NUMERIC(:,ID.angleDegreesPos); % postprocessed positive cyclic angle in degrees, can have
% values larger than 360 degrees. Obtained with function angleDegToPos.m.

disp('Excel file read successfully.');

% Plots:
plot(timeabsReal,angleDegrees,'.-b')
hold on;
plot(timeabsReal,angleDegreesPos,'.-r')
% % Plot angle versus time linear fit:
% try % error control, if fit failed, do not plot
%     plot(fit_angle,'k'); % plot linear fit to msd data as a black line.
% catch ME3
% end
xlabel('t_{abs} (s)'); 
ylabel('orientation angle (deg)');
xlim([0 max(timeabsReal)]);
% legend('angle(deg)','cyclicAngle(deg)');
legend('hide');
hold off;

% OUTPUTS:
time = timeabsReal;
angleRaw = angleDegrees;
