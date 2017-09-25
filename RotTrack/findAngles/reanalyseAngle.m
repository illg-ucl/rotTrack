function reanalyseAngle(excelFileName,frameRateReal)
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
% - frameRateReal: Actual frame rate for the data in frames per second. In
% case analysis needs to be corrected. For Sonia+Andrea's data, frame rate
% is 15 frames per second.

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
% angleDegreesPos = NUMERIC(:,ID.angleDegreesPos); % originally postprocessed positive cyclic angle in degrees, can have
% values larger than 360 degrees. Obtained with function angleDegToPos.m, when this was causing some errors.
angleDegreesPos2 = angleDegToPos(angleDegrees); % use corrected postprocessing function.

disp('Excel file read successfully.');


%% Fit orientation angle versus time to a LINE to get angular velocity of
% particle rotation:
linearFunction = fittype('A + B*x','independent','x'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
% CHECK: values of guesses:
guess_A = 0; % offset guess in degrees.
guess_B = 360; % slope guess in degrees/s when correct frame rate has been entered.
options.StartPoint = [guess_A guess_B];
try % error control
    [fit_angle gof] = fit(timeabsReal,angleDegreesPos2,linearFunction,options);
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

% Plot angle versus time:
subplot(1,2,1);
plot(timeabsReal,angleDegrees,'.-b')
hold on;
plot(timeabsReal,angleDegreesPos2,'.-r')
% Plot angle versus time linear fit:
try % error control, if fit failed, do not plot
    plot(fit_angle,'k'); % plot linear fit to msd data as a black line.
catch ME3
end
xlabel('t_{abs} (s)'); 
ylabel('orientation angle (deg)');
xlim([0 max(timeabsReal)]);
% legend('angle(deg)','cyclicAngle(deg)');
legend('hide');
hold off;

% Display fit results:
subplot(1,2,2); 
axis off; 
% Define different lines of text:
str1(1) = {'Plots shown are for TRAJECTORY '};
str1(2) = {excelFileName};
str1(3) = {'Angular velocity, slope from fit:'};
str1(4) = {[num2str(fit_angle_slope),' +- ',num2str(fit_angle_slope_stDev),' degrees/s.']};
str1(5) = {'Offset from fit:'};
str1(6) = {[num2str(fit_angle_offset),' +- ',num2str(fit_angle_offset_stDev),' degrees/s.']};
str1(7) = {'Rsquared of fit:'};
str1(8) = {num2str(fit_angle_rsq)};
str1(9) = {'Frequency of rotation:'};
str1(10) = {[num2str(fit_angle_slope/360),' +- ',num2str(fit_angle_slope_stDev/360),' Hz.']};
str1(11) = {['Frame rate used for analysis (fps):',num2str(frameRateReal)]};
text(0,0.5,str1)


%% Save figure:
pos1 = strfind(excelFileName,'.xls'); % position of the start of the string '.xls' in the xls input file name.
base_name = excelFileName(1:(pos1-1)); % Filename without the extension.
figName = strcat(base_name,'_bis','.png'); % name of figure file to save to. 
% Export the current figure window at screen size as a png into current
% directory:
set(gcf, 'PaperPositionMode', 'auto')  % Use screen size. (gcf=get current figure)
% h = get(gcf);
print('-dpng','-r300',figName)  % add -r300 (to save at 300 dpi, higher resolution) after -dpng to control resolution.
close; % deletes the current figure (many open figures take too much memory and make Matlab crash).


%%  Export corrected analysis to another sheet of excel file:
% Track data corrected (vectors):
% Next element of results cell array contains track data:
n_trackData = 1; % cell element in which to save track info.
processedTraj{n_trackData}.frameNumber = frameNumber; % vectors.
processedTraj{n_trackData}.timeabsReal = timeabsReal;
processedTraj{n_trackData}.angleDegrees = angleDegrees; % orientation angle by fitting ellipse.
processedTraj{n_trackData}.angleDegreesPos2 = angleDegreesPos2; % positive orientation angle in degrees, cyclic, corrected.   

% Orientation results, angular velocity from fits (numbers):
n_angle = 2;
results_angle.fit_angle_offset = fit_angle_offset; 
results_angle.fit_angle_offset_stDev = fit_angle_offset_stDev;
results_angle.fit_angle_slope = fit_angle_slope;
results_angle.fit_angle_slope_stDev = fit_angle_slope_stDev;
results_angle.fit_angle_rsq = fit_angle_rsq;
results_angle.fit_rotationFreqHz = fit_angle_slope/360;
results_angle.fit_rotationFreqHz_stDev = fit_angle_slope_stDev/360;

processedTraj{n_angle} = results_angle;

% Save corrected results in two other sheets called "Track data corrected"
% and "Angular veloc corrected" within the input excel file:
output_filename = excelFileName; % name of excel file to save to.

dataForSheet1 = [fieldnames(processedTraj{n_trackData})'; num2cell(cell2mat(struct2cell(processedTraj{n_trackData})'))]; % vectors
dataForSheet2 = [fieldnames(processedTraj{n_angle}) struct2cell(processedTraj{n_angle})]; % numbers

warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.

xlswrite(output_filename,dataForSheet1,'Track data corrected'); % write data to sheet 'Track data' in excel file.
xlswrite(output_filename,dataForSheet2,'Angular veloc corrected'); % write data to sheet 'Angular veloc results' in excel file.

