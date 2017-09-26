function reanalyseAngle(excelFileName,frameRateReal,thresh_slope)
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
% the angle and fitting the appropriate linear regions. Use input frame
% rate. All valid for clockwise rotation.
%
% INPUTS:
% - excelFileName: string with the name of the Excel file containing the original extracted
% angle and the postprocessed angle (result of showManyParticleTrajAnalysis.m).
% - frameRateReal: Actual frame rate for the data in frames per second. In
% case analysis needs to be corrected. For Sonia+Andrea's data, frame rate
% is 15 frames per second.
% - thresh_slope: threshold slope in degrees/s. Only linear data with a
% positive slope larger than this threshold is fitted to obtain angular
% velocity and rotation frequency. 360deg/s corresponds to 1Hz. 
% E.g., 250-300deg/s is a good threshold for 10Hz field. For 1Hz field, ~130deg/s is better.
%
% OUTPUTS:
% A new png graph file is generated with name "bis" appended to original
% name, containing the new angle vs time plot, fits and fit results for individual linear sections. 
% A new Excel file "bis" with new sheets is created with the
% revised analysis for the postprocessed angle vs time and for the data and
% fit results for the relevant linear sections.

% PARAMETERS
% - minSectionPoints: minimum number of points in a linear angle vs time
% section for it to be fitted to a line to obtain the slope (angular
% velocity):
minSectionPoints = 6;

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


%% Automatically find segments to fit to a line
% Find linear segments with second derivative equal to zero and fit those
% linear segments that have at least 6-10 points in them.

% Calculate numerical first and second derivatives of angle data (using points two frames away):
deltaTime = timeabsReal(2)-timeabsReal(1); 
diff1 = zeros(length(angleDegreesPos2)-2,1); % initialise column vector
for i = 1:length(diff1)
    % first derivative:
    diff1(i) = (angleDegreesPos2(i+2) - angleDegreesPos2(i))/(2*deltaTime);
end
% Smooth out noise by calculating a simple ('s') moving average:
diff1_smooth = tsmovavg(diff1, 's', 3, 1); % Note that first two points are lost.

% % Aid plots:
% plot(timeabsReal,angleDegreesPos2,'.-r')
% hold on;
% plot(timeabsReal(1:length(timeabsReal)-2),diff1,'.-k')
% plot(timeabsReal(1:length(timeabsReal)-2),diff1_smooth,'.-g')
% min(diff1_smooth)
% close;

% Positions with a slope larger than the input threshold slope:
pos_largeSlope = find(diff1_smooth > thresh_slope);
% Find differences:
pos_diffs = zeros(length(pos_largeSlope)-1,1);
for i = 1:length(pos_diffs)
    pos_diffs(i) = pos_largeSlope(i+1) - pos_largeSlope(i);
end
% Positions of jumps:
pos_jumps = find(pos_diffs~=1);
% Separate into separate sections of large enough slope:
k_start = 1; % initialise
for j = 1:length(pos_jumps)+1
    if j == length(pos_jumps)+1 % do differently for last section
        k_end = length(pos_largeSlope);
    else
        k_end = pos_jumps(j);
    end
    section{j} = pos_largeSlope(k_start:k_end); % section with points with large slope
    k_start = k_end + 1;
end

% We only analyse sections with a number of points larger than
% minSectionPoints:
sectionsToAnalyse = [];
for i = 1:size(section,2)
    if size(section{i},1) > minSectionPoints
        sectionsToAnalyse = [sectionsToAnalyse; i];
    end
end

%% Fit orientation angle versus time to a LINE to get angular velocity of
% particle rotation:

for i = 1:length(sectionsToAnalyse)
    vectorPos = section{sectionsToAnalyse(i)}; % positions to analyse in angle vector
    angleToAnalyse{i} = angleDegreesPos2(vectorPos);
    timeToAnalyse{i} = timeabsReal(vectorPos);    
    % Fit to a line:
    results_angle{i} = fitToLine(timeToAnalyse{i},angleToAnalyse{i});  
end

% Plot angle versus time:
figure;
subplot(1,2,1);
plot(timeabsReal,angleDegrees,'.-b')
hold on;
plot(timeabsReal,angleDegreesPos2,'.-r')
% Plot angle versus time linear fits on same plot:
for i = 1:length(sectionsToAnalyse)
    try % error control, if fit failed, do not plot
        xdata = timeToAnalyse{i};
        offset = results_angle{i}.fit_result_offset;
        slope = results_angle{i}.fit_result_slope;
        ydata = offset + slope*xdata;
        plot(xdata,ydata,'k'); % plot linear fit as black line only over fitted data range.
        % plot(results_angle{i}.fit_result,'b'); % plot linear fit over entire data range as a black line.
    catch ME3
    end
end
xlabel('t_{abs} (s)'); 
ylabel('orientation angle (deg)');
xlim([0 max(timeabsReal)]);
% legend('angle(deg)','cyclicAngle(deg)');
legend('hide');
hold off;


% Calculate final frequency of rotation for particle in Hz, as mean and
% stdev of all slopes fitted or as the only slope fitted with its stdev error:
freqRotHz = [];
for i = 1:length(sectionsToAnalyse)
    freqRotHz = [freqRotHz, results_angle{i}.fit_result_slope/360];
end
if length(freqRotHz)>1
    freqRotHz_final = mean(freqRotHz);
    freqRotHz_stDev = std(freqRotHz);
else
    freqRotHz_final = freqRotHz(1);
    freqRotHz_stDev = results_angle{1}.fit_result_slope_stDev/360;
end


% Display fit results:
subplot(1,2,2); 
axis off; 
% Define different lines of text:
str1(1) = {excelFileName};
str1(2) = {['Frame rate (fps) used: ',num2str(frameRateReal)]};
str1(3) = {['Mean final Freq(Hz): ',num2str(freqRotHz_final),' +- ',num2str(freqRotHz_stDev)]};

k = 6; % total number of writen rows below
for i = 1:length(sectionsToAnalyse)
    str1(4+k*(i-1)) = {'Angular velocity, slope from fit:'};
    str1(5+k*(i-1)) = {[num2str(results_angle{i}.fit_result_slope),' +- ',num2str(results_angle{i}.fit_result_slope_stDev),' degrees/s.']};
    str1(6+k*(i-1)) = {[num2str(results_angle{i}.fit_result_slope/360),' +- ',num2str(results_angle{i}.fit_result_slope_stDev/360),' Hz.']};   
    str1(7+k*(i-1)) = {'Offset from fit:'};
    str1(8+k*(i-1)) = {[num2str(results_angle{i}.fit_result_offset),' +- ',num2str(results_angle{i}.fit_result_offset_stDev),' degrees/s.']};
    str1(9+k*(i-1)) = {['R-squared of fit: ',num2str(results_angle{i}.fit_result_rsq)]};
    end
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

% Save corrected results in several new sheets called "Track data corrected"
% and "Angular veloc corrected i" within the input excel file:
warning off MATLAB:xlswrite:AddSheet % turn warning off when new sheet added to excel file.
output_filename = [base_name,'bis','.xls']; % name of excel file to save to.

dataForSheet1 = [fieldnames(processedTraj{n_trackData})'; num2cell(cell2mat(struct2cell(processedTraj{n_trackData})'))]; % vectors
xlswrite(output_filename,dataForSheet1,'Track data corrected'); % write data to sheet 'Track data corrected' in excel file.

for i = 1:length(sectionsToAnalyse)   
    sheetName1 = ['Data section ',num2str(i)];
    columnTitles = {'timeToAnalyse(s)', 'angleToAnalyse(deg)'};
    vectorsToSave =  num2cell(cell2mat({timeToAnalyse{i} angleToAnalyse{i}}));
    dataForSheetA = [columnTitles; vectorsToSave];
    xlswrite(output_filename,dataForSheetA,sheetName1); % write data to sheet in excel file.
    
    sheetName2 = ['AngVelocFit Section ',num2str(i)];
    dataForSheetB = [fieldnames(results_angle{i}) struct2cell(results_angle{i})]; % numbers
    xlswrite(output_filename,dataForSheetB,sheetName2); % write data to sheet in excel file.    
end

% Save final frequency (mean of varios slopes obtained) in another sheet:
dataForSheetFinalFreq = [{'finalFreqHz_mean'; 'finalFreqHz_stDev'}, num2cell(cell2mat({freqRotHz_final; freqRotHz_stDev}))];
xlswrite(output_filename,dataForSheetFinalFreq,'finalFreqHz');
    
% Save re-analysis parameters:
dataReanalysisParams = [{'frameRateReal(fps)'; 'thresh_slope(deg/s)'; 'minSectionPoints(frames)'}, num2cell(cell2mat({frameRateReal; thresh_slope; minSectionPoints}))];
xlswrite(output_filename,dataReanalysisParams,'Params Reanalysis');
        



