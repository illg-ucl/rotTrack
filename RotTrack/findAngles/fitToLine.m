function resultsFit = fitToLine(xvector,yvector)
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
% Fit two vectors of data (y versus x) to a line and give fit results as
% output.
% 
% INPUTS:
% - xvector: vector of x values. This is time in seconds for fit of angle
% vs time.
% - yvector: vector of y values to fit as a function of x. This is angle in degrees for fit of angle
% vs time.
%
% OUTPUTS:
% - resultsFit is the output structure with fields: 
% fit_result - full result of fit, can be used to plot directly.
% fit_result_offset - offset from fit, in degrees.
% fit_result_offset_stDev  - standard deviation error of offset from fit, in degrees.
% fit_result_slope - slope from fit, in degrees/s.
% fit_result_slope_stDev - standard deviation error of slope from fit, in degrees/s.
% fit_result_rsq - Rsquared of fit.
% fit_rotationFreqHz - slope in Hz, i.e., rotation frequency in Hz.
% fit_rotationFreqHz_stDev - standard deviation error of slope (rotation frequency) in Hz.



%% PARAMETERS:
% Values of guesses for fit, set appropriate values:
guess_A = 0; % offset guess. In degrees for fit of angle vs time.
guess_B = 360; % slope guess. In degrees/s for fit of angle vs time. 360 is a good guess for a rotation frequency of ~1Hz.


%% FIT TO A LINE:

linearFunction = fittype('A + B*x','independent','x'); % define linear funtion to fit to, with 'x' as independent variable;
options = fitoptions('Method','NonlinearLeastSquares');  % Creates a structure of fit options with fields StartPoint, Lower, Upper, etc.
options.StartPoint = [guess_A guess_B]; % use guesses defined in PARAMETERS section above.

try % error control
    [fit_result gof] = fit(xvector,yvector,linearFunction,options);
    fit_param_values = coeffvalues(fit_result); % parameter values resulting from fit. First one is 'A', second one is 'B'.
    % Use coeffnames(fit_result) to find out order of parameters.
    fit_result_offset = fit_param_values(1); % offset from fit.
    fit_result_slope = fit_param_values(2); % slope from fit.
    fit_result_rsq = gof.rsquare; % rsquare coefficient of fit.
    errors = confint(fit_result,0.682); % 68.2% confidence interval for each fit parameter (lower and upper bounds as first and second rows).
    errorSTDEV = (errors(2,:)-errors(1,:))/2; % Standard deviation of each fit parameter (probability to be between -STDEV and +STDEV is 68.2%).
    fit_result_offset_stDev = errorSTDEV(1);
    fit_result_slope_stDev = errorSTDEV(2);
catch ME1 % catch error to prevent exiting function. If there was an error:
    fit_result = [0]; % needed for plot later.
    fit_result_offset = []; 
    fit_result_offset_stDev = []; 
    fit_result_slope = []; 
    fit_result_slope_stDev = []; 
    fit_result_rsq = []; 
end

% Offset in degrees, slope in degrees/s, 
resultsFit.fit_result = fit_result;
resultsFit.fit_result_offset = fit_result_offset; 
resultsFit.fit_result_offset_stDev = fit_result_offset_stDev;
resultsFit.fit_result_slope = fit_result_slope;
resultsFit.fit_result_slope_stDev = fit_result_slope_stDev;
resultsFit.fit_result_rsq = fit_result_rsq;
resultsFit.fit_rotationFreqHz = fit_result_slope/360;
resultsFit.fit_rotationFreqHz_stDev = fit_result_slope_stDev/360;

disp(' ') % empty line
disp([' fit_result_offset = ',num2str(fit_result_offset),' +- ',num2str(fit_result_offset_stDev),';   fit_result_slope = ',num2str(fit_result_slope),' +- ',num2str(fit_result_slope_stDev)]) 
disp([' fit_result_rsq = ',num2str(fit_result_rsq)]) % r squared of fit.
