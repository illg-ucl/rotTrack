function scriptToRun
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

%% To analyse an entire video sequence: 
%
% STEP BY STEP GUIDE for RotTrack, rotational tracking of particles in a video sequence:
%
% E.g., to analyse video "5mT-1Hz.m4v":
% 
% - 1. Make sure the folder containing the RotTrack source code is
% added to the Matlab search path. Go to Home tab -> Set Path -> Add with
% subfolders (select appropriate folder) -> save.
%
% - 2. The source code has been tested using MATLAB R2016a (9.0.0.341360),
% 64-bit (Win64), February 11, 2016. It requires some MATLAB packages to
% run. These are usually included in default MATLAB installations but check
% by typing in the command window:
% >> license('inuse')
% The following packages should appear:
% curve_fitting_toolbox
% image_toolbox
% matlab
% signal_toolbox
% statistics_toolbox
% If they do not appear, they should be installed.
%
% - 3. Make sure the Current Folder within Matlab is the folder that contains
% your image video files. 
data_folder = cd;

% - 4. Define image_label: 
image_label = '5mT-1Hz';
% The code works finding paths and files automatically based on a string label 
% that is equal to the file name (without the file extension). For example, for image video file "210217r25.tif", 
% an appropriate label would be the string '210217r25'. This will be used
% throughout the analysis.
% 
% TRICK: If you don't know the number of frames in a given image you can do:
% frame1 = extract1frameB(1);
% and select the image from a folder. The command window will print the
% number of frames on the image sequence, the file path and frame size.
%
% - 5. PARAMETERS: There are a number of important parameters that need to be set right.
% These are within functions FindTrajectsParticles.m and
% linkTrajSegmentsParticles.m. Find these functions and tweak these parameters
% on the PARAMETERS section at the beginning of each .m file. 
% The key parameters in FindTrajectsParticles.m are:
% - subarray_halfwidth (Default: 25 pixels). Halfwidth of image square subarray which includes particle and background around it.
% - inner_circle_radius (Default: 20 pixels); Radius of circular mask that contains the entire particle. 
% For eliminating coincident positions:
% - d_coincid_cand = 10; % distance (in pixels) for eliminating coincidences in particle-position candidates. Default: 3 pixels.
% - d_coincid_found = 3; % distance for eliminating coincidences in found particle centres.
% - d_01_max (Default: 30 pixels); Max distance in pixels between particle
% centres in one frame and the next one, so that they can be linked into a
% trajectory. This depends on the particle size on the image, the frame rate
% and how much particles move between frames. Tweak to match experimental
% conditions.
% - d_02_max (Default: 30 pix). Similar to above but for linking one frame
% and two frames later.
% The key parameters in linkTrajSegmentsParticles are:
% - d_01_max (Default: 30 pixels); Max distance in pixels between particle
% centres for linking different trajectory segments, i.e., for linking end
% particle position in one trajectory with start particle position in another trajectory.
% - Frames_away_max (Default: 5 pix); Max separation in frames (i.e., prop
% to time) for trajectory segments to be linked. Write 1 if you want no
% jumps (no frame jumps) in the segment linking. 
% For magnetic tweezers, consider how many frames might be dark when magnet
% is brought above sample.
% The remaining parameters can be left as they are. All parameters are
% saved in separate tabs to the excel file ending in "fullTrajs.xls".

% - 6a. Find trajectories and output them to one excel file in the current
% directory.
% Use function: 
% particle_results = FindTrajectsParticles(image_label,start_frame,end_frame,excludedRegions)
% For the input excludedRegions, used to to exclude certain regions from
% image, follow the example below.
% Regions to exclude are given by start coordinates (x_start, y_start) for
% top left corner of rectangle on image and end coordinates (x_end, y_end)
% for bottom right corner of rectangle. Coordinates delimit rectangular
% boxes on image and several rectangles can be defined. All rectangles are
% put together in structure input excludedRegions:
% Example, for Sonia's images we exclude the following regions:
excludedRegions.list_xstart = [1 1 1 130];
excludedRegions.list_xend = [820 112 100 215];
excludedRegions.list_ystart = [581 492 1 1];
excludedRegions.list_yend = [614 580 175 54];
% Then run:
t1 = FindTrajectsParticles(image_label,1,78,excludedRegions);
% Save all result structures in a .mat file:
save 'resultStructures' 't*' 

% - 6b. Link trajectory segments found into longer trajectories:
% linkTrajSegmentsParticles(image_label,start_frame,end_frame,particle_results,data_set_label).
% E.g., for frames 1 to 117 in video "210217r25.tif" in the current directory:
% linkTrajSegmentsParticles(image_label,1,117,t25,'tests'); 
% This generates an Excel file with all the trajectory data,
% "tests_25_fullTrajs.xls", in the current directory folder, for further
% analysis. This file contains two tabs with the parameters used in the
% analysis and another tab with the trajectory data. The key columns are
% Xcom, Ycom (x and y particle centre-of-mass positions in pixels) and
% TrajNumber, the Trajectory Number.
% Note: make sure that the start_frame and end_frame values are kept the
% same throughout all functions ran in steps 6a and 6b.
linkTrajSegmentsParticles(image_label,1,78,t1,'test');

% - 7. Plot and save a .png image of the Trajectory Numbers for the found
% particless overlaid on top of first frame of the video.
% Use function  plotParticleTrajNumbers(image_label,minPointsTraj).
% Input "minPointsTraj" (Default: 10) is the minimum number of points in a trajectory for
% it to be considered (we don't want to analyse tracks with only a few
% points that might appear due to suboptimal tracking when two particles are
% very close together):
plotParticleTrajNumbers(image_label,data_set_label,10)
cd(data_folder); % return to data folder.

% Note: the code does not work optimally for particles that are close together
% and have overlapping diffraction rings, or when there is a very uneven background. 

% - 8. Select good tracks. Three alternative methods 8a, 8b and 8c below.
% - 8a) Inspect tracks manually on a video to decide which to accept as
% good.
% Use function goThroughParticleTracksVideo with input showVideos=1:
% good_tracks = goThroughParticleTracksVideo(image_label,data_set_label,n_traj_start,n_traj_end,minPointsTraj,maxMajorAxisLength,showVideos)
good_tracks = goThroughParticleTracksVideo(image_label,'test',1,'end',6,50,1); 
% The above generates the structure:  
% good_tracks = 
%            image_label: '5mT-1Hz'
%           n_traj_start: 1
%             n_traj_end: 36
%          minPointsTraj: 6
%     good_track_numbers: [5 7 8 12]
%
% - 8b) Note, that doing the above is quite slow, particularly for long tracks
% from long videos. So the alternative is to inspect a few tracks as in 8a),
% make a note of the total number of long-enough tracks (printed on the
% command window during the analysis using goThroughParticleTracksVideo), then
% press Control+C to cancel the function evaluation (will generate no output), and 
% then generate the structure good_tracks2 by hand.
% Exclude particle numbers for particles that are close to each other just
% by looking at the png image generated in step 7.
% E.g., to generate by hand, do:
good_tracks.image_label = image_label;
good_tracks.n_traj_start = 1;
good_tracks.n_traj_end = 36;
good_tracks.minPointsTraj = 6;
good_tracks.maxMajorAxisLength = 50;
good_tracks.good_track_numbers = [5 7 8 12]; 
% Save result (as a .mat file, required for further analysis functions):
output_filename = strcat('good_tracks_',image_label);
save(output_filename,'good_tracks') % save variable good_tracks.
% NOTE: make sure you don't change the name 'good_tracks' to anything else.
% It is used later by function showManyParticleTrajAanalysis.m.
% 
% - 8c) When there are many particles found and long tracks, even 8b) is
% too tedious. In this case, use function goThroughParticleTracksVideo but
% set the input showVideos=0. This way, there is no video display or
% request for user input for each track and all tracks with more than
% "minPointsTraj" points and a max length below "maxMajorAxisLength" are 
% accepted as good tracks. The correct output variable "good_tracks"
% needed for the next step is generated.
% E.g. good_tracks = goThroughParticleTracksVideo(image_label,'test',1,'end',6,50,0); 


% - 9. Analyse each track separatedly.
% This is based on functions showParticleTrajAnalysis.m and
% showManyParticleTrajAnalysis.m
% Running the line below produces one analysis excel file and graph per track:
% processedManyTrajs = showManyParticleTrajAnalysis(image_label,data_set_label,n_traj_start,n_traj_end,start_frame,tsamp,pixelsize_nm,showVideo,saveAvi,minPointsTraj) 
processedTrajs = showManyParticleTrajAnalysis(image_label,'test',1,'end',1,0.0333,1,1,1,6);
% The analysis for each individual trajectory/track consists of:
% - Plot orientation angle (degrees) versus time; Calculate a postprocessed 
% angle and fit all values to a line to get angular velocity in degrees/s.
% - plot x and y (referenced to first point in track) versus time, 
% - plot trajectory on the x-y plane.
% - plot mean square displacement (msd) vs Delta t with error bars first
% all points, then only points with relative error < 150%. Fit the latter
% (tries fitting to a line and to a confined curve and saves results).
% - show the first frame with the trajectory overlaid.
% - print the trajectory number and other info, including angular velocity from fit.
% Saves all plots in a .PNG image file and saves also an EXCEL file with
% all results.
% IMPORTANT: note that the postprocessing of the raw angle assummes
% anti-clockwise rotation!! See function angleDegToPos.m.


% - 10. REANALYSE ANGLE AND ROTATION: If the data is noisy, it is likely
% that the rotation of the particle is not continuous and smooth. Perhaps
% there are several linear sections in the postprocessed continous angle
% plotted vs time, separated by noisy sections. It might be the case that
% the frame rate entered in step 9 needs to be corrected. 
% For any of these cases, functions reanalyseAngle.m and
% scriptToReanalyseAngle.m can be used. 
% The later reanalyses a set of excel files (generated in step 9) containing 
% the track data for frame number, time, angle, etc., to obtain angular
% velocity and frequency of rotation. The raw angle is postprocessed
% assuming anti-clockwise rotation (function angleDegToPos.m). Then the linear 
% sections with a large enough slope (above input thresh_slope) in the 
% the angle vs time plot are automatically detected and fitted to a line. 
% The input frame rate is used to convert frame number into time (s).
% IMPORTANT: calculations assume anti-clockwise rotation.
%
% To reanalyse a set of excel files in a given folder, first make sure
% Matlab's current directory is the directory containing the excel files to
% reanalyse. Then set the correct inputs in the script
% "scriptToReanalyseAngle.m" and run in the command window:
scriptToReanalyseAngle
%
% To reanalyse a single file, do e.g.:
reanalyseAngle('analysis_10mT,1Hz(1)__traj10.xls',15,130,5);
% reanalyseAngle(excelFileName,frameRateReal,thresh_slope,minSectionPoints)
% Inputs:
% - frameRateReal: frame rate for video data in frames per second.
% - thresh_slope: threshold slope in degrees/s. Only linear data sections with a
% positive slope larger than this threshold are fitted to obtain angular
% velocity and rotation frequency. 360deg/s corresponds to 1Hz. 
% For 1Hz field, a threshold~130deg/s works well.
% - minSectionPoints: minimum number of points in a single linear section in the angle vs time
% plot for it to be fitted to a line to obtain the slope (angular velocity).
% Outputs:
% A new folder called 'analysis_BIS' is generated, containing the
% following: 
% - A new png graph file is generated with name "bis" appended to original
% name, containing the new angle vs time plot, fits and fit results for individual linear sections. 
% - A new Excel file "bis" with new sheets is created with the
% revised analysis for the postprocessed angle vs time and for the data and
% fit results for the relevant linear sections.





%% To tests the methods on a single frame:

% % E.g., extract frame 1 from an image sequence selected from a browsing
% % dialog:
% frame1 = extract1frameB(1);
% 
% % Find candidate positions for particles:
% [x1,y1] = findCandidateParticlePositions(frame1,1);
% 
% % Eliminate candidate positions closer to each other than 10 pixels:
% [x1_b,y1_b,pos_to_keep1] = eliminateCoincidentPositions(x1,y1,10);
% 
% Plot:
% figure; imshow(frame1,[],'InitialMagnification',150); hold;
% plot(x1_b,y1_b,'o','Color','y','MarkerSize',14); hold off;
%
% Exclude certain regions from image. Regions to exclude are given by start coordinates (x_start, y_start)
% and end coordinates (x_end, y_end) that delimit rectangular boxes on image.
% For Sonia's images:
% list_xstart = [1 1 1 130];
% list_xend = [820 112 100 215];
% list_ystart = [581 492 1 1];
% list_yend = [614 580 175 54];
% [x1_c,y1_c] = excludeRegions(x1_b,y1_b,list_xstart,list_xend,list_ystart,list_yend);
%
% Plot:
% figure; imshow(frame1,[],'InitialMagnification',150); hold;
% plot(x1_c,y1_c,'o','Color','g','MarkerSize',14); hold off;
%
% % Test finding particle angle on single frame:
% r1 = findParticleAngle1frame(frame1,335,464,40,50);
% s1 = findParticleAngle1frame(frame1,x1_c(1),y1_c(1),50,60);