function particle_results = FindTrajectsParticles(image_label,start_frame,end_frame,excludedRegions)
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
% Function to find all particle trajectories and angles in an input image sequence.
%
% INPUTS:
% - image_label: string that labels a given image sequence found in current
% folder. The code finds the path of the image file automatically based on a string label 
% that is equal to the file name without the file extension. For example,
% for image video file "210217r25.tif", an appropriate label would be the
% string '210217r25'. This will be used throughout the entire RotTrack code.
% - start_frame: first frame of the sequence to be analysed.
% - end_frame: last frame of the sequence to be analysed. One can write
% 'end' if this is not known.
% - excludedRegions: to exclude certain regions from image. 
% Enter [] (so that isempty(excludedRegions) = 1) if you
% don't want to exclude any regions.
% Input a structure as follows to exclude one or several regions.
% Regions to exclude are given by start coordinates (x_start, y_start) for
% top left corner of rectangle on image and end coordinates (x_end, y_end)
% for bottom right corner of rectangle. Coordinates delimit rectangular
% boxes on image and several rectangles can be defined. All rectangles are
% put together in structure input excludedRegions:
% Example, for Sonia's images we exclude the following regions:
% excludedRegions.list_xstart = [1 1 1 130];
% excludedRegions.list_xend = [820 112 100 215];
% excludedRegions.list_ystart = [581 492 1 1];
% excludedRegions.list_yend = [614 580 175 54];
% 
% Example of how to call this function: 
% for frames 1 to 10 of image "Heiko_Thu Jun 24 2010_554.sif" in current folder:
% results1 = FindTrajectsParticles('554',1,10); 
% -------------------------------------------------------
% NOTE: before running this function you should move into the directory which
% contains the image sequence data (labelled by 'image_label').
%
% Reads image sequence data and then for each frame:
% finds candidate particle positions, eliminates coincidences for close candidates, 
% finds actual centre-of-mass particle position and particle orientation
% and accepts only particle centres not too close from edge.
% Results are saved in the structure array particle_results.
% Then links particle positions into trajectory segments: For second frame do differently
% and check only particles in 1st frame and compare to particles in second frame.
% For rest of frames (for frame k): A) first check loose particles (TrajNumber=0) two frames
% ago (k-2) and compare to particless in current frame (k).
% B) then check all particles in previous frame (k-1) and compare to particles in
% current frame (k).
% To decide on the best asignment of pairs of particles (link), we build up
% matrices of pair-wise distances, ratio of intensities and ratio of sigmas
% of all pairs of particles in the two frames being compared, and take the
% winning asignment as that with the smallest pairwise distance.
% For the angular tracking, the original angles obtained (between -90 and
% 90 deg) are left unmodified at this stage.
% 
% start_frame and end_frames are the frames through which the loop runs to
% find centres of particles.
%
% OUTPUT:
% The output, particle_results is a cell array with two elements:
% particle_results = {params, particle_final};
% The first element in the cell array, particle_results{1}, contains all parameters used to run
% the function: "params" (this is a structure array itself).
% The second element in the cell array, particle_results{2}, contains the track
% results: "particle_final".
% "particle_final", is itself a structure array with end_frame by L elements,
% each of which is a structure with fields: CentreX, CentreY, ClipFlag, TooCloseToEdge, FrameNumber and ParticleNumber.
% Along the dimension L, we have the different particle centres found on each
% frame (L will usually be the number of particle centres found in frame_start,
% or the largest number of particles ever accepted).
% For example, particle_final(100,8) is a structure with the above fields
% corresponding to the eighth found particle centre on frame 100.
%
% Example of params: particle_results{1}:
%
%             image_label: '25'
%             start_frame: 1
%               end_frame: 10
%      max_num_candidates: 200
%      subarray_halfwidth: 60
%     inner_circle_radius: 50
%                 rsq_min: 0.9000
%          d_coincid_cand: 3
%         d_coincid_found: 1
%                d_01_max: 30
%                d_02_max: 30
%                     rej: 10000
%               file_name: '210217r25.tif'
%               numFrames: 117
%             frame_Ysize: 1024
%             frame_Xsize: 1280
%
% particle_results1{2} is a struct array with fields:
%     'estimateXcentre'
%     'estimateYcentre'
%     'Xcom' % centre of mass
%     'Ycom'
%     'AngleDegrees'
%     'xpoints_ellipse'
%     'ypoints_ellipse'
%     'xpoints_majorAxis'
%     'ypoints_majorAxis'
%     'majorAxisLength'
%     'minorAxisLength'
%     'ClipFlag'
%     'TooCloseToEdge' 
%     'TrajNumber'
%     'FrameNumber'
%     'ParticleNumber'
% 
% See also end of function for more info on OUTPUT format.


%% DEFINITIONS and PARAMETERS:
% 
% 
% % One can save all parameter values in a file, e.g. "paramsForFindTrajects.m"
% % in the current directory. Just calling the name of that file loads the
% % parameter values into the workspace:
% paramsForFindTrajects
% % In this way, one can save different parameter sets for different data sets
% % in an easy manner and run two Matlabs at the same time working with different parameter sets.

% Alternatively, one can set PARAMETER values within this function below:  
%
% % Define data directory:
% % dir_data = 'Z:\Leake\Heiko Data\';
% % dir_data = 'C:\Isabel\ExperimData\HeikoData\';
% % Print data directory on command window to guide user:
% disp(' ') % empty line
% disp(['The data directory (.sif images) is: ',cd]) % display current directory.

% Maximum number of candidate particles (if eg. 260000 candidate positions are
% found, we get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB"), hence, we limit the
% max number of candidates:
max_num_candidates = 200; % should be around 200.
% Save parameters to results as structure "params":
params.max_num_candidates = max_num_candidates;

% PARAMETERS for finding particle centres (see findparticleCentre1frame.m, inputs to the function):
subarray_halfwidth = 25; % (Default: 60 pixels). Halfwidth of image square subarray which includes particle and background around it.
inner_circle_radius = 20; % (Default: 50 pixels). Radius of circular mask around particle. 
% Save parameters to results as structure "params":
params.subarray_halfwidth = subarray_halfwidth;
params.inner_circle_radius = inner_circle_radius;

% PARAMETERS for eliminating coincident positions:
d_coincid_cand = 10; % distance (in pixels) for eliminating coincidences in particle-position candidates. Default: 3 pixels.
d_coincid_found = 3; % distance for eliminating coincidences in found particle centres.
% Save parameters to results as structure "params":
params.d_coincid_cand = d_coincid_cand; % distance for eliminating coincidences in particle candidates.
params.d_coincid_found = d_coincid_found; % distance for eliminating coincidences in found particle centres.

% PARAMETERS for building trajectories:
% For linking particle centres in current and previous frames:
d_01_max = 30; % max distance in pixels between particle centres in current and previous frames, for linking them into a trajectory (Default: 30 pix).
% Save parameters to results as structure "params":
params.d_01_max = d_01_max; 

% For linking loose particle centres in current frame and 2 frames ago (jump of 1 frame in trajectory):
d_02_max = 30; % max distance in pixels between particle centres in current frame and 2 frames ago. (default: 30 pix).
% Save parameters to results as structure "params":
params.d_02_max = d_02_max; 

% Use a very large number (larger than image size in pixels) for rejected asignments:
rej = 10000;
params.rej = rej; % Save parameters to results as structure "params".

% Alternative way of selecting image sequence file:
% uigetfile opens a file dialog box to choose data file:
% [file_data,path_data] = uigetfile({'*.sif'}, 'Chose image data sequence:');
% strcat('data (.sif image):','  ',path_data,file_data)
% open a file dialog box to choose analysis file:
% [file_analysis,path_analysis] = uigetfile({'*.xls'}, 'Chose analysis file (trajectory):');
% strcat('analysis file (.xls trajectory):','  ',path_analysis,file_analysis)
% 
% disp(' ') % empty line
% disp(['The start frame for finding bright particle trajectories will be ',num2str(start_frame)]) % start_frame is an input.
% disp(['The end frame for finding bright particle trajectories will be ',num2str(end_frame)]) % end_frame is an input.
% -----------------------------------------------------


%% Read in the image-sequence data:

% Read image-sequence file: 
[numFrames frame_Ysize frame_Xsize image_data image_path] = extract_image_sequence_data(image_label);
% See "extract_image_sequence_data.m".
% numFrames is the number of frames in the image sequence.
% To get frame number "p" do: image_data(p).frame_data.
% Frame dimensions are frame_Ysize and frame_Xsize.
% --------------------------------------------------------------

if strcmp(end_frame,'end') % if input end_frame is 'end'
    end_frame = numFrames;
end

% Save input parameters to "params" structure (for output):
params.image_label = image_label;
params.start_frame = start_frame;
params.end_frame = end_frame;

% Save to parameters (for output):
params.file_name = image_path;
params.numFrames = numFrames;
params.frame_Ysize = frame_Ysize;
params.frame_Xsize = frame_Xsize;


%% Find candidate particle positions for start_frame (first frame), and find particle centres for those:

frame = image_data(start_frame).frame_data; % extract matrix data for first frame.
frame = double(frame);
disp(['frame number: ',num2str(start_frame)]) % print frame number to Command Window.

% Xpos is a matrix of the same size as frame, containing x values for all
% pixels and similarly for Ypos (used in future sections):
[Xpos,Ypos] = meshgrid(1:frame_Xsize,1:frame_Ysize);

[candidate_X_000,candidate_Y_000] = findCandidateParticlePositions(frame,1); % Second input: use method 1, which seems to work better.
% candidate_X_000 and candidate_Y_000 are two column vectors of the same
% length containing the x and y coordinates of the candidate particle positions found on the image.
% They contain integer numbers: coordinates or pixel numbers which give
% position on image plane.

disp(['no. of initial particle candidate positions: ',num2str(length(candidate_X_000))])

% Exclude certain regions from image (see INPUT excludedRegions):
if isempty(excludedRegions) == 1 % do not exclude any regions of image
    candidate_X_00 = candidate_X_000;
    candidate_Y_00 = candidate_Y_000;
else
    list_xstart = excludedRegions.list_xstart;
    list_xend = excludedRegions.list_xend;
    list_ystart = excludedRegions.list_ystart;
    list_yend = excludedRegions.list_yend;
    [candidate_X_00,candidate_Y_00] = excludeRegions(candidate_X_000,candidate_Y_000,list_xstart,list_xend,list_ystart,list_yend);
end

disp(['no. of total candidate particle positions after excluding certain regions (as input)): ',num2str(length(candidate_X_00))])


% Error control:
% Limit the max number of candidate particles (if eg. 260000 candidate particles are
% found, we will get an error in function pdist: "Distance matrix has more
% elements than the maximum allowed size in MATLAB").
% Select only the first max_num_candidates then.
if length(candidate_X_00) > max_num_candidates
    candidate_X_00 = candidate_X_00(1:max_num_candidates);
    candidate_Y_00 = candidate_Y_00(1:max_num_candidates);
    disp(['NOTE!! no. of candidate particle positions has been limited to ',num2str(max_num_candidates)])
end

% % Check graphically:
% imshow(frame,[]);
% hold on;
% plot(candidate_X_00,candidate_Y_00,'*');
% figure;

% Eliminate candidate positions closer to each other than d_coincid_cand (3 pixels):
[candidate_X_0,candidate_Y_0,pos_to_keep1] = eliminateCoincidentPositions(candidate_X_00,candidate_Y_00,d_coincid_cand);

disp(['no. of total candidate particle positions after eliminating coincidences: ',num2str(length(candidate_X_0))])

% No need to do background subtraction.
% % Subtract background from entire frame (method 1 is faster) considering particle candidate positions:
% % frameNoBgnd = removeBgnd(frame,candidate_X_0,candidate_Y_0,inner_circle_radius,subarray_halfwidth,1);
% % frame_to_search = frameNoBgnd; % Initialise frame to search for particle centres.
frame_to_search = frame; % Initialise frame to search for particle centres.

% Find particle centres and angles for first frame and decide if we accept them or not:
n =1; % Initialise index n (index for accepted particle centre positions):

% Find particle centres (centre of mass) and angles:
for m = 1:size(candidate_X_0,1) % loop through all candidate particle positions.
    % Now find particle centre and orientation using function findParticleAngle1frame:
    % use candidate particle positions as initial estimates and then refine to find particle centre of mass.   
    particle_result = findParticleAngle1frame(frame_to_search,candidate_X_0(m),candidate_Y_0(m),inner_circle_radius,subarray_halfwidth);
    particle_result.FrameNumber = start_frame; % Add new field containing frame number (time) to result structure.
    
    % Accepted particle centres:
    % if (particle_result.ClipFlag == 0 && ... % use this to exclude particles close to image edge.
    if (particle_result.ClipFlag < 2 && ... % use this to accept all particles;
            particle_result.TooCloseToEdge == 0)
            % Only accept and save result of found particle centre if clipping flag = 0 and if values of rsquare of fits are acceptable.
            particle_result.ParticleNumber = n; % Add new field containing particle number to result structure.
            particle_final(start_frame,n) = particle_result; % store "good" found particle centres.
            % This is also saved in the final result particle_final, structure array.
            % first index is for frame number, second index is for particle number.
            n = n+1; % advance index n for accepted particle centres.
    end
end

% % display the number of accepted particle centres for this frame:
disp(['no. of accepted particle centres in first frame: ',num2str(n-1)])

% Convert results of found particle-centre positions to a useful form :

% Error control: if no particle centres were accepted:
if (n-1) == 0 
    found_particle_CentreX = [];
    found_particle_CentreY = [];
    % Need to create the whole particle_final structure with all its fields
    % here (same fields as output from findParticleAngle1frame, plus fields
    % added in this function), just in case the number of accepted particles in
    % the first frame is zero, in order not to get error: "Subscripted
    % assignment between dissimilar structures".
    % Save empty particle (we need this, otherwise if in the last frame the
    % no. of accepted particles is 0, there will be no result
    % particle_final(end_frame,:) and the following functions will fail).
    particle_final(start_frame,n).estimateXcentre = [];
    particle_final(start_frame,n).estimateYcentre = [];
    particle_final(start_frame,n).Xcom = [];
    particle_final(start_frame,n).Ycom = [];
    particle_final(start_frame,n).AngleDegrees = [];
    particle_final(start_frame,n).xpoints_ellipse = [];
    particle_final(start_frame,n).ypoints_ellipse = [];
    particle_final(start_frame,n).xpoints_majorAxis = [];
    particle_final(start_frame,n).ypoints_majorAxis = [];
    particle_final(start_frame,n).majorAxisLength = [];
    particle_final(start_frame,n).minorAxisLength = [];
    particle_final(start_frame,n).ClipFlag = [];
    particle_final(start_frame,n).TooCloseToEdge = []; 
    particle_final(start_frame,n).FrameNumber = [];
    particle_final(start_frame,n).ParticleNumber = []; 
    particle_final(start_frame,n).TrajNumber = [];
else
    % For plotting found centre of mass:
    found_particle_CentreX = [particle_final(start_frame,:).Xcom]'; % column vector with found Xcom (centre of mass) positions of all particles.
    found_particle_CentreY = [particle_final(start_frame,:).Ycom]'; % column vector with found Ycom positions of all particles.
end

% Check graphically:
figure;
imshow(frame,[]);
hold on;
plot(found_particle_CentreX,found_particle_CentreY,'x','Color','g','MarkerSize',7) % plot accepted particle centres in green.
% For plotting 3 points in major axis (as a line):
for ind = 1:(n-1) % loop through all accepted positions.
    xpoints_majorAxis = particle_final(start_frame,ind).xpoints_majorAxis;
    ypoints_majorAxis = particle_final(start_frame,ind).ypoints_majorAxis;
    plot(xpoints_majorAxis,ypoints_majorAxis,'y','LineWidth',1);
end
pause(0.05); % (this pause is needed to give time for the plot to appear.)
hold off;
% -----------------------------------------------------------------------



%% Loop through selected frames:

tr =1; % initialise trajectory index.

for k = (start_frame+1):end_frame
    % to go through all frames do instead: for k = 1:length(sifData)
    
    frame = image_data(k).frame_data; % extract frame data which is stored in field 'frame_data'.
    frame = double(frame);
    
    imshow(frame,[],'Border','tight','InitialMagnification',150); % show image scaled between its min and max values ([]).
    hold on;
    
    disp(['frame number: ',num2str(k)]) % print frame number to Command Window.
    
    %-------------------------------------
    % Find new candidate particle positions for this frame:      
    [candidate_X_000,candidate_Y_000] = findCandidateParticlePositions(frame,1); % Second input: use method 1, which seems to work better.
    % the subindex "_000" in candidate_X_000 indicates newly found particle
    % candidates for the current frame. 
      
    disp(['no. of initial particle candidate positions: ',num2str(length(candidate_X_000))])
    
    % Exclude certain regions from image (see INPUT excludedRegions):
    if isempty(excludedRegions) == 1 % do not exclude any regions of image
        candidate_X_00 = candidate_X_000;
        candidate_Y_00 = candidate_Y_000;
    else
        [candidate_X_00,candidate_Y_00] = excludeRegions(candidate_X_000,candidate_Y_000,list_xstart,list_xend,list_ystart,list_yend);
    end    
    disp(['no. of total candidate particle positions after excluding certain regions (as input)): ',num2str(length(candidate_X_00))])
      
    % Error control:
    % Limit the max number of candidate particle positions (if eg. 260000 candidates are
    % found, we will get an error in function pdist: "Distance matrix has more
    % elements than the maximum allowed size in MATLAB").
    % Select only the first max_num_candidates then.
    if length(candidate_X_00) > max_num_candidates
        candidate_X_00 = candidate_X_00(1:max_num_candidates);
        candidate_Y_00 = candidate_Y_00(1:max_num_candidates);
        disp(['NOTE!! no. of candidate particle positions has been limited to ',num2str(max_num_candidates)])
    end
           
    % Eliminate candidate positions closer to each other than d_coincid_cand (3 pixels):
    [candidate_X_0,candidate_Y_0,pos_to_keep1] = eliminateCoincidentPositions(candidate_X_00,candidate_Y_00,d_coincid_cand);
    
    disp(['no. of total candidate particle positions after eliminating coincidences: ',num2str(length(candidate_X_0))])
   
    % % Check graphically:
    % imshow(frame,[]);
    % hold on;
    % plot(candidate_X_0,candidate_Y_0,'*');
    % figure;
    
    % No need to do background subtraction.
    % % Subtract background from entire frame (method 1 is faster) considering particle candidate positions:
    % % frameNoBgnd = removeBgnd(frame,candidate_X_0,candidate_Y_0,inner_circle_radius,subarray_halfwidth,1);
    % % frame_to_search = frameNoBgnd; % Initialise frame to search for particle centres.
    frame_to_search = frame; % Initialise frame to search for particle centres.
    
    % Find particle centre of mass and angles and decide if we accept them or not:
    n =1; % Initialise index n (index for accepted particle centre positions):
    
    for m = 1:size(candidate_X_0,1) % loop through all candidate particle positions.
        % Now find particle centre and orientation using function findParticleAngle1frame:
        % use candidate particle positions as initial estimates and then refine to find particle centre of mass.        
        particle_result = findParticleAngle1frame(frame_to_search,candidate_X_0(m),candidate_Y_0(m),inner_circle_radius,subarray_halfwidth);              
        % index k is for frame number, index m is for particle number
        particle_result.FrameNumber = k; % Add new field containing frame number (time) to result structure.
        
        % Accepted particle centres:
        % if (particle_result.ClipFlag == 0 && ... % to exclude particles close to edge
        if (particle_result.ClipFlag < 2 && ... % use this to accept all particles;
                particle_result.TooCloseToEdge == 0)
            % Only accept and save result of found particle centre if clipping flag = 0 and if values of rsquare of fits are acceptable.
            
            particle(k,n) = particle_result; % store accepted found particle centres in this preliminary result.
            % first index is for frame number, second index is for particle number.
            %    %--------------------------------------------------
            %    plot(particle(k,n).Xcom,particle(k,n).Ycom,'o','Color','r','MarkerSize',10); % Plot found particle centres in red
            %    %--------------------------------------------------
                      
            n = n+1; % advance index n for accepted particle centres.
        end
    end
    
    % % display the number of accepted particle centres for this frame:
    disp(['no. of accepted particle centres: ',num2str(n-1)])
    
    % % The following two lines are used together with the previous two
    % "plot" and "imshow" (commented off) lines:
    pause(0.1); % this pause is needed to give time for the plot to appear
    %    hold off;
    
    % Convert results of found particle-centre positions to a useful form:
    
    if (n-1) == 0 % Error control: if no particles were accepted.
        found_particle_CentreX = [];
        found_particle_CentreY = [];
        particle_final(k,n).ParticleNumber = []; % Save empty particle (we need this, otherwise if in the last frame the no. of accepted particles is 0, there will be no result particle_final(end_frame,:) and the following functions will fail).
    else
        found_particle_CentreX = [particle(k,:).Xcom]'; % column vector with found CentreX positions of all candidate particle centres.
        found_particle_CentreY = [particle(k,:).Ycom]'; % column vector with found CentreY positions of all candidate particle centres.
        
        %-------------------------------
        % This is not really necessary but I keep it:
        % Eliminate coincidences in result of last found particle centres for a given frame (for distance <1):
        [found_particle_CentreX found_particle_CentreY pos_final] = eliminateCoincidentPositions(found_particle_CentreX,found_particle_CentreY,d_coincid_found);
        % pos_final contains positions of selected, kept particle centres.
        %-------------------------------
        
        % Save final particle centres to variable particle_final:
        n=1; % index for final kept particles.
        for ii = 1:length(pos_final)
            mientras = particle(k,pos_final(ii)); % intermediate result.
            mientras.ParticleNumber = n; % Add new field containing particle number to result structure.
            if k > (start_frame+1)
                mientras.TrajNumber = []; % to avoid error of disimilar structures
            end
            particle_final(k,n)=mientras; % final result structure of accepted particle centres.
            n = n+1;
        end
        % Plot found particle centres:
        pause(0.5);
        plot(found_particle_CentreX,found_particle_CentreY,'x','Color','g','MarkerSize',7) % plot final accepted particle centres as green cross.
        % pause(0.1); % this pause is needed to give time for the plot to appear
        % For plotting 3 points in major axis (as a line):
        for ind = 1:(n-1) % loop through all accepted positions.
            xpoints_majorAxis = particle_final(k,ind).xpoints_majorAxis;
            ypoints_majorAxis = particle_final(k,ind).ypoints_majorAxis;
            plot(xpoints_majorAxis,ypoints_majorAxis,'y','LineWidth',1);
        end
        pause(0.05); % this pause is needed to give time for the plot to appear
        hold off;
        
        disp(['no. of final found particle centres after eliminating coincidences: ',num2str(length(found_particle_CentreX))])
    end
     
    
%% --------------------
    % LINKING PARTICLE POSITIONS INTO TRAJECTORY SEGMENTS:
    
    % Link found and accepted particle positions into trajectory segments:
    
    % Trajectory index tr is initialised to 1 outside the loop through frames (k loop).
    
    % Do differently FOR SECOND FRAME (k == start_frame+1): compare only accepted particles in
    % previous and current frames:
    if k == start_frame+1 && ... % If second frame and
            (n-1)~=0 && ... % if the number of accepted particle centres is not zero and
            isempty([particle_final(k-1,:).ParticleNumber])==0 && ... % at least 1 accepted particle in previous frame and
            isempty([particle_final(k,:).ParticleNumber])==0 % at least 1 accepted particle in current frame.
        % There are no trajectories jet, so compare accepted particles in previous and current frames:
        N0 = max(cat(1,particle_final(k-1,:).ParticleNumber));  % no. of accepted particles in previous frame.
        % Note: cat(1,particle_final(k-1,:).ParticleNumber) gives a column vector with the values of particleNumber for all non-empty accepted particles in frame k-1.
        N1 = max(cat(1,particle_final(k,:).ParticleNumber));  % no. of accepted particles in current frame.
        
        % Create cell arrays with empty elements to pre-asign sizes:
        d01 = cell(N0,N1); % Note: d01 is a cell array (matrix) but d_01 below is a scalar.
        
        for q0 = 1:N0 % loop though accepted particle centres in previous frame.
            for q1 = 1:N1 % loop though accepted particle centres in current frame.
                % d_01: distance between particle centres in previous and current frames:
                d_01 = sqrt((particle_final(k-1,q0).Xcom-particle_final(k,q1).Xcom)^2+(particle_final(k-1,q0).Ycom-particle_final(k,q1).Ycom)^2);
                %                 d_01
              
                % Accept and save trajectory if particle centres in previous and
                % current frames fulfill the following conditions:
                if d_01 < d_01_max   % see PARAMETERS at start of this function.           
                    % Asign accepted values to cell array elements to store them:
                    d01{q0,q1} = d_01; % use {} for cell arrays.
                else % rejected asignments:
                    d01{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                end
            end
            
%                         d01
%                         [d01{q0,:}]
            % Note that [d01{q0,:}] gives only non-empty elements of row q0
            % in the cell array d01 as a row vector, that's why we had to
            % give a numeric value rej to non-accepted asignments.
            
            % Note that if all asignments in previous step are rejected,
            % [d01{q0,:}] will be a list of rej values, and its minimum will
            % be rej.
            % If list of "linkable" particles, [d01{q0,:}], has no accepted
            % asignments (all values are rej):
            if min([d01{q0,:}]) == rej
                % Asign trajectory number 0 to the particle centre in the previous frame only:
                particle_final(k-1,q0).TrajNumber = 0;
            else % if there is at least one accepted asignment for a given particle centre in the previous frame:
                
                % Decide of all possible accepted particle centres (in current frame)
                % that could be linked to particle centre q0 in previous frame, which one is the best:
                % We take the best as the closest one to particle q0:
                q1_chosen = find([d01{q0,:}] == min([d01{q0,:}])); % find position of the minimum pair-wise distance.
                
                % ERROR CONTROL (on rare occasions, two positions are found above):
                if length(q1_chosen) > 1
                    q1_chosen = q1_chosen(1);
                end
%                             q1_chosen
                
                % Check if there is a better competing asignment for a given particle centre q1 in the current
                % frame from another particle centre in the previous frame.
                % Hence, check also column-wise in matrix d01 to avoid asigning a traj
                % number to a particle centre q1 in the current frame that had already
                % had a traj number asigned to it linking it to a different particle centre q0 in
                % the previous frame, which might be at a shorter distance
                % from it than the current one.
                
%                             [d01{:,q1_chosen}] % chosen column of d01 matrix of distances.
                
                % Asign trajectory numbers to structure particle_final:
                % If the found distance in that column is not the minimum one:
                if q0 ~= find([d01{:,q1_chosen}] == min([d01{:,q1_chosen}]));
                    particle_final(k-1,q0).TrajNumber = 0; % asign trajectory number 0 to particle centre in previous frame.
                else
                    particle_final(k-1,q0).TrajNumber = tr; % asign trajectory number to particle centre in previous frame, to particle_final structure.
                    particle_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to particle centre in current frame.
                    tr = tr+1; % advance trajectory-number index.
                end
            end
        end
        
        
        
    else % for FRAMES k >= start_frame+2, from third chosen frame on:
        
        % A) Compare loose particle centres (TrajNumber is 0) two frames ago (k-2)
        % to found particle centres in current frame (TrajNumber is []):
        % XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        % DO maybe!!
        if k ~= start_frame+1 && (n-1)~=0 && ... % If the number of accepted particle centres is not zero and
                isempty([particle_final(k-2,:).ParticleNumber])==0 && ... % at least 1 accepted particle centre 2 frames ago and
                isempty([particle_final(k,:).ParticleNumber])==0 % at least 1 accepted particle centre in current frame.
            
            N0 = max(cat(1,particle_final(k-2,:).ParticleNumber));  % no. of accepted particle centres 2 frames ago.
            N1 = max(cat(1,particle_final(k,:).ParticleNumber));  % no. of accepted particle centres in current frame.
            
            % Create cell arrays with empty elements to pre-asign sizes.
            d02 = cell(N0,N1); % Note: d02 is a cell array (matrix) but d_02 below is a scalar.
            
            for q0 = 1:N0 % loop though loose accepted particle centres 2 frames ago.
                if particle_final(k-2,q0).TrajNumber == 0 % only for loose (unlinked) particle centres (TrajNumber=0) two frames ago (so only rows q0 in matrix d01 which have unlinked particle centres will fill up).
                    for q1 = 1:N1 % loop though accepted particle centres in current frame.
                        % d_01: distance between particle centres in previous and current frames:
                        d_02 = sqrt((particle_final(k-2,q0).Xcom-particle_final(k,q1).Xcom)^2+(particle_final(k-2,q0).Ycom-particle_final(k,q1).Ycom)^2);
                        
                        %                         d_02
                        
                        % Accept and save trajectory if particle centres in previous and
                        % current frames fulfill the following conditions:
                        if d_02 < d_02_max   % see PARAMETERS at start of this function.                                
                            % Asign accepted values to cell array elements to store them:
                            d02{q0,q1} = d_02; % use {} for cell arrays.
                        else % rejected asignments:
                            d02{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                        end
                    end
                    
%                                         d02
%                                         [d02{q0,:}]
                    % Note that [d02{q0,:}] gives only non-empty elements of the cell array d01 as a row vector.
                    
                    % Note that if all asignments in previous step are rejected,
                    % [d02{q0,:}] will be a list of rej values, and its minimum will be rej.
                    % If list of "linkable" particles, [d02{q0,:}], has no accepted asignments (all values are rej):
                    if min([d02{q0,:}]) == rej
                        % Asign trajectory number 0 to the particle centre two frames ago only:
                        particle_final(k-2,q0).TrajNumber = 0; % point stays loose (unlinked).
                        
                    else % if there is at least one accepted asignment for a given particle centre two frames ago:
                        
                        % Decide of all possible accepted/saved particle centres (in current frame)
                        % that could be linked to particle centre q0 in previous frame (of all possible asignments), which one is the best:
                        % We take the best as the closest one to particle centre q0:
                        q1_chosen = find([d02{q0,:}] == min([d02{q0,:}])); % find position of the minimum pair-wise distance.
                        
%                                             q1_chosen
                        
                        % Check if there is a better competing asignment for a given particle centre q1 in the current
                        % frame from another particle centre q0 two frames ago.
                        % Hence, check also column-wise in d01 to avoid asigning a traj
                        % number to a particle centre q1 in the current frame that had already
                        % had a traj number asigned to it linking it to a different particle centre q0 two frames ago which might be at a shorter distance
                        % from it than the current one.
                        
%                                             [d02{:,q1_chosen}] % chosen column of d01 matrix of distances.
                        
                        % Asign trajectory numbers to structure particle_final:
                        % If the found distance in that column is not the minimum one:
                        if q0 ~= find([d02{:,q1_chosen}] == min([d02{:,q1_chosen}]));
                            % Asign trajectory number 0 to the particle centre two frames ago only:
                            particle_final(k-2,q0).TrajNumber = 0; % point stays loose (unlinked).
                        else
                            particle_final(k-2,q0).TrajNumber = tr; % asign trajectory number to particle centre two frames ago, to particle_final structure.
                            particle_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to particle centre in current frame.
                            tr = tr+1; % advance trajectory-number index.
                        end
                    end
                end
            end
        end
        
        
        % B) Compare loose particle centres (TrajNumber is []) and trajectories (TrajNumber is >0)
        % in previous frame (k-1) to found particle centres in current frame:
        if (n-1)~=0 && ... % If the number of accepted particle centres is not zero and
                isempty([particle_final(k-1,:).ParticleNumber])==0 && ... % at least 1 accepted particle centre in previous frame.
                isempty([particle_final(k,:).ParticleNumber])==0 % at least 1 accepted particle centre in current frame.
            
            % zzzzzzzzzzzzzzzzzzzzzzz
            N0 = max(cat(1,particle_final(k-1,:).ParticleNumber));  % no. of accepted particle centres in previous frame.
            N1 = max(cat(1,particle_final(k,:).ParticleNumber));  % no. of accepted particle centres in current frame.
            
            % Create cell arrays with empty elements to pre-asign sizes.
            d01 = cell(N0,N1);
            
            for q0 = 1:N0 % loop though accepted particle centres in previous frame.
                for q1 = 1:N1 % loop though accepted particle centres in current frame.
                    % d_01: distance between particle centres in previous and current frames:
                    d_01 = sqrt((particle_final(k-1,q0).Xcom-particle_final(k,q1).Xcom)^2+(particle_final(k-1,q0).Ycom-particle_final(k,q1).Ycom)^2);
                    
                    %                     d_01
                    
                    % Accept and save trajectory if particle centres in previous and
                    % current frames fulfill the following conditions:
                    if d_01 < d_01_max   % see PARAMETERS at start of this function.
                        % Asign accepted values to cell array elements to store them:
                        d01{q0,q1} = d_01; % use {} for cell arrays.
                    else % rejected asignments:
                        d01{q0,q1} = rej; % Use rej for asignments not accepted (images usually 512x512arrays, so rej pix is an impossibly large distance, this is why it is chosen here).
                    end
                end
                
%                                 d01
%                                 [d01{q0,:}] % last row of d01 matrix of distances.
                % Note that [d01{q0,:}] gives only non-empty elements
                % of row q0 in the cell array d01, as a row vector.
                
                % Note that if all asignments in previous step are rejected,
                % [d01{q0,:}] will be a list of rej values, and its minimum will
                % be rej.
                % If list of "linkable" particles, [d01{q0,:}], has no accepted
                % asignments (all values are rej):
                if min([d01{q0,:}]) == rej
                    
                    if isempty(particle_final(k-1,q0).TrajNumber) % if point in previous frame was not part of a trajectory (TrajNumber=[]):
                        % Asign trajectory number 0 to the particle centre in the previous frame only:
                        particle_final(k-1,q0).TrajNumber = 0;                       
                    end
                    
                else % if there is at least one accepted asignment for a given particle centre two frames ago:
                    
                    % Decide of all possible accepted/saved particle centres (in current frame)
                    % that could be linked to particle centre q0 in previous frame, which one is the best:
                    % We take the best as the closest one to particle q0:
                    q1_chosen = find([d01{q0,:}] == min([d01{q0,:}])); % find position of the minimum pair-wise distance.
                    
                    % ERROR CONTROL (on rare occasions, two positions are found above):
                    if length(q1_chosen) > 1
                        q1_chosen = q1_chosen(1);
                    end
%                                     q1_chosen
                    
                    % Check if there is a better competing asignment for a given particle centre q1 in the current
                    % frame from another particle centre q0 in the previous frame.
                    % Hence, check also column-wise in d01 to avoid asigning a traj
                    % number to a particle centre q1 in the current frame that had already
                    % had a traj number asigned to it linking it to a different particle centre q0 in
                    % the previous frame which might be at a shorter distance
                    % from it than the current one.
                    
%                                     [d01{:,q1_chosen}]  % chosen column of d01 matrix of distances.
                    
                    % If the found distance in that column is not the minimum one:
                    if q0 ~= find([d01{:,q1_chosen}] == min([d01{:,q1_chosen}]));
                        if isempty(particle_final(k-1,q0).TrajNumber) % if point in previous frame was not part of a trajectory (TrajNumber=[]):
                            % Asign trajectory number 0 to the particle centre in the previous frame only:
                            particle_final(k-1,q0).TrajNumber = 0;
                        end
                    else
                        % Asign trajectory numbers to structure particle_final:
                        if particle_final(k-1,q0).TrajNumber > 0 % if point in previous frame was already part of a trajectory:
                            particle_final(k,q1_chosen).TrajNumber = particle_final(k-1,q0).TrajNumber; % asign that trajectory number to particle centre in current frame.
                        else % if point in previous frame was not part of a trajectory:
                            particle_final(k-1,q0).TrajNumber = tr; % asign new trajectory number to particle centre in previous frame.
                            particle_final(k,q1_chosen).TrajNumber = tr; % asign same trajectory number to particle centre in current frame.
                            tr = tr+1; % advance trajectory-number index.
                        end
                    end
                end
            end
            % zzzzzzzzzzzzzzzzzzzzzzz
        end
    end
    %-----------------
    
    
end  % loop through selected frames



%% OUTPUT OF particle-FINDING PROCESS: final output particle_results:
%
particle_results = {params, particle_final};
% params is a structure array containing all parameters used to run the
% function.
% particle_final, is a structure array with end_frame x L elements,
% each of which is a structure with fields:
%     'estimateXcentre'
%     'estimateYcentre'
%     'Xcom' % centre of mass
%     'Ycom'
%     'AngleDegrees'
%     'xpoints_ellipse'
%     'ypoints_ellipse'
%     'xpoints_majorAxis'
%     'ypoints_majorAxis'
%     'majorAxisLength'
%     'minorAxisLength'
%     'ClipFlag'
%     'TooCloseToEdge' 
%     'TrajNumber'
%     'FrameNumber'
%     'ParticleNumber'
%  
% Along the dimension L, we have the different particle centres found on each
% frame (L will often be the number of particle centres found in frame_start,
% it is always the largest number of particles ever accepted on one frame).
% For example, particle_final(100,8) is a structure with the above fields
% corresponding to the eighth found particle centre on frame 100.
%
% Note that even if we only analyse from start_frame to end_frame,
% particle_final is a list containing empty structure arrays from index 1 to
% index start_frame, and then the found particle centres for the analysed frames start_frame to end_frame.
%
% The result is padded to a fixed number of particle centre structures for each
% frame (the maximum no. of accepted found particle centres of all frames), so that for a given
% frame in which less found particle centres have been accepted, the remaining
% elements are padded with empty structures with empty fields [].
% To check if a given particle centre is empty: isempty(particle_final(101,10).Xcom)
% gives 1 if field "CentreX" of tenth particle centre found and accepted in frame 101
% is empty (equal to []).
%
% e.g. cat(1,particle_final(100,:).SigmaFit) gives a vector column with all the
% non-empty SigmaFit values for all particle centres accepted in frame 100.
