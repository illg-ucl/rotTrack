function particle_result = findParticleAngle1frame(frame,x_estimate,y_estimate,inner_radius,subarray_halfsize)
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
% Function to find the centre-of-mass position and orientation angle of a particle on an image frame.
%
% INPUTS:
% - frame is a matrix containing the image data. Can be obtained doing, e.g.:
% frame1 = extract1frameB(1);
%
% - x_estimate, y_estimate are the initial estimated centre positions fed to the
% iterative method. These can be obtained as the output of function 
% findCandidateParticlePositions. Example:
% [x1,y1] = findCandidateBeadPositions(frame1,1);
% x_estimate = x1(1);  y_estimate = y1(1); 
%
% - inner_radius: radius of inner circular mask in pixels, used to carry
% out a second background subtraction if necessary (usually not). Make sure this is smaller than
% subarray_halfsize. Should include entire particle and its diffraction
% rings.
%
% - subarray_halfsize: the algorithm selects this number of pixels above and below and to
% left and right of estimated centre of particle to form the selected image
% subarray (I). This depends on the particle size, must be larger than
% particle size and include some background around it, default~35-70.
%
%
% OUTPUTS: values resulting from the iterative method. The output is a
% structure "particle_result" with the following fields (see end of file).
%
% Example of how to call this function: s1 = findParticleAngle1frame(frame1,271,364,50,60),
% where frame1 is the image array (single-frame data). Then to obtain results
% from the structure s1, do, for example: s1.estimateXcentre, s1.estimateYcentre, etc.

%-----------------
% PARAMETERS:
d_min = 10; % minimum subarray halfsize in pixels. If candidate bead centre too close to edge, return empty result.
% Additional background subtraction: 1 for yes, 0 for no.
additional_bgnd_subtract = 1; % Default = 0 if background-subtracted image is used as input.
% Thresholding parameter:
threshold_factor = 0.6; % Default: 1. Reduce to e.g. 0.6, to avoid overestimating particle size when particles have blurry edges.
%-----------------

% % Aid plot: plot original frame and estimate centre position
% % (x_estimate,y_estimate) overlayed as red circle on top:
% figure;
% imshow(frame,[],'InitialMagnification',150);
% hold;
% plot(x_estimate,y_estimate,'o','Color','r','MarkerSize',12);
% title('current bead with centre estimate')

d = subarray_halfsize; % subarray halfsize. d needs to be at least equal to the radius of the inner mask, inner_radius. Default: d=8, inner_radius=5.


%% Error control:
if subarray_halfsize < inner_radius
    error('findBeadCentre1frame:one','check input parameters: "subarray_halfsize" cannot be smaller than "inner_radius"')
end

% If the particle candidate is at the edge of the image, flag it with a clipping flag and take a subarray around
% it with a smaller size, as large as possible until the edge is reached, so reasign d:
clipping_flag = 0;
tooCloseToEdge = 0;
d_top = subarray_halfsize;
d_bottom = subarray_halfsize;
d_left = subarray_halfsize;
d_right = subarray_halfsize;
% If the particle is at edge of image, take a smaller subarray around
% it but as large as possible until the edge is reached,
% so reasign the d values:
if (round(y_estimate)-d_top)<1
    d_top = round(y_estimate) - 1;
end
if (round(y_estimate)+d_bottom)>size(frame,1)
    d_bottom = size(frame,1)-round(y_estimate);
end
if (round(x_estimate)-d_left)<1
    d_left = round(x_estimate) - 1;
end
if (round(x_estimate)+d_right)>size(frame,2)
    d_right = size(frame,2)-round(x_estimate);
end
% Chose the minimum distance:
d = min([d_top d_bottom d_left d_right]);

if d < inner_radius
    clipping_flag = 1;
end

if d < d_min
    % Error control: return empty result structure as outuput if particle-centre candidate is closer
    % to edge than d_min:
    tooCloseToEdge = 1;    
    
    particle_result.estimateXcentre = x_estimate; % x-centre estimate input.
    particle_result.estimateYcentre = y_estimate; % y-centre estimate input.
    particle_result.Xcom = []; % x centre-of-mass position within entire frame.
    particle_result.Ycom = []; % y centre-of-mass position within entire frame.
    particle_result.AngleDegrees = [];
    particle_result.xpoints_ellipse = [];
    particle_result.ypoints_ellipse = [];
    particle_result.xpoints_majorAxis = [];
    particle_result.ypoints_majorAxis = [];
    particle_result.majorAxisLength = [];
    particle_result.minorAxisLength = [];
    particle_result.ClipFlag = clipping_flag; % 1 if candidate was closer to edge of image than inner_radius.
    particle_result.TooCloseToEdge = tooCloseToEdge; % 1 if particle candidate was closer to edge of image than d_min.
    particle_result.Area = [];
    
else
         
    %% Create image subarray (I) around particle estimated centre:
    % Create image subarray of size (2*d+1)x(2*d+1) centered on (x_estimate,y_estimate) centroid estimate:
    I = frame(round(y_estimate)-d:round(y_estimate)+d,round(x_estimate)-d:round(x_estimate)+d);
    % % Aid plot:
    % figure;
    % imshow(I,[])
    
    % Create matrices containing the x and y positions corresponding to the
    % subarray I. Ys is a matrix of the same size as I, containing y values for
    % all pixels and similarly for Xs. In Ys, y varies along first dimension
    % (vertical) in matrix. In Xx, x varies along second dimension (horiz).
    [Ys,Xs] = ndgrid(round(y_estimate)-d:round(y_estimate)+d,round(x_estimate)-d:round(x_estimate)+d);
    
    
    %% Additional background subtraction:
    if additional_bgnd_subtract == 1 % I have removed this from input list but leave code here for now.
        % Generate inner circular mask for subarray, to calculate average
        % background per pixel and subtract it:
        inner_mask = zeros(2*d+1);
        % Assign value 1 to pixels for which radius <= inner_radius. C = hypot(A,B) returns SQRT(ABS(A).^2+ABS(B).^2).
        inner_mask = inner_mask | hypot(Xs-x_estimate, Ys-y_estimate) <= inner_radius;
        % Background mask is the negative of the inner circle mask, i.e.,
        % not-inner_mask. It has zeros in centre and ones in surrounding pixels:
        bgnd_mask = double(~inner_mask); % make of class double (not binary) for what follows.
        Ibgnd = I.*bgnd_mask; % bgnd-only image.
        pos_bgnd = find(bgnd_mask==1); % positions of bgnd-only intensities in bgnd mask.
        % Get mean (or median) background intensity per pixel in bgnd region, to
        % exclude hot pixels or nearby bright beads, Ibg_avg:
        % Ibg_avg = median(Ibgnd(pos_bgnd)); % use median.
        Ibg_avg = mean(Ibgnd(pos_bgnd)); % use mean.
        % bg_noise_std = std(Ibgnd(pos_bgnd)); % standard deviation of matrix elements in bgnd region.
        % Calculate background-corrected subarray image:
        I2 = I-Ibg_avg;
    else
        % If additional background subtraction not necessary (typically) do:
        I2 = I;
    end
   
    
    %% Obtain particle angle:
          
    % Threshold image:
    % Rescale sub-image to have values between 0 and 1 (required for thresholding):
    I2 = mat2gray(I2); 
    
    % graythresh gives the threshold and im2bw converts image into a
    % thresholded image with only black (zeros) and white (ones) pixels.
    % BgndMask: background mask: 0 at particle regions, 1 in rest.
    % SignalMask: 1 at particle regions, 0 in rest.
    particle_bgnd = im2bw(I2,threshold_factor*graythresh(I2)); % detects dark particle. Logical variable.
    particle = imcomplement(particle_bgnd);

    %the operation 'majority' is used to connect 2 areas of equal size and in
    %close vicinity.
    particle2 = bwmorph(particle,'majority');
    
    conn_regions = bwconncomp(particle2,8); % find connected components in binary image
    % Visual test:
    % To double check the connected regions, we can create a matrix B of the same
    % size as the matrix �particle�, but with all zeros, and then make equal to 1
    % all the pixel positions corresponding to the second connected region found
    % in mov_conn, to plot them:
    % B = zeros([size(particle2,1),size(particle2,2)]);
    % B(conn_regions.PixelIdxList{1}) = 1;
    % imshow(B,[])
    
    % Get properties of connected regions:
    conn_regions_props = regionprops(conn_regions,'Area'); % area, number of pixels in connected region
    partSizeArray = cat(2,conn_regions_props.Area); % array of total number of pixels in each connected area
    
    if isempty(partSizeArray)
        partSizeMax = 1;
    else
        partSizeMax = max(partSizeArray);
    end
    
    % Extract particle skeleton, keep only largest object.
    % BW2 = bwareaopen(BW, P) removes from a binary image all connected components (objects)
    % that have fewer than P pixels, producing another binary image, BW2:
    particle3 = bwareaopen(particle2,partSizeMax);
    % Visual test:
    % figure;
    % subplot(1,2,1),imshow(particle2,[])
    % subplot(1,2,2),imshow(particle3,[])
    
    % particle_out = particle3; % output
    % skeleton = bwmorph(particle3,'skel',Inf); % output
    % Visual test:
    % figure;
    % subplot(1,2,1),imshow(particle3,[])
    % subplot(1,2,2),imshow(skeleton,[])
    
    
    %% Obtain particle orientation:
    
    conn_regions = bwconncomp(particle3,8); % find connected components in final binary image
    conn_regions_props = regionprops(conn_regions,'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid','Area','PixelList');
    % 'Orientation' = angle (degrees) that the major axis of an ellipsoid
    % matching the connected region makes with respect to the horizontal axis.
    % Note that this angle is in degrees ranging from -90 to 90 degrees
    % only!! 
    % See function calAngle.m.
    % 'MajorAxisLength','MinorAxisLenght' are axes lengths from that same ellipsoid.
    % 'Centroid' is position of centre of mass.
    
    % Error control:
    if size(conn_regions_props,1)>1
        conn_regions_props = conn_regions_props(1); % take only first element.
    end
    
    if isempty(conn_regions_props) == false
        
        foundXcentre = conn_regions_props.Centroid(1); % x centre-of-mass position within image subarray.
        foundYcentre = conn_regions_props.Centroid(2); % y centre-of-mass position within image subarray.
        
        allPointsInMask = conn_regions_props.PixelList; % location of all points in the region. List with raws with [x y] for coordinates of each pixel in region.
        
        % Get all points from ellipse contour to plot them and points in
        % major axis (coordinates within image subarray):
        ellipse = getEllipsePts(conn_regions_props); % get all points in the ellipse that matches the connected region 
        xpoints_ellipse = ellipse.contourXpoints;
        ypoints_ellipse = ellipse.contourYpoints;      
        xpoints_majorAxis = ellipse.majorAxisXpoints;
        ypoints_majorAxis = ellipse.majorAxisYpoints;
        xpoints_mask = allPointsInMask(:,1);
        ypoints_mask = allPointsInMask(:,2);
        % [ySkel, xSkel]=find(skeleton==true); % get all skeleton points
        
        % f = getframe; % capture screen shot of figure. F = getframe gets a frame from the current axes.
        % overlap = frame2im(f); % return associated image data
        
        % Get orientation angle:
        % Angle (degrees) that major axis of ellipsoid makes with x axis,
        % ranging from -90 to 90 degrees:
        angleDeg = conn_regions_props.Orientation;
        % The angle will be converted to positive and cyclic later on after
        % tracking, comparing value of angle in previous frames.
                
        %% Output of the function:
        % The output is a structure "particle_result" with the following fields:
        particle_result.estimateXcentre = x_estimate; % x-centre estimate input.
        particle_result.estimateYcentre = y_estimate; % y-centre estimate input.
        particle_result.Xcom = foundXcentre + (round(x_estimate)-d)-1; % x centre-of-mass position within entire frame.
        particle_result.Ycom = foundYcentre + (round(y_estimate)-d)-1; % y centre-of-mass position within entire frame.
        particle_result.AngleDegrees = angleDeg; % angle of major axis with respect to the horizontal, in degrees. 
        particle_result.xpoints_ellipse = xpoints_ellipse + (round(x_estimate)-d)-1; % all points in ellipse, for plotting.
        particle_result.ypoints_ellipse = ypoints_ellipse + (round(y_estimate)-d)-1;
        particle_result.xpoints_majorAxis = xpoints_majorAxis + (round(x_estimate)-d)-1; % for plotting major axis.
        particle_result.ypoints_majorAxis = ypoints_majorAxis + (round(y_estimate)-d)-1;
        particle_result.majorAxisLength = conn_regions_props.MajorAxisLength;
        particle_result.minorAxisLength = conn_regions_props.MinorAxisLength;
        particle_result.ClipFlag = clipping_flag; % 1 if candidate was closer to edge of image than inner_radius.
        particle_result.TooCloseToEdge = tooCloseToEdge; % 1 if bead candidate was closer to edge of image than d_min.
        particle_result.Area = conn_regions_props.Area; % total number of pixels in connected region.
        % -----------------------------------------------------------------------
        % % Auxiliary stuff below:
        % % Plot results:
%         figure;
%         subplot(3,2,1)
%         imshow(I,[])
%         title('Original subarray')
%         subplot(3,2,2)
%         imshow(I2,[])
%         title('Background-subtracted')
%         subplot(3,2,3)
%         imshow(particle3,[])
%         title('Thresholded image')
%         subplot(3,2,4)
%         imshow(frame,[]);
%         hold on;
%         plot(foundXcentre + (round(x_estimate)-d),foundYcentre + (round(y_estimate)-d),'o','Color','g','MarkerSize',5); 
%         title('particle in whole frame')
%         hold off;
%         % plot subarray as well:
%         subplot(3,2,5)
%         imshow(I2,[]);
%         hold on;
%         plot(foundXcentre,foundYcentre,'x','Color','g','MarkerSize',7); 
%         title('fitted ellipse')
%         % Plot overlap of results on frame:
%         plot(xpoints_ellipse,ypoints_ellipse,'y','LineWidth',1);
%         plot(xpoints_majorAxis,ypoints_majorAxis,'y','LineWidth',1);
%         % plot(xSkel, ySkel, '.b', 'MarkerSize',4);
%         hold off;
%         subplot(3,2,6)
%         imshow(I2,[]);
%         hold on;
%         plot(xpoints_mask,ypoints_mask,'x','Color','g','MarkerSize',3)       
%         title('Mask points')
%         hold off;
    else % if no connected regions (no darker particle) found:
        particle_result.estimateXcentre = x_estimate; % x-centre estimate input.
        particle_result.estimateYcentre = y_estimate; % y-centre estimate input.
        particle_result.Xcom = []; % x centre-of-mass position within entire frame.
        particle_result.Ycom = []; % y centre-of-mass position within entire frame.
        particle_result.AngleDegrees = [];
        particle_result.xpoints_ellipse = []; 
        particle_result.ypoints_ellipse = [];
        particle_result.xpoints_majorAxis = [];
        particle_result.ypoints_majorAxis = [];
        particle_result.majorAxisLength = [];
        particle_result.minorAxisLength = [];
        particle_result.ClipFlag = clipping_flag; % 1 if candidate was closer to edge of image than inner_radius.
        particle_result.TooCloseToEdge = tooCloseToEdge; % 1 if bead candidate was closer to edge of image than d_min.
        particle_result.Area = [];
%         % % Plot results:
%         figure;
%         subplot(1,2,1)
%         imshow(frame,[]);
%         % plot subarray as well:
%         subplot(1,2,2)
%         imshow(I2,[]);

        
    end
end