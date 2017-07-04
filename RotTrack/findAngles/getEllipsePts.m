function ellipse = getEllipsePts(s)
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
% Extract all points (x, y) of ellipse given by input s. 
% INPUT: s is the result of using function regionprops on the found connected regions of an image frame.
% Example of how to use this function:
% conn_regions = bwconncomp(particle3,8); % find connected components in final binary image
% conn_regions_props = regionprops(conn_regions,'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');
% [xElli,yElli] = getEllipsePts(conn_regions_props); 
% Extract elipse coordinates from regions given by 'regionprops'.

phi = linspace(0,2*pi,50)'; % column vector
cosphi = cos(phi);
sinphi = sin(phi);

xcentre = s.Centroid(1);
ycentre = s.Centroid(2);

a = s.MajorAxisLength/2;
b = s.MinorAxisLength/2;

theta = pi*s.Orientation/180;
% Rotation matrix:
R = [cos(theta)   sin(theta)
     -sin(theta)   cos(theta)];

xy = [a*cosphi b*sinphi];
xy_rot = R*xy'; % apply rotation.

% points of the ellipse:
x = xy_rot(1,:) + xcentre;
y = xy_rot(2,:) + ycentre;

% Get major axis:
% x and y positions of three extreme points in major axis, with respect
% to centre of ellipse:
majorAxis_xy = [a 0; 0 0; -a 0]';
majorAxis_xy_rot = R*majorAxis_xy; % rotated major axis.
majorAxisXpoints = majorAxis_xy_rot(1,:) + xcentre;
majorAxisYpoints = majorAxis_xy_rot(2,:) + ycentre;

% Output:
ellipse.contourXpoints = x; % row vector with x positions of ellipse contour;
ellipse.contourYpoints = y; % row vector with y positions of ellipse contour;
ellipse.majorAxisXpoints = majorAxisXpoints; % row vector with x of 3 extreme points of major axis with respect to entire frame.
ellipse.majorAxisYpoints = majorAxisYpoints; % row vector with y of 3 extreme points of major axis with respect to entire frame.
