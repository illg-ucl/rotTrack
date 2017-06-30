function scriptToRun

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




%% To carry out tests of the methods on a single frame:

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


% % Subtract background, method 1 is faster:
% frameNoBgnd1 = removeBgnd(frame1,new_x1,new_y1,50,60,1);
% 
% % Test finding bead centre on single frame:
% s1 = findBeadCentre1frame(frameNoBgnd1,new_x1(1),new_y1(1),50,60);