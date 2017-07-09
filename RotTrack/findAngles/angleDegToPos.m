function angle_pos2 = angleDegToPos(angle_list)
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
% Take an input list of angles (angle_list) between -90 and 90
% degreees and convert it first into a list of angles between 0 and 180
% degrees and then into a list of positive cyclic angles that can have
% values larger than 360 degrees.
% The output (angle_pos) is a list with positive angles in degrees. 


%% Convert list of angles between -90 and 90 degrees into list of angles between 0 and 180 degrees:

% Initialise variable:
angle_pos = zeros(length(angle_list),1); % column vector

for i = 1:length(angle_list)
   if angle_list(i)<0
       angle_pos(i) = 180 + angle_list(i);
   else
       angle_pos(i) = angle_list(i);
   end
end

%% Convert into cyclic positive angles that can take values larger than 360 degrees:
angle_jump = zeros(length(angle_pos)-1,1); % initialise column vector

for i = 1:length(angle_jump)
    angle_jump(i) = angle_pos(i+1)-angle_pos(i);
end

pos_jumps = 1+find(angle_jump < -90); % position of jumps larger than -90 degrees in angle_pos vector.


% Initialise output variable:
angle_pos2 = angle_pos; % column vector

% We leave angle as it is for first angles between 0 and 180 degrees and
% correct the rest to get positive cyclic angles.
for i = 1:length(pos_jumps) % loop through jumps
    % Add 180 degrees to all angles after each jump:
    for j = pos_jumps(i):length(angle_pos2)
       angle_pos2(j) = 180 + angle_pos2(j);  
    end
end

