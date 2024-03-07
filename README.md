# rotTrack - Rotational Tracking of particles

Matlab code for the detection and positional and rotational tracking of (non-spherical) particles in microscopy image sequences (videos).

**Related publication**: 
*Polymeric microellipsoids with programmed magnetic anisotropy for controlled rotation using low (~10 mT) magnetic fields*, A. Bonilla-Brunner, I. Llorente-Garcia, B. Jang, M. Amano Patiño, V. Alimchandani, B. J. Nelson, S. Pane and S. Contera, Applied Materials Today 18, 100511 (2020). https://www.sciencedirect.com/science/article/abs/pii/S2352940719306304?via%3Dihub.

Example:

<p align="center">
  <img src="https://github.com/illg-ucl/rotTrack/blob/master/rotatingParticle_bis_small.png" width=50% height=50%>
</p>

# Copyright and License

Copyright (c) 2017. Isabel Llorente-Garcia, Dept. of Physics and Astronomy, University College London, United Kingdom.

Released and licensed under a BSD 2-Clause License:

https://github.com/illg-ucl/rotTrack/blob/master/LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the BSD 2-Clause License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
BSD 2-Clause License for more details. You should have received 
a copy of the BSD 2-Clause License along with this program. 

Citation: If you use this software for your data analysis please acknowledge 
          it in your publications and cite as follows.
          
              -Citation example 1: 
               RotTrack software. (Version). 2017. Isabel Llorente-Garcia, 
               Dept. of Physics and Astronomy, University College London, United Kingdom.
               https://github.com/illg-ucl/rotTrack. (Download date).
               
              -Citation example 2:
               @Manual{... ,
                title  = {RotTrack software. (Version).},
                author       = {{Isabel Llorente-Garcia}},
                organization = { Dept. of Physics and Astronomy, University College London, United Kingdom.},
                address      = {Gower Place, London, UK.},
                year         = 2017,
                url          = {https://github.com/illg-ucl/rotTrack}
                }

# Guide to run code and basic steps of analysis

* Folder **"scriptsThatHelpRunCode/"** contains various Matlab scripts that help run the code. Script **"scriptToRun.m"** contains a detailed **step-by-step guide** about how to analyse a video sequence doing rotational tracking of particles shown in the video. See also scripts to analyse many video files, e.g. **"scriptToAnalyseManyVideosBis.m"**, and a script to re-analyse the angle obtained (when non-smooth or noisy particle rotation is present), **"scriptToReanalyseAngle.m"**.
* The **basic steps** needed for rotational tracking analysis are as follows (more details in the above mentioned step-by-step guide):
    - Set-up Matlab and correct data folders and current directory;
    - Define "image_label" (string to automatically find the image-sequence file to analyse);
    - Set values of relevant parameters;
    - Find particle trajectories in image sequence and output them to an Excel file in the current directory using function **"FindTrajectsParticles.m"**.
    - Link trajectory segments found into longer trajectories using **"linkTrajSegmentsParticles.m"**.
    - Each track/trajectory is identified by a Trajectory Number. Plot and save a useful .png image of the Trajectory Numbers for the found particless overlaid on top of first frame of the video (using function **"plotParticleTrajNumbers.m"**).
    - Inspect tracks visually (on a video) and manually to validate good tracks by deciding which to accept as good using function **"goThroughParticleTracksVideo.m"**. See options for doing this quickly explained in "scriptToRun.m".
    - Analyse each track separatedly and produce one analysis Excel file and graph per track. This is based on functions **"showParticleTrajAnalysis.m"** and **"showManyParticleTrajAnalysis.m"**. The **analysis includes**: trajectory number, particle orientation angle (degrees) versus time, fit to get angular velocity (degrees/s), particle trajectory on x-y plane, first frame with trajectory overlayed, calculation of particle mean square displacement (MSD) versus delta time with error bars, fit of MSD curve. Calculations assume anti-clockwise rotation.
    - Reanalyse particle angle and rotation speed if needed. If the data has noisy periods and/or the particle rotation is not continuous and smooth, functions **"reanalyseAngle.m"** and **"scriptToReanalyseAngle.m"** can be used to reanalyse a set of Excel files containing track data (frame number, time, angle, etc.), generated in the previous step, and correct the angular velocity and frequency of rotation measured.
The linear sections (angle versus time) with smooth rotation and a large enough slope (above an input threshold value) are automatically detected and fitted to a line (excluding noisy/non-smooth sections).

# Matlab function folders

- **"openImageSequences/"**: functions to open different formats of image sequences (.sif, .dv, .tif, .m4v, .mat or .avi data) and return image data in a useful form. See **"extract_image_sequence_data.m"**.

- **"saveDisplayOthers/"**: contains function **"saveFigurePNG.m"**, which can be used to save a current figure as a .png file.

- **"findParticlePosition/"**: functions to detect and find particle positions on a single image frame. Function **"findCandidateParticlePositions.m"** finds the rough centres (x, y positions) of particles on an image frame, as candidate positions to later refine the particle centre detection. The function recognises dark particles on a lighter background (two methods can be chosen: (1) morphological operations to obtain candidates for particle centres, and (2) Centre-of-Mass method). Function **"removeBgnd.m"** removes or substracts the background from an image of particles via two possible methods (useful for non-uniform backgrounds that would otherwise disturb subsequent steps of the analysis). Function **"excludeRegions.m"** can be used to exclude certain regions of the image. It eliminates all x-y
particle positions found that fall within certain input regions (rectangular boxes). Function **"eliminateCoincidentPositions.m"** eliminates near-coincident particle positions (closer than a certain input distance).

- **"findAngles/"**. The main function is **"findParticleAngle1frame.m"**, which finds the centre-of-mass position and orientation angle of a particle on a single image frame. It uses thresholding, finds a connected region for the particle and obtains orientation as the angle of the major axis of the ellipsoid that best matches the connected region. The output angle is the angle of the major axis with respect to the horizontal, in degrees (ranging from -90 to 90 degrees). The angle is later on converted to positive and cyclic after particle tracking has been done, by comparing to the value of the angle in previous frames. Function **"angleDegToPos.m"** takes an input list of angles between -90 and 90 degreees and converts it into a list of positive cyclic angles (in degrees) that can have values larger than 360 degrees. This function assummes anti-clockwise rotation. Function **"reanalyseAngle.m"** can be used to reanalyse a set of Excel files containing track data (frame number, time, angle, etc.) to correct the angular velocity and frequency of rotation measured. The function automatically detects for a particle the smooth rotation sections which show linear angle versus time (and angular speed - slope - above an input threshold value) and fits them to a line (excluding noisy/non-smooth sections).

- **"TrackParticles/" - Find trajectories and ouput them to Excel file**: functions for tracking particle positions and angles on an image sequence, finding their trajectories (x, y coordinates and time) and orientation angles, and outputting them to an Excel files. 
The two functions used are **"FindTrajectsParticles.m"** and **"linkTrajSegmentsParticles.m"**. The first one finds all particle trajectories and angles in an input image sequence, finds particle centres and joins centres in subsequent frames into trajectory segments. The second one joins segments into longer trajectories for each particle and outputs the resulting trajectory data to an Excel file. It also overlays fully connected trajectories onto the original image sequence in a video. Function **"FindTrajectsParticlesROI.m"** is similar to "FindTrajectsParticles.m" but it uses a region of interest (ROI) to search for particles, instead of excluding regions. Function **"plotParticleTrajNumbers.m"** plots Trajectory Numbers (as they appear on the Excel file generated by "linkTrajSegmentsParticles.m") next to each corresponding particle overlaid on the first frame of the chosen image sequence. This allows easy further analysis, exclusion of particles that are too close, etc.

- **"TrackParticles/" - Inspect and validate tracks**: functions for visually inspecting and validating all particle tracks found in an image sequence. Function **"goThroughParticleTracksVideo.m"** allows visual inspection and validation of the found tracks by showing a video of the found and accepted tracks (overlayed on the original particles) so that you can label each track as good/valid (entering 1) or bad/invalid (entering 0). This process generates structures such as "good_track_nums_label.mat". This visual inspection is useful in case there are objects on the image that you want to exclude, tracks outside acceptable regions, bad tracking or other anomalies. 

- **"TrackParticles/" - Analyse each good track**: function **"showManyParticleTrajAnalysis.m"** (which calls **"showParticleTrajAnalysis.m"**, **"getDisplacement.m"**, etc.) produces one analysis Excel file and one graph per analysed particle track. A folder is generated which contains the Excel and graph files for the analysis of each track in the image sequence. It takes as input a list of "good" track numbers (previously generated). It can show the found particle trajectory and orientation overlayed on the image sequence on a video.

- **"Diffusion-MSDcalculation/"**. Contains main function **"getDisplacement.m"** which analyses trajectory data and calculates the 2D mean square displacements (MSD) as a function of time-interval delta-t for each particle trajectory (track) in a given image sequence (video).

- **"PairWiseDifference/"**: functions for the calculation of pair-wise differences (PwD) for an input vector or series of values. Used for the calculation of the Mean Square Displacement (MSD).
