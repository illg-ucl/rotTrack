function FileWrite(fileinfo,framenum,newframe)

%FileWrite(fileinfo,framenum,newframe)
%
%FileFrame.m is a front end which unifies the commands necessary to write
%frames to common file formats. It overwrites an existing frame.
%This should only be performed on files which have first been coppied.
%   
%INCLUDE:   (none)
%
%INPUTS:    FILEINFO:   A structure containing pertinent information about
%                       the file, generated by FileDetails.m  In order for
%                       this code to work, FILEINFO must have a field .copy
%           FRAMENUM:   The index of the frame to replace
%           NEWFRAME:   The data for a new frame, where (1,1) corisponds to
%                       a pixel in the upper left corner
%
%OUTPUTS:   (none)
%               
%Copyright Scott Parker 10/2009 U. Illinois Urbana-Champaign
%Last modified by Scott Parker on 10/08/2009

% Check if the copy field exists.  If not, display a warning and terminate
% the program
if ~strcmp(fileinfo.copy,'Duplicate File')
  error('This program should only be run on copies, never on the original data')
end


%Call the appropriate function dependent upon file type, by comparing the
%file extension in a non-case sensitive fashion
if strcmpi(fileinfo.format,'.sif')
    %Andor Technology MultiChannel File Format
    SifWrite(fileinfo,framenum,newframe);
elseif strcmpi(fileinfo.format,'.spe')
    %SPE File format (WinSpec)
    SpeWrite(fileinfo,framenum,newframe);
elseif strcmpi(fileinfo.format,'.cine')
    %CINE Vision Research File Format
    WriteCine(fileinfo,framenum,newframe);
end