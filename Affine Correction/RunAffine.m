%% Initialization and Parameters
close all,
clear all,
tic,
InitMovieFiles,

nSegments = 4; %Can only be 9, 6 or 4 at this point, but simple to modify code to allow more
nPCs = 200;
PCuse = 1:200;
mu=.15;
nIC = 150;
maxshift = 5; %Maximum LOCAL (within acquisition) shift. Smaller than maximal GLOBAL shift (including btw acq drift)
%% Load Individual Sessions and Calculate Segment Shifts
display('---------------------Tracking Segments-----------------------')
for j=1:num_files
    j,
    fullfilename = correct_filenames{j};
    TrackSegments(fullfilename,movie_file,maxshift,nSegments);
end

MovFile.acqRef = MovFile.acqRef(:,:,1+(1:length(MovFile.acqFrames)));
trackTime = toc,
%% Load and Affine-Correct Each Acquisition
display('---------------------Applying Motion Correction-----------------------')
Affine_Transform_Frames(apply_filenames,movie_file);
transformTime = toc,
%% Calculate 
display('---------------------Computing Principle Components-----------------------')
Movie_PCA(movie_file,nPCs)

pcaTime = toc,

