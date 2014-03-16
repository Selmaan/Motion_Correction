%Load individual acquisitions, calculate pt shifts for affine and reference
%image and save unmodified to a concatenated .mat file, clear from memory.
%Then reload individual acquitisions from .mat and apply frame-wise motion
%correction affine transform, writing concatenated .mat file. Then load
%pixel subsets (for all frames) and calculate frame covariances (on df/F normalized pixels),
%average over all pixel subsets to get covariance matrix. Eigs on covariance
%for coveval/mixed sig (normalized) and projections of data for mixed filters
%(needs to be normalized). Convert both to doubles and pass to ICA
close all,
clear all,
%% Create File Names
num_files = 6;
nameOffset = 0;
mouse_name = 'LD085';
session_name = '131125';
view_name = 'view1';
slice_name = 'slice2';
tiffpath = cd;
target_file = sprintf('%s_%s_%s_%s',mouse_name,session_name,view_name,slice_name);

for j=1:num_files
    filenames{j} = sprintf('%s_%s_%.3d_%s_%s.tif.tif',mouse_name,session_name,j+nameOffset,view_name,slice_name);
end

%% Load individual acquisition, Track Segments 
display('---------------------Tracking Segments-----------------------')
MovFile = matfile(target_file,'Writable',true);
MovFile.chone_mask = [];
MovFile.acqFrames=[];
MovFile.cated_xShift = [];
MovFile.cated_yShift = [];
MovFile.cated_movie=zeros(0,0,0,'single');
MovFile.acqRef = zeros(0,0,0,'single');

for j=1:num_files
    j,
    fullfilename = filenames{j};
    TrackSegments(fullfilename,target_file);
end

MovFile.acqRef = MovFile.acqRef(:,:,1+(1:length(MovFile.acqFrames)));
%% Load and Affine-Correct Each Acquisition
display('---------------------Applying Motion Correction-----------------------')
corrected_file = sprintf('%s_Corrected',target_file);
Affine_Transform_Frames(target_file,corrected_file)

%% Calculate acquisition-wise covariance matrix
nPCs = 1e3;
nSegs = 10;
display('---------------------Computing Principle Components-----------------------')

Movie_PCA(corrected_file,nPCs,nSegs)

%% Run ICA algorithm
CorrFile = matfile(corrected_file,'Writable',true);

PCuse = 1:200;
mu=.2;
nIC = 150;
ica_A_guess = [];
termtol = 1e-6;
maxrounds = 1e3;
smwidth = 3;
thresh = 2;
arealims = [20 400];
plotting = 1;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(...
     CorrFile.mixedsig,CorrFile.mixedfilters, CorrFile.CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);
 [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation...
    (ica_filters, smwidth, thresh, arealims, plotting);
for i=1:size(ica_segments,1)
    normSeg(:,:,i)=100*ica_segments(i,:,:)/norm(reshape(ica_segments(i,:,:),1,[]));
end