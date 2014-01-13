%% Load and Motion Correct

clear all;
close all;

%get all files names
num_files = 11;
nameOffset = 0;
cated_tiff_filename = 'a16';
indv_tiff_filename = 'FAST';
tiffpath = cd;
smallWin = 5;

for j=1:num_files
%     [tifffilename,tiffpath]=uigetfile('*.tif','pick your tiff file');
%     eval(['fullfilename' num2str(j) ' = [tiffpath tifffilename]']);
filenames{j} = sprintf('%s%.3d.tif',indv_tiff_filename,j+nameOffset);
end


%open files and concatenate; scale for same instensities
cated_movie=[];
cated_corThresh=[];
cated_xShift = [];
cated_yShift = [];
acqRef = [];
for j=1:num_files

    %eval(['fullfilename=fullfilename' num2str(j)]);
    fullfilename = filenames{j};
    info=imfinfo(fullfilename);
    numframes(j)=length(info);
    M=info(1).Width;
    N=info(1).Height;

    chone=uint16(zeros(N,M,numframes(j)));
    for i=1:numframes(j)
        if mod(i,1000)==1
            j
            i
        end
        chone(:,:,i)=imread(fullfilename,'tiff',i,'Info',info);
    end

    %scale movie for seamless intensities
    if j==1
        meanlastframes=median(mean(mean(chone(:,:,1:1e3))));
    end

    meanfirstframes=median(mean(mean(chone(:,:,1:1e3))));
    chone=chone*(meanlastframes/meanfirstframes);
    meanlastframes=median(mean(mean(chone(:,:,end-1e3:end))));
    
    chone=single(chone);
    [xshifts,yshifts,corThresh]=subpixel_minibatch(chone,smallWin,3,0.9,100);
    choneWins = AcquisitionCorrect(chone,xshifts,yshifts);
    choneWins = choneWins(3:end-2,3:end-2,:); 
    [xFshifts,yFshifts,FcorThresh] = track_subpixel_wholeframe_motion_varythresh(choneWins,median(choneWins,3),3,0.925,100);
    choneWins = AcquisitionCorrect(choneWins,xFshifts,yFshifts);
       
    cated_movie=cat(3,cated_movie,chone);
    cated_corThresh=cat(2,cated_corThresh,FcorThresh);
    cated_xShift = cat(2,cated_xShift, xshifts+xFshifts);
    cated_yShift = cat(2,cated_yShift, yshifts+yFshifts);
    acqRef = cat(3,acqRef,median(choneWins(3:end-2,3:end-2,:),3));
end

clear info chone choneWins

[xshiftsAcq,yshiftsAcq]=track_subpixel_wholeframe_motion_varythresh(acqRef,acqRef(:,:,ceil(num_files/2)),10,0.99,100);
for j=1:num_files
    ind = sum(numframes(1:j-1)) + (1:numframes(j));
    xshiftAcq(ind) = xshiftsAcq(j);
    yshiftAcq(ind) = yshiftsAcq(j);
end
% 
% xshiftAcq=reshape(repmat(xshiftsAcq,numframes,1),1,[]);
% yshiftAcq=reshape(repmat(yshiftsAcq,numframes,1),1,[]);

xshift = cated_xShift + xshiftAcq;
yshift = cated_yShift + yshiftAcq;
cated_movie = playback_wholeframe_subpix(cated_movie,xshift,yshift);

% save(cated_tiff_filename,'cated_movie','cated_corThresh','-v7.3')
for j=1:num_files
    j,
    ind = sum(numframes(1:j-1)) + (1:numframes(j));
    AcqMov = cated_movie(:,:,ind);
    AcqThreshs = cated_corThresh(ind);
    save(sprintf('%s_Acq%d',cated_tiff_filename,j),'AcqMov','AcqThreshs','-v7.3')
end
clear AcqMov

%% Calculate piecewise covariance, principle components, and feed to ICA algorithm

M=size(cated_movie,1);
N=size(cated_movie,2);
Z=size(cated_movie,3);
nPCs = 500;
cated_movie=reshape(cated_movie,M*N,size(cated_movie,3));

AcqCov = zeros(M*N,M*N);
for j=1:num_files
    j,
ind = sum(numframes(1:j-1)) + (1:numframes(j));
AcqCov = AcqCov + cov(double(cated_movie(:,ind)'));
end
AcqCov = AcqCov / num_files;

[V,D] = eig(AcqCov);
D=diag(D);
V=fliplr(V);
D=flipud(D);
clear AcqCov
covtrace = sum(D);
CovEvals = D(1:nPCs);
clear D;
mixedsig = V(:,1:nPCs)' * double(cated_movie);
mixedfilters = reshape(V(:,1:nPCs),M,N,nPCs);
save(sprintf('%s_PCs',cated_tiff_filename),'mixedsig','mixedfilters','covtrace','CovEvals')


PCuse = 1:300;
mu=.2;
nIC = ceil(length(PCuse)/2);
ica_A_guess = [];
termtol = 1e-6;
maxrounds = 1e3;
smwidth = 3;
thresh = 2;
arealims = [20 400];
plotting = 1;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(...
     mixedsig,mixedfilters, CovEvals, PCuse, mu, nIC, ica_A_guess, termtol, maxrounds);
[ica_segments, segmentlabel, segcentroid] = CellsortSegmentation...
    (ica_filters, smwidth, thresh, arealims, plotting);
for i=1:size(ica_segments,1)
    normSeg(:,:,i)=100*ica_segments(i,:,:)/norm(reshape(ica_segments(i,:,:),1,[]));
end
normSeg = reshape(normSeg,M*N,[]);
SegTraces = normSeg' * cated_movie;
normSeg = reshape(normSeg,M,N,[]);

for i=1:size(ica_segments,1)
    segSize(i) = squeeze(sum(sum(ica_segments(i,:,:)>0)));
    segSkew(i) = skewness(reshape(ica_segments(i,:,:),1,[]));
end
segSTD = sqrt(std(SegTraces,[],2)./mean(SegTraces,2));
segSize = zscore(segSize);
segSkew = zscore(segSkew);

goodSeg = find(segSize-segSkew > 0);
gTrace = SegTraces(goodSeg,:);
save(sprintf('%s_ICs',cated_tiff_filename),...
    'SegTraces','normSeg','ica_sig','ica_filters','ica_A','ica_segments','segmentlabel','segcentroid',...
    'segSize','segSkew','segSTD','goodSeg','gTrace')
figure,scatter(segSize,segSkew,50,segSTD,'filled')
