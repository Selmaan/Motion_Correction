function TrackSegments(fullfilename,target_file)

MovFile = matfile(target_file,'Writable',true);

info=imfinfo(fullfilename);
numframes=length(info);
M=info(1).Width;
N=info(1).Height;
Z=numframes;

%Load Movie
chone=zeros(N,M,Z,'single');
for frame=1:numframes
    if mod(frame,1000)==1
        frame,
    end
    chone(:,:,frame)=imread(fullfilename,'tiff',frame,'Info',info);
end

%Clip bad image region 
if isempty(MovFile.chone_mask)
    refWin=mean(chone,3);
    imshow(histeq(refWin/max(refWin(:)))),
    h=imrect;
    pause;
    chone_mask = round(getPosition(h));
    MovFile.chone_mask = chone_mask;
else
    chone_mask = MovFile.chone_mask;
end
chone = chone(chone_mask(2):chone_mask(2)+chone_mask(4),chone_mask(1):chone_mask(1)+chone_mask(3),:);
M=chone_mask(3)+1;
N=chone_mask(4)+1;

%Construct Movie Segments
segPos = [];
xind = floor(linspace(1,M/2,3));
yind = floor(linspace(1,N/2,3));
for x=1:length(xind)
    for y=1:length(yind)
        segPos(end+1,:) = [xind(x) yind(y)  floor(M/2) floor(N/2)];
    end
end
nSeg = size(segPos,1);

%First order motion correction
clear xshifts yshifts corThresh,
for Seg = 1:nSeg
    Seg,
    tMov = chone(segPos(Seg,2):segPos(Seg,2)+segPos(Seg,4),segPos(Seg,1):segPos(Seg,1)+segPos(Seg,3),:);
     tBase = prctile(tMov(:),1);
     tTop = prctile(tMov(:),99);
     tMov = (tMov - tBase) / (tTop-tBase);
     tMov(tMov<0) = 0; tMov(tMov>1) = 1;
[xshifts(Seg,:),yshifts(Seg,:)]=track_subpixel_wholeframe_motion_varythresh(...
    tMov,median(tMov,3),5,0.9,100);
end

%Calculate correction for reference image
choneWins = AcquisitionCorrect(chone,mean(xshifts),mean(yshifts));

%Save results to disk
acqFrames = MovFile.acqFrames;
startFrame = sum(acqFrames)+1;
endFrame = startFrame+numframes-1;
MovFile.acqFrames = cat(1,acqFrames,numframes);
MovFile.cated_xShift(1:nSeg,startFrame:endFrame) = xshifts;
MovFile.cated_yShift(1:nSeg,startFrame:endFrame) = yshifts;
MovFile.cated_movie(1:N,1:M,startFrame:endFrame)=chone;
MovFile.acqRef(1:N-10,1:M-10,length(MovFile.acqFrames)+1) = median(choneWins(6:end-5,6:end-5,:),3);
