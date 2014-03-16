function Affine_Transform_Frames(target_file,corrected_file)

MovFile = matfile(target_file,'Writable',false);
CorrFile = matfile(corrected_file,'Writable',true);
M=MovFile.chone_mask(1,3)+1;
N=MovFile.chone_mask(1,4)+1;

segPos = [];
xind = floor(linspace(1,M/2,3));
yind = floor(linspace(1,N/2,3));
for x=1:length(xind)
    for y=1:length(yind)
        segPos(end+1,:) = [xind(x) yind(y)  floor(M/2) floor(N/2)];
    end
end
yoff = segPos(:,2) + floor(segPos(:,4)/2);
xoff = segPos(:,1) + floor(segPos(:,3)/2);
rpts = [xoff, yoff];
R=imref2d([N,M]);

acqRef = MovFile.acqRef;
acqFrames = MovFile.acqFrames;
[xshiftsAcq,yshiftsAcq]=track_subpixel_wholeframe_motion_varythresh(...
    acqRef,acqRef(:,:,ceil(size(acqRef,3)/2)),10,0.99,100);


for j=1:length(acqFrames)
    j,
    ind = sum(acqFrames(1:j-1)) + (1:acqFrames(j));
    xshift = -(MovFile.cated_xShift(:,ind) + xshiftsAcq(j));
    yshift = -(MovFile.cated_yShift(:,ind) + yshiftsAcq(j));
    pxNorm = 100/median(median(acqRef(:,:,j)));
    tMov = pxNorm * MovFile.cated_movie(:,:,ind);
    for frame = 1:size(tMov,3)
        if mod(frame,250)==1
            frame,
        end
        xframe = xshift(:,frame) + xoff;
        yframe = yshift(:,frame) + yoff;
        fpts = [xframe, yframe];
        tform=fitgeotrans(fpts,rpts,'affine');
        tMov(:,:,frame)=imwarp(tMov(:,:,frame),tform,'OutputView',R);   
    end   
    CorrFile.Movie(1:N,1:M,ind) = tMov;
    cropMov(:,:,j) = sum(tMov==0,3);
end

cropMov = sum(cropMov,3);
CorrFile.xmin=find(median(cropMov,1)==0,1,'first');
CorrFile.xmax=find(median(cropMov,1)==0,1,'last');
CorrFile.ymin=find(median(cropMov,2)==0,1,'first');
CorrFile.ymax=find(median(cropMov,2)==0,1,'last');
CorrFile.N=length(CorrFile.ymin:CorrFile.ymax);
CorrFile.M=length(CorrFile.xmin:CorrFile.xmax);
CorrFile.Z=max(ind);
CorrFile.Movie = CorrFile.Movie(CorrFile.ymin:CorrFile.ymax,CorrFile.xmin:CorrFile.xmax,:);
