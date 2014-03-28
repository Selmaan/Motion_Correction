function Affine_Transform_Frames(apply_filenames,movie_file)

%% Load variables / reference points
MovFile = matfile([movie_file '.mat'],'Writable',true);
movie_mask = MovFile.movie_mask;
M=movie_mask(1,3)+1;
N=movie_mask(1,4)+1;

segPos = MovFile.segPos;
yoff = segPos(:,2) + floor(segPos(:,4)/2);
xoff = segPos(:,1) + floor(segPos(:,3)/2);
rpts = [xoff, yoff];
R=imref2d([N,M]);

%% Calculate Session-wise shifts
acqRef = MovFile.acqRef;
acqFrames = MovFile.acqFrames;
[xshiftsAcq,yshiftsAcq]=track_subpixel_wholeframe_motion_varythresh(...
    acqRef,acqRef(:,:,ceil(size(acqRef,3)/2)),10,0.99,100);

%% Loop over acquisitions and apply correction
for j=1:length(acqFrames)
    j,
    %load shifts
    ind = sum(acqFrames(1:j-1)) + (1:acqFrames(j));
    xshift = -(MovFile.cated_xShift(:,ind) + xshiftsAcq(j));
    yshift = -(MovFile.cated_yShift(:,ind) + yshiftsAcq(j));
    
    %load Tiffs and Mask
    t=Tiff(apply_filenames{j});
    origN = t.getTag('ImageLength');
    origM = t.getTag('ImageWidth');
    t.setDirectory(1);
    while ~t.lastDirectory
        t.nextDirectory;
    end
    Z = t.currentDirectory;
    mov = zeros(origN,origM,Z,'single');
    for frame = 1:Z
        t.setDirectory(frame);
        mov(:,:,frame) = t.read;  
        if ~mod(frame, 100)
            fprintf('%1.0f frames loaded.\n', frame);
        end
    end
    mov = mov(movie_mask(2):movie_mask(2)+movie_mask(4),movie_mask(1):movie_mask(1)+movie_mask(3),:);
    t.close();
    
    %Scale for seamless intensity transitions
    meanfirstframes=median(mean(mean(mov(:,:,1:200))));
    if j==1
        meanlastframes = meanfirstframes;
    end
    mov=mov*(meanlastframes/meanfirstframes);
    meanlastframes=median(mean(mean(mov(:,:,end-200:end))));
    
    % Apply affine transform
    parfor frame = 1:Z
        if mod(frame,250)==1
            display(sprintf('frame: %d',frame)),
        end
        xframe = xshift(:,frame) + xoff;
        yframe = yshift(:,frame) + yoff;
        fpts = [xframe, yframe];
        tform=fitgeotrans(fpts,rpts,'affine');
        mov(:,:,frame)=imwarp(mov(:,:,frame),tform,'OutputView',R,'FillValues',nan);   
    end   
    
    %Write to Tiff
    tOut=Tiff(sprintf('%s_Acq%d.tif',movie_file,j),'w');
    tagStruct.RowsPerStrip = 16;
    tagStruct.BitsPerSample = 32;
    tagStruct.SamplesPerPixel = 1;
    tagStruct.ImageLength = N;
    tagStruct.ImageWidth = M;
    tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagStruct.Compression = 1;
    tagStruct.Software = 'MATLAB';
    tagStruct.SampleFormat = 3;
    tOut.setTag(tagStruct);
    
    tOut.write(mov(:,:,1));
    for frame = 2:Z
        tOut.writeDirectory();
        tOut.setTag(tagStruct);
        tOut.write(mov(:,:,frame));
        if ~mod(frame, 100)
            fprintf('%1.0f frames written.\n', frame);
        end
    end
    tOut.close();
    blank_image(:,:,j) = isnan(sum(mov,3));
    corrected_reference(:,:,j) = mean(mov,3);
end

MovFile.blank_image = blank_image;
MovFile.corrected_reference = corrected_reference;

blankM = mean(mean(blank_image,3),1);
blankN = mean(mean(blank_image,3),2);
yMin = find(blankM<(3*min(blankM)),1,'first');
yMax = find(blankM<(3*min(blankM)),1,'last');
xMin = find(blankN<(3*min(blankN)),1,'first');
xMax = find(blankN<(3*min(blankN)),1,'last');
MovFile.blank_frame = [xMin,xMax,yMin,yMax];