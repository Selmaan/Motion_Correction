function Movie_PCA(movie_file,nPCs)

MovFile = matfile(movie_file,'Writable',true);
M=MovFile.movie_mask(1,3)+1;
N=MovFile.movie_mask(1,4)+1;
acqFrames = MovFile.acqFrames;
Z = sum(MovFile.acqFrames);
nAcqs = length(MovFile.acqFrames);
acqCov = zeros(Z,Z,'single');
t=Tiff(sprintf('%s_Acq%d.tif',movie_file,1),'r');
nStrips = t.numberOfStrips;
nRows = t.getTag('RowsPerStrip');
blank_frame = MovFile.blank_frame;

parfor strip = 1:nStrips
    display(sprintf('Strip Number: %d',strip)),
    if strip < nStrips
        tMov = zeros(nRows,M,Z,'single');
    else
        tMov = zeros(mod(N,nRows),M,Z,'single');
    end
    
    for j=1:nAcqs
        t=Tiff(sprintf('%s_Acq%d.tif',movie_file,j),'r');
        for frame = 1:acqFrames(j)
            t.setDirectory(frame);
            tMov(:,:,frame+sum(acqFrames(1:j-1))) = t.readEncodedStrip(strip);  
        end
        t.close();
    end  
    % Blank nan pixels and dF/F normalization
    tMov = reshape(tMov,[],Z);
    pxMean = mean(tMov,2);
    tMov = tMov ./ repmat(pxMean,1,Z) - 1;
    tMov(isnan(pxMean),:) = [];
    acqCov = acqCov + cov(tMov);          
end

acqCov = double(acqCov);
acqCov=acqCov / nStrips;

display('------------------Computing Eigenvalues------------------')
[V,D] = eig(acqCov);
D=flipud(diag(D));
V=fliplr(V);
covtrace = sum(D);
CovEvals = D(1:nPCs);
mixedsig = V(:,1:nPCs)';
clear V D acqCov,

display('------------------Computing Filters------------------')
mixedfilters = zeros(N,M,nPCs);
for strip = 1:nStrips
    display(sprintf('Strip Number: %d',strip)),
    if strip < nStrips
        tMov = zeros(nRows,M,Z,'single');
    else
        tMov = zeros(mod(N,nRows),M,Z,'single');
    end
    
    for j=1:nAcqs
        t=Tiff(sprintf('%s_Acq%d.tif',movie_file,j),'r');
        for frame = 1:acqFrames(j)
            t.setDirectory(frame);
            tMov(:,:,frame+sum(acqFrames(1:j-1))) = t.readEncodedStrip(strip);  
        end
        t.close();
    end
    indN = (strip-1)*nRows+1 : (strip-1)*nRows+size(tMov,1);
    tMov = reshape(tMov,[],Z);
    pxMean = mean(tMov,2);
    tMov = tMov ./ repmat(pxMean,1,Z) - 1;
    mixedfilters(indN,:,:) = reshape(tMov*mixedsig',length(indN),M,nPCs);
end
    
    
mixedfilters=reshape(mixedfilters,N*M,nPCs);
filterMean = nanmean(mixedfilters,1);
filterStd = nanstd(mixedfilters,[],1);
mixedfilters = (mixedfilters-repmat(filterMean,N*M,1))./repmat(filterStd,N*M,1);
mixedfilters = reshape(mixedfilters,N,M,nPCs);

MovFile.mixedfilters = mixedfilters;
MovFile.mixedsig = mixedsig;
MovFile.CovEvals = CovEvals;
MovFile.covtrace = covtrace;