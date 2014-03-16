function Movie_PCA(corrected_file,nPCs,nSegs)

CorrFile = matfile(corrected_file,'Writable',true);
M=CorrFile.M;
N=CorrFile.N;
Z=CorrFile.Z;
acqCov = zeros(Z,Z,'double');

segs = ceil(linspace(0,M,nSegs+1));

for seg=1:nSegs
    seg,
    indM = 1+segs(seg) : segs(seg+1);
    indN = 1:N;
    indZ = 1:Z;
    tMov = reshape(CorrFile.Movie(indN,indM,indZ),[],Z);
    tMov = tMov ./ repmat(mean(tMov,2),1,Z) - 1;
    acqCov = acqCov + cov(double(tMov));
end
acqCov = acqCov/nSegs;
clear tMov,

display('------------------Computing Eigenvalues------------------')
[V,D] = eig(acqCov);
D=flipud(diag(D));
V=fliplr(V);
covtrace = sum(D);
CovEvals = D(1:nPCs);
mixedsig = V(:,1:nPCs)';
clear V D acqCov,

display('------------------Computing Filters------------------')

for seg=1:nSegs
    seg,
    indM = 1+segs(seg) : segs(seg+1);
    indN = 1:N;
    indZ = 1:Z;
    tMov = double(reshape(CorrFile.Movie(indN,indM,indZ),[],Z));
    tMov = tMov ./ repmat(mean(tMov,2),1,Z) - 1;
    mixedfilters(indN,indM,1:nPCs) = reshape(tMov*mixedsig',N,[],nPCs);
end

mixedfilters=reshape(mixedfilters,N*M,nPCs);
filterMean = mean(mixedfilters,1);
filterStd = std(mixedfilters,[],1);
mixedfilters = (mixedfilters-repmat(filterMean,N*M,1))./repmat(filterStd,N*M,1);
mixedfilters = reshape(mixedfilters,N,M,nPCs);

CorrFile.mixedsig = mixedsig;
CorrFile.mixedfilters = mixedfilters;
CorrFile.CovEvals = CovEvals;
CorrFile.covtrace = covtrace;
