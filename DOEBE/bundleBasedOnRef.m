function [EBundle,measure,S] = bundleBasedOnRef(img2d,extractedEndmember,ref,abundance)
% endmember bundle extraction evaluation method
% the cluster of the bundle is based on the ref(reference), ref is a cell
% array that each element is one class of endmembers
% the EBundle is a cell array, each dimension has the same class of
% endmembers with the corresonding dimension of ref
% the result is a cell array that contains all evaluation results
P = size(ref,1);
[L,numE] = size(extractedEndmember);
N = size(img2d,2);
EBundle = cell(P,1);

bundleIndex = zeros(numE,1);
sad = zeros(P,1);
rmse = 100*ones(P,1);
adm = zeros(P,1);
ads = zeros(P,1);
numEachClass = zeros(P,1);

% find which class each endmember belongs to based on the ref
sadEachEndmember = zeros(numE,1);
for i = 1:numE
    tempMinSadEachClass = zeros(P,1);
    for j = 1:P
        tempRef = ref{j,1};
        numRefEachClass = size(tempRef,2);
        tempSadOneClass = zeros(numRefEachClass,1);
        for k = 1:numRefEachClass
            tempSadOneClass(k)= real(SAM(extractedEndmember(:,i),tempRef(:,k)));
        end
        tempMinSadEachClass(j) = min(tempSadOneClass);
    end
    [sadEachEndmember(i),bundleIndex(i)] = min(tempMinSadEachClass);
end
%calculate the mean sad value for each class
for i = 1:P
    idx = bundleIndex == i;
    sad(i) = mean(sadEachEndmember(idx),'omitnan');
end

%calculate the number of endmembers for each class
for i = 1:P
    idx = find(bundleIndex == i);
    numEachClass(i) = length(idx);
    EBundle{i} = extractedEndmember(:,idx);
end

%calculate the mean difference of the mean spectral value from each class
%of the estimated bundles and ref
endmemberForCalAbundance = zeros(L,P);
for i = 1:P
    meanEBundle = mean(EBundle{i},2,'omitnan');
    meanRef = mean(ref{i},2);
    adm(i) = sum(abs(meanEBundle-meanRef))/L;
    
    stdEBundle = std(EBundle{i},1,2,'omitnan'); 
    stdRef = std(ref{i},1,2);
    ads(i) = sum(abs(stdEBundle-stdRef))/L;
    
    endmemberForCalAbundance(:,i) = meanEBundle;
end

% use the collaborative sparse unmixing method to calculate the abundance
SBundle = sunsal([extractedEndmember,0.01*ones(L,1)],img2d,'POSITIVITY','yes','VERBOSE','yes','ADDONE','no', ...
    'LAMBDA', 0.001,'AL_ITERS',2000);
SBundle = SBundle(1:numE,:)./repmat(sum(SBundle(1:numE,:))+eps,numE,1);

% SBundle = sunsal(extractedEndmember,img2d,'POSITIVITY','yes','VERBOSE','yes','ADDONE','no', ...
%     'LAMBDA', 0.00001,'AL_ITERS',2000);
S = zeros(P,N);
for i = 1:P
    idx = bundleIndex == i;
    S(i,:) = sum(SBundle(idx,:),1);
end

% use the Fcls to calculate the abundance
% idx = find(isnan(endmemberForCalAbundance(1,:)));
% endmemberForCalAbundance(:,idx) = [];
% S = hyperFcls(img2d,[endmemberForCalAbundance,0.01*ones(L,1)]);
% S = S(1:P-length(idx),:)./repmat(sum(S(1:P-length(idx),:))+eps,P-length(idx),1);

%calculate the rmse
% abundance(idx,:) = [];
for i=1:P
    rmse(i) = sqrt(sum((abundance(i,:)- S(i,:)).^2)/N);
end

measure.bundleIndex = bundleIndex;
measure.sad = sad;
measure.rmse = rmse;
measure.adm = adm;
measure.ads = ads;
measure.numEachClass = numEachClass; 
end

