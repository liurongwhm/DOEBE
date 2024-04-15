clear;clc;
format long;
format compact;

% Use this code to execute DOEBE
% Input:
%   img3d:      3-D matrix of the image data (row×col×L)
%   img2d:      2-D matrix of the image data (L×N)
%   L:          Number of bands
%   N:          Number of pixels
%   P:          Number of endmembers
%   row:        Number of rows for the 3-D image matrix
%   col:        Number of columns for the 3-D image matrix
%   ref:        Cell array (P×1) restoring the reference endmember spectra and ref{i,1} is the spectra matrix of the ith endmember bundle
%   reference:  Reference endmember spectra matrix which can be expressed as [ref{1,1} ref{2,1}...ref{P,1}]
%   abundance:  Refetence abundance maps (P×N)
%   Subnum:     Number of subimage (Default is 2)
%   PopNum1:    Population size of the multi-objective optimization (Default is 20)
%   PopNum2:    Population size of the single-objective optimization (Default is 20)
%   r:          Threshold of the reconstruction error (Default is 0.006)
%   minpixel:   Minimum number of pixels required for the next iteration (Default is 0.01*N)
%   numRun:     Number of algorithm runs
%Output:
%   em: Spectra of extracted endmembers. em{i,1} is the extracted endmembers of the ith run
%   EBundle: Spectra of extracted endmembers after classified. EBundle{i,1} is the extracted spectra of the ith endmember

%%******************************* Read hyperspectral data *******************************%%
% Input the Urban image
load('F:\Urban\Urban_162.mat')
load('F:\Urban\end6_groundTruth.mat')
reference = M;
abundance = A;
ref = cell(6,1);
for i = 1:6
    ref{i} = reference(:,i);
end
img2d = Y./1000;
[L,N] = size(img2d);
row = 307;
col = 307;
img3d = zeros(row,col,L);
for i = 1:L
    img3d(:,:,i) = reshape(img2d(i,:),row,col);
end
P = 6;   % Number of endmembers

% %Input the Berlin image
% [img2d,img3d] = freadenvi('D:\BerlinUrban2009\BerlinSub1');
% % load('C:\Users\DELL\Desktop\代码\real images\BerlinUrban2009Sub2\abundance.mat')
% load('D:\BerlinUrban2009\library.mat')
% ref = library.refLevel2';
% reference = library.reference;
% img2d = img2d';
% [L,N] = size(img2d);
% row = size(img3d,1);
% col = size(img3d,2);
% % img3d = reshape(img2d',row,col,L);
% P = 6;
% abundance = 0.3*ones(P,col*row);

%%******************************* Set parameters *******************************%%
Subnum = 2;
PopNum1 = 20;
PopNum2 = 20;
r = 0.006;
minPixels = round(0.01*N);
numRun = 1;
%%******************************* Main loop *******************************%%
sad = zeros(numRun,P);
rmse = zeros(numRun,P);
adm = zeros(numRun,P);
ads = zeros(numRun,P);
usedTime = zeros(numRun,1);
num = zeros(numRun,P);
em=cell(numRun,1);

for i = 1:numRun
    tic;
    [indicies,record] = DOEBE(img2d,img3d,Subnum,PopNum1,PopNum2,P,r,minPixels);
    usedTime(i) = toc;
    
    indicies = unique(indicies);
    extractedEndmember = img2d(:,indicies);
    em{i,1} = extractedEndmember;
    
    % Accuracy evaluation
    SAM = SAMpipei(reference,extractedEndmember);
    [EBundle,measure,S] = bundleBasedOnRef(img2d,extractedEndmember,ref,abundance);
    sad(i,:) = measure.sad;
    rmse(i,:) = measure.rmse;
    adm(i,:) = measure.adm;
    ads(i,:) = measure.ads;
    num(i,:) = measure.numEachClass;
    %
    fprintf('the %dth run finished\n',i);
end
meanEachRunSad = mean(sad,2,'omitnan');
meanEachRunRmse = mean(rmse,2,'omitnan');
meanEachRunAdm = mean(adm,2,'omitnan');
meanEachRunAds = mean(ads,2,'omitnan');

resultEachClass = {'meanSad','stdSad','meanRmse','stdRmse','meanAdm','stdAdm','meanAds','stdAds',...,
    'meanTime','stdTime','meanNum','stdNum';mean(sad,'omitnan')',std(sad,'omitnan')',...,
    mean(rmse,'omitnan')',std(rmse,'omitnan')',mean(adm,'omitnan')',std(adm,'omitnan')',...,
    mean(ads,'omitnan')',std(ads,'omitnan')',mean(usedTime,'omitnan')',std(usedTime,'omitnan')',...,
    mean(num,'omitnan')',std(num,'omitnan')'};

resultUnite = {'meanSad','stdSad','meanRmse','stdRmse','meanAdm','stdAdm','meanAds','stdAds',...,
    'meanTime','stdTime','meanNum','stdNum';mean(meanEachRunSad),std(meanEachRunSad),...,
    mean(meanEachRunRmse),std(meanEachRunRmse),mean(meanEachRunAdm),std(meanEachRunAdm),...,
    mean(meanEachRunAds),std(meanEachRunAds),mean(usedTime),std(usedTime),mean(num)',std(num)'};





