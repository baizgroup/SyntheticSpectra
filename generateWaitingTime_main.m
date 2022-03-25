%Script to generate synthetic 2D IR spectra using a Redfield response function.
%Load laser shots and include experimental noise
%Baiz Group, February 2022
addpath(genpath('functions'))

%path to blank laser shots
pathName.dataFolder = 'blank_shots_dataset';
pathName.outputFolder = 'onepeak_CLS_tau1500_broad';

%% Input section
in.manual = false; %plot spectra
% set parameters for model inputs
param = initializeModelInputs();
% Set frequency axis
in.freqAx = (1600:0.5882:1750);
% Number of spectra to generate
in.systemsToGenerate = 200; %how many spectra to generate
% List of waiting times in femtosecond to generate
in.listDelays = [150:14.25:3000]; %repmat(150, [in.systemsToGenerate 1]);
in.crossPeaks = 0; %0 for no cross peaks, and 1 for cross peaks. 2 for dynamic cross peaks

% variables for loading blank laser shots
in.numDataSets = 100; %number of files to load
in.numShots = 5000; %number of shots in each file

in.SNR = Inf;% signal (max) to noise (stdev) ratio
in.numNoisetoInclude = 5; %number of blank shot files to average over to generate noise trajectory

%% loading and preprocessing

%load experimental blank shots - CB 2/18/2022
blankShots = dataLoaderAllShots(pathName.dataFolder,in.numDataSets);
data.diffShots = blankShots.MCT(:,1:2:end,:) - blankShots.MCT(:,2:2:end,:);
clear blankShots;

data.numShotsPerAve = size(param.aux.taxis,2);
data.NumRepeatsinShots = floor(size(data.diffShots,2)./data.numShotsPerAve); %how many repeats in X shots

% split noise by average
for NumAve = 1:size(data.diffShots,3)
    CurrentShots = data.diffShots(:,:,NumAve);
    for n=1:data.numShotsPerAve
        data.OutputNoiseShots(:,n,NumAve) = mean(CurrentShots(:,1:data.numShotsPerAve:end),2);
        CurrentShots = circshift(CurrentShots,[1 -1]);
    end
    clear CurrentShots
end

    mkdir(['output/' pathName.outputFolder '_Inf'])
    mkdir(['output/' pathName.outputFolder '_20'])
    mkdir(['output/' pathName.outputFolder '_10'])
    mkdir(['output/' pathName.outputFolder '_5'])
    mkdir(['output/' pathName.outputFolder '_2'])

%%
for n = 1:in.systemsToGenerate
    % run parameter file to generate random inputs
    param  = initializeModelInputs();
    model.randomizedParams(n,:) = param.packed.xo; %replicate initial parameters (for t2 series);
    %model.randomizedParams(n,:) = param.packed.random; %random parameters
    %within bounds
end

%%


for n = 1:in.systemsToGenerate

    %generate a random noise shot series from the experimental data
    tempRandVec = randperm(size(data.OutputNoiseShots,3));
    AveOutputNoiseShots = mean(data.OutputNoiseShots(:,:,tempRandVec(1:in.numNoisetoInclude)),3);
    clear tempRandVec

    [model.interpSpecInf(:,:,n)]=dynamicModel(model.randomizedParams(n,:),param.aux,in.listDelays(n),in.freqAx,AveOutputNoiseShots,Inf);
    [model.interpSpec20(:,:,n)]=dynamicModel(model.randomizedParams(n,:),param.aux,in.listDelays(n),in.freqAx,AveOutputNoiseShots,20);

    [model.interpSpec10(:,:,n)]=dynamicModel(model.randomizedParams(n,:),param.aux,in.listDelays(n),in.freqAx,AveOutputNoiseShots,10);
    [model.interpSpec5(:,:,n)]=dynamicModel(model.randomizedParams(n,:),param.aux,in.listDelays(n),in.freqAx,AveOutputNoiseShots,5);
    [model.interpSpec2(:,:,n)]=dynamicModel(model.randomizedParams(n,:),param.aux,in.listDelays(n),in.freqAx,AveOutputNoiseShots,2);

    if in.manual
        figure(11); clf;
        spectrumMin = min(min(model.interpSpecInf(:,:,n)));
        spectrumMax = max(max(model.interpSpecInf(:,:,n)));
        contourf(in.freqAx,in.freqAx, model.interpSpecInf(:,:,n),...
            [spectrumMin:spectrumMax/10:spectrumMax]);
        caxis([spectrumMin,spectrumMax]); colormap(cmap2d(20));
        axis square
        line([in.freqAx(1) in.freqAx(end)],[in.freqAx(1)...
            in.freqAx(end)],'color',[0 0 0]);
        xlabel('\omega_1 (cm^{-1})')
        ylabel('\omega_3 (cm^{-1})')
        title(['t2 = ',num2str(in.listDelays(n)),'fs'])
        drawnow
        pause(0.2)
    end

    if ~in.manual
    g=imagesc(flipud(model.interpSpecInf(:,:,n)));
    axis square;
    colormap(gray)
    temp1= g.CData-min(min(g.CData));
    temp2=temp1./max(max(temp1));
    imwrite(temp2, ['output/' pathName.outputFolder '_Inf/spec_1_'  num2str(n,'%0.5d') '.png'])
    imwrite(temp2, ['output/' pathName.outputFolder '_Inf/spec_2_'  num2str(n,'%0.5d') '.png'])
    imwrite(temp2, ['output/' pathName.outputFolder '_Inf/spec_3_'  num2str(n,'%0.5d') '.png'])
    imwrite(temp2, ['output/' pathName.outputFolder '_Inf/spec_4_'  num2str(n,'%0.5d') '.png'])

    g=imagesc(flipud(model.interpSpec20(:,:,n)));
    axis square;
    colormap(gray)
    temp1= g.CData-min(min(g.CData));
    temp2=temp1./max(max(temp1));
    imwrite(temp2, ['output/' pathName.outputFolder '_20/spec_1_'  num2str(n,'%0.5d') '.png'])

    g=imagesc(flipud(model.interpSpec10(:,:,n)));
    axis square;
    colormap(gray)
    temp1= g.CData-min(min(g.CData));
    temp2=temp1./max(max(temp1));
    imwrite(temp2, ['output/' pathName.outputFolder '_10/spec_2_'  num2str(n,'%0.5d') '.png'])

    g=imagesc(flipud(model.interpSpec5(:,:,n)));
    axis square;
    colormap(gray)
    temp1= g.CData-min(min(g.CData));
    temp2=temp1./max(max(temp1));
    imwrite(temp2, ['output/' pathName.outputFolder '_5/spec_3_'  num2str(n,'%0.5d') '.png'])

    g=imagesc(flipud(model.interpSpec2(:,:,n)));
    axis square;
    colormap(gray)
    temp1= g.CData-min(min(g.CData));
    temp2=temp1./max(max(temp1));
    imwrite(temp2, ['output/' pathName.outputFolder '_2/spec_4_'  num2str(n,'%0.5d') '.png'])
    end
    
    disp(num2str(n));
end

save(['output/' pathName.outputFolder '.mat'],'param','model','in')