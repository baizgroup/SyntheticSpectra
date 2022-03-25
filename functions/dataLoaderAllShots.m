function data = dataLoaderAllShots(dataFolder, numFilesToLoad)

% This function loads all shots data at a single T2 delay
% inputs: input, exp file
% outputs: data.MCT(px,shots, delayCount) and data.DAQ{} for mask averaging

%      loaded data will look like:
%            NumMasks2: 0
%     allShotsPath: 'C:\Users\baizgroup\Box Sync\BaizGroup_Shared\Data\2DIR\reena2\'
%          counter: 0
%         rawArray: [10000×384 double]
%       rawDAQData: [9936×2 double]
%          savenow: 1
% all shots data in .mat format is saved as (#shots x pixels), like 10000 shots x 384 pixels
dataFiles = dir([dataFolder '/*.mat']);
dataUrl = 'https://baizgroup.org/data/blank_shots_dataset/blank_shots_dataset_2020_07_31.zip';

%download full dataset (102 files, 1.2GB)
if isempty(dataFiles)
    if ~exist(dataFolder, 'dir')
    mkdir(dataFolder)
    end
    
    getZip(dataFolder, dataUrl)
    %reload files after download/unzip
    dataFiles = dir([dataFolder '/*.mat']);
end

% Load 1 data file to get size
datafile1 = load([dataFiles(1).folder '/' dataFiles(1).name]);
[dim1, dim2] = size(datafile1.rawArray); % dim1 = numShots, dim2 = numPixels

% How many times data was collected at that delay
totRuns = min([length(dataFiles) numFilesToLoad]);

disp(['Loaded data file 1/' num2str(totRuns)]);

% Preallocate data array for SPEED
data.MCT = NaN(dim2, dim1, totRuns);
data.DAQ = cell(totRuns,1);

% Insert first data set into the data array since it's already in memory
data.MCT(:,:,1) = datafile1.rawArray';
data.DAQ{1} = datafile1.rawDAQData;

clear datafile1 dim1 dim2;

% Load all the files in the experiment file list. Start at 2 since I already preloaded the first file after the preallocation step
for fileNum = 2:totRuns
    rawdata = load([dataFiles(fileNum).folder '/' dataFiles(fileNum).name]);
    data.MCT(:,:,fileNum) = rawdata.rawArray'; % transpose it to be pixels x shots
    data.DAQ{fileNum,1} = rawdata.rawDAQData; % had to use a cell since
    disp(['Loaded data file ' num2str(fileNum) '/' num2str(totRuns)]);
end

data.MCT = data.MCT(1:128,:,:); %the the first 128 pixels

end

function getZip(dataFolder, dataUrl)
    tempZipFile = 'blank_shots_dataset_2020_07_31.zip';
    fprintf("Downloading dataset from baizgroup.org ...")
    websave([dataFolder '/' tempZipFile], dataUrl);
    fprintf("done.\n")
    
    fprintf("Extracting zip...")
    unzip([dataFolder '/' tempZipFile],dataFolder);
    fprintf("done.\n")
end
