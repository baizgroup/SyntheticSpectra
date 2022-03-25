% 9/10/2020
% Script to train ANNs on gigantic data set
% Data set is 102 files of data, each file is 10000 blank shots
% total training set is 1020000 blank shots! 

%close all; clear; clc;
%% %%%%%%% User Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath C:\Box\BaizGroup_Shared\Projects\Referencing\Shot-to-shot_Autocorrelations\main_scripts

% Separate script to load variables shared across all main_netcorr scripts
sharedVars;

numDataSets = 100;
in.numShots = 5000; %%%%%%%%%%%%%%%

% File path to all shots experiment file list
% pathName.expAS = [pathName.dataFolder 'allShotsList_giantblank_train_paths.txt']; % giant all shots blank data collected at -10000 fs
pathName.expAS = [pathName.dataFolder 'allShotsList_giantblank_train_paths.txt']; % giant all shots blank data collected at -10000 fs

% File name to save the output file as
pathName.outputName = [pathName.outputFolder 'main_netcorr_train_giantblank_output_'...
    num2str(ANNsettings.options.MaxEpochs) '_' num2str(ANNsettings.options.MiniBatchSize) '_' ANNsettings.netArch '_' num2str(numDataSets) '_sets_mOD.mat'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1. Load data 

% Load all shots experiment file (list of data files names)
explist.expAS = importdata(pathName.expAS); % all shots data exp files
% Fix file paths in experiment list if necessary
explist.expAS = changeFilePaths(in, explist.expAS);

% Load all shots data, blank data for training
trainData.arrayAS = dataLoaderAllShotsSingleDelay(in, explist.expAS(1:numDataSets)); 
