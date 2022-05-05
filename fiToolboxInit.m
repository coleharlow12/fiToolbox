function fiToolboxInit()

% fiToolboxInit()
%
% Initializes the Fluorescence Imaging Toolbox (fiToolbox). Adds sub
% directories to MATLAB path.
%
% Copytight, Henryk Blasinski 2016

close all;
clear all;
clc;

warning off;

addpath(fullfile(fiToolboxRootPath));

addpath(fullfile(fiToolboxRootPath,'camera'));
addpath(fullfile(fiToolboxRootPath,'data'));
addpath(fullfile(fiToolboxRootPath,'results'));

addpath(fullfile(fiToolboxRootPath,'estimation'));
addpath(fullfile(fiToolboxRootPath,'estimation','multistep'));
addpath(fullfile(fiToolboxRootPath,'estimation','nucnorm'));


addpath(fullfile(fiToolboxRootPath,'fluorescence'));
addpath(fullfile(fiToolboxRootPath,'io'));

addpath(fullfile(fiToolboxRootPath,'scripts'));
addpath(fullfile(fiToolboxRootPath,'scripts','validation'));
addpath(fullfile(fiToolboxRootPath,'scripts','simulations'));
addpath(fullfile(fiToolboxRootPath,'scripts','experiments'));
addpath(fullfile(fiToolboxRootPath,'scripts','plots'));

%Added by Cole Harlow in order to simulate coral response to light
addpath(fullfile(fiToolboxRootPath,'scripts','CoralSimulations'));
addpath(fullfile(fiToolboxRootPath,'scripts','CoralSimulations','CoralFluorescence'));
addpath(fullfile(fiToolboxRootPath,'scripts','CoralSimulations','CoralReflectance'));
addpath(fullfile(fiToolboxRootPath,'utilities'));

end



