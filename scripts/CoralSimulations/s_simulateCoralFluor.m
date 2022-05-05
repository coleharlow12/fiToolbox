% Use Image Systems Engineering Toolbox to simulate a set of captures using
% a model of the image acquisition system. Create a fluorescent scene with
% the fiToolbox functions. The captured data contains simulated pixel intensities
% for each of the 24 patches illuminated with each of the 14 narrowband 
% lights and as seen through each of the 8 camera filters.
%
% Copyright, Henryk Blasinski 2016.

% The following variables should be saved after running this script
% cameraGain
%dMatRef
%emRef
%exRef
%flValsRef
%measVals
%reflRef
%reflValsRef
%wave

close all;
clear all;
clc;

ieInit;

dataset = 'McNamara-Boswell';   % Database from which fluorophores are selected
flQe = 1;                       % Fluorescence practical quantum efficiency
                                    % To get this value for Corals see
                                    % Readme.txt in CoralRead
height = 4;                     % Number of test patches (height x width)
width = 6;
nFluorophores = 1;              % Number of fluorophores per test patch

% Save data to file if saveFName ~= []
%dirName = fullfile(fiToolboxRootPath,'data','simulations');
% saveFName = fullfile(dirName,sprintf('%s_%ix%ix%i_qe_%0.2f.mat',dataset,height,width,nFluorophores,flQe));
saveFName = [];

% Wave lengths to perform calculations at
wave = 380:4:1000;
nWaves = length(wave);
deltaL = wave(2) - wave(1);

%% Measurement System Parameters

% Load the illuminants light spectra
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

% Load Camera quantum efficiency
fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);

%% Scene Creation
% Create reflective scene *** NEED TO UNDERSTAND ***
% 1) Don't worry about illuminant that is set later
% 2) The filter is also set later

scene = sceneCreate('macbethEE_IR','',wave);
%scene = sceneCreate('')

% Sets field of view to 5 degrees
scene = sceneSet(scene,'fov',5);
% Sets the scene distance to 1m
scene = sceneSet(scene,'distance',1);

% Location of the spectral data to be read
fName = fullfile(isetRootPath,'data','surfaces','reflectances','macbethChart');

% Interpolates the referenced file to the wavelengths examined, this is the
% known reference. Note that the sceneCreate function references the same
% data
reflRef = ieReadSpectra(fName,wave);

% Create fluorescent scene. *** NEED TO UNDERSTAND ***
% We select fluorophores with peak excitation and emission somewhat within
% the minimum and maximum spectral sampling wavelengths.
% Cole: this function will need to be expanded to load in specific
% emission/excitation spectrums
flScene = fluorescentSceneCreate('height',height,'width',width,'wave',wave,...
                                 'qe',flQe,'nFluorophores',nFluorophores,...
                                 'peakEmRange',[wave(5) wave(end-5)],...
                                 'peakExRange',[wave(5) wave(end-5)],...
                                 'dataSet',dataset);
         
dMatRef = fluorescentSceneGet(flScene,'Donaldson reference');        
exRef = fluorescentSceneGet(flScene,'excitation reference');
emRef = fluorescentSceneGet(flScene,'emission reference');
           
%% Run ISET simulations

%Used to record camera exposure, gain, and offset for post processing
%One of each for each camera exposure for each filter/light combo
cameraExposure = zeros(nFilters,nChannels);
cameraGain = zeros(nFilters,nChannels);
cameraOffset = zeros(nFilters,nChannels);

measVals = zeros(nFilters,nChannels,24);

%Loops through each filter
for f=1:nFilters
    
    [sensor, optics] = createCameraModel(f,'wave',wave);
    
    %Loops through each light source
    for ch = 1:nChannels
        
        fprintf('Simulating filter %i channel %i\n',f,ch);
        
        % Synthesize a fluorescent scene
        localScene = sceneAdjustIlluminant(scene,illuminant(:,ch),0);
        localScene = fiSceneAddFluorescence(localScene,flScene);

        localScene = sceneSet(localScene,'name',sprintf('Filter %i, channel %i',f,ch));
        ieAddObject(localScene);
        

        % Compute the optical image
        oi = oiCreate();
        oi = oiSet(oi,'optics',optics);
        oi = oiSet(oi,'Name',sprintf('Filter %i, channel %i',f,ch));
        oi = oiCompute(localScene,oi);
        ieAddObject(oi);
        

        % Compute the sensor image
        FOV = [sceneGet(localScene,'fov horizontal') sceneGet(localScene,'fov vertical')];
        sensor = sensorSetSizeToFOV(sensor,FOV,oi);
        sensor = sensorSet(sensor,'Name',sprintf('Filter %i, channel %i',f,ch));
        sensor = sensorSet(sensor,'noise flag',0);

        % Adjusts the camera exposure for each filter,illuminant combo
        cameraExposure(f,ch) = autoExposure(oi,sensor,0.95,'luminance');
        sensor = sensorSet(sensor,'exposureTime',cameraExposure(f,ch));

        % This can be done using sensorGet.m instead of the current file
        [cameraGain(f,ch), cameraOffset(f,ch)] = sensorGainAndOffset(localScene,oi,sensor);
        %cameraGain(f,ch) = sensorGet(sensor,'analog Gain');
        %cameraOffset(f,ch) = sensorGet(sensor,'analog Offset');
        
        sensor = sensorCompute(sensor,oi);
        ieAddObject(sensor);
        
        
        % Read out pixel intensities
        sSize = sensorGet(sensor,'size');
        cornerPoints = [1 sSize(1); sSize(2) sSize(1); sSize(2) 1; 1 1];
       
        mVals = macbethSelect(sensor,0,1,cornerPoints);
        mVals = cellfun(@(x) nanmean(x)/(2^sensorGet(sensor,'nbits')),mVals);

        
        measVals(f,ch,:) = mVals;  
    end

end

sceneWindow;
oiWindow;
sensorWindow;

%% Compute linear model predictions

illuminantPhotons = Energy2Quanta(wave,illuminant);

[reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminantPhotons,cameraGain*deltaL,reflRef,dMatRef);

predVals = reflValsRef + flValsRef + repmat(cameraOffset,[1 1 24]);

% Create a scatter plot between the simulated data and the linear model
% data. If everything is OK the points should be on the identity line.
figure;
hold all; grid on; box on;
plot(predVals(:),measVals(:),'.');
plot(linspace(0,1,10),linspace(0,1,10),'r','LineWidth',2);
xlabel('Linear model pixel intensity');
ylabel('ISET simulation pixel intensity');

% Save data
if ~isempty(saveFName)
    save(saveFName,'cameraGain','cameraOffset','measVals','reflRef','dMatRef','exRef','emRef','reflValsRef','flValsRef',...
           'wave');
end