close all;
clear all;
clc;

% Evaluate the accuracy on the entire dataset

dataset = 'McNamara-Boswell';
wave = 380:4:1000;
deltaL = wave(2) - wave(1);
nWaves = length(wave);

flQe = 0.1;

alpha = 0.1;
beta = 0.1;
nu = 0.1;

% Create basis function sets
nReflBasis = 5;
nExBasis = 12;
nEmBasis = 12;

[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = createBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = createBasisSet('emission','wave',wave','n',nEmBasis);


% Load the light spectra (in photons)
fName = fullfile(fiToolboxRootPath,'camera','illuminants');
illuminant = ieReadSpectra(fName,wave);
illuminantPhotons = Energy2Quanta(wave,illuminant);
nChannels = size(illuminant,2);

% Load camera spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters');
filters = ieReadSpectra(fName,wave);

fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters;
nFilters = size(camera,2);
       

% Define which fluorophores will be used
setDir = fullfile(fiToolboxRootPath,'data',dataset);
[fluorophores, ids] = fiReadFluorophoreSet(setDir,'wave',wave,...
                           'peakEmRange',[wave(5) wave(end-5)],...
                           'peakExRange',[wave(5) wave(end-5)]);
nCompounds = length(ids);


% Create reflective scene
scene = sceneCreate('macbethEE_IR','',wave);
scene = sceneSet(scene,'fov',5);
scene = sceneSet(scene,'distance',1);

fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
reflRef = ieReadSpectra(fName,wave);




totalPixelErr = zeros(nCompounds,1);
reflPixelErr = zeros(nCompounds,1);
flPixelErr = zeros(nCompounds,1);
reflErr = zeros(nCompounds,1);
dMatErr = zeros(nCompounds,1);
dMatNormErr = zeros(nCompounds,1);

totalPixelStd = zeros(nCompounds,1);
reflPixelStd = zeros(nCompounds,1);
flPixelStd = zeros(nCompounds,1);
reflStd = zeros(nCompounds,1);
dMatStd = zeros(nCompounds,1);
dMatNormStd = zeros(nCompounds,1);

%% The main cross-validation loop

try
    matlabpool open local
catch 
end

for i=1:nCompounds

    % Create fluorescent scene
    flScene = fluorescentSceneCreate('type','fromfluorophore','wave',wave,'qe',flQe,...
                                     'fluorophore',fluorophores(i));
         
    dMatRef = fluorescentSceneGet(flScene,'Donaldson reference','sceneSize',[4 6]);        
    exRef = fluorescentSceneGet(flScene,'excitation reference','sceneSize',[4 6]); 
    emRef = fluorescentSceneGet(flScene,'emission reference','sceneSize',[4 6]); 
    
    
    %% Run ISET simulations
    
    cameraExposure = zeros(nFilters,nChannels);
    cameraGain = zeros(nFilters,nChannels);
    cameraOffset = zeros(nFilters,nChannels);
    
    measVals = zeros(nFilters,nChannels,24);
    
    for f=1:nFilters
        
        [sensor, optics] = createCameraModel(f,'wave',wave);
        
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
            sensor = sensorSetSizeToFOV(sensor,FOV,localScene,oi);
            sensor = sensorSet(sensor,'Name',sprintf('Filter %i, channel %i',f,ch));
            % sensor = sensorSet(sensor,'noise flag',0);
            
            cameraExposure(f,ch) = autoExposure(oi,sensor,0.95,'luminance');
            sensor = sensorSet(sensor,'exposureTime',cameraExposure(f,ch));
            
            [cameraGain(f,ch), cameraOffset(f,ch)] = sensorGainAndOffset(localScene,oi,sensor);
            
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
    
    [reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminantPhotons,cameraGain*deltaL,reflRef,dMatRef);

    
    %% Estimation
    
    cameraGain = repmat(cameraGain,[1 1 24]);
    cameraOffset = repmat(cameraOffset,[1 1 24]);
    
    
    [ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndMultiFl( measVals, camera, illuminantPhotons, cameraGain*deltaL,...
        cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, nu, 'maxIter',250);
    
    
    %% Evaluation
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    [totalPixelErr(i), totalPixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,24]), reshape(measVals,[nChannels*nFilters,24]), '');
    [reflPixelErr(i), reflPixelStd(i)] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,24]), reshape(reflValsRef,[nChannels*nFilters,24]), '');
    [flPixelErr(i), flPixelStd(i)] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,24]), reshape(flValsRef,[nChannels*nFilters,24]), '');
    
    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, '');
    
    [dMatErr(i), dMatStd(i)] = fiComputeError(dMatEst, dMatRef, '');
    [dMatNormErr(i), dMatNormStd(i)] = fiComputeError(dMatEst, dMatRef, 'normalized');
    
    %% ISET cleanup, remove all objects.

    scenes = vcGetObjects('scene');
    vcDeleteSomeObjects('scene',1:length(scenes));
    oi = vcGetObjects('oi');
    vcDeleteSomeObjects('oi',1:length(oi));
    sensor = vcGetObjects('sensor');
    vcDeleteSomeObjects('sensor',1:length(sensor));

end

try
    matlabpool close
catch
end

%% Save results
dirName = fullfile(fiToolboxRootPath,'results');
if ~exist(dirName,'dir'), mkdir(dirName); end

fName = fullfile(dirName,[dataset '_sim_multiFl.mat']);

save(fName,'fluorophores','ids',...
           'totalPixelErr','reflPixelErr','flPixelErr','reflErr','dMatErr','dMatNormErr',...
           'totalPixelStd','reflPixelStd','flPixelStd','reflStd','dMatStd','dMatNormStd');

       
%% Plot the results
% Sort the fluorophores in the increasing pixel prediction error order.

saveDir = fullfile('~','Dropbox','MsVideo','Notes','FluorescencePaperV2','Figures');
lineStyle = {'rs-','gd-','bo-'};
lw = 2;
mkSz = 8;
fs = 8;
figSize = [0 0 4.5 2.75];
[sortedTotalPixelErr, indx] = sort(totalPixelErr,'ascend');
nFluorophores = length(fluorophores);

figure;
hold all;
pl(1) = plot(1:nFluorophores,sortedTotalPixelErr,'lineWidth',lw);
pl(2) = plot(1:nFluorophores,sortedTotalPixelErr - totalPixelStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(3) = plot(1:nFluorophores,sortedTotalPixelErr + totalPixelStd(indx)/sqrt(24),'lineWidth',0.5*lw);

pl(4) = plot(1:nFluorophores,reflErr(indx),'lineWidth',lw);
pl(5) = plot(1:nFluorophores,reflErr(indx) - reflStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(6) = plot(1:nFluorophores,reflErr(indx) + reflStd(indx)/sqrt(24),'lineWidth',0.5*lw);

pl(7) = plot(1:nFluorophores,dMatNormErr(indx),'lineWidth',lw);
pl(8) = plot(1:nFluorophores,dMatNormErr(indx) - dMatNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);
pl(9) = plot(1:nFluorophores,dMatNormErr(indx) + dMatNormStd(indx)/sqrt(24),'lineWidth',0.5*lw);


set(pl(1:3),'markerSize',mkSz,'color','red');
set(pl(4:6),'markerSize',mkSz,'color','green');
set(pl(7:9),'markerSize',mkSz,'color','blue');


set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition',figSize);
xlabel('Fluorophore index');
ylabel('RMSE');
set(gca,'fontSize',fs);
% print('-depsc',fullfile(saveDir,'multiFl_xVal_alpha.eps'));
