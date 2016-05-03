close all;
clear all;
clc;

wave = 380:4:1000;
nWaves = length(wave);
deltaL = wave(2) - wave(1);

alpha = 0.001;
beta = 0.001;

flQe = 0.5;
nFluorophores = 1;
dataset = 'McNamara-Boswell';
height = 4;
width = 6;
nSamples = height*width;

% Reflectance basis
nReflBasis = 5;
[reflBasis, reflScore] = createBasisSet('reflectance','wave',wave','n',nReflBasis);


% Load the light spectra (in photons)
illuminant = eye(nWaves);
nChannels = size(illuminant,2);

% Load camera spectral properties
camera = eye(nWaves);
nFilters = size(camera,2);

cameraGain = ones(nFilters,nWaves);
cameraOffset = zeros(nFilters,nWaves);

% Create reflective scene
scene = sceneCreate('macbethEE_IR','',wave);
scene = sceneSet(scene,'fov',5);
scene = sceneSet(scene,'distance',1);

fName = fullfile(isetRootPath,'data','surfaces','macbethChart');
reflRef = ieReadSpectra(fName,wave);

% Create fluorescent scene
flScene = fluorescentSceneCreate('height',height,'width',width,'wave',wave,'qe',flQe,'nFluorophores',nFluorophores,...
                                 'peakEmRange',[wave(5) wave(end-5)],...
                                 'peakExRange',[wave(5) wave(end-5)],...
                                 'dataSet',dataset);
                             
dMatRef = fluorescentSceneGet(flScene,'Donaldson reference');         
exRef = fluorescentSceneGet(flScene,'excitation reference');
emRef = fluorescentSceneGet(flScene,'emission reference');

[reflValsRef, flValsRef] = fiComputeReflFlContrib(camera,illuminant,cameraGain*deltaL,reflRef,dMatRef);
measVals = reflValsRef + flValsRef;

cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

nF = max(max(measVals,[],1),[],2);
nF = repmat(nF,[nFilters,nChannels,1]);
measVals = measVals./nF;
cameraGain = cameraGain./nF;



maxExBasis = 25;
maxEmBasis = 25;
[emBasisGrid, exBasisGrid] = meshgrid(1:maxEmBasis,1:maxExBasis);



reflErr = zeros(maxEmBasis,maxExBasis);
reflStd = zeros(maxEmBasis,maxExBasis);

exNormErr = zeros(maxEmBasis,maxExBasis);
exNormStd = zeros(maxEmBasis,maxExBasis);

emNormErr = zeros(maxEmBasis,maxExBasis);
emNormStd = zeros(maxEmBasis,maxExBasis);

pixelErr = zeros(maxEmBasis,maxExBasis);
pixelStd = zeros(maxEmBasis,maxExBasis);


try
    cluster = parcluster('local');
    cluster.NumWorkers = min(numel(emBasisGrid),35);
    pool = parpool(cluster,cluster.NumWorkers);
catch
end

parfor i=1:numel(emBasisGrid);


    exBasis = createBasisSet('excitation','wave',wave','n',exBasisGrid(i));
    emBasis = createBasisSet('emission','wave',wave','n',emBasisGrid(i));
    
    
    [ reflEst, ~, emEst, ~, exEst, ~, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndFl( measVals, camera, cameraGain*deltaL, cameraOffset, illuminant,...
        reflBasis, emBasis, exBasis, alpha, beta, beta, 'maxIter', 250 );

    measValsEst = reflValsEst + flValsEst + cameraOffset;

    [pixelErr(i), pixelStd(i)] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'default');

    [reflErr(i), reflStd(i)] = fiComputeError(reflEst, reflRef, '');

    [exNormErr(i), exNormStd(i)] = fiComputeError(exEst, exRef, 'normalized');
    [emNormErr(i), emNormStd(i)] = fiComputeError(emEst, emRef, 'normalized');
            

end

fName = fullfile(fiToolboxRootPath,'results','evaluation','fl_nBasis.mat');
save(fName,'pixelErr','pixelStd','exNormErr','exNormStd','emNormErr','exNormStd','reflErr','reflStd',...
           'exBasisGrid','emBasisGrid','alpha','beta','exRef','emRef','reflRef',...
           'nReflBasis');


try
    delete(pool);
catch
end
