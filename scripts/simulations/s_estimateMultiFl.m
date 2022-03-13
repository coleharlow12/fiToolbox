% Use the multi-fluorophore algorithm to estimate the reflectance and fluorescence
% properties from simulated data.
%
% Copyright, Henryk Blasinski 2016

close all;
clear variables;
clc;

% Load data simulation data (Need to figure out how to do my own
% simulations)
inFName = 'McNamara-Boswell_4x6x1_qe_0.10';
fName = fullfile(fiToolboxRootPath,'data','simulations',[inFName '.mat']);
load(fName);

% The file contents include the following

% cameraGain: Contains the gain of the camera during measurement for each
% of the 8 (filters) by 14 (illuminants)

% cameraOffset: Contains the offset during measurement for each of the 8
% (filters) by 14 (illuminants)

% measVals: numFilters x numIlluminants x numPixels contains the measured
% values

% wave: Contains the wavelengths at which the reference values are recorded

% reflRef: Contains the reflected value of the spectrum at each of the
% wavelengths

% emRef: Contains the emission of the pixels at each of the wavelengths
% specified in wave(Ground Truth)

% exRef: Contains the excitation of the pixels at each of the wavelengths
% (Grount Truth)

% flValsRef: Fluorescence Values Measured for each of the filters and
% illumants tested (Ground Truth)

% reflValsRef: Reflected Values measured for each of the filters and
% illuminants tested (Ground Truth)

% dMatRef: Donaldson matrix ground truth for each of the 24 pixels (in this
% case the pixels are really each of the squares of the macbeth chart

% Wavelength sampling is specified in the above loaded data file.
deltaL = wave(2) - wave(1);
nWaves = length(wave);

maxIter = 100;

% Create basis function sets
nReflBasis = 5; %The number of reflectance basis functions
nExBasis = 12; %The number of excitation basis
nEmBasis = 12; %The number of emission basis


[reflBasis, reflScore] = fiCreateBasisSet('reflectance','wave',wave','n',nReflBasis);
[exBasis, exScore] = fiCreateBasisSet('excitation','wave',wave','n',nExBasis);
[emBasis, emScore] = fiCreateBasisSet('emission','wave',wave','n',nEmBasis);


% Load the light spectra (in photons)
fName = fullfile(fiToolboxRootPath,'camera','illuminants'); %The file illuminants contains the specta of 14 different LEDs
illuminant = ieReadSpectra(fName,wave); %Interpolates the spectrums of the illuminants
illuminant = Energy2Quanta(wave,illuminant); %converts from watts/m^2 to photons/m^2 at each wavelength
nChannels = size(illuminant,2); %Gets the number of illuminants

% Load camera filter spectral properties
fName = fullfile(fiToolboxRootPath,'camera','filters'); 
filters = ieReadSpectra(fName,wave);

% Reads quantum efficiency of the camera
fName = fullfile(fiToolboxRootPath,'camera','qe');
qe = ieReadSpectra(fName,wave);

camera = diag(qe)*filters; % Scales filters by camera quantum efficiency
                           % Basically the camera quantum efficency is in
                           % series with the filters so the ultimate effect
                           % at each frequency is the multiplication of the
                           % filter qe and the camera qe
      
nFilters = size(camera,2); %Gets the number of filters
       
%% Load simulation data 
nSamples = size(measVals,3);
cameraGain = repmat(cameraGain,[1 1 nSamples]);
cameraOffset = repmat(cameraOffset,[1 1 nSamples]);

for i=1:3
    %Used to test different values of the tuning parameters alpha, beta, eta
    switch i
        case 1
            alpha = 0;
            beta = 0;
            eta = 0.01;
        case 2
            alpha = 0.01;
            beta = 0.01;
            eta = 0.0;
        case 3
            alpha = 0.01;
            beta = 0.01;
            eta = 0.01;
    end
    
    % This script does the estimation
    [ reflEst, reflCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMatEst, reflValsEst, flValsEst, hist  ] = ...
        fiRecReflAndMultiFl( measVals, camera, illuminant, cameraGain*deltaL,...
        cameraOffset, reflBasis, emBasis, exBasis, alpha, beta, beta, eta, 'maxIter',maxIter,...
        'dMatRef',dMatRef,'reflRef',reflRef,'pixelRef',true);
    
    
    
    %% Compute errors
    
    measValsEst = reflValsEst + flValsEst + cameraOffset;
    
    fprintf('====== alpha=%.3f beta=%.3f eta=%.3f======\n',alpha, beta,eta);
    
    
    [err, std] = fiComputeError(reshape(measValsEst,[nChannels*nFilters,nSamples]), reshape(measVals,[nChannels*nFilters,nSamples]), 'absolute');
    fprintf('Total pixel error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(reshape(reflValsEst,[nChannels*nFilters,nSamples]), reshape(reflValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    fprintf('Reflected pixel error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(reshape(flValsEst,[nChannels*nFilters,nSamples]), reshape(flValsRef,[nChannels*nFilters,nSamples]), 'absolute');
    fprintf('Fluoresced pixel error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(reflEst, reflRef, 'absolute');
    fprintf('Reflectance error %.3f, std %.3f\n',err,std);
    
    [err, std] = fiComputeError(dMatEst, dMatRef, 'absolute');
    fprintf('Donaldson Matrix error %f, std %f\n',err,std);
    
    [err, std] = fiComputeError(dMatEst, dMatRef, 'normalized');
    fprintf('Donaldson Matrix (normalized) %.3f, std %.3f\n',err,std);
    
end


%% Plot the results

figure(1);
hold all; grid on; box on;
plot(measValsEst(:),measVals(:),'.');
xlim([0 1]);
ylim([0 1]);
xlabel('Model predicted pixel value');
ylabel('ISET pixel value');

% Convergence
figure(2);
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot([hist{sampleID}.prRes, hist{sampleID}.dualRes],'LineWidth',2);
    xlim([0 length(hist{sampleID}.prRes)]);
    ylim([1e-5 10]);
    set(gca,'yscale','log');
end
end

% Estimate convergence
pixelErr = zeros(maxIter,1);
dMatErr = zeros(maxIter,1);
reflErr = zeros(maxIter,1);
for i=1:nSamples
    pixelErr = pixelErr + hist{i}.pixelErr;
    dMatErr = dMatErr + hist{i}.dMatErr;
    reflErr = reflErr + hist{i}.reflErr;
end

figure(3); 
hold on; grid on; box on;
plot([pixelErr, dMatErr, reflErr]/nSamples);
xlim([1 maxIter]);
set(gca,'yscale','log');
set(gca,'xscale','log');
legend({'Pixel','Donaldson','Reflectance'});


% Estimated vs. ground truth reflectance
figure(4);
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;
    plot(wave,reflEst(:,sampleID),'g','LineWidth',2);
    plot(wave,reflRef(:,sampleID),'b--','LineWidth',2);
    xlim([min(wave) max(wave)]);
    ylim([-0.05 1.05]);

    rmse = sqrt(mean((reflEst(:,sampleID) - reflRef(:,sampleID)).^2));
    title(sprintf('RMSE %.2f',rmse));

end
end

% Estimated vs. ground truth Donaldson matrices: scatter plot
figure(5);
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    hold all; grid on; box on;

    plot(dMatEst{sampleID}(:),dMatRef{sampleID}(:),'.');

    rmse = sqrt(mean((dMatEst{sampleID}(:) - dMatRef{sampleID}(:)).^2));
    title(sprintf('RMSE %.2e',rmse));
end
end

% Estimated vs. ground truth Donaldson matrics: shape
figure(6);
for xx=1:6
for yy=1:4

    plotID = (yy-1)*6 + xx;
    sampleID = (xx-1)*4 + yy;

    subplot(4,6,plotID);
    
    data = [dMatEst{sampleID} dMatRef{sampleID}];
    imagesc(data);

end
end
