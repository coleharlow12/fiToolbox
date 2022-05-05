function [ refl, fl ] = fiComputeReflFlContrib( camera, illuminant, cameraGain, reflectance, dMat )

% [ refl, fl ] = fiComputeReflFlContrib( camera, illuminant, cameraGain, reflectance, dMat )
%
% Given a model of a camera defined by the camera spectral responsivities
% and gains, the illuminant, and surface reflectance and fluorescence 
% properties compute the pixel intensities due to reflected and fluoresced
% light components.
%
% Inputs:
%    camera - a (w x f) matrix containing the spectral responsivity
%      functions of the f camera filters sampled at w wavebands.
%    illuminant - a (w x s) matrix containing the spectral power
%      distributions of s illuminants sampled at w wavebands.
%    cameraGain - a (f x c) matrix of linear camera model gains for each
%      filter-channel combination.
%    reflectance - a (w x s) matrix of surface spectral reflectance of s
%      surfaces sampled at w wavebands.
%    dMat - a (s x 1) cell array containing Donaldson matrices for each of
%      the s surfaces.
%
% Outputs:
%    refl - a (f x c x s) matrix containing pixel intensities of the
%       reflected light component of s surfaces under c lights and captured
%       with f filters.
%    fl - a (f x c x s) matrix containing pixel intensities of the
%       fluoresced light component of s surfaces under c lights and captured
%       with f filters.
%
% Copyright, Henryk Blasinski 2016


nSamples = size(reflectance,2);
nFilters = size(camera,2);
nChannels = size(illuminant,2);

refl = zeros(nFilters, nChannels, nSamples);
fl = zeros(nFilters, nChannels, nSamples);

for s=1:nSamples
    refl(:,:,s) = cameraGain.*(camera'*diag(reflectance(:,s))*illuminant);
    fl(:,:,s) = cameraGain.*(camera'*dMat{s}*illuminant);
end


end

