function [ flScene, fluorophoreIDs ] = fluorescentSceneCreate( varargin )

% [ flScene, fluorophoreIDs ] = fluorescentSceneCreate(...)
%
% Create a scene (i.e. a spatial arrangement) of different fluorophores. The
% scene is effectively a (h x w x n) array of fluorophore structures, where
% h and w represent the height and width respectively, n describes the
% number of fluorescent compounds per spatial location. Current
% implementation does not explicitly permit different spatial locations to have
% different number of fluorophores. This restriction can be bypassed by 
% assigning two or more identical fluorophores to the same spatial location
% and reducing their emission intensity accordingly.
%
% Inputs:
%   'type' - describe the fluorescent scene type to be generated. Available
%      values are
%      'fromfluorophore' - create a scene using a single fluorophore, i.e.
%         all spatial locations have the same fluorescent properties. The user
%         is expected to provide the 'fluorophore' input.
%      'one fluorophore' - create a scene with a single fluorophore (same
%         as above) selected from one of the data sets. Additional expected
%         inputs: 'dataset', 'stokesShiftRange', 'peakEmRange','peakExRange',
%         'fluorophoreIDs'
%      'default' - create a (h x w x n) scene of fluorophores. Additional 
%         expected inputs: 'dataset','stokesShiftRange','peakEmRange',
%         'peakExRange','height','width','nFluorophores','qe'
%   'dataset' - select the dataset from which fluorescent properties are
%      selected. Available values are {'McNamara-Boswell','LifeTechnologies'}
%   'wave' - a (w x 1) vector of wavelength sampling (default = 400:10:700)
%   'name' - a string representing the name of the fluorescent scene.
%   'stokesShiftRange' - a (2 x 1) vector providing the minimum and maximum
%      Stokes shift of fluorophores selected from the dataset (default = [0
%      Inf], i.e. there are no constraints on the Stokes shift).
%   'peakEmRange' - a (2 x 1) vector providing the minimum and maximum peak
%      emission wavelength of the fluorophores selected from the dataset,
%      (default = [420 680]).
%   'peakExRange' - a (2 x 1) vector providing the minimum and maximum peak
%      excitation wavelength of the fluorophores selected from the dataset,
%      (default = [420 680]).
%   'height' - the height (in fluorophores) of the scene (default = 1).
%   'width' - the width (in fluorophores) of the scene (default = 1).
%   'nFluorophores' - the number of fluorophores per spatial location
%      (default = 1).*** IMPORTANT
%   'qe' - the quantum efficiency of a particular spatial location. ***If more
%      than one fluorophore per location is present then the qe is divided
%      between number of fluorophores to represent the qe of the 'spatial
%      location'***IMPORTANT
%   'fluorophoreIDs' - a vector of flurophore identifiers that are read from
%      the databse and used to create the scene (rather than using random 
%      sampling).  ***IMPORTANT
%   'fluorophore' - a fluorophore structure used to generate a scene with
%      one fluorophore.
%
% Outputs:
%   'flScene' - the fluorescent scene structure.
%   'fluorophoreIDs' - the database identifiers of the fluorescent
%      compounds used to create the scene.
%
% Copyright, Henryk Blasinski 2016

p = inputParser;

% Sets all values to default
p.addParamValue('type','default',@ischar);
p.addParamValue('dataset','McNamara-Boswell',@(x) strcmp(x,validatestring(x,{'McNamara-Boswell','LifeTechnologies'})))
p.addParamValue('wave',400:10:700,@isvector);
p.addParamValue('name','Default',@ischar);
p.addParamValue('stokesShiftRange',[0 Inf],@isvector);
p.addParamValue('peakEmRange',[420 680],@isvector);
p.addParamValue('peakExRange',[420 680],@isvector);
p.addParamValue('height',1,@isscalar);
p.addParamValue('width',1,@isscalar);
p.addParamValue('nFluorophores',1,@isscalar);
p.addParamValue('qe',1,@isnumeric);
p.addParamValue('fluorophoreIDs',1,@isnumeric);
p.addParamValue('fluorophore',[],@isstruct);

p.parse(varargin{:});
inputs = p.Results;

flScene.type = 'fluorescent scene';
flScene.name = inputs.name;
flScene = initDefaultSpectrum(flScene,'custom',inputs.wave);

switch inputs.type

    %Single Fluorophore for the entire chart, user provides fluorophore
    %name
    case {'fromfluorophore'}
        flScene = fluorescentSceneSet(flScene,'fluorophores',inputs.fluorophore);        
        flScene = fluorescentSceneSet(flScene,'qe',inputs.qe);
        fluorophoreIDs = [];

    %Single 
    case {'onefluorophore','singlefluorophore'}
        setName = fullfile(fiToolboxRootPath,'data',inputs.dataset);
        flSet = fiReadFluorophoreSet(setName,'wave',inputs.wave,...
            'stokesShiftRange',inputs.stokesShiftRange,...
            'peakEmRange',inputs.peakEmRange,...
            'peakExRange',inputs.peakExRange);

        flScene = fluorescentSceneSet(flScene,'fluorophores',flSet(inputs.fluorophoreIDs));        
        flScene = fluorescentSceneSet(flScene,'qe',inputs.qe/inputs.nFluorophores);
        fluorophoreIDs = inputs.fluorophoreIDs;
    
    case 'default'
        % Create a default fluorescent scene
        setName = fullfile(fiToolboxRootPath,'data',inputs.dataset);

        %Reads all fluorophores in the dataset which meet the requirements
        %given
        [flSet, fluorophoreIDs] = fiReadFluorophoreSet(setName,'wave',inputs.wave,...
            'stokesShiftRange',inputs.stokesShiftRange,...
            'peakEmRange',inputs.peakEmRange,...
            'peakExRange',inputs.peakExRange);
        
        %Total number of fluorophores used across the entire scene (not
        %just one spatial location)
        nFluorophores = inputs.height*inputs.width*inputs.nFluorophores;
        
        %Creates random integers to select fluorophores
        ids = randi(length(flSet),nFluorophores,1);
        selFl = flSet(ids);
        fluorophoreIDs = fluorophoreIDs(ids);
        
        flScene = fluorescentSceneSet(flScene,'fluorophores',reshape(selFl,[inputs.height, inputs.width, inputs.nFluorophores]));        
        flScene = fluorescentSceneSet(flScene,'qe',inputs.qe/inputs.nFluorophores);
        
end



end

