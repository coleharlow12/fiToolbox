function [ basis, score ] = fiCreateBasisSet( type, varargin )
% varargin allows for a variable length input argument list

% [ basis, score ] = fiCreateBasisSet( type, varargin )
%
% Creates a set of basis functions for spectral quantities.
%
% Inputs (required):
%    type - a string specifying the type of spectral quantity for which
%      basis functions are generated. Allowed values are {'reflectance',
%      'excitation','emission'}
%
% Inputs (optional):
%    'wave' - spectral quantities wavelength sampling (default =
%      400:10:700)
%    'n' - number of basis functions (default = 5);
%
% Outputs:
%    basis - a (w x n) array of linear basis functions, where w is the
%      number of wavelength samples.
%    score - a (n x 1) vector of variances captured by each linear basis
%      function.
%
% Copytight, Henryk Blasinski 2016


p = inputParser; %Creates an inputParser object
p.addRequired('type',@(x)(strcmp(x,validatestring(x,{'reflectance','excitation','emission'})))); %Ensures type input is either reflectance excitation or emission
p.addParamValue('wave',400:10:700,@isnumeric); %adds a parameter to the parser. If no value is provided the default is 400:10:700
p.addParamValue('n',5,@isscalar); %Adds a parameter n (the number of basis functions) to the parser. If no value is provided default is 5

p.parse(type,varargin{:});
inputs = p.Results; %Creates structure containing 'type' the wavelengths of data 'wave' and number of basis functions 'n'

switch type
    case 'reflectance'
        
        % Reads the data from macbethChart file. This file contains
        % reflectance data of each of the 24 squares on the macbeth chart
        fName = fullfile(isetRootPath,'data','surfaces','reflectances','macbethChart');
        % Uses extrapolation to get reflectance at input wavelengths
        data = ieReadSpectra(fName,inputs.wave);

    case 'excitation'
        
        fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
        fluorophores = fiReadFluorophoreSet(fName,'peakEmRange',[0 max(inputs.wave)],...
                                                  'peakExRange',[min(inputs.wave) Inf],...
                                                  'wave',inputs.wave);

        data = zeros(length(inputs.wave),length(fluorophores));
        for i=1:length(fluorophores)
            data(:,i) = fluorophoreGet(fluorophores(i),'normalized excitation');
        end

    case 'emission'
        
        fName = fullfile(fiToolboxRootPath,'data','McNamara-Boswell');
        fluorophores = fiReadFluorophoreSet(fName,'peakEmRange',[0 max(inputs.wave)],...
                                                  'peakExRange',[min(inputs.wave) Inf],...
                                                  'wave',inputs.wave);

        data = zeros(length(inputs.wave),length(fluorophores));
        for i=1:length(fluorophores)
            data(:,i) = fluorophoreGet(fluorophores(i),'normalized emission');
        end   

end

% Performs PCA. This means it creates 24 new variables to try and tries to
% reduce the dimensionality of the covariance matrix. Essentially what
% happens is this functon finds the covariance of each of the 24 squares
% with one another. It then stores this in a 24 by 24 matrix. It finds the
% eigenvalues of this and orders them by their ability to explain the
% variance in the data (i.e capture the information in the data). 

[basis, ~, score] = pca(data','centered',false); 
% Rows of x correspond observations and columns to variables.

% Choose as a basis the squares which account for the most variance
basis = basis(:,1:inputs.n); %Outputs the basis functions which give highest variance
score = score(1:inputs.n); %Output the score of the basis functions which give the highest variance

end

