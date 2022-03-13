function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMat, predRefl, predFl, hist  ] = fiRecOneReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin )

% function [ reflEst, rfCoeffs, emEst, emCoeffs, exEst, exCoeffs, dMat, predRefl, predFl, hist  ] = fiRecOneReflAndMultiFl( measVals, cameraMat, illuminant, cameraGain, basisRefl, basisEm, basisEx, alpha, beta, gamma, eta, varargin )

% This function implements the multi fluorophore model and 
% estimates the reflectance and the Donaldson matrix 
% from pixel intensities using an ADMM estimation algorithm. This 
% implementation contains tuning parameters to independently adjust the 
% smoothness of the excitation and emission spectra. 
%
% Inputs (required):
%    measVals - a (f x c) matrix containing pixel intensities of a 
%      surface captured with f camera channels and under c
%      different illuminants.
%    cameraMat - a (w x c) matrix containing the spectral responsivity
%      functions of the c camera channels sampled at w wavebands.
%    cameraGain - a (f x c) matrix of linear camera model gains for each
%      filter-channel combination.
%    basisRefl, basisEm, basisEx - a (w x n) matrices of c linear basis
%      functions representing reflectance, excitation and emission spectra
%      respectively
%    alpha - scalar tuning parameters controlling the smoothness of
%      reflectance estimates
%    beta - scalar tuning parameters controlling the smoothness of
%      fluorescence excitation estimates
%    gamma - scalar tuning parameters controlling the smoothness of
%      fluorescence emission estimates
%    eta - scalar tuning parameters controlling the nuclear norm of the
%      Donaldson matrix estimate
%
% Inputs (optional):
%    'reflRef' - a (w x s) matrix of reference surface reflectances. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'dMatRef' - a (s x 1) cell array of reference Donaldson matrices. If
%      provided the algorithm compute the error between the estimate at
%      iteration i, and the reference. This error is stored in the hist
%      structure. 
%    'pixelRef' - a boolean value indicating if the error between predicted
%      and mesured pixel intensities is to be computed at every iteration.
%      (default = false).
%   'maxIter' - the maximal number of ADMM interations (default = 500).
%   'rescaleRho' - implement a heuristic algorithm improving ADMM
%      convergence from Boyd 2011 (default = true).
%   'tauIncr' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 5).
%   'tauDecr' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 5).
%   'rhoInit' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 0.1).
%   'mu' - heuristic algorithm parameter (see Boyd, 2011) 
%      (default = 10).
%   'nFluorophores' - specify the number of principal components to
%      preserve in the Donaldson matrix estimates. If different than zero,
%      instead of soft-thresholding the singular values, onlty nFluorophore
%      largest singular values are preserved. Setting this parameter
%      makes the problem non convex (default = 0).
%   'epsilon' - the convergence threshold for the ADMM algorithm 
%      (default = 1e-6).
%
% Outputs:
%    reflEst - a (w x 1) vector of the estimated surface spectral reflectance.
%    rfCoeffs - a (s x 1) vector expressing the estimated surface spectral 
%      reflectance in terms of the linear basis weights.
%    emEst - a (w x 1) vector of the estimated emission spectrum (the first
%      principal component of the Donaldson matrix estimate).
%    emCoeffs - a (s x 1) vector expressing the estimated surface emission 
%      spectrum in terms of the linear basis weights.
%    exEst - a (w x 1) vector of the estimated surface excitation spectrum.
%      (the first principal component of the Donaldson matrix estimate).
%    exCoeffs - a (s x 1) vector expressing the estimated surface excitation 
%      spectrum in terms of the linear basis weights.
%    dMat - a (w x w) array representing the Donaldson matrix estimate.
%    predRefl - a (f x c) matrix of values representing the reflected
%      light intensities for each filter-illuminant combination.
%    predFl - a (f x c) matrix of values representing the fluoresced
%      light intensities for each filter-illuminant combination.
%      Specifically, for ideal, noiseless measurements the following holds
%              predRefl + predFl = (measVals - cameraOffset)
%    hist - structure containing the objective function values at 
%      successive minimization steps.
%
% Copyright, Henryk Blasinski 2016



p = inputParser;
p.addParamValue('maxIter',500);
p.addParamValue('tauIncr',5,@isscalar);
p.addParamValue('tauDecr',5,@isscalar);
p.addParamValue('rhoInit',0.1,@isscalar);
p.addParamValue('rescaleRho',true);
p.addParamValue('mu',10,@isscalar);
p.addParamValue('epsilon',1e-6);
p.addParamValue('nFluorophores',0);
p.addParamValue('dMatRef',[]);
p.addParamValue('reflRef',[]);
p.addParamValue('pixelRef',false);
p.parse(varargin{:});
inputs = p.Results;

nEmBasis = size(basisEm,2);         % number of emission basis
nExBasis = size(basisEx,2);         % number of excitment basis
nReflBasis = size(basisRefl,2);     % number of reflection basis
nWaves = size(illuminant,1);        % number of wavelengths the data is available at

scaleFac = max(cameraGain(:));      %Normalizes the camera gain 
cameraGain = cameraGain/scaleFac;   %Normalizes the camera gain 
illuminant = illuminant*scaleFac; 

%% Initialize The Algorithm
% Initialize variables to zero
y1 = zeros(nWaves,1);                  % Reflectance
y2 = zeros(nWaves);                    % Fluorescence excitation-emission matrix
y3 = zeros(nEmBasis,nExBasis);         % Fluorescence excitation-emission wieghts

y1minus = y1;
y2minus = y2;
y3minus = y3;

% Dual variables
u1 = zeros(nWaves,1);
u2 = zeros(nWaves);
u3 = zeros(nEmBasis,nExBasis); 

% Estimate
sol = [];

% In general we use PCG algorithm to solve least-squares problem. Sometimes
% this algorithm is faster if it operates on matrices, other times function
% handles perform better. The commented line below allows to generate some
% of the matrices used in the matrix approach.
% [~, ~, AtA1, Atb1] = generateLSMatrices(measVals,cameraMat,illuminant,cameraGain,basisRefl,basisEx,basisEm, alpha, beta, gamma);

lowerTr = logical(tril(ones(nWaves),-1)); %Extracts the lower triangle of a matrix

%Stores data in a histogram
hist.rho = ones(inputs.maxIter+1,1);
hist.prRes = zeros(inputs.maxIter+1,1);
hist.dualRes = zeros(inputs.maxIter+1,1);
hist.pcgIter = zeros(inputs.maxIter,1);
hist.pcgAcc = zeros(inputs.maxIter,1);

if inputs.pixelRef == true
    hist.pixelErr = zeros(inputs.maxIter,1);
end

if ~isempty(inputs.dMatRef)
    hist.dMatErr = zeros(inputs.maxIter,1);
end

if ~isempty(inputs.reflRef)
    hist.reflErr = zeros(inputs.maxIter,1);
end

%% Iterative Solver

tElapsed = 0; t2 = tic;
for i=1:inputs.maxIter
    
    % Solve for reflectance and fluorescence
    
    t1 = tic;
    
    % Uncomment these lines to use matrix version of PCG
    % [~, ~, AtA2, Atb2] = generatePenaltyMatrices(y1-u1,y2-u2,y3-u3,basisRefl,basisEx,basisEm);
    % [sol, ~, hist.pcgAcc(i), hist.pcgIter(i)] = pcg(AtA1 + (hist.rho(i)/2)*AtA2,Atb1 + (hist.rho(i)/2)*Atb2,[],10000,[],[],sol);
    
    
    % Uncomment these lines to use matrix handle version of PCG
    b = applyAtb(measVals, cameraMat, illuminant, cameraGain, basisRefl, basisEx, basisEm, y1-u1, y2-u2, y3-u3, hist.rho(i)/2);
    AtAhndl = @(x) applyAtA(x, cameraMat, illuminant, cameraGain, basisRefl, basisEx, basisEm, hist.rho(i)/2, alpha, beta, gamma);
    
    %PCG Preconditioned Conjugate gradients method. Attempts to solve
    %linear equations A*X=B for Xs
    [sol, ~, hist.pcgAcc(i), hist.pcgIter(i)] = pcg(AtAhndl,b,[],10000,[],[],sol);
    tElapsed = tElapsed + toc(t1);
    
    
    reWEst = sol(1:nReflBasis);
    wEst = reshape(sol(nReflBasis+1:end),nEmBasis,nExBasis);
    
    % Solve for penalties
    % Box penalty on reflectance
    y1 = basisRefl*reWEst + u1;
    y1(y1 < 0) = 0;
    y1(y1 > 1) = 1;
    
    % Nonnegative penalty on fluorescence (i.e on the lower triangular part
    % only, we dont care about the upper part since we zero it anyway)
    y2 = basisEm*wEst*basisEx' + u2;
    y2(lowerTr & y2 < 0) = 0;

     
    % Rank penalty on fluorescence weights
    tmp = wEst+u3;
    [U,S,V] = svd(tmp);
    if ~isvector(S)
        s = diag(S);
    else
        s = S(1);
    end
    if inputs.nFluorophores <= 0
        s = sign(s).*max(abs(s)-(eta)/hist.rho(i),0);
    else
        s(inputs.nFluorophores+1:end) = 0;
    end
    
    % If the number of excitation end emission bases is different, then we
    % need to pad the S matrix accordingly
    y3 = U*padarray(diag(s),max([nEmBasis nExBasis] - length(s),0),0,'post')*V';
    
    % Update the scaled dual variables
    res1 = basisRefl*reWEst - y1;
    res2 = basisEm*wEst*basisEx' - y2;
    res3 = wEst - y3;
    u1 = u1 + res1;
    u2 = u2 + res2;
    u3 = u3 + res3;
    
    
    hist.prRes(i) = sqrt(norm(res1,'fro')^2 + norm(res2,'fro')^2 + norm(res3,'fro')^2);
    hist.dualRes(i) = sqrt(norm(y1-y1minus,'fro')^2 + norm(y2-y2minus,'fro')^2 + norm(y3-y3minus,'fro')^2);
    
    if (mod(i,10) == 0)
        
        tTotal = toc(t2);
        fprintf('Iter %i pr. res %f, dual res %f\n',i,hist.prRes(i),hist.dualRes(i));
        fprintf('     -> Total time: %.3f sec\n',tTotal);
        fprintf('     -> PCG time:   %.3f sec\n',tElapsed);

        t2 = tic;
        tElapsed = 0;
    end
    
    
    
    if max(hist.prRes(i),hist.dualRes(i)) < inputs.epsilon && (inputs.epsilon >=0)
        break;
    end
    
    % Re-scale the parameter rho
    if hist.prRes(i) > inputs.mu*hist.dualRes(i) && inputs.rescaleRho
        hist.rho(i+1) = hist.rho(i)*inputs.tauIncr;
    end
    if hist.dualRes(i) > inputs.mu*hist.prRes(i) && inputs.rescaleRho
        hist.rho(i+1) = hist.rho(i)/inputs.tauDecr;
    end
    
    u1 = u1*hist.rho(i)/hist.rho(i+1);
    u2 = u2*hist.rho(i)/hist.rho(i+1);
    u3 = u3*hist.rho(i)/hist.rho(i+1);

    
    y1minus = y1;
    y2minus = y2;
    y3minus = y3;
    
    hist.computeTime(i) = toc(t1);
    %% Compute approximation error
    % The opitmization is done, let's compute the estimates. All this is in
    % 'if-else' statements to speed up computation if convergence analysis
    % is not needed.

    if ~isempty(inputs.dMatRef)
        dMat = max(tril(basisEm*wEst*basisEx',-1),0);

        hist.dMatErr(i) = fiComputeError(dMat(:),inputs.dMatRef(:),'normalized');
    end
    
    if ~isempty(inputs.reflRef)
        rfCoeffs = reWEst;
        reflEst = basisRefl*rfCoeffs;
        
        hist.reflErr(i) = fiComputeError(reflEst,inputs.reflRef,'');
    end
    
    if inputs.pixelRef == true
        rfCoeffs = reWEst;
        reflEst = basisRefl*rfCoeffs;    
        
        dMat = max(tril(basisEm*wEst*basisEx',-1),0);    
        
        % Compute the prediction
        predRefl = cameraGain.*(cameraMat'*diag(reflEst)*illuminant);
        predFl = cameraGain.*(cameraMat'*dMat*illuminant);
        
        hist.pixelErr(i) = fiComputeError(predRefl(:) + predFl(:),measVals(:),'');
    end
    
end

% The opitmization is done, let's compute the estimates
rfCoeffs = reWEst;
reflEst = basisRefl*rfCoeffs;

dMat = max(tril(basisEm*wEst*basisEx',-1),0);
[U,S,V] = svd(dMat);

emEst = U(:,1)*sqrt(S(1,1))*sign(min(U(:,1)));
exEst = V(:,1)*sqrt(S(1,1))*sign(min(V(:,1)));

nF = max(exEst);
exEst = exEst/nF;
emEst = emEst*nF;

emCoeffs = basisEm'*emEst;
exCoeffs = basisEx'*exEst;

% Compute the prediction
predRefl = cameraGain.*(cameraMat'*diag(reflEst)*illuminant);
predFl = cameraGain.*(cameraMat'*dMat*illuminant);


end

%% Helper functions used to generate matrices or compute linear operations

function [A, b, AtA, Atb] = generateLSMatrices(measVals, camera, illuminant, cameraGain, basisRefl, basisEx, basisEm, alpha, beta, gamma)

nWaves = size(camera,1);
nFilters = size(camera,2);
nChannels = size(illuminant,2);
nReflBasis = size(basisRefl,2);

nExBasis = size(basisEx,2);
nEmBasis = size(basisEm,2);

b = measVals(:);

% Submatrix representing reflectance
rSubMat = zeros(nChannels*nFilters,nReflBasis);
for i=1:nReflBasis
    tmp = cameraGain.*(camera'*diag(basisRefl(:,i))*illuminant);
    rSubMat(:,i) = tmp(:);
end

% Submatrix representing fluorescence
cmEm = camera'*basisEm;
exLi = basisEx'*illuminant;
flSubMat = kron(exLi',cmEm);
flSubMat = diag(cameraGain(:))*flSubMat;

% Now include smoothness penalties
R = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];
rPen = R*basisRefl;

emPen = kron(basisEx,(R*basisEm));
exPen = kron(R*basisEx,basisEm);

A = [rSubMat, flSubMat;
     sqrt(alpha)*rPen, zeros(size(R,1),nExBasis*nEmBasis);
     zeros(size(emPen,1),nReflBasis) sqrt(beta)*emPen;
     zeros(size(exPen,1),nReflBasis) sqrt(gamma)*exPen];
b = [b; zeros(size(R,1) + size(emPen,1) + size(exPen,1),1)];

AtA = A'*A;
Atb = A'*b;

end


function [A, b, AtA, Atb] = generatePenaltyMatrices(reflPen,flMatPen,flWghtsPen,basisRefl,basisEx,basisEm)

A = [];
b = [];

nWaves = size(basisRefl,1);
nReflBasis = size(basisRefl,2);
nExBasis = size(basisEx,2);
nEmBasis = size(basisEm,2);

% Excitation-emission matrix
tmp = kron(basisEx,basisEm);

AtA = [basisRefl'*basisRefl zeros(nReflBasis,nExBasis*nEmBasis);
       zeros(nExBasis*nEmBasis,nReflBasis) tmp'*tmp + eye(nExBasis*nEmBasis)];
Atb = [basisRefl'*reflPen(:);
       tmp'*flMatPen(:) + flWghtsPen(:)];

end



function res = applyAtA(x, cameraMat, illuminant, cameraGain, reflBasis, exBasis, emBasis, rho, alpha, beta, gamma)

    nWaves = size(cameraMat,1);

    nReflBasis = size(reflBasis,2);
    nExBasis = size(exBasis,2);
    nEmBasis = size(emBasis,2);
    
    nabla = [eye(nWaves-1) zeros(nWaves-1,1)] - [zeros(nWaves-1,1) eye(nWaves-1)];

    cameraGain2 = cameraGain.^2;
    
    wr = x(1:nReflBasis);
    W = x(nReflBasis+1:end);
    W = reshape(W,nEmBasis,nExBasis);
    
    tmp1 = reflBasis'*diag(cameraMat*(cameraGain2.*(cameraMat'*diag(reflBasis*wr)*illuminant)*illuminant'));
    tmp2 = reflBasis'*diag(cameraMat*(cameraGain2.*(cameraMat'*tril(emBasis*W*exBasis',-1)*illuminant)*illuminant'));
    
    % OLD line read: res1 = 2*tmp1 + ...
    % I think the new one below is correct...
    res1 = tmp1 + alpha*reflBasis'*(nabla'*nabla)*reflBasis*wr + rho*(reflBasis'*reflBasis)*wr + tmp2;
    
    tmp3 = emBasis'*tril(cameraMat*(cameraGain2.*(cameraMat'*diag(reflBasis*wr)*illuminant)*illuminant'),-1)*exBasis;
    tmp4 = emBasis'*tril(cameraMat*(cameraGain2.*(cameraMat'*tril(emBasis*W*exBasis',-1)*illuminant)*illuminant'),-1)*exBasis;
                   
    res2 = tmp3 + tmp4 + beta*emBasis'*tril((nabla'*nabla)*tril(emBasis*W*exBasis',-1),-1)*exBasis + ...
                       + gamma*emBasis'*tril(tril(emBasis*W*exBasis',-1)*(nabla'*nabla),-1)*exBasis + ...
                       + rho*(emBasis'*emBasis)*W*(exBasis'*exBasis) + ...
                       + rho*W;
                   
    res = [res1(:); res2(:)];
    
end

function res = applyAtb(measVals, cameraMat, illuminant, cameraGain, reflBasis, exBasis, emBasis, y1, y2, y3, rho)
    %measVals is M in the paper and is the measured value of a pixel under
    %the ith filter, jth illuminant, and kth pixel

    %cameraMat is C in the paper and is the spectral response of the
    %filters scaled by the quantum efficiency of the sensor. [c1,c2,..ci]
    %where c are column vectors

    %illuminant is matrix L in paper [l1,l2,...,lj] where l are column
    %vectors

    %cameraGain is G in the paper and is the camera gain for the ith
    %filter and jth illuminant and the kth measurement in this case.

    %reflBasis Br in the paper. Is the reflectance basis functions

    %exBasis Bx in the paper. Is the excitation basis functions

    %emBasis Bm in the paper. Is the emission basis functions. Columns are
    %the basis functions of this matrix

    %Descriptions of y1,y2,y3, and rho come from the appendix

    %diag function takes diagonal 
    res1 = reflBasis'*diag(cameraMat*(cameraGain.*measVals)*illuminant') + rho*reflBasis'*y1;
    
    res2 = emBasis'*tril(cameraMat*(cameraGain.*measVals)*illuminant',-1)*exBasis + ...
           rho*emBasis'*y2*exBasis + rho*y3;
       
    res = [res1(:); res2(:)];

end
