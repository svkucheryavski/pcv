function Xpv = pcvpca(X, nComp, Center, Scale, CV, CVScope)
%PCVPCA Compute matrix with pseudo-validation set for PCA/SIMCA model
%
% Arguments:
% ----------
% X         matrix with calibration set
% nComp     number of components for PCA decomposition
% Center    logical, mean center data columns or not
% Scale     logical, standardize data columns or not
% CV        defines cross-validation splits pattern (see description)
% CVScope   defines if centering/scaling must be done globally ('global')
% or locally ('local'). Default value is 'global'.
%
% Returns:
% --------
% Xpv       matrix with PV-set (same size as X)
%
% The method computes pseudo-validation matrix Xpv, based on PCA decomposition of
% calibration set X and cross-validation. The 'CV' parameter can be:
%
% A number:
% =========
% 'CV = 1' - full cross-validation (leave one out)
% 'CV = n' - random cross-validation with n segments
%
% A cell array:
% =============
% 'CV = {"rand", n}' - random cross-validation with n segments
% 'CV = {"ven", n}' - systematic cross-validation with n segments
% 'CV = {"loo"}' - full cross-validation (leave one out)
%
% A vector:
% =========
% 'CV = [1, 2, 3, 1, 2, 3] - cross-validation of 6 rows with 3 segments
%
% See more details here: https://github.com/svkucheryavski/pcv
%
% See also: PCVPCR, PCVPLS

   nRows = size(X, 1);
   nCols = size(X, 2);

   % default values
   if nargin < 6
      CVScope = "global";
   end
   
   if nargin < 5
       CV = {"ven", 4};
   end

   if nargin < 4
       Scale = false;
   end

   if nargin < 3
       Center = true;
   end

   % compute indices for cross-validation
   ind = crossval(CV, nRows);

   % compute number of segments
   nSeg = max(ind);

   maxNComp = min([nRows - round(nRows / nSeg) - 1, nCols, 30]);

   if nargin < 2 || nComp > maxNComp
      nComp = maxNComp;
   end

   % compute and save mean and standard deviation
   if Center
      mX = mean(X);
   else
      mX = zeros(1, nCols);
   end

   if Scale
      sX = std(X);
   else
      sX = ones(1, nCols);
   end

   % autoscale data globally
   X = (X - mX) ./ sX;
   
   % create a global model
   [~, ~, P] = svds(X, nComp);
   Pi = eye(nCols) - P * P';

   % prepare empty matrix for pseudo-validation set
   Xpv = zeros(nRows, nCols);

   % cv-loop
   for k = 1:nSeg

      % split data to calibration and validation
      indc = ind ~= k;
      indk = ind == k;

      Xk = X(indk, :);
      Xc = X(indc, :);
      
      % compute mean and standard deviation and autoscale locally in case of local scope
      if CVScope == "local"
         if Center, mXl = mean(Xc); else, mXl = zeros(1, nCols); end
         if Scale, sXl = std(Xc); else, sXl = ones(1, nCols); end
         Xc = (Xc - mXl) ./ sXl;
         Xk = (Xk - mXl) ./ sXl;
      end
      
      % get loadings for local model and rotation matrix between global and local models
      [~, ~, Pk] = svds(Xc, nComp);

      % correct direction of loadings for local model
      a = acos(sum(P .* Pk)) < pi / 2;
      Pk = Pk * diag(a * 2 - 1);

      % compute explained part of Xpv
      Tk = Xk * Pk;
      Xpvpar = Tk * P';
            
      % compute orthogonal part of Xpv
      Ek = Xk - Tk * Pk';
      qk = sum(Ek.^2, 2);
      Xpvorth = getxpvorth(qk, Xk, Pi);

      % rotate the local validation set and save as a part of Xpv
      Xpv(indk, :) = Xpvorth + Xpvpar;      
   end

   % uscenter and unscale the data
   Xpv = Xpv .* sX + mX;
end

