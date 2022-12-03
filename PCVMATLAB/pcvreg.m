% Compute matrix with pseudo-validation set for Regression model
%
% Arguments
% ---------
% X         matrix with predictors from calibration set
% Y         matrix with responses from calibration set
% nComp     number of components for decomposition
% Center    logical, mean center data columns or not
% Scale     logical, standardize data columns or not
% CV        defines cross-validation splits pattern (see 'help crossval')
%
% !!! This is a generic function, use 'pcvpcr()' or 'pcvpls()' instead !!!
%
function [Xpv, D] = pcvreg(X, Y, nComp, Center, Scale, CV, funlist)

   if size(Y, 1) < size(Y, 2)
       Y = Y';
   end

   nRows = size(X, 1);
   nPred = size(X, 2);
   nResp = size(Y, 2);

   % default values
   if nargin < 6
       CV = {"ven", 4};
   end

   if nargin < 5
       Scale = false;
   end

   if nargin < 4
       Center = true;
   end

   % compute indices for cross-validation
   ind = crossval(CV, nRows, Y);

   % compute number of segments
   nSeg = max(ind);

   maxNComp = min([nRows - round(nRows / nSeg) - 1, nPred, 30]);

   if nargin < 2 || nComp > maxNComp
      nComp = maxNComp;
   end

   % compute and save mean and standard deviation
   if Center
      mX = mean(X);
      mY = mean(Y);
   else
      mX = zeros(1, nPred);
      mY = zeros(1, nResp);
   end

   if Scale
      sX = std(X);
      sY = std(Y);
   else
      sX = ones(1, nPred);
      sY = ones(1, nResp);
   end

   % autoscale the calibration set
   X = X - mX;
   Y = Y - mY;

   X = X ./ sX;
   Y = Y ./ sY;

   % prepare empty matrix for pseudo-validation set and scaling
   % coefficients
   Xpv = zeros(nRows, nPred);
   D = zeros(nSeg, nComp);

   % compute global model
   m = funlist.getglobalmodel(X, Y, nComp);

   % loop for computing the PV set
   for k = 1:nSeg

      % split data to calibration and validation
      indc = ind ~= k;
      indk = ind == k;

      Xk = X(indk, :);
      Xc = X(indc, :);
      Yc = Y(indc, :);

      % compute local model for current segment
      mk = funlist.getlocalmodel(Xc, Yc, m);

      % compute explained part of PV-set for current segment
      [Xpvhat, Dk] = funlist.getxpv(m, mk, Xk);
      D(k, :) = diag(Dk);

      % compute the orthogonal part of PV-set
      Xpvorth = getxpvorth(funlist.getqk(Xk, mk), Xk, m.Pi);

      % add the orthogonal part
      Xpv(indk, :) = Xpvhat + Xpvorth;
   end

   % uscenter and unscale the data
   Xpv = Xpv .* sX;
   Xpv = Xpv + mX;
end

