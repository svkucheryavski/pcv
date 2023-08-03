function [Xpv, D] = pcvreg(X, Y, nComp, Center, Scale, CV, funlist, CVScope)
%PCVREG Compute matrix with pseudo-validation set for Regression model
%
% Arguments
% ---------
% X         matrix with predictors from calibration set
% Y         matrix with responses from calibration set
% nComp     number of components for decomposition
% Center    logical, mean center data columns or not
% Scale     logical, standardize data columns or not
% CV        defines cross-validation splits pattern (see 'help crossval')
% CVScope   defines if centering/scaling must be done globally ('global')
% or locally ('local'). Default value is 'global'.
%
% !!! This is a generic function, use 'pcvpcr()' or 'pcvpls()' instead !!!
%

   if size(Y, 1) < size(Y, 2)
       Y = Y';
   end

   nRows = size(X, 1);
   nPred = size(X, 2);
   nResp = size(Y, 2);
   

   % compute indices for cross-validation
   ind = crossval(CV, nRows, Y);
   
   % compute number of segments
   nSeg = max(ind);

   % adjust number of components if necessary
   maxNComp = min([nRows - round(nRows / nSeg) - 1, nPred]);
   if nComp > maxNComp
      nComp = maxNComp;
   end

   % compute and save mean and standard deviation
   if Center
      mXg = mean(X);
      mYg = mean(Y);
   else
      mXg = zeros(1, nPred);
      mYg = zeros(1, nResp);
   end

   if Scale
      sXg = std(X);
      sYg = std(Y);
   else
      sXg = ones(1, nPred);
      sYg = ones(1, nResp);
   end

   % autoscale the calibration set and compute global model
   X = (X - mXg) ./ sXg;
   Y = (Y - mYg) ./ sYg;
   m = funlist.getglobalmodel(X, Y, nComp);

   % prepare empty matrix for pseudo-validation set and scaling
   % coefficients
   Xpv = zeros(nRows, nPred);
   D = zeros(nSeg, nComp);
   
   % loop for computing the PV set
   for k = 1:nSeg

      % split data to calibration and validation
      indc = ind ~= k;
      indk = ind == k;
      
      Xk = X(indk, :);
      Xc = X(indc, :);
      Yc = Y(indc, :);
      
      % autoscale data in case of local scope
      if CVScope == "local"
         
         if Center
            mX = mean(Xc);
            mY = mean(Yc);
         else
            mX = zeros(1, nPred);
            mY = zeros(1, nResp);
         end

         if Scale
            sX = std(Xc);
            sY = std(Yc);
         else
            sX = ones(1, nPred);
            sY = ones(1, nResp);
         end
         
         Xc = (Xc - mX) ./ sX;
         Yc = (Yc - mY) ./ sY;
         Xk = (Xk - mX) ./ sX;
      end
            
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

   % uscenter and unscale the data for global scope
   Xpv = Xpv .* sXg + mXg;
end

