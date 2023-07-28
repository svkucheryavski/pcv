%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for PLS version      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = pcvplsTest
   tests = functiontests(localfunctions);
end

function testCornData(testCase)
% runs tests for Corn dataset

   d = load('corn.mat');
   X = d.X;
   Y = d.Y;

   % compare with reference values from R package
   compareWithReferences(testCase, X, Y);
end

function testCornDataPlots(testCase)
% makes visual tests for Corn data

   d = load('corn.mat');
   X = d.X;
   Y = d.Y;

   % distance plot    
   function predictionPlot(a, Y, Yp, Ypv)
      
      y = Y(:, 1);
      yp = Yp(:, a);
      ypv = Ypv(:, a);

      hold on
      scatter(y, yp, 'db')
      scatter(y, ypv, 'xr')
      hold off
      xlabel('Reference values, y')
      ylabel('Predicted values, yp')
   end   

   A = 30;
   center = true;
   scale = true;
   cv = {'ven', 4};

   [Yp, ~] = plspv(X, Y, X, A, center, scale);
   
   close all
   figure

   Xpv = pcvpls(X, Y, A, center, scale, cv);
   [Ypv, ~] = plspv(X, Y, Xpv, A, center, scale);
   
   subplot(2, 2, 1)
   predictionPlot( 4, Y, Yp, Ypv)
   title('Global, a = 4')
   
   subplot(2, 2, 2)
   predictionPlot(20, Y, Yp, Ypv)
   title('Global, a = 20')

   Xpv = pcvpls(X, Y, A, center, scale, cv, "local");
   [Ypv, ~] = plspv(X, Y, Xpv, A, center, scale);
   
   subplot(2, 2, 3)
   predictionPlot( 4, Y, Yp, Ypv)
   title('Local, a = 4')

   subplot(2, 2, 4)
   predictionPlot(20, Y, Yp, Ypv)
   title('Local, a = 20')

end

function compareWithReferences(testCase, X, Y)
% compare outcomes with reference values from R
   
   cv_cases = {{'ven', 4}, {'ven', 10}, {'loo'}};
   ncomp_cases = [1, 10, 20, 30];
   center = true;

   for cv = cv_cases
      for ncomp = ncomp_cases       
         for scale = [true, false]
            % rng() is not needed because all splits are systematic
            [Xpvg, Dg] = pcvpls(X, Y, ncomp, center, scale, cv{1}, 'global');    
            [Xpvl, Dl] = pcvpls(X, Y, ncomp, center, scale, cv{1}, 'local');    

            [Ypvg, ~] = plspv(X, Y, Xpvg, ncomp, center, scale);
            [Ypvl, ~] = plspv(X, Y, Xpvl, ncomp, center, scale);

            cvString = cv{1}{1};
            if numel(cv{1}) == 2
               cvString = sprintf("%s%d", cvString, cv{1}{2});
            end

            caseDir = "../Test outcomes/pcvpls/";
            bString = {"FALSE", "TRUE"};
            fileSuffix = sprintf("-%d-%s-%s.csv", ncomp, bString{scale + 1}, cvString);

            verifyEqual(testCase, Ypvg, csvread(strcat(caseDir, "Ypvg", fileSuffix)), 'RelTol', 10^-6)
            verifyEqual(testCase, Ypvl, csvread(strcat(caseDir, "Ypvl", fileSuffix)), 'RelTol', 10^-6)
            verifyEqual(testCase, Dg, csvread(strcat(caseDir, "Dg", fileSuffix)), 'RelTol', 10^-6)
            verifyEqual(testCase, Dl, csvread(strcat(caseDir, "Dl", fileSuffix)), 'RelTol', 10^-6)
         end
      end
   end
end

function runOne(testCase, X, Y, ncomp, cv, center, scale, scope)
% run tests for a given combination of parameters

   if scope == 'local'
      fcv = @plscvlocal;
   else
      fcv = @plscvglobal;
   end

   rng(42)
   Xpv = pcvpls(X, Y, ncomp, center, scale, cv, scope);    

   rng(42)
   [Ycv, RMSECV] = fcv(X, Y, ncomp, center, scale, cv);       
   [Ypv, RMSEPV] = plspv(X, Y, Xpv, ncomp, center, scale);
      
   if scope == "global"
      RelTol = 10^-5;
      verifyEqual(testCase, Ypv, Ycv, 'RelTol', RelTol);
      verifyEqual(testCase, RMSECV, RMSEPV, 'RelTol', RelTol);   
   else
      RelTol = 0.05;
      verifyEqual(testCase, RMSECV, RMSEPV, 'RelTol', RelTol);   
   end
   
end

function [Ypv, RMSEPV] = plspv(X, Y, Xpv, ncomp, center, scale)
%PCAPCV fit a global PCA model for X dataset and then applies it to Xpv
% data and computes orthogonal (Qpv) and score (Hpv) distances.
%
   nvar = size(X, 2);
   nobj = size(X, 1);
   
   % center and scale 
   if center, mX = mean(X); else, mX = zeros(1, nvar); end
   if scale, sX = std(X); else, sX = ones(1, nvar); end
   if center, mY = mean(Y); else, mY = 0; end
   if scale, sY = std(Y); else, sY = 1; end

   Xs = (X - mX) ./ sX;
   Xpvs = (Xpv - mX) ./ sX;
   Ys = (Y - mY) ./ sY;
   
   % global model
   [~, R, C] = simpls(Xs, Ys, ncomp);

   % predictions for PV-set
   Tpv = Xpvs * R;
   Ypv = zeros(nobj, ncomp);
   for a = 1:ncomp
      aind = 1:a;
      Ca = C(:, aind);
      Tpva = Tpv(:, aind);
      Ypv(:, a) = Tpva * Ca';
   end
   
   Ypv = Ypv * sY + mY;
   RMSEPV = sqrt(mean((Ypv - Y).^2));
end

function [Ycv, RMSECV] = plscvglobal(X, Y, ncomp, center, scale, cv) 
% PLSCVGLOBAL computes matrices with predicted response values (Ycv) and
% root mean squared error (RMSECV) using cross-validation with global centering and scaling.
% The input parameters are the same as for PLSPCV function.
%
   
   nvar = size(X, 2);
   nobj = size(X, 1);
   
   % get vector with indices of samples in CV-segments   
   cvind = crossval(cv, nobj, Y);
   
   % scale data globally
   if center
      mX = mean(X); 
      mY = mean(Y); 
   else
      mX = zeros(1, nvar); 
      mY = 0; 
   end
   
   if scale
      sX = std(X); 
      sY = std(Y);
   else
      sX = ones(1, nvar); 
      sY = 1;
   end
   
   Xs = (X - mX) ./ sX;
   Ys = (Y - mY) ./ sY;

   % compute distances
   nseg = max(cvind);
   Ycv = zeros(nobj, ncomp);

   for k = 1:nseg
      
      % get indices for local sets
      indc = cvind ~= k;
      indk = cvind == k;
      
      % split data to local calibration and validation sets
      Xc = Xs(indc, :);
      Xk = Xs(indk, :);
      Yc = Ys(indc, :);

      % fit the local model
      [~, Rk, Ck] = simpls(Xc, Yc, ncomp);
      
      % make predictions for local validation set
      Tk = Xk * Rk;

      % compute distances
      for a = 1:ncomp
         aind  = 1:a;
         Tka = Tk(:, aind);
         Cka = Ck(:, aind);
         Ycv(indk, a) = Tka * Cka';
      end
   end
   
   Ycv = Ycv * sY + mY;
   RMSECV = sqrt(mean((Ycv - Y).^2));
end

function [Ycv, RMSECV] = plscvlocal(X, Y, ncomp, center, scale, cv) 
% PLSCVLOCAL computes matrices with orthogonal distances (Qcv) and score 
% distances (Hcv) using cross-validation with local centering and scaling.
% The input parameters are the same as for PVAPCV function.
%

   nvar = size(X, 2);
   nobj = size(X, 1);
   
   % get vector with indices of samples in CV-segments
   cvind = crossval(cv, nobj, Y);
   
   % compute distances
   nseg = max(cvind);
   Ycv = zeros(nobj, ncomp);
   for k = 1:nseg
      
      % get indices for local sets
      indc = cvind ~= k;
      indk = cvind == k;
      
      % split data to local calibration and validation sets
      Xc = X(indc, :);
      Xk = X(indk, :);
      Yc = Y(indc, :);

      % autoscale data locally
      if center
         mX = mean(Xc); 
         mY = mean(Yc); 
      else
         mX = zeros(1, nvar); 
         mY = 0; 
      end
   
      if scale
         sX = std(Xc); 
         sY = std(Yc);
      else
         sX = ones(1, nvar); 
         sY = 1;
      end
      
      Xc = (Xc - mX) ./ sX;
      Xk = (Xk - mX) ./ sX;
      Yc = (Yc - mY) ./ sY;
      
      % fit the local model
      [~, Rk, Ck] = simpls(Xc, Yc, ncomp);
      
      % make predictions for local validation set
      Tk = Xk * Rk;

      % compute and save the distances
      for a = 1:ncomp
         aind  = 1:a;
         Tka = Tk(:, aind);
         Cka = Ck(:, aind);
         Ycv(indk, a) = Tka * Cka';
      end
      
      Ycv(indk, :) = Ycv(indk, :) .* sY + mY;
   end
   
   RMSECV = sqrt(mean((Ycv - Y).^2));
end