%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for PCR version      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = pcvpcrTest
   tests = functiontests(localfunctions);
end

function testRandomData(testCase)
% run tests for random dataset

   X = randn(100, 50);
   Y = X * rand(50, 1);
   runAll(testCase, X, Y);
end

function testCornData(testCase)
% runs tests for Corn dataset

   d = load('corn.mat');
   X = d.X;
   Y = d.Y;

   % automatic tests
   runAll(testCase, X, Y);

   % compare with reference values from R package
   compareWithReferences(testCase, X, Y);

   % manual test
   A = 20;
   center = true;
   scale = false;
   [~, ind] = sort(Y(:, 1));
   cvInd = repmat(1:4, 1, size(X, 1)/4);
   cv = cvInd(ind);

   [Xpv, D] = pcvpcr(X, Y, A, center, scale, cv);
   [Ycv, RMSECV] = pcrcvglobal(X, Y, A, center, scale, cv);
   [Ypv, RMSEPV] = pcrpv(X, Y, Xpv, A, center, scale);

   verifyEqual(testCase, size(D), [4, A]);
   verifyEqual(testCase, Ycv, Ypv, 'RelTol', 10^6);
   verifyEqual(testCase, RMSECV, RMSEPV, 'RelTol', 10^6);
   verifyEqual(testCase, RMSEPV, ...
      [0.363, 0.331, 0.298, 0.298, 0.271, 0.251, 0.213, 0.213, 0.214, 0.216, 0.216, 0.209,...
      0.190, 0.181, 0.183, 0.187, 0.187, 0.183, 0.184, 0.193], 'RelTol', 10^-2);   
end

function testCornDataPlots(testCase)
% makes visual tests for Corn data

   d = load('corn.mat');
   X = d.X;
   Y = d.Y;

   % distance plot    
   function predictionPlot(a, Y, Yp, Ycv, Ypv)
      
      y = Y(:, 1);
      yp = Yp(:, a);
      ypv = Ypv(:, a);
      ycv = Ycv(:, a);

      hold on
      scatter(y, yp, 'db')
      scatter(y, ycv, 'or')
      scatter(y, ypv, 'xk')
      hold off
      xlabel('Reference values, y')
      ylabel('Predicted values, yp')
   end   

   A = 30;
   center = true;
   scale = true;
   cv = {'ven', 4};

   [Yp, ~] = pcrpv(X, Y, X, A, center, scale);
   
   close all
   figure

   Xpv = pcvpcr(X, Y, A, center, scale, cv);
   [Ycv, ~] = pcrcvglobal(X, Y, A, center, scale, cv);
   [Ypv, ~] = pcrpv(X, Y, Xpv, A, center, scale);
   
   subplot(2, 2, 1)
   predictionPlot( 4, Y, Yp, Ycv, Ypv)
   title('Global, a = 4')
   
   subplot(2, 2, 2)
   predictionPlot(20, Y, Yp, Ycv, Ypv)
   title('Global, a = 20')

   Xpv = pcvpcr(X, Y, A, center, scale, cv, "local");
   [Ycv, ~] = pcrcvlocal(X, Y, A, center, scale, cv);
   [Ypv, ~] = pcrpv(X, Y, Xpv, A, center, scale);
   
   subplot(2, 2, 3)
   predictionPlot( 4, Y, Yp, Ycv, Ypv)
   title('Local, a = 4')

   subplot(2, 2, 4)
   predictionPlot(20, Y, Yp, Ycv, Ypv)
   title('Local, a = 20')

end

function runAll(testCase, X, Y)
% automatic tests for different combinations of PCV parameters

   cv_cases = {{'ven', 4}, {'ven', 10}, {'loo'}, {'rand', 10}};
   ncomp_cases = [1, 10, 20, 30];

   for cv = cv_cases
      for ncomp = ncomp_cases       
         for center = [true, false]
            for scale = [true, false]
               runOne(testCase, X, Y, ncomp, cv{1}, center, scale, "global");
            end
         end
      end
   end
   
   for cv = cv_cases
      for ncomp = ncomp_cases       
         for scale = [true, false]
            runOne(testCase, X, Y, ncomp, cv{1}, true, scale, "local");
         end
      end
   end
   
end


function compareWithReferences(testCase, X, Y)
% compare outcomes with reference values from R
   
   cv_cases = {{'ven', 4}, {'ven', 10}, {'loo'}};
   ncomp_cases = [1, 10, 20, 30];

   for cv = cv_cases
      for ncomp = ncomp_cases       
         for center = [true, false]
            for scale = [true, false]
               % rng() is not needed because all splits are systematic
               [Xpvg, Dg] = pcvpcr(X, Y, ncomp, center, scale, cv{1}, 'global');    
               [Xpvl, Dl] = pcvpcr(X, Y, ncomp, center, scale, cv{1}, 'local');    

               [Ypvg, ~] = pcrpv(X, Y, Xpvg, ncomp, center, scale);
               [Ypvl, ~] = pcrpv(X, Y, Xpvl, ncomp, center, scale);

               cvString = cv{1}{1};
               if numel(cv{1}) == 2
                  cvString = sprintf("%s%d", cvString, cv{1}{2});
               end
               
               caseDir = "../References/pcvpcr/";
               bString = {"FALSE", "TRUE"};
               fileSuffix = sprintf("-%d-%s-%s-%s.csv", ncomp, bString{center + 1}, bString{scale + 1}, cvString);
               
               verifyEqual(testCase, Ypvg, csvread(strcat(caseDir, "Ypvg", fileSuffix)), 'RelTol', 10^-6)
               verifyEqual(testCase, Ypvl, csvread(strcat(caseDir, "Ypvl", fileSuffix)), 'RelTol', 10^-6)
               verifyEqual(testCase, Dg, csvread(strcat(caseDir, "Dg", fileSuffix)), 'RelTol', 10^-6)
               verifyEqual(testCase, Dl, csvread(strcat(caseDir, "Dl", fileSuffix)), 'RelTol', 10^-6)
            end
         end
      end
   end
end

function runOne(testCase, X, Y, ncomp, cv, center, scale, scope)
% run tests for a given combination of parameters

   if scope == 'local'
      fcv = @pcrcvlocal;
   else
      fcv = @pcrcvglobal;
   end

   rng(42)
   Xpv = pcvpcr(X, Y, ncomp, center, scale, cv, scope);    

   rng(42)
   [Ycv, RMSECV] = fcv(X, Y, ncomp, center, scale, cv);       
   [Ypv, RMSEPV] = pcrpv(X, Y, Xpv, ncomp, center, scale);
      
   if scope == "global"
      RelTol = 10^-6;
      verifyEqual(testCase, Ypv, Ycv, 'RelTol', RelTol);
      verifyEqual(testCase, RMSECV, RMSEPV, 'RelTol', RelTol);   
   else
      RelTol = 0.05;
      verifyEqual(testCase, RMSECV, RMSEPV, 'RelTol', RelTol);   
   end
   
end

function [Ypv, RMSEPV] = pcrpv(X, Y, Xpv, ncomp, center, scale)
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
   [~, ~, P] = svds(Xs, ncomp);
   T = Xs * P;

   % predictions for PV-set
   Tpv = Xpvs * P;
   Ypv = zeros(nobj, ncomp);
   for a = 1:ncomp
      aind = 1:a;
      Ta = T(:, aind);
      Tpva = Tpv(:, aind);
      Ca = (Ta' * Ta)' \ (Ta' * Ys);
      Ypv(:, a) = Tpva * Ca;
   end
   
   Ypv = Ypv * sY + mY;
   RMSEPV = sqrt(mean((Ypv - Y).^2));
end

function [Ycv, RMSECV] = pcrcvglobal(X, Y, ncomp, center, scale, cv) 
% PCRCVGLOBAL computes matrices with predicted response values (Ycv) and
% root mean squared error (RMSECV) using cross-validation with global centering and scaling.
% The input parameters are the same as for PCRPCV function.
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
      [~, ~, Pk] = svds(Xc, ncomp);
      
      % make predictions for local validation set
      Tc = Xc * Pk;
      Tk = Xk * Pk;

      % compute distances
      for a = 1:ncomp

         aind  = 1:a;
         Tka = Tk(:, aind);
         Tca = Tc(:, aind);
         Ca = (Tca' * Tca)' \ (Tca' * Yc);
         Ycv(indk, a) = Tka * Ca;
      end
   end
   
   Ycv = Ycv * sY + mY;
   RMSECV = sqrt(mean((Ycv - Y).^2));
end

function [Ycv, RMSECV] = pcrcvlocal(X, Y, ncomp, center, scale, cv) 
% PCACVLOCAL computes matrices with orthogonal distances (Qcv) and score 
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
      [~, ~, Pk] = svds(Xc, ncomp);
      
      % make predictions for local validation set
      Tk = Xk * Pk;
      Tc = Xc * Pk;

      % compute and save the distances
      for a = 1:ncomp
         aind  = 1:a;
         Tka = Tk(:, aind);
         Tca = Tc(:, aind);
         Ca = (Tca' * Tca)' \ (Tca' * Yc);
         Ycv(indk, a) = Tka * Ca;
      end
      
      Ycv(indk, :) = Ycv(indk, :) .* sY + mY;
   end
   
   RMSECV = sqrt(mean((Ycv - Y).^2));
end