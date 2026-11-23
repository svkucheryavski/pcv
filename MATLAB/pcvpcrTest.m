%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for PCR version      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = pcvpcrTest
   tests = functiontests(localfunctions);
end


function testEdgeCases(testCase)
% edge cases: constant column, small dataset, invalid CVScope

   Y = [32; 35; 36; 37; 42; 43; 43; 44];
   X = [150 41 28000 119; 160 48 31000 129; 166 47 28000 112; 166 49 14000 123; ...
        175 67 38000 105; 180 80 30000 129; 181 75 31000 105; 180 81 42000 113];

   % constant column must not produce NaN values
   Xc = [X, repmat(42, size(X, 1), 1)];
   rng(42);
   [Xpv, D] = pcvpcr(Xc, Y, 3, true, true, {"ven", 4});
   verifySize(testCase, Xpv, size(Xc));
   verifyFalse(testCase, any(isnan(Xpv(:))));
   verifyFalse(testCase, any(isnan(D(:))));

   % small dataset (6 rows): number of components is reduced automatically
   rng(42);
   [Xpv, ~] = pcvpcr(X(1:6, :), Y(1:6), 10, true, false, {"ven", 3});
   verifyEqual(testCase, size(Xpv, 1), 6);
   verifyFalse(testCase, any(isnan(Xpv(:))));

   % wrong CVScope value must raise an error
   verifyError(testCase, @() pcvpcr(X, Y, 3, true, false, {"ven", 4}, "wrong"), ?MException);

   % uneven split: the largest segment leaves 4 rows, so at most 3 components are possible
   rng(1);
   Xu = randn(20, 10);
   Yu = randn(20, 1);
   cvu = [ones(16, 1); 2 * ones(4, 1)];
   rng(2);
   [Xpv9, D9] = pcvpcr(Xu, Yu, 9, true, false, cvu);
   rng(2);
   Xpv3 = pcvpcr(Xu, Yu, 3, true, false, cvu);
   verifyTrue(testCase, all(isfinite(Xpv9(:))));
   verifyEqual(testCase, Xpv9, Xpv3);
   verifyEqual(testCase, size(D9, 2), 3);
end


function testCornData(testCase)
% runs tests for Corn dataset

   d = load('corn.mat');
   X = d.X;
   Y = d.Y;


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
   [~, RMSEPV] = pcrpv(X, Y, Xpv, A, center, scale);

   verifyEqual(testCase, size(D), [4, A]);
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

   [Yp, ~] = pcrpv(X, Y, X, A, center, scale);
   
   close all
   figure

   Xpv = pcvpcr(X, Y, A, center, scale, cv);
   [Ypv, ~] = pcrpv(X, Y, Xpv, A, center, scale);
   
   subplot(2, 2, 1)
   predictionPlot( 4, Y, Yp, Ypv)
   title('Global, a = 4')
   
   subplot(2, 2, 2)
   predictionPlot(20, Y, Yp, Ypv)
   title('Global, a = 20')

   Xpv = pcvpcr(X, Y, A, center, scale, cv, "local");
   [Ypv, ~] = pcrpv(X, Y, Xpv, A, center, scale);
   
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
            [Xpvg, Dg] = pcvpcr(X, Y, ncomp, center, scale, cv{1}, 'global');    
            [Xpvl, Dl] = pcvpcr(X, Y, ncomp, center, scale, cv{1}, 'local');    

            [Ypvg, ~] = pcrpv(X, Y, Xpvg, ncomp, center, scale);
            [Ypvl, ~] = pcrpv(X, Y, Xpvl, ncomp, center, scale);

            cvString = cv{1}{1};
            if numel(cv{1}) == 2
               cvString = sprintf("%s%d", cvString, cv{1}{2});
            end
               
            caseDir = "../.tests/pcvpcr/";
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

