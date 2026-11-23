%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for PCA version      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = pcvpcaTest
   tests = functiontests(localfunctions);
end


function testEdgeCases(testCase)
% edge cases: constant column, small dataset, invalid CVScope

   X = reshape([1 3 17 19 10 14 7 13 9 11 2 15 8 6 4 12 22 18 4 14 13 12 1 15], 8, 3);

   % constant column must not produce NaN values
   Xc = [X, repmat(42, size(X, 1), 1)];
   rng(42);
   Xpv = pcvpca(Xc, 2, true, true, {"ven", 2});
   verifySize(testCase, Xpv, size(Xc));
   verifyFalse(testCase, any(isnan(Xpv(:))));

   % small dataset (6 rows): number of components is reduced automatically
   rng(42);
   Xpv = pcvpca(X(1:6, :), 10, true, false, {"ven", 3});
   verifyEqual(testCase, size(Xpv, 1), 6);
   verifyFalse(testCase, any(isnan(Xpv(:))));

   % wrong CVScope value must raise an error
   verifyError(testCase, @() pcvpca(X, 2, true, false, {"ven", 4}, "wrong"), ?MException);

   % uneven split: the largest segment leaves 4 rows, so at most 3 components are possible
   rng(1);
   Xu = randn(20, 10);
   cvu = [ones(16, 1); 2 * ones(4, 1)];
   rng(2);
   Xpv9 = pcvpca(Xu, 9, true, false, cvu);
   rng(2);
   Xpv3 = pcvpca(Xu, 3, true, false, cvu);
   verifyTrue(testCase, all(isfinite(Xpv9(:))));
   verifyEqual(testCase, Xpv9, Xpv3);
end

function testCornData(testCase)
% runs tests for Corn dataset

   d = load('corn.mat');
   X = d.X;

   % compare with reference values from R package
   compareWithReferences(testCase, X);

   % manual test
   A = 30;
   center = true;
   scale = false;
   cv = {'ven', 4};
   
   Xpv = pcvpca(X, A, center, scale, cv);
   [Qpv, Hpv] = pcapv(X, Xpv, A, center, scale);
   [Q, H]     = pcapv(X, X, A, center, scale);

   a = 2;   
   
   q = Q(:, a);
   qpv = Qpv(:, a);
   h = H(:, a);
   hpv = Hpv(:, a);
   
   verifyEqual(testCase, sum(q/mean(q) > 4), 4)
   verifyEqual(testCase, sum(qpv/mean(q) > 4), 4)
   verifyEqual(testCase, sum(h/mean(h) > 5), 2)
   verifyEqual(testCase, sum(hpv/mean(h) > 5), 2)

   
   a = 20;   
   
   q = Q(:, a);
   qpv = Qpv(:, a);
   h = H(:, a);
   hpv = Hpv(:, a);
   
   verifyEqual(testCase, sum(q/mean(q) > 4), 0)
   verifyEqual(testCase, sum(qpv/mean(q) > 4), 5)
   verifyEqual(testCase, sum(h/mean(h) > 2), 3)
   verifyEqual(testCase, sum(hpv/mean(h) > 2), 9)
   
end

function testCornDataPlots(testCase)
% makes visual tests for Corn data

   d = load('corn.mat');
   X = d.X;

   % distance plot    
   function distancePlot(a, Q, Qpv, H, Hpv)
      
      q = Q(:, a);
      qpv = Qpv(:, a);

      h = H(:, a);
      hpv = Hpv(:, a);

      hold on
      scatter(h / mean(h), q / mean(q), 'db')
      scatter(hpv / mean(h), qpv / mean(q), 'xr')
      hold off
      xlabel('Score distance, h/h0')
      ylabel('Orthogonal distance, q/q0')
   end   

   A = 30;
   center = true;
   scale = true;
   cv = {'ven', 20};

   [Q, H] = pcapv(X, X, A, center, scale);
   
   close all
   figure

   Xpv = pcvpca(X, A, center, scale, cv);
   [Qpv, Hpv] = pcapv(X, Xpv, A, center, scale);
   
   subplot(2, 2, 1)
   distancePlot( 4, Q, Qpv, H, Hpv)
   title('Global, a = 4')
   
   subplot(2, 2, 2)
   distancePlot(20, Q, Qpv, H, Hpv)
   title('Global, a = 20')

   Xpv = pcvpca(X, A, center, scale, cv, "local");
   [Qpv, Hpv] = pcapv(X, Xpv, A, center, scale);
   
   subplot(2, 2, 3)
   distancePlot( 4, Q, Qpv, H, Hpv)
   title('Local, a = 4')

   subplot(2, 2, 4)
   distancePlot(20, Q, Qpv, H, Hpv)
   title('Local, a = 20')

end

function compareWithReferences(testCase, X)
% compare outcomes with reference values from R
   
   cv_cases = {{'ven', 4}, {'ven', 10}, {'loo'}};
   ncomp_cases = [1, 10, 20, 30];
   center = true;
   
   for cv = cv_cases
      for ncomp = ncomp_cases       
         for scale = [true, false]
            % rng() is not needed because all splits are systematic
            Xpvg = pcvpca(X, ncomp, center, scale, cv{1}, 'global');    
            Xpvl = pcvpca(X, ncomp, center, scale, cv{1}, 'local');    

            [Qpvg, Hpvg] = pcapv(X, Xpvg, ncomp, center, scale);
            [Qpvl, Hpvl] = pcapv(X, Xpvl, ncomp, center, scale);

            cvString = cv{1}{1};
            if numel(cv{1}) == 2
               cvString = sprintf("%s%d", cvString, cv{1}{2});
            end
               
            caseDir = "../.tests/pcvpca/";
            bString = {"FALSE", "TRUE"};
            fileSuffix = sprintf("-%d-%s-%s.csv", ncomp, bString{scale + 1}, cvString);
               
            verifyEqual(testCase, Qpvg, csvread(strcat(caseDir, "Qpvg", fileSuffix)), 'RelTol', 10^-6)
            verifyEqual(testCase, Hpvg, csvread(strcat(caseDir, "Hpvg", fileSuffix)), 'RelTol', 10^-6)
            verifyEqual(testCase, Qpvl, csvread(strcat(caseDir, "Qpvl", fileSuffix)), 'RelTol', 10^-6)
            verifyEqual(testCase, Hpvl, csvread(strcat(caseDir, "Hpvl", fileSuffix)), 'RelTol', 10^-6)
                        
         end
      end
   end
end

function [Qpv, Hpv] = pcapv(X, Xpv, ncomp, center, scale)
%PCAPCV fit a global PCA model for X dataset and then applies it to Xpv
% data and computes orthogonal (Qpv) and score (Hpv) distances.
%
   nvar = size(X, 2);
   nobj = size(X, 1);
   
   % center and scale 
   if center, mX = mean(X); else, mX = zeros(1, nvar); end
   if scale, sX = std(X); else, sX = ones(1, nvar); end

   Xs = (X - mX) ./ sX;
   Xpvs = (Xpv - mX) ./ sX;
   
   % global model
   [~, S, P] = svds(Xs, ncomp);
   s = diag(S) / sqrt(nobj - 1);

   % predictions for PV-set
   Tpv = Xpvs * P;
   Upv = Tpv * diag(1 ./ s);

   Qpv = zeros(nobj, ncomp);
   Hpv = zeros(nobj, ncomp);

   for a = 1:ncomp

      aind = 1:a;
      Pa = P(:, aind);
      Tpva = Tpv(:, aind);
      Upva = Upv(:, aind);

      Epva = Xpvs - Tpva * Pa';
      Qpv(:, a) = sum(Epva.^2, 2);
      Hpv(:, a) = sum(Upva.^2, 2);
   end
   
end
