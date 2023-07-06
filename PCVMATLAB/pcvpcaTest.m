%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for PCA version      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = pcvpcaTest
   tests = functiontests(localfunctions);
end

function testRandomData(testCase)
% run tests for random dataset

   X = randn(100, 50);
   runAll(testCase, X);
end

function testCornData(testCase)
% runs tests for Corn dataset

   d = load('corn.mat');
   X = d.X;

   % automatic tests
   runAll(testCase, X);

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
   function distancePlot(a, Q, Qcv, Qpv, H, Hcv, Hpv)
      
      q = Q(:, a);
      qpv = Qpv(:, a);
      qcv = Qcv(:, a);

      h = H(:, a);
      hpv = Hpv(:, a);
      hcv = Hcv(:, a);

      hold on
      scatter(h / mean(h), q / mean(q), 'db')
      scatter(hcv / mean(h), qcv / mean(q), 'or')
      scatter(hpv / mean(h), qpv / mean(q), 'xk')
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
   [Qcv, Hcv] = pcacvglobal(X, A, center, scale, cv);
   [Qpv, Hpv] = pcapv(X, Xpv, A, center, scale);
   
   subplot(2, 2, 1)
   distancePlot( 4, Q, Qcv, Qpv, H, Hcv, Hpv)
   title('Global, a = 4')
   
   subplot(2, 2, 2)
   distancePlot(20, Q, Qcv, Qpv, H, Hcv, Hpv)
   title('Global, a = 20')

   Xpv = pcvpca(X, A, center, scale, cv, "local");
   [Qcv, Hcv] = pcacvlocal(X, A, center, scale, cv);
   [Qpv, Hpv] = pcapv(X, Xpv, A, center, scale);
   
   subplot(2, 2, 3)
   distancePlot( 4, Q, Qcv, Qpv, H, Hcv, Hpv)
   title('Local, a = 4')

   subplot(2, 2, 4)
   distancePlot(20, Q, Qcv, Qpv, H, Hcv, Hpv)
   title('Local, a = 20')

end

function runAll(testCase, X)
% automatic tests for different combinations of PCV parameters

   cv_cases = {{'ven', 4}, {'ven', 10}, {'loo'}, {'rand', 10}};
   ncomp_cases = [1, 10, 20, 30];

   for cv = cv_cases
      for ncomp = ncomp_cases       
         for center = [true, false]
            for scale = [true, false]
               for scope = ['global', 'local']         
                  runOne(testCase, X, ncomp, cv{1}, center, scale, scope);
               end
            end
         end
      end
   end
end


function compareWithReferences(testCase, X)
% compare outcomes with reference values from R
   
   cv_cases = {{'ven', 4}, {'ven', 10}, {'loo'}};
   ncomp_cases = [1, 10, 20, 30];

   for cv = cv_cases
      for ncomp = ncomp_cases       
         for center = [true, false]
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
               
               caseDir = "../References/pcvpca/";
               bString = {"FALSE", "TRUE"};
               fileSuffix = sprintf("-%d-%s-%s-%s.csv", ncomp, bString{center + 1}, bString{scale + 1}, cvString);
               
               verifyEqual(testCase, Qpvg, csvread(strcat(caseDir, "Qpvg", fileSuffix)), 'RelTol', 10^-6)
               verifyEqual(testCase, Hpvg, csvread(strcat(caseDir, "Hpvg", fileSuffix)), 'RelTol', 10^-6)
               verifyEqual(testCase, Qpvl, csvread(strcat(caseDir, "Qpvl", fileSuffix)), 'RelTol', 10^-6)
               verifyEqual(testCase, Hpvl, csvread(strcat(caseDir, "Hpvl", fileSuffix)), 'RelTol', 10^-6)
            end
         end
      end
   end
end

function runOne(testCase, X, ncomp, cv, center, scale, scope)
% run tests for a given combination of parameters
   
   if scope == 'local'
      fcv = @pcacvlocal;
   else
      fcv = @pcacvglobal;
   end

   rng(42)
   Xpv = pcvpca(X, ncomp, center, scale, cv, scope);    

   rng(42)
   [Qcv, Hcv] = fcv(X, ncomp, center, scale, cv);       
   [Qpv, Hpv] = pcapv(X, Xpv, ncomp, center, scale);

   verifyEqual(testCase, Qpv, Qcv, 'RelTol', 10^-6);
   verifyEqual(testCase, Hpv, Hcv, 'RelTol', 10^-6);   
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

function [Qcv, Hcv] = pcacvglobal(X, ncomp, center, scale, cv) 
% PCACVGLOBAL computes matrices with orthogonal distances (Qcv) and score 
% distances (Hcv) using cross-validation with global centering and scaling.
% The input parameters are the same as for PVAPCV function.
%
   
   nvar = size(X, 2);
   nobj = size(X, 1);
   
   % get vector with indices of samples in CV-segments   
   cvind = crossval(cv, nobj);
   
   % scale data globally
   if center, mX = mean(X); else, mX = zeros(1, nvar); end
   if scale, sX = std(X); else, sX = ones(1, nvar); end
   Xs = (X - mX) ./ sX;

   % fit global model
   [~, S, ~] = svds(Xs, ncomp);
   s = diag(S) / sqrt(nobj - 1);

   % compute distances
   nseg = max(cvind);
   Qcv = zeros(nobj, ncomp);
   Hcv = zeros(nobj, ncomp);

   for k = 1:nseg
      
      % get indices for local sets
      indc = cvind ~= k;
      indk = cvind == k;
      
      % split data to local calibration and validation sets
      Xc = Xs(indc, :);
      Xk = Xs(indk, :);

      % fit the local model
      [~, ~, Pk] = svds(Xc, ncomp);
      
      % make predictions for local validation set
      Tk = Xk * Pk;
      Uk = Tk * diag(1 ./ s);

      % compute distances
      for a = 1:ncomp

         aind  = 1:a;
         Tka = Tk(:, aind);
         Uka = Uk(:, aind);
         Pka = Pk(:, aind);

         Eka = Xk - Tka * Pka';
         Qcv(indk, a) = sum(Eka.^2, 2);
         Hcv(indk, a) = sum(Uka.^2, 2);
      end
   end
end

function [Qcv, Hcv] = pcacvlocal(X, ncomp, center, scale, cv) 
% PCACVLOCAL computes matrices with orthogonal distances (Qcv) and score 
% distances (Hcv) using cross-validation with local centering and scaling.
% The input parameters are the same as for PVAPCV function.
%

   nvar = size(X, 2);
   nobj = size(X, 1);
   
   % get vector with indices of samples in CV-segments
   cvind = crossval(cv, nobj, 1:nobj);
   
   % scale data globally
   if center, mX = mean(X); else, mX = zeros(1, nvar); end
   if scale, sX = std(X); else, sX = ones(1, nvar); end
   Xs = (X - mX) ./ sX;

   % compyte global model and get singular values
   [~, S, ~] = svds(Xs, ncomp);
   s = diag(S) / sqrt(nobj - 1);

   % compute distances
   nseg = max(cvind);
   Qcv = zeros(nobj, ncomp);
   Hcv = zeros(nobj, ncomp);

   for k = 1:nseg
      
      % get indices for local sets
      indc = cvind ~= k;
      indk = cvind == k;
      
      % split data to local calibration and validation sets
      Xc = X(indc, :);
      Xk = X(indk, :);

      % autoscale data locally
      if center, mXl = mean(Xc); else, mXl = zeros(1, nvar); end
      if scale, sXl = std(Xc); else, sXl = ones(1, nvar); end
      Xc = (Xc - mXl) ./ sXl;
      Xk = (Xk - mXl) ./ sXl;
      
      % fit the local model
      [~, ~, Pk] = svds(Xc, ncomp);
      
      % make predictions for local validation set
      Tk = Xk * Pk;
      Uk = Tk * diag(1 ./ s);

      % compute and save the distances
      for a = 1:ncomp

         aind  = 1:a;
         Tka = Tk(:, aind);
         Uka = Uk(:, aind);
         Pka = Pk(:, aind);

         Eka = Xk - Tka * Pka';
         Qcv(indk, a) = sum(Eka.^2, 2);
         Hcv(indk, a) = sum(Uka.^2, 2);
      end
   end
end