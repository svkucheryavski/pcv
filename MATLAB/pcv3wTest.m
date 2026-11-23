%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests for three-way PCV    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tests = pcv3wTest
   tests = functiontests(localfunctions);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSV method                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testPreprocessing(testCase)
% HSV PV-set reproduces cross-validation for all preprocessing settings

   X = simdata3w(24, 10, 8, 1);
   q = quietwarnings(); %#ok<NASGU>
   for center = [true, false]
      for scale = ["none", "mode2", "mode3", "both"]
         [Xpv, info] = pcv3whsv(X, [3, 3], center, scale, [], [], [], 1);
         verifySize(testCase, Xpv, size(X));
         verifyEqual(testCase, info.ncomp, [3, 3]);
         checkhsvcv(testCase, X, Xpv, info, center, scale, 1e-8);
      end
   end
end


function testSplitsAndOptions(testCase)
% HSV PV-set reproduces cross-validation for different splits and settings

   X = simdata3w(24, 10, 8, 2);
   cvs = {{"ven", 3}, {"rand", 5}, {"loo"}, repmat(1:4, 1, 6)};
   q = quietwarnings(); %#ok<NASGU>
   for c = 1:numel(cvs)
      for orth = ["random", "canonical"]
         for align = ["optimal", "greedy"]
            [Xpv, info] = pcv3whsv(X, [4, 3], true, "none", cvs{c}, orth, align, 3);
            checkhsvcv(testCase, X, Xpv, info, true, "none", 1e-8);
         end
      end
   end
end


function testSubranks(testCase)
% some smaller numbers of components are not valid: warning, and they really fail

   % find data where some smaller numbers of components are not valid
   found = false;
   for seed = 1:200
      X = simdata3w(24, 10, 8, seed);
      q = quietwarnings(); %#ok<NASGU>
      [~, info] = pcv3whsv(X, [5, 5], [], [], [], [], [], 1);
      clear q;
      if numel(info.validSubranks.r2) < 5 && numel(info.validSubranks.r3) < 5
         found = true;
         break;
      end
   end
   assertTrue(testCase, found, 'no dataset with invalid smaller numbers of components found');

   [Xpv, info] = verifyWarning(testCase, @() pcv3whsv(X, [5, 5], [], [], [], [], [], 1), ...
      'pcv3w:invalidSubranks');
   vs = info.validSubranks;
   verifyTrue(testCase, any(vs.r2 == 5) && any(vs.r3 == 5));
   checkhsvcv(testCase, X, Xpv, info, true, "none", 1e-8);

   % omitted numbers of components really do not reproduce cross-validation
   cvind = info.cvind;
   M = mean(X, 1);
   Xc = X - M;
   Pc = Xpv - M;
   U2 = reffactors(Xc, 2, 5);
   U3 = reffactors(Xc, 3, 5);
   for mode = [2, 3]
      if mode == 2, valid = vs.r2; else, valid = vs.r3; end
      for r = setdiff(1:5, valid)
         if mode == 2, r2 = r; r3 = 5; else, r2 = 5; r3 = r; end
         maxdiff = 0;
         for k = 1:max(cvind)
            l2 = reffactors(Xc(cvind ~= k, :, :), 2, r2);
            l3 = reffactors(Xc(cvind ~= k, :, :), 3, r3);
            for i = find(cvind == k)'
               Xi = slab(Xc, i);
               Pi = slab(Pc, i);
               qcv = sum(sum((Xi - l2 * (l2' * Xi * l3) * l3').^2));
               qpv = sum(Pi(:).^2) - sum(sum((U2(:, 1:r2)' * Pi * U3(:, 1:r3)).^2));
               maxdiff = max(maxdiff, abs(qpv - qcv) / qcv);
            end
         end
         verifyGreaterThan(testCase, maxdiff, 1e-6, sprintf('mode %d, %d components', mode, r));
      end
   end
end


function testInfo(testCase)
% information returned together with the PV-set is correct

   X = simdata3w(24, 10, 8, 1);
   [~, info] = pcv3whsv(X, [3, 2], [], true, [], "canonical", [], 1);
   verifyEqual(testCase, info.cvind, mod((0:23)', 4) + 1);
   verifyEqual(testCase, info.orth, "canonical");
   verifyEqual(testCase, info.align, "optimal");
   verifyEqual(testCase, info.scale, "both");
   verifySize(testCase, info.prep.M, [10, 8]);
   verifySize(testCase, info.alignment, [4 * (3 + 2), 5]);
   verifyTrue(testCase, all(info.alignment(:, 5) <= 1 + 1e-12));
end


function testCanonicalIsDeterministic(testCase)
% canonical orthogonal part does not depend on random numbers

   X = simdata3w(24, 10, 8, 1);
   Xpv1 = pcv3whsv(X, [3, 3], [], [], [], "canonical", [], 1);
   Xpv2 = pcv3whsv(X, [3, 3], [], [], [], "canonical", [], 2);
   verifyEqual(testCase, Xpv1, Xpv2);
end


function testNCompReduced(testCase)
% too large number of components is reduced with a warning

   X = simdata3w(12, 4, 5, 1);
   [Xpv, info] = verifyWarning(testCase, @() pcv3whsv(X, [6, 6], [], [], {"ven", 3}, [], [], 1), ...
      'pcv3w:ncompReduced');
   verifyEqual(testCase, info.ncomp, [4, 5]);
   checkhsvcv(testCase, X, Xpv, info, true, "none", 1e-8);

   % noiseless data of rank one: numerical rank decides
   rng(1);
   b = rand(1, 6);
   c = rand(1, 1, 5);
   X = randn(12, 1) .* b .* c + 5;
   [Xpv, info] = verifyWarning(testCase, @() pcv3whsv(X, [3, 3], [], [], [], [], [], 1), ...
      'pcv3w:ncompReduced');
   verifyEqual(testCase, info.ncomp, [1, 1]);
   checkhsvcv(testCase, X, Xpv, info, true, "none", 1e-8);

   % second component only in sample 1: local model without it has rank one
   X = randn(12, 1) .* b .* c;
   X(1, :, :) = X(1, :, :) + reshape(rand(6, 1) * rand(1, 5), 1, 6, 5);
   [Xpv, info] = verifyWarning(testCase, @() pcv3whsv(X, [2, 2], false, [], [], [], [], 1), ...
      'pcv3w:ncompReduced');
   verifyEqual(testCase, info.ncomp, [1, 1]);
   checkhsvcv(testCase, X, Xpv, info, false, "none", 1e-8);
end


function testAlignment(testCase)
% optimal and greedy matching give expected assignments

   t = pcv3wtools();
   M = [0.90, 0.80, 0.10; 0.85, 0.10, 0.30; 0.20, 0.75, 0.70];
   Ug = eye(3);
   Ul = M .* [1, -1, 1];

   [Ua, perm, cosines] = t.align(Ug, Ul, "optimal");
   verifyEqual(testCase, perm, [2; 1; 3]);
   verifyEqual(testCase, cosines, [0.80; 0.85; 0.70], 'AbsTol', 1e-12);
   verifyTrue(testCase, all(sum(Ug .* Ua, 1) > 0));

   [~, perm] = t.align(Ug, Ul, "greedy");
   verifyEqual(testCase, perm, [1; 3; 2]);

   % optimal matching reaches the maximal sum for random matrices
   rng(7);
   for n = 2:6
      for rep = 1:5
         A = rand(n);
         p = t.lsap(A);
         verifyEqual(testCase, sort(p), (1:n)');
         verifyLessThan(testCase, abs(sum(A(sub2ind([n, n], (1:n)', p))) - bestsum(A)), 1e-12);
      end
   end

   % greedy matching agrees with the reference routine for random matrices
   for n = 2:8
      for rep = 1:5
         [Ul, ~] = qr(randn(n));
         [~, perm] = t.align(eye(n), Ul, "greedy");
         verifyEqual(testCase, perm, refgreedy(abs(Ul)));
      end
   end
end


function testFixedDirection(testCase)
% orthogonal part uses the fixed direction when needed

   t = pcv3wtools();
   rng(5);
   [U2, ~] = qr(randn(6, 2), 0);
   [U3, ~] = qr(randn(5, 2), 0);

   % direction does not depend on signs of the factors
   cdir = t.complementdir(U2, U3);
   verifyEqual(testCase, t.complementdir(-U2, -U3), cdir, 'AbsTol', 1e-12);
   verifyEqual(testCase, sum(cdir(:).^2), 1, 'AbsTol', 1e-12);
   verifyEqual(testCase, sum(sum(t.proj(cdir, U2, U3).^2)), 0, 'AbsTol', 1e-12);

   % residual lying inside the model space: the fixed direction is used with the right norm
   E = U2 * [1, 3; 2, 4] * U3';
   Eo = t.getxpvorth({E}, {E}, U2, U3, cdir, "canonical");
   verifyEqual(testCase, Eo{1}, cdir * norm(E, 'fro'), 'AbsTol', 1e-12);
   Eo = t.getxpvorth({E; 2 * E}, {E; E}, U2, U3, cdir, "random");
   verifyEqual(testCase, Eo{1}, cdir * norm(E, 'fro'), 'AbsTol', 1e-12);

   % full model space in mode 2 only: direction comes from mode 3
   cdir = t.complementdir(eye(3), U3);
   verifySize(testCase, cdir, [3, 5]);
   verifyEqual(testCase, sum(cdir(:).^2), 1, 'AbsTol', 1e-12);
   verifyEqual(testCase, sum(sum(t.proj(cdir, eye(3), U3).^2)), 0, 'AbsTol', 1e-12);

   % full model space in both modes: no direction
   verifyEmpty(testCase, t.complementdir(eye(3), eye(2)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FU and BU methods          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testFU(testCase)
% FU PV-set is the same as PCVPCA on unfolded data

   X = simdata3w(24, 10, 8, 1);
   [I, J, K] = size(X);
   for scale = [false, true]
      Xpv = pcv3wfu(X, 6, true, scale, [], 5);
      rng(5);
      Xpv2 = pcvpca(reshape(X, I, J * K), 6, true, scale, {"ven", 4});
      verifyEqual(testCase, Xpv, reshape(Xpv2, I, J, K));

      % random split and no centering
      Xpv3 = pcv3wfu(X, 6, false, scale, {"rand", 3}, 5);
      rng(5);
      cvind = crossval({"rand", 3}, I);
      verifyEqual(testCase, Xpv3, reshape(pcvpca(reshape(X, I, J * K), 6, false, scale, cvind), I, J, K));

      % column (j, k) of the unfolded data is column j + J * (k - 1)
      for j = [1, J]
         for k = [1, 2, K]
            verifyEqual(testCase, Xpv(:, j, k), Xpv2(:, j + J * (k - 1)));
         end
      end
   end
end


function testBU(testCase)
% BU PV-set is the same as PCVPCA on each block with one split

   X = simdata3w(24, 10, 8, 1);
   [I, J, K] = size(X);
   cv = {"rand", 4};
   for center = [true, false]
      for mode = [2, 3]
         Xpv = pcv3wbu(X, mode, 4, center, false, cv, 5);
         rng(5);
         cvind = crossval(cv, I);
         if mode == 2, nb = J; else, nb = K; end
         for b = 1:nb
            if mode == 2
               Xb = reshape(X(:, b, :), I, K);
               Pb = reshape(Xpv(:, b, :), I, K);
            else
               Xb = X(:, :, b);
               Pb = Xpv(:, :, b);
            end
            verifyEqual(testCase, Pb, pcvpca(Xb, 4, center, false, cvind));
         end
      end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seed and wrong input       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testSeed(testCase)
% seed makes results reproducible and keeps the caller's random state

   X = simdata3w(24, 10, 8, 1);
   funs = {
      @(s) pcv3whsv(X, [3, 3], [], [], [], [], [], s)
      @(s) pcv3wfu(X, 5, [], [], [], s)
      @(s) pcv3wbu(X, 2, 4, [], [], [], s)
   };
   for f = 1:numel(funs)
      rng(99);
      r1 = rand();
      rng(99);
      P1 = funs{f}(10);
      r2 = rand();
      verifyEqual(testCase, r1, r2);
      verifyEqual(testCase, funs{f}(10), P1);
      d = funs{f}(11) - P1;
      verifyGreaterThan(testCase, max(abs(d(:))), 1e-8);
   end
end


function testWrongInput(testCase)
% wrong input gives an error

   X = simdata3w(24, 10, 8, 1);
   verifyError(testCase, @() pcv3whsv(ones(4, 5), [2, 2]), 'pcv3w:notThreeWay');
   verifyError(testCase, @() pcv3whsv(X(:, 1, :), [2, 2]), 'pcv3w:tooSmall');
   verifyError(testCase, @() pcv3whsv(X), 'pcv3w:wrongNComp');
   ncomps = {2, [0, 2], [2.5, 2], [2, 2, 2], 'ab', [], [NaN, 2]};
   for n = 1:numel(ncomps)
      verifyError(testCase, @() pcv3whsv(X, ncomps{n}), 'pcv3w:wrongNComp');
   end

   Xna = X;
   Xna(1, 1, 1) = NaN;
   verifyError(testCase, @() pcv3whsv(Xna, [2, 2]), 'pcv3w:notFinite');

   Xc = X;
   Xc(:, 2, :) = 1;
   verifyError(testCase, @() pcv3whsv(Xc, [2, 2], true, "mode2"), 'pcv3w:zeroVariance');

   % constant over the samples only: centered values are rounding noise, still an error
   Xc = X;
   Xc(:, 2, :) = repmat(reshape(0.1 * (1:8), 1, 1, 8), 24, 1, 1);
   Xc(:, :, 3) = repmat(0.1 * (1:10), 24, 1);
   verifyErrorText(testCase, @() pcv3whsv(Xc, [2, 2], true, "mode2"), 'pcv3w:zeroVariance', ...
      "mode 2, zero variance in variable(s) 2.");
   verifyErrorText(testCase, @() pcv3whsv(Xc, [2, 2], true, "mode3"), 'pcv3w:zeroVariance', ...
      "mode 3, zero variance in variable(s) 3.");
   verifyErrorText(testCase, @() pcv3whsv(Xc, [2, 2], true, "both"), 'pcv3w:zeroVariance', ...
      "mode 2, zero variance");

   % mode 3 stage of "both" on its own
   Xc = X;
   Xc(:, :, 3) = repmat(0.1 * (1:10), 24, 1);
   verifyErrorText(testCase, @() pcv3whsv(Xc, [2, 2], true, "both"), 'pcv3w:zeroVariance', ...
      "mode 3, zero variance in variable(s) 3.");

   % large values with small but real variation are accepted
   Xl = 1e12 + repmat((-1).^(1:1000)', 1, 2, 2);
   verifyWarningFree(testCase, @() pcv3whsv(Xl, [1, 1], true, "both", [], [], [], 1));

   % without centering a variable with all values equal is still an error
   Xc = X;
   Xc(:, 2, :) = 0.7;
   verifyErrorText(testCase, @() pcv3whsv(Xc, [2, 2], false, "mode2"), 'pcv3w:zeroVariance', ...
      "mode 2, zero variance in variable(s) 2.");

   verifyError(testCase, @() pcv3whsv(X, [2, 2], true, "wrong"), 'pcv3w:wrongScale');
   verifyError(testCase, @() pcv3whsv(X, [2, 2], true, 1), 'pcv3w:wrongScale');
   verifyError(testCase, @() pcv3whsv(X, [2, 2], [], [], [], "wrong"), 'pcv3w:wrongOrth');
   verifyError(testCase, @() pcv3whsv(X, [2, 2], [], [], [], [], "wrong"), 'pcv3w:wrongAlign');
   verifyError(testCase, @() pcv3wbu(X, 1), 'pcv3w:wrongMode');
   verifyError(testCase, @() pcv3wfu(ones(4, 5)), 'pcv3w:notThreeWay');
end


function testWrongSegments(testCase)
% wrong segment numbers in CV give an error

   X = simdata3w(24, 10, 8, 1);
   cvs = {repmat([0, 1, 2], 1, 8), repmat([1, 3], 1, 12), ones(1, 24), ...
      repmat([1.5, 2], 1, 12), [1, 2, 1, 2]};
   for c = 1:numel(cvs)
      verifyError(testCase, @() pcv3whsv(X, [2, 2], [], [], cvs{c}), 'pcv3w:wrongSegments');
      verifyError(testCase, @() pcv3wbu(X, 2, 3, [], [], cvs{c}), 'pcv3w:wrongSegments');
      verifyError(testCase, @() pcv3wfu(X, 3, [], [], cvs{c}), 'pcv3w:wrongSegments');
   end

   % vector given as a row or as a column works the same way
   Xpv1 = pcv3wfu(X, 3, [], [], repmat([1, 2], 1, 12), 1);
   Xpv2 = pcv3wfu(X, 3, [], [], repmat([1; 2], 12, 1), 1);
   verifyEqual(testCase, Xpv1, Xpv2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference values from R    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testReferences(testCase)
% HSV results match the reference values computed in R

   root = string(fileparts(fileparts(mfilename('fullpath'))));
   dataFile = fullfile(root, ".tests", "data", "sim3w-20x8x6.csv");
   caseDir = fullfile(root, ".tests", "pcv3whsv") + filesep;
   assertTrue(testCase, isfile(dataFile) && isfolder(caseDir), 'reference files not found');

   % dataset is stored unfolded, I x (J*K) with mode 2 running fastest
   X = reshape(csvread(dataFile), 20, 8, 6);

   cases = {};
   for center = [true, false]
      for scale = ["none", "mode2", "mode3", "both"]
         cases{end + 1} = {[3, 3], center, scale, "optimal", {"ven", 4}, "ven4"}; %#ok<AGROW>
      end
   end
   for align = ["optimal", "greedy"]
      cases{end + 1} = {[4, 2], true, "none", align, {"ven", 4}, "ven4"}; %#ok<AGROW>
      cases{end + 1} = {[4, 2], true, "none", align, {"loo"}, "loo"}; %#ok<AGROW>
   end

   bString = ["FALSE", "TRUE"];
   q = quietwarnings(); %#ok<NASGU>
   for c = 1:numel(cases)
      [ncomp, center, scale, align, cv, cvText] = cases{c}{:};
      sfx = sprintf("-%dx%d-%s-%s-%s-%s.csv", ncomp(1), ncomp(2), bString(center + 1), ...
         scale, align, cvText);

      [Xpv, info] = pcv3whsv(X, ncomp, center, scale, cv, "canonical", align);
      [Xpvr, infor] = pcv3whsv(X, ncomp, center, scale, cv, "random", align, 1);

      verifyEqual(testCase, reshape(Xpv, 20, 48), csvread(caseDir + "Xpv" + sfx), ...
         'RelTol', 1e-6, 'AbsTol', 1e-8, "Xpv" + sfx);
      verifyEqual(testCase, info.alignment, csvread(caseDir + "align" + sfx), ...
         'RelTol', 1e-6, 'AbsTol', 1e-8, "align" + sfx);

      % random orthogonal part: q and h are still determined (equal to cross-validation)
      [qv, hv] = refqh(X, Xpvr, infor.ncomp, center, scale);
      verifyEqual(testCase, qv, csvread(caseDir + "Q" + sfx), 'RelTol', 1e-6, "Q" + sfx);
      verifyEqual(testCase, hv, csvread(caseDir + "H" + sfx), 'RelTol', 1e-6, "H" + sfx);

      % here greedy and optimal matching differ in some segments
      checkhsvcv(testCase, X, Xpv, info, center, scale, 1e-8);
   end
end


function testJuices(testCase)
% HSV, FU and BU work on Juices data (private, set PCV_JUICES to the path of juices.csv)

   path = getenv('PCV_JUICES');
   assumeTrue(testCase, ~isempty(path) && isfile(path), 'Juices data not available (PCV_JUICES not set)');

   % first column has sample names, second class names, the rest are values
   D = readmatrix(path, 'NumHeaderLines', 1);
   D = D(:, 3:end);
   verifyEqual(testCase, size(D, 2), 61 * 81);
   verifyTrue(testCase, all(isfinite(D(:))));
   X = reshape(D, size(D, 1), 61, 81);

   cv = {"ven", 3};
   [Xpv, info] = pcv3whsv(X, [10, 10], [], [], cv, [], [], 1);
   verifyEqual(testCase, info.ncomp, [10, 10]);
   checkhsvcv(testCase, X, Xpv, info, true, "none", 1e-6);

   Xfu = pcv3wfu(X, 20, [], [], cv, 1);
   verifySize(testCase, Xfu, size(X));
   verifyTrue(testCase, all(isfinite(Xfu(:))));

   for mode = [2, 3]
      Xbu = pcv3wbu(X, mode, 20, [], [], cv, 1);
      verifySize(testCase, Xbu, size(X));
      verifyTrue(testCase, all(isfinite(Xbu(:))));
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper functions           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = simdata3w(I, J, K, seed)
% simulated three-way data: three trilinear components plus noise

   s = rng();
   rng(seed);
   A = randn(I, 3);
   B = rand(J, 3);
   C = rand(K, 3);
   X = zeros(I, J, K);
   for f = 1:3
      X = X + A(:, f) .* reshape(B(:, f), 1, J) .* reshape(C(:, f), 1, 1, K);
   end
   X = X + 0.05 * randn(I, J, K) + 5;
   rng(s);
end


function verifyErrorText(testCase, f, id, text)
% check that f gives an error with the given identifier and text in the message

   verifyError(testCase, f, id);
   try
      f();
   catch e
      verifyTrue(testCase, contains(e.message, text), ...
         sprintf('message "%s" does not contain "%s"', e.message, text));
   end
end


function c = quietwarnings()
% switch off the warnings of the method until the returned object is cleared

   s = [warning('off', 'pcv3w:invalidSubranks'), warning('off', 'pcv3w:ncompReduced')];
   c = onCleanup(@() warning(s));
end


function S = slab(X, i)
% slab i as a J x K matrix

   S = reshape(X(i, :, :), size(X, 2), size(X, 3));
end


function [Xp, M, S2, S3] = refprep(X, center, scale)
% reference preprocessing, written independently of the package code

   [~, J, K] = size(X);
   M = zeros(1, J, K);
   if center, M = mean(X, 1); end
   Xp = X - M;
   S2 = ones(1, J);
   if any(scale == ["mode2", "both"])
      for j = 1:J, v = Xp(:, j, :); S2(j) = std(v(:)); end
   end
   Xp = Xp ./ reshape(S2, 1, J);
   S3 = ones(1, K);
   if any(scale == ["mode3", "both"])
      for k = 1:K, v = Xp(:, :, k); S3(k) = std(v(:)); end
   end
   Xp = Xp ./ reshape(S3, 1, 1, K);
end


function Y = refapply(Y, M, S2, S3)
% apply the reference preprocessing to other data

   Y = (Y - M) ./ reshape(S2, 1, []) ./ reshape(S3, 1, 1, []);
end


function U = reffactors(X, mode, R)
% left singular vectors of mode 2 or mode 3 unfolding

   A = [];
   for i = 1:size(X, 1)
      if mode == 2, A = [A, slab(X, i)]; else, A = [A, slab(X, i)']; end %#ok<AGROW>
   end
   [U, ~, ~] = svd(A, 'econ');
   U = U(:, 1:R);
end


function s = bestsum(A)
% maximal sum of A(i, p(i)) over all permutations p

   n = size(A, 1);
   P = perms(1:n);
   s = max(sum(A(sub2ind([n, n], repmat(1:n, size(P, 1), 1), P)), 2));
end


function checkhsvcv(testCase, X, Xpv, info, center, scale, tol)
% check that projections of HSV PV-set on the global model equal cross-validation results

   R2 = info.ncomp(1);
   R3 = info.ncomp(2);
   cvind = info.cvind;
   al = info.alignment;

   [Xp, M, S2, S3] = refprep(X, center, scale);
   Pp = refapply(Xpv, M, S2, S3);
   U2 = reffactors(Xp, 2, R2);
   U3 = reffactors(Xp, 3, R3);

   for k = 1:max(cvind)
      L2 = reffactors(Xp(cvind ~= k, :, :), 2, R2);
      L3 = reffactors(Xp(cvind ~= k, :, :), 3, R3);

      % reported matching must reach the maximal sum of absolute cosines (checked up to 6 x 6)
      p2 = al(al(:, 1) == k & al(:, 2) == 2, 4);
      p3 = al(al(:, 1) == k & al(:, 2) == 3, 4);
      if info.align == "optimal"
         verifyTrue(testCase, isoptimal(abs(U2' * L2), p2));
         verifyTrue(testCase, isoptimal(abs(U3' * L3), p3));
      else
         verifyEqual(testCase, p2, refgreedy(abs(U2' * L2)));
         verifyEqual(testCase, p3, refgreedy(abs(U3' * L3)));
      end

      % local factors in the order reported by the method, signs matched to the global ones
      U2k = L2(:, p2);
      U3k = L3(:, p3);
      U2k = U2k .* sign(sum(U2 .* U2k, 1));
      U3k = U3k .* sign(sum(U3 .* U3k, 1));

      for i = find(cvind == k)'
         Xi = slab(Xp, i);
         Pi = slab(Pp, i);

         % all components: scores and orthogonal distance
         Tcv = U2k' * Xi * U3k;
         Tpv = U2' * Pi * U3;
         verifyEqual(testCase, Tpv, Tcv, 'AbsTol', tol * max(abs(Tcv(:))));
         qcv = sum(sum((Xi - U2k * Tcv * U3k').^2));
         qpv = sum(Pi(:).^2) - sum(Tpv(:).^2);
         verifyEqual(testCase, qpv, qcv, 'AbsTol', tol * sum(Xi(:).^2));

         % valid smaller numbers of components: compare with local models fitted directly
         for r2 = info.validSubranks.r2
            for r3 = info.validSubranks.r3
               l2 = L2(:, 1:r2);
               l3 = L3(:, 1:r3);
               qcv = sum(sum((Xi - l2 * (l2' * Xi * l3) * l3').^2));
               qpv = sum(Pi(:).^2) - sum(sum(Tpv(1:r2, 1:r3).^2));
               verifyEqual(testCase, qpv, qcv, 'AbsTol', tol * sum(Xi(:).^2));
            end
         end
      end
   end
end


function p = refgreedy(M)
% greedy matching: global components in order choose the best of the unused local ones

   n = size(M, 1);
   p = zeros(n, 1);
   for g = 1:n
      M(g, p(1:g - 1)) = -Inf;
      [~, p(g)] = max(M(g, :));
   end
end


function ok = isoptimal(M, p)
% true if assignment p gives the maximal sum of M(i, p(i)) (true without checking above 6 x 6)

   n = size(M, 1);
   if n > 6
      ok = true;
      return;
   end
   ok = abs(sum(M(sub2ind([n, n], (1:n)', p))) - bestsum(M)) < 1e-12;
end


function [qv, hv] = refqh(X, Xpv, R, center, scale)
% q and h of a PV-set on the global model (h uses the pseudo-inverse of the score scatter)

   [Xp, M, S2, S3] = refprep(X, center, scale);
   U2 = reffactors(Xp, 2, R(1));
   U3 = reffactors(Xp, 3, R(2));

   % scores are vectorized column by column, as in R
   I = size(X, 1);
   T0 = zeros(I, R(1) * R(2));
   T = zeros(I, R(1) * R(2));
   Pp = refapply(Xpv, M, S2, S3);
   for i = 1:I
      T0(i, :) = reshape(U2' * slab(Xp, i) * U3, 1, []);
      T(i, :) = reshape(U2' * slab(Pp, i) * U3, 1, []);
   end

   [~, S, V] = svd(T0, 'econ');
   s = diag(S);
   use = s > sqrt(eps) * s(1);
   qv = sum(reshape(Pp, I, []).^2, 2) - sum(T.^2, 2);
   hv = sum(((T * V(:, use)) ./ s(use)').^2, 2);
end
