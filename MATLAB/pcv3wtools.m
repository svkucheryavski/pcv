function t = pcv3wtools()
%PCV3WTOOLS Helper functions for three-way PCV methods (internal)
%
% Returns a structure with handles to the helper functions used by PCV3WFU, PCV3WBU
% and PCV3WHSV. Not intended to be called directly.

   t = struct( ...
      'getdim', @getdim, ...
      'normalizescale', @normalizescale, ...
      'splitcv', @splitcv, ...
      'setseed', @setseed, ...
      'prep', @prep, ...
      'unprep', @unprep, ...
      'numrank', @numrank, ...
      'slab', @slab, ...
      'unfold', @unfold, ...
      'modefactors', @modefactors, ...
      'proj', @proj, ...
      'complvec', @complvec, ...
      'complementdir', @complementdir, ...
      'getxpvorth', @getxpvorth, ...
      'lsap', @lsap, ...
      'align', @align, ...
      'validsubranks', @validsubranks);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checks and settings       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = getdim(X)
% check that X is a numeric 3-D array and return its size

   % an I x J x 1 array is stored as a matrix, so both cases are reported together
   if ~isnumeric(X) || ndims(X) ~= 3
      error('pcv3w:notThreeWay', 'X must be a numeric 3-D array with at least 2 values in each mode.');
   end

   d = size(X);
   if any(d < 2)
      error('pcv3w:tooSmall', 'X must be a numeric 3-D array with at least 2 values in each mode.');
   end
end


function scale = normalizescale(scale)
% convert parameter "Scale" to one of "none", "mode2", "mode3", "both"

   if islogical(scale) && isscalar(scale)
      if scale, scale = "both"; else, scale = "none"; end
      return;
   end

   if (ischar(scale) || isstring(scale))
      scale = string(scale);
      if isscalar(scale) && any(scale == ["none", "mode2", "mode3", "both"])
         return;
      end
   end

   error('pcv3w:wrongScale', ...
      'Parameter "Scale" must be false, true or one of "none", "mode2", "mode3", "both".');
end


function cvind = splitcv(CV, nObj)
% segment number of each sample, checked to be 1, 2, ..., nSeg

   if isnumeric(CV) && ~isscalar(CV)
      cvind = double(CV(:));
   else
      cvind = crossval(CV, nObj);
      cvind = double(cvind(:));
   end
   nSeg = max(cvind);

   if numel(cvind) ~= nObj || any(~isfinite(cvind)) || any(cvind ~= round(cvind)) || ...
         nSeg < 2 || any(cvind < 1) || numel(unique(cvind)) ~= nSeg
      error('pcv3w:wrongSegments', ...
         'Segment numbers in "CV" must be 1, 2, ..., n (n >= 2), each used at least once.');
   end
end


function cleanup = setseed(seed)
% set the random seed and return an object which restores the previous state when cleared

   if isempty(seed)
      cleanup = [];
      return;
   end

   state = rng();
   rng(seed);
   cleanup = onCleanup(@() rng(state));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = prep(X, center, scale)
% center (mean slab) and scale (within modes) a 3-D array

   d = size(X);
   p = struct('M', zeros(d(2), d(3)), 'S2', ones(d(2), 1), 'S3', ones(d(3), 1), ...
      'center', center, 'scale', scale);

   if center
      p.M = reshape(mean(X, 1), d(2), d(3));
   end
   Xp = X - reshape(p.M, [1, d(2), d(3)]);

   % zero variance check: sd below the rounding noise of centering (10 * I * eps times the
   % largest absolute value of the variable before centering) counts as zero
   A = abs(X) * (10 * d(1) * eps);

   if any(scale == ["mode2", "both"])
      p.S2 = std(unfold(Xp, 2), 0, 2);
      checkstd(p.S2, max(unfold(A, 2), [], 2), 2);
      Xp = Xp ./ reshape(p.S2, [1, d(2), 1]);
      A = A ./ reshape(p.S2, [1, d(2), 1]);
   end

   if any(scale == ["mode3", "both"])
      p.S3 = std(unfold(Xp, 3), 0, 2);
      checkstd(p.S3, max(unfold(A, 3), [], 2), 3);
      Xp = Xp ./ reshape(p.S3, [1, 1, d(3)]);
   end

   p.Xp = Xp;
end


function checkstd(s, tol, mode)
% raise an error if some variables of a mode have zero variance (sd not above tol)

   bad = find(~(isfinite(s) & s > tol));
   if ~isempty(bad)
      error('pcv3w:zeroVariance', 'Cannot scale mode %d, zero variance in variable(s) %s.', ...
         mode, strjoin(string(bad(1:min(5, end))'), ', '));
   end
end


function X = unprep(Xp, p)
% undo preprocessing made by prep() (back to original units)

   d = size(Xp);
   X = Xp .* reshape(p.S2, [1, d(2), 1]) .* reshape(p.S3, [1, 1, d(3)]) + ...
      reshape(p.M, [1, d(2), d(3)]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model helpers             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = numrank(s)
% number of singular values above the relative tolerance

   if isempty(s) || s(1) <= 0
      r = 0;
      return;
   end
   r = sum(s > sqrt(eps) * s(1));
end


function S = slab(X, i)
% slab i of a 3-D array as a J x K matrix

   d = size(X);
   S = reshape(X(i, :, :), d(2), d(3));
end


function M = unfold(X, mode)
% mode 2 unfolding (J x (K*I)) or mode 3 unfolding (K x (J*I))

   d = size(X);
   if numel(d) < 3, d(3) = 1; end
   if mode == 2
      M = reshape(permute(X, [2, 3, 1]), d(2), []);
   else
      M = reshape(permute(X, [3, 2, 1]), d(3), []);
   end
end


function [U2, U3, r2, r3] = modefactors(X, R2, R3)
% mode 2 and mode 3 factors (truncated SVD of unfoldings) and numerical ranks

   [U2, S2] = svd(unfold(X, 2), 'econ');
   [U3, S3] = svd(unfold(X, 3), 'econ');
   s2 = diag(S2);
   s3 = diag(S3);
   U2 = U2(:, 1:min(R2, size(U2, 2)));
   U3 = U3(:, 1:min(R3, size(U3, 2)));
   r2 = numrank(s2);
   r3 = numrank(s3);
end


function P = proj(S, U2, U3)
% projection of a slab onto range(U2) x range(U3)

   P = U2 * (U2' * S * U3) * U3';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Orthogonal part           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = complvec(U)
% unit vector orthogonal to the columns of U

   [~, j] = min(sum(U.^2, 2));
   e = zeros(size(U, 1), 1);
   e(j) = 1;
   w = e - U * (U' * e);
   w = w / sqrt(sum(w.^2));
end


function cdir = complementdir(U2, U3)
% fixed unit direction outside the model space (empty if there is none)

   if size(U2, 2) < size(U2, 1)
      cdir = complvec(U2) * firstfactor(U3)';
   elseif size(U3, 2) < size(U3, 1)
      cdir = firstfactor(U2) * complvec(U3)';
   else
      cdir = [];
   end
end


function u = firstfactor(U)
% first factor with sign fixed, so the direction does not depend on the SVD sign convention

   [~, j] = max(abs(U(:, 1)));
   u = U(:, 1) * sign(U(j, 1));
end


function E = getxpvorth(Xk, Eres, U2, U3, cdir, orth)
% orthogonal parts of the PV-set slabs of one segment: outside the global model space,
% with the norm of the local residual, direction from a random mix of the segment's
% slabs ("random") or from the slab's own local residual ("canonical")

   nk = numel(Xk);
   J = size(U2, 1);
   K = size(U3, 1);
   tol = 100 * eps;

   if orth == "random"
      Z = 2 * rand(nk) - 1;
   end
   E = cell(nk, 1);

   for i = 1:nk
      qi = sum(Eres{i}(:).^2);
      if qi <= 0
         E{i} = zeros(J, K);
         continue;
      end

      if orth == "canonical"
         src = Eres{i};
      else
         src = zeros(J, K);
         for j = 1:nk
            src = src + Z(i, j) * Xk{j};
         end
      end

      Eo = src - proj(src, U2, U3);

      % source lies (almost) inside the model space: use the fixed direction
      if sum(Eo(:).^2) <= tol * sum(src(:).^2)
         if isempty(cdir)
            E{i} = zeros(J, K);
            continue;
         end
         Eo = cdir;
      end

      E{i} = Eo * sqrt(qi / sum(Eo(:).^2));
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matching components       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = lsap(M)
% assignment maximizing the sum of M(i, res(i)) (Hungarian method)

   n = size(M, 1);
   C = max(M(:)) - M;

   % potentials and matching use index 0 as a sentinel, stored at position 1
   u = zeros(n + 1, 1);
   v = zeros(n + 1, 1);
   p = zeros(n + 1, 1);
   way = zeros(n + 1, 1);

   for i = 1:n
      p(1) = i;
      j0 = 0;
      minv = inf(n + 1, 1);
      used = false(n + 1, 1);
      while true
         used(j0 + 1) = true;
         i0 = p(j0 + 1);
         delta = inf;
         j1 = 0;
         for j = 1:n
            if used(j + 1), continue; end
            cur = C(i0, j) - u(i0 + 1) - v(j + 1);
            if cur < minv(j + 1)
               minv(j + 1) = cur;
               way(j + 1) = j0;
            end
            if minv(j + 1) < delta
               delta = minv(j + 1);
               j1 = j;
            end
         end
         for j = 0:n
            if used(j + 1)
               u(p(j + 1) + 1) = u(p(j + 1) + 1) + delta;
               v(j + 1) = v(j + 1) - delta;
            else
               minv(j + 1) = minv(j + 1) - delta;
            end
         end
         j0 = j1;
         if p(j0 + 1) == 0, break; end
      end
      while true
         j1 = way(j0 + 1);
         p(j0 + 1) = p(j1 + 1);
         j0 = j1;
         if j0 == 0, break; end
      end
   end

   res = zeros(n, 1);
   for j = 1:n
      res(p(j + 1)) = j;
   end
end


function [Ual, assigned, cosines] = align(Ug, Ul, method)
% match local factors to global ones by permutation and sign; returns the matched local
% factors, the matched local component for each global one and their absolute cosines

   R = size(Ug, 2);
   if size(Ul, 2) ~= R
      error('pcv3w:wrongFactors', 'Factor matrices must have the same number of columns.');
   end

   M = abs(Ug' * Ul);
   if method == "optimal"
      assigned = lsap(M);
   else
      assigned = zeros(R, 1);
      used = false(R, 1);
      for g = 1:R
         cand = find(~used);
         [~, b] = max(M(g, cand));
         assigned(g) = cand(b);
         used(cand(b)) = true;
      end
   end

   Ual = Ul(:, assigned);
   sgn = sign(sum(Ug .* Ual, 1));
   sgn(sgn == 0) = 1;
   Ual = Ual .* sgn;
   cosines = M(sub2ind([R, R], (1:R)', assigned));
end


function valid = validsubranks(alignment, mode, R)
% numbers of components r for which components 1..r are matched within 1..r in every segment

   ok = true(R, 1);
   rows = alignment(alignment(:, 2) == mode, :);
   for k = unique(rows(:, 1))'
      perm = rows(rows(:, 1) == k, 4);
      ok = ok & (cummax(perm) == (1:R)');
   end
   valid = find(ok)';
end
