function [Xpv, info] = pcv3whsv(X, nComp, Center, Scale, CV, Orth, Align, Seed)
%PCV3WHSV Procrustes cross-validation for three-way data based on HO-SVD (HSV)
%
% Arguments:
% ----------
% X         3-D array with calibration set, size I x J x K (samples x mode 2 x mode 3),
%           all values must be finite
% nComp     vector with two positive integers, [R2, R3], number of components in
%           modes 2 and 3
% Center    logical, subtract the mean slab or not (default: true)
% Scale     scaling within modes: "none" (or false), "mode2", "mode3" or "both"
%           (or true); a variable with zero variance (also when it is constant over
%           the samples) cannot be scaled and gives an error (default: "none")
% CV        defines cross-validation splits, see PCVPCA (default: {"ven", 4})
% Orth      how to construct the orthogonal part: "random" (random mix of the
%           segment's slabs, similar to the two-way method) or "canonical" (each
%           slab's own local residual, no randomness) (default: "random")
% Align     how to match local components to the global ones: "optimal" (maximal
%           sum of absolute cosines) or "greedy" (global components in order of
%           variance choose first) (default: "optimal")
% Seed      optional number, makes random splits and random orthogonal parts
%           reproducible; the state of the random generator is restored afterwards
%
% Returns:
% --------
% Xpv       3-D array with PV-set (same size as X)
% info      structure with fields: ncomp (numbers of components actually used),
%           cvind (segment number of each sample), orth, align, scale, prep (global
%           preprocessing parameters), alignment and validSubranks (see below)
%
% The method models X with truncated HO-SVD in modes 2 and 3 (a Tucker2 model) and
% generates PV-set, which, when projected on the global model, gives the same scores
% and orthogonal distances as cross-validation of the model, sample by sample.
%
% Scaling divides mode 2 and/or mode 3 variables by their standard deviation after
% centering. With "both", mode 2 is scaled first and mode 3 is then scaled using the
% mode 2 scaled data. The PV-set is returned in the original units.
%
% info.alignment is a matrix with columns: segment, mode, global component, matched
% local component and absolute cosine between the two.
%
% A PV-set generated with [R2, R3] components also reproduces cross-validation of the
% models with fewer components, r2 <= R2 and r3 <= R3, unless the local components
% are matched to the global ones across the boundary of the first r2 (or r3)
% components. The numbers of components for which this holds are given in
% info.validSubranks, and the method gives a warning if some of them are missing.
%
% Use [] for an argument to keep its default value.
%
% See also: PCV3WFU, PCV3WBU, PCVPCA

   t = pcv3wtools();
   d = t.getdim(X);

   if nargin < 3 || isempty(Center), Center = true; end
   if nargin < 4 || isempty(Scale), Scale = "none"; end
   if nargin < 5 || isempty(CV), CV = {"ven", 4}; end
   if nargin < 6 || isempty(Orth), Orth = "random"; end
   if nargin < 7 || isempty(Align), Align = "optimal"; end
   if nargin < 8, Seed = []; end

   Orth = string(Orth);
   if ~isscalar(Orth) || ~any(Orth == ["random", "canonical"])
      error('pcv3w:wrongOrth', 'pcv3whsv: parameter "Orth" must be "random" or "canonical".');
   end

   Align = string(Align);
   if ~isscalar(Align) || ~any(Align == ["optimal", "greedy"])
      error('pcv3w:wrongAlign', 'pcv3whsv: parameter "Align" must be "optimal" or "greedy".');
   end

   Scale = t.normalizescale(Scale);

   if ~all(isfinite(X(:)))
      error('pcv3w:notFinite', ...
         'pcv3whsv: X contains NaN or Inf values, interpolate or impute them first.');
   end

   if nargin < 2 || ~isnumeric(nComp) || numel(nComp) ~= 2 || ~all(isfinite(nComp)) || ...
         any(nComp < 1) || any(nComp ~= round(nComp))
      error('pcv3w:wrongNComp', 'pcv3whsv: parameter "nComp" must be two positive integers, [R2, R3].');
   end

   cleanup = t.setseed(Seed); %#ok<NASGU>
   cvind = t.splitcv(CV, d(1));
   [Xpv, info] = run(t, double(X), double(nComp(:)'), cvind, logical(Center), Scale, Orth, Align);
end


function [Xpv, info] = run(t, X, nComp, cvind, Center, Scale, Orth, Align)
% main loop of the method (arguments already checked)

   d = size(X);
   nSeg = max(cvind);

   % global preprocessing and global factors
   prep = t.prep(X, Center, Scale);
   [U2g, U3g, r2min, r3min] = t.modefactors(prep.Xp, nComp(1), nComp(2));

   % local factors for every segment
   U2l = cell(nSeg, 1);
   U3l = cell(nSeg, 1);
   for k = 1:nSeg
      [U2l{k}, U3l{k}, r2, r3] = t.modefactors(prep.Xp(cvind ~= k, :, :), nComp(1), nComp(2));
      r2min = min(r2min, r2);
      r3min = min(r3min, r3);
   end

   % ranks attainable by the global and all local models
   R2 = min(nComp(1), r2min);
   R3 = min(nComp(2), r3min);

   if R2 < 1 || R3 < 1
      error('pcv3w:noComponents', ...
         'pcv3whsv: the data cannot supply at least one component in every mode for all segments.');
   end

   if R2 < nComp(1) || R3 < nComp(2)
      warning('pcv3w:ncompReduced', ...
         'pcv3whsv: "nComp" reduced to [%d, %d] (maximum for all segments).', R2, R3);
   end

   U2 = U2g(:, 1:R2);
   U3 = U3g(:, 1:R3);
   cdir = t.complementdir(U2, U3);

   Xpv = zeros(d);
   alignment = zeros(0, 5);

   for k = 1:nSeg
      indk = find(cvind == k);
      nk = numel(indk);
      Xk = cell(nk, 1);
      for m = 1:nk
         Xk{m} = t.slab(prep.Xp, indk(m));
      end

      % local factors matched to the global ones
      [U2k, perm2, cos2] = t.align(U2, U2l{k}(:, 1:R2), Align);
      [U3k, perm3, cos3] = t.align(U3, U3l{k}(:, 1:R3), Align);
      alignment = [alignment; ...
         repmat([k, 2], R2, 1), (1:R2)', perm2, cos2; ...
         repmat([k, 3], R3, 1), (1:R3)', perm3, cos3]; %#ok<AGROW>

      % explained part: local scores with global factors; residuals of the local model
      Xhat = cell(nk, 1);
      Eres = cell(nk, 1);
      for m = 1:nk
         Ti = U2k' * Xk{m} * U3k;
         Eres{m} = Xk{m} - U2k * Ti * U3k';
         Xhat{m} = U2 * Ti * U3';
      end

      Eo = t.getxpvorth(Xk, Eres, U2, U3, cdir, Orth);
      for m = 1:nk
         Xpv(indk(m), :, :) = reshape(Xhat{m} + Eo{m}, [1, d(2), d(3)]);
      end
   end

   Xpv = t.unprep(Xpv, prep);

   % sub-rank r is valid if in every segment components 1..r are matched within 1..r
   valid2 = t.validsubranks(alignment, 2, R2);
   valid3 = t.validsubranks(alignment, 3, R3);

   if numel(valid2) < R2 || numel(valid3) < R3
      warning('pcv3w:invalidSubranks', ['pcv3whsv: PV-set does not reproduce cross-validation ' ...
         'for some smaller numbers of components, see info.validSubranks.']);
   end

   info = struct();
   info.ncomp = [R2, R3];
   info.cvind = cvind;
   info.orth = Orth;
   info.align = Align;
   info.scale = Scale;
   info.prep = rmfield(prep, 'Xp');
   info.alignment = alignment;
   info.validSubranks = struct('r2', valid2, 'r3', valid3);
end
