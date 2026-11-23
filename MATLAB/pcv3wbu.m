function Xpv = pcv3wbu(X, Mode, nComp, Center, Scale, CV, Seed)
%PCV3WBU Procrustes cross-validation for three-way data based on block unfolding (BU)
%
% Arguments:
% ----------
% X         3-D array with calibration set, size I x J x K (samples x mode 2 x mode 3)
% Mode      which mode to split into blocks, 2 or 3 (default: 2)
% nComp     number of PCA components for each block (default: min(I - 1, size of
%           block, 30))
% Center    logical, center each block or not (default: true)
% Scale     logical, scale each block or not (default: false)
% CV        defines cross-validation splits, see PCVPCA (default: {"ven", 4})
% Seed      optional number, makes random splits and random orthogonal parts
%           reproducible; the state of the random generator is restored afterwards
%
% Returns:
% --------
% Xpv       3-D array with PV-set (same size as X)
%
% The method splits X into J matrices of size I x K (if Mode = 2) or K matrices of
% size I x J (if Mode = 3). Then it generates PV-set for each matrix using PCVPCA with
% global CV scope and combines the PV-sets into a 3-D array. All blocks use the same
% split of the samples into segments.
%
% Use [] for an argument to keep its default value.
%
% See also: PCV3WFU, PCV3WHSV, PCVPCA

   t = pcv3wtools();
   d = t.getdim(X);

   if nargin < 2 || isempty(Mode), Mode = 2; end
   if ~isnumeric(Mode) || ~isscalar(Mode) || ~any(Mode == [2, 3])
      error('pcv3w:wrongMode', 'pcv3wbu: parameter "Mode" must be 2 or 3.');
   end

   if Mode == 2, blockSize = d(3); else, blockSize = d(2); end
   if nargin < 3 || isempty(nComp), nComp = min([d(1) - 1, blockSize, 30]); end
   if nargin < 4 || isempty(Center), Center = true; end
   if nargin < 5 || isempty(Scale), Scale = false; end
   if nargin < 6 || isempty(CV), CV = {"ven", 4}; end
   if nargin < 7, Seed = []; end

   cleanup = t.setseed(Seed); %#ok<NASGU>

   X = double(X);
   cvind = t.splitcv(CV, d(1));
   Xpv = zeros(d);
   for b = 1:d(Mode)
      if Mode == 2
         Xpv(:, b, :) = reshape(pcvpca(reshape(X(:, b, :), d(1), d(3)), nComp, Center, Scale, ...
            cvind, "global"), d(1), 1, d(3));
      else
         Xpv(:, :, b) = pcvpca(X(:, :, b), nComp, Center, Scale, cvind, "global");
      end
   end
end
