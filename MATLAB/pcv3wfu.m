function Xpv = pcv3wfu(X, nComp, Center, Scale, CV, Seed)
%PCV3WFU Procrustes cross-validation for three-way data based on full unfolding (FU)
%
% Arguments:
% ----------
% X         3-D array with calibration set, size I x J x K (samples x mode 2 x mode 3)
% nComp     number of PCA components (default: min(I - 1, J * K, 30))
% Center    logical, center the unfolded data or not (default: true)
% Scale     logical, scale the unfolded data or not, each of the J * K columns
%           separately (default: false)
% CV        defines cross-validation splits, see PCVPCA (default: {"ven", 4})
% Seed      optional number, makes random splits and random orthogonal parts
%           reproducible; the state of the random generator is restored afterwards
%
% Returns:
% --------
% Xpv       3-D array with PV-set (same size as X)
%
% The method unfolds X to a matrix with I rows and J * K columns, generates PV-set
% for the matrix using PCVPCA with global CV scope, and folds the result back to a
% 3-D array. The PV-set reproduces cross-validation of the PCA model of the unfolded
% data. Use enough components, so the PCA subspace contains the subspace of the
% three-way models (e.g. PARAFAC or Tucker3) you are going to validate.
%
% Use [] for an argument to keep its default value.
%
% See also: PCV3WBU, PCV3WHSV, PCVPCA

   t = pcv3wtools();
   d = t.getdim(X);

   if nargin < 2 || isempty(nComp), nComp = min([d(1) - 1, d(2) * d(3), 30]); end
   if nargin < 3 || isempty(Center), Center = true; end
   if nargin < 4 || isempty(Scale), Scale = false; end
   if nargin < 5 || isempty(CV), CV = {"ven", 4}; end
   if nargin < 6, Seed = []; end

   cleanup = t.setseed(Seed); %#ok<NASGU>

   % unfold to I x (J*K) with mode 2 running fastest, generate PV-set and fold back
   cvind = t.splitcv(CV, d(1));
   Xpv = pcvpca(reshape(double(X), d(1), d(2) * d(3)), nComp, Center, Scale, cvind, "global");
   Xpv = reshape(Xpv, d);
end
