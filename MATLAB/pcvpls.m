function [Xpv, D] = pcvpls(X, Y, nComp, Center, Scale, CV, CVScope)
%PCVPLS Compute matrix with pseudo-validation set for PLS model
%
% Arguments
% ---------
% X         matrix with predictors from calibration set
% Y         matrix with responses from calibration set
% nComp     number of components for decomposition
% Center    logical, mean center data columns or not
% Scale     logical, standardize data columns or not
% CV        defines cross-validation splits pattern (see description)
% CVScope   defines if centering/scaling must be done globally ('global')
% or locally ('local'). Default value is 'global'.
%
% Returns:
% --------
% Xpv       matrix with PV-set (same size as X)
% D         matrix with coefficients (ck/c) for each PC and each segment
%
% The method computes pseudo-validation matrix Xpv, based on PLS decomposition of
% calibration set (X, Y) and cross-validation. The 'CV' parameter can be:
%
% A number:
% =========
% 'CV = 1' - full cross-validation (leave one out)
% 'CV = n' - random cross-validation with n segments
%
% A cell array:
% =============
% 'CV = {"rand", n}' - random cross-validation with n segments
% 'CV = {"ven", n}' - systematic cross-validation with n segments
% 'CV = {"loo"}' - full cross-validation (leave one out)
%
% A vector:
% =========
% 'CV = [1, 2, 3, 1, 2, 3] - cross-validation of 6 rows with 3 segments
%
% See more details here: https://github.com/svkucheryavski/pcv
%
% See also: PCVPCA, PCVPCR
%

   % default values
   if nargin < 7
      CVScope = "global";
   end

   if nargin < 6
      CV = {"ven", 4};
   end

   if nargin < 5
      Scale = false;
   end

   if nargin < 4
      Center = true;
   end
   
   if nargin < 3
      nRows = size(X, 1);
      nCols = size(X, 2);
      nComp = min([nRows - 1, nCols, 30]);
   end 
   
   
   % computes a global PLS model
   function m = getglobalmodel(X, Y, nComp)
      [P, R, C] = simpls(X, Y, nComp);
      Pi = eye(size(X, 2)) - R * P';
      m = struct('P', P, 'R', R, 'C', C, 'Pi', Pi, 'nComp', nComp);
   end

   % computes a local PLS model
   function mk = getlocalmodel(Xc, Yc, m)
      [Pk, Rk, Ck] = simpls(Xc, Yc, m.nComp);

      normalize = @(X) X .* (1 ./ sqrt(sum(X.^2)));
      aa = acos(sum((normalize(m.R)) .* normalize(Rk))) < (pi / 2);
      Pk = Pk * diag(aa * 2 - 1);
      Rk = Rk * diag(aa * 2 - 1);
      Ck = Ck * diag(aa * 2 - 1);

      mk = struct('P', Pk, 'R', Rk, 'C', Ck);
   end

   % computes PV-set for current segment
   function [Xpv, Dk] = getxpv(m, mk, Xk)
      Tk = Xk * mk.R;
      Dk = eye(m.nComp);
      for a = 1:m.nComp
         Dk(a, a) = (mk.C(:, a)' * m.C(:, a)) ./ (m.C(:, a)' * m.C(:, a));
      end
      Tpv = Tk * Dk;
      Xpv = Tpv * m.P';
   end

   % computes vector with orthogonal distances
   function qk = getqk(Xk, mk)
      qk = sum( (Xk - Xk * mk.R * mk.P').^2, 2);
   end

   funlist = struct(...
      'getglobalmodel', @getglobalmodel, ...
      'getlocalmodel', @getlocalmodel,...
      'getxpv', @getxpv,...
      'getqk', @getqk...
      );

   [Xpv, D] = pcvreg(X, Y, nComp, Center, Scale, CV, funlist, CVScope);
end
