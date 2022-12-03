function [Xpv, D] = pcvpcr(X, Y, nComp, Center, Scale, CV)
%PCVPCR Compute matrix with pseudo-validation set for PCR model
%
% Arguments
% ---------
% X         matrix with predictors from calibration set
% Y         matrix with responses from calibration set
% nComp     number of components for decomposition
% Center    logical, mean center data columns or not
% Scale     logical, standardize data columns or not
% CV        defines cross-validation splits pattern (see description)
%
% Returns:
% --------
% Xpv       matrix with PV-set (same size as X)
% D         matrix with coefficients (ck/c) for each PC and each segment
%
% The method computes pseudo-validation matrix Xpv, based on PCR decomposition of 
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
% See also: PCVPCA, PCVPLS
%
    
    % computes a global PCR model
    function m = getglobalmodel(X, Y, nComp) 
         [~, ~, P] = svd(X, 'econ');
         P = P(:, 1:nComp);
         Pi = eye(size(X, 2)) - P * P';
         T = X * P;
         C = (T' * T)' \ (T' * Y);
         m = struct('P', P, 'C', C, 'Pi',  Pi, 'nComp', nComp);
    end

    % computes a local PCR model
    function mk = getlocalmodel(Xc, Yc, m) 
         [~, ~, Pk] = svd(Xc, 'econ');
         Pk = Pk(:, 1:m.nComp);
         aa = acos(sum(m.P .* Pk)) < (pi / 2);
         Pk = Pk * diag(aa * 2 - 1);
         Tc = Xc * Pk;
         Ck = (Tc' * Tc)' \ (Tc' * Yc);

         mk = struct('P', Pk, 'C', Ck); 
    end

    % computes PV-set for current segment
    function [Xpv, Dk] = getxpv(m, mk, Xk) 

       Tk = Xk * mk.P;
       Dk = diag(mk.C(:) ./ m.C(:));
       Tpv = Tk * Dk;
       Xpv = Tpv * m.P';
    end

    % computes vector with orthogonal distances
    function qk = getqk(Xk, mk) 
         qk = sum( (Xk - Xk * mk.P * mk.P').^2, 2);
    end
    
    funlist = struct(...
        'getglobalmodel', @getglobalmodel, ...
        'getlocalmodel', @getlocalmodel,...
        'getxpv', @getxpv,...
        'getqk', @getqk...
    );

    [Xpv, D] = pcvreg(X, Y, nComp, Center, Scale, CV, funlist);
end

