function xpvorth = getxpvorth(qk, Xk, PRM)
%XPVORTH Computes orthogonal part of PV-set
%
   nObj = numel(qk);

   % compute the orthogonal component of Xpv
   xpvorth = rand(nObj) * Xk;                                % - project Xk to a random vector
   xpvorth = xpvorth * diag(1 ./ sqrt(sum(xpvorth.^2)));     % - normalize columns
   xpvorth = xpvorth * PRM;                                  % - orthogonalize to global component space
   xpvorth = diag(sqrt(qk ./ sum(xpvorth.^2, 2))) * xpvorth; % - rescale rows
end
