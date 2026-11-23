function xpvorth = getxpvorth(qk, Xk, PRM)
%XPVORTH Computes orthogonal part of PV-set
%
   nObj = numel(qk);

   % compute the orthogonal component of Xpv
   xpvorth = rand(nObj) * Xk;                                % - project Xk to a random vector
   colnorms = sqrt(sum(xpvorth.^2));
   colnorms(colnorms < eps) = 1;
   xpvorth = xpvorth * diag(1 ./ colnorms);                  % - normalize columns
   xpvorth = xpvorth * PRM;                                  % - orthogonalize to global component space
   rownorms = sqrt(sum(xpvorth.^2, 2));
   rownorms(rownorms < eps) = 1;
   xpvorth = diag(sqrt(qk) ./ rownorms) * xpvorth;           % - rescale rows
end
