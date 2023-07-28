function [P, R, C] = simpls(X, Y, nComp)
%SIMPLS Implementation of SIMPLS algorithm for PLS regression
%  X - 
   if size(Y, 1) < size(Y, 2)
      Y = Y';
   end

   nObj = size(X, 1);
   nPred = size(X, 2);
   nResp = size(Y, 2);

   S = X' * Y;
   M = X' * X;

   % y-loadings   
   C = zeros(nResp, nComp);
   
   % x-loadings and weights
   R = zeros(nPred, nComp);
   V = zeros(nPred, nComp);
   P = zeros(nPred, nComp);
   
   % scores
   T = zeros(nObj, nComp);
   U = zeros(nObj, nComp);
   
   for a = 1:nComp
      [SU, ~, ~] = svd(S, 'econ');
      r = SU(:, 1);
      t = X * r;

      tnorm = sqrt(sum(t .* t));
      t = t / tnorm;
      r = r / tnorm;

      p = X' * t;
      c = Y' * t;
      u = Y * c;
      v = p;

      if a > 1
         v = v - V * (V' * p);
         u = u - T * (T' * u);
      end

      v = v / sqrt(sum(v .* v));

      R(:, a) = r;
      V(:, a) = v;
      P(:, a) = p;
      T(:, a) = t;
      U(:, a) = u;
      C(:, a) = c;

      M = M - p' * p;
      S = S - v * (v' * S);
   end   
end
