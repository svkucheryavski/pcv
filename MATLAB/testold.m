% Simple tests for PCV implementation
clear
clc

X = randn(100, 50);

% different sets of arguments
args = {{1}, {10}, {20}, {10, 2}, {10, 4}, {10, 100}, {10, 2, true}, {10, 4, true}, {10, 100, true}};

fprintf("Running: ")
for i = 1:numel(args)
   fprintf("%d ", i)
   Xpv = pcv(X, args{i}{:});
   if sum((Xpv(:) - X(:)).^2) < 10^-6; error("Sets seem to be identical."); end
   if size(Xpv) ~= size(X); error("Wrong dimension."); end
   for j = 1:size(X, 2)
      [h, p] = kstest2(Xpv(:, j), X(:, j));
      if p < 0.01; error("Populations are not the same."); end
   end
end
fprintf("\n")