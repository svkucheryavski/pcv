% This is demo code for Procrustes cross-validation method
%
% The code creates a plot similar to what you can see on Figure 5
% of the original paper (bottom plot in the middle column)
%

% define how many components to use
A = 4;

% load original spectra and compute the pseudo-validation set (A = 6 and K = 4)
load nirsim
spectra_pv = pcv(spectra, 6, 4);

% show the both sets as line plots
figure
subplot 211
plot(wavenumbers, spectra)
subplot 212
plot(wavenumbers, spectra_pv)

% mean center the original spectra and the pv-set
mX = mean(spectra);
Xc = bsxfun(@minus, spectra, mX);
Xpvc = bsxfun(@minus, spectra_pv, mX);

% create a PCA model using mean centered original spectra
[~, ~, P] = svd(Xc, 'econ');
P = P(:, 1:A);

% apply the model to get scores for both sets
Tc = Xc * P;
Tpv = Xpvc * P;

% compute orthogonal distance for given number of components
Ec = Xc - Tc * P';
Epv = Xpvc - Tpv * P';
qc = sum(Ec.^2, 2);
qpv = sum(Epv.^2, 2);

% compute score distance for given number of components
L = sum(Tc.^2) / (size(Tc, 1) - 1);
Uc = Tc * diag(1./sqrt(L));
Upv = Tpv * diag(1./sqrt(L));
hc = sum(Uc.^2, 2);
hpv = sum(Upv.^2, 2);

% compute parameters for the distance distribution
q0 = mean(qc);
h0 = mean(hc);
Nq = round(2 * (q0 / std(qc))^2);
Nh = round(2 * (h0 / std(hc))^2);

% compute total distance, f, and corresponding p-values
fc = qc / q0 * Nq + hc / h0 * Nh;
Nf = Nq + Nh;
pc = cdf('chi2', fc, Nf);
fpv = qpv / q0 * Nq + hpv / h0 * Nh;
ppv = cdf('chi2', fpv, Nf);

% compute number of extremes for different alpha
N = numel(qc);
nextrc = zeros(N, 1);
nextrpv = zeros(N, 1);
nseq = 1:N;
alpha = nseq / N;

for n = nseq
   nextrc(N - n + 1) = sum(pc >= alpha(n));
   nextrpv(N - n + 1) = sum(ppv >= alpha(n));
end

% compute coordinates of the tolerance ellipse
D = 2 * sqrt(nseq .* (1 - alpha));
Nm = [0, nseq - D];
Np = [0, nseq + D];

% now we create the extreme plot
figure
hold on

% show the tolerance ellipse
plot([0, N], [0 N], 'Color', [0.5 0.6 0.7])
plot([0, nseq], Nm, 'Color', [0.5 0.6 0.7])
plot([0, nseq], Np, 'Color', [0.5 0.6 0.7])
plot([nseq; nseq], [Nm(2:end); Np(2:end)], 'Color', [0.6 0.7 0.8])

% show the extremes as points
h1 = scatter(nseq, nextrc, 'b');
h2 = scatter(nseq, nextrpv, 'r');
hold off

% legend and other stuff
box on
grid on
xlabel("Expected")
ylabel("Observed")
title(sprintf("Extremes, A = %d",A))
xlim([0 N])
ylim([0 N])
legend([h1, h2], {'cal', 'pv'}, 'Location', 'SouthEast')