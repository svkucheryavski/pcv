close all
clear
clc

% load original spectra and compute pseudo-validation set
load nirsim
spectra_pv = pcv(spectra, 6, 4);

% show the both sets as line plots
figure
subplot 211
plot(wavenumbers, spectra)
subplot 212
plot(wavenumbers, spectra_pv)

% make PCA model for original spectra


mX = mean(spectra);
Xc = bsxfun(@minus, spectra, mX);
Xpvc = bsxfun(@minus, spectra_pv, mX);

P = pca(Xc);
Tc = Xc * P;
Tpv = Xpvc * P;

A = 4;
Ec = Xc - Tc * P(:, 1:A)';
Epv = Xpvc - Tpv * P(:, 1:A)';
Qc = sum(Ec.^2, 2);
Qpv = sum(Epv.^2, 2);



