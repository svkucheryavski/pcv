% This is demo code for the new version of the Procrustes
% cross-validation method (2022)
%
% In order to use the code you need to download and install the "Procrustes cross-validation"
% toolbox for MATLAB. You can do it from the public repository:
% https://se.mathworks.com/matlabcentral/fileexchange/121468-procrustes-cross-validation
% or by donwloading the .mtlbx file from GitHub (see Releases)
%
% The code below creates plots similar to what you can see on Figures with Corn example
% in the paper describing the new version, submitted in 2022.
%

% load Corn data
close all
clear
clc
load("corn.mat");

%% 1. PCA based examples

A = 20;
I = size(X, 1);

% create pseudo-validation set
Xpv = pcvpca(X, A, true, false, {"ven", 4});

% show plot with original and generated spectra
figure
subplot 221
plot(X')
title("Original data")

subplot 222
plot(Xpv')
title("PV-set")

subplot 223
plot((X - mean(X, 1))')
title("Original data (mean centered)")

subplot 224
plot((Xpv - mean(Xpv, 1))')
title("PV-set (mean centered)")


% mean center calibration and PV sets
mX = mean(X);
Xmc = X - mX;
Xpvmc = Xpv - mX;

% get PCA loadings for the calibration set
[~, ~, P] = svd(Xmc);
P = P(:, 1:A);

% compute eigenvalues - we need them for score distances
T = Xmc * P;
lambda = sum(T.^2) / (size(X, 1) - 1);

% create distance plot for A = 2 and A = 20 similar to shown in the paper
% but only for the new method (see the code for function 'plotdistances' at the
% end of this script)

figure
subplot 121
plotdistances(P, lambda, Xmc, Xpvmc, 2)
xlim([0, 20])
ylim([0, 12])
subplot 122
plotdistances(P, lambda, Xmc, Xpvmc, 20)
xlim([0,  6])
ylim([0, 12])


%% PCR based example

A = 20;

% create PV-set for PCR
Xpv = pcvpcr(X, Y, A, true, false, {"ven", 4});

% show plot with original and generated spectra
figure
subplot 221
plot(X')
title("Original data")

subplot 222
plot(Xpv')
title("PV-set")

subplot 223
plot((X - mean(X, 1))')
title("Original data (mean centered)")

subplot 224
plot((Xpv - mean(Xpv, 1))')
title("PV-set (mean centered)")

%% PLS based example

A = 20;

% create PV-set and also get matrix with scaling coefficients c.k/c (D)
[Xpv, D] = pcvpls(X, Y, A, true, false, {"ven", 4});

% show plot with original and generated spectra
figure
subplot 221
plot(X')
title("Original data")

subplot 222
plot(Xpv')
title("PV-set")

subplot 223
plot((X - mean(X, 1))')
title("Original data (mean centered)")

subplot 224
plot((Xpv - mean(Xpv, 1))')
title("PV-set (mean centered)")

% mean center X-values
mX = mean(X);
Xmc = X - mX;
Xpvmc = Xpv - mX;

% mean center Y-values
mY = mean(Y);
Ymc = Y - mY;

% compute global PLS model for mean centered X and Y
[P, R, C] = simpls(Xmc, Ymc, A);

% get scores for calibration and PV-sets
T = Xmc * R;
Tpv = Xpvmc * R;

% computed predicted values for each set using different number of LVs.
ypc = zeros(size(X, 1), A);
yppv = zeros(size(X, 1), A);
for a = 1:A
   ypc(:, a) = T(:, 1:a) * C(:, 1:a)';
   yppv(:, a) = Tpv(:, 1:a) * C(:, 1:a)';
end

% uncenter the predicted values
ypc = ypc + mY;
yppv = yppv + mY;

% compute errors of prediction
ec = Y - ypc;
epv = Y - yppv;

% compute RMSE
rmsec = sqrt(sum(ec.^2) / size(X, 1));
rmsepv = sqrt(sum(epv.^2) / size(X, 1));

% show plots with performance
figure

% predicted vs. measured values for A = 10
subplot 121
hold on
scatter(Y, ypc(:, 10))
scatter(Y, yppv(:, 10))
hold off
xlim([9.0, 11.5])
ylim([9.0, 11.5])
xlabel("Y, reference")
ylabel("Y, predicted")
title("Predictions (A = 10)")
grid()
box()
legend(["cal", "pv"])

% RMSE plot
subplot 122
hold on
plot(1:A, rmsec(:), '.-')
plot(1:A, rmsepv(:), '.-')
hold off
xlabel("Components")
ylabel("RMSE")
title("RMSE")
grid()
box()
legend(["cal", "pv"])

% show plot with scaling coefficients for K = 4
figure

heatmap(D, 'Colormap', parula, 'ColorLimits', [-1, 2])
xlabel("Components")
ylabel("Segments, k")


%% Helper functions

function [h, q] = getdistances(P, lambda, X, A)
%GETDISTANCES Compute score and orthogonal distance for PCA model and dataset X
   T = X * P(:, 1:A);
   E = X - T * P(:, 1:A)';
   q = sum(E.^2, 2);
   h = sum(T.^2 ./ lambda(1:A), 2);
end

function plotdistances(P, lambda, X, Xpv, A)
%PLOTDISTANCES Show plot with score and orthogonal distance for calibration and PV-set
   [h, q] = getdistances(P, lambda, X, A);
   [hpv, qpv] = getdistances(P, lambda, Xpv, A);

   hold on
   scatter(h / mean(h), q / mean(q), 's')
   scatter(hpv / mean(h), qpv / mean(q), 'x')

   hold off
   title(sprintf("Distance plot (A = %d)", A))
   xlabel("Score distance, h/h0")
   ylabel("Orthogonal distance, q/q0")
   grid()
   box()
   legend(["cal", "pv"])
end