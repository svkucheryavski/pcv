%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                   #
%  This file contains demo code for implementation of Procrustes cross-validation   #
%  method in MATLAB. The code creates plots similar to what you can see on figures  #
%  with Corn example in this paper: https://doi.org/10.1016/j.aca.2023.341096       #
%                                                                                   #
%  Check this for more details: https://github.com/svkucheryavski/pcv/              #
%                                                                                   #
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% !!! Make sure you install the toolbox before running the code !!! %


% load Corn data
close all
clear
clc
load("corn.mat");


%% 1. PCA based examples

A = 20;
I = size(X, 1);

% create pseudo-validation set using PCA based algorithm
Xpv = pcvpca(X, A, true, false, {"ven", 4});


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


% show distance plots for A = 2 and A = 20
% similar to ones shown in Supplementary materials
% see the code for function 'plotdistances' at the end of this script

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
% similar to plots shown in Figure 2
figure

subplot 121
plot((X - mean(X, 1))')
title("Original data (mean centered)")

subplot 122
plot((Xpv - mean(Xpv, 1))')
title("PV-set (mean centered)")


%% PLS based example

A = 20;

% create PV-set and also get matrix with scaling coefficients c.k/c (D)
[Xpv, D] = pcvpls(X, Y, A, true, false, {"ven", 4});

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

% show RMSE and predicted vs. measured plots
% similar to plots shown in Figure 5
figure

% RMSE plot
subplot 131
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

% predicted vs. measured values for A = 10
subplot 132
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

% predicted vs. measured values for A = 20
subplot 133
hold on
scatter(Y, ypc(:, 20))
scatter(Y, yppv(:, 20))
hold off
xlim([9.0, 11.5])
ylim([9.0, 11.5])
xlabel("Y, reference")
ylabel("Y, predicted")
title("Predictions (A = 20)")
grid()
box()
legend(["cal", "pv"])


% show heatmap and boxplot for elements of D
% similar to the first plot shown in Figure 6
figure
heatmap(D, 'Colormap', parula, 'ColorLimits', [-1, 2]);
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