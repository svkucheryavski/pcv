% This is demo code for the new version of the Procrustes 
% cross-validation method (2022)
%
% The code creates plots similar to what you can see on Figures with Corn example
% in the paper describing the new version, submitted in 2022. The Corn dataset was
% downloaded from here: https://eigenvector.com/resources/data-sets/
%

% load Corn data
close all
clear
clc
load("corn.mat");

%% 1. PCA based examples

A = 30;
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

% make PCA model for the calibration set
[~, ~, P] = svd(Xmc);
P = P(:, 1:A);
T = Xmc * P;
lambda = sum(T.^2) / (I - 1);

% get scores for the PV-set
Tpv = Xpvmc * P;

% compute scores and distances for both calibration set and PV-set 
% for a = 1 to A

q = zeros(I, A);
h = zeros(I, A);

qpv = zeros(I, A);
hpv = zeros(I, A);

for a = [2, 20]
    E = Xmc - T(:, 1:a) * P(:, 1:a)'; 
    Epv = Xpvmc - Tpv(:, 1:a) * P(:, 1:a)'; 

    q(:, a) = sum(E.^2, 2);
    qpv(:, a) = sum(Epv.^2, 2);

    h(:, a) = sum(T(:, 1:a).^2 ./ lambda(1:a), 2);
    hpv(:, a) = sum(Tpv(:, 1:a).^2 ./ lambda(1:a), 2);
end

% create distance plot for A = 2 and A = 20 similar to shown in the paper
% but only for the new method
figure

a = 2;
subplot 121
hold on
scatter(h(:, a) / mean(h(:, a)), q(:, a) / mean(q(:, a)), 's')
scatter(hpv(:, a) / mean(h(:, a)), qpv(:, a) / mean(q(:, a)), 'x')
hold off
xlim([0, 20])
ylim([0, 12])
title("Distance plot (A = 2)")
xlabel("Score distance, h/h0")
ylabel("Orthogonal distance, q/q0")
grid()
box()
legend(["cal", "pv"])

a = 20;
subplot 122
hold on
scatter(h(:, a) / mean(h(:, a)), q(:, a) / mean(q(:, a)), 's')
scatter(hpv(:, a) / mean(h(:, a)), qpv(:, a) / mean(q(:, a)), 'x')
hold off
xlim([0,  6])
ylim([0, 12])
title("Distance plot (A = 20)")
xlabel("Score distance, h/h0")
ylabel("Orthogonal distance, q/q0")
grid()
box()
legend(["cal", "pv"])

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
