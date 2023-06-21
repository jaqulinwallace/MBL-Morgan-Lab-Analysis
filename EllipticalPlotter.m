% -------------------------------------------------------------------------
% Elliptical plotting for declustering phenotypes.
% [JW 2023]
% -------------------------------------------------------------------------

clear all; clc; close all;

% -------------------------------------------------------------------------
control_all = readmatrix('');
low_all     = readmatrix('');
high_all    = readmatrix('');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Code starts here:

% Grab data:
AZcontrol = control_all(1,:)';
NNcontrol = control_all(2,:)';
AZlow = low_all(1,:)';
NNlow = low_all(2,:)';
AZhigh = high_all(1,:)';
NNhigh = high_all(2,:)';

% Take log of each vector:
lAZcontrol = log(AZcontrol);
lNNcontrol = log(NNcontrol);
lAZlow     = log(AZlow);
lNNlow     = log(NNlow);
lAZhigh    = log(AZhigh);
lNNhigh    = log(NNhigh);

% Plotting:
figure
subplot(1,3,1)
FitEllipse2Cloud([lAZcontrol,lNNcontrol])
hold on
xlim([2 8])
ylim([2 8])
title('Control')
xlabel('log Distance to AZ')
ylabel('log Distance to NN')
axis square

subplot(1,3,2)
FitEllipse2Cloud([lAZlow,lNNlow])
xlim([2 8])
ylim([2 8])
title('Low')
xlabel('log Distance to AZ')
ylabel('log Distance to NN')
axis square

subplot(1,3,3)
FitEllipse2Cloud([lAZhigh,lNNhigh])
xlim([2 8])
ylim([2 8])
title('High')
xlabel('log Distance to AZ')
ylabel('log Distance to NN')
axis square
hold off
% -------------------------------------------------------------------------