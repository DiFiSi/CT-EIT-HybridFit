%% Prepare workspace, load models and load EIDORS
cd(fileparts(matlab.desktop.editor.getActiveFilename));
addpath(genpath(".\fitting"));
dataFolder = ".\data";

%% Load data
load(fullfile(dataFolder,"sampleEIT.mat"));
fsEIT = 50; % Hz
tEIT = y(:,1);
yEIT = y(:,2);

load(fullfile(dataFolder,"sampleCT.mat"));
fsCT = 0.7; % Hz
tCT = y(:,1);
yCT = y(:,2);

%% Visualize data
figure;

subplot(1,2,1)

xStartEIT = 820;
xMaxEIT = 971;

plot(tEIT,yEIT); hold on;
plot(tEIT(xStartEIT),yEIT(xStartEIT),'r>');
plot(tEIT(xMaxEIT),yEIT(xMaxEIT),'ro');
xlabel("Time [s]");
ylabel("Amplitude [a.u.]");
title("EIT Sample");

subplot(1,2,2)

xStartCT = 3;
xMaxCT = 7;

plot(tCT,yCT); hold on;
plot(tCT(xStartCT),yCT(xStartCT),'r>');
plot(tCT(xMaxCT),yCT(xMaxCT),'ro');
xlabel("Time [s]");
ylabel("Amplitude [HU]");
title("CT Sample");
legend(["Data", "Start", "Peak"]);

%% Run nonlinear fit
disp = true; % turn on display of results

% EIT
resultNonlinEIT = nonlinFit(yEIT(xStartEIT:end), fsEIT, xStartEIT, xMaxEIT, disp);

% CT
resultNonlinCT = nonlinFit(yCT(xStartCT:end), fsCT, xStartCT, xMaxCT, disp);

%% Run hybrid fit
disp = true; % turn on display of results

% EIT
resultHybridEIT = hybridFit(yEIT(xStartEIT:end), fsEIT, xStartEIT, xMaxEIT, disp);

% CT
resultHybridCT = hybridFit(yCT(xStartCT:end), fsCT, xStartCT, xMaxCT, disp);
