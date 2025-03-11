% ==========================================================
% Moment-Curvature Relationship Analysis and Visualization
% ==========================================================
% Author: Victor Calderón (November 2020)
% Updated: Jefferson De la Cuba (February 2025)
% ----------------------------------------------------------------------------
% This script processes structural engineering data to analyze and visualize 
% the moment-curvature relationship of a beam cross-section. It includes unit 
% conversions, plotting with publication-quality formatting, and figure export.
%
% Inputs: 
%   - Dataset File: 'Curvature_Moment.txt' in ../datasets/
%     Columns: [Curvature (1/cm), Moment (kN-m)]
% Outputs:
%   - Figure: 'Moment_Curvature_Plot.eps' in ../outputs/
%     Formatted plot showing curvature (1/cm) vs. moment (kN-m)
%
% -------------------------------------------------------------------------

clear;
close all;
clc;

%% Define Input and Output Directories
currentFolder = fileparts(mfilename('fullpath'));
datasetFolder = fullfile(currentFolder, '..', 'datasets');
outputFolder = fullfile(currentFolder, '..', 'outputs');

% Create outputs directory if non-existent
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Load Data
inputFile = fullfile(datasetFolder, 'Momento_Curvatura.txt');
data = load(inputFile);

% Extract curvature and moment data
curvature = data(:,1) / 100; % Convert from 1/cm
moment = data(:,2); % Moment in kN-m

%% Plot Curvature vs. Moment
figure;
plot(curvature, moment, 'k', 'LineWidth', 2);

% Set axis labels
xlabel('Curvature, \Phi (1/cm)', 'FontSize', 13, 'FontName', 'Times New Roman');
ylabel('Moment, M (kN-m)', 'FontSize', 13, 'FontName', 'Times New Roman');

% Get current axes handle
ax = gca;

% Customize axis properties
set(ax, 'LineWidth', 1.5);
box off;
axis([-0.1e-4 2e-4 -20 400]);
xticks(0:1e-4:2e-4);
yticks(0:100:400);

% Set axis exponents
ax.YAxis.Exponent = 2;
ax.XAxis.Exponent = -4;
set(ax, 'LineWidth', 1.3, 'FontSize', 12, 'FontName', 'Times New Roman');

%% Set Figure Size for Export
width = 12; % cm
height = 10; % cm
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [width height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 width height]);

%% Save Output Figure
outputFile = fullfile(outputFolder, 'Moment_Curvature_Plot.eps');
print(outputFile, '-depsc', '-tiff');
disp(['Plot saved at: ', outputFile]);
close(gcf)