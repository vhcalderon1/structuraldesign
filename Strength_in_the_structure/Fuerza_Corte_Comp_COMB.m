% ===========================================================================
% Structural Shear Force Comparison Analysis
% ===========================================================================
% Author: Victor Calderón (September 2022)
% Updated: Jefferson De la Cuba (February 2025)
% --------------------------------------------------------------------------
% Compares normalized shear distributions between conventional
% and isolated base structures. Processes .txt datasets from
% external directory and generates publication-quality figures.
%
% INPUTS:
% - Conventional_Court.txt   : Conventional base shear data
% - Court_test.txt         : Isolated base shear data
%   (Located in ../datasets relative to execution path)
%
% OUTPUTS:
% - DesignShearComparison.png: Output figure (saved to ../outputs)

%% Clear Workspace and Command Window
clear
close all
clc

%% Parameters
totalWeight = 3800.18; % Total weight in kN

%% Configure Paths
currentFolder = pwd;
datasetFolder = fullfile(currentFolder, '..', 'datasets');
outputFolder = fullfile(currentFolder, '..', 'outputs');

% Ensure output directory exists
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
end

%% Load and Process Data
% Conventional base system
conventionalData = load(fullfile(datasetFolder, 'Conventional_Court.txt'));
heightLevelsConventional = conventionalData(:,1); % Elevation [m]
normalizedShearConventional = conventionalData(:,2)/totalWeight * 100; % V/W [%]

% Isolated base system
isolatedData = load(fullfile(datasetFolder, 'Court_test.txt'));
heightLevelsIsolated = isolatedData(:,1); % Elevation [m]
normalizedShearIsolated = isolatedData(:,2)/totalWeight * 100; % V/W [%]

%% Generate Comparison Figure
figureHandle = figure('Units','centimeters','Position',[0 0 10 12]);

% Create stairs plots
stairs(normalizedShearConventional, heightLevelsConventional,...
    'k--', 'LineWidth', 2);
hold on
stairs(normalizedShearIsolated, heightLevelsIsolated,...
    'k-', 'LineWidth', 2);

%% Figure Formatting
% Axis configuration
xlabel('Normalized Design Shear, V/W (%)',...
    'FontSize',13, 'FontName','Times New Roman');
ylabel('Height Level (m)',...
    'FontSize',13, 'FontName','Times New Roman');
xlim([0 25]); ylim([0 8]);
xticks(0:5:25); yticks(0:1:8);

% Style adjustments
ax = gca;
ax.FontName = 'Times New Roman';
ax.FontSize = 13;
ax.LineWidth = 1.3;
box off

% Legend configuration
legend({'Fixed base', 'Isolated base'},...
    'Location','best',...
    'FontName','Times New Roman',...
    'Box','off');

%% Save and Close Figure
outputFile = fullfile(outputFolder, 'DesignShearComparison.png');
print(figureHandle, outputFile, '-dpng', '-r600');
close(gcf); % Close figure after saving
disp(['Figure saved to: ' outputFile]);