% ==============================================================
% Bilinear Concrete Compression Stress-Strain Model Generator
% ==============================================================
% Author: Victor Calderón (November 2020)
% Updated: Jefferson De la Cuba (February 2025)
%------------------------------------------------------------------------------
% This code generates a bilinear stress-strain model for concrete under compression
% based on ASTM standards. The model consists of an ascending parabolic branch
% and a linear descending branch. The inputs include concrete compressive strength
% in kg/cm², and the outputs are stress and strain vectors along with visualization.
%
% Inputs:
%   - Compressive strength (kg/cm²)
%   - Peak stress (90% of compressive strength)
%   - Concrete modulus (calculated from compressive strength)
%   - material_params.txt : Contains:
%                            1. Compressive strength [kg/cm²] (e.g., 210)
%                            2. Peak stress factor (optional, default = 0.9)
%                            3. Ultimate strain (optional, default = 0.0038)
%   - strain_control.txt  : (Optional) Custom strain range/resolution for curve generation
%
% Outputs:
%   - stress_curve, strain_curve : Data vectors in MATLAB workspace
%   - ../outputs/Concrete_Stress_Strain_Model.eps : Vector graphic (300 DPI)
%   - ../outputs/Concrete_Stress_Strain_Model.tif  : Raster graphic (300 DPI)
%   - ../datasets/ : Directory for input/output data archival
%
%% Clean Workspace and Verify File Structure
clear
close all
clc

% Configure output directory using path traversal
outputDir = '../outputs';
datasetsDir = '../datasets';

% Verify/Create required directories
if ~exist(datasetsDir, 'dir')
    mkdir(datasetsDir)
    fprintf('Created datasets directory: %s\n', datasetsDir)
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
    fprintf('Created outputs directory: %s\n', outputDir)
end

%% ==================== Material Properties ===========================
compressive_strength_kgcm2 = 210;          % Concrete strength [kg/cm²]
peak_stress_kgcm2 = compressive_strength_kgcm2 * 0.9; % 90% of f'_c
concrete_modulus_kgcm2 = round(15100*sqrt(compressive_strength_kgcm2), -2);

%% ==================== Stress-Strain Calculations ====================
[total_stress, total_strain] = calculate_concrete_model(...
    compressive_strength_kgcm2, peak_stress_kgcm2, concrete_modulus_kgcm2);

%% ==================== Visualization & Export ========================
configure_plot(total_strain, total_stress);
save_figure(fullfile(outputDir, 'Concrete_Stress_Strain_Model'));
close(gcf); % Close figure after saving

%% ====================================================================
% Supporting Functions
% =====================================================================

function [total_stress, total_strain] = calculate_concrete_model(~, fp, Ec)
    % Ascending branch calculations
    strain_peak = 1.8 * fp / Ec;
    strain_asc = linspace(0, strain_peak, 100)';
    stress_asc = fp * (2*(strain_asc/strain_peak) - (strain_asc/strain_peak).^2);
    
    % Descending branch calculations
    strain_ult = 0.0038;
    slope = -0.15 * fp / (strain_ult - strain_peak);
    strain_dsc = linspace(strain_peak, strain_ult, 100)';
    stress_dsc = fp + slope * (strain_dsc - strain_peak);
    
    % Combine results
    total_stress = [stress_asc; stress_dsc];
    total_strain = [strain_asc; strain_dsc];
end

function configure_plot(strain, stress)
    figure('Units','centimeters','Position',[0 0 12 10])
    
    % Main plot elements
    plot(strain, stress, 'k', 'LineWidth', 2)
    hold on
    plot([max(strain) max(strain)], [0 max(stress)], '--k', 'LineWidth', 1.5)
    
    % Axis configuration
    ax = gca;
    ax.XAxis.Exponent = -3;
    ax.YLim = [-10 200];
    ax.XLim = [-0.1e-3 4e-3];
    ax.XTick = 0:2e-3:4e-3;
    ax.YTick = 0:50:200;
    
    % Latex-style labels
    xlabel('Concrete Strain, $\varepsilon_c$', 'Interpreter', 'latex', ...
        'FontSize', 13, 'FontName', 'Times New Roman')
    ylabel('Stress, $\sigma_c$ (kg/cm$^2$)', 'Interpreter', 'latex', ...
        'FontSize', 13, 'FontName', 'Times New Roman')
    
    % Legend and grid
    legend({'Concrete Model', 'Ultimate Strain'}, ...
        'Location', 'northeast', 'FontName', 'Times New Roman')
    legend('boxoff')
    set(ax, 'LineWidth',1.3, 'FontSize',12, 'FontName','Times New Roman')
end

function save_figure(fname)
    set(gcf, 'PaperUnits', 'centimeters', ...
             'PaperSize', [12 10], ...
             'PaperPositionMode', 'manual', ...
             'PaperPosition', [0 0 12 10])
    
    % Save as EPS and TIFF
    print(fname, '-depsc', '-r300');
    print(fname, '-dtiff', '-r300');
    
    fprintf('Saved figures:\n  %s.eps\n  %s.tif\n', fname, fname)
end