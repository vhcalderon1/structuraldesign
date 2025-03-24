% ==============================================================
% Bilinear Concrete Compression Stress-Strain Model Generator
% ==============================================================
% Author: Victor Calderón (November 2020)
% Updated: Jefferson De la Cuba (February 2025)
%------------------------------------------------------------------------------
% This code generates a bilinear stress-strain model for concrete under compression
% based on ASTM standards using SI units. The model consists of an ascending parabolic
% branch and linear descending branch, with all calculations in kN/m².
%
% Inputs:
%   - Compressive strength (kN/m²)
%   - Peak stress (90% of compressive strength)
%   - Concrete modulus (calculated from compressive strength)
%   - Compressive strength [kN/m²] (e.g., 21000)
%   - Peak stress factor (optional, default = 0.9)
%   - Ultimate strain (optional, default = 0.0038)n
%
% Outputs:
%   - stress_curve (kN/m²), strain_curve : Data vectors
%   - ../outputs/Concrete_Stress_Strain_Model.png : High-res plot (600 DPI)
%   - ../datasets/ : Input/output data directory
%
%% Clean Workspace and Verify File Structure
clear
close all
clc

% Configure directories
outputDir = '../outputs';
datasetsDir = '../datasets';

% Create directories if missing
if ~exist(datasetsDir, 'dir')
    mkdir(datasetsDir)
    fprintf('Created datasets directory: %s\n', datasetsDir)
end

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
    fprintf('Created outputs directory: %s\n', outputDir)
end

%% ==================== Material Properties ===========================
compressive_strength = 21000;          % Concrete strength [kN/m²]
peak_stress = compressive_strength * 0.9; % 90% of f'_c
concrete_modulus = 149535 * sqrt(compressive_strength); % Modulus [kN/m²]

%% ==================== Stress-Strain Calculations ====================
[total_stress, total_strain] = calculate_concrete_model(...
    compressive_strength, peak_stress, concrete_modulus);

%% ==================== Visualization & Export ========================
configure_plot(total_strain, total_stress);
save_figure(fullfile(outputDir, 'Concrete_Stress_Strain_Model'));
close(gcf);

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
    
    % Main plot configuration
    plot(strain, stress, 'k', 'LineWidth', 2)
    hold on
    plot([max(strain) max(strain)], [0 max(stress)], '--k', 'LineWidth', 1.5)
    
    % Axis settings
    ax = gca;
    ax.XAxis.Exponent = -3;
    ax.YLim = [-1000 25000];
    ax.XLim = [-0.1e-3 4e-3];
    ax.XTick = 0:2e-3:4e-3;
    ax.YTick = 0:5000:25000;
    
    % Labels with LaTeX formatting
    xlabel('Concrete Strain, $\varepsilon_c$', 'Interpreter', 'latex', ...
        'FontSize', 13, 'FontName', 'Times New Roman')
    ylabel('Stress, $\sigma_c$ (kN/m$^2$)', 'Interpreter', 'latex', ...
        'FontSize', 13, 'FontName', 'Times New Roman')
    
    % Legend configuration
    legend({'Concrete Model', 'Ultimate Strain'}, ...
        'Location', 'northeast', 'FontName', 'Times New Roman')
    legend('boxoff')
    set(ax, 'LineWidth',1.3, 'FontSize',12, 'FontName','Times New Roman')
end

function save_figure(fname)
    % Save figure as PNG with 600 DPI
    set(gcf, 'PaperUnits', 'centimeters', ...
             'PaperSize', [12 10], ...
             'PaperPositionMode', 'manual', ...
             'PaperPosition', [0 0 12 10])
    
    print([fname '.png'], '-dpng', '-r600');
    fprintf('Saved figure: %s.png\n', fname);
end