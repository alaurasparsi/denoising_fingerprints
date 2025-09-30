%!/usr/bin/env matlab
% Author: Laura Belli, Sara Stampacchia

%% behaviour's prediction with FC edges selected from >99.9th percentile ICC
clc;
clear;
close all;

%% Initialize paths and set constants
Path2Data = fullfile('..', 'results\');
Path2Results = fullfile('..', 'results\');
Path2Fig = fullfile('..', 'figures\');
Path2Yeo = fullfile('..', 'parcs\shen');
load(fullfile(Path2Yeo, 'shen_yeo_RS7.mat'));
Path2functions = fullfile('..', 'functions\');
load(fullfile(Path2Parc, '\shen_yeo_RS7.mat'));
T = load(fullfile(Path2results,'matrix_Simon.mat')); % ses x sub
matrix_Simon = T.matrix_Simon;
FC_3D = FC_3D_mat_all;
N = size(parc.yeoROIs, 1);
mask_ut = triu(true(N), 1);
numEdges = 38503;
n_sub=10;
n_ses=10;


%% CASE 1: within-session ICC prediction
% Load data and initialize paths
load(fullfile(Path2results, "ICC_2D_within.mat"));
ICC_2D = ICC_within_2D_all;
percentile_threshold = prctile(ICC_2D, 99.9);
temp_ICC_mat_highFG = NaN(N, N);
ICC_mask_highFG_2D = NaN(numEdges, n_ses);
ICC_mat_highFG = zeros(N, N);
FC_sumabs_within = NaN(n_sub, n_ses);

% Build predictor: sum|FC values| per session
for ses = 1:n_ses  
    temp_ICC = ICC_2D(:, ses);
    mask_highFG_2D = temp_ICC >= percentile_threshold(:, ses);
    
    for sub = 1:n_sub
        subject_FC = FC_3D(:, ses, sub);
        masked_FC = subject_FC .* mask_highFG_2D; 
        FC_sum = sum(abs(masked_FC), 'omitnan');
        FC_sumabs_within(sub, ses) = FC_sum;
    end

    % plot ICC top edges 
    ICC_vec_highFG = temp_ICC.* mask_highFG_2D;
    ICC_mat_highFG(mask_ut) = ICC_vec_highFG;
    ICC_mat_within = ICC_mat_highFG + ICC_mat_highFG';
    figure; imagesc(ICC_mat_within(parc.yeoOrder, parc.yeoOrder)); 
    top_title = sprintf("Within-session ICC edges > 99.9th percentile, session %d", ses);
    axis square; colorbar; title(top_title); 

    % Check if NaN and save figure
    nnz(isnan(ICC_mat_highFG));
    save_filename = sprintf('ICC_within_highFG_ses-0%d.png', ses);
    saveas(gcf, fullfile(Path2Fig, save_filename));
    ICC_mask_highFG_2D(:, ses) = ICC_vec_highFG;
    close;
end 

% Save
save(fullfile(Path2results, "ICC_mask_highFG_2D_within.mat"), "ICC_mask_highFG_2D");
save(fullfile(Path2results, 'FC_sumabs_within.mat'), 'FC_sumabs_within');
writematrix(FC_sumabs_within, 'FC_sumabs_within.csv');


%% Plot within-session ICC edges > 95th percentile for visualisation
percentile_threshold = prctile(ICC_2D, 95);
temp_ICC_mat_highFG = NaN(N, N);
ICC_mask_highFG_2D = NaN(numEdges, n_ses);
ICC_mat_highFG = zeros(N, N);

figure('Color', 'w'); set(gcf,'Position',[100 100 1100 800]);
sgtitle('Within-session ICC edges > 95th percentile');
% Get the mask for every session
for ses = 1:n_ses
    temp_ICC_2D = ICC_2D(:, ses);
    mask_highFG_2D = temp_ICC_2D >= percentile_threshold(:, ses);

    % Display top delta ICC edges
    ICC_vec_highFG = temp_ICC_2D.* mask_highFG_2D;
    ICC_mat_highFG(mask_ut) = ICC_vec_highFG';
    ICC_mat_within = ICC_mat_highFG + ICC_mat_highFG';
 
    subplot(4, 3, ses);
    imagesc(ICC_mat_within(parc.yeoOrder, parc.yeoOrder));
    axis square; colormap;
    set(gca, 'XTick', [], 'YTick', []);
    ylabel('');
    xlabel(sprintf(' Session %d', ses), FontSize=8); 
    set(h, 'AlphaData', ~isnan(maskedMatrix(parc.yeoOrder, parc.yeoOrder)));
end
saveas(gcf, fullfile(Path2Fig, 'ICC_within_figure.svg'))
saveas(gcf, fullfile(Path2Fig, 'ICC_within_figure.pdf'))
close all;

%% Prediction: Delta 1/RT ~ Sum|FC edges| 
Sub_vec = repelem(1:10, 1);
mdl = struct();

cmap = [
    0.2670, 0.0049, 0.3294;
    0.2823, 0.1408, 0.4571;
    0.2539, 0.2650, 0.5292;
    0.2068, 0.3718, 0.5533;
    0.1636, 0.4714, 0.5581;
    0.1276, 0.5669, 0.5504;
    0.1340, 0.6586, 0.5176;
    0.2669, 0.7487, 0.4402;
    0.4780, 0.8211, 0.3181;
    0.7414, 0.8735, 0.1505];

for ses = 1:n_ses
    Cog = matrix_Simon(ses, :);
    FC = FC_sumabs_within(:, ses);
    data = table(Cog', FC, Sub_vec', 'VariableNames', {'Behav', 'FC', 'Subject'});
    mdl = fitlm(data, 'Behav ~ FC');

    % Check if FC predictor is significant
    R_squared = mdl.Rsquared.Adjusted;
    p_value = mdl.Coefficients.pValue(2);
    estimate = mdl.Coefficients.Estimate(2);

    % save R and p values
    all_R(ses) = R_squared;
    all_p(ses) = p_value;
    all_estimate(ses) = estimate;

    x = data.FC;
    y = data.Behav;

    % Plot data and fitted line
    figure('Color', 'w'); set(gcf,'Position',[200 200 900 600]);
    hold on; grid on;
    scatter(x, y, 36, cmap(data.Subject, :), 'filled', 'MarkerFaceAlpha', 0.6);

    % Fixed effect line in dark gray
    intercept = mdl.Coefficients.Estimate(strcmp(mdl.CoefficientNames, '(Intercept)'));
    slope = mdl.Coefficients.Estimate(strcmp(mdl.CoefficientNames, 'FC'));
    x_fit = linspace(min(x), max(x), 100);
    y_fit = intercept + slope * x_fit;
    plot(x_fit, y_fit, 'Color', [0.3 0.3 0.3], 'LineWidth', 2);

    xlabel('sum of absolute FC values', 'FontSize', 12);
    ylabel('Δ1/RT', 'FontSize', 12);
    title(sprintf('Within-session prediction of the Simon task performance, session %d', ses), 'FontSize', 14);

    % Annotate with fit stats
    dim = [0.85 0.7 0.2 0.1];
    str = ['Adj. R^2 = ', num2str(R_squared, '%.3f'), newline, ...
       'p-value = ', num2str(p_value, '%.3f'), newline, ...
       '\beta-value = ', num2str(estimate, '%.3f')];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'FontSize', 10, 'Interpreter', 'tex');

    % Save figure
    saveas(gcf, fullfile(Path2Fig, (sprintf('fitlm_pred_case1_ses-0%d.png', ses))));
end
close all;


%% Heatmap with fit indices
all_infos = [round(all_R, 3); round(all_p, 3); round(all_estimate, 3)];
figure('Color', 'w'); set(gcf,'Position',[200 200 900 400]);
h = heatmap(all_infos);
h.YDisplayLabels = {'R^2 Adjusted', 'p-value', '\beta-value'};
h.XDisplayLabels = arrayfun(@(x) sprintf('S%d', x), 1:size(all_infos, 2), 'UniformOutput', false);
h.GridVisible = 'off';
h.FontSize = 10;
title('Fit indices of within-session Simon task prediction');
xlabel('Session');
saveas(gcf, fullfile(Path2Fig, 'heatmap_prediction_case1.png'));
close all;
save(fullfile(Path2newfig, "SD_case1.mat"), "all_infos");

% heatmap sum of absolute FC values
figure('Color', 'w'); set(gcf,'Position',[200 200 900 600]);
h = heatmap(FC_sumabs_within);
h.GridVisible = 'off';
h.FontSize = 12;
viridis_like = [
     68,   1,  84;
     71,  19, 117;
     59,  81, 139;
     33, 145, 140;
     94, 201,  98;
    253, 231,  37
] / 255;
colormap(interp1(linspace(0, 1, size(viridis_like, 1)), viridis_like, linspace(0, 1, 256)));
title('sum of absolute FC values');
xlabel('Session');
ylabel('Subject');
saveas(gcf, fullfile(Path2Fig, 'heatmap_sumabsFC_within.png'));
close all;


%% CASE 2: prediction with dense-sampling ICC
% Load data and initialize variables
load(fullfile(Path2results, 'ICC_2D_ds.mat'), 'ICC_2D_all');
ICC_2D = ICC_2D_all;
load(fullfile(Path2results, 'FC_3D_between.mat'), 'FC_3D_mat_all'); % edge x ses x sub
FC_3D = FC_3D_mat_all;
temp_ICC_mat_highFG = NaN(N, N);
ICC_mask_highFG_2D = NaN(numEdges, n_ses);
ICC_mat_highFG = zeros(N, N);
FC_sumabs_ds = nan(n_sub, n_ses, n_ses);  %subject x mask x model

% Build predictor: sum|FC values| per session and model
for ses = 2:n_ses
    temp_ICC = ICC_2D(:,ses);
    mask_highFG_2D = temp_ICC >= prctile(temp_ICC, 99.9);

    for s = 1:n_sub
        % Extract FC data for all sessions up to 'ses'
        temp_FC_subject = FC_3D(:, 1:ses, s);
        masked_FC = temp_FC_subject .* repmat(mask_highFG_2D, 1, ses);
        % Sum absolute values per model
        temp_FC_sum = sum(abs(masked_FC), 1, 'omitnan');
        FC_sumabs_ds(s, 1:ses, ses) = temp_FC_sum;
    end

    % plot ICC top edges
    top_title = sprintf("Dense-sampling ICC edges > 99.9th percentile, session 1 to %d", ses);
    ICC_vec_highFG = temp_ICC.* mask_highFG_2D;
    ICC_mat_highFG(mask_ut) = ICC_vec_highFG';
    ICC_mat_plot = ICC_mat_highFG + ICC_mat_highFG';
    figure; 
    imagesc(ICC_mat_plot(parc.yeoOrder, parc.yeoOrder)); 
    axis square; 
    colorbar;
    title(top_title); 

    % Check if NaN and save figure
    nnz(isnan(ICC_mat_highFG));
    save_filename = sprintf('ICC_ds_highFG_ses1-0%d.png', ses);
    saveas(gcf, fullfile(Path2Fig, save_filename));
    close;
    ICC_mask_highFG_2D(:, ses) = ICC_vec_highFG;
end

save(fullfile(Path2results, "ICC_mask_highFG_2D_ds.mat"), "ICC_mask_highFG_2D");
save(fullfile(Path2results, "FC_sumabs_ds.mat"), "FC_sumabs_ds");
writematrix(FC_sumabs_ds, 'FC_sumabs_ds.csv');


%% Plot dense-sampling ICC edges > 95th percentile for visualisation
ICC_mask_highFG_2D = NaN(numEdges, n_ses);
temp_ICC_mat_highFG = NaN(N, N);
ICC_mat_highFG = zeros(N, N);
ICC_2D_correct = ICC_2D(:, 2:10);
percentile_threshold = prctile(ICC_2D, 95);

figure('Color', 'w'); set(gcf,'Position',[400 400 1100 800]);
title('Cumulative ICC edges > 95th percentile');
for ses = 1:9
    temp_ICC_2D = ICC_2D_correct(:,ses);
    mask_highFG_2D = temp_ICC_2D >= percentile_threshold(:, ses+1);

    % Display top delta ICC edges
    ICC_vec_highFG = temp_ICC_2D.* mask_highFG_2D;
    ICC_mat_highFG(mask_ut) = ICC_vec_highFG';
    ICC_mat_between = ICC_mat_highFG + ICC_mat_highFG';
 
    subplot(3, 3, ses);
    imagesc(ICC_mat_between(parc.yeoOrder, parc.yeoOrder));
    axis square;
    set(gca, 'XTick', [], 'YTick', []);
    ylabel('');
    xlabel(sprintf(' Session 1 to %d', ses+1), FontSize=8); 
end
print(gcf, '-bestfit')
saveas(gcf, fullfile(Path2Fig, 'ICC_cumulative_figure95.svg'))
close all;


%% Prediction: Delta 1/RT ~ Sum|FC edges| + (1|Subject) + (1|Session) 
load(fullfile(Path2results, 'FC_sumabs_ds.mat'));
Sub = repelem(1:10, 10)';
Ses = repmat(1:10, 1, 10)';
all_R = NaN(1, n_ses);
all_p = NaN(1, n_ses);
all_estimate = NaN(1, n_ses);
all_infos = NaN(3, 10);
for ses = 2:n_ses
    Cog = matrix_Simon(1:ses, :);
    Cog_vec = Cog(:);
    FC_vec = FC_sumabs_ds(:, 1:ses, ses);
    FC_vec = FC_vec(:);
   
    idx = Ses <= ses;
    Sub_vec = Sub(idx);
    Ses_vec = Ses(idx);

subjects = unique(Sub_vec);
cmap = [
    0.2670, 0.0049, 0.3294;
    0.2823, 0.1408, 0.4571;
    0.2539, 0.2650, 0.5292;
    0.2068, 0.3718, 0.5533;
    0.1636, 0.4714, 0.5581;
    0.1276, 0.5669, 0.5504;
    0.1340, 0.6586, 0.5176;
    0.2669, 0.7487, 0.4402;
    0.4780, 0.8211, 0.3181;
    0.7414, 0.8735, 0.1505];

    % Build table
    data = table(Cog_vec, FC_vec, Sub_vec, Ses_vec, ...
        'VariableNames', {'Cog', 'FC', 'Subject', 'Session'});

    % Fit linear mixed-effects model with random intercepts for Subject and Session
    mdl = fitlme(data, 'Cog ~ FC + (1|Subject) + (1|Session)');

    % Extract fixed effect stats
    R_squared = mdl.Rsquared.Adjusted;
    p_value = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Name, 'FC'));
    estimate = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name, 'FC'));

    % Save all R, p and beta values
    all_R(1, ses) = R_squared;
    all_p(1, ses) = p_value;
    all_estimate(1, ses) = estimate;

    % Extract data for plotting
    x = data.FC;
    y = data.Cog;

    % Set the figure
    figure('Color', 'w'); set(gcf,'Position',[200 200 900 600]);
    hold on; grid on;

     % Map subjects to indices and colors
    [~, sub_idx] = ismember(data.Subject, subjects);
    unique_subjects = unique(data.Subject);

    % Scatter plot with all points colored by subject
    scatter(x, y, 36, cmap(sub_idx, :), 'filled', 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');

    % Plot fixed effect regression line
    intercept = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name, '(Intercept)'));
    slope = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name, 'FC'));
    x_fit = linspace(min(x), max(x), 100);
    y_fit = intercept + slope * x_fit;
    plot(x_fit, y_fit, 'Color', [0.3 0.3 0.3], 'LineWidth', 2, 'HandleVisibility', 'off');

    % Add dummy scatter objects per subject for legend
    for i = 1:numel(unique_subjects)
       scatter(nan, nan, 36, cmap(i, :), 'filled', ...
           'DisplayName', sprintf('Subject %s', string(unique_subjects(i))));
    end

    % Plot subject-specific fit slopes
    for i = 1:numel(unique_subjects)
        subj = unique_subjects(i);
        subj_idx = data.Subject == subj;
        x_sub = data.FC(subj_idx);
        y_sub = data.Cog(subj_idx);

        %Only fit line if there's variation in x
         if numel(unique(x_sub)) > 1
             p = polyfit(x_sub, y_sub, 1);
             x_fit = linspace(min(x_sub), max(x_sub), 50);
             y_fit = polyval(p, x_fit);
             plot(x_fit, y_fit, '-', 'Color', cmap(i, :), 'LineWidth', 1.5, ...
                 'HandleVisibility', 'off');
         end
    end

    % Axis labels and title
    xlabel('Sum of absolute FC values', 'FontSize', 12);
    ylabel('Δ1/RT', 'FontSize', 12);
    title(sprintf('Prediction of the Simon task based on dense-sampling ICC, session 1 to %d', ses), ...
        'FontSize', 14);

    % Annotate figure with adjusted R², p-value and beta-value
    dim = [0.78 0.7 0.2 0.1];
    str = ['Adj. R^2 = ', num2str(R_squared, '%.4f'), newline, ...
       'p-value = ', num2str(p_value, '%.4f'), newline, ...
       '\beta-value = ', num2str(estimate, '%.4f')];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', 'FontSize', 10, 'Interpreter', 'tex');

    % Show legend with subject-color mapping
    legend('Location', 'eastoutside');

    % Save figure
    saveas(gcf, fullfile(Path2Fig, sprintf('fitlme_pred_case2_ses1to%d.png', ses)));
end
close all


%% heatmap with fit indices
all_infos = [round(all_R, 4);round(all_p, 4);round(all_estimate, 4)];
figure('Color', 'w');
set(gcf,'Position',[200 200 900 400]);
h = heatmap(all_infos);
h.YDisplayLabels = {'R^2 Adjusted', 'p-value', '\beta-value'};
h.GridVisible = 'off';
h.FontSize = 10;
title('Fit indices of Simon task prediction with dense-sampling ICC');
xlabel('Session');
saveas(gcf, fullfile(Path2Fig,'heatmap_fitlme_case2.png'));
close all;
save(fullfile(Path2newfig, "SD_case2.mat"), "all_infos");

% heatmap sum of absolute FC values
for model = 2:10
    figure('Color', 'w');
    set(gcf,'Position',[200 200 900 600]);
    h = heatmap(FC_sumabs_ds(:, 1:model, model));
    h.GridVisible = 'off';
    h.FontSize = 12;
    viridis_like = [
     68,   1,  84;
     71,  19, 117;
     59,  81, 139;
     33, 145, 140;
     94, 201,  98;
    253, 231,  37
    ] / 255;
    colormap(interp1(linspace(0, 1, size(viridis_like, 1)), viridis_like, linspace(0, 1, 256)));
    title(sprintf('Sum of absolute FC values, model %d', model-1));
    xlabel('Session');
    ylabel('Subject');
    temp = sprintf('heatmap_sumabsFC_case2_model_%d.png', model);
    saveas(gcf, fullfile(Path2Fig, temp));
end
close all;



%% CASE 3: within-session dense-sampling ICC prediction
load(fullfile(Path2results, 'FC_sumabs_selected_edges_within.mat'));
load(fullfile(Path2results, 'matrix_Simon.mat'));
FC_vec = FC_sumabs_within(:);
Sub_vec = repelem(1:10, 10)';
Ses_vec = repmat(1:10, 1, 10)';
Cog = matrix_Simon(:);
subjects = unique(Sub_vec);
N_sub = numel(subjects);
all_R = NaN(1, 10);
all_p = NaN(1, 10);
all_estimate = NaN(1,10);
cmap = [
    0.2670, 0.0049, 0.3294;
    0.2823, 0.1408, 0.4571;
    0.2539, 0.2650, 0.5292;
    0.2068, 0.3718, 0.5533;
    0.1636, 0.4714, 0.5581;
    0.1276, 0.5669, 0.5504;
    0.1340, 0.6586, 0.5176;
    0.2669, 0.7487, 0.4402;
    0.4780, 0.8211, 0.3181;
    0.7414, 0.8735, 0.1505];
n_ses = 10;
for ses = 2:n_ses
    % Select all data where session is ≤ current ses
    idx = Ses_vec <= ses;

    % Build table with selected data
    data = table(Cog(idx), FC_vec(idx), Sub_vec(idx), Ses_vec(idx), ...
        'VariableNames', {'Cog', 'FC', 'Subject', 'Session'});

    % Fit linear mixed-effects model with random intercepts for Subject and Session
    mdl = fitlme(data, 'Cog ~ FC + (1|Subject) + (1|Session)');

    % Extract fixed effect stats
    R_squared = mdl.Rsquared.Adjusted;
    p_value = mdl.Coefficients.pValue(strcmp(mdl.Coefficients.Name, 'FC'));
    estimate = mdl.Coefficients.Estimate(strcmp(mdl.Coefficients.Name, 'FC'));

    % Save all R, p and estimate values
    all_R(1, ses) = R_squared;
    all_p(1, ses) = p_value;
    all_estimate(1, ses) = estimate;

    % Extract data for plotting
    x = data.FC;
    y = data.Cog;
    
    % Plot data and fitted lines
    figure; set(gcf,'Position',[200 200 900 600]);
    hold on; grid on;
    
    % Map subjects to indices and colors
    [~, sub_idx] = ismember(data.Subject, subjects);
    unique_subjects = unique(data.Subject);
    
    % Main scatter plot with all points colored by subject
    scatter(x, y, 36, cmap(sub_idx, :), 'filled', 'MarkerFaceAlpha', 0.8, 'HandleVisibility', 'off');
    
    % Plot subject-specific fit lines
     for i = 1:numel(unique_subjects)
         subj = unique_subjects(i);
         subj_idx = data.Subject == subj;
         x_sub = data.FC(subj_idx);
         y_sub = data.Cog(subj_idx);
     
         % Only fit line if there's variation in x
        if numel(unique(x_sub)) > 1
             p = polyfit(x_sub, y_sub, 1);
             x_fit = linspace(min(x_sub), max(x_sub), 50);
             y_fit = polyval(p, x_fit);
             plot(x_fit, y_fit, '-', 'Color', cmap(i, :), 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('Subject %d', subj));
        end
    end
    
    xlabel('sum of absolute FC values', 'FontSize', 12);
    ylabel('Δ1/RT', 'FontSize', 12);
    title(sprintf('Prediction of the Simon task based on within-session dense-sampling ICC, session 1 to %d', ses), 'FontSize', 14);
    
    % Annotate with global model stats
    dim = [0.15 0.75 0.2 0.1];
    str = ['Adj. R^2 = ', num2str(R_squared, '%.4f'), newline, ...
         'p-value = ', num2str(p_value, '%.4f'), newline, ...
         '\beta-value = ', num2str(estimate, '%.4f')];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
         'BackgroundColor', 'white', 'FontSize', 10, 'Interpreter', 'tex');
    
    % Show legend with subject-color mapping
    legend('Location', 'eastoutside');
    
    % Save figure
    saveas(gcf, fullfile(Path2Fig, (sprintf('fitlme_pred_within_ses1to%d.png', ses))));
end
close all


%% Heatmap with fit indices
all_infos = NaN(3, 10);
all_infos = [round(all_R, 4); round(all_p, 4); round(all_estimate, 4)];
figure('Color', 'w'); set(gcf,'Position',[200 200 900 400]);
h = heatmap(all_infos);
h.YDisplayLabels = {'R^2 Adjusted', 'p-value', '\beta-value'};
h.GridVisible = 'off';
h.FontSize = 10;
title('Fit indices of Simon task prediction with within-session dense-sampling ICC');
xlabel('Session');
saveas(gcf, fullfile(Path2Fig, 'heatmap_fitlme_case3.png'));
close all;
save(fullfile(Path2newfig, "SD_case3.mat"), "all_infos");

% heatmap sum of absolute FC values
figure('Color', 'w'); set(gcf,'Position',[200 200 900 600]);
h = heatmap(FC_sumabs_within);
h.GridVisible = 'off';
h.FontSize = 16;
viridis_like = [
     68,   1,  84;
     71,  19, 117;
     59,  81, 139;
     33, 145, 140;
     94, 201,  98;
    253, 231,  37
] / 255;
colormap(interp1(linspace(0, 1, size(viridis_like, 1)), viridis_like, linspace(0, 1, 256)));
title('sum of absolute FC values');
xlabel('Session');
ylabel('Subject');
saveas(gcf, fullfile(Path2Fig, 'heatmap_sumabsFC_case3.png'));
close all;