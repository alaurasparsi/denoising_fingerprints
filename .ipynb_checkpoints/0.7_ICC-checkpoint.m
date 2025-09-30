%!/usr/bin/env matlab
% Author: Laura Belli
%% edgewise Intra-class Correlation Coefficient (1,1)
clc
clear
close

%% 1.0 set paths and constants
Path2data = 'mnt\data';
Path2Parc = 'mnt\parcs';
Path2results = 'mnt\derivatives';
Path2plots = 'mnt\derivatives\plots';
load(fullfile(Path2Parc, '\shen_yeo_RS7.mat'));
Path2functions = 'mnt\functions';
addpath('mnt\functions\');
load(fullfile(Path2results,"FC_3D_between.mat"));
numEdges = 38503;
N = 278;
n = 10;
n_ses = 10;
n_sub = 10;
mask_ut = triu(true(N), 1);
ICC_2D_all = NaN(numEdges, n_ses);
table_MSR = NaN(numEdges, n_ses);
table_MSW = NaN(numEdges, n_ses);
table_MSC = NaN(numEdges, n_ses);

%% dense-sampling edgewise ICC(1,1)
for ses = 2:n_ses
    ICC_matrix = zeros(N, N);
    ICC_vec = nan(1, numEdges);

    for comp = 1:numEdges
        data_ICC = [];  % Initialize for each edge
        temp_data_ICC = FC_3D_mat_all(comp, 1:ses, :);
        tem_data_ICC = squeeze(temp_data_ICC);

        % Remove rows with NaNs
        rows2delete = isnan(sum(tem_data_ICC, 1));
        tem_data_ICC(:,rows2delete) = [];

        % Compute ICC(1,1) across sessions
        data_ICC = tem_data_ICC';
        ICC_vec(1, comp) = ICC(data_ICC, '1-1');
        [p, table] = anova_rm(data_ICC, 'off');
        SSR = table{3,2};
        SSE = table{4,2};
        SSC = table{2,2};
        SSW = SSE + SSC;
        MSR = SSR / (n-1);
        MSC = SSC / (ses-1);
        MSW = SSW / (n*(ses-1));
        table_MSR(comp, ses) = MSR;
        table_MSW(comp, ses) = MSW;
        table_MSC(comp, ses) = MSC;
    end

    % Build and save vector and matrix
    ICC_matrix(mask_ut) = ICC_vec;
    ICC_matrix = ICC_matrix + ICC_matrix';
    save_filename = sprintf('ICC_ses1-%d_rest_run-01.mat', ses);
    save(fullfile(Path2results, save_filename), 'ICC_matrix', 'ICC_vec');
    ICC_2D_all(:, ses) = ICC_vec;
    
    figure;
    imagesc(ICC_matrix(parc.yeoOrder, parc.yeoOrder));
    caxis([-0.6 1]);
    colorbar;
    axis square;
    title(sprintf('ICC on sessions 1 to %d, Stampacchia denoised', ses), 'Interpreter', 'none');

    %save
    save_ICC_plot = sprintf('ICC_ses1-%d_task-rest_run-01.png', ses);
    saveas(gcf, fullfile(Path2plots, save_ICC_plot));
    close;
end

mean(mean(table_MSR(:, 2:10)));
mean(std(table_MSR(:,2:10)));
mean(mean(table_MSW(:, 2:10)));
mean(std(table_MSW(:,2:10)));

save(fullfile(Path2results, 'ICC_2D_ds.mat'), 'ICC_2D_all');
save(fullfile(Path2results, 'table_MSR_ds.mat'), 'table_MSR');
save(fullfile(Path2results, 'table_MSW_ds.mat'), 'table_MSW');


%% Compute within-session edgewise ICC(1,1)
clc
clear
close all
n = 10;
N = 278;
n_ses = 10;
numEdges = 38503;
table_MSR = NaN(numEdges, n_ses);
table_MSW = NaN(numEdges, n_ses);
table_MSC = NaN(numEdges, n_ses)
ICC_within_2D_all = NaN(numEdges, n_ses);
mask_ut = triu(true(N),1);
load(fullfile(Path2results, 'FC_split4within.mat'));

for ses = 1:n_ses
    ICC_matrix = zeros(N, N);
    ICC_vec = nan(1, numEdges);
    
    for comp=1:numEdges
        current_datA(comp,:) = FC_A_all(comp,:,ses);
        current_datB(comp,:) = FC_B_all(comp,:,ses);

        data_ICC = [];
        temp_data_ICC = [current_datA(comp,:)' current_datB(comp,:)'];
        rows2delete_A = isnan(sum(temp_data_ICC,2));
        temp_data_ICC(rows2delete_A,:) = [];
        ICC_vec(1,comp) = ICC(temp_data_ICC,'1-1') ;
        [p, table] = anova_rm(temp_data_ICC, 'off');
        SSR = table{3,2};
        SSE = table{4,2};
        SSC = table{2,2};
        SSW = SSE + SSC;
        MSR = SSR / (n-1);
        MSC = SSC / (2-1);
        MSW = SSW / (n*(2-1));
        table_MSR(comp, ses) = MSR;
        table_MSW(comp, ses) = MSW;
        table_MSC(comp, ses) = MSC;
    end

    % build symmetric matrix
    ICC_matrix(mask_ut) = ICC_vec;
    ICC_matrix = ICC_matrix + ICC_matrix';

    % Save both matrix and vector for each session
    save_filename = sprintf('ICC_ses-%02d_rest_run-01.mat', ses);
    save(fullfile(Path2results, save_filename), 'ICC_matrix', 'ICC_vec');
    ICC_within_2D_all(:, ses) = ICC_vec;

    % plot and save figures
    figure;
    imagesc(ICC_matrix(parc.yeoOrder, parc.yeoOrder));
    caxis([-0.6 1]);
    axis square;
    colorbar;
    title(sprintf('Within-session %02d ICC(1,1) Stampacchia denoised', ses), 'Interpreter', 'none');
    save_ICC_plot = sprintf('ICC_ses-%02d_task-rest_run-01.png', ses);
    saveas(gcf, fullfile(Path2plots, save_ICC_plot));
    close;
end

mean(mean(table_MSR));
mean(std(table_MSR));
mean(mean(table_MSW));
mean(std(table_MSW));

save(fullfile(Path2results, 'ICC_2D_within.mat'), 'ICC_within_2D_all');
save(fullfile(Path2results, 'table_MSR.mat'), 'table_MSR');
save(fullfile(Path2results, 'table_MSW.mat'), 'table_MSW');