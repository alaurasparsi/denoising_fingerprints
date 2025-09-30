%!/usr/bin/env matlab
% Author: Laura Belli, Sara Stampacchia
%% Functional Connectivity
clc
clear
close
%% 1.0 set paths and constants
Path2results = 'mnt\derivatives';
Path2Parc = 'mnt\parcs';
Path2data = 'mnt\data';
Path2plots = 'mnt\derivatives\plots';
load(fullfile(Path2Parc, '\shen_yeo_RS7.mat'));
n_sub=10;
n_ses=10;
ROIs=278;
N = size(parc.yeoROIs, 1);
mask_ut = triu(true(N), 1);
numEdges = 38503;

% save 2 and 98 percentile values for each FC matrix, in order to build FC plots with the same color scale limits
FC_4D_mat_all = NaN(ROIs, ROIs, n_ses, n_sub);
all_margins = [];
for s = 1:n_sub
    for ses = 1:n_ses
        % build dynamic filename
        data = sprintf('sub-%03d_ses-%02d_task-rest_run-01_Shen_dil_timeseries_Stampacchia_denoised.1D', s, ses);
        temp_TS_ROI = load(fullfile(Path2data, data));

        % compute correlation matrix
        FC_3D = corr(temp_TS_ROI);

        % save FC matrix
        save_FC_matrix = sprintf('FC_matrix_sub_%03d_ses_%02d_task-rest_run-01.mat', s, ses);
        save(fullfile(Path2results, save_FC_matrix), 'FC_3D');
        FC_3D_mat_all(:, ses, s) = FC_3D(mask_ut);

        % compute and store margins
        temp_margins = prctile(FC_3D(:), [2 98]);
        all_margins = [all_margins; temp_margins];  

        % plot connectomes
        figure;
        imagesc(FC_3D(parc.yeoOrder, parc.yeoOrder));
        caxis([-0.6 0.8]);
        axis square;
        colorbar;
        title(sprintf('FC - sub %03d, session %02d, Stampacchia denoised', s, ses), 'Interpreter', 'none');

        % save figure
        save_FC_plot = sprintf('FC_plot_sub_%03d_ses_%02d_task-rest_run-01.png', s, ses);
        saveas(gcf, fullfile(Path2plots, save_FC_plot));
        close;
    end
end
% MISSING ROIs: 158 AND 239 SUB-008, 158 IN SUB-010
save(fullfile(Path2results, 'FC_3D_between'), 'FC_3D_mat_all');



%% Functional Connectivity splitted
% Preallocate storage per session, so that you have all subject together
FC_A_tot = nan(numEdges, n_sub);
FC_B_tot = nan(numEdges, n_sub);

for ses = 1:n_ses
   for s = 1:n_sub

        data = sprintf('sub-%03d_ses-%02d_task-rest_run-01_Shen_dil_timeseries_Sara_denoised.1D', s, ses);
        TS = load(fullfile(Path2data, data));   
        split = length(TS)/2;

        % split dataset
        A = TS(1:split,:);
        B = TS(split+1:390,:);

        % compute and save correlation matrix
        FC_A = corr(A);
        FC_A_mask = FC_A(mask_ut);
        FC_B = corr(B);
        FC_B_mask = FC_B(mask_ut);
        FC_matrix_split = sprintf('FC_matrix_sub_%03d_ses_%02d_split_task-rest_run-01.mat', s, ses);
        save(fullfile(Path2results, FC_matrix_split), 'FC_A', 'FC_B');

        % store FC vectors
        FC_A_tot(:, s) = FC_A_mask(:);
        FC_B_tot(:, s) = FC_B_mask(:);

        % plot splitted connectomes
        figure;
        imagesc(FC_A(parc.yeoOrder, parc.yeoOrder));
        caxis([-0.6 0.8]);
        axis square;
        colorbar;
        title(sprintf('FC - sub %03d, ses %02d split (A), Stampacchia denoised', s, ses), 'Interpreter', 'none');
        save_FC_plot = sprintf('FC_plot_sub_%03d_ses_%02d_splitA_task-rest_run-01.png', s, ses);
        saveas(gcf, fullfile(Path2plots, save_FC_plot));
        close;

        figure;
        imagesc(FC_B(parc.yeoOrder, parc.yeoOrder));
        caxis([-0.6 0.8]);
        axis square;
        colorbar;
        title(sprintf('FC - sub %03d, ses %02d split (B), Stampacchia denoised', s, ses), 'Interpreter', 'none');
        save_FC_plot = sprintf('FC_plot_sub_%03d_ses_%02d_splitB_task-rest_run-01.png', s, ses);
        saveas(gcf, fullfile(Path2plots, save_FC_plot));
        close;
    end
    FC_A_all(:,:,ses) = FC_A_tot(:,:);
    FC_B_all(:,:,ses) = FC_B_tot(:,:);
end

save(fullfile(Path2results, 'FC_split4within'), 'FC_A_all', 'FC_B_all'); % edges, subject, session