%!/usr/bin/env matlab
% Author: Laura Belli

%% thesis figures
Path2newfig = 'mnt\derivatives\prediction';
MD_case1 = load(fullfile(Path2newfig, "MD_case1.mat"));
MD_case2 = load(fullfile(Path2newfig, "MD_case2.mat"));
MD_case3 = load(fullfile(Path2newfig, "MD_case3.mat"));
MD_case1 = MD_case1.all_infos;
MD_case2 = MD_case2.all_infos;
MD_case3 = MD_case3.all_infos;

conditions = {'MD', 'PD', 'SD', 'MDbp', 'PDbp'};
cases = {'case1', 'case2', 'case3'};

% Initialize a struct to store all loaded data
all_data = struct();

for i = 1:length(conditions)
    condition = conditions{i};
    for j = 1:length(cases)
        case_name = cases{j};
        filename = sprintf('%s_%s.mat', condition, case_name);
        filepath = fullfile(Path2newfig, filename);
        data = load(filepath);
        % Store 'all_infos' in struct: all_data.MD.case1, etc.
        all_data.(condition).(case_name) = data.all_infos;
    end
end
%% data manipulation for figures
first_rows = zeros(length(cases), length(conditions), 10);  % 3 x 5 x 10

for i = 1:length(conditions)
    cond = conditions{i};
    for j = 1:length(cases)
        cs = cases{j};
        first_rows(j, i, :) = all_data.(cond).(cs)(1, :);  % extract adjusted R squared
    end
end

%% case 1
data_heatmap_c1 = squeeze(first_rows(1, :, :));
violet_map = [linspace(0.95, 0.9, 64)', linspace(0.9, 0.3, 64)', linspace(1.0, 0.5, 64)'];
figure('Color', 'w'); set(gcf,'Position',[100 100 1100 800]);
h = heatmap(round(data_heatmap_c1', 4));
h.FontSize = 14;
title('Case 1: adusted R2 of within-session ICC predictive models');
h.Colormap = violet_map;
xlabel('denoising model');
ylabel('session');
saveas(gcf, fullfile(Path2newfig, 'R2_case1.svg'));

%% case 2
data_heatmap_c2 = squeeze(first_rows(2, :, 2:10));
violet_map = [linspace(0.95, 0.9, 64)', linspace(0.9, 0.3, 64)', linspace(1.0, 0.5, 64)'];
figure('Color', 'w'); set(gcf,'Position',[100 100 1100 800]);
h = heatmap(round(data_heatmap_c2', 4));
h.FontSize = 14;
title('Case 2: adjusted R² of dense-sampling ICC predictive models');
h.Colormap = violet_map;
xlabel('denoising model');
ylabel('session');
saveas(gcf, fullfile(Path2newfig, 'R2_case2.svg'));

%% case 3
data_heatmap_c3 = squeeze(first_rows(3, :, 2:10));
violet_map = [linspace(0.95, 0.9, 64)', linspace(0.9, 0.3, 64)', linspace(1.0, 0.5, 64)'];
figure('Color', 'w'); set(gcf,'Position',[100 100 1100 800]);
h = heatmap(round(data_heatmap_c3', 4));
h.FontSize = 14;
title('Case 3: adjusted R² of within-session dense-sampling ICC predictive models');
h.Colormap = violet_map;
xlabel('denoising model');
ylabel('session');
saveas(gcf, fullfile(Path2newfig, 'R2_case3.svg'));

%% figure presentation with all SD model results
SD_data = all_data.SD;
violet_map = [linspace(0.95, 0.9, 64)', linspace(0.9, 0.3, 64)', linspace(1.0, 0.5, 64)'];
figure('Color', 'w'); set(gcf,'Position',[100 100 1100 800]);
h = heatmap(SD_data.case2(:,2:10));
h.YDisplayLabels = {'adjusted R^2', 'p-value', '\beta-value'};
h.GridVisible = 'off';
h.FontSize = 14;
title('Case 2: dense-sampling ICC predictive models');
h.Colormap = violet_map;
xlabel('session');