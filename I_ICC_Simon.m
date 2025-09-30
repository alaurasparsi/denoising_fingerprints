%% ICC on Simon's task
% last modification: 13.07.25, Laura Belli

%% 1.0 set paths and variables
addpath('functions\')
N_ses = 10;
N_sub = 10;
ICC_vec = nan(1, N_ses);
mask_ut = triu(true(N_ses), 1);     % mask upper triangle
ICC_matrix = nan(N_ses, N_sub);     % row = session, col = subject
Path2Results = fullfile('..', 'results\');
matrix_Simon = load('..\..\ICC_within\1.prediction\matrix_Simon.mat', 'matrix_Simon');
matrix_Simon = matrix_Simon.matrix_Simon;

% compute ICC(1,1)
for ses = 2:N_ses  % Iterate over sessions
    data_ICC = []; 
    temp_data_ICC = matrix_Simon(1:ses, :);
    ICC_vec(1, ses) = ICC(temp_data_ICC', '1-1');
end
save(fullfile(Path2Results, 'ICC_Simon.mat'), 'ICC_vec');


%% 2.0 plot ICC results
% histogram
Simon_hist = histogram(ICC_vec, "BinLimits", [0, 1], "NumBins", 10);
title("Simon's task reliability: ICC(1,1) on 2 to 10 sessions")

% heatmap
h = heatmap(ICC_vec);
h.MissingDataColor = [1 1 1]; % white for NaNs
h.GridVisible = 'off';
title('Simon task reliability: ICC(1,1) on 2 to 10 sessions')

