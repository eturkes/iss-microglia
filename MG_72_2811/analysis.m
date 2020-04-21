%% Load Data
load(fullfile(pwd, 'latest-software-pixelbased.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'Prob');
set(gca, 'YDir', 'reverse')
set(gca, 'XDir', 'reverse')

%% MG genes
% Read in codebook.
%%%%%%%%%%%%%%%%%%%
fp = fopen(fullfile('codebook_Seppe.txt'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName = tmp{1};
fclose(fp);
%%%%%%%%%%%%%%%%%%%

% Remove unwanted genes.
%%%%%%%%%%%%%%%%%%%%%%%%
GeneName(find(strcmp(GeneName, 'Ptk2b'))) = [];
%%%%%%%%%%%%%%%%%%%%%%%%

iss_change_plot_MG(o,'Prob', GeneName)