%% Load Data
load(fullfile(pwd, 'o.mat'));
BigDapiImage = imread(o.BigDapiFile);

%% Plot All Spots
o.plot(BigDapiImage);
iss_change_plot(o, 'Prob');
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