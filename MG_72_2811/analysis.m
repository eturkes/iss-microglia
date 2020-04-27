%% Load Data
load(fullfile(pwd, 'o.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'Pixel');
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

iss_change_plot_MG(o,'Pixel', GeneName)

%% Individual Genes
clear variables

load(fullfile(pwd, 'o.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

o.plot(o.BigDapiFile, Roi, 'Pixel');
set(gca, 'YDir', 'reverse')
set(gca, 'XDir', 'reverse')

% For all genes.
%%%%%%%%%%%%%%%%
fp = fopen(fullfile('codebook_comb.txt'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName = tmp{1};
fclose(fp);

for i = 1:length(GeneName)
    iss_change_plot_individual(o,'Pixel', GeneName(i))

    saveas(gcf, fullfile('figures', 'latest-software-pixelbased', ...
        'allgenes-individual', GeneName{i}), 'fig')
    saveas(gcf, fullfile('figures', 'latest-software-pixelbased', ...
        'allgenes-individual', GeneName{i}), 'svg')
end
%%%%%%%%%%%%%%%%