%% Prep Data
load(fullfile(pwd, 'o.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

mkdir(fullfile('figures', 'latest-software-pixelbased'))

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'Pixel');
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');

% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'allgenes'), 'fig');
% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'allgenes'), 'svg');

%% All Spots Clusters
SpotSetClustered = get_gene_clusters(o, 'Pixel');
iss_change_plot(o, 'Pixel', o.GeneNames, SpotSetClustered);

% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'allgenes_clusters'), 'fig');
% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'allgenes_clusters'), 'svg');

%% MG genes
% Read in codebook.
%%%%%%%%%%%%%%%%%%%
fp = fopen(fullfile('codebook_Seppe.txt'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneNamesMG = tmp{1};
fclose(fp);
%%%%%%%%%%%%%%%%%%%

% Remove unwanted genes.
%%%%%%%%%%%%%%%%%%%%%%%%
GeneNamesMGFilt = GeneNamesMG;
GeneNamesMGFilt(find(strcmp(GeneNamesMGFilt, 'Ptk2b'))) = [];
%%%%%%%%%%%%%%%%%%%%%%%%

iss_change_plot_MG(o,'Pixel', GeneNamesMGFilt);

% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG'), 'fig');
% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG'), 'svg');

%% MG Clusters
% Must reload figure after calling iss_change_plot_MG().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
o.plot(o.BigDapiFile, Roi, 'Pixel');
set(gca, 'YDir', 'reverse')
set(gca, 'XDir', 'reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iss_change_plot_MG(o, 'Pixel', GeneNamesMGFilt, SpotSetClustered);

% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG_clusters'), 'fig');
% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG_clusters'), 'svg');

%% Individual Genes
% Must reload figure after calling iss_change_plot_MG().
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf;
o.plot(o.BigDapiFile, Roi, 'Pixel');
set(gca, 'YDir', 'reverse')
set(gca, 'XDir', 'reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For all genes.
%%%%%%%%%%%%%%%%
fp = fopen(fullfile('codebook_comb.txt'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneNamesAll = tmp{1};
fclose(fp);

mkdir(fullfile('figures', 'latest-software-pixelbased', 'allgenes-individual'));

for i = 1:length(GeneNamesAll)
    iss_change_plot_individual(o,'Pixel', GeneNamesAll(i));

%     saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'allgenes-individual', 'fig', ...
%         GeneName{i}), 'fig');
%     saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'allgenes-individual', 'svg', ...
%         GeneName{i}), 'svg');
end
%%%%%%%%%%%%%%%%

% For MG genes.
%%%%%%%%%%%%%%%
mkdir(fullfile('figures', 'latest-software-pixelbased', 'MG-individual'));

for i = 1:length(GeneNamesMG)
    iss_change_plot_individual(o,'Pixel', GeneNamesMG(i));

%     saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG-individual', 'fig', ...
%         GeneName{i}), 'fig');
%     saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG-individual', 'svg', ...
%         GeneName{i}), 'svg');
end
%%%%%%%%%%%%%%%