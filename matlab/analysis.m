%    This file is part of iss-microglia.
%    Copyright (C) 2020  Emir Turkes
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%    Emir Turkes can be contacted at emir.turkes@eturkes.com

%% Prep Data
load(fullfile(pwd, 'o.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

mkdir(fullfile('figures', 'latest-software-pixelbased'))

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'Pixel');
caxis([0,25000]);
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
spots_to_cluster = find(ismember(o.GeneNames, GeneNamesMGFilt));
spots_to_cluster = o.quality_threshold('Pixel') & ismember(o.pxSpotCodeNo, spots_to_cluster);
SpotSetClustered = get_gene_clusters(o, 'Pixel', 150, 25, spots_to_cluster);

iss_change_plot_MG(o, 'Pixel', GeneNamesMGFilt, SpotSetClustered);

% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG_clusters'), 'fig');
% saveas(gcf, fullfile('figures', 'latest-software-pixelbased', 'MG_clusters'), 'svg');

%% Individual Genes
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
%         GeneNamesAll{i}), 'svg');
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
%         GeneNamesMG{i}), 'svg');
end
%%%%%%%%%%%%%%%