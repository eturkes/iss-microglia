%    This file is part of iss-microglia.
%    Copyright (C) 2020  Emir Turkes, Sebastiaan De Schepper, UK DRI at UCL
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
load(fullfile('results', 'data', 'pilot', 'o.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);
o.BigDapiFile = fullfile('results', 'data', 'pilot', 'background_image.tif');

o.pIntensityThresh = max(o.pSpotIntensity) + 1;
o.pScoreThresh = 7.5;

% Make codebook subsets.
%%%%%%%%%%%%%%%%%%%%%%%%
fp = fopen(fullfile('codebook_comb.txt'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneNamesAll = tmp{1};
fclose(fp);

fp = fopen(fullfile('assets', 'codebooks', 'pilot', 'codebook_Seppe.txt'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneNamesMG = tmp{1};
fclose(fp);

GeneNamesMGFilt = GeneNamesMG([4:6,8,9,13,14,20]); % MG specific genes.
GeneNamesMGFilt2 = setdiff(GeneNamesMG, GeneNamesMGFilt); % Non-specific MG genes.
%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'Pixel');
caxis([0,25000]);
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');

iss_change_plot(o,'Pixel');
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'all'), 'png');
iss_change_plot_allgenes_MG(o,'Pixel');
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'MG'), 'png');
iss_change_plot_allgenes_MG2(o,'Pixel');
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'MG2'), 'png');
iss_change_plot_allgenes_MG3(o,'Pixel');
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'MG3'), 'png');

%% All Spots Clusters
SpotSetClustered = get_gene_clusters(o, 'Pixel');

iss_change_plot(o, 'Pixel', o.GeneNames, SpotSetClustered);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'clusters', 'all'), 'png');
iss_change_plot_allgenes_MG(o, 'Pixel', o.GeneNames, SpotSetClustered);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'clusters', 'MG'), 'png');
iss_change_plot_allgenes_MG2(o, 'Pixel', o.GeneNames, SpotSetClustered);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'clusters', 'MG2'), 'png');
iss_change_plot_allgenes_MG3(o, 'Pixel', o.GeneNames, SpotSetClustered);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'clusters', 'MG3'), 'png');

%% MG genes
iss_change_plot_allgenes_MG3(o,'Pixel', GeneNamesMG);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'all'), 'svg');

iss_change_plot_MG(o, 'Pixel', GeneNamesMG);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG', 'all'), 'svg');
iss_change_plot_MG(o, 'Pixel', GeneNamesMGFilt);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG', 'specific'), 'svg');
iss_change_plot_MG(o, 'Pixel', GeneNamesMGFilt2);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG', 'non-specific'), 'svg');

iss_change_plot_MG2(o, 'Pixel', GeneNamesMG);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG2', 'all'), 'svg');
iss_change_plot_MG2(o, 'Pixel', GeneNamesMGFilt);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG2', 'specific'), 'svg');
iss_change_plot_MG2(o, 'Pixel', GeneNamesMGFilt2);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG2', 'non-specific'), 'svg');

iss_change_plot_MG3(o, 'Pixel', GeneNamesMG);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG3', 'all'), 'svg');
iss_change_plot_MG3(o, 'Pixel', GeneNamesMGFilt);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG3', 'specific'), 'svg');
iss_change_plot_MG3(o, 'Pixel', GeneNamesMGFilt2);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG3', 'non-specific'), 'svg');

%% MG Clusters
spot_no = find(ismember(o.GeneNames, GeneNamesMG));
spots_to_cluster = o.quality_threshold('Pixel') & ismember(o.pxSpotCodeNo, spot_no);
k = [2;3;5;10;15;25;50];

for i = 1:length(k)
    SpotSetClustered = get_gene_clusters(o, 'Pixel', 75, k(i), spots_to_cluster);

    iss_change_plot_MG(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG', 'all', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
    iss_change_plot_MG2(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG2', 'all', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
    iss_change_plot_MG3(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG3', 'all', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
end

spots_to_cluster = find(ismember(o.GeneNames, GeneNamesMGFilt));
spots_to_cluster = o.quality_threshold('Pixel') & ismember(o.pxSpotCodeNo, spots_to_cluster);
k = [2;3;5;10;15];

for i = 1:length(k)
    SpotSetClustered = get_gene_clusters(o, 'Pixel', 75, k(i), spots_to_cluster);

    iss_change_plot_MG(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG', 'specific', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
    iss_change_plot_MG2(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG2', 'specific', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
    iss_change_plot_MG3(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG3', 'specific', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
end

spots_to_cluster = find(ismember(o.GeneNames, GeneNamesMGFilt2));
spots_to_cluster = o.quality_threshold('Pixel') & ismember(o.pxSpotCodeNo, spots_to_cluster);
k = [2;3;5;10;15;25;50];

for i = 1:length(k)
    SpotSetClustered = get_gene_clusters(o, 'Pixel', 75, k(i), spots_to_cluster);

    iss_change_plot_MG(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG', 'non-specific', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
    iss_change_plot_MG2(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG2', 'non-specific', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
    iss_change_plot_MG3(o, 'Pixel', GeneNamesMG, SpotSetClustered);
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'MG3', 'non-specific', 'clusters', ...
        strcat('k', num2str(k(i)))), 'svg');
end

%% Individual Genes
% All genes.
%%%%%%%%%%%%
for i = 1:length(GeneNamesAll)
    iss_change_plot_individual(o, 'Pixel', GeneNamesAll(i));
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'allgenes', 'individual', ...
        GeneNamesAll{i}), 'svg');
end
%%%%%%%%%%%%

% MG genes.
%%%%%%%%%%%
for i = 1:length(GeneNamesMG)
    iss_change_plot_individual(o, 'Pixel', GeneNamesMG(i));
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'MG', 'individual', GeneNamesMG{i}), 'svg');
end
%%%%%%%%%%%

%% Spatial Cross-correlation
thresh = o.quality_threshold('Pixel');
spots = o.pxSpotGlobalYX(thresh,:);
codes = o.pxSpotCodeNo(thresh,:);
subset = size(spots);
rng(1);
subset = randsample(subset(1), subset(1));
sub_spots = spots(subset,:);
sub_codes = codes(subset,:);
[ccg_out, Pairs, gs, cum_dens, rel_cum_dens, ripley, local_density] = ...
    CCG_2d(sub_spots, sub_codes, 200, 30);

%% Contour Plots
% All genes individual.
%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(GeneNamesAll)
    iss_change_plot_contour(o, 'Pixel', GeneNamesAll(i));
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'allgenes', 'individual', ...
        GeneNamesAll{i}), 'svg');
end
%%%%%%%%%%%%%%%%%%%%%%%

% MG genes individual.
%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(GeneNamesMG)
    iss_change_plot_contour(o, 'Pixel', GeneNamesMG(i));
    saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'individual', ...
        GeneNamesMG{i}), 'svg');
end
%%%%%%%%%%%%%%%%%%%%%%

% MG all.
%%%%%%%%%
iss_change_plot_MG2_contour(o, 'Pixel', GeneNamesMG);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'all'), 'svg');
iss_change_plot_MG2_contour(o, 'Pixel', GeneNamesMGFilt);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'specific'), 'svg');
iss_change_plot_MG2_contour(o, 'Pixel', GeneNamesMGFilt2);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'non-specific'), 'svg');
%%%%%%%%%

% MG clusters all.
%%%%%%%%%%%%%%%%%%
genes = {'Laptm5', 'Sparc', 'Csf1r', 'Ptk2b'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'all', 'c1'), ...
    'svg');

genes = {'Grn', 'Tmem119', 'Bin1', 'Cyfip1', 'Plcg2', 'P2ry12', 'Ccr5', 'Cx3cr1'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'all', 'c2'), ...
    'svg');

genes = {'C1qB', 'C1qC', 'C1qa', 'Olfml3'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'all', 'c3'), ...
    'svg');

genes = {'Pld3', 'Pld4', 'Bin2', 'Atp6v0d2'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'all', 'c4'), ...
    'svg');
%%%%%%%%%%%%%%%%%%

% MG clusters specific.
%%%%%%%%%%%%%%%%%%%%%%%
genes = {'Csf1r'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'specific', ...
    'c1'), 'svg');

genes = {'Tmem119', 'P2ry12', 'Cx3cr1'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'specific', ...
    'c2'), 'svg');

genes = {'C1qB', 'C1qC', 'C1qa', 'Olfml3'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'specific', ...
    'c3'), 'svg');
%%%%%%%%%%%%%%%%%%%%%%%

% MG clusters non-specific.
%%%%%%%%%%%%%%%%%%%%%%%%%%%
genes = {'Laptm5', 'Sparc','Ptk2b'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'non-specific', ...
    'c1'), 'svg');

genes = {'Grn', 'Bin1', 'Cyfip1', 'Plcg2', 'Ccr5'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'non-specific', ...
    'c2'), 'svg');

genes = {'Pld3', 'Pld4', 'Bin2', 'Atp6v0d2'};
iss_change_plot_MG2_contour(o, 'Pixel', genes);
saveas(gcf, fullfile('results', 'figures', 'pilot', 'contour', 'MG', 'clusters', 'non-specific', ...
    'c4'), 'svg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Export Spot Codes
Thresh = o.quality_threshold('Pixel');
SpotCode = o.pxSpotCodeNo(Thresh,:);
writematrix(SpotCode, fullfile('assets', 'codebooks', 'pilot', 'pilot-SpotCode.txt'));
