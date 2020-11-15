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
DataName = 'NLF-150genes';
Sample = 'CTRL-oldRound4-noChan1';
SampleShort = 'CTRL';
Pipeline = 'same-pipeline';

load(fullfile('results', 'data', DataName, Pipeline, Sample, 'o.mat'));
o.BigDapiFile = fullfile('results', 'data', DataName, Pipeline, Sample, ...
    'background_image.tif');
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

o.pIntensityThresh = max(o.pSpotIntensity) + 1;
o.pScoreThresh = 7.5;

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'Pixel');
caxis([0,750]);
set(gca, 'XTick', [], 'YTick', [], 'XDir', 'reverse', 'YDir', 'reverse')
iss_change_plot_MG150(o, 'Pixel');
% saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, ...
%     strcat(SampleShort, '-allgenes')), 'svg');

%% Individual Genes
for i = 1:length(o.GeneNames)
    if i ~= 106
        iss_change_plot_individual_MG150(o, 'Pixel', o.GeneNames(i));
        saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'individual', ...
            strcat(SampleShort, '-', o.GeneNames{i})), 'svg');
    end
end

%% Contour Plots
for i = 1:length(o.GeneNames)
    if i ~= 106 && i ~= 62 && i ~= 118 && i ~= 139
        iss_change_plot_contour_individual_MG150(o, 'Pixel', o.GeneNames(i));
        saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'contour', ...
            strcat(SampleShort, '-', o.GeneNames{i})), 'svg');
    end
end

%% Plots By Gene Group
for i = 1:6
    % Get genes from each gene group.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = fopen(fullfile('results', 'R', 'new-clusters', strcat('group', num2str(i), '.txt')), 'r');
    tmp = textscan(fp, '%s %s', inf);
    genes = tmp{1};
    fclose(fp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iss_change_plot_group_MG150(o, 'Pixel', genes);
    saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, ...
        'by-gene-group', strcat(SampleShort, '-', 'group', num2str(i))), 'svg');
end

%% Contour Plots By Gene Group
for i = 1:6
    % Get genes from each gene group.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = fopen(fullfile('results', 'R', 'new-clusters', strcat('group', num2str(i), '.txt')), 'r');
    tmp = textscan(fp, '%s %s', inf);
    genes = tmp{1};
    fclose(fp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iss_change_plot_contour_MG150(o, 'Pixel', genes);
    saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, ...
        'by-gene-group', 'contour', strcat(SampleShort, '-', 'group', num2str(i))), 'svg');
end

%% Export Spot Codes
Thresh = o.quality_threshold('Pixel');
SpotCode = o.pxSpotCodeNo(Thresh,:);
writematrix(SpotCode, fullfile('assets', 'codebooks', DataName, ...
    strcat(SampleShort, '-', 'SpotCode.txt')));

%% Subset to SLM
clear all;

DataName = 'NLF-150genes';
Sample = 'CTRL-oldRound4-noChan1';
SampleShort = 'CTRL';
Pipeline = 'same-pipeline';

load(fullfile('results', 'data', DataName, Pipeline, Sample, 'o.mat'));

theta = 20;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
o.SpotGlobalYX = o.SpotGlobalYX * R';
o.SpotGlobalYX(:,2) = o.SpotGlobalYX(:,2) - (min(o.SpotGlobalYX(:,2)) - 1);
o.SpotGlobalYX(:,1) = o.SpotGlobalYX(:,1) - (min(o.SpotGlobalYX(:,1)) - 1);
o.pxSpotGlobalYX = o.pxSpotGlobalYX * R';
o.pxSpotGlobalYX(:,2) = o.pxSpotGlobalYX(:,2) - (min(o.pxSpotGlobalYX(:,2)) - 1);
o.pxSpotGlobalYX(:,1) = o.pxSpotGlobalYX(:,1) - (min(o.pxSpotGlobalYX(:,1)) - 1);

o.BigDapiFile = fullfile('results', 'data', DataName, Pipeline, Sample, ...
    'background_image.tif');
o.pIntensityThresh = max(o.pSpotIntensity) + 1;
o.pScoreThresh = 7.5;

% Roi = round([min(o.SpotGlobalYX(:,2)), max(o.SpotGlobalYX(:,2)), min(o.SpotGlobalYX(:,1)), ...
%     max(o.SpotGlobalYX(:,1))]);
% o.plot(o.BigDapiFile, Roi, 'Pixel');
% caxis([0,750]);
% set(gca, 'XTick', [], 'YTick', [], 'XDir', 'reverse', 'YDir', 'reverse')
% iss_change_plot_MG150(o, 'Pixel');

load(fullfile('results', 'data', DataName, Pipeline, Sample, 'xs.mat'));
load(fullfile('results', 'data', DataName, Pipeline, Sample, 'ys.mat'));

xs = vertcat(xs{:});
ys = vertcat(ys{:});
xSub = ismember(o.pxSpotGlobalYX(:,2), xs);
ySub = ismember(o.pxSpotGlobalYX(:,1), ys);
SpotSet = o.quality_threshold('Pixel') & xSub > 0 & ySub > 0;

o.pxSpotGlobalYX(:,2) = o.pxSpotGlobalYX(:,2) - (min(xs) - 1);
o.pxSpotGlobalYX(:,1) = o.pxSpotGlobalYX(:,1) - (min(ys) - 1);
Spots = o.pxSpotGlobalYX(SpotSet,:);
Roi = round([min(Spots(:,2)), max(Spots(:,2)), min(Spots(:,1)), max(Spots(:,1))]);

o.plot(o.BigDapiFile, Roi, 'Pixel');
caxis([0,750]);
set(gca, 'XTick', [], 'YTick', [], 'XDir', 'reverse', 'YDir', 'reverse')
iss_change_plot_MG150(o, 'Pixel', o.GeneNames, SpotSet);
saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'SLM', ...
    strcat(SampleShort, '-allgenes')), 'svg');

%% Individual Genes
for i = 1:length(o.GeneNames)
    if i ~= 106
        iss_change_plot_individual_MG150(o, 'Pixel', o.GeneNames(i), SpotSet);
        saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'SLM', ...
            'individual', strcat(SampleShort, '-', o.GeneNames{i})), 'svg');
    end
end

%% Contour Plots
for i = 1:length(o.GeneNames)
    if i ~= 106 && i ~= 62 && i ~= 118 && i ~= 139
        iss_change_plot_contour_individual_MG150(o, 'Pixel', o.GeneNames(i), SpotSet);
        saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'SLM', ...
            'contour', strcat(SampleShort, '-', o.GeneNames{i})), 'svg');
    end
end

%% Plots By Gene Group
for i = 1:6
    % Get genes from each gene group.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = fopen(fullfile('results', 'R', 'new-clusters', strcat('group', num2str(i), '.txt')), 'r');
    tmp = textscan(fp, '%s %s', inf);
    genes = tmp{1};
    fclose(fp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iss_change_plot_group_MG150(o, 'Pixel', genes, SpotSet);
    saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'SLM', ...
        'by-gene-group', strcat(SampleShort, '-', 'group', num2str(i))), 'svg');
end

%% Contour Plots By Gene Group
for i = 1:6
    % Get genes from each gene group.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fp = fopen(fullfile('results', 'R', 'new-clusters', strcat('group', num2str(i), '.txt')), 'r');
    tmp = textscan(fp, '%s %s', inf);
    genes = tmp{1};
    fclose(fp);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    iss_change_plot_contour_MG150(o, 'Pixel', genes, SpotSet);
    saveas(gcf, fullfile('results', 'figures', DataName, Pipeline, Sample, 'SLM', ...
        'by-gene-group', 'contour', strcat(SampleShort, '-', 'group', num2str(i))), 'svg');
end

%% Export Spot Codes
Thresh = o.quality_threshold('Pixel') & xSub > 0 & ySub > 0;
SpotCode = o.pxSpotCodeNo(Thresh,:);
writematrix(SpotCode, fullfile('assets', 'codebooks', DataName, ...
    strcat(SampleShort, '-', 'SpotCodeSLM.txt')));

SpotSetClustered = get_gene_clusters(o, 'Pixel', 75, 20, SpotSet);
iss_change_plot_MG150(o, 'Pixel', o.GeneNames, SpotSetClustered);
