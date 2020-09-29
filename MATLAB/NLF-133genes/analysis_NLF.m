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
DataName = 'NLF-133genes';
Sample = 'NLF';

load(fullfile('results', 'data', DataName, Sample, 'o.mat'));
Roi = round([1, max(o.SpotGlobalYX(:,2)), 1, max(o.SpotGlobalYX(:,1))]);

% Get genes from codebook.
%%%%%%%%%%%%%%%%%%%%%%%%%%
fp = fopen(fullfile('assets', 'codebooks', DataName, 'codebook-iss-cleaned.tsv'), 'r');
tmp = textscan(fp, '%s %s', inf);
GeneNames = tmp{1};
fclose(fp);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot All Spots
o.plot(o.BigDapiFile, Roi, 'DotProduct');
set(gca, 'YDir', 'reverse');
set(gca, 'XDir', 'reverse');
iss_change_plot_MG133(o, 'DotProduct');
saveas(gcf, fullfile('results', 'figures', DataName, Sample, 'all'), 'svg');

%% Individual Genes
for i = 1:length(GeneNames)
    iss_change_plot_individual_MG133(o, 'DotProduct', GeneNames(i));
    saveas(gcf, fullfile('results', 'figures', DataName, Sample, 'individual', GeneNames{i}), ...
        'svg');
end
