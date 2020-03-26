%    This file is part of iss-scripts.
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

%% extract and filter

%parameters
o = iss;
o.nRounds = 7;
o.nExtraRounds = 1;         %Treat Anchor channel as extra round
o.InputDirectory = pwd;     %Folder path of raw data

SliceNb = 'Round0-6_SplitAnchor00';
mkdir(fullfile(pwd, SliceNb));

%FileBase{r} is the file name of the raw data of round r in o.InputDirectory
o.FileBase = cell(o.nRounds+o.nExtraRounds,1);
o.FileBase{1} = strcat(SliceNb, '1');
o.FileBase{2} = strcat(SliceNb, '2');
o.FileBase{3} = strcat(SliceNb, '3');
o.FileBase{4} = strcat(SliceNb, '4');
o.FileBase{5} = strcat(SliceNb, '5');
o.FileBase{6} = strcat(SliceNb, '6');
o.FileBase{7} = strcat(SliceNb, '7');
o.FileBase{8} = strcat(SliceNb, '8');    %Make sure the last round is the anchor

o.RawFileExtension = '.nd2';
o.TileDirectory = fullfile(pwd, SliceNb, 'tiles');
mkdir(o.TileDirectory);
o.DapiChannel = 1;
o.AnchorChannel =  7;    %Channel that has most spots in anchor round
o.ReferenceRound = 8;
o.FirstBaseChannel = 1;
o.OutputDirectory = fullfile(pwd, SliceNb, 'output');
mkdir(o.OutputDirectory);
o.bpLabels = {'0', '1', '2', '3','4','5','6'}; %order of bases

%These specify the dimensions of the filter. R1 should be approximately the
%radius of the spot and R2 should be double this.
o.ExtractR1 = 'auto';
o.ExtractR2 = 'auto';

o.ExtractScale = 'auto';
o.TilePixelValueShift = 15000;

%Max time (seconds) to wait for raw .nd2 files to be obtained
o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
o.MaxWaitTime = 21600;  

%run code
try
    o = o.extract_and_filter;
catch
    o = o.extract_and_filter_NoGPU;
end
save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');

%% register
o.AutoThresh(:,o.AnchorChannel,o.ReferenceRound) = o.AutoThresh(:,o.AnchorChannel,o.ReferenceRound)*0.25;     %As Anchor Threshold seemed too high
%parameters
o.TileSz = 2048;

%Anchor spots are detected in register2
o.DetectionRadius = 2;
o.SmoothSize = 0;     
o.IsolationRadius1 = 4;
o.IsolationRadius2 = 14;

o.DetectionThresh = 'auto';
o.ThreshParam = 5;
o.MinThresh = 10;
o.minPeaks = 1;
o.InitalShiftAutoMinScoreParam=3;   %a lower value will make it quicker but more likely to fail

%paramaters to find shifts between overlapping tiles
o.RegMinScore = 'auto';     
o.RegStep = [5,5];
o.RegSearch.South.Y = -1900:o.RegStep(1):-1700;
o.RegSearch.South.X = -50:o.RegStep(2):50;
o.RegSearch.East.Y = -50:o.RegStep(1):50;
o.RegSearch.East.X = -1900:o.RegStep(2):-1700;
o.RegWidenSearch = [50,50]; 

%run code
o = o.register2;
save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');

%% find spots

%parameters
o.nBP = 7;

%If a channel or round is faulty, you can ignore it by selecting only the
%good ones in o.UseChannels and o.UseRounds.
o.UseChannels = 1:o.nBP;
o.UseRounds = 1:o.nRounds;

%Search paramaters
o.InitialShiftChannel = 7;      %Channel to use to find initial shifts between rounds
o.FindSpotsMinScore = 'auto';
o.FindSpotsStep = [5,5];
%FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
%structs - have a different range for each round: 
%o.FindSpotsSearch = cell(o.nRounds,1);
o.FindSpotsSearch.Y = -100:o.FindSpotsStep(1):100;
o.FindSpotsSearch.X = -100:o.FindSpotsStep(2):100;
%Make WidenSearch larger if you think you have a large shift between rounds
o.FindSpotsWidenSearch = [50,50]; 

o.PcDist = 3; 
o.MinPCMatches = 1;    %HACK SO IT GETS TO THE END

%run code
o = o.find_spots2;
save(fullfile(o.OutputDirectory, 'oFind_spots'), 'o', '-v7.3');

%% call spots

%parameters
%Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
%the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
o.CodeFile = fullfile('codebook_comb.txt');

%run code
o.CallSpotsCodeNorm = 'WholeCode';      %Other alternative is 'Round'
o = o.call_spots;
o = o.call_spots_prob;
save(fullfile(o.OutputDirectory, 'oCall_spots'), 'o', '-v7.3');

%% plot results

o.CombiQualThresh = 0.7;
BigDapiImage = imread(o.BigDapiFile);
o.plot(BigDapiImage);

%iss_view_codes(o,234321,1);
%o.pIntensityThresh = 100;
%o.pScoreThresh = 10;
%iss_change_plot(o,'Prob');
%iss_view_prob(o,234321,1);
%iss_change_plot(o,'DotProduct');

