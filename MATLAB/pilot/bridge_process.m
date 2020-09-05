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

%% Template of a pipeline that links together all the different processes 
%in the correct order.
%The first two sections below will be different for each experiment so
%should be checked before each run.

%% Parameters that should be checked before each run
%CHECK BEFORE EACH RUN

o = iss;
o.AnchorRound = 8;              %Round that contains Dapi image
o.AnchorChannel =  7;            %Channel that has most spots in anchor round (o.ReferenceRound)
o.DapiChannel = 1;              %Channel in anchor round that contains Dapi images
o.InitialShiftChannel = 7;      %Channel to use to find initial shifts between rounds
o.ReferenceRound = 8;           %Index of anchor round
o.ReferenceChannel = 7;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
o.RawFileExtension = '.nd2';    %Format of raw data
o.LogToFile = 1;                %Set to 1 if you want to save command window to txt file, else set to 0.

%% File Names
%CHECK BEFORE EACH RUN
o.InputDirectory = fullfile('assets', 'data', 'pilot');     %Folder path of raw data

%FileBase{r} is the file name of the raw data of round r in o.InputDirectory
SliceNb = 'Round0-6_SplitAnchor00';
o.FileBase = cell(o.nRounds+o.nExtraRounds,1);
o.FileBase{1} = strcat(SliceNb, '1');
o.FileBase{2} = strcat(SliceNb, '2');
o.FileBase{3} = strcat(SliceNb, '3');
o.FileBase{4} = strcat(SliceNb, '4');
o.FileBase{5} = strcat(SliceNb, '5');
o.FileBase{6} = strcat(SliceNb, '6');
o.FileBase{7} = strcat(SliceNb, '7');
o.FileBase{8} = strcat(SliceNb, '8');    %Make sure the last round is the anchor

o.TileDirectory = fullfile('cache', 'pilot', 'tiles');
o.OutputDirectory = fullfile('cache', 'pilot', 'output');
%Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
%the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
o.CodeFile = fullfile('assets', 'codebooks', 'pilot', 'codebook_comb.txt');

%% Logging
if o.LogToFile
    if isempty(o.LogFile)
        o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
    end
end

%% extract and filter

%parameters
o.nRounds = 7;
o.nExtraRounds = 1;         %Treat Anchor channel as extra round
o.FirstBaseChannel = 1;
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
o = o.extract_and_filter;
save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');

%% register
%Specify tiles to run
%o.EmptyTiles(:) = 1;
%UseTiles = [1,2];
%o.EmptyTiles(UseTiles) = 0;

%o.AutoThresh(:,o.AnchorChannel,o.AnchorRound) = o.AutoThresh(:,o.AnchorChannel,o.AnchorRound)*0.25;     %As Anchor Threshold seemed too high
%parameters
%o.RegMethod = 'PointBased';
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

%If a channel or round is faulty, you can ignore it by selecting only the
%good ones in o.UseChannels and o.UseRounds.
o.nBP = 7;
o.UseChannels = 1:o.nBP;
o.UseRounds = 1:o.nRounds;

%run code
o = o.register2;
save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');

%% find spots

%Search paramaters
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

%run code
if isprop(o,'FindSpotsMethod') && strcmpi(o.FindSpotsMethod, 'Fft')
    o = o.find_spots_FFt;
else
    o = o.find_spots2;
end
save(fullfile(o.OutputDirectory, 'oFind_spots'), 'o', '-v7.3');

%% call spots
%run code
o.CallSpotsCodeNorm = 'WholeCode';      %Other alternative is 'Round'
o = o.call_spots;
[o,LookupTable] = o.call_spots_prob;
save(fullfile(o.OutputDirectory, 'oCall_spots'), 'o', '-v7.3');

%Pixel based
o.PixelFileMaxTiles = 3;
o = o.call_spots_pixel(LookupTable);
save(fullfile(o.OutputDirectory, 'oCall_spots_pixel'), 'o', '-v7.3');

%% plot results

o.CombiQualThresh = 0.7;
Roi = round([1, max(o.SpotGlobalYX(:,2)), ...
1, max(o.SpotGlobalYX(:,1))]);
o.plot(o.BigAnchorFile,Roi,'Prob');
%o.plot(o.BigAnchorFile,Roi,'Pixel');

%iss_view_codes(o,234321,1);
%o.pIntensityThresh = 100;
%o.pScoreThresh = 10;
%iss_change_plot(o,'Prob');
%iss_change_plot(o,'Pixel');
%iss_view_prob(o,234321,1);
%iss_change_plot(o,'DotProduct');
