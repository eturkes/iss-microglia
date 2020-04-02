%% Parameters that should be checked before each run
%CHECK BEFORE EACH RUN

o = iss;
o.AnchorChannel =  7;            %Channel that has most spots in anchor round (o.ReferenceRound)
o.DapiChannel = 1;              %Channel in anchor round that contains Dapi images
o.InitialShiftChannel = 7;      %Channel to use to find initial shifts between rounds
o.ReferenceRound = 8;           %Index of anchor round
o.RawFileExtension = '.nd2';    %Format of raw data

%% File Names
%CHECK BEFORE EACH RUN
o.InputDirectory = pwd;     %Folder path of raw data

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

o.TileDirectory = fullfile(pwd, SliceNb, 'tiles');
o.OutputDirectory = fullfile(pwd, SliceNb, 'output');
%Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
%the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
o.CodeFile = fullfile('codebook_comb.txt');

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
try
    o = o.extract_and_filter_NoGPU;
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
%run code
o.CallSpotsCodeNorm = 'WholeCode';      %Other alternative is 'Round'
o = o.call_spots;
o = o.call_spots_prob;
save(fullfile(o.OutputDirectory, 'oCall_spots'), 'o', '-v7.3');

%% plot results

o.CombiQualThresh = 0.7;
BigDapiImage = imread(o.BigDapiFile);
o.plot(BigDapiImage);
set(gca, 'YDir', 'reverse')
set(gca, 'XDir', 'reverse')

%iss_view_codes(o,234321,1);
o.pIntensityThresh = 100;
o.pScoreThresh = 10;
iss_change_plot(o,'Prob');
%iss_view_prob(o,234321,1);
%iss_change_plot(o,'DotProduct');

