clc
clear all

%% ################################################## %%
%  INPUT DATA: Seismic Volume and Surface of interest  %
% ###################################################  %

% Definition of time window, above and below the horizon of interest.
TIME_ABOVE_HORIZON = 0;
TIME_BELOW_HORIZON = 20;

% Array of strings containing the names of the seismic attribute files.
files = {'Data/Inline 900-903/seismic_amplitude_IL_900-903.segy', ...
         'Data/Inline 900-903/instantaneous_amplitude_IL_900-903.segy', ...
         'Data/Inline 900-903/instantaneous_frecuency_IL_900-903.segy'};

% Array of strings containing the titles of the input seismic attributes.
titlesFiles = {'Seismic amplitude', 'Instantaneous amplitude', 'Instantaneous frequency'};

% Definition of the weights for each attributes. It must be the same
% quantity as the input attributes and their sum must be equal to 1.
% These are used to control the contribution of each seismic attribute
% in the calculation of Fractal Dimension.
ALPHA = [0.33, 0.33, 0.33];

% File to save the calculation of Fractal Dimension.
fileFractalDim = {'Data/EM_MER-U2_Fractal_Dimension.prn'};

% Number of traces to perform Quality Control (QC )of the calculations.
NUMTRACESQC = 5;

% Minimum and Maximum values of input seismic data. Extracted from Textual Header.
MINMAX = {[-7469.42 6192.14], [0.0 5228.60], [0.0 125.0]};
%% ######### %%
%  Load data  %
%  #########  %

% ------------------------------ SURFACE ------------------------------- %

% Opening of the time surface of interest. The format must be X - Y - Time.
% Otherwise, it would generate an error.
% X-Y are coordinates.
assy = fopen('Data/EM_MER-U2.prn');

switch assy
    case -1
        errordlg('The surface file is not in the current directory', ...
            'Error!');
    otherwise
        try
            reg = fscanf(assy,'%10f');
            surfaceData = reshape(reg,3,[])';
        catch id
            errordlg('The file does not have the format X-Y-Time',...
                'Error!');
        end 
end

% ---------------------------- SESMIC DATA ----------------------------- %

% Number of attributes to use in the algorithm
nAttributes = size(files, 2);

% Raise error dialog if the sum of the weights is different from one,
% or they differ from the number of input attributes.
if (sum(ALPHA) < 0.99)
    error('The sum of weights are not equal to 1');
    
elseif (nAttributes ~= size(ALPHA,2))
    error('The number of weights is not equal to input attributes.'); 
    
end

% The following data must be the same for each seismic volume.
% dtSample: Sampling interval (microseconds).
% nSamples: Number of samples.
% nTraces: Number of traces.

% If there was selected only one seismic volume, no verification need it.
if nAttributes >= 2
    [dtSample, nSamples, nTraces] = checkSegyData(files);

else
    [segyHeader] = ReadSegyHeader(files{1});
    dtSample= segyHeader.dt/1000;
    nSamples= segyHeader.ns;
    nTraces = segyHeader.ntraces;

end

%% ################## %%
%  Data preprocessing  %
%  ##################  %

% ----------------------------- TIME COLUMN ---------------------------- %

% Reade segy header
[segyHeader] = ReadSegyHeader(files{1});

% Read segy trace header
% 'SkipData': Read only the header values (Data will return empty)
[~ , segyTrace1Header] = ReadSegy(files{1}, 'traces', 1, 'SkipData', 1);

% From "segyHeader" read time column and sampling interval.
% Multiply and divide by 1000 to get miliseconds.
% From "segyTraceHeader" read "Delay Recording Time".
timeColumn = segyHeader.time*1000 - segyHeader.dt/1000 + segyTrace1Header.DelayRecordingTime;

% nTop and nBottom. Number of samples that defines the top and the bottom
% of time window analysis. extract seismic data for Fractal 
% Dimension analysis. It depends on the time window previously chosen.
nTop = ceil(TIME_ABOVE_HORIZON/dtSample);
nBottom = ceil(TIME_BELOW_HORIZON/dtSample);

% ----------------------------- SEISMIC DATA ---------------------------- %

% X-Y coordinates of all traces in the seismic volume. These coordinates
% will be the same in all input seismic volumes, because the sampling
% interval, the number of samples and the number of traces are the same.
% Also, they are divided by 100 to bring them to metric scale.

% Open and read all trace headers from segy data
[~ , allTraceHeaders] = ReadSegy(files{1}, 'SkipData',1);

% Coordinates of seismic traces and Number of traces
coordinatesSeismicTraces = zeros(nTraces, 3);

for iTrace = 1:nTraces
    % Divide each coordinate by 100 to bring them to metric scale
    coordinatesSeismicTraces(iTrace, :) = [round(allTraceHeaders(iTrace).SourceX / 100), ...
                                           round(allTraceHeaders(iTrace).SourceY / 100), ...
                                           allTraceHeaders(iTrace).TraceNumber];
end


% Select X-Y coordinates and the number of seismic traces where surface
% data was interpolated. 
selectedCoordinates = selectEqualCoord(nTraces,...
                                       surfaceData,...
                                       coordinatesSeismicTraces);

% ----------------------------- SURFACE DATA ---------------------------- %

% Select surface data with coordinates equal to seismic data.
selectedSurface = selectEqualCoord(nTraces, selectedCoordinates, surfaceData);

% ------------------------ PREALLOCATING MATRICES ----------------------- %

%%%% nPointsSurface changed by nTracesSelected
% Number of points on the surface of interest.
% nPointsSurface = size(surfaceData, 1);

% Number of samples selected to Fractal dimension calculations.
nTracesSelected = size(selectedCoordinates, 1);

% This matrix will store the results of Fractal Dimension calculations.
resultFractalDim = zeros(nTracesSelected, 3);

% This vector will store the values ??of R2 in order to make a histogram.
allR2 = zeros(nTracesSelected, 1);

% Definition of random numbers between 1 and the number of traces in
% seismic volumes. This will represent random traces, in order to perform
% Quality Control of Fractal Dimension calculations.
rng(3, 'twister');
tracesQC = round(rand(1, NUMTRACESQC) * nTraces);

% 1. Save seismic attributes traces.
seismicTraces = zeros(nTop+nBottom+1, nAttributes, NUMTRACESQC);

% 2. Save the number of traces, slope, Fractal Dimension and the
% coefficient of determination R2
params = zeros(NUMTRACESQC, 4);

% 3. Save traces coordinates for quality control.
coordinatesQC = zeros(NUMTRACESQC, 2);

% 4. Time interval analysis.
timeAnalysis =  zeros(nTop+nBottom+1, NUMTRACESQC);

% 5. Save divisor (r) and total length (L). The length of both vectors is
% unknown, so pre-allocation is made with a single column.
divisorLength = zeros(2, 1, NUMTRACESQC);

% Number of the traces used to calculations (the coordinates
% of the surface is equal to the coordinate of the seismic trace).
% traceNumbers = zeros(1, nTracesSelected) ;

%% ############################################ %%
%  Multi-attributes Fractal Dimension Algorithm  %
%  ############################################  %

% Initialize waitbar
wb = waitbar(0, 'Analyzing data');

% traceNum = index that locates the trace number read.
traceNum = 1;

% conditionQC: Quality Control Condition. Sets the selection of
% the random trace to save its corresponding data.
conditionQC = false;

% traceQC: index that indicates the trace number to QC.
traceQC = 1;

for iPoint = 1 : nTracesSelected
    
    % Plot waitbar
    waitbar(iPoint/nTracesSelected, wb)
    
%     % Find the position where the surface coordinate is equal to seismic
%     % trace coordinate.
%     condition_time_suface = false;
%     while condition_time_suface == false && traceNum < nTraces
%         condition_time_suface = isequal(surfaceData(iPoint,1:2), coordinatesSeismicTraces(traceNum,1:2));
%         traceNum = traceNum + 1;
%     end
% 
%     % When traceNum>nTraces the search has reached the end of seismic data. This
%     % can also happen if the coordinates of the area of ??interest are not
%     % within the seismic data. In either case, program execution must stop.
%     
%     if traceNum >= nTraces
%         errordlg('Finalizó la búsqueda dentro de los datos sísmicos',...
%                 'Error!');
%         break
%     end
%     
%     traceNum = traceNum - 1;
%     
%     traceNumbers(iPoint) = traceNum;
    
    fprintf('Seismic Trace number: %d', selectedCoordinates(iPoint, 3));

    % Read trace number selected in seismic data.
    seismicData = zeros(nSamples, nAttributes);
    for iAttribute = 1 : nAttributes
        % Read Segy data of each input
        [Data, ~] = ReadSegy(files{iAttribute}, 'traces', selectedCoordinates(iPoint, 3));
        seismicData(:, iAttribute) = Data;
    end
    
    % Check if any of the seismic attributes data are empty. If this case,
    % it will not calculate Fractal Dimension.    
    isEmptyData= any(all(seismicData==0));

    % Continue calculations is seismic attributes data is not empty.
    if ~isEmptyData
        % Get time sample of surface data. "floor": Round towards -inf.
        timeSample = floor(selectedSurface(iPoint,3));

        % Check if "timeSample" is a multiple of sampling interval of seismic
        % data (4 miliseconds in this case).
        % Otherwise, take the next higher multiple.
        modulus = mod(timeSample, dtSample);
        while modulus~=0
            timeSample = timeSample - 1;
            modulus = mod(timeSample, dtSample);
        end

        % Find the position of "timeSample" in "timeColumn" vector.
        indexTimeSample = find(timeSample == timeColumn);
        
        % Select seismic attributes data between time window previously
        % defined. This will be the data under analysis.
        analyzedData = zeros(nTop+nBottom+1,nAttributes);

        for i = 1:nAttributes
            analyzedData(:,i) = seismicData(indexTimeSample - nTop : indexTimeSample + nBottom, i);
        end
        
        % Normalization of data
        normalizedData = normalization(analyzedData , analyzedData);
    
        % Multiplication of the normalized traces by previously defined 
        % weights.
        for i=1:nAttributes
            normalizedData(:,i)= normalizedData(:,i) .* ALPHA(i);
        end

        % The time column in the window of interest does not use for
        % Fractal Dimension calculations. It will be use only for QC.
        % Because the input data was normalized, the time (in the last 
        % column of the array), must also range from 0 to 1.
        
        normalizedData(:,end+1) = zeros(nTop+nBottom+1,1);

        div = 1/(size(normalizedData,1)-1);
        h=0;
        for i=1:size(normalizedData,1)
            normalizedData(i,end)=h*div;
            h=h+1;
        end
        
        % Fractal Dimension calculation
        [r,L,P,R2]=divplot2_EM(normalizedData);

        % Fractal Dimension of the curve in space N-Dimension.
        D = 1 + abs(P);
        
        % This matrix will keep the results of each FRACTAL DIMENSION
        % calculations. The order is X-Y-D.
        resultFractalDim(iPoint,:) = [surfaceData(iPoint,1) surfaceData(iPoint,2) D]; 

        % Save R2 values, in order to perform a histogram with the results.
        allR2(iPoint)=R2;
        
        % For QC, from 5 random traces, the following data will be saved:
        % Input traces, trace number, D, coordinates, r-L
        try
            conditionQC = any(tracesQC >= traceNumbers(iPoint-1) & tracesQC <= traceNumbers(iPoint));
        catch
            continue
        end
        
        if conditionQC
            % Input Trace: Each 3D from the matrix Tr_Analysis_norm_QC has
            % data traces.
            seismicTraces(:,:,traceQC) = analyzedData;
            
            % Each row will correspond to the data of each trace:
            % Trace Number - P - D - R2.
            params(traceQC,1) = traceNum; 
            params(traceQC,2) = P; 
            params(traceQC,3) = D; 
            params(traceQC,4) = R2;
            
            % Coordinates X/Y - r/L
            coordinatesQC(traceQC,:) = surfaceData(iPoint,1:2); 

            % Time column in window analysis (FD_QC_TIME), according to
            % time found in the surface of interest.            
            t_ini_analisis= timeSample - (nTop * dtSample);            
            for i = 1 : (nTop+nBottom+1)
                timeAnalysis(i,traceQC) = t_ini_analisis + dtSample*(i-1);
            end
            
            % Divider size "r" and total length "L". The 3rd dimension
            % of the hyper matrix corresponds to each trace.
            divisorLength(1,1:size(r,2),traceQC) = r; 
            divisorLength(2,1:size(r,2),traceQC) = L; 

            traceQC = traceQC + 1;
        end
        
        st = fclose('all');
    else  
        % Close all open segy files.
        st = fclose('all');
        continue
    end
    
    % Perhaps it would be necessary to establish a condition if it does not
    % obtain any coordinate, which mean that the surface and seismic volume
    % wouldn't be in the same order or do not have the same coordinates.

end 

%% Writing the results file and final graphics.

% Delete empty traces from Fractal Dimension matrix (coordinates equal to
% zero). It is possible that they generate some error in surface loading.

FD_erased = reshape(resultFractalDim,nPointsSurface,[])';
empty_coord = find(FD_erased==0);
FD_erased(empty_coord) = [];
FD_erased = reshape(FD_erased,3,[])';

% Fractal Dimension surface file.
dlmwrite(fileFractalDim{1},FD_erased,'delimiter',' ','precision','%15.5f')

% Graphics.
% It will be both plots in the same window. 
for i=1:size(params,1)
    graph_QC(seismicTraces , params , divisorLength , coordinatesQC , timeAnalysis,...
        titlesFiles,i);
end

% R2 Histogram. All possible empty cells are erased.
figure('Name','Histograma de R^2 en cálculo de Dimensión Fractal')
allR2(allR2==0)=[];
hist(allR2,(0.1:0.1:0.9));
ylabel('Frecuencia'), xlabel('R^2');
title('Histograma de Frecuencia Absoluta. Coeficiente de Determinación.');
