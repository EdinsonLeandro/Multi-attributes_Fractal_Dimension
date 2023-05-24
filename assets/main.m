clc
clear all

%% ################################################## %%
%  INPUT DATA: Seismic Volume and Surface of interest  %
% ###################################################  %

% Definition of time window, above and below the horizon of interest.
TIME_ABOVE_HORIZON = 10;
TIME_BELOW_HORIZON = 50;

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

% Minimum and Maximum values of input seismic data.
% Extracted from Textual Header.
MINMAX = [[-7469.42, 6192.14]; [0.0, 5228.60]; [0.0, 125.0]];

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
    errordlg('The sum of weights are not equal to 1', 'Error!'); 
    
elseif (nAttributes ~= size(ALPHA,2))
    errordlg('The number of weights is not equal to the input attributes.',...
        'Error!'); 
    
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

% Coordinates of seismic traces and trace numbers
coordinatesSeismicTraces = zeros(nTraces, 3);

% Read X and Y coordinates.
% <B = A.'> is equal to <B = transpose(A)>
coordinatesSeismicTraces(:,1) = ReadSegyTraceHeaderValue(files{1}, 'key', 'cdpX').';
coordinatesSeismicTraces(:,2) = ReadSegyTraceHeaderValue(files{1}, 'key', 'cdpY').';

% Divide each coordinate by 100 (metric scale)
coordinatesSeismicTraces = round(coordinatesSeismicTraces./100);

% Linear sequence of trace numbers.
coordinatesSeismicTraces(:, 3) = 1:nTraces;

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
results = zeros(nTracesSelected, 3);

% This vector will store the values of R2 in order to make a histogram.
allR2 = zeros(nTracesSelected, 1);

% Definition of random numbers between 1 and the number of traces in
% seismic volumes. This will represent random traces, in order to perform
% Quality Control of Fractal Dimension calculations.
rng(3, 'twister');
nTracesQc = round(rand(1, NUMTRACESQC) * nTraces);

% Matrices for Quality control.
% 1. Seismic attributes traces.
seismicTracesQc = zeros(nTop+nBottom+1, nAttributes, NUMTRACESQC);

% 2. Number of traces, slope, Fractal Dimension and the coefficient of
% determination R2.
paramsQc = zeros(NUMTRACESQC, 4);

% 3. Traces coordinates.
coordinatesQc = zeros(NUMTRACESQC, 2);

% 4. Time interval analysis.
timeAnalysisQc =  zeros(nTop+nBottom+1, NUMTRACESQC);

% 5. "Dividers" and "Length". Number of columns of both vectors is unknown,
% so pre-allocation is made with one column.
dividerLengthQc = zeros(2, 1, NUMTRACESQC);

%% ############################################ %%
%  Multi-attributes Fractal Dimension Algorithm  %
%  ############################################  %

% Close all open files
fclose('all');

% Initialize waitbar
wb = waitbar(0, 'Analyzing data');

% indexQc: index that indicates the trace number to Quality contol.
indexQc = 1;

for iPoint = 1 : nTracesSelected
    
    % Plot waitbar
    waitbar(iPoint/nTracesSelected, wb)

    % Read trace number selected in seismic data.
    seismicData = zeros(nSamples, nAttributes);
    for jAttribute = 1 : nAttributes
        % Read Segy data of each input
        [Data, ~] = ReadSegy(files{jAttribute}, 'traces', selectedCoordinates(iPoint, 3));
        seismicData(:, jAttribute) = Data;
    end
    
    % Check if any of the seismic attributes data are empty. If this case,
    % it will not calculate Fractal Dimension.    
    isEmptyData= any(all(seismicData==0));

    % Continue calculations is seismic attributes data is not empty.
    if ~isEmptyData
        % Get time sample of surface data. "floor": Round towards -Inf.
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
        analyzedData = zeros(nTop + nBottom + 1, nAttributes);

        for i = 1:nAttributes
            analyzedData(:,i) = seismicData(indexTimeSample - nTop : indexTimeSample + nBottom, i);
        end
        
        % Normalization of data
        normalizedData = minMaxScaler(analyzedData, MINMAX);
    
        % Weighed normalized data
        normalizedData = normalizedData .* ALPHA;

        % Fractal Dimension calculation
        [dividers, length, slope, R2] = divplot2_EM(normalizedData);

        % Fractal Dimension of the curve in space N-Dimension.
        fractalDim = 1 + abs(slope);
        
        % Output. The order is X-Y-fractalDim.
        results(iPoint,:) = [selectedCoordinates(iPoint,1), ...
                             selectedCoordinates(iPoint,2), ...
                             fractalDim]; 

        % Save R2 values, in order to perform a histogram with the results.
        allR2(iPoint) = R2;
        
        % Check if the trace under analysis is the selected for quality control
        if ~isempty(intersect(nTracesQc, selectedCoordinates(iPoint, 3)))
            % Input Trace: Each 3D from the matrix Tr_Analysis_norm_QC has
            % data traces.
            seismicTracesQc(:, :, indexQc) = analyzedData;
            
            % Each column will correspond to the data of each trace:
            % Number of trace - slope - fractalDim - R2.
            paramsQc(indexQc, 1) = selectedCoordinates(iPoint, 3); 
            paramsQc(indexQc, 2) = slope; 
            paramsQc(indexQc, 3) = fractalDim; 
            paramsQc(indexQc, 4) = R2;
            
            % Coordinates X/Y - dividers/length
            coordinatesQc(indexQc, :) = surfaceData(iPoint, 1:2); 

            % Time column under analysis, according to time found in the
            % surface of interest.            
            initialTime = timeSample - (nTop * dtSample);            
            for i = 1 : (nTop+nBottom+1)
                timeAnalysisQc(i,indexQc) = initialTime + dtSample*(i-1);
            end % End for
            
            % "Dividers" size and total "length". The 3rd dimension
            % of the hyper matrix corresponds to each trace.
            dividerLengthQc(1, 1:size(dividers,2), indexQc) = dividers; 
            dividerLengthQc(2, 1:size(dividers,2), indexQc) = length; 

            indexQc = indexQc + 1;
        end % End if

    % Close all open files.
    fclose('all');
    
    end % End if

end % End for

% Close waitbar
close(wb)

%% Results.

% Remove zero rows (empty traces) from Fractal Dimension matrix.
results( ~any(results, 2), : ) = [];

% Fractal Dimension surface file.
dlmwrite(fileFractalDim{1}, results, 'delimiter',' ','precision','%15.5f')

% Close all open files
fclose('all');

% Quality Control Charts.
% It will be both plots in the same window. 
for iParam = 1:size(paramsQc,1)
    qualityControlCharts(seismicTracesQc, paramsQc, dividerLengthQc, ...
        coordinatesQc, timeAnalysisQc, titlesFiles, iParam);
end

% R2 Histogram.
qualityControlHist(allR2)