function [dt,nsamples,ntraces] = checkSegyData(filenames)
% Check if the following data are the same for all input seismic volumes:
% Sampling interval.
% Number of samples.
% Number of traces in the seismic volume.
% Coordinates of the first and last trace.
% 
% Otherwise, the execution must be interrupted.
% "filenames" Matrix is a cell array that contains the names of the input files.

% Description of variables:
% segyHeader:       Segy Header of seismic volumes.
% segyTraceHeader:  Segy Trace Header of seismic traces.
% params:           Matrix containing the parameters to compare.
% checkSampleTrace: Check if the sampling interval, number of samples
%                   and number of traces are the same.
% checkCoordinate:  Check if the coordinates are the same.

nFiles=size(filenames,2);

params = zeros(nFiles,3);
for iFiles = 1:nFiles
    [segyHeader] = ReadSegyHeader(filenames{iFiles});
    params(iFiles,:) = [segyHeader.dt/1000 ; segyHeader.ns ; segyHeader.ntraces];
end
%---------------------------------------------------
%---------------------------------------------------
%---------------------------------------------------
%---------------------------------------------------
% Initialize "checkSampleTrace" matrix.
checkSampleTrace = zeros(nFiles-1,1);
for i = 1 : nFiles-1
    checkSampleTrace(i)=isequal(params(i,:),params(i+1,:));
end

if all(checkSampleTrace)==false;
    error('Error en el SegyHeader de los volúmenes sísmicos');
end


% The following data must be the same for each seismic volume.

% Coordinates of the first trace.
cdp_xy=zeros(nFiles,2,2);
for i=1:nFiles 
    [Data , segyTraceHeader] = ReadSegy(filenames{i},'traces',1,'SkipData',1);
    cdp_xy(i,:,1) = [segyTraceHeader.cdpX , segyTraceHeader.cdpY];
end

% Coordinates of the last trace.
for i=1:nFiles 
    [Data,segyTraceHeader]=ReadSegy(filenames{i},'traces',params(1,3),'SkipData',1);
    cdp_xy(i,:,2) = [segyTraceHeader.cdpX , segyTraceHeader.cdpY];
end

if nFiles >= 2
    checkCoordinate=zeros(nFiles-1,2);
    for j=1:2    
        for i = 1 : nFiles-1
            checkCoordinate(j,i)=isequal(cdp_xy(i,:,j),cdp_xy(i+1,:,j));
        end
    end
        
    if all(checkCoordinate)==false;
        error('Error en las coordenadas de los volúmenes sísmicos');
    end
end

dt= params(1,1);
nsamples= params(1,2);
ntraces= params(1,3);