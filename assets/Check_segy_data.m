function [dt,nsamples,ntraces]=Check_segy_data(Files)
% Check if the following data are the same for all input seismic volumes:
% Sampling interval.
% Number of samples.
% Number of traces in the seismic volume.
% Coordinates of the first and last trace.
% 
% Otherwise, the execution must be interrupted. It doesn't make sense if 
% it is one file, it will generate error.
% "FILE" Matrix is a cell array that contains the names of the input files.
n_file=size(Files,2);

% Segy_H: Segy Header of seismic volumes.
% Segy_TH: Segy Trace Header of seismic traces.
% DaTa: Matrix containing the parameters to compare.
% "condition1": Check if the sampling interval, number of samples
% and number of traces are the same.
% "condition2": Check if the coordinates are the same.

DaTa = zeros(n_file,3);
for i=1:n_file
    [Segy_H] = ReadSegyHeader(Files{i});
    DaTa(i,:) = [Segy_H.dt/1000 ; Segy_H.ns ; Segy_H.ntraces];
end

if n_file >= 2
    condition1=zeros(n_file-1,1);
    for i = 1 : n_file-1
        condition1(i)=isequal(DaTa(i,:),DaTa(i+1,:));
    end

    if all(condition1)==false;
        error('Error en el SegyHeader de los volúmenes sísmicos');
    end
end

% The following data must be the same for each seismic volume.

% Coordinates of the first trace.
cdp_xy=zeros(n_file,2,2);
for i=1:n_file 
    [Data , Segy_TH] = ReadSegy(Files{i},'traces',1,'SkipData',1);
    cdp_xy(i,:,1) = [Segy_TH.cdpX , Segy_TH.cdpY];
end

% Coordinates of the last trace.
for i=1:n_file 
    [Data,Segy_TH]=ReadSegy(Files{i},'traces',DaTa(1,3),'SkipData',1);
    cdp_xy(i,:,2) = [Segy_TH.cdpX , Segy_TH.cdpY];
end

if n_file >= 2
    condition2=zeros(n_file-1,2);
    for j=1:2    
        for i = 1 : n_file-1
            condition2(j,i)=isequal(cdp_xy(i,:,j),cdp_xy(i+1,:,j));
        end
    end
        
    if all(condition2)==false;
        error('Error en las coordenadas de los volúmenes sísmicos');
    end
end

dt= DaTa(1,1);
nsamples= DaTa(1,2);
ntraces= DaTa(1,3);