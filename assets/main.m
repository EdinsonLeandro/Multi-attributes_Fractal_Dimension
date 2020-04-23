clc
clear all

%% INPUT DATA: Seismic Volume and Surface of interest %%
% #################################################### %

% Definition of time window, above and below the horizon of interest.
t1_window=0;
t2_window=20;

% Array of strings containing the names of the seismic attribute files.
Files_AT={'seismic_amplitude_crop_trace1-10000.segy',...
    'instantaneous_amplitude_crop_trace1-10000.segy'};

% Array of strings containing the titles of the input seismic attributes.
% Titulos_AT={'Amp. Sísmica', 'Amp. Instantánea', 'Frec. Instantánea'};
Titulos_AT={'Amp. Sísmica', 'Amp. Instantánea', 'Frec. Instantánea'};

% Definition of the weights for each attributes. It must be the same
% quantity as the input attributes and their sum must be equal to 1.
% These are used to control the contribution of each seismic attribute
% in the calculation of Fractal Dimension.
alpha=[0.5 0.5];

% Opening of the time surface of interest. The format must be X - Y - Time.
% Otherwise, it would generate an error.
% X-Y are coordinates.
assy = fopen('EM_MER-U2.prn');

% File to save the calculation of Fractal Dimension.
File_FD={'EM_MER-U2_Fractal_Dimension.prn'};

% Number of traces to perform Quality Control (QC )of the calculations.
n_tracesQC = 5;

%% PRE-ELIMINARY CALCULATIONS %%
% ############################ %

% ------------------------------ SURFACE ------------------------------- %
% ---------------------------------------------------------------------- %

% The file format of the file must be X-Y-Time. Otherwise it would
% generate an error. 

switch assy
    case -1
        warndlg('The surface file is not in the current directory',...
            'Error!');
        break
    otherwise
        try
            reg = fscanf(assy,'%10f');
            SURFACE=reshape(reg,3,[])';
        catch id
            errordlg('The file does not have the format X-Y-Time',...
                'Error!'); 
            break
        end 
end

% Number of points contained within the surface.
n_surface=size(SURFACE,1);

% ---------------------------- SESMIC DATA ----------------------------- %
% ---------------------------------------------------------------------- %

% Number of attributes to use in the algorithm
n_attributes = size(Files_AT,2);

% If the sum of the weights is different from one, or they differ from the
% number of input attributes, code execution stops.
if (sum(alpha) < 0.98)
    error('La sumatoria de los pesos suministrados no son iguales a 1'); 
elseif (n_attributes ~= size(alpha,2))
    error('El número de pesos no es igual a los atributos de entrada'); 
end

% The following data must be the same for each seismic volume.
% The sampling interval in a SegyHeader is in microseconds.
% dt_sample: Sampling interval.
% n_samples: Number of samples.
% n_traces: Number of traces.

if n_attributes >= 2
    [dt_sample,n_samples,n_traces]=Check_segy_data(Files_AT);
else
    % If there was selected only one seismic volume, the verification 
    % is not done
    [Segy_H] = ReadSegyHeader(Files_AT{1});
    dt_sample= Segy_H.dt/1000;
    n_samples= Segy_H.ns;
    n_traces = Segy_H.ntraces;
end

% X-Y coordinates of each seismic volume. These are divided by 100 to take
% them to metric scale. The coordinates are going to be the
% same in all seismic volumes, because it has already been verified that 
% the seismic volumes have the same parameters.

COORD = zeros(n_attributes*2,n_traces);
for i=1:n_attributes
    COORD((2*i-1),:) = ReadSegyTraceHeaderValue(Files_AT{i},'key','cdpX');
    COORD(2*i,:) = ReadSegyTraceHeaderValue(Files_AT{i},'key','cdpY');
    disp(['Volumen Sísmico:  ',  Files_AT{i}])
end
COORD = reshape(COORD,n_attributes*2,[])';
COORD = round(COORD./100);

% ----------------------------- TIME COLUMN ---------------------------- %
% ---------------------------------------------------------------------- %

% Time column. Time is assumed to start at zero.

t_inicio=0;
t = zeros(n_samples,1);
for i=1:n_samples
    t(i) = t_inicio + dt_sample*(i-1);
end

% k1 and k2. Number of samples to extract seismic data for Fractal 
% Dimension analysis. It depends on the time window previously chosen.

k1 = ceil(t1_window/dt_sample);
k2 = ceil(t2_window/dt_sample);

% --------------------- RESULTS AND QUALITY CONTROL --------------------- %
% ----------------------------------------------------------------------- %

% This matrix will store the results of Fractal Dimension calculations.
FD = zeros(n_surface,3);

% This vector will store the values ??of R2 in order to make a histogram.
R2_ALL = zeros(n_surface,1);

% Definition of 5 five random numbers between 1 and the number of traces in
% seismic volumes. This will represent 5 random traces, in order
% to perform QA of Fractal Dimension calculations.

% QC = rand(1,n_tracesQC);
% QC = round(QC .* n_traces);

% Only for test of attribute volumes with 10,000.
QC = [6557 357 600 9340 6787];

% 1.This matrix will keep the seismic attributes traces.
Tr_Analysis_QC = zeros(k1+k2+1,n_attributes,n_tracesQC);

% 2. This matrix will save the number of traces, slope, Fractal Dimension
% and the coefficient of determination R2
QC_nPD_R2 = zeros(n_tracesQC,4);

% 3. This matrix will retain traces coordinates.
QC_COORD = zeros(n_tracesQC,2);

% 4. This matrix will have time interval analysis.
QC_TIME =  zeros(k1+k2+1,n_tracesQC);

% 5. This matrix will have divisor (r) and total length (L). The length of
% both vectors is unknown, so pre-allocation is made with a single column.
QC_rL = zeros(2,1,n_tracesQC);

%% START %%
% ######## %

wb=waitbar(0,'Analizando Superficie de Interés...');
n=1;

% N: numbers of the traces where the coordinates of the surface of ??interest
% is equal to the coordinate of the read trace.
N = zeros(1,n_surface) ;

% condition_QC: Quality Control Condition. Sets the selection of
% the random trace to save its corresponding data.
condition_QC = false;

% tr: index that indicates the trace number to QC.
tr=1;

for j=1:n_surface

    waitbar(j/n_surface,wb)
    
    % Find the position where the surface coordinate is equal to seismic
    % trace coordinate. Comparison is only with the first trace, because
    % it has been verified that they are all the same.
    
    % n= index that locates the trace number read.
    
    condition_time_suface = false;
    while condition_time_suface == false && n < n_traces
        condition_time_suface = isequal(SURFACE(j,1:2),COORD(n,1:2));
        n = n + 1;
    end
    
    % When n>n_traces the search has reached the end of seismic data. This
    % can also happen if the coordinates of the area of ??interest are not
    % within the seismic data. In either case, program execution must stop.
    
    if n >= n_traces
        errordlg('Finalizó la búsqueda dentro de los datos sísmicos',...
                'Error!');
        break
    end
    
    n = n - 1;
    
    N(j) = n;
    
    disp(['Traza Sísmica nro:  ',  num2str(n)])
    
    % Reading of the "n-th" trace found in seismic data.
    Seismic_Trace = zeros(n_samples,n_attributes);
    for i=1:n_attributes 
        [Data,Segy_TH]=ReadSegy(Files_AT{i},'traces',n);
        Seismic_Trace(:,i) = Data;
    end
    
    % Check if any of the traces are empty. If this case, it will not
    % calculate Fractal Dimension.    
    condition_tr_empty = any(all(Seismic_Trace==0)==1);

    % Both conditions must be met to continue calculations. Otherwise, 
    % it will continue to the next trace.

    if (condition_tr_empty == false) && (condition_time_suface == true)
        % Save the time found
        % ------------------------------------------------------------
        
        % Find the position of the time from the surface in the
        % time column of seismic data, to take data from
        % amplitudes of the corresponding attributes.
        t_surface = floor(SURFACE(j,3));

        % Check if it is a multiple of 4. Otherwise, take the multiple
        % immediate higher
        C = mod(t_surface,dt_sample);
        while C~=0
            t_surface = t_surface - 1;
            C=mod(t_surface,dt_sample);
        end

        % Find when surface time is equal to time vector.
        k=find(t_surface == t);
        
        % According to the position located, the values ??of attributes
        % are extracted using the time window previously defined.
        Tr_Analysis = zeros(k1+k2+1,n_attributes);

        for i = 1:n_attributes;
            Tr_Analysis(:,i) = Seismic_Trace(k-k1:k+k2,i);
        end
        
        % Input data normalization
        Tr_Analysis_norm = normalization(Tr_Analysis , Tr_Analysis);
    
        % Multiplication of the normalized traces by previously defined 
        % weights.
        for i=1:n_attributes
            Tr_Analysis_norm(:,i)= Tr_Analysis_norm(:,i) .* alpha(i);
        end

        % The time column in the window of interest does not use for
        % Fractal Dimension calculations. It will be use only for QC.
        % Because the input data was normalized, the time (in the last 
        % column of the array), must also range from 0 to 1.
        
        Tr_Analysis_norm(:,end+1) = zeros(k1+k2+1,1);

        div = 1/(size(Tr_Analysis_norm,1)-1);
        h=0;
        for i=1:size(Tr_Analysis_norm,1)
            Tr_Analysis_norm(i,end)=h*div;
            h=h+1;
        end
        
        % Fractal Dimension calculation
        [r,L,P,R2]=divplot2_EM(Tr_Analysis_norm);

        % Fractal Dimension of the curve in space N-Dimension.
        D = 1 + abs(P);
        
        % This matrix will keep the results of each FRACTAL DIMENSION
        % calculations. The order is X-Y-D.
        FD(j,:) = [SURFACE(j,1) SURFACE(j,2) D]; 

        % Save R2 values, in order to perform a histogram with the results.
        R2_ALL(j)=R2;
        
        % For QC, from 5 random traces, the following data will be saved:
        % Input traces, trace number, D, coordinates, r-L
        try
            condition_QC = any(QC >= N(j-1) & QC <= N(j));
        catch
            continue
        end
        
        if condition_QC;
            % Input Trace: Each 3D from the matrix Tr_Analysis_norm_QC has
            % data traces.
            Tr_Analysis_QC(:,:,tr) = Tr_Analysis;
            
            % Each row will correspond to the data of each trace:
            % Trace Number - P - D - R2.
            QC_nPD_R2(tr,1) = n; 
            QC_nPD_R2(tr,2) = P; 
            QC_nPD_R2(tr,3) = D; 
            QC_nPD_R2(tr,4) = R2;
            
            % Coordinates X/Y - r/L
            QC_COORD(tr,:) = SURFACE(j,1:2); 

            % Time column in window analysis (FD_QC_TIME), according to
            % time found in the surface of interest.            
            t_ini_analisis= t_surface - (k1 * dt_sample);            
            for i = 1 : (k1+k2+1)
                QC_TIME(i,tr) = t_ini_analisis + dt_sample*(i-1);
            end
            
            % Divider size "r" and total length "L". The 3rd dimension
            % of the hyper matrix corresponds to each trace.
            QC_rL(1,1:size(r,2),tr) = r; 
            QC_rL(2,1:size(r,2),tr) = L; 

            tr = tr + 1;
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

FD_erased = reshape(FD,n_surface,[])';
empty_coord = find(FD_erased==0);
FD_erased(empty_coord) = [];
FD_erased = reshape(FD_erased,3,[])';

% Fractal Dimension surface file.
dlmwrite(File_FD{1},FD_erased,'delimiter',' ','precision','%15.5f')

% Graphics.
% It will be both plots in the same window. 
for i=1:size(QC_nPD_R2,1)
    graph_QC(Tr_Analysis_QC , QC_nPD_R2 , QC_rL , QC_COORD , QC_TIME,...
        Titulos_AT,i);
end

% R2 Histogram. All possible empty cells are erased.
figure('Name','Histograma de R^2 en cálculo de Dimensión Fractal')
R2_ALL(R2_ALL==0)=[];
hist(R2_ALL,(0.1:0.1:0.9));
ylabel('Frecuencia'), xlabel('R^2');
title('Histograma de Frecuencia Absoluta. Coeficiente de Determinación.');
