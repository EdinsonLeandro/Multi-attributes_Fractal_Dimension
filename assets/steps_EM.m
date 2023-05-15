function nSteps = steps_EM(X, divider)
% The original function: n = steps3(X, d)
% calculates the number of steps to walk along the set of points X
% using a divider length d. 
% Written by Gerry Middleton, McMaster Univ., June 1996
% EMail: middleto@mcmaster.ca
% calls nextpt3.m To increase speed at the risk of accuracy
% change this to nextpt.
%
% This is a modified version to calculate the steps number
% needed to approximate the trace X with N-Dimension with "divider".

% Number of samples in data.
[nSamples, ~] = size(X);

% Initialize number of steps and index
nSteps = 0;
index = 2;

% sampleRef = Measure de distance from this point to the next sample.
% The algotithm starts with the first point.
sampleRef = X(1,:);

while index <= nSamples
    % It will compare with the given divisor.

    % Euclidean distance to the next sample.
    diff = sampleRef - X(index,:);
    distanceNext = sqrt(sum(diff.^2,2));
    
    % Euclidean distance to previous sample.
    diff = sampleRef - X(index-1, :);
    distancePrevious = sqrt(sum(diff.^2,2));
    
    % di = distanceNext
    % dim = distancePrevious

    if distanceNext >divider
        if distancePrevious == 0
            ratio = distanceNext/divider;
            nSteps = nSteps+ratio;
            
            % Select next sample
            sampleRef = X(index,:);
            index = index + 1;
        else
            % distancePrevious ~= 0
            dk = distanceNext - divider;
            dj = divider - distancePrevious;
            
            if dk > dj
                % Select same sample
                sampleRef= X(index-1,:);
                nSteps= nSteps+1;
                
            elseif dk < dj
                % Select next sample
                sampleRef = X(index,:);
                nSteps = nSteps + 1;
                index = index + 1;
                
            else
                % dk == dj
                % Select same sample
                sampleRef = X(index-1,:);
                ratio = distancePrevious/divider;
                nSteps = nSteps+ratio;
            end % End if
        end % End if
    
    elseif distanceNext < divider
        if index == nSamples
            ratio = distanceNext/divider;
            nSteps = nSteps+ratio;
            break
        end % End if
        %         
        % For distances less than "divider", the fraction "ratio" is calculated
        %         
        index=index+1;
    
    else
        % distanceNext == divider
        sampleRef = X(index,:);
        nSteps = nSteps+1;
        index = index+1;
    end % End if
    
end
