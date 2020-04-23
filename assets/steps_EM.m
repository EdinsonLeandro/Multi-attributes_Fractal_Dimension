function n = steps_EM(X,d)
% The original function: n = steps3(X, d)
% calculates the number n of steps to walk along the set of points X
% using a divider length d. 
% Written by Gerry Middleton, McMaster Univ., June 1996
% EMail: middleto@mcmaster.ca
% calls nextpt3.m To increase speed at the risk of accuracy
% change this to nextpt.
%
% This is a modified version to calculate the steps number (length "d")
% needed to approximate the trace X with N-Dimension.

n=0;
h=2;

% First Point
Xn=X(1,:);
% Xn=[X(1,1) X(1,2)];

[r c]=size(X);

while h<=r
    % Xn = From this point the distance is measured to the next sample.
    % It will compare with the given divisor.
    % DI  /  DIM = square distance.

    % Calculation of the distance to next sample.
    DI=zeros(1,c);
    for j=1:c
        DI(:,j)=Xn(:,j)-X(h,j);
    end
    DIsum=sum(DI.^2,2);
    di=sqrt(DIsum);
    
    % Calculation of the distance to previous sample.
    DIM=zeros(1,c);
    for j=1:c
        DIM(:,j)=Xn(:,j)-X(h-1,j);
    end
    DIMsum=sum(DIM.^2,2);
    dim=sqrt(DIMsum);
    
%     di = sqrt((Xn(:,1)-X(h,1)).^2 + (Xn(:,2)-X(h,2)).^2);
%     dim= sqrt((Xn(:,1)-X(h-1,1)).^2 + (Xn(:,2)-X(h-1,2)).^2);

    if di>d
        if dim == 0
            nf=di/d;
            n=n+nf;
            Xn=X(h,:);
            h=h+1;
        end
        
        if dim ~= 0
            dk=di-d;
            dj=d-dim;
            if dk>dj
                Xn=X(h-1,:);
                n=n+1;
            else
                Xn=X(h,:);
                n=n+1;
                h=h+1;
            end
            
            if dk==dj
                nf=dim/d;
                n=n+nf;
                Xn=X(h-1,:);
            end
        end
    end
    
    if d>di
        if h==r
            nf=di/d;
            n=n+nf;
            break
        end
        %         
        % For distances less than "d", the fraction "nf" is calculated
        %         
        h=h+1;
    end
    
    if di==d
        Xn=X(h,:);
        n=n+1;
        h=h+1;
    end

% % Counts the steps number n (length "d") that approximate all trace.
    
end
