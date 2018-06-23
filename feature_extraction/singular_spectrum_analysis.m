function out = singular_spectrum_analysis(Y,L,fs,plotFlag)
%%      Calculates the singular spectrum analysis (ssa) for a 1D signal Y with lag L
%       Modifications: 
%       Jonathan Foldager 
%       Stanford University
%       2018
%%      Inputs
%       Y:          Single channel signal
%       L:          Lag in samples / window length, must be 2 < L < T. Theoretically, 
%                   L < T/2.
%       fs:         Sampling frequency
%       plotFlag:   Should plot boolean
%%      Outputs
%       out.C:                  Covariance matrix of trajectory matrix
%       out.PC:                 Principal components of signal;
%       out.RC:                 Reconstruction elements of signal (sum(out.RC,2)~Y);
%       out.Eig.Vec = eigVec;   Eigenvectors of covariance matrix out.C
%       out.Eig.Val = eigVal;   Eigenvalues of covariance matrix  out.C
%%  Algorithm
    N = length(Y);                      % Number of samples
    t = linspace(0,N,N)*1/fs;           % Time vector (for plotting)
    M = N - L + 1;                      % Number of lagged-vectors that X should consists of.
    X = zeros(L,M);                     % Trajectory matrix
    for m=1:M
      X(:,m) = Y((1:N-M+1)+m-1);        % Filling X with time delayed Y's
    end    
    C               = X'*X / (N-M+1);   % Unbiased covariance matrix of X
    [eigVec,eigVal] = eig(C);           % Temporal empirical orthogonal functions
    [eigVal,ind]    = sort(diag(eigVal),'descend'); 
    eigVec          = eigVec(:,ind);                    
%%  Principal components
    PC = X*eigVec;
%%  Reconstruction
    RC = zeros(N,M);
    for m=1:M
        buf=PC(:,m)*eigVec(:,m)'; % invert projection
        buf=buf(end:-1:1,:);
        for n=1:N % anti-diagonal averaging
            RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
        end
    end
%%  Output
    out.C = C;
    out.PC = PC;
    out.RC = RC;
    out.Eig.Vec = eigVec;
    out.Eig.Val = eigVal;
%%  Plot
    if plotFlag
        figure;
        for m=1:4
          subplot(4,1,m);
          plot(t(1:N-M+1),PC(:,m),'k-');
          ylabel(sprintf('PC %d',m));
        end

        figure;
        for m=1:4
          subplot(4,1,m);
          plot(t,RC(:,m),'r-');
          ylabel(sprintf('RC %d',m));
        end    

        figure;
        subplot(2,1,1)
        plot(t,Y,'b-',t,sum(RC(:,:),2),'r-');
        legend('Original','Complete reconstruction');

        subplot(2,1,2)
        plot(t,Y,'b','LineWidth',2);
        plot(t,Y,'b-',t,sum(RC(:,1:2),2),'r-');
        legend('Original','Reconstruction with RCs 1-2');
    end
end