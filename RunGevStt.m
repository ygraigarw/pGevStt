%% Simple code to fit GEV model, all the parameters of which are assumed to BE CONSTANT in time
% GEV(xi, sgm, mu) with xi=xi_0, sigma=sigma_0, mu=mu_0, 
% Will run as is for toy data to check 
% Need to input data as structure (see occurrences of USER INPUT) below
%
% This is a simplified version of pGevNonStt.

%% Set up
clc; clear; clf; pLtx;
VrbNms={'$\xi$';'$\sigma$';'$\mu$'};

%% Simulate a sample of data
if 1; %for testing

    % Time variable (assumed to be defined on 0,1)
    X.nT=10000;
    X.Tim=linspace(0,1,X.nT)';
    
    % True parameters P0=[xi0;sgm0;mu0;]
    X.Prm0=[-0.3;1;5];       
    
    % True parameter estimates in time
    X.XSM0=[ones(X.nT,1)*X.Prm0(1) ones(X.nT,1)*X.Prm0(2) ones(X.nT,1)*X.Prm0(3)];
    
    % Generate data from GEV
    X.Dat=gevrnd(X.XSM0(:,1),X.XSM0(:,2),X.XSM0(:,3));
    
    X, % See the structure
    
end;

%% ***USER INPUT*** Read in your data here
if 0; %set to zero if you want to use simulated data from above
     %X.nT ;  %  1 x 1 number of time points
     %X.Tim ; % nT x 1 times on [0,1]
     %X.Dat ; % nT x 1 data 
     load userInput.mat;
end;

%% Find starting solution by GEV fitting to independent blocks of data
if 1; 
    
    Y.nT=X.nT;
    Y.Tim=X.Tim;
    Y.Dat=X.Dat;
    Y.nT=size(X.Dat,1);
    Y.nB=10;
    Y.Blc=pMakCV(Y.nT,Y.nB,Y.Tim);
    
    clf;
    subplot(2,1,1); plot(Y.Tim, Y.Dat,'k.');
    xlabel 'Time';
    ylabel 'Value';
    pAxsLmt;pDflBig;
    
    tRgr=nan(Y.nB,4);
    for iB=1:Y.nB;
        tPrm=gevfit(Y.Dat(Y.Blc==iB));
        tTim=mean(Y.Tim(Y.Blc==iB));
        tRgr(iB,:)=[tTim tPrm];
    end;
    [jnk,tOrd]=sort(tRgr(:,1));
    Y.Rgr=tRgr(tOrd,:);
    
    Y.XSMStart=mean(Y.Rgr(:,2:4))';
    for j=1:3;
        subplot(2,3,3+j); plot(Y.Rgr(:,1),Y.Rgr(:,j+1),'b.-');
        subplot(2,3,3+j); hold on; plot(Y.Rgr(:,1),ones(Y.nB,1)*Y.XSMStart(j),'k.-');
        ylabel(VrbNms{j},'interpreter','latex'); 
        if j==1; xlabel 'Time'; end;
        pAxsLmt;pDflBig;
    end;
    
end;

%% ***USER INPUT*** Find better solution using MCMC
if 1;
    
    C.nI=10000;       % Number of MCMC iterations - 1e4 minimum when used in anger
    C.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
    
    C.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
    C.AdpItr=1000;    % Number of warm up iterations - don't change
    C.AdpBet=0.05;    % Adaptive MC - don't change
    
    clf;
    C=GevSttMCMC(Y,C);   % Run MCMC algorithm
    
end;

%% Plot MCMC results
if 1;
    
    clf;
    VrbNms={'$\xi$';'$\sigma$';'$\mu$'};
    load MCMC;
    for j=1:3;
        subplot(2,3,j);
        pHst(C.Prm(C.nI-C.n2Plt+1:end,j));
        title(VrbNms{j},'interpreter','latex');
        pAxsLmt; pDflBig;
        fprintf(1,'Median %s = %g\n',VrbNms{j},quantile(C.Prm(C.nI-C.n2Plt+1:end,j),0.5));
    end;
    subplot(2,3,4); plot(C.Nll,'k-'); title 'NLL'; pAxsLmt; pDflBig;
    subplot(2,3,5); plot(C.AccRat); title 'Acceptance rates'; pAxsLmt; pDflBig;
    pGI('GevStt-McmcTracePlot',2);
    
end;


%% Return value plot
if 1;

    %% Calculate return values
    C.RV.RtrPrd=100; % Return period of interest
    C.RV.nRls=1000;  % Number of realisations to use (1000 is good)
    
    % Parameter estimates
    t=randi(C.nI-C.n2Plt,C.RV.nRls,1)+C.n2Plt;
    tXi=C.Prm(t,1);
    tSgm=C.Prm(t,2);
    tMu=C.Prm(t,3);
    C.RV.Est(:,1)=(tSgm./tXi).*( (-log(1-1/C.RV.RtrPrd)).^(-tXi) - 1 ) + tMu; % Return value
    C.Dgn.Prm=[mean(tXi);mean(tSgm);mean(tMu)];
            
    % Summary statistics of differences
    C.RV.Cdf=[sort(C.RV.Est(:,1))];
    C.RV.Qnt=[quantile(C.RV.Est(:,1),[0.025 0.5 0.975])]; % Quantiles
        
    %% Figure
    clf; 
    C.RV.PrbVls=((1:C.RV.nRls)'-0.5)/C.RV.nRls;
    subplot(2,2,1); hold on;
    plot(C.RV.Cdf(:,1),C.RV.PrbVls,'k');
    if isfield(X,'Prm0')==1; % True values are known
        tXi=X.Prm0(1);
        tSgm=X.Prm0(2);
        tMu=X.Prm0(3);
        tRVTru=(tSgm./tXi).*( (-log(1-1/C.RV.RtrPrd)).^(-tXi) - 1 ) + tMu;
        plot(tRVTru*ones(2,1),[0 1]','k--','linewidth',2);
    end;
    pAxsLmt; pDflHug;
    title 'Distribution of RV [True - - -]'; 
    xlabel 'Annual maximum value';
    ylabel('$F_{{Annual Maximum}}$','interpreter','latex')
    pAxsLmt; pDflHug;
    %
    subplot(2,2,2); hold on;
    plot(C.RV.Cdf(:,1),log10(1-C.RV.PrbVls),'k','linewidth',2);
    pAxsLmt; pDflHug;
    title 'Distribution of RV [True - - -]'; 
    xlabel 'Annual maximum value';
    ylabel('$\log_{10}(1-F_{{Annual Maximum}})$')
    pAxsLmt; pDflHug;
    %
    x=(min(X.Dat):0.01:max(X.Dat))';
    z=gevcdf(x,C.Dgn.Prm(1),C.Dgn.Prm(2),C.Dgn.Prm(3));
    subplot(2,1,2); hold on;
    plot(x, log10(1-z),'color','r');
    plot(sort(X.Dat),log10(1-(((1:X.nT)'-0.5)/X.nT)),'ko','linewidth',2);
    title 'Diagnostic for tail (k=data, r=fit)';
    xlabel 'Annual maximum value';
    ylabel('$\log_{10}(1-F_{{Annual Maximum}})$')
    pAxsLmt; pDflHug;
    
    pDatStm; pGI('GevStt-100YearReturnValue',2);
    
    %% Summary statistics to screen
    fprintf(1,'SUMMARY\n');
    fprintf(1,'Quantiles 2.5pc 50pc 97.5pc: %g %g %g\n',C.RV.Qnt(1:3));

end;






























