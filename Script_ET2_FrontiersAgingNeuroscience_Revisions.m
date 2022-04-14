%% This script performs graph theoretical analyses on essential tremor
% patients before and after surgical intervention, and also compares these
% patients to healthy controls matched for age. Several types of analysis
% are implemented:
%
% 1. Analysis of graph theoretical metrics (degree, clustering coefficient
% and eigenvector centrality) from structural covariance graphs generated
% from cortical thickness (CT), surface area (SA) and mean curvature (MC).
% Healthy controls (HCs) are compared to ET patients before thalamotomy
% (ETpre), and post-intervention ET patients (ETpost) to the 
% pre-thalamotomy state.
%
% 2. Analysis of cross-property (CT/SA, SA/MC and CT/MC) relationships, by
% the quantiication of in-degree and out-degree from directional graphs.
% The same contrasts as in 1. are examined.
%
% 3. Supplementary assessment in which mixed modelling is leveraged to
% account for within-subject variance in the data.



%% 0A. Loading of the data

addpath(genpath('Data'));
addpath(genpath('../Utilities'));

% Load the data from all 3 groups
load ET_Clinical
load ET_Freesurfer
load ET_Freesurfer_SC
load Age_HC
load Gender_HC
load ET_Covariates
load CodeBook

% Decides the type of metric to probe
Metric = 'Thickness';
Metric2 = 'Area';
Metric3 = 'MeanCurvature';

% The metrics are loaded, with also subcortical and cerebellar volume
% concatenated to the grey matter regional data

% Cortical thickness
X_HC = [(eval(['HC_',Metric]))';HC_Volume_SC'];
X_BASE = [(eval(['BASE_',Metric]))';BASE_Volume_SC'];
X_YEAR = [(eval(['YEAR_',Metric]))';YEAR_Volume_SC'];

% Surface area
Y_HC = [(eval(['HC_',Metric2]))';HC_Volume_SC'];
Y_BASE = [(eval(['BASE_',Metric2]))';BASE_Volume_SC'];
Y_YEAR = [(eval(['YEAR_',Metric2]))';YEAR_Volume_SC'];

% Surface area
Z_HC = [(eval(['HC_',Metric3]))';HC_Volume_SC'];
Z_BASE = [(eval(['BASE_',Metric3]))';BASE_Volume_SC'];
Z_YEAR = [(eval(['YEAR_',Metric3]))';YEAR_Volume_SC'];

% Number of subjects in total
n_HC = size(X_HC,2);
n_BASE = size(X_BASE,2);
n_YEAR = size(X_YEAR,2);

% Number of regions in the parcellation at play
n_regions = size(X_HC,1);
n_regions_cortical = 68;
n_regions_subcortical = 19;



%% 0B. Definition of analytical parameters and initialization of variables

%%%%%% PARAMETERS

% Will we plot the information or not?
is_plot = 1;

% Will we compute a non-parametric null distribution or not?
is_null = 1;

% Type of graph to generate: choose between Positive (only retains positive
% SCM edges), Negative (only retains and sign flips the negative edges), or
% Both (keeps both edge types with their signs)
Type = 'Positive';

% Number of edges in total with our number of regions
n_edges = n_regions*(n_regions-1)/2;

% Covariates to regress out for "traditional" analyses: note that we demean
% the data for each group in addition to the rest
load('ET_Covariates','BASE_TGV','YEAR_TGV','HC_TGV');
load('ET_Clinical','Age','Gender');
Cov_BASE = [ones(n_BASE,1),Age,Gender,BASE_TGV];
Cov_YEAR = [ones(n_YEAR,1),Age+1,Gender,YEAR_TGV];
Cov_HC = [ones(n_HC,1),Age_HC,Gender_HC,HC_TGV];

% Significance level at which we want to test (uncorrected) in the
% analyses, and type of tail for the tests ('both' = two-tailed)
Alpha = 0.05;
tail = 'both';

% Number of null realizations to go for
n_null_regions = 8000;

% Range of Rho values to probe
Rho_range = 20:10:60;

%%%%%%%% COLORMAPS

% Colormap to use for some representations (red-blue)
CM_RB = cbrewer('div','RdBu',1001);
CM_RB(CM_RB<0) = 0;

CM_RYG = cbrewer('div','RdYlGn',1001);
CM_RYG(CM_RYG < 0) = 0;

CM_Paired = cbrewer('qual','Paired',15);
CM_Paired(CM_Paired < 0) = 0;

% Mapping to lobes (including subcortical areas and cerebellum)
LobesMapping = [[5,9,1,7,5,5,3,5,9,7,1,7,1,5,5,1,1,1,1,7,3,9,1,3,9,1,1,3,5,3,1,5,5,11],1+[5,9,1,7,5,5,3,5,9,7,1,7,1,5,5,1,1,1,1,7,3,9,1,3,9,1,1,3,5,3,1,5,5,11]]';
LobesMapping(69:77) = 13;
LobesMapping(78:86) = 14;
LobesMapping(87) = 15;

CM_YG = cbrewer('seq','YlGn',1001);
CM_YG(CM_YG < 0) = 0;

CM_Subjects = colorcube(n_BASE);
CM_Subjects2 = flipud(cbrewer('div','RdYlGn',n_BASE));
CM_Subjects2(CM_Subjects2 < 0) = 0;

CM_SCM = cbrewer('seq','Greys',1000);
tmp = cbrewer('div','RdYlBu',1000);



%% 0C. Regressing out covariates from the data

% Cortical thickness
for r = 1:n_regions
    [~,X_BASE_res(r,:)] = y_regress_ss(X_BASE(r,:)',Cov_BASE);
    [~,X_YEAR_res(r,:)] = y_regress_ss(X_YEAR(r,:)',Cov_YEAR);
    [~,X_HC_res(r,:)] = y_regress_ss(X_HC(r,:)',Cov_HC);
end

% Surface area
for r = 1:n_regions
    [~,Y_BASE_res(r,:)] = y_regress_ss(Y_BASE(r,:)',Cov_BASE);
    [~,Y_YEAR_res(r,:)] = y_regress_ss(Y_YEAR(r,:)',Cov_YEAR);
    [~,Y_HC_res(r,:)] = y_regress_ss(Y_HC(r,:)',Cov_HC);
end


% Mean curvature
for r = 1:n_regions
    [~,Z_BASE_res(r,:)] = y_regress_ss(Z_BASE(r,:)',Cov_BASE);
    [~,Z_YEAR_res(r,:)] = y_regress_ss(Z_YEAR(r,:)',Cov_YEAR);
    [~,Z_HC_res(r,:)] = y_regress_ss(Z_HC(r,:)',Cov_HC);
end



%% 1. Regional investigations in terms of graph properties

% Number of times we run the null process
idx_rho = 1;

% We run the same across all values for Rho, computing the actual group
% differences for the graph theoretical metrics of interest. Summing into
% area under the curve is done below (Assess_AUC function)
for Rho = Rho_range

    Rho
    
    Delta_Deg_YB(:,idx_rho) = sum(Generate_SCM(Z_YEAR_res, Rho, Type,'Pearson'),2)-sum(Generate_SCM(Z_BASE_res, Rho, Type,'Pearson'),2);
    Delta_CC_YB(:,idx_rho) = clustering_coef_wu(Generate_SCM(Z_YEAR_res, Rho, Type,'Pearson'))-clustering_coef_wu(Generate_SCM(Z_BASE_res, Rho, Type,'Pearson'));
    Delta_EC_YB(:,idx_rho) = eigenvector_centrality_und(Generate_SCM(Z_YEAR_res, Rho, Type,'Pearson'))-eigenvector_centrality_und(Generate_SCM(Z_BASE_res, Rho, Type,'Pearson'));
    
    Delta_Deg_HB(:,idx_rho) = sum(Generate_SCM(Z_HC_res, Rho, Type,'Pearson'),2)-sum(Generate_SCM(Z_BASE_res, Rho, Type,'Pearson'),2);
    Delta_CC_HB(:,idx_rho) = clustering_coef_wu(Generate_SCM(Z_HC_res, Rho, Type,'Pearson'))-clustering_coef_wu(Generate_SCM(Z_BASE_res, Rho, Type,'Pearson'));
    Delta_EC_HB(:,idx_rho) = eigenvector_centrality_und(Generate_SCM(Z_HC_res, Rho, Type,'Pearson'))-eigenvector_centrality_und(Generate_SCM(Z_BASE_res, Rho, Type,'Pearson'));
    
    idx_rho = idx_rho + 1;
end

% For loop for null data generation
for n = 1:n_null_regions
    
    n

    % Shuffled data for null ETpre and ETpost groups
    tmp_X = [Z_BASE_res';Z_YEAR_res'];

    idx = randperm(n_BASE+n_YEAR);

    tmp_X = tmp_X(idx,:);
    tmp_BASE = tmp_X(1:n_BASE,:);
    tmp_YEAR = tmp_X(n_BASE+1:end,:);
    
    % Same for the HC - ETpre case
    tmp_X = [Z_BASE_res';Z_HC_res'];

    idx = randperm(n_BASE+n_HC);

    tmp_X = tmp_X(idx,:);
    tmp_BASE2 = tmp_X(1:n_BASE,:);
    tmp_HC = tmp_X(n_BASE+1:end,:);

    idx_rho = 1;
    
    for Rho = Rho_range
 
        Delta_Deg_YB_null(:,idx_rho,n) = sum(Generate_SCM(tmp_YEAR', Rho, Type,'Pearson'),2)-sum(Generate_SCM(tmp_BASE', Rho, Type,'Pearson'),2);
        Delta_CC_YB_null(:,idx_rho,n) = clustering_coef_wu(Generate_SCM(tmp_YEAR', Rho, Type,'Pearson'))-clustering_coef_wu(Generate_SCM(tmp_BASE', Rho, Type,'Pearson'));
        Delta_BTW_YB_null(:,idx_rho,n) = betweenness_wei(Generate_SCM(tmp_YEAR', Rho, Type,'Pearson'))-betweenness_wei(Generate_SCM(tmp_BASE', Rho, Type,'Pearson'));
        Delta_EC_YB_null(:,idx_rho,n) = eigenvector_centrality_und(Generate_SCM(tmp_YEAR', Rho, Type,'Pearson'))-eigenvector_centrality_und(Generate_SCM(tmp_BASE', Rho, Type,'Pearson'));
        Delta_Eff_YB_null(idx_rho,n) = efficiency_wei(Generate_SCM(tmp_YEAR', Rho, Type,'Pearson'))-efficiency_wei(Generate_SCM(tmp_BASE', Rho, Type,'Pearson'));

        Delta_Deg_HB_null(:,idx_rho,n) = sum(Generate_SCM(tmp_HC', Rho, Type,'Pearson'),2)-sum(Generate_SCM(tmp_BASE2', Rho, Type,'Pearson'),2);
        Delta_CC_HB_null(:,idx_rho,n) = clustering_coef_wu(Generate_SCM(tmp_HC', Rho, Type,'Pearson'))-clustering_coef_wu(Generate_SCM(tmp_BASE2', Rho, Type,'Pearson'));
        Delta_BTW_HB_null(:,idx_rho,n) = betweenness_wei(Generate_SCM(tmp_HC', Rho, Type,'Pearson'))-betweenness_wei(Generate_SCM(tmp_BASE2', Rho, Type,'Pearson'));
        Delta_EC_HB_null(:,idx_rho,n) = eigenvector_centrality_und(Generate_SCM(tmp_HC', Rho, Type,'Pearson'))-eigenvector_centrality_und(Generate_SCM(tmp_BASE2', Rho, Type,'Pearson'));
        Delta_Eff_HB_null(idx_rho,n) = efficiency_wei(Generate_SCM(tmp_HC', Rho, Type,'Pearson'))-efficiency_wei(Generate_SCM(tmp_BASE2', Rho, Type,'Pearson'));
        
        idx_rho = idx_rho + 1;
    end
end

% Significance assessment for the ETpost - ETpre contrast; the returned
% p-values are FDR-corrected for the number of regions at hand. Bonferroni
% correction for the amount of parallel assessments is then addressed by
% setting a significance level at alpha=0.01 or alpha = 0.001 at which to
% assess FDR-corrected p-values
[p_BON_Deg_YB,p_FDR_Deg_YB] = Assess_AUC(Delta_Deg_YB,Delta_Deg_YB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));
[p_BON_CC_YB,p_FDR_CC_YB] = Assess_AUC(Delta_CC_YB,Delta_CC_YB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));
[p_BON_EC_YB,p_FDR_EC_YB] = Assess_AUC(Delta_EC_YB,Delta_EC_YB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));

% Same for the HC - ETpre case
[p_BON_Deg_HB,p_FDR_Deg_HB] = Assess_AUC(Delta_Deg_HB,Delta_Deg_HB_null,tail,Alpha,Rho_range(1:5),CM_Paired(LobesMapping,:));
[p_BON_CC_HB,p_FDR_CC_HB] = Assess_AUC(Delta_CC_HB,Delta_CC_HB_null,tail,Alpha,Rho_range(1:5),CM_Paired(LobesMapping,:));
[p_BON_EC_HB,p_FDR_EC_HB] = Assess_AUC(Delta_EC_HB,Delta_EC_HB_null,tail,Alpha,Rho_range(1:5),CM_Paired(LobesMapping,:));



%% 2. Cross-property graph analysis

idx_rho = 1;

% Actual differences are computed
for Rho = Rho_range
    
    % Correlation across modalities
    tmp_BASE = abs(Threshold_SCM(corr(Y_BASE_res(1:n_regions_cortical,:)',Z_BASE_res(1:n_regions_cortical,:)'),Rho,'Both'));
    tmp_HC = abs(Threshold_SCM(corr(Y_HC_res(1:n_regions_cortical,:)',Z_HC_res(1:n_regions_cortical,:)'),Rho,'Both'));
    tmp_YEAR = abs(Threshold_SCM(corr(Y_YEAR_res(1:n_regions_cortical,:)',Z_YEAR_res(1:n_regions_cortical,:)'),Rho,'Both'));
    
    % Computation of graph metrics
    [kin_BASE,kout_BASE] = degrees_dir(tmp_BASE);
    [kin_HC,kout_HC] = degrees_dir(tmp_HC);
    [kin_YEAR,kout_YEAR] = degrees_dir(tmp_YEAR);

    % Differences of interest
    kin_YB(:,idx_rho) = kin_YEAR - kin_BASE;
    kin_HB(:,idx_rho) = kin_HC - kin_BASE;

    kout_YB(:,idx_rho) = kout_YEAR - kout_BASE;
    kout_HB(:,idx_rho) = kout_HC - kout_BASE;
    
    idx_rho = idx_rho + 1;
end

% Null data scheme
if is_null
    for n = 1:n_null_regions

        n
        
        % Null data groups for the ETpost - ETpre contrast
        tmp_X = [Y_BASE_res(1:n_regions_cortical,:)';Y_YEAR_res(1:n_regions_cortical,:)'];
        tmp_Y = [Z_BASE_res(1:n_regions_cortical,:)';Z_YEAR_res(1:n_regions_cortical,:)'];

        idx = randperm(n_BASE+n_YEAR);

        tmp_X = tmp_X(idx,:);
        tmp_Y = tmp_Y(idx,:);
        
        tmp_BASE_X = tmp_X(1:n_BASE,:);
        tmp_BASE_Y = tmp_Y(1:n_BASE,:);
        
        tmp_YEAR_X = tmp_X(n_BASE+1:end,:);
        tmp_YEAR_Y = tmp_Y(n_BASE+1:end,:);

        % Same for the HC - ETpre case
        tmp_X = [Y_BASE_res(1:n_regions_cortical,:)';Y_HC_res(1:n_regions_cortical,:)'];
        tmp_Y = [Z_BASE_res(1:n_regions_cortical,:)';Z_HC_res(1:n_regions_cortical,:)'];

        idx = randperm(n_BASE+n_HC);

        tmp_X = tmp_X(idx,:);
        tmp_Y = tmp_Y(idx,:);
        
        tmp_BASE2_X = tmp_X(1:n_BASE,:);
        tmp_BASE2_Y = tmp_Y(1:n_BASE,:);
        
        tmp_HC_X = tmp_X(n_BASE+1:end,:);
        tmp_HC_Y = tmp_Y(n_BASE+1:end,:);
        
        idx_rho = 1;
        
        for Rho = Rho_range
            
            % Null covariances
            tmp_null_YEAR = abs(Threshold_SCM(corr(tmp_YEAR_X,tmp_YEAR_Y),Rho,'Both'));
            tmp_null_BASE = abs(Threshold_SCM(corr(tmp_BASE_X,tmp_BASE_Y),Rho,'Both'));
            tmp_null_BASE2 = abs(Threshold_SCM(corr(tmp_BASE2_X,tmp_BASE2_Y),Rho,'Both'));
            tmp_null_HC = abs(Threshold_SCM(corr(tmp_HC_X,tmp_HC_Y),Rho,'Both'));

            % Null graph metrics
            [kin_YEAR,kout_YEAR] = degrees_dir(tmp_null_YEAR);
            [kin_BASE,kout_BASE] = degrees_dir(tmp_null_BASE);
            [kin_BASE2,kout_BASE2] = degrees_dir(tmp_null_BASE2);
            [kin_HC,kout_HC] = degrees_dir(tmp_null_HC);

            % Null differences
            kin_YB_null(:,idx_rho,n) = kin_YEAR-kin_BASE;
            kin_HB_null(:,idx_rho,n) = kin_HC-kin_BASE2;
            kout_YB_null(:,idx_rho,n) = kout_YEAR-kout_BASE;
            kout_HB_null(:,idx_rho,n) = kout_HC-kout_BASE2;
            
            idx_rho = idx_rho + 1;
        end
    end
end

% Significance assessent, done similarly as in section 1.
[p_BON_kin_YB,p_FDR_kin_YB] = Assess_AUC(kin_YB,kin_YB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));
[p_BON_kout_YB,p_FDR_kout_YB] = Assess_AUC(kout_YB,kout_YB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));
[p_BON_kin_HB,p_FDR_kin_HB] = Assess_AUC(kin_HB,kin_HB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));
[p_BON_kout_HB,p_FDR_kout_HB] = Assess_AUC(kout_HB,kout_HB_null,tail,Alpha,Rho_range,CM_Paired(LobesMapping,:));



%% 3. Supplementary analysis (mixed model approach)

% Creation of data for generating a table with everything, with ordering 
% HC -> pre -> post (needed for summoning the MATLAB mixel model utilities)
Age_MM = [Age;Age+1];
Gender_MM = [Gender;Gender];
Subject_MM = nominal([((1:n_BASE)');((1:n_YEAR)')]);
CT_MM = [X_BASE';X_YEAR'];
SA_MM = [Y_BASE';Y_YEAR'];
MC_MM = [Z_BASE';Z_YEAR'];
TGV_MM = [BASE_TGV;YEAR_TGV];
Effects_MM = nominal([ones(n_BASE,1);2*ones(n_YEAR,1)]);

% Initializing the output matrices for all cases at hand
MM_CT = zeros(n_regions,n_regions);
MM_SA = zeros(n_regions,n_regions);
MM_MC = zeros(n_regions,n_regions);

% Looping over individual region pairs
for r1 = 1:n_regions
    for r2 = 1:n_regions
        if r1 ~= r2
            
            %%% For cortical thickness
            
            % Creation of the input mprohometric features (normalized by
            % std) and of the output one (also normalized)
            Input_MM = [X_BASE(r1,:)'/std(X_BASE(r1,:)');X_YEAR(r1,:)'/std(X_YEAR(r1,:)')];
            Output_MM = [X_BASE(r2,:)'/std(X_BASE(r2,:)');X_YEAR(r2,:)'/std(X_YEAR(r2,:)')];
            
            % Creation of the data table
            TBL = table(Input_MM,Output_MM,Effects_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM);

            % Mixed model
            lme = fitlme(TBL,'Output_MM ~ Effects_MM*Input_MM + Age_MM + Gender_MM + TGV_MM + (1|Subject_MM)',...
                'FitMethod','REML','DummyVarCoding','effects');
            
            % Sampling the structural covariance coefficient
            MM_CT(r1,r2) = lme.Coefficients.Estimate(7);
           
            
            
            %%% For surface area
            
            Input_MM = [Y_BASE(r1,:)'/std(Y_BASE(r1,:)');Y_YEAR(r1,:)'/std(Y_YEAR(r1,:)')];
            Output_MM = [Y_BASE(r2,:)'/std(Y_BASE(r2,:)');Y_YEAR(r2,:)'/std(Y_YEAR(r2,:)')];
            
            TBL = table(Input_MM,Output_MM,Effects_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM);

            lme = fitlme(TBL,'Output_MM ~ Effects_MM*Input_MM + Age_MM + Gender_MM + TGV_MM + (1|Subject_MM)',...
                'FitMethod','REML','DummyVarCoding','effects');
  
            MM_SA(r1,r2) = lme.Coefficients.Estimate(7);
           
            
            
            %%% For mean curvature
            Input_MM = [Z_BASE(r1,:)'/std(Z_BASE(r1,:)');Z_YEAR(r1,:)'/std(Z_YEAR(r1,:)')];
            Output_MM = [Z_BASE(r2,:)'/std(Z_BASE(r2,:)');Z_YEAR(r2,:)'/std(Z_YEAR(r2,:)')];
            
            TBL = table(Input_MM,Output_MM,Effects_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM);

            lme = fitlme(TBL,'Output_MM ~ Effects_MM*Input_MM + Age_MM + Gender_MM + TGV_MM + (1|Subject_MM)',...
                'FitMethod','REML','DummyVarCoding','effects');
            
            MM_MC(r1,r2) = lme.Coefficients.Estimate(7);
        end
    end
end