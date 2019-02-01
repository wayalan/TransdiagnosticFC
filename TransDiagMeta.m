% fMRIDataSet is a structure for fMRI data and covariate data at each site
% together with the site name and covariates names at each site
% fMRIDataSet{i}.normal, fMRIDataSet{i}.id_n, fMRIDataSet{i}.patient, fMRIDataSet{i}.id_p, fMRIDataSet{i}.site_name,...
%             fMRIDataSet{i}.behaviours(:,fMRIDataSet{i}.index_beh), fMRIDataSet{i}.covariate_n, fMRIDataSet{i}.covarate_p, ...
%             [fMRIDataSet{i}.covarate_p, fMRIDataSet{i}.behaviours(:,fMRIDataSet{i}.index_cov)]

%% This is the example for transdiagnostic analysis between three groups

%% 1) 
%% identify the those FCs with survived the Bonferonni correction by meta-analysis on many datasets
load fMRIDataSet.mat

%% Multiple Comparison Correction
mcc=2; % 1:Bonferroni, 2:FDR, 3:uncorrected
SigP=0.05; % p-value threshold

%% calcluate FC on each site (dataset), here we have three sets: DEP:HC/BIP:HC/SCZ:HC
rawFC = CalculateFCOnEachDataSite(fMRIDataSet);
% rawFC{1}.patient and rawFC{1}.normal have been vectorized
% number of FC's calculated
numFC = size(rawFC{1}.normal,2);
% number of ROI's considered 
numROI = size(fMRIDataSet{1}.normal{1},2);
for i = 1 : length(fMRIDataSet)
    % Total sample size at each site
    sub1 = sum(fMRIDataSet{i}.id_n);
    sub2 = sum(fMRIDataSet{i}.id_p);
    sampleSize(i) =  sub1 + sub2;
end

%% 2)
%% meta-analysis to identify those FC's survived the Bonferoni correction  (iniFC)
iniMetaResult = []; % Comparison between patient and HC
iniMetaResult_pt = []; % Comparison among patient groups only
for i = 1 : numFC
    % calculate beta and variance and p value at each site for each link
    % between control and patients
    [beta, p_beta, var_beta, ~, effSize] = CalculateBetaOnEachDataSite(fMRIDataSet,rawFC, i);
    % calculate beta and variance and p value at each site for each link
    % between patient groups
    [beta_pt, p_beta_pt, var_beta_pt, ~, effSize_pt] = CalBetaBTPatients(fMRIDataSet,rawFC, i);
    % meta combined effect size for each link
    outputMeta = MetaAnalysis(beta, var_beta, p_beta, sampleSize, 3); % we need to report the direction of the group difference for each link, so we use fixed effect model   
    iniMetaResult(i,:) = [beta, p_beta, outputMeta.EffectSizeMeta, outputMeta.pMeta];   
    iniMetaResult_pt(i,:)=[beta_pt, p_beta_pt, effSize_pt];
    outputMeta2 = MetaAnalysis(effSize, var_beta, p_beta, sampleSize, 3);
    iniMetaEffSize(i,:) = effSize;
    iniESMeta(i,:)=outputMeta2.EffectSizeMeta;
end

%% bonferonni correction to establish the initial feature set
Site=[];
n=ones(1,7);

if mcc==1 %% Bonferroni Correction
    mcP=SigP/numFC;
    mcP(1:4)=mcP;
elseif mcc==2 %% FDR Correction
    % FDR for the p array of meta-analysis
    mcP(4)=FDR(iniMetaResult(:,8),SigP);
    % FDR for the p values between each of the two groups
    mcP(1:3)=FDR(iniMetaResult(:,4:6),SigP);
elseif mcc==3 %% Uncorrected-P
    mcP=SigP;
    mcP(1:4)=mcP;
end

%% 3)
%% Indentify the FC that satisfied the criteria of transdiagnostic and illness-specific model (mcP(i): p-value of group pair i)
for i = 1: numFC
    
    % Transdiagnostic Result: differences found in each of the PT vs HC, as well as having same direction (iniMetaResult(i,8): meta-p)
    if iniMetaResult(i,4) < mcP(1) && iniMetaResult(i,5) < mcP(2) && iniMetaResult(i,6) < mcP(3) && iniMetaResult(i,8) < mcP(4) 
        Site.Trans(n(1,1),1) = i;
        n(1,1)=n(1,1)+1;

    % Illness-specific Result: differences found in the first set of the PT vs HC, no significant difference can be found in any other group set.
    elseif iniMetaResult(i,4) < mcP(1) && iniMetaResult(i,5) > mcP(2) && iniMetaResult(i,6) > mcP(3) && iniMetaResult(i,8) > mcP(4) 
        Site.one(n(1,2),1)=i;
        n(1,2)=n(1,2)+1;

    % Illness-specific Result: differences found in the second set of the PT vs HC, no significant difference can be found in any other group set.
    elseif iniMetaResult(i,4) > mcP(1) && iniMetaResult(i,5) < mcP(2) && iniMetaResult(i,6) > mcP(3) && iniMetaResult(i,8) > mcP(4) 
        Site.two(n(1,3),1)=i;
        n(1,3)=n(1,3)+1;

    % Illness-specific Result: differences found in the third set of the PT vs HC, no significant difference can be found in any other group set.
    elseif iniMetaResult(i,4) > mcP(1) && iniMetaResult(i,5) > mcP(2) && iniMetaResult(i,6) < mcP(3) && iniMetaResult(i,8) > mcP(4) 
        Site.three(n(1,4),1)=i;
        n(1,4)=n(1,4)+1;

    % Conjunction of any two diagnoses: PT1 & PT2
    elseif iniMetaResult(i,4) < mcP(1) && iniMetaResult(i,5) < mcP(2) && iniMetaResult(i,6) > mcP(3) 
        Site.onetwo(n(1,5),1)=i;
        n(1,5)=n(1,5)+1;
    % Conjunction of any two diagnoses: PT2 & PT3
    elseif iniMetaResult(i,4) < mcP(1) && iniMetaResult(i,5) > mcP(2) && iniMetaResult(i,6) < mcP(3) 
        Site.onethree(n(1,6),1)=i;
        n(1,6)=n(1,6)+1;
    % Conjunction of any two diagnoses: PT1 & PT3
    elseif iniMetaResult(i,4) > mcP(1) && iniMetaResult(i,5) < mcP(2) && iniMetaResult(i,6) < mcP(3)
        Site.twothree(n(1,7),1)=i;
        n(1,7)=n(1,7)+1;
    end
end    

% Re-identify the transdiagnostic and illness-specific dysconnectivity 
CheckBtPt=[]; %iniMetaResult_pt(:,4): DEP vs. BIP; iniMetaResult_pt(:,5): DEP vs. SCZ; iniMetaResult_pt(:,4): BIP vs. SCZ
CheckBtPt.CoD(:,1)=Site.CoD(sum(iniMetaResult_pt(Site.Trans,4:6)<0.05,2)==0); %no difference between three patient groups
CheckBtPt.DEP(:,1)=Site.one(sum(iniMetaResult_pt(Site.one,4:5)<0.05,2)==2); %significant difference between: DEP/BIP and DEP/SCZ
CheckBtPt.BIP(:,1)=Site.two(sum(iniMetaResult_pt(Site.two,4:2:6)<0.05,2)==2); %significant difference between: DEP/BIP and BIP/SCZ
CheckBtPt.SCZ(:,1)=Site.three(sum(iniMetaResult_pt(Site.three,5:6)<0.05,2)==2); %significant difference between: DEP/SCZ and BIP/SCZ

Site.CoD=CheckBtPt.CoD;
Site.one=CheckBtPt.DEP;
Site.two=CheckBtPt.BIP;
Site.three=CheckBtPt.SCZ;      
