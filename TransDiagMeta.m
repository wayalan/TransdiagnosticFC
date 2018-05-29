% fMRIDataSet is a structure for fMRI data and covariate data at each site
% together with the site name and covariates names at each site
% fMRIDataSet{i}.normal, fMRIDataSet{i}.id_n, fMRIDataSet{i}.patient, fMRIDataSet{i}.id_p, fMRIDataSet{i}.site_name,...
%             fMRIDataSet{i}.behaviours(:,fMRIDataSet{i}.index_beh), fMRIDataSet{i}.covariate_n, fMRIDataSet{i}.covarate_p, ...
%             [fMRIDataSet{i}.covarate_p, fMRIDataSet{i}.behaviours(:,fMRIDataSet{i}.index_cov)]

%% This is the example for transdiagnostic analysis between three groups

%% 1) 
%% identify the those FCs with survived the Bonferonni correction by meta-analysis on many datasets
load fMRIDataSet_GPreBp_FDok_20170709_Yeo.mat

%% Multiple Comparison Correction
mcc=2; % 1:Bonferroni, 2:FDR, 3:uncorrected
SigP=0.05; % p-value threshold

%% calcluate FC on each site 
rawFC = CalculateFCOnEachDataSite(fMRIDataSet);
% rawFC{1}.patient and rawFC{1}.normal have been vectorized
% number of FC's calculated
numFC = size(rawFC{1}.normal,2);
% number of ROI's considered 
numROI = size(fMRIDataSet{1}.normal{1},2);
for i = 1 : length(fMRIDataSet)
    % sample Size at each site
    sub1 = sum(fMRIDataSet{i}.id_n);
    sub2 = sum(fMRIDataSet{i}.id_p);
    sampleSize(i) =  sub1 + sub2;
end

%% 2)
%% meta-analysis to identify those FC's survived the Bonferoni correction  (iniFC)
iniMetaResult = [];
iniMetaResult_pt = [];
for i = 1 : numFC
    % calculate beta and variance and p value at each site for each link
    % between control and patients
    [beta, p_beta, var_beta, ~, effSize] = CalculateBetaOnEachDataSite(fMRIDataSet,rawFC, i);
    % calculate beta and variance and p value at each site for each link
    % between patients
    [beta_pt, p_beta_pt, var_beta_pt, ~, effSize_pt] = CalBetaBTPatients(fMRIDataSet,rawFC, i);
    % meta for each link
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

if mcc==1 %% Bonferroni
    mcP=SigP/numFC;
    mcP(1:4)=mcP;
elseif mcc==2 %% FDR
	% FDR for the p array of meta-analysis
	mcP(4)=FDR(iniMetaResult(:,8),SigP);
    % FDR for the p values between each of the two groups
	mcP(1:3)=FDR(iniMetaResult(:,4:6),SigP);
        
elseif mcc==3 %% Uncorrected-P
    mcP=SigP;
    mcP(1:4)=mcP;
end

%% 3)
%% Indentify the FC that satisfied the criteria of transdiagnostic and illness-specific model
for i = 1: numFC

	% Transdiagnostic Result: differences found in each of the PT vs HC, as well as having same direction (iniMetaResult(i,8): meta-p)
    if iniMetaResult(i,4) < mcP(1) && iniMetaResult(i,5) < mcP(2) && iniMetaResult(i,6) < mcP(3) && iniMetaResult(i,8) < mcP(4) 
        Site.CoD(n(1,1),1) = i;
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

% PT1 vs HC
    Site.All_1=find(iniMetaResult(:,4)<mcP(1));
% PT2 vs HC
    Site.All_2=find(iniMetaResult(:,5)<mcP(2));
% PT3 vs HC
    Site.All_3=find(iniMetaResult(:,6)<mcP(3));
       

% FCr.CoD=ResFC(:,1:size(FC.CoD,2));
% FCr.uniSite1=ResFC(:,size(FC.CoD,2)+1:size(FC.CoD,2)+size(FC.uniSite1,2));
% FCr.uniSite2=ResFC(:,size(FC.CoD,2)+size(FC.uniSite1,2)+1:size(FC.CoD,2)+size(FC.uniSite1,2)+size(FC.uniSite2,2));
% FCr.uniSite3=ResFC(:,size(FC.CoD,2)+size(FC.uniSite1,2)+size(FC.uniSite2,2)+1:end);

