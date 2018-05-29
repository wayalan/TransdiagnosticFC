%% 0-1 optimization to find best match across datasets
%%  4 years ---early stage
uplimit = 1;

% % %%%%%%%all
% 
% t_index = (Taiwan_Pat_doi > 0 );
% x_index = (Xiangya_Pat_doi > 0 );
% n_index = (Nottingham_Pat_doi > 0);
% h_index = (Huaxi_Pat_doi > 0 );
% cb_index = (Cobre_Pat_doi > 0);

% 
% % % %%%%% first-episode drug-naive
% t_index = (Taiwan_Pat_doi <= -1 );
% x_index = (Xiangya_Pat_doi <= 1 & Xiangya_Pat_naive == 3 );
% n_index = (Nottingham_Pat_doi <= -1);
% h_index = (Huaxi_Pat_doi <= 1 );
% cb_index = (Cobre_Pat_doi <= -1);
% z_index = (Zhongyi_Pat_doi <= 1 &  Zhongyi_Pat_med == 0 );
% % 

% %%%% unmedicated chronic  failed to detect any significant links
% t_index = ( Taiwan_Pat_doi <= -1 );
% x_index = ( Xiangya_Pat_doi > 1 & Xiangya_Pat_naive == 3);
% n_index = ( Nottingham_Pat_doi <= -1);
% h_index = ( Huaxi_Pat_doi > 1 );
% cb_index = (Cobre_Pat_doi<= -1 );


% % %%% medicated early-stage
% t_index = ( Taiwan_Pat_doi <= uplimit );
% x_index = ( Xiangya_Pat_doi <= uplimit & Xiangya_Pat_naive == 2);
% n_index = ( Nottingham_Pat_doi <= uplimit);
% h_index = ( Huaxi_Pat_doi < -1 & Huaxi_Pat_doi <= uplimit );
% cb_index = (Cobre_Pat_doi <= uplimit );
% z_index = (Zhongyi_Pat_doi <= uplimit & Zhongyi_Pat_med == 1);
% % 

% %%%% medicated chronic
t_index = (Taiwan_Pat_doi > uplimit );
x_index = (Xiangya_Pat_doi > uplimit  & Xiangya_Pat_naive == 2);
n_index = (Nottingham_Pat_doi > uplimit );
h_index = (Huaxi_Pat_doi > uplimit );
cb_index = (Cobre_Pat_doi > uplimit );
z_index = (Zhongyi_Pat_doi > 0 & Zhongyi_Pat_med == 1 & Zhongyi_Pat_stage == 2);

% % %%%% all  
% t_index = (Taiwan_Pat_doi > -1 );
% x_index = (Xiangya_Pat_doi > -1 );
% n_index = (Nottingham_Pat_doi > -1 );
% h_index = (Huaxi_Pat_doi > -1 );
% cb_index = (Cobre_Pat_doi > -1 );

disp('tw')
% t_index_nor = ones(1,length(Taiwan_Nor_a_s_e_h_m));
t_index_nor = zeros(1,length(Taiwan_Nor_a_s_e_h_m));
[t_index, t_index_nor] = myMatch4(Taiwan_Pat_a_s_e_h_m(:,[1,2,3,5,6]), ...
    Taiwan_Nor_a_s_e_h_m(:,[1,2,3,5,6]), t_index, t_index_nor);


disp('xy')
% x_index_nor = ones(1,length(Xiangya_Nor_a_s_e_h_m));
x_index_nor = zeros(1,length(Xiangya_Nor_a_s_e_h_m));
[x_index, x_index_nor] = myMatch4(Xiangya_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    Xiangya_Nor_a_s_e_h_m(:,[1,2,3,5,6]), x_index, x_index_nor);

disp('nt')
% n_index_nor =  ones(1,length(Nottingham_Nor_a_s_e_h_m));
n_index_nor =  zeros(1,length(Nottingham_Nor_a_s_e_h_m));
[n_index, n_index_nor] = myMatch4(Nottingham_Pat_a_s_e_h_m(:,[1,2,5,6]),...
    Nottingham_Nor_a_s_e_h_m(:,[1,2,5,6]), n_index, n_index_nor);



disp('hx')
% h_index_nor = ones(1,length(Huaxi_Nor_a_s_e_h_m));
h_index_nor = zeros(1,length(Huaxi_Nor_a_s_e_h_m));
[h_index, h_index_nor] = myMatch4(Huaxi_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    Huaxi_Nor_a_s_e_h_m(:,[1,2,3,5,6]), h_index, h_index_nor);

disp('cb')
% cb_index_nor = ones(1,length(Cobre_Nor_a_s_e_h_m));
cb_index_nor = zeros(1,length(Cobre_Nor_a_s_e_h_m));
[cb_index, cb_index_nor] = myMatch4(Cobre_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    Cobre_Nor_a_s_e_h_m(:,[1,2,3,5,6]), cb_index, cb_index_nor);


% t_index_nor = (t_index_nor > -1);
% x_index_nor = (x_index_nor > -1);
% n_index_nor = (n_index_nor > -1);
% h_index_nor = (h_index_nor > -1);
% cb_index_nor = (cb_index_nor > -1);

disp([sum(t_index), sum(t_index_nor)])
disp([sum(x_index), sum(x_index_nor)])
disp([sum(n_index), sum(n_index_nor)])
disp([sum(h_index), sum(h_index_nor)])
disp([sum(cb_index), sum(cb_index_nor)])

    
%% summerize
% site = 'Taiwan'; sites = 't';
% site = 'Xiangya'; rsites = 'x';
% site = 'Nottingham'; sites = 'n';
% site = 'Huaxi'; sites = 'h';
% site = 'Cobre'; sites = 'cb';

patient = eval([site,'_Pat_a_s_e_h_m']);   p_index = eval([sites, '_index']);
control = eval([site,'_Nor_a_s_e_h_m']);  c_index =  eval([sites, '_index_nor']);
scores_pat =  [eval([site,'_Pat_positive']), eval([site, '_Pat_negative']),...
    eval([site,'_Pat_doi']), eval([site,'_Pat_anhedonia']), eval([site,'_Pat_general'])];

% first line healthy: gender, age, hand, education, p, n , DOI,anhedoina
% second line patient
% gender
numpatient = sum(p_index==1); numcontrol = sum(c_index==1);
sumtable(2,1) = sum(patient(p_index==1, 2) == 1);   sumtable(2,2) = numpatient - sumtable(2,1);
sumtable(1,1) = sum(control(c_index==1, 2) == 1);   sumtable(1,2) = numcontrol - sumtable(1,1);

% age
age = [patient(p_index==1, 1); control(c_index==1, 1)];
id = [ones(1,numpatient), 2*ones(1,numcontrol)];
p = anova1(age,id, 'off');
if p < 0.05
    'age different'
end
sumtable(1,3) = nanmean(age(id==2));   sumtable(1,4) = nanstd(age(id==2));
sumtable(2,3) = nanmean(age(id==1));   sumtable(2,4) = nanstd(age(id==1));

% handedness
sumtable(2,5) = sum(patient(p_index==1, 4) == 1);   sumtable(2,6) = numpatient - sumtable(2,5);
sumtable(1,5) = sum(control(c_index==1, 4) == 1);   sumtable(1,6) = numcontrol - sumtable(1,5);

% education
edu = [patient(p_index==1, 3); control(c_index==1, 3)];
id = [ones(1,numpatient), 2*ones(1,numcontrol)];
p = anova1(edu,id, 'off');
if p < 0.05
    'edu different'
end
sumtable(1,7) = nanmean(edu(id==2));   sumtable(1,8) = nanstd(edu(id==2));
sumtable(2,7) = nanmean(edu(id==1));   sumtable(2,8) = nanstd(edu(id==1));

% motion
motion1 = [patient(p_index==1, 5); control(c_index==1, 5)];
id = [ones(1,numpatient), 2*ones(1,numcontrol)];
p = anova1(motion1,id, 'off');
if p < 0.05
    'motion1 different'
end
sumtable(1,17) = nanmean(motion1(id==2));   sumtable(1,18) = nanstd(motion1(id==2));
sumtable(2,17) = nanmean(motion1(id==1));   sumtable(2,18) = nanstd(motion1(id==1));

motion2 = [patient(p_index==1, 6); control(c_index==1, 6)];
id = [ones(1,numpatient), 2*ones(1,numcontrol)];
p = anova1(motion2, id, 'off');
if p < 0.05
    'motion2 different'
end
sumtable(1,19) = nanmean(motion2(id==2));   sumtable(1,20) = nanstd(motion2(id==2));
sumtable(2,19) = nanmean(motion2(id==1));   sumtable(2,20) = nanstd(motion2(id==1));

% scores: P,N,DOI,Anhedonia,general
scores = scores_pat(p_index==1, :);
sumtable(2,[9,11,13,15,21]) = nanmean(scores);   sumtable(2,[10,12,14,16,22]) = nanstd(scores);


%%  analysis
clear study

tw_id_n = (t_index_nor == 1); tw_id_p = (t_index == 1);
nt_id_n = (n_index_nor == 1); nt_id_p = (n_index == 1);
xy_id_n = (x_index_nor == 1); xy_id_p = (x_index == 1);
hx_id_n = (h_index_nor == 1); hx_id_p = (h_index == 1);
cb_id_n = (cb_index_nor == 1); cb_id_p = (cb_index == 1);

disp('updated')
%% distribution of panss
iOi = 1:14; %index_of_interst
data = [Cobre_Pat_psisub(cb_id_p,iOi); Taiwan_Pat_psisub(tw_id_p,iOi);];% Xiangya_Pat_psisub(xy_id_p,iOi)];
%Huaxi_Pat_psisub(hx_id_p,iOi); %[Xiangya_Pat_psisub(xy_id_p,iOi)];%; Huaxi_Pat_psisub(hx_id_p,iOi)];
meanval = nanmean(data);
varval = nanvar(data);
bar(1:14,meanval, 'facecolor', 'none')
hold on
errorbar(1:14,meanval, sqrt(varval)/2)
hold off
ylabel('severity')
set(gca,'xTick',1:14, 'xTickLabel', xlabels);
xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
xt = get(gca,'XTick'); % 获取横坐标轴刻度句柄
yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄
xtextp=xt; %每个标签放置位置的横坐标，这个自然应该和原来的一样了。
ytextp=(yt(1))*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标
% rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
% 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
% 不同的角度对应不同的旋转位置了，依自己的需求而定了。
text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',90,'fontsize',10);
set(gca,'xticklabel','');% 将原有的标签隐去

%% 
study{1} = link_cal_cov(tw_normal, tw_id_n, tw_patient, tw_id_p, 'Taiwan',...
    Taiwan_Pat_psisub, Taiwan_Nor_a_s_e_h_m(:,[1,2,3,5,6]),...
    Taiwan_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    [Taiwan_Pat_a_s_e_h_m(:,[1,2,3]), Taiwan_Pat_doi]);

study{2} = link_cal_cov(xy_normal, xy_id_n, xy_patient, xy_id_p, 'Xiangya',...
    Xiangya_Pat_psisub, Xiangya_Nor_a_s_e_h_m(:,[1,2,3,5,6]),...
    Xiangya_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    [Xiangya_Pat_a_s_e_h_m(:,[1,2,3]), Xiangya_Pat_medicodosage, Xiangya_Pat_doi]); % Xiangya_Pat_psisub,

study{3} = link_cal_cov(nt_normal, nt_id_n, nt_patient, nt_id_p, 'Nottingham',...
    Nottingham_Pat_doi, Nottingham_Nor_a_s_e_h_m(:,[1,2,5,6]),...
    Nottingham_Pat_a_s_e_h_m(:,[1,2,5,6]),...
    [Nottingham_Pat_a_s_e_h_m(:,[1,2]),Nottingham_Pat_medicodosage]);

study{4} = link_cal_cov(hx_normal, hx_id_n, hx_patient, hx_id_p, 'Huaxi',...
    Huaxi_Pat_psisub, Huaxi_Nor_a_s_e_h_m(:,[1,2,3,5,6]),...
    Huaxi_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    [Huaxi_Pat_a_s_e_h_m(:,[1,2,3]), Huaxi_Pat_doi]);

study{5} = link_cal_cov(cb_normal, cb_id_n, cb_patient, cb_id_p, 'COBRE',...
    Cobre_Pat_psisub, Cobre_Nor_a_s_e_h_m(:,[1,2,3,5,6]),...
    Cobre_Pat_a_s_e_h_m(:,[1,2,3,5,6]),...
    [Cobre_Pat_a_s_e_h_m(:,[1,2,3]),Cobre_Pat_medicodosage, Cobre_Pat_doi]);%Cobre_Pat_psisub, 


disp('Done!')
%% meta-analysis
id = [1,2,3,5,6];
for i = 1 : 4005
    effectSize = [];
    variance = [];
    for j = 1 : length(id)
        effectSize = [effectSize, study{id(j)}{2}(i)];
        variance = [variance, study{id(j)}{2}(4005+i)];
    end
    output = MetaAnalysis(effectSize, variance, [], [], 2);    
    metab(i) = output.EffectSizeMeta;
    metap(i) = output.pMeta;
    output = MetaAnalysis(effectSize, variance, [], [], 3);
    metab_rm(i) = output.EffectSizeMeta;
    metap_rm(i) = output.pMeta;    
    pHeterogeneity(i) = output.pHeterogeneity;
end

%% multiple comparison
ret.prb = metap';
[~, q_fdr] = cca_findsignificance(ret,0.05,2);
PR = (metap < 0.05/4005);


% % ret.prb = metapFE';
% % [~, q_fdr] = cca_findsignificance(ret,0.05,2);
% % PR = (metapFE < q_fdr);

% metabME = metab;
% metapME = metap;
% PRME = PR;
% q_fdrMR = q_fdr;

% metabLO = metab;
% metapLO = metap;
% PRLO = PR;
% q_fdrLO = q_fdr;


%% behaivoural correlation   --  不用Nottingham (3)
id = [1,2,5];
for idk = 1:34  % factors and clusters % subscales of p,n,g, and then, p, n, g, t   
    for i = 1 : 4005
        variance = [];
        effectSize = [];
        for j = 1 : length(id)
            variance = [variance, 1/(study{id(j)}{3}(idk,1)-3)];
            effectSize = [effectSize, gretna_fishertrans(study{id(j)}{3}(idk,i+1))];  % ZCOR
        end
        output = MetaAnalysis(effectSize, variance, [], [], 2);
        metab_corr(idk,i) = (exp(2*output.EffectSizeMeta)-1) / (exp(2*output.EffectSizeMeta)+1);
        metap_corr(idk,i) = output.pMeta;
        output = MetaAnalysis(effectSize, variance, [], [], 3);
        metab_corr_rm(idk,i) = (exp(2*output.EffectSizeMeta)-1) / (exp(2*output.EffectSizeMeta)+1);
        metap_corr_rm(idk,i) = output.pMeta;
        metapHeterogeneity_corr(idk,i) = output.pHeterogeneity;
    end
end   
%    
% save('studiesFirstEpisod_moredata.mat','study', 'metab', 'metap', ...
%   'metap_rm', 'metab_rm', 'pHeterogeneity', 'PR', 'q_fdr', ...
%   'metab_corr', 'metap_corr', 'metab_corr_rm', 'metap_corr_rm', 'metapHeterogeneity_corr')
% save('studiesShortTerm_moredata.mat','study', 'metab', 'metap', ...
%   'metap_rm', 'metab_rm', 'pHeterogeneity', 'PR', 'q_fdr', ...
%   'metab_corr', 'metap_corr', 'metab_corr_rm', 'metap_corr_rm', 'metapHeterogeneity_corr')
% save('studiesLongTerm_moredata.mat','study', 'metab', 'metap', ...
%   'metap_rm', 'metab_rm', 'pHeterogeneity', 'PR', 'q_fdr', ...
%   'metab_corr', 'metap_corr', 'metab_corr_rm', 'metap_corr_rm', 'metapHeterogeneity_corr')

% save('medicatedstages.mat','study', 'metab', 'metap', ...
%   'metap_rm', 'metab_rm', 'pHeterogeneity', 'PR', 'q_fdr', ...
%   'metab_corr', 'metap_corr', 'metab_corr_rm', 'metap_corr_rm', 'metapHeterogeneity_corr')
% save('allstages.mat','study', 'metab', 'metap', ...
%   'metap_rm', 'metab_rm', 'pHeterogeneity', 'PR', 'q_fdr', ...
%   'metab_corr', 'metap_corr', 'metab_corr_rm', 'metap_corr_rm', 'metapHeterogeneity_corr')


%% display the result
clear stage
stage{1} = load('studiesFirstEpisod.mat');  % _factor.mat, _factor_v2.mat
stage{2} = load('studiesShortTerm.mat');
stage{3} = load('studiesLongTerm.mat');
% stage{4} = load('medicatedstages.mat');
% stage{5} = load('allstages.mat');
% %%
% clear study;
% study = stage{1}.study;
% id = [2,4];
% for i = 1 : 4005
%     effectSize = [];
%     variance = [];
%     for j = 1 : length(id)
%         effectSize = [effectSize, study{id(j)}{2}(i)];
%         variance = [variance, study{id(j)}{2}(4005+i)];
%     end
%     output = MetaAnalysis(effectSize, variance, [], [], 2);    
%     metab(i) = output.EffectSizeMeta;
%     metap(i) = output.pMeta;
%     output = MetaAnalysis(effectSize, variance, [], [], 3);
%     metab_rm(i) = output.EffectSizeMeta;
%     metap_rm(i) = output.pMeta;    
%     pHeterogeneity(i) = output.pHeterogeneity;
% end



%% if heterogeneity  <  0.05, fixed effect model failed; the significance holds only
%% when the meta pvalue given by random effect model survives the correction
for i = 1 : 3
    PR = stage{i}.PR;
    finalpvalue = stage{i}.metap;
    finalbvalue = stage{i}.metab;
    rmpvalue = stage{i}.metap_rm;   
    rmbvalue = stage{i}.metab_rm;
    pHeterogeneity = stage{i}.pHeterogeneity;
    for j = 1 : 4005
        if PR(j) == 1 
            if pHeterogeneity(j) < 0.05 
                finalpvalue(j) = rmpvalue(j);
                finalbvalue(j) = rmbvalue(j);
                if rmpvalue(j) < 0.05/4005
                    PR(j) = 1;
                else
                    PR(j) = 0;
                end
            end
        end                
    end
    stage{i}.rePR = PR; 
    ret.prb = finalpvalue';
    [~, q_fdr] = cca_findsignificance(ret,0.05,2);
    disp([q_fdr, sum(finalpvalue<q_fdr)])
    stage{i}.q_fdr = q_fdr;
%     stage{i}.rePR = (finalpvalue<q_fdr);    
    stage{i}.re_pvalue = finalpvalue;
    stage{i}.re_bvalue = finalbvalue;
    disp(sum(PR))
end

%% symptome correlation
for i = 1 : 3
    PR = stage{i}.rePR;
    rmpvalue = stage{i}.metap_corr_rm;    
    pHeterogeneity = stage{i}.metapHeterogeneity_corr;
    fepvalue = stage{i}.metap_corr;
    fervalue = stage{i}.metab_corr;
    for j = 1 : 4005      
       % if PR(j) == 1 
            for k = 1 : size(fepvalue,1)
                if pHeterogeneity(k,j) < 0.05
                    fepvalue(k,j) = rmpvalue(k,j);   
                    fervalue(k,j) = stage{i}.metab_corr_rm(k,j);
                end
            end
       % end                
    end
    stage{i}.recorr_pvalue = fepvalue;
    stage{i}.recorr_rvalue = fervalue;
    ret.prb = fepvalue(:,PR==1); 
    ret.prb = ret.prb(:);
    [~,q_fdr] = cca_findsignificance(ret,0.05,2);
    stage{i}.corrq_fdr = q_fdr;
    stage{i}.corrq_num = sum(ret.prb <= q_fdr);   
    disp([q_fdr, stage{i}.corrq_num])
end

%% display the result
clear resultTable
for i = 1 : 3    
    clear FCname
    clear metap; clear metab; clear metap_rm; clear metab_rm;
    clear pHeterogeneity
    clear metap_corr; clear metab_corr; clear metap_corr_rm; clear metab_corr_rm;
    clear metapHeterogeneity_corr
    clear SiteBValue; clear SitePValue; clear SiteCorr; clear SitePCorr;
    clear SiteGD; clear rePR
    rePR = stage{i}.rePR;
    ind = find(rePR);
    study = stage{i}.study;    
    numbeh = size(stage{i}.metap_corr,1);
    for j = 1 : length(ind)
        FCname{j,1} = edgelabel.edge_roi{ind(j)};
        metap(j,1) = stage{i}.re_pvalue(ind(j));
        metab(j,1) = stage{i}.re_bvalue(ind(j));
        pHeterogeneity(j,1) = stage{i}.pHeterogeneity(ind(j));
        for qq = 1 : numbeh        
            metap_corr(j,qq) = stage{i}.recorr_pvalue(qq,ind(j));
            metab_corr(j,qq) = stage{i}.recorr_rvalue(qq,ind(j));
            metapHeterogeneity_corr(j,qq)  = stage{i}.metapHeterogeneity_corr(qq,ind(j));
        end
        for k = 1 : 5
            if size(study{k},1) > 0                
                SiteBValue(j,k) = study{k}{2}(ind(j));
                SitePValue(j,k) = study{k}{4}(ind(j));
                if k ~= 3
                    for qq = 1 : numbeh
                        SiteCorr(j,(k-1)*34+qq) = study{k}{3}(qq,ind(j)+1);
                        SitePCorr(j,(k-1)*34+qq) = study{k}{5}(qq,ind(j));
                    end
                else
                    for qq = 1 : numbeh
                        SiteCorr(j,(k-1)*34+qq) = -99;
                        SitePCorr(j,(k-1)*34+qq) = -99;
                    end
                end
                SiteGD(j,k) = abs(study{k}{6}(ind(j),2)) - abs(study{k}{6}(ind(j),1));
            else
                SiteBValue(j,k) = -99;
                SitePValue(j,k) = -99;                
                for qq = 1 : numbeh
                    SiteCorr(j,(k-1)*34+qq) = -99;
                    SitePCorr(j,(k-1)*34+qq) = -99;
                end
                SiteGD(j,k) = -99;
            end            
        end        
    end
    
    resultTable{i} = table( FCname,  SiteBValue,  SitePValue, SiteGD, ...
      metab,metap,  pHeterogeneity, SiteCorr,SitePCorr, metapHeterogeneity_corr,...
       metab_corr, metap_corr,  'VariableNames', {'Link','SiteBValue', ...
       'SitePValue', 'SiteGrpDiff', 'metab',  'metap', 'pHeterogeneity',...
       'SiteCorr', 'SitePCorr', 'metapHeterogeneity_corr', 'metab_corr',  'metap_corr'});
   resultTable2{i} = table(FCname, metab, SitePValue, metap, metap_corr, metab_corr);
  
end
clear FCname
clear metap; clear metab; clear metap_rm; clear metab_rm;
clear pHeterogeneity
clear metap_corr; clear metab_corr; clear metap_corr_rm; clear metab_corr_rm;
clear metapHeterogeneity_corr
clear SiteBValue; clear SitePValue; clear SiteCorr; clear SitePCorr;

%%
for i = 1 : 3    
    writetable(resultTable2{i}, ['stage', num2str(i), '_factor_alldata.csv'])
end

%%
alloverlap = find(stage{1}.rePR & stage{2}.rePR & stage{3}.rePR);
[stage{1}.re_bvalue(alloverlap)', stage{1}.re_pvalue(alloverlap)',...
    stage{2}.re_bvalue(alloverlap)', stage{2}.re_pvalue(alloverlap)',...
    stage{3}.re_bvalue(alloverlap)', stage{3}.re_pvalue(alloverlap)']

clc
for i = 1 : 3
    disp(['stage ', num2str(i)])
    for j = 1 : length(alloverlap)
        indx =  find(stage{i}.recorr_pvalue(:,alloverlap(j)) < 0.05);
        disp(edgelabel.edge_roi(alloverlap(j)))
        for k = 1 : length(indx)
            disp([x(indx(k)), ...
                num2str(stage{i}.recorr_pvalue(indx(k), alloverlap(j)))]);
        end
    end
end
%% write link.txt for circos
clear edgeid
t = 1;
for i = 1 : 89
    for j = i+1 : 90
        edgeid(t,:) = [i,j];
        t = t + 1;
    end
end


%%
% clear linkTable
for i = 1 : 3
%     clear stagingedges; clear stagingConnType; clear stagingConnScor; clear ConnType; clear ConnScor;
%     clear hemisphere1; clear parcename1; clear parcename2; clear hemisphere2;
%     stagingedges = edgeid(stage{i}.rePR==1,:);
%     stagingConnType = stage{i}.metab(stage{i}.rePR==1); % connection type; reduced -1, enhanced +1
%     stagingConnScor = -log10(stage{i}.metap(stage{i}.rePR==1)); % connection score: -log10(p)
%     for j = 1 : length(stagingedges)
%         clear tempname
%         tempname = aal90{stagingedges(j,1),1};
%         hemisphere1{j} = lower(tempname(end));
%         parcename1{j} = tempname(1:end-2);
%         clear tempname
%         tempname = aal90{stagingedges(j,2),1};
%         hemisphere2{j} = lower(tempname(end));
%         parcename2{j} = tempname(1:end-2);        
%         
%         if stagingConnType(j) < 0
%             ConnType(j) = -1;
%         else
%             ConnType(j) = 1;
%         end
%         ConnScor(j) = stagingConnScor(j);
%     end   
%     linkTable{i} = table(hemisphere1', parcename1', hemisphere2', parcename2', ConnType', ConnScor');
    writetable(linkTable{i}, ['link', num2str(i), '.xls'])
end

%% display behaivoural correlation
clear xlabels
% for i = 1 : 30
%     xlabels{i} = COBRE_panss.Properties.VariableNames{i+1};
% end
% xlabels{31} = 'Positive Symptome';
% xlabels{32} = 'Negative Symptome';
% xlabels{33} = 'General Symptome';
% xlabels{34} = 'Total Symptome';

for i = 1 : 14
    xlabels{i} = COBRE_panss.Properties.VariableNames{i+1};
end
xlabels{15} = 'Positive Symptome';
xlabels{16} = 'Negative Symptome';
xlabels{17} = 'General Symptome';
xlabels{18} = 'Total Symptome';

stageunion = (stage{1}.rePR | stage{2}.rePR | stage{3}.rePR);
ylabels = [];
t = 1;
for i = 1 : 4005
    if stageunion(i) == 1
        ylabels{t} = edgelabel.edge_roi{i};
        t = t + 1;
    end
end
s1=1;
s2=1;
for i = 1 : length(ylabels)
    if i  <= length(ylabels)/2
        halfylabels1{s1} = ylabels{i};
        s1=s1+1;
    else
        halfylabels2{s2} = ylabels{i};
        s1=s1+1;
    end
end
    
nedge = sum(stageunion);
rmatrixtoshow = zeros(18,nedge,3);
for i = 1 : 3    
    for j = 1 : 18
        clear indd
        clear tempr
        indd = find(stage{i}.recorr_pvalue(j,stageunion) < 0.05);
        tempr = stage{i}.recorr_rvalue(j,stageunion);       
        rmatrixtoshow(j, indd,i) =  tempr(indd);
    end
end

%% progressive change of behavioural correlation
figure
index = [71,72,75,76];% [61,62,64,65]; %19,20,24,25]; %[71,72,75,76]; % 'PoCG.L-THA.L'    'PoCG.L-THA.R'   'PoCG.R-THA.L'    'PoCG.R-THA.R' 
for i = 1 : length(index) 
    subplot(1,4,i)
    for j = 1 : 3
        corrchange(:,j) = rmatrixtoshow(:,index(i),j);
    end
    imagesc(corrchange)
    MyColorMap=[
        0,0,0
        0,0,1
        0,1,1        
        1,1,1        
        1,0,1
        1,1,0
        1,0,0
        ];
    colormap(MyColorMap)    
    caxis([-0.5,0.5])
    if i == 4
        colorbar
    end
    set(gca,'yTick',1:length(xlabels), 'yTickLabel', xlabels, 'xTick', 1:3, 'xTickLabel', {'untreated first-episode', 'treated early-stage', 'trated chronic'});
    xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄
    ytb = get(gca, 'YTickLabel');
    xt = get(gca,'XTick'); % 获取横坐标轴刻度句柄
    yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄
    xtextp=xt; %每个标签放置位置的横坐标，这个自然应该和原来的一样了。
    ytextp=(yt(end)+1)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标
    yxtextp=yt;
    yytextp=(xt(1)-0.5)*ones(1,length(yt));
    % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，
    % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
    % 不同的角度对应不同的旋转位置了，依自己的需求而定了。
    text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',90,'fontsize',10);
    set(gca,'xticklabel','');% 将原有的标签隐去
    if i == 1
    text(yytextp, yxtextp, ytb,'HorizontalAlignment','right','rotation',0,'fontsize',8);
    end
    set(gca,'yticklabel','');% 将原有的标签隐去 
    title(ylabels{index(i)})
end
%%
figurename = {'untreated first-episode', 'treated early stage', 'treated chronic'};
for i = 1 : 3
    figure(i)
    title(figurename{i})
%     ntimes = 1;
%     tempmatrix = zeros(size(rmatrixtoshow,1),ntimes*size(rmatrixtoshow,2));
%     for j = 1 : size(rmatrixtoshow,2)
%         tempmatrix(:,ntimes*j-1) = rmatrixtoshow(:,j,i);
%     end
    h = imagesc(rmatrixtoshow(:,:,i)');
    set(gca,'xTick',1:length(xlabels), 'xTickLabel', xlabels, 'yTick', 1:length(ylabels), 'yTickLabel', ylabels)
    
    xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄    
    ytb = get(gca, 'YTickLabel');
    xt = get(gca,'XTick'); % 获取横坐标轴刻度句柄    
    yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄    
    xtextp=xt; %每个标签放置位置的横坐标，这个自然应该和原来的一样了。                    
    ytextp=(yt(end)+2)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标   
    yxtextp=yt;
    yytextp=(xt(1)-2)*ones(1,length(yt));
    % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，  
    % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
    % 不同的角度对应不同的旋转位置了，依自己的需求而定了。
    text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',90,'fontsize',10);
    set(gca,'xticklabel','');% 将原有的标签隐去
    text(yytextp, yxtextp, ytb,'HorizontalAlignment','right','rotation',0,'fontsize',8);
    set(gca,'yticklabel','');% 将原有的标签隐去
end

%% three stages
%% display behaivoural correlation
clear xlabels
for i = 1 : 7
    xlabels{i} = ['P',num2str(i)]; %[COBRE_panss.Properties.VariableNames{i+1};
end
for i = 8 : 14
    xlabels{i} = ['N', num2str(i-7)];
end
% xlabels{15} = 'Positive Symptome';
% xlabels{16} = 'Negative Symptome';
% xlabels{17} = 'General Symptome';
% xlabels{18} = 'Total Symptome';
clear ylabels
for j = 1 : 3
    t = 1;
    for i = 1 : 4005
        if stage{j}.rePR(i) == 1
            ylabels{j}{t} = edgelabel.edge_roi{i};
            t = t + 1;
        end
    end
    nedge(j) = sum(stage{j}.rePR);
end

clear rmatrixtoshow
for i = 1 : 3    
    rmatrixtoshow{i} = zeros(14,nedge(i));
    for j = 1 : 14
        clear indd
        clear tempr
        indd = find(stage{i}.recorr_pvalue(j,stage{i}.rePR) < 0.05);
        tempr = stage{i}.recorr_rvalue(j,stage{i}.rePR); 
        tempp = stage{i}.recorr_pvalue(j,stage{i}.rePR);
        for ii = 1 : length(indd)
            if tempp(indd(ii)) < 0.05/14
                tempp(indd(ii)) = 2;
            else
                tempp(indd(ii)) = 1;
            end            
        end
        rmatrixtoshow{i}(j, indd) =  sign(tempr(indd)) .* tempp(indd);
    end
end
%%
figurename = {'untreated first-episode', 'treated early stage', 'treated chronic'};
for i = 1 : 3
    figure(i)
%     ntimes = 1;
%     tempmatrix = zeros(size(rmatrixtoshow,1),ntimes*size(rmatrixtoshow,2));
%     for j = 1 : size(rmatrixtoshow,2)
%         tempmatrix(:,ntimes*j-1) = rmatrixtoshow(:,j,i);
%     end
    h = imagesc(rmatrixtoshow{i}');
    MyColorMap=[
        0,0,0
        0,0,1        
        1,1,1        
        1,1,0
        1,0,0
        ];
    colormap(MyColorMap)    
    caxis([-2,2])
    set(gca,'xTick',1:length(xlabels), 'xTickLabel', xlabels, 'yTick', 1:length(ylabels{i}), 'yTickLabel', ylabels{i})
    
    xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄    
    ytb = get(gca, 'YTickLabel');
    xt = get(gca,'XTick'); % 获取横坐标轴刻度句柄    
    yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄    
    xtextp=xt; %每个标签放置位置的横坐标，这个自然应该和原来的一样了。                    
    ytextp=(yt(end)+2)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标   
    yxtextp=yt;
    yytextp=(xt(1)-1)*ones(1,length(yt));
    % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，  
    % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
    % 不同的角度对应不同的旋转位置了，依自己的需求而定了。
    %text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',0,'fontsize',12);
    %set(gca,'xticklabel','');% 将原有的标签隐去
    text(yytextp, yxtextp, ytb,'HorizontalAlignment','right','rotation',0,'fontsize',10);
    set(gca,'yticklabel','');% 将原有的标签隐去
    title(figurename{i})

end

%% biclustering 
rawcluster{1} = [2 2 0 0 2 0 0 2 2 2 0 2 0 0];
colcluster{1} = [1 1 1 0 2 1 1 1 2 0 1 1 0 2 2 1 2 1 1 1 1 2 1 1 1 2 1 2 1];

rawcluster{2} = [2 2 2 2 2 2 2 2 2 2 2 1 2 1];
colcluster{2} = [0 1 0 1 0 0 2 0 0 1 1 1 2 2 2 2 0 0];

rawcluster{3} = [2 2 2 2 0 2 2 2 0 2 2 2 2 2];
colcluster{3} = [2 2 1 1 1 1 2 0 0 0 1 2 1 1 1 2 1 0 1 0 1 0 2 0 2 2 0 2 2 1 0 0 1 2 2 2 2 0 2 2 2 1 2 0 0 0 2 0 0 0 0 0 2];

for i = 1 : 3
    [~,b_sym{i}] = sort(rawcluster{i});  % for x 
    [~,b_fc{i}] = sort(colcluster{i});   % for y
end
clear xlabels2
clear ylabels2
clear rmatrixtoshow2
for i = 1 : 3
    for j = 1 : length(xlabels)
        xlabels2{i}{j} = xlabels{b_sym{i}(j)};
    end
    for j = 1 : length(ylabels{i})
        ylabels2{i}{j} = ylabels{i}{b_fc{i}(j)};
    end
    rmatrixtoshow2{i} = rmatrixtoshow{i}(b_sym{i},b_fc{i});
end


%%
figurename = {'untreated first-episode', 'treated early stage', 'treated chronic'};
for i = 1 : 3
    figure(i)

    h = imagesc(rmatrixtoshow2{i}');
    MyColorMap=[
        0,0,0
        0,0,1        
        1,1,1        
        1,1,0
        1,0,0
        ];
    colormap(MyColorMap)    
    caxis([-2,2])
    set(gca,'xTick',1:length(xlabels2{i}), 'xTickLabel', xlabels2{i}, 'yTick', 1:length(ylabels2{i}), 'yTickLabel', ylabels2{i})
    
    xtb = get(gca,'XTickLabel');% 获取横坐标轴标签句柄    
    ytb = get(gca, 'YTickLabel');
    xt = get(gca,'XTick'); % 获取横坐标轴刻度句柄    
    yt = get(gca,'YTick'); % 获取纵坐标轴刻度句柄    
    xtextp=xt; %每个标签放置位置的横坐标，这个自然应该和原来的一样了。                    
    ytextp=(yt(end)+2)*ones(1,length(xt)); % 设置显示标签的位置，写法不唯一，这里其实是在为每个标签找放置位置的纵坐标   
    yxtextp=yt;
    yytextp=(xt(1)-1)*ones(1,length(yt));
    % rotation，正的旋转角度代表逆时针旋转，旋转轴可以由HorizontalAlignment属性来设定，  
    % 有3个属性值：left，right，center，这里可以改这三个值，以及rotation后的角度，这里写的是45
    % 不同的角度对应不同的旋转位置了，依自己的需求而定了。
    %text(xtextp,ytextp,xtb,'HorizontalAlignment','right','rotation',0,'fontsize',12);
    %set(gca,'xticklabel','');% 将原有的标签隐去
    text(yytextp, yxtextp, ytb,'HorizontalAlignment','right','rotation',0,'fontsize',10);
    set(gca,'yticklabel','');% 将原有的标签隐去
    title(figurename{i})
    grid on
end

