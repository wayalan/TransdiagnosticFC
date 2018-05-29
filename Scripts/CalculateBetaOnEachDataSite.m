function [beta, p_beta, var_beta, lambda, effSize] = CalculateBetaOnEachDataSite(fMRIDataSet, RawFC, targetFC, targetPh, coeff0)
% input: the structure of fMRI data set with TS data and Covariate data
%        the structure of raw functional connectivity on each data set
%        the selected FC id (targetFC) in RawFC to calculate beta value of
%        group
%        the selected behaivoural id (targtPh) in fMRIDataSet{1}.index_beh(id)
% output: beta, significance of beta, variance of beta, lambda (multivariate case only)
% the FC as the response, group and covariantes as the predictors, using
% regressing model to fit for group values, as well as the
% corresponding p-value and variance.

for i = 1 : length(fMRIDataSet)
    dim = size(targetFC,2);
    id_n = fMRIDataSet{i}.id_n; sub1 = sum(id_n);
    id_p = fMRIDataSet{i}.id_p; sub2 = sum(id_p);
    covariate_Nor = fMRIDataSet{i}.covariate_n(id_n,:);
    covariate_Pat = fMRIDataSet{i}.covariate_p(id_p,:);
    % creating X
    x1 = [zeros(sub1,1); ones(sub2,1)]; % group
    % creating covariantes
    x2 = [covariate_Nor; covariate_Pat]; % covariate
    X = [x1,x2];
    if dim == 1
        if nargin < 4
            % one dimention response
            % For meta analysis
            % the effect size is the beta coefficient (b1) of the regression model
            % FC = b0 + b1*group + b2*covariant + e
            % we need to calculate the estimation of b1 and its sample variance
            % creating Y
            data = [RawFC{i}.normal(:,targetFC); RawFC{i}.patient(:,targetFC)];  % FC
            % do the regression
            [b,~,stats] = glmfit(X, data); % the glmfit adds one column of ones by default
            beta(i) = b(2);
            var_beta(i) = stats.covb(2,2);
%             effSize(i) = (2*stats.t(2))/sqrt(stats.dfe); % with equal n
            effSize(i) = (stats.t(2)*(sub1+sub2))/(sqrt(stats.dfe)*sqrt(sub1*sub2)); % with nonequal n
            p_beta(i) = stats.p(2);
            lambda(i) = beta(i);
        else
            % behavioural correlation in patient group
            % For meta analysis
            % the effect size is the beta coefficient (b1) of the regression model
            % Behaviour = b0 + b1*FC + b2*covariant + e
            % we need to calculate the estimation of b1 and its sample variance 
            % creating Y
            behaviour = fMRIDataSet{i}.behaviours(id_p,:);
            index_beh = fMRIDataSet{i}.index_beh;
            index_cov = fMRIDataSet{i}.index_cov;
            data_corr = behaviour(:,index_beh(targetPh));
            data_cov =  [covariate_Pat, behaviour(:,index_cov)];
            %             % creating X
            %             X_corr = [RawFC{i}.patient(:,targetFC), covariate_Pat, behaviour(:,index_cov)];
            %             % do the regression
            %             [b,~,stats] = glmfit(X_corr, data_corr); % the glmfit adds one column of ones by default
            %             beta(i) = b(2);
            %             var_beta(i) = stats.covb(2,2);
            %             p_beta(i) = stats.p(2);
            %             lambda(i) = beta(i);            
            [beta(i), p_beta(i)] = partialcorr(RawFC{i}.patient(:,targetFC), data_corr, ...
               data_cov, 'row', 'pairwise');
           
            lambda(i) = beta(i);
            var_beta(i) = 1 / (length(data_corr)-3);            
        end
    else
        if nargin < 4
            % multivariate resposes
            t1 = table(X(:,1), 'VariableNames', {'group'});
            predictors = 'group';
            for k = 2 : size(X,2)
                t1_temp = table(X(:,k), 'VariableNames', fMRIDataSet{i}.covariate_name(k-1));
                t1 = [t1, t1_temp];
                predictors = [predictors,'+', fMRIDataSet{i}.covariate_name{k-1}];
            end
            % making data table for manova
            data = [RawFC{i}.normal(:,targetFC); RawFC{i}.patient(:,targetFC)];
            for k = 1 : dim
                t1_temp = table(data(:,k), 'VariableNames', {['Y',num2str(k)]});
                t1 = [t1,t1_temp];
            end
            % fit the multivariate model
            rm = fitrm(t1, ['Y1-Y',num2str(dim),'~ ', predictors]);
            tbl = manova(rm);
            lambda(i) = tbl{6,4}; % wilks' lambda
            beta(i) = tbl{6,5}; % F statistic
            var_beta(i) = Fvar(tbl{6,7},tbl{6,8}); % variance of F statistic
            p_beta(i) = tbl{6,9};
        else
            behaviour = fMRIDataSet{i}.behaviours(id_p,:);
            index_beh = fMRIDataSet{i}.index_beh;
            index_cov = fMRIDataSet{i}.index_cov;
            data_corr = behaviour(:,index_beh(targetPh));
            data_cov = [behaviour(:,index_cov),covariate_Pat];
            data = RawFC{i}.patient(:,targetFC);
            if nargin<5
                % cannonical correlation  ���������һ��Ľ��
                % ע��ý��ķ��ز���Ϊ��meta ���ϵ���pֵ������Ϊ��meta���ϵ��                
                % model 1
                [b,~, stats] = glmfit([data, data_cov], data_corr);
                yhat = glmval(b, [data, data_cov], 'identity');
                RSS1 = nansum((data_corr - yhat).^2);
                df1 = size(data,1) - length(b);
                nvar = size(data,2);
                coeff = b(1+1:1+nvar);
                for kk = 1 : nvar
                    var_beta(kk,i) = stats.covb(1+kk,1+kk);
                end              
                [beta(i), p_beta(i)] = partialcorr(data_corr, data * coeff, data_cov, 'row', 'pairwise');
                lambda(:,i) = coeff;
            else                
                [beta(i), p_beta(i)] = partialcorr(data_corr, data * coeff0, data_cov, 'row', 'pairwise');
                lambda(:,i) = coeff0;
                var_beta(i) = 1/(size(data,1)-3);
            end
            
            
            
%             % model 2
%             b2 = glmfit(data_cov, data_corr);
%             yhat = glmval(b2, data_cov,  'identity');
%             RSS0 = nansum((data_corr-yhat).^2);
%             df0 = size(data,1) - length(b2);            
%             beta(i) = RSS0*df1/(RSS1*df0); % F statistic
%             var_beta(i) = Fvar(df0,df1); % variance of F statistic
%             p_beta(i) = 1 - fcdf(beta(i), df0, df1);  
%             lambda(i) = 0;
            
            %     FC = score + covariate   % Ч��ã���������ص�
%             % (multiple) behaivoural correlation of multivariate resposes             
%             behaviour = fMRIDataSet{i}.behaviours(id_p,:);
%             index_beh = fMRIDataSet{i}.index_beh;
%             index_cov = fMRIDataSet{i}.index_cov;
%             data_corr = behaviour(:,index_beh(targetPh));
%             
%             % discretized 
%             
%             low = quantile(data_corr,0.33);
%             med = quantile(data_corr,0.66);
%             index_low = (data_corr <= low);
%             index_hig = (data_corr > med);
%             index_med = (data_corr >low & data_corr <=med);
%             data_corr(index_low) = 0;
%             data_corr(index_med) = 1;
%             data_corr(index_hig) = 2;  
%             
%             data_cov = behaviour(:,index_cov);
%             
%             % forming predictors
%             % behaivours
%             t1 = table(data_corr(:,1), 'VariableNames', {'b1'});
%             predictors = 'b1';
%             for k = 2 : size(data_corr,2)
%                 t1_temp = table(data_corr(:,k), 'VariableNames', {['b',num2str(k)]});
%                 t1 = [t1, t1_temp];
%                 predictors = [predictors,'+', 'b',num2str(k)];
%             end
%             % covariates -- doi or medication
%             for k = 1 : size(data_cov,2)
%                 t1_temp = table(data_cov(:,k), 'VariableNames', {['c',num2str(k)]});
%                 t1 = [t1, t1_temp];
%                 predictors = [predictors,'+', 'c',num2str(k)];
%             end
%             % covariates -- age or edu etc
%             for k = 1 : size(covariate_Pat,2)
%                 t1_temp = table(covariate_Pat(:,k), 'VariableNames', fMRIDataSet{i}.covariate_name(k));
%                 t1 = [t1, t1_temp];
%                 predictors = [predictors,'+', fMRIDataSet{i}.covariate_name{k}];
%             end
%             
%             % forming the data
%             data = RawFC{i}.patient(:,targetFC);
%             for k = 1 : dim
%                 t1_temp = table(data(:,k), 'VariableNames', {['Y',num2str(k)]});
%                 t1 = [t1,t1_temp];
%             end
%             
%             % fit the multivariate model
%             rm = fitrm(t1, ['Y1-Y',num2str(dim),'~ ', predictors]);
%             tbl = manova(rm);
%             lambda(i) = tbl{6,4}; % wilks' lambda
%             beta(i) = tbl{6,5}; % F statistic
%             var_beta(i) = Fvar(tbl{6,7},tbl{6,8}); % variance of F statistic
%             p_beta(i) = tbl{6,9};
        end
        
    end
end


function variance = Fvar(df1, df2)
if df2 > 4
   variance = 2 * df2 * (df1+df2-2) / (df1*(df2-2)^2*(df2-4));
else
    'the variance of F statsitics cannot be estimated'
    variance = -9999;
end






% %% write to a file where R can read
% fid = fopen(FileName{1},'wt');
% % study source sub1, mean_N, std_N, sub2, mean_P, std_P
% for i = 1 : size(study,2)
%     fprintf(fid, '%d\t%s\t', i, study{i}{1});
%     for j = 1 : size(study{i}{2},2)
%         fprintf(fid, '%f', study{i}{2}(j)); 
%         if j < size(study{i}{2},2)
%            fprintf(fid, '\t');
%         end
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);
% %%
% for k = 2 : length(FileName)   
%     fid = fopen(FileName{k},'wt');
%     % study source sub1, corr
%     for i = 1 : size(study,2)
%         fprintf(fid, '%d\t%s\t%d\t', i, study{i}{1}, study{i}{3}(k,1));
%         for j = 2 : size(study{i}{3},2)
%             fprintf(fid, '%f', study{i}{3}(k,j));
%             if j < size(study{i}{3},2)
%                 fprintf(fid, '\t');
%             end
%         end
%         fprintf(fid,'\n');
%     end
%     fclose(fid);
% end
% %%
% disp('written')


