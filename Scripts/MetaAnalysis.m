function output = MetaAnalysis(effectSize, variance, pValue, sampleSize, type)
% three types of meta analysis: type = 1,2,3
    % 1. Z score meta for p value  
    % 2. fixed effect model for beta   
    % 3. DerSimonian and Laird's method
% input: effectSize (Zscore meta requires the sign for effect size on each 
%        data set only), variance(for beta meta only), pValue(for Zsocre
%        meta only), sampleSize(Zscore meta only)
% output: EffectSizeMeta, pMeta, ZScore(Zscore meta only),
%         pHeterogeneity(random effect model only)
% reference: Evangelou and Ioannidis, Meta-analysis methods for genome-wide
% association studies and beyond, Nature Reviews Genetics, 4,379-389, 2013.
k  = length(effectSize); % number of studies
if type == 1
    % 1. Z score meta for p value  
    for i = 1 : k
       Z(i) = norminv(1-pValue(i));
       if isinf(Z(i))
           % the pavlue is very very small
           Z(i) = -norminv(pValue(i));
       end
       w(i) = sqrt(sampleSize(i));
    end
    output.ZScore = sum(Z.*w) / sqrt(sum(w.^2));
    output.pMeta = normcdf(output.ZScore, 'upper');
elseif type == 2
    % 2. fixed effect model for beta   
    for i = 1 : k
        w(i) = 1 / variance(i);
    end
    EffectSizeMeta = sum(w.*effectSize)/sum(w);
    VarianceMeta = 1 / sum(w);
    pMeta = normcdf(abs(EffectSizeMeta)/sqrt(VarianceMeta), 'upper');
    output.EffectSizeMeta = EffectSizeMeta;
    output.pMeta = 2*pMeta;   
elseif type == 3
    % 3. DerSimonian and Laird's method
    for i = 1 : k
        w(i) = 1 / variance(i);
    end
    fixedMeta = sum(w.*effectSize)/sum(w);
    Q = max(sum(w.*(effectSize - fixedMeta).^2),0);
    pHeterogeneity = 1 - chi2cdf(Q, k-1);
    tau2 = max((Q-k+1) / (sum(w) - sum(w.^2)/sum(w)),0);
    for i = 1 : k
        wR(i) = 1 / (1/w(i)+tau2);
    end
    EffectSizeMeta = sum(wR.*effectSize)/sum(wR);
    VarianceMeta = 1 / sum(wR);
    pMeta = normcdf(abs(EffectSizeMeta)/sqrt(VarianceMeta),'upper');
    output.EffectSizeMeta = EffectSizeMeta;
    output.pMeta = 2*pMeta; % two sides
    output.pHeterogeneity = pHeterogeneity;
else
    'not applicable -- MetaAnalysis.m'
end