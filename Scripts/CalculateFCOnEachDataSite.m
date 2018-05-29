function study = CalculateFCOnEachDataSite(fMRIDataSet)
% input: the structure of fMRI data set with TS data and Covariate data
% output: result.normal and result.patient    nsubj X 4005 FC's
for i = 1 : length(fMRIDataSet)    
    % single FC
    study{i}.normal = cal_FC(fMRIDataSet{i}.normal, fMRIDataSet{i}.id_n);
    study{i}.patient = cal_FC(fMRIDataSet{i}.patient, fMRIDataSet{i}.id_p);
    disp(['Dataset',num2str(i), ' is Done!'])
end

function corrZ_N = cal_FC(data, id)
% id must be boolean
tc = data(id);
sub1 = length(tc);
for i=1:sub1
    [corrN(:,:,i), ~] = corrcoef(tc{i});
end
corrZ_N = gretna_fishertrans(corrN);
corrZ_N = mytri(corrZ_N);
