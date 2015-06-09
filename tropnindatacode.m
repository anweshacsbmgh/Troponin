clear

%% Separate the individual patient data from the list
load('tropdata.mat');
% S=table2struct(Troponinexcel);
Troponinexcel(679:682,:)=[];
Troponinexcel(663:666,:)=[];
Troponinexcel(651:656,:)=[];
Troponinexcel(626:628,:)=[];
Troponinexcel(618:621,:)=[];
Troponinexcel(558:564,:)=[];
Troponinexcel(545:547,:)=[];
Troponinexcel(537:539,:)=[];
Troponinexcel(529:531,:)=[];
Troponinexcel(514:523,:)=[];
Troponinexcel(485:489,:)=[];
Troponinexcel(457:459,:)=[];
Troponinexcel(408:410,:)=[];
Troponinexcel(395:396,:)=[];
Troponinexcel(381:382,:)=[];
Troponinexcel(331,:)=[];
Troponinexcel(328,:)=[];
Troponinexcel(325,:)=[];
Troponinexcel(321,:)=[];
Troponinexcel(304:309,:)=[];
Troponinexcel(300,:)=[];
Troponinexcel(269,:)=[];
Troponinexcel(250:251,:)=[];
Troponinexcel(197,:)=[];
Troponinexcel(142:145,:)=[];
Troponinexcel(81:82,:)=[];
Troponinexcel(58:62,:)=[];
Troponinexcel(23:24,:)=[];

complete_feature_set = [10,11,36:63,113:128,158:167,176:179];
% Troponinexcel(isnan(table2array(Troponinexcel(:,complete_feature_set(:,2)))),complete_feature_set(:,2)) = num2cell(5e4);
numOfRows=size(Troponinexcel(:,1));
patientnum = Troponinexcel(:,12);
a = [-2; find(~isnan(table2array(patientnum)))]';
count=1;
for i=1:numOfRows
if count<numOfRows(1,1)
    Troponinexcel(count,158:167) = array2table(max(table2array(Troponinexcel(count,158:167)),table2array(Troponinexcel(count+1,158:167))));
    Troponinexcel(count,176:179) = array2table(max(table2array(Troponinexcel(count,176:179)),table2array(Troponinexcel(count+1,176:179))));
    
    if ismember(count,a)==1
        count=count+3;
    else
        count=count+2;
    end
end
end
binout=0;
binhsout=0;
gend=0;
for i=1:length(a)-1
    Troponinexcel.VarName12(a(i)+3:a(i+1)+1)=i;
    filename=strcat('pat_data',num2str(i),'.mat');
    S=Troponinexcel(a(i)+3:a(i+1)+1,:);
    S.bintnt=max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,10)))>100.*ones(size(Troponinexcel(a(i)+3:a(i+1)+1,10)));
    S.binhstnt=(max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,11)))>34.2).*strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Male')+(max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,11)))>15.4).*strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Female');
    save(filename,'S');
    binout(1,end+1:end+length(S.bintnt)+1)=[S.bintnt;0];
    binhsout(1,end+1:end+length(S.binhstnt)+1)=[S.binhstnt;0];
    gend(1,end+1:end+length(S.Sex)+1)=[strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Female');0];
end
binout(1)=[];
binhsout(1)=[];
gend(1)=[];
binout(end)=[];
binhsout(end)=[];
gend(end)=[];
%
Troponinexcel.binTNT=binout';
Troponinexcel.binhsTNT=binhsout';
Troponinexcel.binGender=gend';
%% Rearranging the tables into matrix to perform logistic regression
Safe_pat_xls = [12 13 18 22 23 29 31 33 37 40 41 54 55 56 59 60 61 68 71 75 77 79 82 88 89 93 95 100 101 104 108]; %green markers on excel sheet

complete_feature_set = [10,11,36:39,41:52,54:63,113:128,158:167,176:179];
for i=1:length(a)-1
    filename=strcat('pat_data',num2str(i),'.mat');
    load (filename)
    [sec_datapoint pos]=max(S.VarName10(3:end));
    patientData(i,:) = [table2array(S(1,complete_feature_set)),table2array(S(pos,complete_feature_set)),~ismember(i,Safe_pat_xls)];
%     patientData(i,47:56)=max(table2array(S(1:2,158:167)));
%     patientData(i,107:116)=max(table2array(S(pos:pos+1,158:167)));
%     patientData(i,57:60)=max(table2array(S(1:2,176:179)));
%     patientData(i,end-4:end-1)=max(table2array(S(pos:pos+1,176:179)));
    gender(i)=S.Sex(1);
%     nondimTNT(i)=log(max(S.VarName10)-min(S.VarName10))*max(S.VarName10);
end
    binaryGender=strcmp(gender,'Female');

% patientData(isnan(patientData)) = 5e4;
outcomeTNT=(max(patientData(:,1),patientData(:,61))>100)'; %110 ng/l TNT as the threshold for heartattack
outcomeHsTNT=(max(patientData(:,2),patientData(:,62))>34.2)'.*strcmp(gender,'Male')+(max(patientData(:,2),patientData(:,62))>15.4)'.*strcmp(gender,'Female'); %34.2 ng/l for male, 15.4 ng/l for female HsTNT based on uptodate.com
troponinOnly=[3:60,63:120]; %considering only troponin assays
HstroponinOnly=[3:60,63:120]; %considering only Hstroponin assays

% dummy=iszero(patientData(:,end));
% patientData(dummy,end)=2;
%% Processing of individual patient data-Exploratory analysis
%{
% B = mnrfit(patientData(:,tropninonly),patientData(:,end))
train_sample = randsample(118,100);
test_sample = find(~ismember(1:118,train_sample));
X = patientData(train_sample,troponinOnly);
% Initialize fitting parameters
initial_theta = zeros(length(troponinOnly),1);
y=outcomeTNT(train_sample)';

% Set regularization parameter lambda to 1 (you should vary this)
lambda = 1;

% Set Options
% options = optimset('Maxiter',3000);
[theta, J, exit_flag] = ...
	fminsearch(@(t)(costFunctionReg(t, X, y, lambda)), initial_theta);
prediction=predict(theta,patientData(test_sample,troponinOnly));

[prediction,outcomeTNT(test_sample)']
%}

%{
%% Initial plotting for problem description
Troponinexcel(isnan(table2array(Troponinexcel(:,complete_feature_set(:,1)))),:)=[];
for i=3:60
    figure
    mark=num2str(Troponinexcel.VarName12);
    
    subplot(2,1,1); % semi-log plot
    semilogy(table2array(Troponinexcel(:,complete_feature_set(:,i))),table2array(Troponinexcel(:,complete_feature_set(:,2))),'w*');
    ln = findobj('type','line');
    set(ln,'marker','.','markers',14,'markerfa','w')
    text(table2array(Troponinexcel(:,complete_feature_set(:,i))),(table2array(Troponinexcel(:,complete_feature_set(:,2)))),mark);
    
    subplot(2,1,2); % normal plot
    plot(table2array(Troponinexcel(:,complete_feature_set(:,i))),table2array(Troponinexcel(:,complete_feature_set(:,2))),'w*');
    ln = findobj('type','line');
    set(ln,'marker','.','markers',14,'markerfa','w')
    text(table2array(Troponinexcel(:,complete_feature_set(:,i))),(table2array(Troponinexcel(:,complete_feature_set(:,2)))),mark);
 
    %     semilogy(table2array(Troponinexcel(:,complete_feature_set(:,i))),table2array(Troponinexcel(:,complete_feature_set(:,1))));
    % ylim([0,1e4]);
    xlabel(Troponinexcel.Properties.VariableNames(complete_feature_set(:,i)));
    ylabel(Troponinexcel.Properties.VariableNames(complete_feature_set(:,2)));
end
for i=3:60
    mdlTNT=fitlm(table2array(Troponinexcel(:,complete_feature_set(:,i))),log(table2array(Troponinexcel(:,complete_feature_set(:,2)))),'quadratic');
    coeffs(:,i-2)=mdlTNT.Coefficients.Estimate;
    rsq(i-2,:)=[mdlTNT.Rsquared.Ordinary mdlTNT.Rsquared.Adjusted];
end
%}
%% Logistic regression using MATLAB inbuilt functions
X = [patientData(:,1),patientData(:,troponinOnly)];
Y=outcomeTNT(:)';
maxdev = chi2inv(.95,1);
rand('twister',0);
opt = statset('display','iter',...
    'TolFun',maxdev,...
    'TolTypeFun','abs');
% fs1 = featureIdxSortbyP(1:150);
holdoutCVP = cvpartition(Y,'holdout',10);
dataTrain = X(:,holdoutCVP.training)';
grpTrain = Y(holdoutCVP.training);
fivefoldCVP = cvpartition(grpTrain,'kfold',5);
inmodel = sequentialfs(@critfun,X,Y',...
    'cv','none',...
    'nullmodel',true,...
    'options',opt,...
    'direction','forward','nfeatures',8);
newFeatures=find(inmodel);
Xnew=patientData(:,newFeatures);
b = glmfit(Xnew,Y','binomial');

%% Stepwise regression using TNT/HsTNT values (1,61-TNT; 2,62-HsTNT)
X = [patientData(:,HstroponinOnly),binaryGender'];
Y=patientData(:,2);
idx=find(patientData(:,2)>2000);
patientData(idx,:)=[];
idx1=isnan(patientData(:,end-1));
patientData(idx1,:)=[];
[b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,Y);

%% Classification learner using y=max(TNT) and x=max(CBC)-min(CBC)

for i=1:length(a)-1
   
    filename=strcat('pat_data',num2str(i),'.mat');
    load (filename)
    if isnan(max(S.VarName11))==1
        continue
    end
    TNT(i)=max(binTNT);
    modTNT(i)=log(max(S.VarName10)-min(S.VarName10))+max(S.VarName10);
    minidx=find(S.VarName10==min(S.VarName10),1);
    maxidx=find(S.VarName10==max(S.VarName10),1);
%     minCBC=min(table2array(S(:,complete_feature_set)));
%     maxCBC=max(table2array(S(:,complete_feature_set)));
minCBC=table2array(S(minidx,complete_feature_set));
maxCBC=table2array(S(maxidx,complete_feature_set));
%     predictors(i,:)=sign(maxCBC-minCBC).*log(abs((maxCBC-minCBC)))+maxCBC;
predictors(i,:)=(maxCBC-minCBC);
end
binoutcome=(TNT>100);
reducedPredictors=[3,5,6,9,12,13,15,17,21,24,25,27,30,31,33,34,35,37,...
    40,43,44,45,46,48,49,51,52,53,54,55,56,57];
reducedPredictorsnew=[3,4,5,7,8,11,14,15,17,19,20,23,24,26,27,29,32,33,34,35,36,37,42,44,45,...
    46,48,49,50,51,52,53,54,56,57,58];
Xpred=[predictors(:,reducedPredictorsnew),binaryGender'];
yHSTNT=(TNT>34.2).*strcmp(gender,'Male')+(TNT>15.4).*strcmp(gender,'Female');
yTNT=TNT>100;
data=[Xpred,TNT'];
idx=find(modTNT>2000);
data(idx,:)=[];
idx1=isnan(data(:,end));
data(idx1,:)=[];
idx2=isinf(data(:,end));
data(idx2,end)=0;

%% identify correlation between predictors
xx1=data(:,1:10);
xx2=data(:,1:10);
plotmatrix(xx1,xx2);

%% new
iii=find(isnan(Troponinexcel.VarName10));
aaa=Troponinexcel;
aaa(iii,:)=[];
colsInpredictor=[36,113,115,116,117,120,121,122,123 259];
dataa=[aaa(:,colsInpredictor),aaa(:,end-1)];

%% logistic regression
% rand('twister',0);
train_sample = randsample(281,200);
test_sample = find(~ismember(1:281,train_sample));
% colsInpredictor=[36,113,115,116,117, 121 123 259];
X = table2array(aaa(train_sample,colsInpredictor));
% Initialize fitting parameters
initial_theta = zeros(length(colsInpredictor),1);
y=table2array(aaa(train_sample,end-1));

[bstep,se,pval,inmodel,statsstep,nextstep,history] = stepwisefit(table2array(aaa(:,colsInpredictor)),...
    table2array(aaa(:,end-1)));
% Set regularization parameter lambda to 1 (you should vary this)
lambda = 1000;

% Set Options
 options = optimset('Maxiter',3000);
[theta, J, exit_flag] = ...
	fminsearch(@(t)(costFunctionReg(t, X, y, lambda)), initial_theta);
prediction=predict(theta,table2array(aaa(test_sample,colsInpredictor)));
sum(abs((prediction-table2array(aaa(test_sample,end-1)))))
% [prediction,table2array(aaa(test_sample,end-2))]

%  inbuilt functions from MATLAB
for j = 1:5
    indices = crossvalind('Kfold',y,5);
    for i = 1:5
        test = (indices == i); train = ~test;
        [blog(:,j,i),dev,stats] = glmfit(X(train,:),y(train),'binomial','logit'); % Logistic regression
%         Fit = glmval(b,X(test,:),'logit');
%          filename=strcat('fit',num2str(i),num2str(j),'.mat');
%     save(filename,'stats')
    end
end