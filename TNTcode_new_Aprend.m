clear
close all
load('tropdata.mat');
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
Troponinexcel(345:346,:)=[]; %because S/n is identical to 343,344
Troponinexcel(331,:)=[];
Troponinexcel(328,:)=[];
Troponinexcel(325,:)=[];
Troponinexcel(321,:)=[];
Troponinexcel(304:309,:)=[];
Troponinexcel(300,:)=[];
Troponinexcel(269,:)=[];
Troponinexcel(250:251,:)=[];
Troponinexcel(197,:)=[];
Troponinexcel(176:177,:)=[]; %because S/n is identical to 72,73
Troponinexcel(142:145,:)=[];
Troponinexcel(81:82,:)=[];
Troponinexcel(58:62,:)=[];
Troponinexcel(23:24,:)=[];
Troponinexcel(19:20,:)=[]; %because S/n is identical to 22,23
numOfRows=size(Troponinexcel(:,1));
patientnum = Troponinexcel(:,12);
a = [-2; find(~isnan(table2array(patientnum)))]';
count=1;

RBCParams = xlsread('tnt.xls');
[flag,patRef] = ismember(RBCParams(:,end),Troponinexcel.VarName14);

for i = 1:length(RBCParams(:,1))
    k=find(patRef(i)<a);
    loc(i) = k(1)-2;
end

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
Safe_pat_xls = [12 13 18 22 23 29 31 33 37 40 41 54 55 56 59 60 61 68 71 75 77 79 82 88 89 93 95 100 101 104 108]'; %green markers on excel sheet

binout=0;
binhsout=0;
xlsbinout=0;
gend=0;
% figure
for i=1:length(a)-1
    Troponinexcel.VarName12(a(i)+3:a(i+1)+1)=i;
    S=Troponinexcel(a(i)+3:a(i+1)+1,:);
    S.bintnt=max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,10)))>10.*ones(size(Troponinexcel(a(i)+3:a(i+1)+1,10)));
    S.binhstnt=(max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,11)))>34.2).*strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Male')+(max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,11)))>15.4).*strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Female');
    maxidx(i)=find(S.VarName10==max(S.VarName10),1);
    %     [flag,idd(i)]=ismember(max(S.VarName10),Troponinexcel.VarName10);
    binout(1,end+1:end+length(S.bintnt)+1)=[S.bintnt;0];
    xlsbinout(1,end+1:end+length(S.bintnt)+1)=[ismember(i,Safe_pat_xls)*ones(length(S.bintnt),1);0];
    binhsout(1,end+1:end+length(S.binhstnt)+1)=[S.binhstnt;0];
    gend(1,end+1:end+length(S.Sex)+1)=[strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Female');0];
    %%%%%%%Plotting TNT to see trends for each patient
%     mark=num2str(i);
    S(find(isnan(S.VarName10)),:)=[];
    S.VarName11(find(isnan(S.VarName11)))=5e4;
%     plot(1:length(S.VarName11),S.VarName11,'k');
%     ln = findobj('type','line');
%     set(ln,'marker','.','markers',16)
%     text(1:length(S.VarName11),S.VarName11,mark);
%     hold on
    %%%%%%%%%%%%%%
    p(i,:) = polyfit(1:length(S.VarName11),S.VarName11',1);
    clear S
end
binout(1)=[];
binhsout(1)=[];
gend(1)=[];
xlsbinout(1)=[];
binout(end)=[];
binhsout(end)=[];
gend(end)=[];
xlsbinout(end)=[];

%
Troponinexcel.binTNT=binout';
Troponinexcel.binhsTNT=xlsbinout';
Troponinexcel.binGender=gend';

[flag,binaryColNo1] = ismember('binGender',Troponinexcel.Properties.VariableNames);
[flag,binaryColNo2] = ismember('binhsTNT',Troponinexcel.Properties.VariableNames);
colsInpredictor=[36,113,115,116,117,120,121,122,123,binaryColNo1,binaryColNo2];

for i=1:length(a)-1
    maxCBC(i,:)=Troponinexcel(a(i)+2+maxidx(i),:);
    maxTNT(i)=Troponinexcel.VarName11(a(i)+2+maxidx(i)); %predictors are CBCs at max(HsTNT)
    %     aaa(i,:)=Troponinexcel(a(i)+3,:);
    aaa(i,:)=Troponinexcel(a(i)+3,:); %predictors are CBC_0s
end
maxTNT(find(isnan(maxTNT)))=5e4; %replacing NaN's with 5e4
data=aaa(:,colsInpredictor(1:end));

ContTNT=aaa(:,colsInpredictor(1:end-1));
ContTNT.contTNT=(maxTNT)';
subContTNT=ContTNT;
ContTNT(find(maxTNT>500),:)=[]; %removing very high HsTNT values for continuous values
% safe patients are marked based on the color coding in the excel sheet

%% logistic regression and stepwise for binary

[bstep,se,pval,inmodel,statsstep,nextstep,history] = stepwisefit(table2array(data(:,1:end-1)),...
    table2array(data(:,end)));

%{
    % rand('twister',0);
train_sample = randsample(118,100);
test_sample = find(~ismember(1:118,train_sample));
X = table2array(data(train_sample,1:end-1));
% Initialize fitting parameters
initial_theta = zeros(length(colsInpredictor)-1,1);
y=table2array(data(train_sample,end));
% Set regularization parameter lambda to 1 (you should vary this)
lambda = 10;

% Set Options
 options = optimset('Maxiter',3000);
[theta, J, exit_flag] = ...
	fminsearch(@(t)(costFunctionReg(t, X, y, lambda)), initial_theta);
prediction=predict(theta,table2array(data(test_sample,1:end-1)));
difference = sum(abs((prediction-table2array(data(test_sample,end)))))
[prediction,table2array(data(test_sample,end))]
%}

%  k-fold crossvalidation using classification
x=table2array(data(:,1:end-1));
y=table2array(data(:,end));
cv = cvpartition(y,'k',10);
classf = @(XTRAIN, ytrain,XTEST)(classify(XTEST,XTRAIN,...
ytrain,'mahalanobis'));

cvMCR = crossval('mcr',x,y,'predfun',classf,'partition',cv)

% Decision tree
view(fitctree(x,y),'Mode','graph')

% k-fold crossvalidation using a different approach
indices = crossvalind('Kfold',y,10);
cp = classperf(y);
for i = 1:10
    test = (indices == i); train = ~test;
    class = classify(x(test,:),x(train,:),y(train,:),'mahalanobis');
    classperf(cp,class,test);
    Err(i)=cp.ErrorRate;
end
[B,FitInfo] = lassoglm(x,y,'binomial','CV',10);
% Training set
cvx = cvpartition(height(data),'Holdout',0.2);
Xtrain = x(training(cvx),1:end-1);
Ytrain = y(training(cvx),:);
% Test set
Xtest = x(cvx.test,1:end-1);
Ytest = y(cvx.test,:);
 glm = fitglm(Xtrain,Ytrain,'Distribution','binomial');
 grouphat=(glm.Fitted.Probability>0.7)+1;
group=Ytrain+1;
C = confusionmat(group,grouphat)
y_glm=glm.predict(Xtest)>0.7;
Ctest = confusionmat(y_glm+1,Ytest+1)
%% Regression using continuous values for hsTNT

[bstepc,sec,pvalc,inmodelc,statsstepc,nextstepc,historyc]=stepwisefit(table2array(ContTNT(:,1:end-1)),table2array(ContTNT(:,end)));
mdl1=stepwiselm(table2array(ContTNT(:,1:end-1)),...
    table2array(ContTNT(:,end)));

ypred=predict(mdl1,table2array(ContTNT(:,1:end-1)));
figure
plot(table2array(ContTNT(:,end)),ypred,'o',1:2000,1:2000);

[B,FitInfo] = lasso(table2array(ContTNT(:,1:end-1)),table2array(ContTNT(:,end)),'Alpha',1,'CV',10);
lassoPlot(B,FitInfo,'PlotType','CV');
%% Analyzing HsTNT with respect to subsets based on increasing of decreasing HsTNT
posSlopes=find(p(:,1)>=0);
subsetData=data(posSlopes,:);
subsetData(find(ismember([2,30,110],posSlopes)),:)=[];

% Training set
cvx1 = cvpartition(height(subsetData),'Holdout',0.2);
Xtrain1 = table2array(subsetData(training(cvx1),1:end-1));
Ytrain1 = table2array(subsetData(training(cvx1),end));
% Test set
Xtest1 = table2array(subsetData(cvx1.test,1:end-1));
Ytest1 = table2array(subsetData(cvx1.test,end));
 glmsub = fitglm(Xtrain1,Ytrain1,'Distribution','binomial');
 grouphat1=(glmsub.Fitted.Probability>0.7)+1;
group1=Ytrain1+1;
Csub = confusionmat(group1,grouphat1)
y_glmsub=glmsub.predict(Xtest1)>0.7;
Csubtest = confusionmat(y_glmsub+1,Ytest1+1)
% for continuous HsTNT regression
subsetDatac=subContTNT(posSlopes,:);
subsetDatac(find(ismember([2,30,110],posSlopes)),:)=[];
subsetDatac(find(table2array(subsetDatac(:,end))>2000),:)=[];

[bstepc,sec,pvalc,inmodelc,statsstepc,nextstepc,historyc]=stepwisefit(table2array(subsetDatac(:,1:end-1)),table2array(subsetDatac(:,end)));
mdl1=stepwiselm(table2array(subsetDatac(:,1:end-1)),...
    table2array(subsetDatac(:,end)));

ypred=predict(mdl1,table2array(subsetDatac(:,1:end-1)));
figure
plot(table2array(subsetDatac(:,end)),ypred,'o',1:2000,1:2000);