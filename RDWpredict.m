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
bbb=Troponinexcel;
for i=1:length(a)-1
    Troponinexcel.VarName12(a(i)+3:a(i+1)+1)=i;
    S=Troponinexcel(a(i)+3:a(i+1)+1,:);
    S.bintnt=max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,10)))>100.*ones(size(Troponinexcel(a(i)+3:a(i+1)+1,10)));
    S.binhstnt=(max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,11)))>34.2).*strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Male')+(max(table2array(Troponinexcel(a(i)+3:a(i+1)+1,11)))>15.4).*strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Female');
    
    binout(1,end+1:end+length(S.bintnt)+1)=[S.bintnt;0];
    binhsout(1,end+1:end+length(S.binhstnt)+1)=[S.binhstnt;0];
    gend(1,end+1:end+length(S.Sex)+1)=[strcmp(Troponinexcel.Sex(a(i)+3:a(i+1)+1),'Female');0];
    aaa(i,:)=Troponinexcel(a(i)+3,:);
    bbb.VarName10(a(i)+3)=NaN;
end
bbb(isnan(bbb.VarName10),:)=[];
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
colsInpredictor=[36,113,115,116,120,121,122,123,10,117];

dataa=table2array(aaa(:,colsInpredictor));
Xpred=table2array(bbb(:,colsInpredictor(1:end-1)));
[b,se,pval,inmodel,stats,nextstep,history]=stepwisefit(dataa(:,1:end-1),dataa(:,end));
mdl1=stepwiselm(dataa(:,1:end-1),dataa(:,end));
[beta,Sigma,E,CovB,logL] = mvregress(dataa(:,1:end-1),dataa(:,end));
idx=find(pval<0.1);
% plotAdjustedResponse(mdl1,2)
% figure
% plotAdjustedResponse(mdl1,3)
% figure
% plotAdjustedResponse(mdl1,5)
% figure
% plotAdjustedResponse(mdl1,6)
% figure
% plotAdjustedResponse(mdl1,7)
% % figure
% % plotAdjustedResponse(mdl1,8)
ypred=predict(mdl1,Xpred);
figure
plot(table2array(bbb(:,117)),ypred,'o',8:26,8:26);
