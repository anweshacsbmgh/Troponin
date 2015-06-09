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

RBCParams = xlsread('TNT.xls');
% [flag,patRefRBC] = ismember(RBCParams(:,end),Troponinexcel.VarName14);
% 
% for i = 1:length(RBCParams(:,1))
%     k=find(patRefRBC(i)<a);
%     loc(i) = k(1)-2;
% end

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
    %     p(i,:) = polyfit(1:length(S.VarName11),S.VarName11',1);
    clear S
end
patientnum = Troponinexcel(:,12);

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
% [flag,binaryColNo2] = ismember('binhsTNT',Troponinexcel.Properties.VariableNames);
colsInpredictor=[36,113,115,116,117,120,121,122,123,binaryColNo1];
CBC=zeros(1,length(colsInpredictor));
TNT=0;
hsTNT=0;
rtcidx=0;
Index = 0;
patref1= 0;
for i=1:length(a)-1
    CBC(end+1:end+length(a(i)+3:a(i+1)+1),:)=table2array(Troponinexcel(a(i)+3:a(i+1)+1,colsInpredictor));
    hsTNT(end+1:end+length(a(i)+3:a(i+1)+1))=Troponinexcel.VarName11(a(i)+3:a(i+1)+1); %predictors are CBCs at max(HsTNT)
    TNT(end+1:end+length(a(i)+3:a(i+1)+1))=Troponinexcel.VarName10(a(i)+3:a(i+1)+1);
    rtcidx(end+1:end+length(a(i)+3:a(i+1)+1))=Troponinexcel.rtcTest(a(i)+3:a(i+1)+1);
    Index(end+1:end+length(a(i)+3:a(i+1)+1))=Troponinexcel.VarName14(a(i)+3:a(i+1)+1);
    patref1(end+1:end+length(a(i)+3:a(i+1)+1),:) = table2array(Troponinexcel(a(i)+3:a(i+1)+1,12));
end
patref1(1)=[];
patref=patref1;
CBC(1,:)=[];
hsTNT(1)=[];
TNT(1)=[];
rtcidx(1)=[];
Index(1)=[];
idxnan = find(isnan(TNT));
idxhsnan = find(isnan(hsTNT));
uniquePat = unique(patref);
%identifying the NaN for hsTNT because the excel had the entry ">5e4"
for i = 1:length(idxhsnan)
    notNaN(i) = idxhsnan(i).*~ismember(idxhsnan(i),idxnan);
end
hsTNT(notNaN(find(notNaN))) = 5e4;
%two CBC per TNT, so the TNT is assigned accordingly
for i=1:length(idxnan)
    TNT(idxnan(i))=TNT(idxnan(i)-1);
    hsTNT(idxnan(i))=hsTNT(idxnan(i)-1);
    CBC(idxnan(i),10) = CBC(idxnan(i)-1,10);
end

data = [CBC(find(rtcidx),:),hsTNT(find(rtcidx))']; %identifying the rows with rtc measurements and creating the dataset

patref = patref(find(rtcidx));
idd = find(hsTNT(find(rtcidx))>50000);
data(idd,:)=[];
patref(idd)=[];
dataCBC = zeros(1,numel(data(1,:)));
%identifying the first and following instances for each patient
%CBCs are future-current, TNT is current
for i = 1:length(uniquePat)
    firstidx(i) = find(uniquePat(i)==patref,1);
  lengthPat(i) = sum(ismember(patref,uniquePat(i)));
  for j=2:lengthPat(i)
  dataCBC(end+1,1:end-1) = data(firstidx(i)+j-1,1:end-1)-data(firstidx(i)+j-2,1:end-1);
  dataCBC(end,end) = data(firstidx(i)+j-2,end);
  end
end
dataCBC(1,:)=[];

for i = 1:numel(dataCBC(1,:))-2
    yy2 = smooth(dataCBC(:,i),dataCBC(:,end),0.1,'rlowess');
    yy21 = smooth(dataCBC(:,i),dataCBC(:,end),0.3,'rlowess');
    yy22=smooth(dataCBC(:,i),dataCBC(:,end),15,'moving');
    smoothY(i,:)=yy2;
    smoothY1(i,:)=yy21;
    smoothY2(i,:)=yy22;
        [bCBC(:,i),bint,r,rint,stats] = regress(yy2,[ones(size(dataCBC(:,i))) dataCBC(:,i)]);
p(i) = stats(3);
end
x=dataCBC(:,1:numel(dataCBC(1,:))-2);
y=dataCBC(:,end);
Label = [cellstr('hsTNT vs WBC'), cellstr('hsTNT vs RBC_i'),cellstr('hsTNT vs HGB'),cellstr('hsTNT vs MCV'),cellstr('hsTNT vs RDW'),...
    cellstr('hsTNT vs MCH'),cellstr('hsTNT vs MCHC'),cellstr('hsTNT vs HCT'),cellstr('hsTNT vs PLT')];
for i=1:numel(dataCBC(1,:))-2
    [xx,ind] = sort(x);
    subplot(3,3,i)
    plot(xx(:,i),y(ind),'w.',xx(:,i),smoothY(i,ind(:,i)),'r-',xx(:,i),smoothY1(i,ind(:,i)),'g-')
    ylim([0 200])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end



%% Smoothing with RBC parameters
for i = 1:length(RBCParams)
    rbcID(i) = ismember(RBCParams(i,end),Index);
    if rbcID(i)>0
        ID(i)=find(Index==RBCParams(i,end),1);
    else
        ID(i)=0;
    end
end

delRow = find(ID==0);
ID(delRow)=[];
RBCParams(delRow,:)=[];
dataRBC1 = [RBCParams(:,1:6),hsTNT(ID)'];
patrefRBC = patref1(ID);
NaNdelRows = find(isnan(RBCParams(:,1)));
dataRBC1(NaNdelRows,:)=[];
patrefRBC(NaNdelRows)=[];
RDW = CBC(ID,5);
RDW(NaNdelRows) = [];
iddRBC = find(dataRBC1(:,end)>50000);
dataRBC1(iddRBC,:)=[];
patrefRBC(iddRBC)=[];
uniquePatRBC = unique(patrefRBC);
[sortPat,patid]=sort(patrefRBC);
dataRBC = dataRBC1(patid,:); 
dataRBC_ = zeros(1,numel(dataRBC(1,:)));

clear lengthPat
clear firstidx
for i = 1:length(uniquePatRBC)
    firstidx(i) = find(uniquePatRBC(i)==sortPat,1);
    lengthPat(i) = sum(ismember(sortPat,uniquePatRBC(i)));
    for j=2:lengthPat(i)
        dataRBC_(end+1,1:end-1) = dataRBC(firstidx(i)+j-1,1:end-1)-dataRBC(firstidx(i)+j-2,1:end-1);
        dataRBC_(end,end) = dataRBC(firstidx(i)+j-2,end);
    end
end
dataRBC_(1,:)=[];

for i = 1:numel(dataRBC_(1,:))-1
    yy2 = smooth(dataRBC_(:,i),dataRBC_(:,end),0.1,'rlowess');
    yy21 = smooth(dataRBC_(:,i),dataRBC_(:,end),0.3,'rlowess');
    yy22=smooth(dataRBC_(:,i),dataRBC_(:,end),15,'moving');
    
    smoothYRBC(i,:)=yy2;
    smoothYRBC1(i,:)=yy21;
    smoothYRBC2(i,:)=yy22;
     [bRBC(:,i),bint,r,rint,stats] = regress(yy2,[ones(size(dataRBC_(:,i))) dataRBC_(:,i)]);
pRBC(i) = stats(3);
end
% figure
% plot(dataRBC(:,end-1),RDW,'.')
Label = [cellstr('hsTNT vs \alpha'), cellstr('hsTNT vs \beta_v'),cellstr('hsTNT vs \beta_h'),cellstr('hsTNT vs D_v'),cellstr('hsTNT vs D_h'),...
    cellstr('hsTNT vs v_c')];

x=dataRBC_(:,1:numel(dataRBC_(1,:))-1);
y=dataRBC_(:,end);
figure
for i=1:numel(dataRBC_(1,:))-1
    [xx,ind] = sort(x);
    subplot((numel(dataRBC_(1,:))-1)/2,2,i)
    plot(xx(:,i),y(ind),'w.',xx(:,i),smoothYRBC(i,ind(:,i)),'r-',xx(:,i),smoothYRBC1(i,ind(:,i)),'g-')
    ylim([0 300])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end


%% Smoothing with PLT distribution
pltinfo=xlsread('/Users/anweshachaudhury/Desktop/Anwesha research/Analyzers/ChangiCompareWBC.xlsx','PLTchangi');
for i = 1:length(pltinfo)
    pltID(i) = ismember(pltinfo(i,1),Index);
    if pltID(i)>0
        ID(i)=find(Index==pltinfo(i,1),1);
    else
        ID(i)=0;
    end
end

delRow = find(ID==0);
ID(delRow)=[];
pltinfo(delRow,:)=[];
patrefplt = patref1(ID);

dataplt1 = [pltinfo(:,[5:19,22]),hsTNT(ID)'];

iddplt = find(dataplt1(:,end)>50000);
dataplt1(iddplt,:)=[];
patrefplt(iddplt,:) = [];

uniquePatplt = unique(patrefplt);
[sortPat,patid]=sort(patrefplt);
dataplt = dataplt1(patid,:); 
dataplt_ = zeros(1,numel(dataplt(1,:)));

clear lengthPat
clear firstidx
for i = 1:length(uniquePatplt)
    firstidx(i) = find(uniquePatplt(i)==sortPat,1);
    lengthPat(i) = sum(ismember(sortPat,uniquePatplt(i)));
    for j=2:lengthPat(i)
        dataplt_(end+1,1:end-1) = dataplt(firstidx(i)+j-1,1:end-1)-dataplt(firstidx(i)+j-2,1:end-1);
        dataplt_(end,end) = dataplt(firstidx(i)+j-2,end);
    end
end
dataplt_(1,:)=[];


for i = 1:numel(dataplt_(1,:))-1
    yy2 = smooth(dataplt_(:,i),dataplt_(:,end),0.1,'rlowess');
    yy21 = smooth(dataplt_(:,i),dataplt_(:,end),0.3,'rlowess');
    yy22=smooth(dataplt_(:,i),dataplt_(:,end),15,'moving');
    
    smoothYplt(i,:)=yy2;
    smoothYplt1(i,:)=yy21;
    smoothYplt2(i,:)=yy22;
    
     [bPLT(:,i),bint,r,rint,stats] = regress(yy2,[ones(size(dataplt_(:,i))) dataplt_(:,i)]);
pPLT(i) = stats(3);
end
Label = [cellstr('hsTNT vs MPV'), cellstr('hsTNT vs \mu_1'),cellstr('hsTNT vs \mu_2'),cellstr('hsTNT vs \sigma_1'),cellstr('hsTNT vs \sigma_2'),cellstr('hsTNT vs Kurtosis_1'),...
    cellstr('hsTNT vs Kurtosis_2'),cellstr('hsTNT vs Skewness_1'),cellstr('hsTNT vs Skewness_2'),cellstr('hsTNT vs Mode_1'),cellstr('hsTNT vs Mode_2'),...
    cellstr('hsTNT vs Median_1'),cellstr('hsTNT vs Median_2'),cellstr('hsTNT vs Interquart_1'),cellstr('hsTNT vs Interquart_2'),cellstr('hsTNT vs corr_{12}')];
x=dataplt_(:,1:numel(dataplt_(1,:))-1);
y=dataplt_(:,end);
figure
for i=1:numel(dataplt_(1,:))-1
    [xx,ind] = sort(x);
    subplot((numel(dataplt_(1,:))-1)/4,4,i)
    plot(xx(:,i),y(ind),'w.',xx(:,i),smoothYplt(i,ind(:,i)),'r-',xx(:,i),smoothYplt1(i,ind(:,i)),'g-')
    ylim([0 300])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end

