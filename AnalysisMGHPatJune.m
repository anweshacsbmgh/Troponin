clear
close all
clc
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/mghPatData.mat');

%Data Smoothing
Label = [cellstr('TNT vs WBC'), cellstr('TNT vs RBC_i'),cellstr('TNT vs HGB'),cellstr('TNT vs MCV'),cellstr('TNT vs RDW'),...
    cellstr('TNT vs MCH'),cellstr('TNT vs MCHC'),cellstr('TNT vs HCT'),cellstr('TNT vs PLT')];
data(data(:,end)<=0.09,:)=[];
for i = 1:numel(data(1,:))-1
    data(:,i) = (data(:,i) - mean(data(:,i)))./std(data(:,i));
end
x=data(:,1:numel(data(1,:))-1);
y=data(:,end);
interestingWBC = data(data(:,1)>=1,:);
restWBC = data(data(:,1)<1,:);
interestingPLT = data(data(:,end-1)<=-1,:);
restPLT = data(data(:,end-1)>-1,:);
for i = 1:numel(data(1,:))-1
    yy2 = smooth(data(:,i),data(:,end),0.04,'rlowess');
    yy21 = smooth(data(:,i),data(:,end),0.3,'rlowess');
    yy22=smooth(data(:,i),data(:,end),15,'moving');
    smoothY(i,:)=yy2;
    smoothY1(i,:)=yy21;
    smoothY2(i,:)=yy22;
    [bCBC(:,i),bint,r,rint,stats] = regress(yy2,[ones(size(data(:,i))) data(:,i)]);
    p(i) = stats(3);
    ynorm = (yy2-min(yy2))./(max(yy2)-min(yy2));
    xnorm = (data(:,i)-min(data(:,i)))./(max(data(:,i))-min(data(:,i)));
    [rCBC,cp] = corrcoef(data(:,i),yy2);
    corrpCBC(:,:,i)=cp;
    corrCBC(:,:,i) = rCBC;
end


for i=1:numel(data(1,:))-1
    [xx,ind] = sort(x);
    subplot(3,3,i)
    plot(xx(:,i),y(ind),'w.',xx(:,i),smoothY(i,ind(:,i)),'r-',xx(:,i),smoothY1(i,ind(:,i)),'g-')
    ylim([0 1])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end
figure
for i=1:numel(data(1,:))-1
    [xx,ind] = sort(x);
    subplot(3,3,i)
    plot(xx(:,i),y(ind),'b.',xx(:,i),smoothY(i,ind(:,i)),'r-',xx(:,i),smoothY1(i,ind(:,i)),'g-')
    ylim([0 1])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end
figure
boxplot(data)
% keyboard
clear smoothY
clear smoothY1
clear smoothY2

%% RBC parameter analysis
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/tDyn.mat')
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/list.mat')

[indi,ia,ic]=unique(tDyn.ix);
dataRBC = [(tDyn(ia,1:6)),cell2table(list.results(indi(1:length(ia)),:))];
delRowsRBC =  find(isnan(dataRBC.alpha));
findindi = find(strcmp(dataRBC.Var1,'<0.01')); %finding rows where result reads "<0.01"
dataRBC.Var1(findindi) = cellstr('0.005');
dataRBC (delRowsRBC,:)=[];

% ResultArray = dataRBC;
dataRBC(strcmp(dataRBC.Var1,'Credit'),:)=[];
delIND1= find(any(cellfun(@isempty,dataRBC.Var1),2));
dataRBC(delIND1,:)=[];
dataRBC(strcmp(dataRBC.Var1,'Refused'),:)=[];
dataRBC(strcmp(dataRBC.Var1,'Cancelled'),:)=[];
dataRBC.Var1=str2num(char(dataRBC.Var1));
dataRBC = table2array(dataRBC);

Label = [cellstr('TNT vs \alpha'), cellstr('TNT vs \beta_v'),cellstr('TNT vs \beta_h'),cellstr('TNT vs D_v'),cellstr('TNT vs D_h'),...
    cellstr('TNT vs v_c')];
for i = 1:numel(dataRBC(1,:))-1
    dataRBC(:,i) = (dataRBC(:,i) - mean(dataRBC(:,i)))./std(dataRBC(:,i));
end
dataRBC(dataRBC(:,end)<=0.09,:)=[];

x=dataRBC(:,1:numel(dataRBC(1,:))-1);
y=dataRBC(:,end);

figure
for i = 1:numel(dataRBC(1,:))-1
    yy2 = smooth(dataRBC(:,i),dataRBC(:,end),0.04,'rlowess');
    yy21 = smooth(dataRBC(:,i),dataRBC(:,end),0.3,'rlowess');
    yy22=smooth(dataRBC(:,i),dataRBC(:,end),15,'moving');
    smoothY(i,:)=yy2;
    smoothY1(i,:)=yy21;
    smoothY2(i,:)=yy22;
    [bRBC(:,i),bint,r,rint,stats] = regress(yy2,[ones(size(dataRBC(:,i))) dataRBC(:,i)]);
    p(i) = stats(3);
    ynorm = (yy2-min(yy2))./(max(yy2)-min(yy2));
    xnorm = (dataRBC(:,i)-min(dataRBC(:,i)))./(max(dataRBC(:,i))-min(dataRBC(:,i)));
    [rRBC,cp] = corrcoef(dataRBC(:,i),yy2);
    corrpRBC(:,:,i)=cp;
    corrRBC(:,:,i) = rRBC;
end


for i=1:numel(dataRBC(1,:))-1
    [xx,ind] = sort(x);
    subplot(3,2,i)
    plot(xx(:,i),y(ind),'w.',xx(:,i),smoothY(i,ind(:,i)),'r-',xx(:,i),smoothY1(i,ind(:,i)),'g-')
    ylim([0 1])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end
figure
for i=1:numel(dataRBC(1,:))-1
    [xx,ind] = sort(x);
    subplot(3,2,i)
    plot(xx(:,i),y(ind),'b.',xx(:,i),smoothY(i,ind(:,i)),'r-',xx(:,i),smoothY1(i,ind(:,i)),'g-')
    ylim([0 1])
    xlim([quantile(xx(:,i),0.005) quantile(xx(:,i),0.995)])
    title(Label(i))
end
figure
boxplot(dataRBC)
