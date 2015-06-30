clear 
close all
clc
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/mghPatData.mat');

%Data Smoothing
Label = [cellstr('TNI vs WBC'), cellstr('TNI vs RBC_i'),cellstr('TNI vs HGB'),cellstr('TNI vs MCV'),cellstr('TNI vs RDW'),...
    cellstr('TNI vs MCH'),cellstr('TNI vs MCHC'),cellstr('TNI vs HCT'),cellstr('TNI vs PLT')];

x=data(:,1:numel(data(1,:))-1);
y=data(:,end);


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
    [rCBC,cp] = corrcoef(xnorm,ynorm);
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