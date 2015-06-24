clear 
close all
clc
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/mghPatData.mat');

%Data Smoothing
Label = [cellstr('TNI vs WBC'), cellstr('TNI vs RBC_i'),cellstr('TNI vs HGB'),cellstr('TNI vs MCV'),cellstr('TNI vs RDW'),...
    cellstr('TNI vs MCH'),cellstr('TNI vs MCHC'),cellstr('TNI vs HCT'),cellstr('TNI vs PLT')];

x=data(:,1:numel(data(1,:))-1);
y=data(:,end);
for i=1:numel(data(1,:))-1
    [xx,ind] = sort(x);
    subplot(3,3,i)
    plot(xx(:,i),y(ind),'w.',xx(:,i),smoothY(i,ind(:,i)),'r-',xx(:,i),smoothY1(i,ind(:,i)),'g-')
    ylim([0 10])
    xlim([quantile(xx(:,i),0.05) quantile(xx(:,i),0.95)])
    title(Label(i))
end