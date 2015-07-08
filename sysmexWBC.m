clear
close all
clc
addpath('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/fcs_read')
fcsdata = fcsread('/Users/anweshachaudhury/Desktop/Anwesha research/MATLABfuns/fcs_read/[XE-5000^A2741][00-10][06142015_101428][L2700284135][DIFF].fcs');
figure
plot(fcsdata(:,1),fcsdata(:,2),'.');
figure
h1 = histfit(fcsdata(:,2),150,'kernel');
[pks,idx,w,p] = findpeaks(-h1(2).YData,'MinPeakProminence',2);
ghost1 = fcsdata(fcsdata(:,2)<=h1(2).XData(idx(1)),:);
rest1=fcsdata(fcsdata(:,2)>h1(2).XData(idx(1)),:);
figure
h2 = histfit(rest1(:,1),150,'kernel');
[pks1,idx1,w1,p1] = findpeaks(-h2(2).YData,'MinPeakProminence',2);
mononuc = rest1(rest1(:,1) <=h2(2).XData(idx1(1)),:);
middle = rest1(rest1(:,1) >h2(2).XData(idx1(1)) & rest1(:,1) <h2(2).XData(idx1(2)),:);
eosino = rest1(rest1(:,1) >=h2(2).XData(idx1(2)),:);
figure
h3 = histfit(mononuc(:,1),50,'kernel');
[pks2,idx2,w2,p2] = findpeaks(-h3(2).YData,'MinPeakProminence',2);
lymph = mononuc(mononuc(:,1)<h3(2).XData(idx2(1)),:);
moreanalysis = [middle;mononuc(mononuc(:,1)>=h3(2).XData(idx2(1)),:)];
arctanv = atand((moreanalysis(:,1))./(moreanalysis(:,2))); %slope of each point from origin
figure
h4 = histfit(arctanv,50,'kernel');
[pks3,idx3,w3,p3] = findpeaks(-h4(2).YData,'MinPeakProminence',3);
mono = moreanalysis(arctanv<=h4(2).XData(idx3(2)),:);
neutro = moreanalysis(arctanv>h4(2).XData(idx3(2)),:);

figure
h5 = histfit(mono(:,1),100,'kernel');
[pks4,idx4,w4,p4] = findpeaks(-h5(2).YData,'MinPeakProminence',2);
moreghost = mono(mono(:,1)>h5(2).XData(idx4(1)),:);
mono(mono(:,1)>h5(2).XData(idx4(1)),:)=[];
figure
h6 = histfit(eosino(:,2),50,'kernel');
[pks5,idx5,w5,p5] = findpeaks(-h6(2).YData,'MinPeakProminence',3);
moreghost1 = eosino(eosino(:,2)>=h6(2).XData(idx5(1)),:);
eosino(eosino(:,2)>=h6(2).XData(idx5(1)),:)=[];
ghost = [ghost1;moreghost;moreghost1];
% arctanneu = atand((neutro(:,1))./(neutro(:,2)));
% figure
% h7 = histogram(arctanneu,70);
% [pks6,idx6,w6,p6] = findpeaks(-h7.Values,'MinPeakProminence',5);
% eosino1 = neutro(arctanneu<h7.BinEdges(idx6(1)),:);
% neutro(arctanneu<h7.BinEdges(idx6(1)),:)=[];
% eosin = [eosino;eosino1];
figure
plot(ghost(:,1),ghost(:,2),'b.',lymph(:,1),lymph(:,2),'m.',mono(:,1),mono(:,2),'g.',...
neutro(:,1),neutro(:,2),'c.',eosino(:,1),eosino(:,2),'r.')