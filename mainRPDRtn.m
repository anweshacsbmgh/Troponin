clear
close all
clc
%% Matching MRN from RPDR with the Maps and ignoring all cases before 2012
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/mghTNJune2015.mat');
load('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/tCBCAbbottAll20150427.mat');
uniqueMRN = unique(tMRN.MRN);
mghPatients = sortrows(mghPatients,{'mrn','lab_date'},{'ascend','ascend'});
tMRN = sortrows(tMRN,{'MRN','dDate'},{'ascend','ascend'});
dtm(:,1:20)=datestr(mghPatients.lab_date);
delrows = str2num(dtm(:,8:11))<2012;
mghPatients(delrows,:)=[];
dtm(delrows,:)=[];
dtt(:,:) = datestr(tMRN.dDate);
for i =1:length(mghPatients.mrn)
    if isempty(find(str2num(mghPatients.mrn{i})==uniqueMRN))==0
        idxmghu(i) = find(str2num(mghPatients.mrn{i})==uniqueMRN);
    else
        idxmghu(i) = 0;
    end
end
mghPatients(idxmghu==0,:)=[];
dtm(idxmghu==0,:)=[];
idxmghu(idxmghu==0)=[];
for i =1:length(tMRN.MRN)
    if isnan(tMRN.MRN(i))==0
        idxut(i) = find(tMRN.MRN(i)==uniqueMRN);
    else
        idxut(i) =0;
    end
end

% diffdate = zeros(length(tMRN.MRN),20); %for mapping tMRN against mghPatients
%     uniqueMRN = unique(tMRN.MRN);
% count = 0;
% flagid=1;
% for i = 1:length(tMRN.MRN)
%     i
%     if tMRN.MRN(i) == uniqueMRN(flagid) && sum(str2num(cell2mat(mghPatients.mrn))==uniqueMRN(flagid))>0
%                     count = find(str2num(cell2mat(mghPatients.mrn))==uniqueMRN(flagid),1);
%
%         for j=1:sum(str2num(cell2mat(mghPatients.mrn))==uniqueMRN(flagid))
%         diffdate(j) = days(datetime(dtt(i,:))-datetime(dtm(count+j-1,:)));
%
%         end
%         if abs(min(diffdate)) < 3
%         match(i) = count-1+find(diffdate==min(diffdate));
%         keyboard
%         else
%             match(i)=0;
%         end
%     else
%         flagid = flagid+1;
%
%     end
%             clear diffdate
%
% end
match=zeros(10,length(tMRN.MRN));
mrns=str2num(cell2mat(mghPatients.mrn));
mghPatients.mrn = mrns;
for i = 1:length(tMRN.MRN)
    if isempty(find(tMRN.MRN(i)==mghPatients.mrn(idxall), 1))==0
        idxall = find(abs(days(datetime(dtt(i,:))-datetime(dtm(:,:))))<3);
        
        match(1:length(find(tMRN.MRN(i)==mghPatients.mrn(idxall))),i) = idxall(tMRN.MRN(i)==mghPatients.mrn(idxall));
    else
        match(1,i)=0;
    end
end