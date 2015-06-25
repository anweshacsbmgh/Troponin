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
    
    if isempty(find(tMRN.MRN(i)==mghPatients.mrn, 1))==0
        idxmrn = find(tMRN.MRN(i)==mghPatients.mrn);
        match(1:length(find(abs(days(datetime(dtt(i,:))-datetime(dtm(idxmrn,:))))<2)),i) = idxmrn(1)-1+find(abs(days(datetime(dtt(i,:))-datetime(dtm(idxmrn,:))))<2);
        % if isempty(find(tMRN.MRN(i)==mghPatients.mrn(idxall), 1))==0
        %         match(1:length(find(tMRN.MRN(i)==mghPatients.mrn(idxall))),i) = idxall(tMRN.MRN(i)==mghPatients.mrn(idxall));
        
    else
        match(1,i)=0;
    end
end
%% matching accession number
clear tMRN
load '/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/tMRN_Accession.mat'
matchid = find(match(1,:));
match(match==0)=1;
accession = tMRN.sAccession;
mghaccession = char(mghPatients.accession);
mghaccession(:,1:4) =[];
for i = 1:length(tMRN.sAccession)
    if isempty(accession{i})==0
        idxmrn = find(tMRN.MRN(i)==mghPatients.mrn);
        if ismember(accession{i},cellstr(mghaccession(idxmrn,:)))==1
            idx(i)=idxmrn(1)-1+find(strcmp(accession{i},cellstr(mghaccession(idxmrn,:))));
        else
            idx(i)=0;
        end
    else
        idx(i)=0;
        
    end
end

%% putting together the list
list(2,:)=[mghPatients(idx(2),:),tMRN(2,:)];

for i = 3:length(idx)
    if idx(i)>0
        list(end+1,:) = [mghPatients(idx(i),:),tMRN(i,:)];
    end
end
list(1,:)=[];
save('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/list.mat','list')
%% Reading from the Datalog extracts
clear
load ('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/list.mat','list')
fnames = dir('/Volumes/MGH-CSB/higgins/data/sapphire/DataLogExtracts-42318az96');
fnames(1:2,:)=[];
string = repmat('%s ', 1, 263);
for i=1:length(fnames)
    filename = strcat('/Volumes/MGH-CSB/higgins/data/sapphire/DataLogExtracts-42318az96/',fnames(i).name);
    fid = fopen(filename);
    excel = textscan(fid,string,2,'delimiter',',');
    exceldate(i,:) = excel{1,2}(2);
    fclose(fid);
    %     clear excel
end
allFiles=struct2table(fnames);
allFiles.date = datestr(exceldate);
allFiles.year = allFiles.date(:,8:11);
% allFile.ddDatnum = datenum(allFiles.date);
allFiles = sortrows(allFiles,{'year','datenum'});
importantColumns = [1,2,8,22,100,101,103,102,106,107,108,109];

for i = 1:height(list)
    k =find(datetime(cell2mat(list.lab_date(i,:)))<=datetime(allFiles.date(:,:)),1)-1;
    if k==0
        k=1;
    end
    
    filename = cell2mat(strcat('/Volumes/MGH-CSB/higgins/data/sapphire/DataLogExtracts-42318az96/',allFiles.name(k)));
    [status, result] = system( ['wc -l ', filename] );
    numlines = str2num(result(1,1:9));
    fid = fopen(filename);
    excel = textscan(fid,string,numlines,'delimiter',',');
    fclose(fid);
    dummy = excel{1,8};
    lineNum = find(strcmp(list.sContainerID(i),dummy),1);
    if isempty(lineNum)==1 
        k1=k+1;
        filename = cell2mat(strcat('/Volumes/MGH-CSB/higgins/data/sapphire/DataLogExtracts-42318az96/',allFiles.name(k1)));
        
        [status, result] = system( ['wc -l ', filename] );
        numlines = str2num(result(1,1:9));
        fid = fopen(filename);
        excel = textscan(fid,string,numlines,'delimiter',',');
        fclose(fid);
        dummy = excel{1,8};
        lineNum = find(strcmp(list.sContainerID(i),dummy),1);
    end
    
    if isempty(lineNum)==1 && k>1
        k1=k-1;
        filename = cell2mat(strcat('/Volumes/MGH-CSB/higgins/data/sapphire/DataLogExtracts-42318az96/',allFiles.name(k1)));
        
        [status, result] = system( ['wc -l ', filename] );
        numlines = str2num(result(1,1:9));
        fid = fopen(filename);
        excel = textscan(fid,string,numlines,'delimiter',',');
        fclose(fid);
        dummy = excel{1,8};
        lineNum = find(strcmp(list.sContainerID(i),dummy),1);
    end
    if isempty(lineNum)==1
            dummy1 = zeros(1,length(importantColumns));
            allReducedLogs (i,j) = mat2cell(dummy1,1);
    else
        
        for j=1:length(importantColumns)
            allReducedLogs (i,j) =  excel{1,importantColumns(j)}(lineNum);
        end
    end
    clear excel
    clear numlines
    clear lineNum
end

delIND = find(any(cellfun(@isempty,allReducedLogs(:,1)),2));
allReducedLogs(delIND,:)=[];
list(delIND,:) = [];
CompList = [list,cell2table(allReducedLogs)];

CompList.Properties.VariableNames{17}='ContainerID';
CompList.Properties.VariableNames{18}='WBC';
CompList.Properties.VariableNames{19}='RBC';
CompList.Properties.VariableNames{20}='HGB';
CompList.Properties.VariableNames{21}='RDW';
CompList.Properties.VariableNames{22}='MCV';
CompList.Properties.VariableNames{23}='MCH';
CompList.Properties.VariableNames{24}='MCHC';
CompList.Properties.VariableNames{25}='HCT';
CompList.Properties.VariableNames{26}='PLT';
findindi = find(strcmp(CompList.results,'<0.01')); %finding rows where result reads "<0.01"
CompList.results(findindi) = cellstr('0.005');

AllResults  = [CompList(:,18:26), CompList(:,5)];

ResultArray = table2array(AllResults);
ResultArray(strcmp(ResultArray(:,end),'Credit'),:)=[];
delIND1= find(any(cellfun(@isempty,ResultArray(:,end)),2));
ResultArray(delIND1,:)=[];
ResultArray(strcmp(ResultArray(:,end),'Refused'),:)=[];
ResultArray(strcmp(ResultArray(:,end),'Cancelled'),:)=[];

data = cellfun(@str2num,ResultArray);
save('/Users/anweshachaudhury/Desktop/mghTN/mghTNpatients/mghPatData.mat');
% % string = ['%d ', '%D ', '%s ', '%s ', '%d ', '%d ', '%d ', '%s ',repmat('%d ', 1, 255)];
% string = repmat('%s ', 1, 263);
% % for i = 1:254
% %     string = strcat(string, '%d ');
% % end
% fnames(1:2,:) =[];
%
% for i=1:length(fnames)
%     file = fnames(i).name;
%     filename = strcat('/Volumes/MGH-CSB/higgins/data/sapphire/DataLogExtracts-42318az96/',file);
%     [status, result] = system( ['wc -l ', filename] );
%     numlines = str2num(result(1,1:9));
%     fid = fopen(filename,'r');
% %     fnamesSub = dir('/Volumes/MGH-CSB/higgins/data/sapphire/file');
%     excel = textscan(fid,string,numlines,'delimiter',',');
%     allDatalogs(i,1:263) = excel;
%     clear excel
%     clear numlines
%     fclose(fid);
% end
% allReducedLogs = cell(1,length(importantColumns));
% importantColumns = [1,2,8,22,100,101,103,102,106,107,108,109];
% for i = 1:length(fnames)
%     for j=1:length(importantColumns)
%     allReducedLogs (end+1:end+numel(allDatalogs{i,1}),j) =  allDatalogs{i,importantColumns(j)};
%     end
% end
%

%% Searching for the corresponding filepath and the name of the FCS files
namesf =  dir('/Volumes/MGH-CSB/higgins/data/sapphire');
namesf([1:3,237:244],:)=[];
namesf = struct2table(namesf);
namesf = sortrows(namesf,'datenum','ascend');

fileDate = datestr(namesf.date(:));
for i = 1:length(ResultArray)
    i
indic(i) = find(datetime(cell2mat(CompList.allReducedLogs2(i,:)))<=datetime(cell2mat(namesf.date(:,:))),1);
comp = datestr(namesf.date(indic(i)));
allindic = find(strcmp(cellstr(comp(1:11)),cellstr(fileDate(:,1:11))));
for j = 1:length(allindic)
    nameSubF = strcat('/Volumes/MGH-CSB/higgins/data/sapphire/',namesf.name(allindic(j)));
    allsubF = struct2table(dir(cell2mat(nameSubF)));
names = char(allsubF.name);
dummyid = find(~ismember(cellstr(names(1:length(names),end-2:end)),'typ'));
datevec = [str2num(names(dummyid,11:14)),str2num(names(dummyid,15:16)),str2num(names(dummyid,17:18)),...
    str2num(names(dummyid,20:21)),str2num(names(dummyid,22:23)),str2num(names(dummyid,24:25))];
 rundate = datestr(datevec);   
 findFile = find(strcmp(cell2mat(CompList.allReducedLogs1(i)),cellstr(num2str(str2num(names(dummyid,27:30)))))); %strcmp(cellstr(rundate(:,1:11)),comp(1:11)) & 
 if isempty(findFile)==0
     assignFile(i,:) = names(findFile,:);
     pathFile(1,i) = nameSubF;
     clear findFile
%      keyboard
     break
 end
 
end
   
end
