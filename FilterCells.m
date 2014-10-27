%Display all tracking results
lst = rdir('/Volumes/quantcellbio$/SHERIFF/sorted/*/**/TrackQuatYFP.mat');
slifet = [];
ThresholdLife = 60; %in frames
NucSizeChnB = 20; %in percentage the nucluer size change after formation (B- begin) in first 3 frames
NucSizChgThdB =  1 - NucSizeChnB/100; %i.e. the least change is 20% -> 0.8
NucSizeChnE = 20; %in percentage the nucluer size change before division (E -end) in last 3 frames
NucSizChgThdE =  1 + NucSizeChnE/100; %i.e. the least change is 20% -> 0.8

for ia = 1:size(lst,1)
     in(ia) = load(lst(ia).name);
    k = in(ia).qparms;
    k(k==0) = nan;
    cellArea = squeeze(k(:,4,:)); %Area -> 4th col
    sel = ~isnan(cellArea);
    
    lifet = sum(sel,2); %Lentth of cell
%     Nffcell = isnan(k(:,2,1)); %not in first frame
%     Nlfcell = isnan(k(:,2,end)); %not in first frame
%     lifet = lifet(Nffcell & Nlfcell); %select cells not in first and last frame
%     slifet = [slifet; lifet];
    
    %idx = Nffcell & Nlfcell & lifet>60;
    
    %select Cells above ThresholdLife
    cellArea = cellArea(lifet > ThresholdLife,:);
    sel = ~isnan(cellArea);
    
    %Looking at the area of the cell at the start and end of the cell
%     sel = squeeze(~isnan(k(:,4,:)));
    lstf = max(bsxfun(@times, sel~=0, 1:size(sel,2)).');
    [~,fstf] = max(sel,[],2);
    
    %Select cells whos area increase from 1st to 3rd frame of its birth
     idx1 = cell2mat(arrayfun(@(x)((cellArea(x,fstf(x)) ./ cellArea(x,fstf(x)+3)) < NucSizChgThdB),1:size(fstf,1),'UniformOutput',0));
%     idx1 = cell2mat(arrayfun(@(x)(cellArea(x,fstf(x)) < cellArea(x,fstf(x)+3)),1:size(fst,1),'UniformOutput',0));
    
    %Select cells whos area increased in the last 3 frames before its birth
    %i.e. NEB
%      idx2 = cell2mat(arrayfun(@(x)((cellArea(x,lstf(x)) >  NucSizChgThdE*cellArea(x,lstf(x)-3))),1:size(fstf,1),'UniformOutput',0));   
    idx2 = false(size(idx1));
     for ib = 1:size(cellArea,1)
         [pk,pkid] = max(diff(cellArea(ib,:),2));
         idx2(ib) = pkid > lstf(ib)-7;
     end
     
     
     
    idx = idx1 & idx2;
    
end