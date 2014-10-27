%Display all tracking results
    lst = rdir('/Volumes/quantcellbio$/SHERIFF/sortedHeLa/*/**/TrackQuatYFP.mat');
slifet = [];scellArea = []; sCellMeanInt = []; sCellSTD = [];
ThresholdLife = 50; %in frames
NucSizeChnB = 20; %in percentage the nucluer size change after formation (B- begin) in first 3 frames
NucSizChgThdB =  1 - NucSizeChnB/100; %i.e. the least change is 20% -> 0.8
NucSizeChnE = 10; %in percentage the nucluer size change before division (E -end) in last 3 frames
NucSizChgThdE =  1 + NucSizeChnE/100; %i.e. the least change is 20% -> 0.8

for ia = 1:size(lst,1)
       in(ia) = load(lst(ia).name);
    k = in(ia).qparms;
    k(k==0) = nan;
    cellArea = squeeze(k(:,4,:)); %Area -> 4th col
    sel = ~isnan(cellArea);
    
    lifet = sum(sel,2); %Lentth of cell
%     slifet = [slifet; lifet];
    
    %idx = Nffcell & Nlfcell & lifet>60;
    
    %select Cells above ThresholdLife
    idx1 = lifet > ThresholdLife;
    cellArea = cellArea(idx1,:);
    k = k(idx1,:,:);
    
    
    %Looking at the area of the cell at the start and end of the cell
%     sel = squeeze(~isnan(k(:,4,:)));
    sel = ~isnan(cellArea);
    lstf = max(bsxfun(@times, sel~=0, 1:size(sel,2)).');
    [~,fstf] = max(sel,[],2);
    
    %Select cells whos area increased in the last 3 frames before cell
    %division  i.e. NEB
%      idx2 = cell2mat(arrayfun(@(x)((cellArea(x,lstf(x)) >  NucSizChgThdE*cellArea(x,lstf(x)-3))),1:size(fstf,1),'UniformOutput',0));   
    idx2 = false(size(sel,1),1);
     for ib = 1:size(cellArea,1)
         [pk,pkid] = max(diff(cellArea(ib,:),2));
         idx2(ib) = pkid > lstf(ib)-7;
     end
     
      %Select cells whos area increase from 1st to 3rd frame of its birth
     idx3 = cell2mat(arrayfun(@(x)((cellArea(x,fstf(x)) ./ cellArea(x,fstf(x)+3)) < NucSizChgThdB),1:size(fstf,1),'UniformOutput',0));
%     idx1 = cell2mat(arrayfun(@(x)(cellArea(x,fstf(x)) < cellArea(x,fstf(x)+3)),1:size(fst,1),'UniformOutput',0));
    
    if(any(idx2) & any(idx3))
    idx = idx2 & idx3'; 
    cellArea = cellArea(idx,:);
    k = k(idx,:,:);
    sel = ~isnan(cellArea);
    lifet = sum(sel,2); %Length of selected cell
    slifet = [slifet; lifet];
    scellArea = [scellArea; cellArea];
    end
    
    
    %Plot SD of the cells
    if(~isempty(k))
    CellMeanInt = squeeze(k(:,5,:)); if(size(CellMeanInt,2) == 1), CellMeanInt = CellMeanInt'; end
    CellSTD = squeeze(k(:,7,:));  if(size(CellSTD,2) == 1), CellSTD = CellSTD'; end
%     try
    sCellMeanInt = [sCellMeanInt; CellMeanInt];
%     catch
%         'wait'
%     end
    
    sCellSTD = [sCellSTD;CellSTD];
    end
end