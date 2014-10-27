%% Init
clear
close all
clc

%% parameters
RegularizationFactor=0.05;
weakThreshFactor=0.9;
cytoRingWidth=3;
cleanSegmentationFactor = strel('disk',5);
minCellSize = 100; % change this according to magnification used

%% USER INPUT
readpth='/Users/rsheriff/Documents/Projects/CellCycleStructure/Data/DATA Timing of cell cycle phases/MCF10A & RPE cells/16May14-RPE-MCF10A-timing_001/sorted';
writepth = readpth; % change if you want to write reasults in a separate temp folder

montagepth = 'montages';
trajpth = 'celltraj';
anExampleWell='C5';
sites={'P00001'}; % {'s1', 's2', 's3'} etc
nuc_channel='YFP';
% cyc_channel='w1';

%% "LOAD Images" where files are
%%pth='/Users/rwollman/Desktop/SilviaData';
%%pth='/Users/silviasantos/Desktop/Test_Matlab_script';
%%pth='/volumes/SS-18Aug09-Cdc25C-fine-intervals_Plate_12324/s2';

WellRow='ABCDEFGH0123456789';

% open file for output
% fid=fopen(fullfile(writepth,'Results.csv'),'w');
% fprintf(fid,'ImgFileName,Site,Well,X,Y,T,Cyc_cyto_int,cyc_nuc_int(r),cyc_nuc_int(r)./cyc_cyto_int(r),histonestd,avgInt,temp_id\n');

%% find all positions and wells
sites=dir(fullfile(readpth,anExampleWell,'P*'));
sites={sites.name};

wells=dir(fullfile(readpth));
wells={wells.name};
wells = wells(cellfun(@(w) all(ismember(w,WellRow)),wells));

%% check what an image size is
tmp_img_lst=dir(fullfile(readpth,wells{1},sites{1},['*' nuc_channel '*.tif']));
% %%info = imfinfo(fullfile(pth,'s1',wells(3).name,tmp_img_lst(1).name));
info = imfinfo(fullfile(readpth,wells{1},sites{1},tmp_img_lst(1).name));
sz=[info.Width info.Height];

%% fill in paramters that are based on image size
[mb,nb] = bestblk(sz,140);
[X,Y]=meshgrid(1:sz(2),1:sz(1));

%% main loop
for s=1:length(wells)

    %% find all wells and loop around them
    for w=1:length(sites)
        %% read directory for image names
        nuc_lst=dir(fullfile(readpth,wells{w},sites{s},['*' nuc_channel '*.tif']));
%        cyc_lst=dir(fullfile(readpth,wells{w},sites{s},['*' cyc_channel '*.tif']));
        fullpth=fullfile(readpth,wells{w},sites{s});
        t0=now;


        %% main for loop
        for i=1:length(nuc_lst)


            %% read images
            nuc=imread(fullfile(fullpth,nuc_lst(i).name));
%             cyc=imread(fullfile(fullpth,cyc_lst(i).name));

            %% rescale intensity
            nuc=mat2gray(nuc);
%             cyc=mat2gray(cyc);

            %% find nuc masks (Primary)
            % we use an "Adaptive" e.g. per block thresholding
            nuc_bw = blkproc(nuc,[mb,nb],@(im) robustThreshold(im));
            nuc_bw=imerode(nuc_bw,ones(3));
            nuc_bw=bwareaopen(nuc_bw,minCellSize);
            lbl=bwlabel(nuc_bw);
            figure(1029), subplot(2,2,1);imagesc(lbl);
            
            
            %% segment cyc using robust threshold - non adaptive
            cyc = nuc;
            cyc_bw=robustThreshold(cyc);
            cyc_bw=imopen(cyc_bw,cleanSegmentationFactor);
            cyc_bw=bwareaopen(cyc_bw,minCellSize);
            lbl=bwlabel(cyc_bw);    
            figure(1029), subplot(2,2,2);imagesc(lbl);
            
            cyc = nuc;
            fl = fspecial('average',100);
            cyc = cyc-imfilter(cyc, fl);
            cyc_bw=robustThreshold(cyc);
            cyc_bw=imopen(cyc_bw,cleanSegmentationFactor);
            cyc_bw=bwareaopen(cyc_bw,minCellSize);
            lbl=bwlabel(cyc_bw);    
            figure(1029), subplot(2,2,3);imagesc(lbl);
            
            cyc = nuc;
            fl = fspecial('average',100);
            cyc = cyc-imfilter(cyc, fl);
            figure(1029), subplot(2,2,4);
            testfunction(cyc);
            
            
% % 
% %            %% Secondary (cytoplasm) segmentation using
% %             % R Jones, AE Carpenter, P Golland (2005) Voronoi-Based Segmentation of Cells on Image Manifolds,
% %             % ICCV Workshop on Computer Vision for Biomedical Image Applications, pp.
% %             % 535-543 (CellProfiler)
% % %             cyc_bw = im2bw(cyc,graythresh(cyc));
% %             cell_bw = IdentifySecPropagateSubfunction(lbl,double(cyc),cyc_bw,RegularizationFactor);
% %             cell_lbl = bwlabel(cell_bw);
% % 
% %             %% "Tertiary" e.g. cytoplasm and nuc seperatly
% %             cyto_lbl = cell_lbl;
% %             cyto_lbl (imerode(nuc_bw,ones(cytoRingWidth))) = 0;
% % %             cell_lbl(~ismember(cell_lbl,cyto_lbl))=0;
% %             nuc_lbl = cell_lbl;
% %             nuc_lbl (~nuc_bw) = 0;
% % 
% %             %% measaurements
% %             [bla,bla,cyc_cyto_int] = grp2cell(cyc(:) ,cyto_lbl(:));
% %             [bla,bla,cyc_nuc_int ] = grp2cell(cyc(:) ,nuc_lbl(:));
% %             [bla,bla,x] = grp2cell(X(:) ,nuc_lbl(:));
% %             [bla,cid,y] = grp2cell(Y(:) ,nuc_lbl(:));
% %             histone = cellfun(@std,grp2cell(nuc(:),nuc_lbl(:)));
% %             avgInt=mean(cyc(:));
% %             Results = [   histone avgInt*ones(length(cyc_cyto_int),1)];
% % 
% %             fprintf('image: %g time: %s\n',i,datestr(now-t0,13));
% % 
% %             %% write segmented image & output file
% %            
% %             imwrite(uint16(cell_lbl),fullfile(fullpth,regexprep(cyc_lst(i).name,[cyc_channel '.TIF'],'segmented.png')));
% %             
% %             %% write to file 
% %             for r=1:length(cyc_cyto_int)
% %                 fprintf(fid,'%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',...
% %                             cyc_lst(i).name,...
% %                             sites{s},...
% %                             wells{w},...
% %                             x(r),y(r),i,...
% %                             cyc_cyto_int(r),...
% %                             cyc_nuc_int(r),...
% %                             cyc_nuc_int(r)./cyc_cyto_int(r),...
% %                             histone(r),...
% %                             avgInt,...
% %                             cid(r));
% %             end
        end %% end analysis loop
    end %well loop

end % end site loop
fclose(fid);
