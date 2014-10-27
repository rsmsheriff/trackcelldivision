function qparms = quanttrack(fullpth,img_channel)

img_lst = dir(fullfile(fullpth,['*' img_channel '*.tif']));
%        cyc_lst=dir(fullfile(readpth,wells{w},sites{s},['*' cyc_channel '*.tif']));
%         fullpth=fullfile(readpth,wells{w},sites{s});
%segt = load(fullfile(fullpth,'Analysis','SegQuantDataH2B_Er6_WS_0.15.mat'));
pimg = imreadstack(fullfile(fullpth,'Analysis','Tracking75_Er12.tif'))-1;
mx = max(pimg(:));
qparms = nan(mx,9,size(pimg,3));
%% main for loop
for ia=1:length(img_lst)
    %% read images
    nuc = imread(fullfile(fullpth,img_lst(ia).name));
    opars = extractpar(nuc,pimg(:,:,ia));
    qparms(1:size(opars,1),:,ia) = opars;
    ia
end
save(fullfile(fullpth,'Analysis',['TrackQuatH2B75_Er6' img_channel]),'qparms');

function ptch = extractpar(img,pimg)
%Creates Patch parameters for the given image(img) with the patch image
%(pimg)
tptch = regionprops(pimg,double(img),'MeanIntensity','Centroid','Area','PixelValues','MinIntensity','MaxIntensity');
vec = 1:size(tptch,1);
nofa = max(pimg(:));

ptch = nan(nofa,9);
ptch(1:nofa,1) = (1:nofa)';
ptch(vec,[2,3]) = reshape([tptch.Centroid],2,[])';  %This is not intensity weighted Centroid
ptch(vec,4) = [tptch.Area]';
ptch(vec,5) = [tptch.MeanIntensity]';
ptch(:,6) = ptch(:,2).*ptch(:,4);

for ia = 1:nofa
    ptch(ia,7) = std(tptch(ia).PixelValues);
end
maskEmptyId = isnan(ptch(:,2));
[tptch(maskEmptyId).MinIntensity] = deal(nan);
[tptch(maskEmptyId).MaxIntensity] = deal(nan);

ptch(vec,8) = [tptch.MinIntensity]';
ptch(vec,9) = [tptch.MaxIntensity]';

% ptch(1:nofa,9) = tptch(1:nofa).MaxIntensity;
% iFA = size(ptch,1);
% [m,n] = size(img);
% timg = img .* repmat(1:n ,m,1);
% cmX_FA =  cell2mat(arrayfun(@(x) sum(timg(pimg==x))/ptch(x,3), 1:iFA,'UniformOutput', false)');
% ptch(:,5) = cmX_FA;
% timg = img .* repmat((1:m)',1,n);
% cmY_FA =  cell2mat(arrayfun(@(x) sum(timg(pimg==x))/ptch(x,3), 1:iFA,'UniformOutput', false)');
% ptch(:,6) = cmY_FA;
%ptch(vec,7) = 1;