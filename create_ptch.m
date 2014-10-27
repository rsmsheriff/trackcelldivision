function ptch = create_ptch(img,pimg)
%Creates Patch parameters for the given image(img) with the patch image
%(pimg)
tptch = regionprops(pimg,img,'MeanIntensity','Centroid','Area');

vec = 1:size(tptch,1);

nofa = max(pimg(:));

ptch(1:nofa,1) = (1:nofa)';



ptch(vec,2) = [tptch.Area]';
ptch(vec,4) = [tptch.MeanIntensity]';

ptch(:,3) = ptch(:,2).*ptch(:,4);




ptch(vec,[5,6]) = reshape([tptch.Centroid],2,[])';  %This is not intensity weighted Centroid

% iFA = size(ptch,1);
% 
% [m,n] = size(img);
% 
% timg = img .* repmat(1:n ,m,1);
% 
% cmX_FA =  cell2mat(arrayfun(@(x) sum(timg(pimg==x))/ptch(x,3), 1:iFA,'UniformOutput', false)');
% 
% ptch(:,5) = cmX_FA;
%  
% timg = img .* repmat((1:m)',1,n);
% 
% cmY_FA =  cell2mat(arrayfun(@(x) sum(timg(pimg==x))/ptch(x,3), 1:iFA,'UniformOutput', false)');
% 
% ptch(:,6) = cmY_FA;


% 
% 
%  cmX_FA =  cell2mat(arrayfun(@(x) sum((tablePxls(find(tablePxls(:,3)==indFA(x)),1).*tablePxls(find(tablePxls(:,3)==indFA(x)),17)))/sum(tablePxls(find(tablePxls(:,3)==indFA(x)),17)), ...
%      [1:iFA],'UniformOutput', false)');


ptch(vec,7) = 1;


end