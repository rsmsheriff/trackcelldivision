function S = segmentimgV4(I_cropped)
I_eq = adapthisteq(I_cropped);
%imagesc(I_eq);
bw = robustThreshold(I_eq);
%imagesc(bw)
bw2 = imfill(bw,'holes');
bw3 = imopen(bw2, ones(5,5));
bw4 = bwareaopen(bw3, 40);
L1 = bwlabel(bw4);
[un, fr]= uniquencount(L1(:));

u = un(fr < 3000);
if(isequal(u,un(2:end)))
    S = L1;
else
bw5 = bw4;
nL1 = zeros(size(L1));
for ia = 1:length(u)
    idx = L1 == u(ia);
    bw5(idx) =  0;
    nL1(idx) = ia;
end
% bw5_perim = bwperim(bw5);
% overlay1 = imoverlay(I_eq, bw5_perim, [.3 1 .3]);
% imagesc(overlay1)
%imagesc(bw5);
mask_em = imextendedmax(I_eq, 0.15);
% imagesc(mask_em)
mask_em = imclose(mask_em, ones(5,5));
mask_em = imfill(mask_em, 'holes');
mask_em = bwareaopen(mask_em, 40);
% overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
% imagesc(overlay2)
I_eq_c = imcomplement(I_eq);  
I_mod = imimposemin(I_eq_c, ~bw5 | mask_em);
L2 = watershed(I_mod);
% imagesc(label2rgb(L))
L2 = double(L2);
md = mode(L2(:));
L2(L2 == md) = 0;
L2 = imdilate(L2,ones(3));
% % bw5 = double(bw4) - double(imdilate(L>=1,ones(3)));
bw6 = double(bw5) - double(L2>=1);
%bw6 = bwareaopen(bw6, 40);
L3 = bwlabel(bw6);
% 
L4 = L2 + (L3 + double(max(L2(:)))* double(bw6>0));

S = nL1 + (L4 + double(max(nL1(:))) * double(L4>0));
stats = regionprops(S,'Area');
%  maskEmptyId = arrayfun(  @(a)isempty(a.Area),  stats  );
% [stats(maskEmptyId).Area] = deal(nan);
S = double(S) .* double(ismember(S,  find([stats.Area]> 100)));
end
%imagesc(S);