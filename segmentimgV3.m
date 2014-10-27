function S = watershedSegV2(I_cropped)
I_eq = adapthisteq(I_cropped);
%imagesc(I_eq);
bw = robustThreshold(I_eq);
imagesc(bw)
bw2 = imfill(bw,'holes');
bw3 = imopen(bw2, ones(5,5));
bw4 = bwareaopen(bw3, 40);
% bw4_perim = bwperim(bw4);
% overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
% imagesc(overlay1)
% mask_em = imextendedmax(I_eq, 0.10);
% imagesc(mask_em)
% mask_em = imclose(mask_em, ones(5,5));
% mask_em = imfill(mask_em, 'holes');
% mask_em = bwareaopen(mask_em, 40);
% overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
% imagesc(overlay2)
% I_eq_c = imcomplement(I_eq);  
% I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);
% L = watershed(I_mod);
% imagesc(label2rgb(L))
% L = double(L);
% md = mode(L(:));
% L(L == md) = 0;
% L = imdilate(L,ones(3));
% % bw5 = double(bw4) - double(imdilate(L>=1,ones(3)));
% bw5 = double(bw4) - double(L>=1);
% K = bwlabel(bw5);
% 
% S = L + (K + double(max(L(:)))* double(bw5>0));

S = bwlabel(bw4);
imagesc(S);