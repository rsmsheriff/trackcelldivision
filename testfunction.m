clear zmtch;
npimg1 = zeros(size(pimg1)); npimg2 = npimg1;
if ~isempty(nlst1)
zmtch(1:length(nlst1),1) = nlst1;
for ia = 1:size(zmtch,1), npimg1(pimg1==zmtch(ia,1)) =  ia;  end
end
if ~isempty(nlst2)
zmtch(1:length(nlst2),2) = nlst2;
for ia = 1:size(zmtch,1), npimg2(pimg2==zmtch(ia,2)) =  ia; end
end
%zmtch = mtch;


figure(1), imagesc(npimg1), figure(2), imagesc(npimg2);

% function testfunction(I_cropped)
% I_eq = adapthisteq(I_cropped);
% %imagesc(I_eq);
% bw = robustThreshold(I_eq);
% %imagesc(bw)
% bw2 = imfill(bw,'holes');
% bw3 = imopen(bw2, ones(5,5));
% bw4 = bwareaopen(bw3, 40);
% bw4_perim = bwperim(bw4);
% overlay1 = imoverlay(I_eq, bw4_perim, [.3 1 .3]);
% %imagesc(overlay1)
% mask_em = imextendedmax(I_eq, 0.15);
% %imagesc(mask_em)
% mask_em = imclose(mask_em, ones(5,5));
% mask_em = imfill(mask_em, 'holes');
% mask_em = bwareaopen(mask_em, 40);
% overlay2 = imoverlay(I_eq, bw4_perim | mask_em, [.3 1 .3]);
% %imagesc(overlay2)
% I_eq_c = imcomplement(I_eq);
% I_mod = imimposemin(I_eq_c, ~bw4 | mask_em);
% L = watershed(I_mod);
% imagesc(label2rgb(L))