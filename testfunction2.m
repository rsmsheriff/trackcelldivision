zmtch = mtch;
npimg1 = zeros(size(pimg1)); npimg2 = npimg1;
for ia = 1:size(zmtch,1), npimg1(pimg1==zmtch(ia,1)) =  ia;  npimg2(pimg2==zmtch(ia,2)) =  ia; end
figure(1), imagesc(npimg1), figure(2), imagesc(npimg2);
% function testfunction2(I)
% hy = fspecial('sobel');
% hx = hy';
% Iy = imfilter(double(I), hy, 'replicate');
% Ix = imfilter(double(I), hx, 'replicate');
% gradmag = sqrt(Ix.^2 + Iy.^2);
% figure, image(gradmag), title('Gradient magnitude (gradmag)')
% L = watershed(gradmag);
% Lrgb = label2rgb(L);
% figure, image(Lrgb), title('Watershed transform of gradient magnitude (Lrgb)')
% se = strel('disk', 20);
% Io = imopen(I, se);
% figure, image(Io), title('Opening (Io)')
% Ie = imerode(I, se);
% Iobr = imreconstruct(Ie, I);
% figure, image(Iobr), title('Opening-by-reconstruction (Iobr)')
% Ioc = imclose(Io, se);
% figure, image(Ioc), title('Opening-closing (Ioc)')
% Iobrd = imdilate(Iobr, se);
% Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
% Iobrcbr = imcomplement(Iobrcbr);
% figure, image(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
% fgm = imregionalmax(Iobrcbr);
% figure, image(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
% I2 = I;
% I2(fgm) = 255;
% figure, image(I2), title('Regional maxima superimposed on original image (I2)')
% se2 = strel(ones(5,5));
% fgm2 = imclose(fgm, se2);
% fgm3 = imerode(fgm2, se2);
% fgm4 = bwareaopen(fgm3, 20);
% I3 = I;
% I3(fgm4) = 255;
% figure, image(I3)
% title('Modified regional maxima superimposed on original image (fgm4)')
% bw = im2bw(Iobrcbr, graythresh(Iobrcbr));
% figure, image(bw), title('Thresholded opening-closing by reconstruction (bw)')
% D = bwdist(bw);
% DL = watershed(D);
% bgm = DL == 0;
% figure, image(bgm), title('Watershed ridge lines (bgm)')
% gradmag2 = imimposemin(gradmag, bgm | fgm4);
% L = watershed(gradmag2);
% I4 = I;
% I4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
% figure, image(I4)
% title('Markers and object boundaries superimposed on original image (I4)')
% Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
% figure, image(Lrgb)
% title('Colored watershed label matrix (Lrgb)')
% figure, image(I), hold on
% himage = image(Lrgb);
% set(himage, 'AlphaData', 0.3);
% title('Lrgb superimposed transparently on original image')