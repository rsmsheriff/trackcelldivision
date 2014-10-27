function [bw,T]=robustThreshold(im)
% im(im<=0) = NaN;
q=prctile(im(:),[5 95]);
indx=im>q(1) & im<q(2);
T=mean(im(indx))+1*std(im(indx));
bw=im>T;