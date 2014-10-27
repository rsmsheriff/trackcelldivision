function T=robustThreshold_value(im)
% im(im<=0) = NaN;
q=prctile(im(:),[5 95]);
indx=im>q(1) & im<q(2);
T=mean(im(indx))+2*std(im(indx));
