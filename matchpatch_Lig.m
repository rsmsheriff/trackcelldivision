function [mtch, n, npimg1, npimg2] = matchpatch(pimg1, pimg2,ptch1,ptch2)
%code to match the patches between two images
% *------------------------------------------------------------------------
% *  Created by Sheriff
% *------------------------------------------------------------------------
% addpath('/Users/mrsheriff/Documents/MATLAB/Clustering/Clustering Interface/'); %for uniquencount.m
%remove the patches which are flagged zero or -1 from the pimg i.e remove
%those which are not with flag 1.
lst = ptch1(ptch1(:,7)~=1,1);
for ia = 1:size(lst)
    pimg1(pimg1 == lst(ia)) = 0;
end
lst = ptch2(ptch2(:,7)~=1,1);
for ia = 1:size(lst)
    pimg2(pimg2 == lst(ia)) = 0;
end
IntTol = [0.025 0.05 0.1 0.15 0.2 0.25 0.5 ];
AreaTol = [0.25 0.5];
maxdisp = [25 50 100 150 200];  %Maximum displacement of cell i.e.radius(in pixels) to search for neighboring focal adhesions
ligase = 6;
climg = logical(logical(pimg1).*logical(pimg2));
mimg1 = pimg1(climg);
mimg2 = pimg2(climg);
mtch = [ mimg1(:) , mimg2(:)];
% [tmtch, freq] = uniquencount(tmtch); %tmtch is total matching (all matches)
mtch = unique(mtch,'rows');
%finding st rict matched and mismatched patches
smtch = mtch;
for ia= 1:2
    [un,ind,~] = unique(sort(mtch(:,ia)));
    ind = ind - [ 0; ind(1:end-1)];    %ind is the frequency of the element
    un = un(ind >1);  %remove element which are repeated more than once
    [~,ind]= setdiff(smtch(:,ia),un);  %find index of non-repeated elements
    smtch = smtch(ind,:); %strict match is smatch
end
smtch1 = crosscheck(smtch,ptch1,ptch2,IntTol(2),AreaTol(2));
%smtch2 = crosscheck(smtch,ptch1,ptch2,ligase*IntTol(end),ligase*AreaTol(end),ligase); %match for ligase
%mtch = [smtch1; smtch2];
mtch = smtch1;
n = size(mtch,1); %n is the total number of strict matches

percentOverlap = 0.5;
tmtch = [pimg1(:), pimg2(:)];
% tmtch = setdiff(tmtch,mtch,'rows');
idx = bsxfun(@and, pimg1(:)>0 , pimg2(:)>0);
tmtch = tmtch(idx,:);
if tmtch
[umtch freq] = uniquencount(tmtch);
else 
    umtch = [];
    freq = [];
end

if mtch
    [rmtch idx] = setdiff(umtch,mtch,'rows');
else
    rmtch = umtch;
    idx = true(size(freq));
end
if(any(rmtch))
    freq = freq(idx);
    % for ia = 1:size(rmtch,1);
    overlap =  [ freq./ ptch1(rmtch(1:end,1),2) , freq./ ptch2(rmtch(1:end,2),2)];
    Bool_overlap = overlap> percentOverlap;
    if(any(Bool_overlap))
        idxoverlap = bsxfun(@or,Bool_overlap(:,1),Bool_overlap(:,2));
        overlapmtch = rmtch(idxoverlap,:);
        omtch = overlapmtch;
        for ia= 1:2
            [un, freq] = uniquencount(overlapmtch(:,ia));
            un = un(freq >2);  %remove element which are repeated more than once
            for ib = 1:size(un,1)
                omtch(omtch(:,ia)== un(ib),:) = [];  %repeating the number of the matched FA twice should check?????????????????
            end
            % [~,ind1]= setdiff(omtch(:,ia),un);  %find index of non-repeated elements
            % omtch = omtch(ind1,:); %non repeated ovelap match is omatch
        end
        omtch = crosscheck(omtch,ptch1,ptch2,IntTol(2),AreaTol(1));
        mtch = [mtch ; omtch];
    end
end

% testfunction2;
%finding lost patches
lmtch = [pimg1(:), pimg2(:)];  %lmtch  is lost matches i.e [x 0] or [0 x] patches
lmtch = unique(lmtch,'rows');
nlst1 = setdiff(lmtch(:,1),[mtch(:,1);0]);
nlst2 = setdiff(lmtch(:,2),[mtch(:,2);0]);


% testfunction;
for ia = 1:length(maxdisp)
    for ib = 1:length(IntTol) - (ia-1)
        nlmtch = nearestneighbourDI(nlst1,nlst2,ptch1,ptch2,maxdisp(ia),IntTol(ib));
        
        if nlmtch
            % nlmtch = crosscheck(nlmtch,ptch1,ptch2,2*IntTol,2*AreaTol);
            mtch = [mtch;nlmtch];
            nlst1 = setdiff(nlst1, [mtch(:,1);0]);
            nlst2 = setdiff(nlst2, [mtch(:,2);0]);
            
        end
    end
end

%Capture NEB matching
 nlmtch = nearestneighbourDI(nlst1,nlst2,ptch1,ptch2,maxdisp(3),2); 
 if nlmtch
     %nlmtch = crosscheck(nlmtch,ptch1,ptch2,ligase*IntTol(end),ligase*AreaTol(end),ligase);
     nlmtch = crosscheck(nlmtch,ptch1,ptch2,ligase*IntTol(end),ligase*AreaTol(end),ligase);
     mtch = [mtch;nlmtch];
     nlst1 = setdiff(nlst1, [mtch(:,1);0]);
     nlst2 = setdiff(nlst2, [mtch(:,2);0]);
     
 end
% testfunction;    
% nlmtch = nearestneighbourDI(nlst1,nlst2,ptch1,ptch2,maxdisp, IntTol*10);
% if nlmtch
% % nlmtch = crosscheck(nlmtch,ptch1,ptch2,5*IntTol,5*AreaTol,ligase);
% mtch = [mtch;nlmtch];
% nlst1 = setdiff(nlst1, mtch(:,1));
% nlst2 = setdiff(nlst2, mtch(:,2));
% end
% end
if nlst1
    mtch = [mtch;nlst1, zeros(size(nlst1))]; % figure(123), imagesc(ismember(pimg1,nlst1));
end
if nlst2
    mtch = [mtch; zeros(size(nlst2)),nlst2]; % figure(124), imagesc(ismember(pimg2,nlst2));
end
% else
%
%     if(nlst1)
%     nlmtch1 = [nlst1(:,1) , zeros(size(nlst1,1),1)];
%     mtch = [mtch;nlmtch1];  %should check if I should unique (nlmtch1,'rows');
%     end
%     if(nlst2)
%     nlmtch2 = [zeros(size(nlst2,1),1), nlst2(:,1) ];
% %     nlmtch = unique([ nlmtch1; nlmtch2],'rows');
%     mtch = [mtch;nlmtch2];
%     end
%
% end
if nargout == 4
    %make similar color to the FA in both
    npimg1 = zeros(size(pimg1));
    npimg2 = zeros(size(pimg2));
    for ia = 1:size(mtch,1)
        if(mtch(ia,1))
            npimg1(pimg1== mtch(ia,1)) = ia;
        end
        if(mtch(ia,2))
            npimg2(pimg2 == mtch(ia,2)) = ia;
        end
    end
    % pimg1 = npimg1;
    % pimg2 = npimg2;
    cmp = colormap(jet(512));
    cmp = cmp(101:512,:);
    cmp = cmp(randperm(412),:);
    cmp(1,:) = [ 1 1 1];
    figure(1253);
    image(pimg1);
    figure(1254)
    image(pimg2);
    figure(1),image(npimg1+1);
    for ia = 1: size(mtch,1)
        if(mtch(ia,1))
            text(round(ptch1( mtch(ia,1) ,6)),round(ptch1( mtch(ia,1),5)), num2str(ia), 'FontSize',8, 'Color','k');
        end
    end
    colormap(cmp);
    figure(2), image(npimg2+1);
    for ia = 1: size(mtch,1)
        if(mtch(ia,2))
            text(round(ptch2( mtch(ia,2) ,6)),round(ptch2( mtch(ia,2),5)), num2str(ia), 'FontSize',8, 'Color','k');
        end
    end
    colormap(cmp);
end
%Convert the mtch variable from double to cell array and put split/merges
%in single cell
% if(~exist('ischn','var'))
if any(mtch(:))
[un1,fr] = uniquencount(mtch(:,1));
un1 = un1(fr>1);
un1 = setdiff(un1,0);
[un2,fr] = uniquencount(mtch(:,2));
un2=  un2(fr>1);
un2 = setdiff(un2,0);
tmtch = mtch;
bolrm = false(size(mtch,1),1);
mtch = num2cell(mtch);
mtch(tmtch == 0) = {[]};
for ia = 1:length(un1)
    lst = find(tmtch(:,1) == un1(ia));
    mtch(lst(1),2) = {tmtch(lst,2)'};
    bolrm(lst(2:end)) = true;
end
for ia = 1:length(un2)
    lst = find(tmtch(:,2) == un2(ia));
    mtch(lst(1),1) = {tmtch(lst,1)'};
    bolrm(lst(2:end)) = true;
end
mtch(bolrm,:) = [];
else
    mtch = cell(1,2);
end
% else
%     mtch = num2cell(mtch);
% end
% function [ind mdist] = mindist(x,y, xv, yv)
% dist = sqrt(((x - xv).^2)   +  ((y - yv).^2));
% mdist = min(dist);
% ind = find(dist == mdist);

function mtch = crosscheck(mtch,ptch1,ptch2,IntTol, AreaTol,ligase)
%check if the integrated intensity matches as well
dIngInt = 2* (ptch1(mtch(:,1),3) - ptch2(mtch(:,2),3)) ./ (ptch1(mtch(:,1),3) + ptch2(mtch(:,2),3));
dArea = 2* (ptch1(mtch(:,1),2) - ptch2(mtch(:,2),2)) ./ (ptch1(mtch(:,1),2) + ptch2(mtch(:,2),2));
if nargin == 6
adsl = (dIngInt./ abs(dIngInt)) < 0 & (dArea ./ abs(dArea)) < 0;
else
    adsl = true(size(dIngInt,1),1);
end

sl1  = abs(dIngInt) < IntTol;
sl2  = abs(dArea) < AreaTol;
mtch = mtch( (sl1 & sl2) & adsl,:);
% function nlmtch = nearestneighbour(nlst1,nlst2,ptch1,ptch2,rad)
% %Find the matching patched based on the distance, scan through a defined
% %radius
% x1 = ptch1(nlst1,6);
% y1 = ptch1(nlst1,5);
% x2 = ptch2(nlst2,6);
% y2 = ptch2(nlst2,5);
% sz1 = size(x1,1);
% sz2 = size(x2,1);
% lmt1 = zeros(sz1,2);
% lmt2 = zeros(sz2,2);
% nlmtch = [];
% if(~isempty(x2) && ~isempty(x1))
%     for ia = 1:sz1
%         [ind mdist] = mindist(x1(ia),y1(ia),x2,y2);
%         lmt1(ia,:) = [ind mdist];
%     end
%     for ia = 1:sz2
%         [ind mdist] = mindist(x2(ia),y2(ia),x1,y1);
%         lmt2(ia,:) = [ind mdist];
%     end
%     nlmtch1 = [];
%     nlmtch2 = [];
%     for ia = sz1:-1:1
%         if(lmt1(ia,2) <= lmt2( lmt1(ia,1),2)  && lmt1(ia,2) < rad)
%             nlmtch1(ia,:) = [ nlst1(ia,1) , nlst2( lmt1(ia,1),1)];
%         else
%             nlmtch1(ia,:) = [ nlst1(ia,1) , 0];
%         end
%     end
%     for ia = sz2:-1:1
%         if(lmt2(ia,2) <= lmt1( lmt2(ia,1),2)  && lmt2(ia,2) < rad)
%             nlmtch2(ia,:) = [nlst1( lmt2(ia,1),1) ,  nlst2(ia,1)];
%         else
%             nlmtch2(ia,:) = [ 0 , nlst2(ia,1) ];
%         end
%     end
%     nlmtch = unique([ nlmtch1; nlmtch2],'rows');
%     nlmtch = nlmtch(nlmtch(:,1)~=0 & nlmtch(:,2)~=0,:);
% end

function nlmtch = nearestneighbourDI(nlst1,nlst2,ptch1,ptch2,rad,Inttol)
%Find the matching patched based on the distance i.e. scan through a defined
%radius, Pick up the right object matching the intensity

xydist = pdist2(ptch1(nlst1,[6,5]),ptch2(nlst2,[6,5]),'euclidean');
Intdist = pdist2(ptch1(nlst1,3),ptch2(nlst2,3),'euclidean');
Intdist(xydist>rad) = Inf;
Intdist = Intdist./ repmat(ptch1(nlst1,3),1,size(nlst2,1));
Intdist(Intdist > Inttol) = Inf;
sz = size(Intdist);
idxlst = zeros(sz(1),2);
ia = 0;
while(~all(isinf(Intdist(:))))
    ia = ia + 1;
    [m,idx] = min(Intdist(:));
    [I,J] = ind2sub(sz,idx);
    idxlst(ia,:) =  [I,J];
    Intdist(I,:) = Inf;
    Intdist(:,J) = Inf;
    
    
end
idxlst = idxlst(idxlst(:,1)>0,:);
nlmtch = [nlst1(idxlst(:,1)), nlst2(idxlst(:,2))];


% function nlmtch = nearestneighbourDI2(nlst1,nlst2,ptch1,ptch2,rad,tol)
% %Find the matching patched based on the distance i.e. scan through a defined
% %radius, Pick up the right object matching the intensity
% 
% xydist = pdist2(ptch1(nlst1,[6,5]),ptch2(nlst2,[6,5]),'euclidean');
% Intdist = pdist2(ptch1(nlst1,3),ptch2(nlst2,3),'euclidean');
% Intdist(xydist>rad) = Inf;
% Intdist = Intdist./ repmat(ptch1(nlst1,3),1,size(nlst2,1));
% %Intdist(Intdist > Inttol) = Inf;
% 
% [a,idx] = min(Intdist,[],2);
% ltol = a < tol;
% nlmtch = [nlst1(ltol), nlst2(idx(ltol))];
