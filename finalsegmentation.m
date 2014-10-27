function [npimg,nptch] = finalsegmentation(pimg,ptch,table)
% pimg  => patch image
% ptch  => patch parameter table
% table => cell array which is column of the finally matched table

nlst = [];
mlst = {};

nindex = [];
mindex = [];


for ia = 1:size(table,1)
    p = (table{ia});
    
    %    try
    
    if(sum(isempty(p))  || sum(isnan(p)))
        p =0;
    end
    
    %    catch
    %        p
    %    end
    
    p = p(:)';
    
%     if((size(p,2) > 1 ) || isequal(p,0))
    if((size(p,2) > 1 ))
        mlst(end+1,1) = {p};   % to be merged
        mindex(end+1) = ia;
    elseif(~(p==0))
        nlst(end+1,1) = p; %chosen normal
        nindex(end+1) = ia;
        
        
    end
end


%nlst(nlst == []) = size(ptch,1) +1;


% Creat new patch parameter table

% ptch(end+1,1) = 0;

nptch = zeros(size(table,1),size(ptch,2));

nptch(nindex',:) = ptch(nlst,:);

%include parameters for combined or merged matches

if(~isequal(mlst,cell(0)))
    for ia = 1:size(mlst,1)
        p = cell2mat(mlst(ia,1));
        nptch = addpatch_S(ptch,nptch,p,mindex(ia));
    end
end


nptch(1:size(table,1),1) = 1:size(table,1);


%Generate new patch image

npimg = newpimg(pimg,nlst,mlst,nindex,mindex);


function nptch = addpatch_S(ptch,nptch,p,n)
%adds multiple patches to one
%ptch -> original patch parameter table
%nptch _> new patch parameter table under which the merged patch will be
%added
%p -> is a vector of patch number that are supposed to be merged
%n is the index or position of this patch in the table


% n = size(nptch,1)+1;

nptch(n,1) = n;

if(p(1) ==0)
    nptch(n,2:7) =0;
else
    
    nptch(n,2) = sum(ptch(p,2));
    nptch(n,3) = sum(ptch(p,3));
    nptch(n,4) = nptch(n,3)/nptch(n,2);
    nptch(n,5) = sum(ptch(p,3).* ptch(p,5))/nptch(n,3);
    nptch(n,6) = sum(ptch(p,3).* ptch(p,6))/nptch(n,3);
    nptch(n,7) = 1;
end

function npimg = newpimg(pimg,nlst,mlst,nindex,mindex)

npimg = zeros(size(pimg));

%     ind = 1;
for ia = 1:size(nlst,1)
    
    npimg(pimg == nlst(ia,1)) = nindex(ia);
    
    %         ind = ind + 1;
    
    
end

for ib = 1:size(mlst,1)
    p = cell2mat(mlst(ib));
    
    for ic = 1:size(p,2)
        if(p(ic))
            %num = handles.mtch{n1}(p1(ic),jj);
            npimg(pimg == p(ic)) = mindex(ib);
        end
    end
    %         ind = ind + 1;
    
end
