function table = bigtableindex(t)


sz = size(t,2);

table = t(1).m;
table(cellfun(@(x)any(x == 0),table)) = {[]};

table(sum(cellfun(@isempty,table),2) == 2,:) = [];

for ia = 2:sz
    
    lst = t(ia).m;
   
 try
    lst(cellfun(@(x)any(x == 0),lst)) = {[]};
    lst(sum(cellfun(@isempty,lst),2) == 2,:) = [];
 catch
     'wait'
 end
    
   

   
        cur_tp = ia;
   
    
    table = edittable(table,cur_tp,lst);
    
   ia 
  
    
end

table(cellfun(@isempty,table)) = {[]};



function table = edittable(table,cur_tp,lst)


%cur_tp is index for comparison image 1
%next_tp is the index for comparison image 2

next_tp = cur_tp +1;

ib = 1;

idx4miss = true(size(lst,1),1);
while(ib<=size(table,1))
    
    
    
    
    key = table{ib,cur_tp};
    
    
    
    %         ind = 1;
    
    table{ib,next_tp} = [];
    
%     if(isequal(key,0))
%         key = NaN;
%     end

    if(~isempty(key))
    
    key = key(:);
    
    for ic = 1:size(key,1)
        
        [pos , st] = findkey(key(ic),lst(:,1));
        
        if(pos)
            table{ib,next_tp} = unique([table{ib,next_tp}; lst{pos,2}(:) ]); %can make unique at the end too ??? to increase speed
            
            if(st)
                
                
                table = shrinktable(table,lst{pos,1},cur_tp);
                
                
            end
            
        
            idx4miss(pos) = false;
            
            
        end
        
        
        
        %             catch
        %                 pos
        %                 key
        %             end
    end
    %         if(~pos) should collect all pos array, find the missing and
    %         should add it at the end to get those which are in one
    
    end
    
    ib = ib+1;
    
end

% misslst = lst(idx4miss,2);

idx4miss(cellfun(@isempty,lst(:,2))) =  false;

tidx = cellfun(@isempty,table(:,end));
table(tidx,end) = {[]};

sz=size(table,1);


if(any(idx4miss))
    table(sz+1:sz+sum(idx4miss),next_tp) = lst(idx4miss,2);
    
end

function [pos st] = findkey(key,lst)


lstel = cellfun(@numel,lst);
for ia = size(lst,1):-1:1, 
    for ib = 1:lstel(ia),
        nlst(ia,ib)  = lst{ia}(ib);       
    end
end

[pos,y] =  ind2sub(size(nlst),find(nlst == key));

st = lstel(pos) > 1; %shrink table

if(isempty(pos))
    pos = [];
    st = [];
else
    pos = pos(1);
end

% tpos = pos;
% tst = st;
% %---------------check
% 
% pos = [];
% st = 0;
% for ia = 1:size(lst,1)
%     
%     gkey = lst{ia};
%     
%     tf = ismember(key,gkey);
%     
%     if(tf)
%         pos = ia;
%         gkey = gkey(:);
%         
%         if(size(gkey,1)>1)
%             st = 1;   %shrink table
%         end
%     end
%     
%     
% end
% 
% if(~isequal(pos,tpos) || ~isequal(st,  tst))
%     'wait'
% end




function stable  = shrinktable(table,keys,ind)

keys = keys(:);

keys = sort(keys);



inpos = findkey(keys(1),table(:,ind));
stable = table;
%continue here

for ia = 2:size(keys,1)
    pos = findkey(keys(ia),stable(:,ind));
    
    if(inpos ~= pos)
        [stable inpos] = mergerows(stable,[inpos,pos]);
%     else
%         stable = table;
    end
    
end

function [table r1] = mergerows(table,rows2mrg)

rows2mrg = sort(rows2mrg);
r1 = rows2mrg(1);
r2 = rows2mrg(2);
for ia = 1:size(table,2)

    table{r1,ia} = [table{r1,ia}(:); table{r2,ia}(:)];
end
table(r2,:) = [];
