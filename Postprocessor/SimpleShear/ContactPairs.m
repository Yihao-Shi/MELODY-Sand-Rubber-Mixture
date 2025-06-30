function [contact_info,RR_num,DR_num,DD_num]=ContactPairs(array,soft_num)
%% Initialization 
contact_info=zeros(0,14);
flag=0;

%% loop for each contact
for i=1:size(array,1)
    temp=array(i,:);
    temp(find(temp(1,:)==-1))=[];
    end1=i-1;
    end2=unique(temp);
    
    row=size(contact_info,1);
    for j=row+1:row+length(end2)
        for m=1:row
            if min(end1,end2(j-row))==contact_info(m,1)&max(end1,end2(j-row))==contact_info(m,2)
                flag=1;break;
            else
                flag=0;
            end
        end
        if flag==0
            contact_info(size(contact_info,1)+1,1)=min(end1,end2(j-row));
            contact_info(size(contact_info,1),2)=max(end1,end2(j-row));
        end
    end
end
[~,idx]=sort(contact_info(:,1))
contact_info=contact_info(idx,:)
contact_info(find(contact_info(:,2)>=797),:)=[];

RR_num=0;DR_num=0;DD_num=0;
for i=1:size(contact_info,1)
   if contact_info(i,1)<soft_num&contact_info(i,2)<soft_num
       DD_num=DD_num+1;
   elseif contact_info(i,1)<soft_num&contact_info(i,2)>=soft_num
       DR_num=DR_num+1;
   elseif contact_info(i,1)>=soft_num&contact_info(i,2)>=soft_num
       RR_num=RR_num+1;
   end
end
end