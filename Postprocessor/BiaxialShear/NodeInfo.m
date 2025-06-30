function [softFieldNode,countorNode]=NodeInfo(FRACTION,total_num,soft_num)
data=fopen([FRACTION 'STATIC_CONTROL.asc']);
scp=0;cp=0;
softFieldNode=zeros(soft_num,1);countorNode=zeros(total_num,1);
while cp<total_num
    tline=fgets(data);
    if strfind(tline,'NODES')~=0&scp<soft_num
        scp=scp+1;
        tline=fgets(data);
        softFieldNode(scp,1)=str2double(tline);
    end
    if strfind(tline,'Closed Linear')~=0
        cp=cp+1;
        tline=fgets(data);
        countorNode(cp,1)=str2double(tline(1:3))-1;
    end
end
fclose(data);
end