clear
proportion='40';
timestep=241;
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
xpos=dlmread([FRACTION 'xpos_info.asc'],' ');
ypos=dlmread([FRACTION 'ypos_info.asc'],' ');
fx=dlmread([FRACTION 'fx_info.asc'],' ');
fy=dlmread([FRACTION 'fy_info.asc'],' ');
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
xpos=xpos(:,3:end-1);
ypos=ypos(:,3:end-1);
fx=fx(:,3:end-1);
fy=fy(:,3:end-1);
master=master(:,3:end-1);
len=len(:,3:end-1);

%% extract data
nodes=64;
particle_num=size(master,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
info=textread('../../Simulation_DATA/DirectShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
rhi=info(6);
x_min=info(1);
x_max=info(2);
x_min=-6*ceil(-x_min/(6*rhi))*rhi;
x_max=6*ceil(x_max/(6*rhi))*rhi;
%% 
[masterPerSec,~]=master_process(particle_num,nodes,master,len,timestep);
[contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);

fxPerSec=zeros(particle_num,nodes);
xposPerSec=zeros(particle_num,1);
fyPerSec=zeros(particle_num,nodes);
yposPerSec=zeros(particle_num,1);
for i=1:particle_num
    xposPerSec(i)=xpos(timestep,i);
    yposPerSec(i)=ypos(timestep,i);
    for j=1:nodes
        fxPerSec(i,j)=fx(timestep,(i-1)*nodes+j);
        fyPerSec(i,j)=fy(timestep,(i-1)*nodes+j);
    end
    while xposPerSec(i)>x_max
        xposPerSec(i)=xposPerSec(i)-(x_max-x_min);
    end
end
    
contact_info(find(contact_info(:,2)==797|contact_info(:,2)==798),:)=[];
for i=1:size(masterPerSec,1)
   for j=1:size(masterPerSec,2)
       end1=min(i-1,masterPerSec(i,j));
       end2=max(i-1,masterPerSec(i,j));
       for k=1:size(contact_info,1)
           if end1==contact_info(k,1)&end2==contact_info(k,2)
               if i-1<masterPerSec(i,j)
                   contact_info(k,3)=contact_info(k,3)+abs(fxPerSec(i,j));
                   contact_info(k,6)=contact_info(k,6)+abs(fyPerSec(i,j));
               else
                   contact_info(k,4)=contact_info(k,4)+abs(fxPerSec(i,j));
                   contact_info(k,7)=contact_info(k,7)+abs(fyPerSec(i,j));
               end
               contact_info(k,5)=contact_info(k,3)+contact_info(k,4);
               contact_info(k,8)=contact_info(k,6)+contact_info(k,7);
               contact_info(k,9)=xposPerSec(end1+1);
               contact_info(k,10)=yposPerSec(end1+1);
               contact_info(k,11)=xposPerSec(end2+1);
               contact_info(k,12)=yposPerSec(end2+1);
               contact_info(k,13)=sqrt(contact_info(k,5)^2+contact_info(k,8)^2);
           end
       end
   end
end

contact_info(find(contact_info(:,13)==0),:)=[];

copy_line=find(abs(contact_info(:,9)-contact_info(:,11))>0.7*(x_max-x_min));
for i=1:size(copy_line,1)
    contact_info(size(contact_info,1)+1,:)=contact_info(copy_line(i),:);
    if contact_info(copy_line(i),9)>contact_info(copy_line(i),11)
        contact_info(copy_line(i),11)=contact_info(copy_line(i),11)+(x_max-x_min);
        contact_info(size(contact_info,1),9)=contact_info(size(contact_info,1),9)-(x_max-x_min);
    else
        contact_info(copy_line(i),11)=contact_info(copy_line(i),11)-(x_max-x_min);
        contact_info(size(contact_info,1),9)=contact_info(size(contact_info,1),9)+(x_max-x_min);
    end
end

for i=1:size(contact_info,1)
    if contact_info(i,1)<soft_num&contact_info(i,2)<soft_num
        contact_info(i,14)=2;
    elseif contact_info(i,1)<soft_num&contact_info(i,2)>=soft_num
        contact_info(i,14)=1;
    else
        contact_info(i,14)=0;
    end
end

%% Output
fid = fopen([FRACTION 'linewidth' num2str(timestep) '.txt'],'w');
for i=1:size(contact_info,1)
    fprintf(fid,'%f %f %f %f %f %f\r\n',contact_info(i,9),contact_info(i,10),contact_info(i,11),contact_info(i,12),contact_info(i,13),contact_info(i,14));
end
fclose(fid); 
