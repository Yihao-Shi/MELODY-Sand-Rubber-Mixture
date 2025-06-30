clear
proportion='0';
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
x_min=info(1);
x_max=info(2);
y_min=info(3);
y_max=info(4);
rhi=info(6);
Force=dlmread([FRACTION 'force_info.asc'],' ');
Disp=dlmread([FRACTION 'displacement_info.asc'],' ');
VerticalDisp=Disp(:,6)-Disp(:,4);
ShearStress=0.5*(Force(:,3)-Force(:,5))/(x_max-x_min);
NormalStress=0.5*(Force(:,6)-Force(:,4))/(x_max-x_min);
ShearStress=ShearStress([1,10:10:2400]);
NormalStress=NormalStress([1,10:10:2400]);
S=(x_max-x_min)*(y_max-y_min+VerticalDisp);
%% 
sigma12=zeros(size(master,1),1);
sigma22=zeros(size(master,1),1);
for timestep=1:1%size(master,1)
    [masterPerSec,~]=master_process(particle_num,nodes,master,len,timestep);
    [contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);

    fxPerSec=zeros(particle_num,nodes);
    xposPerSec=zeros(particle_num,nodes);
    fyPerSec=zeros(particle_num,nodes);
    yposPerSec=zeros(particle_num,nodes);
    for i=1:particle_num
        for j=1:nodes
            fxPerSec(i,j)=fx(timestep,(i-1)*nodes+j);
            xposPerSec(i,j)=xpos(timestep,i);
            fyPerSec(i,j)=fy(timestep,(i-1)*nodes+j);
            yposPerSec(i,j)=ypos(timestep,i);
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
                   if abs(xposPerSec(end1+1,j)-xposPerSec(end2+1,j))>x_max
                       contact_info(k,5)=(x_max-x_min)-abs(xposPerSec(end1+1,j)-xposPerSec(end2+1,j));
                   else
                       contact_info(k,5)=abs(xposPerSec(end1+1,j)-xposPerSec(end2+1,j));
                   end
                   contact_info(k,8)=abs(yposPerSec(end1+1,j)-yposPerSec(end2+1,j));
                   break;
               end
           end
       end
    end
    
    for i=1:size(contact_info,1)
        sigma12(timestep)=sigma12(timestep)+(contact_info(i,3)+contact_info(i,4))*contact_info(i,8)/S(timestep);
        sigma22(timestep)=sigma22(timestep)+(contact_info(i,6)+contact_info(i,7))*contact_info(i,8)/S(timestep);
    end
end
figure,plot(0:241,[0;sigma12./sigma22],1:241,ShearStress./NormalStress,'p');