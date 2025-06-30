clear;
proportion='20';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
locx=dlmread([FRACTION 'locx_info.asc'],' ');
locy=dlmread([FRACTION 'locy_info.asc'],' ');
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
fx=dlmread([FRACTION 'fx_info.asc'],' ');
fy=dlmread([FRACTION 'fy_info.asc'],' ');
xnorm=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorm=dlmread([FRACTION 'ynorm_info.asc'],' ');
id=[];
for i=3:size(locx,2)-1
    if mod(i-2,64)==0
        continue;
    else
        id=[id,i];
    end
end
master=master(:,id);
len=len(:,id);
locx=locx(:,id);
locy=locy(:,id);
xnorm=xnorm(:,id);
ynorm=ynorm(:,id);
fx=fx(:,id);
fy=fy(:,id);
%% extract data
nodes=63;
particle_num=length(id)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
info=textread('../../Simulation_DATA/DirectShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
x_min=info(1);
x_max=info(2);
rhi=info(6);
x_min=-6*ceil(-x_min/(6*rhi))*rhi;
x_max=6*ceil(x_max/(6*rhi))*rhi;
%%
cir=zeros(size(locx,1),1);
iso=zeros(size(locx,1),1);
for timestep=[1:40:240,241]
    locxPerSec=zeros(soft_num,nodes);
    locyPerSec=zeros(soft_num,nodes);
    for i=1:soft_num
        for j=1:nodes
            locxPerSec(i,j)=locx(timestep,(i-1)*nodes+j);
            locyPerSec(i,j)=locy(timestep,(i-1)*nodes+j);
        end
        while mean(locxPerSec(i,:))>x_max
            locxPerSec(i,:)=locxPerSec(i,:)-(x_max-x_min);
        end
    end
    
    equi=zeros(soft_num,1);peri=zeros(soft_num,1);
    for i=1:soft_num
        x=locxPerSec(i,:);y=locyPerSec(i,:);area=0;
        for j=1:nodes-2
            area=area+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
        end
        equi(i)=2*sqrt(area*pi);
        for k=1:nodes
            if k==nodes
                peri(i)=peri(i)+sqrt((x(k)-x(1))^2+(y(k)-y(1))^2);
            else
                peri(i)=peri(i)+sqrt((x(k+1)-x(k))^2+(y(k+1)-y(k))^2);
            end
        end
    end
    factor=equi./peri; 
    cir(timestep)=mean(factor);
    
    fxPerSec=zeros(soft_num,nodes);
    fyPerSec=zeros(soft_num,nodes);
    xnormPerSec=zeros(soft_num,nodes);
    ynormPerSec=zeros(soft_num,nodes);
    for i=1:soft_num
        for j=1:nodes
            fxPerSec(i,j)=fx(timestep,(i-1)*nodes+j);
            fyPerSec(i,j)=fy(timestep,(i-1)*nodes+j);
            xnormPerSec(i,j)=xnorm(timestep,(i-1)*nodes+j);
            ynormPerSec(i,j)=ynorm(timestep,(i-1)*nodes+j);
        end
    end
    fnormPerSec=fxPerSec.*xnormPerSec+fyPerSec.*ynormPerSec;
    
    isofac=zeros(0,2);
    [masterPerSec,lengthPerSec]=master_process(particle_num,nodes,master,len,timestep);
    [contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);
    for i=1:soft_num
        for j=1:nodes
            end1=min(i-1,masterPerSec(i,j));
            end2=max(i-1,masterPerSec(i,j));
            for k=1:size(contact_info,1)
                if end1==contact_info(k,1)&end2==contact_info(k,2)
                    contact_info(k,3)=contact_info(k,3)+lengthPerSec(i,j);
                    contact_info(k,4)=contact_info(k,4)+factor(i);
                    contact_info(k,5)=contact_info(k,5)+1;
                    break;
                end
            end
        end
    end
    for i=1:soft_num
        if contact_info(find(contact_info(:,1)==i-1),2)>=soft_num
            isofac(i,1)=sum(contact_info(find(contact_info(:,1)==i-1),4))/sum(contact_info(find(contact_info(:,1)==i-1),5));
            isofac(i,2)=sum(contact_info(find(contact_info(:,1)==i-1),3));
        end
    end
    ind=find(min(isofac(:,2)));
    isofac(find(isofac(:,1)==0),:)=[];isofac=sortrows(isofac,1);
    
    bins=5;pernum=ceil(length(isofac)/bins);
    fitting=zeros(bins,2);dbin=[isofac(1,1)];
    for i=1:bins-1
        dbin=[dbin,isofac(i*pernum,1)];
    end
    dbin=[dbin,isofac(end,1)];
    for i=1:length(dbin)-1
        temp=find(isofac(:,1)>=dbin(i)&isofac(:,1)<dbin(i+1));
        fitting(i,1)=0.5*(dbin(i)+dbin(i+1));
        fitting(i,2)=mean(isofac(temp,2));
        fitting(i,3)=var(isofac(temp,2),1);
    end
    errorbar(fitting(:,1),fitting(:,2),fitting(:,3));hold on;
end
% [ShearStrain]=Shearstrain(proportion);
% figure,plot(ShearStrain,cir,ShearStrain,iso,'o');