clear;
proportion='20';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
xpos=dlmread([FRACTION 'xpos_info.asc'],' ');
ypos=dlmread([FRACTION 'ypos_info.asc'],' ');
master=masters(:,3:end-1);len=lens(:,3:end-1);
fx=fxs(:,3:end-1);fy=fys(:,3:end-1);
xnorm=xnorms(:,3:end-1);ynorm=ynorms(:,3:end-1);
xpos=xpos(:,3:end-1);ypos=ypos(:,3:end-1);

%% extract data
nodes=64;
particle_num=796;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%%
a_c=zeros(size(master,1),1);a_cDD=zeros(size(master,1),1);a_cDR=zeros(size(master,1),1);a_cRR=zeros(size(master,1),1);
a_n=zeros(size(master,1),1);a_nDD=zeros(size(master,1),1);a_nDR=zeros(size(master,1),1);a_nRR=zeros(size(master,1),1);
bins=20;
dbin=(-90:180/bins:90)*pi/180;
for timestep=1:1
    [masterPerSec,lengthPerSec]=master_process(particle_num,nodes,master,len,timestep);

    fxPerSec=zeros(particle_num,nodes);
    fyPerSec=zeros(particle_num,nodes);
    xnormPerSec=zeros(particle_num,nodes);
    ynormPerSec=zeros(particle_num,nodes);
    xposPerSec=zeros(particle_num,nodes);
    yposPerSec=zeros(particle_num,nodes);
    for i=1:particle_num
        for j=1:nodes
            fxPerSec(i,j)=fx(timestep,(i-1)*nodes+j);
            fyPerSec(i,j)=fy(timestep,(i-1)*nodes+j);
            xnormPerSec(i,j)=xnorm(timestep,(i-1)*nodes+j);
            ynormPerSec(i,j)=ynorm(timestep,(i-1)*nodes+j);
            xposPerSec(i,j)=xpos(timestep,i);
            yposPerSec(i,j)=ypos(timestep,i);
        end
    end
    fnPerSec=abs(fxPerSec.*xnormPerSec+fyPerSec.*ynormPerSec);
    ftPerSec=abs(-fxPerSec.*ynormPerSec+fyPerSec.*xnormPerSec);

    fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),7);
    for i=1:size(masterPerSec,1)
        for j=1:size(masterPerSec,2)
            if masterPerSec(i,j)>=0&masterPerSec(i,j)<796
                fc((i-1)*nodes+j,1)=i-1;fc((i-1)*nodes+j,2)=masterPerSec(i,j);
                fc((i-1)*nodes+j,3)=xnormPerSec(i,j);
                fc((i-1)*nodes+j,4)=ynormPerSec(i,j);
                fc((i-1)*nodes+j,5)=fnPerSec(i,j);
            end
        end
    end
    fc(find(fc(:,5)==0),:)=[];

    RR=[];DR=[];DD=[];
    for i=1:size(fc,1)
        if fc(i,1)<soft_num&fc(i,2)<soft_num
            DD=[DD;fc(i,3),fc(i,4),fc(i,5)];
        elseif fc(i,1)<soft_num&fc(i,2)>=soft_num
            DR=[DR;fc(i,3),fc(i,4),fc(i,5)];
        else
            RR=[RR;fc(i,3),fc(i,4),fc(i,5)];
        end
    end
    Total=[DD;DR;RR];

    [acDD,acDR,acRR,ac,devia_acDD,devia_acDR,devia_acRR,devia_ac]=Fabric_ac(DD,DR,RR);
    [anDD,anDR,anRR,an]=Fabric_an(DD,DR,RR,devia_acDD,devia_acDR,devia_acRR,devia_ac);
    a_c(timestep)=ac;a_cDD(timestep)=acDD;a_cDR(timestep)=acDR;a_cRR(timestep)=acRR;
    a_n(timestep)=an;a_nDD(timestep)=anDD;a_nDR(timestep)=anDR;a_nRR(timestep)=anRR;
end
[ShearStrain]=Shearstrain(proportion);
figure,plot(ShearStrain,a_c,'blue',ShearStrain,a_cDD,'g',ShearStrain,a_cDR,'r',ShearStrain,a_cRR,'black');
figure,plot(ShearStrain,a_n,'blue',ShearStrain,a_nDD,'g',ShearStrain,a_nDR,'r',ShearStrain,a_nRR,'black');