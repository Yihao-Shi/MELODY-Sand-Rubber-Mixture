clear;
proportion='20';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
fx=dlmread([FRACTION 'fx_info.asc'],' ');
fy=dlmread([FRACTION 'fy_info.asc'],' ');
xnorm=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorm=dlmread([FRACTION 'ynorm_info.asc'],' ');
master=master(:,3:end-1);
len=len(:,3:end-1);
xnorm=xnorm(:,3:end-1);
ynorm=ynorm(:,3:end-1);
fx=fx(:,3:end-1);
fy=fy(:,3:end-1);

%% extract data
nodes=64;
particle_num=size(master,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%% 
strong_RR=zeros(size(master,1),1);DD=zeros(size(master,1),1);
strong_DR=zeros(size(master,1),1);DR=zeros(size(master,1),1);
strong_DD=zeros(size(master,1),1);RR=zeros(size(master,1),1);
Aver=zeros(size(master,1),1);
mforce0=0;
for timestep=241:241
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnormPerSec]=master_process(particle_num,nodes,master,fn,timestep);
    [contact_info,DD_num,DR_num,RR_num]=ContactPairs(masterPerSec,soft_num);

    fxPerSec=zeros(particle_num,nodes);
    fyPerSec=zeros(particle_num,nodes);
    xnormPerSec=zeros(particle_num,nodes);
    ynormPerSec=zeros(particle_num,nodes);
    for i=1:particle_num
        for j=1:nodes
            fxPerSec(i,j)=fx(timestep,(i-1)*nodes+j);
            fyPerSec(i,j)=fy(timestep,(i-1)*nodes+j);
            xnormPerSec(i,j)=xnorm(timestep,(i-1)*nodes+j);
            ynormPerSec(i,j)=ynorm(timestep,(i-1)*nodes+j);
        end
    end 
    
    for i=1:size(masterPerSec,1)
       for j=1:size(masterPerSec,2)
           end1=min(i-1,masterPerSec(i,j));
           end2=max(i-1,masterPerSec(i,j));
           for k=1:size(contact_info,1)
               if end1==contact_info(k,1)&end2==contact_info(k,2)
                   if i-1<masterPerSec(i,j)
                        contact_info(k,3)=contact_info(k,3)+fxPerSec(i,j);
                        contact_info(k,4)=contact_info(k,4)+fyPerSec(i,j);
                        contact_info(k,5)=contact_info(k,5)+xnormPerSec(i,j);
                        contact_info(k,6)=contact_info(k,6)+ynormPerSec(i,j);
                        normalize=sqrt(contact_info(k,5)^2+contact_info(k,6)^2);
                        contact_info(k,5)=contact_info(k,5)/normalize;
                        contact_info(k,6)=contact_info(k,6)/normalize;
                    else
                        contact_info(k,7)=contact_info(k,7)+fxPerSec(i,j);
                        contact_info(k,8)=contact_info(k,8)+fyPerSec(i,j);
                        contact_info(k,9)=contact_info(k,9)+xnormPerSec(i,j);
                        contact_info(k,10)=contact_info(k,10)+ynormPerSec(i,j);
                        normalize=sqrt(contact_info(k,9)^2+contact_info(k,10)^2);
                        contact_info(k,9)=contact_info(k,9)/normalize;
                        contact_info(k,10)=contact_info(k,10)/normalize;
                    end
                   break;
               end
           end
       end
    end
    contact_info(:,11)=contact_info(:,3).*contact_info(:,5)+contact_info(:,4).*contact_info(:,6);
    contact_info(:,12)=contact_info(:,7).*contact_info(:,9)+contact_info(:,8).*contact_info(:,10);
    
    RR_force=zeros(0,1);DR_force=zeros(0,1);DD_force=zeros(0,1);
    for i=1:size(contact_info,1)
        if contact_info(i,1)<soft_num&contact_info(i,2)<soft_num
            DD_force(length(DD_force)+1,1)=0.5*(contact_info(i,11)+contact_info(i,12));
        elseif contact_info(i,1)<soft_num&contact_info(i,2)>=soft_num
            DR_force(length(DR_force)+1,1)=0.5*(contact_info(i,11)+contact_info(i,12));
        else
            RR_force(length(RR_force)+1,1)=0.5*(contact_info(i,11)+contact_info(i,12));
        end
    end
    
    DD_force(DD_force==0)=[];DR_force(DR_force==0)=[];RR_force(RR_force==0)=[];
    mForce=(sum(DD_force)+sum(DR_force)+sum(RR_force))/(length(DD_force)+length(DR_force)+length(RR_force));
    DDstrong_num=length(find(DD_force>mForce));DRstrong_num=length(find(DR_force>mForce));RRstrong_num=length(find(RR_force>mForce));
    strong_DD(timestep)=DDstrong_num/(DDstrong_num+DRstrong_num+RRstrong_num);DD(timestep)=mean(DD_force);
    strong_DR(timestep)=DRstrong_num/(DDstrong_num+DRstrong_num+RRstrong_num);DR(timestep)=mean(DR_force);
    strong_RR(timestep)=RRstrong_num/(DDstrong_num+DRstrong_num+RRstrong_num);RR(timestep)=mean(RR_force);
    Aver(timestep)=mForce;
    
    mforce0=6267.5;%the average contact force in specimen V0
    if timestep==241
        force=[DD_force;DR_force;RR_force]./mforce0;
        result=HistRate(force);
    end
end
std = sqrt(sum((force-mean(force)).^2));
[ShearStrain]=Shearstrain(proportion);
figure,plot(ShearStrain,strong_DD,'b',ShearStrain,strong_DR,'g',ShearStrain,strong_RR,'r');
figure,plot(ShearStrain,DD,'b',ShearStrain,DR,'g',ShearStrain,RR,'r');
figure,semilogx(result(:,1),result(:,2),'g');hold on;