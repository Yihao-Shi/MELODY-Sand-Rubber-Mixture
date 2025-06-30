clear;
proportion='40';
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
RR=zeros(size(master,1),1);
DR=zeros(size(master,1),1);
DD=zeros(size(master,1),1);
Aver=zeros(size(master,1),1);
for timestep=1:size(master,1)
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnormPerSec]=master_process(particle_num,nodes,master,fn,timestep);
    [contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);
    
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
    contact_info(:,13)=abs(-contact_info(:,3).*contact_info(:,6)+contact_info(:,4).*contact_info(:,5));
    contact_info(:,14)=abs(-contact_info(:,7).*contact_info(:,10)+contact_info(:,8).*contact_info(:,9));
    
    RR_slid=zeros(0,1);DR_slid=zeros(0,1);DD_slid=zeros(0,1);
    for i=1:size(contact_info,1)
        if abs(contact_info(i,13)/contact_info(i,11))>0.499||abs(contact_info(i,14)/contact_info(i,12))>0.499
            flag=1;
        else
            flag=0;
        end
        if contact_info(i,1)<soft_num&contact_info(i,2)<soft_num
            DD_slid(length(DD_slid)+1,1)=flag;
        elseif contact_info(i,1)<soft_num&contact_info(i,2)>=soft_num
            DR_slid(length(DR_slid)+1,1)=flag;
        else
            RR_slid(length(RR_slid)+1,1)=flag;
        end
    end
    
    DD(timestep)=mean(DD_slid)%/(length(DD_slid)+length(DR_slid)+length(RR_slid));
    DR(timestep)=mean(DR_slid)%/(length(DD_slid)+length(DR_slid)+length(RR_slid));
    RR(timestep)=mean(RR_slid)%/(length(DD_slid)+length(DR_slid)+length(RR_slid));
    Aver(timestep)=(sum(DD_slid)+sum(DR_slid)+sum(RR_slid))/(length(DD_slid)+length(DR_slid)+length(RR_slid));
end
figure,plot(1:241,DD,'r',1:241,DR,'b',1:241,RR,'g',1:241,Aver,'p');