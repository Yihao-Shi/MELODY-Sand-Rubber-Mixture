clear;
proportion='20';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
master=master(:,3:end-1);
len=len(:,3:end-1);
fx=dlmread([FRACTION 'fx_info.asc'],' ');
fy=dlmread([FRACTION 'fy_info.asc'],' ');
xnorm=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorm=dlmread([FRACTION 'ynorm_info.asc'],' ');
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
    [masterPerSec,fnPerSec]=master_process(particle_num,nodes,master,fn,timestep);
    [contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);
    
    lengthPerSec=zeros(particle_num,nodes);
    for i=1:particle_num
        for j=1:nodes
            lengthPerSec(i,j)=len(timestep,(i-1)*nodes+j);
        end
    end
    
    for i=1:size(masterPerSec,1) %% 需要区分每个接触点
       for j=1:size(masterPerSec,2)
           end1=min(i-1,masterPerSec(i,j));
           end2=max(i-1,masterPerSec(i,j));
           for k=1:size(contact_info,1)
               if end1==contact_info(k,1)&end2==contact_info(k,2)
                   contact_info(k,3)=contact_info(k,3)+abs(lengthPerSec(i,j));
                   break;
               end
           end
       end
    end
    contact_info(contact_info(:,3)==0,:)=[];
    
    RR_len=[];DR_len=[];DD_len=[];
    for i=1:size(contact_info,1)
        if contact_info(i,1)<soft_num&contact_info(i,2)<soft_num
            DD_len=[DD_len;contact_info(i,3)];
        elseif contact_info(i,1)<soft_num&contact_info(i,2)>=soft_num
            DR_len=[DR_len;contact_info(i,3)];
        elseif contact_info(i,1)>=soft_num&contact_info(i,2)>=soft_num
            RR_len=[RR_len;contact_info(i,3)];
        end
    end
    
    c_50=2*pi*0.0125;
    DD(timestep)=mean(DD_len)/c_50;
    DR(timestep)=mean(DR_len)/c_50;
    RR(timestep)=mean(RR_len)/c_50;
end
plot(1:241,RR,1:241,DR,1:241,DD)