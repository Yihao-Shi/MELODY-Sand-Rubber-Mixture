clear;
proportion='20';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
mises=dlmread([FRACTION 'mises_info.asc'],' ');
master=master(:,3:end-1);
len=len(:,3:end-1);
mises=mises(:,3:end-1);

%% extract data
nodes=64;
particle_num=size(master,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%%
DR=zeros(size(master,1),1);
DD=zeros(size(master,1),1);
Total=zeros(size(master,1),1);
for timestep=1:size(master,1)
    [masterPerSec,~]=master_process(particle_num,nodes,master,len,timestep);
    [contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);
    
    misesSec=zeros(particle_num,nodes);
    for i=1:soft_num
        for j=1:nodes
            misesSec(i,j)=mises(timestep,(i-1)*nodes+j);
        end
    end
    
    for i=1:size(masterPerSec,1)
       for j=1:size(masterPerSec,2)
           end1=min(i-1,masterPerSec(i,j));
           end2=max(i-1,masterPerSec(i,j));
           for k=1:size(contact_info,1)
               if end1==contact_info(k,1)&end2==contact_info(k,2)
                   if i-1<masterPerSec(i,j)
                       contact_info(k,3)=contact_info(k,3)+misesSec(i,j);
                   else
                       contact_info(k,4)=contact_info(k,4)+misesSec(i,j);
                   end
                   break;
               end
           end
       end
    end
    
    RR_mises=0;DR_mises=0;DD_mises=0;
    for i=1:size(contact_info,1)
        if contact_info(i,1)<soft_num&contact_info(i,2)<soft_num
            DD_mises(length(DD_mises)+1,1)=0.5*(contact_info(i,3)+contact_info(i,4));
        elseif contact_info(i,1)<soft_num&contact_info(i,2)>=soft_num
            DR_mises(length(DR_mises)+1,1)=contact_info(i,3);
        end
    end
    
    DD(timestep)=mean(DD_mises);
    DR(timestep)=mean(DR_mises);
    Total(timestep)=0.5*(DR(timestep)+DD(timestep))/mean(DR_mises);
end
%%
[ShearStrain]=Shearstrain(proportion);
plot(ShearStrain,Total,'b');