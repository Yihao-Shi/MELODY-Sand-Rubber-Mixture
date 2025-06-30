function [masterPerSec,infoPerSec]=master_process(particle_num,countorNode,master,info,timestep)
%% Data modified
masterPerSec=zeros(particle_num,max(countorNode));
infoPerSec=zeros(particle_num,max(countorNode));
for i=1:particle_num
    for j=1:countorNode(i)
        if info(timestep,sum(countorNode(1:i-1))+j)<1e-16
            masterPerSec(i,j)=-1;
        else
            masterPerSec(i,j)=master(timestep,sum(countorNode(1:i-1))+j);
        end
        infoPerSec(i,j)=info(timestep,sum(countorNode(1:i-1))+j);
    end
end
end