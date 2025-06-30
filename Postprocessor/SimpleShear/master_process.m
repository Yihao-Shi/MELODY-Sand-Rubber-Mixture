function [masterPerSec,infoPerSec]=master_process(particle_num,nodes,master,info,timestep)
%% Data modified
masterPerSec=zeros(particle_num,nodes);
infoPerSec=zeros(particle_num,nodes);
for i=1:particle_num
    for j=1:nodes
        if info(timestep,(i-1)*nodes+j)<=1e-5
            masterPerSec(i,j)=-1;
        else
            masterPerSec(i,j)=master(timestep,(i-1)*nodes+j);
        end
        infoPerSec(i,j)=info(timestep,(i-1)*nodes+j);
    end
end
end