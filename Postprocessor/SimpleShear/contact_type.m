clear
proportion='40';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
master=dlmread([FRACTION 'master_info.asc'],' ');
len=dlmread([FRACTION 'length_info.asc'],' ');
master=master(:,3:end-1);
len=len(:,3:end-1);

%% extract data
nodes=64;
particle_num=size(master,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%% 
RR=zeros(size(master,1),1);
DR=zeros(size(master,1),1);
DD=zeros(size(master,1),1);
for timestep=1:size(master,1)
    [masterPerSec,~]=master_process(particle_num,nodes,master,len,timestep);
    [contact_info,RR_num,DR_num,DD_num]=ContactPairs(masterPerSec,soft_num);
    
    RR(timestep)=RR_num;
    DR(timestep)=DR_num;
    DD(timestep)=DD_num;
end
plot(1:241,RR,1:241,DR,1:241,DD);