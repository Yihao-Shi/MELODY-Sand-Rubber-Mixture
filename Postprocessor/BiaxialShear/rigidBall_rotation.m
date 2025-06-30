clear;
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];

%% extract data
particle_num=5025;
total_num=5029;
rigid_num=ceil((1-str2double(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
source=load(['source' proportion '.mat']);
[~,countorNode]=NodeInfo(FRACTION,total_num,soft_num);

%%
timestep=36;
disprPerSec=zeros(particle_num,1);
for i=soft_num+1:total_num
    kinematic=source.kinematics{timestep,i};
    disprPerSec(i,1)=kinematic(1,3);
end

fid = fopen([FRACTION 'rotation' num2str(timestep) '.txt'],'w');
for i=soft_num+1:particle_num
    fprintf(fid,'%f\r\n',disprPerSec(i,1));
end
fclose(fid); 