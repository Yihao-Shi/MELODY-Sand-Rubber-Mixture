clear;
proportion='0';
timestep=241;
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
locx=dlmread([FRACTION 'locx_info.asc'],' ');
locy=dlmread([FRACTION 'locy_info.asc'],' ');
locx=locx(timestep,3:end-1);
locy=locy(timestep,3:end-1);

%% extract data
nodes=64;
particle_num=size(locx,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
info=textread('../../Simulation_DATA/DirectShearTest/ball_pos_PFC/Ball&Wall_info.txt','%f');
x_min=info(1);
x_max=info(2);
rhi=info(6);
x_min=-6*ceil(-x_min/(6*rhi))*rhi;
x_max=6*ceil(x_max/(6*rhi))*rhi;
%%
locxPerSec=zeros(particle_num,nodes);
locyPerSec=zeros(particle_num,nodes);
for i=1:particle_num
    for j=1:nodes
        locxPerSec(i,j)=locx((i-1)*nodes+j);
        locyPerSec(i,j)=locy((i-1)*nodes+j);
    end
    while mean(locxPerSec(i,:))>x_max
        locxPerSec(i,:)=locxPerSec(i,:)-(x_max-x_min);
    end
end
%% Output
fid = fopen([FRACTION 'location' num2str(timestep) '.txt'],'w');
for i=1:size(locxPerSec,1)
    for j=1:2*size(locxPerSec,2)
        if j<=size(locxPerSec,2)
            fprintf(fid,'%f ',locxPerSec(i,j));
        else
            fprintf(fid,'%f ',locyPerSec(i,j-size(locxPerSec,2)));
        end
        if j==2*size(locxPerSec,2)
            fprintf(fid,'\r\n');
        end
    end
end
fclose(fid); 

