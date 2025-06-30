clear;
proportion='0';
FRACTION=['../../Simulation_DATA/DirectShearTest/' proportion 'percent/'];
xpos=dlmread([FRACTION 'xpos_info.asc'],' ');
ypos=dlmread([FRACTION 'ypos_info.asc'],' ');
locx=dlmread([FRACTION 'locx_info.asc'],' ');
locy=dlmread([FRACTION 'locy_info.asc'],' ');
xpos=xpos(:,3:end-1);
ypos=ypos(:,3:end-1);
locx=locx(:,3:end-1);
locy=locy(:,3:end-1);

%% extract data
nodes=64;
particle_num=size(locx,2)/nodes;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

%%
rotate=zeros(size(locx,1),1);
for timestep=1:241
    pos=[];
    for i=soft_num+1:particle_num
        for j=1:nodes
            pos=[pos;xpos(timestep,i),ypos(timestep,i)];
        end
    end
    loc=[locx(timestep,nodes*soft_num+1:nodes*particle_num).',locy(timestep,nodes*soft_num+1:nodes*particle_num).'];
    
    if timestep==1
        orient=[loc(:,1)-pos(:,1),loc(:,2)-pos(:,2)];
    else
        orient_ini=orient;
        orient=[loc(:,1)-pos(:,1),loc(:,2)-pos(:,2)];
        angle=acos((orient(:,1).*orient_ini(:,1)+orient(:,2).*orient_ini(:,2))./(sqrt(orient(:,1).^2+orient(:,2).^2).*sqrt(orient_ini(:,1).^2+orient_ini(:,2).^2)))*180/pi;
        rotate(timestep)=mean(angle);
    end
end
for i=2:length(rotate)
    rotate(i)=rotate(i-1)+rotate(i);
end

fid = fopen([FRACTION 'rotation' num2str(timestep) '.txt'],'w');
for i=1:64:size(angle,1)
    fprintf(fid,'%f\r\n',mean(angle(i:i+62)));
end
fclose(fid); 