clear;
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end
for timestep=1:36
    %%
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
        end
    end

    %% Output
    fid = fopen([FRACTION 'location' num2str(timestep) '.txt'],'w');
    for i=1:particle_num
        for j=1:2*max(countorNode)+1
            if j<=countorNode(i)
                fprintf(fid,'%f ',locxPerSec(i,j));
            elseif j>countorNode(i)&&j<=2*countorNode(i)
                fprintf(fid,'%f ',locyPerSec(i,j-countorNode(i)));
            elseif j>2*countorNode(i)&&j<=2*max(countorNode)
                fprintf(fid,'0 ');
            elseif j==2*max(countorNode)+1
                fprintf(fid,'%f\r\n',countorNode(i));
            end
        end
    end
    fclose(fid); 
end

