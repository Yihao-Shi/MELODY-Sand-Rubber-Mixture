clear;
proportion='10';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   len(i,:)=lens((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).';
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
timestep=1;
fn=fx.*xnorm+fy.*ynorm;
[masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);

locxPerSec=zeros(particle_num,max(countorNode));
locyPerSec=zeros(particle_num,max(countorNode));
for i=1:particle_num
    for j=1:countorNode(i)
        locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
        locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
    end
end

area=zeros(particle_num,1);
for i=1:particle_num
   x=locxPerSec(i,1:countorNode(i));y=locyPerSec(i,1:countorNode(i));
   for j=1:countorNode(i)-2
       area(i)=area(i)+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
   end
end
    
ratio=sum(area(1:soft_num))/(sum(area));    