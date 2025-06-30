clear;
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
masses=dlmread([FRACTION 'mass_info.asc'],' ');

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
   mass(i,:)=masses((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
xstart=-0.625;xend=0.625;
ystart=-1.25;yend=1.25;
dbiny=0.25;dbinx=0.125;
nbiny=(yend-ystart)/dbiny;nbinx=(xend-xstart)/dbinx;
    
x_number_fraction=zeros(36,nbinx);x_volume_fraction=zeros(36,nbinx);
y_number_fraction=zeros(36,nbiny);y_volume_fraction=zeros(36,nbiny);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);

    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    massPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            if i<=soft_num
                massPerSec(i,j)=mass(timestep,sum(countorNode(1:i-1))+j);
            else
                massPerSec(i,j)=1;
            end
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
        end
    end

    area=zeros(particle_num,1);
    center=zeros(particle_num,2);
    for i=1:particle_num
       x=locxPerSec(i,1:countorNode(i));y=locyPerSec(i,1:countorNode(i));m=massPerSec(i,1:countorNode(i));
       for j=1:countorNode(i)-2
           area(i)=area(i)+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
       end

       mass_center=zeros(2);
       for j=1:countorNode(i)
           mass_center(1)=mass_center(1)+x(j)*m(j);
           mass_center(2)=mass_center(2)+y(j)*m(j);
       end
       center(i,1)=mass_center(1)/sum(m);
       center(i,2)=mass_center(2)/sum(m);
    end
    
    rubber_ydistribution=zeros(nbiny,3);rubber_xdistribution=zeros(nbinx,3);
    sand_ydistribution=zeros(nbiny,3);sand_xdistribution=zeros(nbinx,3);
    for i=1:soft_num
        xindex=floor((center(i,1)-xstart)/dbinx)+1;
        yindex=floor((center(i,2)-ystart)/dbiny)+1;
        rubber_xdistribution(xindex,1)=rubber_xdistribution(xindex,1)+area(i);
        rubber_xdistribution(xindex,2)=rubber_xdistribution(xindex,2)+1;
        rubber_ydistribution(yindex,1)=rubber_ydistribution(yindex,1)+area(i);
        rubber_ydistribution(yindex,2)=rubber_ydistribution(yindex,2)+1;
    end

    for i=soft_num+1:particle_num
        xindex=floor((center(i,1)-xstart)/dbinx)+1;
        yindex=floor((center(i,2)-ystart)/dbiny)+1;
        sand_xdistribution(xindex,1)=sand_xdistribution(xindex,1)+area(i);
        sand_xdistribution(xindex,2)=sand_xdistribution(xindex,2)+1;
        sand_ydistribution(yindex,1)=sand_ydistribution(yindex,1)+area(i);
        sand_ydistribution(yindex,2)=sand_ydistribution(yindex,2)+1;
    end
    rubber_ydistribution(:,3)=rubber_ydistribution(:,1)./rubber_ydistribution(:,2);
    sand_ydistribution(:,3)=sand_ydistribution(:,1)./sand_ydistribution(:,2);
    rubber_xdistribution(:,3)=rubber_xdistribution(:,1)./rubber_xdistribution(:,2);
    sand_xdistribution(:,3)=sand_xdistribution(:,1)./sand_xdistribution(:,2);
    xnumber_fraction=rubber_xdistribution(:,2)./(rubber_xdistribution(:,2)+sand_xdistribution(:,2));
    ynumber_fraction=rubber_ydistribution(:,2)./(rubber_ydistribution(:,2)+sand_ydistribution(:,2));
    xvolume_fraction=rubber_xdistribution(:,3)./(rubber_xdistribution(:,3)+sand_xdistribution(:,3));
    yvolume_fraction=rubber_ydistribution(:,3)./(rubber_ydistribution(:,3)+sand_ydistribution(:,3));
    
    x_number_fraction(timestep,1:nbinx)=xnumber_fraction;
    x_volume_fraction(timestep,1:nbinx)=xvolume_fraction;
    y_number_fraction(timestep,1:nbiny)=ynumber_fraction;
    y_volume_fraction(timestep,1:nbiny)=yvolume_fraction;
end

x_number_std=sqrt(sum((x_number_fraction-0.3).^2,2)/size(x_number_fraction,2));
y_number_std=sqrt(sum((y_number_fraction-0.3).^2,2)/size(y_number_fraction,2));