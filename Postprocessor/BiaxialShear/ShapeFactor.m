clear;
proportion='10';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

for i=1:36
   master(i,:)=masters((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locx(i,:)=locxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   locy(i,:)=locys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
circularity=zeros(36,1);
c3=zeros(36,5);c6=zeros(36,5);c8=zeros(36,5);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,~]=master_process(particle_num,countorNode,master,fn,timestep);
    
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    
    area=zeros(particle_num,1);
    perimeter=zeros(particle_num,1);
    for i=1:particle_num
       x=locxPerSec(i,1:countorNode(i));y=locyPerSec(i,1:countorNode(i));
       for j=1:countorNode(i)-2
           area(i)=area(i)+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
       end
       xp=[locxPerSec(i,1:countorNode(i)),locxPerSec(i,1)];yp=[locyPerSec(i,1:countorNode(i)),locyPerSec(i,1)];
       for j=2:countorNode(i)+1
           perimeter(i)=perimeter(i)+sqrt((xp(j)-xp(j-1))^2+(yp(j)-yp(j-1))^2);
       end
    end
    peri_circle=2*sqrt(pi*area);
    cir=peri_circle./perimeter;
    
    fc=zeros(size(masterPerSec,1)*size(masterPerSec,2),2);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
            end
        end
    end
    fc(fc(:,1)==0&fc(:,2)==0,:)=[];
    
    temp=zeros(size(fc,1),14);row=1;
    [slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
    for i=1:length(slavePar)
        test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
        [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
        temp(row,1)=slavePar(i);colume=2;
        for j=1:length(masterPar)
            temp(row,colume)=masterPar(j);
            colume=colume+1;
        end
        temp(row,10)=cir(temp(row,1)+1);
        row=row+1;
    end
    temp(row-1:end,:)=[];
    circularity(timestep)=mean(temp(1:soft_num,10));
end