clear;
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
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
   fx(i,:)=fxs((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   fy(i,:)=fys((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   xnorm(i,:)=xnorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
   ynorm(i,:)=ynorms((i-1)*(sum(countorNode)+2)+3:i*(sum(countorNode)+2),1).'; 
end

%%
rigidNeighbor=zeros(1,36);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,~]=master_process(particle_num,countorNode,master,fn,timestep);
    
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
            temp(row,colume)=masterPar(j)+1;
            colume=colume+1;
        end
        row=row+1;
    end
    temp(row-1:end,:)=[];
    
    t=1;vec=find(temp(:,1)<soft_num);
    for row=1:length(vec)
        c=0;a=0;i=vec(row);
        for j=2:size(temp,2)
            if temp(i,j)>0
                a=a+1;
            end
            if temp(i,j)>soft_num 
                c=c+1;
            end
        end
        if c==a
            rigidNeighbor(t,timestep)=temp(i,1);
            t=t+1;
        end
    end
end