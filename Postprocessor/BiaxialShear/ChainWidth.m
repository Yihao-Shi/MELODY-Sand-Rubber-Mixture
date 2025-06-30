clear
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
masters=dlmread([FRACTION 'master_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
xposs=dlmread([FRACTION 'posx_info.asc'],' ');
yposs=dlmread([FRACTION 'posy_info.asc'],' ');


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
for timestep=1:36
    xpos=xposs(timestep,3:end-1).';ypos=yposs(timestep,3:end-1).';
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);

    fxPerSec=zeros(particle_num,max(countorNode));
    xposPerSec=zeros(particle_num,1);
    fyPerSec=zeros(particle_num,max(countorNode));
    yposPerSec=zeros(particle_num,1);
    xnormPerSec=zeros(particle_num,max(countorNode));
    ynormPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        xposPerSec(i)=xpos(i);
        yposPerSec(i)=ypos(i);
        for j=1:countorNode(i)
            fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
            fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
            xnormPerSec(i,j)=xnorm(timestep,sum(countorNode(1:i-1))+j);
            ynormPerSec(i,j)=ynorm(timestep,sum(countorNode(1:i-1))+j);
        end
    end


    fc=zeros(size(fxPerSec,1)*size(fxPerSec,2),9);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>=0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=min(i-1,masterPerSec(i,j));
                fc(sum(countorNode(1:i-1))+j,2)=max(i-1,masterPerSec(i,j));
                fc(sum(countorNode(1:i-1))+j,3)=xnormPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,4)=ynormPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,5)=fxPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,6)=fyPerSec(i,j);
                if i-1>masterPerSec(i,j)
                    fc(sum(countorNode(1:i-1))+j,3)=-fc(sum(countorNode(1:i-1))+j,3);
                    fc(sum(countorNode(1:i-1))+j,4)=-fc(sum(countorNode(1:i-1))+j,4);
                    fc(sum(countorNode(1:i-1))+j,5)=-fc(sum(countorNode(1:i-1))+j,5);
                    fc(sum(countorNode(1:i-1))+j,6)=-fc(sum(countorNode(1:i-1))+j,6);
                end
            end
        end
    end
    fc(fc(:,1)==0&fc(:,2)==0,:)=[];fc=sortrows(fc,1);

    temp=zeros(size(fc,1),5);count=1;
    [slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
    for i=1:length(slavePar)
        test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
        [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
        for j=1:length(masterPar)
            temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
            temp(count,3)=sum(test(mrows(j):mrows(j+1)-1,3));
            temp(count,4)=sum(test(mrows(j):mrows(j+1)-1,4));
            normalize=sqrt(temp(count,3)^2+temp(count,4)^2);
            temp(count,3)=temp(count,3)/normalize;
            temp(count,4)=temp(count,4)/normalize;
            temp(count,5)=sum(test(mrows(j):mrows(j+1)-1,5));
            temp(count,6)=sum(test(mrows(j):mrows(j+1)-1,6));
            count=count+1;
        end
    end
    temp(count-1:end,:)=[];
    temp(:,7)=temp(:,3).*temp(:,5)+temp(:,4).*temp(:,6);

    for i=1:size(temp,1)
        temp(i,8)=xpos(temp(i,1)+1);
        temp(i,9)=ypos(temp(i,1)+1);
        temp(i,10)=xpos(temp(i,2)+1);
        temp(i,11)=ypos(temp(i,2)+1);
    end

    for i=1:size(temp,1)
        if temp(i,1)<soft_num&&temp(i,2)<soft_num
            temp(i,12)=2;
        elseif temp(i,1)<soft_num&&temp(i,2)>=soft_num||temp(i,1)>=soft_num&&temp(i,2)<soft_num
            temp(i,12)=1;
        else
            temp(i,12)=0;
        end
    end

    %% Output
    fid = fopen([FRACTION 'linewidth' num2str(timestep) '.txt'],'w');
    for i=1:size(temp,1)
        fprintf(fid,'%f %f %f %f %f %f\r\n',temp(i,8),temp(i,9),temp(i,10),temp(i,11),temp(i,7),temp(i,12));
    end
    fclose(fid); 
end
