clear;
proportion='0';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];                      
masters=dlmread([FRACTION 'master_info.asc'],' ');
lens=dlmread([FRACTION 'length_info.asc'],' ');
fxs=dlmread([FRACTION 'fx_info.asc'],' ');
fys=dlmread([FRACTION 'fy_info.asc'],' ');
xnorms=dlmread([FRACTION 'xnorm_info.asc'],' ');
ynorms=dlmread([FRACTION 'ynorm_info.asc'],' ');
locxs=dlmread([FRACTION 'locx_info.asc'],' ');
locys=dlmread([FRACTION 'locy_info.asc'],' ');
xpos=dlmread([FRACTION 'posx_info.asc'],' ');
ypos=dlmread([FRACTION 'posy_info.asc'],' ');
xpos=xpos(:,3:end-1);ypos=ypos(:,3:end-1);

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
bucklingList=cell(1,1);particle=cell(1,1);partInfo=cell(1,1);
for timestep=1:36
    fn=fx.*xnorm+fy.*ynorm;
    [masterPerSec,fnPerSec]=master_process(particle_num,countorNode,master,fn,timestep);
    
    fxPerSec=zeros(particle_num,max(countorNode));
    fyPerSec=zeros(particle_num,max(countorNode));
    xposPerSec=zeros(particle_num,max(countorNode));
    yposPerSec=zeros(particle_num,max(countorNode));
    locxPerSec=zeros(particle_num,max(countorNode));
    locyPerSec=zeros(particle_num,max(countorNode));
    for i=1:particle_num
        for j=1:countorNode(i)
            fxPerSec(i,j)=fx(timestep,sum(countorNode(1:i-1))+j);
            fyPerSec(i,j)=fy(timestep,sum(countorNode(1:i-1))+j);
            xposPerSec(i,j)=xpos(timestep,i);
            yposPerSec(i,j)=ypos(timestep,i);
            locxPerSec(i,j)=locx(timestep,sum(countorNode(1:i-1))+j);
            locyPerSec(i,j)=locy(timestep,sum(countorNode(1:i-1))+j);
        end
    end
    
    area=zeros(particle_num,1);
    for i=1:particle_num
       x=locxPerSec(i,:);y=locyPerSec(i,:);
       for j=1:countorNode(i)-2
           area(i)=area(i)+0.5*(x(1)*y(j+1)-y(1)*x(j+1)+x(j+1)*y(j+2)-y(j+1)*x(j+2)+x(j+2)*y(1)-y(j+2)*x(1));
       end
    end
    
    fc=zeros(size(fxPerSec,1)*size(fxPerSec,2),10);
    for i=1:particle_num
        for j=1:countorNode(i)
            if masterPerSec(i,j)>=0&&masterPerSec(i,j)<5025
                fc(sum(countorNode(1:i-1))+j,1)=i-1;fc(sum(countorNode(1:i-1))+j,2)=masterPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,3)=fxPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,4)=fyPerSec(i,j);
                fc(sum(countorNode(1:i-1))+j,5)=locxPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,j)-xposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1);
                fc(sum(countorNode(1:i-1))+j,6)=locyPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,j)-yposPerSec(fc(sum(countorNode(1:i-1))+j,1)+1,1);
                dist=sqrt(fc(sum(countorNode(1:i-1))+j,5)^2+fc(sum(countorNode(1:i-1))+j,6)^2);
                fc(sum(countorNode(1:i-1))+j,7)=fc(sum(countorNode(1:i-1))+j,3)*fc(sum(countorNode(1:i-1))+j,5)/area(i);
                fc(sum(countorNode(1:i-1))+j,8)=fc(sum(countorNode(1:i-1))+j,3)*fc(sum(countorNode(1:i-1))+j,6)/area(i);
                fc(sum(countorNode(1:i-1))+j,9)=fc(sum(countorNode(1:i-1))+j,4)*fc(sum(countorNode(1:i-1))+j,6)/area(i);
                fc(sum(countorNode(1:i-1))+j,10)=fc(sum(countorNode(1:i-1))+j,4)*fc(sum(countorNode(1:i-1))+j,5)/area(i);
            end
        end
    end
    fc(fc(:,1)==0&fc(:,2)==0,:)=[];
    
    temp=zeros(size(fc,1),5);count=1;
    [slavePar,srows]=unique(fc(:,1));srows(length(srows)+1)=size(fc,1)+1;
    for i=1:length(slavePar)
        test=sortrows(fc(srows(i):srows(i+1)-1,:),2);
        [masterPar,mrows]=unique(test(:,2));mrows(length(mrows)+1)=srows(i+1)-srows(i)+1;
        for j=1:length(masterPar)
            temp(count,1)=slavePar(i);temp(count,2)=masterPar(j);
            temp(count,3)=sum(test(mrows(j):mrows(j+1)-1,7));
            temp(count,4)=(sum(test(mrows(j):mrows(j+1)-1,8))+sum(test(mrows(j):mrows(j+1)-1,9)))/2;
            temp(count,5)=sum(test(mrows(j):mrows(j+1)-1,10));
            count=count+1;
        end
    end
    temp(count-1:end,:)=[];
    
    [partNum,prows]=unique(temp(:,1));prows(length(prows)+1)=size(temp,1)+1;
    partList=zeros(size(partNum,1),9);
    for i=1:length(partNum)
        partList(i,1)=partNum(i);
        partList(i,2)=sum(temp(prows(i):prows(i+1)-1,3));
        partList(i,3)=sum(temp(prows(i):prows(i+1)-1,4));
        partList(i,4)=sum(temp(prows(i):prows(i+1)-1,5));
        partList(i,5)=0.5*(partList(i,2)+partList(i,4)+sqrt((partList(i,2)-partList(i,4))^2+4*partList(i,3)^2));
        partList(i,6)=0.5*(partList(i,2)+partList(i,4)-sqrt((partList(i,2)-partList(i,4))^2+4*partList(i,3)^2));
        partList(i,7)=0.5*atan((2*partList(i,3))/(partList(i,2)-partList(i,4)))+pi/2;
        partList(i,8)=cos(partList(i,7));
        partList(i,9)=sin(partList(i,7));
    end
    partInfo{timestep,1}=partList;
    partList(partList(:,5)<=mean(partList(:,5)),:)=[]; % Filter out all those particles less than or equal to the average
    
    for i=1:length(partList)
        for j=prows(partNum==partList(i,1)):prows(find(partNum==partList(i,1))+1)-1
            flag=ismember(temp(j,2),partList(:,1));
            partList(i,10)=partList(i,10)+flag;
        end
    end
    
    c=find(partList(:,10)==1);
    for i=1:length(c)
        d=temp(partList(c(i),1)==temp(:,1),2);
        for j=1:length(d)
            if partList(d(j)==partList(:,1),10)<=1
                partList(c(i),10)=0;
            end
        end
    end
    partList(partList(:,10)==0,:)=[];

    contList=zeros(sum(partList(:,10)),10);
    for i=1:size(partList,1)
        actCount=1;
        for j=prows(partNum==partList(i,1)):prows(find(partNum==partList(i,1))+1)-1
            flag=ismember(temp(j,2),partList(:,1));
            if flag==1
                contList(sum(partList(1:i-1,10))+actCount,1:2)=temp(j,1:2);
                xSlave=xposPerSec(temp(j,1)+1,1);ySlave=yposPerSec(temp(j,1)+1,1);
                xMaster=xposPerSec(temp(j,2)+1,1);yMaster=yposPerSec(temp(j,2)+1,1);
                contList(sum(partList(1:i-1,10))+actCount,3)=(xMaster-xSlave)/sqrt((xMaster-xSlave)^2+(yMaster-ySlave)^2);
                contList(sum(partList(1:i-1,10))+actCount,4)=(yMaster-ySlave)/sqrt((xMaster-xSlave)^2+(yMaster-ySlave)^2);
                contList(sum(partList(1:i-1,10))+actCount,5:6)=partList(i,8:9);
                contList(sum(partList(1:i-1,10))+actCount,7)=partList(partList(:,1)==contList(sum(partList(1:i-1,10))+actCount,2),5);
                actCount=actCount+1;
            end
        end
    end
    
    row=1;
    while row<=size(contList,1)
        if ~ismember(contList(i,1),contList(:,2))|~ismember(contList(i,2),contList(:,1))
            contList(row,:)=[];
        else
            row=row+1;
        end
    end
    
    for i=1:size(contList,1)
        masterRow=find(contList(:,1)==contList(i,2)&contList(:,2)==contList(i,1));
        if isempty(masterRow)==0
            branVec=[contList(i,3),contList(i,4)];branVecTrans=[contList(masterRow,3),contList(masterRow,4)];
            sigma=[contList(i,5),contList(i,6)];sigmaTrans=[contList(masterRow,5),contList(masterRow,6)];
            alpha=dot(sigma,branVec);alphaTrans=dot(sigmaTrans,branVecTrans);
            contList(i,8)=alpha;
            if abs(alpha)>cos(pi/4)&&abs(alpha)<=1&&abs(alphaTrans)>cos(pi/4)&&abs(alphaTrans)<=1
                contList(i,9)=contList(i,1)+contList(i,2);
            end
        end
    end
    contList(contList(:,9)==0,:)=[];contList=sortrows(contList,1);
    
    pid=unique(contList(:,1));
    for i=1:length(pid)
       if length(contList(pid(i)==contList(:,1),1))==1
           contList(pid(i)==contList(:,1),:)=[];
       end
    end
    
    p=unique(contList(:,1:2));
    for i=1:length(p)
        particle{i,timestep}=[p(i),partList(p(i)==partList(:,1),5)];
    end
        
%     for i=1:size(p,1)
%        X=locxPerSec(p(i,1)+1,1:countorNode(p(i,1)+1)); 
%        Y=locyPerSec(p(i,1)+1,1:countorNode(p(i,1)+1)); 
%        plot(X,Y);hold on;axis equal;
%     end
%     
%     for i=1:size(contList,1)
%         x1=mean(locxPerSec(contList(i,1)+1,1:countorNode(contList(i,1)+1)));
%         x2=mean(locxPerSec(contList(i,2)+1,1:countorNode(contList(i,2)+1)));
%         y1=mean(locyPerSec(contList(i,1)+1,1:countorNode(contList(i,1)+1)));
%         y2=mean(locyPerSec(contList(i,2)+1,1:countorNode(contList(i,2)+1)));
%         plot([x1,x2],[y1,y2]);hold on;
%     end
    
    
    numBucklingUnit=1;
    for i=1:size(contList)
        searchList=[contList(i,1),contList(i,2)];
        rowBuckling=find(contList(:,1)==contList(i,2));
        for j=1:length(rowBuckling)
            if contList(rowBuckling(j),2)~=contList(i,1)
                bucklingList{numBucklingUnit,2*timestep-1}=[searchList,contList(rowBuckling(j),2)];
                bucklingList{numBucklingUnit,2*timestep}=acos(contList(i,3)*contList(rowBuckling(j),3)+contList(i,4)*contList(rowBuckling(j),4))/pi*180;
                numBucklingUnit=numBucklingUnit+1;
            end
        end
    end
end