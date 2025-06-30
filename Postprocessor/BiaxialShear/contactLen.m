clear;
proportion='10';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/' proportion 'percent/'];
contact=load(['source' proportion '.mat']).contacts;

%% extract data
particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;
[~,countorNode]=NodeInfo(FRACTION,particle_num,soft_num);

lengthContact=zeros(36,1);lengthSS=zeros(36,1);lengthSR=zeros(36,1);lengthRR=zeros(36,1);
for timestep=36%1:36
    SS=[];SR=[];RR=[];
    for i=3498%1:particle_num
        particle_contact=contact{timestep,i};
        if size(particle_contact,1)>0
            [masterPar,mrows]=unique(sort(particle_contact(:,4)));mrows=[mrows;size(particle_contact,1)+1];
            particle_contact=sortrows(particle_contact,4);
            for j=1:length(masterPar)
                clen=0;
                for k=mrows(j):mrows(j+1)-1
                   if sqrt(particle_contact(k,16)^2+particle_contact(k,17)^2)>0
                      clen=clen+particle_contact(k,15);
                   end
                end
                if clen>0
                    if i-1<soft_num&&masterPar(j)<soft_num
                        SS=[SS;i-1,masterPar(j),clen];
                    elseif (i-1<soft_num&&masterPar(j)>=soft_num)||(i-1>=soft_num&&masterPar(j)<soft_num)
                        SR=[SR;i-1,masterPar(j),clen];
                    else
                        RR=[RR;i-1,masterPar(j),clen];
                    end
                end
            end
        end
    end
    
    total=[];
    if isempty(SS)==0
        total=[total;SS(:,3)];
        lengthSS(timestep)=mean(SS(:,3));
    end
    if isempty(SR)==0
        total=[total;SR(:,3)];
        lengthSR(timestep)=mean(SR(:,3));
    end
    if isempty(RR)==0
        total=[total;RR(:,3)];
        lengthRR(timestep)=mean(RR(:,3));
    end
    lengthContact(timestep)=mean(total);
end