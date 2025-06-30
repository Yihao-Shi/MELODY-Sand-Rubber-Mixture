clear
proportion='30';
FRACTION=['../../Simulation_DATA/BiaxialShearTest/StrainEnergy/' proportion 'percent/'];
output_path=[FRACTION 'GreenStrain_info/'];
f_ini=70;
f_num=35;
id=[f_ini:f_ini+f_num];
for k=1:length(id)
    if id(k)<10
        file{k}=[FRACTION 'GRAINS_0000' num2str(id(k)) '.vtk'];
    elseif id(k)<100&&id(k)>=10
        file{k}=[FRACTION 'GRAINS_000' num2str(id(k)) '.vtk'];
    elseif id(k)>=100
        file{k}=[FRACTION 'GRAINS_00' num2str(id(k)) '.vtk'];
    end
end

particle_num=5025;
rigid_num=ceil((1-str2num(proportion)/100)*particle_num);
soft_num=particle_num-rigid_num;

[softFieldNode,~]=NodeInfo(FRACTION,particle_num,soft_num);
%% Green-Lagrange Strain
for k=1:length(id)
    data=fopen(file{k});
    while ~feof(data)
        tline=fgetl(data);
        if strfind(tline,'SCALARS Green-Lagrange_XX_Strain float 1')~=0 
            strainxx=zeros(sum(softFieldNode),1);
            tline=fgetl(data);
            for i=1:size(strainxx)
                strainxx(i)=str2double(fgetl(data));
            end
            fid = fopen([output_path 'GreenXX' num2str(id(k)-f_ini) '.txt'],'w');
            fprintf(fid,'%f\r\n',strainxx(:,1));
            fclose(fid); 
        end
        
        if strfind(tline,'SCALARS Green-Lagrange_YY_Strain float 1')~=0 
            strainyy=zeros(sum(softFieldNode),1);
            tline=fgetl(data);
            for i=1:size(strainyy)
                strainyy(i)=str2double(fgetl(data));
            end
            fid = fopen([output_path 'GreenYY' num2str(id(k)-f_ini) '.txt'],'w');
            fprintf(fid,'%f\r\n',strainyy(:,1));
            fclose(fid); 
        end
        
        if strfind(tline,'SCALARS Green-Lagrange_XY_Strain float 1')~=0 
            strainxy=zeros(sum(softFieldNode),1);
            tline=fgetl(data);
            for i=1:size(strainxy)
                strainxy(i)=str2double(fgetl(data));
            end
            fid = fopen([output_path 'GreenXY' num2str(id(k)-f_ini) '.txt'],'w');
            fprintf(fid,'%f\r\n',strainxy(:,1));
            fclose(fid); 
            break;
        end
    end
    fclose(data);
end