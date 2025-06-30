function [CELLS,VERTICES]=LoadCVoroDataFile(FileName)
%clear all
%FileName='C:\CODES\Cvoro2D\bin\Release\data_end.asc';
DataFile=fopen(FileName);
Data = textscan(DataFile,'%s');
Ncells=sscanf(Data{1,1}{3,1},'%f');
Nvertices=sscanf(Data{1,1}{5,1},'%f');

CELLS=cell(Ncells,1);
VERTICES=zeros(Nvertices,2);

for Line=1:size(Data{1,1},1)
    if strcmp(Data{1,1}{Line,1},'CELLS')
        break
    end
end
Line=Line+2;
for i=1:Ncells
    nvert=sscanf(Data{1,1}{Line,1},'%f');
    cel=zeros(nvert,1);
    for j=1:nvert+1
        Line=Line+1;
        cel(j,1)=sscanf(Data{1,1}{Line,1},'%f')+1;
    end
    CELLS{i,1}=cel;
    Line=Line+nvert+2;
end

Line=Line+1;
for i=1:Nvertices
    VERTICES(i,1)=sscanf(Data{1,1}{Line,1},'%f');
    Line=Line+1;
    VERTICES(i,2)=sscanf(Data{1,1}{Line,1},'%f');
    Line=Line+2;
end

fclose(DataFile);

figure;hold on;
for i=1:Ncells
    plot(VERTICES(CELLS{i,1},1),VERTICES(CELLS{i,1},2),'-b')
end
axis equal