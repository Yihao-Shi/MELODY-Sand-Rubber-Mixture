clear all
CellsFileName='C:\Users\ncasas\Documents\MELODY_RESU\MODELE_2\M2-matrices\matrice seule\4x2\creation\data_matrix_4x2.asc';
% GrainsFileName='C:\CODES\Cvoro2D\Executables\Gouge_Matrix_Sample50_01.mat';
Criterion='Centre';
Window=[0,4,0,2];

[Cells,Vertices]=LoadCVoroDataFile(CellsFileName);
axis(Window);

%  plot(VERTICES(CELLS{i,1},1),VERTICES(CELLS{i,1},2),'-b')
for i=1:size(Cells)
    Contours{i,1}(:,1)=Vertices(Cells{i,1},1);
    Contours{i,1}(:,2)=Vertices(Cells{i,1},2);
end
print('Cells.png','-dpng','-r300')
save('Matrix.mat','Cells','Vertices','Contours')


% load(GrainsFileName,'Contours');
% figure;hold on
% for i=1:size(Contours,1)
%     plot(Contours{i,1}(:,1),Contours{i,1}(:,2),'-b')
% end
% axis equal;axis(Window);
% print('Grains.png','-dpng','-r300')

% Status=zeros(size(Cells,1),1);
% for i=1:size(Cells,1)
%     xcel=Vertices(Cells{i,1},1);
%     ycel=Vertices(Cells{i,1},2);
%     if strcmp(Criterion,'Centre')
%         xc=mean(xcel(1:end-1));
%         yc=mean(ycel(1:end-1));
% 
%     end
% end


% MergedCells={};
% for i=1:size(Contours,1)
%     liscells=find(Status==i);
%     segments=[];
%     for j=1:length(liscells)
%         cel=Cells{liscells(j),1};
%         for k=1:size(cel,1)-1
%             segments=[segments;[cel(k),cel(k+1)]];
%         end
%     end
%     s=sort(segments')';
%     segments(:,3)=1e6*s(:,1)+s(:,2);
%     for j=1:size(segments,1)
%         segments(j,4)=size(find(segments(:,3)==segments(j,3)),1);
%     end
%     s=segments(find(segments(:,4)==1),1:2);
%     segments=s(1,:);
%     while size(segments,1)<size(s,1)
%         segments=[segments;s(find(s(:,1)==segments(end,2)),:)];
%     end
%     MergedCells{i,1}=Vertices([segments(:,1);segments(1,1)],:);
% end
% n=size(MergedCells,1);
% for i=1:size(Status)
%     if Status(i)==0
%         n=n+1;
%         MergedCells{n,1}=Vertices(Cells{i,1},:);
%     end
% end
%  save('Matrix.mat','Cells')

% figure;hold on
% for i=1:size(MergedCells,1)
%     plot(MergedCells{i,1}(:,1),MergedCells{i,1}(:,2),'-b')
% end
% axis equal;axis(Window);
% print('Merged.png','-dpng','-r300')
