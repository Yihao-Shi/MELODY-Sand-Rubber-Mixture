function [X,Y,Surfaces,Elongations,Angles,Historique,Cellules,Vertices]=IMC(X,Y,Domain,COVCible,AngleMoy,Rapport_Ani,Niter,ErreurCible,Nom)
global x y
Stot=polyarea(Domain(:,1),Domain(:,2));
Npoints=size(X,1);
Ntranches=30;
Largeur=3;
MoyCible=Stot/Npoints;
StdCible=COVCible*MoyCible;
[MuN,SigN]=Gauss_to_Log(MoyCible,StdCible);
Width=(2*Largeur*SigN)/Ntranches;
Edges=transpose(MuN-Largeur*SigN:Width:MuN+Largeur*SigN);
Centres=transpose(MuN-Largeur*SigN+Width/2:Width:MuN+Largeur*SigN-Width/2);
Cumul=normcdf(Edges,MuN,SigN);
HistoLogCible=zeros(size(Centres,1),1);
for i=1:size(Centres,1)
    HistoLogCible(i,1)=Cumul(i+1,1)-Cumul(i,1);
end
HistoAngleCible=zeros(18,1);
for i=1:18
    angle=((i-0.5)*10-90)/180*3.14159;
    HistoAngleCible(i,1)=(1+Rapport_Ani)/(1-Rapport_Ani)+cos((angle-AngleMoy/180*3.14159)*2);
end
HistoAngleCible=HistoAngleCible/sum(HistoAngleCible);

[X,Y,Vertices,Cellules]=VoronoiDomain2D(X,Y,Domain);
[Surfaces,Angles,Elongations,Voisins]=PropVoronoi2D(X,Y,Vertices,Cellules);
Npoints=size(X,1);
Ncellules=Npoints;

HistoLogExp=transpose(hist(log(Surfaces),Centres))/Ncellules;
HistoLogIni=HistoLogExp;
Erreur1=sqrt(mean((HistoLogCible-HistoLogExp).^2)); %Erreur sur les surfaces
%Erreur1=max(abs(HistoLogCible-HistoLogExp));
HistoAngleExp=transpose(hist(Angles,[-85:10:85]))/Ncellules;
HistoAngleIni=HistoAngleExp;
Erreur2=sqrt(mean((HistoAngleCible-HistoAngleExp).^2)); %Erreur sur les angles
%Erreur2=max(abs(HistoAngleCible-HistoAngleExp));
%Erreur=(Erreur1+Erreur2)/2
Erreur=max([Erreur1,Erreur2])

Historique=[1,Erreur1,Erreur2,Erreur];
Compteur=2;
while Erreur>ErreurCible & Compteur<Niter
    NumPoint=randint(1,1,[1,Npoints]);
    xmi=min(Vertices(transpose(Cellules{NumPoint,1}),1));
    ymi=min(Vertices(transpose(Cellules{NumPoint,1}),2));
    xma=max(Vertices(transpose(Cellules{NumPoint,1}),1));
    yma=max(Vertices(transpose(Cellules{NumPoint,1}),2));
    nvois=size(Voisins{NumPoint,1},2);
    for i=1:nvois
        if min(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),1))<xmi
            xmi=min(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),1));
        end
        if min(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),2))<ymi
            ymi=min(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),2));
        end
        if max(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),1))>xma
            xma=max(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),1));
        end
        if max(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),2))>yma
            yma=max(Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),2));
        end
    end
    Xco=Vertices(transpose(Cellules{NumPoint,1}),1);
    Xco=cat(1,Xco,Xco(1,1));
    Yco=Vertices(transpose(Cellules{NumPoint,1}),2);
    Yco=cat(1,Yco,Yco(1,1));
    OK=0;
    while OK==0
        xtir=xmi+rand*(xma-xmi);
        ytir=ymi+rand*(yma-ymi);
        if inpolygon(xtir,ytir,Xco,Yco)==1
            OK=1;
            break
        end
        for i=1:nvois
            Xvois=Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),1);
            Xvois=cat(1,Xvois,Xvois(1,1));
            Yvois=Vertices(transpose(Cellules{Voisins{NumPoint,1}(1,i),1}),2);
            Yvois=cat(1,Yvois,Yvois(1,1));
            if inpolygon(xtir,ytir,Xvois,Yvois)==1
                OK=1;
                break
            end
        end
    end
    Angle=AnglePolaire(xtir-X(NumPoint,1),ytir-Y(NumPoint,1));
    Distance=sqrt((xtir-X(NumPoint,1))^2+(ytir-Y(NumPoint,1))^2);
    [Xn,Yn,VoisinsN,VerticesN,CellulesN,SurfacesN,AnglesN,ElongationsN]=UpdateVoronoi(X,Y,Domain,NumPoint,Distance,Angle,Voisins,Vertices,Cellules,Surfaces,Angles,Elongations);
        
    HistoLogExp=transpose(hist(log(SurfacesN),Centres))/Ncellules;
    ErreurN1=sqrt(mean((HistoLogCible-HistoLogExp).^2)); %Erreur sur les surfaces
    %ErreurN1=max(abs(HistoLogCible-HistoLogExp));
    HistoAngleExp=transpose(hist(AnglesN,[-85:10:85]))/Ncellules;
    ErreurN2=sqrt(mean((HistoAngleCible-HistoAngleExp).^2)); %Erreur sur les angles
    %ErreurN2=max(abs(HistoAngleCible-HistoAngleExp));
    %ErreurN=(ErreurN1+ErreurN2)/2;
    ErreurN=max([ErreurN1,ErreurN2]);
    
    if ErreurN<Erreur
        Compteur
        X=Xn;
        Y=Yn;
        Erreur=ErreurN
        Voisins=VoisinsN;
        Vertices=VerticesN;
        Cellules=CellulesN;
        Surfaces=SurfacesN;
        Angles=AnglesN;
        Elongations=ElongationsN;
        Historique(size(Historique,1)+1,1)=Compteur;
        Historique(size(Historique,1),2)=ErreurN1;
        Historique(size(Historique,1),3)=ErreurN2;
        Historique(size(Historique,1),4)=ErreurN;
        %save(Nom,'Historique')
        %save([Nom ' ' int2str(Compteur) '.mat'])
        %save(Nom)
    end
    Compteur=Compteur+1;
end

[X,Y,Vertices,Cellules]=VoronoiDomain2D(X,Y,Domain);
[Surfaces,Angles,Elongations,Voisins]=PropVoronoi2D(X,Y,Vertices,Cellules);
