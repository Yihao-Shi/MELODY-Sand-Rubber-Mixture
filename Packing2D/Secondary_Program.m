function [Cellules,Vertices,Historique,ODECS,Proprietes,SolidFraction,Surfaces,Angles,D50,Cu,PDFAnglesExp,PDFSurfacesExp,PDFElongsExp,PDFRoundExp,PDFCircExp,PDFRegulExp]=Secondary_Program(FileName,PFCFileCreation,PFCFileName,Nparticles,TargetSurfaceCOV,TargetMainOrientation,TargetAnisotropy,DomainPoints,NiterMax,TargetError,TargetSolidFraction,NvarOptim,COVSpectrum,TypeCOVSpectrum,TypeSpectrum,FileSpectrum,DescriptorD2,DescriptorD3,SpectrumDecay1,DescriptorD8,SpectrumDecay2,Dmax,Rmin,Pmax)
global D2 D3 Decay1 D8 Decay2 Rmax Ndesc Dtail NVAR FracCible Surface TypeSpectre FileSpectre COVSpectre Cn TypeCOV

Nom=FileName;
Npoints=Nparticles;
COVCible=TargetSurfaceCOV;
AngleMoy=TargetMainOrientation;
Rapport_Ani=TargetAnisotropy;
Domain=DomainPoints;
Niter=NiterMax;
ErreurCible=TargetError;
FracCible=TargetSolidFraction;
NVAR=NvarOptim;
COVSpectre=COVSpectrum;
TypeCOV=TypeCOVSpectrum;
TypeSpectre=TypeSpectrum;
FileSpectre=FileSpectrum;
D2=DescriptorD2;
D3=DescriptorD3;
Decay1=SpectrumDecay1;
D8=DescriptorD8;
Decay2=SpectrumDecay2;
dmax=Dmax;
rmin=Rmin;
pmax=Pmax;

rand('state',sum(100*clock))
if Rapport_Ani==1
    Rapport_Ani=0.9999;
end
if D2==0
    D2=10^-10;
end
if D3==0
    D3=10^-10;
end
if D8==0
    D8=10^-10;
end

xmin=min(Domain(:,1));
xmax=max(Domain(:,1));
ymin=min(Domain(:,2));
ymax=max(Domain(:,2));
Npoints=floor(Npoints/polyarea(Domain(:,1),Domain(:,2))*((xmax-xmin)*(ymax-ymin)));

if COVCible>0.4
    Sample=rand(Npoints,2);
else
    Sample=Halton1(Npoints,2);
end

Xini=xmin+Sample(:,1)*(xmax-xmin);
Yini=ymin+Sample(:,2)*(ymax-ymin);
X=[];
Y=[];
for i=1:size(Xini,1)
    if inpolygon(Xini(i,1),Yini(i,1),Domain(:,1),Domain(:,2))==1
        X=cat(1,X,[Xini(i,1)]);
        Y=cat(1,Y,[Yini(i,1)]);
    end
end
Npoints=size(X,1);   

[Xn,Yn,Surfaces,Elongations,Angles,Historique,Cellules,Vertices]=IMC(X,Y,Domain,COVCible,AngleMoy,Rapport_Ani,Niter,ErreurCible,Nom);
save(Nom)
Couleurs=rand(size(Cellules,1),3);
Disposition=zeros(Npoints,3);
Disposition(:,3)=ones(Npoints,1);
ODECS=cell(Npoints,1);
Proprietes=zeros(Npoints,23);
for i=1:Npoints
    i
    Cell=Vertices(transpose(Cellules{i,1}),:);
    x=Cell(:,1);y=Cell(:,2);
    K=convhull(Cell(:,1),Cell(:,2));Xcell=x(K,1);Ycell=y(K,1);
    [Cercles,Xconv,Yconv]=FillCell(D2,D3,Decay1,D8,Decay2,Xcell,Ycell,NVAR,FracCible,dmax,rmin,pmax);
    ODECS{i,1}=Cercles;
    
    [Angle,L,S,Excentricite,Surface,Xc,Yc,Inertie,Req]=ProprietesODEC(Cercles);
    [TetaContourfin,RContourfin,XContourfin,YContourfin,Perimetre]=ContourFin(Cercles,Xc,Yc,500);
    [Xconv,Yconv,PerimetreConv]=ContourConvexe(XContourfin-Xc,YContourfin-Yc);
    [Xinsc,Yinsc,Rinsc]=RayonInscrit(XContourfin-Xc,YContourfin-Yc);
    [Xcirc,Ycirc,Rcirc]=RayonCirconscrit(Xconv,Yconv);
    [Inertie1,Xc1,Yc1,Inertie2,Xc2,Yc2]=DistMass(Cercles,Surface);
    
    Proprietes(i,1)=size(Cercles,1);
    Proprietes(i,2)=Angle;
    Proprietes(i,3)=L;
    Proprietes(i,4)=S;
    Proprietes(i,5)=Excentricite; %Elongation
    Proprietes(i,6)=Surface;
    Proprietes(i,7)=Perimetre;
    Proprietes(i,8)=Xc;
    Proprietes(i,9)=Yc;
    Proprietes(i,10)=Inertie;
    Proprietes(i,11)=Req;
    Proprietes(i,12)=PerimetreConv;
    Proprietes(i,13)=Rinsc;
    Proprietes(i,14)=Rcirc;
    Proprietes(i,15)=Xc1;
    Proprietes(i,16)=Yc1;
    Proprietes(i,17)=Inertie1;
    Proprietes(i,18)=Xc2;
    Proprietes(i,19)=Yc2;
    Proprietes(i,20)=Inertie2;
    Proprietes(i,21)=mean(Cercles(1,3))/Rinsc; %Roundness
    Proprietes(i,22)=(Rinsc/Rcirc)^0.5; %Circularity
    Proprietes(i,23)=-log10((Perimetre-PerimetreConv)./Perimetre); %Regularity
end
[PDFAnglesExp,PDFSurfacesExp,PDFElongsExp,PDFRoundExp,PDFCircExp,PDFRegulExp,SolidFraction,VoidRatio,D10,D50,Cu]=Statistics(Proprietes,ODECS,Domain);

if PFCFileCreation==1
    Walls=zeros(size(Domain,1)-1,4);
    for i=1:size(Walls,1)
        Walls(i,1:2)=Domain(i,1:2);
        Walls(i,3:4)=Domain(i+1,1:2);
    end
    TextFilePFC(ODECS,Proprietes,Walls,PFCFileName)
end

TraceIMC(X,Y,Xn,Yn,COVCible,AngleMoy,Rapport_Ani,Surfaces,Angles,Elongations,Historique,Cellules,Vertices,Domain);
TraceurODECS(ODECS,Disposition,Cellules,Vertices,Couleurs,[xmin xmax ymin ymax],1);

