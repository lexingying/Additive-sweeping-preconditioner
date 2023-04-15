K=8; %number of waves
NPW=8; % numer of grids per wave
NPML=9; % number of PMLs at the boundary of the domain
CPML=10; % PML parameter for boundary PMLs
NPAD=5; % number of PMLs for each auxiliary subdomain
CPAD=4; % PML parameter for auxiliary PMLs
NLPD=NPAD-1; % number of layers per domain
leafsize=4; % leaf size in nested dissection algorithm

NC = K*NPW; %number of cells
h = 1/NC;
N=NC-1;N1=N;N2=N;N3=N;

[xx,yy,zz]=ndgrid(h*(1:N));

% c is the velocity field.
% Here we initialize a converging lens as an example.
c=ones(size(xx))-0.4*exp(-32*((xx-1/2).^2+(yy-1/2).^2+(zz-1/2).^2));
cmid=(min(c(:))+max(c(:)))/2;c=c/cmid;

ksq=(2*pi*K./c).^2;

gs=(1/NPML)*(0.5:0.5:NPML-0.5);eta=NPML*h;
sigR=CPML/eta*gs.^2;
sR=1./(1+1i*sigR/K);sL=sR(end:-1:1);
s1=[sL,ones(1,2*(N1-2*(NPML-1))-1),sR];
s2=[sL,ones(1,2*(N2-2*(NPML-1))-1),sR];
s3=[sL,ones(1,2*(N3-2*(NPML-1))-1),sR];

gs=(1/NPAD)*(0.5:0.5:NPAD-0.5);eta=NPAD*h;
sigR=CPAD/eta*gs.^2;
pR=1./(1+1i*sigR/K);pL=pR(end:-1:1);

A=setupA3D(h,ksq,s1,s2,s3);

P=setup1(NPML,NLPD,NPAD,pL,pR,h,ksq,s1,s2,s3,leafsize);

funcM=@(inc)apply1(NPAD,P,leafsize,inc);

% f is the source.
% Here we use a Guassian point source at the center of the domain
% as an example
f=exp(-(8*K).^2*((xx-1/2).^2+(yy-1/2).^2+(zz-1/2).^2));

% GMRES solver
[u,flag,relres,iter,resvec]=gmres(A,f(:),40,1e-3,2,funcM);