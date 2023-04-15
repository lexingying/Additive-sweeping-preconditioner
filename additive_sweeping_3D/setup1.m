function[P]=setup1(NPML,NLPD,NPAD,pL,pR,h,ksq,s1,s2,s3,chunksize)
[N1,N2,N3]=size(ksq);

%generate the partitions
NP=ceil((N3-2*(NPML-1))/NLPD);
pttn=cell(NP,1);
pttn{1}=[1,(NPML-1)+NLPD];
for g=2:NP-1
    pttn{g}=(NPML-1)+(g-1)*NLPD+[1,NLPD];
end
pttn{NP}=[(NPML-1)+(NP-1)*NLPD+1,N3];

    %construct matrices
    PM=cell(NP,1);
    if(1)
        b=1;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm:to+(NPAD-1));
        s3new=[s3(2*fm-1:2*to),pR];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PM{b}=db_setup(N1,to-fm+(NPAD-1)+1,chunksize,tmp);
        fprintf('complete: PM fm=%d to=%d\n',fm,to);
    end
    for b=2:NP-1
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm-(NPAD-1):to+(NPAD-1));
        s3new=[pL,s3(2*fm:2*to),pR];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PM{b}=db_setup(N1,to-fm+2*(NPAD-1)+1,chunksize,tmp);
        fprintf('complete: PM fm=%d to=%d\n',fm,to);
    end
    if(1)
        b=NP;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm-(NPAD-1):to);
        s3new=[pL,s3(2*fm:2*to+1)];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PM{b}=db_setup(N1,to-fm+(NPAD-1)+1,chunksize,tmp);
        fprintf('complete: PM fm=%d to=%d\n',fm,to);
    end
    
    %construct forward sweeping matrices
    PF=cell(NP,1);
    for b=2:NP-1
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm:to+(NPAD-1));
        s3new=[s3(2*fm-1:2*to),pR];
        tmp = setupA3D(h,ksqnew,s1,s2,s3new);
        PF{b}=db_setup(N1,to-fm+(NPAD-1)+1,chunksize,tmp);
        fprintf('complete: PF fm=%d to=%d\n',fm,to);
    end
    if(1)
        b=NP;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm:to);
        s3new=s3(2*fm-1:2*to+1);
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PF{b}=db_setup(N1,to-fm+1,chunksize,tmp);
        fprintf('complete: PF fm=%d to=%d\n',fm,to);
    end

    %construct backward sweeping matrices
    PB=cell(NP,1);
    for b=NP-1:-1:2
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm-(NPAD-1):to);
        s3new=[pL,s3(2*fm:2*to+1)];
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PB{b}=db_setup(N1,to-fm+(NPAD-1)+1,chunksize,tmp);
        fprintf('complete: PB fm=%d to=%d\n',fm,to);
    end
    if(1)
        b=1;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(:,:,fm:to);
        s3new=s3(2*fm-1:2*to+1);
        tmp=setupA3D(h,ksqnew,s1,s2,s3new);
        PB{b}=db_setup(N1,to-fm+1,chunksize,tmp);
        fprintf('complete: PB fm=%d to=%d\n',fm,to);
    end
    
    P={PM,PF,PB,pttn,size(ksq),h};
end

