function P = setup1(NPML,NLPD,NPAD,pL,pR,h,ksq,s1,s2)
[N1,N2]=size(ksq);

%generate the partitions
NP=ceil((N1-2*(NPML-1))/NLPD);
pttn=cell(NP,1);
pttn{1}=[1,(NPML-1)+NLPD];
for g=2:NP-1
    pttn{g}=(NPML-1)+(g-1)*NLPD+[1,NLPD];
end
pttn{NP}=[(NPML-1)+(NP-1)*NLPD+1,N1];
    
    %construct matrices
    PM=cell(NP,1);
    if(1)
        b=1;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm:to+(NPAD-1),:);
        s1new=[s1(2*fm-1:2*to),pR];
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PM{b}={L,U};
    end
    for b=2:NP-1
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm-(NPAD-1):to+(NPAD-1),:);
        s1new=[pL,s1(2*fm:2*to),pR];
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PM{b}={L,U};
    end
    if(1)
        b=NP;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm-(NPAD-1):to,:);
        s1new=[pL,s1(2*fm:2*to+1)];
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PM{b}={L,U};
    end
    
    %construct forward sweeping matrices
    PF=cell(NP,1);
    for b=2:NP-1
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm:to+(NPAD-1),:);
        s1new=[s1(2*fm-1:2*to),pR];
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PF{b}={L,U};
    end
    if(1)
        b=NP;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm:to,:);
        s1new=s1(2*fm-1:2*to+1);
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PF{b}={L,U};
    end

    %construct backward sweeping matrices
    PB=cell(NP,1);
    for b=NP-1:-1:2
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm-(NPAD-1):to,:);
        s1new=[pL,s1(2*fm:2*to+1)];
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PB{b}={L,U};
    end
    if(1)
        b=1;
        fm=pttn{b}(1);to=pttn{b}(2);
        ksqnew=ksq(fm:to,:);
        s1new=s1(2*fm-1:2*to+1);
        tmp=setupA2D(h,ksqnew,s1new,s2);
        [L,U]=lu(tmp,0);
        PB{b}={L,U};
    end
    
    P = {PM,PF,PB,pttn,size(ksq),NPAD,h};
end

