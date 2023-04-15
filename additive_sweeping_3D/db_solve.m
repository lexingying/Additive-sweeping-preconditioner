function u = db_solve(N,NLS,P,Tree,f)
  
  NG = N+1;
  L = round(log(NG/P)/log(2)) + 1;
  
  u = f;
  for ell=1:L
    W = 2^(ell-1)*P;
    nck = NG/W;
    %fprintf(1,'%d\n',ell);
    %tic;
    for i=1:nck
      for j=1:nck
        TT = Tree{ell}{i,j};
        GINS = TT.GINS;        GBDS = TT.GBDS;
        invA = TT.invA;        B = TT.B;
        %load
        in = u(GINS);
        bd = u(GBDS);
        %solve
        in = invA*in;
        bd = bd - B*in;
        %save
        u(GINS) = in;
        u(GBDS) = bd;
      end
    end
    %toc;
  end
  %-----------------
  for ell=L:-1:1
    W = 2^(ell-1)*P;
    nck = NG/W;
    %fprintf(1,'%d\n',ell);
    %tic;
    for i=1:nck
      for j=1:nck
        TT = Tree{ell}{i,j};
        GINS = TT.GINS;        GBDS = TT.GBDS;
        C = TT.C;
        %load
        in = u(GINS);
        bd = u(GBDS);
        %solve
        in = in - C*bd;
        %save
        u(GINS) = in;
        u(GBDS) = bd;
      end
    end
    %toc;
  end
  %-----------------
  
  