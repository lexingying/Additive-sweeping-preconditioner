function Tree = db_setup(N,NLS,P,MAT)
  
  NG = N+1;
  L = round(log(NG/P)/log(2)) + 1;
  
  [ii,jj,kk] = ndgrid(1:N,1:N,1:NLS);
  T = zeros(N,N,NLS);
  edge = find(mod(ii,P)==0 | mod(jj,P)==0);
  T(edge) = 1;
  vert = find(mod(ii,P)==0 & mod(jj,P)==0);
  T(vert) = 2;
  
  [I,J,S] = find(MAT);
  tmp = find(T(I)>=1 & T(J)>=1 & T(I)+T(J)<4);  S(tmp) = S(tmp)/2;
  tmp = find(T(I)>=2 & T(J)>=2);  S(tmp) = S(tmp)/4;
  NEW = sparse(I,J,S); %matrix weighted corrected
  
  GRD = reshape(1:N*N*NLS,[N,N,NLS]);
  Tree = cell(L,1);
  if(1)
    ell = 1;
    W = 2^(ell-1)*P;
    nck = NG/W;
    Tree{ell} = cell(nck,nck);
    for i=1:nck
      for j=1:nck
        ia = (i-1)*W;        ib = i*W;        im = (ia+ib)/2;        is = [max(ia,1):min(ib,N)];
        ja = (j-1)*W;        jb = j*W;        jm = (ja+jb)/2;        js = [max(ja,1):min(jb,N)];
        ks = 1:NLS;
        glbidx = GRD(is,js,ks);
        [ii,jj,kk] = ndgrid(is,js,ks);
        gud = ones(size(ii));
        lclidx = zeros(size(ii));
        lclidx(find(gud==1)) = 1:numel(find(gud==1));
        isin = (ii~=ia & ii~=ib & jj~=ja & jj~=jb);
        LINS = lclidx(find(gud==1 & isin==1));
        GINS = glbidx(find(gud==1 & isin==1));
        LBDS = lclidx(find(gud==1 & isin==0));
        GBDS = glbidx(find(gud==1 & isin==0));
        ALL = full(NEW(glbidx,glbidx));
        invA = inv(ALL(LINS,LINS));
        B = ALL(LBDS,LINS);
        C=invA*ALL(LINS,LBDS);
        S = ALL(LBDS,LBDS) - B*C;
        Tree{ell}{i,j} = struct('LINS', LINS, 'LBDS', LBDS, 'GINS', GINS, 'GBDS', GBDS, 'invA', invA, 'B', B,'C',C, 'S', S); clear invA B C S;
      end
    end
  end
  %---
  for ell=2:L
    W = 2^(ell-1)*P;    Q = W/2;
    nck = NG/W;
    Tree{ell} = cell(nck,nck);
    for i=1:nck
      for j=1:nck
        ia = (i-1)*W;        ib = i*W;        im = (ia+ib)/2;        is = [max(ia,1):min(ib,N)];
        ja = (j-1)*W;        jb = j*W;        jm = (ja+jb)/2;        js = [max(ja,1):min(jb,N)];
        ks = 1:NLS;
        glbidx = GRD(is,js,ks);
        [ii,jj,kk] = ndgrid(is,js,ks);
        gud = (ii==ia | ii==im | ii==ib | jj==ja | jj==jm | jj==jb);
        lclidx = zeros(size(ii));
        lclidx(find(gud==1)) = 1:numel(find(gud==1));
        isin = (ii~=ia & ii~=ib & jj~=ja & jj~=jb);
        LINS = lclidx(find(gud==1 & isin==1));
        GINS = glbidx(find(gud==1 & isin==1));
        LBDS = lclidx(find(gud==1 & isin==0));
        GBDS = glbidx(find(gud==1 & isin==0));
        %construct mat
        nb = numel(find(gud==1));
        ALL = zeros(nb,nb);
        cut = ( (ia<=ii&ii<=im) & (ja<=jj&jj<=jm) );        idx = lclidx(find(gud==1 & cut==1));
        ALL(idx,idx) = ALL(idx,idx) + Tree{ell-1}{2*i-1,2*j-1}.S; Tree{ell-1}{2*i-1,2*j-1}.S=[];
        cut = ( (im<=ii|ii==ib) & (ja<=jj&jj<=jm) );        idx = lclidx(find(gud==1 & cut==1));
        ALL(idx,idx) = ALL(idx,idx) + Tree{ell-1}{2*i  ,2*j-1}.S; Tree{ell-1}{2*i  ,2*j-1}.S=[];
        cut = ( (ia<=ii&ii<=im) & (jm<=jj|jj==jb) );        idx = lclidx(find(gud==1 & cut==1));
        ALL(idx,idx) = ALL(idx,idx) + Tree{ell-1}{2*i-1,2*j  }.S; Tree{ell-1}{2*i-1,2*j  }.S=[];
        cut = ( (im<=ii|ii==ib) & (jm<=jj|jj==jb) );        idx = lclidx(find(gud==1 & cut==1));
        ALL(idx,idx) = ALL(idx,idx) + Tree{ell-1}{2*i  ,2*j  }.S; Tree{ell-1}{2*i  ,2*j  }.S=[];
        %Schur complement
        invA = inv(ALL(LINS,LINS));
        B = ALL(LBDS,LINS);
        C=invA*ALL(LINS,LBDS);
        S = ALL(LBDS,LBDS) - B*C;
        Tree{ell}{i,j} = struct('LINS', LINS, 'LBDS', LBDS, 'GINS', GINS, 'GBDS', GBDS, 'invA', invA, 'B', B,'C',C, 'S', S); clear invA B C S;
      end
    end
  end
  