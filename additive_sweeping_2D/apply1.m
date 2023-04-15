function[res]= apply1(P,inc)
PM=P{1};PF=P{2};PB=P{3};pttn=P{4};sz=P{5};NPAD=P{6};h=P{7};
N1=sz(1);N2=sz(2);
inc=reshape(inc,sz);
NP=numel(pttn);
res=zeros(sz);

    AUX = cell(1,NP);
    %self
    if(1)
        b=1;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L = PM{b}{1};        U = PM{b}{2};
        tmp = zeros(cnt+(NPAD-1),N2);
        tmp(1:cnt,:) = inc(fm:to,:);
        bak = U\(L\tmp(:));        bak = reshape(bak,size(tmp));
        AUX{b} = bak;
        res(fm:to,:) = AUX{b}(1:cnt,:);
    end
    for b=2:NP-1
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L = PM{b}{1};        U = PM{b}{2};
        tmp = zeros(cnt+2*(NPAD-1),N2);
        tmp((NPAD-1)+(1:cnt),:) = inc(fm:to,:);
        bak = U\(L\tmp(:));        bak = reshape(bak,size(tmp));
        AUX{b} = bak;
        res(fm:to,:) = AUX{b}((NPAD-1)+(1:cnt),:);
    end
    if(1)
        b = NP;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L = PM{b}{1};        U = PM{b}{2};
        tmp = zeros(cnt+(NPAD-1),N2);
        tmp((NPAD-1)+(1:cnt),:) = inc(fm:to,:);
        bak = U\(L\tmp(:));        bak = reshape(bak,size(tmp));
        AUX{b} = bak;
        res(fm:to,:) = AUX{b}((NPAD-1)+(1:cnt),:);
    end
    
    %forward sweeping
    if(1)
        b = 1;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        cur = AUX{b}(cnt,:);
    end
    for b=2:NP-1
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L=PF{b}{1};U=PF{b}{2};
        rhs=[(-1/(h*h))*cur;zeros(cnt+(NPAD-1)-1,N2)];
        tmp = U\(L\rhs(:));        tmp = reshape(tmp,cnt+(NPAD-1),N2);
        res(fm:to,:) = res(fm:to,:) + tmp(1:cnt,:);
        %update new
        cur = tmp(cnt,:) + AUX{b}((NPAD-1)+cnt,:);
    end
    if(1)
        b = NP;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L=PF{b}{1};U=PF{b}{2};
        rhs=[(-1/(h*h))*cur;zeros(cnt-1,N2)];
        tmp = U\(L\rhs(:));tmp = reshape(tmp,cnt,N2);
        res(fm:to,:) = res(fm:to,:) + tmp;
    end
    
    %backward sweeping
    if(1)
        b = NP;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        cur = AUX{b}((NPAD-1)+1,:); 
    end
    for b=NP-1:-1:2
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L=PB{b}{1};U=PB{b}{2};
        rhs=[zeros(cnt+(NPAD-1)-1,N2);(-1/(h*h))*cur];
        tmp = U\(L\rhs(:));        tmp = reshape(tmp,cnt+(NPAD-1),N2);
        res(fm:to,:) = res(fm:to,:)+tmp((NPAD-1)+(1:cnt),:);
        %update new
        cur = tmp((NPAD-1)+1,:) + AUX{b}((NPAD-1)+1,:);
    end
    if(1)
        b = 1;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        L=PB{b}{1};U=PB{b}{2};
        rhs=[zeros(cnt-1,N2);(-1/(h*h))*cur];
        tmp = U\(L\rhs(:));tmp = reshape(tmp,cnt,N2);
        res(fm:to,:) = res(fm:to,:) + tmp;
    end
    res = res(:);
end

