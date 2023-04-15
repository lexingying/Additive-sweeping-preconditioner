function[res]=apply1(NPAD,P,leafsize,inc)
PM=P{1};PF=P{2};PB=P{3};
pttn=P{4};
sz=P{5};
h=P{6};
N1=sz(1);N2=sz(2);N3=sz(3);
inc = reshape(inc,sz);
NP = numel(pttn);
res = zeros(sz);    
    
    AUX = cell(NP,1);
    %self
    if(1)
        b=1;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        tmp = zeros(N1,N2,cnt+(NPAD-1));
        tmp(:,:,1:cnt) = inc(:,:,fm:to);
        bak=db_solve(N1,cnt+(NPAD-1),leafsize,PM{b},tmp);
        AUX{b} = bak;
        res(:,:,fm:to) = AUX{b}(:,:,1:cnt);
    end
    for b=2:NP-1
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        tmp = zeros(N1,N2,cnt+2*(NPAD-1));
        tmp(:,:,(NPAD-1)+(1:cnt)) = inc(:,:,fm:to);
        bak=db_solve(N1,cnt+2*(NPAD-1),leafsize,PM{b},tmp);
        AUX{b} = bak;
        res(:,:,fm:to) = AUX{b}(:,:,(NPAD-1)+(1:cnt));
    end
    if(1)
        b = NP;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        tmp = zeros(N1,N2,cnt+(NPAD-1));
        tmp(:,:,(NPAD-1)+(1:cnt)) = inc(:,:,fm:to);
        bak=db_solve(N1,cnt+(NPAD-1),leafsize,PM{b},tmp);
        AUX{b} = bak;
        res(:,:,fm:to) = AUX{b}(:,:,(NPAD-1)+(1:cnt));
    end
    
    %forward sweeping
    if(1)
        b = 1;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        cur = AUX{b}(:,:,cnt);
    end
    for b=2:NP-1
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        rhs=zeros(N1,N2,cnt+(NPAD-1));rhs(:,:,1)=(-1/(h*h))*cur;
        tmp =db_solve(N1,cnt+(NPAD-1),leafsize,PF{b},rhs);
        res(:,:,fm:to) = res(:,:,fm:to) + tmp(:,:,1:cnt);
        %update new
        cur = tmp(:,:,cnt) + AUX{b}(:,:,(NPAD-1)+cnt);
    end
    if(1)
        b = NP;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        rhs=zeros(N1,N2,cnt);rhs(:,:,1)=(-1/(h*h))*cur;
        tmp =db_solve(N1,cnt,leafsize,PF{b},rhs);
        res(:,:,fm:to) = res(:,:,fm:to) + tmp;
    end
    
    %backward sweeping
    if(1)
        b = NP;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        cur = AUX{b}(:,:,(NPAD-1)+1); 
    end
    for b=NP-1:-1:2
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        rhs=zeros(N1,N2,cnt+(NPAD-1));rhs(:,:,end)=(-1/(h*h))*cur;
        tmp =db_solve(N1,cnt+(NPAD-1),leafsize,PB{b},rhs);
        res(:,:,fm:to) = res(:,:,fm:to) + tmp(:,:,(NPAD-1)+(1:cnt));
        %update new
        cur = tmp(:,:,(NPAD-1)+1) + AUX{b}(:,:,(NPAD-1)+1);
    end
    if(1)
        b = 1;
        fm = pttn{b}(1);        to = pttn{b}(2);        cnt = to-fm+1;
        rhs=zeros(N1,N2,cnt);rhs(:,:,end)=(-1/(h*h))*cur;
        tmp =db_solve(N1,cnt,leafsize,PB{b},rhs);
        res(:,:,fm:to) = res(:,:,fm:to) + tmp;
    end
    res = res(:);
end

