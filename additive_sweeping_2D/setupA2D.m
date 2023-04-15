function[A]=setupA2D(h,ksq,s1,s2)
    [N1,N2]=size(ksq);P1=N1+2;P2=N2+2;
    
    %form matrix
    idx = zeros(P1,P2);
    idx(2:P1-1,2:P2-1) = reshape(1:N1*N2,N1,N2);
    
    MD1 = 2:P1-1;LF1 = 1:P1-2;RT1 = 3:P1;
    MD2 = 2:P2-1;LF2 = 1:P2-2;RT2 = 3:P2;
    
    if(length(s1)~=N1*2+1||length(s2)~=N2*2+1)
        fprintf('error, s1 length %d, s2 length %d, c size %d * %d\n',length(s1),length(s2),N1,N2);
    end
    
    s1=reshape(s1,[length(s1),1]);
    s2=reshape(s2,[1,length(s2)]);
    
    Il1 = idx(MD1,MD2);
    Jl1 = idx(LF1,MD2);
    Sl1 = repmat(1/(h*h)*(s1(2:2:2*N1).*s1(1:2:2*N1-1)),[1,N2]);
    
    Ir1 = idx(MD1,MD2);
    Jr1 = idx(RT1,MD2);
    Sr1 = repmat(1/(h*h)*(s1(2:2:2*N1).*s1(3:2:2*N1+1)),[1,N2]);
    
    Il2 = idx(MD1,MD2);
    Jl2 = idx(MD1,LF2);
    Sl2 = repmat(1/(h*h)*(s2(2:2:2*N2).*s2(1:2:2*N2-1)),[N1,1]);
    
    Ir2 = idx(MD1,MD2);
    Jr2 = idx(MD1,RT2);
    Sr2 = repmat(1/(h*h)*(s2(2:2:2*N2).*s2(3:2:2*N2+1)),[N1,1]);
    
    Im = idx(MD1,MD2);
    Jm = idx(MD1,MD2);
    Sm = -(Sl1+Sr1+Sl2+Sr2) + ksq;
    
    Is = [Il1(:); Ir1(:); Ir2(:); Il2(:); Im(:)];
    Js = [Jl1(:); Jr1(:); Jr2(:); Jl2(:); Jm(:)];
    Ss = [Sl1(:); Sr1(:); Sr2(:); Sl2(:); Sm(:)];
    
    gd = find(Js>0);
    Is = Is(gd); %keep the good and shift
    Js = Js(gd);
    Ss = Ss(gd);
    A = sparse(Is,Js,Ss);
end
