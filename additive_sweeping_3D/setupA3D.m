function[A]=setupA3D(h,ksq,s1,s2,s3)
    [N1,N2,N3]=size(ksq);P1=N1+2;P2=N2+2;P3=N3+2;
    
    %form matrix
    idx = zeros(P1,P2,P3);
    idx(2:P1-1,2:P2-1,2:P3-1) = reshape(1:N1*N2*N3,N1,N2,N3);
    
    MD1=2:P1-1;LF1=1:P1-2;RT1=3:P1;
    MD2=2:P2-1;LF2=1:P2-2;RT2=3:P2;
    MD3=2:P3-1;LF3=1:P3-2;RT3=3:P3;
    
    if(length(s1)~=N1*2+1||length(s2)~=N2*2+1||length(s3)~=N3*2+1)
        fprintf('error, s1 length %d, s2 length %d, s3 length %d, ksq size %d * %d * %d\n',length(s1),length(s2),length(s3),N1,N2,N3);
    end
    
    s1=reshape(s1,[length(s1),1,1]);
    s2=reshape(s2,[1,length(s2),1]);
    s3=reshape(s3,[1,1,length(s3)]);
    
    Il1 = idx(MD1,MD2,MD3);
    Jl1 = idx(LF1,MD2,MD3);
    Sl1 = repmat(1/(h*h)*(s1(2:2:2*N1).*s1(1:2:2*N1-1)),[1,N2,N3]);
    
    Ir1 = idx(MD1,MD2,MD3);
    Jr1 = idx(RT1,MD2,MD3);
    Sr1 = repmat(1/(h*h)*(s1(2:2:2*N1).*s1(3:2:2*N1+1)),[1,N2,N3]);
    
    Il2 = idx(MD1,MD2,MD3);
    Jl2 = idx(MD1,LF2,MD3);
    Sl2 = repmat(1/(h*h)*(s2(2:2:2*N2).*s2(1:2:2*N2-1)),[N1,1,N3]);
    
    Ir2 = idx(MD1,MD2,MD3);
    Jr2 = idx(MD1,RT2,MD3);
    Sr2 = repmat(1/(h*h)*(s2(2:2:2*N2).*s2(3:2:2*N2+1)),[N1,1,N3]);
        
    Il3 = idx(MD1,MD2,MD3);
    Jl3 = idx(MD1,MD2,LF3);
    Sl3 = repmat(1/(h*h)*(s3(2:2:2*N3).*s3(1:2:2*N3-1)),[N1,N2,1]);
    
    Ir3 = idx(MD1,MD2,MD3);
    Jr3 = idx(MD1,MD2,RT3);
    Sr3 = repmat(1/(h*h)*(s3(2:2:2*N3).*s3(3:2:2*N3+1)),[N1,N2,1]);
    
    Im = idx(MD1,MD2,MD3);
    Jm = idx(MD1,MD2,MD3);
    Sm = -(Sl1+Sl2+Sl3+Sr1+Sr2+Sr3) + ksq;
    
    Is = [Il1(:);Ir1(:);Il2(:);Ir2(:);Il3(:);Ir3(:);Im(:)];
    Js = [Jl1(:);Jr1(:);Jl2(:);Jr2(:);Jl3(:);Jr3(:);Jm(:)];
    Ss = [Sl1(:);Sr1(:);Sl2(:);Sr2(:);Sl3(:);Sr3(:);Sm(:)];
    
    gd = find(Js>0);
    Is = Is(gd); %keep the good and shift
    Js = Js(gd);
    Ss = Ss(gd);
    A = sparse(Is,Js,Ss);
end
