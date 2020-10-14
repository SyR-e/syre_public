
hsEq
wsEq
teq
Acu=Aslot*Kcu;
w=(wsEq-2*teq)/3;
h2=hsEq-teq/2-w;
h1=hsEq-teq-w;
Acu1=w*h1;
Ps11=2*hsEq-w;
Ains1=Ps11*teq/2;
Acu2=(2*hsEq+wsEq-2*teq)*w;
Acu2=2*w*(hsEq-teq/2)+(w+teq)*(hsEq-h2-teq/2);

xv=[0,0;
    0,wsEq;
    0+teq/2,0+teq/2;
    0+teq/2,wsEq-teq/2;
    teq/2+w,teq/2+w;
    teq+w,teq+w;
    teq+2*w,teq+2*w;
    wsEq-w-teq/2,wsEq-w-teq/2];
yv=[0,hsEq;
    hsEq,hsEq;
    0,hsEq-teq/2;
    hsEq-teq/2,hsEq-teq/2;
    0,h2;
    0,h1;
    0,h1;
    0,h2];

figure(1);plot(xv',yv','b'); axis  equal