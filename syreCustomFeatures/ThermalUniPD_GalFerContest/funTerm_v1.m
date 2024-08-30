function [Twind_mean,Twind_min,Twind_max] = funTerm_v1(geo,mat,opt,T,P,M)
%%
n_repetition_stack=opt.n_repetition_stack;
n_repetition_head=opt.n_repetition_head;
lstack=(geo.l*1e-3)/2;
lhead=geo.l*geo.lfact*1e-3;
r_min_head=opt.r_min_head;
r_max_head=opt.r_max_head;
Nt=opt.Nt;
Time=opt.Time;
Tamb=opt.Tamb;
times=linspace(0,Time,Nt);
q=opt.q;
%%
deltaT=times(2)-times(1);
%%
if opt.plot
    figure
    hold on
    patch('Faces',T,'Vertices',P,'CData',M,'FaceColor','flat','FaceAlpha',1)
    axis equal
    colorbar
    drawnow
end
%%
barTri=zeros(size(T,1),3);
for ii=1:size(T,1)
    barTri(ii,1:3)=sum(P(T(ii,:),:))/3;
end
[PN,PM,PL,PB] =funTri2Prism_2l(T,M,P(:,[1,3,2]),[n_repetition_head,n_repetition_stack],[lhead,lstack],barTri(:,[1,3,2]));
% delete elements 
rr=sqrt(PB(:,1).^2+PB(:,3).^2);
% todelete=find(PB(:,2)<lhead & (rr>r_max_head | rr<r_min_head));
todelete=find(PB(:,2)<lhead &  rr<r_min_head);
PM(todelete,:)=[];
PL(todelete,:)=[];
PB(todelete,:)=[];
rr(todelete)=[];
% assign new flag for head windings
headwind=(PB(:,2)<lhead & rr<r_max_head);
exthead=(PB(:,2)<lhead & rr>r_max_head);
PL(headwind)=max(PL)+1; 
PL(exthead)=max(PL)+1; 
% delete nodes
[keep_nodes]=intersect(1:size(PN,1),unique(PM(:)));
map=zeros(size(PN,1),1);
map(keep_nodes)=1:length(keep_nodes);
PN=PN(keep_nodes,:);
PM(:)=map(PM(:));
%
maxrr=max(rr);
%
[TTM,TTL,TTN]=funPrism2Tetra(PN,PM,PL);
nTet=size(TTM,1);
%% assign material
labels=unique(TTL);
l_V=zeros(nTet,1);
c_V=zeros(nTet,1);
r_V=zeros(nTet,1);
q_V=zeros(nTet,1);
for kk = 1:length(labels)
    ii=labels(kk);
    ww=(TTL==ii); % tetra elements with tag ii
    if ii<=length(mat.CondTerm)
        if mat.CondTerm(ii)==0 || mat.kgm3(ii)==0 || mat.HeatCap(ii)==0
            error(['missing material properties ii=',num2str(ii)])
        end
    else
        error(['missing material properties ii=',num2str(ii)])
    end
    l_V(ww)=mat.CondTerm(ii); % conducibilità termica
    r_V(ww)=mat.kgm3(ii);     % densità
    c_V(ww)=mat.HeatCap(ii);  % calore specifico

%     if ii==25 % LayerAir 
%         l_V(ww)=mat.LayerAir.CondTerm; % conducibilità termica
%         r_V(ww)=mat.LayerAir.kgm3;     % densità
%         c_V(ww)=mat.LayerAir.HeatCap;  % calore specifico
% %         h_w(ww)=0;                     % convective coefficient
% %         q_V(ww)=0;                     % power 
%     elseif ii == 27 || ii==1 % extra iron head || iron
%         l_V(ww)=mat.Stator.CondTerm; % conducibilità termica
%         r_V(ww)=mat.Stator.kgm3;     % densità
%         c_V(ww)=mat.Stator.HeatCap;  % calore specifico
%     elseif ii == 26 || ii==3 % winding head || winding
%         l_V(ww)=mat.SlotCond.CondTerm; % conducibilità termica
%         r_V(ww)=mat.SlotCond.kgm3;     % densità
%         c_V(ww)=mat.SlotCond.HeatCap;  % calore specifico
%     elseif ii == 9 || ii==8 % magnet
%         l_V(ww)=mat.LayerMag.CondTerm; % conducibilità termica
%         r_V(ww)=mat.LayerMag.kgm3;     % densità
%         c_V(ww)=mat.LayerMag.HeatCap;  % calore specifico
%     elseif ii == 7 % shaft and insulation
%         l_V(ww)=mat.Shaft.CondTerm; % conducibilità termica
%         r_V(ww)=mat.Shaft.kgm3;     % densità
%         c_V(ww)=mat.Shaft.HeatCap;  % calore specifico        
%     else
%        error('missing material')
% %         h_V(ww)=0;
% %         q_V(ww)=0;
%     end
end
% q_V((TTL==3 | TTL==26))=1;

if isempty(mat.q)
    error('mat.q is empty')
else
    for ii = 1:length(mat.q)
        loc=TTL==mat.q(ii);
        if isempty(loc)
             error(['wrong assignement in mat.q, mat.q(',num2str(ii),') exceed max flag number'])
        else
            q_V(loc)=1;
        end
    end
end
h=opt.h;
%%
% tic
[~,~,~,F1] = gcd_mexed(TTM.');
% toc
% nv=size(D1,2);
% nf=size(c,2);
% D volumes x faces
% Dpos = (D1+abs(D1))/2;
% Dneg = (D1-abs(D1))/2;
% [r_pos_d,c_pos_d,val_pos_d] = find(Dpos);
% [r_neg_d,c_neg_d,val_neg_d] = find(Dneg);
% Matrix_D = sparse(c_pos_d,val_pos_d,ones(size(c_pos_d,1),1),nv,nf);
% Matrix_D = Matrix_D+sparse(c_neg_d,abs(val_neg_d),-ones(size(c_neg_d,1),1),nv,nf);
%  Dnew = full(Dnew);
ind_face_free=find(F1(5,:)==0);
% [vffree,J]=find(Matrix_D(:,ind_face_free)); % volumi attaccati alle facce libere
%% find faces connected to fluid
rrn=sqrt(TTN(:,1).^2+TTN(:,3).^2);
rrf=rrn(F1(1:3,ind_face_free));
rrfl=find(prod(rrf>maxrr));
ind_face_fl=ind_face_free(rrfl);
%%
if opt.plot
    ccol=zeros(length(ind_face_free),1);
    ccol(rrfl)=1;
    figure
    hold on
    patch('Faces',F1(1:3,ind_face_free).','Vertices',TTN,'CData',ccol,'FaceColor','flat','FaceAlpha',1)
    axis equal
    view(3)
    colorbar
    drawnow
    figure
    patch('Faces',TTM(:,1:3),'Vertices',TTN,'CData',TTL,'FaceColor','flat','FaceAlpha',1)
    hold on
    patch('Faces',TTM(:,[1,2,4]),'Vertices',TTN,'CData',TTL,'FaceColor','flat','FaceAlpha',1)
    patch('Faces',TTM(:,[2,3,4]),'Vertices',TTN,'CData',TTL,'FaceColor','flat','FaceAlpha',1)
    patch('Faces',TTM(:,[1,3,4]),'Vertices',TTN,'CData',TTL,'FaceColor','flat','FaceAlpha',1)
    view(3)
    axis equal
    colorbar
    drawnow
    %
    for kk = 1:length(labels)
        ii=labels(kk);
        ww=TTL==ii;
        figure
        patch('Faces',TTM(:,1:3),'Vertices',TTN,'CData',TTL(:),'FaceColor','flat','FaceAlpha',0)
        hold on
        patch('Faces',TTM(:,[1,2,4]),'Vertices',TTN,'CData',TTL(:),'FaceColor','flat','FaceAlpha',0)
        patch('Faces',TTM(:,[2,3,4]),'Vertices',TTN,'CData',TTL(:),'FaceColor','flat','FaceAlpha',0)
        patch('Faces',TTM(:,[1,3,4]),'Vertices',TTN,'CData',TTL(:),'FaceColor','flat','FaceAlpha',0)
        %
        patch('Faces',TTM(ww,1:3),'Vertices',TTN,'CData',TTL(ww),'FaceColor','flat','FaceAlpha',1)
        hold on
        patch('Faces',TTM(ww,[1,2,4]),'Vertices',TTN,'CData',TTL(ww),'FaceColor','flat','FaceAlpha',1)
        patch('Faces',TTM(ww,[2,3,4]),'Vertices',TTN,'CData',TTL(ww),'FaceColor','flat','FaceAlpha',1)
        patch('Faces',TTM(ww,[1,3,4]),'Vertices',TTN,'CData',TTL(ww),'FaceColor','flat','FaceAlpha',1)
        view(3)
        axis equal
        colorbar
        title(ii)
        drawnow
    end
end
%%
% tic
[A,A_,rhs,N2T,Vols] = funFEM_loc(F1(:,ind_face_fl),TTN,TTM,c_V.*r_V,l_V,h,deltaT,q_V,q,Tamb);
% toc
% [A,A_,rhs,N2T,Vols] = funFEM_loc(F1(:,ind_face_fl),TTN,TTM,0*c_V.*r_V,l_V,h,1,q_V,q,Tamb);
% xx=A\rhs
% return
%%
% tic
x=rhs*0+Tamb;
for ii = 1:Nt
    x=A\(A_*x+rhs);
end
% toc
%% RICCC
y=N2T*x;
ww=find(sum(TTL==mat.Tout,2)); % windings
Twind_mean=sum(Vols(ww).*y(ww))/sum(Vols(ww));
Twind_min=min(y(ww));%sum(Vols(ww).*y(ww))/sum(Vols(ww));
Twind_max=max(y(ww));
%%
if opt.plot_final
    figure
    patch('Faces',TTM(:,1:3),'Vertices',TTN,'CData',x,'Facecolor','interp','FaceAlpha',1.0)
    patch('Faces',TTM(:,[1,2,4]),'Vertices',TTN,'CData',x,'Facecolor','interp','FaceAlpha',1.0)
    patch('Faces',TTM(:,[2,3,4]),'Vertices',TTN,'CData',x,'Facecolor','interp','FaceAlpha',1.0)
    patch('Faces',TTM(:,[3,1,4]),'Vertices',TTN,'CData',x,'Facecolor','interp','FaceAlpha',1.0)
    axis equal
    view(3)
    colorbar
    colormap hot
    drawnow
end
end
%%
function [A,A_,rhs,N2T,Vols] = funFEM_loc(F,P,VP,rc,l,h,deltaT,q_V,q,Tamb)
rhs=zeros(max(VP(:)),1);
Nn = size(P,1);
Ne = size(VP,1);
inds = 0;
A_Si = zeros(16*Ne,1); 
A_Sj = zeros(16*Ne,1); 
A_Sk = zeros(16*Ne,1);
A_Si_ = zeros(16*Ne,1); 
A_Sj_ = zeros(16*Ne,1); 
A_Sk_ = zeros(16*Ne,1); 
N2T_ = zeros(4,Ne); 
Vols=zeros(Ne,1);
ii = 1:4; ji = 1:4; [ii, ji] = ndgrid(ii,ji);
ii = ii(:); ji = ji(:);
%%
wg=[0.0091694;0.016027;0.021157;0.03698;0.0091694;0.016027;0.021157;0.03698];
p1=[0.54415;0.54415;0.12251;0.12251;0.54415;0.54415;0.12251;0.12251];
p2=[0.294;0.07068;0.56593;0.13605;0.294;0.07068;0.56593;0.13605];
p3=[0.034203;0.081396;0.065839;0.15668;0.12765;0.30377;0.24571;0.58475];
p4=[0.12765;0.30377;0.24571;0.58475;0.034203;0.081396;0.065839;0.15668];
%%
for jj = 1:Ne
    eleNodes = VP(jj,1:4);
    pp=P(eleNodes,:);
    PG = p1*pp(1,:)+p2*pp(2,:)+p3*pp(3,:)+p4*pp(4,:);
    [Vol] = fun_my_vol_tetra_loc(pp);
    N2T_(:,jj)=eleNodes;
    Vols(jj)=Vol;
    WG = wg*Vol*6;
    B = [[1, 1, 1, 1];pp'];
    invB=inv(B);
    Kl2_K = invB(:,2:4)*invB(:,2:4)'.*abs(Vol);
    Kl2 = zeros(16,1);
    % 
%     if abs(rc(jj))>0
    for kk = 1:8 
        xsi=invB*[1;PG(kk,1);PG(kk,2);PG(kk,3)]; 
        ss=1;
        for ee=1:4
            rhs(eleNodes(ee),1)=rhs(eleNodes(ee),1)+xsi(ee)*q_V(jj)*WG(kk)*deltaT;
            for qq = 1:4
                Kl2(ss,1)=Kl2(ss,1)+xsi(ee)*xsi(qq)*WG(kk);
                ss=ss+1;
            end
        end
    end
%     end
    inds = inds(end) + (1:16);
    A_Si(inds) = eleNodes(ii);
    A_Sj(inds) = eleNodes(ji);
    A_Sk(inds) = Kl2(:)*rc(jj)+Kl2_K(:)*l(jj)*deltaT;
    A_Si_(inds) = eleNodes(ii);
    A_Sj_(inds) = eleNodes(ji);
    A_Sk_(inds) = Kl2(:)*rc(jj);    
end
%%
denq=q/sum(Vols(q_V==1));
rhs=rhs*denq;
%%
for jj=1:size(F,2) 
   p=F(1:4,jj);
   ptri=P(p,:);
   [PG,WG,~]=GaussGraglia7_loc(ptri.');
   PG=PG.';
   %
%    v=vffree(jj);
   eleNodes = F(1:4,jj);
   ptet=P(eleNodes,:);
   B = [[1, 1, 1, 1];ptet'];
   Kl2 = zeros(16,1);
%    if abs(h(v))>0
   for hh = 1:7
       xsi=B\[1;PG(hh,1);PG(hh,2);PG(hh,3)]; 
        ss=1;
        for ee=1:4
            rhs(p(ee),1)=rhs(p(ee),1)+xsi(ee)*WG(hh)*h*deltaT*Tamb;
            for qq = 1:4
                Kl2(ss,1)=Kl2(ss,1)+xsi(ee)*xsi(qq)*WG(hh);
                ss=ss+1;
            end
        end
   end
%    end
    inds = inds(end) + (1:16);
    A_Si(inds) = eleNodes(ii);
    A_Sj(inds) = eleNodes(ji);
    A_Sk(inds) = Kl2(:)*h*deltaT;
end
%%
A = sparse(A_Si,A_Sj,A_Sk,Nn,Nn);
A_ = sparse(A_Si_,A_Sj_,A_Sk_,Nn,Nn);
N2T=sparse(1:Ne,N2T_(1,:),0.25*ones(Ne,1),Ne,Nn)+...
    sparse(1:Ne,N2T_(2,:),0.25*ones(Ne,1),Ne,Nn)+...
    sparse(1:Ne,N2T_(3,:),0.25*ones(Ne,1),Ne,Nn)+...
    sparse(1:Ne,N2T_(4,:),0.25*ones(Ne,1),Ne,Nn);
end
%%
function [PP,w,area]=GaussGraglia7_loc(NN)
r1=NN(:,1);
r2=NN(:,2);
r3=NN(:,3);
e1=r3-r2;
e2=r1-r3;
area=norm(cross(e1,e2))/2;
%points    
v1=1/3;
v2=(6-sqrt(15))/21;
v3=(9+2*sqrt(15))/21;
v4=(6+sqrt(15))/21;
v5=(9-2*sqrt(15))/21;
eta=[v1;v2;v2;v3;v4;v4;v5];
csi=[v1;v2;v3;v2;v4;v5;v4];
%weights
wa=9/40*area;
wb=(155-sqrt(15))/1200*area;
wc=(155+sqrt(15))/1200*area;
w=[wa,wb,wb,wb,wc,wc,wc];
N1=1-eta-csi;
N2=eta;
N3=csi;
PP=([N1,N2,N3]*[r1';r2';r3'])';
end
%
function [V] = fun_my_vol_tetra_loc(tetra)
v1 =  tetra(2,:) - tetra(1,:);
v2 =  tetra(3,:) - tetra(1,:);
v3 =  tetra(4,:) - tetra(1,:);
V = det([v1;v2;v3])/6;
end
%%
function [PN,PM,PL,PB] = funTri2Prism_2l(TM,TL,TN,n_repetitions,lung,barTri)
%% To prism 
n_repetition=sum(n_repetitions);
ll1=linspace(0,lung(1),n_repetitions(1)+1);
ll2=linspace(lung(1),lung(2),n_repetitions(2)+1);
ll=[ll1,ll2(2:end)];
n_points_tri=size(TN,1);
n_triangles=size(TM,1);
n_points_pri=n_points_tri*(n_repetition+1);
n_prisms=n_triangles*n_repetition;
PN=zeros(n_points_pri,3);
PM=zeros(n_prisms,6);
PL=zeros(n_prisms,1);
PB=zeros(n_prisms,3);
for ii = 1:n_repetition
    PN([1:n_points_tri]+(ii-1)*n_points_tri,:)=...
        TN+repmat([0 ll(ii) 0],n_points_tri,1);
    %
    PM([1:n_triangles]+(ii-1)*n_triangles,:)=...
    [TM+n_points_tri*(ii-1),TM+n_points_tri*(ii)];
    %
    PL(1+(ii-1)*n_triangles:(ii)*n_triangles,1)=TL;
    PB(1+(ii-1)*n_triangles:(ii)*n_triangles,:)=barTri+repmat([0 0.5*(ll(ii)+ll(ii+1)) 0],n_triangles,1);
end
ii=n_repetition+1;
PN([1:n_points_tri]+(ii-1)*n_points_tri,:)=...
        TN+repmat([0 ll(ii) 0],n_points_tri,1);
end
%%
function [TM,TL,TN] = funPrism2Tetra(PN,PM,PL)
%%
n_prisms=size(PM,1);
n_points_pri=size(PN,1);
TN=[PN;zeros(n_prisms,3)]; % tetra nodes
n_points_tet=n_points_pri+n_prisms;
n_tets=n_prisms*8;
TM=zeros(n_tets,4);
TL=zeros(n_tets,1);
for ii = 1:n_prisms
    node=sum(PN(PM(ii,:),1:3))/6; % central node
    TN(n_points_pri+ii,:)=node; % add central node to points
    tet1=[PM(ii,[3,2,1]),n_points_pri+ii]; % tet 1 
    tet2=[PM(ii,[4,5,6]),n_points_pri+ii]; % tet 2
    edg1=sort(PM(ii,[1,2]));
    edg2=sort(PM(ii,[2,3]));
    edg3=sort(PM(ii,[3,1]));
    edg1p=sort(PM(ii,3+[1,2]));
    edg2p=sort(PM(ii,3+[2,3]));
    edg3p=sort(PM(ii,3+[3,1]));    
    tet3=[edg1,edg1p(1),n_points_pri+ii];
    tet4=[edg1p([2,1]),edg1(2),n_points_pri+ii];
    tet5=[edg2,edg2p(1),n_points_pri+ii];
    tet6=[edg2p([2,1]),edg2(2),n_points_pri+ii];
    tet7=[edg3([2,1]),edg3p(1),n_points_pri+ii];
    tet8=[edg3p,edg3(2),n_points_pri+ii];
    TM(1+(ii-1)*8,:)=tet1;
    TM(2+(ii-1)*8,:)=tet2;
    TM(3+(ii-1)*8,:)=tet3;
    TM(4+(ii-1)*8,:)=tet4;
    TM(5+(ii-1)*8,:)=tet5;
    TM(6+(ii-1)*8,:)=tet6;
    TM(7+(ii-1)*8,:)=tet7;
    TM(8+(ii-1)*8,:)=tet8;
    TL([1:8]+(ii-1)*8)=PL(ii);
end
for ii = 1:n_tets
    tetra=TN(TM(ii,:),:);
    [V] = fun_vol_tetra(tetra);
    if V<0
        TM(ii,:)=TM(ii,[3,2,1,4]);
    end
end
end
function [V] = fun_vol_tetra(tetra)
    v1 =  tetra(2,1:3) - tetra(1,1:3);
    v2 =  tetra(3,1:3) - tetra(1,1:3);
    v3 =  tetra(4,1:3) - tetra(1,1:3);
    V = (v1(1)*(v2(2)*v3(3)-v2(3)*v3(2))-v1(2)*(v2(1)*v3(3)...
        -v2(3)*v3(1))+v1(3)*(v2(1)*v3(2)-v2(2)*v3(1)))/6;
end