
XpontRadBarSx = zeros(nlay);
YpontRadBarSx = zeros(nlay);
XpontRadBarDx = zeros(nlay);
YpontRadBarDx = zeros(nlay);
XpontRadDx    = zeros(nlay);
YpontRadDx    = zeros(nlay);
XpontRadSx    = zeros(nlay);
YpontRadSx    = zeros(nlay);

switch geo.RotType
    case 'Vtype' % rev. Gallo 07/05/2018
        %% Vtype
        % NOTA: Ho anticipato il calcolo spessore ponticello radiale in funzione della velocità di rotazione (Massa appesa ai ponticelli)
        % in "nodes_rotor_Vtype"
        %Nota: il calcolo è stato rivisto per tenere conto della massa dei magneti
        %per calcolare la massa appesa ai ponticelli - rev.Gallo
        
        %Inizializzazione delle variabili
        XcRaccpontRadSx=zeros(1,1);
        YcRaccpontRadSx=zeros(1,1);
        RcRaccpontRadSx=zeros(1,1);
        XcRaccpontRadDx=zeros(1,1);
        YcRaccpontRadDx=zeros(1,1);
        RcRaccpontRadDx=zeros(1,1);
        
        %%% DISEGNO PONTICELLI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jj=1:nlay
            m=tan(angle);
            %Definisco il raccordo ponticelli radiali mediante racc_pont
            %formule trigonometriche necessarie:
            cos_sottr=cos(pi/2)*cos(angle)+sin(pi/2)*sin(angle); %cos(alpha-beta)
            sin_sottr=sin(pi/2)*cos(angle)-cos(pi/2)*sin(angle); %sin(alpha-beta)
            
            hpont=pont(1)/2;
            %a1=-m;
            %b1=1;
            %c1=m*B1k;
            [a1,b1,c1]=retta_per_2pti(B1k,0,xxD1k,yyD1k);
            a2=0;
            b2=1;
            c2=-hpont;
            [x_temp1,y_temp1]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
            x10=x_temp1+racc_pont/cos_sottr; %devo definirlo cosi per avere condione di tangenza con retta lato barriera
            y10=y_temp1;
            [a1,b1,c1]=retta_per_2pti(B2k,0,xxD2k,yyD2k);
            %a1=-m;
            %b1=1;
            %c1=m*B2k;
            [x_temp2,y_temp2]=intersezione_tra_rette(a1,b1,c1,a2,b2,c2);
            x11=x_temp2-racc_pont/cos_sottr; %devo definirlo cosi per avere condione di tangenza con retta lato barriera
            y11=y_temp2;
            
            %Definisco arco di raccordo ponticelli radiali barriera di flusso
            
            % lato interno
            [a1,b1,c1]=retta_per_2pti(B1k,0,xxD1k,yyD1k);   % retta lato barriera
            a2=0;
            b2=1;
            c2=-pont/2; % retta limite ponticello
            % calcolo delle rette parallele per il calcolo del centro del
            % raccordo
            [a1p,b1p,c1p]=retta_per_2pti(B1k+racc_pont*sin(angle),0-racc_pont*cos(angle),xxD1k+racc_pont*sin(angle),yyD1k-racc_pont*cos(angle));
            a2p=0;
            b2p=1;
            c2p=-(pont/2+racc_pont);
            [xc1,yc1]=intersezione_tra_rette(a1p,b1p,c1p,a2p,b2p,c2p);  % centro dell'arco di raccordo
            xA=xc1;
            yA=pont/2;
            xB=xc1-racc_pont*sin(angle);
            yB=yc1+racc_pont*cos(angle);
            %[xB,yB]=intersezione_tra_rette(a1,b1,c1,a2p,b2p,c2p);       % punto su lato barriera
            %[xA,yA]=intersezione_tra_rette(a1p,b1p,c1p,a2,b2,c2);       % punto su lato ponticello
            % lato esterno
            [a1,b1,c1]=retta_per_2pti(B2k,0,xxD2k,yyD2k);   % retta lato barriera
            a2=0;
            b2=1;
            c2=-pont/2; % retta limite ponticello
            % calcolo delle rette parallele per il calcolo del centro del
            % raccordo
            [a1p,b1p,c1p]=retta_per_2pti(B2k-racc_pont*sin(angle),0+racc_pont*cos(angle),xxD2k-racc_pont*sin(angle),yyD2k+racc_pont*cos(angle));
            a2p=0;
            b2p=1;
            c2p=-(pont/2+racc_pont);
            [xc2,yc2]=intersezione_tra_rette(a1p,b1p,c1p,a2p,b2p,c2p);  % centro dell'arco di raccordo
            xC=xc2;
            yC=pont/2;
            xD=xc2+racc_pont*sin(angle);
            yD=yc2-racc_pont*cos(angle);
            %[xD,yD]=intersezione_tra_rette(a1,b1,c1,a2p,b2p,c2p);       % punto su lato barriera
            %[xC,yC]=intersezione_tra_rette(a1p,b1p,c1p,a2,b2,c2);       % punto su lato ponticello
            rc1=racc_pont;
            rc2=racc_pont;
            
            
            %[xc1,yc1,rc1,xA,yA,xB,yB]=racc_retteoblique(B1k,0,xxD1k,yyD1k,x10,y10,angle,racc_pont);
            %[xc2,yc2,rc2,xC,yC,xD,yD]=racc_retteoblique(B2k,0,xxD2k,yyD2k,x11,y11,angle,racc_pont);
            
            if (xC(jj) > xA(jj)) && (pont(jj)>0) %if (x11 > x10) && (pont>0)  %Se i due punti non si incrociano (raccordo non troppo grande rispetto alla
                % larghezza della barriera), e se lo spessore del ponticello non è troppo piccolo, allora procedo al disegno del ponticello
                XpontRadSx=xA;
                YpontRadSx=yA;
                
                %caso particolare#1: condizione yt<0 ordinata punto di tangenza barriera esterna è negativa (alpha< 0.25)
                %if (isnan(xD) && isnan(yD))%condizione in cui xc2,yc2,xD,YD sono NaN (caso yt<0 ordinata punto di tangenza barriera è esterna è negativa, non ho tratto lineare
                %                 if yt_old<0                           % ma solo parte raccordata come lato barriera
                %                     XpontRadDx(jj)=xxD2k(jj); %siccome il raccordo non si può disegnare, prolungo tratto rettilineo fino a intersecare
                %                     YpontRadDx(jj)=yyD2k(jj); %il lato della barriera,che nel caso yt<0 coincide con xxD2k,yyD2k, ho uno spigolo VIVO.
                %                     %raccordo destro su ponticello radiale è sostituito con spigolo vivo
                %                 else
                XpontRadDx=xC;
                YpontRadDx=yC;
                %                 end
                %caso particolare#2: yt>0 ancora positivo ma valido anche per yt<0, lato esterno barriera è composto da solo parte raccordata (alpha<<)
                if yt_old(1) < hpont(1) || yD(1) > yyD2k(1) %devo ridefinire i punti caratteristi che definiscono il raccordo destro
                    
                    %Definisco centro e raggio arco di raccordo destro su ponticello tangenziale
                    xc3=XcRibTraf2;
                    yc3=YcRibTraf2;
                    rc3=RcRibTraf2;
                    %                     dverif=sqrt((xxD2k-xc3)^2+(yyD2k-yc3)^2);
                    
                    [xE,yE]=intersezione_cerchi(xc3,yc3,rc3,xc2,yc2,rc2); %punto di intersezione tra i due archi di raccordo
                    %circonferenze secanti, ci sono due punti di intersezione, devo scegliere quello con ordinata più bassa
                    if xE(1)<xE(2)
                        XpontRadBarDx=xE(1); %i due punti vanno a coincidere con il punto di intersezione
                        YpontRadBarDx=yE(1);
                        xxD2k=xE(1); %modifico anche coordinate punto caratteristico arco di raccordo superiore destro
                        yyD2k=yE(1);
                    else
                        XpontRadBarDx=xE(2); %i due punti vanno a coincidere con il punto di intersezione
                        YpontRadBarDx=yE(2);
                        xxD2k=xE(2); %modifico anche coordinate punto caratteristico arco di raccordo superiore destro
                        yyD2k=yE(2);
                    end
                    
                else
                    XpontRadBarDx=xD;
                    YpontRadBarDx=yD;
                end
                
                XpontRadBarSx=xB;
                YpontRadBarSx=yB;
                
                %Coordinate centro e raggio raccordi ponticelli radiali
                XcRaccpontRadSx=xc1;
                YcRaccpontRadSx=yc1;
                RcRaccpontRadSx=rc1;
                XcRaccpontRadDx=xc2;
                YcRaccpontRadDx=yc2;
                RcRaccpontRadDx=rc2;
                
                
            else % Se invece i punti (x10,y10) e (x11,y11) si incrociano (ovvero raccordo troppo grande rispetto ...
                % allo spessore barriera), non disegno nessun ponticello
                pont=0;
                XpontRadBarSx=B1k;
                YpontRadBarSx=0;
                XpontRadBarDx=B2k;
                YpontRadBarDx=0;
                XpontRadSx=NaN;
                YpontRadSx=0;
                XpontRadDx=NaN;
                YpontRadDx=0;
                
                %Coordinate centro e raggio raccordi ponticelli radiali
                XcRaccpontRadSx=NaN;
                YcRaccpontRadSx=NaN;
                RcRaccpontRadSx=NaN;
                XcRaccpontRadDx=NaN;
                YcRaccpontRadDx=NaN;
                RcRaccpontRadDx=NaN;
                
            end
        end
        
    case 'Circular'
        % Valutazione ponticelli radiali:
        r_all = [];                                                                 % Nel vettore r_all si memorizzano invece le distanze, sempre rispetto al
        % Starting and ending point of the flux barrier...
        for jj = 1:nlay                                                             % che individuano l'inizio e la fine delle nlay barriere.
            r_all = [r_all x0-B2k(jj) x0-B1k(jj)];
        end
        % median point of flux barrier...
        rbeta=(r_all(1:2:end)+r_all(2:2:end))./2;
        % Posizione banane su asse q (convenzione assi VAGATI)
        XBanqdx=x0-r_all(1:2:end);
        XBanqsx=x0-r_all(2:2:end);
        
        A = zeros(1,nlay);  % carrier Area (Fe)
        Ab = A;             % barrier Area (Air or PM)
        tmp_bary = A;       % carrier center of gravity
        tmp_bary_b = A;     % barrier center of gravity
        rG = A;             % aggregate center of gravity (nlay equivalent masses, Fe + PM)
        
        % carriers (Fe): Area and Center of Gravity 
        for j = 1 : nlay
            if j == 1
                % first carrier
                [x,y] = calc_intersezione_cerchi(r,r_all(1),x0);
                dist = x0 - x;
                theta = atan2(y,dist);
                theta2 = atan2(y,x);
                A1 = 0.5*(r^2*theta2 - abs(0.5*det([x y 1; 0 0 1; x -y 1])));
                A2 = 0.5*(r_all(j)^2*theta - abs(0.5*det([x y 1; x0 0 1; x -y 1])));
                A(j) = A1+A2;
                bary1 = (2*r*(sin(theta2))^3)/(3*(theta2-sin(theta2)*cos(theta2)));
                bary2 = x0-(2*r_all(1)*(sin(theta))^3)/(3*(theta-sin(theta)*cos(theta)));
                tmp_bary(j) = (2*A1*bary1 + 2*A2*bary2)/(2*(A1+A2));
            else
                [x,y] = calc_intersezione_cerchi(r,r_all(2*j-1),x0);
                dist = x0 - x;
                theta = atan2(y,dist);
                A(j) = (r_all(2*j-1)^2 - r_all(2*j-2)^2)*theta*0.5;
                tmp_bary(j) = x0 - (2/3*sin(theta)/theta*(r_all(2*j-1)^3 - r_all(2*j-2)^3)/(r_all(2*j-1)^2 - r_all(2*j-2)^2));
            end
        end
        
        % barriers (Air or PM): Area and Center of Gravity
        for j = 1 : nlay
            [x,y] = calc_intersezione_cerchi(r,r_all(2*j),x0);
            dist = x0 - x;
            theta = atan2(y,dist);
            Ab(j) = (r_all(2*j)^2 - r_all(2*j-1)^2)*theta*0.5;
            tmp_bary_b(j) = x0 - (2/3*sin(theta)/theta*(r_all(2*j)^3 - r_all(2*j-1)^3)/(r_all(2*j)^2 - r_all(2*j-1)^2));
        end
        
        % islands as seen by the bridges
        Afe = cumsum(A);
        Abarr = cumsum(Ab);
        
        % mass and center of gravity of interest for bridges
        if (mean(unique(mat.LayerMag.Br)) > 0)
            % Fe + PM
            mass = Afe*rhoFE+Abarr*rhoPM*geo.BarFillFac;
            rG(1) = (tmp_bary(1)*Afe(1)*rhoFE+tmp_bary_b(1)*Abarr(1)*rhoPM*geo.BarFillFac)/(mass(1));
            for ii = 2 : nlay
                rG(ii)=(mass(ii-1)*rG(ii-1)+tmp_bary(ii)*A(ii)*rhoFE+tmp_bary_b(ii)*Ab(ii)*rhoPM*geo.BarFillFac)/mass(ii);  % baricentri delle regioni di ferro sorrette dalle due U
            end
        else
            % Fe + Air
            mass = Afe*rhoFE;
            rG(1) = tmp_bary(1);     
            for ii = 2 : nlay
                rG(ii)=(Afe(ii-1)*rG(ii-1)+A(ii)*tmp_bary(ii))/Afe(ii);  % baricentri delle regioni di ferro sorrette dalle due U
            end
        end
        
        M_tot = mass * 2 * l * 1e-9;   % kg
        F_centrifuga = M_tot .* rG/1000 * (nmax * pi/30)^2;
        
        if geo.radial_ribs_eval == 0
            pont = F_centrifuga/(sigma_max * l);    % mm
        else
            pont = geo.pontR;                        % input from GUI
        end
        
        pont(pont<pont0) = 0;   % pont dimension >= pont0
        
        %%% Radial Ribs Node Coordinates %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for jj = 1 : 2 : length(r_all)
            x10 = x0 - r_all(jj) - racc_pont;   % (Ponticello radiale_vertice in alto a dx)=(Punto in basso a dx del raccordo))
            x11 = x0 - r_all(jj+1) + racc_pont; % (Ponticello radiale_vertice in alto a sx)=(Punto in basso a sx del raccordo))
            if (x11 < x10) && (pont(ceil(jj/2))>0)
                % Se i due punti non si incrociano (raccordo non troppo
                % grande rispetto alla larghezza barriera)
                hpont = pont(ceil(jj/2));
                y10 = hpont/2;                     
                YpontRadBarDx(ceil(jj/2)) = y10 + racc_pont;
                [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj),0,YpontRadBarDx(ceil(jj/2)));
                XpontRadBarDx(ceil(jj/2)) = 2*x0 - x_temp;
                YpontRadBarSx(ceil(jj/2)) = y10 + racc_pont;
                [x_temp,~] = intersezione_retta_circonferenza(x0,0,r_all(jj+1),0,YpontRadBarSx(ceil(jj/2)));
                XpontRadBarSx(ceil(jj/2)) = 2*x0 - x_temp;
                XpontRadDx(ceil(jj/2))  = XpontRadBarDx(ceil(jj/2))  - racc_pont;
                YpontRadDx(ceil(jj/2))  = y10;
                XpontRadSx(ceil(jj/2))  = XpontRadBarSx(ceil(jj/2))  + racc_pont;
                YpontRadSx(ceil(jj/2))  = y10;
                
            else
                % Se invece i punti 10 e 11 si incrociano (ovvero raccordo troppo grande rispetto)
                % alla larghezza barriera), non faccio il ponticello e disegno solo una linea
                hpont = 0;
                XpontRadBarSx(ceil(jj/2)) = XBanqsx(ceil(jj/2));
                YpontRadBarSx(ceil(jj/2)) = 0;
                XpontRadBarDx(ceil(jj/2)) = XBanqdx(ceil(jj/2));
                YpontRadBarDx(ceil(jj/2)) = 0;
                XpontRadDx(ceil(jj/2)) = NaN;
                YpontRadDx(ceil(jj/2)) = 0;
                XpontRadSx(ceil(jj/2)) = NaN;
                YpontRadSx(ceil(jj/2)) = 0;
            end
        end
        
    otherwise
        %%  Seg - calcolo area Ferro
        [xrot_traf,yrot_traf]=calc_intersezione_cerchi(r, rTemp, x0);
        Dx=(r-xrot_traf(1))/5;  %per avere almeno 5 divisioni;
        xcir_plot=[r:-Dx:r*cos(pi/2/p)];
        ycir_plot=sqrt(r^2-xcir_plot.^2);
        VectCir=find(xcir_plot>=xrot_traf(1));
        x_ext_rot=xcir_plot(VectCir);
        y_ext_rot=ycir_plot(VectCir);
        
        A=[];
        for ii=1:nlay
            if ii==1
                X=[B2k(ii), XpBar2(ii), xxD2k(ii),xpont(ii),fliplr(x_ext_rot)];
                Y=[0, YpBar2(ii), yyD2k(ii),ypont(ii),fliplr(y_ext_rot)];
                %         figure(100);hold on;fill(X,Y,'r');hold off;
                A(ii)=polyarea(X,Y);
                tmp_bary(ii,:)=centroid(X',Y');         %MarcoP
                clear X Y;
            else
                X=[B1k(ii-1), XpBar1(ii-1), xxD1k(ii-1),xpont(ii-1),xrot_traf(ii-1),xrot_traf(ii),xpont(ii),xxD2k(ii),XpBar2(ii),B2k(ii)];
                Y=[0, YpBar1(ii-1), yyD1k(ii-1),ypont(ii-1),yrot_traf(ii-1),yrot_traf(ii),ypont(ii),yyD2k(ii),YpBar2(ii),0];
                %         figure(100);hold on;fill(X,Y,'r'); hold off;
                A(ii)=polyarea(X,Y);
                tmp_bary(ii,:)=centroid(X',Y');         %MarcoP
                clear X Y;
            end
        end
        Afe=cumsum(A);
        
        %MarcoP
        rG(1)=tmp_bary(1,1);     % baricentro della regione di ferro sorretta dal ponticello a I
        for ii=2:nlay
            rG(ii)=(Afe(ii-1)*rG(ii-1)+A(ii)*tmp_bary(ii,1))/Afe(ii);  % baricentri delle regioni di ferro sorrette dalle due U
        end
        
        M_Fe = 2*Afe*l * 1e-9 * rhoFE ;   % massa ferro appeso ai ponticelli
        
        %% Seg - calcolo area magnete
        % barriera ad I
        X_Ibarr=[xxD1k(1) xxD2k(1) B2k(1) B1k(1)];
        Y_Ibarr=[yyD1k(1) yyD2k(1) 0 0];
        areaI=2*polyarea(X_Ibarr,Y_Ibarr);  % area della barriera a I in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e gli eventuali ponticelli)
        % % barriera ad U centrale
        if nlay == 1
            area_barr_withoutPM = [areaI];
        else
            % barriere ad U
            for kk=2:nlay
                areaUbars(kk-1) = (2*YpBar2(kk)*hc(kk))+(2*hc(kk)*abs((xxD2k(kk)+j*yyD2k(kk))-(XpBar2(kk)+j*YpBar2(kk)))); %area della barriera a U in cui inserire i magneti (escludo i due semicerchi all'estremità della barriera e le due "gambe" della "U")
            end
            area_barr_withoutPM = [areaI areaUbars];
        end
        rG_PM = (B1k+B2k)/2;    % vettore con i baricentri dei magneti
        area_barr_withPM=area_barr_withoutPM;
        A_PM=cumsum(area_barr_withPM);
        M_PM=A_PM*rhoPM*1e-9*2*l;
        
        %% Seg - calcolo e disegno ponticelli
        F_centrifuga = (M_Fe+M_PM) .* rG/1000 *  (nmax * pi/30)^2;
        
        if geo.radial_ribs_eval == 0
            pont = F_centrifuga/(sigma_max * l);    % mm
        else
            pont = geo.pontR;
        end
        
        XpontRadBarSx = zeros(1,nlay);
        YpontRadBarSx = zeros(1,nlay);
        XpontRadBarDx = zeros(1,nlay);
        YpontRadBarDx = zeros(1,nlay);
        XpontRadDx    = zeros(1,nlay);
        YpontRadDx    = zeros(1,nlay);
        XpontRadSx    = zeros(1,nlay);
        YpontRadSx    = zeros(1,nlay);
        
        XpontSplitBarSx = nan(2,nlay);
        YpontSplitBarSx = nan(2,nlay);
        XpontSplitBarDx = nan(2,nlay);
        YpontSplitBarDx = nan(2,nlay);
        XpontSplitDx    = nan(2,nlay);
        YpontSplitDx    = nan(2,nlay);
        XpontSplitSx    = nan(2,nlay);
        YpontSplitSx    = nan(2,nlay);
        
        
        if geo.radial_ribs_split
            % Controllo spessore ponticello:
            % -        pont<pont0   --> no ponticello
            % - pont0<=pont<2*pont0 --> disegno 2 ponticelli spessi pont0
            % - pont <=2*pont0      --> disegno 2 ponticelli spessi pont/2
            for jj=1:nlay
                if pont(jj)<pont0
                    pont(jj)=0;
                elseif pont(jj)<2*pont0
                    pont(jj)=2*pont0;
                end
            end
            
            hpont=pont/2;
            
            
            for ii=1:nlay
                if hpont(ii)>0
                    XpontSplitBarSx(1,ii)=B1k(ii);
                    XpontSplitBarSx(2,ii)=B1k(ii);
                    YpontSplitBarSx(1,ii)=YpBar2(ii);
                    YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii)-2*racc_pont-hpont(ii);
                    XpontSplitBarDx(1,ii)=B2k(ii);
                    XpontSplitBarDx(2,ii)=B2k(ii);
                    YpontSplitBarDx(1,ii)=YpBar2(ii);
                    YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii)-2*racc_pont-hpont(ii);
                    XpontSplitSx(1,ii)=B1k(ii)+racc_pont;
                    XpontSplitSx(2,ii)=B1k(ii)+racc_pont;
                    YpontSplitSx(1,ii)=YpBar2(ii)-racc_pont;
                    YpontSplitSx(2,ii)=YpontSplitSx(1,ii)-hpont(ii);
                    XpontSplitDx(1,ii)=B2k(ii)-racc_pont;
                    XpontSplitDx(2,ii)=B2k(ii)-racc_pont;
                    YpontSplitDx(1,ii)=YpBar2(ii)-racc_pont;
                    YpontSplitDx(2,ii)=YpontSplitDx(1,ii)-hpont(ii);
                else
                    XpontSplitBarSx(1,ii)=B1k(ii);
                    XpontSplitBarSx(2,ii)=B1k(ii);
                    YpontSplitBarSx(1,ii)=YpBar2(ii);
                    YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii);
                    XpontSplitBarDx(1,ii)=B2k(ii);
                    XpontSplitBarDx(2,ii)=B2k(ii);
                    YpontSplitBarDx(1,ii)=YpBar2(ii);
                    YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii);
                    XpontSplitSx(1,ii)=NaN;
                    XpontSplitSx(2,ii)=NaN;
                    YpontSplitSx(1,ii)=NaN;
                    YpontSplitSx(2,ii)=NaN;
                    XpontSplitDx(1,ii)=NaN;
                    XpontSplitDx(2,ii)=NaN;
                    YpontSplitDx(1,ii)=NaN;
                    YpontSplitDx(2,ii)=NaN;
                end
                
                XpontRadBarSx(ii)=B1k(ii);
                YpontRadBarSx(ii)=0;
                XpontRadBarDx(ii)=B2k(ii);
                YpontRadBarDx(ii)=0;
                XpontRadDx(ii)=NaN;
                YpontRadDx(ii)=0;
                XpontRadSx(ii)=NaN;
                YpontRadSx(ii)=0;
            end
            
        else
            for jj=1:nlay
                if (pont(jj) < pont0) % non disegno i ponticelli radiali il cui spessore è minore della tolleranza di lavorazione per gli altri tipi di rotore
                    pont(jj)=0;
                end
            end
            
            hpont=pont/2;
            rac_pont=abs(B1k-B2k)/4;
            rac_pont=pont0;
            
            for ii=1:nlay
                if hpont(ii)>0
                    XpontRadBarSx(ii)=B1k(ii);
                    YpontRadBarSx(ii)=hpont(ii)+racc_pont;
                    XpontRadBarDx(ii)=B2k(ii);
                    YpontRadBarDx(ii)=hpont(ii)+racc_pont;
                    XpontRadDx(ii)=B2k(ii)-racc_pont;
                    YpontRadDx(ii)=hpont(ii);
                    XpontRadSx(ii)=B1k(ii)+racc_pont;
                    YpontRadSx(ii)=hpont(ii);
                else
                    XpontRadBarSx(ii)=B1k(ii);
                    YpontRadBarSx(ii)=0;
                    XpontRadBarDx(ii)=B2k(ii);
                    YpontRadBarDx(ii)=0;
                    XpontRadDx(ii)=NaN;
                    YpontRadDx(ii)=0;
                    XpontRadSx(ii)=NaN;
                    YpontRadSx(ii)=0;
                end
                
                XpontSplitBarSx(1,ii)=B1k(ii);
                XpontSplitBarSx(2,ii)=B1k(ii);
                YpontSplitBarSx(1,ii)=YpBar2(ii);
                YpontSplitBarSx(2,ii)=YpontSplitBarSx(1,ii);
                XpontSplitBarDx(1,ii)=B2k(ii);
                XpontSplitBarDx(2,ii)=B2k(ii);
                YpontSplitBarDx(1,ii)=YpBar2(ii);
                YpontSplitBarDx(2,ii)=YpontSplitBarDx(1,ii);
                XpontSplitSx(1,ii)=NaN;
                XpontSplitSx(2,ii)=NaN;
                YpontSplitSx(1,ii)=NaN;
                YpontSplitSx(2,ii)=NaN;
                XpontSplitDx(1,ii)=NaN;
                XpontSplitDx(2,ii)=NaN;
                YpontSplitDx(1,ii)=NaN;
                YpontSplitDx(2,ii)=NaN;
                
            end
        end
end

