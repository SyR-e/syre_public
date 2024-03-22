function draw_rotor_in_MCAD(mcad,geo,~,~)

switch geo.RotType
    case 'Circular'
        invoke(mcad,'SetVariable','MotorType_MotorLAB','SYNCREL');
        invoke(mcad,'SetVariable','BPMRotor','13'); %U
        % rotor parameters
        invoke(mcad,'SetVariable','Shaft_Dia',geo.Ar*2);
        invoke(mcad,'SetVariable','Pole_number',geo.p*2);

        invoke(mcad,'SetVariable','Magnet_Layers',int2str(geo.nlay));  %number of magnetic layers

        %centre posts
        tmp=fliplr(geo.pontR);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_CentrePost_Array',tmp);

        nlay=geo.nlay;
        mg_leng=zeros(1,nlay);
        mg_leng=mat2str(mg_leng);
        mg_leng=mg_leng(2:end-1);
        mg_leng(mg_leng==' ')=':';
        invoke(mcad,'SetVariable','UMagnet_Length_Inner_Array',mg_leng);
        invoke(mcad,'SetVariable','UMagnet_Length_Outer_Array',mg_leng);

        alpha=cumsum(geo.dalpha);
        delta = [alpha(1) diff(alpha) (90/geo.p)-alpha(end)];

        if geo.p==2
            tmp=(-delta(2)).*ones(1,nlay);
        elseif geo.p==3
            tmp=(-delta(2)*2).*ones(1,nlay);
        elseif geo.p==4
            tmp=(-delta(2)*3.8).*ones(1,nlay);
        else
            tmp=(-delta(2)).*ones(1,nlay);
        end
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_OuterAngleOffset_Array',tmp);

        tmp=fliplr(geo.hc);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_Thickness_Inner_Array',tmp);
        invoke(mcad,'SetVariable','UShape_Thickness_Outer_Array',tmp);

        %         hc=geo.hc(1);
        tmp=0;
        for i=1:1:nlay
            tmp(i)=2*geo.B1k(nlay-i+1);
        end
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_InnerDiameter_Array',tmp);

        tmp=fliplr(geo.pontT);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_BridgeThickness_Array',tmp);

        if geo.p==2
            tmp=sqrt(2)/2*geo.xxD1k-sqrt(2)/2*geo.yyD1k;
        else
            m=tan(pi/(2*geo.p));
            tmp=(abs(geo.yyD1k-m*geo.xxD1k))/(sqrt(1+m^2));
        end
        tmp=2*fliplr(tmp);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_WebThickness_Array',tmp);

    case 'Seg'
        invoke(mcad,'SetVariable','MotorType_MotorLAB','BPM');
        if any(geo.hcShrink<1)
            invoke(mcad,'SetVariable','BPMRotor','13'); %U
        else
            invoke(mcad,'SetVariable','BPMRotor','11'); %interior V (web)
        end

        tmp=geo.Ar*2;
        invoke(mcad,'SetVariable','Shaft_Dia',tmp);
        tmp=geo.p*2;
        invoke(mcad,'SetVariable','Pole_number',tmp);

        invoke(mcad,'SetVariable','Magnet_Layers',int2str(geo.nlay));  %number of magnetic layers

        %centre posts
        tmp=fliplr(geo.pontR);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_CentrePost_Array',tmp);

        nlay=geo.nlay;
        mg_leng=zeros(1,nlay);
        mg_leng=mat2str(mg_leng);
        mg_leng=mg_leng(2:end-1);
        mg_leng(mg_leng==' ')=':';
        invoke(mcad,'SetVariable','UMagnet_Length_Inner_Array',mg_leng);
        invoke(mcad,'SetVariable','UMagnet_Length_Outer_Array',mg_leng);

        tmp=zeros(1,nlay);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_OuterAngleOffset_Array',tmp);

        tmp=fliplr(geo.hc);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_Thickness_Inner_Array',tmp);
        invoke(mcad,'SetVariable','UShape_Thickness_Outer_Array',tmp);

        %         hc=geo.hc(1);
        tmp=0;
        for i=1:1:nlay
            tmp(i)=2*geo.B1k(nlay-i+1);
        end
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_InnerDiameter_Array',tmp);

        %%%%%%
        tmp=fliplr(geo.pontT+0.5);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_BridgeThickness_Array',tmp);

        if geo.p==2
            tmp=sqrt(2)*(geo.B1k-geo.YpBar1);
        else
            m=tan(pi/(2*geo.p));
            tmp=2*(abs(geo.YpBar1-m*geo.B1k))/(sqrt(1+m^2));
        end
        tmp=fliplr(tmp);
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UShape_WebThickness_Array',tmp);

        %radial ribs
        if geo.radial_ribs_split==1
            tmp=fliplr(geo.pontR/2);
            tmp=mat2str(tmp);
            tmp=tmp(2:end-1);
            tmp(tmp==' ')=':';
            invoke(mcad,'SetVariable','UShape_Post_Inner_Array',tmp);

            tmp=zeros(1,geo.nlay);
            tmp=mat2str(tmp);
            tmp=tmp(2:end-1);
            tmp(tmp==' ')=':';
            invoke(mcad,'SetVariable','UShape_CentrePost_Array',tmp);
        else
            tmp=fliplr(geo.pontR);
            tmp=mat2str(tmp);
            tmp=tmp(2:end-1);
            tmp(tmp==' ')=':';
            invoke(mcad,'SetVariable','UShape_CentrePost_Array',tmp);
        end

        %magnet
        tmp=fliplr(geo.PMdim(1,:));
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UMagnet_Length_Inner_Array',tmp);

        tmp=fliplr(geo.PMdim(2,:));
        tmp=mat2str(tmp);
        tmp=tmp(2:end-1);
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','UMagnet_Length_Outer_Array',tmp);

    case 'SPM'
        invoke(mcad,'SetVariable','BPMRotor','0'); %surface radial

        tmp=geo.Ar*2;
        invoke(mcad,'SetVariable','Shaft_Dia',tmp);

        tmp=geo.hc-geo.betaPMshape*geo.hc;   %bombatura magneti
        invoke(mcad,'SetVariable','MagnetReduction',tmp);

        tmp=geo.p*2;
        invoke(mcad,'SetVariable','Pole_number',tmp);

        invoke(mcad,'SetVariable','Airgap',geo.g);
        invoke(mcad,'SetVariable','Magnet_Thickness',geo.hc_pu);
        invoke(mcad,'SetVariable','Magnet_Arc_[ED]',geo.dalpha);
        invoke(mcad,'SetVariable','Magnet_Thickness',geo.hc);
        invoke(mcad,'SetVariable','CircumferentialSegments','2');

    case 'Vtype'
        invoke(mcad,'SetVariable','MotorType_MotorLAB','BPM');
        invoke(mcad,'SetVariable','BPMRotor','11'); %interior V (web)

        tmp=geo.Ar*2;
        invoke(mcad,'SetVariable','Shaft_Dia',tmp);
        tmp=geo.p*2;
        invoke(mcad,'SetVariable','Pole_number',tmp);

        tmp=geo.nlay;
        invoke(mcad,'SetVariable','VMagnet_Layers',tmp);

        tmp=fliplr(geo.hc);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','MagnetThickness_Array',tmp);

        tmp=fliplr(geo.pontT);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','BridgeThickness_Array',tmp);

        tmp=2*atan(geo.yyD1k./(geo.xxD1k-geo.B1k))*180/pi;
        tmp=fliplr(tmp);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','PoleVAngle_Array',tmp);

        tmp=fliplr(geo.pontR);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','VShapeMagnetPost_Array',tmp);

        tmp=fliplr(2*geo.yPMC2b);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','MagnetSeparation_Array',tmp);

        tmp=zeros(1,geo.nlay);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','WebLength_Array',tmp);

        %web thickness (at the barrier end point)
        tmp=atan(geo.yyD1k./(geo.xxD1k-geo.B1k));
        add_x=geo.hc.*cos(tmp);
        add_y=geo.hc.*sin(tmp);
        add_x(end)=add_x(1)*0.7;
        add_y(end)=add_y(1)*0.7;
        xxD1k_web=geo.xxD1k+add_x.*1.1;
        yyD1k_web=geo.yyD1k+add_y.*1.1;
        m=tan(pi/(2*geo.p));
        tmp=2*(abs(yyD1k_web-m*xxD1k_web))./(sqrt(1+m^2));
        tmp=fliplr(tmp);
        if geo.p==3
            if geo.nlay>2
                tmp(1)=tmp(1)*1.2;
            else
                tmp(1)=tmp(1)*1.1;
            end
        end
        if geo.p==4
            if geo.nlay>2
                tmp(1)=tmp(1)*1.47;
            else
                tmp(1)=tmp(1)*1.1;
            end
        end
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','WebThickness_Array',tmp);

        %tangential ribs
        m=yyD1k_web./xxD1k_web;
        ypont=(geo.r.*sin(atan(m)));
        xpont=(geo.r.*cos(atan(m)));
        pontT_mcad=sqrt((xpont-xxD1k_web).*(xpont-xxD1k_web)+(ypont-yyD1k_web).*(ypont-yyD1k_web));
        tmp=fliplr(pontT_mcad);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','BridgeThickness_Array',tmp);

        m_c=-geo.xxD2k./geo.yyD2k;
        alpha=atan(-1./m_c);
        x_c=geo.r*cos(alpha);
        y_c=geo.r*sin(alpha);
        q_c=y_c-x_c.*m_c;
        m_b=(geo.yyD1k-geo.yPMC1b)./(geo.xxD1k-geo.xPMC1b);
        q_b=geo.yyD2k-m_b.*geo.xxD2k;
        xPoleArc=(-pontT_mcad.*sqrt(1+m_c.*m_c)-q_b+q_c)./(m_b-m_c);
        yPoleArc=m_b.*xPoleArc+q_b;
        tmp=fliplr(2*geo.p.*atan(yPoleArc./xPoleArc)*180/pi);
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','PoleArc_Array',tmp);

        %magnet
        tmp=fliplr(geo.PMdim(1,:));
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','MagnetBarWidth_Array',tmp);

        tmp=fliplr(geo.PMclear(1,:));
        tmp=mat2str(tmp);
        if geo.nlay>1
            tmp=tmp(2:end-1);
        end
        tmp(tmp==' ')=':';
        invoke(mcad,'SetVariable','VShapeMagnetClearance_Array',tmp);

        %corner rounding
        invoke(mcad,'SetVariable','CornerRounding_Rotor','1');
        tmp=geo.hc(1)/2;
        invoke(mcad,'SetVariable','CornerRoundingRadius_Rotor',tmp);
end