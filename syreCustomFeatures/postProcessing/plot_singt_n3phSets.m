function [out,wf] = plot_singt_n3phSets(pathname,filename)

if nargin()<2
    [filename, pathname, ~] = uigetfile([cd '\.mat'], 'LOAD DATA');
end

load([pathname filename])


% multi-3-phase computation
% if geo.win.n3phase>1
for ii=1:geo.win.n3phase
    th = out.SOL.th;

    fa = out.SOL.fa(ii,:);
    fb = out.SOL.fb(ii,:);
    fc = out.SOL.fc(ii,:);
    fdq = abc2dq(fa,fb,fc,(th+geo.th0(ii)-geo.th0(1))*pi/180);
    out.SOL.sets(ii).fa = fa;
    out.SOL.sets(ii).fb = fb;
    out.SOL.sets(ii).fc = fc;
    out.SOL.sets(ii).fd = fdq(1,:);
    out.SOL.sets(ii).fq = fdq(2,:);
    out.SOL.sets(ii).f0 = (fa+fb+fc)/3;

    ia = out.SOL.ia(ii,:);
    ib = out.SOL.ib(ii,:);
    ic = out.SOL.ic(ii,:);
    idq = abc2dq(ia,ib,ic,(th+geo.th0(ii)-geo.th0(1))*pi/180);
    out.SOL.sets(ii).ia = ia;
    out.SOL.sets(ii).ib = ib;
    out.SOL.sets(ii).ic = ic;
    out.SOL.sets(ii).id = idq(1,:);
    out.SOL.sets(ii).iq = idq(2,:);
    out.SOL.sets(ii).i0 = (ia+ib+ic)/3;
end
% end



% average computations
if isfield(out.SOL,'sets')
    for ii=1:length(out.SOL.sets)
        out.(['fd' int2str(ii)]) = mean(out.SOL.sets(ii).fd);
        out.(['fq' int2str(ii)]) = mean(out.SOL.sets(ii).fq);
    end
end

% wf - refers to complete waveform
nRep = 360/per.delta_sim_singt; % number of repetition needed

wf.fd       = repmat(out.SOL.fd,1,nRep);
wf.fq       = repmat(out.SOL.fq,1,nRep);
wf.T        = repmat(out.SOL.T,1,nRep);
wf.th       = linspace(0,360,length(wf.fd)+1);
wf.th       = wf.th(1:end-1);

fph = phaseQuantityDecoding(out.SOL.fa,out.SOL.fb,out.SOL.fc,per.delta_sim_singt);
wf.fa = [fph.a fph.a(:,1)];
wf.fb = [fph.b fph.b(:,1)];
wf.fc = [fph.c fph.c(:,1)];

if isfield(out.SOL,'sets')
    for ii=1:length(out.SOL.sets)
        wf.sets(ii).fd = repmat(out.SOL.sets(ii).fd,1,nRep);
        wf.sets(ii).fq = repmat(out.SOL.sets(ii).fq,1,nRep);
        wf.sets(ii).id = repmat(out.SOL.sets(ii).id,1,nRep);
        wf.sets(ii).iq = repmat(out.SOL.sets(ii).iq,1,nRep);
        
        tmp = phaseQuantityDecoding(out.SOL.sets(ii).fa,out.SOL.sets(ii).fb,out.SOL.sets(ii).fc,per.delta_sim_singt);
        wf.sets(ii).fa = tmp.a;
        wf.sets(ii).fb = tmp.b;
        wf.sets(ii).fc = tmp.c;
        wf.sets(ii).f0 = (tmp.a+tmp.b+tmp.c)/3;

        tmp = phaseQuantityDecoding(out.SOL.sets(ii).ia,out.SOL.sets(ii).ib,out.SOL.sets(ii).ic,per.delta_sim_singt);
        wf.sets(ii).ia = tmp.a;
        wf.sets(ii).ib = tmp.b;
        wf.sets(ii).ic = tmp.c;
        wf.sets(ii).i0 = (tmp.a+tmp.b+tmp.c)/3;
    end
end

% plots

figNames{1} = 'n3ph_wf_ABCcurrentsVSpos';
figNames{2} = 'n3ph_wf_ABCfluxLinkagesVSpos';
figNames{3} = 'n3ph_wf_DQcurrentsVSpos';
figNames{4} = 'n3ph_wf_DQfluxLinkagesVSpos';
figNames{5} = 'n3ph_wf_ABScurrentsVSpos';
figNames{6} = 'n3ph_wf_ABSfluxLinkagesVSpos';
figNames{7} = 'n3ph_wf_torqueVSpos';
figNames{8} = 'n3ph_vector_currents';
figNames{9} = 'n3ph_vector_fluxLinkages';


for ii=1:length(figNames)
    hfig(ii) = figure();
    figSetting(14,14)
    hax(ii) = axes('OuterPosition',[0 0 1 1]);
    set(hfig(ii),'FileName',[pathname figNames{ii} '.fig'],'Name',figNames{ii})
    switch ii
        case 1
            xlabel('$\theta_r$ (elt deg)')
            ylabel('$i_{abc}$ (A)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 2
            xlabel('$\theta_r$ (elt deg)')
            ylabel('$\lambda_{abc}$ (Vs)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 3
            xlabel('$\theta_r$ (elt deg)')
            ylabel('$i_{dq}$ (A)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 4
            xlabel('$\theta_r$ (elt deg)')
            ylabel('$\lambda_{dq}$ (Vs)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 5
            xlabel('$\theta_r$ (elt deg)')
            ylabel('|$i_{dq}|$ (A)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 6
            xlabel('$\theta_r$ (elt deg)')
            ylabel('$|\lambda_{dq}|$ (Vs)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 7
            xlabel('$\theta_r$ (elt deg)')
            ylabel('$T$ (Nm)')
            set(gca,'XLim',[0 360],'XTick',0:60:360);
        case 8
            set(gca,'DataAspectRatio',[1 1 1]);
            xlabel('$i_d$ (A)')
            ylabel('$i_q$ (A)')
        case 9
            set(gca,'DataAspectRatio',[1 1 1]);
            xlabel('$\lambda_d$ (A)')
            ylabel('$\lambda_q$ (A)')
    end
    hleg(ii) = legend(hax(ii),'show');
end

colors = get(hax(1),'ColorOrder');

for ii=1:length(wf.sets)
    % fig 1 - phase currents
    plot(hax(1),wf.th,wf.sets(ii).ia,'Color',colors(ii,:),'LineStyle','-','DisplayName',['$i_{a,' int2str(ii) '}$'])
    plot(hax(1),wf.th,wf.sets(ii).ib,'Color',colors(ii,:),'LineStyle','--','DisplayName',['$i_{b,' int2str(ii) '}$'])
    plot(hax(1),wf.th,wf.sets(ii).ic,'Color',colors(ii,:),'LineStyle',':','DisplayName',['$i_{c,' int2str(ii) '}$'])
    % fig2 - phase flux linkages
    plot(hax(2),wf.th,wf.sets(ii).fa,'Color',colors(ii,:),'LineStyle','-','DisplayName',['$\lambda_{a,' int2str(ii) '}$'])
    plot(hax(2),wf.th,wf.sets(ii).fb,'Color',colors(ii,:),'LineStyle','--','DisplayName',['$\lambda_{b,' int2str(ii) '}$'])
    plot(hax(2),wf.th,wf.sets(ii).fc,'Color',colors(ii,:),'LineStyle',':','DisplayName',['$\lambda_{c,' int2str(ii) '}$'])
    % fig3 - dq0 currents
    plot(hax(3),wf.th,wf.sets(ii).id,'Color',colors(ii,:),'LineStyle','-','DisplayName',['$i_{d,' int2str(ii) '}$'])
    plot(hax(3),wf.th,wf.sets(ii).iq,'Color',colors(ii,:),'LineStyle','--','DisplayName',['$i_{q,' int2str(ii) '}$'])
    plot(hax(3),wf.th,wf.sets(ii).i0,'Color',colors(ii,:),'LineStyle',':','DisplayName',['$i_{0,' int2str(ii) '}$'])
    % fig4 - dq0 flux linkages
    plot(hax(4),wf.th,wf.sets(ii).fd,'Color',colors(ii,:),'LineStyle','-','DisplayName',['$\lambda_{d,' int2str(ii) '}$'])
    plot(hax(4),wf.th,wf.sets(ii).fq,'Color',colors(ii,:),'LineStyle','--','DisplayName',['$\lambda_{q,' int2str(ii) '}$'])
    plot(hax(4),wf.th,wf.sets(ii).f0,'Color',colors(ii,:),'LineStyle',':','DisplayName',['$\lambda_{0,' int2str(ii) '}$'])
    % fig5 - currents amplitude
    plot(hax(5),wf.th,abs(wf.sets(ii).id+j*wf.sets(ii).iq),'-','Color',colors(ii,:),'DisplayName',['set ' int2str(ii)])
    % fig6 - flux linkages amplitude
    plot(hax(6),wf.th,abs(wf.sets(ii).fd+j*wf.sets(ii).fq),'-','Color',colors(ii,:),'DisplayName',['set ' int2str(ii)])
    % fig7 - torque
    % plot(hax(7),wf.th,3/2*geo.p*(wf.sets(ii).fd.*wf.sets(ii).iq-wf.sets(ii).fq.*wf.sets(ii).id),'-','Color',colors(ii,:),'DisplayName',['set ' int2str(ii)])
    % fig8 - current vector
    plot(hax(8),[0 mean(wf.sets(ii).id)],[0 mean(wf.sets(ii).iq)],'Color',colors(ii,:),'Marker','o','MarkerIndices',2,'DisplayName',['set ' int2str(ii)])
    % fig9 - flux linkage vector
    plot(hax(9),[0 mean(wf.sets(ii).fd)],[0 mean(wf.sets(ii).fq)],'Color',colors(ii,:),'Marker','o','MarkerIndices',2,'DisplayName',['set ' int2str(ii)])
end

plot(hax(7),wf.th,wf.T,'-k','DisplayName','total')

save([pathname 'output_singleSet.mat'],'out','geo','per','mat','wf')


if nargout()==0
    clear out
end


