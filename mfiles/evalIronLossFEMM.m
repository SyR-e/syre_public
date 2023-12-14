% Copyright 2020
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function [SOL] = evalIronLossFEMM(geo,per,mat,SOL,method)
%
% [SOL] = evalIronLossFEMM(geo,mat,per,SOL,method)
%
% Compute iron loss from FEMM mesh data, according to the variable method.
% method = 0 --> FEMM method, as the tutorial on the website. Flux density
%                is divided in all the harmonic component and the Steinmetz
%                equation is applied for each frequency, both for
%                hysteresis and eddy current loss
% method = 1 --> iGSE is used to compute the hysteresis loss, by dividing
%                the minor loops from the major hysteresis loop. Eddy
%                current loss are computed as the previous method. For the
%                hystereis loops detection, just the main flux density
%                direction for each elements is considered
% method = 2 --> iGSE, as method=1, but for each element, the loss are
%                computed for the main flux density direction and the
%                ortogonal direction
% method = 3 --> iGSE, as method=1, but for each element, the loss are
%                computed for the main flux density direction and the
%                ortogonal direction, using the R factor introduced in
%                2004-Bottauscio-"Additional loss in IM under synchronous
%                no-load condition", with the curve R=f(B) reported in
%                2012-Dlala-"Numerical Investigation of the Effects of
%                Loading and Short Harmonics on the Core Loss of IMs".
%
% Default method: 2
% IT IS STRONGLY suggested to use the case 2 since it embeds the DC bias
% effect on the iron loss as well as the reaction field effect and the
% 3D effect on the magnet loss!


% load data from input structures
p       = geo.p;
if strcmp(geo.RotType,'Seg')
    PMdim0 = geo.PMdim;
    PMdim0(1,:) = PMdim0(1,:)*2;
    PMdim1   = PMdim0(PMdim0>0);
    
    PMNc0      = geo.PMNc;
    PMNc0(1,:) = PMNc0(1,:)/2;
    PMNc1      = PMNc0(geo.PMdim>0);
    PMNc1(PMNc1<1) = 1;
    PMclear1   = geo.PMclear(geo.PMdim>0);
    hc1        = [geo.hc(geo.PMdim(1,:)>0)'; geo.hc(geo.PMdim(2,:)>0)'];

    PMdim   =[];
    PMclear =[];
    hc      =[];
    PMNc    =[];
    for ii=1:length(PMNc1)
        PMdim   = [PMdim; repmat(PMdim1(ii),[PMNc1(ii)*2 1])];
        PMclear = [PMclear; repmat(PMclear1(ii),[PMNc1(ii)*2 1])];
        PMNc    = [PMNc; repmat(PMNc1(ii),[PMNc1(ii)*2 1])];
        hc      = [hc; repmat(hc1(ii),[PMNc1(ii)*2 1])];
    end
%     PMdim   = [PMdim; PMdim];
%     PMNc    = [PMNc; PMNc];
%     PMclear = [PMclear; PMclear];
%     hc      = [hc; hc];
end
l       = geo.l/geo.PMNa;


speed   = per.EvalSpeed;

khS     = mat.Stator.kh;
alphaS  = mat.Stator.alpha;
betaS   = mat.Stator.beta;
keS     = mat.Stator.ke;
rhoS    = mat.Stator.kgm3;
khR     = mat.Rotor.kh;
alphaR  = mat.Rotor.alpha;
betaR   = mat.Rotor.beta;
keR     = mat.Rotor.ke;
rhoR    = mat.Rotor.kgm3;
sigmaPM = mat.LayerMag.sigmaPM;
mu      = mat.LayerMag.mu*4*pi*1e-7;

vol     = SOL.vol;
bs      = SOL.bs;
br      = SOL.br;
groNo   = SOL.groNo;
am      = SOL.am;

% from 180° data to 360°
switch per.delta_sim_singt
    case 360
        th = SOL.th;
        bs = SOL.bs;
        br = SOL.br;
        am = SOL.am;
    case 180
        th = [SOL.th SOL.th+180];
        bs = [bs; -bs];
        br = [br; br];
        am = [am; am];
    otherwise
        error('FEA simulation span not correct!!!')
end

nsim    = length(th);

EleNo   = length(groNo);
numPM   = max(groNo)-200;
f0      = speed*p/60; % fundamental frequency in Hz

%% Cases
switch method
    case 0        % FEMM method: harmonic decomposition and sum

        % elaboration of the flux densities for the iron loss computation
        bsx = abs(fft(real(bs)))*2/nsim;
        bsy = abs(fft(imag(bs)))*2/nsim;
        brx = abs(fft(real(br)))*2/nsim;
        bry = abs(fft(imag(br)))*2/nsim;
        Bs = (bsx.^2+bsy.^2).^0.5;
        Br = (brx.^2+bry.^2).^0.5;

        % computation of of the fft of A, at the center of each element, and
        % subtracting the average value of each PM piece
        Jm = (fft(am))*2/nsim;
        Jo = zeros(size(Jm));
        for mm=1:numPM
            index = 1:1:EleNo;
            index = index(groNo==(200+mm));
            vmag = sum(vol(index));
            tmp = (Jm(:,index).*vol(index))/vmag;
            tmp = sum(tmp,2);
            Jo(:,index) = repmat(tmp,1,length(index));
            Jm(:,index) = Jm(:,index)-Jo(:,index);
        end

        % loss computation
        freq = (0:1:(nsim-1))*f0;   % frequency associated to each harmonic component
        freq(floor(nsim/2):end) = 0;
        freq = repmat(freq.',1,EleNo);
        vol  = repmat(vol,nsim,1);

        index = 1:1:EleNo;
        psh = (khS*freq(:,index).^alphaS.*Bs(:,index).^betaS)*rhoS.*vol(:,index);
        psc = (keS*freq(:,index).^2.*Bs(:,index).^2)*rhoS.*vol(:,index);

        index = 1:1:EleNo;
        prh = (khR*freq(:,index).^alphaR.*Br(:,index).^betaR)*rhoR.*vol(:,index);
        prc = (keR*freq(:,index).^2.*Br(:,index).^2)*rhoR.*vol(:,index);

        index = 1:1:EleNo;
        ppm = 1/2*sigmaPM*(2*pi*freq(:,index)).^2.*abs(Jm(:,index)).^2.*vol(:,index);

    case 1    % iGSE+eddy current, just the main axis
        % hysteresis with minor loops

        th = linspace(0,2*pi,1001);

        kiS = khS/((2*pi)^(alphaS-1)*trapz(th,abs(cos(th)).^alphaS*2^(betaS-alphaS)));
        kiR = khR/((2*pi)^(alphaR-1)*trapz(th,abs(cos(th)).^alphaR*2^(betaR-alphaR)));

        % stator
        index = 1:1:EleNo;
        index = index(groNo==12);
        nS = numel(index);

        psh = zeros(1,EleNo);

        for ii=1:nS
            x = th;
            y = bs(:,index(ii));
            [major,minor] = minorLoopDetection(x,y);
            if ~isempty(major.y)
                dt = (major.x(2)-major.x(1))*(pi/180)/(2*pi*f0);
                dB = diff([major.y major.y(1)]);
                DB = max(major.y)-min(major.y);
                Pmajor = 1/(dt*length(major.x))*trapz(dt*(0:1:numel(major.x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                Tmajor = dt*length(major.x);
                if ~isempty(minor(1).y)
                    Pminor = nan(1,length(minor));
                    Tminor = nan(1,length(minor));
                    for mm=1:length(minor)
                        dB = diff([minor(mm).y minor(mm).y(1)]);
                        DB = max(minor(mm).y)-min(minor(mm).y);
                        Pminor(mm) = 1/(dt*length(minor(mm).x))*trapz(dt*(0:1:numel(minor(mm).x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                        Tminor(mm) = dt*length(minor(mm).x);
                    end
                else
                    Pminor = 0;
                    Tminor = 0;
                end
                psh(index(ii)) = sum([Pmajor Pminor].*[Tmajor Tminor])/sum([Tmajor Tminor]);
            end
        end

        % rotor
        index = 1:1:EleNo;
        index = index(groNo==22);
        nR = numel(index);

        prh = zeros(1,EleNo);

        for ii=1:nR
            x = th;
            y = br(:,index(ii));
            [major,minor] = minorLoopDetection(x,y);
            if ~isempty(major.x)
                dt = (major.x(2)-major.x(1))*(pi/180)/(2*pi*f0);
                dB = diff([major.y major.y(1)]);
                DB = max(major.y)-min(major.y);
                Pmajor = 1/(dt*length(major.x))*trapz(dt*(0:1:numel(major.x)-1),kiR*abs(dB/dt).^alphaR*DB^(betaR-alphaR));
                Tmajor = dt*length(major.x);
                if ~isempty(minor(1).y)
                    Pminor = nan(1,length(minor));
                    Tminor = nan(1,length(minor));
                    for mm=1:length(minor)
                        dB = diff([minor(mm).y minor(mm).y(1)]);
                        DB = max(minor(mm).y)-min(minor(mm).y);
                        Pminor(mm) = 1/(dt*length(minor(mm).x))*trapz(dt*(0:1:numel(minor(mm).x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                        Tminor(mm) = dt*length(minor(mm).x);
                    end
                else
                    Pminor = 0;
                    Tminor = 0;
                end
                prh(index(ii)) = sum([Pmajor Pminor].*[Tmajor Tminor])/sum([Tmajor Tminor]);
            end
        end

        psh = psh.*vol.*rhoS;
        prh = prh.*vol.*rhoR;


        % eddy current (harmonic composition)
        bsx = abs(fft(real(bs)))*2/nsim;
        bsy = abs(fft(imag(bs)))*2/nsim;
        brx = abs(fft(real(br)))*2/nsim;
        bry = abs(fft(imag(br)))*2/nsim;
        Bs = (bsx.^2+bsy.^2).^0.5;
        Br = (brx.^2+bry.^2).^0.5;

        freq = (0:1:(nsim-1))*f0;   % frequency associated to each harmonic component
        freq(floor(nsim/2):end) = 0;
        freq = repmat(freq.',1,EleNo);
        vol  = repmat(vol,nsim,1);

        %     index = 1:1:EleNo;
        %     index = index(groNo==12);
        index = 1:1:EleNo;
        psc = (mat.Stator.ke*freq(:,index).^2.*Bs(:,index).^2)*mat.Stator.kgm3.*vol(:,index);
        psc = sum(psc,1);

        %     index = 1:1:EleNo;
        %     index = index(groNo==22);
        index = 1:1:EleNo;
        prc = (mat.Rotor.ke*freq(:,index).^2.*Br(:,index).^2)*mat.Rotor.kgm3.*vol(:,index);
        prc = sum(prc,1);

        % permanent magnet (as FEMM)
        Jm = (fft(am))*2/nsim;
        Jo = zeros(size(Jm));
        for mm=1:numPM
            index = 1:1:EleNo;
            index = index(groNo==(200+mm));
            vmag = sum(vol(index));
            tmp = (Jm(:,index).*vol(index))/vmag;
            tmp = sum(tmp,2);
            Jo(:,index) = repmat(tmp,1,length(index));
            Jm(:,index) = Jm(:,index)-Jo(:,index);
        end

        index = 1:1:EleNo;
        ppm = 1/2*mat.LayerMag.sigmaPM*(2*pi*freq(:,index)).^2.*abs(Jm(:,index)).^2.*vol(:,index);

    case 2 % iGSE+eddy current, both axis, no correction factor
        % hysteresis with minor loops

        th1 = linspace(0,2*pi,1001);

        kiS = khS/((2*pi)^(alphaS-1)*trapz(th1,abs(cos(th1)).^alphaS*2^(betaS-alphaS)));
        kiR = khR/((2*pi)^(alphaR-1)*trapz(th1,abs(cos(th1)).^alphaR*2^(betaR-alphaR)));

        % stator
        index = 1:1:EleNo;
        index = index(groNo==12);
        nS = numel(index);

        psh = zeros(1,EleNo);

        sigma = -100*0; %% for kmech (mech stress impact on stator) from Yamazaki

        for ii=1:nS
            for dd=1:2
                x = th;
                y = bs(:,index(ii));
                % detection of the main/quadrature flux density axis
                ang = angle(y);
                ang(ang<0) = ang(ang<0)+pi;
                ang = mean(ang);
                y = y*exp(-j*ang);
                if dd==1
                    y = real(y);
                else
                    y = imag(y);
                end
                [major,minor] = minorLoopDetection(x,y);
             
                if ~isempty(major.y)
                    dt = (major.x(2)-major.x(1))*(pi/180)/(2*pi*f0);
                    dB = diff([major.y major.y(1)]);
                    DB = max(major.y)-min(major.y);
                    Pmajor = 1/(dt*length(major.x))*trapz(dt*(0:1:numel(major.x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                    Tmajor = dt*length(major.x);
                    if ~isempty(minor(1).y)
                        Pminor = nan(1,length(minor));
                        Tminor = nan(1,length(minor));
                        for mm=1:length(minor)
                            dB = diff([minor(mm).y minor(mm).y(1)]);
                            DB = max(minor(mm).y)-min(minor(mm).y);
                            Pminor(mm) = 1/(dt*length(minor(mm).x))*trapz(dt*(0:1:numel(minor(mm).x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                            kdc_min   = 0.65*mean(abs(minor(mm).y))^2.1+1;
                            kmech_min = 1 + (4.9-1)*exp(-mean(abs(minor(mm).y))/0.7)*(1-exp(sigma/100));
                            Pminor(mm) = Pminor(mm)*kdc_min*kmech_min;
                            Tminor(mm) = dt*length(minor(mm).x);
                        end
                    else
                        Pminor = 0;
                        Tminor = 0;
                    end
                    %%dc bias
                    kdc_maj = 0.65*mean(abs(bs(:,index(ii))))^2.1+1;
                    kmech_maj = 1 + (4.9-1)*exp(-mean(abs(bs(:,index(ii))))/0.7)*(1-exp(sigma/100));
                    psh(index(ii)) = psh(index(ii))+sum([Pmajor*kdc_maj*kmech_maj Pminor].*[Tmajor Tminor])/sum([Tmajor Tminor]);
                end
            end
        end


        % rotor
        index = 1:1:EleNo;
        index = index(groNo==22);
        nR = numel(index);

        prh = zeros(1,EleNo);

        for ii=1:nR
            for dd=1:2
                x = th;
                y = br(:,index(ii));
                % detection of the main/quadrature flux density axis
                ang = angle(y);
                ang(ang<0) = ang(ang<0)+pi;
                ang = mean(ang);
                y = y*exp(-j*ang);
                if dd==1
                    y = real(y);
                else
                    y = imag(y);
                end
                [major,minor] = minorLoopDetection(x,y);
                if ~isempty(major.y)
                    dt = (major.x(2)-major.x(1))*(pi/180)/(2*pi*f0);
                    dB = diff([major.y major.y(1)]);
                    DB = max(major.y)-min(major.y);
                    Pmajor = 1/(dt*length(major.x))*trapz(dt*(0:1:numel(major.x)-1),kiR*abs(dB/dt).^alphaR*DB^(betaR-alphaR));
                    Tmajor = dt*length(major.x);
                    if ~isempty(minor(1).y)
                        Pminor = nan(1,length(minor));
                        Tminor = nan(1,length(minor));
                        for mm=1:length(minor)
                            dB = diff([minor(mm).y minor(mm).y(1)]);
                            DB = max(minor(mm).y)-min(minor(mm).y);
                            Pminor(mm) = 1/(dt*length(minor(mm).x))*trapz(dt*(0:1:numel(minor(mm).x)-1),kiR*abs(dB/dt).^alphaR*DB^(betaR-alphaR));
                            kdc_min(mm) = 0.65*mean(abs(minor(mm).y))^2.1+1;
                            Pminor(mm) = Pminor(mm)*kdc_min(mm);
                            Tminor(mm) = dt*length(minor(mm).x);
                        end
                    else
                        Pminor = 0;
                        Tminor = 0;
                    end
                    kdc_maj = 0.65*mean(abs(bs(:,index(ii))))^2.1+1;
                    prh(index(ii)) = prh(index(ii))+sum([Pmajor*kdc_maj Pminor].*[Tmajor Tminor])/sum([Tmajor Tminor]);
                end
            end
        end

        psh = psh.*vol.*rhoS;
        prh = prh.*vol.*rhoR;


        % eddy current (harmonic decomposition)
        bsx = abs(fft(real(bs)))*2/nsim;
        bsy = abs(fft(imag(bs)))*2/nsim;
        brx = abs(fft(real(br)))*2/nsim;
        bry = abs(fft(imag(br)))*2/nsim;
        Bs = (bsx.^2+bsy.^2).^0.5;
        Br = (brx.^2+bry.^2).^0.5;

        freq = (0:1:(nsim-1))*f0;   % frequency associated to each harmonic component
        freq(floor(nsim/2):end) = 0;
        freq = repmat(freq.',1,EleNo);
        vol  = repmat(vol,nsim,1);

        index = 1:1:EleNo;
        psc = (mat.Stator.ke*freq(:,index).^2.*Bs(:,index).^2)*mat.Stator.kgm3.*vol(:,index);
        psc = sum(psc,1);

        index = 1:1:EleNo;
        prc = (mat.Rotor.ke*freq(:,index).^2.*Br(:,index).^2)*mat.Rotor.kgm3.*vol(:,index);
        prc = sum(prc,1);

        % permanent magnet (as FEMM)
        Jm = fft(am)*2/nsim;
        Jo = zeros(size(Jm));

        Crf = ones(nsim,EleNo);
        Cef = ones(nsim,EleNo);
        ppm_PM = zeros(numPM,1);
        Crf_array = ones(nsim,numPM);
        Cef_array = ones(nsim,numPM);

        for mm=1:numPM
            index = 1:1:EleNo;
            index = index(groNo==(200+mm));
            vmag = sum(vol(index));
            tmp = (Jm(:,index).*vol(index))/vmag;
            tmp = sum(tmp,2);
            Jo(:,index) = repmat(tmp,1,length(index));
            Jm(:,index) = Jm(:,index)-Jo(:,index);

            if strcmp(geo.RotType,'Seg')
                %%%CORRECTION FACTOR - REACTION FIELD  
                %M. Hullmann and B. Ponick, "General Analytical Description of the Effects of Segmentation on Eddy Current Losses in Rectangular Magnets," 
                %2022 International Conference on Electrical Machines (ICEM), Valencia, Spain, 2022, pp. 1757-1762, doi: 10.1109/ICEM51905.2022.9910629.
                delta = 1000*sqrt(1./(pi*mu*mat.LayerMag.sigmaPM*freq(:,1)).*(1+repmat(PMclear(mm)/hc(mm),[nsim 1])));
                PMratio = (PMdim(mm)/PMNc(mm))./delta;
                Crf_array(:,mm) = 6./(PMratio).^3.*(sinh(PMratio)-sin(PMratio))./(cosh(PMratio)+cos(PMratio));
                Crf_array(isnan(Crf_array)) = 1;
                Crf(:,index) = repmat(Crf_array(:,mm),1,length(index));

                %%%CORRECTION FACTOR - END SIDE EFFECT
                n = 1:1:10^3;
                n = [10e-3 n];
                lambdaN = (2*n+1)*pi/(PMdim(mm)/PMNc(mm));
                betaN = sqrt(lambdaN.^2+2*1i./delta.^2);
                Cef_sigma = ((lambdaN.^2-2*imag(betaN).^2).*real(betaN).*lambdaN.^3.*sinh(real(betaN)*l)+(lambdaN.^2+2*real(betaN).^2).*imag(betaN).*lambdaN.^3.*sin(imag(betaN)*l))./((2*n+1).^5.*abs(betaN).^6.*(cosh(real(betaN)*l)+cos(imag(betaN)*l)));
                Cef_sigma(isnan(Cef_sigma)) = 0;
                Cef_array(:,mm) = 1-32*(PMdim(mm)/PMNc(mm))/pi^5/l.*PMratio.^3.*(cosh(PMratio)+cos(PMratio))./(sinh(PMratio)-sin(PMratio)).*sum(Cef_sigma')';
                Cef_array(isnan(Cef_array)) = 1;
                Cef(:,index) = repmat(Cef_array(:,mm),1,length(index));

                tmp = 1/2*mat.LayerMag.sigmaPM*(2*pi*freq).^2.*(abs(Jm)).^2.*vol.*Crf.*Cef;
                ppm_PM(mm) = sum(sum(tmp(:,index)));
            end
        end

        ppm_noRFno3D = 1/2*mat.LayerMag.sigmaPM*(2*pi*freq).^2.*(abs(Jm)).^2.*vol;
        ppm_no3D = 1/2*mat.LayerMag.sigmaPM*(2*pi*freq).^2.*(abs(Jm)).^2.*vol.*Crf;
        ppm = 1/2*mat.LayerMag.sigmaPM*(2*pi*freq).^2.*(abs(Jm)).^2.*vol.*Crf.*Cef;
        
%         tmp = sum(ppm_noRFno3D,2);
%         ppmold_plot = tmp((1:(nsim/2)-1));
%         
%         tmp = sum(ppm,2);
%         ppm_plot = tmp((1:(nsim/2)-1));
% 
%         figure
%         figSetting(14,6,10)
%         yyaxis left
%         ylabel('$P_{PM,stat}$ [W]')
%         xlabel('$f$ [kHz]')
%         %bar(freq(1:(nsim/2)-1,5990)/1000,abs(Jm(1:(nsim/2)-1,5980)/max(Jm(1:(nsim/2)-1,5980))),0.5)
%         bar(freq(1:(nsim/2)-1,1)/1000,ppmold_plot*(2*geo.p/geo.ps),'HandleVisibility','off')
% %         bar(freq(1:(nsim/2)-1,1)/1000,ppm_plot*(2*geo.p/geo.ps),'HandleVisibility','off')
%         yyaxis right
%         ylabel('$k$')
%         plot(freq(1:(nsim/2)-1,1)/1000,Crf_array(1:(nsim/2)-1,1),'DisplayName','$k_{rf}$')
%         plot(freq(2:(nsim/2)-1,1)/1000,Cef_array(2:(nsim/2)-1,1),'DisplayName','$k_{ef}$')
%  xlim([0 35])
% 
%         figure
%         figSetting
%         plot(PMratio(1:(nsim/2)-1),Crf_array(1:(nsim/2)-1,1))
%         plot(PMratio(1:(nsim/2)-1),Cef_array(1:(nsim/2)-1,1))
%         xlim([0 6])
%         yticks([0:0.1:1])
% 
%                 figure
%         figSetting
%          plot(freq(2:(nsim/2)-1,1)/1000,Cef_array(2:(nsim/2)-1,1))
%         plot(freq(2:(nsim/2)-1,1)/1000,Crf_array(2:(nsim/2)-1,1))
%         xlim([0 100])
%         ylabel('$k_{ef}$')
%         xlabel('$f$ [kHz]')
%         ylim([0.8 1.2])
%         yticks([0:0.1:2])

    case 3 % iGSE+eddy current, both axis, correction factor on hysteresis loss from 2004-Bottauscio-Additional Loss
        % hysteresis with minor loops

        Rh.x = [0.0000 0.0929 0.1908 0.2887 0.3866 0.4845 0.5824 0.6803 0.7782 0.8760 0.9739 1.0718 1.1697 1.2676 1.3655 1.4634 1.5613 1.6592 1.7571 1.8550 2.0000];
        Rh.y = [2.3750 2.2325 2.0900 1.9475 1.8525 1.7575 1.6625 1.5675 1.5200 1.5010 1.4725 1.4630 1.4250 1.3775 1.2825 1.2065 1.0735 0.9111 0.7144 0.4807 0.0000];

        th1 = linspace(0,2*pi,1001);

        kiS = khS/((2*pi)^(alphaS-1)*trapz(th1,abs(cos(th1)).^alphaS*2^(betaS-alphaS)));
        kiR = khR/((2*pi)^(alphaR-1)*trapz(th1,abs(cos(th1)).^alphaR*2^(betaR-alphaR)));

        % stator
        index = 1:1:EleNo;
        index = index(groNo==12);
        nS = numel(index);

        psh = zeros(1,EleNo);

        for ii=1:nS
            for dd=1:2
                x = th;
                y = bs(:,index(ii));
                % detection of the main/quadrature flux density axis
                ang = angle(y);
                ang(ang<0) = ang(ang<0)+pi;
                ang = mean(ang);
                y = y*exp(-j*ang);
                if dd==1
                    y = real(y);
                    Bmax = max(abs(y));
                    kRot = 1;
                else
                    y = imag(y);
                    R = interp1(Rh.x,Rh.y,Bmax,'linear',0);
                    kRot = R-1;
                end
                [major,minor] = minorLoopDetection(x,y);
                if ~isempty(major.y)
                    dt = (major.x(2)-major.x(1))*(pi/180)/(2*pi*f0);
                    dB = diff([major.y major.y(1)]);
                    DB = max(major.y)-min(major.y);
                    Pmajor = 1/(dt*length(major.x))*trapz(dt*(0:1:numel(major.x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                    Tmajor = dt*length(major.x);
                    if ~isempty(minor(1).y)
                        Pminor = nan(1,length(minor));
                        Tminor = nan(1,length(minor));
                        for mm=1:length(minor)
                            dB = diff([minor(mm).y minor(mm).y(1)]);
                            DB = max(minor(mm).y)-min(minor(mm).y);
                            Pminor(mm) = 1/(dt*length(minor(mm).x))*trapz(dt*(0:1:numel(minor(mm).x)-1),kiS*abs(dB/dt).^alphaS*DB^(betaS-alphaS));
                            Tminor(mm) = dt*length(minor(mm).x);
                        end
                    else
                        Pminor = 0;
                        Tminor = 0;
                    end
                    psh(index(ii)) = psh(index(ii))+sum([Pmajor Pminor].*[Tmajor Tminor])/sum([Tmajor Tminor])*kRot;
                end
            end
        end

        % rotor
        index = 1:1:EleNo;
        index = index(groNo==22);
        nR = numel(index);

        prh = zeros(1,EleNo);

        for ii=1:nR
            for dd=1:2
                x = th;
                y = br(:,index(ii));
                % detection of the main/quadrature flux density axis
                ang = angle(y);
                ang(ang<0) = ang(ang<0)+pi;
                ang = mean(ang);
                y = y*exp(-j*ang);
                if dd==1
                    y = real(y);
                    Bmax = max(abs(y));
                    kRot = 1;
                else
                    y = imag(y);
                    R = interp1(Rh.x,Rh.y,Bmax,'linear',0);
                    kRot = R-1;
                end
                [major,minor] = minorLoopDetection(x,y);
                if ~isempty(major.y)
                    dt = (major.x(2)-major.x(1))*(pi/180)/(2*pi*f0);
                    dB = diff([major.y major.y(1)]);
                    DB = max(major.y)-min(major.y);
                    Pmajor = 1/(dt*length(major.x))*trapz(dt*(0:1:numel(major.x)-1),kiR*abs(dB/dt).^alphaR*DB^(betaR-alphaR));
                    Tmajor = dt*length(major.x);
                    if ~isempty(minor(1).y)
                        Pminor = nan(1,length(minor));
                        Tminor = nan(1,length(minor));
                        for mm=1:length(minor)
                            dB = diff([minor(mm).y minor(mm).y(1)]);
                            DB = max(minor(mm).y)-min(minor(mm).y);
                            Pminor(mm) = 1/(dt*length(minor(mm).x))*trapz(dt*(0:1:numel(minor(mm).x)-1),kiR*abs(dB/dt).^alphaR*DB^(betaR-alphaR));
                            Tminor(mm) = dt*length(minor(mm).x);
                        end
                    else
                        Pminor = 0;
                        Tminor = 0;
                    end
                    prh(index(ii)) = prh(index(ii))+sum([Pmajor Pminor].*[Tmajor Tminor])/sum([Tmajor Tminor])*kRot;
                end
            end
        end

        psh = psh.*vol.*rhoS;
        prh = prh.*vol.*rhoR;


        % eddy current (harmonic decomposition)
        bsx = abs(fft(real(bs)))*2/nsim;
        bsy = abs(fft(imag(bs)))*2/nsim;
        brx = abs(fft(real(br)))*2/nsim;
        bry = abs(fft(imag(br)))*2/nsim;
        Bs = (bsx.^2+bsy.^2).^0.5;
        Br = (brx.^2+bry.^2).^0.5;

        freq = (0:1:(nsim-1))*f0;   % frequency associated to each harmonic component
        freq(floor(nsim/2):end) = 0;
        freq = repmat(freq.',1,EleNo);
        vol  = repmat(vol,nsim,1);

        %     index = 1:1:EleNo;
        %     index = index(groNo==12);
        index = 1:1:EleNo;
        psc = (mat.Stator.ke*freq(:,index).^2.*Bs(:,index).^2)*mat.Stator.kgm3.*vol(:,index);
        psc = sum(psc,1);

        %     index = 1:1:EleNo;
        %     index = index(groNo==22);
        index = 1:1:EleNo;
        prc = (mat.Rotor.ke*freq(:,index).^2.*Br(:,index).^2)*mat.Rotor.kgm3.*vol(:,index);
        prc = sum(prc,1);

        % permanent magnet (as FEMM)
        Jm = (fft(am))*2/nsim;
        Jo = zeros(size(Jm));
        for mm=1:numPM
            index = 1:1:EleNo;
            index = index(groNo==(200+mm));
            vmag = sum(vol(index));
            tmp = (Jm(:,index).*vol(index))/vmag;
            tmp = sum(tmp,2);
            Jo(:,index) = repmat(tmp,1,length(index));
            Jm(:,index) = Jm(:,index)-Jo(:,index);
        end

        index = 1:1:EleNo;
        ppm = 1/2*mat.LayerMag.sigmaPM*(2*pi*freq(:,index)).^2.*abs(Jm(:,index)).^2.*vol(:,index);
end


% Output data
SOL.psh   = psh;
SOL.psc   = psc;
SOL.prh   = prh;
SOL.prc   = prc;
SOL.ppm           = ppm;
SOL.ppm_no3D      = ppm_no3D;  
SOL.ppm_noRFno3D  = ppm_noRFno3D;
SOL.ppm_PM        = ppm_PM;
SOL.freq  = freq;
% SOL.bs    = bs;
% SOL.br    = br;
% SOL.am    = am;
SOL.Jm    = Jm;
