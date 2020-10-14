% Copyright 2014
%
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
%
%        http://www.apache.org/licenses/LICENSE-2.0
%
%    Unless required by applicable law or agreed to in writing, dx
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

function out=spettro_v2(ValoreNelTempo,n)

dime = max(size(ValoreNelTempo));
% non si possono richiedere più armoniche del numero di campioni-1
if (n > dime)
    n = dime-1;
end

ValoreNelTempo=reshape(ValoreNelTempo,dime,1);

a=fft(ValoreNelTempo);
Continua=a(1)/dime;
Armoniche=2*abs(a(2:dime))/dime;  % la procedura seguente serve per eliminare il valore medio plottando solo le armoniche
% /dme restituisce i valori pari al
% valore medio Ma sono in percentuale
% rispetto alla coppia media?????????????
%Armoniche=abs(a(1:dime))/dime;
%keyboard
numarm = 1:n;
% keyboard
figure
%bar([0 numarm],[Continua;Armoniche(numarm)],0.5,'r')
bar(numarm,Armoniche(numarm),0.5,'r')
xlabel('harmonic order'),grid,
set(gca,'xTickLabel',numarm),
set(gca,'xTick',1:1:length(numarm))

%     figure(fig+1)
%     plot(ValoreNelTempo,'b-');
%     grid on

% out=[Continua;Armoniche(numarm)];
out=Armoniche(numarm);