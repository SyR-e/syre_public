% Copyright 2014
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

function objects(fid)
%% Generazione sezione OBJECTS
% Scrittura intestazione OBJECTS
testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'SECTION';
fprintf(fid,'%s\n',testo);
testo = '2';
fprintf(fid,'%s\n',testo);
testo = 'OBJECTS';
fprintf(fid,'%s\n',testo);

%scrivere testo sezione

testo = '0';
fprintf(fid,'%s\n',testo);
testo = 'ENDSEC';
fprintf(fid,'%s\n',testo);

