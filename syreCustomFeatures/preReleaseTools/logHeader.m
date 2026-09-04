% Copyright 2023
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

function logHeader(logfile,logTitle)

fid = fopen(logfile,'w');
fprintf(fid,[logTitle '\r\n']);
fprintf(fid,['Date           : ' char(datetime('now')) '\r\n']);
fprintf(fid,['Computer name  : ' getenv('computername') '\r\n']);
fprintf(fid,['CPU            : ' feature('getcpu') '\r\n']);
fprintf(fid,['# of cores     : ' int2str(feature('numcores')) '\r\n']);
if ispc
    [~,tmp] = memory;
    fprintf(fid,['RAM            : ' int2str(tmp.PhysicalMemory.Total/1024/1024/1024) ' GB\r\n' ]);
end
fprintf(fid,['OS             : ' feature('getos') '\r\n']);
tmp = ver('matlab');
fprintf(fid,['Matlab release : ' tmp.Release(2:end-1) '\r\n' ]);

fprintf(fid,'\r\n\r\n');

fclose(fid);
edit(logfile);