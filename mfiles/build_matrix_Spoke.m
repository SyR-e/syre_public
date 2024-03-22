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

function rotore = build_matrix_Spoke(temp,geo)


% Thesis Automotive Engineering
% Matrix created by Franco Porras

%% Geometry for line

xxB1k  = temp.xxB1k;    
yyB1k  = temp.yyB1k;    
xxB2k  = temp.xxB2k;    
yyB2k  = temp.yyB2k;        
xxD1k  = temp.xxD1k;    
yyD1k  = temp.yyD1k;    
xxD2k  = temp.xxD2k;    
yyD2k  = temp.yyD2k;    

rotore = [];
Mag    = [];

materialCodes;
% codMatFeRot   = 5 Rotor
% codMatBar     = 6   PM
indexEle = 1;

rotore = [rotore
        xxB1k yyB1k xxD1k yyD1k NaN NaN 0 codMatFeRot indexEle
        xxD1k yyD1k xxD2k yyD2k NaN NaN 0 codMatFeRot indexEle
        xxD2k yyD2k xxB2k yyB2k NaN NaN 0 codMatFeRot indexEle
        ];

indexEle = 2;

Mag = [Mag
        xxB1k yyB1k xxD1k yyD1k NaN NaN 0 codMatBar indexEle
        xxD1k yyD1k xxD2k yyD2k NaN NaN 0 codMatBar indexEle
        xxD2k yyD2k xxB2k yyB2k NaN NaN 0 codMatBar indexEle
        ];

rotore = [rotore;Mag];

