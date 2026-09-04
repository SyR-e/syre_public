% Copyright 2026
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

function [a2, b2, c2] = retta_ruotata_attorno_punto(a, b, c, x0, y0, theta)
%Ruota una retta 2D implicita attorno a un punto.
%
%   [a2, b2, c2] = ruotaRettaAttornoPunto(a, b, c, x0, y0, theta)
%
%   Ruota la retta:
%
%       a*x + b*y + c = 0
%
%   di un angolo theta attorno al punto P = (x0, y0) non necessariamente
%   appartenente alla retta.
%
%   Input:
%       a, b, c : coefficienti della retta iniziale
%       x0, y0  : coordinate del punto di rotazione
%       theta   : angolo di rotazione in gradi
%
%   Output:
%       a2, b2, c2 : coefficienti della retta ruotata
%
%   Convenzione:
%       theta > 0 indica una rotazione antioraria.

    % Ruota il vettore normale della retta
    a2 = a*cosd(theta) - b*sind(theta);
    b2 = a*sind(theta) + b*cosd(theta);

    % Aggiorna il termine noto per ruotare la retta attorno a (x0, y0)
    c2 = c + (a - a2)*x0 + (b - b2)*y0;
end