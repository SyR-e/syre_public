% stampa un vettore double m x n in formato tabella per il C

function StampaVarg(f1,Xvar,m,n,var,commento,formato)

% m numero righe
% n numero colonne
fprintf(f1,'\r\n');
fprintf(f1,commento);
fprintf(f1,'\r\n');
fprintf(f1,['float  ']);
if(m>1)
    fprintf(f1,[var '[' num2str(m) ']' '[' num2str(n) ']' ' = {\n']);
else
    fprintf(f1,[var '[' num2str(n) ']' ' = \n']);
end
for indexFound = 1:1:m*n,
    Xprint=(Xvar(indexFound));
    if(mod(indexFound-1,n)==0)
        fprintf(f1,['{' formato],Xprint);
    else
        fprintf(f1,[',' formato],Xprint);
    end
    if((mod(indexFound,n)==0))
        fprintf(f1,'},\n');
    end
end
if(m>1)
    fprintf(f1,['};\r\n']);
else
    fprintf(f1,[';\r\n']);
end

