function PrintClarke(motorModel,n_set)


UserMacrosH_path = [motorModel.data.pathname motorModel.data.motorName '_ctrl_INST\User_functions\Inc\User_Macros.h'];


l = n_set;
k=3;
n = k*l;

%% Define theta_n MATRIX

theta_n = zeros(l,k);

for j=1:l
    for i=1:k
        theta_n(j,i) = (180/n)*(2*l*(i-1)+j-1);
    end
end

% Compute Clake transformation

for j=1:l

    tmp_cos = (2/3)*[cosd(theta_n(j,1)) cosd(theta_n(j,2)) cosd(theta_n(j,3))];
    tmp_sin = (2/3)*[sind(theta_n(j,1)) sind(theta_n(j,2)) sind(theta_n(j,3))];
    eval(['clarke_' num2str(j) '=[ tmp_cos; tmp_sin];']);

end

for j=1:l

    tmp_cos = [cosd(theta_n(j,1)) ; cosd(theta_n(j,2)); cosd(theta_n(j,3))];
    tmp_sin = [sind(theta_n(j,1)) ; sind(theta_n(j,2)) ; sind(theta_n(j,3))];
    eval(['clarke_inv_' num2str(j) '=[ tmp_cos tmp_sin];']);

end


fid = fopen(UserMacrosH_path,'r');
i = 1;
tline = fgetl(fid);
readData{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    readData{i} = tline;
end
fclose(fid);

readData(end) = [];
User_Macros = string(readData)';

User_Macros = [ User_Macros; strings(2,1)];
User_Macros = [User_Macros;sprintf("//Direct Clarke transformation (a,b,c)--> (alpha,beta)")];
User_Macros = [ User_Macros; strings(2,1)];

for i=1:l
    clarke_tmp = eval(['clarke_' num2str(i)]);
    User_Macros = [User_Macros;sprintf("#define _clarke%d(abc,ab); \\",i)];
    User_Macros = [ User_Macros;sprintf(blanks(4)+"ab.alpha = %.4f*abc.a+(%.4f)*abc.b+(%.4f)*abc.c;  \\",clarke_tmp(1,1),clarke_tmp(1,2),clarke_tmp(1,3))];
    User_Macros = [ User_Macros;sprintf(blanks(4)+"ab.beta = %.4f*abc.a+(%.4f)*abc.b+(%.4f)*abc.c;  \\",clarke_tmp(2,1),clarke_tmp(2,2),clarke_tmp(2,3))];
    User_Macros = [ User_Macros; strings(1,1)];
end

User_Macros = [User_Macros;sprintf("//Inverse Clarke transformation (alpha,beta)--> (a,b,c)")];
User_Macros = [ User_Macros; strings(2,1)];

for i=1:l
    clarke_inv_tmp = eval(['clarke_inv_' num2str(i)]);
    User_Macros = [User_Macros;sprintf("#define _invclarke%d(ab,abc);  \\",i)];
    User_Macros = [ User_Macros;sprintf(blanks(4)+"abc.a = %.4f*ab.alpha+(%.4f)*ab.beta;  \\",clarke_inv_tmp(1,1),clarke_inv_tmp(1,2))];
    User_Macros = [ User_Macros;sprintf(blanks(4)+"abc.b = %.4f*ab.alpha+(%.4f)*ab.beta;  \\",clarke_inv_tmp(2,1),clarke_inv_tmp(2,2))];
    User_Macros = [ User_Macros;sprintf(blanks(4)+"abc.c = %.4f*ab.alpha+(%.4f)*ab.beta;  \\",clarke_inv_tmp(3,1),clarke_inv_tmp(3,2))];
    User_Macros = [ User_Macros; strings(1,1)];
end


%% -------------------------Stampa del nuovo file------------------------%%
fid = fopen(UserMacrosH_path, 'w');
for i = 1:numel(User_Macros)
    fprintf(fid,'%s\n', User_Macros{i});
end
fclose(fid);










end