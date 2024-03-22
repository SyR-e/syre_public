function ComputeDecouplingMatrix(motorModel,n_set)

UserMacrosH_path = [motorModel.data.pathname motorModel.data.motorName '_ctrl_INST\User_functions\Inc\User_Macros.h'];


n = n_set;  % number of sets

TD = zeros(n,n);

TD(1,:) = 1;

x=1;
row = 2;
for k=1:(n-1)    
    
    wk = sqrt((n*n-n*k)/(n-k+1));
    qk = -sqrt(n/((n-k)*(n-k+1)));
    
    TD(row,x) = wk;
    TD(row,x+1:end) = qk;
    x=x+1;
    row = row+1;

end


TD = (1/n)*TD;
inv_TD = inv(TD);

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



% Direct Decoupling Algorithm

User_Macros = [ User_Macros; strings(2,1)];
User_Macros = [User_Macros;sprintf("//Decoupling algorithm dq1 dq2..dqn -> dq_cm dq_dm1 dq_dm2...dq_dm(n-1)")];
User_Macros = [ User_Macros; strings(1,1)];

% Define Macro

tmp_s = sprintf('#define _Decoupling(');

for i=1:n
tmp_s = [tmp_s sprintf('dq%d,',i)];
end

tmp_s = [tmp_s sprintf('dq_cm')];

for i=1:(n-1)
tmp_s = [tmp_s sprintf(',dq_dm%d',i)];
end

User_Macros = [User_Macros; sprintf("%s); \\",tmp_s)];


% Common mode subspace

tmp_d = sprintf('dq_cm.d =');
tmp_q = sprintf('dq_cm.q =');


for i=1:n
    tmp_d = [tmp_d sprintf('+(%.4f)*dq%d.d',TD(1,i),i)];
    tmp_q = [tmp_q sprintf('+(%.4f)*dq%d.q',TD(1,i),i)];
end
User_Macros = [User_Macros; sprintf(blanks(4)+"%s; \\",tmp_d)];
User_Macros = [User_Macros; sprintf(blanks(4)+"%s; \\",tmp_q)];

% Differential mode subspace

for i=1:(n-1)
        
    tmp_d = sprintf('dq_dm%d.d =',i);
    tmp_q = sprintf('dq_dm%d.q =',i);

    for j=1:n
        tmp_d = [tmp_d sprintf('+(%.4f)*dq%d.d',TD(i+1,j),j)];
        tmp_q = [tmp_q sprintf('+(%.4f)*dq%d.q',TD(i+1,j),j)];
    end
    User_Macros = [User_Macros; sprintf(blanks(4)+"%s; \\",tmp_d)];
    User_Macros = [User_Macros; sprintf(blanks(4)+"%s; \\",tmp_q)];

end


% Inverse Decoupling Algorithm

User_Macros = [ User_Macros; strings(2,1)];
User_Macros = [User_Macros;sprintf("//Decoupling algorithm dq_cm dq_dm1 dq_dm2...dq_dm(n-1)-> dq1 dq2..dqn  ")];
User_Macros = [ User_Macros; strings(1,1)];

% Define Macro

tmp_s = sprintf('#define _InvDecoupling(');
tmp_s = [tmp_s sprintf('dq_cm')];

for i=1:(n-1)
    tmp_s = [tmp_s sprintf(',dq_dm%d',i)];
end

for i=1:n
    tmp_s = [tmp_s sprintf(',dq%d',i)];
end

User_Macros = [User_Macros; sprintf("%s); \\",tmp_s)];


for i=1:n
    tmp_d = sprintf('dq%d.d=(%.4f)*dq_cm.d',i,inv_TD(i,1));
    tmp_q = sprintf('dq%d.q=(%.4f)*dq_cm.q',i,inv_TD(i,1));
    for j=1:(n-1)
        tmp_d = [tmp_d sprintf('+(%.4f)*dq_dm%d.d',inv_TD(i,j+1),j)]; 
        tmp_q = [tmp_q sprintf('+(%.4f)*dq_dm%d.q',inv_TD(i,j+1),j)];
    end
    User_Macros = [User_Macros; sprintf(blanks(4)+"%s; \\",tmp_d)];
    User_Macros = [User_Macros; sprintf(blanks(4)+"%s; \\",tmp_q)];


end    

fid = fopen(UserMacrosH_path, 'w');
for i = 1:numel(User_Macros)
    fprintf(fid,'%s\n', User_Macros{i});
end
fclose(fid);








end