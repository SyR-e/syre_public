defipypath='C:\Program Files\AnsysEM\AnsysEM20.1\Win64\common\IronPython\ipy64.exe';
[ipy64exe, ipypath] = uigetfile(defipypath,'Select ipy64.exe (Ansys) Directory');
ipypath=strcat('"',ipypath,ipy64exe,'" "');
close_project(ipypath);

function close_project(ipypath)
command=strcat(ipypath,pwd,'\close_projectAM.py"');
dos(command);
end
