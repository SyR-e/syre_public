import sys
sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64")
sys.path.append(r"C:/Program Files/AnsysEM/AnsysEM20.1/Win64/PythonFiles/DesktopPlugin")
import ScriptEnv
ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")
oDesktop.RestoreWindow()
oDesktop.CloseProject("VGB_HS")
