%#create the FEM Coil components
function FEMCoil_Creation_JMAG (circuit, phaseLabel, position)
        circuit.CreateComponent('Coil', strcat(phaseLabel, '-Phase Coil'));
        circuit.CreateInstance(strcat(phaseLabel, '-Phase Coil'), position(1), position(2));
        circuit.GetComponent(strcat(phaseLabel, '-Phase Coil')).SetValue('Turn', 'coilTurn');
        circuit.GetComponent(strcat(phaseLabel, '-Phase Coil')).SetValue('Resistance', 'coilRes');
        circuit.GetComponent(strcat(phaseLabel, '-Phase Coil')).SetValue('LeakageInductance', 'coilLind');
end

