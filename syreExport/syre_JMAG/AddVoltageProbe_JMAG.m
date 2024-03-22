function AddVoltageProbe (circuit, index, position, phaseLabel)
        circuit.CreateComponent('VoltageProbe', strcat('VP',num2str(index)));
        circuit.CreateInstance(strcat('VP',num2str(index)), position(1), position(2));
        circuit.GetComponent(strcat('VP',num2str(index))).SetName(phaseLabel);
end
