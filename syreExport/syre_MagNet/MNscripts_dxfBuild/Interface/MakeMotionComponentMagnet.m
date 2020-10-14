function magDataUpt = MakeMotionComponentMagnet(magData,compMot,compMotName,per_mec,pos_mec)



mh = magData.magnetHandler;

permec = num2str(per_mec);
posmec = num2str(pos_mec);
nComp = length(compMot);

invoke(mh, 'processCommand',['ReDim ArrayOfValues(',num2str(nComp-1),')']);
for k = 1 : nComp
    invoke(mh, 'processCommand',['ArrayOfValues(',num2str(k-1),')="',compMot{k},'"']);
end

invoke(mh, 'processCommand', 'CALL getDocument().makeMotionComponent(ArrayOfValues)');
invoke(mh, 'processCommand', ['CALL getDocument().setMotionSourceType("',compMotName,'", infoVelocityDriven)']);
invoke(mh, 'processCommand', ['CALL getDocument().setParameter("',compMotName,'", "PositionVsTime", "[0%ms,0%deg,',permec,'%ms,',posmec,'%deg]", infoArrayParameter)']);
invoke(mh, 'processCommand', ['CALL getDocument().setMotionRotaryCenter("',compMotName,'", Array(0, 0, 0))']);
invoke(mh, 'processCommand', ['CALL getDocument().setMotionRotaryAxis("',compMotName,'", Array(0, 0, 1))']);


magDataUpt = magData;