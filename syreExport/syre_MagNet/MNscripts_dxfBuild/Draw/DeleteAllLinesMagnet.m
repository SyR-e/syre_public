function DeleteAllLinesMagnet(magData,x1,y1,x2,y2)


mh = magData.magnetHandler;

invoke(mh,'processCommand',['CALL getDocument().getView().selectIn(',num2str(x1),', ',num2str(y1),', ',num2str(x2),', ',num2str(y2),', infoSetSelection, Array(infoSliceLine, infoSliceArc))']);
invoke(mh,'processCommand','CALL getDocument().getView().deleteSelection()');
