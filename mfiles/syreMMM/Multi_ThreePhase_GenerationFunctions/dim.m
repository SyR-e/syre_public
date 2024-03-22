function dim(blockpath,d1,d2,d3,d4)
    pos = get_param(blockpath,'Position');
    set_param(blockpath,'Position',[pos(1)+d1 pos(2)+d2 pos(3)+d3 pos(4)+d4]);
end