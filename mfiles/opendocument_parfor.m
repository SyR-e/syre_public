% ActiveFEMM (C)2006 David Meeker, dmeeker@ieee.org

function opendocument_parfor(fn,h)
callfemm_parfor([ 'open(' , quote(fn) , ')' ],h);

