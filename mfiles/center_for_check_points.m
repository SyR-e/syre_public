function [beta1,beta2]=center_for_check_points(temp,geo)

                if isfield(temp,'xxD1k')
                    Rbeta=geo.x0-temp.Bx0;
                    beta1=atan(temp.yyD1k./(x0-temp.xxD1k));
                    beta2=atan(temp.yyD2k./(x0-temp.xxD2k));
                end
end