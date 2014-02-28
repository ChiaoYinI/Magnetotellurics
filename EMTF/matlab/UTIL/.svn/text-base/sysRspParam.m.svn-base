function param = sysRspParam(type,str)

% Replicates EMTF/D/fcorsu.f for some filters:
% given a filter type and a string as recorded
% in SP files, creates a parameter structure.
%
% Example filters:
% 'PZ', '0 3 1984.31   (-6.28319,10.8825) (-6.28319,-10.8825) (-12.5664,0.0)'
% 'DT', '.2455'
%
% Last modified by A. Kelbert, 25 Apr 2008

switch type
    
    case 'PZ'
        % First read the info from the string
        [val,count,errmsg,nchar] = sscanf(str,'%d',2);
        if count == 2
            nzeros = val(1);
            npoles = val(2);
        else
            error('Unable to read the number of poles and zeros from PZ filter');
        end
        ind = nchar;
        [val,count,errmsg,nchar] = sscanf(str(ind:end),'%f',1);
        if count == 1
            param.a0 = val;
        else
            error('Unable to read a0 from PZ filter');
        end
        ind = ind+nchar;
        [val,count,errmsg,nchar] = sscanf(str(ind:end),' (%f,%f)',2*nzeros);
        if count == 2*nzeros
            param.zeros = [];
            for j = 1:nzeros
                param.zeros(j) = val(2*j-1)+i*val(2*j);
            end
        else
            error('Unable to read the zeros PZ filter');
        end
        ind = ind+nchar;
        [val,count,errmsg,nchar] = sscanf(str(ind:end),' (%f,%f)',2*npoles);
        if count == 2*npoles
            param.poles = [];
            for j = 1:npoles
                param.poles(j) = val(2*j-1)+i*val(2*j);
            end
        else
            error('Unable to read the poles from PZ filter');
        end
        
    case 'DT'
        % Read the time offset from the string
        [val,count,errmsg,nchar] = sscanf(str,'%f',2);
        if count == 1
            param = val;
        else
            error('Unable to read the time offset from DT filter');
        end

    otherwise
        error('Sorry, this filter not yet implemented');
end