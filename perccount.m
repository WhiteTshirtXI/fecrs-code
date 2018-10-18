% Reports the percentage of the job done to the screen
%
% perccount(jj,maxjj)
%
% ========================================================================
% Inputs
% ========================================================================
% I = current iteration between 1 and Imax
% Imax = maximum number of iterations

function  perccount(jj,maxjj)

persistent lastCall;
if(nargin == 2)
    if isempty(lastCall)
        lastCall = -1;
    end
    if(lastCall ~= floor(((jj-1)/maxjj) * 100))
        if(jj ~= 1)
            fprintf(1,'\b\b\b');
        else
            fprintf(1,'\n\tPercentage complete: ');
        end
        pc_done = num2str(floor(((jj-1)/maxjj) * 100));
        if(length(pc_done) == 1)
            pc_done(2) = pc_done(1);
            pc_done(1) = '0';
        end
        fprintf(1,'%s%%',pc_done);
    end
    
    lastCall = floor(((jj-1)/maxjj) * 100);
    if(jj == maxjj)
        fprintf(1,'\b\b\b100%%\n\n');
    end
else
    error('Error: PERCCOUNT needs two input arguments...');
end