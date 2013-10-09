% set xlabel
if(isempty(char(datalab.t_unit)))
    xlab = 't';
else
    xlab = strcat('t [',datalab.t_unit,']');
end

%  
indjump = find(diff(datalab.t) == 0);
jumps = size(indjump,1);

win = 1
% open new window
figure(win);
% set subplot properties
nrows = double(datalab.nysub(win));
ncols = double(datalab.nxsub(win));
nsubplt = double(datalab.nplt(win));

i = 0;

% draw objective function
if(not(isempty(datalab.obj)))
    if(i > nsubplt)
        win = win+1;
        % open new window
        figure(win);
        % set subplot properties
        nrows = double(datalab.nysub(win));
        ncols = double(datalab.nxsub(win));
        nsubplt = double(datalab.nplt(win));
        i = -1;
    end
    subplot(nrows,ncols,1); 
    bar(datalab.obj(1:datalab.iter_idx))
    xlabel('SQP Iteration')
    % set ylabel
    if(isempty(char(datalab.obj_unit)))
        ylabel(strcat('obj'))
    else
        ylabel(strcat('obj [',datalab.obj_unit,']'))
    end
    % set title
    title(datalab.obj_title,'FontSize',12)
    i = i+1;
end

% draw model state duration
if(not(isempty(datalab.h)))
    [m,n] = size(datalab.h);
    for k = 1:n
        if(i+k > nsubplt)
            win = win+1;
            % open new window
            figure(win);
            % set subplot properties
            nrows = double(datalab.nysub(win));
            ncols = double(datalab.nxsub(win));
            nsubplt = double(datalab.nplt(win));
            i = -k+1;
        end
        subplot(nrows,ncols,k+i); 
        plot(datalab.h(1:datalab.iter_idx,k),'+k','MarkerSize',5)
        xlabel('SQP Iteration')
        % set ylabel
        if(isempty(char(datalab.h_unit(k))))
            ylabel(strcat('h_',int2str(k-1)))
        else
            ylabel(strcat('h_{',int2str(k-1),'} [',datalab.h_unit(k),']'))
        end
        % set title
        title(datalab.h_title(k),'FontSize',12)
    end
    i = i+n;
end

% draw global model parameter
if(not(isempty(datalab.p)))
    [m,n] = size(datalab.p);
    for k = 1:n
        if(i+k > nsubplt)
            win = win+1;
            % open new window
            figure(win);
            % set subplot properties
            nrows = double(datalab.nysub(win));
            ncols = double(datalab.nxsub(win));
            nsubplt = double(datalab.nplt(win));
            i = -k+1;
        end
        subplot(nrows,ncols,k+i); 
        plot(datalab.p(1:datalab.iter_idx,k),'+k','MarkerSize',5)
        xlabel('SQP Iteration')
        % set ylabel
        if(isempty(char(datalab.p_unit(k))))
            ylabel(strcat('p_{',int2str(k-1),'}'))
        else
            ylabel(strcat('p_{',int2str(k-1),'} [',datalab.p_unit(k),']'))
        end
        % set title
        title(datalab.p_title(k),'FontSize',12)
    end
    i = i+n;
end

% draw control function
if(not(isempty(datalab.u)))
    [m,n]=size(datalab.u);
    for k=1:n
        if(i+k > nsubplt)
            win = win+1;
            % open new window
            figure(win);
            % set subplot properties
            nrows = double(datalab.nysub(win));
            ncols = double(datalab.nxsub(win));
            nsubplt = double(datalab.nplt(win));
            i = -k+1;
        end
        subplot(nrows,ncols,k+i);  
        plot(datalab.t,datalab.u(:,k),'-g')
        xlabel(xlab)
        % set ylabel
        if(isempty(char(datalab.u_unit(k))))
            ylabel(strcat('u_{',int2str(k-1),'}(t)'))
        else
            ylabel(strcat('u_{',int2str(k-1),'}(t)',' [',datalab.u_unit(k),']'))
        end
        % set title
        title(datalab.u_title(k),'FontSize',12)
        grid on
    end
    i = i+n;
end

% draw differential state function
if(not(isempty(datalab.xd)))
    [m,n]=size(datalab.xd);
    for k=1:n
        if(i+k > nsubplt)
            win = win+1;
            % open new window
            figure(win);
            % set subplot properties
            nrows = double(datalab.nysub(win));
            ncols = double(datalab.nxsub(win));
            nsubplt = double(datalab.nplt(win));
            i=-k+1;
        end
        subplot(nrows,ncols,k+i);
        hold off
        plot(datalab.t(1:indjump(1)),datalab.xd(1:indjump(1),k),'-b')
        hold on
        for j = 2:jumps
            plot(datalab.t(indjump(j-1)+1:indjump(j)),datalab.xd(indjump(j-1)+1:indjump(j),k),'-b')
        end
        xlabel(xlab)
        % set ylabel
        if(isempty(char(datalab.xd_unit(k))))
            ylabel(strcat('x_{',int2str(k-1),'}(t)'))
        else
            ylabel(strcat('x_{',int2str(k-1),'}(t)',' [',datalab.xd_unit(k),']'))
        end
        % set title
        title(datalab.xd_title(k),'FontSize',12)
        grid on
    end
    i = i+n;
end

% draw algebraic state function
if(not(isempty(datalab.xa)))
    [m,n] = size(datalab.xa);
    for k = 1:n
        if(i+k > nsubplt)
            win = win+1;
            % open new window
            figure(win);
            % set subplot properties
            nrows = double(datalab.nysub(win));
            ncols = double(datalab.nxsub(win));
            nsubplt = double(datalab.nplt(win));
            i = -k+1;
        end
        subplot(nrows,ncols,k+i); 
        hold off
        plot(datalab.t(1:indjump(1)),datalab.xa(1:indjump(1),k),'-r')
        hold on
        for j = 2:jumps
            plot(datalab.t(indjump(j-1)+1:indjump(j)),datalab.xa(indjump(j-1)+1:indjump(j),k),'-r')
        end
        xlabel(xlab)
        % set ylabel
        if(isempty(char(datalab.xa_unit(k))))
            ylabel(strcat('z_{',int2str(k-1),'}(t)'))
        else
            ylabel(strcat('z_{',int2str(k-1),'}(t)',' [',datalab.xa_unit(k),']'))
        end
        % set title
        title(datalab.xa_title(k),'FontSize',12)
        grid on
    end
end

drawnow

