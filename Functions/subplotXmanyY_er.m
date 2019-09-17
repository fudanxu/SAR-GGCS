function [ax] = subplotXmanyY_er(m,p,Listener)
%SUBPLOTXMANYY_ER Creates subplot axes with common x-axis and minimizes the
%white space between the different axes. The white space removal is
%performed when the figure's 'SizeChangedFcn' callback is
%initiated. Initiates after first figure resize or first time figure
%becomes visible.
%
%   ax = subplotXmanyY_er(m,p)
%   ax = subplotXmanyY_er(m,p,'listner')    - Note: Reduces performance
%
% Example code:
%
% t = (-10:0.1:10)';
% f = sin(t) + 0.2*randn(numel(t),1);
% dfdt = [0;diff(f)];
% figure('name','subplotXmanyY_er','Visible','off');
% ax(1) = subplotXmanyY_er(2,1);
% plot(t,f);ylabel('Dispalcement [m]');title('f(t) = sin(t) + noise');grid on
% set(gca,'XTickLabel',{})
% ax(2) = subplotXmanyY_er(2,2);
% plot(t,dfdt);legend('f''');ylabel('Velocity [m/s]');grid on
% xlabel('Time [s]');grid on
% set(gcf,'Visible','on')
% linkaxes(ax,'x');   % Link axes for zooming
%
% Author:           Eduard Reitmann
% Version:          1.0 (2017-04-21)
% Original date:    2017-04-21

n = 1;
if nargin == 0
    m = 1;p = 1;
end
% Position matrices [left bottom width height]
left = repmat(0:1/n:1-1/n,m,1);
bottom = repmat(fliplr(0:1/m:1-1/m)',1,n);
% Position
rows = repmat(1:m,n,1);rows = rows(:);
cols = repmat((1:n)',m,1);
row = rows(p);
col = cols(p);
% Perform checks
if max(p) > m*n
    error('Outside of range')
end
nposcomb = sum((rows >= min(row) & rows <= max(row)) & (cols >= min(col) & cols <= max(col)));
if numel(p) ~= nposcomb
    error('Not valid selection')
end
l = left(1,min(col));
b = bottom(max(row),1);
w = (1/n)*numel(unique(col));
h = (1/m)*numel(unique(row));
ax = axes('OuterPosition',[l b w h],...
    'Units','normalized',...
    'NextPlot','replacechildren');
setappdata(ax,'FixedOuterPosition',p);
set(gcf,'PaperPositionMode','manual','Units','normalized','SizeChangedFcn',@fillaxes)
%% Add additional listener
if nargin == 3
    % Add lsiterner to fillaxes when adding labels and titles -> Will reduce
    % performance !!! (Alternative: Set figure visible 'off' then 'on')
    if strcmpi(Listener,'listner')
        addlistener(ax,'ChildAdded',@fillaxes);
        addlistener(ax,'MarkedDirty',@fillaxes);
        warning('off','MATLAB:callback:error')
    end
end
end
function fillaxes(~,b)
f = b.Source;
% f = a.Parent;
% Default margin
margin = 0.005;
% Fill all axes
ax = findobj(f,'Type','axes');
n = numel(ax);
l = 0;
for i = 1:n
    P{i,1} = getappdata(ax(i),'FixedOuterPosition');
    TI = ax(i).TightInset;
    topbot(i,1) = sum(TI([2 4]));
    % Find maximum length
    if TI(1)+margin > l
        l = TI(1)+margin;
    end
end
% Find maximum width
w = 1;
for i = 1:n
    TI = ax(i).TightInset;
    if 1 - l - TI(3) - margin < w
        w = 1 - l - TI(3) - margin;
    end
end
m = numel(cell2mat(P'));
h = (1-sum(topbot)-(n+1)*margin)/m;
%% Pack axes from bottom
maxP = cellfun(@max, P);
[~,packorder] = sort(maxP,'descend');
bot = 0;
for i = 1:n
    ii = packorder(i);
    p = getappdata(ax(ii),'FixedOuterPosition');
    TI = ax(ii).TightInset;
    b = bot + margin + TI(2);    
    NewPosition = [l b w h*numel(p)];
    NewPosition(NewPosition<0) = 0;    
    if any(ax(ii).Position ~= NewPosition)
        ax(ii).Position = NewPosition;
    end    
    bot = bot + margin + h*numel(p) + topbot(ii);    
end
end