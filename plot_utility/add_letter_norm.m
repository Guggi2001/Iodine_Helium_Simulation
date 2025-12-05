function t = add_letter_norm(letter, position, size, ax, z)
%ADD_LETTER 
%places a letter in the bottom left corner of the current axis
if nargin <1
    letter = 'a';
end

if nargin<2
    position = 'bottomleft';
end

if nargin<3
    size=14;
end

 color = [0.1, 0.1, 0.1];


% if nargin<4;
%    color = [0.1, 0.1, 0.1];
% end
if nargin<4
ax = gca;
end

 if nargin<5
  z = ax.ZLim(2);
 end



xmin = 0;
ymin = 0;
xwidth = 1;

ywidth = 1;

reverse = strcmp(ax.XDir, 'reverse');


 
if strcmp('bottomleft', position)
    t = text(ax, xmin + (0.05 + 0.9*reverse)*xwidth,ymin + 0.1*ywidth,z,letter,  'color', color, 'fontsize', size,'Units','normalized');
end

if strcmp('topleft', position)
    t = text(ax,xmin + (0.05 + 0.90*reverse)*xwidth,ymin + 0.9*ywidth,z,letter,'color', color, 'fontsize', size, 'Units','normalized');
end

if strcmp('bottomright', position) 
    t = text(ax,xmin + (0.8 - 0.7*reverse)*xwidth,ymin + 0.1*ywidth,z,letter, 'color', color, 'fontsize', size,'Units','normalized');
end


if strcmp('topright', position) 
    t = text(ax, xmin + (0.8 - 0.7*reverse)*xwidth,ymin + 0.9*ywidth,z,letter, 'color', color, 'fontsize', size, 'Units','normalized');
enda

end

