function load_fig_into_tiledlayout(figpath, tl,tile)
%LOAD_FIG_INTO_TILEDLAYOUT loads a figure from figpath, and copies its
%contents into a tiledlayout tile
%nexttile

original_children = tl.Children;

f = openfig(figpath);
copyobj(f.Children,tl);

close(f)


for child = tl.Children'
    if ~ismember(child, original_children) & isa(child, 'matlab.graphics.axis.Axes')
        child.Layout.Tile = tile;
        child.Title.String = sprintf('tile %.0f', tile);
        
    end

end

end

