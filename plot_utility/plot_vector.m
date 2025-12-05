function plot_vector(base, v, scale)
if nargin<3
    scale = 1;
end

quiver3(base(1),base(2),base(3), v(1), v(2),v(3), scale);

end

