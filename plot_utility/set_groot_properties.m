objectType = 'Line';
propertyName = 'LineWidth';
set(groot, ['Default', objectType, propertyName], 1.6)

objectType = 'Axes';
propertyName = 'FontSize';
set(groot, ['Default', objectType, propertyName], 15)

objectType = 'Scatter';
propertyName = 'Marker';
set(groot, ['Default', objectType, propertyName], 'o')
propertyName = 'Linewidth';
set(groot, ['Default', objectType, propertyName], 1.3)

objectType = 'Surface';
propertyName = 'EdgeColor';
set(groot, ['Default', objectType, propertyName], 'none')

global screensize
screensize = get( groot, 'Screensize' );
