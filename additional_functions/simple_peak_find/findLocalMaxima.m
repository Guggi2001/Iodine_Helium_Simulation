function localMaxima_filtered = findLocalMaxima(y, m, f, mindiff)






    % Ensure m is odd to have a center point
    if mod(m, 2) == 0
        error('Window size m must be odd.');
    end

    % Calculate the offset for the center point
    offset = floor(m / 2);

    % Initialize the list of local maxima
    localMaxima = [];

    % Iterate over the data vector with the sliding window
    for i = (offset + 1):(length(y) - offset)
        % Extract the current window
        window = y((i - offset):(i + offset));

        % Check if the center point is larger than f% of the entries around it
        if sum(y(i)>window)/numel(window)>f/100
            localMaxima = [localMaxima, i];
        end
    end


    ids = localMaxima;
    
   localMaxima_filtered = [];

   while length(ids)>0


        diff = ids - ids(1)

        b_diff = diff<mindiff;
        
        ids_in_cluster = ids(b_diff);


        id_max = find(y(ids_in_cluster) ==max(y(ids_in_cluster)), 1)
        

       
        localMaxima_filtered = [localMaxima_filtered, ids_in_cluster(id_max)];

        ids(b_diff) = [];
        

   end





end
