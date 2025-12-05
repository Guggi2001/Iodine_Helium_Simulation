       global plot_processed_with_ROI
    plot_processed_with_ROI = false;
    % from 26.8.24, 600 mW probe, 40 bar, 14 K
    res = plot_processed_VMI(32910, 1, [528.4614  382.9408], 1);
    close all
    global velocity_factor

    vf = velocity_factor;
    vf =  6.4149*0.9675;
    
    v = res.r*vf;

    v_centers = 100:200:(1861-100);

    parameter_array = [];

    error_array = [];

       figure(2)
           
    leg = {};

    for vc = v_centers
    b_v= v>vc-100 & v<vc+100;

 
    y2 = mean(res.image_polar(:, b_v),2);
    y2 = movmean(y2,4);
    %y2 = y2 - min(y2);
    %y2 = y2/max(y2);



    fun2 = @(beta, phi) [real(beta(1)*cos(phi).^(beta(2))+ beta(3)) ; imag(real(beta(1)*cos(phi).^(beta(2))) + beta(3))] ;
    [beta2,R,J,COVB,MSE] = nlinfit(res.phi', [y2; 0*y2], fun2, [1,2, 0.1]);
    parci = nlparci(beta, R,"covariance",COVB);
    error = parci(:,2) - parci(:,1);

    hold on
    fitresult = fun2(beta2, res.phi')
    %plot(res.phi, fitresult(1:length(y2)));
    plot(res.phi,y2 );
    parameter_array = [parameter_array; beta2];
    error_array =[error_array; error'];

    leg{end+1} = sprintf("%.f m/s ", vc);

    end
    legend(leg);

    recolor_lines
    xlabel('v / m/s');
    

    % anisotropy
    %https://pyabel.readthedocs.io/en/latest/anisotropy_parameter.html
    %I = sigma*(1 + beta cos(theta))
    % my fit function is of the form  beta(1) cos(phi)^(beta(2)) + beta(3)
    % this means that the beta isotropy as it is usually defined is:
    % beta(1) / beta(3)

    anisotropy = parameter_array(:,1)./parameter_array(:,3);
    anisotropy_error = error_array(:,1)./parameter_array(:,3) + parameter_array(:,1)./parameter_array(:,3).^2 .* error_array(:,3);
    
    figure
    errorbar(v_centers, anisotropy, anisotropy_error);
    xlabel('v / m/s');
    ylabel('anisotropy');
