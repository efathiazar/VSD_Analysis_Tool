function r_ellipse = fit_ellipse2(rois)% input ellipse parameters
    theta_grid = linspace(0,2*pi);
    phi = (rois(5))*180/pi;
    X0=rois(1);
    Y0=rois(2);
    a=rois(3);
    b=rois(4);

    % the ellipse in x and y coordinates 
    ellipse_x_r  = X0 + a*cos( theta_grid );
    ellipse_y_r  = Y0 + b*sin( theta_grid );

    %Define a rotation matrix
    R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

    %let's rotate the ellipse to some angle phii
    r_ellipse = (R * [ellipse_x_r;ellipse_y_r])';
end