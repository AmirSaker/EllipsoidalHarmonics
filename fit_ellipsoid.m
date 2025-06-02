function [a, b, c] = fit_ellipsoid(v, plot_fig)
    % v : N x 3 array of points (x, y, z)
    % plot_fig : boolean to plot fitted ellipsoid

    x = v(:,1);
    y = v(:,2);
    z = v(:,3);

    % Least squares fitting to quadratic surface
    A = [x.^2, y.^2, z.^2];
    b_vec = ones(size(x));
    coeffs = A \ b_vec;

    % Extract semi-axes from fitted inverse-squares
    a = sqrt(1 / coeffs(1));
    b = sqrt(1 / coeffs(2));
    c = sqrt(1 / coeffs(3));

    % Optional plotting
    if plot_fig
        [X,Y,Z] = ellipsoid(0, 0, 0, a, b, c, 50);
        figure;
        surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;
        scatter3(x, y, z, 10, 'filled');
        axis equal;
        title('Fitted Ellipsoid');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        legend('Ellipsoid', 'Point Cloud');
        hold off;
    end
end
