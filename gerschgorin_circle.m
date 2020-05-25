% Draw Gerschgorin Circle ¸Ç¶ûÔ²
A = [0 3 2;2 7 1;1 2 15]; % original matrix
T = [1 -0.3 -0.1;-0.2 -1 0.2;0 0.2 1];    % transform matrix

N = length(A);
A_trans = T*A/T;              % A after transformed
centers       = zeros(1,N);
centers_trans = zeros(1,N);
radius        = zeros(1,N);
radius_trans  = zeros(1,N);

% find center of circles and their radius
for i=1:N
    centers(i)       = A(i,i);
    centers_trans(i) = A_trans(i,i);
    radius(i)        = sum(abs(A(i,:))) - abs(centers(i));
    radius_trans(i)  = sum(abs(A_trans(i,:))) - abs(centers_trans(i));
end

% prepare plot data
theta = linspace(0,2*pi,100);
sin_theta = sin(theta);
cos_theta = cos(theta);
plot_x1 = zeros(3,100); % ori
plot_x2 = zeros(3,100); % trans
plot_y1 = zeros(3,100);
plot_y2 = zeros(3,100);
for i = 1:N
    plot_x1(i,:) = radius(i)*sin_theta + real(centers(i));
    plot_y1(i,:) = radius(i)*cos_theta + imag(centers(i));
    plot_x2(i,:) = radius_trans(i)*sin_theta + real(centers_trans(i));
    plot_y2(i,:) = radius_trans(i)*cos_theta + imag(centers_trans(i));
end

% plot
figure(1)
scatter(real(centers),imag(centers),'b+');hold on;   % center
scatter(real(centers_trans),imag(centers_trans),'r*');hold on;
for i = 1:N
    plot(plot_x1(i,:),plot_y1(i,:),'b');hold on;
    plot(plot_x2(i,:),plot_y2(i,:),'r--');hold on;
end

