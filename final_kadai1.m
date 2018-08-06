clear;
close;
hold on;
%generate training data
n = 100;
x = 3 * (rand(n, 2) - 0.5);
radius = x(:, 1).^2 + x(:, 2).^2;
y = (radius > 0.7 + 0.1 * randn(n, 1)) & (radius < 2.2 + 0.1 * randn(n, 1));
y = 2 * y -1;

%plot(x(:,1),x(:,2),'x');
w = [1,1];
x = x';
lambda = 1;
alpha = 0.01;

N = 100;
Jarr = zeros(N,1);
for training_itr = 1:N
    round_J = 0;
    J = 0;
    for i = 1 : n
        round_J = round_J - ( y(i) * exp(-y(i) * w * x(:,i)) ) .* x(:,i)' ./ ( 1 + exp(-y(i) * w * x(:,i)));
        J = J + log(1 + exp(- y(i) * w * x(:,i)) ) ;
    end
    round_J = round_J + 2 * lambda * w;
    J = J + lambda * (w * w');
    Jarr(training_itr) = J;
    w = w - alpha * round_J;
end

plot(Jarr);
alpha = 1;
w = [1,1];
for training_itr = 1:N
    round_J = 0;
    roundround_J = 0;
    J = 0;
    for i = 1 : n
        ye_scalar = - y(i) * exp(-y(i) * w * x(:,i));
        bunshi = 1 + exp(-y(i) * w * x(:,i));
        
        round_J = round_J + ye_scalar .* x(:,i)' ./ bunshi;
        roundround_J = roundround_J + ( bunshi .* ( - y(i) * ye_scalar .* (x(:,i) * x(:,i)') ) - ye_scalar^2 .* ( x(:,i) * x(:,i)' )) ./ bunshi^2;
        
        J = J + log(1 + exp(- y(i) * w * x(:,i)) );
    end
    round_J = round_J + 2 * lambda * w;
    roundround_J = roundround_J + 2 * lambda * eye(2);
    J = J + lambda * (w * w');
    Jarr(training_itr) = J;
    w = w - alpha .* (roundround_J^(-1) * round_J')';
end

plot(Jarr);

legend('Steepest Descent','Newton');

