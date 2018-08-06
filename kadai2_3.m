clear;
%close;

A = [3 0.5;0.5 1];
mu = [1;2];
eta = max(eig(2 * A))^(-1) * 0.01;
lambda = [2,4,6];

N = 1000;
ht = 0;
epsiron = 0.00001;
results = zeros(N,3);
w_log = zeros(N,2,3);
for lambda_index = 1:3
    w = [10;12];
    for i = 1:N
                
        ht = ht + (2 .* A * (w - mu))' * (2 .* A * (w - mu));
        w = w - eta / (sqrt(ht) + epsiron) .* (2 .* A * (w - mu));
        w_log(i,:,lambda_index) = w';
        results(i,lambda_index) = (w - mu)' * A * (w - mu) + lambda(lambda_index) * norm(w,1);
    end
end

hold on;
w_results = zeros(N,3);
w_log = w_log - w_log(N,:,:);
for i = 1:3
    for j = 1:N
        w_results(j,i) = norm(w_log(j,:,i));
    end
end
for i = 1:3
    plot(1:N,log10(w_results(:,i)));
end
legend('PG:\lambda = 2','PG:\lambda = 4','PG:\lambda = 6','APG:\lambda = 2','APG:\lambda = 4','APG:\lambda = 6','AdaG:\lambda = 2','AdaG:\lambda = 4','AdaG:\lambda = 6');

