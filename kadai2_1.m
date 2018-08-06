clear;
close;
A = [3 0.5;0.5 1];
mu = [1;2];
eta = max(eig(2 * A))^(-1) * 0.1;
lambda = [2,4,6];

N = 1000;

results = zeros(N,3);
w_log = zeros(N,2,3);
for lambda_index = 1:3
    w = [11;12];
    for i = 1:N
        temp = w - eta .* (2 * A * (w - mu) ) ;
        
        for index = 1:2
            if temp(index) > lambda(lambda_index) * eta
                w(index) = temp(index) - lambda(lambda_index) *eta;
                
            elseif abs(temp(index) ) < lambda(lambda_index) * eta
                w(index) = 0;
                
            elseif temp(index) < - lambda(lambda_index) * eta
                w(index) = temp(index) + lambda(lambda_index) * eta;
                
            end
        end
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
    semilogx(1:N,log10(w_results(:,i)));
end
xlabel('t');
ylabel('log10( $$\| \omega(t) - \hat{\omega} \|$$)','Fontsize',12,'Interpreter','latex'); 
%legend('PG:\lambda = 2','PG:\lambda = 4','PG:\lambda = 6');


