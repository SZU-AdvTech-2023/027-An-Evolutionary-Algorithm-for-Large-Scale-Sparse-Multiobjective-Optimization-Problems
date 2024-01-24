function [Offspring, out_index_list, chosen_groups] = Group_Operate(Problem, Parent1, Parent2, number_of_groups)

    [proC, disC, ~, disM] = deal(1, 20, 1, 20);
    [N, D] = size(Parent1);
    
    %% Genetic operators for real encoding
    %% 对实数编码的决策变量进行遗传算子操作
    beta = zeros(N, D);
    mu   = rand(N, D);
    beta(mu <= 0.5) = (2 * mu(mu <= 0.5)).^(1 / (disC + 1));
    beta(mu > 0.5)  = (2 - 2 * mu(mu > 0.5)).^(-1 / (disC + 1));
    beta = beta .* (-1).^randi([0, 1], N, D);
    beta(rand(N, D) < 0.5) = 1;
    beta(repmat(rand(N, 1) > proC, 1, D)) = 1;
    Offspring = (Parent1 + Parent2) / 2 + beta .* (Parent1 - Parent2) / 2;
    Lower = repmat(Problem.lower, N, 1);
    Upper = repmat(Problem.upper, N, 1);
    [out_index_list, ~] = CreateGroups(number_of_groups, Offspring, D); 
    chosen_groups = randi(number_of_groups, size(out_index_list, 1), 1);
    Site = out_index_list == chosen_groups;
    mu = rand(N, 1);
    mu = repmat(mu, 1, D);
    temp = Site & mu <= 0.5;
    Offspring = min(max(Offspring, Lower), Upper);
    Offspring(temp) = Offspring(temp) + (Upper(temp) - Lower(temp)) .* ((2 .* mu(temp) + (1 - 2 .* mu(temp)) .*...
                      (1 - (Offspring(temp) - Lower(temp)) ./ (Upper(temp) - Lower(temp))).^(disM + 1)).^(1 / (disM + 1)) - 1);
    temp = Site & mu > 0.5;
    Offspring(temp) = Offspring(temp) + (Upper(temp) - Lower(temp)) .* (1 - (2 .* (1 - mu(temp)) + 2 .* (mu(temp) - 0.5) .*...
                      (1 - (Upper(temp) - Offspring(temp)) ./ (Upper(temp) - Lower(temp))).^(disM + 1)).^(1 / (disM + 1)));    
end


function [out_index_array, number_of_groups_array] = CreateGroups(number_of_groups, x_prime, number_of_variables)
    % 使用按序分组

    out_index_array = [];
    number_of_groups_array = [];
    no_of_solutions = size(x_prime, 1);
    for sol = 1:no_of_solutions        
        vars_per_group = floor(number_of_variables / number_of_groups);
        vars = x_prime(sol, :);
        [~, I] = sort(vars);
        out_index_list = ones(1, number_of_variables);
        for i = 1:number_of_groups-1
            out_index_list(I(((i-1) * vars_per_group) + 1:i * vars_per_group)) = i;
        end
        out_index_list(I(((number_of_groups - 1) * vars_per_group) + 1:end)) = number_of_groups;    
        out_index_array = [out_index_array; out_index_list];
        number_of_groups_array = [number_of_groups_array; number_of_groups];    
    end
end
