function [off_dec, off_mask] = GeneOperator(Problem, parent_dec, parent_mask, fitness)
    
    [N, ~]       = size(parent_dec);
    p1_dec      = parent_dec(1:floor(end/2), :);
    p2_dec      = parent_dec(floor(end/2)+1:floor(end/2)*2, :);
    p1_mask     = parent_mask(1:floor(end/2), :);
    p2_mask     = parent_mask(floor(end/2)+1:floor(end/2)*2, :);
    
    %% dec的杂交和变异
    if any(Problem.encoding ~= 4)
        [off_dec, groupIndex, chosengroups] = Group_Operate(Problem, p1_dec, p2_dec, 4); % 4为组数
        off_dec(:, Problem.encoding == 4) = 1;
    else
        off_dec = ones(size(p1_dec));
    end

    %% mask的杂交
    off_mask = p1_mask;
    for i = 1 : N/2
        if rand < 0.5
            temp = find(p1_mask(i, :) & ~p2_mask(i, :));
            temp = temp(TS(-fitness(temp)));
            off_mask(i, temp) = 0;
        else
            temp = find(~p1_mask(i, :) & p2_mask(i, :));
            temp = temp(TS(fitness(temp)));
            off_mask(i, temp) = p2_mask(i, temp);
        end
    end

    %% mask的变异
    if any(Problem.encoding ~= 4)
        chose_temp = groupIndex == chosengroups;
        for i = 1 : N/2
            if rand < 0.5
                temp = find(off_mask(i, :) & chose_temp(i, :));
                temp = temp(TS(-fitness(temp)));
                off_mask(i, temp) = 0;
            else
                temp = find(~off_mask(i, :) & chose_temp(i, :));
                temp = temp(TS(fitness(temp)));
                off_mask(i, temp) = 1;
            end
        end
    end
end

function temp = TS(fitness)
    % 二元锦标赛
    if isempty(fitness)
        temp = [];
    else
        temp = TournamentSelection(2, 1, fitness);
    end
end
