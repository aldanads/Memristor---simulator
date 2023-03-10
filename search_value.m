

function data = search_value(root,target)

    % Base case --> We go all the way down and the node is null
    % ATTENTION! -> Take care, some of them are CELLS and DOUBLES
    if iscell(root.data)
        aux = root.data{1}(1);
    else
        aux = root.data;
    end

    % If it is None, is the leaf
    if ~isempty(root.left) && ~iscell(root.left.data)
        aux_l = root.left.data;
    elseif ~isempty(root.left)
        aux_l = root.left.data{1}(1);
    end
    
    if isempty(root)
        data = false;
        return

    % We find the value among the leafs - Only interested in leafs, so we set the two extra conditions:
    % root.left and root.right is None.
    % The target should be smaller than the leaf we are checking to select that node
    elseif isempty(root.left) && isempty(root.right) && target <= aux
        data = root.data;
        return
    else

        if target <= aux_l
            data = search_value(root.left,target);
        else
            % Remember: nodes are the sum of their children. 
            % Because of that, if the target is greater than the left side, 
            % we operate target - root.left.data
            data = search_value(root.right,target-aux_l);
        end
    end



end