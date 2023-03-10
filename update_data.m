% Update the balanced binary tree which node as a sum of the value of the
% children.
% Starting from the leaf to the top, so the root will be the sum of all the
% leafs.
% For sorted trees, it improves binary search

function root = update_data(root)
    if isempty(root)
        return;
    end

    % We start at the leafs
    root.left = update_data(root.left);
    root.right = update_data(root.right);

    % Superior nodes are the sum of their children
    if ~isempty(root.left) && ~isempty(root.right)
        if iscell(root.left.data)
            aux_l = root.left.data{1}(1);
        else
            aux_l = root.left.data;
        end

        if iscell(root.right.data)
            aux_r = root.right.data{1}(1);
        else
            aux_r = root.right.data;
        end
        root.data = aux_l + aux_r;

    elseif ~isempty(root.left)
        if iscell(root.left.data)
            root.data = root.left.data{1}(1);
        else
            root.data = root.left.data;
        end
    elseif ~isempty(root.right)
        if iscell(root.right.data)
            root.data = root.right.data{1}(1);
        else
            root.data = root.right.data;
        end
    end

end