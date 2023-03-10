% Build balanced tree

function root = build_tree(arr)
    if isempty(arr)
        root = [];
        return
    end
    if length(arr) == 1
        root = Node(arr);
        return
    end
    
        mid = floor(length(arr)/2);
        root = Node([]);
        root.left = build_tree(arr(1:mid));
        root.right = build_tree(arr(mid+1:end));
    
end



