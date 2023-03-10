%% Balanced binary tree 

classdef Node
    properties
        data
        left
        right
    end

    methods
        function obj = Node(data)
            obj.data = data;
            obj.left = [];
            obj.right = [];

        end
    end

end


