n1 = 4:5;
n2 = 0:4;
f1 = 10:2:24;
f2 = 10:2:24;
matrix_of_index = ones(length(n1), length(f1), length(n2), length(f2));
As = @(n, f) n*pi*f^2/4;
pos = struct('n1', nan, 'f1', nan, 'n2', nan, 'f2', nan, 'area', nan);
pos = repmat(pos, length(matrix_of_index),1);

for r = 1:numel(matrix_of_index)
    index = cell(length(size(matrix_of_index)), 1);
    [index{:}] = ind2sub(size(matrix_of_index), r);
    a = As(n1(index{1}), f1(index{2})) + As(n2(index{3}), f2(index{4}));
    pos(r).n1 = n1(index{1});
    pos(r).f1 = f1(index{2});
    pos(r).n2 = n2(index{3});
    pos(r).f2 = f2(index{4});
    pos(r).area = a;
end

pos = struct2table(pos);
[~, sortingIndex] = sort(pos.area);
posSorted = pos(sortingIndex, :);


        