%filename: side_check2.m
%Input:     Matrix of three point: (2x3 or 3x2)
%           P(p_x,p_y),Q(q_x,q_y),R(r_x,r_y)
%Output:    ���P Q R���㹹��һ���ҹգ�����True
%           �ҹ���ζ�ŵ�R��PQ�������Ҳ�.��ʱ

function result = side_check2(D)
if all(size(D) ~= [3,2])
    if all(size(D)==[2,3])
        D = D';
    else
        error('error dimension')
    end
end
result = (det([[1;1;1], D]) < 0 );