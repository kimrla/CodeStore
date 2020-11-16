% �����������͵�B�������ߣ���Ҫǰ������������.m�ļ�
clear all;
%���ƶ���
P = [9.036145, 21.084337, 37.607573, 51.893287, 61.187608;
    51.779661, 70.084746, 50.254237, 69.745763, 49.576271];
 
n = 4; k = 2;
 
flag = 2;
% flag = 1�����ƾ���B��������
% flag = 2, ����׼����B��������
% flag = 3, ���Ʒֶ�Bezier����
 
switch flag
    case 1
        NodeVector = linspace(0, 1, n+k+2); % ����B�����Ľڵ�ʸ��
 
        % ������������
        plot(P(1, 1:n+1), P(2, 1:n+1),...
                        'o','LineWidth',1,...
                        'MarkerEdgeColor','k',...
                        'MarkerFaceColor','g',...
                        'MarkerSize',6);
        line(P(1, 1:n+1), P(2, 1:n+1));
        Nik = zeros(n+1, 1);
        for u = k/(n+k+1) : 0.001 : (n+1)/(n+k+1)
            % for u = 0 : 0.005 : 1
            for i = 0 : 1 : n
                Nik(i+1, 1) = BaseFunction(i, k , u, NodeVector);
            end
        p_u = P * Nik;
        line(p_u(1,1), p_u(2,1), 'Marker','.','LineStyle','-', 'Color',[.3 .6 .9]);
        end
    case 2
        NodeVector = U_quasi_uniform(n, k); % ׼����B�����Ľڵ�ʸ��
        DrawSpline(n, k, P, NodeVector);
    case 3
        NodeVector = U_piecewise_Bezier(n, k);  % �ֶ�Bezier���ߵĽڵ�ʸ��
        DrawSpline(n, k, P, NodeVector);
    otherwise
        fprintf('error!\n');
end