function u = new_knots (kc, kf)
% Find the new knots with their multiplicity//寻找kf比kc中多的节点，重复的元素也重复输出
[valc, multc] = unique (kc, 'last');
multc = diff ([0 multc(:)']);
[valf, multf] = unique (kf, 'last');
multf = diff ([0 multf(:)']);

unew = setdiff (kf, kc);
[~,posf] = ismember (unew, valf);
mult_new = multf(posf);

[urep, indc, indf] = intersect (valc, valf);
mult_rep = multf(indf) - multc(indc);
urep = urep(mult_rep>0);
mult_rep = mult_rep(mult_rep>0);

mult = [mult_new mult_rep];
u = [unew, urep];

ind = zeros (numel(kf)-numel(kc), 1);
ind(cumsum([1 mult(:)'])) = 1;
u = sort (u(cumsum(ind(1:end-1))));

end