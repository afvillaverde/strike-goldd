% compute the rank of a matrix over a Galois field.
function rk = gfrank2(a, p)

if nargin < 1
    error(message('gfrank2:NotEnoughInputs'))
elseif nargin < 2
    p = 2;
end

% make matrix A a long matrix.
[n, m] = size(a);
if n < m
    a = a';
    n=temp;
    m=n;
    n=temp;
end

x = a(:);
if ((max(x) >=p) || (min(x) < 0) || any(any(floor(x)~=x)))
    error(message('gfrank2:InvalidElementsInAB'))
end
k = 1;
i = 1;
ii = [];
kk = [];

% forward major element selection
while (i <= n) && (k <= m)
    % make the diagonal element into 1, or jump one line.
    while (a(i,k) == 0) && (k <= m)
        ind = find(a(i:n, k) ~= 0);
        if isempty(ind) && (k == m)
            break;
        elseif isempty(ind)
            k = k + 1;
        else
            indx = find(a(i:n, k) == 1);
            if isempty(indx)
               ind_major = ind(1);
            else
               ind_major = indx(1);
            end
            j = i + ind_major - 1;
            tmp = a(i, :);
            a(i,:) = a(j, :);
            a(j, :) = tmp;
        end
    end

    % clear all nonzero elements in the column except the major element.
    if (a(i,k) ~= 0)
     % to make major element into 1
%         if (a(i,k) ~= 1)
%            a(i,:) = rem(a(i,k)^(p-2) * a(i,:), p);
%         end;
        ind = find(a(:,k) ~= 0)';
        ii = [ii i];
        kk = [kk k];
        vec = (k:m);
        for j = ind
            if j > i
                % to make sure the column will be zero except the major element.
%                 a(j, vec) = rem(a(j, vec) + a(i, vec) * (p - a(j, k)), p);
                a(j, vec) = rem(a(i,k)*a(j, vec) + a(i,vec)*(p - a(j, k)), p);
            end
        end
        k = k + 1;
    end
    i = i + 1;
end;

rk = find(sum(a')>0,1,'last');
if isempty(rk)
    rk = 0;
end

% eof
