function result = LDPC(H1,bit_piece,maxit)
% % Test:
% clear all; close all; clc;
% H1 = [0 1 0 1 1 0 0 1; 1 1 1 0 0 1 0 0; 0 0 1 0 0 1 1 1; 1 0 0 1 1 0 1 0];
% bit_piece = [1 1 0 1 1 1 0 1];
% maxit = 3;
[rows,col] = find(H1);
[r1,c1] = size(H1);
result = bit_piece;
syndrome = mod(result*H1',2);
iter = 0;
while (norm(syndrome) ~= 0 && iter < maxit)
    iter = iter + 1;  
    for i = 1:c1
        c_nodes_fs = rows(col==i);   % f nodes connected to the c node.
        H_nodes =  H1(c_nodes_fs,:);   % Rows of H that represents the f nodes.
        [r2,c2] = size(H_nodes);
        [rows1,col1] = find(H_nodes);
        A = repmat(result,r2,1);   % Matrix of encoded bits of the same size as the H_nodes.
        v_nodes_c = H_nodes.*A;
        nodes_val = mod(sum(v_nodes_c,2),2);
        if(norm(nodes_val)>0)
            diff_val = nodes_val(nodes_val>0);
            same_val = nodes_val(nodes_val ==0);
            [r3,c3] =  size(diff_val);
            [r4,c4] = size(same_val);
            B = repmat(result,r3,1);
            vote_fn = bitxor(diff_val,B(:,i));  % This gives us the values given by the f nodes.
            vote_f = [vote_fn;repmat(result(i),r4,1)];
        else
            vote_f = A(:,i);
        end %if
        
        vote = (sum(vote_f)+result(i));
        if(vote>(length(vote_f)+1)/2)
            res = 1;
        elseif (vote == (length(vote_f)+1)/2)
            res = result(i);
        else
            res = 0;
        end %if
        result(i) = res;
        syndrome = mod(result*H1',2);
        
        if (iter < maxit && norm(syndrome) ~= 0)
            A = repmat(result,r2,1);   % Matrix of encoded bits of the same size as the H_nodes.
            v_nodes_c = H_nodes.*A;
            v_nodes_c(:,i) = [];
            nodes_val = mod(sum(v_nodes_c,2),2);
            if(norm(nodes_val)>0)
                diff_val = nodes_val(nodes_val>0);
                same_val = nodes_val(nodes_val ==0);
                [r3,c3] =  size(diff_val);
                [r4,c4] = size(same_val);
                B = repmat(result,r3,1);
                vote_fn = bitxor(diff_val,B(:,i));  % This gives us the values given by the f nodes.
                vote_f = [vote_fn;repmat(result(i),r4,1)];
            else
                vote_f = A(:,i);
            end %if
            
            for m = 1:length(vote_f)
                ignore_f = vote_f;
                ignore_f(m) = [];
                res = zeros(length(vote_f),1);
                vote = (sum(ignore_f)+result(i));
                if(vote>length(vote_f)/2)
                    res = 1;
                elseif (vote == length(vote_f)/2)
                    res = result(i);
                else
                 	res = 0;
                end %if
            end 
            result(i) = res;
        end
    end % for
    syndrome = mod(result*H1',2);
end
result = result.';