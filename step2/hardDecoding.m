function result = hardDecoding(H,bit_rx,maxit)
    length_f=size(H,1); % number of f-nodes
    length_c=size(H,2); % number of c-nodes
    result = []; % corrected received vector
    decision = bit_rx; % The c-nodes first receive the received vector
    final_decision = decision;
    condition = mod(decision'*H',2); % syndrom
    iter = 0; % number of iterations
    decision=decision*ones(1,length_f);

    while(norm(condition) ~= 0 && iter < maxit) % while the syndrom is not zero 
        c=decision;
        c_next=zeros(length_c,length_f); % every f-node sends changes to all the c-nodes. Changes of f1 are stored in the first column of c_next for example
        for i=1:length_f % for every f-node
            f=H(i,:); 
            indexC=find(f==1);   % find the indexes of the c's connect to fi
            for indC=1:length(indexC)
                remInd=indexC;
                remInd(indC)=[]; % remove the contribution from c(k)
                check=mod(sum(c(remInd,i)),2); % modulo-2 sum of all the c's connected to fi, except c(k)
                if(check==0)
                    c_next(indexC(indC),i)=0; % if the sum is 0, c(k) should be 0
                else
                    c_next(indexC(indC),i)=1; % if the sum is 1, c(k) should be 1 such that the sum of all c's (c(k) included) gives 0
                end
            end
        end
        decision=zeros(length_c,1); % Necessary ?
        for j=1:length_c % for every c-node
            c_node=H(:,j); 
            indexF=find(c_node==1);   % check the f_nodes that are connected to c(j)
            f_decision = c_next(j,indexF); % the decision is made taking into account the responses of the f's connected to c(j)
            c_decision = [bit_rx(j) f_decision]; % add the decision of the received vector
            % QUESTION : Should be the decision of bit_rx or of
            % previous decision ?
             % Always the initial bit_rx received vector !
            nbVotes=length(c_decision); % majority voting
            if(sum(c_decision) == nbVotes/2)
                % if there are 50% vote for 1 and 50% for 0, the decision
                % is the original received bit
                decision(j) = bit_rx(j);
            elseif(sum(c_decision)>nbVotes/2)
                % There are more than 50% 1s --> the decision is 1
                decision(j) =1;
            else 
                decision(j) = 0;
            end
        end
        final_decision = decision;
        condition = mod(decision'*H',2); % syndrom, algorithm stops when syndrom is zero
        iter = iter+1;
            
        % For more than 1 iteration :
        % If the syndrom is not zero, the c for the next iteration should
        % be recomputed WITHOUT taking into account that node in the voting
        % If the syndrom is zero, the final decision is made with all nodes
        % voting
            
        if(norm(condition) ~= 0 && iter < maxit)
            decision=zeros(length_c,length_f); % Necessary ?
                
            for j=1:length_c % for every c-node
                c_node=H(:,j); 
                indexF=find(c_node==1);   % check the f_nodes that are connected to c(j)
                    
                for i=1:length(indexF) % for every f-node connected to c(j)
                    remIndF = indexF;
                    remIndF(i) = []; % remove contribution from f(i)

                    f_decision = c_next(j,remIndF); % the decision is made taking into account the responses of the f's connected to c(j), without f(j)
                    c_decision = [bit_rx(j) f_decision]; % add the decision of the received vector
                    nbVotes=length(c_decision); % majority voting
                    if(sum(c_decision) == nbVotes/2)
                        % if there are 50% vote for 1 and 50% for 0, the decision
                        % is the original received bit
                        decision(j,indexF(i)) = bit_rx(j);
                    elseif(sum(c_decision)>nbVotes/2)
                        % There are more than 50% 1s --> the decision is 1
                        decision(j, indexF(i)) =1;
                    else 
                        decision(j,indexF(i)) = 0;
                    end
                end
            end       
        end
    end
    % final_decision = decision;
    result = [result final_decision];
    %u = [u decision];
end