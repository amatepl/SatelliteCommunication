function corrected_bit = tanner(piece_bit_rx, H, maxit)
% Hard decoding of LDPC with Tanner diagram
% INPUTS :
% - piece_bit_rx : frame of bit received 
% - H            : parity check matrix reconstructed with permutation
% - maxit        : max number of iteration
% OUTPUT :
% - corrected_bit: LDPC decoded bit

[row col] = size(H);
ci = piece_bit_rx';
% Initialization
r = zeros(row, col);
% Asscociate the ci matrix with non-zero elements of H
c = H.*repmat(ci, row, 1);
 
% Iteration
for n = 1:maxit
   % ----- First step -----
   for i = 1:row
      % Find non-zeros in the column to know where is the connection edges
      edges = find(H(i, :));
      % sum of all ci associated of non-zero element of H with the ci and 
      % take te modulo 2 (even = 0 or odd = 1). 
      for k = 1:length(edges)
         r(i, edges(k)) = mod(sum(c(i, edges)) + c(i, edges(k)), 2);
      end 
   end
   
   % ------ Second step ------
   for j = 1:col
      % Find non-zeros in the rows to compare the bits received at one node
      comp = find(H(:, j));
      
      % Number of 1s in a row
      numOfOnes = length(find(r(comp, j)));
      
      for k = 1:length(comp)        
         % c is corrected if there are more ones than 0 with the comparison
         % of all bits received at one node
         if numOfOnes + c(j) >= length(comp) - numOfOnes + r(comp(k), j)
            c(comp(k), j) = 1;
         else
            c(comp(k), j) = 0;
         end    
      end
      
      % Bit decoding
      if numOfOnes + c(j) >= length(comp) - numOfOnes
         corrected_bit(j) = 1;
      else
         corrected_bit(j) = 0;
      end
             
   end 
   
end