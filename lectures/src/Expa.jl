function expa(stoichiometric_matrix::Array{Float64,2})::Array{Float64,2}

    # initialize -
    (ℳ,ℛ) = size(stoichiometric_matrix) 
    T = [Matrix{Float64}(I,ℛ,ℛ) transpose(stoichiometric_matrix)];

    # what is the dimension of the tables?
    (NRT,NCT) = size(T)

    # outer loop -
    for col_index = 1:ℳ

        # adjust the col index (from the org matrix -)
	    local_col_index = col_index + ℛ;
        
        # get the current RHS -
        RHS = T[:,(ℛ + 1):end] # this gets the cols belonging to the stoichiometric matrix -

        # Pick a col -
	    TMP_COL = RHS[:,col_index]; # this is grabbing a metabolite -

        # how many zeros, + and - in this col -
        idx_zero_array = findall(x->x==0.0,TMP_COL)
        idx_positive_array = findall(x->x>0.0,TMP_COL)
        idx_negative_array = findall(x->x<0.0,TMP_COL)

        # make new T -
        TNEW = T[idx_zero_array,:]

        # zero out non-zero elements -
        TMP_ARR = Array{Array{Float64,1},1}()
        for idx_positive ∈ idx_positive_array
            for idx_negative ∈ idx_negative_array
            
                # setup α and β -
			    α = abs(T[idx_positive,local_col_index]);
			    β = abs(T[idx_negative,local_col_index]); 

                # compute a tmp row -
                TMP_ROW = α*T[idx_negative,:] .+ β*T[idx_positive,:];

                # grab the tmp row -
                push!(TMP_ARR,TMP_ROW)
            end
        end
    
        # ok, so we have a new set of possible pathways, check for independent pathways -
        P = transpose(hcat(TMP_ARR...))
		(NRP,_) = size(P)
        OK_ARR = Matrix{Float64}(I,NRP,NRP);
        for outer_row_index ∈ 1:NRP
            for inner_row_index ∈ 1:NRP
                
                if (outer_row_index != inner_row_index)
                    
                    # grab the index of zeros for the outer and inner row -
                    idx_zero_outer = findall(x->x==0.0,P[outer_row_index,:])
                    idx_zero_inner = findall(x->x==0.0,P[inner_row_index,:])

                    # Ok, so now what?
				    Z = setdiff(idx_zero_outer,idx_zero_inner);
				    if (isempty(Z))
					    OK_ARR[outer_row_index, inner_row_index] = 0.0;
				    end
                end

            end # inner 
        end # outer

		# which rows should we keep?
		idx_rows_to_keep = findall(x->x==1.0,diag(OK_ARR))
		P_IND = P[idx_rows_to_keep,:]

		if (isempty(P_IND) == false)
			# Update T -
			T = vcat(TNEW,P_IND);
		end
    end # main -

    # return - 
    return T

end # function -