function transition(time::Int64, system::Array{Int64,2}, parameters::Dict{String,Any})::Array{Float64,2}
    
    # initialize -
    T = zeros(2,2)
    A = system[time,1]
    B = system[time,2]
        
    # get parameters -
    K1 = parameters["K1"]
    n1 = parameters["n1"]
    K2 = parameters["K2"]
    n2 = parameters["n2"]

    # compute -
    p1 = A^n1/(K1^n1+A^n1)
    p2 = B^n2/(K2^n2+B^n2)
    
    # package -
    T[1,1] = 1 - p1
    T[1,2] = p1
    T[2,1] = p2
    T[2,2] = 1 - p2
    
    # return -
    return T
end

function machine(time::Int64, state::Int64, system::Array{Int64,2}, parameters::Dict{String,Any})::Int64
       
    # compute the probabilities -
    T = transition(time, system, parameters)
    
    # generate a random number -
    r = rand() # uniform from 0 -> 1
    
    # what are my choice?
    choices = 1.0 .- T[state,:]
    
    # what is the index of the max?
    idx_max = argmax(choices)
    p = choices[idx_max]
    
    if (r >= p)
        return argmax(choices)
    else
        return argmin(choices)
    end
end