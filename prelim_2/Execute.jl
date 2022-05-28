include("Include.jl")


function main()

    # setup parameters -
    parameters = Dict{String,Any}()
    parameters["K1"] = 20
    parameters["n1"] = 4
    parameters["K2"] = 40
    parameters["n2"] = 8

    # build the Poisson -
    λ = 10
    d = Poisson(λ)

    # S -
    S = [-1 0 ; 0 1]

    # set the system dimensions -
    number_of_chemical_species = 2
    number_of_time_steps = 100
    number_of_sample_paths = 100

    # chemical state -
    A = 100
    B = 0
    chemical_state = Array{Int64,3}(undef, number_of_time_steps, number_of_chemical_species, number_of_sample_paths)

    # main simulation loop -
    for s ∈ 1:number_of_sample_paths
        
        # machine state -
        machine_state = Vector{Int64}(undef, number_of_time_steps)
        machine_state[1] = 1
        
        # setup the chemical state -
        chemical_state[1,1,s] = A
        chemical_state[1,2,s] = B
        
        for t ∈ 2:number_of_time_steps
        
            # run the machine -
            machine_state[t] = machine(t, machine_state[t-1], chemical_state[:,:,s], parameters)

            # what state are we in?
            current_machine_state = machine_state[t]
            if (current_machine_state == 2)

                # the enzyme has activity if machine in state 2 -
                Δ = ones(2)*rand(d)
                rV = S*Δ

                # do we have enough A?
                if ((chemical_state[t-1,1,s] + rV[1]) > 0)

                    # update the chemical state -
                    chemical_state[t,1,s] = chemical_state[t-1,1,s] + rV[1]
                    chemical_state[t,2,s] = chemical_state[t-1,2,s] + rV[2]
                    
                else
                    
                    # state stays the same -
                    chemical_state[t,1,s] = chemical_state[t-1,1,s]
                    chemical_state[t,2,s] = chemical_state[t-1,2,s]
                end
            else

                # state stays the same -
                chemical_state[t,1,s] = chemical_state[t-1,1,s]
                chemical_state[t,2,s] = chemical_state[t-1,2,s]
            end
        end 
    end

    return (chemical_state, number_of_sample_paths)
end

# run the sim -
(chemical_state, number_of_sample_paths) = main();

# make the plots -
for s ∈ 1:number_of_sample_paths
    
    if (s == 1)
        plot(chemical_state[:, 1,s], legend=false, c=:red)
        plot!(chemical_state[:,2,s], c=:black)
    else
        plot!(chemical_state[:,1,s], legend=false, c=:red)
        plot!(chemical_state[:,2,s], c=:black)
    end
end
xlabel!("Time step index", fontsize=18)
ylabel!("Molecule number", fontsize=18)
current()