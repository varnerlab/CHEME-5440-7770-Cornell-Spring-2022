# include packages and other codes -
include("Include.jl")

# make an executable script -
function main()

    # ask the user for thier net id -
    print("What's your netid? ")

    # Calling rdeadline() function
    netid = readline()
    uuid_value = uuid4()

    # create a message, then display it -
    message = "Hello World from $(netid)! Your secret (random) 5440/7770 code is: $(uuid_value)"
    println(message)
end

# execute - 
main()