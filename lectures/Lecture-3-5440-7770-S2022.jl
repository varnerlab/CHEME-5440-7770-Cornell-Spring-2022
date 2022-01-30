### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 9e2566fd-4fc9-4828-913f-f755c5e5cf50
idx_ep = 4

# ╔═╡ 907af081-4a01-4a51-bda9-3a34295f9208
function ingredients(path::String)
	
	# this is from the Julia source code (evalfile in base/loading.jl)
	# but with the modification that it returns the module instead of the last object
	name = Symbol("lib")
	m = Module(name)
	Core.eval(m,
        Expr(:toplevel,
             :(eval(x) = $(Expr(:core, :eval))($name, x)),
             :(include(x) = $(Expr(:top, :include))($name, x)),
             :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
             :(include($path))))
	m
end

# ╔═╡ c51f8430-8147-11ec-3501-498c9bf67b57
begin

	# import some packages -
	using PlutoUI
	using LinearAlgebra
	
	# setup paths -
	const _PATH_TO_NOTEBOOK = pwd()
	const _PATH_TO_DATA = joinpath(_PATH_TO_NOTEBOOK,"data")
	const _PATH_TO_FIGS = joinpath(_PATH_TO_NOTEBOOK,"figs")
	const _PATH_TO_SRC = joinpath(_PATH_TO_NOTEBOOK,"src")

	# load the class lib as LectureLib -
	lib = ingredients(joinpath(_PATH_TO_SRC, "Include.jl"));

	# return -
	nothing
end

# ╔═╡ be4b1854-1a14-4452-ad7e-22d614740a10
begin

	# Setup a collection of reaction strings -
	reaction_array = Array{String,1}()

	# Setup a collection of reaction strings -
	reaction_array = Array{String,1}()

	# encode the reactions -
	push!(reaction_array,"v₁,A+ATP,B+ADP,false")
	push!(reaction_array,"v₂,B+NAD,C+NADH,false")
	push!(reaction_array,"v₃,C+ADP,D+ATP,false")
	push!(reaction_array,"v₄,D+NADH,E+NAD,false")
	push!(reaction_array,"v₅,D+NAD,F+NADH,false")
	push!(reaction_array,"v₆,F+NADH,G+NAD,false") # change
	push!(reaction_array,"v₇,H,F,false") # change
	push!(reaction_array,"v₈,G+ATP,I+ADP,false") # change
	push!(reaction_array,"v₉,I+NADH,H+NAD,false") # change
	push!(reaction_array,"v₁₀,G,J,false")
	push!(reaction_array,"vₐ,∅,A,false")
	push!(reaction_array,"vₑ,E,∅,false")
	push!(reaction_array,"vₕ,H,∅,false") # change
	push!(reaction_array,"vⱼ,J,∅,false") # change
	push!(reaction_array,"vATP,ATP,∅,false")
	push!(reaction_array,"vADP,ADP,∅,false")
	push!(reaction_array,"vNADH,NADH,∅,false")
	push!(reaction_array,"vNAD,NAD,∅,false")
	
	# compute the stoichiometric matrix -
	(S, species_array, reaction_name_array) = lib.build_stoichiometric_matrix(reaction_array; 
		expand=true);

	# show -
	nothing
end

# ╔═╡ 8965f69d-3014-46b3-816d-7d6f7fd57adf
(ℳ,ℛ) = size(S)

# ╔═╡ f305c659-76d6-4a89-b231-80f5fbf1517a
S

# ╔═╡ 276a5591-bf84-4750-95ea-9fc48ea7bacb
function expa(stoichiometric_matrix::Array{Float64,2})

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
				    Z = setdiff(idx_zero_inner,idx_zero_outer);
				    if (isempty(Z))
					    OK_ARR[outer_row_index, inner_row_index] = 0.0;
					else
						OK_ARR[outer_row_index, inner_row_index] = 1.0;				
					end
                end

            end # inner 
        end # outer


		idx_rows_to_keep = Array{Int64,1}()
		for row_index ∈ 1:NRP
			TMP_OK_ROW = OK_ARR[row_index,:]
			N = length(findall(x->x==1.0,TMP_OK_ROW))
			if (N == NRP)
				push!(idx_rows_to_keep,row_index)
			end
		end

		# which rows should we keep?
		P_IND = P[idx_rows_to_keep,:]

		if (isempty(P_IND) == false)
			# Update T -
			T = vcat(TNEW,P_IND);
		end
    end # main -

	
    # return - 
    return T

end # function -

# ╔═╡ 8be1e489-7382-4315-8c4c-111abdead290
PM = expa(S)

# ╔═╡ c425a427-071c-45dc-a417-27a96f3da32e
N = PM[:,1:ℛ]

# ╔═╡ 95bb3270-8487-43ee-9795-766dd36ad092
N[3,:]

# ╔═╡ b25cf5a8-483a-418a-960c-ae6a9deedf8b
idxnz = findall(x->x!=0,N[idx_ep,:])

# ╔═╡ 57b93c34-6528-41d9-8af3-b5b7fce5d9b6
reaction_name_array[idxnz]

# ╔═╡ bbfeb47e-8bbd-4c4b-9d29-65cb99ca26c1
N[3,idxnz]

# ╔═╡ d201abce-202c-44f6-98a3-67e47c2a99f4
html"""
<style>
main {
    max-width: 900px;
    width: 75%;
    margin: auto;
    font-family: "Roboto, monospace";
}

a {
    color: blue;
    text-decoration: none;
}

.H1 {
    padding: 0px 30px;
}
</style>"""

# ╔═╡ 647315d2-e0d7-4a55-bab8-e8ff9d784b97
html"""
<script>
	// initialize -
	var section = 0;
	var subsection = 0;
	var subsubsection = 0;
	var headers = document.querySelectorAll('h3, h5, h6');
	
	// main loop -
	for (var i=0; i < headers.length; i++) {
	    
		var header = headers[i];
	    var text = header.innerText;
	    var original = header.getAttribute("text-original");
	    if (original === null) {
	        
			// Save original header text
	        header.setAttribute("text-original", text);
	    } else {
	        
			// Replace with original text before adding section number
	        text = header.getAttribute("text-original");
	    }
	
	    var numbering = "";
	    switch (header.tagName) {
	        case 'H3':
	            section += 1;
	            numbering = section + ".";
	            subsection = 0;
	            break;
	        case 'H5':
	            subsection += 1;
	            numbering = section + "." + subsection;
	            break;
			case 'H6':
	            subsubsection += 1;
	            numbering = section + "." + subsection + "." + subsubsection;
	            break;
	    }
		// update the header text 
		header.innerText = numbering + " " + text;
	};
</script>"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.32"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0b5cfbb704034b5b4c1869e36634438a047df065"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "ae6145ca68947569058866e443df69587acc1806"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.32"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═be4b1854-1a14-4452-ad7e-22d614740a10
# ╠═8965f69d-3014-46b3-816d-7d6f7fd57adf
# ╠═f305c659-76d6-4a89-b231-80f5fbf1517a
# ╠═8be1e489-7382-4315-8c4c-111abdead290
# ╠═c425a427-071c-45dc-a417-27a96f3da32e
# ╠═9e2566fd-4fc9-4828-913f-f755c5e5cf50
# ╠═95bb3270-8487-43ee-9795-766dd36ad092
# ╠═b25cf5a8-483a-418a-960c-ae6a9deedf8b
# ╠═bbfeb47e-8bbd-4c4b-9d29-65cb99ca26c1
# ╠═57b93c34-6528-41d9-8af3-b5b7fce5d9b6
# ╠═276a5591-bf84-4750-95ea-9fc48ea7bacb
# ╠═c51f8430-8147-11ec-3501-498c9bf67b57
# ╠═907af081-4a01-4a51-bda9-3a34295f9208
# ╠═d201abce-202c-44f6-98a3-67e47c2a99f4
# ╠═647315d2-e0d7-4a55-bab8-e8ff9d784b97
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
