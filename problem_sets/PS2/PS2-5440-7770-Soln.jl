### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 3d4c79d0-10ec-4ba1-be3e-ce6cd8f83e9e
md"""
### CHEME 5440/7770: Problem Set 2 (PS2) solution
"""

# ╔═╡ 594f9194-91ef-4a41-999b-a553cca2c180
md"""
##### Build the Stoichiometric Matrix
"""

# ╔═╡ d93fad69-7f2f-420c-ba78-0133fc6e1572
md"""
##### Compute the extreme pathway matrix using expa
"""

# ╔═╡ 82ddb8c4-d286-47e4-9aa4-02f95b6792c9
idx_pathway = 13

# ╔═╡ 6389894e-e9b4-4530-b10c-feedb44862f9
md"""
##### Compute binary stoichiometric array
"""

# ╔═╡ 1e1670b6-f87d-463f-bd38-dbb2e6103653
md"""
##### Metabolite Connectivity Array (MCA)
"""

# ╔═╡ 641b3e7d-e0c8-44f8-bbcd-044bff755046
md"""
##### Reaction Connectivity Array (RCA)
"""

# ╔═╡ 5115eed1-546a-4d87-bfb3-d58ecebe653a
md"""
### Is there a relationship between reaction connectivity and usage?
"""

# ╔═╡ d0f77fde-88d5-11ec-3620-fb46e6aebbca
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

# ╔═╡ f7e07c69-f8b3-4956-90bc-bca6dbd41ac5
begin
	# import some packages -
	using PlutoUI
	using PrettyTables
	using LinearAlgebra
	
	# setup paths -
	const _PATH_TO_NOTEBOOK = pwd()
	const _PATH_TO_SRC = joinpath(_PATH_TO_NOTEBOOK,"src")
	const _PATH_TO_CONFIG = joinpath(_PATH_TO_NOTEBOOK,"config")
	const _PATH_TO_FIGS = joinpath(_PATH_TO_NOTEBOOK, "figs")

	# load the PS2 code lib -
	lib = ingredients(joinpath(_PATH_TO_SRC, "Include.jl"));

	# return -
	nothing
end

# ╔═╡ a7cec06e-8a50-4c34-81a3-0f303099cfd5
PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS,"Fig-Urea-cycle.png"))

# ╔═╡ 32f498cf-a379-4178-9ea2-b022c571227d
begin

	# load/parse the network file -
	_PATH_TO_NETWORK_FILE = joinpath(_PATH_TO_CONFIG,"Network_v2.net")
	list_of_reactions = lib.read_reaction_file(_PATH_TO_NETWORK_FILE);

	# build the stoichiometric array -
	(S, mna, rna) = lib.build_stoichiometric_matrix(list_of_reactions; expand=true);

	# get expanded reaction strings -
	expanded_reaction_array = lib.expand_reversible_reactions(list_of_reactions)

	# return -
	nothing
end

# ╔═╡ db048167-7ca6-4845-889c-6960b61a77dd
(ℳ,ℛ) = size(S)

# ╔═╡ 902d9a8b-2f11-425a-91cd-ccb2345503f1
begin

	# compute the extreme pathways -
	T = lib.expa(S)
	P = T[:,1:ℛ]
	C = T[:,(ℛ + 1):end]

	# show -
	nothing
end

# ╔═╡ 1ac2f945-868c-4262-b903-a30db71ce21d
P

# ╔═╡ b95e014c-e486-4608-a645-b62bb5ac2fb6
begin

	# compute the frequency of reaction usage in P -
	(RP,CP) = size(P)
	reaction_usage_array = Array{Int64,1}(undef,CP)
	for j ∈ 1:CP
		usage = 0
		for i ∈ 1:RP
			if (P[i,j] != 0.0)
				usage += 1 
			end
		end
		reaction_usage_array[j] = usage
	end
	
end

# ╔═╡ a5fb1e84-dc81-49b3-aa60-552e2bebc547
let

	# visualize pathways -
	idx_nz = findall(x->x!=0.0,P[idx_pathway,:])

	# initialize storage for the table -
	state_array = Array{Any,2}(undef, length(idx_nz),3)
	for (i,index) ∈ enumerate(idx_nz)
		state_array[i,1] = i
		state_array[i,2] = P[idx_pathway,index]
		state_array[i,3] = expanded_reaction_array[index]
	end

	# setup the header -
	header_data = (["index","coeff","reaction"])

	with_terminal() do
		pretty_table(state_array; header=header_data, alignment=:l)
	end
		
end

# ╔═╡ 84c09869-1893-450c-8ff7-e5033077a48c
B = lib.binary_stoichiometric_matrix(S)

# ╔═╡ 9786e10f-b1eb-43f4-8e9f-cce5b8f3af1f
# compute the metabolite connectivity array -
MCA = B*transpose(B)

# ╔═╡ e42e400f-0447-40d5-a0bd-7f86f4a1abfa
begin

	# build a ranking table -
	mca_diagonal = diag(MCA)

	# sort -
	idx_sort = sortperm(mca_diagonal; rev=true)

	# build the mca table -
	state_array = Array{Any,2}(undef,ℳ,3)
	for i ∈ 1:ℳ
		state_array[i,1] = i
		state_array[i,2] = mna[idx_sort[i]]
		state_array[i,3] = mca_diagonal[idx_sort[i]]
	end

	# setup header -
	header_data = (["rank","metabolite","connections"])
	
	with_terminal() do
		pretty_table(state_array; header = header_data)
	end
end

# ╔═╡ ade866b0-b419-45a5-83f7-4a0058adcc55
# compute the reaction connectivity array -
RCA = transpose(B)*B

# ╔═╡ 273d535a-c28e-451f-9a03-5b5397d13798
let
	# build a ranking table -
	rca_diagonal = diag(RCA)

	# sort -
	idx_sort = sortperm(rca_diagonal; rev=true)

	# build the mca table -
	state_array = Array{Any,2}(undef,ℳ,3)
	for i ∈ 1:ℳ
		state_array[i,1] = i
		state_array[i,2] = rna[idx_sort[i]]
		state_array[i,3] = rca_diagonal[idx_sort[i]]
	end

	# setup header -
	header_data = (["rank","reaction","connections"])
	
	with_terminal() do
		pretty_table(state_array; header = header_data)
	end
end

# ╔═╡ 250b0ed2-e465-4e4c-b2eb-4d5ae80387cc
let

	# build usage table -
	state_table = Array{Any,2}(undef, CP, 5)
	idx_sort = sortperm(reaction_usage_array; rev=true)
	for i ∈ 1:CP
		state_table[i,1] = i
		state_table[i,2] = idx_sort[i]
		state_table[i,3] = reaction_usage_array[idx_sort[i]]
		state_table[i,4] = (1/RP)*reaction_usage_array[idx_sort[i]]
		state_table[i,5] = expanded_reaction_array[idx_sort[i]]
	end

	# setup header -
	header_data = (["rank","index","usage","freq", "reaction"])
	
	with_terminal() do
		pretty_table(state_table; header=header_data, alignment=:l)
	end
end

# ╔═╡ b0f24d05-af37-477d-a709-b011cba825fc
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
</style>"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"

[compat]
PlutoUI = "~0.7.34"
PrettyTables = "~1.3.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
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

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

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

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

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

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

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
# ╟─3d4c79d0-10ec-4ba1-be3e-ce6cd8f83e9e
# ╠═a7cec06e-8a50-4c34-81a3-0f303099cfd5
# ╟─594f9194-91ef-4a41-999b-a553cca2c180
# ╠═32f498cf-a379-4178-9ea2-b022c571227d
# ╠═db048167-7ca6-4845-889c-6960b61a77dd
# ╟─d93fad69-7f2f-420c-ba78-0133fc6e1572
# ╠═902d9a8b-2f11-425a-91cd-ccb2345503f1
# ╠═1ac2f945-868c-4262-b903-a30db71ce21d
# ╠═82ddb8c4-d286-47e4-9aa4-02f95b6792c9
# ╠═a5fb1e84-dc81-49b3-aa60-552e2bebc547
# ╟─6389894e-e9b4-4530-b10c-feedb44862f9
# ╠═84c09869-1893-450c-8ff7-e5033077a48c
# ╟─1e1670b6-f87d-463f-bd38-dbb2e6103653
# ╠═9786e10f-b1eb-43f4-8e9f-cce5b8f3af1f
# ╟─e42e400f-0447-40d5-a0bd-7f86f4a1abfa
# ╟─641b3e7d-e0c8-44f8-bbcd-044bff755046
# ╠═ade866b0-b419-45a5-83f7-4a0058adcc55
# ╟─273d535a-c28e-451f-9a03-5b5397d13798
# ╟─5115eed1-546a-4d87-bfb3-d58ecebe653a
# ╠═b95e014c-e486-4608-a645-b62bb5ac2fb6
# ╠═250b0ed2-e465-4e4c-b2eb-4d5ae80387cc
# ╠═f7e07c69-f8b3-4956-90bc-bca6dbd41ac5
# ╠═d0f77fde-88d5-11ec-3620-fb46e6aebbca
# ╠═b0f24d05-af37-477d-a709-b011cba825fc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
