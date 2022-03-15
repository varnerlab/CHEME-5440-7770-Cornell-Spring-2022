### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ c835ff26-9df0-416c-a64f-9e04db57543b
md"""
### Example: Analysis of GFP expression using FBA
"""

# ╔═╡ d052ff3c-2a8d-4897-b3ce-aa08d06f8124
md"""
##### Setup: objective coefficient array
"""

# ╔═╡ d4c47ae1-d670-48f6-a2f1-61c8ff5620db
md"""
##### Setup: species bounds array
"""

# ╔═╡ 213041c0-2bfe-4c10-81da-5301b4d46b58
md"""
##### Setup: flux bounds array
"""

# ╔═╡ 02756d50-5881-49c8-be5a-d228b90b3912
md"""
##### Solve
"""

# ╔═╡ f12fdff4-0ce5-4d07-a278-16d4f41ede2f
md"""
__Table__: Open extent table $\dot{\epsilon}_{j}$ (units: $\star$mol/time) for the transcription and translation of green fluorescent protein (GFP).
"""

# ╔═╡ 6f852067-8c2e-45d0-baa9-73875ca49e12
md"""
__Table__: Input/output mol transfer rate $\dot{n}_{i,\star}$ (units: ⋆mol/time) for the transcription and translation of green fluorescent protein (GFP).
"""

# ╔═╡ 34d2f753-d74a-457a-83bc-44e08d6323ca
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

# ╔═╡ 656bc222-a43f-11ec-1f09-6b332c3b672f
begin

	# load some packages -
	using PlutoUI
	using BSON
	using DataFrames
	using TOML
	using GLPK
	using PrettyTables
	
	# setup paths -
	_PATH_TO_ROOT = pwd()
	_PATH_TO_DATA = joinpath(_PATH_TO_ROOT, "data")
	_PATH_TO_FIGS = joinpath(_PATH_TO_ROOT, "figs")
	_PATH_TO_SRC = joinpath(_PATH_TO_ROOT, "src")

	# load the 5440/7770 code lib -
    lib = ingredients(joinpath(_PATH_TO_SRC, "Include.jl"))

	# show -
	nothing
end

# ╔═╡ c7321b9e-b741-4d2d-8fd8-f3b8954d368f
begin

	# load the GFP model file as a dictionary -
	model_file_path = joinpath(_PATH_TO_DATA,"GFP.bson")
	model = BSON.load(model_file_path, @__MODULE__)[:system]

	# show - 
	nothing
end

# ╔═╡ e97e290a-d920-4797-8315-c65adf94c8d5
model

# ╔═╡ 18137921-25a7-4572-ba8c-39203031c24c
S = model["stoichiometric_matrix"]

# ╔═╡ 75394d6e-f185-4d3c-8b27-4d27c1a7c94a
(ℳ,ℛ) = size(S)

# ╔═╡ 0562d302-30dd-4c88-a9bb-5ec4100b3060
begin
	# setup the objective coefficient array -
	c_vector = model["objective_coefficient_array"]
	reaction_table = model["reaction_table"]
	
	# find the index of the reaction we want to maximize -
	target_reaction_id = "mRNA_EGFP_translation"
	idx_target = findfirst(x->x==target_reaction_id, reaction_table[!,:id])

	# update -
	if (isnothing(c_vector) == true)
		c_vector = zeros(ℛ)
	end
	c_vector[idx_target] = -1.0
	model["objective_coefficient_array"] = c_vector

	# show -
	nothing
end

# ╔═╡ 29a92990-1e95-49cd-aac0-7b68e1099a1d
model["species_table"][!,:symbol]

# ╔═╡ e4538ec6-678d-43dd-a01e-18a92dfe8e18
begin

	# get the default flux bounds array -
	flux_bounds_array = model["flux_bounds_array"]

	# show -
	nothing
end

# ╔═╡ 758c8e5e-7a41-40c6-945c-bf68e93947e9
aa_map = TOML.parsefile(joinpath(_PATH_TO_DATA,"AAMap.toml"))

# ╔═╡ a1f07b91-e911-4811-b02b-64ddc487e6af
begin

	# species_bounds_array 
	species_bounds_array = model["species_bounds_array"]
	
	# update the upper bound -
	species_bounds_array[:,2] .= Inf # everything can go out of the box at whatever rate -

	# update the lower bound -
	species_bounds_array[5:9,1] .= -1000.0

	# we need to update the bounds of all of the AA's -
	species_symbol_array = model["species_table"][!,:symbol]
	for (key, aa_symbol) ∈ aa_map

		# find index of aa_symbol -
		aa_index = findfirst(x->x==aa_symbol, species_symbol_array)

		# update the bound -
		species_bounds_array[aa_index,1] = -100.0
	end

	# show -
	nothing
end

# ╔═╡ fbd77bc8-9388-4e79-9036-bdca3f03fd4c
begin

	# estimate the optimal flux -
	result = lib.flux(S,flux_bounds_array,species_bounds_array,c_vector);
	
end

# ╔═╡ ea038677-c71d-4e64-992d-c36497fe52f8
with_terminal() do

	v = result.calculated_flux_array
	
	# let's build a table -
	state_array = Array{Any,2}(undef, ℛ, 5)
	for i ∈ 1:ℛ
		state_array[i,1] = i
		state_array[i,2] = v[i]
		state_array[i,3] = reaction_table[i,:id]
		state_array[i,4] = reaction_table[i,:forward_reaction]
		state_array[i,5] = reaction_table[i,:reverse_reaction]
	end

	# setup the header -
	header_data = (
		["i", "ϵ", "reaction name", "foward", "reverse"], 
		["", "⋆mol/time", "", "", ""]
	);
	
	pretty_table(state_array; alignment=:l, header=header_data)
end

# ╔═╡ ada9b431-b449-40c8-9a65-77a88259615e
with_terminal() do

	# get the "flux" -
	ϵ_dot = result.calculated_flux_array
	
	# get compute the 
	n_dot_in = abs.(species_bounds_array[:,1])
	n_dot_out = n_dot_in + S*ϵ_dot
	
	# setup species table -
	state_array = Array{Any, 2}(undef, ℳ, 5)
	for i ∈ 1:ℳ
		state_array[i,1] = i
		state_array[i,2] = species_symbol_array[i]
		state_array[i,3] = round(n_dot_in[i]; digits=2)
		state_array[i,4] = round(n_dot_out[i]; digits=2)
		state_array[i,5] = round((n_dot_out[i] - n_dot_in[i]); digits=2)
	end

	# setup the header -
	header_data = (
		["i", "species", "nᵢ_in", "nᵢ_out", "Δn"], 
		["", "", "⋆mol/time", "⋆mol/time", "⋆mol/time"]
	);
	
	pretty_table(state_array; header=header_data)
end

# ╔═╡ 446ef6c4-b3c5-4b17-9cec-098dac93774b
html"""
<style>
main {
    max-width: 900px;
    width: 70%;
    margin: auto;
    font-family: "Roboto, monospace";
}

a {
    color: blue;
    text-decoration: none;
}
</style>"""

# ╔═╡ cb4f47b8-e63a-4638-a832-38a1b38a329c
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
BSON = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GLPK = "60bf3e95-4087-53dc-ae20-288a0d20c6a6"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
TOML = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[compat]
BSON = "~0.3.5"
DataFrames = "~1.3.2"
GLPK = "~1.0.0"
PlutoUI = "~0.7.37"
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

[[deps.BSON]]
git-tree-sha1 = "306bb5574b0c1c56d7e1207581516c557d105cad"
uuid = "fbb218c0-5317-5bc6-957e-2ee96dd4b1f0"
version = "0.3.5"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

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

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "ae02104e835f219b8930c7664b8012c93475c340"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.2"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLPK]]
deps = ["GLPK_jll", "MathOptInterface"]
git-tree-sha1 = "7971e2ce3715a873b539174137bd8c4e19ac7a8f"
uuid = "60bf3e95-4087-53dc-ae20-288a0d20c6a6"
version = "1.0.0"

[[deps.GLPK_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "fe68622f32828aa92275895fdb324a85894a5b1b"
uuid = "e8aa6df9-e6ca-548a-97ff-1f85fc5b8b98"
version = "5.0.1+0"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"

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

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

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

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "JSON", "LinearAlgebra", "MutableArithmetics", "OrderedCollections", "Printf", "SparseArrays", "Test", "Unicode"]
git-tree-sha1 = "a62df301482a41cb7b1db095a4e6949ba7eb3349"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.1.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "ba8c0f8732a24facba709388c74ba99dcbfdda1e"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "db3a23166af8aebf4db5ef87ac5b00d36eb771e2"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

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

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

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
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

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
# ╠═c7321b9e-b741-4d2d-8fd8-f3b8954d368f
# ╟─c835ff26-9df0-416c-a64f-9e04db57543b
# ╠═e97e290a-d920-4797-8315-c65adf94c8d5
# ╠═18137921-25a7-4572-ba8c-39203031c24c
# ╠═75394d6e-f185-4d3c-8b27-4d27c1a7c94a
# ╟─d052ff3c-2a8d-4897-b3ce-aa08d06f8124
# ╠═0562d302-30dd-4c88-a9bb-5ec4100b3060
# ╟─d4c47ae1-d670-48f6-a2f1-61c8ff5620db
# ╠═29a92990-1e95-49cd-aac0-7b68e1099a1d
# ╠═758c8e5e-7a41-40c6-945c-bf68e93947e9
# ╠═a1f07b91-e911-4811-b02b-64ddc487e6af
# ╠═213041c0-2bfe-4c10-81da-5301b4d46b58
# ╠═e4538ec6-678d-43dd-a01e-18a92dfe8e18
# ╟─02756d50-5881-49c8-be5a-d228b90b3912
# ╠═fbd77bc8-9388-4e79-9036-bdca3f03fd4c
# ╠═f12fdff4-0ce5-4d07-a278-16d4f41ede2f
# ╠═ea038677-c71d-4e64-992d-c36497fe52f8
# ╟─6f852067-8c2e-45d0-baa9-73875ca49e12
# ╠═ada9b431-b449-40c8-9a65-77a88259615e
# ╠═656bc222-a43f-11ec-1f09-6b332c3b672f
# ╠═34d2f753-d74a-457a-83bc-44e08d6323ca
# ╠═446ef6c4-b3c5-4b17-9cec-098dac93774b
# ╠═cb4f47b8-e63a-4638-a832-38a1b38a329c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
