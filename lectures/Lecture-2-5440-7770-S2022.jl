### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ fe33473f-084c-4d42-b37a-b3f2cb8ff1f0
begin

	# import some packages -
	using PlutoUI

	# setup paths -
	const _PATH_TO_NOTEBOOK = pwd()
	const _PATH_TO_DATA = joinpath(_PATH_TO_NOTEBOOK,"data")
	const _PATH_TO_FIGS = joinpath(_PATH_TO_NOTEBOOK,"figs")
	const _PATH_TO_SRC = joinpath(_PATH_TO_NOTEBOOK,"src")
	
	# return -
	nothing
end

# ╔═╡ 6627255a-b788-4aaa-bf8e-c98a39553dea
md"""
### Introduction to the Constraint-Based Perspective
"""

# ╔═╡ a5da8cc5-7ec8-4072-bf36-9e54e77c2a3f
md"""
### Dynamic species material balances in the continuum limit

Let's consider the general case in which we have $\mathcal{R}$ chemical reactions and $\mathcal{M}$ metabolites (chemical species) operating in some well-mixed physical (or logical) control volume $\Omega$. In the general continuum well-mixed limit, the species mole balance around each chemical species $i$ takes the form:

$$\frac{d}{dt}\left(C_{i}\beta\right) = \sum_{s=1}^{\mathcal{S}}d_{is}\dot{n}_{is} + \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right)\beta \qquad{i=1,2,\dots,\mathcal{M}}$$

where $C_{i}$ denotes the concentration (or number density) of species $i$ calculated with respect to a _system basis_ $\beta$. 
The first set of terms on the right-hand side denotes the rate of physical or _logical_ transport (units: mol or number per unit time) of chemical species $i$ into or from control volume $\Omega$, where $d_{is}$ denotes a direction parameter. The second set of terms denote the chemical reaction terms which describe the rate of production (consumption) of species $i$ by enzyme-catalyzed chemical reaction(s) $j$. The  net rate per unit $\beta$ is denoted by $\hat{r}_{j}$, while $\sigma_{ij}$ denotes the stoichiometric coefficient relating species $i$ with reaction $j$:

* A stoichiometric coefficient $\sigma_{ij}$ > 0 implies that metabolite $i$ is __produced__ by reaction $j$
* A stoichiometric coefficient $\sigma_{ij}$  = 0 implies that metabolite $i$ is __not connected__ to reaction $j$
* A stoichiometric coefficient $\sigma_{ij}$ < 0 implies that metabolite $i$ is __consumed__ by reaction $j$

"""

# ╔═╡ cf04d58a-230d-47ea-ba4a-072719cc9aa2
md"""
##### General intracellular species material balances

Let's use the general balance equations above, and develop a general model that we can use to write material balance equations for both extracellular and intracellular models. In particular, let's expand the accumulation term (chain rule):

$$\frac{d}{dt}\left(C_{i}\beta\right) = \beta\frac{dC_{i}}{dt}+C_{i}\frac{d\beta}{dt}\qquad{i=1,\dots,\mathcal{M}}$$

which can be substituted into the general balance equation to give the General continuum species material balance model:

$$\frac{dC_{i}}{dt} = \frac{1}{\beta}\left(\sum_{s=1}^{\mathcal{S}}d_{is}\dot{n}_{is}\right) + \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - \frac{C_{i}}{\beta}\frac{d\beta}{dt}\qquad{i=1,\dots,\mathcal{M}}$$

When building models of the _intracellular_ state of a cell, there are no physical transport terms (material is not flowing or from the cell). Intracellular balances in industrial metabolic engineering applications for liquid cultures are typically written in _biomass specific units_, e.g., 

$$\beta = XV_{R}$$

where $X$ denotes the _dry_ cell mass per liter in the culture (units: gDW/L) and $V_{R}$ denotes the volume of the culture. Substituting $\beta$ into the general case gives:

$$\frac{dC_{i}}{dt} = \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - 
\frac{C_{i}}{XV_{R}}\left(X\frac{dV_{R}}{dt}+V_{R}\frac{dX}{dt}\right)\qquad{i=1,\dots,\mathcal{M}}$$

Expanding the dilution terms gives:

$$\frac{dC_{i}}{dt} = \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - 
\frac{C_{i}}{V_{R}}\frac{dV_{R}}{dt}-\frac{C_{i}}{X}\frac{dX}{dt}\qquad{i=1,\dots,\mathcal{M}}$$

"""

# ╔═╡ 0bc7cc05-d2da-4b41-9c5a-06c21e3caf70
md"""
###### Batch culture 

In a batch culture, there is no liquid volume change ($dV_{R}/dt$ = 0). However, cell mass accumulates over time in the culture. Thus, the dilution terms become:

$$\frac{dC_{i}}{dt} = \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - {\mu}C_{i}\qquad{i=1,\dots,\mathcal{M}}$$

where $\mu$ denotes the _specific growth rate_ of the cells in the batch culture $\mu = X^{-1}\dot{X}$ (units: hr$^{-1}$). 
"""

# ╔═╡ d31d5e9c-a5d7-4f27-899c-bacfe350cd64
md"""
###### Continuous culture 

In continuous culture, there is also no volume change (volume flow rate = volume flow rate out). However, there are transport terms into and from the reactor for _extracellular_ variables. In particular, the cell mass balance equation for a continuous culture (assuming a sterile feed stream) is given by:

$$\frac{dX}{dt} = \left(\mu-D\right)X$$

where $\mu$ denotes the specific growth rate (units: hr$^{-1}$), and $D$ denotes the _dilution rate_ (units: hr$^{-1}$). Substituting the cell mass time rate of change into the general species material balance gives:

$$\frac{dC_{i}}{dt} = \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - 
C_{i}\left(\mu - D\right)\qquad{i=1,\dots,\mathcal{M}}$$


"""

# ╔═╡ 13d7ca39-f956-479e-b560-d6338c1d8b12
md"""
###### Fed-batch culture 

In continuous culture, there is also no volume change. However, there are transport terms into and from the reactor. 

"""

# ╔═╡ 84d7abb6-38ea-48b8-b598-e658d4c52544
md"""
### The constraint based perspective
"""

# ╔═╡ 38df7ed4-ad9c-4d72-9d77-3aa08f9eec12
md"""
### Formulation and properties of the stoichiometric matrix S
The stoichiometric matrix $\mathbf{S}$ is a $\mathcal{M}\times\mathcal{R}$ array holding the 
coefficients $\sigma_{ij}$; each row of $\mathbf{S}$ describes a metabolite, while each column represents a possible chemical reaction. Thus, the stoichiometric matrix is a _digital_ representation of the biochemical capabilities of a metabolic network. However, when formulating the stoichiometric matrix, there is always the question of what should be included to describe known phenomena. For example, should every possible reaction or just some _subset_ of all possible reactions be included in the stoichiometric matrix when exploring the features of a metabolic dataset, e.g., measurements of metabolites coming into and from the network? The practical significance of this question is two-fold:

* Fewer material balances can describe the behavior of a minimal network, and
* Fewer chemical reactions result in less uncertainty when estimating the kinetic rates $\hat{r}_{j}$.

To explore this question, let's look at an example from the paper:

* [Bordbar A, Nagarajan H, Lewis NE, et al. Minimal metabolic pathway structure is consistent with associated biomolecular interactions. Mol Syst Biol. 2014;10(7):737. Published 2014 Jul 1. doi:10.15252/msb.20145243](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299494/)
"""

# ╔═╡ d9abf0a6-c968-4d2d-afad-58694eb287c0
md"""
$(PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS,"Fig-ToyNetwork-MinSpan.png")))
"""

# ╔═╡ b941f57c-f1a2-4363-9a37-22cae5f2e825
begin

	# Setup a collection of reaction strings -
	reaction_array = Array{String,1}()

	# encode the reactions -
	push!(reaction_array,"v₁,A+ATP,B+ADP,false")
	push!(reaction_array,"v₂,B+NAD,C+NADH,false")
	push!(reaction_array,"v₃,C+ADP,D+ATP,false")
	push!(reaction_array,"v₄,D+NADH,E+NAD,false")
	push!(reaction_array,"v₅,D+NAD,F+NADH,false")
	push!(reaction_array,"v₆,F+NADH,G+NAD,true")
	push!(reaction_array,"v₇,F,H,true")
	push!(reaction_array,"v₈,G+ATP,I+ADP,true")
	push!(reaction_array,"v₉,H+NAD,I+NADH,true")
	push!(reaction_array,"v₁₀,G,J,true")
	push!(reaction_array,"vₐ,Aₓ,A,false")
	push!(reaction_array,"vₑ,E,Eₓ,false")
	push!(reaction_array,"vₕ,H,Hₓ,true")
	push!(reaction_array,"vⱼ,J,Jₓ,true")
	push!(reaction_array,"vATP,ATP,ATPₓ,true")
	push!(reaction_array,"vADP,ADP,ADPₓ,true")
	push!(reaction_array,"vNADH,NADH,NADHₓ,true")
	push!(reaction_array,"vNAD,NAD,NADₓ,true")

	
	
end

# ╔═╡ 000109db-2cf2-439f-8803-306bb92b18c2
split("H+NAD",'+')

# ╔═╡ 2840b124-c1bd-4abf-ae37-5cc0b909a75c


# ╔═╡ 84c245b4-21b9-430b-8859-7a348c003141
function extract_species_dictionary(reaction_phrase::String;
	direction::Float64 = -1.0)::Dict{String,Float64}

	# initialize -
	species_symbol_dictionary = Dict{String,Float64}()
	
	# ok, do we hve a +?
	component_array = split(reaction_phrase,'+');
	for component ∈ component_array

		if (contains(component,'*') == true)
			
			tmp_array = split(component,'*')
			st_coeff = direction*parse(Float64,tmp_array[1])
			species_symbol = String(tmp_array[2])
			species_symbol_dictionary[species_symbol] = st_coeff
		end
		
		# strip any spaces -
		species_symbol = component |> lstrip |> rstrip
		species_symbol_dictionary[species_symbol] = direction*1.0
	end

	# return -
	return species_symbol_dictionary
end

# ╔═╡ 7fd25bfd-bea0-4656-a51a-6f10a1850287
function build_stoichiometric_matrix(reactions::Array{String,1})

	# initialize -
	species_array = Array{String,1}()
	reaction_dictionary_array = Array{Dict{String,Float64},1}()
	number_of_reactions = length(reactions)
	
	# first: discover the species list -
	for reaction_string ∈ reactions

		# initialize tmp storage -
		tmp_dictionary = Dict{String,Float64}()
		
		# split the reaction into its components -
		component_array = split(reaction_string,',');

		# reactant phrase => 2, and product phrase => 3
		reactant_phrase = String.(component_array[2]);
		product_phrase = String.(component_array[3]);

		# generate species lists for the reactants and products, then merge -
		tmp_dict_reactant = extract_species_dictionary(reactant_phrase; direction = -1.0);
		tmp_dict_product = extract_species_dictionary(product_phrase; direction = 1.0);
		merge!(tmp_dictionary, tmp_dict_reactant)
		merge!(tmp_dictionary, tmp_dict_product)

		# grab the tmp_dictionary for later -
		push!(reaction_dictionary_array, tmp_dictionary)

		# the species that we need to look at are the keys of the tmp_dictionary -
		tmp_species_list = keys(tmp_dictionary)
		
		# we need a unique species list, so check to see if we have already discovered this species -
		for tmp_species ∈ tmp_species_list

			if (in(tmp_species, species_array) == false)

				# ok, we have *not* seen this species before, let's grab it -
				push!(species_array, tmp_species)
			end
		end
	end

	# sort alphabetically -
	sort!(species_array)
	
	# now that we have a *unique* species array, let's initialize some storage for the stoichiometric array
	S = zeros(length(species_array), number_of_reactions);

	# last: fill in the values for stoichiometric coefficents -
	for (row_index, species_symbol) ∈ enumerate(species_array)
		for (col_index, reaction_dictionary) ∈ enumerate(reaction_dictionary_array)

			# ok: is this species symbol in my reaction dictionary?
			if (haskey(reaction_dictionary, species_symbol) == true)
				S[row_index,col_index] = reaction_dictionary[species_symbol]
			end
		end
	end

	# return -
	return S
end

# ╔═╡ 6495c567-7756-48cf-8c70-797c7e00f420
tmp = build_stoichiometric_matrix(reaction_array)

# ╔═╡ 8862ce06-cf11-4654-a471-6aa135a86917


# ╔═╡ b1b251d8-7e23-11ec-09d1-97ace15b3bd3
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

# ╔═╡ 9833473a-3fc3-4b61-bb07-050d6e159cb2
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
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
git-tree-sha1 = "92f91ba9e5941fc781fecf5494ac1da87bdac775"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.0"

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
# ╟─6627255a-b788-4aaa-bf8e-c98a39553dea
# ╟─a5da8cc5-7ec8-4072-bf36-9e54e77c2a3f
# ╟─cf04d58a-230d-47ea-ba4a-072719cc9aa2
# ╟─0bc7cc05-d2da-4b41-9c5a-06c21e3caf70
# ╟─d31d5e9c-a5d7-4f27-899c-bacfe350cd64
# ╟─13d7ca39-f956-479e-b560-d6338c1d8b12
# ╟─84d7abb6-38ea-48b8-b598-e658d4c52544
# ╟─38df7ed4-ad9c-4d72-9d77-3aa08f9eec12
# ╟─d9abf0a6-c968-4d2d-afad-58694eb287c0
# ╠═b941f57c-f1a2-4363-9a37-22cae5f2e825
# ╠═6495c567-7756-48cf-8c70-797c7e00f420
# ╠═000109db-2cf2-439f-8803-306bb92b18c2
# ╠═7fd25bfd-bea0-4656-a51a-6f10a1850287
# ╠═2840b124-c1bd-4abf-ae37-5cc0b909a75c
# ╠═84c245b4-21b9-430b-8859-7a348c003141
# ╠═8862ce06-cf11-4654-a471-6aa135a86917
# ╠═fe33473f-084c-4d42-b37a-b3f2cb8ff1f0
# ╠═b1b251d8-7e23-11ec-09d1-97ace15b3bd3
# ╟─9833473a-3fc3-4b61-bb07-050d6e159cb2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
