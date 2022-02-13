### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 219d899e-8e54-437e-b20b-840e923ed0c7
begin

	# load -
	using PlutoUI
	
end

# ╔═╡ 83a9dd26-8ab0-11ec-05a7-97fba528df76
md"""
### Simple and Complex Models of Enzyme Kinetics (Part 2)
"""

# ╔═╡ 967500aa-664a-4bca-81b6-1ee95bec595a
md"""
### Kinetics of multiple substrate reactions

##### Ping–pong mechanisms
Enzymes with ping–pong mechanisms include [oxidoreductases](https://en.wikipedia.org/wiki/Oxidoreductase), [transferases](https://en.wikipedia.org/wiki/Transferase), and [serine proteases](https://en.wikipedia.org/wiki/Serine_protease) such as [trypsin](https://en.wikipedia.org/wiki/Trypsin), [chymotrypsin](https://en.wikipedia.org/wiki/Chymotrypsin) and several enzymes of the blood clotting cascade. 

##### Random-order mechanisms

##### Power-law kinetics and biochemical systems theory (BST)
Power-law kinetics are a flexible tool to describe multiple substrate kinetics. Suppose reaction $v_{i}$ depends upon $j=1,2,\dots,\mathcal{F}$ factors. These factors can be concentration e.g., substrates or products well as other type of data e.g., cate:

$$v_{i} = \alpha_{i}\prod_{j=1}^{\mathcal{F}}X_{j}^{f_{ij}}\qquad{i=1,2,\dots,\mathcal{R}}$$

where $\alpha_{j}$ denotes the rate constant for reaction $j$, $X_{j}$ denotes the abundance of
factor $j$ and $f_{ij}$ denotes the kinetic order of factor $j$ in reaction $j$. Power-law kinetics are a prominent feature of Biochemical Systems Theory (BST). Biochemical systems theory has been developed since the 1960s by Michael Savageau, Eberhard Voit, and others for the systems analysis of biochemical processes:

* [Atkinson MR, Savageau MA, Myers JT, Ninfa AJ. Development of genetic circuitry exhibiting toggle switch or oscillatory behavior in Escherichia coli. Cell. 2003 May 30;113(5):597-607. doi: 10.1016/s0092-8674(03)00346-5. PMID: 12787501.](https://pubmed.ncbi.nlm.nih.gov/12787501/)
* [Alvarez-Vasquez F, Sims KJ, Cowart LA, Okamoto Y, Voit EO, Hannun YA. Simulation and validation of modeled sphingolipid metabolism in Saccharomyces cerevisiae. Nature. 2005 Jan 27;433(7024):425-30. doi: 10.1038/nature03232. PMID: 15674294.](https://pubmed.ncbi.nlm.nih.gov/15674294/)
* [Goel G, Chou IC, Voit EO. Biological systems modeling and analysis: a biomolecular technique of the twenty-first century. J Biomol Tech. 2006;17(4):252-269.](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC2291792/?report=classic)

##### General multisubstrate kinetics
Suppose the irreversible rate $v_{i}$ is dependent upon susbtrates $S_{j},j=1,2,\dots,\mathcal{S}$, then the multiple saturation kinetic form is given by:

$$v_{i} = V_{max,i}\left[\frac{\prod_{j}\frac{S_{j}}{K_{j}}}{\prod_{j}\left(1+\frac{S_{j}}{K_{j}}\right) - 1}\right]\qquad{i=1,2,\dots,\mathcal{R}}$$

where $V_{max,i}$ denote the maximum reaction rate (units: concentration/time), $S_{j}$ denotes the
substrate concentration (units: concentration) and $K_{j}$ denotes the saturation constant for substrate $j$.

* [Liebermeister W, Klipp E. Bringing metabolic networks to life: convenience rate law and thermodynamic constraints. Theor Biol Med Model. 2006;3:41. Published 2006 Dec 15. doi:10.1186/1742-4682-3-41](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC1781438/)

"""

# ╔═╡ 7c136883-cff7-4104-a952-ffc5929b871f
md"""
### Discrete State Model (DSM) Enzyme Regulation model
Suppose we model the rate $v_{j}$ as the product of a kinetic limit (a simple model of the rate) and a correction term that accounts for the missing regulation:

$$v_{j} = r_{j}\theta\left(...\right)_{j}$$

where $v_{j}$ denotes the overall rate (units: $\mu$M/time), $r_{j}$ denotes the kinetic limit i.e., the maximum rate of conversion (units: $\mu$M/time) and 
$0\leq \theta\left(...\right)_{j}\leq 1$ (units: dimensionless) is a control function that describes the influence of effector molecules. 

Suppose an enzyme $E$ can exits in one of $s=1,2\dots,\mathcal{S}$ possible microstates, where each microstate $s$ has some pseudo energy $\epsilon_{s}$. Some microstates will lead to activity (the ability to carry out the chemical reactions, while others will not). For each microstate $s$, let's assign a pseudo energy $\epsilon_{s}$, where by definition $\epsilon_{1}=0$; we assume the base state has the lowest energy. Next, suppose the probability that enzyme $E$ is in microstate $s$ follows a [Boltzmann distribution](https://en.wikipedia.org/wiki/Boltzmann_distribution) which says:

$$p_{i} = \frac{1}{Z} \times f_{i}\exp\left(-\beta\epsilon_{i}\right)\qquad{i=1,2,\dots,\mathcal{S}}$$

where $p_{i}$ denotes the probability that enzyme $E$ is in microstate $i=1,2,\dots,\mathcal{S}$, $f_{i}$ denotes a state-specific factor $f_{i}\in\left[0,1\right]$, $\beta$ denotes the [thermodynamic beta](https://en.wikipedia.org/wiki/Thermodynamic_beta) and $Z$ denotes a normalization factor (called the [Partiton function](https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)) in the statistical physics community). We can find $Z$ using the summation law of discrete probolity e.g.,  $\sum_{s}p_{s} = 1$ which gives:

$$Z = \sum_{s=1}^{\mathcal{S}}f_{i}\exp\left(-\beta\epsilon_{i}\right)$$

which gives:

$$p_{i} = \frac{f_{i}\exp\left(-\beta\epsilon_{i}\right)}{\displaystyle \sum_{s=1}^{\mathcal{S}}f_{i}\exp\left(-\beta\epsilon_{i}\right)}\qquad{i=1,2,\dots,\mathcal{S}}$$.

Finally, we relate the probability that enzyme $E$ is in microstate $s$ back to the $\theta$ control function by computing the overall probability that the desired event happens, e.g., enzyme $E$ catalyzes the reaction of interests. We know if $\Omega = \left\{1,2,\dots,\mathcal{S}\right\}$, then we can define the subset $\mathcal{A}\subseteq\Omega$ in which the desired event happens. Given $\mathcal{A}$, the $\theta$ function becomes:

$$\theta=\sum_{s\in{\mathcal{A}}}p_{s}$$
"""

# ╔═╡ b722071c-4256-4fea-84c2-3a209822a1c1
md"""
##### DSM model redux
"""

# ╔═╡ aaa1cfb4-3ba2-43dc-9ed7-214664cb8029
md"""
##### Does the DSM describe experimental data?
"""

# ╔═╡ 3f4ae7c0-f425-49c9-9df7-7424eaca7aa2
md"""
### Kinetics and flux balance analysis calculations
"""

# ╔═╡ f7a88d98-7149-485f-ade1-b1d39f28daaa
md"""
### Summary and Conclusions
"""

# ╔═╡ 5391cecc-bbd5-44e4-9a92-2184b4252f20
md"""
### Next Time
"""

# ╔═╡ 518e00f4-16a2-4ca1-a365-d9ef17437d97
TableOfContents(title="📚 Table of Contents", indent=true, depth=5, aside=true)

# ╔═╡ c375b63a-c383-4476-84fc-4b3078de8c64
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

# ╔═╡ 9780c35f-f044-410c-a16d-a38378c6feed
html"""
<style>
main {
    max-width: 860px;
    width: 70%;
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
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.34"
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
git-tree-sha1 = "8979e9802b4ac3d58c503a20f2824ad67f9074dd"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.34"

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
# ╟─83a9dd26-8ab0-11ec-05a7-97fba528df76
# ╟─967500aa-664a-4bca-81b6-1ee95bec595a
# ╟─7c136883-cff7-4104-a952-ffc5929b871f
# ╟─b722071c-4256-4fea-84c2-3a209822a1c1
# ╟─aaa1cfb4-3ba2-43dc-9ed7-214664cb8029
# ╟─3f4ae7c0-f425-49c9-9df7-7424eaca7aa2
# ╟─f7a88d98-7149-485f-ade1-b1d39f28daaa
# ╟─5391cecc-bbd5-44e4-9a92-2184b4252f20
# ╠═518e00f4-16a2-4ca1-a365-d9ef17437d97
# ╠═219d899e-8e54-437e-b20b-840e923ed0c7
# ╠═c375b63a-c383-4476-84fc-4b3078de8c64
# ╠═9780c35f-f044-410c-a16d-a38378c6feed
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
