### A Pluto.jl notebook ###
# v0.18.1

using Markdown
using InteractiveUtils

# ╔═╡ fdb00d02-cb5b-49e0-aa5b-af108ab078f4
begin
	using PlutoUI
end

# ╔═╡ deb978d1-54a6-4d1d-9722-0a44778c4a94
md"""
### Discrete State Models for Transcriptional Regulation

In this lecture, we'll take a deeper dive into the modeling and analysis of transcription (TX) and translation (TL). In particular, we will:

* Revisit the mRNA and protein balance equations
* Introduce the discrete state promoter model formulation
* Explore gene expression dynamics in a population of single cells

"""

# ╔═╡ ffddbade-a162-4ebc-a81c-3ff672ae3a2c
md"""

### Revisit: mRNA and protein balance equations

A few lectures ago, we introduced balance equations for messenger RNA (mRNA) and protein abundance.
The material balance equations governing the abundance of the ith mRNA species ($m_{i}$) and the ith protein species ($p_{i}$) are given by (in cellmass specific units in a growing population):

$$\begin{eqnarray}
\frac{dm_{i}}{dt} & = & r_{X,i}u_{i} - (\theta_{m,i}+\mu)m_{i}+\lambda_{i}\\
\frac{dp_{i}}{dt} & = & r_{L,i}w_{i} - (\theta_{p,i}+\mu)p_{i}
\end{eqnarray}$$

The $r_{X, i}$ denote the _kinetic limit_ of transcription (if a gene were `ON`, this is the rate of production that you would expect), while the $r_{L, i}$ denote the _kinetic limit_ of translation. The second set of terms represents the degradation and dilution due to growth. The term $\lambda_{i}$ denotes the rate of unregulated or background expression (units: conc/time). Finally, the $u_{i}\in\left[0,1\right]$ and $w_{i}\in\left[0,1\right]$ terms describe the input controlling the transcription of gene $i$, and the translation of mRNA $i$, respectively. There is a __huge__ variety of approaches to formulate these control functions (many are super interesting)! 

However, these equations hold a dark secret; they actually assume a well-mixed population of cells. As well shall see, the well-mixed assumption is not valid (but does it matter?)

##### Derivation

Let's consider the integral form of the mol balance around the ith mRNA $m_{i}$ written with respect to some volume basis $\beta$ in a _batch_ system (no convective or diffusive transport terms):

$$\frac{d}{dt}\left(\int_{\beta}m_{i}d\beta\right) = \int_{\beta}\left(
\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right)d\beta\qquad{i=1,2,\dots,\mathcal{M}}$$

The term on the left-hand side is the accumulation term, while the right-hand side is the reaction term where $\sigma_{ij}$ denotes the stoichiometric coefficient for species $i$ in reaction $j$, 
$\mathcal{R}$ represents the number of reactions, and $\hat{r}_{j}$ denotes the kinetics of reaction $j$ per unit $\beta$. 

###### Well mixed assumption
If the system has no $\beta$ dependence, i.e., the system is well mixed such that there is no variation of $m_{i}$ with $\beta$, then we can pull $m_{i}$ out of the integral:

$$\frac{d}{dt}\left(\int_{\beta}m_{i}d\beta\right) \simeq \frac{d}{dt}\left(m_{i}\beta\right)$$

Similarly, we can pull the out the reaction terms to give:

$$\int_{\beta}\left(
\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right)d\beta \simeq \beta\left(
\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right)$$

###### Alternative explanation
Another way to think about these arguments is to assume there _is_ variation with $\beta$, but we approximate the integrals as:

$$\int_{\beta}m_{i}d\beta \simeq \langle m_{i} \rangle\beta$$

and

$$\int_{\beta}\left(
\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right)d\beta \simeq 
\langle \sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j} \rangle\beta$$

where $\langle\cdot\rangle$ denotes an averaged quantity over $\beta$. 

###### Putting it all together
Either argument gives a balance equation of the form (shown for the alternative perspective):

$$\frac{d}{dt}\left(\langle m_{i} \rangle\beta\right) = \langle \sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j} \rangle\beta\qquad{i=1,2,\dots,\mathcal{M}}$$

When we use cell mass specific units, we assume $\beta = XV_{R}$, where $X$ denotes the cell mass abundance (units: gDW/L) in the culture, while $V_{R}$ (units: L) represents the culture volume.

"""

# ╔═╡ 9aff53c6-0d8e-431f-b690-72c9b97be740
md"""
### Discrete state model for promoter functions 

Suppose a promoter $P$ can exits in one of $s=1,2\dots,\mathcal{S}$ possible _disjoint_ microstates, where each microstate $s$ has some pseudo energy $\epsilon_{s}$. Some microstates will lead to expression (the ability to produce an mRNA molecule by transcription), while others do not. For each microstate $s$, let's assign a pseudo energy $\epsilon_{s}$, where by definition $\epsilon_{1}=0$; we assume the base state has zero energy. Next, suppose the probability that promoter $P$ is in microstate $s$ follows a [Boltzmann distribution](https://en.wikipedia.org/wiki/Boltzmann_distribution) which says:

$$p_{i} = \frac{1}{Z} \times f_{i}\exp\left(-\beta\epsilon_{i}\right)\qquad{i=1,2,\dots,\mathcal{S}}$$

where $p_{i}$ denotes the probability that promoter $P$ is in microstate $i=1,2,\dots,\mathcal{S}$, $f_{i}$ denotes a state-specific factor $f_{i}\in\left[0,1\right]$, $\beta$ denotes the [thermodynamic beta](https://en.wikipedia.org/wiki/Thermodynamic_beta) and $Z$ denotes a normalization factor (called the [Partiton function](https://en.wikipedia.org/wiki/Partition_function_(statistical_mechanics)) in the statistical physics community). We can find $Z$ using the summation law of discrete probolity e.g.,  $\sum_{s}p_{s} = 1$ which gives:

$$Z = \sum_{s=1}^{\mathcal{S}}f_{i}\exp\left(-\beta\epsilon_{i}\right)$$

which gives:

$$p_{i} = \frac{f_{i}\exp\left(-\beta\epsilon_{i}\right)}{\displaystyle \sum_{s=1}^{\mathcal{S}}f_{i}\exp\left(-\beta\epsilon_{i}\right)}\qquad{i=1,2,\dots,\mathcal{S}}$$.

Finally, we relate the probability that promoter $P$ is in microstate $s$ back to the $u\left(\dots\right)$ control function by computing the overall probability that the desired event happens, e.g., promoter $P$ undergoes transcription. We know if $\Omega = \left\{1,2,\dots,\mathcal{S}\right\}$, then we can define the subset $\mathcal{A}\subseteq\Omega$ in which the desired event happens (in this case transcription). Given $\mathcal{A}$, the $u\left(\dots\right)$ function becomes:

$$u=\sum_{s\in{\mathcal{A}}}p_{s}$$

##### Example

* [Moon TS, Lou C, Tamsir A, Stanton BC, Voigt CA. Genetic programs constructed from layered logic gates in single cells. Nature. 2012;491(7423):249-253. doi:10.1038/nature11516](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC3904217/)
"""

# ╔═╡ 8494b59a-bcb4-4836-a3dd-5a72d4060d7b
md"""
### A microscopic view of gene expression

Typical transcription and translation models (such as those we'll use in this class) are population-averaged descriptions of gene expression dynamics. However, is this what is happening inside a single cell? For example, is there significant cell to cell variation in the _Escherichia coli_ cells we are using to make the product `XYZ`? Yes, there is. However, can we still use population-averaged models? It depends. 

##### Example:

* [Golding I, Paulsson J, Zawilski SM, Cox EC. Real-time kinetics of gene activity in individual bacteria. Cell. 2005 Dec 16;123(6):1025-36. doi: 10.1016/j.cell.2005.09.031. PMID: 16360033.](https://pubmed.ncbi.nlm.nih.gov/16360033/)
"""

# ╔═╡ ba43f178-01c8-4801-b750-32b84b923707
md"""
### Summary and conclusions

In this lecture we:

* Derived the mRNA and protein balance equations for a population of cells
* Introduced discrete state promoter models used them to capture averaged expression dynamics for logic circuits 
* Explored gene expression dynamics in a population of single cells
"""

# ╔═╡ f25cda13-83fd-470a-b755-ec7b4877b19d
md"""
### Next time
"""

# ╔═╡ 9a85e011-6188-41f0-b309-b0a869ec5e9d
TableOfContents(title="📚 Table of Contents", indent=true, depth=5, aside=true)

# ╔═╡ f5cafa0c-4173-4354-a873-0953e9bcfcb1
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

# ╔═╡ 248cbb68-9ec6-11ec-05b5-abe80f6c2b53
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
PlutoUI = "~0.7.37"
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
# ╟─deb978d1-54a6-4d1d-9722-0a44778c4a94
# ╠═ffddbade-a162-4ebc-a81c-3ff672ae3a2c
# ╟─9aff53c6-0d8e-431f-b690-72c9b97be740
# ╟─8494b59a-bcb4-4836-a3dd-5a72d4060d7b
# ╠═ba43f178-01c8-4801-b750-32b84b923707
# ╟─f25cda13-83fd-470a-b755-ec7b4877b19d
# ╠═9a85e011-6188-41f0-b309-b0a869ec5e9d
# ╠═fdb00d02-cb5b-49e0-aa5b-af108ab078f4
# ╠═f5cafa0c-4173-4354-a873-0953e9bcfcb1
# ╠═248cbb68-9ec6-11ec-05b5-abe80f6c2b53
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
