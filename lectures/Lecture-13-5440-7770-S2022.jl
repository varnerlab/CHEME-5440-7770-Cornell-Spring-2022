### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ c3530026-5a87-4447-88f5-fc744cb48ff2
begin
	
	using PlutoUI 
end

# ╔═╡ 2c895e14-e71a-435d-86d4-c1b101877b3c
md"""
### The End of the Beginning: Whole-Cell Metabolic Engineering

In the first lecture, we introduced metabolic engineering. Today, we'll close the loop at put together all we learned. Metabolic engineering is the practice of optimizing genetic and regulatory processes within cells to increase the cell's production of a desired small molecule or protein product of interest. 

Metabolic engineers manipulate the biochemical networks cells use to convert raw materials into molecules necessary for the cell's survival. Metabolic engineering specifically seeks to:

1. Mathematically model biochemical networks, calculate the yield (product divided by substrate) of valuable products and identify parts of the network that constrain the production of these products of interest. 
1. Use genetic engineering techniques to modify the biochemical network to relieve constraints limiting production. Metabolic engineers can then model the modified network to calculate the new product yield and identify new constraints (back to 1).

In this (the final lecture of Part I of the course), we'll:

* Incorporate the costs of transcription and translation into the flux balance analysis problem
* Introduce a workaround (often used in practice) to account for the metabolic cost of gene expression
* Close the loop: using a metabolic model and flux balance to optimize glycan production in _E.coli_. 

"""

# ╔═╡ ce5c8e9c-7652-4c3e-98a7-98e3567dc695
md"""
### Integrating the Cost and Logic of Gene Expression With Flux Balance Analysis

The [Allen and Palsson study](https://pubmed.ncbi.nlm.nih.gov/12453446/) gave us a roadmap of integrating the cost of gene expression with the metabolic operation of the cell. Further, [Vilkhovoy et al.](https://pubmed.ncbi.nlm.nih.gov/29944340/) showed how we could do this in a cell-free system. However, how do you account for gene expression in metabolic models at the whole-genome scale?

* [Thiele I, Jamshidi N, Fleming RM, Palsson BØ. Genome-scale reconstruction of Escherichia coli's transcriptional and translational machinery: a knowledge base, its mathematical formulation, and its functional characterization. PLoS Comput Biol. 2009;5(3):e1000312. doi:10.1371/journal.pcbi.1000312](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC2648898/?report=classic)
* [Lerman JA, Hyduke DR, Latif H, Portnoy VA, Lewis NE, Orth JD, Schrimpe-Rutledge AC, Smith RD, Adkins JN, Zengler K, Palsson BO. In silico method for modeling metabolism and gene product expression at genome scale. Nat Commun. 2012 Jul 3;3:929. doi: 10.1038/ncomms1928. PMID: 22760628; PMCID: PMC3827721.](https://pubmed.ncbi.nlm.nih.gov/22760628/)
"""

# ╔═╡ bfd134fe-b311-4768-a53b-2c5c4c1b121a
md"""
##### Aside: A question of timescales
One open question is whether we need to solve the gene expression and flux estimation problem simultaneously. The time scale of enzyme activity (and hence metabolic flux) is much faster than gene expression. Why? 

* The $k_{cat}$ for [phosphofructokinase-2 (pfk) in _E.coli_ is 9240 s$^{-1}$](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104955&ver=7&trm=kcat+pfk+in+E.+coli&org=). Thus, pfk has a characteristic time scale of $\tau\sim{k_{cat}}^{-1}$ or about 1$\times$10$^{-4}$ s.
* The $k_{cat}$ for RNAP polymerase is $k_{cat}$=$e_{x}L^{-1}$. If we assume [$e_{x}\sim{35}$ nt/s](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111871&ver=1&trm=elongation+rate+RNAP+E.coli+&org=) and average read length [$L\sim{924}$ nt](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111922&ver=3&trm=avarge+gene+length+in+E.+coli&org=), this gives $k_{cat}\sim$ 0.038 s$^{-1}$ or a timescale of $\tau\sim$26.4 s.

Thus, the time scale of metabolic flux (e.g., the timescale of enzyme activity) is approximately 5-orders of magnitude faster than gene expression. Therefore, gene expression seems constant  from the perspective of metabolic flux; hence, we treat these scales separately in traditional models. 
"""

# ╔═╡ 8041cb6c-b80f-4e81-8d1d-5a2e9942d809
md"""
##### Aside: How do you account for macromolecular synthesis without the E/ME formulation?
The E/ME formulation of Palsson and coworkers is detailed and requires more time investment than perhaps you're willing to give in a metabolic engineering application. How do you account for the energy and resource cost of producing enzymes (and all the other machinery) required to make a cell?

The traditional way metabolic engineers address this question is to treat cell mass formation as just another reaction, e.g., the growth reaction $\mu$: 

$$\left\{precursors\right\}~\rightarrow~cell$$

where each cellmass precursor has some stoichiometric coefficient $\sigma_{i,\mu}$ which describes how much of precursor $i$ is consumed (or produced) during cellmass formation. This formulation leads to species material balances of the form (assuming cellmass specific units, in a constant volume batch culture):

$$\frac{dx_{i}}{dt} = \sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j} - \left(x_{i}-\sigma_{i,\mu}\right)\mu\qquad{i=1,2,\dots,\mathcal{M}}$$

At steady state, for example in a flux balance analysis problem, these balances form the species constraints:

$$\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j} - \left(x_{i}-\sigma_{i,\mu}\right)\mu = 0\qquad{i=1,2,\dots,\mathcal{M}}$$

Palsson and coworkers gave an example cell mass reaction for _E.coli_ in:

* [Orth JD, Fleming RM, Palsson BØ. Reconstruction and Use of Microbial Metabolic Networks: the Core Escherichia coli Metabolic Model as an Educational Guide. EcoSal Plus. 2010 Sep;4(1). doi: 10.1128/ecosalplus.10.2.1. PMID: 26443778.](https://pubmed.ncbi.nlm.nih.gov/26443778/)

A recent paper on sampling uncertain biomass stoichiometric coefficients:

* [Dinh et al. (2022) Quantifying the propagation of parametric uncertainty on flux balance analysis. Metabolic Engineering, 69:26-39; https://doi.org/10.1016/j.ymben.2021.10.012](https://www.sciencedirect.com/science/article/pii/S1096717621001634)

"""

# ╔═╡ 06c3eafa-10b3-4a47-9a06-39f6abd29db3
md"""
### Using flux balance analysis to improve metabolic performance
Ultimately a metabolic engineer uses a model to improve the performance of a production host or cell-free system. There is a huge variety of design approaches to improve metabolic function; let's look at a recent one from Cornell:

* [Wayman JA, Glasscock C, Mansell TJ, DeLisa MP, Varner JD. Improving designer glycan production in Escherichia coli through model-guided metabolic engineering. Metab Eng Commun. 2019;9:e00088. Published 2019 Mar 29. doi:10.1016/j.mec.2019.e00088](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6454127/)

"""

# ╔═╡ b15b540b-8f13-4f21-b2da-510d34143b33
md"""
### Summary and conclusions

Today, we:

* Incorporated the costs of transcription and translation into the flux balance analysis problem
* Introduced a workaround (often used in practice) to account for the metabolic cost of gene expression
* Closed the loop: using a metabolic model and flux balance to optimize glycan production in _E.coli_. 
"""

# ╔═╡ a818128b-d642-49e0-8c7c-db36b37ea882
md"""
### Alas, if we just had some more time, it would have been cool to talk about:

##### Better constraints give better estimates of flux:

* [Buescher JM, Antoniewicz MR, Boros LG, et al. A roadmap for interpreting (13)C metabolite labeling patterns from cells. Curr Opin Biotechnol. 2015;34:189-201. doi:10.1016/j.copbio.2015.02.003](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC4552607/)

##### Flux balance analysis prediction of mutant response is not correct:

* [Segrè D, Vitkup D, Church GM. Analysis of optimality in natural and perturbed metabolic networks. Proc Natl Acad Sci U S A. 2002 Nov 12;99(23):15112-7. doi: 10.1073/pnas.232349399. Epub 2002 Nov 1. PMID: 12415116; PMCID: PMC137552.](https://pubmed.ncbi.nlm.nih.gov/12415116/)

##### Molecular crowding plays a role in metabolism (non-representative random sample):


* [Beg QK, Vazquez A, Ernst J, de Menezes MA, Bar-Joseph Z, Barabási AL, Oltvai ZN. Intracellular crowding defines the mode and sequence of substrate uptake by Escherichia coli and constrains its metabolic activity. Proc Natl Acad Sci U S A. 2007 Jul 31;104(31):12663-8. doi: 10.1073/pnas.0609845104. Epub 2007 Jul 24. PMID: 17652176; PMCID: PMC1937523.](https://pubmed.ncbi.nlm.nih.gov/17652176/)'
* [Conrado RJ, Mansell TJ, Varner JD, DeLisa MP. Stochastic reaction-diffusion simulation of enzyme compartmentalization reveals improved catalytic efficiency for a synthetic metabolic pathway. Metab Eng. 2007 Jul;9(4):355-63. doi: 10.1016/j.ymben.2007.05.002. Epub 2007 May 26. PMID: 17601761.](https://pubmed.ncbi.nlm.nih.gov/17601761/)
* [Klumpp S, Scott M, Pedersen S, Hwa T. Molecular crowding limits translation and cell growth. Proc Natl Acad Sci U S A. 2013 Oct 15;110(42):16754-9. doi: 10.1073/pnas.1310377110. Epub 2013 Sep 30. PMID: 24082144; PMCID: PMC3801028.](https://pubmed.ncbi.nlm.nih.gov/22760628/)

"""

# ╔═╡ 84d8a1c4-483c-4cc2-bdd0-7c88d26fd319
TableOfContents(title="📚 Lecture Outline", indent=true, depth=5, aside=true)

# ╔═╡ 7988263a-a5f5-11ec-3836-150105d04264
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

# ╔═╡ dee24b55-a5c5-4a47-9e1c-88ea0cf1cd61
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
# ╟─2c895e14-e71a-435d-86d4-c1b101877b3c
# ╟─ce5c8e9c-7652-4c3e-98a7-98e3567dc695
# ╟─bfd134fe-b311-4768-a53b-2c5c4c1b121a
# ╟─8041cb6c-b80f-4e81-8d1d-5a2e9942d809
# ╟─06c3eafa-10b3-4a47-9a06-39f6abd29db3
# ╟─b15b540b-8f13-4f21-b2da-510d34143b33
# ╟─a818128b-d642-49e0-8c7c-db36b37ea882
# ╠═84d8a1c4-483c-4cc2-bdd0-7c88d26fd319
# ╠═c3530026-5a87-4447-88f5-fc744cb48ff2
# ╟─7988263a-a5f5-11ec-3836-150105d04264
# ╟─dee24b55-a5c5-4a47-9e1c-88ea0cf1cd61
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
