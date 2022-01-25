### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ a3101ae7-2699-43e7-aa28-332b0447bbc5
begin

	# Setup the Julia environment -
	using PlutoUI
	
	# setup paths -
	const _PATH_TO_NOTEBOOK = pwd()
	const _PATH_TO_DATA = joinpath(_PATH_TO_NOTEBOOK,"data")
	const _PATH_TO_FIGS = joinpath(_PATH_TO_NOTEBOOK,"figs")
	const _PATH_TO_SRC = joinpath(_PATH_TO_NOTEBOOK,"src")
	
	# return -
	nothing
end

# ╔═╡ 45e60528-7b8d-11ec-3aa3-65625ebbbfb2
md"""
### Introduction to Metabolic Engineering

Metabolic engineering is the practice of optimizing genetic and regulatory processes within cells to increase the cell's production of a desired small molecule or protein product of interest. 

Metabolic engineers manipulate the biochemical networks used by cells to convert raw materials into molecules necessary for the cell's survival. Metabolic engineering specifically seeks to:

1. Mathematically model biochemical networks, calculate the yield (product divided substrate) of useful products and identify parts of the network that constrain the production of the products of interest. 
1. Use genetic engineering techniques to modify the biochemical network in order to relieve constraints limiting production. The modified network can then be modeled to calculate the new product yield, and to identify new constraints (back to 1).

Resources for metabolic engineering:

* [Bailey JE. Toward a science of metabolic engineering. Science. 1991 Jun 21;252(5013):1668-75. doi: 10.1126/science.2047876. PMID: 2047876.](https://pubmed.ncbi.nlm.nih.gov/2047876/)

* [Stephanopoulos G, Vallino JJ. Network rigidity and metabolic engineering in metabolite overproduction. Science. 1991 Jun 21;252(5013):1675-81. doi: 10.1126/science.1904627. PMID: 1904627.](https://pubmed.ncbi.nlm.nih.gov/1904627/)

Resources for biochemical network information

* [Kanehisa M, Goto S. KEGG: kyoto encyclopedia of genes and genomes. Nucleic Acids Res. 2000 Jan 1;28(1):27-30. doi: 10.1093/nar/28.1.27. PMID: 10592173; PMCID: PMC102409.](https://www.genome.jp/kegg/)

* [Karp, Peter D et al. “The BioCyc collection of microbial genomes and metabolic pathways.” Briefings in bioinformatics vol. 20,4 (2019): 1085-1093. doi:10.1093/bib/bbx085](https://pubmed.ncbi.nlm.nih.gov/29447345/)

* [Gama-Castro, Socorro et al. “RegulonDB version 9.0: high-level integration of gene regulation, coexpression, motif clustering and beyond.” Nucleic acids research vol. 44,D1 (2016): D133-43. doi:10.1093/nar/gkv1156](https://pubmed.ncbi.nlm.nih.gov/26527724/)
"""

# ╔═╡ 96d697e1-d190-4d9b-aeb8-53466cc1cce2
md"""
__Fig 1.__ The overall metabolic map from the KEGG database. Each dot (_node_) is a metabolite, each line (_edge_) is a metabolic reaction. 
"""

# ╔═╡ 87c121a0-91ed-436f-a981-a8a5198c3226
PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS,"KEGG-map01100.png"))

# ╔═╡ 2e67cdc9-9cd4-4a20-9ac7-1ecfd4b46967
md"""
### Example biological products

Biotechnology can be divided into (roughly) two categories, industrial biotechnology, and medical biotechnology. Metabolic engineering plays a key role in both sectors. Industrial biotechnology is typically consumer-focused e.g., components of consumer products such as detergents, food products, and small-molecule chemical feedstocks. On the other hand, medical biotechnology develops molecules for human (and animal) health applications, e.g., antibodies, therapeutic proteins, vaccines, etc. 

"""

# ╔═╡ 21f2c20b-0386-4532-a81d-18b05df6c3cf
PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS,"Figs-BiotechMarket-2020.png"))

# ╔═╡ 226b8bdd-13aa-4ab6-b4cb-73dd9cba5cb2
md"""

##### Monoclonal Antibodies (mAbs) and therapeutic proteins 
[Monoclonal antibodies (mAbs)](https://en.wikipedia.org/wiki/Monoclonal_antibody) are important molecules for human health e.g., cancer treatments such as [Herceptin](https://www.herceptin.com) or everyday laboratory uses such as affinity reagents used [Western blotting](https://www.nature.com/scitable/definition/western-blot-288/). In addition to mAbs, there are a huge variety of therapeutic proteins e.g., clotting factors or recombinant human insulin [Humulin R](https://www.accessdata.fda.gov/drugsatfda_docs/label/2015/018780s150lbl.pdf) products in the `biologics` space.

* [Walsh G. Biopharmaceutical benchmarks 2018. Nat Biotechnol. 2018 Dec 6;36(12):1136-1145. doi: 10.1038/nbt.4305. PMID: 30520869](https://pubmed.ncbi.nlm.nih.gov/30520869/)
"""

# ╔═╡ 46a13c8a-1acf-4e22-9f29-44ac8519ce26
md"""
##### Origin story: The first biologic Humulin-R

Herbert Boyer et al. developed a portable expression system (circular DNA construct called a plasmid) to express the human insulin protein in _Escherichia coli_ (Humulin R). He formed a company called [Genentech](https://www.gene.com) and licensed the technology to [Eli Lilly & Company](https://www.lilly.com) in 1982. The rest is biotech history. 
"""

# ╔═╡ 0a830a04-cfdd-4086-9f4d-8af6deed7e8f
PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS,"Figs-HumulinR.png"))

# ╔═╡ 419db086-ce0c-4e47-9f0e-63dd1a5150d3
md"""

__Insulin references__:
* [Cohen SN, Chang AC, Boyer HW, Helling RB. Construction of biologically functional bacterial plasmids in vitro. Proc Natl Acad Sci U S A. 1973 Nov;70(11):3240-4. doi: 10.1073/pnas.70.11.3240. PMID: 4594039; PMCID: PMC427208.](https://pubmed.ncbi.nlm.nih.gov/4594039/)
* [Riggs AD. Making, Cloning, and the Expression of Human Insulin Genes in Bacteria: The Path to Humulin. Endocr Rev. 2021;42(3):374-380. doi:10.1210/endrev/bnaa029](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8152450/)
* [Hirsch IB, Juneja R, Beals JM, Antalis CJ, Wright EE. The Evolution of Insulin and How it Informs Therapy and Treatment Choices. Endocr Rev. 2020;41(5):733-755. doi:10.1210/endrev/bnaa015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7366348/)

__Insulin Litigation__ (since is always _just_ science):
* [Fox JL. Insulin patent dispute revisits old biotechnology battleground. Nat Biotechnol. 1997 Apr;15(4):307. doi: 10.1038/nbt0497-307. PMID: 9094114.](https://pubmed.ncbi.nlm.nih.gov/9094114/)

"""

# ╔═╡ cf05dc5a-b772-45e6-9f16-64448ab0f42b
md"""
### What makes biology so cool? Choices. Many possible choices.
I became interested in biology, metabolism, models, etc because, unlike traditional chemical systems, biological systems are controlled. A cell can sense the world around them, take stock of their internal state and make choices, in essence, they can reprogram themselves to meet a changing world. There are fast choices that operate on a milli- or μ-second time scale (regulation of enzyme activity, called [allosteric regulation](https://en.wikipedia.org/wiki/Allosteric_regulation)) and slow choices (gene expression) which operate on a time scale of tens of minutes. 

##### The origin story of omics. Yeast reprogramming
The promise of `omics` technologies was that we could measure everything, all at the same time, inside a population of cells. Promise and reality turned out to be a little different. The next wave of omics technology (occurring now) is that we can measure everything inside a single cell. But we have the same basic problem as before: what can we do with data that is noisy, and often has no direct physical interpretation e.g., has no units or is scaled, etc. From a traditional modeling perspective: not much. But still: it's cool.

* [DeRisi JL, Iyer VR, Brown PO. Exploring the metabolic and genetic control of gene expression on a genomic scale. Science. 1997 Oct 24;278(5338):680-6. doi: 10.1126/science.278.5338.680. PMID: 9381177.](https://pubmed.ncbi.nlm.nih.gov/9381177/)
"""

# ╔═╡ 092b0b15-8693-4280-a1b3-e23604db94f7
md"""
Yeast cells make a choice, by how and why? Fig. 3 Reproduced from DeRisi et al, Science 278, 1997.
"""

# ╔═╡ 9994bb96-e76d-416a-913b-8f53cbb4fbcc
PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS, "Fig-3-DeRisi-Brown-Science-278-1997.png"))

# ╔═╡ 731ff161-434b-4ae1-89a5-37d6d6f165dc
md"""
### Conclusions

In this lecture we:

* Introduced Metabolic Engineering: using models and engineering principles to design the production of metabolic products
* Disucssed two classes of biotechology: industrial and medical biotechnology. Industrial biotechnology is primarly focused on consumer products, while medical biotechnology focusses on products humart to human health.
* Discussed modes of biological regulation (one of the reasons this problem is hard). regulation can occur on a fast time scale (regulation of enzyme activity) and on a slow scale (regulation of gene expression).
"""

# ╔═╡ 76eb92dc-600f-49c2-8ca3-54422cff4f00
md"""
### Next time

* We'll build our first mathematical model of coupled enzyme-caytalyzed reactions and take stocxk of what we know, and what we don't know. 
* We'll introduce a life changing way of thinking: Constraints based analysis
"""

# ╔═╡ 42299ff5-93ac-4447-9b93-bf9222ef821b
TableOfContents(title="📚 Table of Contents", indent=true, depth=5, aside=true)

# ╔═╡ 3f5396b0-4b17-4a20-9f88-4e7dadbcbe82
html"""
<script>
	// initialize -
	var section = 0;
	var subsection = 0;
	var headers = document.querySelectorAll('h3, h5');
	
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
	    }
		// update the header text 
		header.innerText = numbering + " " + text;
	};
</script>"""

# ╔═╡ 64329198-85ce-47ea-a8d9-e664481a9658
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
PlutoUI = "~0.7.30"
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
git-tree-sha1 = "5c0eb9099596090bb3215260ceca687b888a1575"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.30"

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
# ╟─45e60528-7b8d-11ec-3aa3-65625ebbbfb2
# ╟─96d697e1-d190-4d9b-aeb8-53466cc1cce2
# ╟─87c121a0-91ed-436f-a981-a8a5198c3226
# ╟─2e67cdc9-9cd4-4a20-9ac7-1ecfd4b46967
# ╟─21f2c20b-0386-4532-a81d-18b05df6c3cf
# ╟─226b8bdd-13aa-4ab6-b4cb-73dd9cba5cb2
# ╟─46a13c8a-1acf-4e22-9f29-44ac8519ce26
# ╟─0a830a04-cfdd-4086-9f4d-8af6deed7e8f
# ╟─419db086-ce0c-4e47-9f0e-63dd1a5150d3
# ╟─cf05dc5a-b772-45e6-9f16-64448ab0f42b
# ╟─092b0b15-8693-4280-a1b3-e23604db94f7
# ╟─9994bb96-e76d-416a-913b-8f53cbb4fbcc
# ╟─731ff161-434b-4ae1-89a5-37d6d6f165dc
# ╠═76eb92dc-600f-49c2-8ca3-54422cff4f00
# ╠═42299ff5-93ac-4447-9b93-bf9222ef821b
# ╠═a3101ae7-2699-43e7-aa28-332b0447bbc5
# ╠═3f5396b0-4b17-4a20-9f88-4e7dadbcbe82
# ╠═64329198-85ce-47ea-a8d9-e664481a9658
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
