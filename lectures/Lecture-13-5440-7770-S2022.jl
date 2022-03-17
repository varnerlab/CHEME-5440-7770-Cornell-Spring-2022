### A Pluto.jl notebook ###
# v0.18.2

using Markdown
using InteractiveUtils

# ╔═╡ 2c895e14-e71a-435d-86d4-c1b101877b3c
md"""
### The End of the Beginning: Whole-Cell Metabolic Engineering

In the first lecture, we introduced metabolic engineering. Today, we'll close the loop at put together all we learned. Metabolic engineering is the practice of optimizing genetic and regulatory processes within cells to increase the cell's production of a desired small molecule or protein product of interest. 

Metabolic engineers manipulate the biochemical networks cells use to convert raw materials into molecules necessary for the cell's survival. Metabolic engineering specifically seeks to:

1. Mathematically model biochemical networks, calculate the yield (product divided by substrate) of valuable products and identify parts of the network that constrain the production of these products of interest. 
1. Use genetic engineering techniques to modify the biochemical network to relieve constraints limiting production. Metabolic engineers can then model the modified network to calculate the new product yield and identify new constraints (back to 1).

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
##### Aside: A Question of Timescales
One open question is whether we need to solve the gene expression and flux estimation problem at the same time. The time scale of enzyme activity (and hence metabolic flux) is much faster than gene expression. Why? 

* The $k_{cat}$ for [phosphofructokinase-2 (pfk) in _E.coli_ is 9240 s$^{-1}$](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104955&ver=7&trm=kcat+pfk+in+E.+coli&org=). Thus, pfk has a characteristic time scale of $\tau\sim{k_{cat}}^{-1}$ or about 1$\times$10$^{-4}$ s.
* The $k_{cat}$ for RNAP polymerase is $k_{cat}$=$e_{x}L^{-1}$. If we assume [$e_{x}\sim{35}$ nt/s](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111871&ver=1&trm=elongation+rate+RNAP+E.coli+&org=) and average read length [$L\sim{924}$ nt](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111922&ver=3&trm=avarge+gene+length+in+E.+coli&org=), this gives $k_{cat}\sim$ 0.038 s$^{-1}$ or a timescale of $\tau\sim$26.4 s.

Thus, the time scale of metabolic flux (the timescale of enzyme activity) is 5-orders of magnitude faster than gene expression. Therefore, from the perspective of metabolic flux, gene expression seems constant; hence, we treat these scales separately in traditional models. 
"""

# ╔═╡ 06c3eafa-10b3-4a47-9a06-39f6abd29db3
md"""
### Using metabolic models to improve metabolic performance
Ultimately a metabolic engineer uses a model to improve the performance of a production host or cell-free system. There is a huge variety of design approaches to improve metabolic function; let's look at a recent one from Cornell:

* [Wayman JA, Glasscock C, Mansell TJ, DeLisa MP, Varner JD. Improving designer glycan production in Escherichia coli through model-guided metabolic engineering. Metab Eng Commun. 2019;9:e00088. Published 2019 Mar 29. doi:10.1016/j.mec.2019.e00088](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC6454127/)

"""

# ╔═╡ b15b540b-8f13-4f21-b2da-510d34143b33
md"""
### Summary and conclusions
"""

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
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─2c895e14-e71a-435d-86d4-c1b101877b3c
# ╟─ce5c8e9c-7652-4c3e-98a7-98e3567dc695
# ╟─bfd134fe-b311-4768-a53b-2c5c4c1b121a
# ╟─06c3eafa-10b3-4a47-9a06-39f6abd29db3
# ╠═b15b540b-8f13-4f21-b2da-510d34143b33
# ╟─7988263a-a5f5-11ec-3836-150105d04264
# ╟─dee24b55-a5c5-4a47-9e1c-88ea0cf1cd61
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
