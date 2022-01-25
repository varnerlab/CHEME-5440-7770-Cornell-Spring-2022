### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

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
##### General continuum species material balance model

Let's use the general balance equations above, and develop a general model that we can use to write material balance equations for both extracellular and intracellular models. In particular, let's expand the accumulation term (chain rule):

$$\frac{d}{dt}\left(C_{i}\beta\right) = \beta\frac{dC_{i}}{dt}+C_{i}\frac{d\beta}{dt}\qquad{i=1,\dots,\mathcal{M}}$$

which can be substituted into the general balance equation to give the General continuum species material balance model:

$$\frac{dC_{i}}{dt} = \frac{1}{\beta}\left(\sum_{s=1}^{\mathcal{S}}d_{is}\dot{n}_{is}\right) + \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - \frac{C_{i}}{\beta}\frac{d\beta}{dt}\qquad{i=1,\dots,\mathcal{M}}$$

"""

# ╔═╡ 053889fa-1895-4a8a-b2b1-6df85b275c0f
md"""
##### Dynamic intracellular material balances equations
When building models of the _intracellular_ state of a cell, there are no physical transport terms (material is not flowing or from the cell). Thus, the transport terms for the intracellular balances are equal zero:

$$\frac{1}{\beta}\left(\sum_{s=1}^{\mathcal{S}}d_{is}\dot{n}_{is}\right) = 0$$

Intracellular balances in industrial metabolic engineering applications for liquid cultures are typically written in _biomass specific units_, e.g., 

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
###### Batch case 

In a batch culture, there is no in/out flow e.g., no volume change. In this case the dilution terms becomes:

$$\frac{dC_{i}}{dt} = \left(\sum_{j=1}^{\mathcal{R}}\sigma_{ij}\hat{r}_{j}\right) - {\mu}C_{i}\qquad{i=1,\dots,\mathcal{M}}$$

where $\mu$ denotes the specific growth rate of the cells in the batch culture $\mu = X^{-1}\dot{X}$.
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
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.1"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─6627255a-b788-4aaa-bf8e-c98a39553dea
# ╟─a5da8cc5-7ec8-4072-bf36-9e54e77c2a3f
# ╟─cf04d58a-230d-47ea-ba4a-072719cc9aa2
# ╟─053889fa-1895-4a8a-b2b1-6df85b275c0f
# ╟─0bc7cc05-d2da-4b41-9c5a-06c21e3caf70
# ╟─d31d5e9c-a5d7-4f27-899c-bacfe350cd64
# ╟─13d7ca39-f956-479e-b560-d6338c1d8b12
# ╠═b1b251d8-7e23-11ec-09d1-97ace15b3bd3
# ╟─9833473a-3fc3-4b61-bb07-050d6e159cb2
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
