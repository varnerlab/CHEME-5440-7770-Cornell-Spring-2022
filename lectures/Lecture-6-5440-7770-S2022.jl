### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 8ba6a00a-7a53-4dcc-af4d-cb85aa9b4d68
md"""
### Simple and Complex Models of Enzyme Kinetics

The flux bounds are essential constraints in flux balance analysis calculations and the convex decomposition of the stoichiometric array. Beyond their role in the flux estimation problem, the flux bounds are _integrative_, i.e., these constraints integrate many types of genetic and biochemical information into the flux estimation problem. 

Flux bounds constrain the values that each reaction in a metabolic network can take. A general model for these bounds is given by:

$$-\delta_{j}\left[{V_{max,j}^{\circ}}\left(\frac{e}{e^{\circ}}\right)\theta_{j}\left(\dots\right){f_{j}\left(\dots\right)}\right]\leq{v_{j}}\leq{V_{max,j}^{\circ}}\left(\frac{e}{e^{\circ}}\right)\theta_{j}\left(\dots\right){f_{j}\left(\dots\right)}$$

where $V_{max,j}^{\circ}$ denotes the maximum reaction velocity (units: flux) computed at some characteristic enzyme abundance (units: concentration), the ratio $e/e^{\circ}$ is a correction for enzyme abundance (units: dimensioness), $\theta_{j}\left(\dots\right)\in\left[0,1\right]$ is the fraction of maximial enzyme activity (a function or measurement producing units: dimensionless), and $f_{j}\left(\dots\right)$ is a function describing the substrate dependence of the reaction rate $j$ (units: dimensionless). Both $\theta_{j}\left(\dots\right)$ and $f_{j}\left(\dots\right)$ could have associated parameters, e.g., saturation or binding constants, etc.  Finally, the quanity $\delta_{j}\in\left\{0,1\right\}$ is a _binary_ variable: 
* If reaction $j$ is __reversible__ $\delta_{j}=1$ or,
* If reaction $j$ is __irreversible__ $\delta_{j}=0$

Today, let's focus on approaches for computing the value of these bounds, i.e., the form and value of the terms on the brackets. 

In this lecture, we will:

1. Develop simple models of enzyme kinetics using the Michaelis–Menten approach
1. Introduce the fundamental concepts underlying more complex models of allosteric regulation (the fast control mechanisms that we mentioned previously)
1. Introduce effective discrete regulation models to capture allosteric regulation
"""

# ╔═╡ 17724757-604e-4c77-a42b-238ba121e88f
md"""
### Simple: Michaelis–Menten kinetics
Let's assume we have a well-mixed test tube containing an enzyme $E$ (a protein that catalyzes chemical reactions), which converts substrate $S$ (the starting compound) 
into product $P$ according to three elementary reactions:

$$\begin{eqnarray}
	E+S&\rightleftharpoons&{E:S}\\
	{E:S}&\longrightarrow&E+P
\end{eqnarray}$$

The kinetics of each elementary step can be written using mass-action kinetics, i.e.,

$$\begin{eqnarray}
	r_{1} & = & k_{1}\left[E\right]\left[S\right]\\
	r_{2} & = & k_{2}\left[E:S\right]\\
	r_{3} & = & k_{3}\left[E:S\right]
\end{eqnarray}$$

where $\left[\cdot\right]$ denotes a species concentration, and $k_{j}$ denotes the rate constant governing the $jth$ elementary reaction:

* The rate $r_{1}$ describes the _association_ rate between the enzyme and substrate,
* The rate $r_{2}$ represents the rate of _dissociation_ of the enzyme-substrate complex, and
* The $r_{3}$ denotes the rate of _chemical conversion_ of the bound substrate into the product (we assume the dissociation of the product from the enzyme is fast).

The enzyme must obey the relationship:

$$\left[E_{T}\right] = \left[E\right] + \left[E:S\right]$$

where $\left[E_{T}\right]$ denotes the total enzyme concentration in the tube, 
$\left[E\right]$ denotes the free enzyme concentration (not bound to substrate) while
$\left[E:S\right]$ denotes the enzyme substrate complex. 

To estimate the _overall_ rate of enzymatic conversion ($v$) of $S$ to $P$, we need to stipulate a single rate-limiting step out of the set of elementary reactions describing the transformation.
Let's assume that the rate of chemical conversion ($r_{3}$) is the \emph{slowest} step, i.e., the substrate bounces on/off the enzyme quickly with only a fraction of these binding events resulting
in a successful chemical transformation. Thus, the overall rate is then given by:

$$v = k_{3}\left[E:S\right]$$

Let's also assume that we already know (or can estimate) the rate constants $k_{1},k_{2}$ and $k_{3}$. When this is true, the only unknown is $\left[E:S\right]$.
However, we can relate $\left[E:S\right]$ to variables we know ($E_{T}$ and at least initially $S$) through the enzyme balance, and the _pseudo-steady-state assumption_ for the reaction intermediate 
$\left[E:S\right]$:

$$\frac{d\left[E:S\right]}{dt} = k_{1}\left[E\right]\left[S\right] - k_{2}\left[E:S\right] - k_{3}\left[E:S\right]\simeq{0}$$

Rearranging Eqn. \eqref{eqn:pssa} and solving for $\left[E:S\right]$ gives the relationship:

$$\left[E:S\right]\simeq\frac{k_{1}}{k_{2}+k_{3}}\left[E\right]\left[S\right]$$

where the ratio of constants is defined as the Michaelis-Menten saturation coefficient or $K_{M}$:

$$\frac{1}{K_{M}}\equiv\frac{k_{1}}{k_{2}+k_{3}}$$

Substituting Eqn. \eqref{eqn:es} into the overall rate yields:

$$v = k_{3}\frac{\left[E\right]\left[S\right]}{K_{M}}$$

However, we do not know the free enzyme concentration of $\left[E\right]$. 
We use the total enzyme balance to get $\left[E\right]$. 
Substituting Eqn. \eqref{eqn:es} into the enzyme balance Eqn. \eqref{eqn:total-enzyme-balance} and solving for $\left[E\right]$ yields:

$$\left[E\right] = \frac{\left[E_T\right]K_{M}}{K_{M}+\left[S\right]}$$

Lastly, we can substitute Eqn. \eqref{eqn:eqn-finally} into Eqn. \eqref{eqn:final-v} to arrive at the final expression for $v$:

$$v = V_{max}\frac{\left[S\right]}{K_{M}+\left[S\right]}$$

where $V_{max}\equiv{k_{3}}\left[E_{T}\right]$. 

__Limiting cases:__
* when $S\gg{K}_{M}$, the rate becomes close to $V_{max}$.
* when $S\ll{K}_{M}$ the rate appears to be linear with respect to substrate concentration. Lastly, it is easy to show that when $K_{M}\simeq S$ the reaction rate equals $v\simeq 1/2V_{max}$. 


"""

# ╔═╡ 7efaf27c-8e0a-40f9-ac28-af90357450a3
begin
	# run kinetic sims -
end

# ╔═╡ 55ea8324-f36d-40bd-99e3-92e7fcbafca9
md"""
### Complex: MWC and Sequential kinetic models

"""

# ╔═╡ 3502da52-7d2b-4c79-82ce-7424d756cd9b
md"""
### Medium: Effective discrete state kinetic models

Let the model take the form $\hat{r}_{j} = r_{j}v\left(...\right)_{j}$
where $\hat{r}_{j}$, the overall rate of the PFK reaction ($\mu$M h$^{-1}$), 
is the product of a kinetic limit $r_{j}$ ($\mu$M h$^{-1}$), i.e., the maximum rate of conversion, 
and a control variable $0\leq v\left(...\right)_{j}\leq 1$ (dimensionless) that describes the influence
of effector molecules. Since we have only a single enzyme let j = 1. 

The model for the kinetic limit was given in the problem:

$r_{1} = k_{cat}E_{1}\left(\frac{F6P}{K_{F6P}+F6P}\right)\left(\frac{ATP}{K_{ATP}+ATP}\right)$

along with all the parameters associated with this model and the concentration of enzyme ($E_{1}$), F6P, and ATP. Thus, we can calculate the kinetic limit. However, we need to postulate a form for the v-variable and estimate the parameters in the v-variable expression. We'll do this from the data.

###### v-model
In the problem we proposed a three state model for $v\left(...\right)_{j}$: State 0: no effector+no substrate, State 1: no effector+substrate and State 2: effector+substrate. Following the procedure discussed in class, given this three state model, the v-variable takes the form:

$v_{1} = \frac{W_{1}+W_{2}f_{2}}{1+W_{1}+W_{2}f_{2}}$

where $W_{i}$ denotes the weight of state i (dimensionless), and $f_{2}$ denote the fraction of bound effector given by: $f_{2} = (x/K_{2})^{n_{2}}/(1+(x/K_{2})^{n_{2}})$. The quanity x denotes the concentration of the effector 3$^{\prime}$-5$^{\prime}$-AMP, $K_{2}$ denotes a binding constant (units of concentration) and $n_{2}$ denotes a binding cooperativity parameter (dimensionless). State 0 has no activity, while State 1 (with weight $W_{1}$) has basal enzyme activity, and State 2 describes the enhanced activity following from effector binding (with weight $W_{2}f_{2}$). 

 - The probability of state 0 is $p_{0}$ = $1/Z$ where $Z$=$1+W_{1}+W_{2}f_{2}$.
 - The probability of state 1 is $p_{1}$ = $W_{1}/Z$ where $Z$=$1+W_{1}+W_{2}f_{2}$.
 - The probability of state 2 is $p_{2}$ = $W_{2}f_{2}/Z$ where $Z$=$1+W_{1}+W_{2}f_{2}$.

"""

# ╔═╡ d3568c5d-16bf-4698-9891-0be65d62b36c
md"""
### Summary and Conclusions

In this lecture we:

1. Developed simple models of enzyme kinetics using the Michaelis–Menten approach
1. Introduced the fundamental concepts underlying more complex models of allosteric regulation (the fast control mechanisms that we mentioned previously)
1. Introduced effective discrete regulation models to simulate allosteric regulation

"""

# ╔═╡ 55d4ab48-9589-4fd9-ac06-3338fdb418c4
md"""
### Next Time

"""

# ╔═╡ 06a90381-b4cb-4d32-af0d-2f084364129c
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

# ╔═╡ b6890de0-8923-11ec-3552-3113bdd53f86
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
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[deps]
"""

# ╔═╡ Cell order:
# ╟─8ba6a00a-7a53-4dcc-af4d-cb85aa9b4d68
# ╟─17724757-604e-4c77-a42b-238ba121e88f
# ╠═7efaf27c-8e0a-40f9-ac28-af90357450a3
# ╟─55ea8324-f36d-40bd-99e3-92e7fcbafca9
# ╟─3502da52-7d2b-4c79-82ce-7424d756cd9b
# ╟─d3568c5d-16bf-4698-9891-0be65d62b36c
# ╠═55d4ab48-9589-4fd9-ac06-3338fdb418c4
# ╠═06a90381-b4cb-4d32-af0d-2f084364129c
# ╠═b6890de0-8923-11ec-3552-3113bdd53f86
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
