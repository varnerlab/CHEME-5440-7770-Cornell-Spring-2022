# CHEME-5440-7770-Cornell-Spring-2022

Repository for CHEME 5440 and 7770 lectures, examples, problem sets, codes, etc.

This repository holds lecture notes and example problems discussed in lecture. The problems are structured in [Pluto](https://github.com/fonsp/Pluto.jl) notebooks and use the [Julia](https://julialang.org) programming language.

### Installing Julia and Pluto

[Julia](https://julialang.org) is open source, free and runs on all major operating systems and platforms. To install
[Julia](https://julialang.org) and [Pluto](https://github.com/fonsp/Pluto.jl) please check out the tutorial for
[MIT 18.S191/6.S083/22.S092 course from Fall 2020](https://computationalthinking.mit.edu/Fall20/installation/).

1. [Install Julia (we are using v1.7.x, newer versions of Julia should also work)](https://julialang.org/downloads/)
1. [Install Pluto.jl](https://github.com/fonsp/Pluto.jl#installation)
1. Clone this repo:

    ```bash
    git clone https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022.git
    ```

1. From the Julia REPL (`julia`), run Pluto (a web browser window will pop-up):

    ```julia
    julia> using Pluto
    julia> Pluto.run()
    ```

    _Or you can simply type the following in a terminal:_

    ```bash
    julia -E "using Pluto; Pluto.run()"
    ```

1. From Pluto, open one of the `.jl` lecture notebook files located in the `CHEME-5440-7770-Cornell-Spring-2022/lectures` directory—enjoy!

### IDE

We use [Visual Studio Code](https://code.visualstudio.com) and GitHub desktop because they are awesome!

### Lectures

* [Lecture 1: Introduction to Metabolic Engineering](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-1-5440-7770-S2022.jl.html)
* [Lecture 2: Constraint-Based Perspective and Tools](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-2-5440-7770-S2022.jl.html)
* [Lecture 3: Our first metabolic flux calculation](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-3-5440-7770-S2022.jl.html)
* [Lecture 4: Computing the reversibility of a metabolic reaction](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-4-5440-7770-S2022.jl.html)
* [Lecture 5: Minimum Energy Metabolic Flux Distributions](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-5-5440-7770-S2022.jl.html)
* [Lecture 6: Introduction to Enzyme Kinetics](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-6-5440-7770-S2022.jl.html)
* [Lecture 7: Multisubstrate enzyme kinetics and inhibitors](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-7-5440-7770-S2022.jl.html)
* [Lecture 8: Carbon catabolite repression and central dogma](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-8-5440-7770-S2022.jl.html)
* [Lecture 9: Boolean and Fuzzy Logic in Metabolic Models](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-9-5440-7770-S2022.jl.html)
* [Lecture 10: Discrete State Models for Transcriptional Regulation](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-10-5440-7770-S2022.jl.html)
* [Lecture 11: Stochastic simulation of single cell gene expression](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-11-5440-7770-S2022.jl.html)
* [Lecture 12: TX/TL and flux balance analysis](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-12-5440-7770-S2022.jl.html)
* [Lecture 13: The end of the beginning](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/Lecture-13-5440-7770-S2022.jl.html)

### Problem set solutions

* [Problem set 2 solution](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/PS2-5440-7770-Soln.jl.html)
* [Problem set 3 solution](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/html/PS3-5440-7770-S2022-C2-Soln.jl.html)

### Prelim 1 solutions
* [Prelim 1 Q2 solution](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/prelim_1/html/P1-Q2-Soln-S2022.jl.html)
* [Prelim 1 Q3 solution](https://htmlview.glitch.me/?https://github.com/varnerlab/CHEME-5440-7770-Cornell-Spring-2022/blob/main/prelim_1/html/P1-Q3-Soln-S2022.jl.html)
