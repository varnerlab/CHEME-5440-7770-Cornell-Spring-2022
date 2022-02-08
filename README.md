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
