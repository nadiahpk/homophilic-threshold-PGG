# Evolutionary dynamics of agents playing a nonlinear public goods game under genetically homophilic group formation

[![DOI](https://zenodo.org/badge/443937530.svg)](https://zenodo.org/badge/latestdoi/443937530)

## About the code

This code is used to study how genetically homophilic group formation influences the evolutionary dynamics of agents playing a nonlinear public goods game, primarily the threshold game. It uses the higher-order genetic association approach of Ohtsuki (2014; Phil Trans Roy Soc B), and features three different group-formation models, where groups are formed sequentially and current members tend to recruit or attract kin.

## About the manuscript

Kristensen, N.P., Ohtsuki, H., Chilsholm, R.A. (2022) Ancestral social environments plus nonlinear benefits can explain cooperation in human societies, *Scientific Reports*

**Abstract**

Human cooperation (paying a cost to benefit others) is puzzling from a Darwinian perspective,
particularly in groups with strangers who cannot repay nor are family members. The beneficial effects
of cooperation typically increase nonlinearly with the number of cooperators, e.g., increasing returns
when cooperation is low and diminishing returns when cooperation is high. Such nonlinearity can
allow cooperation between strangers to persist evolutionarily if a large enough proportion of the
population are already cooperators. However, if a lone cooperator faces a conflict between the group's
and its own interests (a social dilemma), that raises the question of how cooperation arose in the first
place. We use a mathematically tractable evolutionary model to formalise a chronological narrative
that has previously only been investigated verbally: given that ancient humans interacted mostly with
family members (genetic homophily), cooperation evolved first by kin selection, and then persisted
in situations with nonlinear benefits as homophily declined or even if interactions with strangers
became the norm. The model also predicts the coexistence of cooperators and defectors observed in
the human population (polymorphism), and may explain why cooperators in behavioural experiments
prefer to condition their contribution on the contributions of others (conditional cooperation in public
goods games).

## Code quickstart

A quickstart tutorial to calculate the change in the proportion of Cooperators, $\Delta p$, can be found in: `/tutorials/calculate_deltap.pdf`

Matrix $M$ used to calculate $\boldsymbol{\theta}$ from $\mathbf{F}$ is calculated using the the script:
`/scripts/matrix_M/save_matrix_Ms.py`.

Precalculated $M$ matrices up to size $n=24$ are stored in `/results/matrix_M/`
and can be read using `read_matrix_M()` in `/functions/my_functions.py`.

A quickstart tutorial for how to use `read_matrix_M()` is provided in:
`/tutorials/matrix_M.pdf`.

The combinatorial term needed to calculate family partition structure probabilities for the members-recruit
group-formation model (i.e., $\sum_{\vec{\mathbf{n}}} \sum_{\mathbf{m}} C(\mathbf{m}, \vec{\mathbf{n}}) \prod_k m_k$)
can be calculated using the script: 
`scripts/members_recruit/sum_product_mistakes/save_sum_prod_mistakes.py`.

Precalculated terms up to size $n=18$ and are stored in:
`/results/members_recruit/sum_product_mistakes/`.

## License

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or distribute this software, either in source code form or as a compiled binary, for any purpose, commercial or non-commercial, and by any means.

In jurisdictions that recognize copyright laws, the author or authors of this software dedicate any and all copyright interest in the software to the public domain. We make this dedication for the benefit of the public at large and to the detriment of our heirs and successors. We intend this dedication to be an overt act of relinquishment in perpetuity of all present and future rights to this software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <http://unlicense.org/>
