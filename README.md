# Evolutionary dynamics of agents playing a nonlinear public goods game under genetically homophilic group formation

## About the code

This code is used to study how genetically homophilic group formation influences the evolutionary dynamics of agents playing a nonlinear public goods game, primarily the threshold game. It uses the higher-order genetic association approach of Ohtsuki (2014; Phil Trans Roy Soc B), and features three different group-formation models, where groups are formed sequentially and current members tend to recruit or attract kin.

## About the manuscript

Kristensen, N.P., Ohtsuki, H., Chilsholm, R.A. Ancestral social environments and nonlinear payoffs can explain cooperation in human societies

**Abstract**

Human cooperation (paying a cost to benefit others) is puzzling from a Darwinian perspective,
particularly in groups with strangers who cannot repay nor are family members. Real-world
scenarios typically involve a nonlinear relationship between the payoff and the number of
cooperators in a group, and nonlinear payoffs can allow cooperation between strangers to persist
evolutionarily provided that a large enough proportion of the population are already cooperators.
However, if a lone cooperator faces a conflict between the group’s and their own interests (a social
dilemma), that raises the question of how cooperation arose in the first place. We use a
mathematically tractable evolutionary model to formalise a chronological narrative that has
previously only been investigated verbally: given that ancient humans interacted mostly with family
members (genetic homophily), cooperation evolved first by kin selection, and then persisted in
situations with nonlinear payoffs as homophily declined or even if interactions with strangers
became the norm. The model also predicts the coexistence of cooperators and defectors observed in
the human population (polymorphism), and may explain why cooperators in behavioural
experiments prefer to condition their contribution on the contributions of others (conditional
cooperation in public goods games).

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
can be calculated using the script: \\
`scripts/members_recruit/sum_product_mistakes/save_sum_prod_mistakes.py`.

Precalculated terms up to size $n=18$ and are stored in:
`/results/members_recruit/sum_product_mistakes/`.
