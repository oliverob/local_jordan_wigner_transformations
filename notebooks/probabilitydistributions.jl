### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 82eb253a-2378-4074-9dab-4cb3c7eea868
using PlutoUI

# ╔═╡ e3f5a522-6051-489e-9f72-89d9ecfc9b2c
using StatsPlots

# ╔═╡ bf81ee3b-dd04-4983-9f2e-b52badad1da1
 using HypertextLiteral

# ╔═╡ 8261efb7-0e4f-4866-92cb-50b0d1052d5c
using LaTeXStrings

# ╔═╡ d8ed59fe-a83e-11eb-07f8-d7123a34a1c5
md"# Probability Distributions"

# ╔═╡ 4fd5c54c-4c43-454a-81dd-5ac7245feb85
TableOfContents(; depth = 2)

# ╔═╡ cf85ab67-9d99-4551-93b6-db3494a3bc73
md"## Introduction"

# ╔═╡ 46f5142e-369f-4c78-8d91-cd267933d8be
md"Probability is a conserved quantity distributed over a space ``S``.  A probability distribution is then a mathematical description of how much of this conserved quantity, named probability, is assigned to most[^1] subsets ``A`` of ``S``, ``A \subseteq S``.  The set ``S`` is called the **sample space**."

# ╔═╡ 57d26450-bb69-46b2-93a8-bbee93c91436
md"**Example: Gas Station Keypads**.  Imagine you are at a gas station, about to enter the zip code associated with your credit card.  Have you ever noticed that the digits of the local zip codes are often more worn than the other digits?  In fact, the wear on the keypad describes a probability distribution over the sample space ``S = \{0, 1, 2, 3, 4, 5, 6, 7, 8, 9\}``.  Probability is higher for, or distributed more to, the numbers of the digits of the local zip codes, and is lower for, or distributed less to, the other numbers.  Further, we should expect repeated numbers within the local zip codes to have higher probability than less common numbers of the local zip codes."

# ╔═╡ 18f3e3ce-9fd2-4a55-bcb2-5224f4bb943c
md"**Example: Fair Die**.  A fair die is a more common choice for a first example.  The sample space ``S = \{1, 2, 3, 4, 5, 6\}`` and probability is evenly, or uniformly, distributed over the singletons (single element subsets) of ``S``."

# ╔═╡ 0dfb6258-4451-4d64-a704-9da3ec92ca07
md"**Example: Fair Coin**.  A fair coin with sides ``H`` for heads and ``T`` for tails is another common first example choice. A fair coin is another example of a uniform distribution of probability over its sample space ``S = \{H, T\}``.  A nice feature of this example is that the sample space is not numbers, thus highlighting the generality of probability distributions.  Probability distributions need not apply to sample spaces of only numbers.  The sample space ``S`` can really be almost any set."

# ╔═╡ 8aa81756-2bc5-462b-bbc2-84708fab3a70
md"The gas station keypad example above contrasts nicely against the next two examples.  The gas station keypad's distribution is a non-uniform distribution of probability (over singletons of the sample space), whereas the other two examples represent a uniform distribution (of probability over singletons of the sample space).  These three examples will thus help us generalize our concept of probability."

# ╔═╡ 227003ce-5470-4f04-93a5-6cbe84a2358d
md"A probability distribution is any function that assigns a number, also referred to as probability, to (most[^1]) subsets of a sample space and also satisfies the axioms of probability, provided below.  It's common to think of probability as a number, but in fact it is a function that assigns probability to subsets of ``S``.  After discussing the Axioms of Probability, we'll discuss some more complicated examples of probability distributions, to help ensure that we come away with a general concept of probability that can be more easily adapted to new situations."

# ╔═╡ febc72ae-2ef6-484a-bc90-5844eeedae76
md"## Probability, a Set Function"

# ╔═╡ 2b0c948d-2eb1-461a-8999-0802b7e9e819
md"A probability distribution is a function that applies sets, specifically a function that applies to (most) subsets of the sample space ``S``.  Hence it is a set function.  For an arbitrary subset ``A`` of ``S``, ``A \subseteq S``, we write

$$\mathbb{P}[A] = p$$

for some ``p \in [0, 1]``. The specific value of ``p`` is determined by the specific distribution being used.  

* If ``\mathbb{P}`` represents a fair coin and ``A = \{H\}``, then we tend to believe ``\mathbb{P}[A] = \mathbb{P}[\{H\}] = 0.5``.
* If ``\mathbb{P}`` represents a fair die and ``A = \{1, 2\}``, then we would have ``\mathbb{P}[\{1, 2\}] = 2/6``.
* If ``\mathbb{P}`` represents a gas station keypad in Chico, CA, where the five local zip codes are ``(95926, 95928, 95929, 95973, 95976)``, then we might expect to see wear on local gas station keypads that, relatively closely, aligns with the probabilities ``\mathbb{P}[\{9\}] = 11/25`` and ``\mathbb{P}[\{5\}] = 5/25``, and ``\mathbb{P}[\{2\}] = 3 / 25``."

# ╔═╡ 6ddca240-5797-4564-adf3-0c788c8d5a80
md"## A Note on Notation"

# ╔═╡ d7da0c9d-8f19-49d0-b093-1d0816b6c8ed
md"Statisticians have a bad habit of reusing the same name, namely the letter P which I write in bold as ``\mathbb{P}``, for different functions. In fact, we'll see this same issue again in the notebook on density functions."

# ╔═╡ 0774a8f8-e6d1-44dc-b186-d481b713ba00
md"There are an uncountable number of probability distributions.  It is common in the world of statistics to reference each one with ``\mathbb{P}``.  One hopes that the specific probability distribution of interest is clear from the context, but one also imagines that navigating such distinctions is more difficult when first learning this material."

# ╔═╡ 33743bd0-0d79-4999-b93b-4d27d448b408
md"## Axioms of Probability"

# ╔═╡ 410b2fff-0ed5-4ff4-82c1-bc3409ffd0af
md"The three axioms of probability are

1. ``\mathbb{P}[A] \geq 0`` for (most) ``A \subseteq S``,
2. ``\mathbb{P}[S] = 1``, and
3. ``\mathbb{P}[\cup_n A_n] = \sum_{n} \mathbb{P}[A_n]`` where ``A_i \cap A_j = \emptyset`` and ``i \not= j``.
"

# ╔═╡ 1ad18652-cebf-4dfe-aa01-19834ba4b14a
md"The first axiom says that probability is bounded below by zero.  This should make intuitive sense, since we reasonably believe that probability can not be negative.  The second axiom suggests that probability is bounded above by 1, but technically says that the probability assigned to the entire sample space is 1.  Since any given probability distribution only assigns probability to subsets of the sample space ``S``, of which ``S`` is a valid subset, we can not distribute more than 100% of our conserved quantity, known as probability.  The third axiom says that probability is assigned to a union of disjoint sets as the sum of the components of the union."

# ╔═╡ 427adede-fe09-4c59-b8ca-df898c4a7924
md"The third axiom of probability is relatively easy to visualize.  The entire box, let's call it ``A``, below is written as a union of disjoint sets ``A_n`` for ``n \in \{1, 2, 3, 4\}``,

$$A = \cup_{n \in \{1, 2, 3, 4, 5\}} A_n.$$"

# ╔═╡ 8358fab7-4b42-45ec-bcdf-661747dde8f7
plot(0:5, t = :vline, xlimit = (0, 5), 
	annotations = [
		(0.5, 1.5, L"A_1"), 
		(1.5, 1.5, L"A_2"), 
		(2.5, 1.5, L"A_3"),
		(3.5, 1.5, L"A_4"),
		(4.5, 1.5, L"A_5")],
	grid = false, ticks = false, label = false,
	framestyle = :box)

# ╔═╡ 9612e871-4b10-4cff-9215-4da903d33d06
md"Since the ``A_n`` are disjoint, we can apply the third axiom of probability to find the probability of ``A``,

$$\mathbb{P}[A] = \mathbb{P}[\cup_{n \in \{1, 2, 3, 4, 5\}} A_n] = \mathbb{P}[A_1] + \mathbb{P}[A_2] + \mathbb{P}[A_3] + \mathbb{P}[A_4] + \mathbb{P}[A_5] = \sum_{n \in \{1, 2, 3, 4, 5\}} \mathbb{P}[A_n].$$"

# ╔═╡ 59588307-aaa6-45ab-862b-759b9dfc6810
md"## Generalizing Probability"

# ╔═╡ a91d6013-be1a-406d-a193-906e0c8d5a8c
md"Probability distributions dictate the assignment of probability to arbitrary subsets of the sample space.  Only in the simplest of cases do probability calculations become calculations of the number of elements in the subset ``A`` relative to the number of elements in an enclosing space, which is often but not always the sample space ``S``, such as

$$\mathbb{P}[A] = |A| / |S|.$$

For instance, this definition immediately won't apply to uncountable sample spaces."

# ╔═╡ 272ce517-d9ec-4afe-b5c1-90d6c4c89899
md"In general, probability is defined to be area under a probability distribution's  density function.  See the notebook on [TODO add link to] Density Functions for more details.  For now, we will assume that a density function of a probability distribution is the function that enables the assignment of probability to subsets of the sample space."

# ╔═╡ f10cecb7-cb04-43ee-8db0-f623d99115c1
md"Luckily, defining probability as area under a density function conforms to 

$$\mathbb{P}[A] = |A| / |S|$$

whent ``\mathbb{P}`` represents a uniform distribution (of probability) over a finite sample space ``S``, and enables more complex definitions of probability.  Let's start with the fair die example from above."

# ╔═╡ 26cdcadc-f212-4cd2-a159-c85fc1a6574f
md"The plot below represents the density function for a fair die, which we will see in Section Examples is written in statistics notation as ``\text{Uniform}(1, 6)``."

# ╔═╡ 6f0da015-f70a-490c-b116-9b9204371f76
md"We use this canonical image of area under a function as our starting point.  We will adapt this idea in a few specific ways so as to 

1. satisfy the Axioms of Probability, and 
2. expand the problem domain of probability distributions and thus statistics."

# ╔═╡ 28dce3a3-63b4-4496-bf99-2779d29881d9


# ╔═╡ 203f0879-875e-4a25-a7df-3c54f72ff41f
begin
	cols = get_color_palette(:auto, 17)
	gaskeypad = (0 => 0.0, 1 => 0.0, 2 => 3/25, 3 => 1/25, 4 => 0.0, 5 => 5/25,
		6 => 2/25, 7 => 2/25, 8 => 1/25, 9 => 11/25)	
	y = 0:10
	yfill = 6:9
	p2 = plot(y, [last(gkp) for gkp in gaskeypad], t = :scatter, 
		xlabel = L"x", ylabel = L"g(x)", label = false)
	for yf in yfill
	 plot!(p2, [yf, yf], [0.0, gaskeypad[yf+1][2]], t = :line, 
			c = cols[1], linewidth = 2, label = false)
	end
	p2
end

# ╔═╡ 483c1814-c485-4ac2-8940-8c0b3b68de59
begin
	x = 1:6
	xfill = 2:4
	p = plot(x, repeat([1/6], 6), t = :scatter,
		ylimits = (0, 0.25),
		xlabel = L"x", ylabel = "Density", label = false)
	for xf in xfill
		plot!(p, [xf, xf], [0, 1/6],
			t = :line, c = cols[1], label = false)
	end
	p
end


# ╔═╡ f54cb70e-b6f1-44f8-b36c-01c6f373fa90


# ╔═╡ b9c057ef-bff2-479d-a38f-473b7ca4c309


# ╔═╡ f52bc4fd-1a07-4963-8efc-f67152f76152
md"Alternatively, consider the density function ``g(x)`` from above.  There's not simple way to write this function, even a piece-wise function would be messy.  Often in cases like this, a table is used."

# ╔═╡ a1b61767-02c5-4678-8cde-3b86cc30bbd7


# ╔═╡ d2bca69d-353f-48d9-bde0-677892b3fa24
md"## Examples"

# ╔═╡ b284dddb-2da3-4e93-9ecb-b6592dcab087
md"Many common probability distributions are grouped into named families.  Although, not all distributions are so classified.  For instance, the gas station keypad example above is an example of an unnamed probability distribution."

# ╔═╡ 38f19a20-e002-41b5-976c-22c6af542e89
md"A point of clarity is in order.  Statistics authors often fail to mention the family component of probability distributions at all.  It is hard to blame such authors for neglecting this naunce.  As you'll see here, typing/saying, for instance, the Binomial family of probability distributions quickly becomes exhausting.  As such, I too will likely fail to mention the family aspect in all but these opening paragraphs.  Further, I'll probably even refer to probability distributions as, simply, distributions."

# ╔═╡ 7f773008-ac15-44bf-a276-15ea69c7e4b1
md"**Example: Binomial Distribution**.  "

# ╔═╡ e76c82e9-711c-4450-89fe-d795344d2cdb
md"### Gamma Distribution"

# ╔═╡ 07a816ec-d1bb-4720-b381-a493643e747d
md"## Footnotes"

# ╔═╡ 72b053b8-9327-4128-931e-95b3e9305571
md"
[^1]: Technically probability can only be assigned to measurable sets.  Unfortunately, I will probably not mention this technicality at all times when I should.  TODO see reference instead.
"

# ╔═╡ 0d89f93b-3a3a-4c1b-8e63-552ed41c6d4c
render_table(gkp) = @htl("""
	<table><caption><h4>Gas Keypad Distribution</h4></caption>
	<tbody>
	<tr><td>Digit$(map(gaskeypad) do gkp @htl("<td>$(string(first(gkp)))") end)
	<tr><td>Density$(map(gaskeypad) do gkp @htl("<td>$(string(last(gkp)))") end)</tbody></table>""")

# ╔═╡ 3ad2e477-ba5b-4676-b8d3-c99fa50b91d5
render_table(gaskeypad)

# ╔═╡ Cell order:
# ╠═d8ed59fe-a83e-11eb-07f8-d7123a34a1c5
# ╠═4fd5c54c-4c43-454a-81dd-5ac7245feb85
# ╠═cf85ab67-9d99-4551-93b6-db3494a3bc73
# ╠═46f5142e-369f-4c78-8d91-cd267933d8be
# ╠═57d26450-bb69-46b2-93a8-bbee93c91436
# ╠═18f3e3ce-9fd2-4a55-bcb2-5224f4bb943c
# ╠═0dfb6258-4451-4d64-a704-9da3ec92ca07
# ╠═8aa81756-2bc5-462b-bbc2-84708fab3a70
# ╠═227003ce-5470-4f04-93a5-6cbe84a2358d
# ╠═febc72ae-2ef6-484a-bc90-5844eeedae76
# ╠═2b0c948d-2eb1-461a-8999-0802b7e9e819
# ╠═6ddca240-5797-4564-adf3-0c788c8d5a80
# ╠═d7da0c9d-8f19-49d0-b093-1d0816b6c8ed
# ╠═0774a8f8-e6d1-44dc-b186-d481b713ba00
# ╠═33743bd0-0d79-4999-b93b-4d27d448b408
# ╠═410b2fff-0ed5-4ff4-82c1-bc3409ffd0af
# ╠═1ad18652-cebf-4dfe-aa01-19834ba4b14a
# ╠═427adede-fe09-4c59-b8ca-df898c4a7924
# ╠═8358fab7-4b42-45ec-bcdf-661747dde8f7
# ╠═9612e871-4b10-4cff-9215-4da903d33d06
# ╠═59588307-aaa6-45ab-862b-759b9dfc6810
# ╠═a91d6013-be1a-406d-a193-906e0c8d5a8c
# ╠═272ce517-d9ec-4afe-b5c1-90d6c4c89899
# ╠═f10cecb7-cb04-43ee-8db0-f623d99115c1
# ╠═26cdcadc-f212-4cd2-a159-c85fc1a6574f
# ╠═483c1814-c485-4ac2-8940-8c0b3b68de59
# ╠═6f0da015-f70a-490c-b116-9b9204371f76
# ╠═28dce3a3-63b4-4496-bf99-2779d29881d9
# ╠═203f0879-875e-4a25-a7df-3c54f72ff41f
# ╠═f54cb70e-b6f1-44f8-b36c-01c6f373fa90
# ╠═b9c057ef-bff2-479d-a38f-473b7ca4c309
# ╠═f52bc4fd-1a07-4963-8efc-f67152f76152
# ╠═3ad2e477-ba5b-4676-b8d3-c99fa50b91d5
# ╠═a1b61767-02c5-4678-8cde-3b86cc30bbd7
# ╠═d2bca69d-353f-48d9-bde0-677892b3fa24
# ╠═b284dddb-2da3-4e93-9ecb-b6592dcab087
# ╠═38f19a20-e002-41b5-976c-22c6af542e89
# ╠═7f773008-ac15-44bf-a276-15ea69c7e4b1
# ╠═e76c82e9-711c-4450-89fe-d795344d2cdb
# ╠═07a816ec-d1bb-4720-b381-a493643e747d
# ╠═72b053b8-9327-4128-931e-95b3e9305571
# ╠═0d89f93b-3a3a-4c1b-8e63-552ed41c6d4c
# ╠═82eb253a-2378-4074-9dab-4cb3c7eea868
# ╠═e3f5a522-6051-489e-9f72-89d9ecfc9b2c
# ╠═bf81ee3b-dd04-4983-9f2e-b52badad1da1
# ╠═8261efb7-0e4f-4866-92cb-50b0d1052d5c
