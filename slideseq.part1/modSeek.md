**CHAPTER 4**

**AN ACCURATE SIMILARITY MEASURE AND AGGREGATION APPROACH BASED ON
CORRELATION RANK**

**ABSTRACT**

In this chapter, we develop RBP (rank-based Pearson), a novel rank-based
framework for aggregating thousands of coexpression rankings based on
rank-transformed similarities of Pearson correlations. Previously, gene
correlation search has used Pearson correlation values or derived
z-scores (value-based) as basis of integration. We compare RBP with
value-based approaches. We also investigate and compare two classes of
commonly used meta-analysis strategies: concatenation-based and
transformation-based approaches.

**INTRODUCTION**

The choice of similarity measure is an important consideration in the
design of similarity search system. The most common measures used by
biologists are Pearson correlation and Spearman correlation. These have
the desirable characteristics of being easy to calculate and leading to
interpretable results. More advanced measures have been proposed
recently, such as the maximal information coefficient (Reshef et al.
2011), distance correlation (Székely et al. 2007), which can detect
complex correlation beyond just linear relationships, but the
disadvantage with these measures is that they can be difficult to
interpret (de Siqueira Santos et al. 2013), since non-linear
relationships may not always find an explanation in biological context.
These measures (MIC, MI) however are interesting in a theoretical point
of view (de Siqueira Santos et al. 2013).

Recently, there have been several efforts (Obayashi & Kinoshita 2009;
Ruan et al. 2010; Kolde et al. 2012) to utilize a rank-transformed
similarities of Pearson correlations (rank-based measures) for
clustering genes from expression datasets. These efforts are motivated
by the fact that the distribution of Pearson correlation values,
especially in gene clustering applications, can often be a skewed
distribution owing to hubby genes (Han et al. 2004) in the network, or
well-connected genes. Hubby genes cause the node degrees in a
coexpression network to be uneven with hubby genes possessing very large
node degrees in yeast coexpression network (Ruan et al. 2010). Hubby
genes may gather to form large dense clusters, which prevent the
formation of smaller clusters and produce many singletons (Ruan et al.
2010). Rank-based approach sought to solve this bias by considering the
rank of Pearson correlation in each gene's correlation vector. As the
correlation rank vector will be normalized to the same range (rank=1...
*N*), all genes' correlation rank vectors will be comparable to each
other and thus normalized. Some approach performs an additional
processing step, such as converting the edges of top ranked X neighbors
of a gene to 1 and ignoring the rest of edges (or assigning 0 to rest of
genes) (Ruan et al. 2010), which is used for constructing a binary,
unweighted coexpression network. The rank-based approach has been
successfully used for module finding in coexpression network (Ruan et
al. 2010). It was noted that the use of rank of correlation leads to a
better fit to the property of a scale-free network when the degree of
nodes is analyzed, than a Pearson value-based network (Ruan et al.
2010). Additionally, the rank of Pearson has been demonstrated to work
better than Pearson on single expression datasets in retrieval of
functionally co-annotated genes (Obayashi & Kinoshita 2009). Within
rank-based measures, the mutual rank similarity measure is a popular
approach first mentioned in CoxpressDB (Obayashi & Kinoshita 2011),
which computes the relative rank of two genes in the correlation vector.
Mutual rank has been frequently used in searching for functionally
co-annotated pairs in both array and RNA-seq datasets (Obayashi et al.
2013; Van Dam et al. 2015).

Though the above rank-based measures have been able to allow users
analyze the correlations on a single dataset, biologists often are faced
with hundreds of datasets. Combined analyses of multiple datasets' gene
correlations, also termed as meta-analysis, would lead to more
statistically robust biological correlations that can be tested in labs.
Meta-analysis approaches for correlation search problem are largely
divided into two classes, based on the terminologies coined in a recent
review (Ritchie et al. 2015): 1) *concatenation-based* approaches, where
one concatenates multiple datasets into a large dataset then analyzes
meta-dataset, and 2) *transformation-based* approaches, where one
analyzes datasets individually and combines the results of many analyses
(Ritchie et al. 2015). Unfortunately, although methods applying each of
two meta-analytic principles are routinely used, no comparisons have
been made so far between them especially for a large number of datasets
and for coexpression mining problem.

By combining the advantage of rank-based measures with meta-analysis, we
can deliver superior correlation search performance. Thus, we will
perform a comprehensive evaluation of the two types of approaches using
the large model organism and human gene expression compendia. We will
also evaluate the performance of rank-based correlation measures, such
as mutual rank, in large-scale correlation search involving thousands of
datasets, as this has not been done previously.

At the same time, storing the correlation matrices for thousands of
expression datasets is both costly to maintain and inefficient to
integrate. Our estimate suggests that the storage requirement for 1000
human datasets' correlation matrices is approximately 600GB! Efficient
time and space saving strategies must be developed.

In this work, we explore rank-based integration solutions for the
problem of large-scale gene correlation search. This is different from
the previously developed Pearson z-score value-based integration, SEEK.
The change is motivated by our promising results delivered by rank-based
measure in clustering and small-scale functional studies. Here we
propose a new rank-based Pearson (RBP) similarity measure and
aggregation approach for the problem of integrating thousands of
correlation gene rankings to arrive at a consensus ranking. We
illustrate the advantage of this approach over prior value-based
approach in retrieving functionally co-annotated genes for expression
datasets deposited for human and diverse model organisms. We also
illustrate the space-efficient property of RBP that achieves many folds
of size reduction over full rank correlation matrix. Importantly, the
RBP approach adopts a rank position based weighting system in the
aggregation of rank lists, where this weighting is dependent on
geometric distribution parameter. This makes flexible integration of
rank lists possible, and large space reduction possible. The efficiency
makes RBP a very appealing approach for very large-scale analysis and
on-the-fly aggregation. Additionally, we compare and investigate the two
widely adopted approaches to meta-analysis: analyze-and-merge, and
merge-and-analyze approaches.

**RESULTS**

**Large-compendia meta-analysis: concatenation-based versus
transformation-based strategies**

We first survey the different rank-based aggregations that work in
large-scale setting (with thousands of datasets). These approaches are
divided into two types:

1\) **Merge and analyze**: concatenate multiple datasets into a large
dataset then analyze once. Let
$D\  = \ \left\{ D_{1},\ D_{2},\ D_{3}\ldots \right\}$ be a collection
of datasets with $n_{1}$, $n_{2}$, $n_{3}$... number of conditions. Let
$N_{x} = n_{1} + n_{2} + n_{3} + \ldots$. Combined dataset will have
$N_{x}$ conditions. Combined analysis of this type involves performing
the analysis once on the merged dataset with $N_{x}$ conditions. In the
example of coexpression analysis, this involves computation of Pearson
correlations on the merged dataset. This approach sidesteps the need to
aggregate multiple lists since only one rank-list is generated (the
original data collection is concatenated to be seen as one dataset as a
preparation step)

2\) **Analyze and merge**: analyze datasets individually then combine
the analysis results. This is typically the scenario of most rank list
aggregations, where each dataset is analyzed to generate one ranking (or
intermediate result). Then the rankings from all datasets are aggregated
together to produce a master ranking. In the case of coexpression
analysis, Pearson correlations are calculated for each dataset in the
collection. These correlations are used to generate a gene ranking (of
similar genes with respect to a query of interest). Then the rankings
are combined across all datasets to generate a consensus ranking of
genes, describing the consensus correlation ranking.

These two types of approaches have been termed in other literature as
*concatenation-based* and *transformation-based* approaches (Ritchie et
al. 2015), and they are equally relevant in our correlation search
context. Below I name several examples of these two approaches in the
context of coexpression analysis, and mention some of their advantages
and disadvantages.

Used in prior coexpression works, the first approach has the potential
to increase statistical power of the analysis, though its shortcomings
and challenges are significant. Specifically, this approach is most
often hampered by batch effects and experimental variations, because it
is not trivial to concatenate different datasets together to form a
merged matrix. The reasons for this are several: technology differences,
experimental variations, and biases. These in turn create incomparable
expression values that bias the downstream analysis. Thus, correctly
merging matrices almost always requires prior correction of batch
effects between datasets, which is often not trivial to resolve
completely. Secondly, since the first approach merges datasets,
individual constituent datasets are assumed to have equal weight. For
example, Pearson correlations calculated on a concatenated dataset will
typically treat all conditions equally (the boundary of datasets
disappear in merged dataset). This is done because it is not possible to
know ahead of time which biological context users will query, or know
which conditions or datasets one should weight as a result. These
"unweighted" Pearson correlations may not be ideal for addressing
certain context-specific queries where specific conditions need to be
emphasized and weighted. Despite this inflexibility, the first approach
has the practical advantage of being efficient, since Pearson of a
merged dataset will be fixed and can be determined in advance. Example
of this type of approach includes CoxpressDB (Obayashi & Kinoshita
2011), GeneFriends (Van Dam et al. 2015) and aggregated dataset
correlation (compared in SPELL (Hibbs et al. 2007)).

The second type of approach solves many problems of the first type, and
allows flexible query-dependent weighting of datasets. This approach
performs analysis on individual datasets since conditions from within a
dataset are most comparable to each other. This avoids batch effect
issues in merging datasets directly, but rather the results of analyses
are combined. Examples of this approach included Pearson z-score-based
methods SEEK (Zhu et al. 2015), SPELL (Hibbs et al. 2007), rank-based
methods that aggregate rank-lists: MEM (Adler et al. 2009), Robust Rank
Aggregation (Kolde et al. 2012), RankProd (Hong et al. 2006). In
agreement with our expectations, our prior work on SEEK shows that SEEK
outperforms the aggregated dataset correlation (i.e., the first
approach) by many fold of magnitude in retrieval of functionally
co-annotated gene.

**Pearson correlation value versus rank**

The second comparison that will be made is between Pearson correlation
values and ranks. Pearson correlation values, though widely used, can be
severely affected by the problem of hubby genes. These hubby genes
usually attract a large number of other genes to form large, dense
clusters of genes (or sometimes called hubs) in the network (Ruan et al.
2010), since these hubby genes have the highest Pearson correlation
values and other genes' correlation to the hubby genes tend to be very
high. As a result, the large clusters formed from Pearson correlation
values very often prevent weakly correlated nonhubby genes from forming
their own smaller clusters. To resolve this problem, rank-based approach
presents an appealing alternative to value-based Pearson correlation.
Rank-based approach normalizes Pearson correlation values on a per gene
basis: for each gene *g*, it ranks all other genes based on correlation
with *g*, assigning genes with rank of 1 up to *N*, the number of genes.

We propose RBP (Rank-based Pearson), a novel, simple rank-based
similarity measure and aggregation method for coexpression rank-lists.
RBP differs from other rank-based approach in that RBP relies on the
exponential transformation on the rank of Pearson correlation values,
followed by gene hubbiness normalization (**Methods**). The geometric
parameter has the desirable effect of emphasizing the items at the top
of the list, thereby achieving rank position-based weighting
(**Methods**). This method has substantial benefit in large-scale
coexpression rank-lists aggregation as we will show in later sections.
We compare RBP to mutual rank integrated with RankProd algorithm and,
separately also, to value-based approaches (**Table 1**). **Table 1**
summarizes all of the comparisons that we will make.

Employing exponential weighting on ranks has been previously used as an
effective feature selection method in machine learning applications
(Haury et al. 2011). Haury et al has used an exponentially decreasing
function of the rank when aggregating ranks of features across multiple
samples (this area of application is called ensemble feature selection).
In that study, a similar problem was encountered, where small
perturbations of the training data result in unstable signatures (Haury
et al. 2011). Rank aggregation was proposed to resolve this problem by
"stabilizing" variable selection (Haury et al. 2011). In another area of
application, Rank-Biased Precision uses a similar
geometric-parameterized version for robust measurement of retrieval
effectiveness (Moffat & Zobel 2008).

**Table 1**: Overview of search methods compared in this paper and in
our previous paper

| Method                          | Description of similarity measure and integration method                                                                                                 | Type of integration strategy | Is value-based? | Is rank-based? |
|---------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------|-----------------|----------------|
| Agg dset MR (COXPRESdb \[1\])   | Concatenated dataset's mutual rank                                                                                                                       | Merge and analyze            |                 | Yes            |
| Agg dset cor                    | Concatenated dataset's Pearson correlation                                                                                                               | Merge and analyze            | Yes             |                |
| MR \[2\] + RankProd \[3\]       | Mutual rank for individual dataset. Aggregation of rank-lists across datasets by computing rank geometric mean (RankProd). Rank-list is *weighted*\*.    | Analyze and merge            |                 | Yes            |
| SEEK \[4\]                      | Pearson correlation z-score values per dataset. Integrate all datasets' Pearson correlation values by computing *weighted*\* mean.                       | Analyze and merge            | Yes             |                |
| SEEK with RBP (proposed)        | Pearson correlation ranks per dataset. Aggregation of rank lists across datasets with an exponential emphasis on rank position. Rank-lists are weighted. | Analyze and merge            |                 | Yes            |
| MEM \[5\](order statistics)\*\* | Pearson correlation ranks. Aggregation of rank lists by binomial probability distribution.                                                               | Analyze and merge            |                 | Yes            |

\[1\] COXPRESdb: a database to compare gene coexpression in seven model
animals. Nucleic Acids Res. 2011 Jan; 39: D1016-22.

\[2\] Rank of Correlation Coefficient as a Comparable Measure for
Biological Significance of Gene Coexpression. DNA Research. 2009.

\[3\] RankProd: a bioconductor package for detecting differentially
expressed genes in meta-analysis. Bioinformatics. 2006.

\[4\] Targeted exploration and analysis of large cross-platform human
transcriptomic compendia. Nature Methods. 2015 Mar; 12(3):211-4.

\[5\] Mining for coexpression across hundreds of datasets using novel
rank aggregation and visualization methods. Genome Biology. 2009.

\*Under weighted integration, we calculate weighted average correlation,
$
\frac{1}{\sum_{d}^{}w_{d}}\sum_{d}^{}{w_{d}r}$ (for value-based), or
weighted average rank,
$\left( \prod_{d}^{}\text{Rank}^{w_{d}} \right)^{\frac{1}{\sum_{d}^{}w_{d}}}$
(for rank-based). Weight is assigned on the basis of relevance of each
dataset to the query genes, as used in SEEK.

\*\*Compared in Zhu et al Nat. Methods paper, so comparison is skipped
in this paper.

**Large-scale functional evaluations illustrate the advantage of RBP**

Our evaluation set-up tests each approach's ability to retrieve
co-annotated genes from the GO biological process annotations. To
evaluate the effect of the rank-list based correlation search, we tested
the search approach on diverse queries that are drawn from GO slim
biological process annotations. Each GO slim gene set defines a set of
co-regulated genes in a slim biological process or pathway (these are
called slim because they represent a level in the GO hierarchy that is
neither too general nor too specific). For each GO slim term, we divide
the genes into two portions, one which is used to choose our queries,
and the other which is used for evaluating the results of retrieving the
queries. The queries are randomly drawn to ensure they are diverse and
unbiased.

The new approach RBP has been implemented in the SEEK system and
conveniently labeled as "SEEK with RBP" in our comparisons. First, we
compared the new measure SEEK with RBP against the previous SEEK with
Pearson z-scores for each of the 5 model organisms and human, keeping
other factors such as the compendium composition, search algorithms
constant. The results suggest that RBP consistently outperforms Pearson
z-scores in the SEEK framework with relative performance improvement at
10% recall ranging from 112% (fly, yeast, mouse) to 130% (human), and
even at larger 20% recall the fold improvement of RBP remains at
117-140% (**Figure 1**). The improvement is not the result of few
"outlier" GO slim terms that achieved massive gain, because when we
switched to using log-average as performance aggregator, the fold
difference between the two is similarly large ranging from 113% (mouse),
to 121% (fly, yeast) and to 130% (worm). Such a result shows the
superiority of rank-based method such as RBP over value based method
especially in the context of integrating thousands of coexpression
rankings.

<img src="media/image1.png" style="width:5.99981in;height:4.06479in" />

**Figure 1**: rank-based vs. value-based approach.

Next, we compared the RBP approach against the previously widely adopted
merge-and-analyze meta-analysis (i.e. COXPRESSdb (Obayashi & Kinoshita
2011)), still keeping the comparison rank-based. We used RBP to first
compute coexpression ranking of each dataset and combined them via the
query-dependent weighting. For the merge-and-analyze, we calculate
Pearson mutual rank on the merged dataset (i.e., all datasets combined
together) as is done previously (see COXPRESSdb, this is named "agg dset
MR" in **Figure 2**). In such a way, it treats all arrays in the merged
dataset equally. Thus, the comparison would also demonstrate the
difference between a weighted (the former) and an unweighted Pearson
approach (the latter). From **Figure 2**, we clearly witness a rather
large performance gap (up to 200%) between the two approaches in terms
of functional evaluation. One important factor playing to the
disadvantage of the merge-and-analyze approach is the unweighted nature,
which permits a large number of irrelevant, bad quality arrays within
the merged dataset to affect the Pearson correlation values. Of course,
unweighted correlations are quite informative for general, housekeeping
biological processes which are expected to be present in all arrays,
such as rNRA processing, zinc ion transport, fatty acid oxidation,
meiosis I, cell cell recognition, cilium morphogenesis --the
fold-over-random values are at 649, 258, 211, 133, 124, and 121
respectively at 1% recall (Supplementary Figure X). However, by using
flexible query-based correlation weighting, the analyze-and-merge
approach is allowed to deliver superior functional search performance
for a much more diverse set of GO terms -- in 289 of the 407 GO slim
biological processes (in mouse) the weighted approach wins
(Supplementary Figure X).

<img src="media/image2.png" style="width:6in;height:4.29217in" />

**Figure 2**: merge-and-analyze (or concatenation-based, Agg Dset MR)
vs. analyze-and-merge (or transformation-based, SEEK with RBP) approach.

**Robust over various compendium sizes**

The model organism expression collections offer us a good opportunity to
assess the robustness of the performance over different number of
datasets, different number of genes, and datasets generated by various
array technologies (as listed in Table X). In terms of the robustness,
we note that the performance gain of RBP over value-based z-score
integration holds consistently over diverse organisms coming from
various clades and genome size (**Figure 2**). This shows the general
applicability of RBP. It should also be noted that some organisms'
compendium may require special tuning of parameter p to achieve optimal
performance. For example, C. elegans' RBP performance is optimal at
p=0.999, while at p=0.99 it is within only 1% difference from SEEK. For
mouse, fly, human, yeast, we chose the p=0.99. For D. rerio, we chose
p=0.999 similar to C. elegans. Though the parameter p seems to depend on
the size of the genome, this is not always the case (as it may also
depend on the judgment standards that seem to prefer results of certain
depth in the ranking). Our recommendation is to train the system with
sample queries to see what parameter value works best within the
possible choices (p=0.99, 0.995, 0.999).

**Space-efficient property**

Due to the rank-position weighting, RBP can utilize a sparse correlation
rank matrix for each dataset which avoids the need to store lowly ranked
genes (as their rank-position weights are negligible). By doing so, we
were able to dramatically reduce the storage size of the correlation
matrix (**Table 2**). The fold reduction in storage space compared with
full correlation rank matrices ranges from 3.2-fold in D. rerio to as
much as 5.3-fold in M. musculus, given U=3000 for all organisms (U
represents the depth of the ranking to preserve). Importantly, this
reduction does not prevent the search algorithm from achieving accurate
search performance. We compared with the ideal case where complete
rank-lists have been integrated using the state-of-the-art RankProd
approach. When we integrate coexpressed gene rankings derived from
sparse rank matrix across datasets using RBP, the performance is between
-3% and +8.1% away from the state-of-the-art RankProd approach which
integrates in the presence of full gene rank information (**Figure 3,
Table 2**). In fact, in all but one case, we even observe a positive
gain in the performance from using RBP. These results together
illustrate that storage saving does not need to sacrifice quality of the
consensus ranking, because RBP is able to effectively utilize partial
gene rankings to obtain good performance.

**Table 2**: Space-efficient property of RBP. Uncompromised search
performance.

<img src="media/image3.png" style="width:6.32374in;height:2.18858in" />

**CONCLUSIONS**

With thousands of genome-wide expression datasets now deposited,
meta-analysis becomes more important now than ever. The proposed RBP
rank-list integration provides a general framework for integrating
rank-based results from hundreds, or thousands of studies. RBP is
expected to have usage beyond coexpression rankings, such as in
differentially expressed gene integration. The space-saving advantage of
RBP will be important in search applications where the number of
features (such as genes) make it no longer practical to calculate and
store the full correlation matrix, such as in exon or transcript-based
correlation search.

<img src="media/image4.png" style="width:6.41573in;height:4.43473in" />

**Figure 3**: Mutual Rank (MR) + RankProd vs. SEEK with RBP,
illustrating the small difference between using full rank matrix
(former) and partial or shaved matrix with rank-position emphasis
(latter).

**Methods**

**Defining similarity measure**

A correlation matrix is defined as a $G \times G$ matrix, where each
entry in the matrix is the correlation, $r$, of gene $g_{x}$ and
$g_{y}$. First, r is symmetric ($r\left( x,\ y \right) = r(y,x)$). We
first perform row-wise rank transformation of matrix to generate a
Pearson rank-based matrix where $rank(x,\ y)$ describes the rank of y
with respect to x, when the entire set of correlation values in the row
x is ranked (i.e., the values \[$r(x,\ g_{1})$, $r(x,\ g_{2})$,
$r(x,\ g_{3})$...\]. Note that this matrix is asymmetric. Next, we
define:

Symmetric rank measure:

$$\text{minimal\ rank}\left( x,\ y \right) = \min\left( \text{rank}\left( x,y \right),\ rank\left( y,\ x \right) \right)$$

We introduce a rank position weighting system, which converts the
minimal rank matrix to a geometrically-based rank matrix (RBP):

$$\text{RBP}\left( x,y \right) = \left( 1 - p \right)p^{\text{minimal\ rank}\left( x,y \right) - 1}$$

We define approximated RBP, $\text{RB}P^{a}$, whose purpose is to
"shave" the matrix RBP according to the top depth U:

$$\text{RB}P^{a}\left( x,y \right) = \left\{ \begin{matrix}
\text{RBP}\left( x,\ y \right)\text{\ \ if\ minimal\ rank}\left( x,\ y \right) \leq U \\
0\ \ if\ minimal\ rank\left( x,y \right) > U \\
\end{matrix} \right.\ $$

Where $U \approx 0.10 \times \left| L \right|$, the size of the list;
this $U$ represents the top depth of the rank list a user would like to
preserve, here 10% of the list size. The idea is that only the top U
genes are assumed to be important. But the weight on the rank position
should decrease gracefully. This $U$ and the parameter p both control
the weight on the rank position. If $U$ is determined using the
aforementioned guideline, then p should be appropriately set so that for
rank position $U$ or beyond, $\left( 1 - p \right)p^{U - 1}$ would be
infinitesimally small (i.e., zero), so that it fits the description of
$\text{RB}P^{a}\left( x,y \right) = 0$ in equation X. One reasonable
recommendation is to set the parameter p such that
$\sum_{t \in \left\lbrack 1,\ U \right\rbrack}^{}{\left( 1 - p \right)p^{t - 1}}$
is reasonably close to 1 (loosely defined as 0.95-0.99). For our
organisms, we found that the following parameters satisfy this
recommendation: (U=500, p=0.99), (U=3000, p=0.999). To represent
$\text{RB}P^{a}$ in terms of data structure, we choose a sparse matrix
representation since most entries are zero.

Finally, we adopt the following normalization procedure (hubbiness
correction):

$$\text{RB}P^{*} = D^{- \frac{1}{2}}\left( \text{RB}P^{a} \right)D^{- \frac{1}{2}}$$

Where $D$ is a diagonal matrix with (i, i) element equal to the sum of
i-th row in $\text{RB}P^{a}$.

The above procedure and the selection of $U$ and $p$ parameters allow us
to create a sparse, rank-based representation of the original Pearson
correlation matrix. This sparse matrix $\text{RB}P^{*}(x,\ y)$ (defined
in equation X) importantly has between $U$ and $2U$ number of entries,
which compared to the original full dimension (with $G^{2}$ entries),
has achieved an enormous reduction in space, since zero-values in
$\text{RB}P^{*}$ are not stored.

**Search system with rank-based integration: RBP**

We first obtain rank-based RBP\*(x, y) matrix for each dataset d,
RBP\*d(x, y). The goal is to retrieve the query's surrounding genes
using rank-based similarity matrix RBP\* that is defined per dataset.
Thus given query gene-set Q, we calculate an aggregated gene rank per
dataset as:

$$R_{d}\left( Q,\ g \right) = \sum_{q \in Q}^{}{\frac{1}{|Q|}\text{RB}{P_{d}}^{*}(q,g)}$$

We linearly combine $R_{d}\left( Q,\ g \right)$ over all datasets by
multiplying individual dataset's ranking with its corresponding dataset
weight w:

$$R_{C}\left( Q,\ g \right) = \frac{1}{\sum_{d \in D}^{}w_{d}}\sum_{d \in D}^{}{w_{d}R_{d}(Q,\ g)}$$

The resulting rank-list R~C~ is the consensus gene-ranking describing
query's correlated genes. The weight w~d~ is the query-based dataset
weighting (published in SEEK, the idea is to use coexpression of query
genes to weight how informative each dataset is).

**Comparison with value-based approach: Pearson correlation z-scores**

We perform the same query-based search algorithm as used for RBP, but
with similarity measure changed to Pearson correlation z-scores (with
hubbiness correction), defined as in our previous paper.

**Comparison with rank-based integration: mutual rank and RankProd**

For the purpose of comparing against the integration when complete gene
ranking is available, we used Mutual rank coupled with RankProd, which
is a rank-list aggregation approach that works with full ranking
information available. First, mutual rank matrix is defined by:

$$MR(x,y) = \sqrt{\text{rank}\left( x,y \right) \times rank(y,x)}$$

$$\text{RP}_{d}\left( Q,\ g \right) = \left( \prod_{q \in Q}^{}{MR_{d}\left( q,g \right)} \right)^{\frac{1}{\left| Q \right|}}$$

$$\text{RP}_{C}\left( Q,\ g \right) = \left( \prod_{d \in D}^{}{\text{RP}_{d}\left( Q,\ g \right)^{w_{d}}} \right)^{\frac{1}{\sum_{d \in D}^{}w_{d}}}$$

RP~C~ describes the weighted consensus ranking. Similar to above, the
weight w~d~ is the query-based dataset weighting used in SEEK.
