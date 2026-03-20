# How to choose your index parameters

## The size of the partitions

### Constraints and impact of parameters

A partition represents one or several columns of REINDEER2's index matrix-like representation. In practive, its size is affected by the number of datasets indexed, the number of abundance levels, the size of the Bloom filters and the number of partitions.

All these parameters are linked by the following formula:

$
partition\_size = (nb\_datasets*levels*filters\_size)/nb\_partitions
$

However, due to implementation restriction, **the maximum size of a partition is always $2^{32}$ bits**.


### Indexing with REINDEER2

When indexing a collection, REINDEER2 starts by checking if the number of partitions is high enough to handle the size of the index.

If not, REINDEER2 will set by default the smallest number of partitions that will work with the parameters. However, since the number of partitions have been chosen to fit the parameters (including the number of datasets), adding new datasets won't be possible.


### Merging indexes with REINDEER2

When merging indexes, REINDEER2 checks that the parameters used for both the indexes were the same. It also checks that the number of partitions is high enough to handle the new size of the future merged index.

If not, the merge won't be available. It is not yet possible to split the indexes into a larger number of partitions in order to allow the merge.



## How to choose your parameters?

### Expected number of datasets

For small indexes (under 500 datasets), high values (bf=32, levels=255) of parameters won't be an issue.

In future upgrades, a parameter will be added to indicate to REINDEER2 the expected number of datasets that will be handled by the index (now or with future index merges), allowing it to propose an appropriate number of partitions suiting the parameters.

Until then, it is important to determine beforehand the expected number of datasets in order to choose the right number of partitions.



### Bloom filters size

The Bloom filters size is highly correlated with the false positive rate at query time. Depending on the average number of distinct k-mers in your datasets, you can find the right Bloom filters size in log2 (between 26 and 34) suiting your expectations in terms of average false positive rate.


| Bloom filters size |   0.10 |   0.25 |   0.50 |   0.75 |      1 |   1.25 |   1.50 |
| ------------------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| 26                 | 77.47% | 97.59% | 99.94% |   100% |   100% |   100% |   100% |
| 28                 | 31.10% | 60.60% | 84.47% | 93.88% | 97.59% | 99.05% | 99.63% |
| 30                 |  8.89% | 20.77% | 37.23% | 50.27% | 60.60% | 68.78% | 75.27% |
| 32                 |  2.30% |  5.65% | 10.99% | 16.02% | 20.77% | 25.25% | 29.48% |
| 34                 |  0.58% |  1.44% |  2.87% |  4.27% |  5.65% |  7.02% |  8.36% |



### Number of levels

The number of abundance levels determine directly the average approximation ratio between discretized abundances and real abundances values.

You will find in the following figure the average approximation ratio associated with a range of numbers of abundance levels.


![Abundance levels and approximation ratio](https://github.com/Yohan-HernandezCourbevoie/REINDEER2/tree/dev/doc/levels_and_approximation.pdf)



## How to determine the right number of partitions?

From the previous relations between our parameters, we can calculate the appropriate number of partitions using the following formula :

$
nb\_partitions >= (nb\_datasets*levels*filters\_size)/ 2^{32}
$

Keep in mind that the number of partitions calculated this way is the minimal number of partitions for these parameters. While the number of levels and the Bloom filters size are fixed after the indexing step, the number of datasets can still increase and these potential changes must be taken into account in the calculation.



