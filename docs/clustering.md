Clustering:
- always use sizes, not copy numbers (more separation in numbers)
- GMM will try `--max_num_clusters` + 4 components (higher than target number of clusters sometimes needed for optimum results)
- find best GMM cluster using AIC
- remove clusters of size <= 5 (`--max_bad_cluster_size`) where difference between mean and median differ by 10% of either so that large alleles with less supports are kept
- merge clusters separated by 10 (`--min_cluster_d`)
- merge singleton to clusters when the singleton is within 10% of closest member of cluster
- remove singletons
- merge clusters so that final clustering has `--max_num_clusters` - closest clusters merged first
- remove clusters less than `--min_cluster_size` unless median cluster size >= 5kb (`--max_check_size`)
-- `--min_cluster_size` can be fractional (<1.0) or integers

