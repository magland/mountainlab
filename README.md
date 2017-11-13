# MountainLab Spike Sorting Software

## About

MountainSort (a component of MountainLab) is spike sorting software developed by Jeremy Magland, Alex Barnett, and Leslie Greengard at the Center for Computational Biology, Flatiron Institute in close collaboration with Jason Chung and Loren Frank at UCSF department of Physiology. It is part of MountainLab, a general framework for data analysis and visualization.

MountainLab software is being developed by Jeremy Magland and Witold Wysota.

The software comprises tools for processing electrophysiological recordings and for visualizing and validating the results.

Contact the authors for information on the slack team for users and developers.

**Announcement:** We have a planned new release of the software for Thursday, November 16th. See below for more information.

## Installation

[Installation instructions](old/doc/installation.md)

## How to run spike sorting

[The first sort](old/doc/the_first_sort.md)

## Working branches

* Announcement: We have a planned new release of the software for Thursday, November 16th. You can explore a preview at https://github.com/flatironinstitute/mountainsort

* [ms3](https://github.com/magland/mountainlab/tree/ms3) - development branch with the ms3 processing pipeline (preferred)

* [2017_06 branch](https://github.com/magland/mountainlab/tree/2017_06) - snapshot with only critical bug fix updates

## Some demo videos

* [Installation demo and introduction to Mountainlab](https://www.youtube.com/watch?v=P-WqvIvmx84) -- skip to the end of the video to see the example sort and the GUI. Also describes how the software is organized and some of its philosophy.

* [Demo of additional WIP GUI for viewing very large datasets -- spikeview](https://www.youtube.com/watch?v=z1V1di8sQOI)

## Tests

[Repo of unit (and not so unit) tests](https://github.com/magland/mountainlab_tests)

## Data formats used in MountainLab

[The .mda file format](old/doc/mda_format.md)

## Data management

[The .prv data management system](old/doc/prv_system.md)

## Automated curation/annotation

Documentation on using annotation scripts will be forthcoming.

Because one of the goals of mountainsort is to enable reproducible spike sorting, we strongly advise against manual modifications that go beyond merging bursting clusters and perhaps rejecting certain noise clusters.  Instead, we suggest that you export the cluster metrics along with the sorted clusters and then set cutoffs for inclusion of data in analyses based on those metrics.  This will make it easy to describe your subsequent analyses and easy to determine how those analyses do or do not depend on cluster quality.

The isolation and noise overlap metrics we provide work well for the situations we have focused on, but they can be supplemented or replaced by other objective metrics as needed. Such metric processors may be included in the pipeline as post-processing plugins as C++, matlab, or python modules. Contact us if you you would like to contribute additional cluster metrics, or need help with integration.

## Other documentation

[A guide to using MountainSort with snippets, rather than continuous data acquisition](https://github.com/mari-sosa/Mountainsort_for_snippets/blob/master/mountainsort_for_snippets.md)

An old guide: [Cluster metrics and automated curation](old/doc/metrics_automated_curation.md)

## Acknowledgements

Thanks to all the users on the slack team for ongoing testing and feedback.

## References

[Barnett, Alex H., Jeremy F. Magland, and Leslie F. Greengard. "Validation of Neural Spike Sorting Algorithms without Ground-truth Information." Journal of Neuroscience Methods 264 (2016): 65-77.](http://www.ncbi.nlm.nih.gov/pubmed/26930629) [Link to arXiv](http://arxiv.org/abs/1508.06936)

Magland, Jeremy F., and Alex H. Barnett. Unimodal clustering using isotonic regression: ISO-SPLIT. [Link to arXiv](http://arxiv.org/abs/1508.04841)

