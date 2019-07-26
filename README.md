---
author: Nicola F. Müller
level: Intermediate
title: Estimating Reassortment Networks with CoalRe
subtitle: Estimating Reassortment Networks with CoalRe
beastversion: 2.5.2
---


# Background

Phylogenetic trees are often used to describe the history of genetic sequences. 
There are however several processes that recombine genetic material.
Different parts of genome can then code for different histories and the shared history of the full genome can not be represented anymore by a tree.
Reassortment is such a process that can reshuffle segments when people are infected by multiple viruses.
In order to allow inference in the presence of reassortment, we introduced the coalescent with reassortment and a MCMC framework to infer ressortment networks, the embedding of segment trees and evolutionary parameters.

----

# Programs used in this Exercise

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

BEAST2 ([http://www.beast2.org](http://www.beast2.org)) is a free software package for Bayesian evolutionary analysis of molecular sequences using MCMC and strictly oriented toward inference using rooted, time-measured phylogenetic trees. This tutorial is written for BEAST v{{ page.beastversion }} {% cite BEAST2book2014 --file CoupledMCMC-Tutorial/master-refs %}.


### BEAUti2 - Bayesian Evolutionary Analysis Utility

BEAUti2 is a graphical user interface tool for generating BEAST2 XML configuration files.

Both BEAST2 and BEAUti2 are Java programs, which means that the exact same code runs on all platforms. For us it simply means that the interface will be the same on all platforms. The screenshots used in this tutorial are taken on a Mac OS X computer; however, both programs will have the same layout and functionality on both Windows and Linux. BEAUti2 is provided as a part of the BEAST2 package so you do not need to install it separately.

### TreeAnnotator

TreeAnnotator is used to summarise the posterior sample of trees to produce a maximum clade credibility tree. It can also be used to summarise and visualise the posterior estimates of other tree parameters (e.g. node height).

TreeAnnotator is provided as a part of the BEAST2 package so you do not need to install it separately.


----

# Practical: Setting up an coalescent with reassortment analysis

In this tutorial, we will describe how to set up a coalecent with reassortment analysis  

## The Data
The dataset consists of 25 HA and NA sequences of pandemic 2009 like human influenza A/H1N1 sampled between 2012 and 2017.

## Download CoalRe
First, we have to download the CoalRe package by using the BEAUTi package manager. Go to `File >> Manage Packages` and download the CoalRe package .

<figure>
	<a id="fig:example1"></a>
	<img style="width:50%;" src="figures/Package.png" alt="">
	<figcaption>Figure 1: Download the CoalRe package</figcaption>
</figure>

After the package is installed, re-start BEUAti

## Load in the sequences
In order to load in the sequences into BEAUti, they can be dragged and dropped into the paritions window.
Then a window will pop up, where we'll have to specify that all the sequence files just loaded in are nucleotide sequences

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/LoadSequences.png" alt="">
	<figcaption>Figure 2: Drag and drop the sequence files into the oartitions window.</figcaption>
</figure>

In order to account for rate variations across the different nucleotide sites, we next have to split both alignments into codon positions.
To do so, select one of the segments and the press split and select the 1,2 +3 option.
Then, repeat the same thing for the other segment.

<figure>
	<a id="fig:example1"></a>
	<img style="width:30%;" src="figures/SplitAlignments.png" alt="">
	<figcaption>Figure 2: Split alignments into codon positions.</figcaption>
</figure>

## Linking the site models and clock models 
Next, we'll have to select all partitions, and then press `Link Site Models`, `Link Clock Models`.
 
## Setting up the sampling times
Now, we'll have to set up the sampling times. 
To do so, go to the `Tip Dates` partition and select `Use tip dates`
Specify the dates format to be `yyyy-M-dd` and press the `Auto-configure` button.
In the new window, select split on character 

<figure>
	<a id="fig:example1"></a>
	<img style="width:30%;" src="figures/SplitCharacter.png" alt="">
	<figcaption>Figure 2: Split character and tke the second group to get the sampling times .</figcaption>
</figure>

In order to set the tip dates for both segments, select both partitions an press `OK`

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/CloneTipDates.png" alt="">
	<figcaption>Figure 2: Clone tip dates .</figcaption>
</figure>

## Setting up the site model

As a site model, we will use an $HKY+\Gamma_4$ model. To do so, set the site model from `JC69` to `HKY`.
Also, set the `Gamma Category Count` to 4 and make sure to click `estimate` for the `Substitution Rate`. The last part allows each segment, as well as the first two and the third codon position to have different relative rates of evolution.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/SiteModel.png" alt="">
	<figcaption>Figure 2: Setting up the site models .</figcaption>
</figure>

We next have to go back to the `Partitions` ta, select all Partitions and then press `Unlink Site Models`. This allows each parition to have its own site model with its own parameters. 

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/UnlinkSiteModels.png" alt="">
	<figcaption>Figure 2: Setting up the site models.</figcaption>
</figure>

## Setting up the Priors

We can leave the clock model as is, which means that we use a Strict Clock Model and directly go the the `Priors` tab.
The first and most important thing we have to do here is to to change the `Yule Model` to `Coalescent With Reassortment Constant Population`. 
This has to be done for both (!) trees.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/CoalescentWithReassortment.png" alt="">
	<figcaption>Figure 2: Changing the Yule model to the coalescent with reassortment.</figcaption>
</figure>

We now have to set the prior distribution on the parameters.
The prior distribution on the effective population size can be left as is, but we have to change the prior on the reassortment rate.
Set the prior distribution on the reassortment rate to be an exponential distribution with mean 0.25. This means that we assume on average one reassortment event per lineage every 4 years.
 
<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/ReassortmentPrior.png" alt="">
	<figcaption>Figure 2: Setting the prior distribution on the reassortment rate.</figcaption>
</figure>

## Setting up the Chain Length
Next, we can go to the `MCMC` panel to set the chain length. For this analysis with just a few sequences, a chain length of 5 million should be enough.
This was the last step of setting up the xml and we can now save it by going to `File > Save as`

<figure>
	<a id="fig:example1"></a>
	<img style="width:20%;" src="figures/SaveXml.png" alt="">
	<figcaption>Figure 2: Save xml.</figcaption>
</figure>

## Run the xml	
Next, open `BEAST` and run the xml. This should take about 15-20 minutes.
Alternative, the folder `precooked_runs` contains the log files of the run.

## Inspect the run in Tracer

Next, we have to check whether everything converged.
To do so, we can opern the program Tracer and load the `*.log` file.
All ESS values should optimally be above 200.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/Tracer.png" alt="">
	<figcaption>Figure 2: Check convergence in Tracer.</figcaption>
</figure>

Next, we can check what rates were inferred.
The `clockRate` denotes the average rate of evolution across all segments.
The 'reassortmentRateCwR.alltrees` denotes the rate of reassortment (from present to past) per lineage and year.


## Summarize the posterior distribution of networks
Next, we can summarize the distribution of networks by maximizing the clade credibilities.
To do so, open `BEAUti` and select `File > Launch Apps`.
Then, select `Reassortment Network Annotator`.

Next, choose the `networks.trees` file as input for the `Reassortment Network log file` and choose the file where the mcc network should be saved to and press analyse.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/Summary.png" alt="">
	<figcaption>Figure 2: Produce the maximum clade credibility network.</figcaption>
</figure>


## Visualize the network using icytree.org
Next, open your browser and go to the webpage [icytree.org](icytree.org)
The resulting mcc network file can now be drag and dropped into icytree to visualize the network.

<figure>
	<a id="fig:example1"></a>
	<img style="width:70%;" src="figures/IcyTree.png" alt="">
	<figcaption>Figure 2: Produce the maximum clade credibility network.</figcaption>
</figure>


----

# Useful Links

- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)
- Join the BEAST user discussion: [http://groups.google.com/group/beast-users](http://groups.google.com/group/beast-users)

----

# Relevant References

{% bibliography --cited --file CoupledMCMC-Tutorial/master-refs %}
