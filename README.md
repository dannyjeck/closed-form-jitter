# Closed-form-jitter

Code for performing the interval jitter hypothesis test in a closed-form fashion.

## Publication

The code here was was developed for the following [technical report](https://arxiv.org/abs/1502.07907). 

```
@TechReport{Jeck_Niebur15c,
  Title                    = {Closed Form Jitter Analysis of Neuronal Spike Trains},
  Author                   = {Daniel Jeck and Ernst Niebur},
  Institution              = {Zanvyl Krieger Mind/Brain Institute, Johns Hopkins University},
  Year                     = {2015},
  Note                     = {arXiv:1502.07907 [q-bio.NC]},
  Number                   = {DJEN-2015.1},

  Owner                    = {niebur},
  Timestamp                = {2015.02.27}
}
```

If you use this code in your own work please cite it. 

This work was also presented at CISS 2015

D. Jeck and E. Niebur, "Closed form jitter methods for neuronal spike train analysis," 2015 49th Annual Conference on Information Sciences and Systems (CISS), Baltimore, MD, 2015, pp. 1-3.
doi: 10.1109/CISS.2015.7086908

## Files

The function to perform the interval jitter hypothesis test is [jitter_closed_form.m](./jitter_closed_form.m), which calls [generate_table.m](./generate_table.m).


[speed_test.m](./speed_test.m) and [speed_test2.m](./speed_test2.m) were used to generate the figures in the technical report.
