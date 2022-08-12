# To Do: Expected
In this file I will list anything that needs to be changed in the Expected
files that I haven't already done.

# spline_model_monthly.R

* Change directories of files being read in to point to the appropriate directory
* Make sure I am using the correct directories for saving output files
* I still have in the script code that computes expecteds with nb variation from the gam model. If I take this code out, since I set the seed at the beginning of the script, the output will be slightly changed. Is it alright to make this change?
* I do not compute Sweden and Germany expecteds with linear annual trends. Should I?

# spline_model_annual.R

* Change directories of files being read in to point to the appropriate directory
* Make sure I am using the correct directories for saving output files
* I still have in the script code that computes expecteds with nb variation from the gam model. If I take this code out, since I set the seed at the beginning of the script, the output will be slightly changed. Is it alright to make this change?

# multinomial_model.R

* Change directories of files being read in to point to the appropriate directory
* Make sure I am using the correct directories for saving output files
* I still have in the script code that computes expecteds with nb variation from the gam model. If I take this code out, since I set the seed at the beginning of the script, the output will be slightly changed. Is it alright to make this change?

# Temperature data

The temperature data that I use in the multinomial_model.R script is the output
of temperature_clean.R. This processes 
WilliamData/multinomial_model/temperatures.csv and adds in data for Niue (as
described in temperature_clean.R). The output temperature_data_clean.Rdata needs
to be in another directory.