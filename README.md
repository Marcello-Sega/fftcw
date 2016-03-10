### Welcome to the GitHub Page of fftcw.
I am using this code mostly for my personal purposes, and it is still full of tentative/suboptimal code, but I thought that anyway somebody might find it useful. I will soon strip it down to a more usable format and, provide some usage guidelines. If you happen to be interested in this code, drop me an email, that might be the driving force to bring this project to a "publishable" state.

### Simple usage example

If you have just one scalar value of which, you want to compute the time autocorrelation, and the time series is stored as a single column in a file (```series.dat```) at intervals of time of 0.001 units, then you can call

    fftcw -t 0.001 -n 1 -d 1 -T series.dat -o correlation.dat

The result will be stored in ```correlation.dat```, with the first two columns specifying the indices of the quantities being correlated (here, 0 and 0), followed by the time lag in the third column, and eventually then correlation function.


### Go back to [my Homepage on GitHub](http://marcello-sega.github.io/)
