COMP3456 - Assignment 1: Gibbs sampler for Motif Finding
========================================================

Hi Mike.

It's me, Alex's Gibbs sampler.

Usage
-----
You can run me like this:

mike@mikebox:$ (python|pypy) sampler.py file l k n

         file = <input file containing DNA matrix> 
         l = <length of motif to search for> 
         k = <number of identical results needed to acquire "convergence">
         n = <number of times to re-run the algorithm from a random starting point>

Don't forget those parameters. It feels all weird when someone runs me without them. Or worse, with extras. Ew.
