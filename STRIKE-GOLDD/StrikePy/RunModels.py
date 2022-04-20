"""
This file is used to run the different models. 
Simply call the function 'strike_goldd()' with no arguments (the default options file 'options.py' will be used) or by passing
the name of a custom options file as the only argument. 

An example of both calls is given below, just run this file to perform the test.
 ***Only the function call without arguments is analysed (delete the # before the second call to analyse the HIV model)
"""

from strike_goldd import strike_goldd

# Without arguments...
strike_goldd()

# Passing to it the options of the predefined model HIV...
#strike_goldd('options_HIV')

