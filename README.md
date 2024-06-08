# hodgkin-huxley
Some code to simulate an axon's response to an applied stimulus current 
based on the Hodgkin-Huxley equations as described in the 1952 paper 
(cited below). Allows for visualization of how the gating variables m, h, 
and n influence ionic currents and membrane voltage, how the gating 
variables change based on their respective rate constants a and B, and 
how those rate constants are dependent on the membrane voltage.

Provides code to simulate an action potential using the moden convention 
(resting potential at -60mV and positive ions into the cell = voltage gets 
more positive) instead of the classic HH convention because that's what we're used to. 
(-60mV is used to be consistent with what is shown in "The Annotated Hodgkin and Huxley")

I then replicate the graphs from the Hodgkin-Huxley paper using both the 
modern convention and the original Hodgkin-Huxley convention (resting 
potential at 0mV and sign of voltage is flipped)

In "The Annotated Hodgkin and Huxley" (linked below), the modern versions 
of the graphs should match what is simulated here using the modern convention. 
I noted the page numbers in "The Annotated Hodgkin and Huxley" for each of the graphs 
so you can easily compare them to the text and mess with how changing different parameters would affect the output.

HODGKIN AL, HUXLEY AF. A quantitative description of membrane current and 
its application to conduction and excitation in nerve. 
J Physiol. 1952 Aug;117(4):500-44. doi: 10.1113/jphysiol.1952.sp004764. 
PMID: 12991237; PMCID: PMC1392413.

The Annotated Hodgkin and Huxley: A Reader's Guide 
https://press.princeton.edu/books/paperback/9780691220635/the-annotated-hodgkin-and-huxley
