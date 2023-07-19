# step-imp
Small collection of code used to generate stepped-impedance filter designs for microwave applications.
### Current features
- Microstrip filter designs
- Supports Chebyshev type 1 and Butterworth filter designs, upto the 9th order
- Plots expected transfer function
- Calculates g values
- Calculates inductance and capacitance of each element
- Calculates lengths and widths of elements
- Exports dimensions to Ansys HFSS for electromagnetic simulations
### Planned features
- Better Ansys HFSS integration
- Better generalization of length calculations for higher orders to account for shunt susceptance and series reactance.
- Three-dimensional filter designs not using microstrip?