*****************
Terminology (WIP)
*****************

Protein *

Sidechain *

Residue

    * Each residue has many attributes
        * Conformation Angles (Angles of Rotation)
            * phi, psi, ome (omega), omep (omega prime)
        * Bond Angles
            * a1, a2, a3, a4, a5, a6, a7
        * Bond Lengths
            * L1, L2, L3, L4, L5
        * Secondary Structure
            * ss
        * Chi Angles, properties vary per Residue
            * chi1, chi2, chi3, chi4, chi5
        * bm
        * bs
        * bg
        * h_bond_energy
        * zeta
        * terminal_flag
        * xpr

Plots

    * In the plot each hit is sorted into a box based on the values of its attributes listed as "X Axis" and "Y Axis". Any attribute can be chosen as X or Y. The "Plotted Attribute" of all the hits that end up in the same box are averaged (and the standard deviation is calculated). If no hits fall in a box there is no meaningful average so the box is left blank. The occupied boxes are colored based on the average value in that box.

    * The exception is the default, when the "Plotted Attribute" is "Observations". Then the color is simply calculated from the number of hits that end up in the box, with zero being black. Regardless of the "Plotted Attribute" the same boxes should always be black.

Search

    * Returns a set of residues
    * When the "Plotted Attribute" is "Observations" the color is simply calculated from the number of hits that end up in the box, with zero being black. Regardless of the "Plotted Attribute" the same boxes should always be black.
