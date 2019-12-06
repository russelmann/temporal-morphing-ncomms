//////////////////////////////////////////////////////////////////
// Configuration //////////////////////////////////////////////////

// Main structure parameters
spring_th = 0.4;  // spring thickness, mm
d         = 8.0;  // distance between the bases, mm (use 0 to print bulk)

// Optional enhancements
thrd_h    = 33;

// Other structure parameters
thickness = 5.0;  // thickness of the structure, mm (for 2 parts)
h         = 1.0;  // spring height, mm
alpha     = 37;   // spring angle, deg
base_d    = 2.0;  // base diameter, mm
conn_a    = 2;    // inner connector from base center
conn_b    = 4;    // outer connector from facet
thrd_len  = 4;    // length of threshold

// Handle dimensions, mm (height depends on the structure)
total_length = 40;
handle_width =  5;

// Alignment parameters, mm
male_radius   = 1.00;

// Printing gaps
female_radius = 1.04;
gap = 0.05; // for adaptor

///////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

include <inputs.scad>

include <inputs.scad>;
