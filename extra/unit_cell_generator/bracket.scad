include <params.scad>

handle_length = (total_length - d) / 2;

// Measurement block output
l = (d/2) / cos(alpha) + spring_th * tan(alpha)/2;

module sbar()
{
    difference()
    {
        union()
        {
            // Brackets
            for (i = [-1:2:1])
                for (j = [0:1])
                {
                    rotate([0,0,180*j])
                    translate([-d/2,0,0])
                    rotate([0,0,i*alpha])
                    translate([0,-spring_th/2,0])
                    cube([l,spring_th,h]);
                }
                
            // Handles
            for (k = [0:1])
            {
                rotate([0,0,180*k])
                {
                    difference() {
                        union() {
                    translate([d/2,0,0])
                    cylinder(thickness/2,base_d,base_d,$fs=0.3);
                    translate([d/2,-handle_width/2,0])
                    cube([handle_length,handle_width,thickness/2]);
                        }
                    // Slot
                    translate([
                        thrd_h - (handle_length + d/2) - thrd_len,
                        -handle_width/2,thickness/2-1.5])
                    cube([thrd_len, handle_width, 1.5]);
                    }
                }
            }
            
            // Male connectors
            translate([d/2 + conn_a, 0, 0])
            cylinder(thickness, male_radius, male_radius,$fs=0.3);
            translate([d/2 + handle_length - conn_b, 0, 0])
            cylinder(thickness, male_radius, male_radius,$fs=0.3);
        }
        
        // Female connectors
        translate([-d/2 - conn_a, 0, -0.1])
        cylinder(thickness, female_radius, female_radius, $fs=0.3);
        translate([-d/2 - handle_length + conn_b, 0, -0.1])
        cylinder(thickness, female_radius, female_radius, $fs=0.3);
    }
}

// Add caption
module sbar_capt()
{
    difference()
    {
        sbar();
        translate([handle_length + d/2 - 1, -0.5, 0.4])
        rotate([180, 0, 180])
        scale([0.2, 0.2, 1])
        linear_extrude(height = 0.5)
        text(str(spring_th, "/", d));
    }
}

translate([0, -handle_width/2 - 2, 0])
sbar_capt();
translate([0, handle_width/2 + 2, 0])
sbar_capt();

// Connectors between halves
translate([total_length/2-3, -handle_width/2, 0])
cube([1,handle_width,0.5]);
translate([total_length/2-6, -handle_width/2, 0])
cube([1,handle_width,0.5]);
