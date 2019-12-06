include <params.scad>

module copy_mirror(vec=[1,0,0]) 
{ 
    children(); 
    mirror(vec) children(); 
} 

extrude = 10;
depth = total_length - thrd_h + 6;

//rotate([90,0,0])
copy_mirror()
{
    base_thck = thickness/2 + 6;
    
    difference() {
        // Base
        cube([base_thck, extrude, total_length - thrd_h + 11]);
        // Slot
        cube([thickness/2 + gap, extrude, depth + gap]);
    }
    // Pin
    translate([0, 0, 2+gap])
    cube([1.5-gap, extrude, thrd_len-gap*2]);

    // Back wall
    cube([thickness, extrude - handle_width - gap,
        total_length - thrd_h + 2]);
    translate([0, 0, total_length - thrd_h + 5])
    cube([thickness, extrude - handle_width - gap, 2]);
    
    // Bottom platform
    translate([0, 0, total_length - thrd_h + 11 - 4])
    cube([base_thck, extrude * 2 - handle_width, 4]);
}
