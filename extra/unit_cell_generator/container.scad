include <params.scad>

num_br = 20; // number of brackets to pack
fth = 5; // frame thickness
fhg = 2; // frame height

fgap = 0.5; // gap

handle_length = (total_length - d) / 2;

difference()
{
    fr_len = (thickness + fgap * 2) * num_br;
    fr_wid = total_length + fgap * 2;
    translate([-fth, -fth - fr_wid/2, -fhg])
    cube([fr_len + 2 * fth, fr_wid + 2 * fth,
        handle_width/2 + fhg]);
    {
        translate([0, -fr_wid/2, 0])
        cube([fr_len, fr_wid, handle_width]);
        translate([0, -5, -fhg])
        cube([fr_len, 10, 10]);
        for (k = [0:1])
        {
            mirror([0,k,0])
            translate([0, 1-total_length/2, -fhg])
            cube([fr_len, 6, 10]);
        }
        translate([-fth/2-1, -10, 0])
        cube([2, 20, handle_width/2]);
    }
}
translate([num_br*(thickness+1)+fth/2-1+gap,
    -10+2*gap, handle_width/2])
cube([2-gap*2, 20-gap, handle_width/2-gap]);


for (i = [0:num_br-1])
{
    translate([i * (thickness + fgap * 2), 0, 0])
    for (k = [0:1])
    {
        mirror([0,k,0])
        translate([gap, gap, 0])
        translate([thickness/2 - 1.5 + fgap,
            thrd_h - (handle_length + d/2) - thrd_len, 0])
        cube([1.5 * 2 - gap*2, thrd_len - gap*2, 1]);
    }
}