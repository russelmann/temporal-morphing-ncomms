% Holder used for shells with a matching interface

length = 70;
width = 10;
thickness = 10;
slot_d = 4.5;
bite_d = 2.3*2;
bite_l = 20;

eps = 1e-3;

difference() {
    translate([0,-width/2, -thickness/2])
    cube([length, width, thickness]);
    translate([bite_l + 10, -slot_d/2, -thickness/2-1])
    cube([length - bite_l - 20, slot_d, thickness+2]);
    translate([-eps, -width/2 -eps, -bite_d/2])
    cube([bite_l, width + eps*2, bite_d]);
}
for (s = [-1:2:1])
{
    translate([2, -width/2, -1+s*slot_d/2])
    cube([1.8, width, 2]);
}
