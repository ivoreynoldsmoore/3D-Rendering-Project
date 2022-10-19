Controls:
w - forwards relative to the camera
a - left
s - backwards
d - right
spacebar - upwards
control - downwards
q - rotate about y axis
e - rotate about y axis other direction
up arrow - rotate about x axis
down arrow - rotate about x axis other direction
1 - wireframe rendering
2 - rasterise rendering
3 - raytracing rendering
p - save image
o - enable orbit
l - disable orbit
click - debug information for raytracing at point clicked

To switch model, comment out the current "importOBJ" line, and uncomment the line which loads in the desired model.

The maximum number of reflections is bound by MAXDEPTH defined at the top of the file.

The number of lights for soft shading is set in main, but will always round down to the nearest square number.

The default renderer type is set in main.

The models have specifications that are set in main.

I changed "abs()" to "fabs()", but I haven't tested if my program works on Lab machines, sorry.

