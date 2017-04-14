# Structure of white dwarf stars
White dwarf stars are incredibly dense objects composed of heavy nuclei and their electrons. The most stable nucleus,  56 Fe usually dominates. The density and composition of the white dwarf depends on the size of the original star and when it collapsed - the larger the star, the greater the central density of the resulting white dwarf. If the star collapses into a white dwarf prematurely, before the fusion process has had a chance to run its course, there may also be some  12C present.

Finding the mass and radius of a white dwarf involves solving two coupled differential equations simultaneously. Computationally, two of the most important methods used to solve ordinary differential equations are those of Euler (low accuracy but simple algorithm) and Runge-Kutta (high accuracy, more complicated algorithm).

These are both implemented in this project.

## Build and run
Make sure you have cpgplot installed, the C-callable version of pgplot <http://www.astro.caltech.edu/~tjp/pgplot/>

### Command line
Build:
```
make proj_2
```

Then run:
```
./proj_2
```

### Within CLion, JetBrains' IDE for C projects
The `CMakeLists.txt` file should give you everything you need to specify build targets and required libraries, allowing you to build and run within the IDE.