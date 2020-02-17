%module(directors="1", threads="1", allprotected="1") NucleusNucleusInteraction

%include "exception.i"

%{
#include "CRPropa.h"
#include "NucleusNucleusInteraction.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "NucleusNucleusInteraction.h"




