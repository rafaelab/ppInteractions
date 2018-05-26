%module(directors="1", threads="1", allprotected="1") ProtonProtonInteraction

%include "exception.i"

%{
#include "CRPropa.h"
#include "ProtonProtonInteraction.h"
%}

/* import crpropa in wrapper */
%import (module="crpropa") "crpropa.i"

/* include plugin parts to generate wrappers for */
%include "ProtonProtonInteraction.h"




