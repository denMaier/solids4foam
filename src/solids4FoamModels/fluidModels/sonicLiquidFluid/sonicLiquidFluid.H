/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    sonicLiquidFluid

Description
    sonicLiquidFoam fluid model.

Author
    Based on sonicLiquidFoam.
    Ported by Philip Cardiff and Iago Oliveira.

SourceFiles
    sonicLiquidFluid.C

\*---------------------------------------------------------------------------*/

// There are significant differences between OpenFOAM/FOAM versions so to
// improve readability and maintainability, the ported versions will be given
// in separate files
#ifdef OPENFOAM_ORG
    #include "sonicLiquidFluid.foundation.H"
#elif OPENFOAM_COM
    #include "sonicLiquidFluid.esi.H"
#else
    #include "sonicLiquidFluid.foamextend.H"
#endif

// ************************************************************************* //
