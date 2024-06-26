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
    BlockBiCGStabSolver

Description
    Typedefs Preconditioned Bi-Conjugate Gradient stabilised solver.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    blockICBiCGStabSolvers.C

\*---------------------------------------------------------------------------*/

#ifndef blockBiCGStabSolvers_H
#define blockBiCGStabSolvers_H

#include "BlockBiCGStabSolver.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef BlockBiCGStabSolver<scalar> blockBiCGStabSolverScalar;
typedef BlockBiCGStabSolver<vector> blockBiCGStabSolverVector;
typedef BlockBiCGStabSolver<tensor> blockBiCGStabSolverTensor;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
