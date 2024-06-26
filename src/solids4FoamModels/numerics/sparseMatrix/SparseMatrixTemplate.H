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
    SparseMatrixTemplate

Description
    A simple but convenient matrix for storing sparse matrices.
    Values (templated) are stored using a HashTable where the Hash is based on
    the (rowI, columns) key.

    Example usage:

        SparseMatrixTemplate<tensor> matT;
        matT(1, 0) = tensor(1,2,3,4,5,6,7,8,9);
        matT(0, 2) += 3*I;

        SparseMatrixTemplate<scalar> matS;
        matS(1, 0) = 1;
        matS(0, 2) += 3;

Author
    Philip Cardiff, UCD.

SourceFiles
    SparseMatrixTemplate.C

\*---------------------------------------------------------------------------*/

#ifndef SparseMatrixTemplate_H
#define SparseMatrixTemplate_H

#include "HashTable.H"
#include "FixedList.H"
#include "tensor.H"
#include "refCount.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------* \
                       Class SparseMatrixTemplate Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SparseMatrixTemplate
:
    public refCount
{
    typedef HashTable
    <
        Type, FixedList<label, 2>, FixedList<label, 2>::Hash<>
    > SparseMatrixTemplateData;

    // Private data

        //- Hashtable (map) where the values are stored
        SparseMatrixTemplateData data_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const SparseMatrixTemplate&);

public:

    //- Runtime type information
    TypeName("SparseMatrixTemplate");


    // Constructors

        //- Construct given approximate size
        SparseMatrixTemplate(const label size = 128);

        //- Copy constructer
        SparseMatrixTemplate(const SparseMatrixTemplate& mat);


    // Destructor

        virtual ~SparseMatrixTemplate()
        {}


    // Member Functions

        // Access

            //- Const access to the data
            const SparseMatrixTemplateData& data() const
            {
                return data_;
            }

            //- Non-const access to the data
            SparseMatrixTemplateData& data()
            {
                return data_;
            }

            //- Count and return the number of block rows in the matrix
            label nBlockRows() const;

            //- Print out the matrix coefficients
            void print() const;

        // Modifiers

            //- Clear the matrix but do not change its capacity
            void clear()
            {
                data_.clear();
            }

        // Operators

            //- Non-const access to a value entry. If the entry does not exist
            //  then create it and initialise it to zero
            Type& operator()(const label rowI, const label colI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    typedef SparseMatrixTemplate<scalar> sparseScalarMatrix;
    typedef SparseMatrixTemplate<vector> sparseVectorMatrix;
    typedef SparseMatrixTemplate<tensor> sparseTensorMatrix ;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SparseMatrixTemplate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
