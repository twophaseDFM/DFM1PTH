// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup FacetTests
 * \brief The spatial parameters for the single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULK_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULK_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/io/grid/griddata.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the single-phase facet coupling test.
 */
template< class GridGeometry, class Scalar >
class BulkSpatialParams
: public FVSpatialParamsOneP< GridGeometry, Scalar, BulkSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = BulkSpatialParams< GridGeometry, Scalar >;
    using ParentType = FVSpatialParamsOneP< GridGeometry, Scalar, ThisType >;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    //! Export the type used for permeabilities
    using PermeabilityType = Scalar;

    BulkSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
    				  std::shared_ptr<const Dumux::GridData<Grid>> gridData,
					  const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    , gridDataPtr_(gridData)
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
//    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
//    { return permeability_; }

//    template< class ElementSolution, class FluidState >
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol)const
//								  const FluidState& fs) const
    {
        	return permeability_;
    }

    //! Returns the porosity
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
    	return porosity_;
    }

    //! returns the domain marker for an element
    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

private:
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
    PermeabilityType permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
