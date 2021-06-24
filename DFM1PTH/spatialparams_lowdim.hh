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
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIM_SPATIALPARAMS_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIM_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/io/grid/griddata.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the single-phase facet coupling test.
 */
template< class GridGeometry, class Scalar >
class LowDimSpatialParams
: public FVSpatialParamsOneP< GridGeometry, Scalar, LowDimSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = LowDimSpatialParams< GridGeometry, Scalar >;
    using ParentType = FVSpatialParamsOneP< GridGeometry, Scalar, ThisType >;
    using SubControlVolume = typename GridGeometry::SubControlVolume;

    using GridView = typename GridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    //! Export the type used for permeabilities
    using PermeabilityType = Scalar;

    LowDimSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
    				  std::shared_ptr<const Dumux::GridData<Grid>> gridData,
					  const std::string& paramGroup = "")
    : ParentType(gridGeometry)
    , gridDataPtr_(gridData)
    {
        porosity_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Porosity");
        a1_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture1"));
        a2_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture2"));
        a3_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture3"));
        a4_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture4"));
        a5_ = (getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture5"));
        kn_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.Stiffness"));
        wte_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.WettingPhaseThermalExpansionCoefficient"));
        a_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.a"));
        b_ = (getParamFromGroup<Scalar>(paramGroup, "Problem.b"));
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
//    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
//    { return permeability_; }

//    template< class ElementSolution, class FluidState>
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol)const
//								  const FluidState& fs) const
    {
    	int pressureIdx = 0;
    	int temperatureIdx = 1;
        const auto& priVars = elemSol[scv.localDofIndex()];

		const GlobalPosition& globalPos = scv.center();
		const auto domainHeight = 100.0 + 6000;
		const auto initialTemperature = 283.0 + (domainHeight - globalPos[dimWorld-1])*0.03;
		const auto ThermalExpan = (priVars[temperatureIdx] - initialTemperature) * wte_;
		Scalar densityW = 1000.0;
		const auto g = -9.81;
		const auto initialPressure = 1e5 - (domainHeight - globalPos[dimWorld-1])*densityW*g;
		const auto deltaP = priVars[pressureIdx] - initialPressure;

        Scalar a1 = a1_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a1_;
        Scalar a2 = a2_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a2_;
        Scalar a3 = a3_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a3_;
        Scalar a4 = a4_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a4_;
        Scalar a5 = a5_ + a_ * deltaP/kn_ + b_ * ThermalExpan * a5_;

        if (getElementDomainMarker(element) == 1)
        	return a1*a1/12;
        else if(getElementDomainMarker(element) == 2)
        	return a2*a2/12;
        else if(getElementDomainMarker(element) == 3)
        	return a3*a3/12;
        else if(getElementDomainMarker(element) == 4)
        	return a4*a4/12;
        else
        	return a5*a5/12;

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
    Scalar porosity_;
    Scalar a1_,a2_,a3_,a4_,a5_;
    Scalar kn_, wte_;
    Scalar a_, b_ ;
};

} // end namespace Dumux

#endif
