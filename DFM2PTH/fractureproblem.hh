// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
  * \ingroup MultiDomain
  * \ingroup MultiDomainFacet
  * \ingroup TwoPTests
  * \brief The sub-problem for the fracture domain in the exercise on two-phase flow in fractured porous media.
  */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_PROBLEM_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_PROBLEM_HH

// include the base problem and properties we inherit from
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/io/grid/griddata.hh>
#include <dune/common/indices.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
  * \brief The sub-problem for the fracture domain in the exercise on two-phase flow in fractured porous media.
 */
template<class TypeTag>
class FractureSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;


    // some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx,
    };

public:
    //! The constructor
    FractureSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                       std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
					   std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                       const std::string& paramGroup)
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , gridDataPtr_(gridData)
    , aperture1_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture1"))
    , aperture2_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture2"))
    , aperture3_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture3"))
    , aperture4_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture4"))
    , aperture5_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture5"))
    , IsPureCO2_(getParamFromGroup<Scalar>(paramGroup, "Problem.IsPureCO2"))
    {
        // initialize the fluid system, i.e. the tabulation
        // of water properties. Use the default p/T ranges.
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        FluidSystem::init();
    }

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        // We only use no-flow boundary conditions for all immersed fractures
        // in the domain (fracture tips that do not touch the domain boundary)
        // Otherwise, we would lose mass leaving across the fracture tips.
        values.setAllNeumann();

        // However, there is one fracture reaching the top boundary. For this
        // fracture tip we set Dirichlet Bcs as in the matrix domain
        // TODO dumux-course-task A
        // Change boundary conditions
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6)
            values.setAllDirichlet();

        return values;
    }

    //! Evaluate the source term in a sub-control volume of an element
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain using the function in the coupling manager
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);

        // these sources are in kg/s, divide by volume and extrusion to have it in kg/s/m³
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Set the aperture as extrusion factor.
//    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
//    {
//        // We treat the fractures as lower-dimensional in the grid,
//        // but we have to give it the aperture as extrusion factor
//        // such that the dimensions are correct in the end.
//    	return aperture3_;
//    }
    template< class ElementSolution >
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scvf,
                           const ElementSolution& elemSol) const
    {
//        static constexpr auto fractureDomainId = Dune::index_constant<1>();
//        const auto& fractureProblem = couplingManager().problem(fractureDomainId);
//
//        // use helper function in coupling manager to obtain the element this face couples to
//        const auto fractureElement = couplingManager().getLowDimElement(element, scvf);
//
//        // obtain marker from the spatial params of the fracture problem
//        const auto fractureElementMarker = fractureProblem.spatialParams().getElementDomainMarker(fractureElement);
//
		if (getElementDomainMarker(element) == 1)
			return aperture1_;
		else if (getElementDomainMarker(element) == 2)
			return aperture2_;
		else if (getElementDomainMarker(element) == 3)
			return aperture3_;
		else if (getElementDomainMarker(element) == 4)
			return aperture4_;
		else
        	return aperture5_;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // For the grid used here, the height of the domain is equal
        // to the maximum y-coordinate
        const auto domainHeight = this->gridGeometry().bBoxMax()[1] + 6000;

        // we assume a constant water density of 1000 for initial conditions!
        const auto& g = this->spatialParams().gravity(globalPos);
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 - (domainHeight - globalPos[1])*densityW*g[1];
        values[temperatureIdx] = 283.0 + (domainHeight - globalPos[1])*0.03;
        if (IsPureCO2_)
        	{values[saturationIdx] = 1.0 - eps_;}
        else
        	{values[saturationIdx] = 0.0;}

        return values;
    }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
//    Scalar temperature() const
//    { return 383.15; /*10°*/ }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
    Scalar aperture1_,aperture2_,aperture3_,aperture4_,aperture5_;
    bool IsPureCO2_;
    static constexpr Scalar eps_ = 1e-7;
};

} // end namespace Dumux

#endif
