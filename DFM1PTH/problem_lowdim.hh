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
 * \brief The problem for the lower-dimensional domain in the single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIMPROBLEM_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_LOWDIMPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/io/grid/griddata.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The lower-dimensional test problem for the incompressible
 *        one-phase model with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePLowDimProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Grid = typename GridView::Grid;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    static constexpr int dimWorld = GridView::dimensionworld;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum {
        // index of the transport equation
        conti0EqIdx = Indices::conti0EqIdx,
        energyEqIdx = Indices::energyEqIdx
    };

public:
    OnePLowDimProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                      std::shared_ptr<CouplingManager> couplingManager,
					  std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                      const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , gridDataPtr_(gridData)
    , couplingManagerPtr_(couplingManager)
    , aperture_(getParam<Scalar>("Problem.FractureAperture"))
    , aperture1_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture1"))
    , aperture2_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture2"))
    , aperture3_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture3"))
    , aperture4_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture4"))
    , aperture5_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture5"))
    , kn_(getParamFromGroup<Scalar>(paramGroup, "Problem.Stiffness"))
    , wte_(getParamFromGroup<Scalar>(paramGroup, "Problem.WettingPhaseThermalExpansionCoefficient"))
    , nte_(getParamFromGroup<Scalar>(paramGroup, "Problem.NonWettingPhaseThermalExpansionCoefficient"))
    , a_(getParamFromGroup<Scalar>(paramGroup, "Problem.a"))
    , b_(getParamFromGroup<Scalar>(paramGroup, "Problem.b"))
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" +
                         getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6)
            values.setAllDirichlet();
        return values;
    }

    //! Evaluates the source term at a given position.
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! Sets the aperture as extrusion factor.
//    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
//    { return aperture_; }
    template< class ElementSolution >
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume& scv,
                           const ElementSolution& elemSol) const
    {
		const auto& priVars = elemSol[scv.localDofIndex()];

		const GlobalPosition& globalPos = scv.center();
		const auto initialValues = initialAtPos(globalPos);
		const auto ThermalExpan = (priVars[temperatureIdx] - initialValues[temperatureIdx]) * wte_;
		const auto deltaP = priVars[pressureIdx] - initialValues[pressureIdx];

		if (getElementDomainMarker(element) == 1)
			return aperture1_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture1_ ;
		else if (getElementDomainMarker(element) == 2)
			return aperture2_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture2_;
		else if (getElementDomainMarker(element) == 3)
			return aperture3_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture3_;
		else if (getElementDomainMarker(element) == 4)
			return aperture4_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture4_;
		else
        	return aperture5_ + a_ * deltaP/kn_ + b_ * ThermalExpan * aperture5_;
    }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const auto domainHeight = this->gridGeometry().bBoxMax()[dimWorld-1] + 6000;

        // we assume a constant water density of 1000 for initial conditions!
        const auto& g = this->spatialParams().gravity(globalPos);
        PrimaryVariables values(0.0);
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 - (domainHeight - globalPos[dimWorld-1])*densityW*g[dimWorld-1];
        values[temperatureIdx] = 283.0 + (domainHeight - globalPos[dimWorld-1])*0.03;

        return values;
    }

    //! Returns the temperature in \f$\mathrm{[K]}\f$ in the domain.
//    Scalar temperature() const
//    { return 283.15; /*10°*/ }

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    int getElementDomainMarker(const Element& element) const
    { return gridDataPtr_->getElementDomainMarker(element); }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
    Scalar aperture_;
    Scalar aperture1_,aperture2_,aperture3_,aperture4_,aperture5_;
    std::string problemName_;
    static constexpr Scalar eps_ = 1e-7;
    Scalar kn_, nte_, wte_;
    Scalar a_, b_;
};

} // end namespace Dumux

#endif
