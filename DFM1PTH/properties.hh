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
 * \brief The properties for the bulk domain in the single-phase facet coupling test.
 */
#ifndef DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULK_PROPERTIES_HH
#define DUMUX_TEST_TPFAFACETCOUPLING_ONEP_BULK_PROPERTIES_HH

#include <dune/alugrid/grid.hh>
#include <dune/foamgrid/foamgrid.hh>

//#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/material/components/co2.hh>
#include <dumux/material/components/h2o.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>

#include "problem_bulk.hh"
#include "problem_lowdim.hh"
#include "co2tables.hh"
#include "spatialparams_lowdim.hh"
#include "spatialparams_bulk.hh"


// default for the bulk grid type
//#ifndef BULKGRIDTYPE
//#define BULKGRIDTYPE Dune::ALUGrid<3, 3, Dune::simplex, Dune::nonconforming>
//#endif

//#ifndef LOWDIMGRIDTYPE
//#define LOWDIMGRIDTYPE Dune::FoamGrid<2, 3>
//#endif

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct OnePBulk { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePBulkTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePBulk>; };
struct OnePBulkMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePBulk>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
//struct Grid<TypeTag, TTag::OnePBulk> { using type = BULKGRIDTYPE; };
struct Grid<TypeTag, TTag::OnePBulk> { using type = Dune::ALUGrid<2, 2, Dune::simplex, Dune::conforming>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePBulk> { using type = OnePBulkProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePBulk>
{
    using type = BulkSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid systemHeterogeneousCO2Tables::CO2Tables
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePBulk>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
//    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
//    using type = FluidSystems::OnePLiquid< Scalar, Components::CO2<Scalar, HeterogeneousCO2Tables::CO2Tables> >;
    using type = FluidSystems::OnePLiquid< Scalar, Components::H2O<Scalar> >;
};

// create the type tag nodes
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OnePNI>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };

// we need an additional type tag for the test using mpfa in the bulk domain
struct OnePLowDimMpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
//struct Grid<TypeTag, TTag::OnePLowDim> { using type = LOWDIMGRIDTYPE; };
struct Grid<TypeTag, TTag::OnePLowDim> { using type = Dune::FoamGrid<1, 2>; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLowDim> { using type = OnePLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLowDim>
{
    using type = LowDimSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLowDim>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
//    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
//    using type = FluidSystems::OnePLiquid< Scalar, Components::CO2<Scalar, HeterogeneousCO2Tables::CO2Tables> >;
    using type = FluidSystems::OnePLiquid< Scalar, Components::H2O<Scalar> >;
};

// obtain/define some types to be used below in the property definitions and in main
template< class BulkTypeTag, class LowDimTypeTag >
class TestTraits
{
    using BulkFVGridGeometry = GetPropType<BulkTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<BulkTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set cm property in the sub-problems
using TpfaTraits = TestTraits<TTag::OnePBulkTpfa, TTag::OnePLowDimTpfa>;
using MpfaTraits = TestTraits<TTag::OnePBulkMpfa, TTag::OnePLowDimMpfa>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePBulkMpfa> { using type = typename MpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimMpfa> { using type = typename MpfaTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
