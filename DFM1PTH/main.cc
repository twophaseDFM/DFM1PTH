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
 * \brief Test for the one-phase facet coupling model.
 */

#include <config.h>

#include <iostream>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/dumuxmessage.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>

#include <dumux/multidomain/facet/gridmanager.hh>

#include <dumux/io/vtkoutputmodule.hh>

#include "computevelocities.hh"

#include "properties.hh"


// main program
int main(int argc, char** argv)
{
    using namespace Dumux;

    //////////////////////////////////////////////////////
    //////////////////////////////////////////////////////

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // initialize parameter tree
    Parameters::init(argc, argv);

    //////////////////////////////////////////////////////
    // try to create the grids (from the given grid file)
    //////////////////////////////////////////////////////
    using BulkProblemTypeTag = Properties::TTag::BULKTYPETAG;
    using LowDimProblemTypeTag = Properties::TTag::LOWDIMTYPETAG;
    using BulkGrid = GetPropType<BulkProblemTypeTag, Properties::Grid>;
    using LowDimGrid = GetPropType<LowDimProblemTypeTag, Properties::Grid>;
    using BulkFVGridGeometry = GetPropType<BulkProblemTypeTag, Properties::GridGeometry>;
	using LowDimFVGridGeometry = GetPropType<LowDimProblemTypeTag, Properties::GridGeometry>;
	using TheMultiDomainTraits = Dumux::MultiDomainTraits<BulkProblemTypeTag, LowDimProblemTypeTag>;
	using TheCouplingMapper = Dumux::FacetCouplingMapper<BulkFVGridGeometry, LowDimFVGridGeometry>;
    using TheCouplingManager = Dumux::FacetCouplingManager<TheMultiDomainTraits, TheCouplingMapper>;


    using GridManager = FacetCouplingGridManager<BulkGrid, LowDimGrid>;
    GridManager gridManager;
    gridManager.init();
    gridManager.loadBalance();

    ////////////////////////////////////////////////////////////
    // run stationary, non-linear problem on this grid
    ////////////////////////////////////////////////////////////
    static constexpr std::size_t matrixGridId = 0;
    static constexpr std::size_t fractureGridId = 1;

    // we compute on the leaf grid views
    const auto& bulkGridView = gridManager.template grid<0>().leafGridView();
    const auto& lowDimGridView = gridManager.template grid<1>().leafGridView();

    // create the finite volume grid geometries
    auto bulkFvGridGeometry = std::make_shared<BulkFVGridGeometry>(bulkGridView);
    auto lowDimFvGridGeometry = std::make_shared<LowDimFVGridGeometry>(lowDimGridView);
    bulkFvGridGeometry->update();
    lowDimFvGridGeometry->update();

    // the coupling mapper
    using TestTraits = Properties::TestTraits<BulkProblemTypeTag, LowDimProblemTypeTag>;
    auto couplingMapper = std::make_shared<typename TestTraits::CouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, gridManager.getEmbeddings());

    // the coupling manager
    using CouplingManager = typename TestTraits::CouplingManager;
    auto couplingManager = std::make_shared<CouplingManager>();

    // the problems (boundary conditions)
    using BulkProblem = GetPropType<BulkProblemTypeTag, Properties::Problem>;
    using LowDimProblem = GetPropType<LowDimProblemTypeTag, Properties::Problem>;
    auto bulkGridData = gridManager.getGridData()->template getSubDomainGridData<matrixGridId>();
    auto bulkSpatialParams = std::make_shared<typename BulkProblem::SpatialParams>(bulkFvGridGeometry, bulkGridData, "Bulk");
    auto bulkProblem = std::make_shared<BulkProblem>(bulkFvGridGeometry, bulkSpatialParams, couplingManager, "Bulk");
    auto LowDimGridData = gridManager.getGridData()->template getSubDomainGridData<fractureGridId>();
    auto lowDimSpatialParams = std::make_shared<typename LowDimProblem::SpatialParams>(lowDimFvGridGeometry, LowDimGridData, "LowDim");
    auto lowDimProblem = std::make_shared<LowDimProblem>(lowDimFvGridGeometry, lowDimSpatialParams, couplingManager, LowDimGridData, "LowDim");

    // the solution vector
//    using MDTraits = typename TestTraits::MDTraits;
    using SolutionVector = typename TheMultiDomainTraits::SolutionVector;
    SolutionVector x, xOld;

    static const auto bulkId = typename TheMultiDomainTraits::template SubDomain<0>::Index();
    static const auto lowDimId = typename TheMultiDomainTraits::template SubDomain<1>::Index();
    x[bulkId].resize(bulkFvGridGeometry->numDofs());
    x[lowDimId].resize(lowDimFvGridGeometry->numDofs());
    bulkProblem->applyInitialSolution(x[bulkId]);
    lowDimProblem->applyInitialSolution(x[lowDimId]);

    // instantiate the class holding the coupling maps between the domains
    // this needs the information on embeddings (connectivity between matrix
    // and fracture domain). This information is extracted directly from the
    // grid during file read and can therefore be obtained from the grid manager.
    const auto embeddings = gridManager.getEmbeddings();
//    auto couplingMapper = std::make_shared<TheCouplingMapper>();
    couplingMapper->update(*bulkFvGridGeometry, *lowDimFvGridGeometry, embeddings);

    // the grid variables
    using BulkGridVariables = GetPropType<BulkProblemTypeTag, Properties::GridVariables>;
    using LowDimGridVariables = GetPropType<LowDimProblemTypeTag, Properties::GridVariables>;
    auto bulkGridVariables = std::make_shared<BulkGridVariables>(bulkProblem, bulkFvGridGeometry);
    auto lowDimGridVariables = std::make_shared<LowDimGridVariables>(lowDimProblem, lowDimFvGridGeometry);
    bulkGridVariables->init(x[bulkId]);
    lowDimGridVariables->init(x[lowDimId]);
    xOld = x;

    // initialize coupling manager
//    auto couplingManager = std::make_shared<TheCouplingManager>();
    couplingManager->init(bulkProblem, lowDimProblem, couplingMapper, x);

    // intialize the vtk output module
    using BulkSolutionVector = std::decay_t<decltype(x[bulkId])>;
    using LowDimSolutionVector = std::decay_t<decltype(x[lowDimId])>;
    VtkOutputModule<BulkGridVariables, BulkSolutionVector> bulkVtkWriter(*bulkGridVariables, x[bulkId], bulkProblem->name(), "Bulk");
    VtkOutputModule<LowDimGridVariables, LowDimSolutionVector> lowDimVtkWriter(*lowDimGridVariables, x[lowDimId], lowDimProblem->name(), "LowDim");

    // Add model specific output fields
    using BulkIOFields = GetPropType<BulkProblemTypeTag, Properties::IOFields>;
    using LowIOFields = GetPropType<LowDimProblemTypeTag, Properties::IOFields>;
//    using lowDimVelocityOutput = GetPropType<LowDimProblemTypeTag, Properties::VelocityOutput>;
//    lowDimVtkWriter.addVelocityOutput(std::make_shared<lowDimVelocityOutput>(*lowDimGridVariables));
//    using bulkVelocityOutput = GetPropType<BulkProblemTypeTag, Properties::VelocityOutput>;
//    bulkVtkWriter.addVelocityOutput(std::make_shared<bulkVelocityOutput>(*bulkGridVariables));
    BulkIOFields::initOutputModule(bulkVtkWriter);
    LowIOFields::initOutputModule(lowDimVtkWriter);

    std::vector<int> bulkDomainMarkers(bulkFvGridGeometry->gridView().size(0));
    for (const auto& element : elements(bulkFvGridGeometry->gridView()))
        bulkDomainMarkers[bulkFvGridGeometry->elementMapper().index(element)] = bulkProblem->spatialParams().getElementDomainMarker(element);
    bulkVtkWriter.addField(bulkDomainMarkers, "domainMarker");

    std::vector<int> lowDimDomainMarkers(lowDimFvGridGeometry->gridView().size(0));
    for (const auto& element : elements(lowDimFvGridGeometry->gridView()))
        lowDimDomainMarkers[lowDimFvGridGeometry->elementMapper().index(element)] = lowDimProblem->spatialParams().getElementDomainMarker(element);
    lowDimVtkWriter.addField(lowDimDomainMarkers, "domainMarker");

    // write initial solution
    bulkVtkWriter.write(0.0);
    lowDimVtkWriter.write(0.0);

    using Velocity = Dune::FieldVector<double, 2>;
    using VelocityVector = std::vector<Velocity>;

	VelocityVector velocityBulk(bulkFvGridGeometry->gridView().size(0));
	VelocityVector velocityLowDim(lowDimFvGridGeometry->gridView().size(0));

	bulkVtkWriter.addField(velocityBulk, "velocity");
	lowDimVtkWriter.addField(velocityLowDim, "velocity");


    // get some time loop parameters
    const auto tEnd = getParam<double>("TimeLoop.TEnd");
    const auto maxDt = getParam<double>("TimeLoop.MaxTimeStepSize");
    auto dt = getParam<double>("TimeLoop.DtInitial");

    // instantiate time loop
    auto timeLoop = std::make_shared< TimeLoop<double> >(/*startTime*/0.0, dt, tEnd);
    timeLoop->setMaxTimeStepSize(maxDt);

    // the assembler
    using Assembler = MultiDomainFVAssembler<TheMultiDomainTraits, CouplingManager, DiffMethod::numeric, /*implicit?*/true>;
    auto assembler = std::make_shared<Assembler>( std::make_tuple(bulkProblem, lowDimProblem),
                                                  std::make_tuple(bulkFvGridGeometry, lowDimFvGridGeometry),
                                                  std::make_tuple(bulkGridVariables, lowDimGridVariables),
                                                  couplingManager,
												  timeLoop,xOld);

    // the linear solver
    using LinearSolver = ILU0BiCGSTABBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto newtonSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    // time loop
    timeLoop->start(); do
    {
        // set previous solution for storage evaluations
        assembler->setPreviousSolution(xOld);

        // solve the non-linear system with time step control
        newtonSolver->solve(x, *timeLoop);

        computeVelocities<BulkProblemTypeTag>(bulkId, *assembler, *couplingManager, *bulkFvGridGeometry, *bulkGridVariables, x[bulkId], velocityBulk);
		computeVelocities<LowDimProblemTypeTag>(lowDimId, *assembler, *couplingManager, *lowDimFvGridGeometry, *lowDimGridVariables, x[lowDimId], velocityLowDim);

        // make the new solution the old solution
        xOld = x;
        bulkGridVariables->advanceTimeStep();
        lowDimGridVariables->advanceTimeStep();

        // advance to the time loop to the next step
        timeLoop->advanceTimeStep();

        // write vtk output
        bulkVtkWriter.write(timeLoop->time());
        lowDimVtkWriter.write(timeLoop->time());

        // report statistics of this time step
        timeLoop->reportTimeStep();

        // set new dt as suggested by the Newton solver
        timeLoop->setTimeStepSize(newtonSolver->suggestTimeStepSize(timeLoop->timeStepSize()));

    } while (!timeLoop->finished());

    // output some Newton statistics
    newtonSolver->report();

    // report time loop statistics
    timeLoop->finalize();

//    // linearize & solve
//    newtonSolver->solve(x);

//    // update grid variables for output
//    bulkGridVariables->update(x[bulkId]);
//    lowDimGridVariables->update(x[lowDimId]);

//    // write vtk output
//    bulkVtkWriter.write(1.0);
//    lowDimVtkWriter.write(1.0);

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/false);

    return 0;
}
