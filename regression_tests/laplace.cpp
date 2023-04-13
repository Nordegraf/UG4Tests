/*
 * Copyright (c) 2023:  G-CSC, Goethe University Frankfurt
 * Author: Niklas Conen
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "ug.h"
#include "ugbase.h"
#include "../../ConvectionDiffusion/convection_diffusion_base.h"
#include "../../ConvectionDiffusion/fv1/convection_diffusion_fv1.h"
#include "lib_algebra/operator/linear_solver/bicgstab.h"
#include "lib_algebra/operator/preconditioner/preconditioners.h"
#include "lib_algebra/operator/linear_solver/agglomerating_solver.h"
#include "../../SuperLU/super_lu.h"

#include "testcase.h"


namespace ug
{
    namespace test
    {
        /**
         * \brief Laplace testcase
         */
        class Laplace : public Testcase<3>
        {
            typedef Testcase<3> base_type;
            typedef ug::ConvectionDiffusionPlugin::ConvectionDiffusionFV1<TDomain> TConvDiff;
            typedef ug::AssembledMultiGridCycle<TDomain, TAlgebra> GMG;

            using Testcase::Testcase;

        public:
            /**
             * Runs the Laplace testcase
             */
            void run()
            {
                AlgebraType algebra("CPU", 1);
                ug::bridge::InitUG(3, algebra);

                // Domain
                m_spDomain = make_sp(new TDomain());
                LoadDomain(*m_spDomain, m_gridname.c_str());
                refine(4);

                // Approximation Space
                m_spApproxSpace = make_sp(new TApproxSpace(m_spDomain));
                m_spApproxSpace->add("c", "Lagrange", 1);
                m_spApproxSpace->init_top_surface();

                // Element Discretization
                SmartPtr<TConvDiff> cd = make_sp(new TConvDiff("c", "Inner"));
                cd->set_diffusion(1.0);
                cd->set_reaction(0.0);
                m_spElemDisc = cd;

                // Dirichlet Boundary Conditions
                SmartPtr<TDirichletBoundary> boundary = make_sp(new TDirichletBoundary());
                boundary->add(-1, "c", "bndNegative");
                boundary->add(1, "c", "bndPositive");
                m_spDirichlet = boundary;
                
                // Domain Discretization
                m_spDomainDisc = make_sp(new TDomainDiscretization(m_spApproxSpace));
                m_spDomainDisc->add(m_spElemDisc);
                m_spDomainDisc->add(m_spDirichlet);

                // Solver Configuration
                // Jacobi smoother with default damping of 0.66
                SmartPtr<Jacobi<TAlgebra>> smoother = make_sp(new Jacobi<TAlgebra>(0.66));

                // Base Solver: SuperLU
                SmartPtr<ILinearOperatorInverse<vector_type, vector_type>> superlu = make_sp(new SuperLUSolver<TAlgebra>());
                SmartPtr<AgglomeratingSolver<TAlgebra>> baseSolver = make_sp(new AgglomeratingSolver<TAlgebra>(superlu));

                // Transfer
                SmartPtr<StdTransfer<TDomain, TAlgebra>> transfer = make_sp(new StdTransfer<TDomain, TAlgebra>());
                transfer->enable_p1_lagrange_optimization(true);

                // Geometric Multigrid Preconditioner
                SmartPtr<GMG> gmg = make_sp(new GMG(m_spApproxSpace));
                gmg->set_base_solver(baseSolver);
                gmg->set_smoother(smoother);
                gmg->set_base_level(0);
                gmg->set_cycle_type("V");
                gmg->set_num_presmooth(3);
                gmg->set_num_postsmooth(3);
                gmg->set_rap(false);
                gmg->set_smooth_on_surface_rim(false);
                gmg->set_emulate_full_refined_grid(false);
                gmg->set_gathered_base_solver_if_ambiguous(false);
                gmg->set_transfer(transfer);

                // Convergence Check
                SmartPtr<StdConvCheck<vector_type>> ConvCheck = make_sp(new StdConvCheck<vector_type>(100, 1e-12, 1e-6, true));

                // BiCGStab Solver
                m_spSolver = make_sp(new BiCGStab<TAlgebra::vector_type>());
                m_spSolver->set_preconditioner(gmg);
                m_spSolver->set_convergence_check(ConvCheck);

                // Assemble Linear Operator
                m_spOp = make_sp(new AssembledLinearOperator<TAlgebra>(m_spDomainDisc));
                m_spU = make_sp(new TGridFunction(m_spApproxSpace));
                m_spB = make_sp(new TGridFunction(m_spApproxSpace));

                m_spU->set(0.0);
                m_spDomainDisc->adjust_solution(*m_spU);
                m_spDomainDisc->assemble_linear(*m_spOp, *m_spB);

                // Solve
                SmartPtr<vector_type> u = m_spU;
                m_spSolver->init(m_spOp, *m_spU);
                m_spSolver->apply(*m_spU, *m_spB);

                // Save Solution
                SmartPtr<std::vector<double>> sol = make_sp(new std::vector<double>);
                m_spOp->get_values(*sol);
                m_spSolution = sol;

                /*SaveMatrixForConnectionViewer(*m_spU, *m_spOp, "laplace_matrix.mat");
                SaveVectorForConnectionViewer(*m_spB, "laplace_rhs.vec");
                VTKOutput<TGridFunction::dim> out;
                out.print("laplace3d.vtk", *m_spU, true);*/
            }

        protected:
            SmartPtr<TDirichletBoundaryBase> m_spDirichlet;
            SmartPtr<AssembledLinearOperator<TAlgebra>> m_spOp;
            SmartPtr<TGridFunction> m_spU;
            SmartPtr<TGridFunction> m_spB;
            SmartPtr<BiCGStab<TAlgebra::vector_type>> m_spSolver;
        };

    } // namespace RegressionTest
} // namespace ug