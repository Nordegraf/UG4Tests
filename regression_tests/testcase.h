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

#include <string>

#include "ug.h"
#include "ugbase.h"
#include "lib_disc/domain.h"
#include "lib_grid/refinement/global_multi_grid_refiner.h"

namespace ug
{
    namespace test
    {
        /**
         * \brief Base class for all testcases for regression testing
         * 
         * \tparam dim Dimension of the problem
         */
        template <int dim>
        class Testcase
        {

        protected:
            typedef ug::CPUAlgebra TAlgebra;
            typedef TAlgebra::vector_type vector_type;
            typedef TAlgebra::matrix_type matrix_type;
            typedef Domain<dim> TDomain;
            typedef ApproximationSpace<TDomain> TApproxSpace;
            typedef DirichletBoundary<TDomain, TAlgebra> TDirichletBoundary;
            typedef IDomainConstraint<TDomain, TAlgebra> TDirichletBoundaryBase;
            typedef DomainDiscretization<TDomain, TAlgebra> TDomainDiscretization;
            typedef ug::GridFunction<TDomain, TAlgebra> TGridFunction;
            typedef ug::IElemDisc<TDomain> TElemDisc;
            typedef std::string string;

        public:
            /**
             * Constructor
             *
             * \param[in]    gridname    Name of the grid file
             * \param[in]    reference   Name of the reference file
             */
            Testcase(string grid, string reference)
            {
                m_gridname = grid;
                m_reference = reference;
            }

            void run()
            {
                UG_THROW("run() function of testcase not implemented.");
            }

            /**
             * Compares the solution with the reference solution
             * 
             * \return true if the solution is equal to the reference solution
             */
            bool compare()
            {
                read_reference();

                for (size_t i = 0; i < m_spSolution->size(); i++)
                {
                    if (!isEqual((*m_spSolution)[i], (*m_spReference)[i]))
                    {
                        std::cout << "Not equal at " << i << std::endl;
                        return false;
                    }
                }

                return true;
            }

        protected:
            /**
             * \brief Refines the grid
             * 
             * \param[in] numRefs 
             */
            void refine(int numRefs)
            {
                GlobalMultiGridRefiner ref(*m_spDomain->grid(), m_spDomain->refinement_projector());

                for (int i = 0; i < numRefs; i++)
                    ref.refine();
            }

            /**
             * writes a vector containing the reference solution to a file
             * 
             * \param[in] vec  vector to write
             */
            void write_reference(std::vector<double> &vec)
            {
                std::ofstream output_file(m_reference);

                std::ostream_iterator<double> output_iterator(output_file, "\n");
                std::copy(std::begin(vec), std::end(vec), output_iterator);
            }

            /**
             * reads the reference solution file and stores it in m_spReference
             * the path to the reference solution is saved in m_reference when the testcase is created
             */
            void read_reference()
            {
                std::ifstream is(m_reference);
                std::istream_iterator<double> start(is), end;
                m_spReference = make_sp(new std::vector<double>(start, end));
            }

            /**
             * checks if two doubles are equal within a tolerance of 0.000001
             * 
             * \param a first number
             * \param b second number
             * \return true if a and b are equal within the tolerance
             */
            bool isEqual(double a, double b)
            {
                return abs(a - b) < 0.000001;
            }

            SmartPtr<TDomain> m_spDomain;
            SmartPtr<TApproxSpace> m_spApproxSpace;
            SmartPtr<TElemDisc> m_spElemDisc;
            SmartPtr<TDomainDiscretization> m_spDomainDisc;
            SmartPtr<std::vector<double>> m_spReference;
            SmartPtr<std::vector<double>> m_spSolution;
            string m_gridname;
            string m_reference;
        };

    } // namespace RegressionTest
} // namespace ug
