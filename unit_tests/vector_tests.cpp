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

#include <gtest/gtest.h>
#include <cmath>

#include "common/common.h"
#include "common/math/ugmath.h"

namespace ug
{
    namespace test
    {

        template <typename ValueType>
        class VectorTests : public ::testing::Test
        {
        protected:
            static const int dim = 3;
            typedef MathVector<dim, ValueType> VectorType;

            VectorTests()
            {
                reset();
            };

            VectorType a, b, c, d, e;

            void reset()
            {
                for (int i = 0; i < dim; ++i)
                {
                    a[i] = urand(0, 10);
                    b[i] = urand(0, 10);
                }

                c = 0;
                d = 0;
                e = 0;
            }
        };

        using TestTypes = ::testing::Types<float, double>;
        TYPED_TEST_SUITE(VectorTests, TestTypes);

        TYPED_TEST(VectorTests, VecAppend)
        {
            using VectorType = typename TestFixture::VectorType;

            this->reset();
            VectorType &a = this->a;
            VectorType &b = this->b;
            VectorType &c = this->c;
            VectorType &d = this->d;
            VectorType &e = this->e;

            VectorType acopy;
            VecCopy(acopy, a, 0);

            VecAppend(a, b);
            for (int i = 0; i < this->dim; ++i)
                EXPECT_EQ(a[i], acopy[i] + b[i]);

            VecAppend(c, a, b);
            VecAppend(d, a, b, c);
            VecAppend(e, a, b, c, d);

            for (int i = 0; i < this->dim; ++i)
            {
                EXPECT_EQ(c[i], a[i] + b[i]);
                EXPECT_EQ(d[i], a[i] + b[i] + c[i]);
                EXPECT_EQ(e[i], a[i] + b[i] + c[i] + d[i]);
            }
        }

        TYPED_TEST(VectorTests, VecScaleAppend)
        {
            using VectorType = typename TestFixture::VectorType;

            this->reset();
            VectorType &a = this->a;
            VectorType &b = this->b;
            VectorType &c = this->c;
            VectorType &d = this->d;
            VectorType &e = this->e;

            VecScaleAppend(c, 2, b);
            for (int i = 0; i < this->dim; ++i)
                EXPECT_EQ(c[i], 2 * b[i]);
            
            c = 0;

            VecScaleAppend(c, 2, a, 3, b);
            VecScaleAppend(d, 2, a, 3, b, 4, c);
            VecScaleAppend(e, 2, a, 3, b, 4, c, 5, d);

            for (int i = 0; i < this->dim; ++i)
            {
                EXPECT_EQ(c[i], 2 * a[i] + 3 * b[i]);
                EXPECT_EQ(d[i], 2 * a[i] + 3 * b[i] + 4 * c[i]);
                EXPECT_EQ(e[i], 2 * a[i] + 3 * b[i] + 4 * c[i] + 5 * d[i]);
            }
        }

        TYPED_TEST(VectorTests, VecAdd)
        {
            using VectorType = typename TestFixture::VectorType;

            this->reset();
            VectorType &a = this->a;
            VectorType &b = this->b;
            VectorType &c = this->c;
            VectorType &d = this->d;
            VectorType &e = this->e;

            VecAdd(c, a, b);
            VecAdd(d, a, b, c);
            VecAdd(e, a, b, c, d);

            for (int i = 0; i < this->dim; ++i)
            {
                EXPECT_EQ(c[i], a[i] + b[i]);
                EXPECT_EQ(d[i], a[i] + b[i] + c[i]);
                EXPECT_EQ(e[i], a[i] + b[i] + c[i] + d[i]);
            }
        }

        TYPED_TEST(VectorTests, VecSub)
        {
            using VectorType = typename TestFixture::VectorType;

            this->reset();
            VectorType &a = this->a;
            VectorType &b = this->b;
            VectorType &c = this->c;

            VecSubtract(c, a, b);

            for (int i = 0; i < this->dim; ++i)
                EXPECT_EQ(c[i], a[i] - b[i]);
        }

        TYPED_TEST(VectorTests, VecPow)
        {
            using VectorType = typename TestFixture::VectorType;

            this->reset();
            VectorType &a = this->a;
            VectorType &c = this->c;

            VecPow(c, a, 2);

            for (int i = 0; i < this->dim; ++i)
                EXPECT_EQ(c[i], std::pow(a[i], 2));
        }


    } // namespace test
} // namespace ug
