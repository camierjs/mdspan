// Copyright (c) 2010-2023, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifdef _WIN32
#define _USE_MATH_DEFINES
#include <cmath>
#endif

#include <list>
#include <type_traits>

#include "unit_tests.hpp"

#include "fem/mdgridfunc.hpp"
#include "general/forall.hpp"
#include "general/mdarray.hpp"
#include "linalg/mdvector.hpp"

#if 1
#include <typeinfo>
#ifndef _MSC_VER
#include <cxxabi.h>
#endif
#include <memory>
#include <string>
#include <cstdlib>

template <class T>
std::string type_name()
{
   typedef typename std::remove_reference<T>::type TR;
   std::unique_ptr<char, void(*)(void*)> own(
#ifndef _MSC_VER
      abi::__cxa_demangle(typeid(TR).name(), nullptr,
                          nullptr, nullptr),
#else
      nullptr,
#endif
      std::free
   );
   std::string r = own != nullptr ? own.get() : typeid(TR).name();
   if (std::is_const<TR>::value) { r += " const"; }
   if (std::is_volatile<TR>::value) { r += " volatile"; }
   if (std::is_lvalue_reference<T>::value) { r += "&"; }
   else if (std::is_rvalue_reference<T>::value) { r += "&&"; }
   return r;
}
#endif

using namespace mfem;

TEST_CASE("MDArray", "[MDSpan], [MDArray]")
{
   SECTION("Types")
   {
      MDArray<int,3> mda;
      REQUIRE(mda.Size() == 0);
      REQUIRE(std::is_same<decltype(mda.HostRead()), int const*>());
      REQUIRE(std::is_same<decltype(mda.MDHostRead()), mfem::MDTensor<3, int const> const>());
   }

   SECTION("SetSize")
   {
      constexpr int NA = 11, NB = 22, NC = 33;
      MDArray<int,3> abc;
      abc.SetSize(NA, NB, NC);
      REQUIRE(abc.Size() == NA*NB*NC);
      abc.HostRead();
      abc.MDHostRead();
   }

   SECTION("SetLayout")
   {
      constexpr int NA = 11, NB = 22, NC = 33;
      constexpr int na =  0, nb =  1, nc =  2;
      MDLayout<3> layout_012({na,nb,nc});

      MDArray<int,3> abc(NA, NB, NC);
      MDArray<int,3, MDLayout<3>> abc_ini(NA, NB, NC);
      MDArray<int,3, MDLayout<3>> abc_set(NA, NB, NC);
      abc_set.SetLayout(layout_012);

      REQUIRE(abc_set.Offset(na,nb,nc) == abc.Offset(na,nb,nc));
      REQUIRE(abc_set.Offset(na,nb,nc) == abc_ini.Offset(na,nb,nc));
   }
}

TEST_CASE("MDVector", "[MDSpan], [MDVector]")
{
   SECTION("Types")
   {
      MDVector<3> mdv;
      REQUIRE(mdv.Size() == 0);
      REQUIRE(std::is_same<decltype(mdv.HostRead()), double const*>());
      REQUIRE(std::is_same<decltype(mdv.MDHostRead()), mfem::MDTensor<3, double const> const>());
   }

   SECTION("SetSize")
   {
      constexpr int NA = 11, NB = 22, NC = 33;
      MDVector<3> abc;
      abc.SetSize(NA, NB, NC);
      REQUIRE(abc.Size() == NA*NB*NC);
   }

   SECTION("SetLayout")
   {
      constexpr int NA = 11, NB = 22, NC = 33;
      constexpr int na =  0, nb =  1, nc =  2;
      MDLayout<3> layout_012({na,nb,nc});

      MDVector<3> abc(NA, NB, NC);
      MDVector<3, MDLayout<3>> abc_ini(NA, NB, NC);
      MDVector<3, MDLayout<3>> abc_set(NA, NB, NC);
      abc_set.SetLayout(layout_012);

      REQUIRE(abc_set.Offset(na,nb,nc) == abc.Offset(na,nb,nc));
      REQUIRE(abc_set.Offset(na,nb,nc) == abc_ini.Offset(na,nb,nc));
   }

   SECTION("Offset")
   {
      constexpr int NA = 18, NB = 2, NC = 36, ND = 32;
      // Fortran col major: (N1:18, 2, 36, Nd:32)
      //                    (    0, 1,  2,     3)
      // = 0 + 18( 1 + 2( 2 + 36( 3))) = 3978
      MDVector<4> left(NA,NB,NC,ND); // default layout is LayoutLeft
      REQUIRE(left.Offset(0,1,2,3) == 3978);

      // C/C++ row major: (N1:18, 2, 36, Nd:32)
      //                  (    0, 1,  2,     3)
      // = 3 + 32( 2 + 36( 1 + 2(0))) = 1219
      // = (((0)*2 + 1) * 36 + 2) * 32 + 3
      MDVector<4> right(NA, NB, NC, ND);
      right.SetLayout(MDLayout<4>({3,2,1,0}));
      REQUIRE(right.Offset(0,1,2,3) == 1219);

      MDVector<4,MDLayoutRight<4>> right4(NA, NB, NC, ND);
      right4.SetLayout(MDLayoutRight<4>({3,2,1,0}));
      REQUIRE(right4.Offset(0,1,2,3) == 1219);
   }
}

TEST_CASE("MDGridFunction layouts", "[MDSpan], [MDGridFunction]")
{
   constexpr int NE = 18, NG = 2, NA = 32;

   const bool all = launch_all_non_regression_tests;
   auto p = all ? GENERATE(1,2) : 2;
   auto nx = all ? GENERATE(3,5) : 2;
   auto dim = all ? GENERATE(1,2,3) : 2;

   auto MakeCartesian = [](int dim, int nx)
   {
      return dim == 2 ? Mesh::MakeCartesian2D(nx, nx, Element::QUADRILATERAL):
             dim == 3 ? Mesh::MakeCartesian3D(nx, nx, nx, Element::HEXAHEDRON):
             Mesh::MakeCartesian1D(nx);
   };

   Mesh mesh = MakeCartesian(dim, nx);

   H1_FECollection fec(p, dim);
   FiniteElementSpace fes(&mesh, &fec);
   const int ND = fes.GetNDofs();

   SECTION("Types")
   {
      MDGridFunction<4> mdgf(NE, NG, &fes, NA);
      REQUIRE(mdgf.Size() == (NE * NG * fes.GetVSize() * NA));
      REQUIRE(std::is_same<decltype(mdgf.HostRead()), double const*>());
      REQUIRE(std::is_same<decltype(mdgf.MDHostRead()), mfem::MDTensor<4, double const> const>());
   }

   SECTION("EGDA, left")
   {
      MDGridFunction<4> gsa(NE, NG, &fes, NA);
      const int gsa_0123 = gsa.Offset(0, 1, 2, 3);
      REQUIRE(gsa_0123 == 0 + 1*(NE) + 2*(NE*NG) + 3*(NE*NG*ND));
      mfem::Array<int> vdofs(ND);
      gsa.GetVDofs(0, 1, vdofs, 3);
      REQUIRE(0 <= vdofs.Min());
      REQUIRE(vdofs[2] == gsa_0123);
      REQUIRE(vdofs.Max() < NE*NG*NA*ND);

   }

   SECTION("EGDA, right")
   {
      MDGridFunction<4, MDLayoutRight<4>> gsa(NE, NG, &fes, NA);
      const int gsa_0123 = gsa.Offset(0,1,2,3);
      REQUIRE(gsa_0123 == 0*(NG*ND*NA) + 1*(ND*NA) + 2*(NA) + 3);
      mfem::Array<int> vdofs(ND);
      gsa.GetVDofs(0, 1, vdofs, 3);
      REQUIRE(0 <= vdofs.Min());
      REQUIRE(vdofs[2] == gsa_0123);
      REQUIRE(vdofs.Max() < NE*NG*NA*ND);
   }
}

TEST_CASE("MDGridFunction reshapes", "[MDSpan], [MDGridFunction]")
{
   auto is_equal = [&](const Vector &a, const Vector &b)
   {
      REQUIRE(a.Size() == b.Size());
      for (int i = 0; i < a.Size(); i++)
      {
         const double va = a.GetData()[i], vb = b.GetData()[i];
         REQUIRE(va == MFEM_Approx(vb));
      };
      return true;
   };

   SECTION("MDReshapes")
   {
      constexpr int p = 2;
      constexpr int dim = 3;
      constexpr int nx = 5, ny = 3, nz = 2;
      Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, mfem::Element::HEXAHEDRON);

      H1_FECollection fec_mesh(p, dim);
      FiniteElementSpace fes_mesh(&mesh, &fec_mesh, dim);
      mesh.SetNodalFESpace(&fes_mesh);

      L2_FECollection fec(p, dim);
      FiniteElementSpace fes(&mesh, &fec);

      const std::list<MDLayout<3>> layouts =
      { {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,1,0}, {2,0,1} };

      for (auto &layout: layouts)
      {
         constexpr int numGroups = 4, numAngles = 7;

         MDGridFunction<3> psi(&fes, numGroups, numAngles);
         psi.SetLayout(layout);

         const GridFunction *nodes = mesh.GetNodes();
         const FiniteElementSpace *mfes = mesh.GetNodalFESpace();
         const int ng = numGroups, na = numAngles, ne = mfes->GetNE();
         const ElementDofOrdering e_ordering = ElementDofOrdering::LEXICOGRAPHIC;
         const Operator *R = mfes->GetElementRestriction(e_ordering);
         REQUIRE(R);
         const FiniteElement *mfe = mfes->GetFE(0);
         const int nd = mfe->GetDof(), vdim = mfes->GetVDim();
         Vector nodes_e(vdim*nd*ne); nodes_e.UseDevice(true);
         constexpr int D1D = p + 1;
         REQUIRE(fes.GetVSize() == D1D*D1D*D1D*ne);
         nodes_e.Read();
         REQUIRE(nodes);
         R->Mult(*nodes, nodes_e);
         const auto X = mfem::Reshape(nodes_e.Read(), D1D, D1D, D1D, vdim, ne);
         auto dY = psi.MDWrite();

         MDGridFunction<3> rY1(&fes, numGroups, numAngles);
         rY1.SetLayout(MDLayout<3>(layout));
         REQUIRE(rY1.Size() == psi.Size());
         auto drY1 = rY1.MDWrite();

         MDGridFunction<3> rY2(&fes, numGroups, numAngles);
         rY2.SetLayout(MDLayout<3>(layout));
         REQUIRE(ng%2 == 0);
         auto drY2 = rY2.MDReshape<4>(rY2.Write(),
                                      D1D*D1D*D1D*ne,
                                      std::array<int,2> {2, ng/2},
                                      na);

         MDGridFunction<3> rY3(&fes, numGroups, numAngles);
         rY3.SetLayout(MDLayout<3>(layout));
         auto drY3 = rY3.MDReshape<6>(rY3.Write(),
                                      std::array<int,4> {D1D, D1D, D1D, ne},
                                      ng, na);

         MDGridFunction<3> rY4(&fes, numGroups, numAngles);
         rY4.SetLayout(MDLayout<3>(layout));
         auto drY4 = rY4.MDReshape<7>(rY4.Write(),
                                      std::array<int,4> {D1D, D1D, D1D, ne},
                                      std::array<int,2> {1, ng},
                                      na);

         MDGridFunction<3> rY5(&fes, numGroups, numAngles);
         rY5.SetLayout(MDLayout<3>(layout));
         auto drY5 = rY5.MDReshape<7>(rY5.Write(),
                                      std::array<int,4> {D1D, D1D, D1D, ne},
                                      std::array<int,2> {2, ng/2},
                                      na);

         const double exp_m08 = exp(-0.8);

         MFEM_FORALL_3D(ega, ne*ng*na, D1D,D1D,D1D,
         {
            const int e = ega/(ng*na), ga = ega%(ng*na), g = ga/na, a = ga%na;
            MFEM_FOREACH_THREAD(dz,z,D1D)
            {
               MFEM_FOREACH_THREAD(dy,y,D1D)
               {
                  MFEM_FOREACH_THREAD(dx,x,D1D)
                  {
                     const int xyze = dx + D1D*(dy + D1D*(dz + D1D*(e)));

                     const double p0 = X(dx,dy,dz,0,e), p1 = X(dx,dy,dz,1,e);
                     const double value = 1.0 - exp_m08*cos(M_PI*p0)*cos(M_PI*p1);

                     dY(xyze,g,a) = value;
                     drY1(xyze,g,a) = value;
                     drY2(xyze,g%2,g/2,a) = value;
                     drY3(dx,dy,dz,e, g, a) = value;
                     drY4(dx,dy,dz,e, 0,g, a) = value;
                     drY5(dx,dy,dz,e, g%2,g/2, a) = value;
                  }
               }
            }
         });
         psi.MDHostRead(); rY1.HostRead();
         REQUIRE(is_equal((Vector&)rY1, (Vector&)psi));
         REQUIRE(is_equal((Vector&)rY2, (Vector&)psi));
         REQUIRE(is_equal((Vector&)rY3, (Vector&)psi));
         REQUIRE(is_equal((Vector&)rY4, (Vector&)psi));
         REQUIRE(is_equal((Vector&)rY5, (Vector&)psi));
      }
   }
}
