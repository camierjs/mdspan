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

#ifndef MFEM_MDGRIDFUNC
#define MFEM_MDGRIDFUNC

#include "../config/config.hpp"

#include "fem/gridfunc.hpp"
#include "general/mdspan.hpp"

namespace mfem
{

template<int N, class Layout = MDLayoutLeft<N>>
class MDGridFunction : public MDSpan<GridFunction, N, Layout>
{
   using base_t = MDSpan<GridFunction, N, Layout>;
   using base_t::Nd;
   using base_t::Sd;

public:

   using GridFunction::data;

   using GridFunction::Read;
   using GridFunction::Write;
   using GridFunction::HostRead;
   using GridFunction::HostWrite;
   using GridFunction::operator=;

   MDGridFunction(): base_t() { }

   MDGridFunction(MDGridFunction&&) = delete;
   MDGridFunction(const MDGridFunction&) = delete;
   MDGridFunction& operator=(MDGridFunction&&) = delete;

   template <typename... Ts>
   MDGridFunction(FiniteElementSpace *fes, Ts... args): MDGridFunction(args...)
   {
      this->GridFunction::fes = fes;
      base_t::Setup(fes->GetNDofs(), args...);
   }

   template <typename... Ts>
   MDGridFunction(int n, Ts... args): MDGridFunction(args...)
   {
      base_t::Setup(n, args...);
   }

   template <int n = 1, typename... Ts>
   void GetScalarGridFunction(GridFunction &gf, Ts... args) const
   {
      FiniteElementSpace *fes = GridFunction::fes;
      MFEM_VERIFY(fes->GetNDofs() == Nd[n-1], "Error in dofs size!");
      gf.SetSpace(fes);
      for (int s = 0; s < Nd[n-1]; s++)
      {
         gf[s] = data[get_vdofs_offset +
                      MDOffset<n,N,int,Ts...>::offset(Sd, s, args...)];
      }
      get_vdofs_offset = 0; // re-init for next calls
   }

   template <int n = 1, typename... Ts>
   void GetScalarGridFunction(int dim, Ts&&... args) const
   {
      get_vdofs_offset += dim * Sd[n-1];
      MDGridFunction::GetScalarGridFunction<n+1>(std::forward<Ts>(args)...);
   }

   template <int n = 1, typename... Ts>
   void SetScalarGridFunction(const GridFunction &gf, Ts... args)
   {
      MFEM_VERIFY(GridFunction::fes->GetNDofs() == Nd[n-1], "Error in fespace size!");
      for (int s = 0; s < Nd[n-1]; s++)
      {
         data[get_vdofs_offset +
              MDOffset<n,N,int,Ts...>::offset(Sd, s, args...)] = gf[s];
      }
      get_vdofs_offset = 0; // re-init for next calls
   }

   template <int n = 1, typename... Ts>
   void SetScalarGridFunction(int dim, Ts... args)
   {
      get_vdofs_offset += dim * Sd[n-1];
      MDGridFunction::SetScalarGridFunction<n+1>(args...);
   }

private:
   mutable int get_vdofs_offset = 0;
};

} // namespace mfem

#endif // MFEM_MDGRIDFUNC
